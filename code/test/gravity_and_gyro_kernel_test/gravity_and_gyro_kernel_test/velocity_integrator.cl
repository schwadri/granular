typedef unsigned int index_t;

#define BUILD_FOR_NV_GT330M

//hardware specific constants and utilities

inline index_t bank_line_float(index_t i);

//Nvidia GT330M specific sizes
#ifdef BUILD_FOR_NV_GT330M
  #define WORK_GROUP_SIZE               512
  #define HALF_WORK_GROUP_SIZE          (WORK_GROUP_SIZE / 2)
  #define QUARTER_WORK_GROUP_SIZE       (WORK_GROUP_SIZE / 4)
  #define QUARTER_HALF_WORK_GROUP_SIZE  (HALF_WORK_GROUP_SIZE / 4)
  #define LOCAL_MEM_BANK_COUNT    16
  //how many consecutive bytes belong to the same local memory bank
  #define LOCAL_MEM_BANK_WIDTH    4
  #define LOCAL_MEM_BANK_SIZE     1024

/*this method returns for a given offset the local memory bank line
  where the offset points to. the resulting value can be used to introduce
  index shifts to avoid memory bank conflicts*/
inline index_t bank_line_float(index_t i) {
  //i / 16
  return i >> 4;
}
inline index_t bank_line_double(index_t i) {
  //i / 8
  return i >> 3;
}
#endif

#ifdef BUILD_FOR_NV_FERMI
  #define WORK_GROUP_SIZE               1536
  #define HALF_WORK_GROUP_SIZE          (WORK_GROUP_SIZE / 2)
  #define QUARTER_WORK_GROUP_SIZE       (WORK_GROUP_SIZE / 4)
  #define QUARTER_HALF_WORK_GROUP_SIZE  (HALF_WORK_GROUP_SIZE / 4)
  #define LOCAL_MEM_BANK_WIDTH -1
  #define LOCAL_MEM_BANK_COUNT 32
  #define LOCAL_MEM_BANK_SIZE  2048

inline index_t bank_line_float(index_t i) {
  return i >> 5;
}
#endif




//rigid body problem specific constants
//rigid body with quaternion specific state sizes
#define BODY_STATE_SIZE 8
#define BODY_INERTIA_SIZE 4
#define BODY_Q_SIZE 8
#define BODY_U_SIZE 8
#define BODY_X_SIZE 4
#define BODY_P_SIZE 4
#define BODY_DX_SIZE 4
#define BODY_DP_SIZE 4
#define BODY_Q_U_SIZE 16
#define HALF_BODY_STATE_SIZE (BODY_STATE_SIZE / 2)

/* the opencl cross function performs the cross
   product between the first three components of
   the input vectors, but we need it for the last
   three
 */
inline float4 back_cross(float4 a, float4 b);

inline float4 back_cross(float4 a, float4 b) {
  return (float4)(
    (float)0.0,
    a.s2 * b.s3 - a.s3 * b.s2,
    a.s3 * b.s1 - a.s1 * b.s3,
    a.s1 * b.s2 - a.s2 * b.s1
  );
}

/* velocity integration on float4 vectors, better readability but
   not optimal for the nvidia gpu architecture
 */
__kernel void integrate_velocities_v4(
  __global float4 * g_dx, __global float4 * g_dp,
  __global float4 * g_inertia, __global float4 * g_inv_inertia,
  float dt,
  index_t nbodies
) {
  //get global thread id
  index_t g_id    = get_global_id(0);
  index_t body_id = g_id;

  //if the global thread id is bigger than the number of bodies
  if(body_id >= nbodies) return;
  
  float4 theta      = g_inertia[body_id];
  float4 inv_theta  = g_inv_inertia[body_id];
  
  //step 1: integrate translational velocities
  float4 vx_B = g_dx[body_id];
  //for now only apply gravity
  float4 vx_E = vx_B + dt * (float4)(0.0f, 0.0f, -9.81f, 0.0f);
  g_dx[g_id] = vx_E;
  
  //step 2: integrate angular velocities
  float4 omega_B = g_dp[body_id];
  
  float4 gyro    = inv_theta * back_cross(omega_B, theta * omega_B);
  //multiply gyroscopic term with inverse inertia matrix
  float4 omega_E = omega_B - dt * gyro;
  g_dp[body_id]  = omega_E;
}

/* integrate_velocities performs the same operations as
   integrate_velocities_v4 but it is better optimized for the nvidia arch.
 */
#undef BANK_SHIFT_SIZE
#define BANK_SHIFT_SIZE (WORK_GROUP_SIZE / LOCAL_MEM_BANK_COUNT)
__kernel
void integrate_velocities(
  __global float * g_dx, __global float * g_dp,
  __global float * g_inertia, __global float * g_inv_inertia,
  float dt,
  index_t nbodies
) {
  //step 0. get global and local thread ids
  index_t g_id  = get_global_id(0);
  index_t l_id  = get_local_id(0);
  
  //if the global thread id is bigger than the number of bodies times their
  //half state size kill the thread.
  if(g_id >= nbodies * BODY_DX_SIZE) return;
  
  //step 1: integrate translational velocities
  
  {
    float fexti;
    index_t dx_row_id = l_id & 0x3; //l_id mod 4
    
    //FIXME: thread divergence fixing? not worthwhile here as it destroys
    //coalesced global memory access
    if(dx_row_id == 2)
      fexti = -9.81f;
    else
      fexti = 0.0f;
    g_dx[g_id] = g_dx[g_id] + dt * fexti;
  }
  
  //step 2: integrate angular velocities
  {
    //load coalesced
    //NOTE: check, that there is enough local memory for these
    __local float l_inertia[WORK_GROUP_SIZE + BANK_SHIFT_SIZE];
    __local float l_inv_inertia[WORK_GROUP_SIZE + BANK_SHIFT_SIZE];
    __local float l_omega[WORK_GROUP_SIZE + BANK_SHIFT_SIZE];
    
    //coalesced load of inertia and omega values
    float a = g_inertia[g_id];
    float b = g_inv_inertia[g_id];
    float c = g_dp[g_id];
    
    //determine shift to counter bank conflicts
    index_t l_shift       = bank_line_float(l_id);
    index_t shifted_l_id  = l_id + l_shift;

    //write to shifted local memory location
    l_inertia[shifted_l_id] = a;
    l_inv_inertia[shifted_l_id] = b;
    l_omega[shifted_l_id] = c;
    
    //reorganize threads into groups which work on the same rows of omega but for different bodies
    index_t dp_row_id           = l_id >> 7;                            // l_id / QUARTER_WORK_GROUP_SIZE
    index_t body_id             = l_id & (QUARTER_WORK_GROUP_SIZE - 1); // l_id mod QUARTER_WORK_GROUP_SIZE
    index_t body_offset         = body_id << 2;                         // body_id * 4
    index_t body_l_shift        = bank_line_float(body_offset);
    index_t shifted_body_offset = body_offset + body_l_shift;
    index_t shifted_row_offset  = shifted_body_offset + dp_row_id;
    
    //make sure that all writes to local memory are visible by all threads in the same work group 
    mem_fence(CLK_LOCAL_MEM_FENCE);
    
    //load inertia and omega for this thread
    float4 inv_inertia  = *(__local float4 *)(l_inv_inertia + shifted_body_offset);
    float4 omega_B      = *(__local float4 *)(l_omega       + shifted_body_offset);
    float4 inertia      = *(__local float4 *)(l_inertia     + shifted_body_offset);
    
    //determine component of omega_E
    float omega_E_i;
    switch(dp_row_id) {
      case 0:
        omega_E_i = omega_B[0];
        break;
      case 1:
        omega_E_i = omega_B[1] - dt * inv_inertia[1] * (inertia[3] - inertia[2]) * omega_B[2] * omega_B[3];
        break;
      case 2:
        omega_E_i = omega_B[2] - dt * inv_inertia[2] * (inertia[1] - inertia[3]) * omega_B[3] * omega_B[1];
        break;
      case 3:
        omega_E_i = omega_B[3] - dt * inv_inertia[3] * (inertia[2] - inertia[1]) * omega_B[1] * omega_B[2];
        break;
    }
    //write back to local memory
    l_omega[shifted_row_offset] = omega_E_i;
    
    mem_fence(CLK_LOCAL_MEM_FENCE);
    //FIXME: might be improved
    //coalesced write back new omegas to global
    g_dp[g_id] = l_omega[shifted_l_id];
  }
}


