typedef unsigned int index_t;
typedef unsigned int bool_t;
typedef float real_t;
typedef float2 real2_t;
typedef float4 real4_t;
typedef float8 real8_t;


#define BUILD_FOR_NV_GT330M

//hardware specific constants and utilities

inline index_t bank_line_float(index_t i);

//Nvidia GT330M specific sizes
#ifdef BUILD_FOR_NV_GT330M
  #define BLOCK_SIZE              512
  #define WORK_GROUP_SIZE         512
  #define HALF_BLOCK_SIZE         (BLOCK_SIZE / 2)
  #define QUARTER_BLOCK_SIZE      (BLOCK_SIZE / 4)
  #define QUARTER_HALF_BLOCK_SIZE (HALF_BLOCK_SIZE / 4)
  #define LOCAL_MEM_BANK_COUNT    16
  //how many consecutive bytes belong to the same local memory bank
  #define LOCAL_MEM_BANK_WIDTH    4
  #define LOCAL_MEM_BANK_SIZE     1024
  //FIXME:problem specific constant. move elsewhere
  #define SHIFT_SIZE              32

/*this method returns for a given offset the local memory bank line
  where the offset points to. the resulting value can be used to introduce
  index shifts to avoid memory bank conflicts*/
inline index_t bank_line_float(index_t i) {
  return i >> 4;
}
#endif

#ifdef BUILD_FOR_NV_FERMI
  #define BLOCK_SIZE           1536
  #define WORK_GROUP_SIZE      1536
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
#define BODY_Q_SIZE 8
#define BODY_U_SIZE 8
#define BODY_X_SIZE 4
#define BODY_P_SIZE 4
#define BODY_DX_SIZE 4
#define BODY_DP_SIZE 4
#define BODY_Q_U_SIZE 16
#define HALF_BODY_STATE_SIZE (BODY_STATE_SIZE / 2)

/* NOTE: -on nvidia geforce gt 330m the euler kernel is much faster than euler_v4 or euler_v8.
          it seems that the vector expressions are handled serially as the nv architecture is in essence a scalar
          architecture, wherease amd gpu's have real simd units.
 
  experiences so far: The aim was to implement the euler integration step of the generalized coordinates for
  a rigid body with quaternions on the gpu and optimize it to gain first experiences with the constraints of this
  computing model and the opencl framework. All experiences were made on the Nvidia the GT 330M, which supports only
  single precision floating point.
  Most important Specs:
  The GT 330M has 6 SM's, each of which has 8 SP's, 16kb local memory divided in 16 memory banks à 1024bytes.
  The warp size is 32 and the maximum number of warps allowed on a SM is 16 e.g. 512 threads per SP which is also the
  maximum work group size.
 
  I started out with a simple vector addition and extended this to the euler_v8 kernel,
  because the rigid body with quaternion can be represented by 7 coordinates. Using the vector builtin vector primitives of
  the opencl C kernel programing language and implementing the quaternion update with it was straightforward. This lead to the
  kernel euler_v8_w_quat. But this kernel is not optimal for Nvidia cards because they handle vector expressions sequentially,
  whereas amd cards have simd/vector lanes that process them in parallel. To figure out how bad my implementation is, i
  wrote the euler kernel and compared it to euler_v8. euler_v8 ran ~5 times slower. So in the euler_v8_quat kernel every thread
  handles one single rigid body and this is not optimal. To improve on that kernel euli was written. Here every thread handles
  a single coordinate update of a rigid body. To minimize global memory transfers the coordinates are loaded into local memory
  at the beginning of the kernel. Then the threads are divided into two groups to avoid thread divergence problems. The first
  group is performing the position update and the second group the quaternion update. Now I was facing several new problems.
  The threads in the second group suffer from bank conflicts because the rigid body state sizes are properly aligned with the
  memory bank count and threads in the same warp are updating different quaternion entries which have to be calculated differently
  so thread divergence is still present. Fixing these resulted in euli2, 3 and euli7. euli7 performs quite reasonably. 
  To fix the thread divergence problem two different approaches were taken. The first which is implemented in euli2, 3 and 7 is to split the 
  quaternion handling thread into four groups where each group of consecutive threads works on the same quaternion component but for different bodies.
  The second was to find a formulation where each quaternion component can be calculated the same way, by using index magic and lookup tables. 
  This was implemented in euli4,5 and 6.
  The second approach turned out to perform very poorly, e.g. ~10 times slower than the other as too many index operations are necessary to
  hide the differences of the different quaternion components. Maybe there might be a better way to do this, but it seems not worthwhile
  to continue on that road.

  Finally a different approach was chosen which resulted in the euler_split* kernels. As the position and orientation updates are structurally 
  very different it should be better to perform them in separate kernels or sequentially. Combining the approach of split state arrays with
  the quaternion integration part of euli7 resulted in euler_split4 which is the fastest kernel so far, e.g. ~80% speed of the euler kernel
  which does a simple mad operation. 

  TODO: -properly handle the case very the number of bodies is not properly aligned with the body_per_work_group size
        -quaternion normalization!

  -implementation:
    - efficient synchronization only possible within a work group 
      - single work item is load store consistent (thank god)
        mem_fence
        barrier ? when to use it instead of mem_fence?
        read_mem_fence ? what is it good for?
        write_mem_fence

  -optimization:
    -memory
      minimize global memory transfers
       -> is it better to have separate arrays for particle states,
          or should a single particle have all of its state in a contiguous
          memory block?
      coalescing global memory access
      avoid local memory bank conflicts

 
    -control
      avoid thread divergence
      register blocking
      register dependency latency
 
    -synchronization
      within a warp (nv specific) synchronization is not necessary as it runs in lock step and
      behaves like a single work item which is load store consistent
 */

__kernel void euler_v8(
  __global float8 * s, __global float8 * ds,
  real_t dt,
  index_t size
) {
  index_t i = get_global_id(0);

  if(i >= size) return;

  float8 x = s[i];
  float8 v = ds[i];

  x = x + dt * v;

  s[i]  = x;
}

__kernel void euler_v4(
  __global float4 * s, __global float4 * ds,
  real_t dt,
  index_t size
) {
  index_t i = get_global_id(0);

  if(i >= size * 2) return;

  float4 x = s[i];
  float4 v = ds[i];

  x = x + dt * v;

  s[i]  = x;
}

__kernel void euler(
  __global float * s, __global float * ds,
  real_t dt,
  index_t size
) {
  index_t i = get_global_id(0);

  if(i >= size * 8) return;

  float x = s[i];
  float v = ds[i];

  x = x + dt * v;

  s[i]  = x;
}



//define prototypes
float4 quat_mul(float4 a, float4 b);
float4 omega_to_dp(float4 p, float4 o);

/*multiply two quaternions*/
float4 quat_mul(float4 a, float4 b) {
  float4 r;
  r[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
  r[1] = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2];            
  r[2] = a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1];
  r[3] = a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0];
  return r;
}

/*map angular velocity \a o represented in the body fixed frame
 to a quaternion differential*/
float4 omega_to_dp(float4 p, float4 o) {
  float4 dp;
  dp[0] = (real_t)(0.5) * (p[0] * o[0] - p[1] * o[0] - p[2] * o[1] - p[3] * o[2]);
  dp[1] = (real_t)(0.5) * (p[1] * o[0] + p[0] * o[0] - p[3] * o[1] + p[2] * o[2]);
  dp[2] = (real_t)(0.5) * (p[2] * o[0] + p[3] * o[0] + p[0] * o[1] - p[1] * o[2]);
  dp[3] = (real_t)(0.5) * (p[3] * o[0] - p[2] * o[0] + p[1] * o[1] + p[0] * o[2]);
  return dp;
}

/* $\dot p = F(q) u$
 */
__kernel void euler_v8_w_quat(
  __global float8 * s, __global float8 * ds,
  real_t dt,
  index_t size
) {
  //get global work item / thread id
  index_t i = get_global_id(0);
  //if the id is greater than the number of bodies, kill this thread
  if(i >= size) return;
  
  //load generalized coordinates from global memory
  float8 x      = s[i];
  //load generalized velocities from global memory
  float8 u      = ds[i];
  //extract the angular velocity from the generalized velocity vector
  float4 omega  = u.s4567;
  //extract the quaternion from the generalized coordinates vector
  float4 p      = x.s4567;
  //determine the quaternion time differential from the angular velocity
  //and the current quaternion
  float4 dp     = omega_to_dp(p, omega);
  //combine the translational velocity with the quaternion time differential
  float8 dx_dt  = (float8)(u.s0123, dp);
  //do euler integration step
  x = x + dt * dx_dt;
  
  //write new generalized coordinates back to global memory
  s[i]  = x;
}



/* FIXME: this is not finished
 the idea is to load the rigid body phase state from an interleaved
 array, but i don't think that this has any advantages, because
 we only update the positions and leave the velocities unchanged
 -> break coalescing during writeback to global memory
 */
__kernel void euler_single_state(
  __global float * g_s,
  real_t dt,
  index_t nbodies
) {
  __local float l_s[2 * BLOCK_SIZE];
  index_t id  = get_global_id(0);
  index_t gsize  = get_global_size(0);
  index_t lid = get_local_id(0);
  //if the id is greater than the number of bodies times their state size, kill this thread
  if(id >= nbodies * BODY_Q_U_SIZE) return;
  
  //step 1
  //copy first half of bodies into local memory BLOCK_SIZE / (BODY_Q_U_SIZE)
  l_s[lid]  = g_s[id];
  mem_fence(CLK_LOCAL_MEM_FENCE);
  
  //step 2
  //copy second half of bodies into local memory
  l_s[lid + BLOCK_SIZE]  = g_s[id + gsize];
  
  
  //step 3 do something with the data locally
  
  //step 4 write back to global
  
  //first half
  g_s[id]  = l_s[lid];
  //second half
  //g_s[id + gsize]  = l_s[lid + BLOCK_SIZE];
}

/** combined variants: here we use 8 threads per body and perform the integration of the position and orientation
    in parallel.
    NOTE: it turned out not to be worthwhile to do the integration of the position and orientation at the same time
          as they are structurally very different.
 */


/*
  - first step to fix thread divergence. group threads into two groups. 
    the first group handles the position update and the second group
    computes the quaternion update, which is more involved
 */
__kernel void euli(
  __global float * g_s, __global float * g_ds,
  real_t dt,
  index_t nbodies
) {
  __local float l_s[BLOCK_SIZE];
  __local float l_ds[BLOCK_SIZE];
  index_t id  = get_global_id(0);
  index_t lid = get_local_id(0);
  //if the id is greater than the number of bodies time their state size, kill this thread
  if(id >= nbodies * BODY_STATE_SIZE) return;
  
  //step 1
  //copy global body state for BLOCK_SIZE / BODY_STATE_SIZE bodies into local memory
  //all threads read a linear piece of memory to enable coalesced memory access
  l_s[lid]  = g_s[id];
  l_ds[lid] = g_ds[id];
  //make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);
  
  //step 2
  //reorganize threads
  if(lid < (BLOCK_SIZE / 2)) {
    //the first half performs the translational velocity integration
    //per position update four threads are needed per body
    //the next four threads have to be shifted to the next multipled of BODY_STATE_SIZE
    // 8 * ceil(lid/4) + mod(lid, 4)
    index_t stride =  ((lid & 0xfffffffc) << 1) | (lid & 0x03);
    //index_t stride = BODY_STATE_SIZE * (lid / 4) + (lid % 4);
    //perform the actual operation
    l_s[stride] = l_s[stride] + dt * l_ds[stride];
  } else {
    //the other half has to do the quaternion update
    //calculate a local half block thread index which runs from 0 to BLOCK_SIZE / 2 - 1
    index_t s_lid = lid - (BLOCK_SIZE / 2);
    //this half block is further grouped into groups of four threads per body
    //determine this threads local index in a 4 group
    //index_t s_dpid = (s_lid % 4);
    index_t s_dpid = (s_lid & 0x03);
    //
    //index_t stride = BODY_STATE_SIZE * (s_lid / 4) + 4;
    index_t stride = ((s_lid & 0xfffffffc) << 1) + 4;
    index_t row_stride = stride + s_dpid;
    
    //load the angular velocity and quaternion from local memory
    float4 p = *(__local float4 *)(l_s + stride);
    float4 o = *(__local float4 *)(l_ds + stride);
    
    //calculate a component of the quaternion time differential where the 4-group local
    //thread id determines which component to calculate
    float dpi;
    switch(s_dpid) {
      case 0:
        dpi = (real_t)(0.5) * (p[0] * o[0] - p[1] * o[1] - p[2] * o[2] - p[3] * o[3]);
        break;
      case 1:
        dpi = (real_t)(0.5) * (p[1] * o[0] + p[0] * o[1] - p[3] * o[2] + p[2] * o[3]);
        break;
      case 2:
        dpi = (real_t)(0.5) * (p[2] * o[0] + p[3] * o[1] + p[0] * o[2] - p[1] * o[3]);
        break;
      case 3:
        dpi = (real_t)(0.5) * (p[3] * o[0] - p[2] * o[1] + p[1] * o[2] + p[0] * o[3]);
        break;
    }
    //integrate the quaternion component and write back to local memory
    l_s[row_stride] = l_s[row_stride] + dt * dpi;
  }
  //make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);
  
  //step 3
  //write back to global
  //again consecutive threads write consecutive data for coalesced access
  g_s[id]  = l_s[lid];
  //g_ds[id] = l_ds[id];
}

/*
  - fixed thread divergence by grouping threads into warps that handle
    the same row of the body state update.
  - suffers from bank conflicts
 */
__kernel void euli2(
  __global float * g_s, __global float * g_ds,
  real_t dt,
  index_t nbodies
) {
  __local float l_s[BLOCK_SIZE];
  __local float l_ds[BLOCK_SIZE];
  index_t id  = get_global_id(0);
  index_t lid = get_local_id(0);
  //if the id is greater than the number of bodies time their state size, kill this thread
  if(id >= nbodies * BODY_STATE_SIZE) return;
  
  //step 1
  //copy global body state for BLOCK_SIZE / BODY_STATE_SIZE bodies into local memory
  //all threads read a linear piece of memory to enable coalesced memory access
  l_s[lid]  = g_s[id];
  l_ds[lid] = g_ds[id];
  //make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);
  
  //step 2
  //reorganize threads
  if(lid < HALF_BLOCK_SIZE) {
    //the first half performs the translational velocity integration
    //per position update four threads are needed per body
    //the next four threads have to be shifted to the next multipled of BODY_STATE_SIZE
    // 8 * ceil(lid/4) + mod(lid, 4)
    index_t stride =  ((lid & 0xfffffffc) << 1) | (lid & 0x03);
    //index_t stride = BODY_STATE_SIZE * (lid / 4) + (lid % 4);
    //perform the actual operation
    l_s[stride] = l_s[stride] + dt * l_ds[stride];
  } else {
    //the other half has to do the quaternion update
    //calculate a local half block thread index which runs from 0 to BLOCK_SIZE / 2 - 1
    index_t s_lid = lid - HALF_BLOCK_SIZE;
    //this half block is further grouped into groups of four threads per body
    //determine this threads local index in a 4 group
    //index_t s_dpid = s_lid / QUARTER_HALF_BLOCK_SIZE; s_lid / 64 = s_lid / (2^6)
    index_t s_dpid = s_lid >> 6;
    
    //index_t body_id = s_lid % QUARTER_HALF_BLOCK_SIZE;
    index_t body_id = s_lid & (QUARTER_HALF_BLOCK_SIZE - 1);
    index_t stride = BODY_STATE_SIZE * body_id + HALF_BODY_STATE_SIZE;
    index_t row_stride = stride + s_dpid;
    
    //load the angular velocity and quaternion from local memory
    float4 p = *(__local float4 *)(l_s + stride);
    float4 o = *(__local float4 *)(l_ds + stride);
    
    //calculate a component of the quaternion time differential where the 4-group local
    //thread id determines which component to calculate
    float dpi;
    switch(s_dpid) {
      case 0:
        dpi = (real_t)(0.5f) * (p[0] * o[0] - p[1] * o[1] - p[2] * o[2] - p[3] * o[3]);
        break;
      case 1:
        dpi = (real_t)(0.5f) * (p[1] * o[0] + p[0] * o[1] - p[3] * o[2] + p[2] * o[3]);
        break;
      case 2:
        dpi = (real_t)(0.5f) * (p[2] * o[0] + p[3] * o[1] + p[0] * o[2] - p[1] * o[3]);
        break;
      case 3:
        dpi = (real_t)(0.5f) * (p[3] * o[0] - p[2] * o[1] + p[1] * o[2] + p[0] * o[3]);
        break;
    }
    //integrate the quaternion component and write back to local memory
    l_s[row_stride] = l_s[row_stride] + dt * dpi;
  }
  //make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);
  
  //step 3
  //write back to global
  //again consecutive threads write consecutive data for coalesced access
  g_s[id]  = l_s[lid];
}


/*
  - fixed thread divergence
  - fixed bank conflicts
 */
__kernel void euli3(
  __global float * g_s, __global float * g_ds,
  real_t dt,
  index_t nbodies
) {
  __local float l_s[BLOCK_SIZE + SHIFT_SIZE];
  __local float l_ds[BLOCK_SIZE + SHIFT_SIZE];
  index_t g_id  = get_global_id(0);
  index_t l_id = get_local_id(0);
  //if the id is greater than the number of bodies time their state size, kill this thread
  if(g_id >= nbodies * BODY_STATE_SIZE) return;
  
  //step 1
  //copy global body state for BLOCK_SIZE / BODY_STATE_SIZE bodies into local memory
  //all threads read a linear piece of memory to enable coalesced memory access
  //and introduce a shift
  //index_t shift = lid / GT330M_BANK_COUNT;
  index_t shift = bank_line_float(l_id);
  l_s[l_id  + shift]  = g_s[g_id];
  l_ds[l_id + shift]  = g_ds[g_id];
  
  
  //step 2
  //reorganize threads
  if(l_id < HALF_BLOCK_SIZE) {
    //the first half performs the translational velocity integration
    //per position update four threads are needed per body
    //the next four threads have to be shifted to the next multipled of BODY_STATE_SIZE
    // 8 * ceil(lid/4) + mod(lid, 4)
    index_t body_stride = ((l_id & 0xfffffffc) << 1);
    //index_t l_shift = body_stride / GT330M_BANK_COUNT;
    index_t l_shift = body_stride >> 4;
    index_t x_id  = (l_id & 0x03);
    index_t stride =  body_stride + x_id + l_shift;
    //index_t stride = BODY_STATE_SIZE * (lid / 4) + (lid % 4);
    
    //make sure that all the data in the local memory is ready
    mem_fence(CLK_LOCAL_MEM_FENCE);
    //perform the actual operation
    l_s[stride] = l_s[stride] + dt * l_ds[stride];
  } else {
    //the other half has to do the quaternion update
    //calculate a local half block thread index which runs from 0 to BLOCK_SIZE / 2 - 1
    index_t s_lid = l_id - HALF_BLOCK_SIZE;
    //this half block is further grouped into groups of four threads per body
    //determine this threads local index in a 4 group
    //index_t s_dpid = s_lid / QUARTER_HALF_BLOCK_SIZE; s_lid / 64 = s_lid / (2^6)
    index_t s_dpid = s_lid >> 6;
    
    //index_t body_id = s_lid % QUARTER_HALF_BLOCK_SIZE;
    index_t body_id     = s_lid & (QUARTER_HALF_BLOCK_SIZE - 1);
    index_t body_stride = BODY_STATE_SIZE * body_id;
    index_t l_shift = body_stride >> 4;
    index_t stride = body_stride + HALF_BODY_STATE_SIZE + l_shift;
    index_t row_stride = stride + s_dpid;
    
    //make sure that all the data in the local memory is ready
    mem_fence(CLK_LOCAL_MEM_FENCE);
    
    //load the angular velocity and quaternion from local memory
    float4 p = *(__local float4 *)(l_s + stride);
    float4 o = *(__local float4 *)(l_ds + stride);
    
    //calculate a component of the quaternion time differential where the 4-group local
    //thread id determines which component to calculate
    float dpi;
    switch(s_dpid) {
      case 0:
        dpi = (real_t)(0.5f) * (p[0] * o[0] - p[1] * o[1] - p[2] * o[2] - p[3] * o[3]);
        break;
      case 1:
        dpi = (real_t)(0.5f) * (p[1] * o[0] + p[0] * o[1] - p[3] * o[2] + p[2] * o[3]);
        break;
      case 2:
        dpi = (real_t)(0.5f) * (p[2] * o[0] + p[3] * o[1] + p[0] * o[2] - p[1] * o[3]);
        break;
      case 3:
        dpi = (real_t)(0.5f) * (p[3] * o[0] - p[2] * o[1] + p[1] * o[2] + p[0] * o[3]);
        break;
    }
    //integrate the quaternion component and write back to local memory
    l_s[row_stride] = l_s[row_stride] + dt * dpi;
  }
  //make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);
  
  //step 3
  //write back to global
  //again consecutive threads write consecutive data for coalesced access
  g_s[g_id]  = l_s[l_id + shift];
}

/*
 - fixed thread divergence
 - fixed bank conflicts
 - try to improve quaternion mul, but trying to fix thread divergence
   through fancy index magic needs too many additional operations and
   lookup tables. so there is no benefit to circumvent the problem this way
 */
__kernel void euli4(
  __global float * g_s, __global float * g_ds,
  real_t dt,
  index_t nbodies
) {
  __local float l_s[BLOCK_SIZE + SHIFT_SIZE];
  __local float l_ds[BLOCK_SIZE + SHIFT_SIZE];
  index_t id  = get_global_id(0);
  index_t lid = get_local_id(0);
  //if the id is greater than the number of bodies time their state size, kill this thread
  if(id >= nbodies * BODY_STATE_SIZE) return;
  
  //step 1
  //copy global body state for BLOCK_SIZE / BODY_STATE_SIZE bodies into local memory
  //all threads read a linear piece of memory to enable coalesced memory access
  //and introduce a shift to avoid banking conflicts
  //index_t shift = lid / GT330M_BANK_COUNT;
  index_t shift = lid >> 4;
  l_s[lid + shift]  = g_s[id];
  l_ds[lid + shift] = g_ds[id];
  //make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);
  
  //step 2
  //reorganize threads
  if(lid < HALF_BLOCK_SIZE) {
    //the first half performs the translational velocity integration
    //per position update four threads are needed per body
    //the next four threads have to be shifted to the next multipled of BODY_STATE_SIZE
    // 8 * ceil(lid/4) + mod(lid, 4)
    index_t body_stride = ((lid & 0xfffffffc) << 1);
    //index_t l_shift = body_stride / GT330M_BANK_COUNT;
    index_t l_shift = body_stride >> 4;
    index_t x_id  = (lid & 0x03);
    index_t stride =  body_stride + x_id + l_shift;
    //index_t stride = BODY_STATE_SIZE * (lid / 4) + (lid % 4);
    //perform the actual operation
    l_s[stride] = l_s[stride] + dt * l_ds[stride];
  } else {
    //the other half has to do the quaternion update
    //calculate a local half block thread index which runs from 0 to BLOCK_SIZE / 2 - 1
    index_t s_lid = lid - HALF_BLOCK_SIZE;
    //this half block is further grouped into groups of four threads per body
    //determine this threads local index in a 4 group
    //index_t s_dpid = s_lid / QUARTER_HALF_BLOCK_SIZE; s_lid / 64 = s_lid / (2^6)
    index_t s_dpid = s_lid >> 6;
    
    //index_t body_id = s_lid % QUARTER_HALF_BLOCK_SIZE;
    index_t body_id     = s_lid & (QUARTER_HALF_BLOCK_SIZE - 1);
    index_t body_stride = BODY_STATE_SIZE * body_id;
    index_t l_shift = body_stride >> 4;
    index_t stride = body_stride + HALF_BODY_STATE_SIZE + l_shift;
    index_t row_stride = stride + s_dpid;
    
    //load the angular velocity and quaternion from local memory
    float4 p = *(__local float4 *)(l_s + stride);
    float4 o = *(__local float4 *)(l_ds + stride);
    
    //calculate a component of the quaternion time differential where the 4-group local
    //s_dpid determines which component to calculate
    
    //use lookup tables and index magic to make the quaternion multiplication
    //rows look exactly the same from a computational viewpoint so we
    //don't have thread divergence
    __local int idx_l_of[4] = {
      0, 3, 1, 2
    };
    __local int idx_l_inc[4] = {
      1, -1, 1, -1
    };
    __local int idx_r_of[4] = {
      0, 2, 3, 1
    };
    __local float sgn[4] = {-0.5f, 0.5f, 0.5f, 0.5f};
    
    int ai0 = idx_l_of[s_dpid];
    int bi0 = idx_r_of[s_dpid];
    int ai  = idx_l_inc[s_dpid];
    int bi1 = (bi0 + 1) & 0x3;
    int bi2 = (bi0 + 2) & 0x3;
    int bi3 = (bi0 + 3) & 0x3;
    int ai1 = (ai0 + ai) & 0x3;
    int ai2 = (ai1 + ai) & 0x3;
    int ai3 = (ai2 + ai) & 0x3;
    float dpi = sgn[s_dpid] * (- p[ai0] * o[bi0] + p[ai1] * o[bi1] + p[ai2] * o[bi2] + p[ai3] * o[bi3]);
    
    //integrate the quaternion component and write back to local memory
    l_s[row_stride] = l_s[row_stride] + dt * dpi;
  }
  //make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);
  
  //step 3
  //write back to global
  //again consecutive threads write consecutive data for coalesced access
  g_s[id]  = l_s[lid + shift];
}

/* this is based on euli. instead of grouping all the quaternion threads by component into warps
   we use fancy index magic to avoid thread divergence. but this uses too much time and we still
   have banking conflicts.
 */
__kernel void euli5(
  __global float * g_s, __global float * g_ds,
  real_t dt,
  index_t nbodies
) {
  __local float l_s[BLOCK_SIZE];
  __local float l_ds[BLOCK_SIZE];
  index_t id  = get_global_id(0);
  index_t lid = get_local_id(0);
  //if the id is greater than the number of bodies time their state size, kill this thread
  if(id >= nbodies * BODY_STATE_SIZE) return;
  
  //step 1
  //copy global body state for BLOCK_SIZE / BODY_STATE_SIZE bodies into local memory
  //all threads read a linear piece of memory to enable coalesced memory access
  l_s[lid]  = g_s[id];
  l_ds[lid] = g_ds[id];
  //make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);
  
  //step 2
  //reorganize threads
  if(lid < (BLOCK_SIZE / 2)) {
    //the first half performs the translational velocity integration
    //per position update four threads are needed per body
    //the next four threads have to be shifted to the next multipled of BODY_STATE_SIZE
    // 8 * ceil(lid/4) + mod(lid, 4)
    index_t stride =  ((lid & 0xfffffffc) << 1) | (lid & 0x03);
    //index_t stride = BODY_STATE_SIZE * (lid / 4) + (lid % 4);
    //perform the actual operation
    l_s[stride] = l_s[stride] + dt * l_ds[stride];
  } else {
    //the other half has to do the quaternion update
    //calculate a local half block thread index which runs from 0 to BLOCK_SIZE / 2 - 1
    index_t s_lid = lid - (BLOCK_SIZE / 2);
    //this half block is further grouped into groups of four threads per body
    //where each thread in a four group computes one of the four components
    //of the quaternion product 0.5 * p * (s, omega)
    
    //compute the quaternion component to handle by this thread
    //index_t s_dpid = (s_lid % 4);
    index_t s_dpid = (s_lid & 0x03);

    //compute the beginning of the quaternion for the body this thread is working on
    //in local memory
    //index_t stride = BODY_STATE_SIZE * (s_lid / 4) + 4;
    index_t stride = ((s_lid & 0xfffffffc) << 1) + 4;
    index_t row_stride = stride + s_dpid;
    
    //load the angular velocity and quaternion for this body from local memory
    //all threads inside a four group load the same values -> broadcast operation
    float4 p = *(__local float4 *)(l_s + stride);
    float4 o = *(__local float4 *)(l_ds + stride);
    
    //calculate a component of the quaternion time differential where the 4-group local
    //s_dpid determines which component to calculate
    
    //use lookup tables and index magic to make the quaternion multiplication
    //rows look exactly the same from a computational viewpoint so we
    //don't have thread divergence
    __local int idx_l_of[4] = {
      0, 3, 1, 2
    };
    __local int idx_l_inc[4] = {
      1, -1, 1, -1
    };
    __local int idx_r_of[4] = {
      0, 2, 3, 1
    };
    float sgn[4] = {-0.5f, 0.5f, 0.5f, 0.5f};
    
    int ai0 = idx_l_of[s_dpid];
    int bi0 = idx_r_of[s_dpid];
    int ai  = idx_l_inc[s_dpid];
    int bi1 = (bi0 + 1) & 0x3;
    int bi2 = (bi0 + 2) & 0x3;
    int bi3 = (bi0 + 3) & 0x3;
    int ai1 = (ai0 + ai) & 0x3;
    int ai2 = (ai1 + ai) & 0x3;
    int ai3 = (ai2 + ai) & 0x3;
    float dpi = sgn[s_dpid] * (- p[ai0] * o[bi0] + p[ai1] * o[bi1] + p[ai2] * o[bi2] + p[ai3] * o[bi3]);    
    
    //integrate the quaternion component and write back to local memory
    l_s[row_stride] = l_s[row_stride] + dt * dpi;
  }
  //make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);
  
  //step 3
  //write back to global
  //again consecutive threads write consecutive data for coalesced access
  g_s[id]  = l_s[lid];
  //g_ds[id] = l_ds[id];
}

/* this is based on euli5. instead of grouping all the quaternion threads by component into warps
 we use fancy index magic to avoid thread divergence.
  - fixed banking conflicts by introducing a bank shift after every two bodies
 */
__kernel void euli6(
  __global float * g_s, __global float * g_ds,
  real_t dt,
  index_t nbodies
) {
  __local float l_s[BLOCK_SIZE + SHIFT_SIZE];
  __local float l_ds[BLOCK_SIZE + SHIFT_SIZE];
  index_t id  = get_global_id(0);
  index_t lid = get_local_id(0);
  //if the id is greater than the number of bodies time their state size, kill this thread
  if(id >= nbodies * BODY_STATE_SIZE) return;
  
  //step 1
  //copy global body state for BLOCK_SIZE / BODY_STATE_SIZE bodies into local memory
  //all threads read a linear piece of memory to enable coalesced memory access
  //and introduce a shift to avoid banking conflicts
  //index_t shift = lid / GT330M_BANK_COUNT;
  index_t shift = lid >> 4;
  l_s[lid + shift]  = g_s[id];
  l_ds[lid + shift] = g_ds[id];
  //make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);
  
  //step 2
  //reorganize threads
  if(lid < HALF_BLOCK_SIZE) {
    //the first half performs the translational velocity integration
    //per position update four threads are needed per body
    //the next four threads have to be shifted to the next multipled of BODY_STATE_SIZE
    // 8 * ceil(lid/4) + mod(lid, 4)
    index_t body_stride = ((lid & 0xfffffffc) << 1);
    //index_t l_shift = body_stride / GT330M_BANK_COUNT;
    index_t l_shift = bank_line_float(body_stride);
    index_t x_id  = (lid & 0x03);
    index_t stride =  body_stride + x_id + l_shift;
    //index_t stride = BODY_STATE_SIZE * (lid / 4) + (lid % 4);
    //perform the actual operation
    l_s[stride] = l_s[stride] + dt * l_ds[stride];
  } else {
    //the other half has to do the quaternion update
    //calculate a local half block thread index which runs from 0 to BLOCK_SIZE / 2 - 1
    index_t s_lid = lid - HALF_BLOCK_SIZE;
    //this half block is further grouped into groups of four threads per body
    //where each thread in a four group computes one of the four components
    //of the quaternion product 0.5 * p * (s, omega)
    
    //compute the quaternion component to handle by this thread
    //index_t s_dpid = (s_lid % 4);
    index_t s_dpid = (s_lid & 0x03);
    
    //compute the beginning of the quaternion for the body this thread is working on
    //in local memory
    
    index_t body_stride = ((s_lid & 0xfffffffc) << 1);
    index_t l_shift = bank_line_float(body_stride);
    index_t stride = body_stride + HALF_BODY_STATE_SIZE + l_shift;
    index_t row_stride = stride + s_dpid;
    
    //load the angular velocity and quaternion for this body from local memory
    //all threads inside a four group load the same values -> broadcast operation
    float4 p = *(__local float4 *)(l_s + stride);
    float4 o = *(__local float4 *)(l_ds + stride);
    
    //calculate a component of the quaternion time differential where the 4-group local
    //s_dpid determines which component to calculate
    
    //use lookup tables and index magic to make the quaternion multiplication
    //rows look exactly the same from a computational viewpoint so we
    //don't have thread divergence
    __local int idx_l_of[4] = {
      0, 3, 1, 2
    };
    __local int idx_l_inc[4] = {
      1, -1, 1, -1
    };
    __local int idx_r_of[4] = {
      0, 2, 3, 1
    };
    float sgn[4] = {-0.5f, 0.5f, 0.5f, 0.5f};
    
    int ai0 = idx_l_of[s_dpid];
    int bi0 = idx_r_of[s_dpid];
    int ai  = idx_l_inc[s_dpid];
    int bi1 = (bi0 + 1) & 0x3;
    int bi2 = (bi0 + 2) & 0x3;
    int bi3 = (bi0 + 3) & 0x3;
    int ai1 = (ai0 + ai) & 0x3;
    int ai2 = (ai1 + ai) & 0x3;
    int ai3 = (ai2 + ai) & 0x3;
    float dpi = sgn[s_dpid] * (- p[ai0] * o[bi0] + p[ai1] * o[bi1] + p[ai2] * o[bi2] + p[ai3] * o[bi3]);    
    
    //integrate the quaternion component and write back to local memory
    l_s[row_stride] = l_s[row_stride] + dt * dpi;
  }
  //make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);
  
  //step 3
  //write back to global
  //again consecutive threads write consecutive data for coalesced access
  g_s[id]  = l_s[lid + shift];
}

/** split variants: here we use 4 threads per body and perform the integration of the position and orientation
    sequentially
 */

/* this kernel first integrates the positions and then the orientations
   as these integration steps are structurally quite different and hard
   to parallelize
*/
__kernel void euler_split(
  __global float * g_x, __global float * g_dx,
  __global float * g_p, __global float * g_dp,
  real_t dt,
  index_t nbodies
) {
  //get global thread id
  index_t g_id  = get_global_id(0);
  
  //if the global thread id is bigger than the number of bodies times their
  //half state size kill the thread.
  if(g_id >= nbodies * 4) return;
  
  //step 1: integrate positions
  {
    g_x[g_id] = g_x[g_id] + dt * g_dx[g_id];
  }
  
  //step 2: integrate quaternions
  {
    index_t body_id     = g_id / 4;
    index_t p_row_id    = g_id & 0x3; // g_id modulo 4
    index_t body_offset = body_id * 4;
    float4 p = *(__global float4 *)(g_p+body_offset);
    float4 o = *(__global float4 *)(g_dp+body_offset);
    float dpi;
    switch(p_row_id) {
      case 0:
        dpi = (real_t)(0.5) * (p[0] * o[0] - p[1] * o[1] - p[2] * o[2] - p[3] * o[3]);
        break;
      case 1:
        dpi = (real_t)(0.5) * (p[1] * o[0] + p[0] * o[1] - p[3] * o[2] + p[2] * o[3]);
        break;
      case 2:
        dpi = (real_t)(0.5) * (p[2] * o[0] + p[3] * o[1] + p[0] * o[2] - p[1] * o[3]);
        break;
      case 3:
        dpi = (real_t)(0.5) * (p[3] * o[0] - p[2] * o[1] + p[1] * o[2] + p[0] * o[3]);
        break;
    }

    g_p[g_id] = p[p_row_id] + dt * dpi;
  }
}

/* this kernel first integrates the positions and then the orientations
 as these integration steps are structurally quite different and hard
 to parallelize
 - try to fix thread divergence with index magic -> FAIL
 */
__kernel void euler_split2(
                          __global float * g_x, __global float * g_dx,
                          __global float * g_p, __global float * g_dp,
                          real_t dt,
                          index_t nbodies
                          ) {
  //get global thread id
  index_t g_id  = get_global_id(0);
  
  //if the global thread id is bigger than the number of bodies times their
  //half state size kill the thread.
  if(g_id >= nbodies * 4) return;
  
  //step 1: integrate positions
  {
    g_x[g_id] = g_x[g_id] + dt * g_dx[g_id];
  }
  
  //step 2: integrate quaternions
  {
    index_t body_id     = g_id / 4;
    index_t p_row_id    = g_id & 0x3; // g_id modulo 4
    index_t body_offset = body_id * 4;
    float4 p = *(__global float4 *)(g_p+body_offset);
    float4 o = *(__global float4 *)(g_dp+body_offset);
    
    //use lookup tables and index magic to make the quaternion multiplication
    //rows look exactly the same from a computational viewpoint so we
    //don't have thread divergence
    __local int idx_l_of[4] = {
      0, 3, 1, 2
    };
    __local int idx_l_inc[4] = {
      1, -1, 1, -1
    };
    __local int idx_r_of[4] = {
      0, 2, 3, 1
    };
    __local float sgn[4] = {-0.5f, 0.5f, 0.5f, 0.5f};
    
    int ai0 = idx_l_of[p_row_id];
    int bi0 = idx_r_of[p_row_id];
    int ai  = idx_l_inc[p_row_id];
    int bi1 = (bi0 + 1) & 0x3;
    int bi2 = (bi0 + 2) & 0x3;
    int bi3 = (bi0 + 3) & 0x3;
    int ai1 = (ai0 + ai) & 0x3;
    int ai2 = (ai1 + ai) & 0x3;
    int ai3 = (ai2 + ai) & 0x3;
    float dpi = sgn[p_row_id] * (- p[ai0] * o[bi0] + p[ai1] * o[bi1] + p[ai2] * o[bi2] + p[ai3] * o[bi3]);
    
    g_p[g_id] = p[p_row_id] + dt * dpi;
  }
}

/* this kernel first integrates the positions and then the orientations
 as these integration steps are structurally quite different and hard
 to parallelize
 - fix thread divergence by grouping threads based on quat rows -> destroys good coalescing
   properties
 */
__kernel void euler_split3(
                          __global float * g_x, __global float * g_dx,
                          __global float * g_p, __global float * g_dp,
                          real_t dt,
                          index_t nbodies
                          ) {
  //get global thread id
  index_t g_id  = get_global_id(0);
  
  //if the global thread id is bigger than the number of bodies times their
  //half state size kill the thread.
  if(g_id >= nbodies * 4) return;
  
  //step 1: integrate positions
  {
    g_x[g_id] = g_x[g_id] + dt * g_dx[g_id];
  }
  
  //step 2: integrate quaternions
  {
    index_t l_id  = get_local_id(0);
    //regroup threads
    //index_t p_row_id = l_id / QUARTER_BLOCK_SIZE;
    index_t p_row_id = l_id >> 7;
    
    //index_t body_id     = g_id % QUARTER_BLOCK_SIZE;
    index_t body_id     = g_id & (QUARTER_BLOCK_SIZE - 1); //works because power of two

    //index_t body_offset = body_id * 4;
    index_t body_offset = body_id << 2;
    float4 p = *(__global float4 *)(g_p  + body_offset);
    float4 o = *(__global float4 *)(g_dp + body_offset);
    float dpi;
    switch(p_row_id) {
      case 0:
        dpi = (real_t)(0.5) * (p[0] * o[0] - p[1] * o[1] - p[2] * o[2] - p[3] * o[3]);
        break;
      case 1:
        dpi = (real_t)(0.5) * (p[1] * o[0] + p[0] * o[1] - p[3] * o[2] + p[2] * o[3]);
        break;
      case 2:
        dpi = (real_t)(0.5) * (p[2] * o[0] + p[3] * o[1] + p[0] * o[2] - p[1] * o[3]);
        break;
      case 3:
        dpi = (real_t)(0.5) * (p[3] * o[0] - p[2] * o[1] + p[1] * o[2] + p[0] * o[3]);
        break;
    }
    
    g_p[g_id] = p[p_row_id] + dt * dpi;
  }
}

/* PERFORMANCE: faster than euli7 / euler_split5
 
 this kernel first integrates the positions and then the orientations
 as these integration steps are structurally quite different and hard
 to parallelize
 - fix thread divergence by grouping threads based on quat rows
 - use local memory to keep coalescing behaviour
 */
#undef SHIFT_SIZE
#define SHIFT_SIZE WORK_GROUP_SIZE / LOCAL_MEM_BANK_SIZE
__kernel void euler_split4(
  __global float * g_x, __global float * g_dx,
  __global float * g_p, __global float * g_dp,
  real_t dt,
  index_t nbodies
) {
  //get global thread id
  index_t g_id  = get_global_id(0);
  
  //if the global thread id is bigger than the number of bodies times their
  //half state size kill the thread.
  if(g_id >= nbodies * BODY_X_SIZE) return;
  
  //step 1: integrate positions
  {
    g_x[g_id] = g_x[g_id] + dt * g_dx[g_id];
  }
  
  //step 2: integrate quaternions
  {
    __local float l_p[BLOCK_SIZE + SHIFT_SIZE];
    __local float l_dp[BLOCK_SIZE + SHIFT_SIZE];

    index_t l_id  = get_local_id(0);
    
    //step 2.1.1 load quaternion data from global into local memory for good
    //coalescing. keep in registers until we have calculated to local indices
    float a = g_p[g_id];
    float b = g_dp[g_id];
    
    //step 2.1.2 introduce for every 16 floats a shift of one float to 
    //avoid banking errors during the local operations
    index_t shift = l_id >> 4;
    index_t l_id_shifted = l_id + shift;
    
    //step 2.1.3 safe values from global memory into local memory
    l_p[l_id_shifted]  = a;
    l_dp[l_id_shifted] = b;
    
    
    //step 2.2 regroup threads to minimize thread divergence. threads in the same warp
    //should evaluate the same quaternion row
    //this is done by dividing the whole work group into four groups
    
    //determine the row index of the thread [0, 4)
    //index_t p_row_id = l_id / QUARTER_BLOCK_SIZE;
    index_t p_row_id = l_id >> 7;
    
    //determine the body index of the thread [0, QUARTER_BLOCK_SIZE)
    //index_t body_id     = g_id % QUARTER_BLOCK_SIZE;
    index_t body_id  = g_id & (QUARTER_BLOCK_SIZE - 1); //works because power of two
    
    //determine the offset of the body in the local memory
    //index_t body_offset = body_id * 4;
    index_t body_offset = body_id << 2;
    
    //determine the bank shift for this offset
    index_t offset_shift = body_offset >> 4;
    index_t shifted_body_offset = body_offset + offset_shift;
    
    
    //now make sure that all the data in the local memory is ready
    mem_fence(CLK_LOCAL_MEM_FENCE);
    
    //step 2.3 load the quaternion and the angular velocity vector from local memory
    float4 p = *(__local float4 *)(l_p  + shifted_body_offset);
    float4 o = *(__local float4 *)(l_dp + shifted_body_offset);
    
    //compute the row offset with banking shift for this thread
    index_t shifted_row_offset = shifted_body_offset + p_row_id;
    
    //step 2.4 compute row
    float dpi;
    switch(p_row_id) {
      case 0:
        dpi = (real_t)(0.5) * (p[0] * o[0] - p[1] * o[1] - p[2] * o[2] - p[3] * o[3]);
        break;
      case 1:
        dpi = (real_t)(0.5) * (p[1] * o[0] + p[0] * o[1] - p[3] * o[2] + p[2] * o[3]);
        break;
      case 2:
        dpi = (real_t)(0.5) * (p[2] * o[0] + p[3] * o[1] + p[0] * o[2] - p[1] * o[3]);
        break;
      case 3:
        dpi = (real_t)(0.5) * (p[3] * o[0] - p[2] * o[1] + p[1] * o[2] + p[0] * o[3]);
        break;
    }
    //step 2.5 write back to local memory
    l_dp[shifted_row_offset] = dpi;
  
    //now make sure that all the data in the local memory is ready
    mem_fence(CLK_LOCAL_MEM_FENCE);
    
    //step 2.6 write back to global memory with old thread ordering for coalescing
    g_p[g_id] = a + dt * l_dp[l_id_shifted];
  }
}

/* base on euli3 new favourite
 - fixed thread divergence
 - fixed banking conflicts
 */
__kernel void euli7(
  __global float * g_s, __global float * g_ds,
  real_t dt,
  index_t nbodies
) {
  __local float l_s[BLOCK_SIZE + SHIFT_SIZE];
  __local float l_ds[BLOCK_SIZE + SHIFT_SIZE];
  index_t g_id  = get_global_id(0);
  index_t l_id = get_local_id(0);
  //if the id is greater than the number of bodies time their state size, kill this thread
  if(g_id >= nbodies * BODY_STATE_SIZE) return;
  
  //step 1
  //copy global body state for BLOCK_SIZE / BODY_STATE_SIZE bodies into local memory
  //all threads read a linear piece of memory to enable coalesced memory access
  //and introduce a shift
  //index_t shift = lid / GT330M_BANK_COUNT;
  index_t shift = bank_line_float(l_id);
  index_t shifted_l_id = l_id + shift;
  float state = g_s[g_id];
  l_s[shifted_l_id]  = state;
  l_ds[shifted_l_id] = g_ds[g_id];
  
  //step 2
  //reorganize threads
  if(l_id < HALF_BLOCK_SIZE) {
    //the first half of threads has to do the quaternion update
    //this half block is further grouped into groups of four threads per body
    //determine this threads local index in a 4 group
    //index_t row_id = lid / QUARTER_HALF_BLOCK_SIZE; lid / 64 = lid / (2^6)
    index_t row_id = l_id >> 6;
    
    //index_t body_id = l_id % QUARTER_HALF_BLOCK_SIZE;
    index_t body_id     = l_id & (QUARTER_HALF_BLOCK_SIZE - 1);
    index_t body_stride = BODY_STATE_SIZE * body_id;
    index_t l_shift = bank_line_float(body_stride);
    index_t stride = body_stride + HALF_BODY_STATE_SIZE + l_shift;
    index_t row_stride = stride + row_id;
    
    //make sure that all the data in the local memory is ready
    mem_fence(CLK_LOCAL_MEM_FENCE);
    
    //load the angular velocity and quaternion from local memory
    float4 p = *(__local float4 *)(l_s + stride);
    float4 o = *(__local float4 *)(l_ds + stride);
    
    //calculate a component of the quaternion time differential where the 4-group local
    //thread id determines which component to calculate
    float dpi;
    switch(row_id) {
      case 0:
        dpi = (real_t)(0.5f) * (p[0] * o[0] - p[1] * o[1] - p[2] * o[2] - p[3] * o[3]);
        break;
      case 1:
        dpi = (real_t)(0.5f) * (p[1] * o[0] + p[0] * o[1] - p[3] * o[2] + p[2] * o[3]);
        break;
      case 2:
        dpi = (real_t)(0.5f) * (p[2] * o[0] + p[3] * o[1] + p[0] * o[2] - p[1] * o[3]);
        break;
      case 3:
        dpi = (real_t)(0.5f) * (p[3] * o[0] - p[2] * o[1] + p[1] * o[2] + p[0] * o[3]);
        break;
    }
    //write quaternion time differential to local memory
    l_ds[row_stride] = dpi;
  }
  //make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);
  
  //load state differential from local memory
  float dsi = dt * l_ds[shifted_l_id];
  //step 3
  //perform euler step and write back to global
  //write back to global
  //again consecutive threads write consecutive data for coalesced access
  g_s[g_id]  = state + dsi;
}

/* PERFORMANCE SIMILAR TO euli7
 
 this kernel first integrates the positions and then the orientations
 as these integration steps are structurally quite different and hard
 to parallelize
 - fix thread divergence by grouping threads based on quat rows
 - use local memory to keep coalescing behaviour
 */
__kernel void euler_split5(
  __global float * g_x, __global float * g_dx,
  __global float * g_p, __global float * g_dp,
  real_t dt,
  index_t nbodies
) {
  //get global thread id
  index_t g_id  = get_global_id(0);
  index_t g_gid = get_group_id(0);
  index_t l_id  = get_local_id(0);
  //fix global id
  g_id = g_id - g_gid * HALF_BLOCK_SIZE;
  
  //if the global thread id is bigger than the number of bodies times their
  //half state size kill the thread.
  if(g_id >= nbodies * 4) return;
  
  //step 1: integrate positions
  if(l_id < HALF_BLOCK_SIZE) {
    g_x[g_id] = g_x[g_id] + dt * g_dx[g_id];
  }
  else {
    l_id = l_id & (HALF_BLOCK_SIZE - 1);
    g_id = g_id - HALF_BLOCK_SIZE;
    float a = g_p[g_id];
    float b = g_dp[g_id];
    //step 2: integrate quaternions
    __local float l_p[HALF_BLOCK_SIZE + SHIFT_SIZE];
    __local float l_dp[HALF_BLOCK_SIZE + SHIFT_SIZE];
    
    //step 2.1.1 load quaternion data from global into local memory for good
    //coalescing. keep in registers until we have calculated to local indices
    
    
    //step 2.1.2 introduce for every 16 floats a shift of one float to 
    //avoid banking errors during the local operations
    index_t shift = bank_line_float(l_id);
    index_t l_id_shifted = l_id + shift;
    
    //step 2.1.3 safe values from global memory into local memory
    l_p[l_id_shifted]  = a;
    l_dp[l_id_shifted] = b;
    
    
    //step 2.2 regroup threads to minimize thread divergence. threads in the same warp
    //should evaluate the same quaternion row
    //this is done by dividing the whole work group into four groups
    
    //determine the row index of the thread [0, 4)
    //index_t p_row_id = l_id / QUARTER_HALF_BLOCK_SIZE;
    index_t p_row_id = l_id >> 6;
    
    //determine the body index of the thread [0, QUARTER_HALF_BLOCK_SIZE)
    //index_t body_id     = g_id % QUARTER_HALF_BLOCK_SIZE;
    index_t body_id  = g_id & (QUARTER_HALF_BLOCK_SIZE - 1); //works because power of two
    
    //determine the offset of the body in the local memory
    //index_t body_offset = body_id * 4;
    index_t body_offset = body_id << 2;
    
    //determine the bank shift for this offset
    index_t offset_shift = bank_line_float(body_offset);
    index_t shifted_body_offset = body_offset + offset_shift;
    
    //now make sure that all the data in the local memory is ready
    mem_fence(CLK_LOCAL_MEM_FENCE);
    
    //step 2.3 load the quaternion and the angular velocity vector from local memory
    float4 p = *(__local float4 *)(l_p  + shifted_body_offset);
    float4 o = *(__local float4 *)(l_dp + shifted_body_offset);
    
    //compute the row offset with banking shift for this thread
    index_t shifted_row_offset = shifted_body_offset + p_row_id;
    
    //step 2.4 compute row
    float dpi;
    switch(p_row_id) {
      case 0:
        dpi = (real_t)(0.5) * (p[0] * o[0] - p[1] * o[1] - p[2] * o[2] - p[3] * o[3]);
        break;
      case 1:
        dpi = (real_t)(0.5) * (p[1] * o[0] + p[0] * o[1] - p[3] * o[2] + p[2] * o[3]);
        break;
      case 2:
        dpi = (real_t)(0.5) * (p[2] * o[0] + p[3] * o[1] + p[0] * o[2] - p[1] * o[3]);
        break;
      case 3:
        dpi = (real_t)(0.5) * (p[3] * o[0] - p[2] * o[1] + p[1] * o[2] + p[0] * o[3]);
        break;
    }
    //step 2.5 write back to local memory
    l_dp[shifted_row_offset] = dpi;
    
    //now make sure that all the data in the local memory is ready
    mem_fence(CLK_LOCAL_MEM_FENCE);
    
    //step 2.6 write back to global memory with old thread ordering for coalescing
    g_p[g_id] = a + dt * l_dp[l_id_shifted];
  }
}

__kernel void euler_position(
  __global float * g_x, __global float * g_dx,
  real_t dt,
  index_t nbodies
) {
  //get global thread id
  index_t g_id  = get_global_id(0);
  
  //if the global thread id is bigger than the number of bodies times their
  //half state size kill the thread.
  if(g_id >= nbodies * BODY_X_SIZE) return;
  
  //step 1: integrate positions
  g_x[g_id] = g_x[g_id] + dt * g_dx[g_id];
}

#undef SHIFT_SIZE
#define SHIFT_SIZE WORK_GROUP_SIZE / LOCAL_MEM_BANK_SIZE
__kernel void euler_orientation(
  __global float * g_p, __global float * g_dp,
  real_t dt,
  index_t nbodies
) {
  __local float l_p[WORK_GROUP_SIZE + SHIFT_SIZE];
  __local float l_dp[WORK_GROUP_SIZE + SHIFT_SIZE];
  //get global thread id
  index_t g_id  = get_global_id(0);
  index_t l_id  = get_local_id(0);
  
  //if the global thread id is bigger than the number of bodies times their
  //half state size kill the thread.
  if(g_id >= nbodies * BODY_P_SIZE) return;
  
  //step 2.1.1 load quaternion data from global into local memory for good
  //coalescing. keep in registers until we have calculated to local indices
  float a = g_p[g_id];
  float b = g_dp[g_id];

  //step 2.1.2 introduce for every 16 floats a shift of one float to 
  //avoid banking conflicts during the local operations
  index_t shift = bank_line_float(l_id);
  index_t l_id_shifted = l_id + shift;

  //step 2.1.3 safe values from global memory into local memory
  l_p[l_id_shifted]  = a;
  l_dp[l_id_shifted] = b;

  //step 2.2 regroup threads to minimize thread divergence. threads in the same warp
  //should evaluate the same quaternion row
  //this is done by dividing the whole work group into four groups

  //determine the row index of the thread [0, 4)
  //index_t p_row_id = l_id / QUARTER_BLOCK_SIZE;
  index_t p_row_id = l_id >> 7;

  //determine the body index of the thread [0, QUARTER_BLOCK_SIZE)
  //index_t body_id     = g_id % QUARTER_BLOCK_SIZE;
  index_t body_id  = l_id & (QUARTER_BLOCK_SIZE - 1); //works because power of two

  //determine the offset of the body in the local memory
  //index_t body_offset = body_id * 4;
  index_t body_offset = body_id << 2;

  //determine the bank shift for this offset
  index_t offset_shift = bank_line_float(body_offset);
  index_t shifted_body_offset = body_offset + offset_shift;

  //now make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);

  //step 2.3 load the quaternion and the angular velocity vector from local memory
  float4 p = *(__local float4 *)(l_p  + shifted_body_offset);
  float4 o = *(__local float4 *)(l_dp + shifted_body_offset);

  //compute the row offset with banking shift for this thread
  index_t shifted_row_offset = shifted_body_offset + p_row_id;

  //step 2.4 compute row
  float dpi;
  switch(p_row_id) {
    case 0:
      dpi = (real_t)(0.5) * (p[0] * o[0] - p[1] * o[1] - p[2] * o[2] - p[3] * o[3]);
      break;
    case 1:
      dpi = (real_t)(0.5) * (p[1] * o[0] + p[0] * o[1] - p[3] * o[2] + p[2] * o[3]);
      break;
    case 2:
      dpi = (real_t)(0.5) * (p[2] * o[0] + p[3] * o[1] + p[0] * o[2] - p[1] * o[3]);
      break;
    case 3:
      dpi = (real_t)(0.5) * (p[3] * o[0] - p[2] * o[1] + p[1] * o[2] + p[0] * o[3]);
      break;
  }
  //step 2.5 write back to local memory
  l_dp[shifted_row_offset] = dpi;

  //now make sure that all the data in the local memory is ready
  mem_fence(CLK_LOCAL_MEM_FENCE);

  //step 2.6 write back to global memory with old thread ordering for coalescing
  g_p[g_id] = a + dt * l_dp[l_id_shifted];
}



