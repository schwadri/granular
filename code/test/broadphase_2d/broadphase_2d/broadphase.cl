typedef unsigned int index_t;

#define BUILD_FOR_NV_GT330M

//hardware specific constants and utilities

inline index_t bank_line_float(index_t i);
inline index_t bank_line_double(index_t i);
inline index_t bank_line_uint(index_t i);

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
inline index_t bank_line_uint(index_t i) {
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
  //i / 32
  return i >> 5;
}
inline index_t bank_line_double(index_t i) {
  //i / 16
  return i >> 4;
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

/* this kernel takes the coordinates of the com of a single rigid body
   and returns an cell index for a regular 2d grid, which contains the com
   of the body. Every body is handled by a separate thread. It returns an array
   of (body_id, cell_id) pairs. This list is subsequently sorted by cell_id.
 */
__kernel void map_body_to_cell_v4(
  __global float4 * g_x, __global float4 * g_p,
  __global uint2 * g_c,
  float4 l_x, uint4 n_c,
  index_t nbodies
) {
  index_t g_id = get_global_id(0);
  index_t body_id = g_id;
  if(body_id >= nbodies) return;
  float4 x  = g_x[g_id];
  uint cx   = clamp((index_t)floor(x.x / l_x[0]), (index_t)0, (index_t)(n_c[0] - 1));
  uint cy   = clamp((index_t)floor(x.y / l_x[1]), (index_t)0, (index_t)(n_c[1] - 1));
  uint cell_id = cx + n_c[0] * cy;
  g_c[body_id] = (uint2)(body_id, cell_id);
}

/* way slower than the v4 version
 */
#define BANK_SHIFT (WORK_GROUP_SIZE / LOCAL_MEM_BANK_COUNT)
__kernel void map_body_to_cell(
  __global float * g_x, __global float * g_p,
  __global uint  * g_c,
  float4 l_x, uint4 n_c,
  index_t nbodies
) {
  index_t g_id = get_global_id(0);
  if(g_id >= nbodies * 4) return;
  index_t l_id = get_local_id(0);

  __local uint l_cell_ids[WORK_GROUP_SIZE + BANK_SHIFT];
  
  index_t shifted_l_id = l_id + bank_line_uint(l_id);
  index_t component_id  = g_id & 0x3;
  float   x_i  = g_x[g_id];
  index_t cx_i  = clamp(
    (index_t)floor(x_i / l_x[component_id]),
    (index_t)0,
    (index_t)(n_c[component_id] - 1)
  );
  l_cell_ids[shifted_l_id] = cx_i;
  
  if(l_id < QUARTER_WORK_GROUP_SIZE) {

    index_t offset = l_id << 2;
    index_t shift = bank_line_uint(offset);
    mem_fence(CLK_LOCAL_MEM_FENCE);
    index_t cx_0 = l_cell_ids[offset + shift];
    index_t cx_1 = l_cell_ids[offset + shift + 1];
    index_t cell_id = cx_0 + n_c[0] * cx_1;
    index_t body_id       = g_id - (3 * QUARTER_WORK_GROUP_SIZE) * get_group_id(0);
    g_c[body_id] = cell_id;
  }
}

