typedef unsigned int index_t;

__kernel void barrier_test(
  __global float * g_s,
  index_t sizes
) {
  index_t id  = get_global_id(0);
  //thread 0 continues, thread 1 stops immediately
  if(id >= 4) return;

  float val  = g_s[id];
  float4 p;
  p[0] = val;
  float4 o;
  
  //lookup tables to make all rows of quaternion multiplication look equal
  int const idx_l_of[4] = {
    0, 3, 1, 2
  };
  int const  idx_l_inc[4] = {
    1, -1, 1, -1
  };
  int const  idx_r_of[4] = {
    0, 2, 3, 1
  };
  float const sgn[4] = {-1.0f, 1.0f, 1.0f, 1.0f};
  
    int ai0 = idx_l_of[id];
    int bi0 = idx_r_of[id];
    int ai  = idx_l_inc[id];
    int bi1 = (bi0 + 1) & 0x3;
    int bi2 = (bi0 + 2) & 0x3;
    int bi3 = (bi0 + 3) & 0x3;
    int ai1 = (ai0 + ai) & 0x3;
    int ai2 = (ai1 + ai) & 0x3;
    int ai3 = (ai2 + ai) & 0x3;
    float r = sgn[id] * (
                     - p[ai0] * o[bi0] 
                     + p[ai1] * o[bi1]
                     + p[ai2] * o[bi2]
                     + p[ai3] * o[bi3]
                     );

  //this barrier is never going to be reached by thread 1... what happens?
  barrier(CLK_GLOBAL_MEM_FENCE);
 
  g_s[id]  = r;

}

float signum(int i);

float signum(int i) {
  return (i > 0) * 2.0f - 1.0f;
}

__kernel void quaternion_test(
                           __global float * g_s,
                           index_t sizes
                           ) {
  index_t id  = get_global_id(0);
  //thread 0 continues, thread 1 stops immediately
  if(id >= 4) return;
  
  //float val  = g_s[id];
  float4 p = (float4)(1.0f, 2.0f, 3.0f, 5.0f);
  float4 o = (float4)(7.0f, 11.0f, 13.0f, 17.0f);
  
  //lookup tables to make all rows of quaternion multiplication look equal
  __private int idx_l_of[4] = {
    0, 3, 1, 2
  };
  __private int idx_l_inc[4] = {
    1, -1, 1, -1
  };
  __private int idx_r_of[4] = {
    0, 2, 3, 1
  };
  //float const sgn[4] = {-1.0f, 1.0f, 1.0f, 1.0f};
  
  int ai0 = idx_l_of[id];
  int bi0 = idx_r_of[id];
  int ai  = idx_l_inc[id];
  int bi1 = bi0 + 1;
  int bi2 = bi0 + 2;
  int bi3 = bi0 + 3;
  int ai1 = (ai0 + ai) & 0x3;
  int ai2 = (ai1 + ai) & 0x3;
  int ai3 = (ai2 + ai) & 0x3;
  float r = signum(id) * (
                       - p[ai0] * o[bi0] 
                       + p[ai1] * o[bi1]
                       + p[ai2] * o[bi2]
                       + p[ai3] * o[bi3]
                       );
  
  //this barrier is never going to be reached by thread 1... what happens?
  
  g_s[id]  = r;
  
}



