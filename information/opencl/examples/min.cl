/*
 *  min.cl
 *  gpgpu
 *
 *  Created by Adrian Schweizer on 11/22/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */


#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable
 
// 9. The source buffer is accessed as 4-vectors.

__kernel void minp(
  __global uint4 * src,
  __global  uint * gmin,
  __local   uint * lmin,
  __global  uint * dbg,
  size_t nitems,
  uint dev
) {
 // 10. Set up __global memory access pattern.

  uint count = ( nitems / 4 ) / get_global_size(0);
  uint idx = (dev == 0) ? get_global_id(0) * count : get_global_id(0); 
  uint stride = (dev == 0) ? 1 : get_global_size(0); 
  uint pmin = (uint) -1;

  // 11. First, compute private min, for this work-item. 

  for(int n=0; n < count; n++, idx += stride) {
    pmin = min(pmin, src[idx].x ); 
    pmin = min(pmin, src[idx].y ); 
    pmin = min(pmin, src[idx].z ); 
    pmin = min(pmin, src[idx].w ); 
  }

  // 12. Reduce min values inside work-group. 

  if(get_local_id(0) == 0) 
    lmin[0] = (uint) -1; 

  barrier(CLK_LOCAL_MEM_FENCE); 

  lmin = min(lmin, pmin); //(void) atom_min( lmin, pmin );

  barrier(CLK_LOCAL_MEM_FENCE);

  // Write out to __global.

  if(get_local_id(0) == 0)
    gmin[get_group_id(0)] = lmin[0];

  // Dump some debug information.

  if(get_global_id(0) == 0) { 
    dbg[0] = get_num_groups(0);
    dbg[1] = get_global_size(0);
    dbg[2] = count;
    dbg[3] = stride;
  }
}

// 13. Reduce work-group min values from __global to __global.
 
__kernel void reduce(__global uint4 *src,
 __global uint *gmin) { 
 (void) atom_min( gmin, gmin[get_global_id(0)] ) ;
} 