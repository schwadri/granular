/*
 *  csr_jor_prox.cl
 *  gpgpu
 *
 *  Created by Adrian Schweizer on 11/22/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

typedef unsigned int  index_t;
typedef float         real_t;
typedef unsigned char bool_t;

/** sparse JORprox for unilateral contacts without friction
  unoptimized!
*/
__kernel void csr_jor_prox(
  __global index_t * rows,
  __global index_t * columns,
  __global real_t * tij,
  __global real_t * b,
  __global real_t * la0,
  __global real_t * la1,
  __global bool_t * converged,
  real_t tol_rel,
  index_t size
) {

	index_t i = get_global_id(0);
  if(i == 0)
		converged[0] = true;
	if(i < size) {
		__global index_t * jptr     = columns + rows[i];
		__global index_t * jptr_end = columns + rows[i + 1];
		__global real_t * tptr      = tij + rows[i];
		
		real_t la = b[i];
		for(;jptr < jptr_end; ++jptr, ++tptr) {
			la = la + (*tptr) * la0[*jptr];
		}
		//prox r^+
		la1[i] = max((real_t)0, la);
		//check convergence criterion
		bool convgd = fabs(la1[i] - la0[i]) <= tol_rel * fabs(la0[i]);
		if(!convgd)
			converged[0] = false;
	}
}

__kernel void all_converged(
  __global bool_t * converged,
  index_t size
) {

	index_t i = get_global_id(0);
	if(i < size) {
		if(!converged[i])
      converged[0] = false;
	}
}

