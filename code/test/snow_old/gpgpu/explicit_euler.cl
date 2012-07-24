/*
 *  explicit_euler.c
 *  gpgpu
 *
 *  Created by Adrian Schweizer on 11/22/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */


__kernel void euler(
  __global float4 * s, __global float4 * ds, float dt
) {
  unsigned int i = get_global_id(0);
      
  float4 x = s[i];
  float4 v = ds[i];

  x = x + dt * v;

  s[i]  = x;
}
