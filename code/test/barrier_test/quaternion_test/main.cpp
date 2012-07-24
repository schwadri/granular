//
//  main.cpp
//  quaternion_test
//
//  Created by Adrian Schweizer on 7/3/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#include <iostream>

float signum(int i) {
  return (i > 0) * 2.0f - 1.0f;
}

int main(int argc, const char * argv[])
{
  
  /* this is a test to see whether it is possible to
   * code the quaternion multiplication in such a way
   * that each component of the quaternion is calculated
   * in exactly the same way so they could be done in parallel
   */
  int a[4] = {1, 2, 3, 5};
  int b[4] = {7, 11, 13, 17};
  int c[4];
  int d[4];
  
  //quaternion multiplication hardcoded
  c[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
  c[1] = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2];
  c[2] = a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1];
  c[3] = a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0];
  
  //with lookup table to make all lines look equal
  int idx_l_of[4] = {
    0, 3, 1, 2
  };
  int idx_l_inc[4] = {
    1, -1, 1, -1
  };
  int idx_r_of[4] = {
    0, 2, 3, 1
  };
  
  float const sgn[4] = {-1.0f, 1.0f, 1.0f, 1.0f};
  
  for(int i = 0; i < 4; ++i) {
    int ai0 = idx_l_of[i];
    int bi0 = idx_r_of[i];
    int ai  = idx_l_inc[i];
    int bi1 = (bi0 + 1) & 0x3;
    int bi2 = (bi0 + 2) & 0x3;
    int bi3 = (bi0 + 3) & 0x3;
    int ai1 = (ai0 + ai) & 0x3;
    int ai2 = (ai1 + ai) & 0x3;
    int ai3 = (ai2 + ai) & 0x3;
    d[i] = sgn[i] * (
      - a[ai0] * b[bi0] 
      + a[ai1] * b[bi1]
      + a[ai2] * b[bi2]
      + a[ai3] * b[bi3]
    );
  }

  std::cout << " c = " << c[0] << ", " << c[1] << ", " << c[2] << ", " << c[3] << '\n';
  std::cout << " d = " << d[0] << ", " << d[1] << ", " << d[2] << ", " << d[3] << '\n';
  
  for(int i = 0; i < 4; ++i)
  std::cout << " sgn(" << i << ") = " << signum(i) << std::endl;
  
  return 0;
}

