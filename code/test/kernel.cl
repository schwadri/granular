typedef float real;
//#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable
//#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
//#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
//#pragma OPENCL EXTENSION cl_amd_fp64 : enable
//#pragma OPENCL EXTENSION cl_khr_fp64 : enable

real prox_r(real const p) {
  return max((real)0.0, p);
}

struct row {
  unsigned int boffset;
  unsigned int bcount;
  real r_n, r_t1, r_t2;
};

__kernel void jprox(
  __global real * pnew,
  __global real const * pold,
  __global real const * c,
  __global struct row const * rows,
  __global unsigned int const * bj,
  __global real const * g,
  real const mu
)
{
  size_t ri  = get_global_id(0);
  size_t i   = ri * 3;
  real XN  = c[i];
  real XT1 = c[i + 1],
       XT2 = c[i + 2];
  __global struct row const * r = &rows[ri];
  __global real const * rg = g + r->boffset * 9;
  __global unsigned int const * rj  = bj + r->boffset;

  for(int bi = 0; bi < r->bcount; ++bi) {
    __global real const * bp = &pold[rj[bi]];
    for(int j = 0; j < 3; ++j) {
      XN  += rg[0 + j] * bp[j];
      XT1 += rg[3 + j] * bp[j];
      XT2 += rg[6 + j] * bp[j];
    }
    rg += 9;
  }

  real PN   = prox_r(pold[i] - r->r_n * XN);
  real PT1  = pold[i + 1] - r->r_t1 * XT1;
  real PT2  = pold[i + 2] - r->r_t2 * XT2;
  real PT_norm = PT1 * PT1 + PT2 * PT2;
  real a = mu * PN;
  if(PT_norm > a * a) {
    PT_norm = a / sqrt(PT_norm);
    PT1 *= PT_norm;
    PT2 *= PT_norm;
  }
  pnew[i] = PN;
  pnew[i + 1] = PT1;
  pnew[i + 2] = PT2;
}
