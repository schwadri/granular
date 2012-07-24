#include <iostream>
#include <xmmintrin.h>

#include "../core/util/time.hpp"
#include "../core/util/real.hpp"

typedef float real;
#include "../core/util/progress_bar.hpp"

struct vec4 {
  __m128 values;
};

int main() {
  //vec4 * bla = new vec4;
  //bla->values = _mm_setzero_ps();

  int n = 1000000;

  float *a, *b;

  a = new float[4 * n];
  b = new float[4 * n];
  //a = __builtin_assume_aligned(a, 16);
  //b = __builtin_assume_aligned(b, 16);

  util::time_pos t;
  std::cout << t.seconds << " s " << t.micro_seconds << " us " << std::endl;

  for(int i = 0; i < n; ++i)
    a[i] += b[i];
  //*((__m128*)(&a[i])) =  *((__m128*)&a[i]) + *((__m128*)&b[i]);
  //*((__m128*)(&a[i])) = __builtin_ia32_addps(*(__m128*)&a[i], *(__m128*)&b[i]);

  util::duration dt = util::time_pos() - t;
  std::cout << dt << std::endl;


  for(int i = 0; i < n; ++i)
    a[i];

  delete[] a;
  delete[] b;

  return EXIT_SUCCESS;
}