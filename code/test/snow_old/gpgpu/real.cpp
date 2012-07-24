#include <iostream>   //using cout, cerr
#include <limits>     //using numeric_limits
#include <cmath>
#include <tr1/cstdint>

using std::tr1::uint64_t;
using std::tr1::int64_t;
using std::tr1::uint32_t;
using std::tr1::int32_t;

bool is_nan(float f) {
  uint32_t const & b = reinterpret_cast<uint32_t const &>(f);
  return (b & 0x7f800000) == 0x7f800000 && (b & 0x007fffff);
}

bool is_nan(double f) {
  uint64_t const & b = reinterpret_cast<uint64_t const &>(f);
  return (b & 0x7ff0000000000000l) == 0x7ff0000000000000l && (b & 0x000fffffffffffffl);
}

bool is_snan(double f) {
  uint64_t const & b = reinterpret_cast<uint64_t const &>(f);
  return (b & 0x7ff0000000000000l) == 0x7ff0000000000000l && (b & 0x000fffffffffffffl) && !(b & 0x0008000000000000l);
}

bool is_qnan(double f) {
  uint64_t const & b = reinterpret_cast<uint64_t const &>(f);
  return (b & 0x7ff0000000000000l) == 0x7ff0000000000000l && (b & 0x000fffffffffffffl) && (b & 0x0008000000000000l);
}

bool is_snan(float f) {
  uint32_t const & b = reinterpret_cast<uint32_t const &>(f);
  return (b & 0x7f800000) == 0x7f800000 && (b & 0x007fffff) && !(b & 0x00400000);
}

bool is_qnan(float f) {
  uint32_t const & b = reinterpret_cast<uint32_t const &>(f);
  return (b & 0x7f800000) == 0x7f800000 && (b & 0x007fffff) && (b & 0x00400000);
}

bool is_inf(double f) {
  uint64_t const & b = reinterpret_cast<uint64_t const &>(f);
  return (b & 0x7fffffffffffffffl) == 0x7ff0000000000000l;
}

bool is_inf(float f) {
  uint32_t const & b = reinterpret_cast<uint32_t const &>(f);
  return (b & 0x7fffffff) == 0x7f800000;
}

uint32_t payload(float f) {
  uint32_t const & b = reinterpret_cast<uint32_t const &>(f);
  return (b & 0x003fffff);
}

uint64_t payload(double f) {
  uint64_t const & b = reinterpret_cast<uint64_t const &>(f);
  return (b & 0x0007ffffffffffffl);
}

uint64_t mantissa(double f) {
  uint64_t const & b = reinterpret_cast<uint64_t const &>(f);
  return (b & 0x000fffffffffffffl);
}

int64_t exponent(double f) {
  int64_t const & b = reinterpret_cast<int64_t const &>(f);
  return ((b & 0x7ff0000000000000l) >> 52) - 1023;
}

uint32_t mantissa(float v) {
  uint32_t bits = reinterpret_cast<uint32_t const &>(v);
  return bits & 0x007fffff;
}

int32_t exponent(float v) {
  int32_t bits = reinterpret_cast<int32_t const &>(v);
  return ((bits & 0x7f800000) >> 23) - 127;
}

double sign(double v) {
  uint64_t bits = reinterpret_cast<uint64_t const &>(v);
  return bits & 0x8000000000000000l ? -1.0 : 1.0;
}

float sign(float v) {
  uint32_t bits = reinterpret_cast<uint32_t const &>(v);
  return bits & 0x80000000 ? -1.0 : 1.0;
}


double next(double prev) {
  int64_t bits = reinterpret_cast<int64_t const &>(prev);

  if(bits < 0) {
    --bits;
    if(bits > 0)
      bits = 0;
  } else
    ++bits;
  return reinterpret_cast<double const &>(bits);
}

double prev(double v) {
  int64_t bits = reinterpret_cast<int64_t const &>(v);

  if(bits < 0)
    ++bits;
  else {
    --bits;
    if(bits < 0)
      bits = 0x8000000000000000l;
  }
  return reinterpret_cast<double const &>(bits);
}

float next(float prev) {
  int32_t bits = reinterpret_cast<int32_t const &>(prev);

  if(bits < 0) {
    --bits;
    if(bits > 0)
      bits = 0;
  } else
    ++bits;

  return reinterpret_cast<float const &>(bits);
}

float prev(float v) {
  int32_t bits = reinterpret_cast<int32_t const &>(v);

  if(bits < 0)
    ++bits;
  else {
    --bits;
    if(bits < 0)
      bits = 0x80000000;
  }
  return reinterpret_cast<float const &>(bits);
}

uint32_t ulp_distance(float a, float b) {
  int32_t bits1 = reinterpret_cast<int32_t const &>(a);
  int32_t bits2 = reinterpret_cast<int32_t const &>(b);

  //if float a is negative build two's complement
  if(bits1 < 0)
    bits1 = 0x80000000 - bits1;

  //if float b is negative build two's complement
  if(bits2 < 0)
    bits2 = 0x80000000 - bits2;

  return std::abs(bits1 - bits2);
}

uint64_t ulp_distance(double a, double b) {
  int64_t bits1 = reinterpret_cast<int64_t const &>(a);
  int64_t bits2 = reinterpret_cast<int64_t const &>(b);

  if(bits1 < 0)
    bits1 = 0x8000000000000000l - bits1;

  if(bits2 < 0)
    bits2 = 0x8000000000000000l - bits2;

  return std::abs(bits1 - bits2);
}

using namespace std;

int main() {

  double a = 0;
  cout << "start float iteration...\n";
  while(!is_inf(a))
    a = next(a);
  cout << "reached inf\n";

  return EXIT_SUCCESS;
}
