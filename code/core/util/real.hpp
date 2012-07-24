#ifndef _core_util_real_hpp_
#define _core_util_real_hpp_

#include <iostream>   //using cout, cerr
#include <limits>     //using numeric_limits
#include <cmath>
#include <tr1/cstdint>

using std::tr1::uint64_t;
using std::tr1::int64_t;
using std::tr1::uint32_t;
using std::tr1::int32_t;

inline bool is_nan(float f) {
  uint32_t const & b = reinterpret_cast<uint32_t const &>(f);
  return (b & 0x7f800000) == 0x7f800000 && (b & 0x007fffff);
}

inline bool is_nan(double f) {
  uint64_t const & b = reinterpret_cast<uint64_t const &>(f);
  return (b & 0x7ff0000000000000l) == 0x7ff0000000000000l && (b & 0x000fffffffffffffl);
}

inline bool is_snan(double f) {
  uint64_t const & b = reinterpret_cast<uint64_t const &>(f);
  return (b & 0x7ff0000000000000l) == 0x7ff0000000000000l && (b & 0x000fffffffffffffl) && !(b & 0x0008000000000000l);
}

inline bool is_qnan(double f) {
  uint64_t const & b = reinterpret_cast<uint64_t const &>(f);
  return (b & 0x7ff0000000000000l) == 0x7ff0000000000000l && (b & 0x000fffffffffffffl) && (b & 0x0008000000000000l);
}

inline bool is_snan(float f) {
  uint32_t const & b = reinterpret_cast<uint32_t const &>(f);
  return (b & 0x7f800000) == 0x7f800000 && (b & 0x007fffff) && !(b & 0x00400000);
}

inline bool is_qnan(float f) {
  uint32_t const & b = reinterpret_cast<uint32_t const &>(f);
  return (b & 0x7f800000) == 0x7f800000 && (b & 0x007fffff) && (b & 0x00400000);
}

inline bool is_inf(double f) {
  uint64_t const & b = reinterpret_cast<uint64_t const &>(f);
  return (b & 0x7fffffffffffffffl) == 0x7ff0000000000000l;
}

inline bool is_inf(float f) {
  uint32_t const & b = reinterpret_cast<uint32_t const &>(f);
  return (b & 0x7fffffff) == 0x7f800000;
}

inline uint32_t payload(float f) {
  uint32_t const & b = reinterpret_cast<uint32_t const &>(f);
  return (b & 0x003fffff);
}

inline uint64_t payload(double f) {
  uint64_t const & b = reinterpret_cast<uint64_t const &>(f);
  return (b & 0x0007ffffffffffffl);
}

inline uint64_t mantissa(double f) {
  uint64_t const & b = reinterpret_cast<uint64_t const &>(f);
  return (b & 0x000fffffffffffffl);
}

inline int64_t exponent(double f) {
  int64_t const & b = reinterpret_cast<int64_t const &>(f);
  return ((b & 0x7ff0000000000000l) >> 52) - 1023;
}

inline uint32_t mantissa(float v) {
  uint32_t bits = reinterpret_cast<uint32_t const &>(v);
  return bits & 0x007fffff;
}

inline int32_t exponent(float v) {
  int32_t bits = reinterpret_cast<int32_t const &>(v);
  return ((bits & 0x7f800000) >> 23) - 127;
}

inline double sign(double v) {
  uint64_t bits = reinterpret_cast<uint64_t const &>(v);
  return bits & 0x8000000000000000l ? -1.0 : 1.0;
}

inline float sign(float v) {
  uint32_t bits = reinterpret_cast<uint32_t const &>(v);
  return bits & 0x80000000 ? -1.0 : 1.0;
}


inline double next(double prev) {
  int64_t bits = reinterpret_cast<int64_t const &>(prev);

  if(bits < 0) {
    --bits;
    if(bits > 0)
      bits = 0;
  } else
    ++bits;
  return reinterpret_cast<double const &>(bits);
}

inline double prev(double v) {
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

inline float next(float prev) {
  int32_t bits = reinterpret_cast<int32_t const &>(prev);

  if(bits < 0) {
    --bits;
    if(bits > 0)
      bits = 0;
  } else
    ++bits;

  return reinterpret_cast<float const &>(bits);
}

inline float prev(float v) {
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

inline uint32_t ulp_distance(float a, float b) {
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

inline uint64_t ulp_distance(double a, double b) {
  int64_t bits1 = reinterpret_cast<int64_t const &>(a);
  int64_t bits2 = reinterpret_cast<int64_t const &>(b);

  if(bits1 < 0)
    bits1 = 0x8000000000000000l - bits1;

  if(bits2 < 0)
    bits2 = 0x8000000000000000l - bits2;

  return std::abs(bits1 - bits2);
}

/**
 * set the fpu control word to change rounding and precision behaviour
 */
inline void set_fcw(uint16_t cw) {
  asm(
    "fldcw %0\n"
    : : "m" (cw)
  );
}

/**
 * get the current fpu control word content
 */
inline uint16_t get_fcw() {
  uint16_t r;
  asm(
    "fstcw %0\n"
    : "=m" (r) : "m" (r) : "memory"
  );
  return r;
}

void show_rounding_mode() {
  uint16_t cw = get_fcw();
  switch((cw >> 10) & 0x3) {
    case 0:
      std::cout << "round to nearest, break ties with even\n";
      break;
    case 1:
      std::cout << "round to -infinity\n";
      break;
    case 2:
      std::cout << "round to +infinity\n";
      break;
    case 3:
      std::cout << "truncate (round to zero)\n";
      break;
  }
}

void show_precision_mode() {
  uint16_t cw = get_fcw();
  switch((cw >> 8) & 0x3) {
    case 0:
      std::cout << "24 bits mantissa\n";
      break;
    case 1:
      std::cout << "reserved\n";
      break;
    case 2:
      std::cout << "53 bits mantissa\n";
      break;
    case 3:
      std::cout << "64 bits mantissa\n";
      break;
  }
}

#endif // _core_util_real_hpp_
