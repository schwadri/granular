#ifndef _linalg_hpp_
#define _linalg_hpp_

template <int I>
struct int_ {
  typedef int_ type;
  static int const value = I;
};

#include <core/math/linalg/vector/vector.hpp>
#include <core/math/linalg/vector/ops.hpp>
#include <core/math/linalg/vector/io.hpp>

#define CORE_MATH_LIMIT_MATRIX_DIMENSION 6

#include <core/math/linalg/matrix/matrix.hpp>
#include <core/math/linalg/matrix/ops.hpp>
#include <core/math/linalg/matrix/util.hpp>
#include <core/math/linalg/matrix/io.hpp>
#include <core/math/linalg/utility/quaternion.hpp>
#include <core/math/clamp.hpp>

namespace lin = core::math::linalg;
using core::math::clamp;
using lin::to_matrix;

typedef lin::vector<real, 2>        vector2r;
typedef lin::vector<real, 3>        vector3r;
typedef lin::vector<real, 4>        vector4r;
typedef lin::vector<real, 6>        vector6r;
typedef lin::zero_vector<real>      zero_vector;
typedef lin::scalar_vector<real>    scalar_vector;

#define MATRIX_LAYOUT lin::layout::row_major

template <typename T, int N>
    struct mat
    {
        typedef lin::matrix<
            T,
            N, N,
            MATRIX_LAYOUT
        > type;
    };


typedef mat<real, 3>::type                          matrix3r;
typedef mat<real, 4>::type                          matrix4r;
typedef mat<real, 2>::type                          matrix2x2r;
typedef mat<real, 3>::type                          matrix3x3r;
typedef mat<real, 6>::type                          matrix6x6r;
typedef lin::matrix<real, 6, 3, MATRIX_LAYOUT>      matrix6x3r;
typedef lin::matrix<real, 2, 3, MATRIX_LAYOUT>      matrix2x3r;
typedef lin::util<real, 4, lin::layout::row_major>  matrix4_gen;
typedef lin::util<real, 3, lin::layout::row_major>  matrix3_gen;
typedef lin::zero_matrix<real>                      zero_matrix;
typedef lin::identity_matrix<real>                  identity_matrix;
typedef lin::scalar_diag_matrix<real>               scalar_diag_matrix;

#endif // _linalg_hpp_

