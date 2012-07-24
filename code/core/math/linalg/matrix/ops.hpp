#ifndef _core_math_linalg_matrix_ops_hpp_
#define _core_math_linalg_matrix_ops_hpp_

/** \file ops.hpp
 *  \author Adrian Schweizer
 *  \created  $Di 08 Jul 03:47:20 pm CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 *  \modified $Di 29 Jul 01:10:14 pm CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 */

#include <algorithm>
#include <cmath>

#include "fwd.hpp"
#include "../vector/fwd.hpp"

namespace core {

    namespace math {

        namespace linalg {

            template <typename T, template <int, int> class Layout>
                matrix<T, 3, 3, Layout> inverse(matrix<T, 3, 3, Layout> const & m)
                {
                    T m00 = (m[1][1]*m[2][2]-m[2][1]*m[1][2]),
                      m10 = (m[0][1]*m[2][2]-m[2][1]*m[0][2]),
                      m20 = (m[0][1]*m[1][2]-m[1][1]*m[0][2]);

                    T det =   m[0][0]*m00
                            - m[1][0]*m10
                            + m[2][0]*m20;

                    T inv_det = 1.0 / det;

                    return matrix<T, 3, 3, Layout>(
                                                       m00 * inv_det,                               -m10 * inv_det,                                m20 * inv_det,
                        -(m[1][0]*m[2][2]-m[2][0]*m[1][2]) * inv_det,  (m[0][0]*m[2][2]-m[2][0]*m[0][2]) * inv_det, -(m[0][0]*m[1][2]-m[1][0]*m[0][2]) * inv_det,
                         (m[1][0]*m[2][1]-m[2][0]*m[1][1]) * inv_det, -(m[0][0]*m[2][1]-m[2][0]*m[0][1]) * inv_det,  (m[0][0]*m[1][1]-m[1][0]*m[0][1]) * inv_det
                    );
                }


            template <typename T, template <int, int> class Layout>
                matrix<T, 2, 2, Layout> inverse(matrix<T, 2, 2, Layout> const & m)
                {
                    T det =   m[0][0]*m[1][1] - m[1][0]*m[0][1];
                    T inv_det = 1.0 / det;

                    return matrix<T, 2, 2, Layout>(
                         m[1][1] * inv_det, -m[0][1] * inv_det,
                        -m[1][0] * inv_det,  m[0][0] * inv_det
                    );
                }

            //fast inversion for homogeneous 4x4 matrices representing only rotations and translations
            template <typename T, template <int, int> class Layout>
                matrix<T, 4, 4, Layout> hom_inverse(matrix<T, 4, 4, Layout> const & m)
                {
                    using std::swap;
                    matrix<T, 4, 4, Layout> r(m);
                    swap(r[0][1], r[1][0]);
                    swap(r[0][2], r[2][0]);
                    swap(r[1][2], r[2][1]);

                    T   t1  =   r[0][0]*r[0][3] + r[0][1]*r[1][3] + r[0][2]*r[2][3],
                        t2  =   r[1][0]*r[0][3] + r[1][1]*r[1][3] + r[1][2]*r[2][3],
                        t3  =   r[2][0]*r[0][3] + r[2][1]*r[1][3] + r[3][2]*r[2][3];
                    r[0][3]  =   -t1;
                    r[1][3]  =   -t2;
                    r[2][3]  =   -t3;
                    r[3][0]  =   0;
                    r[3][1]  =   0;
                    r[3][2]  =   0;
                    r[3][3]  =   1.0;

                    return r;
                }

            template <typename T, int M, int N, template <int, int> class Layout>
                matrix<T, N, M, Layout> transpose(matrix<T, M, N, Layout> const & m)
                {
                    matrix<T, N, M, Layout> r;

                    for(int i = 0; i < M; ++i)
                        for(int j = 0; j < N; ++j)
                            r(j, i) = m(i, j);
                    return r;
                }

            //fast inversion for orthogonal matrices
            template <typename T, int M, template <int, int> class Layout>
                matrix<T, M, M, Layout> ortho_inverse(matrix<T, M, M, Layout> const & m)
                {
                    return transpose(m);
                }

            template <typename T, int M, int N, int O, template <int, int> class Layout>
                matrix<T, M, O, Layout> operator * (matrix<T, M, N, Layout> const & a, matrix<T, N, O, Layout> const & b)
                {
                    matrix<T, M, O, Layout> r;

                    for(int i = 0; i < M; ++i)
                        for(int j = 0; j < O; ++j)
                        {
                            r[i][j] = 0.0;

                            for(int k = 0; k < N; ++k)
                                r[i][j] += a[i][k] * b[k][j];
                        }
                    return r;
                }

            //FIXME: properly handle binary promotion of types
            //for now vector and matrix have to have the same element type
            template <typename T, int M, int N, template <int, int> class Layout>//, typename T1>
                vector<T, M> operator * (matrix<T, M, N, Layout> const & m, vector<T, N> const & v)
                {
                    vector<T, M> result;

                    for(int i = 0; i < M; ++i)
                    {
                        result[i] = 0.0;

                        for(int j = 0; j < N; ++j)
                            result[i] += m(i, j) * v[j];
                    }
                    return result;
                }

            template <typename T, int M, int N, template <int, int> class Layout, typename Vector, typename Range>
                vector<T, M> operator * (matrix<T, M, N, Layout> const & m, vector_static_range_proxy<Vector, Range> const & v)
                {
                    vector<T, M> result;

                    for(int i = 0; i < M; ++i)
                    {
                        result[i] = 0.0;

                        for(int j = 0; j < N; ++j)
                            result[j] += m(i, j) * v[j];
                    }
                    return result;
                }

            template <typename T, int M, int N, template <int, int> class Layout>
                matrix<T, M, N, Layout> operator + (matrix<T, M, N, Layout> const & m1, matrix<T, M, N, Layout> const & m2)
                {
                    matrix<T, M, N, Layout> r(m1);
                    r += m2;
                    return r;
                }

            template <typename T, int M, int N, template <int, int> class Layout>
                matrix<T, M, N, Layout> operator - (matrix<T, M, N, Layout> const & m1, matrix<T, M, N, Layout> const & m2)
                {
                    matrix<T, M, N, Layout> r(m1);
                    r -= m2;
                    return r;
                }

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_matrix_ops_hpp_
