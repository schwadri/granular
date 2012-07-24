#ifndef _core_math_linalg_matrix_util_hpp_
#define _core_math_linalg_matrix_util_hpp_

/** \file util.hpp Contains the template class \a vector, that represents a linear algebraic vector.
 *  \author Adrian Schweizer
 *  \created  $Di 08 Jul 10:51:00 pm CEST 2008 schwadri@PwnererMachine.local$
 *  \modified $Di 29 Jul 09:48:40 am CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 */

#include <cmath>

#include "layout.hpp"
#include "fwd.hpp"
#include "../vector/fwd.hpp"

namespace core {

    namespace math {

        namespace linalg {

            template <typename T>
                struct zero_matrix
                {
                    T const & operator()(int i, int j) const
                    {
                        return 0;
                    }
                };

            template <typename T>
                struct identity_matrix
                {
                    T const & operator()(int i, int j) const
                    {
                        return (i == j) ? 1 : 0;
                    }
                };

            template <typename T>
                struct scalar_diag_matrix
                {
                    scalar_diag_matrix(T const & value)
                    :   m_value(value)
                    { }

                    T const & operator()(int i, int j) const
                    {
                        return (i == j) ? m_value : 0;
                    }

                    T const & value() const
                    {
                        return m_value;
                    }

                private:
                    T m_value;
                };

            //FIXME: helper for homogenous transformations. transform a 2 M vector with a 3x3 basic_matrix, or a 3 M vector with a 4x4 basic_matrix
            /*vector<value_type,M-1> operator *(vector<value_type,M-1>& v)
            {
                vector<value_type,M-1> result;
                for(int i1=0;i1<M;++i1)
                {
                    result[i1]=0;
                    for(int i2=0;i2<(M-1);++i2)
                        result[i1]+= m_data[i1+i2*M]*v[i2];
                    result[i1]+=m_data[i1+(M-1)*M];
                }
                return result;
            }
            */

            template <typename T, int M, template <int, int> class Layout = layout::row_major>
                struct util;

            template <typename T, template <int, int> class Layout>
                struct util<T, 4, Layout>
                {
                    typedef matrix<T, 4, 4, Layout> type;

                    static type identity()
                    {
                        return type(
                            1, 0, 0, 0,
                            0, 1, 0, 0,
                            0, 0, 1, 0,
                            0, 0, 0, 1
                        );
                    }

                    static type rotatex(T const & angle)
                    {
                        using std::cos;
                        using std::sin;
                        T c = static_cast<T>(cos(angle));
                        T s = static_cast<T>(sin(angle));

                        return type(
                            1, 0, 0, 0,
                            0, c,-s, 0,
                            0, s, c, 0,
                            0, 0, 0, 1
                        );
                    }

                    static type rotatey(T const & angle)
                    {
                        using std::cos;
                        using std::sin;
                        T c = static_cast<T>(cos(angle));
                        T s = static_cast<T>(sin(angle));

                        return type(
                            c, 0, s, 0,
                            0, 1, 0, 0,
                            -s, 0, c, 0,
                            0, 0, 0, 1
                        );
                    }

                    static type rotatez(T const & angle)
                    {
                        using std::cos;
                        using std::sin;
                        T c = static_cast<T>(cos(angle));
                        T s = static_cast<T>(sin(angle));

                        return type(
                            c,-s, 0, 0,
                            s, c, 0, 0,
                            0, 0, 1, 0,
                            0, 0, 0, 1
                        );
                    }

                    static type mirrorx()
                    {
                        return type(
                            -1, 0, 0, 0,
                            0, 1, 0, 0,
                            0, 0, 1, 0,
                            0, 0, 0, 1
                        );
                    }

                    static type mirrory()
                    {
                        return type(
                            1, 0, 0, 0,
                            0,-1, 0, 0,
                            0, 0, 1, 0,
                            0, 0, 0, 1
                        );
                    }

                    static type mirrorz()
                    {
                        return type(
                            1, 0, 0, 0,
                            0, 1, 0, 0,
                            0, 0,-1, 0,
                            0, 0, 0, 1
                        );
                    }

                    static type scale(vector<T, 3> const & scale)
                    {
                        return type(
                            scale[0],        0,        0, 0,
                                   0, scale[1],        0, 0,
                                   0,        0, scale[2], 0,
                                   0,        0,        0, 1
                        );
                    }

                    static type scale(T const & scale)
                    {
                        return type(
                            scale,     0,     0, 0,
                                0, scale,     0, 0,
                                0,     0, scale, 0,
                                0,     0,     0, 1
                        );
                    }

                    static type translate(vector<T, 3> const & v)
                    {
                        return type(
                            1, 0, 0, v[0],
                            0, 1, 0, v[1],
                            0, 0, 1, v[2],
                            0, 0, 0, 1
                        );
                    }

                    struct projection
                    {

                        static type frustum(vector<T, 2> const & bl, vector<T, 2> const & ur, T const & near, T const & far)
                        {
                            return type(
                                2.0 * near/(ur[0] - bl[0]), 0, (ur[0]+bl[0])/(ur[0]-bl[0]), 0,
                                0, 2.0 *near/(ur[1]-bl[1]), (ur[1]+bl[1])/(ur[1]-bl[1]),0,
                                0, 0, -(far+near)/(far-near), -2.0*far*near/(far-near),
                                0, 0, -1, 0
                            );
                        }

                        static type frustumsym(vector<T, 2> const & size, T const & near, T const & far)
                        {
                            return type(
                                2*near/(size[0]), 0, 0, 0,
                                0, 2*near/(size[1]), 0, 0,
                                0, 0, -(far+near)/(far-near), -2*far*near/(far-near),
                                0, 0, -1, 0
                            );
                        }

                        static type ortho(vector<T, 2> const & bl, vector<T, 2> const & ur, T const & near, T const & far)
                        {
                            return type(
                                2/(ur[0]-bl[0]), 0, 0, -(ur[0]+bl[0]) / (ur[0]-bl[0]),
                                0, 2/(ur[1]-bl[1]), 0, -(ur[1]+bl[1]) / (ur[1]-bl[1]),
                                0, 0, -2/(far-near), -2*(far+near)/(far-near),
                                0, 0, 0, 1
                            );
                        }

                        static type perspective(T const & fov, T const & aspect, T const & near, T const & far)
                        {
                            using std::tan;
                            return type(
                                1.0/(tan(fov * 0.5)*aspect), 0, 0, 0,
                                0, 1.0/(tan(fov * 0.5)), 0, 0,
                                0, 0, -(far+near)/(far-near), -2.0 *far*near/(far-near),
                                0, 0, -1, 0
                            );
                        }

                    };
                };

            template <typename T, template <int, int> class Layout>
                struct util<T, 3, Layout>
                {
                    typedef matrix<T, 3, 3, Layout> type;

                    static type identity()
                    {
                        return type(
                            1, 0, 0,
                            0, 1, 0,
                            0, 0, 1
                        );
                    }

                    static type rotatex(T const & angle)
                    {
                        using std::cos;
                        using std::sin;
                        T c = static_cast<T>(cos(angle));
                        T s = static_cast<T>(sin(angle));

                        return type(
                            1, 0, 0,
                            0, c,-s,
                            0, s, c
                        );
                    }

                    static type rotatey(T const & angle)
                    {
                        using std::cos;
                        using std::sin;
                        T c = static_cast<T>(cos(angle));
                        T s = static_cast<T>(sin(angle));

                        return type(
                            c, 0, s,
                            0, 1, 0,
                            -s, 0, c
                        );
                    }

                    static type rotatez(T const & angle)
                    {
                        using std::cos;
                        using std::sin;
                        T c = static_cast<T>(cos(angle));
                        T s = static_cast<T>(sin(angle));

                        return type(
                            c,-s, 0,
                            s, c, 0,
                            0, 0, 1
                        );
                    }

                    static type mirrorx()
                    {
                        return type(
                            -1, 0, 0,
                            0, 1, 0,
                            0, 0, 1
                        );
                    }

                    static type mirrory()
                    {
                        return type(
                            1, 0, 0,
                            0,-1, 0,
                            0, 0, 1
                        );
                    }

                    static type mirrorz()
                    {
                        return type(
                            1, 0, 0,
                            0, 1, 0,
                            0, 0,-1
                        );
                    }

                    static type scale(vector<T, 3> const & scale)
                    {
                        return type(
                            scale[0],        0,        0,
                                   0, scale[1],        0,
                                   0,        0, scale[2]
                        );
                    }

                    static type scale(T const & scale)
                    {
                        return type(
                            scale,     0,     0,
                                0, scale,     0,
                                0,     0, scale
                        );
                    }
                };
        } // namespace linalg
    } //namespace math
} // namespace core

#endif // _core_math_linalg_matrix_util_hpp_
