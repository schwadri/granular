#ifndef _core_math_linalg_vector_ops_hpp_
#define _core_math_linalg_vector_ops_hpp_

/** \file ops.hpp
 *  \author Adrian Schweizer
 *  \created  $Di 08 Jul 12:21:34 am CEST 2008 schwadri@PwnererMachine.local$
 *  \modified $Di 08 Jul 01:38:10 pm CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 */

#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>
#include <iosfwd>

#include <boost/preprocessor/control/expr_if.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/preprocessor/comparison/less.hpp>

#include "norm.hpp"
#include "transform.hpp"

namespace core {

    namespace math {

        namespace linalg {

            template <typename T, int D>
                vector<T, D>  operator + (vector<T,D> const & v0, vector<T, D> const & v1)
                {
                    vector<T, D> temp(v0);
                    temp += v1;
                    return temp;
                }

            template <typename T, int D>
                vector<T, D> operator - (vector<T, D> const & v0, vector<T, D> const & v1)
                {
                    vector<T, D> temp(v0);
                    temp -= v1;
                    return temp;
                }

            template <typename T, int D, typename S>
                vector<T,D>  operator + (vector<T, D> const & v, S const & s)
                {
                    vector<T, D> temp(v);
                    temp += s;
                    return temp;
                }

            template <typename T, int D, typename S>
                vector<T,D> operator - (vector<T, D> const & v, S const & s)
                {
                    vector<T, D> temp(v);
                    temp -= s;
                    return temp;
                }

            template <typename T, int D, typename S>
                vector<T, D>  operator + (S const & s, vector<T, D> const & v)
                {
                    vector<T, D> temp(v);
                    temp += s;
                    return temp;
                }

            template <typename T, int D, typename S>
                vector<T, D> operator - (S const & s, vector<T, D> const & v)
                {
                    vector<T, D> temp(s);
                    temp -= v;
                    return temp;
                }
            template <typename T, int D>
                vector<T, D>  operator * (vector<T, D> const & v0, vector<T, D> const & v1)
                {
                    vector<T, D> temp(v0);
                    temp *= v1;
                    return temp;
                }
            
            template <typename T, int D, typename S>
                vector<T, D>  operator * (vector<T, D> const & v, S const & s)
                {
                    vector<T, D> temp(v);
                    temp *= s;
                    return temp;
                }

            template <typename T, int D, typename S>
                vector<T, D>  operator / (vector<T, D> const & v, S const & s)
                {
                    vector<T, D> temp(v);
                    temp /= s;
                    return temp;
                }

            template <typename T, int D, typename S>
                vector<T, D>  operator * (S const & s, vector<T, D> const & v)
                {
                    vector<T, D> temp(v);
                    temp *= s;
                    return temp;
                }

            namespace aux {

                template<class T>
                    struct divides_value
                    {
                        divides_value(const T& t)
                        :   value(t)
                        { }

                        T operator()(const T& v0)
                        {
                            return value/v0;
                        }

                        T value;
                    };
            }

            template <typename T, int D, class S>
                vector<T, D>  operator / (S const & s, vector<T, D> const & v)
                {
                    vector<T, D> temp(v);
                    transform(v, temp, aux::divides_value<T>(s));
                    return temp;
                }

            template <typename T>
                vector<T, 3> cross(vector<T, 3> const & v1, vector<T, 3> const & v2)
                {
                    return vector<T, 3>(
                        v1[1]*v2[2] - v1[2]*v2[1],
                        v1[2]*v2[0] - v1[0]*v2[2],
                        v1[0]*v2[1] - v1[1]*v2[0]
                    );
                }

            template <typename T>
                T cross(vector<T, 2> const & v1, vector<T, 2> const & v2)
                {
                    return v1[0]*v2[1] - v1[0]*v2[1];
                }

            //scalar product shortcut
            template <typename T, int D>
                T dot(vector<T, D> const & v0, vector<T, D> const & v1)
                {
                    using std::inner_product;
                    using std::plus;
                    using std::multiplies;
                    return inner_product(
                        v0.begin(),
                        v0.end(),
                        v1.begin(),
                        T(0),
                        plus<T>(),
                        multiplies<T>()
                    );
                }

            //norm_2 shortcut
            template <typename T, int D>
                T norm(vector<T, D> const & v)
                {
                    return norm_2(v);
                }

            template <typename T, int Size>
                T norm_sqr(vector<T, Size> const & v)
                {
                    return norm_2_sqr(v);
                }

            //relational operators for vector<T,n> for n from 2-GENERATE_MAX_PARAMS

            #define GENERATE_PARAM_OP_ADDER(z,n,OP)  && lhs[BOOST_PP_ADD(n,1)] OP rhs[BOOST_PP_ADD(n,1)]

            #define CORE_GENERATE_OBJECT_FUNCS(n)                                           \
                    BOOST_PP_EXPR_IF(   BOOST_PP_LESS(n,CORE_MATH_LIMIT_VECTOR_DIMENSION),  \
                                        GENERATE_OBJECT_FUNC(n,==)                          \
                                        GENERATE_OBJECT_FUNC(n,<)                           \
                                        GENERATE_OBJECT_FUNC(n,>)                           \
                                        GENERATE_OBJECT_FUNC(n,<=)                          \
                                        GENERATE_OBJECT_FUNC(n,>=))

            #define GENERATE_OBJECT_FUNC(n,OP)                                                      \
                template <typename T>                                                               \
                    inline bool operator OP(                                                        \
                        vector<T,BOOST_PP_ADD(n,1)> const & lhs,                                    \
                        vector<T,BOOST_PP_ADD(n, 1)> const & rhs                                    \
                    )                                                                               \
                    {                                                                               \
                        return lhs[0] OP rhs[0] BOOST_PP_REPEAT(n, GENERATE_PARAM_OP_ADDER, OP);    \
                    }                                                                               \

            template <typename T, int D>
                inline bool operator!=(vector<T, D> const & lhs, vector<T, D> const & rhs)
                {
                    return !(lhs == rhs);
                }



            #define BOOST_PP_LOCAL_MACRO   CORE_GENERATE_OBJECT_FUNCS
            #define BOOST_PP_LOCAL_LIMITS (1, CORE_MATH_LIMIT_VECTOR_DIMENSION)
            //DO NOT REMOVE WHITESPACE   ^
            #include BOOST_PP_LOCAL_ITERATE()

            #undef GENERATE_PARAM_OP_ADDER
            #undef GENERATE_OBJECT_FUNC
            #undef CORE_GENERATE_OBJECT_FUNCS

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_vector_hpp_
