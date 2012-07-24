#ifndef _core_math_linalg_traits_spec_hpp_
#define _core_math_linalg_traits_spec_hpp_

/* \file scalar_vector_spec.hpp
 *  \author Adrian Schweizer
 *  \created  $Di 08 Jul 11:04:03 am CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 *  \modified $Di 08 Jul 11:04:03 am CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 */

#include "fwd_decls.hpp"

namespace core {

    namespace math {

        namespace linalg {

            //fwd decl
            template <typename T>
                struct is_vector;

            template <typename T>
                struct is_scalar;

            //specialization for standard types
            #define CORE_MATH_IS_VECTOR(T)              \
            template <>                                 \
                struct is_vector<T> : mpl::true_        \
                { };

            #define CORE_MATH_IS_NOT_VECTOR(T)          \
            template <>                                 \
                struct is_vector<T> : mpl::false_       \
            { };

            #define CORE_MATH_IS_SCALAR(T)              \
            template <>                                 \
                struct is_scalar<T> : mpl::true_        \
                { };

            #define CORE_MATH_INTEGRAL_SPECIALIZATION(T)\
                CORE_MATH_IS_SCALAR(T)                  \
                CORE_MATH_IS_NOT_VECTOR(T)

            //specializations for builtin types
            CORE_MATH_INTEGRAL_SPECIALIZATION(char)
            CORE_MATH_INTEGRAL_SPECIALIZATION(unsigned char)
            CORE_MATH_INTEGRAL_SPECIALIZATION(short)
            CORE_MATH_INTEGRAL_SPECIALIZATION(unsigned short)
            CORE_MATH_INTEGRAL_SPECIALIZATION(int)
            CORE_MATH_INTEGRAL_SPECIALIZATION(unsigned int)
            CORE_MATH_INTEGRAL_SPECIALIZATION(float)
            CORE_MATH_INTEGRAL_SPECIALIZATION(double)
            CORE_MATH_INTEGRAL_SPECIALIZATION(bool)


            #undef CORE_MATH_INTEGRAL_SPECIALIZATION
            #undef CORE_MATH_IS_VECTOR
            #undef CORE_MATH_IS_NOT_VECTOR
            #undef CORE_MATH_IS_SCALAR

            //specialization for known classes
            template <typename T, std::size_t N>
                struct is_vector<boost::array<T, N> >
                :   mpl::if_<
                        is_scalar<T>,
                        mpl::true_,
                        mpl::false_
                    >::type
                { };

            template <typename T, std::size_t N>
                struct is_vector<T[N]>
                :   mpl::if_<
                        is_scalar<T>,
                        mpl::true_,
                        mpl::false_
                    >::type
                { };

            template <typename T, typename A>
                struct is_vector<std::vector<T, A> >
                :   mpl::if_<
                        is_scalar<T>,
                        mpl::true_,
                        mpl::false_
                    >::type
                { };

            template <typename T>
                struct is_vector<std::valarray<T> >
                :   mpl::if_<
                        is_scalar<T>,
                        mpl::true_,
                        mpl::false_
                    >::type
                { };

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_traits_spec_hpp_
