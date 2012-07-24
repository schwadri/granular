#ifndef _core_math_linalg_matrix_matrix_hpp_
#define _core_math_linalg_matrix_matrix_hpp_

/** \file basic_matrix.hpp Contains the template class basic_matrix, that represents a mathematical basic_matrix.
 *  \author Adrian Schweizer
 *  \created  $Di 08 Jul 09:17:59 pm CEST 2008 schwadri@PwnererMachine.local$
 *  \modified $Mi 10 Sep 05:13:53 pm CEST 2008 schwadri@guest-docking-hg-1-154.ethz.ch$
 */

#include <algorithm>
#include <cmath>

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/iteration/local.hpp>

#include "fwd.hpp"
#include "init.hpp"
#include "layout.hpp"
#include "basic_matrix.hpp"
#include "util.hpp"

#ifndef CORE_MATH_LIMIT_MATRIX_DIMENSION
    #define CORE_MATH_LIMIT_MATRIX_DIMENSION 4
#endif

namespace core {

    namespace math {

        namespace linalg {

            #define GENERATE_TEMPLATE_PARAM_INNER(z, n, m) \
                typename BOOST_PP_CAT(T, BOOST_PP_CAT(m, n))

            #define GENERATE_TEMPLATE_PARAM(z, m, n) \
                BOOST_PP_ENUM(n, GENERATE_TEMPLATE_PARAM_INNER, m )

            #define GENERATE_ARG_INNER(z, n, m) \
                BOOST_PP_CAT(T, BOOST_PP_CAT(m, n)) BOOST_PP_CAT(value, BOOST_PP_CAT(m, n))

            #define GENERATE_ARG(z, m, n) \
                BOOST_PP_ENUM(n, GENERATE_ARG_INNER, m)

            #define GENERATE_ARG_ASSIGN_INNER(z, n, m) \
                (*this)(m, n) = static_cast<T>(BOOST_PP_CAT(value, BOOST_PP_CAT(m, n)));

            #define GENERATE_ARG_ASSIGN(z, m, n)  \
                BOOST_PP_REPEAT(n, GENERATE_ARG_ASSIGN_INNER, m)

            #define CORE_GENERATE_OBJECT_TEMPLATE_SPEC(m) \
                BOOST_PP_REPEAT(                                        \
                    BOOST_PP_SUB(CORE_MATH_LIMIT_MATRIX_DIMENSION, 1),  \
                    CORE_GENERATE_TEMPLATE_SPEC_ADDER,                  \
                    m                                                   \
                )

            #define CORE_GENERATE_TEMPLATE_SPEC_ADDER(z, n, m) \
                CORE_GENERATE_TEMPLATE_SPEC(m, BOOST_PP_ADD(n, 2))

            #define CORE_GENERATE_TEMPLATE_SPEC(m, n)                                                       \
                template <typename T, template <int, int> class Layout>                                     \
                    class matrix<T, m, n, Layout>                                                           \
                    :   public basic_matrix<T, m, n, Layout, matrix<T, m, n, Layout> >                      \
                    {                                                                                       \
                        typedef basic_matrix<T, m, n, Layout, matrix<T, m, n, Layout> > base_type;          \
                        typedef matrix                                                  type;               \
                    public:                                                                                 \
                                                                                                            \
                        template <int M, int N>                                                             \
                            struct rebind                                                                   \
                            {                                                                               \
                                typedef matrix<T, M, N, Layout> type;                                       \
                            };                                                                              \
                                                                                                            \
                        matrix() { }                                                                        \
                                                                                                            \
                        matrix(T const & value)                                                             \
                        :   base_type(value)                                                                \
                        { }                                                                                 \
                                                                                                            \
                        template <typename S>                                                               \
                            matrix(identity_matrix<S> const & im)                                           \
                            {                                                                               \
                                for(int i = 0; i < m; ++i)                                                  \
                                    for(int j = 0; j < n; ++j)                                              \
                                        if(i == j)                                                          \
                                            (*this)(i, j) = 1.0;                                            \
                                        else                                                                \
                                            (*this)(i, j) = 0.0;                                            \
                                                                                                            \
                            }                                                                               \
                                                                                                            \
                        template <typename S>                                                               \
                            matrix(scalar_diag_matrix<S> const & sm)                                        \
                            {                                                                               \
                                for(int i = 0; i < m; ++i)                                                  \
                                    for(int j = 0; j < n; ++j)                                              \
                                        if(i == j)                                                          \
                                            (*this)(i, j) = sm.value();                                     \
                                        else                                                                \
                                            (*this)(i, j) = 0.0;                                            \
                                                                                                            \
                            }                                                                               \
                                                                                                            \
                        template <BOOST_PP_ENUM(m, GENERATE_TEMPLATE_PARAM, n)>                             \
                            matrix(BOOST_PP_ENUM(m, GENERATE_ARG, n))                                       \
                            {                                                                               \
                                BOOST_PP_REPEAT(m, GENERATE_ARG_ASSIGN, n)                                  \
                            }                                                                               \
                                                                                                            \
                        typedef matrix_init<                                                                \
                            matrix,                                                                         \
                            1 < n ? 0 : 1,                                                                  \
                            1 < n ? 1 : 0                                                                   \
                        > next_init;                                                                        \
                                                                                                            \
                        template <typename Scalar>                                                          \
                            next_init operator = (Scalar const & value)                                     \
                            {                                                                               \
                                (*this)(0, 0) = value;                                                      \
                                                                                                            \
                                return next_init(*this);                                                    \
                            }                                                                               \
                                                                                                            \
                        template <typename S>                                                               \
                            type & operator=(zero_matrix<S> const & zm)                                     \
                            {                                                                               \
                                using std::fill;                                                            \
                                fill(this->begin(), this->end(), 0.0);                                      \
                                return *this;                                                               \
                            }                                                                               \
                                                                                                            \
                        template <typename S>                                                               \
                            type & operator=(identity_matrix<S> const & im)                                 \
                            {                                                                               \
                                for(int i = 0; i < m; ++i)                                                  \
                                    for(int j = 0; j < n; ++j)                                              \
                                        if(i == j)                                                          \
                                            (*this)(i, j) = 1.0;                                            \
                                        else                                                                \
                                            (*this)(i, j) = 0.0;                                            \
                                                                                                            \
                                return *this;                                                               \
                            }                                                                               \
                                                                                                            \
                        template <typename S>                                                               \
                            type & operator=(scalar_diag_matrix<S> const & sm)                              \
                            {                                                                               \
                                for(int i = 0; i < m; ++i)                                                  \
                                    for(int j = 0; j < n; ++j)                                              \
                                        if(i == j)                                                          \
                                            (*this)(i, j) = sm.value();                                     \
                                        else                                                                \
                                            (*this)(i, j) = 0.0;                                            \
                                                                                                            \
                                return *this;                                                               \
                            }                                                                               \
                    };

            #define BOOST_PP_LOCAL_MACRO   CORE_GENERATE_OBJECT_TEMPLATE_SPEC
            #define BOOST_PP_LOCAL_LIMITS (2, CORE_MATH_LIMIT_MATRIX_DIMENSION)
            //DO NOT REMOVE WHITESPACE   ^
            #include BOOST_PP_LOCAL_ITERATE()


            #undef CORE_GENERATE_OBJECT_BASE_SPEC
            #undef CORE_GENERATE_OBJECT_TEMPLATE_SPEC
            #undef CORE_GENERATE_OBJECT_TEMPLATE_SPEC2

            #undef TEMPLATE_PARAM_ADDER
            #undef TEMPLATE_PARAM_INNER
            #undef PARAM_ARG_INNER
            #undef PARAM_ARG2_INNER
            #undef GENERATE_PARAM_ARG_ADDER
            #undef GENERATE_PARAM_ARG2_ADDER

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_matrix_matrix_hpp_
