#ifndef _core_math_linalg_vector_vector_hpp_
#define _core_math_linalg_vector_vector_hpp_

/** \file vector.hpp Contains the template class \a vector, that represents a linear algebraic vector.
 *  \author Adrian Schweizer
 *  \created  $Do 23 Aug 10:32:12 pm CEST 2007 schwadri@SchwadriComp.local$
 *  \modified $Di 29 Jul 01:31:12 pm CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 */

#include <algorithm>

#include <boost/preprocessor/control/expr_if.hpp>
#include <boost/preprocessor/empty.hpp>
#include <boost/preprocessor/identity.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/preprocessor/facilities/expand.hpp>
#include <boost/preprocessor/comparison/greater_equal.hpp>
#include <boost/preprocessor/comparison/less.hpp>
#include <boost/mpl/bool.hpp>

#ifndef CORE_MATH_LIMIT_VECTOR_DIMENSION
    #define CORE_MATH_LIMIT_VECTOR_DIMENSION   6
#endif

#include "fwd.hpp"
#include "basic_vector.hpp"
#include "../traits/field.hpp"
#include "../traits/is_vector.hpp"
#include "../traits/dimension.hpp"
#include "init.hpp"
#include "util.hpp"

namespace core {

    namespace math {

        namespace linalg {

            //fwd decl
            template <typename Vector, typename Range>
                struct vector_static_range_proxy;

            namespace mpl = boost::mpl;

            #define GENERATE_FORMAL_PARAM_DECL(z, n, arg) \
                T2 BOOST_PP_CAT(arg, n)

            #define INIT_VALUE_FROM_ARG(z, n, arg) \
                this->elems[n] = static_cast<T>(BOOST_PP_CAT(arg, n));

            #define CORE_GENERATE_OBJECT_FUNC(n)                                                \
                template <typename T>                                                           \
                    class vector<T, n>                                                          \
                    :   public basic_vector<T, n, vector<T, n> >                                \
                    {                                                                           \
                        typedef basic_vector<T, n, vector<T, n> >   base_type;                  \
                        typedef vector<T, n>                        type;                       \
                                                                                                \
                    public:                                                                     \
                                                                                                \
                        template <int S>                                                        \
                            struct rebind                                                       \
                            {                                                                   \
                                typedef vector<T, S> type;                                      \
                            };                                                                  \
                                                                                                \
                        struct initializer                                                      \
                        {                                                                       \
                            T values[n];                                                        \
                        };                                                                      \
                                                                                                \
                        vector() { }                                                            \
                                                                                                \
                        explicit vector(T const & value)                                        \
                        :   base_type(value)                                                    \
                        { }                                                                     \
                                                                                                \
                        vector(initializer const & init)                                        \
                        {                                                                       \
                            using std::copy;                                                    \
                            copy(init.values, init.values + n, this->begin());                  \
                        }                                                                       \
                                                                                                \
                        template <typename S>                                                   \
                            vector(zero_vector<S> const & zv)                                   \
                            {                                                                   \
                                using std::fill;                                                \
                                fill(this->begin(), this->end(), 0.0);                          \
                            }                                                                   \
                                                                                                \
                        template <typename S>                                                   \
                            vector(scalar_vector<S> const & sv)                                 \
                            {                                                                   \
                                using std::fill;                                                \
                                fill(this->begin(), this->end(), sv.value());                   \
                            }                                                                   \
                                                                                                \
                        template <typename Vector, typename Range>                              \
                            vector(vector_static_range_proxy<Vector, Range> const & vr)         \
                            {                                                                   \
                                using std::copy;                                                \
                                copy(vr.begin(), vr.end(), this->begin());                      \
                            }                                                                   \
                                                                                                \
                                                                                                \
                        template <typename T2>                                                  \
                            vector(BOOST_PP_ENUM(n, GENERATE_FORMAL_PARAM_DECL, arg))           \
                            {                                                                   \
                                BOOST_PP_REPEAT(n, INIT_VALUE_FROM_ARG, arg)                    \
                            }                                                                   \
                                                                                                \
                        template <typename Scalar>                                              \
                            vector_init<type, 1> operator=(Scalar const & value)                \
                            {                                                                   \
                                this->elems[0] = value;                                         \
                                return vector_init<type, 1>(*this);                             \
                            }                                                                   \
                                                                                                \
                        template <typename S>                                                   \
                            type & operator=(zero_vector<S> const & zv)                         \
                            {                                                                   \
                                using std::fill;                                                \
                                fill(this->begin(), this->end(), 0.0);                          \
                                return *this;                                                   \
                            }                                                                   \
                                                                                                \
                        template <typename S>                                                   \
                            type & operator=(scalar_vector<S> const & sv)                       \
                            {                                                                   \
                                using std::fill;                                                \
                                fill(this->begin(), this->end(), sv.value());                   \
                                return *this;                                                   \
                            }                                                                   \
                                                                                                \
                    };

            #define BOOST_PP_LOCAL_MACRO   CORE_GENERATE_OBJECT_FUNC
            #define BOOST_PP_LOCAL_LIMITS (2, CORE_MATH_LIMIT_VECTOR_DIMENSION)
            //DO NOT REMOVE WHITESPACE   ^
            #include BOOST_PP_LOCAL_ITERATE()

            #undef GENERATE_PARAM_ARG_ADDER
            #undef GENERATE_PARAM_ARG2_ADDER
            #undef GENERATE_OBJECT_FUNCS

            //overload math and vector traits
            template <int D, typename T> struct is_vector<vector<T, D> > : mpl::true_ { };
            template <int D, typename T> struct dimension<vector<T, D> > : detail::static_dimension_impl<vector<T, D>, D> { };
            template <int D, typename T> struct field<vector<T, D> > { typedef T type; };

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_vector_vector_hpp_
