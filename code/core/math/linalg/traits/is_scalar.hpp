#ifndef _core_math_linalg_traits_is_scalar_hpp_
#define _core_math_linalg_traits_is_scalar_hpp_

#include <boost/mpl/bool.hpp>

namespace core {

    namespace math {

        namespace linalg {

        namespace mpl = boost::mpl;

        //mathematical type traits

        /** \brief template metafunction to query whether a type is a scalar
         *
         *  Example Usage:
         *  \code 
         *  typedef math::is_scalar<int>::type                  a;  // returns true_
         *  typedef math::is_scalar<std::vector<int> >::type    b;  // returns false_
         *  typedef math::is_scalar<double[10]>::type           c;  // false_
         *  bool really = typedef math::is_scalar<Point>::value;    // false
         *  \endcode
         *  this is true for types that implement all the operators a normal real 
         *  has (+,-,/,*,+=,-=,/=,*=, ==,<,>,<=,>=)
         */
        template<class T>
            struct is_scalar : mpl::false_
            { };

        } // namespace linalg
    } // namespace math
} // namespace core

#include "detail/scalar_vector_spec.hpp"

#endif // _core_math_linalg_traits_is_scalar_hpp_
