#ifndef _core_math_linalg_traits_is_vector_hpp_
#define _core_math_linalg_traits_is_vector_hpp_

#include <boost/mpl/bool.hpp>

namespace core {

    namespace math {

        namespace linalg {

            namespace mpl = boost::mpl;

            /** \brief template metafunction to check whether a type is a vector or not
             *
             *  Example Usage:
             *  \code
             *  typedef math::is_vector<int>::type                  a;  // returns false_
             *  typedef math::is_vector<std::vector<int> >::type    b;  // returns true_
             *  bool really = typedef math::is_vector<std::vector<int> >::value;  // returns true
             *  typedef math::is_vector<double[10]>::type           c;  // true_
             *  bool really1 = typedef math::is_vector<Point>::value;   // true
             *  \endcode
             *  In general this metafunction evaluates to \a true_ if the type in question implements operator[](int index)
             *  and the underlying field type is a scalar.
             */
            template<class T>
                struct is_vector : mpl::false_
                { };

        } // namespace linalg
    } // namespace math
} // namespace core

#include "detail/scalar_vector_spec.hpp"

#endif // _core_math_linalg_traits_is_vector_hpp_
