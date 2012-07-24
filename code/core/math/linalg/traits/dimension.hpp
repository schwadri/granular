#ifndef _core_math_linalg_traits_dimension_hpp_
#define _core_math_linalg_traits_dimension_hpp_

/**\file dimension.hpp Contains a metafunction to query the size of a vector type.
 */

#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/size_t.hpp>

#include "detail/unsupported_type.hpp"
#include "is_scalar.hpp"

namespace core {

    namespace math {

        namespace linalg {

            namespace mpl = boost::mpl;

            //if T is a not supported vector type there is no template specialization
            //and template instanciation results in a compiler error
            /** \brief metafunction to query the dimension of a vector type
             *
             *  Usage: \code dimension<VectorType>::value \endcode
             *  Depending on whether the vector is of static (compile time defined) dimension
             *  or dynamic vector \a value is N or 0. If The vector has a dynamic storage, its 
             *  dimension can be queried at runtime like this:
             *  \code dimension<VectorType>::get(v) \endcode
             */
            template <typename T>
                struct dimension
                :   mpl::eval_if<
                        is_scalar<T>,
                        mpl::size_t<1>,
                        unsupported_type<T>
                    >::type
                {
                    static std::size_t get(const T& t)
                    {
                        return 0;
                    }
                };

        } // namespace linalg
    } // namespace math
} // namespace core

#include "detail/vector_traits_spec.hpp"

#endif // _core_math_linalg_traits_dimension_hpp_
