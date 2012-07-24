#ifndef _core_math_linalg_traits_field_hpp_
#define _core_math_linalg_traits_field_hpp_

/**\file field.hpp Contains metafunctions to query the size of a vector type
 *  and the scalar type of the field over which it is defined
 */

#include <boost/mpl/if.hpp>

#include "detail/unsupported_type.hpp"
#include "is_scalar.hpp"

namespace core {

    namespace math {

        namespace linalg {

            /** \brief a metafunction that returns the field type of the vector
             *
             *  Usage: \code field<VectorType>::type \endcode
             */
            template <typename T>
                struct field
                :   mpl::if_<
                        is_scalar<T>,
                        T,
                        unsupported_type<T>
                    >
                { };

        } // namespace linalg
    } // namespace math
} // namespace core

#include "detail/vector_traits_spec.hpp"

#endif // _core_math_linalg_traits_field_hpp_
