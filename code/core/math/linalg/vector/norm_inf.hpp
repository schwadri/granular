#ifndef _core_math_linalg_vector_norm_inf_hpp_
#define _core_math_linalg_vector_norm_inf_hpp_

#include "../traits/field.hpp"
#include "../detail/norm_inf_impl.hpp"

namespace core {

    namespace math {

        namespace linalg {

            /** \brief \f$ || . ||_\infty \f$*/
            template <typename T>
                typename field<T>::type norm_inf(T const & v)
                {
                    typedef typename detail::choose_norm_inf_impl<T>::type impl_type;
                    return impl_type::eval(v);
                }

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_vector_norm_inf_hpp_
