#ifndef _core_math_linalg_vector_norm_1_hpp_
#define _core_math_linalg_vector_norm_1_hpp_

#include "../detail/norm_1_impl.hpp"

namespace core {

    namespace math {

        namespace linalg {

            /** \brief \f$ || . ||_1 \f$*/
            template <typename T>
                typename field<T>::type norm_1(T const & v)
                {
                    typedef typename detail::choose_norm_1_impl<T>::type impl_type;
                    return impl_type::eval(v);
                }

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_vector_norm_1_hpp_
