#ifndef _core_math_linalg_vector_norm_2_hpp_
#define _core_math_linalg_vector_norm_2_hpp_

#include <cmath>

#include "../detail/norm_2_impl.hpp"

namespace core {

    namespace math {

        namespace linalg {

            /** \brief calculates \f$|| . ||_2^2\f$*/
            template <typename T>
                typename field<T>::type norm_2_sqr(T const & v)
                {
                    typedef typename detail::choose_norm_2_sqr_impl<T>::type impl_type;
                    return impl_type::eval(v);
                }

            /** \brief calculates \f$||.||_2\f$*/
            template <typename T>
                typename field<T>::type norm_2(T const & v)
                {
                    using std::sqrt;
                    return sqrt(norm_2_sqr(v));
                }

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_vector_norm_2_hpp_
