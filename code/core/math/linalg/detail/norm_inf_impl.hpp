#ifndef _core_math_linalg_detail_norm_inf_impl_hpp_
#define _core_math_linalg_detail_norm_inf_impl_hpp_

#include <cmath>
#include <boost/mpl/eval_if.hpp>

#include "../traits/is_scalar.hpp"
#include "../traits/field.hpp"

#include "../vector/accumulate.hpp"
#include "../utility/aux_functional.hpp"

namespace core {

    namespace math {

        namespace linalg {

            namespace detail {

                template <typename T>
                    struct norm_inf_scalar_impl
                    {
                        static T eval(const T& v)
                        {
                            using std::abs;
                            return abs(v);
                        }
                    };

                template <typename V>
                    struct norm_inf_vector_impl
                    {
                        typedef typename field<V>::type scalar_type;
                        static scalar_type eval(const V& v)
                        {
                            return accumulate(
                                v,
                                scalar_type(0),
                                aux::max_abs<scalar_type>()
                            );
                        }
                    };

                template <typename T>
                    struct choose_norm_inf_impl
                    :   mpl::eval_if<
                            is_scalar<T>,
                            mpl::identity<norm_inf_scalar_impl<T> >,
                            mpl::identity<norm_inf_vector_impl<T> >
                        >
                    { };

            } // namespace detail
        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_detail_norm_inf_impl_hpp_
