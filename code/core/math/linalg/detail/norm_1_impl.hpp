#ifndef _core_math_linalg_detail_norm_1_impl_hpp_
#define _core_math_linalg_detail_norm_1_impl_hpp_

#include <cmath>

#include "../traits/is_scalar.hpp"
#include "../traits/field.hpp"

#include "../vector/accumulate.hpp"
#include "../utility/aux_functional.hpp"

namespace core {

    namespace math {

        namespace linalg {

            namespace detail {

                template<class T>
                    struct norm_1_scalar_impl
                    {
                        static T eval(const T& v)
                        {
                            using std::abs;
                            return abs(v);
                        }
                    };

                template<class V>
                    struct norm_1_vector_impl
                    {
                        typedef typename field<V>::type scalar_type;
                        static typename field<V>::type eval(const V& v)
                        {
                            return accumulate(
                                v,
                                scalar_type(0),
                                aux::plus_abs<scalar_type>()
                            );
                        }
                    };

                template<class T>
                    struct choose_norm_1_impl
                    :   mpl::eval_if<
                            is_scalar<T>,
                            mpl::identity<norm_1_scalar_impl<T> >,
                            mpl::identity<norm_1_vector_impl<T> >
                        >
                    { };

            } // namespace detail
        } // namespace linalg
    } // namespace vector
} // namespace math

#endif // _core_math_linalg_detail_norm_1_impl_hpp_
