#ifndef _core_math_linalg_detail_norm_2_impl_hpp_
#define _core_math_linalg_detail_norm_2_impl_hpp_

#include "../traits/field.hpp"
#include "../traits/is_scalar.hpp"

#include "../vector/accumulate.hpp"
#include "../utility/aux_functional.hpp"

namespace core {

    namespace math {

        namespace linalg {

            /** \internal */
            namespace detail {

                template<class T>
                    struct norm_2_sqr_scalar_impl
                    {
                        static T eval(const T& v)
                        {
                            return v*v;
                        }
                    };

                template<class V>
                    struct norm_2_sqr_vector_impl
                    {
                        typedef typename field<V>::type scalar_type;
                        static scalar_type eval(const V& v)
                        {
                            return accumulate(v,scalar_type(0),aux::plus_sqr<scalar_type>());
                        }
                    };

                template<class T>
                    struct choose_norm_2_sqr_impl
                    :   mpl::eval_if<
                            is_scalar<T>,
                            mpl::identity<norm_2_sqr_scalar_impl<T> >,
                            mpl::identity<norm_2_sqr_vector_impl<T> >
                        >
                    { };

            } // namespace detail

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_detail_norm_2_impl_hpp_
