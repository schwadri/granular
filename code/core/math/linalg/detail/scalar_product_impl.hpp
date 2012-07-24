#ifndef _scalar_product_impl_hpp_
#define _scalar_product_impl_hpp_

/** \file inner_product_impl.hpp Implementations for math::inner_product
 *  \internal
 *  \author Adrian Schweizer
 *  \date   8.7.2007
 */

#include "../math_traits.hpp"
#include "../vector_traits.hpp"

#include <numeric>
#include <functional>

#include "../inner_product.hpp"

namespace core {

    namespace math {

        namespace detail {

            template<class V>
                struct scalar_product_vector_impl
                {
                    typedef typename field<V>::type scalar_type;
                    static typename field<V>::type eval(const V& v0, const V& v1)
                    {
                        return inner_product(v0,v1,scalar_type(0),std::plus<scalar_type>(),std::multiplies<scalar_type>());
                    }
                };

            template<class T>
                struct scalar_product_scalar_impl
                {
                    static T eval(const T& v0, const T& v1)
                    {
                        return v0*v1;
                    }
                };

            template<class T>
                struct choose_scalar_product_impl
                :   public boost::mpl::eval_if_c<is_scalar<T>::value,
                                                boost::mpl::identity<scalar_product_scalar_impl<T> >,
                                                boost::mpl::identity<scalar_product_vector_impl<T> >
                                                >
                { };

        } // namespace detail

    } // namespace math

} // namespace core

#endif // _scalar_product_impl_hpp_
