#ifndef _vector_traits_spec_hpp_
#define _vector_traits_spec_hpp_

/* \file vector_traits_spec.hpp Contains specializations of the field
 *  and dimension metafunctions
 *  \author Adrian Schweizer
 *  \created  $Di 08 Jul 11:03:32 am CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 *  \modified $Di 08 Jul 11:08:49 am CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 */

#include <boost/mpl/size_t.hpp>

#include "fwd_decls.hpp"

namespace core {

    namespace math {

        namespace linalg {

            //fwd decl
            template <typename V>
                struct dimension;

            template <typename V>
                struct field;

            namespace detail {

                template <typename T, std::size_t N>
                    struct static_dimension_impl : mpl::size_t<N>
                    {
                        static std::size_t get(const T& t)
                        {
                            return N;
                        }
                    };

                template <typename T>
                    struct stl_dynamic_dimension_impl : mpl::size_t<0>
                    {
                        static std::size_t get(const T& t)
                        {
                            return t.size();
                        }
                    };

            } // namespace detail

            template <typename T, std::size_t N>
                struct dimension<boost::array<T, N> >
                : detail::static_dimension_impl<boost::array<T, N>, N>
                { };

            template <typename T, std::size_t N>
                struct dimension<T[N]>
                : detail::static_dimension_impl<T[N], N>
                { };

            template <typename T, typename A>
                struct dimension<std::vector<T, A> >
                : detail::stl_dynamic_dimension_impl<std::vector<T, A> >
                { };

            template <typename T>
                struct dimension<std::valarray<T> >
                : detail::stl_dynamic_dimension_impl<std::valarray<T> >
                { };

            //specialization for Point, Index, std::vector, std::valarray, boost::array and builtin static arrays
            template <typename T, std::size_t N>
                struct field<boost::array<T, N> >
                {
                    typedef T type;
                };

            template <typename T, std::size_t N>
                struct field<T[N] >
                {
                    typedef T type;
                };

            template <typename T, typename A>
                struct field<std::vector<T, A> >
                {
                    typedef T type;
                };

            template <typename T>
                struct field<std::valarray<T> >
                {
                    typedef T type;
                };

            //TODO:queries to determine the behaviour of the different operator implementations of the vectors

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _vector_traits_spec_hpp_
