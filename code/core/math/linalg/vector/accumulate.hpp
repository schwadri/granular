#ifndef _core_math_linalg_vector_accumulate_hpp_
#define _core_math_linalg_vector_accumulate_hpp_

#include <numeric>

#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>

#include "../traits/is_scalar.hpp"
#include "../traits/field.hpp"

namespace core {

    namespace math {

        namespace linalg {

            namespace detail {

                namespace mpl = boost::mpl;

                struct accumulate_impl
                {
                    template <typename Vector, typename BinaryFunctor>
                        inline static typename field<Vector>::type eval(Vector const & v, typename field<Vector>::type const & init_value, BinaryFunctor binary_op)
                        {
                            using std::accumulate;

                            return accumulate(
                                v.begin(),
                                v.end(),
                                init_value,binary_op
                            );
                        }
                };

                struct accumulate_scalar_impl
                {
                    template <typename S, typename F>
                        inline static S eval(S const & value, const S& init_value, F binary_op)
                        {
                            return binary_op(init_value, value);
                        }
                };

                template <typename T, typename F>
                    struct choose_accumulate_impl
                    {
                        //check if this is a vector with static dimension
                        typedef typename mpl::eval_if<
                            is_scalar<T>,
                            mpl::identity<accumulate_scalar_impl>,
                            mpl::identity<accumulate_impl>
                        >::type impl_type;

                        typedef typename field<T>::type scalar_type;
                        inline static scalar_type eval(const T& v, const scalar_type& init_value, F binary_op)
                        {
                            return impl_type::template eval<T, F>(v, init_value, binary_op);
                        }
                    };

            } // namespace detail

            /** \brief std::accumulate-like method with optional loop unrolling for static arrays and vectors
             *
             */
            template <typename Vector, typename BinaryFunctor>
                typename field<Vector>::type accumulate(Vector const & v, typename field<Vector>::type const & init_value, BinaryFunctor & binary_op)
                {
                    typedef typename detail::choose_accumulate_impl<Vector, BinaryFunctor&> impl_type;
                    return impl_type::eval(v,init_value,binary_op);
                }

            template <typename Vector, typename BinaryFunctor>
                typename field<Vector>::type accumulate(Vector const & v, typename field<Vector>::type const & init_value, BinaryFunctor const & i_binary_op)
                {
                    BinaryFunctor binary_op(i_binary_op);
                    typedef typename detail::choose_accumulate_impl<Vector, BinaryFunctor&> impl_type;
                    return impl_type::eval(v, init_value, binary_op);
                }

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_vector_accumulate_hpp_
