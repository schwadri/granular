#ifndef _core_math_linalg_vector_transform_hpp_
#define _core_math_linalg_vector_transform_hpp_

#include <cstddef>

#include <boost/mpl/identity.hpp>
#include <boost/mpl/eval_if.hpp>

#include "../traits/field.hpp"
#include "../traits/is_scalar.hpp"

namespace core {

    namespace math {

        namespace linalg {

            namespace detail {

                namespace mpl = boost::mpl;

                struct transform_impl
                {
                    template <typename V, typename VR, typename F>
                        static void eval(V const & v1, V const & v2, VR & v3, F binary_op)
                        {
                            //check that all 3 vectors have the same dimension
                            std::size_t dim =   dimension<V>::get(v1);

                            assert(
                                dim == dimension<V>::get(v2) &&
                                dim == dimension<V>::get(v3)
                            );

                            using std::transform;
                            transform(
                                v1.begin(), v1.end(),
                                v2.begin(), v3.begin(),
                                binary_op
                            );
                        };

                    template <typename V, typename VR, typename F>
                        static void eval(V const & v1, VR & v2, F unary_op)
                        {
                            //check that both vectors have the same dimension
                            std::size_t dim = dimension<V>::get(v1);
                            assert(dim == dimension<V>::get(v2));

                            using std::transform;
                            transform(
                                v1.begin(), v1.end(),
                                v2.begin(),
                                unary_op
                            );
                        };
                };

                struct transform_scalar_impl
                {
                    template <typename V, typename VR, typename F>
                        static void eval(const V& v1, const V& v2, VR& v3, F binary_op)
                        {
                            v3 = binary_op(v1,v2);
                        };

                    template <typename V, typename VR, typename F>
                        static void eval(const V& v1, VR& v2, F unary_op)
                        {
                            v2 = unary_op(v1);
                        };
                };

                template <typename V, typename VR, typename F>
                    struct choose_transform_impl
                    {
                        //check if this is a vector with static dimension
                        typedef typename mpl::eval_if<
                            is_scalar<V>,
                            mpl::identity<transform_scalar_impl>,
                            mpl::identity<transform_impl>
                        >::type  impl_type;

                        static void eval(V const & v1, V const & v2, VR & v3, F binary_op)
                        {
                            impl_type::template eval<V, VR, F>(v1, v2, v3, binary_op);
                        }

                        static void eval(V const & v1, VR & v2, F unary_op)
                        {
                            impl_type::template eval<V, VR, F>(v1, v2, unary_op);
                        }
                    };
            } // namespace detail

            template <typename Vector, typename VectorR,  typename BinaryFunctor>
                void transform(Vector const & v1, Vector const & v2, VectorR & v3, BinaryFunctor & binary_op)
                {
                    typedef detail::choose_transform_impl<Vector, VectorR, BinaryFunctor&>   impl_type;
                    impl_type::eval(v1, v2, v3, binary_op);
                };

            template <typename Vector, typename VectorR,  typename BinaryFunctor>
                void transform(Vector const & v1, Vector const & v2, VectorR& v3, BinaryFunctor const & i_binary_op)
                {
                    BinaryFunctor binary_op(i_binary_op);
                    typedef detail::choose_transform_impl<Vector, VectorR, BinaryFunctor&>   impl_type;
                    impl_type::eval(v1, v2, v3, binary_op);
                };

            template <typename Vector, typename VectorR,  typename UnaryFunctor>
                void transform(Vector const & v1, VectorR & v2, UnaryFunctor & unary_op)
                {
                    typedef detail::choose_transform_impl<Vector, VectorR, UnaryFunctor&>   impl_type;
                    impl_type::eval(v1, v2, unary_op);
                };

            template <typename Vector, typename VectorR, typename UnaryFunctor>
                void transform(Vector const & v1, VectorR & v2, UnaryFunctor const & i_unary_op)
                {
                    UnaryFunctor unary_op(i_unary_op);
                    typedef detail::choose_transform_impl<Vector, VectorR, UnaryFunctor&>   impl_type;
                    impl_type::eval(v1, v2, unary_op);
                };
        } // namespace linalg

    } // namespace math
} // namespace core

#endif // _core_math_linalg_vector_transform_hpp_
