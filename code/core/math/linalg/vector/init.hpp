#ifndef _core_math_linalg_vector_init_hpp_
#define _core_math_linalg_vector_init_hpp_

/** \file init.hpp
 *  \author Adrian Schweizer
 *  \created  $Mi 29 Aug 03:11:38 pm CEST 2007 schwadri@SchwadriComp.local$
 *  \modified $Fr 18 Jul 05:37:23 pm CEST 2008 schwadri@SchwadriLaptop.local$
 */

#include <boost/static_assert.hpp>

namespace core {

    namespace math {

        namespace linalg {

            //forward decl
            template <typename T, int N>
                class vector;

            template <typename Vector, int I>
                struct vector_init
                {

                    vector_init(Vector & v_)
                    :   v(v_)
                    { }

                    typedef vector_init<
                        Vector,
                        I + 1
                    > next_type;

                    //if this assert fails, you provided too many values
                    BOOST_STATIC_ASSERT(I <= Vector::size);

                    template <typename Scalar>
                        next_type operator,(Scalar const & value)
                        {
                            v[I] = value;
                            return next_type(v);
                        }

                    Vector & v;
                };

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_vector_init_hpp_
