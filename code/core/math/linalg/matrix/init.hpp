#ifndef _core_math_linalg_matrix_init_hpp_
#define _core_math_linalg_matrix_init_hpp_

/** \file init.hpp
 *  \author Adrian Schweizer
 *  \created  $Mi 29 Aug 03:11:38 pm CEST 2007 schwadri@SchwadriComp.local$
 *  \modified $Di 08 Jul 11:28:22 pm CEST 2008 schwadri@PwnererMachine.local$
 */

#include <boost/static_assert.hpp>

#include "fwd.hpp"

namespace core {

    namespace math {

        namespace linalg {

            template <typename Matrix, int I, int J>
                struct matrix_init;

            template <
                typename T,
                int M, int N,
                int I, int J,
                template <int, int> class Layout
            >
                struct matrix_init<matrix<T, M, N, Layout>, I, J>
                {
                    typedef matrix<T, M, N, Layout> matrix_t;

                    matrix_init(matrix_t & m_)
                    :   m(m_)
                    { }

                    //if this assert fails, you provided too many values
                    BOOST_STATIC_ASSERT(I < M && J < N || I == M && J == 0);

                    typedef matrix_init<
                        matrix_t,
                        (J + 1) < N ? I : I + 1,
                        (J + 1) < N ? J + 1 : 0
                    > next_type;

                    template <typename Scalar>
                        next_type  operator,(Scalar const & value)
                        {
                            m(I, J) = value;
                            return next_type(m);
                        }

                    matrix_t & m;
                };

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_matrix_init_hpp_
