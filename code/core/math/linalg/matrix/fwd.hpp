#ifndef _core_math_linalg_matrix_fwd_hpp_
#define _core_math_linalg_matrix_fwd_hpp_

/** \file fwd.hpp
 *  \author Adrian Schweizer
 *  \created  $Di 08 Jul 09:17:59 pm CEST 2008 schwadri@PwnererMachine.local$
 *  \modified $Di 08 Jul 11:33:59 pm CEST 2008 schwadri@PwnererMachine.local$
 */

namespace core {

    namespace math {

        namespace linalg {

            namespace layout {
                //forward decl
                template <int, int>
                    struct row_major;
            }

            template <
                typename T,
                int M,
                int N = M,
                template <int, int> class Layout = layout::row_major
            >
                class matrix;

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_matrix_fwd_hpp_
