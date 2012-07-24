#ifndef _core_math_linalg_matrix_layout_hpp_
#define _core_math_linalg_matrix_layout_hpp_

/** \file layout.hpp
 *  \author Adrian Schweizer
 *  \created  $Di 08 Jul 09:59:09 pm CEST 2008 schwadri@PwnererMachine.local$
 *  \modified $Di 08 Jul 09:59:09 pm CEST 2008 schwadri@PwnererMachine.local$
 */
namespace core {

    namespace math {

        namespace linalg {

            namespace layout {

                template <int M, int N>
                    struct row_major
                    {
                        static int index_to_offset(int i, int j)
                        {
                            return N * i + j;
                        }
                    };

                template <int M, int N>
                    struct column_major
                    {
                        static int index_to_offset(int i, int j)
                        {
                            return i + M * j;
                        }
                    };
            }
        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_matrix_layout_hpp_
