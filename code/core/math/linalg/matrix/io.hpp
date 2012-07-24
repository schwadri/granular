#ifndef _core_math_linalg_matrix_io_hpp_
#define _core_math_linalg_matrix_io_hpp_

/** \file io.hpp
 *  \author Adrian Schweizer
 *  \created  $Di 08 Jul 06:27:33 pm CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 *  \modified $Mi 10 Sep 05:23:01 pm CEST 2008 schwadri@guest-docking-hg-1-154.ethz.ch$
 */

#include <iosfwd>

#include "fwd.hpp"

namespace boost { namespace serialization {

        //serialization methods
        template <typename Archive, typename T, int M, int N, template <int, int> class Layout>
            void serialize(Archive & ar, core::math::linalg::matrix<T, M, N, Layout> & m, const unsigned int version)
            {
                split_free(ar, m, version);
            }
} }

namespace core {

    namespace math {

        namespace linalg {

            template <typename Archive, typename T, int M, int N, template <int, int> class Layout>
                void load(Archive & ar, matrix<T, M, N, Layout> & m, const unsigned int version)
                {
                    for(int i = 0; i < M; ++i)
                        for(int j = 0; j < N; ++j)
                            ar & m(i, j);
                }

            template <typename Archive, typename T, int M, int N, template <int, int> class Layout>
                void save(Archive & ar, matrix<T, M, N, Layout> const & m, const unsigned int version)
                {
                    for(int i = 0; i < M; ++i)
                        for(int j = 0; j < N; ++j)
                        {
                            //FIXME:
                            T t = m(i, j);
                            ar & t;
                        }
                }

            //serialization methods
            template <typename CharT, typename CharTraits, typename T, int M, int N, template <int, int> class Layout>
                std::basic_ostream<CharT, CharTraits> & operator << (std::basic_ostream<CharT, CharTraits> & out, matrix<T, M, N, Layout> const & m)
                {
                    //FIXME: this only works with char

                    out << '[';

                    for(int i = 0; i < M; ++i)
                    {
                        bool fst = true;
                        for(int j = 0; j < N; ++j)
                            if(fst)
                            {
                                out << m(i, j);
                                fst = false;
                            }
                            else
                                out << ", " << m(i, j);
                        if(i + 1 < M)
                            out << "; ";
                    }

                    out << ']';

                    return out;
                }

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_matrix_io_hpp_
