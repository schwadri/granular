#ifndef _core_math_linalg_vector_io_hpp_
#define _core_math_linalg_vector_io_hpp_

/** \file io.hpp
 *  \author Adrian Schweizer
 *  \created  $Di 08 Jul 01:26:27 pm CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 *  \modified $Mi 10 Sep 05:10:29 pm CEST 2008 schwadri@guest-docking-hg-1-154.ethz.ch$
 */

#include <iosfwd>
#include <algorithm>
#include <boost/serialization/split_free.hpp>

#include "fwd.hpp"

namespace boost { namespace serialization {

        //serialization methods
        template <typename Archive, typename T, int Size>
            void serialize(Archive & ar, core::math::linalg::vector<T, Size> & v, const unsigned int version)
            {
                split_free(ar, v, version);
            }
} }

namespace core {

    namespace math {

        namespace linalg {

            template <typename Archive, typename T, int Size>
                void load(Archive & ar, vector<T, Size> & v, const unsigned int version)
                {
                    for(int i = 0; i < Size; ++i)
                        ar & v[i];
                }

            template <typename Archive, typename T, int Size>
                void save(Archive & ar, vector<T, Size> const & v, const unsigned int version)
                {
                    for(int i = 0; i < Size; ++i)
                        ar & v[i];
                }


            template <typename CharT, typename CharTraits, typename T, int Size>
                std::basic_ostream<CharT, CharTraits> & operator<<(std::basic_ostream<CharT, CharTraits> & out, vector<T, Size> const & v)
                {
                    //FIXME: this only works with char

                    out << '[';

                    using std::copy;
                    copy(v.begin(), v.end() - 1, std::ostream_iterator<T>(out, ", "));

                    out << v.back() << ']';

                    return out;
                }

            template <typename CharT, typename CharTraits, typename Vector, typename Range>
                std::basic_ostream<CharT, CharTraits> & operator<<(std::basic_ostream<CharT, CharTraits> & out, vector_static_range_proxy<Vector, Range> const & v)
                {
                    //FIXME: this only works with char

                    out << '[';

                    using std::copy;
                    copy(v.begin(), v.end() - 1, std::ostream_iterator<typename Vector::value_type>(out, ", "));

                    out << v.back() << ']';

                    return out;
                }

            //TODO:for istream

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_vector_io_hpp_
