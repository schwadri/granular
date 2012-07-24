/** \file clamp.hpp
 *  \author Adrian Schweizer
 *  \created  $Do 23 Aug 10:32:12 pm CEST 2007 schwadri@SchwadriComp.local$
 *  \modified $Di 08 Jul 11:08:47 pm CEST 2008 schwadri@PwnererMachine.local$
 */

#ifndef _core_math_clamp_hpp_
#define _core_math_clamp_hpp_

#include <cmath> // using min, max

namespace core {

    namespace math {

        template <typename T>
            T clamp(T t, T t0, T t1) {
                using std::min;
                using std::max;
                return min(max(t, t0), t1);
            }

    } // namespace math
} // namespace core

#endif // _core_math_clamp_hpp_
