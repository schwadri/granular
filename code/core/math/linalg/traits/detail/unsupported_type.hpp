#ifndef _core_math_linalg_traits_detail_unsupported_type_hpp_
#define _core_math_linalg_traits_detail_unsupported_type_hpp_

/** \file unsupported_type.hpp
 *  \author Adrian Schweizer
 *  \created  $Di 08 Jul 11:06:19 am CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 *  \modified $Di 08 Jul 11:06:37 am CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 */

//TODO: move somewhere where it makes sence
//keep msvc from barking
#if defined _MSC_VER && _MSC_VER >= 1300
    #if !defined _CRT_SECURE_NO_WARNINGS
        #define _CRT_SECURE_NO_WARNINGS
    #endif
    #if !defined _SCL_SECURE_NO_WARNINGS
        #define _SCL_SECURE_NO_WARNINGS
    #endif
    #if !defined _CRT_NONSTDC_NO_WARNINGS
        #define _CRT_NONSTDC_NO_WARNINGS
    #endif
#endif

namespace core {

    namespace math {

        namespace linalg {

            template<class T>
                struct unsupported_type;

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_traits_detail_unsupported_type_hpp_
