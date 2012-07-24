#ifndef _core_math_linalg_utility_quaternion_hpp_
#define _core_math_linalg_utility_quaternion_hpp_

/** \file quaternion.hpp Supporting functionality for quaternions.
 *  \author Adrian Schweizer
 *  \created  $Do 30 Aug 07:27:48 pm CEST 2007 schwadri@SchwadriComp.local$
 *  \modified $Di 29 Jul 10:35:41 am CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 */

#include <boost/math/quaternion.hpp>
#include "../matrix/fwd.hpp"
#include "../vector/fwd.hpp"

namespace core {

    namespace math {

        namespace linalg {

            /** \brief convert a \p vector<T, 4> representing a quaternion to a matrix*/
            template <template <int, int> class Layout, typename T>
                linalg::matrix<T, 3, 3, Layout> to_matrix(linalg::vector<T, 4> const & q)
                {
                    return linalg::matrix<T, 3, 3, Layout>(
                            1.0-2*(q[2]*q[2])-2*(q[3]*q[3]),
                                2*(q[1]*q[2])-2*(q[0]*q[3]),
                                2*(q[1]*q[3])+2*(q[0]*q[2]),

                                2*(q[1]*q[2])+2*(q[0]*q[3]),
                            1.0-2*(q[1]*q[1])-2*(q[3]*q[3]),
                                2*(q[2]*q[3])-2*(q[0]*q[1]),

                                2*(q[1]*q[3])-2*(q[0]*q[2]),
                                2*(q[2]*q[3])+2*(q[0]*q[1]),
                            1.0-2*(q[1]*q[1])-2*(q[2]*q[2]));
                }

            /** \brief convert a \p boost::math::quaternion to a matrix*/
            template <template <int, int> class Layout, typename T>
                linalg::matrix<T, 3, 3, Layout> to_matrix(boost::math::quaternion<T> const & q)
                {
                    return linalg::matrix<T, 3, 3, Layout>(
                            1.0-2*(q.R_component_3()*q.R_component_3())-2*(q.R_component_4()*q.R_component_4()),
                                2*(q.R_component_2()*q.R_component_3())-2*(q.R_component_1()*q.R_component_4()),
                                2*(q.R_component_2()*q.R_component_4())+2*(q.R_component_1()*q.R_component_3()),

                                2*(q.R_component_2()*q.R_component_3())+2*(q.R_component_1()*q.R_component_4()),
                            1.0-2*(q.R_component_2()*q.R_component_2())-2*(q.R_component_4()*q.R_component_4()),
                                2*(q.R_component_3()*q.R_component_4())-2*(q.R_component_1()*q.R_component_2()),

                                2*(q.R_component_2()*q.R_component_4())-2*(q.R_component_1()*q.R_component_3()),
                                2*(q.R_component_3()*q.R_component_4())+2*(q.R_component_1()*q.R_component_2()),
                            1.0-2*(q.R_component_2()*q.R_component_2())-2*(q.R_component_3()*q.R_component_3()));
                }

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_utility_quaternion_hpp_
