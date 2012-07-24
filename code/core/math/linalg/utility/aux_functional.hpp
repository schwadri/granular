#ifndef _core_math_linalg_utility_aux_functional_hpp_
#define _core_math_linalg_utility_aux_functional_hpp_

#include <cmath>
#include <algorithm>

namespace core {
    
    namespace math {
        
        namespace linalg {

            /** \brief contains auxiliary binary function objects for use with stl-algorithms and their static counterparts*/
            namespace aux {

                /** \brief binary functor calculating \f$ x_1 + \mathrm{abs} \left(x_2\right) \f$*/
                template <typename T>
                    struct plus_abs
                    {
                        T operator()(T const & val, T const & r) const
                        {
                            using std::abs;
                            return val + abs(r);
                        }
                    };

                /** \brief binary functor that calculates \f$ x_1 + x_2^2 \f$*/
                template<class T>
                    struct plus_sqr
                    {
                        T operator()(T const & val, T const & r) const
                        {
                            return val + r*r;
                        }
                    };

                /** \brief binary functor calculating \f$ \max \left\{x_1 , \mathrm{abs} \left(x_2\right)\right\} \f$*/
                template <typename T>
                    struct max_abs
                    {
                        T operator()(T const & l, T const & r) const
                        {
                            using std::max;
                            using std::abs;
                            return max(l, abs(r));
                        }
                    };

            } // namespace aux
        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_utility_aux_functional_hpp_
