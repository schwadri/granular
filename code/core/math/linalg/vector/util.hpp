#ifndef _core_math_linalg_vector_util_hpp_
#define _core_math_linalg_vector_util_hpp_

/** \file util.hpp Contains the template class \a vector, that represents a linear algebraic vector.
 *  \author Adrian Schweizer
 *  \created  $Di 08 Jul 10:51:00 pm CEST 2008 schwadri@PwnererMachine.local$
 *  \modified $Di 29 Jul 09:32:41 am CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 */


namespace core {

    namespace math {

        namespace linalg {

            template <typename T>
                class zero_vector
                {
                    T operator[](int i) const
                    {
                        return 0;
                    }

                    T operator()(int i) const
                    {
                        return 0;
                    }
                };

            template <typename T>
                struct scalar_vector
                {
                    scalar_vector(T const & value)
                    :   m_value(value)
                    { }

                    T const & operator[](int i) const
                    {
                        return m_value;
                    }

                    T const & operator()(int i) const
                    {
                        return m_value;
                    }

                    T const & value() const
                    {
                        return m_value;
                    }

                private:
                    T m_value;
                };

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_vector_util_hpp_
