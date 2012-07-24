#ifndef _core_math_linalg_vector_basic_vector_hpp_
#define _core_math_linalg_vector_basic_vector_hpp_

/** \file basic_vector.hpp
 *  \author Adrian Schweizer
 *  \created  $Do 23 Aug 10:32:12 pm CEST 2007 schwadri@SchwadriComp.local$
 *  \modified $Fr 18 Jul 06:10:37 pm CEST 2008 schwadri@SchwadriLaptop.local$
 */

#include <algorithm>
#include <functional>

#include <boost/array.hpp>

#include "range.hpp"

namespace core {

    namespace math {

        namespace linalg {

            template <typename F, int Size, typename Derived>
                class basic_vector
                :   public boost::array<F, Size>
                {
                public:
                    typedef basic_vector<
                        F,
                        Size,
                        Derived
                    > type;

                    typedef boost::array<F, Size>               base_type;
                    typedef F                                   value_type;
                    typedef value_type *                        pointer;
                    typedef value_type const *                  const_pointer;
                    typedef typename base_type::iterator        iterator;
                    typedef typename base_type::const_iterator  const_iterator;

                    static const int size = Size;

                    using base_type::begin;
                    using base_type::end;
                    using base_type::operator[];

                    /** \name Structors */      //@{

                    /** Trivial Constructor*/
                    basic_vector()
                    { }

                    /** Copy-Constructor*/
                    template <typename F2>
                        basic_vector(basic_vector<F2, Size, Derived> const & other)
                        {
                            using std::copy;
                            copy(other.begin(), other.end(), begin());
                        }

                    explicit basic_vector(value_type const & value)
                    {
                        using std::fill;
                        fill(begin(), end(), value);
                    }

                    //@}

                    /** \name Iterators */ //@{
                    /*iterator begin()
                    {
                        return m_data;
                    }

                    iterator end()
                    {
                        return m_data + size;
                    }

                    const_iterator begin() const
                    {
                        return m_data;
                    }

                    const_iterator end() const
                    {
                        return m_data + size;
                    }*/

                    //@}

                    /** \name Operators */ //@{

                    template <typename T, T First, T Last>
                        vector_static_range_proxy<Derived, v_range_c<T, First, Last> > operator[](v_range_c<T, First, Last> (*hint)())
                        {
                            return vector_static_range_proxy<Derived, v_range_c<T, First, Last> >(*This());
                        }

                    template <typename T, T First, T Last>
                        vector_static_range_proxy<Derived const, v_range_c<T, First, Last> > const operator[](v_range_c<T, First, Last> (*hint)()) const
                        {
                            return vector_static_range_proxy<Derived const, v_range_c<T, First, Last> >(*This());
                        }

                    template <typename Scalar>
                        Derived & operator *= (Scalar const & s)
                        {
                            using std::transform;
                            using std::bind2nd;
                            using std::multiplies;
                            transform(
                                begin(),
                                end(),
                                begin(),
                                bind2nd(multiplies<Scalar>(), s)
                            );

                            return static_cast<Derived &>(*this);
                        }

                    template <typename Scalar>
                        Derived & operator += (Scalar const & s)
                        {
                            using std::transform;
                            using std::bind2nd;
                            using std::plus;
                            transform(
                                begin(),
                                end(),
                                begin(),
                                bind2nd(plus<Scalar>(), s)
                            );

                            return static_cast<Derived &>(*this);
                        }

                    template <typename Scalar>
                        Derived & operator -= (Scalar const & s)
                        {
                            using std::transform;
                            using std::bind2nd;
                            using std::minus;
                            transform(
                                begin(),
                                end(),
                                begin(),
                                bind2nd(minus<Scalar>(), s)
                            );

                            return static_cast<Derived &>(*this);
                        }

                    template <typename Scalar>
                        Derived & operator /= (Scalar const & s)
                        {
                            using std::transform;
                            using std::bind2nd;
                            using std::divides;
                            transform(
                                begin(),
                                end(),
                                begin(),
                                bind2nd(divides<Scalar>(), s)
                            );

                            return static_cast<Derived &>(*this);
                        }

                    Derived & operator += (Derived const & other)
                    {
                        using std::transform;
                        using std::plus;
                        transform(
                            begin(),
                            end(),
                            other.begin(),
                            begin(),
                            plus<value_type>()
                        );

                        return static_cast<Derived &>(*this);
                    }

                    Derived & operator -= (Derived const & other)
                    {
                        using std::transform;
                        using std::minus;
                        transform(
                            begin(),
                            end(),
                            other.begin(),
                            begin(),
                            minus<value_type>()
                        );

                        return static_cast<Derived &>(*this);
                    }

                    Derived & operator *= (Derived const & other)
                    {
                        using std::transform;
                        using std::multiplies;
                        transform(
                            begin(),
                            end(),
                            other.begin(),
                            begin(),
                            multiplies<value_type>()
                        );

                        return *This();
                    }

                    Derived & operator /= (Derived const & other)
                    {
                        using std::transform;
                        using std::divides;
                        transform(
                            begin(),
                            end(),
                            other.begin(),
                            begin(),
                            divides<value_type>()
                        );

                        return *This();
                    }

                    Derived operator-() const
                    {
                        Derived result;
                        using std::transform;
                        using std::negate;
                        transform(
                            begin(),
                            end(),
                            result.begin(),
                            negate<value_type>()
                        );
                        return result;
                    }

                private:
                    Derived * This()
                    {
                        return static_cast<Derived *>(this);
                    }

                    Derived const * This() const
                    {
                        return static_cast<Derived const *>(this);
                    }

                };

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_vector_basic_vector_hpp_
