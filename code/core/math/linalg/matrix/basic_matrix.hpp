#ifndef _core_math_linalg_matrix_basic_matrix_hpp_
#define _core_math_linalg_matrix_basic_matrix_hpp_

/** \file basic_matrix.hpp
 *  \author Adrian Schweizer
 *  \created  $Mo 23 Aug 10:32:12 pm CNST 2007 schwadri@SchwadriComp.local$
 *  \modified $Mo 28 Jul 01:28:08 pm CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 */

#include <algorithm>
#include <functional>

#include <boost/array.hpp>

#include "layout.hpp"
#include "fwd.hpp"
#include "range.hpp"

namespace core {

    namespace math {

        namespace linalg {

            template <typename T, int M, int N, template <int, int> class Layout, typename Derived>
                class basic_matrix
                :   boost::array<T, M * N>
                {
                    typedef Layout<M, N>    layout;
                public:
                    struct row_proxy
                    {
                        row_proxy(T * data)
                        :   m_data(data)
                        { }

                        T & operator [](int j)
                        {
                            return *(m_data + layout::index_to_offset(0, j));
                        }

                    private:
                        T * m_data;
                    };

                    struct const_row_proxy
                    {
                        const_row_proxy(T const * data)
                        :   m_data(data)
                        { }

                        T operator [](int j) const
                        {
                            return *(m_data + layout::index_to_offset(0, j));
                        }

                    private:
                        T const * m_data;
                    };

                    typedef T                                       value_type;
                    typedef basic_matrix<T, M, N, Layout, Derived>  my_type;
                    typedef boost::array<T, M * N>                  base_type;
                    typedef value_type*                             pointer;
                    typedef value_type const *                      const_pointer;
                    typedef pointer                                 iterator;
                    typedef const_pointer                           const_iterator;
                    typedef value_type&                             reference;
                    typedef value_type const &                      const_reference;

                    using base_type::begin;
                    using base_type::end;

                    basic_matrix()
                    { }

                    explicit basic_matrix(value_type const & value)
                    {
                        using std::fill;
                        fill(begin(), end(), value);
                    }

                    row_proxy operator[](int i)
                    {
                        return row_proxy(this->elems + layout::index_to_offset(i, 0));
                    }

                    const_row_proxy operator[](int i) const
                    {
                        return const_row_proxy(this->elems + layout::index_to_offset(i, 0));
                    }

                    template <typename U, U First0, U Last0, U First1, U Last1>
                        matrix_static_range_proxy<Derived, m_range_c<U, First0, Last0, First1, Last1> > operator[](m_range_c<U, First0, Last0, First1, Last1> (*hint)())
                        {
                            return matrix_static_range_proxy<Derived, m_range_c<U, First0, Last0, First1, Last1> >(*This());
                        }

                    template <typename U, U First0, U Last0, U First1, U Last1>
                        matrix_static_range_proxy<Derived const, m_range_c<U, First0, Last0, First1, Last1> > const operator[](m_range_c<U, First0, Last0, First1, Last1> (*hint)()) const
                        {
                            return matrix_static_range_proxy<Derived const, m_range_c<T, First0, Last0, First1, Last1> >(*This());
                        }

                    value_type & operator()(int i, int j)
                    {
                        return *(this->elems + layout::index_to_offset(i, j));
                    }

                    value_type operator()(int i, int j) const
                    {
                        return *(this->elems + layout::index_to_offset(i, j));
                    }

                    template <typename Scalar>
                        Derived & operator *= (Scalar const & s)
                        {
                            using std::transform;
                            using std::multiplies;
                            using std::bind2nd;

                            transform(
                                begin(),
                                end(),
                                begin(),
                                bind2nd(multiplies<Scalar>(), s)
                            );

                            return *This();
                        }

                    template <typename Scalar>
                        Derived & operator += (Scalar const & s)
                        {
                            using std::transform;
                            using std::plus;
                            using std::bind2nd;

                            transform(
                                begin(),
                                end(),
                                begin(),
                                bind2nd(plus<Scalar>(), s)
                            );

                            return *This();
                        }

                    template <typename Scalar>
                        Derived &  operator -= (Scalar const & s)
                        {
                            using std::transform;
                            using std::minus;
                            using std::bind2nd;

                            transform(
                                begin(),
                                end(),
                                begin(),
                                bind2nd(minus<Scalar>(), s)
                            );

                            return *This();
                        }

                    template <typename Scalar>
                        Derived & operator /= (Scalar const & s)
                        {
                            using std::transform;
                            using std::divides;
                            using std::bind2nd;

                            transform(
                                begin(),
                                end(),
                                begin(),
                                bind2nd(divides<Scalar>(), s)
                            );

                            return *This();
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

                        return *This();
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

                        return *This();
                    }


                    template <template <int, int> class L>
                        Derived & operator *= (matrix<value_type, N, N, L> const & other)
                        {
                            for(int i = 0; i < M; ++i)
                            {
                                value_type temp_row[M];
                                for(int j = 0; j < N; ++j)
                                {
                                    temp_row[j] = 0.0;

                                    for(int u = 0; u < N; ++u)
                                        temp_row[j] += (*this)(i, u) * other(u, j);
                                }

                                for(int j = 0; j < N; ++j)
                                    (*this)(i, j) = temp_row[j];
                            }

                            return *This();
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

#endif // _core_math_linalg_matrix_basic_matrix_hpp_
