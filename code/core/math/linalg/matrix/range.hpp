#ifndef _core_math_linalg_matrix_range_hpp_
#define _core_math_linalg_matrix_range_hpp_

/** \file range.hpp
 *  \author Adrian Schweizer
 *  \created  $Fr 18 Jul 12:26:36 pm CEST 2008 schwadri@SchwadriLaptop.local$
 *  \modified $Di 29 Jul 10:54:55 am CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 */

namespace core {

    namespace math {

        namespace linalg {

            template <typename T, T First0, T Last0, T First1, T Last1>
                struct m_range_c
                {
                    typedef T   value_type;
                    static const T first0    = First0;
                    static const T last0     = Last0;
                    static const T first1    = First1;
                    static const T last1     = Last1;
                };

            template <typename T>
                struct m_range_d
                {
                    m_range_d(T const & first0_, T const & last0_, T const & first1_, T const & last1_)
                    :   first0(first0_),
                        last0(last0_),
                        first1(first1_),
                        last1(last1_)
                    { }

                    T first0,
                      last0,
                      first1,
                      last1;
                };

            template <typename Matrix, typename Range>
                struct matrix_static_range_proxy
                {
                    typedef matrix_static_range_proxy       type;
                    typedef typename Matrix::value_type     value_type;
                    typedef Range                           range;
                    static const int                        size0 = Range::last0 - Range::first0;
                    static const int                        size1 = Range::last1 - Range::first1;
                    typedef typename Matrix::template rebind<size0, size1>::type matrix_type;

                    matrix_static_range_proxy(Matrix & m_)
                    :   m(m_)
                    { }

                    value_type & operator()(int i, int j)
                    {
                        return m(Range::first0 + i, Range::first1 + j);
                    }

                    value_type const & operator()(int i, int j) const
                    {
                        return m(Range::first0 + i, Range::first1 + j);
                    }

                    type & operator=(matrix_type const & other)
                    {
                        for(int i = 0; i < size0; ++i)
                            for(int j = 0; j < size1; ++j)
                                m(i + Range::first0, j + Range::first1) = other(i, j);

                        return *this;
                    }


                    template <typename Scalar>
                        type & operator += (Scalar const & s)
                        {
                            for(int i = 0; i < size0; ++i)
                                for(int j = 0; j < size1; ++j)
                                    m(i + Range::first0, j + Range::first1) += s;

                            return *this;
                        }

                    template <typename Scalar>
                        type & operator -= (Scalar const & s)
                        {
                            for(int i = 0; i < size0; ++i)
                                for(int j = 0; j < size1; ++j)
                                    m(i + Range::first0, j + Range::first1) -= s;

                            return *this;
                        }
                    template <typename Scalar>
                        type & operator *= (Scalar const & s)
                        {
                            for(int i = 0; i < size0; ++i)
                                for(int j = 0; j < size1; ++j)
                                    m(i + Range::first0, j + Range::first1) *= s;

                            return *this;
                        }

                    template <typename Scalar>
                        type & operator /= (Scalar const & s)
                        {
                            for(int i = 0; i < size0; ++i)
                                for(int j = 0; j < size1; ++j)
                                    m(i + Range::first0, j + Range::first1) /= s;

                            return *this;
                        }

                private:
                    Matrix & m;
                };

            template <typename Vector>
                struct matrix_range_proxy
                {
                };

            template <int First0, int Last0, int First1, int Last1>
                m_range_c<int, First0, Last0, First1, Last1> range()
                {
                    return m_range_c<int, First0, Last0, First1, Last1>();
                }

            template <typename T>
                m_range_d<T> range(T const & first0, T const & last0, T const & first1, T const & last1)
                {
                    return m_range_d<T>(first0, last0, first1, last1);
                }

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_matrix_range_hpp_
