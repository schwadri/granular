#ifndef _core_math_linalg_vector_range_hpp_
#define _core_math_linalg_vector_range_hpp_

/** \file range.hpp
 *  \author Adrian Schweizer
 *  \created  $Fr 18 Jul 10:49:53 am CEST 2008 schwadri@SchwadriLaptop.local$
 *  \modified $Di 29 Jul 11:35:52 am CEST 2008 schwadri@guest-docking-hg-1-155.ethz.ch$
 */

#include <algorithm>
#include <functional>

#include <boost/static_assert.hpp>

#include "../traits/field.hpp"

#include "fwd.hpp"
#include "init.hpp"

namespace core {

    namespace math {

        namespace linalg {

            template <typename T, T First, T Last>
                struct v_range_c
                {
                    typedef T   value_type;
                    static const T first    = First;
                    static const T last     = Last;
                };

            template <typename T>
                struct v_range_d
                {
                    v_range_d(T const & first_, T const & last_)
                    :   first(first_),
                        last(last_)
                    { }

                    T first,
                      last;
                };

            template <typename Vector, typename Range>
                struct vector_static_range_proxy
                {
                    typedef vector_static_range_proxy                       type;
                    typedef typename Vector::value_type                     value_type;
                    typedef Range                                           range;
                    static const int                                        size = Range::last - Range::first;
                    typedef typename Vector::iterator                       iterator;
                    typedef typename Vector::const_iterator                 const_iterator;
                    typedef typename Vector::template rebind<size>::type    vector_type;

                    vector_static_range_proxy(Vector & v_)
                    :   v(v_)
                    { }

                    value_type & operator()(int i)
                    {
                        return v[i + Range::first];
                    }

                    value_type & operator[](int i)
                    {
                        return v[i + Range::first];
                    }

                    value_type const & operator()(int i) const
                    {
                        return v[i + Range::first];
                    }

                    value_type const & operator[](int i) const
                    {
                        return v[i + Range::first];
                    }

                    value_type & front()
                    {
                        return v[Range::first];
                    }

                    value_type const & front() const
                    {
                        return v[Range::first];
                    }

                    value_type & back()
                    {
                        return v[Range::last - 1];
                    }

                    value_type const & back() const
                    {
                        return v[Range::last - 1];
                    }

                    template <typename T, T First, T Last>
                        vector_static_range_proxy<
                            Vector,
                            v_range_c<
                                T,
                                First + Range::first,
                                Last + Range::first
                            >
                        > operator[](v_range_c<T, First, Last> (*hint)())
                        {
                            BOOST_STATIC_ASSERT((size > First && size >= Last));
                            return vector_static_range_proxy<
                                Vector,
                                v_range_c<
                                    T,
                                    First + Range::first,
                                    Last + Range::first
                                >
                            >(v);
                        }

                    template <typename T, T First, T Last>
                        vector_static_range_proxy<
                            Vector const,
                            v_range_c<
                                T,
                                First + Range::first,
                                Last + Range::first
                            >
                        > operator[](v_range_c<T, First, Last> (*hint)()) const
                        {
                            BOOST_STATIC_ASSERT((size > First && size >= Last));
                            return vector_static_range_proxy<
                                Vector const,
                                v_range_c<
                                    T,
                                    First + Range::first,
                                    Last + Range::first
                                >
                            >(v);
                        }

                    iterator begin()
                    {
                        return v.begin() + Range::first;
                    }

                    iterator end()
                    {
                        return v.begin() + Range::last;
                    }

                    const_iterator begin() const
                    {
                        return v.begin() + Range::first;
                    }

                    const_iterator end() const
                    {
                        return v.begin() + Range::last;
                    }

                    template <typename Scalar>
                        vector_init<type, 1> operator=(Scalar const & value)
                        {
                            v[Range::first] = value;
                            return vector_init<type, 1>(*this);
                        }

                    template <typename T>
                        type & operator=(vector<T, size> const & v)
                        {
                            using std::copy;
                            copy(
                                v.begin(),
                                v.end(),
                                begin()
                            );

                            return *this;
                        }

                    template <typename Scalar>
                        type & operator *= (Scalar const & s)
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

                            return *this;
                        }

                    template <typename Scalar>
                        type & operator += (Scalar const & s)
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

                            return *this;
                        }

                    template <typename Scalar>
                        type & operator -= (Scalar const & s)
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

                            return *this;
                        }

                    template <typename Scalar>
                        type & operator /= (Scalar const & s)
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

                            return *this;
                        }

                    type & operator += (vector_type const & other)
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

                        return *this;
                    }

                    type & operator -= (vector_type const & other)
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

                        return *this;
                    }

                    type & operator *= (vector_type const & other)
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

                        return *this;
                    }

                    type & operator /= (vector_type const & other)
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

                        return *this;
                    }

                    vector_type operator-() const
                    {
                        vector_type result;
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
                    Vector & v;
                };

            template <typename S, typename Vector, typename Range>
                typename vector_static_range_proxy<Vector, Range>::vector_type operator * (vector_static_range_proxy<Vector, Range> const & vrp, S const & s)
                {
                    typename vector_static_range_proxy<Vector, Range>::vector_type r;
                    using std::copy;

                    copy(vrp.begin(), vrp.end(), r.begin());
                    r *= s;
                    return r;
                }

            template <typename S, typename Vector, typename Range>
                typename vector_static_range_proxy<Vector, Range>::vector_type operator * (S const & s, vector_static_range_proxy<Vector, Range> const & vrp)
                {
                    typename vector_static_range_proxy<Vector, Range>::vector_type r;
                    using std::copy;

                    copy(vrp.begin(), vrp.end(), r.begin());
                    r *= s;
                    return r;
                }

            template <typename T, typename Vector, typename Range>
                vector<T, 3> cross(vector_static_range_proxy<Vector, Range> const & v1, vector<T, 3> const & v2)
                {
                    return vector<T, 3>(
                        v1[1]*v2[2] - v1[2]*v2[1],
                        v1[2]*v2[0] - v1[0]*v2[2],
                        v1[0]*v2[1] - v1[1]*v2[0]
                    );
                }

            //TODO:
            template <typename Vector, typename Range>
                struct vector_range_proxy
                {
                    vector_range_proxy(Vector & v_, Range const & rng_)
                    :   v(v_),
                        rng(rng_)
                    { }

                private:
                    Vector &    v;
                    Range       rng;
                };

            template <int First, int Last>
                v_range_c<int, First, Last> range()
                {
                    return v_range_c<int, First, Last>();
                }

            //TODO:
            template <typename T>
                v_range_d<T> range(T const & first, T const & last)
                {
                    return v_range_d<T>(first, last);
                }

            template <typename V, typename R>
                struct field<vector_static_range_proxy<V, R> >
                :   field<V>
                { };

        } // namespace linalg
    } // namespace math
} // namespace core

#endif // _core_math_linalg_vector_range_hpp_
