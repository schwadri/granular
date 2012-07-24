/** \file tiny_linalg.hpp
 *  \author Adrian Schweizer (mailto: schwadri@ethz.ch)
 *  \created $Tue 10 May 08:37:28 PM CEST 2011 schwadri@schwadrilaptop.local$
 * \modified $Thu 21 Jul 04:57:28 PM CEST 2011 schwadri@vpn-global-dhcp2-55.ethz.ch$
 *  \description:
 *    this file contains a partial implementation
 *    of vector and matrix algebra operations for dense fixed size vectors and
 *    matrices. Furthermore the quaternion algebra is implemented as well.
 *
 *   implemented are:
 *    for vectors:
 *      u(i)
 *
 *      u + v
 *      s * u
 *
 *      -u
 *      u - v
 *      dot(u, v)
 *      norm_2_sqr(u)
 *      norm_2(u)
 *
 *      A * u
 *      u * A
 *
 *      cout << u
 *
 *    for matrices
 *      A(i, j)
 *
 *      -A
 *
 *      A * B
 *
 *      cout << A
 *    for quaternions:
 *      p(i)
 *
 *      p + q
 *      s * p
 *      p * p
 *      p / q
 *      real(p)
 *      imag(p)
 *      conj(p)
 *      inv(p)
 *
 *      p - q
 *      -p
 *
 *      modulus_sqr(p)  = |p|^2
 *      modulus(p)      = |p|
 *
 *      cout << p
 */

#ifndef _tiny_linalg_hpp_
#define _tiny_linalg_hpp_

#include <cassert>  //using assert
#include <cmath>    //using sqrt
#include <ostream>  //using ostream
#include <limits>   //using numeric_limits

namespace linalg {

  /** \name vector algebra*/ //@{
  template <typename F, int Size>
    struct vector {

      F & operator()(int i) {
        assert(i < Size && i >= 0);
        //FIXME: ugly pointer magic.
        return *(&coordinate + i);
      }

      F const & operator()(int i) const {
        assert(i < Size && i >= 0);
        //FIXME: ugly pointer magic.
        return *(&coordinate + i);
      }

      F & operator[](int i) {
        return (*this)(i);
      }

      F const & operator[](int i) const {
        return (*this)(i);
      }

      F                   coordinate;
      vector<F, Size - 1> tail;
    };

  template <typename F>
    struct vector<F, 1> {

      F & operator()(int i) {
        assert(i == 0);
        return coordinate;
      }

      F const & operator()(int i) const {
        assert(i == 0);
        return coordinate;
      }

      F  coordinate;
    };

  /** \brief add two vectors \p x and \p y together*/
  template <typename F, int Size>
    vector<F, Size> operator + (vector<F, Size> const & x, vector<F, Size> const & y) {
      return (vector<F, Size>){x.coordinate + y.coordinate, x.tail + y.tail};
    }

  template <typename F>
    vector<F, 1> operator + (vector<F, 1> const & x, vector<F, 1> const & y) {
      return (vector<F, 1>){x.coordinate + y.coordinate};
    }

  /** \brief multiply a vector \p x with a \p scalar*/
  template <typename F, int Size>
    vector<F, Size> operator * (F const & scalar, vector<F, Size> const & x) {
      return (vector<F, Size>){scalar * x.coordinate, scalar * x.tail};
    }

  template <typename F>
    vector<F, 1> operator * (F const & scalar, vector<F, 1> const & x) {
      return (vector<F, 1>){scalar * x.coordinate};
    }

  template <typename F, int Size>
    vector<F, Size> operator - (vector<F, Size> const & x) {
      return F(-1.0) * x;
    }

  template <typename F, int Size>
    vector<F, Size> operator - (vector<F, Size> const & x, vector<F, Size> const & y) {
      return x + (-y);
    }

  template <typename F, int Size>
    F dot(vector<F, Size> const & x, vector<F, Size> const & y) {
      return x.coordinate * y.coordinate + dot(x.tail, y.tail);
    }

  template <typename F>
    F dot(vector<F, 1> const & x, vector<F, 1> const & y) {
      return x.coordinate * y.coordinate;
    }

  template <typename F, int Size>
    F norm_2_sqr(vector<F, Size> const & x) {
      return dot(x, x);
    }

  template <typename F, int Size>
    F norm_2(vector<F, Size> const & x) {
      using std::sqrt;
      return norm_2_sqr(x);
    }

  template <typename CharT, typename Traits, typename F, int N>
    std::basic_ostream<CharT, Traits> & operator << (
        std::basic_ostream<CharT, Traits> & out,
        vector<F, N> const & v
    ) {
      out << '[';
        for(int i = 0; i < N - 1; ++i)
          out << v(i) << " ";
      return out << v(N - 1) << ']';
    }

  //@}
  /** \name matrix algebra*/ //@{
  template <typename F, int M, int N>
    struct matrix {

      F const & operator()(int i, int j) const {
        assert(i < M && i >= 0);
        //FIXME: ugly pointer magic.
        return (*(&row + i))(j);
      }

      F & operator()(int i, int j) {
        assert(i < M && i >= 0);
        //FIXME: ugly pointer magic.
        return (*(&row + i))(j);
      }

      vector<F, N>         row;
      matrix<F, M - 1, N>  tail;
    };

  template <typename F, int N>
    struct matrix<F, 1, N> {

      F const & operator()(int i, int j) const {
        assert(i == 0);
        return row(j);
      }

      F & operator()(int i, int j) {
        assert(i == 0);
        return row(j);
      }

      vector<F, N>  row;
    };

  template <typename F, int M, int N>
    matrix<F, M, N> operator * (F const & scalar, matrix<F, M, N> const & a) {
      return (matrix<F, M, N>){scalar * a.row, scalar * a.tail};
    }

  template <typename F, int N>
    matrix<F, 1, N> operator * (F const & scalar, matrix<F, 1, N> const & a) {
      return (matrix<F, 1, N>){scalar * a.row};
    }

  template <typename F, int M, int N>
    matrix<F, M, N> operator - (matrix<F, M, N> const & a) {
      return F(-1.0) * a;
    }

  template <typename F, int M, int K, int N>
    matrix<F, M, N> operator * (matrix<F, M, K> const & a, matrix<F, K, N> const & b) {
      return (matrix<F, M, N>){a.row * b, a.tail * b};
    }

  template <typename F, int K, int N>
    matrix<F, 1, N> operator * (matrix<F, 1, K> const & a, matrix<F, K, N> const & b) {
      return (matrix<F, 1, N>){a.row * b};
    }

  template <typename F, int M, int N>
    vector<F, M> operator * (matrix<F, M, N> const & a, vector<F, N> const & x) {
      return (vector<F, M>){dot(a.row, x), a.tail * x};
    }

  template <typename F, int N>
    vector<F, 1> operator * (matrix<F, 1, N> const & a, vector<F, N> const & x) {
      return (vector<F, 1>){dot(a.row, x)};
    }

  template <typename F, int M, int N>
    vector<F, N> operator * (vector<F, M> const & x, matrix<F, M, N> const & y) {
      return x.coordinate * y.row + x.tail * y.tail;
    }

  template <typename F, int N>
    vector<F, N> operator * (vector<F, 1> const & x, matrix<F, 1, N> const & y) {
      return x.coordinate * y.row;
    }

  template <typename CharT, typename Traits, typename F, int M, int N>
    std::basic_ostream<CharT, Traits> & operator << (
        std::basic_ostream<CharT, Traits> & out,
        matrix<F, M, N> const & m
    ) {
      out << '[';
        for(int i = 0; i < M; ++i) {
          for(int j = 0; j < N - 1; ++j)
            out << m(i, j) << ' ';
          out << m(i, N - 1);
         if(i < M - 1) out << "; ";
        }
      return out << ']';
    }
  //@}

  /** \name quaternion algebra*/ //@{
  template <typename F>
    struct quaternion {

      F const & operator()(int const & i) const {
        assert(i >= 0 && i < 4);
        return coords[i];
      }

      F & operator()(int const & i) {
        assert(i >= 0 && i < 4);
        return coords[i];
      }

      F coords[4];
    };

  template <typename F>
    quaternion<F> conj(quaternion<F> const & p) {
      F const (& a)[4] = p.coords;
      return quaternion<F> {a[0], -a[1], -a[2], -a[3]};
    }

  /** extract the imaginary part*/
  template <typename F>
    quaternion<F> imag(quaternion<F> const & p) {
      F const (& a)[4] = p.coords;
      return quaternion<F> {F(0), a[1], a[2], a[3]};
    }

  /** extract the real part*/
  template <typename F>
    F real(quaternion<F> const & p) {
      F const (& a)[4] = p.coords;
      return /*{*/a[0];//, F(0), F(0), F(0)};
    }

  template <typename F>
    F modulus_sqr(quaternion<F> const & p) {
      F const (& a)[4] = p.coords;
      return a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3];
    }

  /** inv
   */
  template <typename F>
    quaternion<F> inv(quaternion<F> const & p) {
      return (F(1.0) / modulus_sqr(p)) * conj(p);
    }

  template <typename F>
    F modulus(quaternion<F> const & p) {
      using std::sqrt;
      return sqrt(modulus_sqr(p));
    }

  /** quaternion times scalar
   */
  template <typename F>
    quaternion<F> operator * (F const & scalar, quaternion<F> const & p) {
      F const (& a)[4] = p.coords;
      return quaternion<F> {scalar * a[0], scalar * a[1], scalar * a[2], scalar * a[3]};
    }

  /** quaternion + quaternion
  */
  template <typename F>
    quaternion<F> operator + (quaternion<F> const & p, quaternion<F> const & q) {
      F const (& a)[4] = p.coords;
      F const (& b)[4] = q.coords;
      return quaternion<F> {a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]};
    }

  /** quaternion - quaternion
  */
  template <typename F>
    quaternion<F> operator - (quaternion<F> const & p, quaternion<F> const & q) {
      F const (& a)[4] = p.coords;
      F const (& b)[4] = q.coords;
      return quaternion<F> {a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]};
    }

  /** quaternion negation
  */
  template <typename F>
    quaternion<F> operator - (quaternion<F> const & p) {
      F const (& a)[4] = p.coords;
      return quaternion<F> {-a[0], -a[1], -a[2], -a[3]};
    }

  /** quaternion multiplication operator
    \f$ \qp \qmul \qq =
    \left(p^0 q^0 - p^i q^i - p^j q^j - p^k q^k \right) \qid
    + \left(p^0 q^i + p^i q^0 + p^j q^k - p^k q^j \right) \qi
    + \left(p^0 q^j + p^j q^0 + p^k q^i - p^i q^k \right) \qj
    + \left(p^0 q^k + p^k q^0 + p^i q^j - p^j q^i \right) \qk \f$
    */
  template <typename F>
    quaternion<F> operator * (quaternion<F> const & p, quaternion<F> const & q) {
      F const (& a)[4] = p.coords;
      F const (& b)[4] = q.coords;
      return quaternion<F> {
          a[0] * b[0] - a[1] * b[1]
        - a[2] * b[2] - a[3] * b[3],

          a[0] * b[1] + a[1] * b[0]
        + a[2] * b[3] - a[3] * b[2],

          a[0] * b[2] - a[1] * b[3]
        + a[2] * b[0] + a[3] * b[1],

          a[0] * b[3] + a[1] * b[2]
        - a[2] * b[1] + a[3] * b[0]
      };
    }

  /** quaternion division
  */
  template <typename F>
    quaternion<F> operator / (quaternion<F> const & p, quaternion<F> const & q) {
      F const (& a)[4] = p.coords;
      F const (& b)[4] = q.coords;
      F const inv_modulus_sqr = F(1) / modulus_sqr(q);
      return quaternion<F> {
        (- a[0] * b[0] + a[1] * b[1]
         + a[2] * b[2] + a[3] * b[3]) * inv_modulus_sqr,

        (- a[0] * b[1] - a[1] * b[0]
         - a[2] * b[3] + a[3] * b[2]) * inv_modulus_sqr,

        (- a[0] * b[2] + a[1] * b[3]
         - a[2] * b[0] - a[3] * b[1]) * inv_modulus_sqr,

        (- a[0] * b[3] - a[1] * b[2]
         + a[2] * b[1] - a[3] * b[0]) * inv_modulus_sqr
      };
    }

  /** convert a quaternion \p p to a rotation matrix times uniform scaling
   *  matrix. if \p p is a unit quaternion the rotation matrix corresponds
   *  to a pure rotation
   */
  template <typename F>
    matrix<F, 3, 3> to_matrix(quaternion<F> const & p) {
      F const (& a)[4] = p.coords;
      /*
       * p |-] (v |-] map_to_V (p^(0,v)p))
       */
      return matrix<F, 3, 3> {
        //alt: (more efficient, or more stable?)
        //Field(1.0) - Field(2.0) * (a[2] * a[2]) - F(2.0) * (a[3] * a[3]),
        (a[0] * a[0]) + (a[1] * a[1]) - (a[2] * a[2]) - (a[3] * a[3]),
        F(2.0) * (a[1] * a[2]) - F(2.0) * (a[0] * a[3]),
        F(2.0) * (a[1] * a[3]) + F(2.0) * (a[0] * a[2]),

        F(2.0) * (a[1] * a[2]) + F(2.0) * (a[0] * a[3]),

        //alt: (more efficient?)
        //Field(1.0) - Field(2.0) * (a[1] * a[1]) - F(2.0) * (a[3] * a[3]),
        (a[0] * a[0]) - (a[1] * a[1]) + (a[2] * a[2]) - (a[3] * a[3]),
        F(2.0) * (a[2] * a[3]) - F(2.0) * (a[0] * a[1]),

        F(2.0) * (a[1] * a[3]) - F(2.0) * (a[0] * a[2]),
        F(2.0) * (a[2] * a[3]) + F(2.0) * (a[0] * a[1]),

        //alt:
        //Field(1.0) - Field(2.0) * (a[1] * a[1]) - F(2.0) * (a[2] * a[2])
        (a[0] * a[0]) - (a[1] * a[1]) - (a[2] * a[2]) + (a[3] * a[3])
      };
    }

  /** recover a unit quaternion from a rotation matrix
   *  if the matrix does not represent a pure rotation,
   *  the result is undefined
   */
  template <typename F>
    quaternion<F> to_quaternion(matrix<F, 3, 3> const & R) {
      using std::abs; using std::sqrt;
      int uvw[] = {0, 1, 2};

      F dmax = abs(R(0, 0));

      for(int i = 1; i < 3; ++i)
        if(dmax < abs(R(i, i))) {
          uvw[0] = i;
          uvw[1] = (i + 1) % 3;
          uvw[2] = (i + 2) % 3;
          dmax = abs(R(i, i));
        }

      F r = sqrt(F(1) + R(uvw[0], uvw[0]) - R(uvw[1], uvw[1]) - R(uvw[2], uvw[2]));
      F q0;
      F q[3];
      if(r < std::numeric_limits<F>::epsilon()) {
        q0 = F(1); q[0] = F(0); q[1] = F(0); q[2] = F(0);
      }
      else {
        q0 = (R(uvw[2], uvw[1]) - R(uvw[1], uvw[2])) / (F(2) * r);
        q[uvw[0]] = r / F(2);
        q[uvw[1]] = (R(uvw[0], uvw[1]) + R(uvw[1], uvw[0])) / (F(2) * r);
        q[uvw[2]] = (R(uvw[2], uvw[0]) + R(uvw[0], uvw[2])) / (F(2) * r);
      }
      return quaternion<F> {q0, q[0], q[1], q[2]};
    }

  template <typename CharT, typename Traits, typename F>
    std::basic_ostream<CharT, Traits> & operator << (
      std::basic_ostream<CharT, Traits> & out,
      quaternion<F> & p
    ) {
      F const (& a)[4] = p.coords;
      return out << '(' << a[0] << ',' << a[1] << ',' << a[2] << ',' << a[3] << ')';
    }

  //@}


} //namespace linalg

#endif // _tiny_linalg_hpp_
