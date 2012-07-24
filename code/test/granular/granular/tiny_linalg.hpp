/** \file tiny_linalg.hpp
 *  \author Adrian Schweizer (mailto: schwadri@ethz.ch)
 *  \created $Tue 10 May 08:37:28 PM CEST 2011 schwadri@schwadrilaptop.local$
 * \modified $Tue 30 Aug 12:31:41 PM CEST 2011 schwadri@SchwadriBook.local$
 *  \description:
 *    this file contains a partial implementation
 *    of vector and matrix algebra operations for dense fixed size vectors and
 *    matrices. Furthermore the quaternion algebra is implemented as well.
 *
 *   implemented are:
 *    for vectors:
 *      u[i]
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
 *      for 3 vectors
 *      cross
 *
 *      A * u
 *
 *      cout << u
 *
 *    for matrices
 *      A[i][j]
 *
 *      A + B
 *      A - B
 *      s * A
 *
 *      -A
 *
 *      A * B
 *
 *      for 3x3 matrices:
 *      det
 *      transpose
 *      inverse
 *
 *      cout << A
 *    for quaternions:
 *      p(i)
 *
 *      p + q
 *      s * p
 *      p * q
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

#include <algorithm>  //using transform
#include <numeric>    //using inner_product
#include <functional> //using plus, times
#include <cassert>    //using assert
#include <cmath>      //using sqrt
#include <ostream>    //using ostream
#include <iterator>   //using ostream_iterator
#include <limits>     //using numeric_limits

namespace linalg {

  /** \name vector algebra*/ //@{
  template <typename F, int Size>
    struct vector {
      F & operator[](int i) {
        assert(i < Size && i >= 0);
        return coords[i];
      }

      F const & operator[](int i) const {
        assert(i < Size && i >= 0);
        return coords[i];
      }
      F coords[Size];
    };

  /** \brief add two vectors \p x and \p y together*/
  template <typename F, int Size>
    vector<F, Size> operator + (vector<F, Size> const & a, vector<F, Size> const & b) {
      vector<F, Size> r;
      using std::transform; using std::plus;
      transform(
        a.coords,
        a.coords + Size,
        b.coords,
        r.coords,
        plus<F>()
      );
      return r;
    }

  /** \brief multiply a vector \p x with a \p scalar*/
  template <typename F, int Size>
    vector<F, Size> operator * (F const & s, vector<F, Size> const & x) {
      vector<F, Size> r;
      using std::transform; using std::multiplies; using std::bind1st;
      transform(
        x.coords,
        x.coords + Size,
        r.coords,
        bind1st(multiplies<F>(), s)
      );
      return r;
    }

  /** \brief negate the vector. This is a shortcut for (-1) * v */
  template <typename F, int Size>
    vector<F, Size> operator - (vector<F, Size> const & x) {
      vector<F, Size> r;
      using std::transform; using std::negate;
      transform(
        x.coords,
        x.coords + Size,
        r.coords,
        negate<F>()
      );
      return r;
    }

  /** \brief substract vector \p b from \p a. This is a shortcut for a + (-1) * b */
  template <typename F, int Size>
    vector<F, Size> operator - (vector<F, Size> const & a, vector<F, Size> const & b) {
      vector<F, Size> r;
      using std::transform; using std::minus;
      transform(
        a.coords,
        a.coords + Size,
        b.coords,
        r.coords,
        minus<F>()
      );
      return r;
    }

  /** \brief determine the inner product of \p a and \p b*/
  template <typename F, int Size>
    F dot(vector<F, Size> const & a, vector<F, Size> const & b) {
      using std::inner_product;
      return inner_product(
        a.coords,
        a.coords + Size,
        b.coords,
        F(0)
      );
    }
  /** \brief determine the cross product of \p a and \p b*/
  template <typename F>
    inline vector<F, 3> cross(vector<F, 3> const & a, vector<F, 3> const & b) {
      vector<F, 3> r = (vector<F, 3>){
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
      };
      return r;
    }

  template <typename F, int Size>
    F norm_2_sqr(vector<F, Size> const & x) {
      return dot(x, x);
    }

  template <typename F, int Size>
    F norm_2(vector<F, Size> const & x) {
      using std::sqrt;
      return sqrt(norm_2_sqr(x));
    }

  template <typename CharT, typename Traits, typename F, int N>
    std::basic_ostream<CharT, Traits> & operator << (
        std::basic_ostream<CharT, Traits> & out,
        vector<F, N> const & v
    ) {
      out << '[';
      using std::copy; using std::ostream_iterator;
      copy(
        v.coords,
        v.coords + (N - 1),
        ostream_iterator<F>(out, " ")
      );
      return out << v.coords[N - 1] << ']';
    }

  //@}
  /** \name matrix algebra*/ //@{
  template <typename F, int M, int N>
    struct matrix {

      vector<F, N> const & operator[](int i) const {
        assert(i < M && i >= 0);
        return rows[i];
      }

      vector<F, N> & operator[](int i) {
        assert(i < M && i >= 0);
        return rows[i];
      }

      vector<F, N>  rows[M];
    };

  /** \brief add to matrices */
  template <typename F, int M, int N>
    matrix<F, M, N> operator + (matrix<F, M, N> const & a, matrix<F, M, N> const & b) {
      matrix<F, M, N> r;
      using std::transform; using std::plus;
      transform(
        a.rows,
        a.rows + M,
        b.rows,
        r.rows,
        plus<vector<F, N> >()
      );
      return r;
    }

  template <typename F, int M, int N>
    matrix<F, M, N> operator - (matrix<F, M, N> const & a, matrix<F, M, N> const & b) {
      matrix<F, M, N> r;
      using std::transform; using std::minus;
      transform(
        a.rows,
        a.rows + M,
        b.rows,
        r.rows,
        minus<vector<F, N> >()
      );
      return r;
    }

  template <typename S, typename V>
    struct scalar_product : std::binary_function<S, V, V> {
      V operator()(S const & s, V const & v) const {
        return s * v;
      }
    };

  template <typename F, int M, int N>
    matrix<F, M, N> operator * (F const & s, matrix<F, M, N> const & a) {
      matrix<F, M, N> r;
      using std::transform; using std::multiplies; using std::bind1st;
      transform(
        a.rows,
        a.rows + M,
        r.rows,
        bind1st(scalar_product<F, vector<F, N> >(), s)
      );
      return r;
    }

  template <typename F, int M, int N>
    matrix<F, M, N> operator - (matrix<F, M, N> const & a) {
      return F(-1.0) * a;
    }

  template <typename F, int M, int K, int N>
    matrix<F, M, N> operator * (matrix<F, M, K> const & a, matrix<F, K, N> const & b) {
      matrix<F, M, N> r;
      for(int i = 0; i < M; ++i)
        for(int j = 0; j < N; ++j) {
          int s = F(0);
          for(int k = 0; k < K; ++k)
            s = s + a[i][k] * b[k][j];
          r[i][j] = s;
        }
      return r;
    }

  template <typename F, int M, int N>
    vector<F, M> operator * (matrix<F, M, N> const & a, vector<F, N> const & x) {
      vector<F, M> r;
      using std::transform; using std::bind2nd;
      for(int i = 0; i < M; ++i) {
          F s(0);
        for(int j = 0; j < N; ++j)
          s = s + a[i][j] * x[j];
        r[i] = s;
      }
      return r;
    }

  template <typename F>
    matrix<F, 3, 3> inverse(matrix<F, 3, 3> const & m) {
      F m00 = (m[1][1] * m[2][2] - m[2][1] * m[1][2]),
        m10 = (m[0][1] * m[2][2] - m[2][1] * m[0][2]),
        m20 = (m[0][1] * m[1][2] - m[1][1] * m[0][2]);

      F det =   m[0][0] * m00
              - m[1][0] * m10
              + m[2][0] * m20;

      F inv_det = F(1) / det;

      matrix<F, 3, 3> r = {
        m00 * inv_det, -m10 * inv_det, m20 * inv_det,
        -(m[1][0] * m[2][2] - m[2][0] * m[1][2]) * inv_det,  (m[0][0] * m[2][2] - m[2][0] * m[0][2]) * inv_det, -(m[0][0] * m[1][2] - m[1][0] * m[0][2]) * inv_det,
         (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * inv_det, -(m[0][0] * m[2][1] - m[2][0] * m[0][1]) * inv_det,  (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * inv_det
      };
      return r;
    }

  template <typename F>
    F det(matrix<F, 3, 3> const & a) {
      return
        + a[0][0] * a[1][1] * a[2][2]
        + a[1][0] * a[2][1] * a[0][2]
        + a[2][0] * a[0][1] * a[1][2]
        - a[0][0] * a[2][1] * a[1][2]
        - a[1][0] * a[0][1] * a[2][2]
        - a[2][0] * a[1][1] * a[0][2]
      ;
    }

  template <typename F>
    matrix<F, 3, 3> transpose(matrix<F, 3, 3> const & a) {
      matrix<F, 3, 3> m = {
        a[0][0], a[1][0], a[2][0],
        a[0][1], a[1][1], a[2][1],
        a[0][2], a[1][2], a[2][2]
      };
      return m;
    }

  template <typename CharT, typename Traits, typename F, int M, int N>
    std::basic_ostream<CharT, Traits> & operator << (
        std::basic_ostream<CharT, Traits> & out,
        matrix<F, M, N> const & m
    ) {
      out << '[';
        for(int i = 0; i < M; ++i) {
          for(int j = 0; j < N - 1; ++j)
            out << m[i][j] << ' ';
          out << m[i][N - 1];
         if(i < M - 1) out << "; ";
        }
      return out << ']';
    }
  //@}

  /** \name quaternion algebra*/ //@{
  template <typename F>
    struct quaternion {

      F const & operator[](int const & i) const {
        assert(i >= 0 && i < 4);
        return coords[i];
      }

      F & operator[](int const & i) {
        assert(i >= 0 && i < 4);
        return coords[i];
      }

      F coords[4];
    };

  template <typename F>
    quaternion<F> conj(quaternion<F> const & p) {
      F const (& a)[4] = p.coords;
      quaternion<F> r = {a[0], -a[1], -a[2], -a[3]};
      return r;
    }

  /** extract the imaginary part*/
  template <typename F>
    quaternion<F> imag(quaternion<F> const & p) {
      F const (& a)[4] = p.coords;
      quaternion<F> r = {F(0), a[1], a[2], a[3]};
      return r;
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
      quaternion<F> r = {scalar * a[0], scalar * a[1], scalar * a[2], scalar * a[3]};
      return r;
    }

  /** quaternion + quaternion
  */
  template <typename F>
    quaternion<F> operator + (quaternion<F> const & p, quaternion<F> const & q) {
      F const (& a)[4] = p.coords;
      F const (& b)[4] = q.coords;
      quaternion<F> r = {a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]};
      return r;
    }

  /** quaternion - quaternion
  */
  template <typename F>
    quaternion<F> operator - (quaternion<F> const & p, quaternion<F> const & q) {
      F const (& a)[4] = p.coords;
      F const (& b)[4] = q.coords;
      quaternion<F> r = {a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]};
      return r;
    }

  /** quaternion negation
  */
  template <typename F>
    quaternion<F> operator - (quaternion<F> const & p) {
      F const (& a)[4] = p.coords;
      quaternion<F> r = {-a[0], -a[1], -a[2], -a[3]};
      return r;
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
      quaternion<F> r = {
          a[0] * b[0] - a[1] * b[1]
        - a[2] * b[2] - a[3] * b[3],

          a[0] * b[1] + a[1] * b[0]
        + a[2] * b[3] - a[3] * b[2],

          a[0] * b[2] - a[1] * b[3]
        + a[2] * b[0] + a[3] * b[1],

          a[0] * b[3] + a[1] * b[2]
        - a[2] * b[1] + a[3] * b[0]
      };
      return r;
    }

  /** quaternion division
  */
  template <typename F>
    quaternion<F> operator / (quaternion<F> const & p, quaternion<F> const & q) {
      F const (& a)[4] = p.coords;
      F const (& b)[4] = q.coords;
      F const inv_modulus_sqr = F(1) / modulus_sqr(q);
      quaternion<F> r = {
        (- a[0] * b[0] + a[1] * b[1]
         + a[2] * b[2] + a[3] * b[3]) * inv_modulus_sqr,

        (- a[0] * b[1] - a[1] * b[0]
         - a[2] * b[3] + a[3] * b[2]) * inv_modulus_sqr,

        (- a[0] * b[2] + a[1] * b[3]
         - a[2] * b[0] - a[3] * b[1]) * inv_modulus_sqr,

        (- a[0] * b[3] - a[1] * b[2]
         + a[2] * b[1] - a[3] * b[0]) * inv_modulus_sqr
      };
      return r;
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
      matrix<F, 3, 3> r = {
        //alt: (more efficient, or more stable?)
        (a[0] * a[0]) + (a[1] * a[1]) - (a[2] * a[2]) - (a[3] * a[3]),
        F(2.0) * (a[1] * a[2]) - F(2.0) * (a[0] * a[3]),
        F(2.0) * (a[1] * a[3]) + F(2.0) * (a[0] * a[2]),

        F(2.0) * (a[1] * a[2]) + F(2.0) * (a[0] * a[3]),

        //alt: (more efficient?)
        (a[0] * a[0]) - (a[1] * a[1]) + (a[2] * a[2]) - (a[3] * a[3]),
        F(2.0) * (a[2] * a[3]) - F(2.0) * (a[0] * a[1]),

        F(2.0) * (a[1] * a[3]) - F(2.0) * (a[0] * a[2]),
        F(2.0) * (a[2] * a[3]) + F(2.0) * (a[0] * a[1]),

        //alt:
        (a[0] * a[0]) - (a[1] * a[1]) - (a[2] * a[2]) + (a[3] * a[3])
      };
      return r;
    }

  /** recover a unit quaternion from a rotation matrix
   *  if the matrix does not represent a pure rotation,
   *  the result is undefined
   */
  template <typename F>
    quaternion<F> to_quaternion(matrix<F, 3, 3> const & R) {
      using std::abs; using std::sqrt;
      int uvw[] = {0, 1, 2};

      F dmax = abs(R[0][0]);

      for(int i = 1; i < 3; ++i)
        if(dmax < abs(R[i][i])) {
          uvw[0] = i;
          uvw[1] = (i + 1) % 3;
          uvw[2] = (i + 2) % 3;
          dmax = abs(R[i][i]);
        }

      F r = sqrt(F(1) + R[uvw[0]][uvw[0]] - R[uvw[1]][uvw[1]] - R[uvw[2]][uvw[2]]);
      F q0;
      F q[3];
      if(r < std::numeric_limits<F>::epsilon()) {
        q0 = F(1); q[0] = F(0); q[1] = F(0); q[2] = F(0);
      }
      else {
        q0 = (R[uvw[2]][uvw[1]] - R[uvw[1]][uvw[2]]) / (F(2) * r);
        q[uvw[0]] = r / F(2);
        q[uvw[1]] = (R[uvw[0]][uvw[1]] + R[uvw[1]][uvw[0]]) / (F(2) * r);
        q[uvw[2]] = (R[uvw[2]][uvw[0]] + R[uvw[0]][uvw[2]]) / (F(2) * r);
      }
      quaternion<F> rs = {q0, q[0], q[1], q[2]};
      return rs;
    }

  template <typename CharT, typename Traits, typename F>
    std::basic_ostream<CharT, Traits> & operator << (
      std::basic_ostream<CharT, Traits> & out,
      quaternion<F> & p
    ) {
      F const (& a)[4] = p.coords;
      return out << '[' << a[0] << ' ' << a[1] << ' ' << a[2] << ' ' << a[3] << ']';
    }

  //@}


} //namespace linalg

#endif // _tiny_linalg_hpp_
