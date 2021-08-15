
#include <iostream>
#include <iomanip>

#include <boost/math/special_functions/factorials.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>

using namespace boost;
using namespace accumulators;

using boost::math::factorial;

template <class Float1>
typename boost::math::tools::promote_args<Float1>::type
LimLi(const int &n, const Float1 &x) {

  typedef typename boost::math::tools::promote_args<Float1>::type Type;

  if(x < 0) return Type(0);
  else if(n != 0) return -pow(x, n) / factorial<Type>(n);
  else if(x > 0) return Type(-1);
  else return Type(-0.5);
}


#include "oldFunctions.hpp"

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
mInt0to1LimLim1(const int &m, const Float1 &a, Float2 d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  if(a == 0) {
    std::cout << "Something is wrong! The function mInt0to1LimLim1 can not be called with a = 0" << std::endl;
    abort();
  }

  d = -d / a;
  if(d > 0 && d < 1) {
    return pow(d, m) / fabs(a);
  }
  else if(d == 0 || d == 1) {
    return (m == 0) ? 0.5 / fabs(a) : 0.5 * pow(d, m) / fabs(a);
  }
  else {
    return Type(0);
  }
}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
mInt0to1LimLi0(const int &m, const Float1 &a, Float2 d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  if(a == 0) {
    return -LimLi(0, d) / Type(m + 1);
  }

  d = -d / a;

  if(a > 0) {
    if(d <= 0) return 1 / Type(m + 1);
    else if(d < 1) return (Type(1) - pow(d, m + 1)) / Type(m + 1);
    else return Type(0);
  }
  else if(a < 0) {
    if(d >= 1) return 1 / Type(m + 1);
    else if(d > 0) pow(d, m + 1) / Type(m + 1);
    else return 0;
  }
}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
mInt0to1LimLi(const int &s, const unsigned &m, Float1 a, Float2 d) {
  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  switch(s) {
    case -1:
      return mInt0to1LimLim1(m, a, d);
      break;

    case 0:
      return mInt0to1LimLi0(m, a, d);
      break;

    default:

      if(a == 0) {
        return -LimLi(s, d) / Type(m + 1);
      }

      Type INT(0);
      Type x(a + d);

      if(x < 0 || d < 0) { // in all these cases no-significant digit cancellation occurs
        if(x >= 0) {
          Type g =  1 / (-a);
          for(unsigned i = 1; i <= m + 1; i++) {
            INT += g * LimLi(s + i, x);
            g *= (m + 1 - i) / (-a);
          }
        }
        else if(d >= 0.) {
          INT += factorial<Type>(m) * LimLi(s + m + 1, d) / (a * pow(-a, m));
        }
      }
      else { //alternative formula to avoid significant digit cancellations when s>1, and (a+d) and d are non-negative and d >> a
        for(int i = 1; i <= s; i++) {
          INT -= pow(-a, s - i) / factorial<Type>(m + 1 + s - i) * LimLi(i, x) ;
        }
        INT += pow(-a, s) / factorial<Type>(m + 1 + s);
        INT *= factorial<Type>(m);
      }
      return INT;
  }
}


template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
LSIm1(const int &m, const Float1 &a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  if(a == 0) {
    std::cout << "Something is wrong! The function LSIm1 can not be called with a = 0" << std::endl;
    abort();
  }

  if(fabs(d / a) < 1) {
    return (m == 0) ? Type(1) / fabs(a) : pow(-d / a, m) / fabs(a);
  }
  else if(fabs(d / a) == 1) {
    return 0.5 * pow(-d / a, m) / fabs(a);
  }
  else {
    return Type(0);
  }
}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
LSI0(const int &m, const Float1 &a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  if(a == 0) {
    return -LimLi(0, d) * (1 + pow(-1, m)) / Type(m + 1);
  }

  if(d <= -fabs(a)) {
    return 0.;
  }
  else if(d >= fabs(a)) {
    return (1. + pow(-1., m)) / Type(m + 1);
  }
  else {
    if(a > 0) {
      return (1. - pow(-d / a, m + 1)) / Type(m + 1);
    }
    else {
      return (pow(-d / a, m + 1) + pow(-1., m)) / Type(m + 1);
    }
  }
}


template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
LSI(const int &s, const unsigned &m, Float1 a, Float2 d) {
  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  switch(s) {
    case -1:
      return LSIm1(m, a, d);
      break;

    case 0:
      return LSI0(m, a, d);
      break;

    default:

      if(a == 0) {
        return -LimLi(s, d) * (1 + pow(-1, m)) / Type(m + 1);
      }

      Type INT(0);

      Type x(a + d);
      Type y(-a + d);

      if(d < fabs(a)) { // in all these cases no-significant digit cancellation occurs
        if(x > 0.) {
          Type xDa = x / a;
          Type c1 = s + 1;
          Type c2 = m + 1;
          Type fx =  pow(x, c1) / (a * factorial<Type>(s + 1));
          for(unsigned i = 0; i < m + 1; i++) {
            INT += fx;
            fx *= - xDa * (--c2) / (++c1);
          }
        }
        else if(y > 0.) {
          Type yDa = y / a;
          Type c1 = s + 1;
          Type c2 = m + 1;
          Type fy =  - pow(-1., m) * pow(y, c1) / (a * factorial<Type>(s + 1));
          for(unsigned i = 0; i < m + 1; i++) {
            INT += fy;
            fy *= yDa * (--c2) / (++c1);
          }
        }
      }
      else { //alternative formula to avoid significant digit cancellations when s>1, and (a+d) and (-a+d) are non-negative and d >> a

        Type px = 1.;
        Type py = pow(-1, m + s);

        Type c1 = m + 1 + s;
        Type c2 = 1;
        Type f1 = pow(-a, s) / factorial<Type>(m + 1 + s);

        for(int i = 0; i < s + 1; i++) {
          INT += f1 * (px + py);
          f1 *= -(c1--) / (a * c2++);
          px *= x;
          py *= -y;
        }
        INT *= factorial<Type>(m);
      }
      return INT;
      break;
  }
}





template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
TriangleReduced(const int &s, const std::vector<unsigned> &m_input, const std::vector <Float1> &a_input, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];

  Type TRI = 0.;
  if(b == 0) { //parallel to right edge
    TRI = (-LimLi(s + 1, a + d) +
           factorial<Type>(m + n + 1) * pow(-1. / a, m + n + 1) * LimLi(s + m + n + 2, d)
          ) / (a * (n + 1.));
  }
  else if(a == 0) { //parallel to bottom edge
    TRI = (factorial<Type>(n) * pow(-1. / b, n) * LimLi(s + n + 1, d) -
           factorial<Type>(m + n + 1) * pow(-1. / b, m + n + 1) * LimLi(s + m + n + 2, d)
          ) / (b * (m + 1.));
  }
  else if(a + b == 0) { // parallel to y=x edge//TODO
    if(a + d  >= 0.) {
      for(unsigned i = 1; i <= m + 1; i++) {
        TRI += pow(-1. / a, i) * LimLi(s + n + 1 + i, a + d) / factorial<Type>(m + 1 - i);
      }
      TRI *= factorial<Type>(n) * factorial<Type>(m) * pow(1. / a, n + 1);
    }
    TRI += LimLi(s + 1, d) / (a * (m + n + 1));
  }
  else { // general case
    if(fabs(b) > fabs(a)) {
      if(d > 0.) {
        for(unsigned i = 1; i <= n + 1; i++) {
          TRI += factorial<Type>(m + n + 1 - i) / factorial<Type>(n + 1 - i) * pow((a + b) / b, i);
        }
        TRI *= LimLi(s + m + n + 2, d) / pow(-a - b, m + n + 2) ;
      }
      TRI += mInt0to1LimLi(s + n + 1, m, a, d) * pow(-1. / b, n + 1);
      TRI *= factorial<Type>(n);
    }
    else {
      if(d > 0.) {
        for(unsigned i = 1; i <= m + 1; i++) {
          TRI += factorial<Type>(m + n + 1 - i) / factorial<Type>(m + 1 - i) * pow((a + b) / a, i);
        }
        TRI *= - LimLi(s + m + n + 2, d) / pow(-a - b, m + n + 2) ;
      }
      if(a + d > 0.) {
        Type TRI2 = 0.;
        for(unsigned i = 1; i <= m + 1; i++) {
          TRI2 +=  LimLi(s + n + i + 1, a + d) / (factorial<Type>(m + 1 - i) * pow(-a, i));
        }
        TRI += TRI2 * factorial<Type>(n) / pow(-b, n + 1) ;
      }
      TRI *= factorial<Type>(m);
    }

  }
  return TRI;
}



template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
TriangleFull(const int &s, const std::vector<unsigned> &m_input, const std::vector <Float1> &a_input, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];

  Type TRI = 0.;
  //std::cout << "FULL " << a << " " << b << " " << d << " ";
  if(b == 0) { //parallel to right edge
    if(a == 0) TRI = -LimLi(s, d) / ((m + n + 2) * (n + 1));
    else TRI = mInt0to1LimLi(s, m + n + 1, a, d) / (n + 1.);
  }
  else if(a == 0) { //parallel to bottom edge
    TRI = (mInt0to1LimLi(s, n, b, d) - mInt0to1LimLi(s, m + n + 1, b, d)) / (m + 1.);
  }
  else if(a + b == 0) { //parallel to y = x edge //TODO
    for(unsigned i = 1; i <= n + 1; i++) {
      TRI += pow(1. / a, i) * LimLi(s + i, d) / (factorial<Type>(n + 1 - i) * (m + n + 2 - i));
    }
    TRI += mInt0to1LimLi(s + n + 1, m, a, d) / pow(a, n + 1);
    TRI *= factorial<Type>(n);
  }
  else { //generic case
    if(fabs(a) < fabs(b)) { // |a| < |b|
      for(unsigned j = 1; j <= n + 1; j++) {
        TRI -= mInt0to1LimLi(s + j, m + n + 1 - j, a + b, d) * pow(-1. / b, j)
               / factorial<Type>(n + 1 - j);
      }
      TRI += mInt0to1LimLi(s + n + 1, m, a, d) * pow(-1. / b, n + 1);
      TRI *= factorial<Type>(n);
    }
    else { // |b| < |a|
      for(unsigned j = 1; j <= m + 1; j++) {
        TRI += (mInt0to1LimLi(s + j, m + n + 1 - j, a + b, d) -
                mInt0to1LimLi(s + j, n, b, a + d)) * pow(-1. / a, j)
               / factorial<Type>(m + 1 - j);
      }
      TRI *= factorial<Type>(m);
    }
  }

  return TRI;
}



template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
TriangleA(const int &s, const std::vector<unsigned> &m, const std::vector <Float1> &a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  if(fabs(a[0]) <= 1.0e-10 && fabs(a[1]) <= 1.0e-10) return -LimLi(s, d) / ((m[0] + m[1] + 2) * (m[1] + 1)); // only from higher dimensions calls n = <0,0,c>


  switch(s) {
    case -1:
      if(a[0] + a[1] + d <= 0) {
        return TriangleReduced(-1, m, a, d);
      }
      else return TriangleReduced(-1, m, std::vector<Type> {-a[0], -a[1]}, -d);
      break;
    case 0:
      if(a[0] + a[1] + d <= 0) {
        return TriangleReduced(0, m, a, d);
      }
      else return Type(1) / ((1 + m[1]) * (2 + m[1] + m[0])) - TriangleReduced(0, m, std::vector<Type> {-a[0], -a[1]}, -d);
      break;
    default:

      if(a[0] + a[1] + d <= 0) {
        return TriangleReduced(s, m, a, d);
      }
      else if(a[0] + a[1] + d > 0 && d > 0 && a[1] > 0 && d / a[1] < 0.1 && fabs((a[0] + a[1]) / (a[0] - a[1])) < 0.015) {
        Type INT = TriangleFull(s, m, a, d);
        //std::cout<< INT << " ";
        return INT;
      }
      else {

        accumulator_set<Type, stats<tag::sum_kahan> > acc;


        Type INT(0.);
        Type m1f = factorial<Type>(m[1]);
        INT += m1f / factorial<Type>(m[1] + s) * pow(-a[1], s)
               * TriangleA(0, std::vector<unsigned> {m[0], m[1] + s}, a, d) ;

        acc(m1f / factorial<Type>(m[1] + s) * pow(-a[1], s)
            * TriangleA(0, std::vector<unsigned> {m[0], m[1] + s}, a, d)) ;

        //std::cout << INT << " ";

        for(int i = s - 1; i >= 0; i--) {
          INT += m1f / factorial<Type>(m[1] + 1 + i) * pow(-a[1], i)
                 * mInt0to1LimLi(s - i, m[0] + m[1] + 1 + i, a[0] + a[1], d);
          acc(m1f / factorial<Type>(m[1] + 1 + i) * pow(-a[1], i)
              * mInt0to1LimLi(s - i, m[0] + m[1] + 1 + i, a[0] + a[1], d));
          //std::cout << INT << " ";
        }
        //std::cout << INT << "\n";
        return sum_kahan(acc);
      }
  }
}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Triangle(const int &s, const std::vector<unsigned> &m, const std::vector <Float1> &a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  return TriangleA(s, m,  std::vector<Type> {-a[0], a[1]}, d + a[0]);
}


template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
TetrahedronB(const int &s, const std::vector<unsigned> &m_input, const std::vector <Float1> &a_input, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const Type &c = a_input[2];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];
  const unsigned &o = m_input[2];

  Type TET = 0.;
  if(c == 0) {  //plane is parallel to z-axis
    for(unsigned i = 0; i <= o + 1; i++)  {
      TET += pow(-1., i) / (factorial<Type>(i) * factorial<Type>(o + 1 - i)) *
             TriangleA(s, std::vector<unsigned> {m + o + 1 - i, n + i}, std::vector<Type> {a, b}, d);
    }
    TET *= factorial<Type>(o);
  }
  else { //all other cases
    for(unsigned i = 1; i <= o + 1; i++)  {
      Type TETi = 0.;
      for(unsigned j = 0; j <= o + 1 - i; j++)  {
        TETi -= pow(-1., j) / (factorial<Type>(j) * factorial<Type>(o + 1 - i - j)) *
                TriangleA(s + i, std::vector<unsigned> {m + o + 1 - i - j, n + j}, std::vector<Type> {a + c, b - c}, d);
        //std::cout << TETi << " ";
      }
      TET += TETi * pow(-1. / c, i);
      //std::cout << TET << std::endl;
    }

    std::cout << -pow(-1. / c, o) * TriangleA(s + o + 1, std::vector<unsigned> {m, n}, std::vector<Type> {a, b}, d) / c << " ";
    std::cout << TET << " ";

    TET -= pow(-1. / c, o) * TriangleA(s + o + 1, std::vector<unsigned> {m, n}, std::vector<Type> {a, b}, d) / c;
    //std::cout << TET << " ";
    TET *= factorial<Type>(o);
  }
  return TET;
}

#include <bits/stdc++.h>

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
TetrahedronA(const int &s, std::vector<unsigned> m, std::vector <Float1> a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  if(fabs(a[2]) < fabs(a[1])) {
    std::swap(a[2], a[1]);
    std::swap(m[2], m[1]);
  }

  return TetrahedronB(s, m, a, d);

  switch(s) {
    case -1:
      if(a[0] + a[1] + d <= 0) return TetrahedronB(-1, m, a, d);
      else return TetrahedronB(-1, m, std::vector<Type> {-a[0], -a[1], -a[2]}, -d);
      break;
    case 0:
      if(a[0] + a[1] + d <= 0) {
        return TetrahedronB(0, m, a, d);
      }
      else return factorial<Type>(m[1]) * factorial<Type>(m[2]) /
                    (factorial<Type>(2. + m[1] + m[2]) * (3. + m[0] + m[1] + m[2]))
                    - TetrahedronB(0, m, std::vector<Type> {-a[0], -a[1], -a[2]}, -d);
      break;
    default:
      //return TetrahedronB(s, m, a, d);


      Type TET =  pow(-a[2], s) / factorial<Type>(m[2] + s) *
                  TetrahedronA(0, std::vector<unsigned> {m[0], m[1], m[2] + s}, a, d) ;
      for(unsigned i = 1; i <= s; i++)  {
        Type TETi = 0.;
        for(unsigned j = 0; j <= m[2] + i; j++)  {
          TETi += pow(-1., j) / (factorial<Type>(j) * factorial<Type>(m[2] + i - j)) *
                  TriangleA(s + 1 - i, std::vector<unsigned> {m[0] + m[2] + i - j, m[1] + j}, std::vector<Type> {a[0] + a[2], a[1] - a[2]}, d);
        }
        TET += TETi * pow(-a[2], i - 1);
      }
      TET *= factorial <Type> (m[2]);

      return TET;

  }
}


template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Tetrahedron(const int &s, const std::vector<unsigned> &m, const std::vector <Float1> &a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;
  if(fabs(a[0]) <= fabs(a[1]) && fabs(a[0]) <= fabs(a[2])) {
    if(fabs(a[1]) <= fabs(a[2]))
      return TetrahedronA(s, std::vector<unsigned> {m[0], m[1], m[2]}, std::vector<Type> {-a[0], a[1], a[2]}, d + a[0]) ;
    else
      return TetrahedronA(s, std::vector<unsigned> {m[0], m[2], m[1]}, std::vector<Type> {-a[0], a[2], a[1]}, d + a[0]) ;
  }
  else if(fabs(a[1]) <= fabs(a[2])) {
    if(fabs(a[0]) <= fabs(a[2]))
      return TetrahedronA(s, std::vector<unsigned> {m[1], m[0], m[2]}, std::vector<Type> {-a[1], a[0], a[2]}, d + a[1]) ;
    else
      return TetrahedronA(s, std::vector<unsigned> {m[1], m[2], m[0]}, std::vector<Type> {-a[1], a[2], a[0]}, d + a[1]) ;
  }
  else {
    if(fabs(a[0]) <= fabs(a[1]))
      return TetrahedronA(s, std::vector<unsigned> {m[2], m[0], m[1]}, std::vector<Type> {-a[2], a[0], a[1]}, d + a[2]) ;
    else
      return TetrahedronA(s, std::vector<unsigned> {m[2], m[1], m[0]}, std::vector<Type> {-a[2], a[1], a[0]}, d + a[2]) ;
  }

}




template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
PrismB(const int &s, const std::vector<unsigned> &m_input, const std::vector <Float1> &a_input, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const Type &c = a_input[2];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];
  const unsigned &o = m_input[2];

  Type PRI = 0.;
  if(c == 0) { //plane is parallel to z-axis
    if(o % 2 == 0)
      PRI = Type(2) / (o + 1) * Triangle(s, std::vector<unsigned> {m, n}, std::vector<Type> {a, b}, d);
    else
      PRI = 0.;
  }
  else {//all other cases
    for(unsigned k = 1; k <= o + 1; k++)  {
      PRI += (pow(-1., k - 1) * Triangle(s + k, std::vector<unsigned> {m, n}, std::vector<Type> {a, b }, c + d) +
              pow(-1., o - 1) * Triangle(s + k, std::vector<unsigned> {m, n}, std::vector<Type> {a, b }, -c + d))
             / (factorial<Type>(o + 1 - k) * pow(c, k));
    }
    PRI *= factorial<Type>(o);
  }
  return PRI;
}



template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Prism(const int &s, const std::vector<unsigned> &m, const std::vector <Float1> &a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  if(a[0] + a[1] + d <= 0) {
    return PrismB(s, m, a, d);
  }
  switch(s) {
    case -1:
      return PrismB(s, m, std::vector<Type> {-a[0], -a[1], -a[2]}, -d);
      break;
    case 0:
      return (Type(1. + pow(-1, m[2])) / Type(m[2] + 1)) / Type((m[0] + m[1] + 2) * (m[1] + 1))
             - PrismB(s, m, std::vector<Type> {-a[0], -a[1], -a[2]}, -d);
      break;
    default:
      return PrismB(s, m, a, d);
  }
}


template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
HyperCubeBOld(const int &s, unsigned i, const std::vector<unsigned> &m, const std::vector <Float1> &a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  Type HPI = 0.;
  if(i > 0) {
    Type aI = 1. / a[i];

    int sl = (m[i] % 2 == 1) ? 1. : -1.; // this is (-1)^(m-1)
    int sr = 1; // this is (-1)^(j-1) for j = 1, ..., m + 1
    Type c = aI; // this is m!/(m+1-j)! 1/a^j for j = 1,...,m + 1

    for(int j = 1; j <= m[i] + 1;  c *= aI * (m[i] + 1 - j), sr *= -1, j++) {
      HPI += c * (sl * HyperCubeBOld(s + j, i - 1, m, a, -a[i] + d) + sr * HyperCubeBOld(s + j, i - 1, m, a, a[i] + d));
    }
  }
  else {
    HPI = LSI(s, m[0], a[0], d);
  }
  return HPI;
}



template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
HyperCubeA(const unsigned &nn, const int &s, std::vector<unsigned> &m,
           const std::vector <Float1> &a, const std::vector <Float1> &ma,
           const Float2 &d, const Float2 &md);

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
HyperCubeB(const unsigned &n, const int &s, std::vector<unsigned> &m,
           const std::vector <Float1> &a, const std::vector <Float1> &ma,
           const Float2 &d, const Float2 &md) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  Type HCI = 0.;
  Type aI = 1. / a[n];

  int sl = (m[n] % 2 == 1) ? 1. : -1.; // this is (-1)^(m-1)
  int sr = 1; // this is (-1)^(j-1) for j = 1, ..., m + 1
  Type c = aI; // this is m!/(m+1-j)! 1/a^j for j = 1,...,m + 1

  for(int j = 1; j <= m[n] + 1;  c *= aI * (m[n] + 1 - j), sr = -sr, j++) {
    HCI += c * (sl * HyperCubeA(n - 1, s + j, m, a, ma, -a[n] + d, a[n] - d) + sr * HyperCubeA(n - 1, s + j, m, a, ma, a[n] + d, -a[n] - d));
  }
  return HCI;
}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
HyperCubeC(const unsigned &n, const int &s, std::vector<unsigned>& m,
           const std::vector <Float1> &a, const std::vector <Float1> &ma,
           const Float2 &d, const Float2 &md) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  Type HCI = 0.;
  Type an = a[n];

  unsigned mn = m[n];
  Type mnf = factorial<Type>(mn);
  m[n] += s;
  HCI += mnf / factorial<Type>(mn + s) * pow(-an, s) * HyperCubeA(n, 0, m, a, ma, d, md) ;
  m[n] -= s;

  Type s1 = 1;
  Type s2 = (mn % 2 == 0) ? 1. : -1.; //this is pow(-1, mn);
  Type c1 = 1 / Type(mn + 1);

  for(unsigned i = 0; i < s; s1 = -s1, c1 *= an / (mn + 2 + i), i++) {
    HCI += c1 * (s1 * HyperCubeA(n - 1, s - i, m, ma, a, an + d, -an - d) +
                 s2 * HyperCubeA(n - 1, s - i, m, a, ma, -an + d, an - d));
  }
  return  HCI;

}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
HyperCubeA(const unsigned &n, const int &s, std::vector<unsigned> &m,
           const std::vector <Float1> &a, const std::vector <Float1> &ma,
           const Float2 &d, const Float2 &md) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  if(n == 0)  return LSI(s, m[0], a[0], d);

  switch(s) {
    case -1: // interface integral
      if(d <= 0) {
        return HyperCubeB(n, s, m, a, ma, d, md);
      }
      else {
        return HyperCubeB(n, s, m, ma, a, md, d);
      }
      break;
    case 0: // step function integral
      if(d <= 0) {
        return HyperCubeB(n, s, m, a, ma, d, md);
      }
      else {
        Type fullIntegral = 1;
        for(unsigned i = 0; i <= n; i++) {
          fullIntegral *= (m[i] % 2 == 0) ? 2 / Type(m[i] + 1.) : 0;
        }
        return (fullIntegral - HyperCubeB(n, s, m, ma, a, md, d));
      }
      break;
    default: // all other cases, whatever they mean
      if(d < fabs(a[n])) {
        return HyperCubeB(n, s, m, a, ma, d, md);
      }
      else {
        return HyperCubeC(n, s, m, a, ma, d, md);
      }
  }
}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
HyperCube(const int &s, std::vector<unsigned> m, std::vector <Float1> a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  Type HCI = 1.;

  for(int i = a.size() - 1; i >= 0; i--) {
    if(a[i] == 0) {
      if(m[i] % 2 == 1) { // odd monomial on symmetric domain
        return 0.;
      }
      else { // integral of the even monomial from -1 to 1
        HCI *= 2 / Type(m[i] + 1);
        a.erase(a.begin() + i);
        m.erase(m.begin() + i);
      }
    }
  }


  // all the left a[i] coefficients \ne 0
  for(unsigned i = 0; i < a.size() - 1; i++) { // reorder iterated integral from the smallest to the largest coefficient a[i]
    for(unsigned j = i + 1; j < a.size(); j++) {
      if(fabs(a[i]) > fabs(a[j])) {
        std::swap(a[i], a[j]);
        std::swap(m[i], m[j]);
      }
    }

  }
  std::vector <Float1> ma(a.size());
  for(unsigned i = 0; i < a.size(); i++)  ma[i] = -a[i];


  return HCI * HyperCubeA(a.size() - 1, s, m, a, ma, d, -d);
}



#include "TestLine.hpp"
#include "TestQuad.hpp"
#include "TestTriangle.hpp"
#include "TestHexahedron.hpp"
#include "TestTetrahedron.hpp"

int main(int, char**) {

  typedef double myType;

  bool line = false;//true;//true;
  bool quad = false;//true;//false;//true;
  bool triangle = false;//false;//true;
  bool hexahedron = false;//true;
  bool tetrahedron = false;//true;//false;//true;//true;

  myType eps = 5.0e-12;
  if(line) {
    clock_t time = clock();
    TestLine(eps);
    std::cout << "Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
  }
  eps = 5.0e-11;
  if(triangle) TestTriangle(eps);

  eps = 1.0e-12;
  if(quad) {
    clock_t time = clock();
    TestQuad(eps);
    std::cout << "Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
  }
  eps = 1.0e-12;
  if(hexahedron) {

    TestHexahedron(eps);
    clock_t time = clock();
    TestHexahedronTime(eps);
    std::cout << "Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;

    std::vector<double>a = {.25, .15, .335, .455, 0.75};
    std::vector<unsigned>m = {0, 0, 0, 0, 0};
    std::cout << HyperCube(0, m, a, 0) << std::endl;

    a = {1., 0., 0., 0.};
    m = {0, 0, 0, 0};
    std::cout << HyperCube(-1, m, a, 0) << std::endl;

  }
  eps = 1.0e-10;
  if(tetrahedron) TestTetrahedron(eps);




  //typedef double Type;
  typedef boost::multiprecision::cpp_bin_float_oct Type;
  //typedef boost::multiprecision::cpp_bin_float_quad Type;
  //typedef boost::multiprecision::cpp_bin_float_double_extended Type;
  //typedef boost::multiprecision::cpp_bin_float_double Type;

  std::cout.precision(40);
  Type a = Type(0.), b = Type(0.), c = Type(-1.), d = Type(0.001);
  int s = 0;
  unsigned m = 7, n = 5, o = 6;

  Type G0 = TetrahedronA(s, std::vector<unsigned> {m, n, o}, std::vector<Type> {a, b, c}, d);

  

  Type S0 = TriangleA(s + o + 1, std::vector<unsigned> {m, n}, std::vector<Type> {a, b}, d) / pow(-c, o + 1) ;

  Type fn = factorial<Type>(n);
  Type cMb = c - b;
  Type cMbOc = (c - b) / c;

  Type S1 = 0;
  Type fi = Type(1 / fn);
  for(unsigned i = 1; i <= o + 1; i++) {

    fi = -fi / (n + i);
    Type fj =  pow(cMbOc, i) * fn;
    Type S1i = fj;
    for(unsigned j = 1; j < i; j++) {
      fj = fj / cMbOc * (n + j) / j;
      S1i += fj;
    }
    S1 -= fi * TriangleA(s, std::vector<unsigned> {m + o + 1 - i, n + i}, std::vector<Type> {a + c, b - c}, d) /
          factorial<Type>(o + 1 - i) * S1i;

  }

  Type S2 = 0;
  for(unsigned i = 1; i <= o + 1; i++) {
    Type S2i = 0.;
    for(unsigned k = 1; k <= i; k++) {
      Type S2ik = 0.;
      Type f =  fn / factorial<Type>(n + k);
      for(unsigned j = 0; j <= o + 1 - i; j++) {
        S2ik += f / factorial<Type>(o + 1 - i - j);
        f = - f * Type(n + j + 1) / Type((j + 1) * (n + k + j + 1)) ;
      }
      S2i += pow(cMb, k - 1) * mInt0to1LimLi(s + i + 1 - k, m + n + o + 1 - i + k, a + b, d) * S2ik;
    }
    S2 -= S2i / pow(-c, i);
  }
  std::cout << S0 << " " << S1 + S2 << " " << S0 + S2 << std::endl;
  std::cout << (S0 + S1 + S2) * factorial<Type>(o) << std::endl;
  std::cout << G0 << std::endl;
}




















