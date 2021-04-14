
#include <iostream>
#include <iomanip>

#include <boost/math/special_functions/factorials.hpp>

using boost::math::factorial;

template <class Float1>
typename boost::math::tools::promote_args<Float1>::type
LimLi(const int &n, const Float1 &x) {

  typedef typename boost::math::tools::promote_args<Float1>::type Type;

  if(x < 0.) return 0.;
  else if(n != 0) return -pow(x, n) / factorial<Type>(n);
  else if(x > 0.) return -1.;
  else return -0.5;
}

#include "oldFunctions.hpp"

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Int0to1LimLi(const int &s, const unsigned &m, Float1 a, Float2 d) {
  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  Type INT = 0.;
  bool unitStepFunctionComplement = false;

  if(d > 0) { //this assures that for s < 1, either limit a + d or d below is negative
    if(s == -1) { //this is the Dirac distribution, the sign of the normal can be changed and INT(delta(ax+d)) = INT(delta(-ax-d))
      a = -a;
      d = -d;
    }
    else if(s == 0) { //this is the unit step function, the sign of the normal can be changed if INT(U(ax+d)) = INT(1) - INT(U(-ax-d))
      a = -a;
      d = -d;
      unitStepFunctionComplement = true;
      INT = -1. / (m + 1);
    }
  }

  Type x(a + d);

  if(x < 0 || d < 0) { // in all these cases no-significant digit cancellation occurs
    if(x >= 0.) {
      Type g =  1. / (-a);
      for(unsigned i = 1; i <= m + 1; i++) {
        INT += g * LimLi(s + i, x);
        g *= (m + 1 - i) / (-a);
      }
    }
    else if(d >= 0.) {
      INT += pow(-1. / a, m) * LimLi(s + m + 1, d) / a * factorial<Type>(m);
    }
  }
  else { //alternative formula to avoid significant digit cancellations when s>1, and (a+d) and d are non-negative and d >> a
    for(unsigned i = 1; i <= s; i++) {
      INT -= pow(-a, s - i) / factorial<Type>(m + 1 + s - i) * LimLi(i, x) ;
    }
    INT += pow(-a, s) / factorial<Type>(m + 1 + s);
    INT *= factorial<Type>(m);
  }
  return (unitStepFunctionComplement) ? -INT : INT;
}



template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Intm1to1LimLi(const int &s, const unsigned &m, Float1 a, Float2 d) {
  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  Type INT = 0.;
  bool unitStepFunctionComplement = false;

  if(d > 0) { //this assures that for s < 1, either limit x or y below is negative
    if(s == -1) { //this is the Dirac distribution, the sign of the normal can be changed and INT(delta(ax+d)) = INT(delta(-ax-d))
      a = -a;
      d = -d;
    }
    else if(s == 0) { //this is the unit step function, the sign of the normal can be changed if INT(U(ax+d)) = INT(1) - INT(U(-ax-d))
      a = -a;
      d = -d;
      unitStepFunctionComplement = true;
      INT = -(1. + pow(-1., m)) / (m + 1);
    }
  }

  Type x(a + d);
  Type y(-a + d);

  if(x < 0 || y < 0) { // in all these cases no-significant digit cancellation occurs
    if(x >= 0.) {
      Type g =  1. / (-a);
      for(unsigned i = 1; i <= m + 1; i++) {
        INT += g * LimLi(s + i, x);
        g *= (m + 1 - i) / (-a);
      }
    }
    else if(y >= 0.) {
      Type g =  pow(-1., m) / a;
      for(unsigned i = 1; i <= m + 1; i++) {
        INT += g * LimLi(s + i, y);
        g *= (m + 1 - i) / a;
      }
    }
  }
  else { //alternative formula to avoid significant digit cancellations when s>1, and (a+d) and (-a+d) are non-negative and d >> a
    for(unsigned i = 1; i <= s; i++) {
      INT += pow(-a, s - i) / factorial<Type>(m + 1 + s - i) *
             (pow(x, i) + pow(-1, m + s - i) * pow(y, i)) / factorial<Type>(i);
    }
    INT += pow(-a, s) * (1 + pow(-1, m + s)) / factorial<Type>(m + 1 + s);
    INT *= factorial<Type>(m);
  }
  return (unitStepFunctionComplement) ? -INT : INT;
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
  //std::cout << "REDC " << a << " " << b << " " << d << " ";
  if(b == 0) { //parallel to right edge
    if(a == 0) TRI = -LimLi(s, d) / ((m + n + 2) * (n + 1));
    else TRI = (-LimLi(s + 1, a + d) +
                  factorial<Type>(m + n + 1) * pow(-1. / a, m + n + 1) * LimLi(s + m + n + 2, d)
                 ) / (a * (n + 1.));
  }
  else if(a == 0) { //parallel to bottom edge
    TRI = (factorial<Type>(n) * pow(-1. / b, n) * LimLi(s + n + 1, d) -
           factorial<Type>(m + n + 1) * pow(-1. / b, m + n + 1) * LimLi(s + m + n + 2, d)
          ) / (b * (m + 1.));
  }
  else if(a + b == 0) { // parallel to y=x edge//TODO
    for(unsigned i = 1; i <= m + 1; i++) {
      TRI += pow(-1. / a, i) * LimLi(s + n + 1 + i, a + d) / factorial<Type>(m + 1 - i);
    }
    TRI *= factorial<Type>(n) * factorial<Type>(m) * pow(1. / a, n + 1);
    TRI += LimLi(s + 1, d) / (a * (m + n + 1));
  }
  else { // general case
    for(unsigned i = 1; i <= n + 1; i++) {
      TRI += factorial<Type>(m + n + 1 - i) / factorial<Type>(n + 1 - i) * pow((a + b) / b, i);
    }
    TRI *= pow(-1., m + n) * LimLi(s + m + n + 2, d) / pow(a + b, m + n + 2) ;

    TRI += Int0to1LimLi(s + n + 1, m, a, d) * pow(-1. / b, n + 1);
    TRI *= factorial<Type>(n);

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
    else TRI = Int0to1LimLi(s, m + n + 1, a, d) / (n + 1.);
  }
  else if(a == 0) { //parallel to bottom edge
    TRI = (Int0to1LimLi(s, n, b, d) - Int0to1LimLi(s, m + n + 1, b, d)) / (m + 1.);
  }
  else if(a + b == 0) { //parallel to y = x edge //TODO
    for(unsigned i = 1; i <= n + 1; i++) {
      TRI += pow(1. / a, i) * LimLi(s + i, d) / (factorial<Type>(n + 1 - i) * (m + n + 2 - i));
    }
    TRI += Int0to1LimLi(s + n + 1, m, a, d) / pow(a, n + 1);
    TRI *= factorial<Type>(n);
  }
  else { //generic case
    if(fabs(a) < fabs(b)) { // |a| < |b|
      for(unsigned j = 1; j <= n + 1; j++) {
        TRI -= Int0to1LimLi(s + j, m + n + 1 - j, a + b, d) * pow(-1. / b, j)
               / factorial<Type>(n + 1 - j);
      }
      TRI += Int0to1LimLi(s + n + 1, m, a, d) * pow(-1. / b, n + 1);
      TRI *= factorial<Type>(n);
    }
    else { // |b| < |a|
      for(unsigned j = 1; j <= m + 1; j++) {
        TRI += (Int0to1LimLi(s + j, m + n + 1 - j, a + b, d) -
                Int0to1LimLi(s + j, n, b, a + d)) * pow(-1. / a, j)
               / factorial<Type>(m + 1 - j);
      }
      TRI *= factorial<Type>(m);
    }
  }

  return TRI;
}



template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Triangle(const int &s, const std::vector<unsigned> &m, const std::vector <Float1> &a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  if(a[0] + a[1] + d <= 0) {
    return TriangleReduced(s, m, a, d);
  }
  switch(s) {
    case -1:
      return TriangleReduced(s, m, std::vector<Type> {-a[0], -a[1]}, -d);
      break;
    case 0:
      return 1. / ((1. + m[1]) * (2. + m[1] + m[0])) - TriangleReduced(s, m, std::vector<Type> {-a[0], -a[1]}, -d);
      break;
    default:
      return TriangleFull(s, m, a, d);
  }
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
             Triangle(s, std::vector<unsigned> {m + o + 1 - i, n + i}, std::vector<Type> {a, b}, d);
    }
    TET *= factorial<Type>(o);
  }
  else { //all other cases
    //std::cout<<"BBB\n";
    for(unsigned i = 1; i <= o + 1; i++)  {
      Type TETi = 0.;
      for(unsigned j = 0; j <= o + 1 - i; j++)  {
        TETi -= pow(-1., j) / (factorial<Type>(j) * factorial<Type>(o + 1 - i - j)) *
                Triangle(s + i, std::vector<unsigned> {m + o + 1 - i - j, n + j}, std::vector<Type> {a + c, b - c}, d);
      }
      TET += TETi * pow(-1. / c, i);
    }
    TET -= pow(-1. / c, o) * Triangle(s + o + 1, std::vector<unsigned> {m, n}, std::vector<Type> {a, b}, d) / c;
    TET *= factorial<Type>(o);
  }
  return TET;
}



template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Tetrahedron(const int &s, const std::vector<unsigned> &m, const std::vector <Float1> &a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  if(a[0] + a[1] + d <= 0) {
    return TetrahedronB(s, m, a, d);
  }
  switch(s) {
    case -1:
      return TetrahedronB(s, m, std::vector<Type> {-a[0], -a[1], -a[2]}, -d);
      break;
    case 0:
      return factorial<Type>(m[1]) * factorial<Type>(m[2]) /
             (factorial<Type>(2. + m[1] + m[2]) * (3. + m[0] + m[1] + m[2]))
             - TetrahedronB(s, m, std::vector<Type> {-a[0], -a[1], -a[2]}, -d);
      break;
    default:
      return TetrahedronB(s, m, a, d);
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
      PRI = 2. / (o + 1) * Triangle(s, std::vector<unsigned> {m, n}, std::vector<Type> {a, b}, d);
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
      return ((1. + pow(-1., m[2])) / (m[2] + 1)) / ((m[0] + m[1] + 2) * (m[1] + 1))
             - PrismB(s, m, std::vector<Type> {-a[0], -a[1], -a[2]}, -d);
      break;
    default:
      return PrismB(s, m, a, d);
  }
}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
HyperCubeB(const int &s, unsigned i, const std::vector<unsigned> &m, const std::vector <Float1> &a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  Type HPI = 0.;
  if(i > 0) {
    Type aI = 1. / a[i];

    int sl = (m[i] % 2 == 1) ? 1. : -1.; // this is (-1)^(m-1)
    int sr = 1; // this is (-1)^(j-1) for j = 1, ..., m + 1
    Type c = aI; // this is m!/(m+1-j)! 1/a^j for j = 1,...,m + 1

    for(int j = 1; j <= m[i] + 1;  c *= aI * (m[i] + 1 - j), sr *= -1, j++) {
      HPI += c * (sl * HyperCubeB(s + j, i - 1, m, a, -a[i] + d) + sr * HyperCubeB(s + j, i - 1, m, a, a[i] + d));
    }
  }
  else {
    HPI = Intm1to1LimLi(s, m[0], a[0], d);
  }
  return HPI;
}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
HyperCube(const int &s, std::vector<unsigned> m, std::vector <Float1> a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  Type HCI = 1.;
  for(int i = a.size() - 1; i >= 0; i--) {
    if(a[i] == 0) {
      if(m[i] % 2 == 1) {
        return 0.;
      }
      else {
        HCI *= 2. / (m[i] + 1.);
        a.erase(a.begin() + i);
        m.erase(m.begin() + i);
      }
    }
  }

  for(unsigned i = 0; i < a.size() - 1; i++) {
    for(unsigned j = 1; j < a.size(); j++) {
      if(fabs(a[i]) > fabs(a[j])) {
        std::swap(a[i], a[j]);
        std::swap(m[i], m[j]);
      }
    }
  }

  if(d <= 0) {
    return HCI * HyperCubeB(s, a.size() - 1, m, a, d);
  }
  Type intOf1 = 1.;
  switch(s) {
    case -1:
      for(unsigned i = 0; i < a.size(); i++) a[i] = -a[i];
      return HCI * HyperCubeB(s, a.size() - 1, m, a, -d);
      break;
    case 0:
      for(unsigned i = 0; i < a.size(); i++) {
        a[i] = -a[i];
        intOf1 *= (1. + pow(-1., m[i])) / (m[i] + 1.);
      }
      return HCI * (intOf1 - HyperCubeB(s, a.size() - 1, m, a, -d));
      break;
    default:
      HCI * HyperCubeB(s, a.size() - 1, m, a, d);
  }
}



#include "TestLine.hpp"
#include "TestQuad.hpp"
#include "TestTriangle.hpp"
#include "TestHexahedron.hpp"
#include "TestTetrahedron.hpp"

int main(int, char**) {

  typedef double myType;

  bool line = true; // false;//true;
  bool quad = true; //false;//true;
  bool triangle = true;
  bool hexahedron =true; //false;//true;
  bool tetrahedron = true; //false;//true;

  myType eps = 5.0e-11;
  if(line) TestLine(eps);

  eps = 5.0e-11;
  if(triangle) TestTriangle(eps);

  eps = 1.0e-12;
  if(quad) TestQuad(eps);

  eps = 1.0e-10;
  if(hexahedron) TestHexahedron(eps);

  eps = 1.0e-10;
  if(tetrahedron) TestTetrahedron(eps);



}









