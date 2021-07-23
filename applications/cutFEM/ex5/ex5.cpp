
//#include "FemusInit.hpp"
//using namespace femus;

#include <boost/math/special_functions/factorials.hpp>
#include <algorithm>
#include <iostream>
#include <stdlib.h>

using boost::math::factorial;

template<typename Type>
Type LimLi(const int &n, const Type &x) {
  if(x < 0.) return 0.;
  else if(n != 0) return -pow(x, n) / factorial<Type>(n);
  else if(x > 0.) return -1.;
  else return -0.5;
}

template<typename Type>
Type HyperCubeB(const int &s, unsigned i, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type &d) {

  Type aiI = 1. / a[i];

  Type dl = -a[i] + d;
  Type dr = a[i] + d;

  Type HPI = 0.;

  int sl = (m[i] % 2 == 1) ? 1. : -1.; // this is (-1)^(m-1)
  int sr = 1; // this is (-1)^(j-1) for j = 1, ..., m + 1
  Type c = aiI; // this is m!/(m+1-j)! 1/a^j for j = 1,...,m + 1
  if(i > 0) {
    for(int j = 1; j <= m[i] + 1;  c *= aiI * (m[i] + 1 - j), sr *= -1, j++) {
      HPI += c * (sl * HyperCubeB(s + j, i - 1, m, a, dl) + sr * HyperCubeB(s + j, i - 1, m, a, dr));
    }
  }
  else {
    for(int j = 1; j <= m[i] + 1;  c *= aiI * (m[i] + 1 - j), sr *= -1, j++) {
      HPI += -c * (sl * LimLi(s + j, dl) + sr * LimLi(s + j, dr));
    }
  }
  return HPI;
}

template<typename Type>
Type HyperCubeA(const int &s, std::vector<unsigned> m, std::vector <Type> a, const Type &d) {
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
  return HCI * HyperCubeB(s, a.size() - 1, m, a, d);
}

template<typename Type>
Type Int0to1LimLi(const int &s, const unsigned &m, const Type &a, const Type &d) {
  Type TRI = 0.;
  for(unsigned i = 1; i <= m + 1; i++) {
    TRI += pow(-1. / a, i) * LimLi(s + i, a + d) / factorial<Type>(m + 1 - i);
  }
  TRI += pow(-1. / a, m) * LimLi(s + m + 1, d) / a;
  TRI *= factorial<Type>(m);
  return TRI;
}


template<typename Type>
Type TriangleReduced(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type &d) {

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];

  Type TRI = 0.;
  //std::cout << "REDC ";
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
  else if(a + b == 0) { // parallel to y=x edge
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


template<typename Type>
Type TriangleFull(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type &d) {

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];

  Type TRI = 0.;
  //std::cout << "FULL ";
  if(b == 0) { //parallel to right edge
    if(a == 0) TRI = -LimLi(s, d) / ((m + n + 2) * (n + 1));
    else TRI = Int0to1LimLi(s, m + n + 1, a, d) / (n + 1.);
  }
  else if(a == 0) { //parallel to bottom edge
    TRI = (Int0to1LimLi(s, n, b, d) - Int0to1LimLi(s, m + n + 1, b, d)) / (m + 1.);
  }
  else if(a + b == 0) { //parallel to y = x edge
    for(unsigned i = 1; i <= n + 1; i++) {
      TRI += pow(1. / a, i) * LimLi(s + i, d) / (factorial<Type>(n + 1 - i) * (m + n + 2 - i));
    }
    TRI += Int0to1LimLi(s + n + 1, m, a, d) / pow(a, n + 1);
    TRI *= factorial<Type>(n);
  }
  else { //generic case
    for(unsigned j = 1; j <= n + 1; j++) {
      TRI -= Int0to1LimLi(s + j, m + n + 1 - j, a + b, d) * pow(-1. / b, j)
             / factorial<Type>(n + 1 - j);
    }
    TRI += Int0to1LimLi(s + n + 1, m, a, d) * pow(-1. / b, n + 1);
    TRI *= factorial<Type>(n);
  }
  return TRI;
}

template<typename Type>
Type Triangle(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type &d) {

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


template<typename Type>
Type TetrahedronB(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type &d) {

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
//   else if(a == 0 && b == 0) {//parallel to the base
//     for(unsigned i = 1; i <= o + 1; i++)  {
//       Type TETi = 0.;
//       for(unsigned j = 0; j <= o + 1 - i; j++)  {
//         TETi -= pow(-1., j) / (factorial<Type>(j) * factorial<Type>(o + 1 - i - j)) *
//                 Triangle(s + i, std::vector<unsigned> {m + o + 1 - i - j, n + j}, std::vector<Type> {c, -c}, d);
//       }
//       TET += TETi * pow(-1. / c, i);
//     }
//     TET += pow(-1. / c, o) * LimLi(s + o + 1, d) / (c * (n + 1) * (n + m + 2));
//     TET *= factorial<Type>(o);
//   }
//   else if(a + c == 0 && b - c == 0) {
//     for(unsigned i = 1; i <= o + 1; i++)  {
//       Type TETi = 0.;
//       for(unsigned j = 0; j <= o + 1 - i; j++)  {
//         TETi += pow(-1., j) / (factorial<Type>(j) * factorial<Type>(o + 1 - i - j)) *
//                 LimLi(s + i, d) / (n + j + 1);
//       }
//       TET += TETi * pow(-1. / c, i) / (m + n + o + 3 - i);
//     }
//     std::cout << "AAA " << TET << " ";
//
//     TET -= pow(-1. / c, o) * Triangle(s + o + 1, std::vector<unsigned> {m, n}, std::vector<Type> {a, b}, d) / c;
//     TET *= factorial<Type>(o);
//   }
  else { //all other cases
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

template<typename Type>
Type Tetrahedron(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type &d) {

  if(a[0] + a[1] + a[2] + d <= 0) {
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


template<typename Type>
Type PrismB(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type &d) {

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

template<typename Type>
Type Prism(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type &d) {

  if(a[0] + a[1] + a[2] + d <= 0) {
    return PrismB(s, m, a, d);
  }
  switch(s) {
    case -1:
      return PrismB(s, m, std::vector<Type> {-a[0], -a[1], -a[2]}, -d);
      break;
    case 0:
      if(m[2] % 2 == 0)
        return 2. / ((m[0] + m[1] + 2) * (m[1] + 1) * (m[2] + 1))
               - PrismB(s, m, std::vector<Type> {-a[0], -a[1], -a[2]}, -d);
      else return 0.;
      break;
    default:
      return PrismB(s, m, a, d);
  }
}


int main(int argc, char** args) {



  long double a1 = cos(M_PI / 4.) * sin(M_PI / 3.);
  long double a2 = sin(M_PI / 4.) * sin(M_PI / 3.);
  long double a3 = sqrt(1 - a1 * a1 - a2 * a2);
  long double d0 = 0.;

  a1 = 0.00001;
  a2 = 0.00001;
  a3 = sqrt(1. - a1 * a1 - a2 * a2);
  d0 = 1. * a3;

  std::cout << "Area = " << 2. * sqrt(1. + pow(a1 / a3, 2.) + pow(a2 / a3, 2.)) << std::endl;


  std::vector <long double> a = {a1, a2, a3};
  long double d = d0;
  std::vector <unsigned> m = {0, 0, 0};

  unsigned degree = 2;
  unsigned N = degree + 1;
  int s = -1;

  std::cout << "Integral Value2 = " << HyperCubeA(0, m, a, d) << std::endl;

  unsigned cnt = N * (N + 1) * (N + 2) / 6;
  std::vector < std::vector < unsigned > > idx(cnt);
  std::vector < double > f(cnt);
  for(unsigned  i = 0; i < idx.size(); i++) idx[i].resize(3);

  cnt = 0u;
  for(unsigned l = 0 ; l < N; l++) {// sum of the monomial powers
    for(int m = l; m >= 0; m--) { // x^m
      for(int n = l - m; n >= 0; n--) { // y^n
        idx[cnt][0] = m;
        idx[cnt][1] = n;
        idx[cnt][2] = l - m - n; // z^{o}
        f[cnt] = HyperCubeA(s, idx[cnt], a, d);
        cnt++;
      }
    }
  }

  std::cout << " 3D element hex, degree " << degree << ", size " << f.size() << std::endl;
  for(unsigned  i = 0; i < idx.size(); i++) {
    std::cout << "F[" << idx[i][0] << "][" << idx[i][1] << "][" << idx[i][2] << "] = " << f[i] << std::endl;
  }

  m = {0, 0};
  a = {-1., 0};
  d =  1;
  std::cout << "Expected value 0.5, value obtained = " << Triangle(-1,  m, a, d) << std::endl;
  std::cout << "Expected value 0.5, value obtained = " << Triangle(0,  m, a, d) << std::endl;

  a = {0., -1};
  d = 0.;
  std::cout << "Expected value 0.5, value obtained = " << Triangle(-1,  m, a, d) << std::endl;
  std::cout << "Expected value 0, value obtained = " << Triangle(0,  m, a, d) << std::endl;

  a = {sqrt(2.) / 2., -sqrt(2.) / 2.};
  d = 0.;
  std::cout << "Expected value " << sqrt(2.) / 2. << ", value obtained = " << Triangle(-1,  m, a, d) << std::endl;
  std::cout << "Expected value 0.5, value obtained = " << Triangle(0, m, a, d) << std::endl;

  a = {sqrt(2.) / 2., -sqrt(2.) / 2.};
  d = -0.5 * a[0];
  std::cout << "Expected value 0.125, value obtained = " << Triangle(0, m, a, d) << std::endl;

  a = {-sqrt(2.) / 2., -sqrt(2.) / 2.};
  d = -a[0];
  std::cout << "Expected value 0.25, value obtained = " << Triangle(0, m, a, d) << std::endl;

  a = {sqrt(2.) / 2., sqrt(2.) / 2.};
  d = -a[0];
  std::cout << "Expected value 0.25, value obtained = " << Triangle(0, m, a, d) << std::endl;

  m = {3, 5};
  a = {0., 1.};
  d =  0.5;
  std::cout << "Expected value 0.5, value obtained = " << Triangle(-1,  m, a, d) << std::endl;
  std::cout << "Expected value 0.5, value obtained = " << Triangle(0,  m, a, d) << std::endl;

  std::cout << "Expected value 0.5, value obtained = " << Triangle(-1,  m, a, d) << std::endl;
  std::cout << "Expected value 0., value obtained = " << Triangle(0,  m, a, d) << std::endl;

  m = {0, 0};
  a = {0.035, .5};
  d =  -0.125;
  std::cout << "Expected value 0.5, value obtained = " << Triangle(-1,  m, a, d) << std::endl;
  std::cout << "Expected value 0.5, value obtained = " << Triangle(0,  m, a, d) << std::endl;

  std::cout << "Expected value 0.5, value obtained = " << Triangle(-1,  m, a, d) << std::endl;
  std::cout << "Expected value 0., value obtained = " << Triangle(0,  m, a, d) << std::endl;


  std::cout << "testing the Thetrahedron \n";
  m = {4, 4, 4};
  a = {0., 0., 0.};
  for(unsigned i = 0; i < 100; i++) {
    a[0] = static_cast <long double>(rand()) / RAND_MAX;
    a[1] = static_cast <long double>(rand()) / RAND_MAX;
    a[2] = -sqrt(1 - a[0] * a[0] + a[1] * a[1]);
    d = static_cast <long double>(rand()) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) std::cout << "Wrong choice of coefficients for the test\n";

    if(fabs(TetrahedronB(-1,  m, a, d) - Tetrahedron(-1,  m, a, d)) > 1.0e-12)
      std::cout << "test surface failed" << std::endl;

    if(fabs(TetrahedronB(0,  m, a, d) - Tetrahedron(0,  m, a, d)) > 1.0e-12)
      std::cout << "test volume failed" << std::endl;
  }

  for(unsigned i = 0; i < 100; i++) {
    a[0] = static_cast <long double>(rand()) / RAND_MAX;
    a[1] = -a[0];
    a[2] = sqrt(1 - a[0] * a[0] + a[1] * a[1]);
    d = static_cast <long double>(rand()) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) std::cout << "Wrong choice of coefficients for the test\n";

    if(fabs(TetrahedronB(-1,  m, a, d) - Tetrahedron(-1,  m, a, d)) > 1.0e-12)
      std::cout << "test surface failed" << std::endl;

    if(fabs(TetrahedronB(0,  m, a, d) - Tetrahedron(0,  m, a, d)) > 1.0e-12)
      std::cout << "test volume failed" << std::endl;
  }

  for(unsigned i = 0; i < 100; i++) {
    a[0] = 0;
    a[1] = 0;
    a[2] = 1;
    d = static_cast <long double>(rand()) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) std::cout << "Wrong choice of coefficients for the test\n";

    if(fabs(TetrahedronB(-1,  m, a, d) - Tetrahedron(-1,  m, a, d)) > 1.0e-12)
      std::cout << "test surface failed" << std::endl;

    if(fabs(TetrahedronB(0,  m, a, d) - Tetrahedron(0,  m, a, d)) > 1.0e-12)
      std::cout << "test volume failed" << std::endl;
  }

  for(unsigned i = 0; i < 100; i++) {
    a[0] = 0.;
    a[1] = 1.;
    a[2] = 0.;
    d = static_cast <long double>(rand()) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) std::cout << "Wrong choice of coefficients for the test\n";

    if(fabs(TetrahedronB(-1,  m, a, d) - Tetrahedron(-1,  m, a, d)) > 1.0e-12)
      std::cout << "test surface failed" << std::endl;

    if(fabs(TetrahedronB(0,  m, a, d) - Tetrahedron(0,  m, a, d)) > 1.0e-12)
      std::cout << "test volume failed" << std::endl;
  }

  for(unsigned i = 0; i < 100; i++) {
    a[0] = 1.;
    a[1] = 0.;
    a[2] = 0.;
    d = static_cast <long double>(rand()) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) std::cout << "Wrong choice of coefficients for the test\n";

    if(fabs(TetrahedronB(-1,  m, a, d) - Tetrahedron(-1,  m, a, d)) > 1.0e-12)
      std::cout << "test surface failed" << std::endl;

    if(fabs(TetrahedronB(0,  m, a, d) - Tetrahedron(0,  m, a, d)) > 1.0e-12)
      std::cout << "test volume failed" << std::endl;
  }

  for(unsigned i = 0; i < 100; i++) {
    a[0] = 1.;
    a[1] = 0.;
    a[2] = 0.;
    d = static_cast <long double>(rand()) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) std::cout << "Wrong choice of coefficients for the test\n";

    if(fabs(TetrahedronB(-1,  m, a, d) - Tetrahedron(-1,  m, a, d)) > 1.0e-12)
      std::cout << "test surface failed" << std::endl;

    if(fabs(TetrahedronB(0,  m, a, d) - Tetrahedron(0,  m, a, d)) > 1.0e-12)
      std::cout << "test volume failed" << std::endl;
  }


  for(unsigned i = 0; i < 100; i++) {
    a[0] = static_cast <long double>(rand()) / RAND_MAX;
    a[1] = sqrt(1. - a[0] * a[0]);
    a[2] = 0.;
    d = static_cast <long double>(rand()) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) std::cout << "Wrong choice of coefficients for the test\n";

    if(fabs(TetrahedronB(-1,  m, a, d) - Tetrahedron(-1,  m, a, d)) > 1.0e-12)
      std::cout << "test surface failed" << std::endl;

    if(fabs(TetrahedronB(0,  m, a, d) - Tetrahedron(0,  m, a, d)) > 1.0e-12)
      std::cout << "test volume failed" << std::endl;
  }

  for(unsigned i = 0; i < 100; i++) {
    a[0] = static_cast <long double>(rand()) / RAND_MAX;
    a[1] = 0.;
    a[2] = sqrt(1. - a[0] * a[0]);
    d = static_cast <long double>(rand()) / RAND_MAX;

    if(fabs(a[0]) < 1.0e-3) a[0] = 0.;
    if(fabs(a[1]) < 1.0e-3) a[1] = 0.;
    if(fabs(a[2]) < 1.0e-3) a[2] = 0.;
    double long det = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= det;
    a[1] /= det;
    a[2] /= det;

    if(a[0] + a[1] + d <= 0) std::cout << "Wrong choice of coefficients for the test\n";

    if(fabs(TetrahedronB(-1,  m, a, d) - Tetrahedron(-1,  m, a, d)) > 1.0e-12) {
      std::cout << "test surface failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << " " << TetrahedronB(-1,  m, a, d) << " " << Tetrahedron(-1,  m, a, d) << std::endl;
    }

    if(fabs(TetrahedronB(0,  m, a, d) - Tetrahedron(0,  m, a, d)) > 1.0e-12) {
      std::cout << "test volume failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << " " << TetrahedronB(0,  m, a, d) << " " << Tetrahedron(0,  m, a, d) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1; i++) {
      
    std::cout<<"AAA\n";
    a[0] = 0.00856861;
    a[1] = 0.;
    a[2] = 0.999963; 
    d = 0.996042;

    if(fabs(a[0]) < 1.0e-3) a[0] = 0.;
    if(fabs(a[1]) < 1.0e-3) a[1] = 0.;
    if(fabs(a[2]) < 1.0e-3) a[2] = 0.;
    double long det = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= det;
    a[1] /= det;
    a[2] /= det;

    if(a[0] + a[1] + d <= 0) std::cout << "Wrong choice of coefficients for the test\n";

    if(fabs(TetrahedronB(-1,  m, a, d) - Tetrahedron(-1,  m, a, d)) > 1.0e-12) {
      std::cout << "test surface failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << " " << TetrahedronB(-1,  m, a, d) << " " << Tetrahedron(-1,  m, a, d) << std::endl;
    }

    if(fabs(TetrahedronB(0,  m, a, d) - Tetrahedron(0,  m, a, d)) > 1.0e-12) {
      std::cout << "test volume failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << " " << TetrahedronB(0,  m, a, d) << " " << Tetrahedron(0,  m, a, d) << std::endl;
    }
  }
  

//   a = {0., 0., 1.};
//   d =  0.;
//   std::cout << "Expected value 0.5, value obtained = " << TetrahedronB(-1,  m, a, d) << std::endl;
//   std::cout << "Expected value 0.5, value obtained = " << TetrahedronB(0,  m, a, d) << std::endl;
//
//   std::cout << "Expected value 0.5, value obtained = " << Tetrahedron(-1,  m, a, d) << std::endl;
//   std::cout << "Expected value 0.5, value obtained = " << Tetrahedron(0,  m, a, d) << std::endl;
//
//   m = {2, 0, 2};
//   a = {0., 1., 0.};
//   d =  0.;
//   std::cout << "Expected value 0.5, value obtained = " << TetrahedronB(-1,  m, a, d) << std::endl;
//   std::cout << "Expected value 0.5, value obtained = " << TetrahedronB(0,  m, a, d) << std::endl;
//
//   std::cout << "Expected value 0.5, value obtained = " << Tetrahedron(-1,  m, a, d) << std::endl;
//   std::cout << "Expected value 0.5, value obtained = " << Tetrahedron(0,  m, a, d) << std::endl;

//   m = {0, 0, 0};
//   a = {1. / sqrt(3.), -1. / sqrt(3.), -1. / sqrt(3.)};
//   d =  0.;
//   std::cout << "Expected value 0.5, value obtained = " << TetrahedronB(-1,  m, a, d) << std::endl;
//   std::cout << "Expected value 0.5, value obtained = " << TetrahedronB(0,  m, a, d) << std::endl;
//
//   std::cout << "Expected value 0.5, value obtained = " << Tetrahedron(-1,  m, a, d) << std::endl;
//   std::cout << "Expected value 0.5, value obtained = " << Tetrahedron(0,  m, a, d) << std::endl;
//
//   m = {0, 0, 0};
//   a = {-1. / sqrt(3.), 1. / sqrt(3.), 1. / sqrt(3.)};
//   d =  0.;
//   std::cout << "Expected value 0.5, value obtained = " << TetrahedronB(-1,  m, a, d) << std::endl;
//   std::cout << "Expected value 0.5, value obtained = " << TetrahedronB(0,  m, a, d) << std::endl;
//
//   std::cout << "Expected value 0.5, value obtained = " << Tetrahedron(-1,  m, a, d) << std::endl;
//   std::cout << "Expected value 0.5, value obtained = " << Tetrahedron(0,  m, a, d) << std::endl;


  return 1;



}



















