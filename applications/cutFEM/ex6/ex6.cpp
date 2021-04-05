
#include <iostream>
#include <iomanip>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/factorials.hpp>

using boost::multiprecision::cpp_dec_float_50;
using boost::multiprecision::cpp_dec_float_100;
using boost::math::factorial;

// namespace boost{ namespace multiprecision{
//
// typedef number<cpp_dec_float<200> > cpp_dec_float_200;
//
// }} // namespaces


template <class Float1>
typename boost::math::tools::promote_args<Float1>::type
LimLi(const int &n, const Float1 &x) {

  typedef typename boost::math::tools::promote_args<Float1>::type Type;

  if(x < 0.) return 0.;
  else if(n != 0) return -pow(x, n) / factorial<Type>(n);
  else if(x > 0.) return -1.;
  else return -0.5;
}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Int0to1LimLi(const int &s, const unsigned &m, const Float1 &a, const Float2 &d) {

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;
  Type TRI(0.);
  Type ma(1. / (-a));
  Type f(factorial<Type>(m));
  for(unsigned i = 1; i <= m + 1; i++) {
    TRI += ma * LimLi(s + i, a + d) / f;
    f /= (m + 1 - i);
    ma /= (-a);
  }
  TRI += pow(-1. / a, m) * LimLi(s + m + 1, d) / a;
  TRI *= factorial<Type>(m);
  return TRI;

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





int main(int, char**) {

  std::cout << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);
  std::cout << LimLi(3, 2.5) << std::endl;
  std::cout << std::scientific << std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10);
  std::cout << LimLi(3, cpp_dec_float_50(2.5)) << std::endl;

  std::cout << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);
  std::cout << Int0to1LimLi(0, 10, 0.1001, 0.999) << std::endl;

  std::cout << std::scientific << std::setprecision(std::numeric_limits<long double>::digits10);
  std::cout << Int0to1LimLi(0, 10, (long double)0.1001, (long double) 0.999) << std::endl;

  std::cout << std::scientific << std::setprecision(std::numeric_limits<cpp_dec_float_50>::digits10);
  std::cout << Int0to1LimLi(0, 10, cpp_dec_float_50(0.00001001), cpp_dec_float_50(0.999)) << std::endl;

  std::cout << std::scientific << std::setprecision(std::numeric_limits<cpp_dec_float_100>::digits10);
  std::cout << Int0to1LimLi(0, 10, cpp_dec_float_100(0.0000001001), cpp_dec_float_100(0.999)) << std::endl;
  {
    std::vector<long double> a{0, 0};
    std::vector<unsigned> m = {4, 4};

    for(unsigned i = 0; i < 1; i++) {

      std::cout << "BBB\n";
      a[0] = 0.0856861;
      a[1] = 0.;

      long double d = 0.996042;

      std::cout << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);
      std::cout << "test volume failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n" << TriangleFull(0,  m, a, d) << "\n" << Triangle(0,  m, a, d) << std::endl;


    }
  }
  {
    std::vector<cpp_dec_float_50> a{0, 0};
    std::vector<unsigned> m = {4, 4};

    for(unsigned i = 0; i < 1; i++) {

      std::cout << "BBB\n";
      a[0] = 0.000000856861;
      a[1] = 0.;

      cpp_dec_float_50 d = 0.996042;

      std::cout << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);
      std::cout << "test volume failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n" << TriangleFull(0,  m, a, d) << "\n" << Triangle(0,  m, a, d) << std::endl;


    }
  }

  {
    std::vector<cpp_dec_float_100> a{0, 0};
    std::vector<unsigned> m = {4, 4};

    for(unsigned i = 0; i < 1; i++) {

      std::cout << "BBB\n";
      a[0] = 0.0000000000856861;
      a[1] = 0.;

      //cpp_dec_float_100 d = 0.996042;
      double d = 0.996042;

      std::cout << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);
      std::cout << "test volume failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n" << TriangleFull(0,  m, a, d) << "\n" << Triangle(0,  m, a, d) << std::endl;
    }
  }
  {
    std::vector<long double> a{0, 0, 0};
    std::vector<unsigned> m = {4, 4, 4};

    a[0] = 0.00000856861;
    a[1] = 0.;
    a[2] = 0.999963;
    long double d = 0.996042;

    if(fabs(a[0]) < 1.0e-8) a[0] = 0.;
    if(fabs(a[1]) < 1.0e-8) a[1] = 0.;
    if(fabs(a[2]) < 1.0e-8) a[2] = 0.;
    long double det = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= det;
    a[1] /= det;
    a[2] /= det;

    std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << " " << TetrahedronB(-1,  m, a, d) << " " << Tetrahedron(-1,  m, a, d) << std::endl;
    std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << " " << TetrahedronB(0,  m, a, d) << " " << Tetrahedron(0,  m, a, d) << std::endl;
  }

  {
    std::vector<boost::multiprecision::cpp_dec_float_100> a{0, 0, 0};
    std::vector<unsigned> m = {7, 7, 7};

    a[0] = 0.000000856861;
    a[1] = 0.;
    a[2] = 0.999963;
    boost::multiprecision::cpp_dec_float_100 d = 0.996042;

    if(fabs(a[0]) < 1.0e-18) a[0] = 0.;
    if(fabs(a[1]) < 1.0e-18) a[1] = 0.;
    if(fabs(a[2]) < 1.0e-18) a[2] = 0.;
    boost::multiprecision::cpp_dec_float_100 det = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= det;
    a[1] /= det;
    a[2] /= det;

    std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << " " << TetrahedronB(-1,  m, a, d) << " " << Tetrahedron(-1,  m, a, d) << std::endl;
    std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << " " << TetrahedronB(0,  m, a, d) << " " << Tetrahedron(0,  m, a, d) << std::endl;
  }


}

