
#include <iostream>
#include <iomanip>

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/beta.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>

using namespace boost;
using namespace accumulators;

using boost::math::factorial;
using boost::math::beta;


#include <boost/multiprecision/cpp_bin_float.hpp>
namespace boost {
  namespace multiprecision {
    typedef number < backends::cpp_bin_float < 24, backends::digit_base_2, void, boost::int16_t, -126, 127 >, et_off >         cpp_bin_float_single;
    typedef number < backends::cpp_bin_float < 53, backends::digit_base_2, void, boost::int16_t, -1022, 1023 >, et_off >       cpp_bin_float_double;
    typedef number < backends::cpp_bin_float < 64, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >     cpp_bin_float_double_extended;
    typedef number < backends::cpp_bin_float < 113, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >    cpp_bin_float_quad;
    typedef number < backends::cpp_bin_float < 237, backends::digit_base_2, void, boost::int32_t, -262142, 262143 >, et_off >  cpp_bin_float_oct;
  }
} // namespaces

using boost::multiprecision::cpp_bin_float_oct;
using boost::multiprecision::cpp_bin_float_quad;

template <class TypeA>
TypeA F(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d, const TypeA &x1, const TypeA &x2) {

  std::cout << "F" << std::endl;

  TypeA sum(0);
  for(unsigned i = 0; i <= s; i++) {
    TypeA c1 = pow(a[1], i) / ((m[1] + i + 1) * factorial<TypeA>(i));
    for(unsigned k = 0; k <= s - i; k++) {
      sum += c1 * pow(a[0], k) * pow(d, s - i - k) * (pow(x2, m[0] + k + 1) - pow(x1, m[0] + k + 1)) /
             (factorial<TypeA>(s - i - k) *  factorial<TypeA>(k) * (m[0] + k + 1));
    }
  }
  return  sum;
}

template <class TypeA>
TypeA G(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d, const TypeA &x1, const TypeA &x2) {

  std::cout << "G" << std::endl;

  TypeA sum1 = ((m[1] % 2 == 0) ? -1 : 1) * (factorial<TypeA>(m[1]) / pow(a[1], m[1] + 1));

  std::cout << sum1 << std::endl;

  TypeA sum = 0;
  for(unsigned j = 0; j <= s + m[1] + 1; j++) {
    sum += (pow(x2, m[0] + j + 1) - pow(x1, m[0] + j + 1)) * (pow(a[0], j) * pow(d, s + m[1] + 1 - j)) / (factorial<TypeA>(s + m[1] + 1 - j) * factorial<TypeA>(j) * (m[0] + j + 1));

  }
  std::cout << sum << std::endl;
  
  return sum * sum1;
}

template <class TypeIO, class TypeA>
TypeIO SquareA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeIO> &a, const TypeIO &d) {

  std::vector<TypeA> aA ={a[0], a[1]};
  TypeA dA(d);
  
  TypeA xf = (-aA[1] - dA) / aA[0];
  TypeA xg = - dA / aA[0];
  

  TypeA INT(0);

  if(a[0] > 0) {
    if(xf < 0) {
      if(xg < 0) {
        INT = F<TypeA>(s, m, aA, dA, 0, 1);
      }
      else if(xg < 1) {
        INT =  F<TypeA>(s, m, aA, dA, 0, 1) - G<TypeA>(s, m, aA, dA, 0, xg) ;
      }
      else {
        INT =  F<TypeA>(s, m, aA, dA, 0, 1) - G<TypeA>(s, m, aA, dA, 0, 1)  ;
      }
    }
    else if(xf < 1) {
      if(xg < 0) {
        INT = G<TypeA>(s, m, aA, dA, 0, xf) +  F<TypeA>(s, m, aA, dA, xf, 1);
      }
      else if(xg < xf) {
        INT = G<TypeA>(s, m, aA, dA, xg, xf) +  F<TypeA>(s, m, aA, dA, xf, 1);
      }
      else if(xg < 1) {
        INT = -G<TypeA>(s, m, aA, dA, xf, xg) +  F<TypeA>(s, m, aA, dA, xf, 1);
      }
      else {
        INT =  F<TypeA>(s, m, aA, dA, xf, 1) - G<TypeA>(s, m, aA, dA, xf, 1)  ; //xg=1,xf=0,d=-1,a[0]=a[1]=1
      }
    }
    else {
      if(xg < 0) {
        INT = G<TypeA>(s, m, aA, dA, 0, 1);
      }
      else if(xg < 1) {
        INT = G<TypeA>(s, m, aA, dA, xg, 1);
      }
    }
  }
  else {
    if(xf > 1) {
      if(xg > 1) {
        INT =  F<TypeA>(s, m, aA, dA, 0, 1);
      }
      else if(xg > 0) {
        INT =  F<TypeA>(s, m, aA, dA, 0, 1) - G<TypeA>(s, m, aA, dA, 0, xg);
      }
      else {
        INT =  F<TypeA>(s, m, aA, dA, 0, 1) - G<TypeA>(s, m, aA, dA, 0, 1); //xf=1,xg=0,a[0]=-1,a[1]=1,d=0
      }
    }
    else if(xf > 0) {
      if(xg > 1) {
        INT =  F<TypeA>(s, m, aA, dA, 0, xf) + G<TypeA>(s, m, aA, dA, xf, 1);
      }
      else if(xg > xf) {
        INT =  F<TypeA>(s, m, aA, dA, 0, xf) + G<TypeA>(s, m, aA, dA, xf, xg);
      }
      else if(xg > 0) {
        INT =  F<TypeA>(s, m, aA, dA, 0, xf) - G<TypeA>(s, m, aA, dA, xg, xf);
      }
      else {
        INT =  F<TypeA>(s, m, aA, dA, 0, xf) - G<TypeA>(s, m, aA, dA, 0, xf)  ;
      }
    }
    else {
      if(xg > 1) {
        INT = G<TypeA>(s, m, aA, dA, 0, 1);
      }
      else if(xg > 0) {
        INT =  G<TypeA>(s, m, aA, dA, 0, xg);
      }
    }
  }

  return static_cast<TypeIO>(INT);
}


int main() {

  std::vector<double> a = {1, -1};
  std::vector<unsigned> m = {0, 0};
  double d = 0.;
  unsigned q = 6;
  double dt = 0.;
  double delta = 0.;
  double temp1 = 0.;

//   for(unsigned k = 0; k < 401; k++) {
//
//     a[0] = -1 / (1. - (2 * delta));
//     d = delta / (1 - (2 * delta));
//     delta += 0.005;
//
//     if(abs(delta - 0.5) > 0.00000000001){
//         dt = SquareA<double, double>(0, m, a, d);
//     }
//     temp1 = dt - 0.5;
//     //std::cout << "dt =  " << dt << std::endl;
//
//     if(!(abs(temp1) < 0.000000000000001)){
//
//         std::cout << "failedmmm!! " << dt << std::endl;
//     }
// //     std::cout << "new = " << SquareA<double, double>(0, m, a, d) << std::endl;
//      std::cout << "a = " << a[0] << std::endl;
// //     std::cout << "delta = " << delta << std::endl;
//
//
// //     for(unsigned j = 0; j < 10; j++) {
// //
// //
// //
// //     }
//   }

  for(unsigned k = 0; k < q; k++) {

    for(unsigned j = 0; j < q; j++) {
      double temp  = 1 / (double)((j + 1) * (k + j + 2));
      double diff = temp - SquareA <double, double> (0, {k, j}, a, d);//TODO Fix ordering...?

      //std::cout << "new = " << SquareA<double, double>(0, {k, j}, a, d) << std::endl;
      //std::cout << "analytical = " << temp << std::endl;
      std::cout << "difference for m = " << k << " and n = " << j << "difference = " << diff << std::endl;
    }
  }

  a[0] = -1;
  a[1] = -1;

  m[0] = 3;
  m[1] = 5;
  d = 1.0e-10;

  double I1 = pow(d, 2 + m[0] + m[1]) * factorial<double> (m[0]) * factorial<double>(m[1]) / factorial<double> (2 + m[0] + m[1]);

  std::cout << I1 << " " << SquareA<double, double>(0, m, a, d) << std::endl;


  double eps = 1.0e-10;

  a[0] = -1;
  a[1] = -eps;

  m[0] = 3;
  m[1] = 5;
  d = eps;

  //I1 = pow(eps, -1 - m[1]) * beta<double> (m[1]+1,  m[0]+2, eps) / (1 + m[0]);
  std::cout << SquareA<double, double>(0, m, a, d) << std::endl;
  std::cout << SquareA<double, cpp_bin_float_oct>(0, m, a, d) << std::endl;
  std::cout << SquareA<double, double>(0, m, a, d) << std::endl;


  return 1;
}
