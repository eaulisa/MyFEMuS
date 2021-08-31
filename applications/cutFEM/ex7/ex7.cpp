
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

using boost::multiprecision::cpp_dec_float_50;
using boost::multiprecision::cpp_dec_float_100;

namespace boost {
  namespace multiprecision {
    typedef number<cpp_dec_float<200> > cpp_dec_float_200;

    typedef number < backends::cpp_bin_float < 24, backends::digit_base_2, void, boost::int16_t, -126, 127 >, et_off >         cpp_bin_float_single;
    typedef number < backends::cpp_bin_float < 53, backends::digit_base_2, void, boost::int16_t, -1022, 1023 >, et_off >       cpp_bin_float_double;
    typedef number < backends::cpp_bin_float < 64, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >     cpp_bin_float_double_extended;
    typedef number < backends::cpp_bin_float < 113, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >    cpp_bin_float_quad;
    typedef number < backends::cpp_bin_float < 237, backends::digit_base_2, void, boost::int32_t, -262142, 262143 >, et_off >  cpp_bin_float_oct;

    typedef number<cpp_dec_float<7> > cpp_dec_float_7;
    typedef number<cpp_dec_float<14> > cpp_dec_float_14;
    typedef number<cpp_dec_float<21> > cpp_dec_float_21;
    typedef number<cpp_dec_float<28> > cpp_dec_float_28;
    typedef number<cpp_dec_float<35> > cpp_dec_float_35;
    typedef number<cpp_dec_float<42> > cpp_dec_float_42;
    typedef number<cpp_dec_float<49> > cpp_dec_float_49;
    typedef number<cpp_dec_float<56> > cpp_dec_float_56;
    typedef number<cpp_dec_float<63> > cpp_dec_float_63;
    typedef number<cpp_dec_float<63> > cpp_dec_float_70;

  }
} // namespaces

#include <iostream>
#include <iomanip>

#include <boost/math/special_functions/factorials.hpp>

using namespace boost;

using boost::math::factorial;

template <class Type>
Type LimLi(const int &n, const Type &x) {

  if(x < 0) return Type(0);
  else if(n != 0) return -pow(x, n) / factorial<Type>(n);
  else if(x > 0) return Type(-1);
  else return Type(-0.5);
}


template <class Type>
Type LSIm1(const int &m, const Type &a, Type d) {

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

template <class Type>
Type LSI0(const int &m, const Type &a, Type d) {

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

template <class Type>
Type LSI(const int &s, const unsigned &m, const Type a, Type d) {
    
  switch(s) {
    case -1:
      return LSIm1(m, a, d);
      break;

    case 0:
      return LSI0(m, a, d);
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


#include <bits/stdc++.h>



template <class Type>
Type HyperCubeA(const unsigned & n, const int &s, std::vector<unsigned> &m,
           const std::vector <Type> &a, const std::vector <Type> &ma,
           const Type & d, const Type & md);

template <class Type>
Type HyperCubeC(const unsigned & n, const int &s, std::vector<unsigned>& m,
           const std::vector <Type> &a, const std::vector <Type> &ma,
           const Type & d, const Type & md);

template <class Type>
Type HyperCubeB(const unsigned & n, const int &s, std::vector<unsigned> &m,
           const std::vector <Type> &a, const std::vector <Type> &ma,
           const Type & d, const Type & md) {
    
  Type aI = 1. / a[n];
  Type c = aI; // this is m!/(m+1-j)! 1/a^j for j = 1,...,m + 1

  Type I1 = 0;

  for(int j = 1; j <= m[n] + 1;  c *= -aI * (m[n] + 1 - j), j++) {
    I1 += c * HyperCubeA(n - 1, s + j, m, a, ma, d + a[n], -d - a[n]);
  }

  Type I2 = factorial<Type>(m[n]) * pow(-aI, m[n] + 1) * HyperCubeA(n - 1, s + m[n] + 1, m, a, ma, d, -d);
  //std::cout << I1 << " " << I2 << "\n";
  return I1 + I2;

//   if(std::max(fabs(I1), fabs(I2)) == 0 || fabs(I1 / I2 + 1)  > 1.0e-3) {
//     return I1 + I2;
//   }
//   else  {
//     std::cout << "?";
//     return HyperCubeC(n, s, m, a, ma, d, md);
//   }


}

template <class Type>
Type HyperCubeC(const unsigned & n, const int &s, std::vector<unsigned>& m,
           const std::vector <Type> &a, const std::vector <Type> &ma,
           const Type & d, const Type & md){

  Type HCI = 0.;
  Type an = a[n];

  unsigned mn = m[n];
  Type mnf = factorial<Type>(mn);
//   m[n] += (s + 1);
//   Type I1 = mnf / factorial<Type>(mn + s + 1) * pow(-an, s + 1) * HyperCubeA(n, -1, m, a, ma, d, md) ;
//   m[n] -= (s + 1);

  m[n] += s;
  Type I1 = mnf / factorial<Type>(mn + s) * pow(-an, s) * HyperCubeA(n, 0, m, a, ma, d, md) ;
  m[n] -= s;

  Type I2 = 0;
  Type c1 = 1 / Type(mn + 1);

  for(unsigned i = 0; i < s; c1 *= -an / (mn + 2 + i), i++) {
    I2 += c1 * HyperCubeA(n - 1, s - i, m, a, ma, d + an, -(d + an));
  }
  //std::cout << I1 << " " << I2 << "\n";
  if(std::max(fabs(I1), fabs(I2)) == 0 || fabs(I1 / I2 + 1)  > 1.0e-3) {
    return I1 + I2;
  }
  else  {
    std::cout << "?";
    return HyperCubeB(n, s, m, a, ma, d, md);
  }

}

template <class Type>
Type HyperCubeA(const unsigned & n, const int &s, std::vector<unsigned> &m,
           const std::vector <Type> &a, const std::vector <Type> &ma,
           const Type & d, const Type & md) {

  
  if(n == 0)  return LSI(s, m[0], a[0], d);

//   Type sum = d;
//   for(unsigned i = 0; i < a.size(); i++) {
//     sum += a[i];
//   }
//
//   Type I1, I2;
//
//   if(sum <= 0) {
//     std::cout << "!";
//     return HyperCubeB(n, s, m, a, ma, d, md);
//   }
//   else {
//     switch(s) {
//       case -1:
//         return HyperCubeB(n, -1, m, ma, a, md, d);
//       case 0:
//         I1 = 1;
//         for(unsigned i = 0; i < a.size(); i++) {
//           I1 *= 1. / Type(m[i] + 1);
//         }
//         I2 = -HyperCubeB(n, 0, m, ma, a, md, d);
//         if(std::max(fabs(I1), fabs(I2)) == 0 || fabs(I1 / I2 + 1)  > 1.0e-3) {
//           return I1 + I2;
//         }
//         else  {
//           std::cout << "@";
//           return HyperCubeB(n, 0, m, a, ma, d, md);
//         }
//       default:
//         std::cout << "*";
//         return HyperCubeC(n, s, m, a, ma, d, md);
//     }
//   }

  switch(s) {
    case -1: // interface integral
      if(d <= 0) {
        return HyperCubeB(n, -1, m, a, ma, d, md);
      }
      else {
        return HyperCubeB(n, -1, m, ma, a, md, d);
      }
      break;
    case 0: // step function integral
      if(d <= 0) {
        return HyperCubeB(n, 0, m, a, ma, d, md);
      }
      else {
        Type I1 = 1;
        for(unsigned i = 0; i <= n; i++) {
          I1 *= 1 / Type(m[i] + 1);
        }
        Type I2 = HyperCubeB(n, 0, m, ma, a, md, d);
        if(std::max(fabs(I1), fabs(I2)) == 0 || fabs(I1 / I2 + 1)  > 1.0e-3) {
          return I1 + I2;
        }
        else  {
          std::cout << "@";
          return HyperCubeB(n, 0, m, a, ma, d, md);
        }

      }
      break;
    default: // all other cases, whatever they mean
      if(d < fabs(a[n])) {
        //std::cout << "!";
        return HyperCubeB(n, s, m, a, ma, d, md);
      }
      else {
        //std::cout << "*";
        return HyperCubeC(n, s, m, a, ma, d, md);
      }
  }



}

template <class Type>
Type HyperCube(const int &s, std::vector<unsigned> m, std::vector <Type> a, Type d) {

//   for(unsigned i = 0; i < a.size(); i++) {
//     d -= a[i];
//     a[i] *= 2;
//   }

  //std::cout<<a.size();
  
  int dim = a.size();
  Type HCI = (dim == 2 ) ? 4 : 8;
  //std::cout << HCI;

  for(int i = a.size() - 1; i >= 0; i--) {
    if(a[i] == 0) {
      HCI *= 1 / Type(m[i] + 1);
      a.erase(a.begin() + i);
      m.erase(m.begin() + i);
    }
  }

  if(a.size() > 0) {
    // all the left a[i] coefficients \ne 0
    for(unsigned i = 0; i < a.size() - 1; i++) { // reorder iterated integral from the smallest to the largest coefficient a[i]
      for(unsigned j = i + 1; j < a.size(); j++) {
        if(fabs(a[i]) > fabs(a[j])) {
          std::swap(a[i], a[j]);
          std::swap(m[i], m[j]);
        }
      }
    }

    std::vector <Type> ma(a.size());
    for(unsigned i = 0; i < a.size(); i++)  ma[i] = -a[i];

    return HCI * HyperCubeA(a.size() - 1, s, m, a, ma, d, -d);
  }
  else {
    return -HCI * LimLi(s, d);
  }
}




#include "TestHyperCube.hpp"

int main(int, char**) {

  typedef double Type;

  std::vector<Type> a = {-1, 1};
  Type norm = sqrt(a[0] * a[0] + a[1] * a[1]);
  a[0] /= norm;
  a[1] /= norm;
  std::vector<Type> ma(2);
  std::vector<unsigned> m = {0, 0};
  Type d = 0. / norm;

  ma[0] = -a[0];
  ma[1] = -a[1];

  std::cout.precision(16);

  std::cout << HyperCube(-1, m, a, d) << std::endl;

  //return 1;

  bool quad = false;//false
  bool hex = true;//false

  Type eps = 5.0e-12;

  if(quad) {
    clock_t time = clock();
    TestQuad(eps);
    std::cout << "Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
  }
  
  if(hex) {
    clock_t time = clock();
    TestHex(eps);
    std::cout << "Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
  }
  return 1;






}























