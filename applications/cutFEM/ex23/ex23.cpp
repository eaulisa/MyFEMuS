#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <algorithm>    // std::sort
#include <ctime>
#include <cstdlib>
#include <climits>
#include <typeinfo>

#include <boost/math/special_functions/factorials.hpp>
//#include <boost/math/special_functions/pow.hpp>
using namespace std;
using boost::math::factorial;


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

template <class Type>
struct Point {
    Type x;
    Type y;
};

template <class Type>
struct Parabola {
    Type k;
    Type b;
    Type d;
};

// Function to calculate the equation of a parabola
template <class Type>
Parabola <Type>  get_parabola_equation( Point <Type> p1, Point <Type> p2, Point <Type> p3 , Type det) {
    Type x1 = p1.x, x2 = p2.x, x3 = p3.x;
    Type y1 = p1.y, y2 = p2.y, y3 = p3.y;
//     Type det = x1 * x1 * (x2 - x3) -x1* (x2*x2 - x3*x3)+ x2*x3*(x2 - x3) ;
//     Type denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
    Type k = (y1 * (x2 - x3) + y2 * (x3 - x1) + y3 * (x1 - x2)) / det;
    Type b = (y1 * (x3*x3 - x2*x2) + y2 * (x1*x1 - x3*x3)+ y3 * ((x2*x2 - x1*x1))) / det;
    Type d = (y1 * x2 * x3 * (x2 - x3) + y2 * x3 * x1 * (x3 - x1) + y3 * x1 * x2 * (x1 - x2)) / det;
    return {k, b, d};
}

template <class Type>
void random_polynomial(std::vector <Type> &a1) {
  a1[0] = ((double(std::rand()) / double(RAND_MAX)) * (20)) - 10;
  a1[1] = ((double(std::rand()) / double(RAND_MAX)) * (20)) - 10;
  a1[2] = ((double(std::rand()) / double(RAND_MAX)) * (20)) - 10;
  a1[3] = ((double(std::rand()) / double(RAND_MAX)) * (20)) - 10;
}
template <class TypeIO, class TypeA>
void GetIntervalall(const std::vector <TypeIO> &a1, const std::vector <TypeIO> &a2, std::vector< std::pair<TypeIO, TypeIO> > &I1, std::vector< std::pair<TypeIO, TypeIO> > &I2, std::vector<std::pair<TypeIO, TypeIO>> &I3) {
  I1.resize(0);
  I2.resize(0);
  I3.resize(0);
  std::vector <TypeA> x(6);
  x[0] = 0;

  unsigned cnt = 1;
  TypeA delta;
  std::vector<TypeA> a = {static_cast<TypeA>(a1[0]), static_cast<TypeA>(a1[1]), static_cast<TypeA>(a1[2])};
  for(unsigned k = 0; k < 2; k++) {
    if( k == 1){
      a = {static_cast<TypeA>(a2[0]), static_cast<TypeA>(a2[1]), static_cast<TypeA>(a2[2])};
    }

    if (a[0] == 0){
      TypeA y = -a[2] / a[1];
      if(y < 1 && y > 0) {
        x[cnt] = y;
        cnt++;
      }
    }
    else {
      delta = static_cast<TypeA>(a[1] * a[1] - 4 * a[0] * a[2]);
      if(delta >0) {
        TypeA sqrtdelta = sqrt(delta);
        int sign = (a[0] > 0) ? 1 : -1;

        for(unsigned i = 0; i < 2; i++) {  //TODO
            TypeA y = (-a[1] - sign * sqrtdelta) / (2 * a[0]);
            if(y >= 1) break;
            else if(y > 0) {
              x[cnt] = y;
              cnt++;
            }
          sign *= -1;
        }
      }
    }
  }
  x[cnt] = 1 ;
  cnt++;

  x.resize(cnt);
  std::sort(x.begin() + 1, x.end() - 1); //TODO
//   for(unsigned i = 0; i < cnt; i++) {
// //     std::cout << "x = " << x[i] << std::endl;
//   }
  for(unsigned i = 0 ; i < cnt - 1 ; i++) {
    TypeA xm = (x[i] + x[i + 1]) / 2;
    TypeA f1 = a1[0] * xm * xm + a1[1] * xm + a1[2] ;
    TypeA f2 = a2[0] * xm * xm + a2[1] * xm + a2[2] ;
    if(f1 > 0) {
      if(f2 > 0) {
        I3.resize(I3.size() + 1, std::pair<TypeIO, TypeIO>(static_cast<TypeIO>(x[i]), static_cast<TypeIO>(x[i + 1])));
      }
      else {
        I1.resize(I1.size() + 1, std::pair<TypeIO, TypeIO>(static_cast<TypeIO>(x[i]), static_cast<TypeIO>(x[i + 1])));
      }
    }
    else if(f2 > 0) {
      I2.resize(I2.size() + 1, std::pair<TypeIO, TypeIO>(static_cast<TypeIO>(x[i]), static_cast<TypeIO>(x[i + 1])));
    }
  }

// Reduce size for adjacent intervals
  for(unsigned i = 1; i < I1.size(); i++) {
    if(I1[i - 1].second == I1[i].first) {
      I1[i - 1].second = I1[i].second;
      I1.erase(I1.begin() + i);
      i--;
    }
  }
  for(unsigned i = 1; i < I2.size(); i++) {
    if(I2[i - 1].second == I2[i].first) {
      I2[i - 1].second = I2[i].second;
      I2.erase(I2.begin() + i);
      i--;
    }
  }
  for(unsigned i = 1; i < I3.size(); i++) {
    if(I3[i - 1].second == I3[i].first) {
      I3[i - 1].second = I3[i].second;
      I3.erase(I3.begin() + i);
      i--;
    }
  }

}

template <class Type>
Type easy_integral_A2(const unsigned &m, const unsigned &n, const int &s, const Type &a, const Type &c, const std::vector <Type> &pol1, const std::vector< std::pair<Type, Type> > &I2) {

  Type A2(0);
  if(a == 0) {
    int rMax = s + n + 1;
    for(unsigned i = 0; i < I2.size(); i++) {
      unsigned r_pm_p1 = m + 1;
      unsigned rMax_mr_pm_p1 = 2 * rMax + m + 1;

//       BEGIN preevaluating
        std::vector <Type> diff_x_pow(2 * rMax + m + 2, 0) ;
//         Type x1pi = I2[i].first;
//         Type x2pi = I2[i].second;
        Type x1pi = pow(I2[i].first,m+1);
        Type x2pi = pow(I2[i].second,m+1);
        for(unsigned pwr = m+1; pwr <= 2 * rMax + m + 1 ; pwr++, x1pi *= I2[i].first, x2pi *= I2[i].second) {
            diff_x_pow[pwr] = (x2pi - x1pi) / (pwr) ;
        }
//         END

    Type pterm = pol1[0] * pol1[2];

      for(int r = 0; r <= rMax; r++) {
        Type term(1);
        Type sum1 = pow(pol1[1], r);
        Type sum2(0);
        unsigned r_p1_m2p =  r + 1;
        unsigned rMax_mr_pp = rMax - r;
        for(int p = 1; p <= r / 2; p++) {
          r_p1_m2p -= 2;
          rMax_mr_pp += 1;
  //           sum += (pow(pol1[0], p) * pow(pol1[1], r - 2 * p) * pow(pol1[2], rMax + p - r)) / (factorial<Type>(p) * factorial<Type>(r - 2 * p) * factorial<Type>(rMax + p - r));
  //           term *= pol1[0] * pol1[2] *(r - 2 * p + 1)*(r - 2 * p + 2) / ( p * (rMax + p - r));
          term *= pterm * r_p1_m2p * (r_p1_m2p + 1) / (p * rMax_mr_pp);
          sum1 += term * pow(pol1[1], r - 2 * p);
        }
        sum1 = sum1 / (factorial<Type>(r) * factorial<Type>(rMax - r));
        sum2 = (r == rMax) ? 0 : sum1 * pow(pol1[0], rMax - r);
        sum1 *= pow(pol1[2], rMax - r);
//         for(unsigned i = 0; i < I2.size(); i++) {

        A2 += sum1 * diff_x_pow[r_pm_p1] +  sum2 * diff_x_pow[rMax_mr_pm_p1] ;
//           A2 += sum1 * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1 +  sum2 * (pow(I2[i].second, rMax_mr_pm_p1) - pow(I2[i].first, rMax_mr_pm_p1)) / rMax_mr_pm_p1 ;
//         }
        r_pm_p1 ++;
        rMax_mr_pm_p1 --;
      }
    }

    A2 *= ((n % 2 == 0) ? -1 : 1) * factorial<Type>(n) / pow(c, n + 1);

    return A2;
  }
  else {
    std::vector <Type> k(3);
    std::cout.precision(20);

    k[0] = pol1[0] / (a * a);
    k[1] = pol1[1] / a;
    k[2] = k[0] * c * c - k[1] * c + pol1[2];
    k[1] -= 2 * c * k[0];

    std::vector <Type> A(s + n + 2, 0);
    std::vector <Type> B(s + n + 2, 0);

    unsigned qMax = s + n + 1;

    //BEGIN pre-evalate A[q] and B[q].
    if(k[1] != 0) {  //regular
      Type kterms = (k[0] * k[2]) / (k[1] * k[1]);
      for(int q = 0; q <= qMax; q++) {
        Type term(1);
        A[q] = term;
        unsigned q_p1_m2r = q + 1;
        unsigned qMax_mq_pr = qMax - q;
        for(int r = 1; r <= q / 2; r++) {
          q_p1_m2r -= 2;
          qMax_mq_pr += 1;
          //term *= k[0] * k[2] * (q - 2 * r + 1) * (q - 2 * r + 2) / (r * (s + n + 1 + r - q) * k[1] * k[1]);
          term *= kterms * q_p1_m2r * (q_p1_m2r + 1) / (r * qMax_mq_pr);
          A[q] += term ;
        }
        A[q] *= pow(k[1], q) / (factorial<Type>(q) * factorial<Type>(qMax - q));
        B[q] = A[q] * (pow(k[0], qMax - q));
        A[q] *= pow(k[2], qMax - q) ;
      }
    }

      //END pre-evalate A[q] and B[q].

    else { // (special case if k[1] =  small )
      for(unsigned w = 0; w < I2.size(); w++)  {
        Type u1 = a * I2[w].first + c;
        Type u2 = a * I2[w].second + c;
//         k[2] = pol1[2] - pol1[1]*c /(2*a);
        for(int p = 0; p <= m; p++) {
          Type sum(0);
          for(int q = 0; q <= qMax; q++) {
            int pwr = 2 * q - n + p ;
            sum += pow(k[2], qMax - q) * pow(k[0], q) / (factorial<Type>(q) * factorial<Type>(qMax - q)) * ((pwr == 0) ? log(u2 / u1) : ((pow(u2, pwr) - pow(u1, pwr)) / (pwr))) ;
          }
          A2 += sum * pow(-c, m - p) / (factorial<Type>(p) * factorial<Type>(m - p)) ;
        }
      }
      A2 *= pow(-1, n + 1) * factorial<Type>(n) * factorial<Type>(m) / pow(a, m + 1) ;

      return A2;
    }


    //integration starts from here.....
    for(unsigned w = 0; w < I2.size(); w++)  {
      Type u1 = a * I2[w].first + c;
      Type u2 = a * I2[w].second + c;

      if(u1 == 0 || u2 == 0) {   // TODO need to fix this. if we have double parts of any region. This is wrong .

        Type A2i(0);
        Type c_0 = (a * pol1[1] - pol1[0] * c) / (a * a);
        int pMax = s + n + 1 ;
        Type p_term(1);
        Type p_sum(1);
        for(int p = 1; p <= pMax; p++) {
          Type q_term(1);
          Type q_sum = q_term;

          for(int q = 1; q <= s; q++) {
            Type r_pm_p1 = p + q + 1;
            q_term *= a * (s - q + 1) / (c * q);
            q_sum += q_term * (pow(I2[w].second, r_pm_p1) - pow(I2[w].first, r_pm_p1)) / r_pm_p1;
          }
          q_sum *= pow(c, s) / factorial<Type>(s) ;
          p_term *= pol1[0] * (pMax - p + 1) / (a * p * c_0);
          p_sum += p_term * q_sum ;
        }

        A2 += p_sum * ((n % 2 == 0) ? -1 : 1) * factorial<Type>(n) * factorial<Type>(s) * pow(c_0, pMax) / factorial<Type>(pMax) ;
      }
      else {

//        {
// //         Type A2i(0);
// //         for(unsigned p = 0; p <= m; p++) {
// //           Type sum1(0);
// //           for(unsigned q = 0; q <= qMax; q++) {
// //             int pwr = p + q - n;
// //             sum1 += A[q] * ((pwr== 0) ? log(u2 / u1) : (pow(u2, pwr) - pow(u1, pwr)) / (pwr));
// //           }
// //           Type sum2(0);
// //           for(unsigned q = 0; q < qMax; q++) {
// //             int pwr= 2 * s + n + 2 + p - q;
// //             sum2 += B[q] * (pow(u2,pwr) - pow(u1,pwr)) / (pwr);
// //           }
// //           A2i += (sum1 + sum2) * pow(-c, m - p) / (factorial<Type>(p) * factorial<Type>(m - p));
// //         }
// //         A2 += A2i * pow(-1, n + 1) * factorial<Type>(n) * factorial<Type>(m) / pow(a, m + 1) ;
//        }

        // BEGIN pre evalution of power of U
        std::vector <Type> diff_u_pow(2 * s + 2 * n + m + 3, 0) ;
        Type u1pi = 1. / pow(u1, n);
        Type u2pi = 1. / pow(u2, n);
        for(int pwr = 0; pwr <= n - 1 ; pwr++, u1pi *= u1, u2pi *= u2) {
          int actual_pwr = pwr - n;
          diff_u_pow[pwr] = (u2pi - u1pi) / actual_pwr ;
        }
//             Type u1pi = 1./u1;
//             Type u2pi = 1./u2;
//             for(int pwr = n-1; pwr >= 0 ; pwr--, u1pi /= u1, u2pi /= u2) {
//               int actual_pwr = pwr - n;
//               diff_u_pow[pwr] = (u2pi - u1pi) / actual_pwr ;
//             }

        diff_u_pow[n] = log(u2 / u1) ;
        u1pi = u1;
        u2pi = u2;
        for(int pwr = n + 1; pwr <= 2 * qMax + m ; pwr++, u1pi *= u1, u2pi *= u2) {
          int actual_pwr = pwr - n;
          diff_u_pow[pwr] = (u2pi - u1pi) / actual_pwr ;
        }
        // END pre evaluation of power


        Type A2i(0);
        for(int p = 0; p <= m; p++) {
          Type sum1(0);
          for(int q = 0; q <= qMax; q++) {
//             int pwr = p + q;                       // added n with original power
//             sum1 += A[q] * diff_u_pow[pwr] ;
            sum1 += A[q] * diff_u_pow[p + q] ;

          }
          Type sum2(0);
          for(int q = 0; q < qMax; q++) {
//             int pwr = 2 * qMax + p - q;            // added n with original power
//             sum2 += B[q] * diff_u_pow[pwr];
            sum2 += B[q] * diff_u_pow[2 * qMax + p - q];
          }
          A2i += (sum1 + sum2) * pow(-c, m - p) / (factorial<Type>(p) * factorial<Type>(m - p));
        }
        A2 += A2i * ((n % 2 == 0) ? -1 : 1) */* pow(-1, n + 1) **/ factorial<Type>(n) * factorial<Type>(m) / pow(a, m + 1) ;
      }

    }
    //total
    //  A2 *= pow(-1, n + 1) * factorial<Type>(n) * factorial<Type>(m) / pow(a, m + 1);
// //           std::cout << "final. A2= " << A2 << std::endl;
    return A2;
  }
}


template <class Type>
Type easy_integral_A3(const unsigned &m, const unsigned &n, const int &s, const Type &a, const Type &c, const std::vector <Type> &pol1, const std::vector< std::pair<Type, Type> > &I3) {
  Type A3(0);
  if(a == 0) {
    for(int i = 0; i <= s; i++) {
      for(unsigned w = 0; w < I3.size(); w++) {
        int pMax = s - i;
        for(int r = 0; r <= pMax; r++) {
          Type sum1 = pow(pol1[2], pMax - r)/factorial<Type>(pMax - r) ;
          Type sum = sum1 * (pow(I3[w].second,2 * r + m + 1) - pow(I3[w].first, 2 * r + m + 1))/ (2 * r + m + 1);
          for(int p = 1; p <= pMax - r; p++) {
            sum1 *= pol1[1] * (pMax - r - p + 1) / (pol1[2] * p)   ;
            sum += sum1 * (pow(I3[w].second,2 * r + m+p + 1) - pow(I3[w].first, 2 * r +p+ m + 1))/ (2 * r + m+p + 1) ;
          }
          A3 += sum * pow(pol1[0], r) / (factorial<Type>(r) * (n + i + 1) * factorial<Type>(i)) ;
        }
      }
    }
  }

//     if(a == 0) { //preevaluate this but it is expensive!!!!
//     for(int i = 0; i <= s; i++) {
//       for(unsigned w = 0; w < I3.size(); w++) {
//         //BEGIN preevaluating
// //         std::vector <Type> diff_x_pow(2 * s + 2, 0) ;
// //         Type x1pi = pow(I3[w].first,m+1);
// //         Type x2pi = pow(I3[w].second,m+1);
// //         for(unsigned pwr = m+1; pwr <= 2 * s + m + 1 ; pwr++, x1pi *= I3[w].first, x2pi *= I3[w].second) {
// //             //               diff_u_pow[pwr] = (pow(u2, pwr) - pow(u1, pwr)) / (pwr) ;
// //             diff_x_pow[pwr-m] = (x2pi - x1pi) / (pwr) ;
// //         }
//         //END
//         int pMax = s - i;
//         // #1
//         for(int r = 0; r <= pMax; r++) {
//
//           Type sum1 = pow(pol1[2], pMax - r)/factorial<Type>(pMax - r) ;
// //           Type sum = sum1 * diff_x_pow[2 * r + 1];
// //           int r_pm_p1 = 2 * r + m + 1;
// //           Type sum = sum1 * (pow(I3[w].second,r_pm_p1) - pow(I3[w].first, r_pm_p1))/ (r_pm_p1);
//           Type sum = sum1 * (pow(I3[w].second,2 * r + m + 1) - pow(I3[w].first, 2 * r + m + 1))/ (2 * r + m + 1);
//
//           for(int p = 1; p <= pMax - r; p++) {
// //             r_pm_p1++;
//             sum1 *= pol1[1] * (pMax - r - p + 1) / (pol1[2] * p)   ;
// //             sum += sum1 * diff_x_pow[2 * r + p + 1] ;
//             sum += sum1 * (pow(I3[w].second,2 * r + m+p + 1) - pow(I3[w].first, 2 * r +p+ m + 1))/ (2 * r + m+p + 1) ;
// //          sum += sum1 * (pow(I3[w].second,r_pm_p1) - pow(I3[w].first, r_pm_p1))/(r_pm_p1) ;
// //             sum += pow(pol1[1], p) * pow(pol1[2], pMax - r - p) * diff_x_pow[r_pm_p1] / (factorial<Type>(p) * factorial<Type>(pMax - r - p) ) ;
//           }
//           A3 += sum * pow(pol1[0], r) / (factorial<Type>(r) * (n + i + 1) * factorial<Type>(i)) ;
//         }
//       }
//     }
//   }
//     if(a == 0) { // TODO optimize
//     for(int i = 0; i <= s; i++) {
//       for(unsigned w = 0; w < I3.size(); w++) {
//         int pMax = s - i;
//         // #1
//         for(int r = 0; r <= pMax; r++) {
//           Type sum = 0;
//
//           for(int p = 0; p <= pMax - r; p++) {
//             Type r_pm_p1 = 2 * r + p + m + 1;
//             sum += pow(pol1[1], p) * pow(pol1[2], s - i - r - p) * (pow(I3[w].second, r_pm_p1) - pow(I3[w].first, r_pm_p1)) / (factorial<Type>(p) * factorial<Type>(s - i - r - p) * r_pm_p1) ;
//           }
//           A3 += sum * pow(pol1[0], r) / (factorial<Type>(r) * (n + i + 1) * factorial<Type>(i)) ;
//         }
//       }
//     }
//   }

  else {
    std::vector <Type> k(3);
    k[0] = pol1[0] / (a * a);
    k[1] = pol1[1] / a;
    k[2] = k[0] * c * c - k[1] * c + pol1[2];
    k[1] -= 2 * c * k[0];

    for(int i = 0; i <= s; i++) {
      unsigned qMax = s - i;
      if(k[1] == 0) {   // if k[1] is small
        for(unsigned w = 0; w < I3.size(); w++) {
          Type u1 = a * I3[w].first + c;
          Type u2 = a * I3[w].second + c;

          // BEGIN pre evalution of power of all
          std::vector <Type> diff_u_pow(m + 2 * s + 2, 0) ;
          Type u1pi = u1;
          Type u2pi = u2;
          for(unsigned pwr = 1; pwr <= m + 2 * s + 1 ; pwr++, u1pi *= u1, u2pi *= u2) {
            //               diff_u_pow[pwr] = (pow(u2, pwr) - pow(u1, pwr)) / (pwr) ;
            diff_u_pow[pwr] = (u2pi - u1pi) / (pwr) ; // TODO TOCHECK
          }
          std::vector <Type> pow_c(m + 1, 0) ;
          pow_c[0] = 1;
          for(unsigned pwr = 1; pwr <= m ; pwr++) {
            pow_c[pwr] = (-c) * pow_c[pwr - 1] ;
          }
          std::vector <Type> pow_k0(s + 1, 0) ;
          std::vector <Type> pow_k2(s + 1, 0);
          pow_k0[0] = 1;
          pow_k2[0] = 1;
          for(unsigned pwr = 1; pwr <= s ; pwr++) {
            pow_k0[pwr] = k[0] * pow_k0[pwr - 1] ;
            pow_k2[pwr] = k[0] * pow_k2[pwr - 1] ;
          }
          // END pre evalution of power of all

          //         k[2] = pol1[2] - pol1[1]*c /(2*a);
          for(int p = 0; p <= m; p++) {
            Type sum(0);
            for(int q = 0; q <= qMax; q++) {
              int pwr = 2 * q + i + p + 1 ;
              sum += pow_k2[qMax - q] * pow_k0[q] * diff_u_pow[pwr] / (factorial<Type>(q) * factorial<Type>(qMax - q))  ;
            }
            A3 += sum * pow_c[m - p] / (factorial<Type>(p) * factorial<Type>(m - p)) ;
          }
        }
      }

      else { // main integral
        // BEGIN pre evaluation A[q] and B[q]
        std::vector <Type> A(s - i + 1, 0);  // size of all this vector changes.
        std::vector <Type> B(s - i + 1, 0);
        if(k[1] != 0) {
          Type kterms = (k[0] * k[2]) / (k[1] * k[1]);
          for(int q = 0; q <= qMax; q++) {
            Type term(1);
            A[q] = term;
            unsigned q_p1_m2r = q + 1;
            unsigned qMax_mq_pr = qMax - q;

            for(int r = 1; r <= q / 2; r++) {
              q_p1_m2r -= 2;
              qMax_mq_pr += 1;
              //term *= k[0] * k[2] * (q - 2 * r + 1) * (q - 2 * r + 2) / (r * (s + n + 1 + r - q) * k[1] * k[1]);
              term *= kterms * q_p1_m2r * (q_p1_m2r + 1) / (r * qMax_mq_pr);
              A[q] += term ;
            }
            //           B[q] = A[q] * (pow(k[1], q) * pow(k[0], qMax - q)) / (factorial<Type>(q) * factorial<Type>(qMax - q));
            //           A[q] *= (pow(k[1], q) * pow(k[2], qMax - q)) / (factorial<Type>(q) * factorial<Type>(qMax - q));

            A[q] *= pow(k[1], q) / (factorial<Type>(q) * factorial<Type>(qMax - q));
            B[q] = A[q] * pow(k[0], qMax - q);
            A[q] *= pow(k[2], qMax - q);

            //         std::cout<<"A["<<q<<"] = " << A[q] <<"  B[] ="<< B[q] << std::endl;
            //         std::cout << "A[" << q << "] = " << A[q] << "  B[] =" << B[q] << std::endl;
          }
        }
        // END  pre evaluation

        for(unsigned w = 0; w < I3.size(); w++) {
          Type u1 = a * I3[w].first + c;
          Type u2 = a * I3[w].second + c;

          // BEGIN pre evalution of power of U
          std::vector <Type> diff_u_pow(m + 2 * s + 2, 0) ;
          Type u1pi = u1;
          Type u2pi = u2;
          for(unsigned pwr = 1; pwr <= m + 2 * s + 1 ; pwr++, u1pi *= u1, u2pi *= u2) {
            //               diff_u_pow[pwr] = (pow(u2, pwr) - pow(u1, pwr)) / (pwr) ;
            diff_u_pow[pwr] = (u2pi - u1pi) / (pwr) ;
          }
          // END
          // BEGIN pre evalution of power of -c
          std::vector <Type> pow_c(m + 1, 0) ;
          pow_c[0] = 1;
          for(unsigned pwr = 1; pwr <= m ; pwr++) {
            pow_c[pwr] = (-c) * pow_c[pwr - 1] ;
          }
          // END pre evalution of power of -c

          //Type A3i(0);
          for(unsigned p = 0; p <= m; p++) {
            Type sum1(0);
            int pwr = p + i + 1;
            for(unsigned q = 0; q <= qMax; q++, pwr++) {
              //             int pwr = p + q + i + 1;
              sum1 += A[q] * diff_u_pow[pwr];
            }
            Type sum2(0);
            pwr = 2 * s - i + p + 1;
            for(unsigned q = 0; q < qMax; q++, pwr--) {
              //int pwr = 2 * s - i + p - q + 1;
              sum2 += B[q] * diff_u_pow[pwr];
            }
            A3 += (sum1 + sum2) * pow_c[m - p] / (factorial<Type>(p) * factorial<Type>(m - p));
          }
        }
      }
      A3 /= ((n + i + 1) * factorial<Type>(i)) ;
    }
    A3 *= factorial<Type>(m) / pow(a, m + 1);
    //     std::cout<< "final. A3= "<< A3 << std::endl;
  }
  return A3;
}



template <class Type>
void creat_parabola_table(std::vector< std::vector< std::vector< std::vector< Type >>>> &parabola_table , const int &partition, const unsigned &m, const unsigned &n, const int &s){

  Type area ;
  Type a(0);
  std::vector <Type> pol1(3, 0);
  std::vector <Type> pol2(3, 0);
  unsigned int count;
  Point <Type> p1, p2, p3 ;
//   cout << " Number | (x1,y1) | (x2,y2) | (x3,y3) |  k  |  b  |  d  |  c  | Area |" <<endl;
  double del_x = 1./partition;
//   cout << "del_x " << del_x << endl;
  for (int table = 0; table <=7; table++){  // BEGIN preevaluation of the table
//     cout << "Table " << table << endl;
    parabola_table.resize(parabola_table.size() + 1);
    parabola_table[table].resize(0);
    count = 0;

    for (double i1=0.;i1<=1.;i1+= del_x){
      for (double i2=0.;i2<=1.;i2+=del_x){
        for (double i3=0.;i3<=1.;i3+=del_x){
//           cout << " i3 = " << i3 << endl;
           switch (table) {
                case 0:
                    p1 = {0, i1};
                    p2 = {i2, 1};
                    break;
                case 1:
                    p1 = {0, i1};
                    p2 = {1,i2};
                    break;
                case 2:
                    p1 = {0, i1};
                    p2 = {i2, 0};
                    break;
                case 3:
                    p1 = {i1,1};
                    p2 = {i2, 1};
                    if(i2 >= i1){p1 = {0, 0};p2 = {0, 0};}
                    break;
                case 4:
                    p1 = {i1,1};
                    p2 = {1, i2};
                    break;
                case 5:
                    p1 = {i1,1};
                    p2 = {i2, 0};
                    break;
                case 6:
                    p1 = {1, i1};
                    p2 = {i2, 0};
                    break;
                case 7:
                    p1 = {i1, 0};
                    p2 = {i2, 0};
                    if(i2 >= i1){p1 = {0, 0};p2 = {0, 0};}
                    break;
           }
           p3 = {0.5 * (p1.x + p2.x) , i3 };

           p1 = {static_cast<Type>(p1.x), static_cast<Type>(p1.y)};
           p2 = {static_cast<Type>(p2.x), static_cast<Type>(p2.y)};
           p3 = {static_cast<Type>(p3.x), static_cast<Type>(p3.y)};
           Type c(-1) ;
           Type det = p1.x * p1.x * (p2.x - p3.x) -p1.x* (p2.x*p2.x - p3.x*p3.x)+ p2.x*p3.x*(p2.x - p3.x) ;// only sort out the points parallel to y axis
           if (det !=0){
            parabola_table[table].resize(parabola_table[table].size() + 1);
            parabola_table[table][count].resize(2);

            for(int normal=0; normal <=1; normal++ ){
                Type A1 (0), A2 (0), A3 (0);
                Parabola <Type> parabola = get_parabola_equation(p1, p2, p3, det);
                if (normal == 1){ // To calculate other side
                    parabola.k *= -1;
                    parabola.b *= -1;
                    parabola.d *= -1;
                    c *= -1;
                }

                pol1[0] = parabola.k;
                pol1[1] = parabola.b;
                pol1[2] = parabola.d + c;
                pol2[0] = parabola.k;
                pol2[1] = parabola.b;
                pol2[2] = parabola.d;


                std::vector< std::pair <Type, Type> > I1, I2, I3 ;
                GetIntervalall<Type, double>(pol1, pol2, I1, I2, I3);
                if(I1.size() > 0) {
                    A1 = easy_integral_A3(m, n, s, a, c, pol2, I1) -  easy_integral_A2(m, n, s, a, c, pol2, I1);
                }
                if(I2.size() > 0) {
                    A2 = easy_integral_A2(m, n, s, a, c, pol2, I2);
                }
                if(I3.size() > 0) {
                    A3 = easy_integral_A3(m, n, s, a, c, pol2, I3);
                }
                area = A1 + A2 + A3;
//                 cout << " \n" << count << ". (" << p1.x << ", " << p1.y << ") ," << "(" << p2.x << ", " << p2.y << ") ," << "(" << p3.x << ", " << p3.y << ")  : "  <<  parabola.k << ", " << parabola.b << ", " << parabola.d << ", " << c << ", " << area<<  endl;

//                 parabola_table[table].resize(parabola_table[table].size() + 1);
//                 parabola_table[table][count].resize(2);
                parabola_table[table][count][normal].resize(15);

                parabola_table[table][count][normal][0] = count;
                parabola_table[table][count][normal][1] = static_cast<Type>(i1);
                parabola_table[table][count][normal][2] = static_cast<Type>(i2);
                parabola_table[table][count][normal][3] = static_cast<Type>(i3);

                parabola_table[table][count][normal][4] = p1.x;
                parabola_table[table][count][normal][5] = p1.y;
                parabola_table[table][count][normal][6] = p2.x;
                parabola_table[table][count][normal][7] = p2.y;
                parabola_table[table][count][normal][8] = p3.x;
                parabola_table[table][count][normal][9] = p3.y;

                parabola_table[table][count][normal][10] = parabola.k;
                parabola_table[table][count][normal][11] = parabola.b;
                parabola_table[table][count][normal][12] = parabola.d;
                parabola_table[table][count][normal][13] = c;
                parabola_table[table][count][normal][14] = area;

//                 for (int nm=0; nm<=11; nm++){
//                 cout << parabola_table[normal][count][nm] << " " ;
//                 }
            }
            count ++ ;
          }
        }
      }
    }
  }

}

template <class Type>
void inverse_parabola(std::vector< std::vector< std::vector< std::vector< Type >>>> &parabola_table ,const std::vector <Type> &given_parabola, vector <Type> &intersection, int &normal, std::vector< Type > &interp_point, std::vector< std::vector< Type >> & interp_table, unsigned int &table_number, const int &partition){

    interp_point.resize(0);
    intersection.resize(0);

    Type k = given_parabola[0];
    Type b = given_parabola[1];
    Type d = given_parabola[2];
    Type c = given_parabola[3];
    bool left=0 , top = 0 , right = 0 , bottom = 0;
     if (normal == 0){ k = k/c; b=b/c; d=d/c; c=-1;}
    else if (normal == 1){ k = -k/c; b=-b/c; d=-d/c; c=1;}



      if (d>=0 && d<=1){ //LEFT
        intersection.resize(intersection.size()+2);
        intersection[intersection.size()-2] = 0;
        intersection[intersection.size()-1] = d;
        left = 1 ;
        interp_point.resize(interp_point.size()+1);
        interp_point[interp_point.size()-1] = d;
      }
      // LEFT-TOP solve kx^2+bx+d-1 =0  ; Table 0
      if (k == 0){
      Type x =  (1-d)/b ;
        if(x < 1 && x> 0) {
          intersection.resize(intersection.size()+2);
          intersection[intersection.size()-2] = x;
          intersection[intersection.size()-1] = 1;
          interp_point.resize(interp_point.size()+1);
          interp_point[interp_point.size()-1] = x;
          top =1;
          if (left ==1) table_number = 0 ;
        }
      }
      else {
        Type delta = b*b - 4*k*(d+c);
        if(delta >0) {
          Type sqrtdelta = sqrt(delta);
          int sign = (k > 0) ? 1 : -1;

          for(unsigned i = 0; i < 2; i++) {
            Type x = (- b - sign * sqrtdelta) / (2 * k);
            cout<< "Top x = "<< x<< endl;
            if(x > 1) break;
            else if(x > 0) {
              intersection.resize(intersection.size()+2);
              intersection[intersection.size()-2] = x;
              intersection[intersection.size()-1] = 1;
              interp_point.resize(interp_point.size()+1);
              interp_point[interp_point.size()-1] = x;
              if (top ==1){table_number = 3 ;}
              top = 1;
              if (left ==1){table_number = 0 ;break;}

            }
            sign *= -1;
          }
        }
      }
      Type y_1=k+b+d; //LEFT-RIGHT x=1 ; Table 1
      if (interp_point.size() < 2 && y_1 >= 0 && y_1 <= 1){ //TODO check sign when normal change
          intersection.resize(intersection.size()+2);
          intersection[intersection.size()-2] = 1;
          intersection[intersection.size()-1] = y_1;
          interp_point.resize(interp_point.size()+1);
          interp_point[interp_point.size()-1] = y_1;
          if (left ==1){table_number = 1 ;}
          if  (top ==1){table_number = 4 ;}
          right = 1 ;
      }

      if (interp_point.size() < 2){ //LEFT-BOTTOM  solve kx^2+bx+d =0 ; Table 2
        if (k == 0){
          Type x =  -d/b ;
          if(x < 1 && x> 0) {
            intersection.resize(intersection.size()+2);
            intersection[intersection.size()-2] = x;
            intersection[intersection.size()-1] = 0;
            interp_point.resize(interp_point.size()+1);
            interp_point[interp_point.size()-1] = x;
            if (left ==1){table_number = 2 ;}
            if (right ==1){table_number = 6 ;}
            if (top ==1){table_number = 5 ;}
          }
        }

        else {
          Type delta = b*b - 4*k*d;
          if(delta >0) {
            Type sqrtdelta = sqrt(delta);
            int sign = (k > 0) ? 1 : -1;

            for(unsigned i = 0; i < 2; i++) {
              Type x = (- b - sign * sqrtdelta) / (2 * k);
              if(x > 1) break;
              else if(x > 0) {
                intersection.resize(intersection.size()+2);
                intersection[intersection.size()-2] = x;
                intersection[intersection.size()-1] = 0;
                interp_point.resize(interp_point.size()+1);
                interp_point[interp_point.size()-1] = x;
                if (bottom ==1){table_number = 7 ;}
                if (left ==1){table_number = 2 ;}
                if (right ==1){table_number = 6 ;}
                if (top ==1){table_number = 5 ;}    // TODO check the table
                bottom = 1;
              }
              sign *= -1;
            }
          }
        }
      }

      if (interp_point.size() == 2){  // finding p3
            Type mid_point = 0.5*(intersection[0]+intersection[2]);
            intersection.resize(intersection.size()+2);
            intersection[intersection.size()-2] = mid_point;
            intersection[intersection.size()-1] = k*mid_point*mid_point + b*mid_point + d;


            interp_point.resize(3);
            interp_point[2] = k*mid_point*mid_point + b*mid_point + d;
            if (interp_point[2]<0 && interp_point[2]>1) cout << " parabola at midpoint is outside : check your intersections" << endl;


//       cout<< " given parabola 55= " << k << " " << b << " "<< d << " "<< c <<endl;
      }
//       else cout<<"some thing is wrong in intersection points : dealing with more or less then two intersection points " << interp_point.size()<< "table " << table_number <<endl;

      cout<< " given parabola : k= " << k << "; b= " << b << "; d= "<< d << "; c = "<< c <<";"<< endl;
      cout<< " inter section line = " << " LEFT " << left << " Top "<< top << " Right "<< right << " bottom " << bottom <<endl;

      //finding closed points
      if (interp_point.size() == 3){
        Type x0,x1,y0,y1,z0,z1;
        x0 = floor(interp_point[0]* partition) / partition;
        x1 = ceil(interp_point[0]* partition) / partition;
        y0 = floor(interp_point[1]* partition) / partition;
        y1 = ceil(interp_point[1]* partition) / partition;
        z0 = floor(interp_point[2]* partition) / partition;
        z1 = ceil(interp_point[2]* partition) / partition;
//         cout<< " intersection points = " << interp_point[0] << " " << interp_point[1] << " "<< interp_point[2] <<endl;
//         cout<< x0 << " " << x1 << " "<< y0 << " " << y1 << " "<< z0 << " " << z1 << " "<< endl;

        interp_table[0] = {x0,y0,z0,-1 };
        interp_table[1] = {x0,y0,z1,-2 };
        interp_table[2] = {x0,y1,z0,-3 };
        interp_table[3] = {x0,y1,z1,-4 };
        interp_table[4] = {x1,y0,z0,-5 };
        interp_table[5] = {x1,y0,z1,-6 };
        interp_table[6] = {x1,y1,z0,-7 };
        interp_table[7] = {x1,y1,z1,-8 };

//               cout<< " parabola_table[table_number].size() = " << parabola_table[table_number].size() << " " <<endl;
  //       interp_table = {x0,y0,z0,-1,
  //                       x0,y0,z1,-2,
  //                       x0,y1,z0,-3,
  //                       x0,y1,z1,-4,
  //                       x1,y0,z0,-5,
  //                       x1,y0,z1,-6,
  //                       x1,y1,z0,-7,
  //                       x1,y1,z1,-8 };
        for(unsigned int i = 0; i <=7 ; i++){
          for (unsigned int count_x = 0; count_x <= parabola_table[table_number].size(); count_x++ ){
            if(fabs(parabola_table[table_number][count_x][normal][1] - interp_table[i][0]) < 0.00000000001){
//                   cout<< "count x =" << count_x << " " << parabola_table[table_number][count_x][normal][1] << " " << interp_table[i][0]  <<endl;
              for (unsigned int count_y = count_x; count_y <= parabola_table[table_number].size(); count_y++ ){
                if(fabs(parabola_table[table_number][count_y][normal][2] - interp_table[i][1]) < 0.000000000001){
//                   cout<< "count y =" << count_y <<endl;
                  for (unsigned int count_z = count_y; count_z <= parabola_table[table_number].size(); count_z++ ){
                    if(fabs(parabola_table[table_number][count_z][normal][3] - interp_table[i][2]) < 0.000000000001 ){
//                       cout<< "count z =" << count_z <<endl;
                      interp_table[i][3] = parabola_table[table_number][count_z][normal][14] ;
//                       cout << " area = "<< parabola_table[table_number][count_z][normal][14];
                      count_y = parabola_table[table_number].size()+1;
                      count_x = parabola_table[table_number].size()+1;
                      break;
                    }
                  }
                }
              }
            }
          }
        }
  //                 parabola_table[table][count][normal][0]


//         cout<< " given parabola = " << k << " " << b << " "<< d << " "<< c <<endl;
        cout<< " intersection points = " << interp_point[0] << " " << interp_point[1] << " "<< interp_point[2] <<endl;
        cout << " interpolation table = " << endl;
        for(unsigned int i = 0; i <=7 ; i++){
          cout << interp_table[i][0] << " " <<interp_table[i][1] << " " <<interp_table[i][2] << " " <<interp_table[i][3]<< endl;
        }
      }


  // END inverse function //TODO I have to be carefull with special cases.

  }



template <class Type>
Type trilinier_interpolation(std::vector< std::vector< Type >> & interp_table , const std::vector< Type > &interp_point, const int &partition ){

  Type x = interp_point[0];
  Type y = interp_point[1];
  Type z = interp_point[2];

  Type x0 = floor(interp_point[0] * partition) / partition;
  Type x1 = ceil(interp_point[0]  * partition) / partition;
  Type y0 = floor(interp_point[1] * partition) / partition;
  Type y1 = ceil(interp_point[1]  * partition) / partition;
  Type z0 = floor(interp_point[2] * partition) / partition;
  Type z1 = ceil(interp_point[2]  * partition) / partition;

  Type x_d = (x-x0)/(x1-x0);
  Type y_d = (y-y0)/(y1-y0);
  Type z_d = (z-z0)/(z1-z0);

  Type c_00 = interp_table[0][3] * (1-x_d) + interp_table[4][3] * x_d ;
  Type c_01 = interp_table[1][3] * (1-x_d) + interp_table[5][3] * x_d ;
  Type c_10 = interp_table[2][3] * (1-x_d) + interp_table[6][3] * x_d ;
  Type c_11 = interp_table[3][3] * (1-x_d) + interp_table[7][3] * x_d ;

  Type c_0 = c_00 * (1-y_d) + c_10 * y_d ;
  Type c_1 = c_01 * (1-y_d) + c_11 * y_d ;

  Type c = c_0 * (1-z_d) + c_1 * z_d ;

  return c;
}




int main() {

  typedef cpp_bin_float_oct Type;      //     typedef double Type;  //  std::cout.precision(20);
  unsigned int m = 0;
  unsigned int n = 0;
  int s = 0;
  unsigned int partition = 15;

  std::vector< std::vector< std::vector< std::vector< Type >>>> parabola_table(0) ;

  std::srand(10);

  creat_parabola_table<Type>(parabola_table, partition, m,n,s);

  clock_t t = clock();

    std::vector <Type> given_parabola(4);
    std:: vector <Type> intersection(0);
    std:: vector <Type> interp_point(0);
    std:: vector <std:: vector <Type>> interp_table(8, std:: vector <Type> (4));
    int normal = 0;
    unsigned int table_number;
    Type interp_area;
    for (unsigned i=0;i<20;i++){
      cout<< " "<<i<<endl;
      random_polynomial(given_parabola);
      inverse_parabola<Type>(parabola_table , given_parabola, intersection, normal, interp_point, interp_table, table_number, partition);
        if ( interp_point.size() == 3) {
           interp_area = trilinier_interpolation<Type>( interp_table , interp_point, partition) ;
           cout << " interpolated area = " << interp_area << "\n"<<endl;
        }
        else if (interp_point.size() == 0){
          interp_area = (given_parabola[2] >0)? 0 : 1 ;
          cout << " NO INTERSECTION : area = " << interp_area << "\n"<<endl;
        }
        else cout<<"something is not right about intersection points. Total intersection : " << interp_point.size() <<endl;
    }


  for (int table = 0; table <=7; table++){
//     cout << "Table " << table << endl;
    cout<< parabola_table[table].size() << " + " ;
  }
  t = clock() - t;
  std::cout << "Time taken " << (Type)(t) / CLOCKS_PER_SEC << std::endl;

    return 0;
}




