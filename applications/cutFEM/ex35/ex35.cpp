
#include <typeinfo>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <fstream>
#include <vector>




// #include <iostream>
// #include <algorithm>
// #include <vector>
// #include <cmath>
// #include <algorithm>    // std::sort
// #include <ctime>
// #include <cstdlib>
// #include <climits>

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
struct PointT {
  Type x;
  Type y;
};

template <class Type>
struct Parabola {
  Type k;
  Type b;
  Type d;
};


#include "cutFemWeightParabola.hpp"
#include "Rebuild.hpp"
#include "PolynomialBases.hpp"
#include "Fem.hpp"


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
    if(k == 1) {
      a = {static_cast<TypeA>(a2[0]), static_cast<TypeA>(a2[1]), static_cast<TypeA>(a2[2])};
    }

    if(a[0] == 0) {
      TypeA y = -a[2] / a[1];
      if(y < 1 && y > 0) {
        x[cnt] = y;
        cnt++;
      }
    }
    else {
      delta = static_cast<TypeA>(a[1] * a[1] - 4 * a[0] * a[2]);
      if(delta > 0) {
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
Parabola <Type>  get_parabola_equation(const PointT <Type> p1, const PointT <Type> p2, const PointT <Type> p3) {
  Type x1 = p1.x, x2 = p2.x, x3 = p3.x;
  Type y1 = p1.y, y2 = p2.y, y3 = p3.y;
  Type k, b, d;
  Type det = p1.x * p1.x * (p2.x - p3.x) - p1.x * (p2.x * p2.x - p3.x * p3.x) + p2.x * p3.x * (p2.x - p3.x);

//     Type det = x1 * x1 * (x2 - x3) -x1* (x2*x2 - x3*x3)+ x2*x3*(x2 - x3) ;
//     Type denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
  if(fabs(det) > 0.0000000000000001) {
    k = -(y1 * (x2 - x3) + y2 * (x3 - x1) + y3 * (x1 - x2)) / det;
    b = -(y1 * (x3 * x3 - x2 * x2) + y2 * (x1 * x1 - x3 * x3) + y3 * ((x2 * x2 - x1 * x1))) / det;
    d = -(y1 * x2 * x3 * (x2 - x3) + y2 * x3 * x1 * (x3 - x1) + y3 * x1 * x2 * (x1 - x2)) / det;
  }
  else {
    Type slope = (y1 - y2) / (x1 - x2) ;
    k = 0;
    b = -slope;
    d = x2 * slope - y2 ;
  }

//     else {This will give us a straight line paralal to y axix.
//        k = -0.;
//        b = -1.;
//        d = p1.x ;
//     }

  if(fabs(k) < 1.e-14) k = 0 ;
  if(fabs(b) < 1.e-14) b = 0 ;
  if(fabs(d) < 1.e-14) d = 0 ;

  return {k, b, d};
}

template <class Type>
Type trig_integral_A2(const unsigned &m, const unsigned &n, const int &s, const Type &a, const Type &c, const std::vector <Type> &pol1, const std::vector< std::pair<Type, Type> > &I2) {
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
      Type x1pi = pow(I2[i].first, m + 1);
      Type x2pi = pow(I2[i].second, m + 1);
      for(unsigned pwr = m + 1; pwr <= 2 * rMax + m + 1 ; pwr++, x1pi *= I2[i].first, x2pi *= I2[i].second) {
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
//change the sign here. Why did I put -1 for even n before? I don't understand. But for now lets change the sign
    A2 *= ((n % 2 == 0) ? -1 : 1) * factorial<Type>(n) / pow(c, n + 1);

    return A2;
  }

  //Did not do any change in here.

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
Type trig_integral_A3(const unsigned &m, const unsigned &n, const int &s, const Type &a, const Type &c, const std::vector <Type> &pol1, const std::vector< std::pair<Type, Type> > &I3) {
  Type A3(0);
  if(a == 0) {
    for(int i = 0; i <= s; i++) {
      for(unsigned w = 0; w < I3.size(); w++) {
        int pMax = s - i;
        for(int r = 0; r <= pMax; r++) {
          Type sum1 = pow(pol1[2], pMax - r) / factorial<Type>(pMax - r) ;
          Type sum = sum1 * (pow(I3[w].second, 2 * r + m + 1 + n + 1 - i) - pow(I3[w].first, 2 * r + m + 1 + n + 1 - i)) / (2 * r + m + 1 + n + 1 - i);
          for(int p = 1; p <= pMax - r; p++) {
            sum1 *= pol1[1] * (pMax - r - p + 1) / (pol1[2] * p)   ;
            sum += sum1 * (pow(I3[w].second, 2 * r + m + p + 1) - pow(I3[w].first, 2 * r + p + m + 1)) / (2 * r + m + p + 1) ;
          }
          A3 += sum * pow(pol1[0], r) / (factorial<Type>(r) * (n + i + 1) * factorial<Type>(i)) ;
        }
      }
    }
  }


  //Did not change anything in here
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
Type find_trig_area_2intersection_formula_second(const unsigned &m, const unsigned &n, const int &s, const Type &a, Type c, const int &table,  PointT <Type> &p1,  PointT <Type> &p2, const PointT <Type> &p3) {
  Type area(0);
  Type A1(0), A2(0), A3(0);
  std::vector< std::pair <Type, Type> > I1, I2, I3 ;

  std::vector <Type> pol1(3, 0);
  std::vector <Type> pol2(3, 0);
  std::vector <Type> intersection;
  std::vector <Type> interp_point;
  Parabola <Type> parabola;
  Type ankor(0) ;
  parabola = get_parabola_equation(p1, p2, p3);

  Type k = parabola.k;
  Type b = parabola.b;
  Type d = parabola.d;
  Type singleintersection;

//     cout << "\n---------------------- \n points = \n("<<p1.x<<","<<p1.y<<")\n"<<"("<<p2.x<<","<<p2.y<<")\n"<<"("<<p3.x<<","<<p3.y<<")\n"<<endl;
//     cout<< "parabola = "<<k<<"x^2 +"<<b<<"x+"<<d<<"+y=0"<<endl;

  bool do_line = 0;

//     if(fabs(k - 0.) < 0.00000000001){
  if(k != 0) {
    if(table == 1) { //we only use modified integrals if it is concave down
      if(k > 0) {
        ankor = p1.x;
        do_line = 1;
        Type delta = (b + 1.) * (b + 1.) - 4 * k * d;
        //         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
        if(delta >= 0) {
          Type sqrtdelta = sqrt(delta);
          int sign = (k > 0) ? 1 : -1;
          for(unsigned i = 0; i < 2; i++) {
            Type x = (- (b + 1) - sign * sqrtdelta) / (2 * k);
            //             cout<< "Top x = "<< x<< endl;
            if(x > 0 && x < 1 && x > p1.x) {
              ankor = x;

            }
            sign *= -1;
          }
        }
      }
      //                           cout<<"\nankor" << ankor <<endl;

    }

    else if(table == 2) { // There are six possible cases we have to use modified integrals
      do_line = 1;
      if(k >= 0) { //concave down (2 possible scenerio)
        ankor = p1.x;
        if(p1.x < p2.x) { //case (a) take highest p1.x
          Type delta = (b + 1.) * (b + 1.) - 4 * k * d;
          if(delta >= 0) {
            Type sqrtdelta = sqrt(delta);
            int sign = 1;    //TODO we can get rid of this if
            for(unsigned i = 0; i < 2; i++) {
              Type x = (- (b + 1) - sign * sqrtdelta) / (2 * k);
              //          cout<< "Top x = "<< x<< endl;
              if(x > 0 && x < 1 && x > p1.x) {
                ankor = x;
              }
              sign *= -1;
            }
          }
        }
        else { //case b and c  take lowest p1.x
          Type delta = (b + 1.) * (b + 1.) - 4 * k * d;
          if(delta >= 0) {
            Type sqrtdelta = sqrt(delta);
            int sign = (k > 0) ? -1 : 1;  //this gives us highest p1.x first then lowest.
            for(unsigned i = 0; i < 2; i++) {
              Type x = (- (b + 1) - sign * sqrtdelta) / (2 * k);
              //          cout<< "Top x = "<< x<< endl;
              if(x > 0 && x < 1 && x < p1.x) {
                ankor = x;
              }
              sign *= -1;
            }
          }
        }

      }
      else {    //concave up k<0 , table 2 case d,e,f
        ankor = p2.x ;
        if(p1.x < p2.x) { //case (d and e) take lowest
          Type delta = b * b - 4 * k * d;
          if(delta >= 0) {
            Type sqrtdelta = sqrt(delta);
            int sign = (k > 0) ? 1 : -1;    //TODO we can get rid of this if
            for(unsigned i = 0; i < 2; i++) {
              Type x = (- b - sign * sqrtdelta) / (2 * k);
              //          cout<< "Top x = "<< x<< endl;
              if(x > 0 && x < 1 && x < p2.x) {
                ankor = x;
              }
              sign *= -1;
            }
          }
        }
        else { //case f,  x1>x2
          Type delta = (b) * (b) - 4 * k * d;
          if(delta >= 0) {
            Type sqrtdelta = sqrt(delta);
            int sign = (k > 0) ? -1 : 1;  //this gives us highest x first then lowest.
            for(unsigned i = 0; i < 2; i++) {
              Type x = (- (b) - sign * sqrtdelta) / (2 * k);
              //          cout<< "Top x = "<< x<< endl;
              if(x > 0 && x < 1 && x > p2.x) { //highest p2.x
                ankor = x;
              }
              sign *= -1;
            }
          }
        }
      }
    }
    else if(table == 3) {
      if(k < 0) {
        ankor = p2.x;
        do_line = 1;
        Type delta = b * b - 4 * k * d;
        //         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
        if(delta >= 0) {
          Type sqrtdelta = sqrt(delta);
          int sign = (k > 0) ? 1 : -1;
          for(unsigned i = 0; i < 2; i++) {
            Type x = (- b - sign * sqrtdelta) / (2 * k);
            //             cout<< "Top x = "<< x<< endl;
            if(x < 1 && x > 0 && x > p2.x) {    //highest p2.x  //TODO should use if ( x>0 &&  x < 1 && x < p2.x )
              ankor = x;
            }
            sign *= -1;
          }
        }
      }
    }

  }
  pol2[0] = parabola.k;
  pol2[1] = parabola.b;
  pol2[2] = parabola.d;

  if(do_line) { // TODO we donot have to dot them separately. We can do this for each table above

    if(table == 1) {
      I1.resize(0);
      I1.resize(1, std::pair<Type, Type>(ankor, static_cast<Type>(1)));  //not sure if it is taking value. Lets do I1 manually.
      area = trig_integral_A3(m, n, s, a, c, pol2, {{ankor, static_cast<Type>(1)}}) -  trig_integral_A2(m, n, s, a, c, pol2, {{ankor, static_cast<Type>(1)}});
    }
    else if(table == 2) { //TODO
      if(k > 0) {
        if(p1.x < p2.x) {   //a  //TODO check all the others. I have just fixed one problem in case a. It was (p1.x,1) . It should be (p2.x,1) . Check other cases.
          I1.resize(0);
          I1.resize(1, std::pair<Type, Type>(static_cast<Type>(ankor), static_cast<Type>(p2.x)));
          I3.resize(0);
          I3.resize(1, std::pair<Type, Type>(static_cast<Type>(p2.x), static_cast<Type>(1)));  //not sure if it is taking value. Lets do I1 manually.
          area = trig_integral_A3(m, n, s, a, c, pol2, {{ankor, p2.x}}) -  trig_integral_A2(m, n, s, a, c, pol2, {{ankor, p2.x}}) + trig_integral_A3(m, n, s, a, c, pol2, {{p2.x, static_cast<Type>(1)}});
        }
        else {  //b,c
//               I1.resize(0);
//               I1.resize(1, std::pair<Type, Type>(static_cast<Type>(p2.x), static_cast<Type>(p1.x)));
//               I3.resize(0);
//               I3.resize(1, std::pair<Type, Type>(static_cast<Type>(0), static_cast<Type>(p2.x)));
//               area = trig_integral_A3(m, n, s, a, c, pol2, {{p2.x, p1.x}}) -  trig_integral_A2(m, n, s, a, c, pol2, {{p2.x, p1.x}}) + trig_integral_A3(m, n, s, a, c, pol2, {{ static_cast<Type>(0), p2.x}});

          I2.resize(0);  //Here we are integrating oposite site. normal -1 .
          I2.resize(1, std::pair<Type, Type>(static_cast<Type>(p2.x), static_cast<Type>(ankor)));
          I3.resize(0);
          I3.resize(1, std::pair<Type, Type>(static_cast<Type>(ankor), static_cast<Type>(1)));
//               area = trig_integral_A2(m, n, s, a, c, pol2, {{p2.x, p1.x}}) + trig_integral_A3(m, n, s, a, c, pol2, {{ p1.x, static_cast<Type>(1)}});
          area = trig_integral_A2(m, n, s, a, c, pol2, I2) + trig_integral_A3(m, n, s, a, c, pol2, I3);
        }
      }
      else {
        if(p1.x < p2.x) { //d,e
          I1.resize(0);
          I1.resize(1, std::pair<Type, Type>(static_cast<Type>(p1.x), static_cast<Type>(ankor)));
          I3.resize(0);
          I3.resize(1, std::pair<Type, Type>(static_cast<Type>(ankor), static_cast<Type>(1)));
          area = trig_integral_A3(m, n, s, a, c, pol2, {{p1.x, ankor}}) -  trig_integral_A2(m, n, s, a, c, pol2, {{p1.x, ankor}}) + trig_integral_A3(m, n, s, a, c, pol2, {{ ankor, static_cast<Type>(1)}});
        }
        else { //f //TODO check
//               I1.resize(0);
//               I1.resize(1, std::pair<Type, Type>(static_cast<Type>(p2.x), static_cast<Type>(p1.x)));
//               I3.resize(0);
//               I3.resize(1, std::pair<Type, Type>(static_cast<Type>(0), static_cast<Type>(p2.x)));
//               area = trig_integral_A3(m, n, s, a, c, pol2, {{p2.x, p1.x}}) -  trig_integral_A2(m, n, s, a, c, pol2, {{p2.x, p1.x}}) + trig_integral_A3(m, n, s, a, c, pol2, {{ static_cast<Type>(0), p2.x}});

          I2.resize(0);  //Here we are integrating oposite site. normal -1 .
          I2.resize(1, std::pair<Type, Type>(static_cast<Type>(ankor), static_cast<Type>(p1.x)));
          I3.resize(0);
          I3.resize(1, std::pair<Type, Type>(static_cast<Type>(p1.x), static_cast<Type>(1)));
//               area = trig_integral_A2(m, n, s, a, c, pol2, {{p2.x, p1.x}}) + trig_integral_A3(m, n, s, a, c, pol2, {{ p1.x, static_cast<Type>(1)}});
          area = trig_integral_A2(m, n, s, a, c, pol2, I2) + trig_integral_A3(m, n, s, a, c, pol2, I3);


        }
      }

    }

    else if(table == 3) {
      I1.resize(0);
      I3.resize(0);
      I1.resize(1, std::pair<Type, Type>(static_cast<Type>(ankor), static_cast<Type>(1)));
      I3.resize(1, std::pair<Type, Type>(static_cast<Type>(0), static_cast<Type>(ankor)));
      area = trig_integral_A3(m, n, s, a, c, pol2, I1) -  trig_integral_A2(m, n, s, a, c, pol2, I1) + trig_integral_A3(m, n, s, a, c, pol2, I3);
    }
  }

  else {

    pol1[0] = k + a;
    pol1[1] = b + c;
    pol1[2] = d;
    GetIntervalall<Type, double>(pol1, pol2, I1, I2, I3);

    if(I1.size() > 0) {
      A1 = trig_integral_A3(m, n, s, a, c, pol2, I1) -  trig_integral_A2(m, n, s, a, c, pol2, I1);
    }
    if(I2.size() > 0) {
      A2 = trig_integral_A2(m, n, s, a, c, pol2, I2);
    }
    if(I3.size() > 0) {
      A3 = trig_integral_A3(m, n, s, a, c, pol2, I3);
    }
    area = A1 + A2 + A3;
  }

  //this is just calculating area forcefuly

//         pol1[0] = k+a; pol1[1] = b + c; pol1[2] = d;
//         GetIntervalall<Type, double>(pol1, pol2, I1, I2, I3);
//
//         if(I1.size() > 0) {
//             A1 = trig_integral_A3(m, n, s, a, c, pol2, I1) -  trig_integral_A2(m, n, s, a, c, pol2, I1);
//         }
//         if(I2.size() > 0) {
//             A2 = trig_integral_A2(m, n, s, a, c, pol2, I2);
//         }
//         if(I3.size() > 0) {
//             A3 = trig_integral_A3(m, n, s, a, c, pol2, I3);
//         }
//         area = A1 + A2 + A3;

  return area ;
}

template <class Type>
void change_points_covert_table(const PointT <Type> &p1, const PointT <Type> &p2, const PointT <Type> &p3, const int &actual_table, int &old_table, PointT <Type> &q1, PointT <Type> &q2, PointT <Type> &q3, bool &vertical) {

  if(actual_table < 3) {
    vertical = true ;
    q1 = {(1. - p1.x), p1.y};
    q2 = {(1. - p2.x), p2.y};
    q3 = {(1. - p3.x), p3.y};

    if(actual_table == 0) { // swap q1 and q2
      old_table = 1 ;
      q1 = {(1. - p2.x), p2.y};
      q2 = {(1. - p1.x), p1.y};
    }
    else if(actual_table == 1) old_table = 3 ;
    else if(actual_table == 2) old_table = 2 ;
  }
  else {
    vertical = false ;
    q1 = {(1. - p1.y), p1.x};
    q2 = {(1. - p2.y), p2.x};
    q3 = {(1. - p3.y), p3.x};

    if(actual_table == 3) old_table = 2 ;
    else if(actual_table == 4) old_table = 3 ;
    else if(actual_table == 5) {
      old_table = 1 ;
      q1 = {(1. - p2.y), p2.x};
      q2 = {(1. - p1.y), p1.x};
    }
  }
}


template <class Type>
Type find_trig_area_2intersection_formula_first(const unsigned &m, const unsigned &n, const int &s, const Type &a, Type c, const int &actual_table,  PointT <Type> &p1,  PointT <Type> &p2, const PointT <Type> &p3) {
  Type area(0);
  PointT <Type> q1, q2, q3 ;
  int old_table;
  bool vertical;

  change_points_covert_table<Type>(p1, p2, p3, actual_table, old_table, q1, q2, q3, vertical);
  area = find_trig_area_2intersection_formula_second<Type>(m, n, s, a, c, old_table,  q1,  q2, q3);

  return area;
}


struct Point3D {
  double x, y, z;

  // Default constructor
  Point3D() : x(0), y(0), z(0) {}

  // Constructor with double values
  Point3D(double x, double y, double z) : x(x), y(y), z(z) {}

  // Constructor with pointers (as in your original definition)
  Point3D(double* x, double* y, double* z) : x(*x), y(*y), z(*z) {}
};

template <class Type>
Type trilinier_interpolation(std::vector< std::vector< Type >> & interp_table, const std::vector< Type > &interp_point) {

  Type x = interp_point[0];
  Type y = interp_point[1];
  Type z = interp_point[2];

  Type c_000 = interp_table[0][3];
  Type c_001 = interp_table[1][3];
  Type c_010 = interp_table[2][3];
  Type c_011 = interp_table[3][3];
  Type c_100 = interp_table[4][3];
  Type c_101 = interp_table[5][3];
  Type c_110 = interp_table[6][3];
  Type c_111 = interp_table[7][3];

  Type x0 = interp_table[0][0];
  Type x1 = interp_table[7][0];
  Type y0 = interp_table[0][1];
  Type y1 = interp_table[7][1];
  Type z0 = interp_table[0][2];
  Type z1 = interp_table[7][2];

  Type x_d = (x - x0) / (x1 - x0);
  Type y_d = (y - y0) / (y1 - y0);
  Type z_d = (z - z0) / (z1 - z0);

  Type c_00 = c_000 * (1 - x_d) + c_100 * x_d ;
  Type c_01 = c_001 * (1 - x_d) + c_101 * x_d ;
  Type c_10 = c_010 * (1 - x_d) + c_110 * x_d ;
  Type c_11 = c_011 * (1 - x_d) + c_111 * x_d ;

  Type c_0 = c_00 * (1 - y_d) + c_10 * y_d ;
  Type c_1 = c_01 * (1 - y_d) + c_11 * y_d ;

  Type cc = c_0 * (1 - z_d) + c_1 * z_d ;

  return cc;
}

template <class Type>
void trilinier_interpolation_vector(const std::vector< std::vector< Type >> & interp_table, const std::vector< std::vector< Type >> & interp_table_values, const std::vector< Type > &interp_point, std::vector< Type > &interp_point_values) {

  interp_point_values.resize(interp_table_values[0].size());
  Type x = interp_point[0];
  Type y = interp_point[1];
  Type z = interp_point[2];

  Type x0 = interp_table[0][0];
  Type x1 = interp_table[7][0];
  Type y0 = interp_table[0][1];
  Type y1 = interp_table[7][1];
  Type z0 = interp_table[0][2];
  Type z1 = interp_table[7][2];

  Type x_d = (x - x0) / (x1 - x0);
  Type y_d = (y - y0) / (y1 - y0);
  Type z_d = (z - z0) / (z1 - z0);

  for(unsigned i = 0; i < interp_table_values[0].size(); i++) {
    Type c_000 = interp_table_values[0][i];
    Type c_001 = interp_table_values[1][i];
    Type c_010 = interp_table_values[2][i];
    Type c_011 = interp_table_values[3][i];
    Type c_100 = interp_table_values[4][i];
    Type c_101 = interp_table_values[5][i];
    Type c_110 = interp_table_values[6][i];
    Type c_111 = interp_table_values[7][i];

    Type c_00 = c_000 * (1 - x_d) + c_100 * x_d ;
    Type c_01 = c_001 * (1 - x_d) + c_101 * x_d ;
    Type c_10 = c_010 * (1 - x_d) + c_110 * x_d ;
    Type c_11 = c_011 * (1 - x_d) + c_111 * x_d ;

    Type c_0 = c_00 * (1 - y_d) + c_10 * y_d ;
    Type c_1 = c_01 * (1 - y_d) + c_11 * y_d ;

    interp_point_values[i] = c_0 * (1 - z_d) + c_1 * z_d ;
  }
}

template <class Type>
void trilinier_interpolation_vector_fixed(const std::vector<std::vector<Type>> & interp_table,
                                          const std::vector<std::vector<Type>> & interp_table_values,
                                          const std::vector<Type> &interp_point,
                                          std::vector<Type> &interp_point_values) {

  interp_point_values.resize(interp_table_values[0].size());

  // Find the minimum and maximum coordinates of the cube
//     Type x_min = interp_table[0][0];
//     Type x_max = interp_table[0][0];
//     Type y_min = interp_table[0][1];
//     Type y_max = interp_table[0][1];
//     Type z_min = interp_table[0][2];
//     Type z_max = interp_table[0][2];
//
//     for (size_t i = 1; i < interp_table.size(); ++i) {
//         x_min = std::min(x_min, interp_table[i][0]);
//         x_max = std::max(x_max, interp_table[i][0]);
//         y_min = std::min(y_min, interp_table[i][1]);
//         y_max = std::max(y_max, interp_table[i][1]);
//         z_min = std::min(z_min, interp_table[i][2]);
//         z_max = std::max(z_max, interp_table[i][2]);
//     }

  Type x_min = interp_table[0][0];
  Type x_max = interp_table[1][0];
  Type y_min = interp_table[0][1];
  Type y_max = interp_table[2][1];
  Type z_min = interp_table[0][2];
  Type z_max = interp_table[4][2];


  std::cout << "Cube bounds: x[" << x_min << ", " << x_max << "] y[" << y_min << ", " << y_max
            << "] z[" << z_min << ", " << z_max << "]" << std::endl;





  // Check if the dimensions are valid (non-zero)
  Type x_dim = x_max - x_min;
  Type y_dim = y_max - y_min;
  Type z_dim = z_max - z_min;

  if(x_dim < std::numeric_limits<Type>::epsilon() ||
      y_dim < std::numeric_limits<Type>::epsilon() ||
      z_dim < std::numeric_limits<Type>::epsilon()) {
    std::cout << "WARNING: Cube has zero or near-zero dimension!" << std::endl;
    // Use small epsilon values to avoid division by zero
    x_dim = (x_dim < std::numeric_limits<Type>::epsilon()) ?
            std::numeric_limits<Type>::epsilon() * 10 : x_dim;
    y_dim = (y_dim < std::numeric_limits<Type>::epsilon()) ?
            std::numeric_limits<Type>::epsilon() * 10 : y_dim;
    z_dim = (z_dim < std::numeric_limits<Type>::epsilon()) ?
            std::numeric_limits<Type>::epsilon() * 10 : z_dim;
  }

  // Calculate normalized coordinates of the query point
  Type x_d = (interp_point[0] - x_min) / x_dim;
  Type y_d = (interp_point[1] - y_min) / y_dim;
  Type z_d = (interp_point[2] - z_min) / z_dim;

  // Optional: Clamp coordinates to [0,1] if slightly outside the bounds
  x_d = std::max(static_cast<Type>(0), std::min(static_cast<Type>(1), x_d));
  y_d = std::max(static_cast<Type>(0), std::min(static_cast<Type>(1), y_d));
  z_d = std::max(static_cast<Type>(0), std::min(static_cast<Type>(1), z_d));

  std::cout << "Normalized coordinates: (" << x_d << ", " << y_d << ", " << z_d << ")" << std::endl;

  // Map vertices to their correct positions in the unit cube
  // We need to identify which vertex corresponds to which corner

  // Define corners of the unit cube
  std::vector<std::vector<int>> corners = {
    {0, 0, 0}, // lower-left-front  (000)
    {1, 0, 0}, // lower-right-front (100)
    {0, 1, 0}, // lower-left-back   (010)
    {1, 1, 0}, // lower-right-back  (110)
    {0, 0, 1}, // upper-left-front  (001)
    {1, 0, 1}, // upper-right-front (101)
    {0, 1, 1}, // upper-left-back   (011)
    {1, 1, 1}  // upper-right-back  (111)
  };

  // Map vertices to corners
  std::vector<int> vertex_map(8, -1);
  for(size_t i = 0; i < interp_table.size(); ++i) {
    // Determine which corner this vertex is closest to
    int closest_corner = 0;
    Type min_distance = std::numeric_limits<Type>::max();

    for(int c = 0; c < 8; ++c) {
      // Calculate the normalized position of this vertex
      Type nx = (interp_table[i][0] - x_min) / x_dim;
      Type ny = (interp_table[i][1] - y_min) / y_dim;
      Type nz = (interp_table[i][2] - z_min) / z_dim;

      // Calculate distance to this corner
      Type dx = nx - corners[c][0];
      Type dy = ny - corners[c][1];
      Type dz = nz - corners[c][2];
      Type distance = dx * dx + dy * dy + dz * dz;

      if(distance < min_distance) {
        min_distance = distance;
        closest_corner = c;
      }
    }

    // Record this vertex as being closest to this corner
    vertex_map[closest_corner] = i;
  }

  // Check that all corners have a vertex assigned
  bool all_corners_mapped = true;
  for(int i = 0; i < 8; ++i) {
    if(vertex_map[i] == -1) {
      all_corners_mapped = false;
      std::cout << "WARNING: Corner " << i << " doesn't have a vertex assigned!" << std::endl;

      // As a fallback, find the first unused vertex
      for(size_t j = 0; j < interp_table.size(); ++j) {
        bool is_used = false;
        for(int k = 0; k < 8; ++k) {
          if(vertex_map[k] == static_cast<int>(j)) {
            is_used = true;
            break;
          }
        }

        if(!is_used) {
          vertex_map[i] = j;
          break;
        }
      }
    }
  }

  std::cout << "Vertex mapping: ";
  for(int i = 0; i < 8; ++i) {
    std::cout << vertex_map[i] << " ";
  }
  std::cout << std::endl;

  // Now perform the trilinear interpolation using the mapped vertices
  for(unsigned i = 0; i < interp_table_values[0].size(); i++) {
    // Get the values at the 8 corners of the cube using the vertex mapping
    Type c_000 = interp_table_values[vertex_map[0]][i]; // (0,0,0)
    Type c_100 = interp_table_values[vertex_map[1]][i]; // (1,0,0)
    Type c_010 = interp_table_values[vertex_map[2]][i]; // (0,1,0)
    Type c_110 = interp_table_values[vertex_map[3]][i]; // (1,1,0)
    Type c_001 = interp_table_values[vertex_map[4]][i]; // (0,0,1)
    Type c_101 = interp_table_values[vertex_map[5]][i]; // (1,0,1)
    Type c_011 = interp_table_values[vertex_map[6]][i]; // (0,1,1)
    Type c_111 = interp_table_values[vertex_map[7]][i]; // (1,1,1)

    // Perform trilinear interpolation
    Type c_00 = c_000 * (1 - x_d) + c_100 * x_d;
    Type c_01 = c_001 * (1 - x_d) + c_101 * x_d;
    Type c_10 = c_010 * (1 - x_d) + c_110 * x_d;
    Type c_11 = c_011 * (1 - x_d) + c_111 * x_d;

    Type c_0 = c_00 * (1 - y_d) + c_10 * y_d;
    Type c_1 = c_01 * (1 - y_d) + c_11 * y_d;

    interp_point_values[i] = c_0 * (1 - z_d) + c_1 * z_d;
  }
}

template <class Type>
Type trilinier_interpolation_FEM_orientation(const int table, std::vector< std::vector< Type >> & interp_table, const std::vector< Type > &interp_point) {

  Type x = interp_point[0];
  Type y = interp_point[1];

  Type denominator(0), denominatorz0(0), denominatorz1(0);

  if(table < 2 || table > 3) denominator = y - 2;
  else denominator = x + y - 2;

  Type z = -(2 * interp_point[2]) / denominator;


  Type c_000 = interp_table[0][3];
  Type c_100 = interp_table[1][3];
  Type c_110 = interp_table[2][3];
  Type c_010 = interp_table[3][3];
  Type c_001 = interp_table[4][3];
  Type c_101 = interp_table[5][3];
  Type c_111 = interp_table[6][3];
  Type c_011 = interp_table[7][3];

  Type x0 = interp_table[0][0];
  Type x1 = interp_table[1][0];
  Type y0 = interp_table[0][1];
  Type y1 = interp_table[7][1];

  if(table < 2 || table > 3) {
    denominatorz0 = interp_table[0][1] - 2.;
    denominatorz1 = interp_table[7][1] - 2.;
  }
  else {
    denominatorz0  = interp_table[0][0] + interp_table[0][1] - 2.;
    denominatorz1 = interp_table[7][0] + interp_table[7][1] - 2.;
  }
  Type z0 = - (2 * interp_table[0][2]) / denominatorz0;
  Type z1 = - (2 * interp_table[7][2]) / denominatorz1;

  Type x_d = (x - x0) / (x1 - x0);
  Type y_d = (y - y0) / (y1 - y0);
  Type z_d = (z - z0) / (z1 - z0);

  Type c_00 = c_000 * (1 - x_d) + c_100 * x_d ;
  Type c_01 = c_001 * (1 - x_d) + c_101 * x_d ;
  Type c_10 = c_010 * (1 - x_d) + c_110 * x_d ;
  Type c_11 = c_011 * (1 - x_d) + c_111 * x_d ;

  Type c_0 = c_00 * (1 - y_d) + c_10 * y_d ;
  Type c_1 = c_01 * (1 - y_d) + c_11 * y_d ;

  Type cc = c_0 * (1 - z_d) + c_1 * z_d ;

  return cc;
}


template <class Type>
void trilinier_interpolation_vector_FEM_orientation(const std::vector< std::vector< Type >> & interp_table, const std::vector< std::vector< Type >> & interp_table_values, const std::vector< Type > &interp_point, std::vector< Type > &interp_point_values) {

  interp_point_values.resize(interp_table_values[0].size());
  Type x = interp_point[0];
  Type y = interp_point[1];
  Type z = interp_point[2];

  Type x0 = interp_table[0][0];
  Type x1 = interp_table[1][0];
  Type y0 = interp_table[0][1];
  Type y1 = interp_table[7][1];
  Type z0 = interp_table[0][2];
  Type z1 = interp_table[7][2];

  Type x_d = (x - x0) / (x1 - x0);
  Type y_d = (y - y0) / (y1 - y0);
  Type z_d = (z - z0) / (z1 - z0);

  for(unsigned i = 0; i < interp_table_values[0].size(); i++) {

    Type c_000 = interp_table_values[0][i];
    Type c_100 = interp_table_values[1][i];
    Type c_110 = interp_table_values[2][i];
    Type c_010 = interp_table_values[3][i];
    Type c_001 = interp_table_values[4][i];
    Type c_101 = interp_table_values[5][i];
    Type c_111 = interp_table_values[6][i];
    Type c_011 = interp_table_values[7][i];

    Type c_00 = c_000 * (1 - x_d) + c_100 * x_d ;
    Type c_01 = c_001 * (1 - x_d) + c_101 * x_d ;
    Type c_10 = c_010 * (1 - x_d) + c_110 * x_d ;
    Type c_11 = c_011 * (1 - x_d) + c_111 * x_d ;

    Type c_0 = c_00 * (1 - y_d) + c_10 * y_d ;
    Type c_1 = c_01 * (1 - y_d) + c_11 * y_d ;

    interp_point_values[i] = c_0 * (1 - z_d) + c_1 * z_d ;
  }
}


template <class Type>
void trilinear_interpolation_vector_deformed(
  const std::vector<std::vector<Type>>& cube_vertices,       // 8 vertices of the deformed cube
  const std::vector<std::vector<Type>>& vertex_values,       // Values at each vertex
  const std::vector<Type>& query_point,                      // Point to interpolate at
  std::vector<Type>& interpolated_values) {                  // Result storage

  interpolated_values.resize(vertex_values[0].size());

  // We need to compute the trilinear coordinates (weights) for the query point
  // First, we'll map the query point to a reference unit cube coordinate system
  // using an inverse mapping approach

  // For the new corner ordering:
  // 0: (0,0,0) - Bottom-front-left
  // 1: (1,0,0) - Bottom-front-right
  // 2: (1,1,0) - Bottom-back-right
  // 3: (0,1,0) - Bottom-back-left
  // 4: (0,0,1) - Top-front-left
  // 5: (1,0,1) - Top-front-right
  // 6: (1,1,1) - Top-back-right
  // 7: (0,1,1) - Top-back-left

  Type u = 0.5, v = 0.5, w = 0.5;  // Initial guess at center of unit cube
  const int MAX_ITERATIONS = 20;
  const Type TOLERANCE = Type(1e-8);

  // Newton-Raphson iterations to find u,v,w coordinates
  for(int iter = 0; iter < MAX_ITERATIONS; ++iter) {
    // Compute the interpolated position at current u,v,w
    std::vector<Type> interpolated_position(3, 0);

    // Basis functions for the corners of a unit cube with counter-clockwise ordering
    Type N[8];
    N[0] = (1 - u) * (1 - v) * (1 - w); // (0,0,0) - Bottom-front-left
    N[1] = u     * (1 - v) * (1 - w); // (1,0,0) - Bottom-front-right
    N[2] = u     * v     * (1 - w); // (1,1,0) - Bottom-back-right
    N[3] = (1 - u) * v     * (1 - w); // (0,1,0) - Bottom-back-left
    N[4] = (1 - u) * (1 - v) * w;  // (0,0,1) - Top-front-left
    N[5] = u     * (1 - v) * w;    // (1,0,1) - Top-front-right
    N[6] = u     * v     * w;      // (1,1,1) - Top-back-right
    N[7] = (1 - u) * v     * w;    // (0,1,1) - Top-back-left

    // P = sum(N_i * P_i)
    for(int i = 0; i < 8; ++i) {
      for(int j = 0; j < 3; ++j) {
        interpolated_position[j] += N[i] * cube_vertices[i][j];
      }
    }

    // Compute the residual (difference between actual and current interpolated position)
    std::vector<Type> residual(3);
    for(int j = 0; j < 3; ++j) {
      residual[j] = query_point[j] - interpolated_position[j];
    }

    // Early exit if we're close enough
    if(std::abs(residual[0]) < TOLERANCE &&
        std::abs(residual[1]) < TOLERANCE &&
        std::abs(residual[2]) < TOLERANCE) {
      break;
    }

    // Compute Jacobian matrix
    // J_ij = dP_i/dj where j is u,v,w and P_i is x,y,z
    std::vector<std::vector<Type>> jacobian(3, std::vector<Type>(3, 0));

    // Derivatives of basis functions with respect to u
    Type dNdu[8];
    dNdu[0] = -(1 - v) * (1 - w);  // d/du[(1-u)(1-v)(1-w)]
    dNdu[1] = (1 - v) * (1 - w);   // d/du[u(1-v)(1-w)]
    dNdu[2] = v * (1 - w);         // d/du[u*v*(1-w)]
    dNdu[3] = -v * (1 - w);        // d/du[(1-u)*v*(1-w)]
    dNdu[4] = -(1 - v) * w;        // d/du[(1-u)(1-v)w]
    dNdu[5] = (1 - v) * w;         // d/du[u(1-v)w]
    dNdu[6] = v * w;               // d/du[u*v*w]
    dNdu[7] = -v * w;              // d/du[(1-u)*v*w]

    // Derivatives of basis functions with respect to v
    Type dNdv[8];
    dNdv[0] = -(1 - u) * (1 - w);  // d/dv[(1-u)(1-v)(1-w)]
    dNdv[1] = -u * (1 - w);        // d/dv[u(1-v)(1-w)]
    dNdv[2] = u * (1 - w);         // d/dv[u*v*(1-w)]
    dNdv[3] = (1 - u) * (1 - w);   // d/dv[(1-u)*v*(1-w)]
    dNdv[4] = -(1 - u) * w;        // d/dv[(1-u)(1-v)w]
    dNdv[5] = -u * w;              // d/dv[u(1-v)w]
    dNdv[6] = u * w;               // d/dv[u*v*w]
    dNdv[7] = (1 - u) * w;         // d/dv[(1-u)*v*w]

    // Derivatives of basis functions with respect to w
    Type dNdw[8];
    dNdw[0] = -(1 - u) * (1 - v);  // d/dw[(1-u)(1-v)(1-w)]
    dNdw[1] = -u * (1 - v);        // d/dw[u(1-v)(1-w)]
    dNdw[2] = -u * v;              // d/dw[u*v*(1-w)]
    dNdw[3] = -(1 - u) * v;        // d/dw[(1-u)*v*(1-w)]
    dNdw[4] = (1 - u) * (1 - v);   // d/dw[(1-u)(1-v)w]
    dNdw[5] = u * (1 - v);         // d/dw[u(1-v)w]
    dNdw[6] = u * v;               // d/dw[u*v*w]
    dNdw[7] = (1 - u) * v;         // d/dw[(1-u)*v*w]

    // Fill the Jacobian
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 8; ++j) {
        jacobian[i][0] += dNdu[j] * cube_vertices[j][i];
        jacobian[i][1] += dNdv[j] * cube_vertices[j][i];
        jacobian[i][2] += dNdw[j] * cube_vertices[j][i];
      }
    }

    // Solve the system J * delta = residual for delta using Cramer's rule
    // (for a 3x3 system this is faster than general matrix inversion)
    Type det = jacobian[0][0] * (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1])
               - jacobian[0][1] * (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0])
               + jacobian[0][2] * (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]);

    if(std::abs(det) < TOLERANCE) {
      // Jacobian is close to singular, we might need a better initial guess
      // or the point might be outside the cube
      break;
    }

    Type delta_u = (residual[0] * (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1])
                    - jacobian[0][1] * (residual[1] * jacobian[2][2] - jacobian[1][2] * residual[2])
                    + jacobian[0][2] * (residual[1] * jacobian[2][1] - jacobian[1][1] * residual[2])) / det;

    Type delta_v = (jacobian[0][0] * (residual[1] * jacobian[2][2] - jacobian[1][2] * residual[2])
                    - residual[0] * (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0])
                    + jacobian[0][2] * (jacobian[1][0] * residual[2] - residual[1] * jacobian[2][0])) / det;

    Type delta_w = (jacobian[0][0] * (jacobian[1][1] * residual[2] - residual[1] * jacobian[2][1])
                    - jacobian[0][1] * (jacobian[1][0] * residual[2] - residual[1] * jacobian[2][0])
                    + residual[0] * (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0])) / det;

    // Update u, v, w
    u += delta_u;
    v += delta_v;
    w += delta_w;

    // Clamp to unit cube
    u = std::max(Type(0), std::min(Type(1), u));
    v = std::max(Type(0), std::min(Type(1), v));
    w = std::max(Type(0), std::min(Type(1), w));
  }

  // Now that we have u, v, w, compute the final interpolated values
  // Recompute basis functions for the final u, v, w based on new corner ordering
  Type N[8];
  N[0] = (1 - u) * (1 - v) * (1 - w); // (0,0,0) - Bottom-front-left
  N[1] = u     * (1 - v) * (1 - w); // (1,0,0) - Bottom-front-right
  N[2] = u     * v     * (1 - w); // (1,1,0) - Bottom-back-right
  N[3] = (1 - u) * v     * (1 - w); // (0,1,0) - Bottom-back-left
  N[4] = (1 - u) * (1 - v) * w;  // (0,0,1) - Top-front-left
  N[5] = u     * (1 - v) * w;    // (1,0,1) - Top-front-right
  N[6] = u     * v     * w;      // (1,1,1) - Top-back-right
  N[7] = (1 - u) * v     * w;    // (0,1,1) - Top-back-left

  // Interpolate all value components
  for(size_t i = 0; i < interpolated_values.size(); ++i) {
    interpolated_values[i] = 0;
    for(int j = 0; j < 8; ++j) {
      interpolated_values[i] += N[j] * vertex_values[j][i];
    }
  }
}



template <class Type>
void trilinear_interpolation_vector_remap_to_unitcube(
  const int table,
  const std::vector<std::vector<Type>>& cube_vertices,       // 8 vertices of the deformed cube
  const std::vector<std::vector<Type>>& vertex_values,       // Values at each vertex
  const std::vector<Type>& query_point,                      // Point to interpolate at
  std::vector<Type>& interpolated_values) {

  std::cout << "DEBUG: Input parameters" << std::endl;
  std::cout << "Table: " << table << std::endl;
  std::cout << "Query point: (";
  for(size_t i = 0; i < query_point.size(); ++i) {
    std::cout << query_point[i];
    if(i < query_point.size() - 1) std::cout << ", ";
  }
  std::cout << ")" << std::endl;

  // Create new cube vertices with the modified mapping
  std::vector<std::vector<Type>> new_cube_vertices(cube_vertices.size());

  // Create the transformed cube vertices with table-dependent mapping
//     std::cout << "\nTransforming cube vertices:" << std::endl;
  for(size_t i = 0; i < cube_vertices.size(); i++) {
    new_cube_vertices[i] = cube_vertices[i];  // Copy all coordinates

    // Extract coordinates for clearer code
    Type x = cube_vertices[i][0];
    Type y = cube_vertices[i][1];
    Type z = cube_vertices[i][2];
    Type denominator = 0;

    // Apply the mapping for z-coordinate based on table number
    if(table < 2 || table > 3) {
      denominator = y - 2;
    }
    else {
      denominator = x + y - 2;
    }

    // Apply the transformation if denominator is not close to zero
    if(std::abs(denominator) > std::numeric_limits<Type>::epsilon()) {
      Type new_z = -(2 * z) / denominator;
//             std::cout << "  Transformed z: " << z << " -> " << new_z << std::endl;
      new_cube_vertices[i][2] = new_z;
    }
    else {
      // Handle case where denominator is close to zero
      std::cout << "  WARNING: Denominator close to zero. Keeping original z value." << std::endl;
      new_cube_vertices[i][2] = z; // Keep original as a fallback
    }
  }

  // Print new cube vertices
  std::cout << "\nTransformed cube vertices:" << std::endl;
  for(size_t i = 0; i < new_cube_vertices.size(); ++i) {
    std::cout << "Vertex " << i << ": (";
    for(size_t j = 0; j < new_cube_vertices[i].size(); ++j) {
      std::cout << new_cube_vertices[i][j];
      if(j < new_cube_vertices[i].size() - 1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
  }

  // Create new query point with the same table-dependent transformation
  std::vector<Type> new_query_point = query_point;

  // Extract coordinates for clearer code
  Type x = query_point[0];
  Type y = query_point[1];
  Type z = query_point[2];
  Type denominator = 0;

//     std::cout << "\nTransforming query point:" << std::endl;
  // Apply the mapping for z-coordinate based on table number
  if(table < 2 || table > 3) {
    denominator = y - 2;
  }
  else {
    denominator = x + y - 2;
  }

  // Apply the transformation if denominator is not close to zero
  if(std::abs(denominator) > std::numeric_limits<Type>::epsilon()) {
    Type new_z = -(2 * z) / denominator;
//         std::cout << "  Transformed z: " << z << " -> " << new_z << std::endl;
    new_query_point[2] = new_z;
  }
  else {
    // Handle case where denominator is close to zero
    std::cout << "  WARNING: Denominator close to zero. Keeping original z value." << std::endl;
    new_query_point[2] = z; // Keep original as a fallback
  }

  std::cout << "\nTransformed query point: (";
  for(size_t i = 0; i < new_query_point.size(); ++i) {
    std::cout << new_query_point[i];
    if(i < new_query_point.size() - 1) std::cout << ", ";
  }
  std::cout << ")" << std::endl;

  std::cout << "\nCalling trilinear interpolation with transformed coordinates..." << std::endl;
  // Everything is in unit_cube. We can use the standard trilinear interpolation now
  trilinier_interpolation_vector_FEM_orientation(new_cube_vertices, vertex_values, new_query_point, interpolated_values);

}


template <class Type>
Type trilinear_interpolation_deformed(
  const std::vector<std::vector<double>>& cube_vertices_and_values,  // Vertices (x,y,z) and value
  const std::vector<double>& query_point) {                          // Point to interpolate at

  // Extract vertices and values
  std::vector<std::vector<Type>> cube_vertices(8, std::vector<Type>(3));
  std::vector<std::vector<Type>> vertex_values(8, std::vector<Type>(1));

  for(int i = 0; i < 8; ++i) {
    // First 3 values are coordinates
    cube_vertices[i][0] = cube_vertices_and_values[i][0];
    cube_vertices[i][1] = cube_vertices_and_values[i][1];
    cube_vertices[i][2] = cube_vertices_and_values[i][2];

    // Last value is the scalar value at this vertex
    vertex_values[i][0] = cube_vertices_and_values[i][3];
  }

  // Prepare query point and result
  std::vector<Type> interp_point(3);
  interp_point[0] = query_point[0];
  interp_point[1] = query_point[1];
  interp_point[2] = query_point[2];

  std::vector<Type> result;

  // Call the vector version
  trilinear_interpolation_vector_deformed(cube_vertices, vertex_values, interp_point, result);

  // Return the single scalar result
  return result[0];
}

template <class Type>   //TODO change this based on 6 table
void get_p1_p2_p3(const int &table, const std::vector<double> &corner, PointT <Type> &p1, PointT <Type> &p2, PointT <Type> &p3) {
  double epsilon = 0.000000000000001;
  Type i1_pm_eps(-1), i2_pm_eps(-1), i3_pm_eps(-1);

  // std::cout << "Corner " << i << ": (" << corner[0] << ", " << corner[1] << ", " << corner[2] << ") - Print Something\n";

  switch(table) {
  case 0:
    i1_pm_eps = static_cast<Type>(corner[0]);
    i2_pm_eps = static_cast<Type>(corner[1] + epsilon);
    i3_pm_eps = static_cast<Type>(corner[2] + epsilon);

    p1 = {static_cast<Type>(0), i1_pm_eps};
    p2 = {i2_pm_eps, static_cast<Type>(1) - i2_pm_eps};
    p3 = {(p1.x + p2.x) * 0.5, i3_pm_eps};
    break;
  case 1:
    i1_pm_eps = static_cast<Type>(corner[0]);
    i2_pm_eps = static_cast<Type>(corner[1] + epsilon);
    i3_pm_eps = static_cast<Type>(corner[2] + epsilon);

    p1 = {static_cast<Type>(0), i1_pm_eps};
    p2 = {i2_pm_eps, static_cast<Type>(0)};
    p3 = {(p1.x + p2.x) * 0.5, i3_pm_eps};
    break;
  case 2:
    //Do we really need epsilon on this table?
    i1_pm_eps = static_cast<Type>(corner[0]);
    i2_pm_eps = static_cast<Type>(corner[1] - epsilon);
    i3_pm_eps = static_cast<Type>(corner[2] + epsilon);

    p1 = {i1_pm_eps, static_cast<Type>(1) - i1_pm_eps};
    p2 = {i2_pm_eps, static_cast<Type>(0)};
    p3 = {(p1.x + p2.x) * 0.5, i3_pm_eps};
    break;
  case 3:
    i1_pm_eps = static_cast<Type>(corner[0]);
    i2_pm_eps = static_cast<Type>(corner[1] - epsilon);
    i3_pm_eps = static_cast<Type>(corner[2] + epsilon);

    p1 = {static_cast<Type>(1) - i1_pm_eps, i1_pm_eps};
    p2 = {static_cast<Type>(0), i2_pm_eps};
    p3 = {i3_pm_eps, (p1.y + p2.y) * 0.5};

//             i1_pm_eps = static_cast<Type>(corner[0] - epsilon);
//             i2_pm_eps = static_cast<Type>(corner[1] );
//             p1 = {static_cast<Type>(0), i1_pm_eps};
//             p2 = {i2_pm_eps, static_cast<Type>(1) - i2_pm_eps};
//             p3 = {static_cast<Type>(corner[2] + epsilon), (p1.y + p2.y)*0.5};
    break;
  case 4:
    i1_pm_eps = static_cast<Type>(corner[0]);
    i2_pm_eps = static_cast<Type>(corner[1] + epsilon);
    i3_pm_eps = static_cast<Type>(corner[2] + epsilon);

    p1 = {i1_pm_eps, static_cast<Type>(0)};
    p2 = {static_cast<Type>(0), i2_pm_eps};
    p3 = {i3_pm_eps, (p1.y + p2.y) * 0.5};

//             i1_pm_eps = static_cast<Type>(corner[0] + epsilon);
//             i2_pm_eps = static_cast<Type>(corner[1] );
//             p1 = {static_cast<Type>(0), i1_pm_eps};
//             p2 = {i2_pm_eps, static_cast<Type>(0)};
//             p3 = {static_cast<Type>(corner[2]), (p1.y + p2.y)*0.5};
    break;

  case 5:
    i1_pm_eps = static_cast<Type>(corner[0]);
    i2_pm_eps = static_cast<Type>(corner[1] + epsilon);
    i3_pm_eps = static_cast<Type>(corner[2] + epsilon);

    p1 = {i1_pm_eps, static_cast<Type>(0)};
    p2 = {static_cast<Type>(1) - i2_pm_eps, i2_pm_eps};
    p3 = {i3_pm_eps, (p1.y + p2.y) * 0.5};

//             i1_pm_eps = static_cast<Type>(corner[0] - epsilon);
//             i2_pm_eps = static_cast<Type>(corner[1] );
//             p1 = {i1_pm_eps, static_cast<Type>(1) - i1_pm_eps};
//             p2 = {i2_pm_eps, static_cast<Type>(0)};
//             p3 = {static_cast<Type>(corner[2]), (p1.y + p2.y)*0.5};
    break;
  }

}

template <class Type>
void find_actual_table_trig(const PointT <Type> &p1, const PointT <Type> &p2,  PointT <Type> &p3, int &actual_table, int &old_table, Point3D &searchP, PointT <Type> &q1, PointT <Type> &q2, PointT <Type> &q3, bool &vertical) {    //TODO re arrange q1 and q2 as this will create a problem in two intersection formula.
  double epsilon = 0.0000000000001;
  vertical = true;
  bool xSpan = false;
  bool ySpan = false;

  if((p1.x < p3.x && p3.x < p2.x) || (p2.x < p3.x && p3.x < p1.x)) xSpan = true ;
  if((p1.y < p3.y && p3.y < p2.y) || (p2.y < p3.y && p3.y < p1.y)) ySpan = true ;

  cout << "xspan = " << xSpan << " ySpan = " << ySpan << endl;

  if(xSpan) {
    if(ySpan) {
      double dx = min(fabs(static_cast<double>(p1.x - p3.x)), fabs(static_cast<double>(p1.x - p3.x)));
      double dy = min(fabs(static_cast<double>(p1.y - p3.y)), fabs(static_cast<double>(p1.y - p3.y)));
      if(dx >= dy) vertical = true;
      else vertical = false;
    }
    else {
      vertical = true;
    }
  }
  else {
    if(ySpan) vertical = false;
    else {
      std::cout << " The parabola formed by this three points is not a function. Use line cuts " << std::endl;

    }
  }

  if(vertical) {

    q1 = {(1. - p1.x), p1.y};
    q2 = {(1. - p2.x), p2.y};
    q3 = {(1. - p3.x), p3.y};

    cout << "q1 = (" << q1.x << "," << q1.y << ")" << endl;
    cout << "q2 = (" << q2.x << "," << q2.y << ")" << endl;
    cout << "q3 = (" << q3.x << "," << q3.y << ")" << endl;

    if(fabs(q3.x - (q1.x + q2.x) / 2.) > epsilon) {

      Parabola <Type> parabola = get_parabola_equation(q1, q2, q3) ;
      cout << " ****** Third point q3 changes position from (" << q3.x << "," << q3.y << ") to (" ;
      q3.x = (q1.x + q2.x) / 2. ;
      q3.y = - parabola.k * q3.x * q3.x - parabola.b * q3.x - parabola.d ;
      cout << q3.x << "," << q3.y << ") .******** " << endl;

      p3.x = 1. - q3.x ;
      p3.y = q3.y ;

      std::cout << "parabola " << parabola.k << "x^2+" << parabola.b << "x+" << parabola.d << " + y = 0 " << std::endl;
    }

    if(fabs(p1.x - 0.) < epsilon) {
      if(fabs(p2.x + p2.y - 1.) < epsilon) {
        actual_table = 0;
        old_table = 1;
        searchP = {static_cast<double>(p1.y), static_cast<double>(p2.x), static_cast<double>(p3.y)};
        //swap
        q1 = {(1. - p2.x), p2.y};
        q2 = {(1. - p1.x), p1.y};

      }
      else if(fabs(p2.y - 0.) < epsilon) {
        actual_table = 1;
        old_table = 3;
        searchP = {static_cast<double>(p1.y), static_cast<double>(p2.x), static_cast<double>(p3.y)};
      }
    }

    else if(fabs(p2.x - 0.) < epsilon) {
      if(fabs(p1.x + p1.y - 1.) < epsilon) {
        actual_table = 0;
        old_table = 1;
        searchP = {static_cast<double>(p2.y), static_cast<double>(p1.x), static_cast<double>(p3.y)};
      }
      else if(fabs(p1.y - 0.) < epsilon) {
        actual_table = 1;
        old_table = 3;
        searchP = {static_cast<double>(p2.y), static_cast<double>(p1.x), static_cast<double>(p3.y)};
        //swap
        q1 = {(1. - p2.x), p2.y};
        q2 = {(1. - p1.x), p1.y};
      }
    }

    else if(fabs(p1.x + p1.y - 1.) < epsilon) {
      if(fabs(p2.y - 0.) < epsilon) {
        actual_table = 2;
        old_table = 2;
        searchP = {static_cast<double>(p1.x), static_cast<double>(p2.x), static_cast<double>(p3.y)};
      }
    }
    else if(fabs(p2.x + p2.y - 1.) < epsilon) {
      if(fabs(p1.y - 0.) < epsilon) {
        actual_table = 2;
        old_table = 2;
        searchP = {static_cast<double>(p2.x), static_cast<double>(p1.x), static_cast<double>(p3.y)};
        //swap
        q1 = {(1. - p2.x), p2.y};
        q2 = {(1. - p1.x), p1.y};
      }
    }
  }
  else { //Horizontal
    q1 = {(1. - p1.y), p1.x};
    q2 = {(1. - p2.y), p2.x};
    q3 = {(1. - p3.y), p3.x};

    cout << "q1 = (" << q1.x << "," << q1.y << ")" << endl;
    cout << "q2 = (" << q2.x << "," << q2.y << ")" << endl;
    cout << "q3 = (" << q3.x << "," << q3.y << ")" << endl;
    if(fabs(q3.x - (q1.x + q2.x) / 2.) > epsilon) {
      Parabola <Type> parabola = get_parabola_equation(q1, q2, q3) ;
      cout << " ****** Third point q3 changes position from (" << q3.x << "," << q3.y << ") to (" ;
      q3.x = (q1.x + q2.x) / 2. ;
      q3.y = - parabola.k * q3.x * q3.x - parabola.b * q3.x - parabola.d ;
      cout << q3.x << "," << q3.y << ") .******** " << endl;
      p3.y = 1. - q3.x ;
      p3.x = q3.y;

      std::cout << "parabola " << parabola.k << "x^2+" << parabola.b << "x+" << parabola.d << " + y = 0 " << std::endl;
    }


    if(fabs(p1.x - 0.) < epsilon) {
      if(fabs(p2.x + p2.y - 1.) < epsilon) {
        actual_table = 3;
        old_table = 2;
        searchP = {static_cast<double>(p2.y), static_cast<double>(p1.y), static_cast<double>(p3.x)};
        //swap
        q1 = {(1. - p2.y), p2.x};
        q2 = {(1. - p1.y), p1.x};
      }
      else if(fabs(p2.y - 0.) < epsilon) {
        actual_table = 4;
        old_table = 3;
        searchP = {static_cast<double>(p2.x), static_cast<double>(p1.y), static_cast<double>(p3.x)};
        //swap
        q1 = {(1. - p2.y), p2.x};
        q2 = {(1. - p1.y), p1.x};
      }
    }
    else if(fabs(p2.x - 0.) < epsilon) {
      if(fabs(p1.x + p1.y - 1.) < epsilon) {
        actual_table = 3;
        old_table = 2;
        searchP = {static_cast<double>(p1.y), static_cast<double>(p2.y), static_cast<double>(p3.x)};
      }
      else if(fabs(p1.y - 0.) < epsilon) {
        actual_table = 4;
        old_table = 3;
        searchP = {static_cast<double>(p1.x), static_cast<double>(p2.y), static_cast<double>(p3.x)};
      }
    }

    else if(fabs(p1.x + p1.y - 1.) < epsilon) {
      if(fabs(p2.y - 0.) < epsilon) {
        actual_table = 5;
        old_table = 1;
        searchP = {static_cast<double>(p2.x), static_cast<double>(p1.y), static_cast<double>(p3.x)};
      }
    }
    else if(fabs(p2.x + p2.y - 1.) < epsilon) {
      if(fabs(p1.y - 0.) < epsilon) {
        actual_table = 5;
        old_table = 1;
        searchP = {static_cast<double>(p1.x), static_cast<double>(p2.y), static_cast<double>(p3.x)};
        //swap
        q1 = {(1. - p2.y), p2.x};
        q2 = {(1. - p1.y), p1.x};
      }
    }
  }
}

double GaussIntegral(const int &xExp, const int &yExp, const double* xg, const double* yg, const std::vector<double> &interp_point_weights, const double* gaussWeight) {
  double Integral = 0;
  for(unsigned ig = 0; ig < interp_point_weights.size(); ig++) {
    Integral += pow(xg[ig], xExp) * pow(yg[ig], yExp) * interp_point_weights[ig] * gaussWeight[ig];
  }
  return Integral;
}

template <class Type>
class OctreeNode {
  public:
    std::vector<Point3D> corners;  // All 8 corners of the node
    bool isLeaf;
    std::vector<OctreeNode> children;
    std::vector<std::vector<double>> cornerAreas;
    std::vector<std::vector<double>> cornerWeights;
    std::vector<std::vector<double>> midWeights;
    int table;
    unsigned depth;
    unsigned qM;
    int s = 0;
    Type a = 0;
    double relative_error = -1;
    double relative_error_opposite = -1;
    CutFemWeightParabola<double, Type>* _Pweights;

    OctreeNode(const std::vector<Point3D>& _corners, const int& _table, const int& _depth, const unsigned& _qM, CutFemWeightParabola<double, Type>* Pweights)
      : corners(_corners), isLeaf(true), table(_table), depth(_depth), qM(_qM), _Pweights(Pweights) {
      if(corners.size() != 8) {
        throw std::invalid_argument("OctreeNode must be initialized with exactly 8 corners");
      }
    }

    ~OctreeNode() {
      // No need to manually delete children
    }

    void getCorners() {
      // Corners are already stored, so we just need to calculate areas and weights
      Type area(0);
      Type c(1);
      cornerAreas.resize(8);
      cornerWeights.resize(8);
      PointT<Type> p1, p2, p3;

      for(size_t i = 0; i < corners.size(); ++i) {
        const auto& corner = corners[i];
        std::vector<double> corner_vec = {corner.x, corner.y, corner.z};
        get_p1_p2_p3(table, corner_vec, p1, p2, p3);

        int count = 0;
        for(unsigned qq = 0; qq <= qM; qq++) {
          for(unsigned jj = 0; jj <= qq; jj++) {
            unsigned ii = qq - jj;
            area = find_trig_area_2intersection_formula_first(ii, jj, s, a, c, table, p1, p2, p3);
            cornerAreas[i].push_back(static_cast<double>(area));
            count++;
          }
        }
        (*_Pweights)(s, a, c, table, p1, p2, p3, cornerWeights[i]);
      }
    }

    void getmiddlepoints_weights() {
      // Calculate the midpoints first
      std::vector<Point3D> midpoints(19);
      calculateMidpoints(midpoints);

      // Resize the midWeights vector to hold weights for all 19 midpoints
      midWeights.resize(19);

      // Calculate weights for each midpoint
      Type c(1);
      PointT<Type> p1, p2, p3;

      for(size_t i = 0; i < midpoints.size(); ++i) {
        const auto& midpoint = midpoints[i];
        std::vector<double> midpoint_vec = {midpoint.x, midpoint.y, midpoint.z};

        // Get p1, p2, p3 for this midpoint
        get_p1_p2_p3(table, midpoint_vec, p1, p2, p3);

        // Use _Pweights to calculate the weights for this midpoint
        (*_Pweights)(s, a, c, table, p1, p2, p3, midWeights[i]);
      }
    }

    void subdivideWithRelativeError_old(int maxDepth, double maxRelativeError, int currentDepth = 0) {
      if(currentDepth >= maxDepth || !isLeaf) {
        getCorners();
        return;
      }

      getCorners();

      // Calculate midpoints for subdivision
      std::vector<Point3D> midpoints(19);
      calculateMidpoints(midpoints);

      std::vector<double> relativeErrors;
      std::vector<double> relativeErrorsOpposite;
      for(const auto& midpoint : midpoints) {
        std::vector<double> interp_point = {midpoint.x, midpoint.y, midpoint.z};
        Type f_area(0);
        Type c = 1;
        PointT<Type> p1, p2, p3;
        get_p1_p2_p3(table, interp_point, p1, p2, p3);
        std::vector<std::vector<double>> interpolation_vector(8);
        int count = 0;
//               for (unsigned qq = 0; qq <= qM; qq++) {
        for(unsigned qq = 0; qq <= 0; qq++) {    //just based on area
          for(unsigned jj = 0; jj <= qq; jj++) {
            unsigned ii = qq - jj;
            for(size_t ic = 0; ic < corners.size(); ++ic) {
              interpolation_vector[ic] = {corners[ic].x, corners[ic].y, corners[ic].z, cornerAreas[ic][count]};
            }

            // Use the new deformed interpolation function


            double interp_area = trilinier_interpolation_FEM_orientation<double>(table, interpolation_vector, interp_point);

            f_area = find_trig_area_2intersection_formula_first(jj, ii, s, a, c, table, p1, p2, p3);
            double formula_area = static_cast<double>(f_area);
            double r_error = fabs(formula_area - interp_area) / formula_area;
            double r_error_opposite = fabs(formula_area - interp_area) / (1. / (ii + jj + 2.) * (jj + 1.) - formula_area);
            relativeErrors.push_back(r_error);
            relativeErrorsOpposite.push_back(r_error_opposite);
//                       count++;
          }
        }
      }
      relative_error = *std::max_element(relativeErrors.begin(), relativeErrors.end());
      relative_error_opposite = *std::max_element(relativeErrorsOpposite.begin(), relativeErrorsOpposite.end());

      bool force_subdevide = false;

      if(table == 0 || table == 5) {
        if(corners[1].x > 0.9999999999) {
          if(corners[1].y < 0.0000000001) {
            force_subdevide = true;
          }
        }
      }
      else if(table == 1 || table == 4) {
        if(corners[0].x < 0.0000000001) {
          if(corners[0].y < 0.0000000001) {
            force_subdevide = true;
          }
        }
      }

      else if(table == 2 || table == 3) {
        if(corners[2].x > 0.9999999999) {
          if(corners[2].y > 0.9999999999) {
            force_subdevide = true;
          }
        }
      }

//         if((table == 0 || table == 5) && corners[1].x > 0.9999 && corners[1].y < 0.0001) force_subdevide = true;
//         else if((table == 1 || table == 4) && corners[0].x < 0.0001 && corners[0].y < 0.0001) force_subdevide = true;
//         else if ((table == 2 || table ==3) && corners[2].x > 0.9999 && corners[2].y > 0.9999) force_subdevide = true;


      if(depth <= 3 || force_subdevide || relative_error > maxRelativeError || relative_error_opposite > maxRelativeError) {
        isLeaf = false;
        children.reserve(children.size() + 8);
        std::vector<std::vector<Point3D>> childCorners = subdivideCorners();
        for(const auto& childCorner : childCorners) {
          children.emplace_back(childCorner, table, depth + 1, qM, _Pweights);
        }

        for(auto& child : children) {
          child.subdivideWithRelativeError(maxDepth, maxRelativeError, currentDepth + 1);
        }

//             if(force_subdevide){
//                for (auto& child : children) {
//                 child.subdivideWithRelativeError(10, maxRelativeError, currentDepth + 1);
//                }
//             }
//             else{
//               for (auto& child : children) {
//                 child.subdivideWithRelativeError(maxDepth, maxRelativeError, currentDepth + 1);
//               }
//             }
      }
    }


    void subdivideWithRelativeError_now(int maxDepth, double maxRelativeError, int currentDepth = 0) {
      if(!isLeaf) {
        getCorners();
        return;  // Already subdivided, no need to process further
      }

      getCorners();

      // Determine if this octant contains the target corner we want to force subdivide to depth 10
      bool containsTargetCorner = false;

      if(table == 0 || table == 5) {
        // Check if this octant contains the corner at x=1, y=0
        for(const auto& corner : corners) {
          if(corner.x > 0.9999999999 && corner.y < 0.0000000001) {
            containsTargetCorner = true;
            break;
          }
        }
      }
      else if(table == 1 || table == 4) {
        // Check if this octant contains the corner at x=0, y=0
        for(const auto& corner : corners) {
          if(corner.x < 0.0000000001 && corner.y < 0.0000000001) {
            containsTargetCorner = true;
            break;
          }
        }
      }
      else if(table == 2 || table == 3) {
        // Check if this octant contains the corner at x=1, y=1
        for(const auto& corner : corners) {
          if(corner.x > 0.9999999999 && corner.y > 0.9999999999) {
            containsTargetCorner = true;
            break;
          }
        }
      }

      // Exit early if we've reached maximum depth and this isn't a target corner
      if(currentDepth >= maxDepth && !containsTargetCorner) {
        return;
      }

      // Exit early if we've reached depth 10 for target corners
      if(currentDepth >= 9 && containsTargetCorner) {
        return;
      }

      getCorners();

      // Calculate midpoints for subdivision
      std::vector<Point3D> midpoints(19);
      calculateMidpoints(midpoints);

      // We already determined containsTargetCorner above
      bool shouldForceSubdivide = containsTargetCorner && currentDepth < 9;

      std::vector<double> relativeErrors;
      std::vector<double> relativeErrorsOpposite;
      for(const auto& midpoint : midpoints) {
        std::vector<double> interp_point = {midpoint.x, midpoint.y, midpoint.z};
        Type f_area(0);
        Type c = 1;
        PointT<Type> p1, p2, p3;
        get_p1_p2_p3(table, interp_point, p1, p2, p3);
        std::vector<std::vector<double>> interpolation_vector(8);
        int count = 0;
        for(unsigned qq = 0; qq <= 0; qq++) {    //just based on area
          for(unsigned jj = 0; jj <= qq; jj++) {
            unsigned ii = qq - jj;
            for(size_t ic = 0; ic < corners.size(); ++ic) {
              interpolation_vector[ic] = {corners[ic].x, corners[ic].y, corners[ic].z, cornerAreas[ic][count]};
            }

            // Use the new deformed interpolation function
            double interp_area = trilinier_interpolation_FEM_orientation<double>(table, interpolation_vector, interp_point);

            f_area = find_trig_area_2intersection_formula_first(jj, ii, s, a, c, table, p1, p2, p3);
            double formula_area = static_cast<double>(f_area);



            double r_error = fabs(formula_area - interp_area) / fabs(formula_area);
//                     if(isnan(r_error)) r_error = 1.0 ;
            double r_error_opposite = fabs(formula_area - interp_area) / (0.5 - fabs(formula_area));
//                     if(isnan(r_error_opposite)) r_error_opposite = 1.0 ;
            relativeErrors.push_back(r_error);
            relativeErrorsOpposite.push_back(r_error_opposite);

            // Handle potential division by zero to avoid NaN values
//                     double r_error = 0.0;
//                     if (fabs(formula_area) > 0.000000000001) {  // Check if denominator is not too close to zero
//                         r_error = fabs(formula_area - interp_area) / fabs(formula_area);
//                     }
//                     else r_error = 1.0 ;
            // For opposite error calculation
//                     double denom_opposite = 1./ ((ii + jj + 2.) * (jj + 1.)) - formula_area;
//                     double r_error_opposite = 0.0;
//                     if (fabs(denom_opposite) > 0.000000000001) {  // Check if denominator is not too close to zero
//                         r_error_opposite = fabs(formula_area - interp_area) / fabs(denom_opposite);
//                     }
//                     else r_error_opposite = 1.0 ;
            /*
                                if(!isnan(r_error)){
                                  relativeErrors.push_back(r_error);
                                }
                                if(!isnan(r_error_opposite)){
                                  relativeErrorsOpposite.push_back(r_error_opposite);
                                }*/
          }
        }
      }

      relative_error = *std::max_element(relativeErrors.begin(), relativeErrors.end());
      relative_error_opposite = *std::max_element(relativeErrorsOpposite.begin(), relativeErrorsOpposite.end());

      // Decide whether to subdivide based on depth, error, or target corner
      bool shouldSubdivide = (currentDepth < 3) ||
                             (relative_error > maxRelativeError) ||
                             (relative_error_opposite > maxRelativeError) ||
                             shouldForceSubdivide;  // Force target corner to subdivide to depth 10

      if(shouldSubdivide) {
        isLeaf = false;
        children.reserve(children.size() + 8);
        std::vector<std::vector<Point3D>> childCorners = subdivideCorners();

        for(const auto& childCorner : childCorners) {
          children.emplace_back(childCorner, table, depth + 1, qM, _Pweights);
        }

        // Recursively subdivide all children
        for(auto& child : children) {
          // Pass the same maxDepth for all children
          // The per-child containsTargetCorner check will handle special treatment for target corners
          child.subdivideWithRelativeError(maxDepth, maxRelativeError, currentDepth + 1);
        }
      }
    }

    void subdivideWithRelativeError(int maxDepth, double maxRelativeError, int currentDepth = 0) {
      if(!isLeaf) {
        getCorners();
        getmiddlepoints_weights();
        return;  // Already subdivided, no need to process further
      }

      getCorners();
      getmiddlepoints_weights();

      // Determine if this octant contains the target corner we want to force subdivide to depth 10
      bool containsTargetCorner = false;

      if(table == 0 || table == 5) {
        // Check if this octant contains the corner at x=1, y=0
        for(const auto& corner : corners) {
          if(corner.x > 0.9999999999 && corner.y < 0.0000000001) {
            containsTargetCorner = true;
            break;
          }
        }
      }
      else if(table == 1 || table == 4) {
        // Check if this octant contains the corner at x=0, y=0
        for(const auto& corner : corners) {
          if(corner.x < 0.0000000001 && corner.y < 0.0000000001) {
            containsTargetCorner = true;
            break;
          }
        }
      }
      else if(table == 2 || table == 3) {
        // Check if this octant contains the corner at x=1, y=1
        for(const auto& corner : corners) {
          if(corner.x > 0.9999999999 && corner.y > 0.9999999999) {
            containsTargetCorner = true;
            break;
          }
        }
      }

      // Exit early if we've reached maximum depth and this isn't a target corner
      if(currentDepth >= maxDepth && !containsTargetCorner) {
        return;
      }

      // Exit early if we've reached depth 10 for target corners
      if(currentDepth >= 9 && containsTargetCorner) {
        return;
      }

//       getCorners();

      // We already determined containsTargetCorner above
      bool shouldForceSubdivide = containsTargetCorner && currentDepth < 9;

      double average_area = 0.0;
      double max_area = cornerAreas[0][0];
      double min_area = cornerAreas[0][0];

      // Calculate average, find max and min
      for(size_t i = 0; i < 8; ++i) {
        average_area += cornerAreas[i][0];
        if(cornerAreas[i][0] > max_area) max_area = cornerAreas[i][0];
        if(cornerAreas[i][0] < min_area) min_area = cornerAreas[i][0];
      }
      average_area /= 8.0;

      relative_error          = std::fabs(max_area - min_area) / average_area;
      relative_error_opposite = std::fabs(max_area - min_area) / (0.5 - average_area);

      double denominator1 = 2.0 * std::atan(100.0 * (average_area + 0.0001)) / M_PI;
      double denominator2 = 2.0 * std::atan(100.0 * ((0.5 - average_area) + 0.0001)) / M_PI;

      // Calculate relative errors with new denominators
      relative_error = std::fabs(max_area - min_area) / denominator1;
      relative_error_opposite = std::fabs(max_area - min_area) / denominator2;



      // Calculate midpoints for subdivision
//         std::vector<Point3D> midpoints(19);
//         calculateMidpoints(midpoints);
//
//         std::vector<double> relativeErrors;
//         std::vector<double> relativeErrorsOpposite;
//         for (const auto& midpoint : midpoints) {
//             std::vector<double> interp_point = {midpoint.x, midpoint.y, midpoint.z};
//             Type f_area(0);
//             Type c = 1;
//             PointT<Type> p1, p2, p3;
//             get_p1_p2_p3(table, interp_point, p1, p2, p3);
//             std::vector<std::vector<double>> interpolation_vector(8);
//             int count = 0;
//             for (unsigned qq = 0; qq <= 0; qq++) {   //just based on area
//                 for (unsigned jj = 0; jj <= qq; jj++) {
//                     unsigned ii = qq - jj;
//                     for (size_t ic = 0; ic < corners.size(); ++ic) {
//                         interpolation_vector[ic] = {corners[ic].x, corners[ic].y, corners[ic].z, cornerAreas[ic][count]};
//                     }
//
//                     // Use the new deformed interpolation function
//                     double interp_area = trilinier_interpolation_FEM_orientation<double>(table, interpolation_vector, interp_point);
//
//                     f_area = find_trig_area_2intersection_formula_first(jj, ii, s, a, c, table, p1, p2, p3);
//                     double formula_area = static_cast<double>(f_area);
//
//
//
//                     double r_error = fabs(formula_area - interp_area) / fabs(formula_area);
// //                     if(isnan(r_error)) r_error = 1.0 ;
//                     double r_error_opposite = fabs(formula_area - interp_area) / (0.5 - fabs(formula_area));
// //                     if(isnan(r_error_opposite)) r_error_opposite = 1.0 ;
//                     relativeErrors.push_back(r_error);
//                     relativeErrorsOpposite.push_back(r_error_opposite);
//
//                     // Handle potential division by zero to avoid NaN values
// //                     double r_error = 0.0;
// //                     if (fabs(formula_area) > 0.000000000001) {  // Check if denominator is not too close to zero
// //                         r_error = fabs(formula_area - interp_area) / fabs(formula_area);
// //                     }
// //                     else r_error = 1.0 ;
//                     // For opposite error calculation
// //                     double denom_opposite = 1./ ((ii + jj + 2.) * (jj + 1.)) - formula_area;
// //                     double r_error_opposite = 0.0;
// //                     if (fabs(denom_opposite) > 0.000000000001) {  // Check if denominator is not too close to zero
// //                         r_error_opposite = fabs(formula_area - interp_area) / fabs(denom_opposite);
// //                     }
// //                     else r_error_opposite = 1.0 ;
// /*
//                     if(!isnan(r_error)){
//                       relativeErrors.push_back(r_error);
//                     }
//                     if(!isnan(r_error_opposite)){
//                       relativeErrorsOpposite.push_back(r_error_opposite);
//                     }*/
//                 }
//             }
//         }

//           relative_error = *std::max_element(relativeErrors.begin(), relativeErrors.end());
//           relative_error_opposite = *std::max_element(relativeErrorsOpposite.begin(), relativeErrorsOpposite.end());

      // Decide whether to subdivide based on depth, error, or target corner
      /*        bool shouldSubdivide = (currentDepth < 3) ||
                                    (box_relative_error > 0.1) ||
                                    (box_relative_error_opposite > 0.1) ||
                                    (relative_error > maxRelativeError) ||
                                    (relative_error_opposite > maxRelativeError) ||
                                    shouldForceSubdivide;*/  // Force target corner to subdivide to depth 10

      bool shouldSubdivide = (currentDepth < 3) ||
                             (relative_error > maxRelativeError) ||
                             (relative_error_opposite > maxRelativeError) ||
                             shouldForceSubdivide;  // Force target corner to subdivide to depth 10

      if(shouldSubdivide) {
        isLeaf = false;
        children.reserve(children.size() + 8);
        std::vector<std::vector<Point3D>> childCorners = subdivideCorners();

        for(const auto& childCorner : childCorners) {
          children.emplace_back(childCorner, table, depth + 1, qM, _Pweights);
        }

        // Recursively subdivide all children
        for(auto& child : children) {
          // Pass the same maxDepth for all children
          // The per-child containsTargetCorner check will handle special treatment for target corners
          child.subdivideWithRelativeError(maxDepth, maxRelativeError, currentDepth + 1);
        }
      }
    }

    OctreeNode* search(const Point3D& point) {
      // First check if point is even in this node
      if(contains(point)) {
        std::cout << "\nFound containing node at depth " << depth << ":\n";
        std::cout << "Point: (" << point.x << ", " << point.y << ", " << point.z << ")\n";
        std::cout << "Rel_error =" << relative_error << ", " << relative_error_opposite << "\n Node corners:\n";
        for(size_t i = 0; i < corners.size(); ++i) {
          std::cout << "Corner " << i << ": ("
                    << corners[i].x << ", "
                    << corners[i].y << ", "
                    << corners[i].z << ") "
                    << cornerAreas[i][0] << "\n";
        }

        // If this is a leaf node, we're done
        if(isLeaf) {
          std::cout << "This is a leaf node - returning\n";
          return this;
        }

        // If not a leaf, check children
//         std::cout << "Checking " << children.size() << " children\n";
        for(auto& child : children) {
          OctreeNode* result = child.search(point);
          if(result != nullptr) {
            return result;
          }
        }

        // If no children contain the point, return this node
        std::cout << "No children contain point - returning this node\n";
        return this;
      }

      // Point not in this node
      return nullptr;
    }

    void saveOctreeToCSV(const std::string& filename) const {
      std::ofstream ofs(filename, std::ios::binary);
      if(!ofs.is_open()) {
        std::cerr << "Error: Unable to open file for writing: " << filename << std::endl;
        return;
      }
      serialize(ofs);
      ofs.close();
    }

    void loadOctreeFromCSV(const std::string& filename) {
      std::ifstream ifs(filename, std::ios::binary);
      if(!ifs.is_open()) {
        std::cerr << "Error: Unable to open file for reading: " << filename << std::endl;
        return;
      }
      deserialize(ifs);
      ifs.close();
    }

  private:

    void serialize(std::ofstream& ofs) const {
      ofs.write(reinterpret_cast<const char*>(&isLeaf), sizeof(isLeaf));
      ofs.write(reinterpret_cast<const char*>(&depth), sizeof(depth));
      ofs.write(reinterpret_cast<const char*>(&relative_error), sizeof(relative_error));
      ofs.write(reinterpret_cast<const char*>(&relative_error_opposite), sizeof(relative_error_opposite));

      serializeVector(ofs, corners);
      serializeVector(ofs, cornerAreas);
      serializeVector(ofs, cornerWeights);

      // Serialize midWeights
      serializeVector(ofs, midWeights);

      size_t childCount = children.size();
      ofs.write(reinterpret_cast<const char*>(&childCount), sizeof(childCount));
      for(const auto& child : children) {
        child.serialize(ofs);
      }
    }

    void deserialize(std::ifstream& ifs) {
      ifs.read(reinterpret_cast<char*>(&isLeaf), sizeof(isLeaf));
      ifs.read(reinterpret_cast<char*>(&depth), sizeof(depth));
      ifs.read(reinterpret_cast<char*>(&relative_error), sizeof(relative_error));
      ifs.read(reinterpret_cast<char*>(&relative_error_opposite), sizeof(relative_error_opposite));

      deserializeVector(ifs, corners);
      deserializeVector(ifs, cornerAreas);
      deserializeVector(ifs, cornerWeights);

      // Deserialize midWeights
      deserializeVector(ifs, midWeights);

      size_t childCount;
      ifs.read(reinterpret_cast<char*>(&childCount), sizeof(childCount));
      children.clear();
      children.reserve(childCount);
      for(size_t i = 0; i < childCount; ++i) {
        children.emplace_back(corners, 0, 0, 0, nullptr);
        children.back().deserialize(ifs);
      }
    }

    void serializeVector(std::ofstream& ofs, const std::vector<std::vector<double>>& vec) const {
      size_t size = vec.size();
      ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));

      for(const auto& innerVec : vec) {
        size_t innerSize = innerVec.size();
        ofs.write(reinterpret_cast<const char*>(&innerSize), sizeof(innerSize));
        ofs.write(reinterpret_cast<const char*>(innerVec.data()), innerSize * sizeof(double));
      }
    }


    void deserializeVector(std::ifstream& ifs, std::vector<std::vector<double>>& vec) {
      size_t size;
      ifs.read(reinterpret_cast<char*>(&size), sizeof(size));

      vec.resize(size);
      for(size_t i = 0; i < size; ++i) {
        size_t innerSize;
        ifs.read(reinterpret_cast<char*>(&innerSize), sizeof(innerSize));
        vec[i].resize(innerSize);
        ifs.read(reinterpret_cast<char*>(vec[i].data()), innerSize * sizeof(double));
      }
    }

    void serializeVector(std::ofstream& ofs, const std::vector<Point3D>& vec) const {
      size_t size = vec.size();
      ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
      for(const auto& point : vec) {
        ofs.write(reinterpret_cast<const char*>(&point.x), sizeof(double));
        ofs.write(reinterpret_cast<const char*>(&point.y), sizeof(double));
        ofs.write(reinterpret_cast<const char*>(&point.z), sizeof(double));
      }
    }

    void deserializeVector(std::ifstream& ifs, std::vector<Point3D>& vec) {
      size_t size;
      ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
      vec.resize(size);
      for(size_t i = 0; i < size; ++i) {
        ifs.read(reinterpret_cast<char*>(&vec[i].x), sizeof(double));
        ifs.read(reinterpret_cast<char*>(&vec[i].y), sizeof(double));
        ifs.read(reinterpret_cast<char*>(&vec[i].z), sizeof(double));
      }
    }

    Point3D middleof(Point3D& point1, Point3D& point2) {
      Point3D midpoint{
        (point1.x + point2.x) / 2.,
        (point1.y + point2.y) / 2.,
        (point1.z + point2.z) / 2.
      };
      return midpoint;
    }

    void calculateMidpoints(std::vector<Point3D>& midpoints) {
      midpoints.resize(19);

      // Midpoints of bottom face edges (counter-clockwise)
      midpoints[0] = middleof(corners[0], corners[1]); // Bottom edge: front
      midpoints[1] = middleof(corners[1], corners[2]); // Bottom edge: right
      midpoints[2] = middleof(corners[2], corners[3]); // Bottom edge: back
      midpoints[3] = middleof(corners[3], corners[0]); // Bottom edge: left

      // Midpoints of top face edges (counter-clockwise)
      midpoints[4] = middleof(corners[4], corners[5]); // Top edge: front
      midpoints[5] = middleof(corners[5], corners[6]); // Top edge: right
      midpoints[6] = middleof(corners[6], corners[7]); // Top edge: back
      midpoints[7] = middleof(corners[7], corners[4]); // Top edge: left

      // Midpoints of vertical edges
      midpoints[8] = middleof(corners[0], corners[4]);  // Vertical edge: front-left
      midpoints[9] = middleof(corners[1], corners[5]);  // Vertical edge: front-right
      midpoints[10] = middleof(corners[2], corners[6]); // Vertical edge: back-right
      midpoints[11] = middleof(corners[3], corners[7]); // Vertical edge: back-left

      // Midpoints of faces
      midpoints[12] = middleof(midpoints[0], midpoints[2]); // Bottom face center
      midpoints[13] = middleof(midpoints[4], midpoints[6]); // Top face center
      midpoints[14] = middleof(midpoints[0], midpoints[4]); // Front face center
      midpoints[15] = middleof(midpoints[1], midpoints[5]); // Right face center
      midpoints[16] = middleof(midpoints[2], midpoints[6]); // Back face center
      midpoints[17] = middleof(midpoints[3], midpoints[7]); // Left face center

      // Center point of the octree node
      midpoints[18] = middleof(midpoints[12], midpoints[13]); // Center of the cube
    }


    std::vector<std::vector<Point3D>> subdivideCorners() {
      std::vector<std::vector<Point3D>> childCorners(8, std::vector<Point3D>(8));
      std::vector<Point3D> midpoints(19);
      calculateMidpoints(midpoints);

      // Child 0: Bottom-front-left (following counter-clockwise convention)
      childCorners[0] = {
        corners[0],      // Bottom-front-left
        midpoints[0],    // Bottom-front-middle
        midpoints[12],   // Bottom-center
        midpoints[3],    // Bottom-left-middle
        midpoints[8],    // Front-left-middle
        midpoints[14],   // Front-center
        midpoints[18],   // Center
        midpoints[17]    // Left-center
      };

      // Child 1: Bottom-front-right
      childCorners[1] = {
        midpoints[0],    // Bottom-front-middle
        corners[1],      // Bottom-front-right
        midpoints[1],    // Bottom-right-middle
        midpoints[12],   // Bottom-center
        midpoints[14],   // Front-center
        midpoints[9],    // Front-right-middle
        midpoints[15],   // Right-center
        midpoints[18]    // Center
      };

      // Child 2: Bottom-back-right
      childCorners[2] = {
        midpoints[12],   // Bottom-center
        midpoints[1],    // Bottom-right-middle
        corners[2],      // Bottom-back-right
        midpoints[2],    // Bottom-back-middle
        midpoints[18],   // Center
        midpoints[15],   // Right-center
        midpoints[10],   // Back-right-middle
        midpoints[16]    // Back-center
      };

      // Child 3: Bottom-back-left
      childCorners[3] = {
        midpoints[3],    // Bottom-left-middle
        midpoints[12],   // Bottom-center
        midpoints[2],    // Bottom-back-middle
        corners[3],      // Bottom-back-left
        midpoints[17],   // Left-center
        midpoints[18],   // Center
        midpoints[16],   // Back-center
        midpoints[11]    // Back-left-middle
      };

      // Child 4: Top-front-left
      childCorners[4] = {
        midpoints[8],    // Front-left-middle
        midpoints[14],   // Front-center
        midpoints[18],   // Center
        midpoints[17],   // Left-center
        corners[4],      // Top-front-left
        midpoints[4],    // Top-front-middle
        midpoints[13],   // Top-center
        midpoints[7]     // Top-left-middle
      };

      // Child 5: Top-front-right
      childCorners[5] = {
        midpoints[14],   // Front-center
        midpoints[9],    // Front-right-middle
        midpoints[15],   // Right-center
        midpoints[18],   // Center
        midpoints[4],    // Top-front-middle
        corners[5],      // Top-front-right
        midpoints[5],    // Top-right-middle
        midpoints[13]    // Top-center
      };

      // Child 6: Top-back-right
      childCorners[6] = {
        midpoints[18],   // Center
        midpoints[15],   // Right-center
        midpoints[10],   // Back-right-middle
        midpoints[16],   // Back-center
        midpoints[13],   // Top-center
        midpoints[5],    // Top-right-middle
        corners[6],      // Top-back-right
        midpoints[6]     // Top-back-middle
      };

      // Child 7: Top-back-left
      childCorners[7] = {
        midpoints[17],   // Left-center
        midpoints[18],   // Center
        midpoints[16],   // Back-center
        midpoints[11],   // Back-left-middle
        midpoints[7],    // Top-left-middle
        midpoints[13],   // Top-center
        midpoints[6],    // Top-back-middle
        corners[7]       // Top-back-left
      };

      return childCorners;
    }

    bool contains(const Point3D& point) const {
      const double EPSILON = 1e-10;

      // Helper function to calculate dot product
      auto dot = [](const Point3D & a, const Point3D & b) -> double {
        return a.x * b.x + a.y * b.y + a.z * b.z;
      };

      // Helper function to calculate cross product
      auto cross = [](const Point3D & a, const Point3D & b) -> Point3D {
        return Point3D{
          a.y * b.z - a.z * b.y,
          a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x
        };
      };

      // Helper function to create vector from two points
      auto makeVector = [](const Point3D & from, const Point3D & to) -> Point3D {
        return Point3D{to.x - from.x, to.y - from.y, to.z - from.z};
      };

      // Define the six faces of the hexahedron
      // Each face is defined by four corners in counter-clockwise order when viewed from outside
      // This ensures the normal is pointing outward
      const std::vector<std::vector<int>> faces = {
        {0, 3, 2, 1},    // bottom face (viewed from below)
        {4, 5, 6, 7},    // top face
        {0, 1, 5, 4},    // front face
        {1, 2, 6, 5},    // right face
        {2, 3, 7, 6},    // back face
        {3, 0, 4, 7}     // left face
      };

      // Check if point is on the correct side of all faces
      for(const auto& face : faces) {
        // Get three points from the face to define the plane
        const Point3D& v0 = corners[face[0]];
        const Point3D& v1 = corners[face[1]];
        const Point3D& v2 = corners[face[2]];

        // Calculate face normal using cross product (pointing outward)
        Point3D edge1 = makeVector(v0, v1);
        Point3D edge2 = makeVector(v0, v2);
        Point3D normal = cross(edge1, edge2);

        // Calculate signed distance from point to plane
        Point3D toPoint = makeVector(v0, point);
        double signedDist = dot(normal, toPoint);

        // For a point inside the octree, the signed distance should be negative
        // (point is on the opposite side of the face normal)
        if(signedDist > EPSILON) {
          return false;
        }
      }

      return true;
    }
};

template <class Type>
void generateAndLoadOctrees(const int &maxDepth, const int &degree, const double &percent, CutFemWeightParabola<double, Type> &Pweights, std::vector<OctreeNode<Type>>& loadedRoots) {


  std::vector<std::vector<Point3D>> initialCorners = {
    // Bottom face: counter-clockwise starting from (0,0,0)
    // Then top face: counter-clockwise starting from point above (0,0,0)
    { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 0.5}, {0.0, 1.0, 0.5}
    },

    { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 0.5}, {0.0, 1.0, 0.5}
    },

    { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}, {1.0, 0.0, 0.5}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.5}
    },

    { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}, {1.0, 0.0, 0.5}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.5}
    },  //table 3 is same as case 2

    { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 0.5}, {0.0, 1.0, 0.5}
    },  //table 4 is same as case 1

    { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 0.5}, {0.0, 1.0, 0.5}
    }   //table 5 is same as case 0
  };
  for(int ttable = 0; ttable < 6; ++ttable) {
    std::string filename = "save/octree_table_" + std::to_string(ttable) + "_maxdepth_" + std::to_string(maxDepth) + "_per_" + std::to_string(percent) + "_degree_" + std::to_string(degree) + ".csv";

    FILE *fp;
    fp = fopen(filename.c_str(), "r");
    if(fp != NULL) {
      fclose(fp);
    }
    else {
      cout << "creating the tables" << ttable << endl;

      OctreeNode<Type> root(initialCorners[ttable], ttable, 0, degree, &Pweights);
      root.subdivideWithRelativeError(maxDepth, percent);

      root.saveOctreeToCSV(filename);
      std::cout << "Octree Structure:\n";
    }
  }

  // Load the octree structure and vectors from the CSV file
  cout << "Loading the tables" << endl;
  loadedRoots.clear();
  loadedRoots.reserve(2);
  for(int ttable = 0; ttable < 6; ++ttable) {

    loadedRoots.emplace_back(initialCorners[ttable], ttable, 0, degree, nullptr);
    loadedRoots.back().loadOctreeFromCSV("save/octree_table_" + std::to_string(ttable) + "_maxdepth_" + std::to_string(maxDepth) + "_per_" + std::to_string(percent) + "_degree_" + std::to_string(degree) + ".csv");
  }
}




template <class Type>
void printOctreeStructure(const OctreeNode<Type>& node, int depth = 0) {
  // Print indentation based on depth
  std::string indent(depth * 2, ' ');

  // Print node information
  std::cout << indent << "Node at depth " << node.depth << ": ";
  std::cout << "Is Leaf: " << (node.isLeaf ? "Yes" : "No") << "\n";

  // Print relative error information for all nodes (not just leaves)
  std::cout << indent << "  Relative Error: " << node.relative_error << "\n";
  std::cout << indent << "  Relative Error Opposite: " << node.relative_error_opposite << "\n";

  // Print corner information
  std::cout << indent << "  Corners:\n";
  for(size_t i = 0; i < node.corners.size(); ++i) {
    std::cout << indent << "    Corner " << i << ": ("
              << node.corners[i].x << ", "
              << node.corners[i].y << ", "
              << node.corners[i].z << ") = " << node.cornerAreas[i][0] << "\n";
  }
  for(size_t i = 0; i < node.midWeights.size(); ++i) {
    std::cout << indent << "    midpoint weights " << i << ": ("
              << node.midWeights[i][0] << node.midWeights[i][1] << node.midWeights[i][2] << "\n";
  }

  // If it's a leaf node, print detailed corner areas and weights
  if(node.isLeaf) {
//         std::cout << indent << "  Corner Areas:\n";
//         for (size_t i = 0; i < node.cornerAreas.size(); ++i) {
//             std::cout << indent << "    Corner " << i << " Areas:\n";
//             for (size_t j = 0; j < node.cornerAreas[i].size(); ++j) {
//                 std::cout << indent << "      [" << j << "]: "
//                           << node.cornerAreas[i][j] << "\n";
//             }
//         }
//
//         std::cout << indent << "  Corner Weights:\n";
//         for (size_t i = 0; i < node.cornerWeights.size(); ++i) {
//             std::cout << indent << "    Corner " << i << " Weights:\n";
//             for (size_t j = 0; j < node.cornerWeights[i].size(); ++j) {
//                 std::cout << indent << "      [" << j << "]: "
//                           << node.cornerWeights[i][j] << "\n";

//             }
//         }
  }
  else {
    // If not a leaf, recursively print children
    for(const auto& child : node.children) {
      printOctreeStructure(child, depth + 1);
    }
  }
}


int checkVectorRelation(const std::vector<int>& vec1, const std::vector<int>& vec2) {
  // Check if the sizes of the vectors are different
  if(vec1.size() != vec2.size()) {
    std::cerr << "Warning: Number of sign do not match" << std::endl;
    return 0;
  }

  int equalCount = 0;
  int negativeCount = 0;

  // Iterate through the vectors to count equal and negative elements
  for(size_t i = 0; i < vec1.size(); ++i) {
    if(vec1[i] == vec2[i]) {
      ++equalCount;
    }
    if(vec1[i] == -vec2[i]) {
      ++negativeCount;
    }
  }

  if(equalCount > negativeCount) {
    return 1;
  }
  else if(negativeCount > equalCount) {
    return -1;
  }
  else {
    std::cerr << "Warning: Equal number of positive and negative sign on the corner" << std::endl;
    return 0;
  }
}


// Utility function to check if a point is inside the circle
bool isPointInCircle(double x, double y, double centerX, double centerY, double radius) {
  double dx = x - centerX;
  double dy = y - centerY;
  return (dx * dx + dy * dy) <= radius * radius;
}

// Function to find circle-line intersection
std::vector<std::pair<double, double>> findCircleLineIntersection(
                                      double x1, double y1, double x2, double y2,
double centerX, double centerY, double radius) {

  std::vector<std::pair<double, double>> intersections;

  // Translate circle center to origin
  x1 -= centerX;
  y1 -= centerY;
  x2 -= centerX;
  y2 -= centerY;

  // Calculate line direction vector
  double dx = x2 - x1;
  double dy = y2 - y1;

  // Calculate quadratic equation coefficients
  // (tx + x1)² + (ty + y1)² = r², where t is parameter along line
  double a = dx * dx + dy * dy;
  double b = 2 * (x1 * dx + y1 * dy);
  double c = x1 * x1 + y1 * y1 - radius * radius;

  // Calculate discriminant
  double discriminant = b * b - 4 * a * c;

  if(discriminant >= 0) {
    // Calculate intersection parameters
    double t1 = (-b + sqrt(discriminant)) / (2 * a);
    double t2 = (-b - sqrt(discriminant)) / (2 * a);

    // Calculate intersection points
    double ix1 = x1 + t1 * dx + centerX;
    double iy1 = y1 + t1 * dy + centerY;
    double ix2 = x1 + t2 * dx + centerX;
    double iy2 = y1 + t2 * dy + centerY;

    // Check if points lie on the line segment
    auto isOnSegment = [](double x, double y, double x1, double y1, double x2, double y2) {
      return x >= std::min(x1, x2) && x <= std::max(x1, x2) &&
             y >= std::min(y1, y2) && y <= std::max(y1, y2);
    };

    if(t1 >= 0 && t1 <= 1 &&
        isOnSegment(ix1, iy1, x1 + centerX, y1 + centerY, x2 + centerX, y2 + centerY)) {
      intersections.push_back({ix1, iy1});
    }

    if(discriminant > 0 && t2 >= 0 && t2 <= 1 &&
        isOnSegment(ix2, iy2, x1 + centerX, y1 + centerY, x2 + centerX, y2 + centerY)) {
      intersections.push_back({ix2, iy2});
    }

    if(intersections.size() == 2) {



    }
  }

  return intersections;
}


// Modified function to find circle-triangle intersection points with vertical flag
std::vector<std::pair<double, double>> findCircleTriangleIntersectionPoints(
                                      double x1, double y1, double x2, double y2, double x3, double y3,
double centerX, double centerY, double radius, bool& vertical) {

  std::vector<std::pair<double, double>> intersections;

  // Check each edge of the triangle
  auto points1 = findCircleLineIntersection(x1, y1, x2, y2, centerX, centerY, radius);
  auto points2 = findCircleLineIntersection(x2, y2, x3, y3, centerX, centerY, radius);
  auto points3 = findCircleLineIntersection(x3, y3, x1, y1, centerX, centerY, radius);

  // Collect all intersection points
  intersections.insert(intersections.end(), points1.begin(), points1.end());
  intersections.insert(intersections.end(), points2.begin(), points2.end());
  intersections.insert(intersections.end(), points3.begin(), points3.end());

  // Sort and remove duplicates
  std::sort(intersections.begin(), intersections.end());
  intersections.erase(std::unique(intersections.begin(), intersections.end()), intersections.end());

  // If we have exactly two intersection points, add a third point
  if(intersections.size() == 2) {
    double x_diff = std::abs(intersections[1].first - intersections[0].first);
    double y_diff = std::abs(intersections[1].second - intersections[0].second);

    vertical = (x_diff >= y_diff);

    if(vertical) {
      // Find third point using x-midpoint
      double x3 = (intersections[0].first + intersections[1].first) / 2.0;
      double y3 = centerY + sqrt(radius * radius - pow(x3 - centerX, 2));
      // Choose the y value that's between the two intersection points
      if(y3 < std::min(intersections[0].second, intersections[1].second) ||
          y3 > std::max(intersections[0].second, intersections[1].second)) {
        y3 = centerY - sqrt(radius * radius - pow(x3 - centerX, 2));
      }
      intersections.push_back({x3, y3});
    }
    else {
      // Find third point using y-midpoint
      double y3 = (intersections[0].second + intersections[1].second) / 2.0;
      double x3 = centerX + sqrt(radius * radius - pow(y3 - centerY, 2));
      // Choose the x value that's between the two intersection points
      if(x3 < std::min(intersections[0].first, intersections[1].first) ||
          x3 > std::max(intersections[0].first, intersections[1].first)) {
        x3 = centerX - sqrt(radius * radius - pow(y3 - centerY, 2));
      }
      intersections.push_back({x3, y3});
    }
  }

  return intersections;
}

// Modified function to check circle-triangle intersection
bool doesCircleIntersectTriangle(double x1, double y1, double x2, double y2, double x3, double y3,
                                 double centerX, double centerY, double radius) {
  // Check if any vertex is inside the circle
  bool v1Inside = isPointInCircle(x1, y1, centerX, centerY, radius);
  bool v2Inside = isPointInCircle(x2, y2, centerX, centerY, radius);
  bool v3Inside = isPointInCircle(x3, y3, centerX, centerY, radius);

  // If all vertices are inside, return false (fully contained)
  if(v1Inside && v2Inside && v3Inside) {
    return false;
  }

  // Check each edge for intersection
  auto edge1 = findCircleLineIntersection(x1, y1, x2, y2, centerX, centerY, radius);
  auto edge2 = findCircleLineIntersection(x2, y2, x3, y3, centerX, centerY, radius);
  auto edge3 = findCircleLineIntersection(x3, y3, x1, y1, centerX, centerY, radius);

  // If any edge intersects, return true
  return !edge1.empty() || !edge2.empty() || !edge3.empty();
}

// Function to transform points from physical to reference triangle
void transformToReferenceTriangle(
  double px1, double py1, double px2, double py2, double px3, double py3,  // Physical triangle
  double x, double y,  // Point to transform
  double& xi, double& eta)  // Reference coordinates
{
  // Calculate transformation matrix components
  double J11 = px2 - px1;
  double J12 = px3 - px1;
  double J21 = py2 - py1;
  double J22 = py3 - py1;

  // Calculate determinant
  double det = J11 * J22 - J12 * J21;

  // Transform point
  double dx = x - px1;
  double dy = y - py1;

  xi = (J22 * dx - J12 * dy) / det;
  eta = (-J21 * dx + J11 * dy) / det;
}

// Helper function to print triangle state
void printTriangleState(int triangleIndex,
                        double x1, double y1, double x2, double y2, double x3, double y3,
                        double centerX, double centerY, double radius,
                        const std::vector<std::pair<double, double>>& intersections,
                        const std::vector<std::pair<double, double>>& refIntersections,
                        double totalarea, double area, bool vertical) {
//     std::cout << "\n=== Triangle " << triangleIndex << " ===" << std::endl;
  std::cout << "Physical coordinates: " << std::endl;
  std::cout << "V1: (" << x1 << ", " << y1 << ")" << std::endl;
  std::cout << "V2: (" << x2 << ", " << y2 << ")" << std::endl;
  std::cout << "V3: (" << x3 << ", " << y3 << ")" << std::endl;




  // Determine triangle state
  bool v1Inside = isPointInCircle(x1, y1, centerX, centerY, radius);
  bool v2Inside = isPointInCircle(x2, y2, centerX, centerY, radius);
  bool v3Inside = isPointInCircle(x3, y3, centerX, centerY, radius);

  if(v1Inside && v2Inside && v3Inside) {
    std::cout << "State: INSIDE" << std::endl;
  }
  else if(!intersections.empty()) {
    std::cout << "State: CUT" << std::endl;
  }
  else {
    std::cout << "State: OUTSIDE" << std::endl;
  }

  // Print intersection points
  if(!intersections.empty()) {
    std::cout << "\nPhysical intersection points:" << std::endl;
    for(const auto& p : intersections) {
      std::cout << "(" << p.first << ", " << p.second << ")" << std::endl;
    }

    std::cout << "\nReference intersection points:" << std::endl;
    for(const auto& p : refIntersections) {
      std::cout << "(" << p.first << ", " << p.second << ")" << std::endl;
    }


    // Calculate and print parabola points
//         std::cout << "\nParabola points (physical space):" << std::endl;
//         for (double t = 0; t <= 1; t += 0.1) {
//             // Linear interpolation between intersection points
//             double x = intersections[0].first + t * (intersections[1].first - intersections[0].first);
//             // Calculate y using circle equation
//             double y = centerY + sqrt(radius * radius - (x - centerX) * (x - centerX));
//             // Check if this y-value is between the intersection points
//             if (y >= std::min(intersections[0].second, intersections[1].second) &&
//                 y <= std::max(intersections[0].second, intersections[1].second)) {
//                 std::cout << "(" << x << ", " << y << ")" << std::endl;
//             }
//         }
  }
  std::cout << "vertical: " << vertical << std::endl;
  std::cout << "area: " << area << std::endl;

  std::cout << "Cumulative area: " << totalarea << std::endl;
}


template <class Type>
void GetTriquadraticInterpolationVector(int table, const std::vector< std::vector< Type >> & vtxCrdn, const std::vector< std::vector< Type >> & vtxValues, const std::vector< std::vector< Type >> & mdlValues, const std::vector< Type > &intPoint, std::vector< Type > &intValues) {

  intValues.resize(vtxValues[0].size());
  const Type & x = intPoint[0];
  const Type & y = intPoint[1];
  const Type z = (table == 2 || table == 3) ? - (2. * intPoint[2]) / (x + y - 2.) : -(2. * intPoint[2]) / (y - 2.);


  // const Type & x = 0.125 * (vtxCrdn[0][0] + vtxCrdn[1][0] + vtxCrdn[2][0] + vtxCrdn[3][0] + vtxCrdn[4][0] + vtxCrdn[5][0] + vtxCrdn[6][0] + vtxCrdn[7][0]);
  // const Type & y = 0.125 * (vtxCrdn[0][1] + vtxCrdn[1][1] + vtxCrdn[2][1] + vtxCrdn[3][1] + vtxCrdn[4][1] + vtxCrdn[5][1] + vtxCrdn[6][1] + vtxCrdn[7][1]);
  //
  // const Type & zz = 0.125 * (vtxCrdn[0][2] + vtxCrdn[1][2] + vtxCrdn[2][2] + vtxCrdn[3][2] + vtxCrdn[4][2] + vtxCrdn[5][2] + vtxCrdn[6][2] + vtxCrdn[7][2]);
  // const Type z = (table == 2 || table == 3) ? - (2. * zz) / (x + y - 2.) : -(2. * zz) / (y - 2.);


  const Type & x0 = vtxCrdn[0][0];
  const Type & y0 = vtxCrdn[0][1];
  const Type z0 = (table == 2 || table == 3) ? - (2. * vtxCrdn[0][2]) / (vtxCrdn[0][0] + vtxCrdn[0][1] - 2.) : -(2. * vtxCrdn[0][2]) / (vtxCrdn[0][1] - 2.);

  const Type & x1 = vtxCrdn[1][0];
  const Type & y1 = vtxCrdn[7][1];
  const Type z1 = (table == 2 || table == 3) ? - (2. * vtxCrdn[7][2]) / (vtxCrdn[7][0] + vtxCrdn[7][1] - 2.) : -(2. * vtxCrdn[7][2]) / (vtxCrdn[7][1] - 2.);

  Type s1 = (x - x0) / (x1 - x0);
  Type s2 = (y - y0) / (y1 - y0);
  Type s3 = (z - z0) / (z1 - z0);



  //abort();

  double phi[3], phj[3], phk[3];
  phi[0] = (1. - s1) * (1. - 2. * s1);
  phi[1] = 4. * s1 * (1. -  s1);
  phi[2] = s1 * (2. * s1 - 1.);

  phj[0] = (1. - s2) * (1. - 2. * s2);
  phj[1] = 4. * s2 * (1. -  s2);
  phj[2] = s2 * (2. * s2 - 1.);

  phk[0] = (1. - s3) * (1. - 2. * s3);
  phk[1] = 4. * s3 * (1. -  s3);
  phk[2] = s3 * (2. * s3 - 1.);

  // std::cerr<<table<<" "<< s1 <<" "<< s2<<" "<<s3<<std::endl;
  // std::cerr<<table<<" "<< phi[0] <<" "<< phi[1]<<" "<<phi[2]<<std::endl;
  // std::cerr<<table<<" "<< phj[0] <<" "<< phj[1]<<" "<<phj[2]<<std::endl;
  // std::cerr<<table<<" "<< phk[0] <<" "<< phk[1]<<" "<<phk[2]<<std::endl;


  for(unsigned k = 0; k < vtxValues[0].size(); k++) {

    double W[27];
    for(unsigned j = 0; j < 8; j++)     {
      W[j] = vtxValues[j][k];
    }
    for(unsigned j = 8; j < 27; j++)      {
      W[j] = mdlValues[j - 8][k] ;
    }

    double b[3][3];
    b[0][0] = W[0] * phi[0] + W[8] * phi[1] + W[1] * phi[2];
    b[0][1] = W[11] * phi[0] + W[20] * phi[1] + W[9] * phi[2];
    b[0][2] = W[3] * phi[0] + W[10] * phi[1] + W[2] * phi[2];

    b[1][0] = W[16] * phi[0] + W[22] * phi[1] + W[17] * phi[2];
    b[1][1] = W[25] * phi[0] + W[26] * phi[1] + W[23] * phi[2];
    b[1][2] = W[19] * phi[0] + W[24] * phi[1] + W[18] * phi[2];

    b[2][0] = W[4] * phi[0] + W[12] * phi[1] + W[5] * phi[2];
    b[2][1] = W[15] * phi[0] + W[21] * phi[1] + W[13] * phi[2];
    b[2][2] = W[7] * phi[0] + W[14] * phi[1] + W[6] * phi[2];

    // std::cerr << b[0][0] - W[8] << " ";
    // std::cerr << b[0][1] - W[20] << " ";
    // std::cerr << b[0][2] - W[10] << std::endl;
    //
    //  std::cerr << b[1][0] - W[22] << " ";
    // std::cerr << b[1][1] - W[26] << " ";
    // std::cerr << b[1][2] - W[24] << std::endl;


    double a[3] = {0., 0., 0.};

    for(unsigned i = 0; i < 3; i++) {
      for(unsigned j = 0; j < 3; j++) {
        a[i] += b[i][j] * phj[j];
      }
    }

    // std::cerr << a[0] - W[20] << " ";
    // std::cerr << a[1] - W[26] << " ";
    // std::cerr << a[2] - W[21] << " ";

    intValues[k] = a[0] * phk[0] + a[1] * phk[1] +  a[2] * phk[2];
    //std::cerr << intValues[k] - W[26] << std::endl;
  }
  //std::cerr<<std::endl;
}

#include "FemusInit.hpp"
using namespace femus;



int main (int argc, char** args) {

  FemusInit mpinit (argc, args, MPI_COMM_WORLD);


  typedef cpp_bin_float_oct Type;
  unsigned int m = 0;
  unsigned int n = 0;
  int s = 0;
  Type k, b, d, a = 0, c = 1;

  std::cout.precision(16);

  PointT <Type> p1, p2, p3;
  p1 = { static_cast<Type>(0.), static_cast<Type>(0.4471) };
  p2 = { static_cast<Type>(0.4471), static_cast<Type>(0.) };
  p3 = { static_cast<Type>((p1.x + p2.x) / 2.0), static_cast<Type>(0.5) };

  std::vector<double>weightCF;
  CutFemWeightParabola <double, Type> Pweights(TRI, 3, "legendre");
  Pweights(s, a, c, 0, p1, p2, p3, weightCF);
  std::vector<std::vector<double>>badcase;
  std::vector<std::vector<double>>case_summary;

  // Circle parameters
  double centerX = 0.5;
  double centerY = 0.5;
  double radius = 0.397;
//  double radius = 0.19;
//  Mesh parameters



  int triangleIndex = 0;
  int fourintersection = 0;

  cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= " << endl;

  Fem fem = Fem(3 * 2, 2);
  unsigned TRI = 4;
  unsigned linear = 0;
  const elem_type *femQuad = fem.GetFiniteElement(TRI, linear);
  unsigned nInt;

  unsigned nPoints = 3;
  unsigned dim = 2;
  short unsigned ielType = 4; //quad = 3 and tri = 4
  unsigned femType = 0; //linear FEM

  std::vector<std::vector<double>> xv ;
  std::vector<std::vector<double>> unitxv = {{0., 1., 0.}, {0., 0., 1.}};
//     std::vector<std::vector<double>> xv = {{0.125, 0.125, 0.}, {0.25, 0.375, 0.375}};
  std::vector<double> A = {1., 0., 1., -1., -1., 0.342391}; //ax^2+bxy+cy^2+dx+ey+f =0
//     std::vector<double> A = {1., 0., 1., -1., -1., 0.46};


  int maxDepth = 5;
  int degree = 1;
  double percent = -1;
  //   std::vector<OctreeNode<Type>> roots;
  std::vector<OctreeNode<Type>>loadedRoots;

  generateAndLoadOctrees<Type>(maxDepth, degree, percent, Pweights, loadedRoots);

//     printOctreeStructure(loadedRoots[0],0);
//     return 1;
  // Loop through the mesh

  int nd = 4;  // Number of divisions per side

  unsigned nlevel = 7;
  std::vector<double> ndLevel(nlevel);
  std::vector<double> absError(nlevel);

  for(unsigned level = 0; level < nlevel; level++, nd *= 2) {
    double totalArea = 0.0;
    double totalformulaArea = 0.0;
    double h = 1.0 / nd;  // Grid spacing
    for(int i = 0; i < nd; i++) {
      for(int j = 0; j < nd; j++) {
        for(int t = 0; t < 2; t++) {
          triangleIndex++;

//           cout << " ======= Triangle ====== " << triangleIndex << endl;
          double x1, y1, x2, y2, x3, y3, area ;
          bool normal = true;

          if(t == 0) {
            // First triangle
            x1 = i * h;
            y1 = j * h;
            x2 = (i + 1) * h;
            y2 = j * h;
            x3 = i * h;
            y3 = (j + 1) * h;
          }
          else {
            // Second triangle
            x1 = (i + 1) * h;
            y1 = j * h;
            x2 = (i + 1) * h;
            y2 = (j + 1) * h;
            x3 = i * h;
            y3 = (j + 1) * h;
          }

//               if (t == 0) {
//                     x2 = i * h;     y2 = j * h;
//                     x3 = (i+1) * h; y3 = j * h;
//                     x1 = i * h;     y1 = (j+1) * h;
//                 } else {
//                     x2 = (i+1) * h; y2 = j * h;
//                     x3 = (i+1) * h; y3 = (j+1) * h;
//                     x1 = i * h;     y1 = (j+1) * h;
//                 }


//               if (t == 0) {
//                     x3 = i * h;     y3 = j * h;
//                     x1 = (i+1) * h; y1 = j * h;
//                     x2 = i * h;     y2 = (j+1) * h;
//                 } else {
//                     x3 = (i+1) * h; y3 = j * h;
//                     x1 = (i+1) * h; y1 = (j+1) * h;
//                     x2 = i * h;     y2 = (j+1) * h;
//                 }

          xv = {{x1, x2, x3}, {y1, y2, y3}};  //physical triangle

          cout << "Physical triangle = (" << x1 << ", " << y1 << "), (" << x2 << ", " << y2 << "), (" << x3 << ", " << y3 << ") " << endl;

          bool vertical = false;
          std::vector<std::pair<double, double>> refIntersections;

          if(doesCircleIntersectTriangle(x1, y1, x2, y2, x3, y3, centerX, centerY, radius)) {
            auto intersections = findCircleTriangleIntersectionPoints(
                                   x1, y1, x2, y2, x3, y3, centerX, centerY, radius, vertical);

            cout << "Vertical (intersection)? " << vertical <<  "intersection = " << intersections.size() << endl;

            if(intersections.size() == 4) {
              std::cout << "Triangle " << triangleIndex << " has four intersections - skipping : We assume the tri is inside. area = " << 0.5 * h*h << std::endl;
              totalArea += 0.5 * h * h;  // Add full triangle area
              continue;
            }

            if(intersections.size() == 3) {
              // Transform to reference triangle points
              PointT<Type> q1, q2, q3;
              std::vector<int> circlesign(3);
              std::vector<int> parabolasign(3);
              int actual_table = 9999;
              int old_table = 999;
              Parabola <Type> parabola;
              Point3D searchP;

              std::pair<std::vector<std::vector<double>>, std::vector<double>> xp = GetCellPointsFromQuadric(xv, A, nPoints, nInt);
              std::vector < std::vector < std::vector <double > > > aP(1);
              ProjectNodalToPolynomialCoefficients(aP[femType], xv, ielType, femType);

              std::vector<int> xvsign(3);
              std::vector<int> unitxvsign(3);
              std::vector<std::vector<double>> xi(nPoints, std::vector<double>(2, 0.));

              for(unsigned i = 0; i < nPoints; i++) {
                bool inverseMapping = GetInverseMapping(femType, ielType, aP, xp.first[i], xi[i], 100);        //This maps the phsical points to {(-1,-1),(1,1)} box for quad. For triangle it maps to (0,0),(1,0),(0,1)
                std::cout << " \nx[i] physical value " << i << " " << xp.first[i][0] << " " << xp.first[i][1] << std::endl;
                std::cout << " x[i] value in reference " << i << " (" << xi[i][0] << ", " << xi[i][1] << ")" << std::endl;

                //       xi[i] = {0.5 * (xi[i][0] + 1.), 0.5 * (xi[i][1] + 1.)};                                        // //This maps the points to unit box. For quad
                //       std::cout << "value in reference triangle : check if it works" << i << " " << xi[i][0] << " " << xi[i][1] << std::endl;
              }

              cout << " xv sign = {" ;
              for(unsigned l = 0; l < 3 ; l++) {
                xvsign[l] = ((A[0] * xv[0][l] * xv[0][l] + A[1] * xv[0][l] * xv[1][l] + A[2] * xv[1][l] * xv[1][l] + A[3] * xv[0][l] + A[4] * xv[1][l] + A[5]) >= 0) ? 1 : -1 ;
                cout << xvsign[l] << ", ";
              }
              cout << "} " << endl;

              std::vector<double> phi, gradPhi;
              std::vector<double> Xg(femQuad->GetGaussPointNumber(), 0);
              std::vector<double> Yg(femQuad->GetGaussPointNumber(), 0);
              std::vector<double> Jg(femQuad->GetGaussPointNumber(), 0);
              for(unsigned ig = 0; ig < femQuad->GetGaussPointNumber(); ig++) {
                // *** get gauss point weight, test function and test function partial derivatives ***
                femQuad->Jacobian(xv, ig, Jg[ig], phi, gradPhi);
                for(unsigned i = 0; i < phi.size(); i++) {
                  Xg[ig] += phi[i] * xv[0][i];
                  Yg[ig] += phi[i] * xv[1][i];
                }
//                           std::cout <<"checking gauss points and jacobian"<<ig<<" "<<Xg[ig]<<" "<<Yg[ig]<<" "<<Jg[ig]<<std::endl;
              }

              //points in reference domain
              p1 = { static_cast<Type>(xi[0][0]), static_cast<Type>(xi[0][1]) };
              p2 = { static_cast<Type>(xi[2][0]), static_cast<Type>(xi[2][1]) };
              p3 = { static_cast<Type>(xi[1][0]), static_cast<Type>(xi[1][1]) };

              find_actual_table_trig<Type>(p1, p2, p3, actual_table, old_table, searchP, q1, q2, q3, vertical);

              cout << "Vertical (table)? = " << vertical << " actual_table_number = " << actual_table << " : points in actual table : ( " << p1.x << "," << p1.y << " )" << " , ( " << p2.x << "," << p2.y << " )" << " , ( " << p3.x << "," << p3.y << " ) &&&&&&  old_table_number" << old_table << " : points in old table : ( " << q1.x << "," << q1.y << " )" << " , ( " << q2.x << "," << q2.y << " )" << " , ( " << q3.x << "," << q3.y << " ) " << endl;

              //Just testing the signsxv
              if(vertical) {
                parabola = get_parabola_equation(p1, p2, p3);
                std::cout << "parabola " << parabola.k << "x^2+" << parabola.b << "x+" << parabola.d << " + y = 0 " << std::endl;
                //check the parabola sign

                cout << " unit triangle sign = {" ;
                for(unsigned l = 0; l < 3 ; l++) {
                  unitxvsign[l] = ((parabola.k * unitxv[0][l] * unitxv[0][l] + parabola.b * unitxv[0][l] + parabola.d + unitxv[1][l]) > 0) ? 1 : -1;
                  cout << unitxvsign[l] << ", ";
                }
                cout << "} " << endl;

//                           p3.x = (p1.x + p2.x) / 2;
//                           p3.y = -parabola.k * p3.x * p3.x - parabola.b * p3.x - parabola.d ;

              }
              else {
                cout << " it is a Horizontal parabola .......@.......@........@........" << endl;
                p1 = { static_cast<Type>(xi[0][1]), static_cast<Type>(xi[0][0]) };
                p2 = { static_cast<Type>(xi[2][1]), static_cast<Type>(xi[2][0]) };
                p3 = { static_cast<Type>(xi[1][1]), static_cast<Type>(xi[1][0]) };
                cout <<  "( " << p1.x << "," << p1.y << " )" << " , ( " << p2.x << "," << p2.y << " )" << " , ( " << p3.x << "," << p3.y << " ) " << endl;

                Parabola <Type> parabola = get_parabola_equation(p1, p2, p3);
                cout << parabola.k << "y^2+ " << parabola.b << "y+ " << parabola.d << "+x =0 " << endl;

                //use horizotal parabola for the normal
                cout << " unit triangle sign = {" ;
                for(unsigned l = 0; l < 3 ; l++) {
                  unitxvsign[l] = ((static_cast<double>(parabola.k) * unitxv[1][l] * unitxv[1][l] + static_cast<double>(parabola.b) * unitxv[1][l] + static_cast<double>(parabola.d) + unitxv[0][l]) > 0) ? 1 : -1;
                  cout << unitxvsign[l] << ", ";
                }
                cout << "} " << endl;
//                           p3.x = (p1.x + p2.x) / 2.;
//                           p3.y = -parabola.k * p3.x * p3.x - parabola.b * p3.x - parabola.d ;
              }

              int checksign = checkVectorRelation(xvsign, unitxvsign);   //TODO checkVectorRelation output int
              if(checksign == 1) {
                normal = false;
              }
              Type ref_formula_area(0);




              if(normal) {   //using second area beacause we want to use the formula on q directly after the transformation. TODO or should we always use first area integral with P
                ref_formula_area = find_trig_area_2intersection_formula_second<Type>(0, 0, 0, 0, 1, old_table,  q1,  q2, q3);
              }

              else {
                ref_formula_area = 0.5 - find_trig_area_2intersection_formula_second<Type>(0, 0, 0, 0, 1, old_table,  q1,  q2, q3);
              }



//                         if (normal){   //using second area beacause we want to use the formula on q directly after the transformation. TODO or should we always use first area integral with P
//                           ref_formula_area = find_trig_area_2intersection_formula_first<Type>(0, 0, 0, 0, 1, actual_table,  p1,  p2, p3);
//                         }
//
//                         else{
//                           ref_formula_area = 0.5 - find_trig_area_2intersection_formula_first<Type>(0, 0, 0, 0, 1, actual_table,  p1,  p2, p3);
//                         }






              totalformulaArea += static_cast<double>(ref_formula_area) * h * h;


              std::vector<double>weightCF, interp_point_weights, interp_point_integrals;

              // Search in octree
              OctreeNode<Type>* result = loadedRoots[actual_table].search(searchP);

              cout << "Search point " << searchP.x << " " << searchP.y << " " << searchP.z << endl;
              cout << "actual table " << actual_table << "result " << result << endl;

              if(result) {
                std::vector<double> interp_point = {searchP.x, searchP.y, searchP.z};
                std::vector<std::vector<double>> corners(8, std::vector<double>(3));

                for(size_t i = 0; i < result->corners.size(); ++i) {
                  corners[i][0] = result->corners[i].x;
                  corners[i][1] = result->corners[i].y;
                  corners[i][2] = result->corners[i].z;
                }


                PointT <Type> pp1, pp2, pp3;     //This is just to print the parabolas on the corner
                for(size_t i = 0; i < result->corners.size(); ++i) {
                  get_p1_p2_p3(actual_table, corners[i], pp1, pp2, pp3);
                  parabola = get_parabola_equation(pp1, pp2, pp3);
                  std::cout << parabola.k << "x^2+" << parabola.b << "x+" << parabola.d << " + y = 0 " << std::endl;
                }


                std::vector<double> interp_point_weights;
//                             trilinear_interpolation_vector_deformed(corners, result->cornerWeights,interp_point, interp_point_weights);


                std::cout << "\n interp Point: (" << interp_point[0] << ", " << interp_point[1] << ", " << interp_point[2] << ")\n";


//                             trilinear_interpolation_vector_deformed(corners, result->cornerAreas, interp_point, interp_point_weights);

                trilinear_interpolation_vector_remap_to_unitcube(actual_table, corners, result->cornerAreas, interp_point, interp_point_weights);

                std::cout << " interpolated integrals = ";
                for(size_t j = 0; j < interp_point_weights.size(); ++j) {
                  std::cout << interp_point_weights[j] << ", ";
                }
                std::cout << " )" << std::endl;

//                 trilinear_interpolation_vector_deformed(corners, result->cornerWeights, interp_point, interp_point_weights);

                /*
                                trilinear_interpolation_vector_remap_to_unitcube(actual_table,corners, result->cornerWeights, interp_point, interp_point_weights);*/

                GetTriquadraticInterpolationVector(actual_table, corners, result->cornerWeights, result->midWeights, interp_point, interp_point_weights);

                std::cout << " interpolated weights = ";
                for(size_t j = 0; j < interp_point_weights.size(); ++j) {
                  std::cout << interp_point_weights[j] << ", ";
                }
                std::cout << " )" << std::endl;

                std::vector<double>modified_weights(interp_point_weights.size());

                if(!normal) {
                  for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {

                    interp_point_weights[aq] = 1 - interp_point_weights[aq];
                  }
                }

                area = GaussIntegral(0, 0, Xg.data(), Yg.data(), interp_point_weights, Jg.data());


                totalArea += area;  // No need to scale it. jacobian does it for us.


//                             printTriangleState(triangleIndex, x1, y1, x2, y2, x3, y3,
//                                              centerX, centerY, radius,
//                                              intersections, refIntersections, totalArea,area,vertical);
//                             std::cout << "scaled area: " << area*h*h << std::endl;

                double reldif = fabs(area / h / h - static_cast<double>((ref_formula_area)));
                cout << "--" << endl;
                cout << "Ref quadrature area = " << area / h / h << endl;
                cout << " Ref formula area = " << ref_formula_area << endl;
                cout << "Physical quadrature area = " << area << endl;
                cout << " Physical formula  area = " << (ref_formula_area)*h*h << endl;
                cout << " difference between formula and quadrature = " << reldif << endl;
                cout << " Relative difference  " << reldif / static_cast<double>((ref_formula_area)) << endl;
                if(reldif / static_cast<double>((ref_formula_area)) > percent / 2.) {
                  std::vector<double> badcaseEntry = {
                    static_cast<double>(triangleIndex),
                    static_cast<double>(actual_table),
                    static_cast<double>(ref_formula_area),
                    static_cast<double>((ref_formula_area)*h * h),
                    static_cast<double>(area),
                    static_cast<double>(reldif),
                    static_cast<double>(reldif / static_cast<double>((ref_formula_area)))
                  };
                  badcase.push_back(badcaseEntry);

                }

                std::vector<double> caseEntry = {
                  static_cast<double>(triangleIndex),
                  static_cast<double>(actual_table),
//                                     static_cast<double>(ref_formula_area),
//                                     static_cast<double>((ref_formula_area)*h*h),
//                                     static_cast<double>(area),
                  reldif,
                  reldif / static_cast<double>((ref_formula_area))
                };
                case_summary.push_back(caseEntry);


                cout << "sign & table : " << checksign << " " << actual_table;
                (vertical) ? cout << "v \n\n" : cout << "h\n\n" << endl;

              }
            }
          }

          // If triangle is completely inside circle
          else if(isPointInCircle(x1, y1, centerX, centerY, radius) &&
                  isPointInCircle(x2, y2, centerX, centerY, radius) &&
                  isPointInCircle(x3, y3, centerX, centerY, radius)) {
            totalArea += 0.5 * h * h;  // Add full triangle area
            totalformulaArea += 0.5 * h * h;

            std::vector<std::pair<double, double>> emptyIntersections;
//                     printTriangleState(triangleIndex, x1, y1, x2, y2, x3, y3,
//                                      centerX, centerY, radius,
//                                      emptyIntersections, emptyIntersections, totalArea,area,vertical);
            std::cout << "It is inside. Scaled area: " << 0.5 * h*h << std::endl;
            cout << " Cumulative area = " << totalArea << endl;
          }
          else {
            std::vector<std::pair<double, double>> emptyIntersections;
//                     printTriangleState(triangleIndex, x1, y1, x2, y2, x3, y3,
//                                      centerX, centerY, radius,
//                                      emptyIntersections, emptyIntersections, totalArea,area,vertical);
            std::cout << "It is outside. Scaled area: " << 0 << std::endl;


          }
        }
      }
    }

    // Calculate errors
    double analyticalArea = M_PI * radius * radius;
    double absoluteError = fabs(totalArea - analyticalArea);
    double relativeError = absoluteError / analyticalArea;
    double formulaAbsError = fabs(totalformulaArea - analyticalArea);
    double formulaRelError = formulaAbsError / analyticalArea;

    cout << "Number of four intersection = " << fourintersection << endl;
    std::cout << "\n=== Final Results ===" << std::endl;
    std::cout << "Formula Area: " << totalformulaArea << std::endl;
    std::cout << "Numerical Area: " << totalArea << std::endl;
    std::cout << "Analytical Area: " << analyticalArea << std::endl;
    std::cout << "Numerical Absolute Error: " << absoluteError << std::endl;
    std::cout << "Relative Error: " << relativeError * 100 << "%" << std::endl;

    std::cout << "Formula Absolute Error: " << formulaAbsError << std::endl;
    std::cout << "formulaRelative Error: " << formulaRelError * 100 << "%" << std::endl;
//     cout<< " badcase = " ;
//     for(int i = 0; i < badcase.size(); i++){
//       std::cout << " " << badcase[i]<< std::endl;
//     }

//             std::cout << "All Cases:" << std::endl;
//         cout<< "index   table   difference            relative     "<<endl;
//         for (size_t i = 0; i < case_summary.size(); ++i) {
//             std::cout << "Case " << i << ": ";
//             for (const auto& value : case_summary[i]) {
//                 std::cout << value << " ";
//             }
//             std::cout << std::endl;
//         }

    std::cout << "Bad Cases:" << std::endl;
    cout << "index table   analytica_area(max1/2)      analytical_scale      interpolation_area         difference         relative     " << endl;
    for(size_t i = 0; i < badcase.size(); ++i) {
      std::cout << "Case " << i << ": ";
      for(const auto& value : badcase[i]) {
        std::cout << value << " ";
      }
      std::cout << std::endl;
    }

    ndLevel[level] = nd;
    absError[level] = absoluteError;

  }

  for(unsigned level = 0; level < nlevel; level++) {
    std::cout << level << " " << ndLevel[level] << " " << absError[level] << " " << std::endl;
  }

  return 0;
}



double GetBiquadraticInterpolation(const std::vector<double> & W, const std::vector<double> & s) {


  const double &s1 = s[0];
  const double &s2 = s[1];
  const double &s3 = s[2];

  double phi[3];
  phi[0] = (1. - s1) * (1. - 2. * s1);
  phi[1] = 4. * s1 * (1. -  s1);
  phi[2] = s1 * (2. * s1 - 1.);

  double b[3][3];
  b[0][0] = W[0] * phi[0] + W[8] * phi[1] + W[1] * phi[2];
  b[0][1] = W[11] * phi[0] + W[20] * phi[1] + W[9] * phi[2];
  b[0][2] = W[3] * phi[0] + W[10] * phi[1] + W[2] * phi[2];

  b[1][0] = W[16] * phi[0] + W[22] * phi[1] + W[17] * phi[2];
  b[1][1] = W[25] * phi[0] + W[26] * phi[1] + W[23] * phi[2];
  b[1][2] = W[19] * phi[0] + W[24] * phi[1] + W[18] * phi[2];

  b[1][0] = W[4] * phi[0] + W[12] * phi[1] + W[5] * phi[2];
  b[1][1] = W[15] * phi[0] + W[21] * phi[1] + W[13] * phi[2];
  b[1][2] = W[7] * phi[0] + W[14] * phi[1] + W[6] * phi[2];

  phi[0] = (1. - s2) * (1. - 2. * s2);
  phi[1] = 4. * s2 * (1. -  s2);
  phi[2] = s2 * (2. * s2 - 1.);
  double a[3] = {0., 0., 0.};
  for(unsigned i = 0; i < 3; i++) {
    for(unsigned j = 0; j < 3; j++) {
      a[i] += b[i][j] * phi[j];
    }
  }

  phi[0] = (1. - s3) * (1. - 2. * s3);
  phi[1] = 4. * s3 * (1. -  s3);
  phi[2] = s3 * (2. * s3 - 1.);
  double weight = a[0] * phi[0] + a[1] * phi[1] +  a[2] * phi[2];

  return weight;

}

template <class Type>
Type GetTriquadraticInterpolation(const int table, std::vector< std::vector< Type >> & interp_table, std::vector< std::vector< Type >> & midpoint_table, const std::vector< Type > &interp_point) {

  Type x = interp_point[0];
  Type y = interp_point[1];

  Type denominator(0), denominatorz0(0), denominatorz1(0);

  if(table < 2 || table > 3) denominator = y - 2;
  else denominator = x + y - 2;

  Type z = -(2 * interp_point[2]) / denominator;


  Type c_000 = interp_table[0][3];
  Type c_100 = interp_table[1][3];
  Type c_110 = interp_table[2][3];
  Type c_010 = interp_table[3][3];
  Type c_001 = interp_table[4][3];
  Type c_101 = interp_table[5][3];
  Type c_111 = interp_table[6][3];
  Type c_011 = interp_table[7][3];

  Type x0 = interp_table[0][0];
  Type x1 = interp_table[1][0];
  Type y0 = interp_table[0][1];
  Type y1 = interp_table[7][1];
  if(table < 2 || table > 3) {
    denominatorz0 = interp_table[0][1] - 2.;
    denominatorz1 = interp_table[7][1] - 2.;
  }
  else {
    denominatorz0  = interp_table[0][0] + interp_table[0][1] - 2.;
    denominatorz1 = interp_table[7][0] + interp_table[7][1] - 2.;
  }
  Type z0 = - (2 * interp_table[0][2]) / denominatorz0;
  Type z1 = - (2 * interp_table[7][2]) / denominatorz1;

  double W[27];
  for(unsigned j = 0; j < 27; j++) {
    if(j < 8) {
      W[j] = interp_table[j][3] ;
    }
    else {
      W[j] = midpoint_table[j][3] ;
    }
  }

  const double &s1 = (x - x0) / (x1 - x0);
  const double &s2 = (y - y0) / (y1 - y0);
  const double &s3 = (z - z0) / (z1 - z0);

  double phi[3];
  phi[0] = (1. - s1) * (1. - 2. * s1);
  phi[1] = 4. * s1 * (1. -  s1);
  phi[2] = s1 * (2. * s1 - 1.);

  double b[3][3];
  b[0][0] = W[0] * phi[0] + W[8] * phi[1] + W[1] * phi[2];
  b[0][1] = W[11] * phi[0] + W[20] * phi[1] + W[9] * phi[2];
  b[0][2] = W[3] * phi[0] + W[10] * phi[1] + W[2] * phi[2];

  b[1][0] = W[16] * phi[0] + W[22] * phi[1] + W[17] * phi[2];
  b[1][1] = W[25] * phi[0] + W[26] * phi[1] + W[23] * phi[2];
  b[1][2] = W[19] * phi[0] + W[24] * phi[1] + W[18] * phi[2];

  b[1][0] = W[4] * phi[0] + W[12] * phi[1] + W[5] * phi[2];
  b[1][1] = W[15] * phi[0] + W[21] * phi[1] + W[13] * phi[2];
  b[1][2] = W[7] * phi[0] + W[14] * phi[1] + W[6] * phi[2];


  phi[0] = (1. - s2) * (1. - 2. * s2);
  phi[1] = 4. * s2 * (1. -  s2);
  phi[2] = s2 * (2. * s2 - 1.);
  double a[3] = {0., 0., 0.};
  for(unsigned i = 0; i < 3; i++) {
    for(unsigned j = 0; j < 3; j++) {
      a[i] += b[i][j] * phi[j];
    }
  }

  phi[0] = (1. - s3) * (1. - 2. * s3);
  phi[1] = 4. * s3 * (1. -  s3);
  phi[2] = s3 * (2. * s3 - 1.);
  Type weight = a[0] * phi[0] + a[1] * phi[1] +  a[2] * phi[2];
  return weight;
}
