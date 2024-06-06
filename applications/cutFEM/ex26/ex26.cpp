 
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <algorithm>    // std::sort
#include <ctime>
#include <cstdlib>
#include <climits>
#include <typeinfo>
#include <fstream>
#include <sstream>
#include <string>


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

// Function to calculate the equation of a parabola
template <class Type>
Parabola <Type>  get_parabola_equation( PointT <Type> p1, PointT <Type> p2, PointT <Type> p3) {
    Type x1 = p1.x, x2 = p2.x, x3 = p3.x;
    Type y1 = p1.y, y2 = p2.y, y3 = p3.y;
    Type k,b,d;
    Type det = p1.x * p1.x * (p2.x - p3.x) -p1.x* (p2.x*p2.x - p3.x*p3.x)+ p2.x*p3.x*(p2.x - p3.x);

//     Type det = x1 * x1 * (x2 - x3) -x1* (x2*x2 - x3*x3)+ x2*x3*(x2 - x3) ;
//     Type denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
    if(fabs(det) >0.00000001){
       k = -(y1 * (x2 - x3) + y2 * (x3 - x1) + y3 * (x1 - x2)) / det;
       b = -(y1 * (x3*x3 - x2*x2) + y2 * (x1*x1 - x3*x3)+ y3 * ((x2*x2 - x1*x1))) / det;
       d = -(y1 * x2 * x3 * (x2 - x3) + y2 * x3 * x1 * (x3 - x1) + y3 * x1 * x2 * (x1 - x2)) / det;
    }
    else{
      Type slope = (y1-y2)/(x1-x2) ;
      k=0;
      b= -slope;
      d= x2*slope - y2 ;
    }

//     else {This will give us a straight line paralal to y axix.
//        k = -0.;
//        b = -1.;
//        d = p1.x ;
//     }

    if (fabs(k) < 1.e-12) k = 0 ;
    if (fabs(b) < 1.e-12) b = 0 ;
    if (fabs(d) < 1.e-12) d = 0 ;

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
void CheckIntersection(int &intersect_number, unsigned int &table_number , std::vector <Type> &intersection, std::vector <Type> &interp_point,  Parabola <Type> &parabola){

  table_number = -1;
  intersect_number=0;
  intersection.resize(0);
  interp_point.resize(0);
  Type k = parabola.k;
  Type b = parabola.b;
  Type d = parabola.d;
  Type c = 1;
  int left =0 , top = 0, right = 0, bottom = 0 ;
//   cout<< " parabola I get from solving system of linear equation :  " << parabola.k <<"x^2 + "<< parabola.b <<"x + "<< parabola.d << "+" <<  c << " y=0"  <<endl;

      if (-d>=0 && -d<=1){ //LEFT
        intersection.resize(intersection.size()+2);
        intersection[intersection.size()-2] = 0;
        intersection[intersection.size()-1] = -d;
        left = 1 ;
        intersect_number += 1;
        interp_point.resize(interp_point.size()+1);
        interp_point[interp_point.size()-1] = -d;
//         cout << " left = " << left ;
      }
      // LEFT-TOP solve kx^2+bx+d-1 =0  ; Table 0
      if (k == 0){
      Type x =  (-1-d)/b ;
        if(x <= 1 && x>= 0) {
          intersection.resize(intersection.size()+2);
          intersection[intersection.size()-2] = x;
          intersection[intersection.size()-1] = 1;
          interp_point.resize(interp_point.size()+1);
          interp_point[interp_point.size()-1] = x;
          top =1;
          intersect_number += 1;
          if (left ==1) table_number = 0 ;
//           cout << " top = " << top ;
        }
      }
      else {
        Type delta = b*b - 4*k*(d+1);
        cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
        if (delta >= 0){
              Type sqrtdelta = sqrt(delta);
              int sign = (k > 0) ? 1 : -1;

              for(unsigned i = 0; i < 2; i++) {
                Type x = (- b - sign * sqrtdelta) / (2 * k);
    //             cout<< "Top x = "<< x<< endl;
                if(x > 1) break;
                else if(x >= 0) {
                  intersection.resize(intersection.size()+2);
                  intersection[intersection.size()-2] = x;
                  intersection[intersection.size()-1] = 1;
                  interp_point.resize(interp_point.size()+1);
                  interp_point[interp_point.size()-1] = x;
                  intersect_number += 1;
                  if (top ==1){table_number = 3 ;}
                  top += 1;
                  if (left ==1){table_number = 0 ;}
    //               cout << " top = " << top ;
                }
                sign *= -1;
              }
            }
      }
      Type y_1=-(k+b+d); //LEFT-RIGHT x=1 ; Table 1
      if (y_1 >= 0 && y_1 <= 1){ //TODO check sign when normal change
          intersection.resize(intersection.size()+2);
          intersection[intersection.size()-2] = 1;
          intersection[intersection.size()-1] = y_1;
          interp_point.resize(interp_point.size()+1);
          interp_point[interp_point.size()-1] = y_1;
          intersect_number += 1;
          if (left ==1){table_number = 1 ;}
          if  (top >=1){table_number = 4 ;}
          right = 1 ;
//           cout << " right = " << right ;
      }

        //LEFT-BOTTOM  solve kx^2+bx+d =0 ; Table 2
      if (k == 0){
          Type x =  -d/b ;
          if(x <= 1 && x>= 0) {
            intersection.resize(intersection.size()+2);
            intersection[intersection.size()-2] = x;
            intersection[intersection.size()-1] = 0;
            interp_point.resize(interp_point.size()+1);
            interp_point[interp_point.size()-1] = x;
            intersect_number += 1;
            if (left == 1){table_number = 2 ;}
            if (right == 1){table_number = 6 ;}
            if (top >= 1){table_number = 5 ;}
            bottom = 1;
//             cout << " bottom = " << bottom ;
          }
      }

      else {
          Type delta = b*b - 4*k*d;
          if(delta >=0) {
            Type sqrtdelta = sqrt(delta);
            int sign = (k > 0) ? 1 : -1;

            for(unsigned i = 0; i < 2; i++) {
              Type x = (- b - sign * sqrtdelta) / (2 * k);
//               cout << " bottom root = " << x ;
              if(x > 1) break;
              else if(x >= 0) {
                intersection.resize(intersection.size()+2);
                intersection[intersection.size()-2] = x;
                intersection[intersection.size()-1] = 0;
                interp_point.resize(interp_point.size()+1);
                interp_point[interp_point.size()-1] = x;
                if (bottom >=1){table_number = 7 ;}
                if (left ==1){table_number = 2 ;}
                if (right ==1){table_number = 6 ;}
                if (top ==1){table_number = 5 ;}    // TODO check the table
                bottom += 1;
                intersect_number += 1;
//                 cout << " bottom = " << bottom ;
              }
              sign *= -1;
            }
          }
      }

      if (intersect_number == 4){
        if(left == 1 && top == 2 && right == 1){
          table_number = 0;
          Type swap;
          swap = interp_point[2];
          interp_point[2] = interp_point[3];
          interp_point[3] = swap;
        }
        else if (left == 1 && top == 2 && bottom == 1){
          table_number = 1;
          Type swap;
          swap = interp_point[1];
          interp_point[1] = interp_point[2];
          interp_point[2] = interp_point[3];
          interp_point[3] = swap;
        }

        else if(left == 1 && bottom == 2 && right ==1){
          table_number = 2;
          Type swap;
          swap = interp_point[2];
          interp_point[2] = interp_point[1];
          interp_point[1] = swap;
        }

        else if (left == 1 && bottom == 2 && top == 1){
          table_number = 3;
          Type swap;
          swap = interp_point[1];
          interp_point[1] = interp_point[2];
          interp_point[2] = interp_point[3];
          interp_point[3] = swap;
        }

        else if (right == 1 && top == 2 && bottom == 1){
          table_number = 4;
          Type swap;
          swap = interp_point[1];
          interp_point[1] = interp_point[0];
          interp_point[0] = interp_point[2];
          interp_point[2] = interp_point[3];
          interp_point[3] = swap;
        }

        else if (right == 1 &&  bottom == 2 && top == 1){
          table_number = 5;
          Type swap;
          swap = interp_point[2];
          interp_point[2] = interp_point[0];
          interp_point[0] = interp_point[1];
          interp_point[1] = interp_point[3];
          interp_point[3] = swap;
        }

        else if( top == 2 && bottom == 2){
//           if(bottom == 1) table_number = 2;
//           else if(top == 1 ) table_number =3;
          if (interp_point[0] > interp_point[3]){
            table_number = 6 ;
            Type swap;
            swap = interp_point[0];
            interp_point[0] = interp_point [3];
            interp_point[3] = swap;
          }
          else{
            table_number = 7;
            Type swap;
            swap = interp_point[1];
            interp_point[1] = interp_point [2];
            interp_point[2] = swap;

          }
        }

      }
/*
cout<< " " << " left " << left << " top "<< top << " right "<< right << " bottom " << bottom  << " table number :"<< table_number << " number of intersection " << intersect_number <<endl;*/

}


template <class Type>
void creat_parabola_table(std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< Type >>>>>> &parabola_table, const int &partition, const unsigned &m, const unsigned &n, const int &s){

  Type area ;
  Type a(0);
  std::vector <Type> pol1(3, 0);
  std::vector <Type> pol2(3, 0);
  unsigned int count;
  PointT <Type> p1, p2, p3 ;
//   cout << " Number | (x1,y1) | (x2,y2) | (x3,y3) |  k  |  b  |  d  |  c  | Area |" <<endl;
  Type del_x = 1/static_cast<Type>(partition);
  Type epsilon (0.001);
//   cout << "del_x " << del_x << endl;
  for (int table = 0; table <=7; table++){  // BEGIN preevaluation of the table
//     cout << "Table " << table << endl;
    parabola_table.resize(parabola_table.size() + 1);
//     parabola_table[table].resize(0);
    count = 0;

    for (unsigned int i1 = 0 ; i1<= partition ; i1++){

      parabola_table[table].resize(parabola_table[table].size() + 1);

      for (unsigned int i2=0 ;i2 <= partition ;i2++){

        parabola_table[table][i1].resize(parabola_table[table][i1].size() + 1);

        for (unsigned int i3=0 ;i3 <= partition ;i3++){

          parabola_table[table][i1][i2].resize(parabola_table[table][i1][i2].size() + 1);
//           cout << " i3 = " << i3 << endl;

          Type i1_pm_eps , i2_pm_eps;

           switch (table) {
                case 0:
                    i1_pm_eps = static_cast<Type>(i1*del_x - epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x + epsilon);
//                    if (i1 == partition ) i1_pm_eps = static_cast<Type>(i1*del_x - epsilon);     //it keeps my i2 in (0,1)
                    p1 = {static_cast<Type>(0), i1_pm_eps};
                    p2 = {i2_pm_eps, static_cast<Type>(1)};

                    break;
                case 1:
                    i1_pm_eps = static_cast<Type>(i1*del_x + epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x - epsilon);
                    p1 = {static_cast<Type>(0), i1_pm_eps};
                    p2 = {static_cast<Type>(1), i2_pm_eps};
                    break;
                case 2:
                    i1_pm_eps = static_cast<Type>(i1*del_x + epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x + epsilon);
                    p1 = {static_cast<Type>(0), i1_pm_eps};
                    p2 = {i2_pm_eps, static_cast<Type>(0)};
                    break;
                case 3:
                    i1_pm_eps = static_cast<Type>(i1*del_x - epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x + epsilon);
                    p1 = {i1_pm_eps,static_cast<Type>(1)};
                    p2 = {i2_pm_eps, static_cast<Type>(1)};
                    break;
                case 4:
                    i1_pm_eps = static_cast<Type>(i1*del_x - epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x - epsilon);
                    p1 = {i1_pm_eps, static_cast<Type>(1)};
                    p2 = {static_cast<Type>(1), i2_pm_eps};
                    break;
                case 5:
                    i1_pm_eps = static_cast<Type>(i1*del_x + epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x - epsilon);
                    p1 = {i1_pm_eps,static_cast<Type>(1)};
                    p2 = {i2_pm_eps, static_cast<Type>(0)};
                    break;
                case 6:
                    i1_pm_eps = static_cast<Type>(i1*del_x + epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x - epsilon);
                    p1 = {static_cast<Type>(1), i1_pm_eps};
                    p2 = {i2_pm_eps, static_cast<Type>(0)};
                    break;
                case 7:
                    i1_pm_eps = static_cast<Type>(i1*del_x - epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x + epsilon);
                    p1 = {i1_pm_eps, static_cast<Type>(0)};
                    p2 = {i2_pm_eps, static_cast<Type>(0)};
                    break;
           }

           p3 = {0.5 * (p1.x + p2.x) , i3*del_x };
           if (i3 == 0 ) p3.y = i3*del_x + epsilon;
           if (i3 == partition ) p3.y = i3*del_x - epsilon;



           p1 = {static_cast<Type>(p1.x), static_cast<Type>(p1.y)};
           p2 = {static_cast<Type>(p2.x), static_cast<Type>(p2.y)};
           p3 = {static_cast<Type>(p3.x), static_cast<Type>(p3.y)};


           Type c(1) ;
           Type det = p1.x * p1.x * (p2.x - p3.x) -p1.x* (p2.x*p2.x - p3.x*p3.x)+ p2.x*p3.x*(p2.x - p3.x) ;// only sort out the points parallel to y axis

            Parabola <Type> parabola;
            int intersect_number;
            std::vector <Type> intersection;
            std::vector <Type> interp_point;    //never used in this function. it was used in interpolation;
            unsigned int table_number = table;

            parabola = get_parabola_equation(p1, p2, p3);

            CheckIntersection <Type> (intersect_number, table_number, intersection, interp_point, parabola);

/*
            if (intersect_number > 2){
              bool convrt2line = 0;
              if (table_number == 5) convrt2line = 1;
              else if(pol2[0] < 0 ){
                if(table_number == 2 || table_number == 6) convrt2line = 1;
              }
              else if ( pol2[0] > 0 ){
                if(table_number == 0 || table_number == 4) convrt2line = 1;
              }
              if (convrt2line == 1){
                Type slope = (p3.y-p1.y)/(p3.x-p1.x);
                c=static_cast<Type>(1);
                pol2[0] = static_cast<Type>(0);    //k=0
                pol2[1] = -slope;
                pol2[2] = slope*p1.x - p1.y ;
                pol1[0] = pol2[0];    //k=0
                pol1[1] = pol2[1];
                pol1[2] = pol2[2] + c ;
//                 cout << pol2[0] << " " <<pol2[1] << " " << pol2[3] << " " << " " << c << endl;
              }
            }*/

          bool special_parabola = 0;
          int concaveUp = 0;
          if(parabola.k < 0) concaveUp = 1; //kx^2+bx+d+cy = 0 ; c is initialy 1 => y = -kx^2 - bx - d .
          else if (parabola.k > 0)concaveUp= -1;


            //TODO creat a new model special integral.

            if (intersect_number > 2 ){
              if (table == 4 && concaveUp == -1){ //TODO if it's concave up we are not going to do anything .. meaning if c and k has different sign.
                special_parabola = 1 ;
                const Type One(1) , Zero (0);

                Type mxArea = static_cast<Type>(1) / ((m + static_cast<Type>(1)) * (n + static_cast<Type>(1)));
                Type smallArea =  (1 - p1.x)*(1 - p2.y) ;
                Type restOfArea = mxArea * (1 - smallArea) ;

                // Mapping

                PointT <Type> q1 = {Zero, One};
                PointT <Type> q2 = {One, Zero};
                PointT <Type> q3 = {(p3.x - p1.x)/(One - p1.x) , (p3.y - p2.y)/(One - p2.y)};
                Parabola<Type> q_parabola = get_parabola_equation(q1, q2, q3);

                pol1[0] = q_parabola.k;
                pol1[1] = q_parabola.b;
                pol1[2] = q_parabola.d + c;
                pol2[0] = q_parabola.k;
                pol2[1] = q_parabola.b;
                pol2[2] = q_parabola.d;

                cout<<endl;

                cout<< "Q = (" << q1.x << ", " << q1.y << ") ," << "(" << q2.x << ", " << q2.y << ") ," << "(" << q3.x << ", " << q3.y << ")  : Q parabola = " << q_parabola.k <<"x^2 + "<< q_parabola.b <<"x + "<<q_parabola.d << endl;


                for(int normal=0; normal <=1; normal++ ){
                  Type A1 (0), A2 (0), A3 (0);

                  if (normal == 1){ // To calculate other side
                      pol1[0] *= -1;
                      pol1[1] *= -1;
                      pol1[2] *= -1;
                      pol2[0] *= -1;
                      pol2[1] *= -1;
                      pol2[2] *= -1;
                      c *= -1;
                  }


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
                  cout<< "special method : " << pol2[0] << "x^2+ " << pol2[1] << "x+ " << pol2[2] << "+ " << c << "y = 0    area ="<< area<<" ";
                  //Now scaling
                  area = area * smallArea ;
                  cout<< "scaling area: " << area<<" ";
                  if (normal == 1 ) area = area + restOfArea ;
                  cout<< "final area: " << area<<" ";


                  cout << " \n" <<"table : " << table << " " << count << ". (" << p1.x << ", " << p1.y << ") ," << "(" << p2.x << ", " << p2.y << ") ," << "(" << p3.x << ", " << p3.y << ")  : "  <<  parabola.k << "x^2+ " << parabola.b << "x+ " << parabola.d << "+ " << 1 << "y = 0     .area=" << area<<  endl;
                  //Just to keep all the values to positive:
                  if (p1.x < 0 ) p1.x = 0;
                  if (p1.y < 0 ) p1.y = 0;
                  if (p2.x < 0 ) p2.x = 0;
                  if (p2.y < 0 ) p2.y = 0;
                  if (p3.x < 0 ) p3.x = 0;
                  if (p3.y < 0 ) p3.y = 0;

                  parabola_table[table][i1][i2][i3].resize(2) ;
                  parabola_table[table][i1][i2][i3][normal].resize(15);

                  parabola_table[table][i1][i2][i3][normal][0] = count;
                  parabola_table[table][i1][i2][i3][normal][1] = i1_pm_eps;
                  parabola_table[table][i1][i2][i3][normal][2] = i2_pm_eps;
                  parabola_table[table][i1][i2][i3][normal][3] = static_cast<Type>(i3*del_x);

                  parabola_table[table][i1][i2][i3][normal][4] = p1.x;
                  parabola_table[table][i1][i2][i3][normal][5] = p1.y;
                  parabola_table[table][i1][i2][i3][normal][6] = p2.x;
                  parabola_table[table][i1][i2][i3][normal][7] = p2.y;
                  parabola_table[table][i1][i2][i3][normal][8] = p3.x;
                  parabola_table[table][i1][i2][i3][normal][9] = p3.y;

                  if(normal == 0){
                    parabola_table[table][i1][i2][i3][normal][10] = parabola.k;
                    parabola_table[table][i1][i2][i3][normal][11] = parabola.b;
                    parabola_table[table][i1][i2][i3][normal][12] = parabola.d;
                  }
                  else{
                    parabola_table[table][i1][i2][i3][normal][10] = -parabola.k;
                    parabola_table[table][i1][i2][i3][normal][11] = -parabola.b;
                    parabola_table[table][i1][i2][i3][normal][12] = -parabola.d;
                  }

                  parabola_table[table][i1][i2][i3][normal][13] = c;
                  parabola_table[table][i1][i2][i3][normal][14] = area;
                }
              }
              else if (table == 0 && concaveUp == -1){ //TODO if it's concave up we are not going to do anything .. meaning if c and k has different sign.
                special_parabola = 1 ;
                const Type One(1) , Zero (0);

                Type mxArea = static_cast<Type>(1) / ((m + static_cast<Type>(1)) * (n + static_cast<Type>(1)));
                Type smallArea =  (p2.x - p1.x) * (p2.y - p1.y) ;
                Type restOfArea = mxArea * (1 - smallArea) ;

                // Mapping

                PointT <Type> q1 = {Zero, Zero};
                PointT <Type> q2 = {One, One};
                PointT <Type> q3 = {(p3.x - p1.x)/(p2.x - p1.x) , (p3.y - p1.y)/(p2.y - p1.y)};
                Parabola<Type> q_parabola = get_parabola_equation(q1, q2, q3);

                pol1[0] = q_parabola.k;
                pol1[1] = q_parabola.b;
                pol1[2] = q_parabola.d + c;
                pol2[0] = q_parabola.k;
                pol2[1] = q_parabola.b;
                pol2[2] = q_parabola.d;

                cout<<endl;

                cout<< "Q = (" << q1.x << ", " << q1.y << ") ," << "(" << q2.x << ", " << q2.y << ") ," << "(" << q3.x << ", " << q3.y << ")  : Q parabola = " << q_parabola.k <<"x^2 + "<< q_parabola.b <<"x + "<<q_parabola.d << endl;


                for(int normal=0; normal <=1; normal++ ){
                  Type A1 (0), A2 (0), A3 (0);

                  if (normal == 1){ // To calculate other side
                      pol1[0] *= -1;
                      pol1[1] *= -1;
                      pol1[2] *= -1;
                      pol2[0] *= -1;
                      pol2[1] *= -1;
                      pol2[2] *= -1;
                      c *= -1;
                  }


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
                  cout<< "special method : " << pol2[0] << "x^2+ " << pol2[1] << "x+ " << pol2[2] << "+ " << c << "y = 0    area ="<< area<<" ";
                  //Now scaling
                  area = area * smallArea ;
                  cout<< "scaling area: " << area<<" ";
                  if (normal == 1 ) area = area + restOfArea ;
                  cout<< "final area: " << area<<" ";


                  cout << " \n" <<"table : " << table << " " << count << ". (" << p1.x << ", " << p1.y << ") ," << "(" << p2.x << ", " << p2.y << ") ," << "(" << p3.x << ", " << p3.y << ")  : "  <<  parabola.k << "x^2+ " << parabola.b << "x+ " << parabola.d << "+ " << 1 << "y = 0     .area=" << area<<  endl;

                                    //Just to keep all the values to positive:
                  if (p1.x < 0 ) p1.x = 0;
                  if (p1.y < 0 ) p1.y = 0;
                  if (p2.x < 0 ) p2.x = 0;
                  if (p2.y < 0 ) p2.y = 0;
                  if (p3.x < 0 ) p3.x = 0;
                  if (p3.y < 0 ) p3.y = 0;

                  parabola_table[table][i1][i2][i3].resize(2) ;
                  parabola_table[table][i1][i2][i3][normal].resize(15);

                  parabola_table[table][i1][i2][i3][normal][0] = count;
                  parabola_table[table][i1][i2][i3][normal][1] = i1_pm_eps;
                  parabola_table[table][i1][i2][i3][normal][2] = i2_pm_eps;
                  parabola_table[table][i1][i2][i3][normal][3] = static_cast<Type>(i3*del_x);

                  parabola_table[table][i1][i2][i3][normal][4] = p1.x;
                  parabola_table[table][i1][i2][i3][normal][5] = p1.y;
                  parabola_table[table][i1][i2][i3][normal][6] = p2.x;
                  parabola_table[table][i1][i2][i3][normal][7] = p2.y;
                  parabola_table[table][i1][i2][i3][normal][8] = p3.x;
                  parabola_table[table][i1][i2][i3][normal][9] = p3.y;

                  if(normal == 0){
                    parabola_table[table][i1][i2][i3][normal][10] = parabola.k;
                    parabola_table[table][i1][i2][i3][normal][11] = parabola.b;
                    parabola_table[table][i1][i2][i3][normal][12] = parabola.d;
                  }
                  else{
                    parabola_table[table][i1][i2][i3][normal][10] = -parabola.k;
                    parabola_table[table][i1][i2][i3][normal][11] = -parabola.b;
                    parabola_table[table][i1][i2][i3][normal][12] = -parabola.d;
                  }

                  parabola_table[table][i1][i2][i3][normal][13] = c;
                  parabola_table[table][i1][i2][i3][normal][14] = area;
                }
              }
              else if (table == 2 && concaveUp == 1){ //TODO if it's concave up we are not going to do anything .. meaning if c and k has different sign.
                special_parabola = 1 ;
                const Type One(1) , Zero (0);

                Type mxArea = static_cast<Type>(1) / ((m + static_cast<Type>(1)) * (n + static_cast<Type>(1)));
                Type smallArea =  (p2.x - p1.x) * (p1.y - p2.y) ;
                Type restOfArea = mxArea * (1 - smallArea) ;

                // Mapping

                PointT <Type> q1 = {Zero, One};
                PointT <Type> q2 = {One, Zero};
                PointT <Type> q3 = {(p3.x - p1.x)/(p2.x - p1.x) , (p3.y - p2.y)/(p1.y - p2.y)};
                Parabola<Type> q_parabola = get_parabola_equation(q1, q2, q3);

                pol1[0] = q_parabola.k;
                pol1[1] = q_parabola.b;
                pol1[2] = q_parabola.d + c;
                pol2[0] = q_parabola.k;
                pol2[1] = q_parabola.b;
                pol2[2] = q_parabola.d;

                cout<<endl;

                cout<< "Q = (" << q1.x << ", " << q1.y << ") ," << "(" << q2.x << ", " << q2.y << ") ," << "(" << q3.x << ", " << q3.y << ")  : Q parabola = " << q_parabola.k <<"x^2 + "<< q_parabola.b <<"x + "<<q_parabola.d << endl;


                for(int normal=0; normal <=1; normal++ ){
                  Type A1 (0), A2 (0), A3 (0);

                  if (normal == 1){ // To calculate other side
                      pol1[0] *= -1;
                      pol1[1] *= -1;
                      pol1[2] *= -1;
                      pol2[0] *= -1;
                      pol2[1] *= -1;
                      pol2[2] *= -1;
                      c *= -1;
                  }


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
                  cout<< "special method : " << pol2[0] << "x^2+ " << pol2[1] << "x+ " << pol2[2] << "+ " << c << "y = 0    area =" << area << " ";
                  //Now scaling
                  area = area * smallArea ;
                  cout<< "scaling area: " << area<<" ";
                  if (normal == 0 ) area = area + restOfArea ;
                  cout<< "final area: " << area<<" ";


                  cout << " \n" <<"table : " << table << " " << count << ". (" << p1.x << ", " << p1.y << ") ," << "(" << p2.x << ", " << p2.y << ") ," << "(" << p3.x << ", " << p3.y << ")  : "  <<  parabola.k << "x^2+ " << parabola.b << "x+ " << parabola.d << "+ " << 1 << "y = 0     .area=" << area<<  endl;

                  //Just to keep all the values to positive:
                  if (p1.x < 0 ) p1.x = 0;
                  if (p1.y < 0 ) p1.y = 0;
                  if (p2.x < 0 ) p2.x = 0;
                  if (p2.y < 0 ) p2.y = 0;
                  if (p3.x < 0 ) p3.x = 0;
                  if (p3.y < 0 ) p3.y = 0;


                  parabola_table[table][i1][i2][i3].resize(2) ;
                  parabola_table[table][i1][i2][i3][normal].resize(15);

                  parabola_table[table][i1][i2][i3][normal][0] = count;
                  parabola_table[table][i1][i2][i3][normal][1] = i1_pm_eps;
                  parabola_table[table][i1][i2][i3][normal][2] = i2_pm_eps;
                  parabola_table[table][i1][i2][i3][normal][3] = static_cast<Type>(i3*del_x);

                  parabola_table[table][i1][i2][i3][normal][4] = p1.x;
                  parabola_table[table][i1][i2][i3][normal][5] = p1.y;
                  parabola_table[table][i1][i2][i3][normal][6] = p2.x;
                  parabola_table[table][i1][i2][i3][normal][7] = p2.y;
                  parabola_table[table][i1][i2][i3][normal][8] = p3.x;
                  parabola_table[table][i1][i2][i3][normal][9] = p3.y;

                  if(normal == 0){
                    parabola_table[table][i1][i2][i3][normal][10] = parabola.k;
                    parabola_table[table][i1][i2][i3][normal][11] = parabola.b;
                    parabola_table[table][i1][i2][i3][normal][12] = parabola.d;
                  }
                  else{
                    parabola_table[table][i1][i2][i3][normal][10] = -parabola.k;
                    parabola_table[table][i1][i2][i3][normal][11] = -parabola.b;
                    parabola_table[table][i1][i2][i3][normal][12] = -parabola.d;
                  }

                  parabola_table[table][i1][i2][i3][normal][13] = c;
                  parabola_table[table][i1][i2][i3][normal][14] = area;
                }
              }

              else if (table == 6 && concaveUp == 1){ //TODO if it's concave up we are not going to do anything .. meaning if c and k has different sign.
                special_parabola = 1 ;
                const Type One(1) , Zero (0);

                Type mxArea = static_cast<Type>(1) / ((m + static_cast<Type>(1)) * (n + static_cast<Type>(1)));
                Type smallArea =  (p1.x - p2.x) * (p1.y - p2.y) ;
                Type restOfArea = mxArea * (1 - smallArea) ;

                // Mapping

                PointT <Type> q1 = {One, One};
                PointT <Type> q2 = {Zero, Zero};
                PointT <Type> q3 = {(p3.x - p2.x)/(p1.x - p2.x) , (p3.y - p2.y)/(p1.y - p2.y)};
                Parabola<Type> q_parabola = get_parabola_equation(q1, q2, q3);

                pol1[0] = q_parabola.k;
                pol1[1] = q_parabola.b;
                pol1[2] = q_parabola.d + c;
                pol2[0] = q_parabola.k;
                pol2[1] = q_parabola.b;
                pol2[2] = q_parabola.d;

                cout<<endl;

                cout<< "Q = (" << q1.x << ", " << q1.y << ") ," << "(" << q2.x << ", " << q2.y << ") ," << "(" << q3.x << ", " << q3.y << ")  : Q parabola = " << q_parabola.k <<"x^2 + "<< q_parabola.b <<"x + "<<q_parabola.d << endl;


                for(int normal=0; normal <=1; normal++ ){
                  Type A1 (0), A2 (0), A3 (0);

                  if (normal == 1){ // To calculate other side
                      pol1[0] *= -1;
                      pol1[1] *= -1;
                      pol1[2] *= -1;
                      pol2[0] *= -1;
                      pol2[1] *= -1;
                      pol2[2] *= -1;
                      c *= -1;
                  }


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
                  cout<< "special method : " << pol2[0] << "x^2+ " << pol2[1] << "x+ " << pol2[2] << "+ " << c << "y = 0    area =" << area << " ";
                  //Now scaling
                  area = area * smallArea ;
                  cout<< "scaling area: " << area<<" ";

                  if (normal == 0 ) area = area + restOfArea ;
                  cout<< "final area: " << area<<" ";


                  cout << " \n" <<"table : " << table << " " << count << ". (" << p1.x << ", " << p1.y << ") ," << "(" << p2.x << ", " << p2.y << ") ," << "(" << p3.x << ", " << p3.y << ")  : "  <<  parabola.k << "x^2+ " << parabola.b << "x+ " << parabola.d << "+ " << 1 << "y = 0     .area=" << area<<  endl;

                  //Just to keep all the values to positive:
                  if (p1.x < 0 ) p1.x = 0;
                  if (p1.y < 0 ) p1.y = 0;
                  if (p2.x < 0 ) p2.x = 0;
                  if (p2.y < 0 ) p2.y = 0;
                  if (p3.x < 0 ) p3.x = 0;
                  if (p3.y < 0 ) p3.y = 0;


                  parabola_table[table][i1][i2][i3].resize(2) ;
                  parabola_table[table][i1][i2][i3][normal].resize(15);

                  parabola_table[table][i1][i2][i3][normal][0] = count;
                  parabola_table[table][i1][i2][i3][normal][1] = i1_pm_eps;
                  parabola_table[table][i1][i2][i3][normal][2] = i2_pm_eps;
                  parabola_table[table][i1][i2][i3][normal][3] = static_cast<Type>(i3*del_x);

                  parabola_table[table][i1][i2][i3][normal][4] = p1.x;
                  parabola_table[table][i1][i2][i3][normal][5] = p1.y;
                  parabola_table[table][i1][i2][i3][normal][6] = p2.x;
                  parabola_table[table][i1][i2][i3][normal][7] = p2.y;
                  parabola_table[table][i1][i2][i3][normal][8] = p3.x;
                  parabola_table[table][i1][i2][i3][normal][9] = p3.y;

                  if(normal == 0){
                    parabola_table[table][i1][i2][i3][normal][10] = parabola.k;
                    parabola_table[table][i1][i2][i3][normal][11] = parabola.b;
                    parabola_table[table][i1][i2][i3][normal][12] = parabola.d;
                  }
                  else{
                    parabola_table[table][i1][i2][i3][normal][10] = -parabola.k;
                    parabola_table[table][i1][i2][i3][normal][11] = -parabola.b;
                    parabola_table[table][i1][i2][i3][normal][12] = -parabola.d;
                  }

                  parabola_table[table][i1][i2][i3][normal][13] = c;
                  parabola_table[table][i1][i2][i3][normal][14] = area;
                }
              }

            }

            //Using line instead of 4 intersection
            //TODO It is still linearizing concave down table 4 ; fix it

            if(!special_parabola){

              if (intersect_number == 2 || table == 1 || table == 3 || table == 7){ // we are not changing 4 intersection when do table-1 left-right
                  if (det !=0){
                      pol1[0] = parabola.k;
                      pol1[1] = parabola.b;
                      pol1[2] = parabola.d + c;
                      pol2[0] = parabola.k;
                      pol2[1] = parabola.b;
                      pol2[2] = parabola.d;
                  }
                  else {  //TODO decode what I did here. Couldn't figure it out what I did here. We are not using this probably after introducing epsilon
                      cout << " STOP : FIGURE OUT WHAT YOU DID : " << endl;
                      pol1[0] = static_cast<Type>(0);
                      pol1[1] = static_cast<Type>(-1);
                      pol1[2] = p1.x;
                      pol2[0] = pol1[0];
                      pol2[1] = pol1[1];
                      pol2[2] = p1.x;
                      c=static_cast<Type>(0);
  //                     cout << "went in straight line : " <<pol2[0] << " " <<pol2[1] << " " << pol2[3] << " " << " " << c << endl;
                  }
              }

//               else if (intersect_number > 2){
//                   Type slope = (p3.y-p1.y)/(p3.x-p1.x);
//                   c=static_cast<Type>(1);
//                   pol2[0] = static_cast<Type>(0);    //k=0
//                   pol2[1] = -slope;
//                   pol2[2] = slope*p1.x - p1.y ;
//                   pol1[0] = pol2[0];    //k=0
//                   pol1[1] = pol2[1];
//                   pol1[2] = pol2[2] + c ;
//   //                 cout << pol2[0] << " " <<pol2[1] << " " << pol2[3] << " " << " " << c << endl;
//                   cout<< " intersection << . used a straight line " << endl;
//               }
              for(int normal=0; normal <=1; normal++ ){
                  Type A1 (0), A2 (0), A3 (0);

                  if (normal == 1){ // To calculate other side
                      pol1[0] *= -1;
                      pol1[1] *= -1;
                      pol1[2] *= -1;
                      pol2[0] *= -1;
                      pol2[1] *= -1;
                      pol2[2] *= -1;
                      c *= -1;
                  }


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

                  cout << " \n" <<"table : " << table << " " << count << ". (" << p1.x << ", " << p1.y << ") ," << "(" << p2.x << ", " << p2.y << ") ," << "(" << p3.x << ", " << p3.y << ")  : "  <<  pol2[0] << "x^2+ " << pol2[1] << "x+ " << pol2[2] << "+ " << c << "y = 0     .area=" << area<<  endl;


                                    //Just to keep all the values to positive:
                  if (p1.x < 0 ) p1.x = 0;
                  if (p1.y < 0 ) p1.y = 0;
                  if (p2.x < 0 ) p2.x = 0;
                  if (p2.y < 0 ) p2.y = 0;
                  if (p3.x < 0 ) p3.x = 0;
                  if (p3.y < 0 ) p3.y = 0;


                  parabola_table[table][i1][i2][i3].resize(2) ;
                  parabola_table[table][i1][i2][i3][normal].resize(15);

                  parabola_table[table][i1][i2][i3][normal][0] = count;
                  parabola_table[table][i1][i2][i3][normal][1] = i1_pm_eps;
                  parabola_table[table][i1][i2][i3][normal][2] = i2_pm_eps;
                  parabola_table[table][i1][i2][i3][normal][3] = static_cast<Type>(i3*del_x);

                  parabola_table[table][i1][i2][i3][normal][4] = p1.x;
                  parabola_table[table][i1][i2][i3][normal][5] = p1.y;
                  parabola_table[table][i1][i2][i3][normal][6] = p2.x;
                  parabola_table[table][i1][i2][i3][normal][7] = p2.y;
                  parabola_table[table][i1][i2][i3][normal][8] = p3.x;
                  parabola_table[table][i1][i2][i3][normal][9] = p3.y;

                  parabola_table[table][i1][i2][i3][normal][10] = pol2[0];
                  parabola_table[table][i1][i2][i3][normal][11] = pol2[1];
                  parabola_table[table][i1][i2][i3][normal][12] = pol2[2];
                  parabola_table[table][i1][i2][i3][normal][13] = c;
                  parabola_table[table][i1][i2][i3][normal][14] = area;

  //                 for (int nm=0; nm<=11; nm++){
  //                 cout << parabola_table[normal][count][nm] << " " ;
  //                 }
              }
            }
            count ++ ;
        }
      }
    }
  }

}


template <class Type>
void creat_parabola_table_4intersection(std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< Type >>>>>> &parabola_table_4intersection , const int &partition, const unsigned &m, const unsigned &n, const int &s){

  Type area ;
  Type a(0);
  std::vector <Type> pol1(3, 0);
  std::vector <Type> pol2(3, 0);
  unsigned int count;

  PointT <Type> p1, p2, p3 ;
//   cout << " Number | (x1,y1) | (x2,y2) | (x3,y3) |  k  |  b  |  d  |  c  | Area |" <<endl;
  Type del_x = static_cast<Type>(1)/partition;
  Type epsilon (0.001);
//   cout << "del_x " << del_x << endl;
  for (int table = 0; table <=7; table++){  // BEGIN preevaluation of the table
//     cout << "Table " << table << endl;
    parabola_table_4intersection.resize(parabola_table_4intersection.size() + 1);
//     parabola_table_4intersection[table].resize(0);
    count = 0;

    for (unsigned int i1 = 0 ; i1<= partition ; i1++){
      parabola_table_4intersection[table].resize(parabola_table_4intersection[table].size() + 1);
      for (unsigned int i2=0 ;i2 <= partition ;i2++){
        parabola_table_4intersection[table][i1].resize(parabola_table_4intersection[table][i1].size() + 1);
        for (unsigned int i3=0 ;i3 <= partition ;i3++){
          parabola_table_4intersection[table][i1][i2].resize(parabola_table_4intersection[table][i1][i2].size() + 1);
//           cout << " i3 = " << i3 << endl;

           Type i1_pm_eps , i2_pm_eps, i3_pm_eps;

           switch (table) {
                case 0:  // Left-Top-Top - right
                    i1_pm_eps = static_cast<Type>(i1*del_x - epsilon);  // this creates all parabola with same normal
                    i2_pm_eps = static_cast<Type>(i2*del_x + epsilon);
                    i3_pm_eps = static_cast<Type>(i3*del_x - epsilon);
                      if (i2 == partition ) i2_pm_eps = static_cast<Type>(i2*del_x - epsilon); //it keeps my i2 in (0,1)

                    p1 = {static_cast<Type>(0), i1_pm_eps};
                    p2 = {i2_pm_eps, static_cast<Type>(1)};
                    p3 = {static_cast<Type>(1),i3_pm_eps};
                    break;

                case 1:   // Left-Top-Top - bottom
                    i1_pm_eps = static_cast<Type>(i1*del_x - epsilon);  //why do I need to use epslon here?
                    i2_pm_eps = static_cast<Type>(i2*del_x + epsilon);
                    i3_pm_eps = static_cast<Type>(i3*del_x + 2*epsilon);
                    p1 = {static_cast<Type>(0), i1_pm_eps};
                    p2 = {i2_pm_eps, static_cast<Type>(1)};
                    p3 = {i3_pm_eps, static_cast<Type>(0)};
                    break;


                case 2:   // Left - Bottom - Bottom - right
                    i1_pm_eps = static_cast<Type>(i1*del_x + epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x + epsilon);
                    i3_pm_eps = static_cast<Type>(i3*del_x + epsilon);
                    if (i2 == partition ) i2_pm_eps = static_cast<Type>(i2*del_x - epsilon);

                    p1 = {static_cast<Type>(0), i1_pm_eps};
                    p2 = {i2_pm_eps, static_cast<Type>(0)};
                    p3 = {static_cast<Type>(1), i3_pm_eps};
                    break;


                case 3:   // Left - Bottom - Bottom - top
                    i1_pm_eps = static_cast<Type>(i1*del_x + epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x + epsilon);
                    i3_pm_eps = static_cast<Type>(i3*del_x + 2*epsilon);
                    p1 = {static_cast<Type>(0), i1_pm_eps};
                    p2 = {i2_pm_eps, static_cast<Type>(0)};
                    p3 = {i3_pm_eps, static_cast<Type>(1)};
                    break;

                case 4:   // Right- Top - Bottom
                    i1_pm_eps = static_cast<Type>(i1*del_x - epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x - epsilon);
                    i3_pm_eps = static_cast<Type>(i3*del_x - 2*epsilon);
                    p1 = {static_cast<Type>(1), i1_pm_eps};
                    p2 = {i2_pm_eps, static_cast<Type>(1)};
                    p3 = {i3_pm_eps, static_cast<Type>(0)};
                    break;


                case 5:  // Right - Bottom - Top
                    i1_pm_eps = static_cast<Type>(i1*del_x + epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x - epsilon);
                    i3_pm_eps = static_cast<Type>(i3*del_x - 2*epsilon);
                    p1 = {static_cast<Type>(1), i1_pm_eps};
                    p2 = {i2_pm_eps, static_cast<Type>(0)};
                    p3 = {i3_pm_eps, static_cast<Type>(1)};
                    break;


                case 6:   //  Bottom - Top - Top - Bottom
                    i1_pm_eps = static_cast<Type>(i1*del_x - epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x + epsilon);
                    i3_pm_eps = static_cast<Type>(i3*del_x + 2*epsilon);
                    p1 = {i1_pm_eps, static_cast<Type>(0)};
                    p2 = {i2_pm_eps, static_cast<Type>(1)};
                    p3 = {i3_pm_eps, static_cast<Type>(0) };
                    break;


                case 7:  // Top -Bottom -Bottom -Top
                    i1_pm_eps = static_cast<Type>(i1*del_x - epsilon);
                    i2_pm_eps = static_cast<Type>(i2*del_x + epsilon);
                    i3_pm_eps = static_cast<Type>(i3*del_x + 2*epsilon);
                    p1 = {i1_pm_eps, static_cast<Type>(1)};
                    p2 = {i2_pm_eps, static_cast<Type>(0)};
                    p3 = {i3_pm_eps, static_cast<Type>(1)};
                    break;
           }
/*
           p1 = {static_cast<Type>(p1.x), static_cast<Type>(p1.y)};
           p2 = {static_cast<Type>(p2.x), static_cast<Type>(p2.y)};
           p3 = {static_cast<Type>(p3.x), static_cast<Type>(p3.y)};*/
           Type c(1) ;
           Type det = p1.x * p1.x * (p2.x - p3.x) -p1.x* (p2.x*p2.x - p3.x*p3.x)+ p2.x*p3.x*(p2.x - p3.x) ;// only sort out the points parallel to y axis

            Parabola <Type> parabola ;
            int intersect_number;
            std::vector <Type> intersection;
            std::vector <Type> interp_point ; //never used in this function. it was used in interpolation;
            unsigned int table_number = table ;

            parabola = get_parabola_equation(p1, p2, p3);

            CheckIntersection <Type> (intersect_number, table_number , intersection, interp_point , parabola);


//             cout << "intersect_number = " << intersect_number << endl;
//             cout << "point intersected = " << intersection.size() ;
//             for (int j =0;j < intersection.size() ; j++){
//              cout <<" " << intersection[j]  ;
//             }



                if (det !=0){
                    pol1[0] = parabola.k;
                    pol1[1] = parabola.b;
                    pol1[2] = parabola.d + c;
                    pol2[0] = parabola.k;
                    pol2[1] = parabola.b;
                    pol2[2] = parabola.d;
                }
                else cout << " determinant 0 . we have to use  previous value : "<< endl;

            if (intersect_number < 4) cout << " ----------- determinant = " << det <<
              " intersection number is something other than 4 check please " << " =/////////////////============///////////////////=======================/////////////////////===============////////////// intersection point = "  << intersect_number<< endl;
//             else cout << " It intersect the unit box in numbers other than 2 or 4 " ;

            for(int normal=0; normal <=1; normal++ ){
                Type A1 (0), A2 (0), A3 (0);

                if (normal == 1){ // To calculate other side
                    pol1[0] *= -1;
                    pol1[1] *= -1;
                    pol1[2] *= -1;
                    pol2[0] *= -1;
                    pol2[1] *= -1;
                    pol2[2] *= -1;
                    c *= -1;
                }


                std::vector< std::pair <Type, Type> > I1(0), I2(0), I3(0) ;
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

//                 cout << " \n" <<"table : " << table << " " << count << ". (" << p1.x << ", " << p1.y << ") ," << "(" << p2.x << ", " << p2.y << ") ," << "(" << p3.x << ", " << p3.y << ")  : "  <<  pol2[0] << ", " << pol2[1] << ", " << pol2[2] << ", " << c << ", " << area<<  endl;

//                 parabola_table_4intersection[table].resize(parabola_table_4intersection[table].size() + 1);
//                 parabola_table_4intersection[table][count].resize(2);

                parabola_table_4intersection[table][i1][i2][i3].resize(2) ;
                parabola_table_4intersection[table][i1][i2][i3][normal].resize(15);

                parabola_table_4intersection[table][i1][i2][i3][normal][0] = count;
                parabola_table_4intersection[table][i1][i2][i3][normal][1] = i1_pm_eps;
                parabola_table_4intersection[table][i1][i2][i3][normal][2] = i2_pm_eps;
                parabola_table_4intersection[table][i1][i2][i3][normal][3] = i3_pm_eps;

                parabola_table_4intersection[table][i1][i2][i3][normal][4] = p1.x;
                parabola_table_4intersection[table][i1][i2][i3][normal][5] = p1.y;
                parabola_table_4intersection[table][i1][i2][i3][normal][6] = p2.x;
                parabola_table_4intersection[table][i1][i2][i3][normal][7] = p2.y;
                parabola_table_4intersection[table][i1][i2][i3][normal][8] = p3.x;
                parabola_table_4intersection[table][i1][i2][i3][normal][9] = p3.y;

                parabola_table_4intersection[table][i1][i2][i3][normal][10] = pol2[0];
                parabola_table_4intersection[table][i1][i2][i3][normal][11] = pol2[1];
                parabola_table_4intersection[table][i1][i2][i3][normal][12] = pol2[2];
                parabola_table_4intersection[table][i1][i2][i3][normal][13] = c;
                parabola_table_4intersection[table][i1][i2][i3][normal][14] = area;


                if(det == 0 && count != 0){      // TODO check if I still using this.
                 cout<< " using previous value : " <<  parabola_table_4intersection[table][i1][i2][i3-1][normal][14] << "instead of using the formula area : " << area << endl;
                 parabola_table_4intersection[table][i1][i2][i3][normal][14] = parabola_table_4intersection[table][i1][i2][i3-1][normal][14] ;
                }


//                 cout << " \n" <<"table : " << table << " " << count << " " <<parabola_table_4intersection[table][count][normal][1]  << " " << parabola_table_4intersection[table][count][normal][2] << " " << parabola_table_4intersection[table][count][normal][3] << " ==> " << parabola_table_4intersection[table][count][normal][10]  << " " << parabola_table_4intersection[table][count][normal][11] << " " << parabola_table_4intersection[table][count][normal][12]<<" " << parabola_table_4intersection[table][count][normal][13]  << " " << parabola_table_4intersection[table][count][normal][14] << endl;
//                 cout << parabola_table_4intersection[normal][count][nm] << " " ;
            }
            count ++ ;
//           }         //det ==0 ends
        }
      }
    }
  }

}



template <class Type>
void inverse_parabola(std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< Type >>>>>> &parabola_table, std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< Type >>>>>> &parabola_table_4intersection, const std::vector <Type> &given_parabola, vector <Type> &intersection, int &normal, std::vector< Type > &interp_point, std::vector< std::vector< Type >> & interp_table, unsigned int &table_number, const int &partition){

    int intersect_number(0);
    interp_point.resize(0);
    intersection.resize(0);

    Type k = given_parabola[0];
    Type b = given_parabola[1];
    Type d = given_parabola[2];
    Type c = given_parabola[3];

    cout << " original parabola " << k <<"x^2 + "<< b <<"x + "<< d << " + " << c <<"y = 0" << endl;

    bool left=0 , top = 0 , right = 0 , bottom = 0;

    cout << "normal before = " << normal ;
    if (c<0) normal = (normal+1) % 2 ;
    cout << " normal after = " << normal <<endl;
    k = k/c; b=b/c; d=d/c; c=static_cast<Type>(1);  //TODO 1 4/24 sign of c is different.



//     else if (normal == 1){ k = -k/c; b=-b/c; d=-d/c; c=-1;}

    cout << " simplified parabola " << k <<"x^2 + "<< b <<"x + "<< d << " + " << c <<"y = 0" << endl;


    {
//       if (-d>=0 && -d<=1){ //LEFT
//         intersection.resize(intersection.size()+2);
//         intersection[intersection.size()-2] = 0;
//         intersection[intersection.size()-1] = -d;
//         left = 1 ;
//         interp_point.resize(interp_point.size()+1);
//         interp_point[interp_point.size()-1] = -d;
//       }
//       // LEFT-TOP solve kx^2+bx+d-1 =0  ; Table 0
//       if (k == 0){
//       Type x =  (-1-d)/b ;
//         if(x < 1 && x> 0) {
//           intersection.resize(intersection.size()+2);
//           intersection[intersection.size()-2] = x;
//           intersection[intersection.size()-1] = 1;
//           interp_point.resize(interp_point.size()+1);
//           interp_point[interp_point.size()-1] = x;
//           top =1;
//           if (left ==1) table_number = 0 ;
//         }
//       }
//       else {
//         Type delta = b*b - 4*k*(d+1);
//         if(delta >0) {
//           Type sqrtdelta = sqrt(delta);
//           int sign = (k > 0) ? 1 : -1;
//
//           for(unsigned i = 0; i < 2; i++) {
//             Type x = (- b - sign * sqrtdelta) / (2 * k);
// //             cout<< "Top x = "<< x<< endl;
//             if(x > 1) break;
//             else if(x > 0) {
//               intersection.resize(intersection.size()+2);
//               intersection[intersection.size()-2] = x;
//               intersection[intersection.size()-1] = 1;
//               interp_point.resize(interp_point.size()+1);
//               interp_point[interp_point.size()-1] = x;
//               if (top ==1){table_number = 3 ;}
//               top = 1;
//               if (left ==1){table_number = 0 ;break;}
//
//             }
//             sign *= -1;
//           }
//         }
//       }
//       Type y_1=-(k+b+d); //LEFT-RIGHT x=1 ; Table 1
//       if (y_1 >= 0 && y_1 <= 1){ //TODO check sign when normal change
//           intersection.resize(intersection.size()+2);
//           intersection[intersection.size()-2] = 1;
//           intersection[intersection.size()-1] = y_1;
//           interp_point.resize(interp_point.size()+1);
//           interp_point[interp_point.size()-1] = y_1;
//           if (left ==1){table_number = 1 ;}
//           if  (top ==1){table_number = 4 ;}
//           right = 1 ;
//       }
//           //LEFT-BOTTOM  solve kx^2+bx+d =0 ; Table 2
//         if (k == 0){
//           Type x =  -d/b ;
//           if(x < 1 && x> 0) {
//             intersection.resize(intersection.size()+2);
//             intersection[intersection.size()-2] = x;
//             intersection[intersection.size()-1] = 0;
//             interp_point.resize(interp_point.size()+1);
//             interp_point[interp_point.size()-1] = x;
//             if (left ==1){table_number = 2 ;}
//             if (right ==1){table_number = 6 ;}
//             if (top ==1){table_number = 5 ;}
//           }
//         }
//
//         else {
//           Type delta = b*b - 4*k*d;
//           if(delta >0) {
//             Type sqrtdelta = sqrt(delta);
//             int sign = (k > 0) ? 1 : -1;
//
//             for(unsigned i = 0; i < 2; i++) {
//               Type x = (- b - sign * sqrtdelta) / (2 * k);
//               if(x > 1) break;
//               else if(x > 0) {
//                 intersection.resize(intersection.size()+2);
//                 intersection[intersection.size()-2] = x;
//                 intersection[intersection.size()-1] = 0;
//                 interp_point.resize(interp_point.size()+1);
//                 interp_point[interp_point.size()-1] = x;
//                 if (bottom ==1){table_number = 7 ;}
//                 if (left ==1){table_number = 2 ;}
//                 if (right ==1){table_number = 6 ;}
//                 if (top ==1){table_number = 5 ;}    // TODO check the table
//                 bottom = 1;
//               }
//               sign *= -1;
//             }
//           }
//         }
    }

      Parabola <Type> parabola{k,b,d};
      CheckIntersection <Type> (intersect_number, table_number , intersection, interp_point , parabola);


      if (interp_point.size() == 2){  // finding p3
            Type mid_point = 0.5*(intersection[0]+intersection[2]);
            intersection.resize(intersection.size()+2);
            intersection[intersection.size()-2] = mid_point;
            intersection[intersection.size()-1] = - (k*mid_point*mid_point + b*mid_point + d);


            interp_point.resize(3);
            interp_point[2] = -(k*mid_point*mid_point + b*mid_point + d);
            if (interp_point[2]<0 && interp_point[2]>1) cout << " parabola at midpoint is outside : check your intersections" << endl;
      }

      cout<< " inter section line = " << " left " << left << " top "<< top << " right "<< right << " bottom " << bottom  << " table number :"<< table_number << " number of intersection " << intersect_number <<endl;

      //finding closed points
      if (interp_point.size() >= 3){
        unsigned i1_0,i1_1,i2_0,i2_1,i3_0,i3_1 ;
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



      if(intersect_number == 2){
        cout << " Table used - 2 intersection " << endl;


//         interp_point[0] = static_cast<double>(interp_point[0]);

        i1_0 = floor(static_cast<double>(interp_point[0] * partition)) ;
        i1_1 = ceil( static_cast<double>(interp_point[0] * partition)) ;
        i2_0 = floor(static_cast<double>(interp_point[1] * partition)) ;
        i2_1 = ceil( static_cast<double>(interp_point[1] * partition)) ;
        i3_0 = floor(static_cast<double>(interp_point[2] * partition)) ;
        i3_1 = ceil( static_cast<double>(interp_point[2] * partition)) ;



        interp_table[0][3] = parabola_table[table_number][i1_0][i2_0][i3_0][normal][14] ;
        interp_table[1][3] = parabola_table[table_number][i1_0][i2_0][i3_1][normal][14] ;
        interp_table[2][3] = parabola_table[table_number][i1_0][i2_1][i3_0][normal][14] ;
        interp_table[3][3] = parabola_table[table_number][i1_0][i2_1][i3_1][normal][14] ;
        interp_table[4][3] = parabola_table[table_number][i1_1][i2_0][i3_0][normal][14] ;
        interp_table[5][3] = parabola_table[table_number][i1_1][i2_0][i3_1][normal][14] ;
        interp_table[6][3] = parabola_table[table_number][i1_1][i2_1][i3_0][normal][14] ;
        interp_table[7][3] = parabola_table[table_number][i1_1][i2_1][i3_1][normal][14] ;


        cout<< " given parabola = " << k << " " << b << " "<< d << " "<< c <<endl;
//         cout<< " intersection points = " << interp_point[0] << " " << interp_point[1] << " "<< interp_point[2] <<endl;
        cout << " interpolation table = " << endl;
        for(unsigned int i = 0; i <=7 ; i++){
          cout << interp_table[i][0] << " " <<interp_table[i][1] << " " <<interp_table[i][2] << " " <<interp_table[i][3]<< endl;
        }

        interp_table[0] = {parabola_table[table_number][i1_0][i2_0][i3_0][normal][1],parabola_table[table_number][i1_0][i2_0][i3_0][normal][2],parabola_table[table_number][i1_0][i2_0][i3_0][normal][3],parabola_table[table_number][i1_0][i2_0][i3_0][normal][14] };
        interp_table[1] = {parabola_table[table_number][i1_0][i2_0][i3_1][normal][1],parabola_table[table_number][i1_0][i2_0][i3_1][normal][2],parabola_table[table_number][i1_0][i2_0][i3_1][normal][3],parabola_table[table_number][i1_0][i2_0][i3_1][normal][14]};
        interp_table[2] = {parabola_table[table_number][i1_0][i2_1][i3_0][normal][1],parabola_table[table_number][i1_0][i2_1][i3_0][normal][2],parabola_table[table_number][i1_0][i2_1][i3_0][normal][3],parabola_table[table_number][i1_0][i2_1][i3_0][normal][14] };
        interp_table[3] = {parabola_table[table_number][i1_0][i2_1][i3_1][normal][1],parabola_table[table_number][i1_0][i2_1][i3_1][normal][2],parabola_table[table_number][i1_0][i2_1][i3_1][normal][3],parabola_table[table_number][i1_0][i2_1][i3_1][normal][14] };
        interp_table[4] = {parabola_table[table_number][i1_1][i2_0][i3_0][normal][1],parabola_table[table_number][i1_1][i2_0][i3_0][normal][2],parabola_table[table_number][i1_1][i2_0][i3_0][normal][3],parabola_table[table_number][i1_1][i2_0][i3_0][normal][14] };
        interp_table[5] = {parabola_table[table_number][i1_1][i2_0][i3_1][normal][1],parabola_table[table_number][i1_1][i2_0][i3_1][normal][2],parabola_table[table_number][i1_1][i2_0][i3_1][normal][3],parabola_table[table_number][i1_1][i2_0][i3_1][normal][14] };
        interp_table[6] = {parabola_table[table_number][i1_1][i2_1][i3_0][normal][1],parabola_table[table_number][i1_1][i2_1][i3_0][normal][2],parabola_table[table_number][i1_1][i2_1][i3_0][normal][3],parabola_table[table_number][i1_1][i2_1][i3_0][normal][14] };
        interp_table[7] = {parabola_table[table_number][i1_1][i2_1][i3_1][normal][1],parabola_table[table_number][i1_1][i2_1][i3_1][normal][2],parabola_table[table_number][i1_1][i2_1][i3_1][normal][3],parabola_table[table_number][i1_1][i2_1][i3_1][normal][14] };



       if(table_number == 5){ // checking normal

        int modified_normal;
        modified_normal = normal;

        cout << " modified_normal = "<< modified_normal << endl;
        if (normal == 1){k*=-1;b*=-1;d*=-1;c*=-1;}

        std::vector< Type > grad_given_par{2*k*intersection[4]+b , c}; // will it work like this?

        cout << " gradient of given parabola = " << grad_given_par[0] << " " << grad_given_par[1] << endl;

        std::vector< Type > grad_interp_par(2);

        grad_interp_par = {2*parabola_table[table_number][i1_0][i2_0][i3_0][normal][10]*  intersection[4]   +  parabola_table[table_number][i1_0][i2_0][i3_0][normal][11] , parabola_table[table_number][i1_0][i2_0][i3_0][normal][13]};
//         if(grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0) modified_normal = (normal+1)%2;
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;
        //TODO check what happens if it is 0.
        cout << " gradient of interp parabola = " << grad_interp_par[0] << " " << grad_interp_par[1] << ". dot product " << grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] << " modified normal = " << modified_normal << endl;

          interp_table[0] = {parabola_table[table_number][i1_0][i2_0][i3_0][modified_normal][1],parabola_table[table_number][i1_0][i2_0][i3_0][modified_normal][2],parabola_table[table_number][i1_0][i2_0][i3_0][modified_normal][3],parabola_table[table_number][i1_0][i2_0][i3_0][modified_normal][14] };

        grad_interp_par = {2*parabola_table[table_number][i1_0][i2_0][i3_1][normal][10]*  intersection[4]   +  parabola_table[table_number][i1_0][i2_0][i3_1][normal][11] , parabola_table[table_number][i1_0][i2_0][i3_1][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;
        cout << " gradient of interp parabola = " << grad_interp_par[0] << " " << grad_interp_par[1] << ". dot product " << grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1]<< " modified normal = " << modified_normal << endl;

          interp_table[1] = {parabola_table[table_number][i1_0][i2_0][i3_1][modified_normal][1],parabola_table[table_number][i1_0][i2_0][i3_1][modified_normal][2],parabola_table[table_number][i1_0][i2_0][i3_1][modified_normal][3],parabola_table[table_number][i1_0][i2_0][i3_1][modified_normal][14]};

        grad_interp_par = {2*parabola_table[table_number][i1_0][i2_1][i3_0][normal][10]*  intersection[4]   +  parabola_table[table_number][i1_0][i2_1][i3_0][normal][11] , parabola_table[table_number][i1_0][i2_1][i3_0][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;
        cout << " gradient of interp parabola = " << grad_interp_par[0] << " " << grad_interp_par[1] << ". dot product " << grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1]<< " modified normal = " << modified_normal << endl;

          interp_table[2] = {parabola_table[table_number][i1_0][i2_1][i3_0][modified_normal][1],parabola_table[table_number][i1_0][i2_1][i3_0][modified_normal][2],parabola_table[table_number][i1_0][i2_1][i3_0][modified_normal][3],parabola_table[table_number][i1_0][i2_1][i3_0][modified_normal][14] };

        grad_interp_par = {2*parabola_table[table_number][i1_0][i2_1][i3_1][normal][10]* intersection[4]   +  parabola_table[table_number][i1_0][i2_1][i3_1][normal][11] , parabola_table[table_number][i1_0][i2_1][i3_1][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;
        cout << " gradient of interp parabola = " << grad_interp_par[0] << " " << grad_interp_par[1] << ". dot product " << grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1]<< " modified normal = " << modified_normal << endl;

          interp_table[3] = {parabola_table[table_number][i1_0][i2_1][i3_1][modified_normal][1],parabola_table[table_number][i1_0][i2_1][i3_1][modified_normal][2],parabola_table[table_number][i1_0][i2_1][i3_1][modified_normal][3],parabola_table[table_number][i1_0][i2_1][i3_1][modified_normal][14] };

        grad_interp_par = {2*parabola_table[table_number][i1_1][i2_0][i3_0][normal][10]* intersection[4]   +  parabola_table[table_number][i1_1][i2_0][i3_0][normal][11] , parabola_table[table_number][i1_1][i2_0][i3_0][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;
        cout << " gradient of interp parabola = " << grad_interp_par[0] << " " << grad_interp_par[1] << ". dot product " << grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1]<< " modified normal = " << modified_normal << endl;

          interp_table[4] = {parabola_table[table_number][i1_1][i2_0][i3_0][modified_normal][1],parabola_table[table_number][i1_1][i2_0][i3_0][modified_normal][2],parabola_table[table_number][i1_1][i2_0][i3_0][modified_normal][3],parabola_table[table_number][i1_1][i2_0][i3_0][modified_normal][14] };

        grad_interp_par = {2*parabola_table[table_number][i1_1][i2_0][i3_1][normal][10]* intersection[4]  +  parabola_table[table_number][i1_1][i2_0][i3_1][normal][11] , parabola_table[table_number][i1_1][i2_0][i3_1][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;
        cout << " gradient of interp parabola = " << grad_interp_par[0] << " " << grad_interp_par[1] << ". dot product " << grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1]<< " modified normal = " << modified_normal << endl;

          interp_table[5] = {parabola_table[table_number][i1_1][i2_0][i3_1][modified_normal][1],parabola_table[table_number][i1_1][i2_0][i3_1][modified_normal][2],parabola_table[table_number][i1_1][i2_0][i3_1][modified_normal][3],parabola_table[table_number][i1_1][i2_0][i3_1][modified_normal][14] };

        grad_interp_par = {2*parabola_table[table_number][i1_1][i2_1][i3_0][normal][10]*  intersection[4]   +  parabola_table[table_number][i1_1][i2_1][i3_0][normal][11] , parabola_table[table_number][i1_1][i2_1][i3_0][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;
        cout << " gradient of interp parabola = " << grad_interp_par[0] << " " << grad_interp_par[1] << ". dot product " << grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1]<< " modified normal = " << modified_normal << endl;

          interp_table[6] = {parabola_table[table_number][i1_1][i2_1][i3_0][modified_normal][1],parabola_table[table_number][i1_1][i2_1][i3_0][modified_normal][2],parabola_table[table_number][i1_1][i2_1][i3_0][modified_normal][3],parabola_table[table_number][i1_1][i2_1][i3_0][modified_normal][14] };

        grad_interp_par = {2*parabola_table[table_number][i1_1][i2_1][i3_1][normal][10]* intersection[4]   +  parabola_table[table_number][i1_1][i2_1][i3_1][normal][11] , parabola_table[table_number][i1_1][i2_1][i3_1][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;
        cout << " gradient of interp parabola = " << grad_interp_par[0] << " " << grad_interp_par[1] << ". dot product " << grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1]<< " modified normal = " << modified_normal << endl;


          interp_table[7] = {parabola_table[table_number][i1_1][i2_1][i3_1][modified_normal][1],parabola_table[table_number][i1_1][i2_1][i3_1][modified_normal][2],parabola_table[table_number][i1_1][i2_1][i3_1][modified_normal][3],parabola_table[table_number][i1_1][i2_1][i3_1][modified_normal][14] };

       }

        for(unsigned int i = 0; i <=7 ; i++){
          cout <<" direct table  : " <<  interp_table[i][0] << " " <<interp_table[i][1] << " " <<interp_table[i][2] << " " <<interp_table[i][3]<< endl;
        }

          {  // all the parabola equation of the inpolation points
            cout<< "\n "<< parabola_table[table_number][i1_0][i2_0][i3_0][normal][10]<<"x^2 + "<<parabola_table[table_number][i1_0][i2_0][i3_0][normal][11]<<"x + "<<parabola_table[table_number][i1_0][i2_0][i3_0][normal][12]<< " + " << parabola_table[table_number][i1_0][i2_0][i3_0][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table[table_number][i1_0][i2_0][i3_1][normal][10]<<"x^2 + "<<parabola_table[table_number][i1_0][i2_0][i3_1][normal][11]<<"x + "<<parabola_table[table_number][i1_0][i2_0][i3_1][normal][12]<< " + " << parabola_table[table_number][i1_0][i2_0][i3_1][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table[table_number][i1_0][i2_1][i3_0][normal][10]<<"x^2 + "<<parabola_table[table_number][i1_0][i2_1][i3_0][normal][11]<<"x + "<<parabola_table[table_number][i1_0][i2_1][i3_0][normal][12]<< " + " << parabola_table[table_number][i1_0][i2_1][i3_0][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table[table_number][i1_0][i2_1][i3_1][normal][10]<<"x^2 + "<<parabola_table[table_number][i1_0][i2_1][i3_1][normal][11]<<"x + "<<parabola_table[table_number][i1_0][i2_1][i3_1][normal][12]<< " + " << parabola_table[table_number][i1_0][i2_1][i3_1][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table[table_number][i1_1][i2_0][i3_0][normal][10]<<"x^2 + "<<parabola_table[table_number][i1_1][i2_0][i3_0][normal][11]<<"x + "<<parabola_table[table_number][i1_1][i2_0][i3_0][normal][12]<< " + " << parabola_table[table_number][i1_1][i2_0][i3_0][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table[table_number][i1_1][i2_0][i3_1][normal][10]<<"x^2 + "<<parabola_table[table_number][i1_1][i2_0][i3_1][normal][11]<<"x + "<<parabola_table[table_number][i1_1][i2_0][i3_1][normal][12]<< " + " << parabola_table[table_number][i1_1][i2_0][i3_1][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table[table_number][i1_1][i2_1][i3_0][normal][10]<<"x^2 + "<<parabola_table[table_number][i1_1][i2_1][i3_0][normal][11]<<"x + "<<parabola_table[table_number][i1_1][i2_1][i3_0][normal][12]<< " + " << parabola_table[table_number][i1_1][i2_1][i3_0][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table[table_number][i1_1][i2_1][i3_1][normal][10]<<"x^2 + "<<parabola_table[table_number][i1_1][i2_1][i3_1][normal][11]<<"x + "<<parabola_table[table_number][i1_1][i2_1][i3_1][normal][12]<< " + " << parabola_table[table_number][i1_1][i2_1][i3_1][normal][13] <<"y = 0"<<endl;

          }





        }

//         for(unsigned int i = 0; i <=7 ; i++){
//           for (unsigned int count_x = 0; count_x < parabola_table[table_number].size(); count_x++ ){
//             if(fabs(parabola_table[table_number][count_x][normal][1] - interp_table[i][0]) < 0.00000000001){
// //                   cout<< "count x =" << count_x << " " << parabola_table[table_number][count_x][normal][1] << " " << interp_table[i][0]  <<endl;
//               for (unsigned int count_y = count_x; count_y < parabola_table[table_number].size(); count_y++ ){
//                 if(fabs(parabola_table[table_number][count_y][normal][2] - interp_table[i][1]) < 0.000000000001){
// //                   cout<< "count y =" << count_y <<endl;
//                   for (unsigned int count_z = count_y; count_z < parabola_table[table_number].size(); count_z++ ){
//                     if(fabs(parabola_table[table_number][count_z][normal][3] - interp_table[i][2]) < 0.000000000001 ){
// //                       cout<< "count z =" << count_z <<endl;
//                       interp_table[i][3] = parabola_table[table_number][count_z][normal][14] ;
// //                       cout << " area = "<< parabola_table[table_number][count_z][normal][14];
//                       count_y = parabola_table[table_number].size()+1;
//                       count_x = parabola_table[table_number].size()+1;
//                       break;
//                     }
//                   }
//                 }
//               }
//             }
//           }
//         }


      else if (intersect_number>2){
        cout << " Table used - 4 intersection " << endl;
        {
//           for(unsigned int i = 0; i <=7 ; i++){
//             for (unsigned int count_x = 0; count_x < parabola_table_4intersection[table_number].size(); count_x++ ){
//               if(fabs(parabola_table_4intersection[table_number][count_x][normal][1] - interp_table[i][0]) < 0.000001){
//   //                   cout<< "count x =" << count_x << " " << parabola_table_4intersection[table_number][count_x][normal][1] << " " << interp_table[i][0]  <<endl;
//                 for (unsigned int count_y = count_x; count_y < parabola_table_4intersection[table_number].size(); count_y++ ){
//                   if(fabs(parabola_table_4intersection[table_number][count_y][normal][2] - interp_table[i][1]) < 0.000001){
//   //                   cout<< "count y =" << count_y <<endl;
//                     for (unsigned int count_z = count_y; count_z < parabola_table_4intersection[table_number].size(); count_z++ ){
//                       if(fabs(parabola_table_4intersection[table_number][count_z][normal][3] - interp_table[i][2]) < 0.000001 ){
//   //                       cout<< "count z =" << count_z <<endl;
//                         interp_table[i][3] = parabola_table_4intersection[table_number][count_z][normal][14] ;
//   //                       cout << " area = "<< parabola_table_4intersection[table_number][count_z][normal][14];
//                         count_y = parabola_table_4intersection[table_number].size()+1;
//                         count_x = parabola_table_4intersection[table_number].size()+1;
//                         break;
//                       }
//                     }
//                   }
//                 }
//               }
//             }
//           }
        }

        i1_0 = floor(static_cast<double>(interp_point[0] * partition)) ;
        i1_1 = ceil( static_cast<double>(interp_point[0] * partition)) ;
        i2_0 = floor(static_cast<double>(interp_point[1] * partition)) ;
        i2_1 = ceil( static_cast<double>(interp_point[1] * partition)) ;
        i3_0 = floor(static_cast<double>(interp_point[2] * partition)) ;
        i3_1 = ceil( static_cast<double>(interp_point[2] * partition)) ;



        interp_table[0][3] = parabola_table_4intersection[table_number][i1_0][i2_0][i3_0][normal][14] ;
        interp_table[1][3] = parabola_table_4intersection[table_number][i1_0][i2_0][i3_1][normal][14] ;
        interp_table[2][3] = parabola_table_4intersection[table_number][i1_0][i2_1][i3_0][normal][14] ;
        interp_table[3][3] = parabola_table_4intersection[table_number][i1_0][i2_1][i3_1][normal][14] ;
        interp_table[4][3] = parabola_table_4intersection[table_number][i1_1][i2_0][i3_0][normal][14] ;
        interp_table[5][3] = parabola_table_4intersection[table_number][i1_1][i2_0][i3_1][normal][14] ;
        interp_table[6][3] = parabola_table_4intersection[table_number][i1_1][i2_1][i3_0][normal][14] ;
        interp_table[7][3] = parabola_table_4intersection[table_number][i1_1][i2_1][i3_1][normal][14] ;

        cout<< " given parabola = " << k << " " << b << " "<< d << " "<< c <<endl;
//         cout<< " intersection points = " << interp_point[0] << " " << interp_point[1] << " "<< interp_point[2] <<endl;
        cout << " interpolation table = " << endl;
        for(unsigned int i = 0; i <=7 ; i++){
          cout << interp_table[i][0] << " " <<interp_table[i][1] << " " <<interp_table[i][2] << " " <<interp_table[i][3]<< endl;
        }

        interp_table[0] = {parabola_table_4intersection[table_number][i1_0][i2_0][i3_0][normal][1],parabola_table_4intersection[table_number][i1_0][i2_0][i3_0][normal][2],parabola_table_4intersection[table_number][i1_0][i2_0][i3_0][normal][3],parabola_table_4intersection[table_number][i1_0][i2_0][i3_0][normal][14] };
        interp_table[1] = {parabola_table_4intersection[table_number][i1_0][i2_0][i3_1][normal][1],parabola_table_4intersection[table_number][i1_0][i2_0][i3_1][normal][2],parabola_table_4intersection[table_number][i1_0][i2_0][i3_1][normal][3],parabola_table_4intersection[table_number][i1_0][i2_0][i3_1][normal][14]};
        interp_table[2] = {parabola_table_4intersection[table_number][i1_0][i2_1][i3_0][normal][1],parabola_table_4intersection[table_number][i1_0][i2_1][i3_0][normal][2],parabola_table_4intersection[table_number][i1_0][i2_1][i3_0][normal][3],parabola_table_4intersection[table_number][i1_0][i2_1][i3_0][normal][14] };
        interp_table[3] = {parabola_table_4intersection[table_number][i1_0][i2_1][i3_1][normal][1],parabola_table_4intersection[table_number][i1_0][i2_1][i3_1][normal][2],parabola_table_4intersection[table_number][i1_0][i2_1][i3_1][normal][3],parabola_table_4intersection[table_number][i1_0][i2_1][i3_1][normal][14] };
        interp_table[4] = {parabola_table_4intersection[table_number][i1_1][i2_0][i3_0][normal][1],parabola_table_4intersection[table_number][i1_1][i2_0][i3_0][normal][2],parabola_table_4intersection[table_number][i1_1][i2_0][i3_0][normal][3],parabola_table_4intersection[table_number][i1_1][i2_0][i3_0][normal][14] };
        interp_table[5] = {parabola_table_4intersection[table_number][i1_1][i2_0][i3_1][normal][1],parabola_table_4intersection[table_number][i1_1][i2_0][i3_1][normal][2],parabola_table_4intersection[table_number][i1_1][i2_0][i3_1][normal][3],parabola_table_4intersection[table_number][i1_1][i2_0][i3_1][normal][14] };
        interp_table[6] = {parabola_table_4intersection[table_number][i1_1][i2_1][i3_0][normal][1],parabola_table_4intersection[table_number][i1_1][i2_1][i3_0][normal][2],parabola_table_4intersection[table_number][i1_1][i2_1][i3_0][normal][3],parabola_table_4intersection[table_number][i1_1][i2_1][i3_0][normal][14] };
        interp_table[7] = {parabola_table_4intersection[table_number][i1_1][i2_1][i3_1][normal][1],parabola_table_4intersection[table_number][i1_1][i2_1][i3_1][normal][2],parabola_table_4intersection[table_number][i1_1][i2_1][i3_1][normal][3],parabola_table_4intersection[table_number][i1_1][i2_1][i3_1][normal][14] };
        for(unsigned int i = 0; i <=7 ; i++){
          cout <<" direct table  : " <<  interp_table[i][0] << " " <<interp_table[i][1] << " " <<interp_table[i][2] << " " <<interp_table[i][3]<< endl;
        }

                  {  // all the parabola equation of the inpolation points
            cout<< "\n "<< parabola_table_4intersection[table_number][i1_0][i2_0][i3_0][normal][10]<<"x^2 + "<<parabola_table_4intersection[table_number][i1_0][i2_0][i3_0][normal][11]<<"x + "<<parabola_table_4intersection[table_number][i1_0][i2_0][i3_0][normal][12]<< " + " << parabola_table_4intersection[table_number][i1_0][i2_0][i3_0][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table_4intersection[table_number][i1_0][i2_0][i3_1][normal][10]<<"x^2 + "<<parabola_table_4intersection[table_number][i1_0][i2_0][i3_1][normal][11]<<"x + "<<parabola_table_4intersection[table_number][i1_0][i2_0][i3_1][normal][12]<< " + " << parabola_table_4intersection[table_number][i1_0][i2_0][i3_1][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table_4intersection[table_number][i1_0][i2_1][i3_0][normal][10]<<"x^2 + "<<parabola_table_4intersection[table_number][i1_0][i2_1][i3_0][normal][11]<<"x + "<<parabola_table_4intersection[table_number][i1_0][i2_1][i3_0][normal][12]<< " + " << parabola_table_4intersection[table_number][i1_0][i2_1][i3_0][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table_4intersection[table_number][i1_0][i2_1][i3_1][normal][10]<<"x^2 + "<<parabola_table_4intersection[table_number][i1_0][i2_1][i3_1][normal][11]<<"x + "<<parabola_table_4intersection[table_number][i1_0][i2_1][i3_1][normal][12]<< " + " << parabola_table_4intersection[table_number][i1_0][i2_1][i3_1][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table_4intersection[table_number][i1_1][i2_0][i3_0][normal][10]<<"x^2 + "<<parabola_table_4intersection[table_number][i1_1][i2_0][i3_0][normal][11]<<"x + "<<parabola_table_4intersection[table_number][i1_1][i2_0][i3_0][normal][12]<< " + " << parabola_table_4intersection[table_number][i1_1][i2_0][i3_0][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table_4intersection[table_number][i1_1][i2_0][i3_1][normal][10]<<"x^2 + "<<parabola_table_4intersection[table_number][i1_1][i2_0][i3_1][normal][11]<<"x + "<<parabola_table_4intersection[table_number][i1_1][i2_0][i3_1][normal][12]<< " + " << parabola_table_4intersection[table_number][i1_1][i2_0][i3_1][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table_4intersection[table_number][i1_1][i2_1][i3_0][normal][10]<<"x^2 + "<<parabola_table_4intersection[table_number][i1_1][i2_1][i3_0][normal][11]<<"x + "<<parabola_table_4intersection[table_number][i1_1][i2_1][i3_0][normal][12]<< " + " << parabola_table_4intersection[table_number][i1_1][i2_1][i3_0][normal][13] <<"y = 0"<<endl;

            cout<< parabola_table_4intersection[table_number][i1_1][i2_1][i3_1][normal][10]<<"x^2 + "<<parabola_table_4intersection[table_number][i1_1][i2_1][i3_1][normal][11]<<"x + "<<parabola_table_4intersection[table_number][i1_1][i2_1][i3_1][normal][12]<< " + " << parabola_table_4intersection[table_number][i1_1][i2_1][i3_1][normal][13] <<"y = 0"<<endl;

          }


      }

  //                 parabola_table[table][count][normal][0]





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

  Type cc = c_0 * (1-z_d) + c_1 * z_d ;

  return cc;
}



template <class Type>
void justIntersectionNumber(int &intersect_number, Parabola <Type> &parabola){

  Type k = parabola.k;
  Type b = parabola.b;
  Type d = parabola.d;
  Type c = 1;
//   cout<< " parabola I get from solving system of linear equation :  " << parabola.k <<"x^2 + "<< parabola.b <<"x + "<< parabola.d << "+" <<  c << " y=0"  <<endl;

      if (-d>=0 && -d<=1){ //LEFT
        intersect_number += 1;
      }
      // LEFT-TOP solve kx^2+bx+d-1 =0  ; Table 0
      if (k == 0){
      Type x =  (-1-d)/b ;
        if(x <= 1 && x>= 0) {
          intersect_number += 1;
        }
      }
      else {
        Type delta = b*b - 4*k*(d+1);
        cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
        if (delta >= 0){
              Type sqrtdelta = sqrt(delta);
              int sign = (k > 0) ? 1 : -1;

              for(unsigned i = 0; i < 2; i++) {
                Type x = (- b - sign * sqrtdelta) / (2 * k);
    //             cout<< "Top x = "<< x<< endl;
                if(x > 1) break;
                else if(x >= 0) {
                  intersect_number += 1;
                }
                sign *= -1;
              }
            }
      }
      Type y_1=-(k+b+d); //LEFT-RIGHT x=1 ; Table 1
      if (y_1 >= 0 && y_1 <= 1){
          intersect_number += 1;
      }

        //LEFT-BOTTOM  solve kx^2+bx+d =0 ; Table 2
//                 cout << " bottom = " << bottom ;
      if (k == 0){
          Type x =  -d/b ;
          if(x <= 1 && x>= 0) {
            intersect_number += 1;
          }
      }

      else {
          Type delta = b*b - 4*k*d;
          if(delta >=0) {
            Type sqrtdelta = sqrt(delta);
            int sign = (k > 0) ? 1 : -1;
            for(unsigned i = 0; i < 2; i++) {
              Type x = (- b - sign * sqrtdelta) / (2 * k);
//               cout << " bottom root = " << x ;
              if(x > 1) break;
              else if(x >= 0) {
                intersect_number += 1;
              }
              sign *= -1;
            }
          }
      }
}




int main() {

  typedef cpp_bin_float_oct Type;      //     typedef double Type;
    std::cout.precision(16);
  unsigned int m = 0;
  unsigned int n = 0;
  int s = 0;
  unsigned int partition = 10;
  Parabola <Type> parabola;


    // Specify the path to your CSV file
    std::string csvFilePath = "data.csv";

    // Open the CSV file
    std::ifstream file(csvFilePath);

    // Check if the file is open
    if (!file.is_open()) {
        // If the file cannot be opened, print an error message and exit
        std::cerr << "Error opening file: " << csvFilePath << std::endl;
        return 1;
    }
    cout<< "file read successfully"<< endl;
    // Define vectors to store data from CSV
    std::vector<std::vector<std::string>> data;

    // Read the header line separately
    std::string header;
    std::getline(file, header);

    // Add the new column name to the header
    header += ",exact value,difference";

    // Store the header in the data vector
    std::vector<std::string> headerRow;
    std::istringstream headerIss(header);
    std::string headerToken;
    while (std::getline(headerIss, headerToken, ',')) {
        // Tokenize the header by comma and store values in the headerRow vector
        headerRow.push_back(headerToken);
    }
    // Add the headerRow vector to the data vector
    data.push_back(headerRow);

    // Assuming header is a vector containing column names

//     for (size_t i = 0; i < header.size(); ++i) {
// //     std::cout << "Header[" << i << "] = " << header[i] << std::endl;
//     }


    // Read the rest of the file
    std::string line;
    while(std::getline(file, line)){
    //for (int lineNumber = 0; std::getline(file, line); ++lineNumber) {
        std::istringstream iss(line);
        std::vector<std::string> row;
        int intersect_number(0) ;

        // Tokenize the line by comma and store values in the row vector
        std::string token;
        while (std::getline(iss, token, ',')) {
            row.push_back(token);
        }

//             for (size_t i = 0; i < row.size(); ++i) {
//               std::cout << "row [" << i << "] = " << row[i] << std::endl;
//             }

        // Extract numbers from the row and add them
        PointT<Type> p1,p2,p3 ;
        double F_values;

//         cout<< "inside line "<< " row size = "<< row.size() << " row = " << row[0] << " " << row[1] << " " << row[2]<<" " <<row[3] << " " <<row[4] << endl;


        for (int i = 1; i < row.size(); ++i) {
            // Convert each element of the row to a double and add to the sum
            double value = std::stod(row[i]);
//             std::cout << "Debug: i = " << i << ", row.size() = " << row.size() << "value " << row[i] << std::endl;
            if      (i==1) p1.x = static_cast<Type>(value);
            else if (i==2) p1.y = static_cast<Type>(value);
            else if (i==3) p2.x = static_cast<Type>(value);
            else if (i==4) p2.y = static_cast<Type>(value);
            else if (i==5) p3.x = static_cast<Type>(value);
            else if (i==6) p3.y = static_cast<Type>(value);
            else if (i==7) F_values = value ;


//             cout <<" " << p1.x << " "<<p1.y << " "<<p2.x << " "<<p2.y << " "<<p3.x << " "<<p3.y << " "<<F_values <<endl;
        }
        //now that extracted these points we have to generate the parabola equations.
        parabola = get_parabola_equation(p1, p2, p3);

        justIntersectionNumber(intersect_number, parabola);

        Type formula_area ;
        double area;


          if (intersect_number == 2){
                std::vector< std::pair <Type, Type> > I1, I2, I3 ;
                std::vector <Type> pol1(3, 0);
                std::vector <Type> pol2(3, 0);
                Type A1,A2,A3, a=0 , c= 1;
                pol1[0] = parabola.k;
                pol1[1] = parabola.b;
                pol1[2] = parabola.d + c;
                pol2[0] = parabola.k;
                pol2[1] = parabola.b;
                pol2[2] = parabola.d;
//                 cout<< " given parabola : =" << pol2[0] << "x^2 + " << pol2[1] << "x +"<< pol2[2] << " + "<< c <<"y=0"<< endl;
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
                formula_area = A1 + A2 + A3;
                area = static_cast<double>(formula_area);
          }

          else{
              area = F_values ;
          }

        double difference = area - F_values;

        // Convert the sum to a string and add it to the row vector as the "exact value"
        row.push_back(std::to_string(area));
        row.push_back(std::to_string(difference));

        // Add the row vector to the data vector
        data.push_back(row);
            cout<<" calculation on going, area = " << area << " difference = "<< difference  <<endl;
    }


    cout<<" file reading done " <<endl;
    // Close the file
    file.close();

    // Specify the path for the new CSV file
    std::string newCsvFilePath = "new_data.csv";


    // Open the new file for writing
    std::ofstream outFile(newCsvFilePath);

    // Check if the new file is open
    if (!outFile.is_open()) {
        // If the new file cannot be opened for writing, print an error message and exit
        std::cerr << "Error opening new file for writing: " << newCsvFilePath << std::endl;
        return 1;
    }

    outFile << std::fixed << std::setprecision(15);

    // Write the modified data to the new CSV file
    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            // Write each element of the row to the new file
            outFile << row[i];

            // Add a comma if it's not the last element in the row
            if (i < row.size() - 1) {
                outFile << ",";
            }
        }
        // Add a newline character after each row
        outFile << "\n";
    }

    // Close the new output file
    outFile.close();

    // Print a success message
    std::cout << "New file successfully created and saved." << std::endl;

    return 0;
}




// int main() {
//
//   typedef cpp_bin_float_oct Type;      //     typedef double Type;
//   std::cout.precision(16);
//   unsigned int m = 0;
//   unsigned int n = 0;
//   int s = 0;
//   unsigned int partition = 10;
//   Parabola <Type> parabola;
//
//
//     // Specify the path to your CSV file
//     std::string csvFilePath = "data.csv";
//
//     // Open the CSV file
//     std::ifstream file(csvFilePath);
//
//     // Check if the file is open
//     if (!file.is_open()) {
//         // If the file cannot be opened, print an error message and exit
//         std::cerr << "Error opening file: " << csvFilePath << std::endl;
//         return 1;
//     }
//     cout<< "file read successfully"<< endl;
//     // Define vectors to store data from CSV
//     std::vector<std::vector<std::string>> data;
//
//     // Read the header line separately
//     std::string header;
//     std::getline(file, header);
//
//     // Add the new column name to the header
//     header += ",exact value,difference";
//
//     // Store the header in the data vector
//     std::vector<std::string> headerRow;
//     std::istringstream headerIss(header);
//     std::string headerToken;
//     while (std::getline(headerIss, headerToken, ',')) {
//         // Tokenize the header by comma and store values in the headerRow vector
//         headerRow.push_back(headerToken);
//     }
//     // Add the headerRow vector to the data vector
//     data.push_back(headerRow);
//
//     // Assuming header is a vector containing column names
//
// //     for (size_t i = 0; i < header.size(); ++i) {
// // //     std::cout << "Header[" << i << "] = " << header[i] << std::endl;
// //     }
//
//
//     // Read the rest of the file
//     std::string line;
//     while(std::getline(file, line)){
//     //for (int lineNumber = 0; std::getline(file, line); ++lineNumber) {
//         std::istringstream iss(line);
//         std::vector<std::string> row;
//
//         // Tokenize the line by comma and store values in the row vector
//         std::string token;
//         while (std::getline(iss, token, ',')) {
//             row.push_back(token);
//         }
//
// //             for (size_t i = 0; i < row.size(); ++i) {
// //               std::cout << "row [" << i << "] = " << row[i] << std::endl;
// //             }
//
//         // Extract numbers from the row and add them
//         PointT<Type> p1,p2,p3 ;
//         double F_values;
//
// //         cout<< "inside line "<< " row size = "<< row.size() << " row = " << row[0] << " " << row[1] << " " << row[2]<<" " <<row[3] << " " <<row[4] << endl;
//
//
//             Parabola <Type> parabola;
//             int intersect_number;
//             std::vector <Type> intersection;
//             std::vector <Type> interp_point;    //never used in this function. it was used in interpolation;
//             unsigned int table;
//             std::vector< std::pair <Type, Type> > I1, I2, I3 ;
//             std::vector <Type> pol1(3, 0);
//             std::vector <Type> pol2(3, 0);
//             Type A1,A2,A3, a=0 , c= 1;
//             Type formula_area ;
//             double area;
//             for (int i = 1; i < row.size(); ++i) {
//                 // Convert each element of the row to a double and add to the sum
//                 double value = std::stod(row[i]);
//     //             std::cout << "Debug: i = " << i << ", row.size() = " << row.size() << "value " << row[i] << std::endl;
//                 if      (i==1) p1.x = static_cast<Type>(value);
//                 else if (i==2) p1.y = static_cast<Type>(value);
//                 else if (i==3) p2.x = static_cast<Type>(value);
//                 else if (i==4) p2.y = static_cast<Type>(value);
//                 else if (i==5) p3.x = static_cast<Type>(value);
//                 else if (i==6) p3.y = static_cast<Type>(value);
//                 else if (i==7) {
//                   cout << " it is 7 " << endl;
//                   F_values = value ;
//                 }
//
//     //             cout <<" " << p1.x << " "<<p1.y << " "<<p2.x << " "<<p2.y << " "<<p3.x << " "<<p3.y << " "<<F_values <<endl;
//             }
//             Type det = p1.x * p1.x * (p2.x - p3.x) -p1.x* (p2.x*p2.x - p3.x*p3.x)+ p2.x*p3.x*(p2.x - p3.x) ;// only sort out the points parallel to y axis
//             parabola = get_parabola_equation(p1, p2, p3);
//             CheckIntersection <Type> (intersect_number, table, intersection, interp_point, parabola);
//
//             bool concaveUp = 0;
//             bool do_line = 0;
//             if(parabola.k < 0) concaveUp = 1; //kx^2+bx+d+cy = 0 ; c is initialy 1 => y = -kx^2 - bx - d .
//             if (intersect_number > 2){
//                 if(!concaveUp){
//                     if (table == 0 || table == 1 || table == 4 || table == 5){
//                         do_line = 1 ;
//                     }
//                 }
//                 else {
//                     if (table == 5 || table == 6 ){
//                         do_line = 1 ;
//                     }
//
//                 }
//             }
//
//             if (intersect_number == 2 || do_line == 0){ // we are not changing 4 intersection when do table-1 left-right
//                 if (det !=0){
//                     pol1[0] = parabola.k;
//                     pol1[1] = parabola.b;
//                     pol1[2] = parabola.d + c;
//                     pol2[0] = parabola.k;
//                     pol2[1] = parabola.b;
//                     pol2[2] = parabola.d;
//                 }
//                 else {  //TODO decode what I did here. Couldn't figure it out what I did here. We are not using this probably after introducing epsilon
//
//                     pol1[0] = static_cast<Type>(0);
//                     pol1[1] = static_cast<Type>(-1);
//                     pol1[2] = p1.x;
//                     pol2[0] = pol1[0];
//                     pol2[1] = pol1[1];
//                     pol2[2] = p1.x;
//                     c=static_cast<Type>(0);
//
//                 }
//             }
//             else if (intersect_number > 2){
//                 Type slope = (p3.y-p1.y)/(p3.x-p1.x);
//                 c=static_cast<Type>(1);
//                 pol2[0] = static_cast<Type>(0);    //k=0
//                 pol2[1] = -slope;
//                 pol2[2] = slope*p1.x - p1.y ;
//                 pol1[0] = pol2[0];    //k=0
//                 pol1[1] = pol2[1];
//                 pol1[2] = pol2[2] + c ;
//             }
//
// //                 cout<< " given parabola : =" << pol2[0] << "x^2 + " << pol2[1] << "x +"<< pol2[2] << " + "<< c <<"y=0"<< endl;
//             GetIntervalall<Type, double>(pol1, pol2, I1, I2, I3);
//             if(I1.size() > 0) {
//                A1 = easy_integral_A3(m, n, s, a, c, pol2, I1) -  easy_integral_A2(m, n, s, a, c, pol2, I1);
//             }
//             if(I2.size() > 0) {
//                A2 = easy_integral_A2(m, n, s, a, c, pol2, I2);
//             }
//             if(I3.size() > 0) {
//                A3 = easy_integral_A3(m, n, s, a, c, pol2, I3);
//             }
//             formula_area = A1 + A2 + A3;
//             area = static_cast<double>(formula_area);
//
//
//         double difference = area - F_values;
//
//         // Convert the sum to a string and add it to the row vector as the "exact value"
//         row.push_back(std::to_string(area));
//         row.push_back(std::to_string(difference));
//
//         // Add the row vector to the data vector
//         data.push_back(row);
//             cout<<" calculation on going, area = " << area << " difference = "<< difference  <<endl;
//     }
//
//
//     cout<<" file reading done " <<endl;
//     // Close the file
//     file.close();
//
//     // Specify the path for the new CSV file
//     std::string newCsvFilePath = "new_data.csv";
//
//
//     // Open the new file for writing
//     std::ofstream outFile(newCsvFilePath);
//
//     // Check if the new file is open
//     if (!outFile.is_open()) {
//         // If the new file cannot be opened for writing, print an error message and exit
//         std::cerr << "Error opening new file for writing: " << newCsvFilePath << std::endl;
//         return 1;
//     }
//
//     outFile << std::fixed << std::setprecision(15);
//
//     // Write the modified data to the new CSV file
//     for (const auto& row : data) {
//         for (size_t i = 0; i < row.size(); ++i) {
//             // Write each element of the row to the new file
//             outFile << row[i];
//
//             // Add a comma if it's not the last element in the row
//             if (i < row.size() - 1) {
//                 outFile << ",";
//             }
//         }
//         // Add a newline character after each row
//         outFile << "\n";
//     }
//
//     // Close the new output file
//     outFile.close();
//
//     // Print a success message
//     std::cout << "New file successfully created and saved." << std::endl;
//
//     return 0;
// }
