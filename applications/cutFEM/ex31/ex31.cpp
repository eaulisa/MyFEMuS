
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
void random_polynomial(std::vector <Type> &a1, std::vector <Type> &a2) {
  a1[0] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
  a1[1] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
  a1[2] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
  a2[0] = a1[0] ;                                               // this gives us a = 0 ;
//   a2[1] = a1[1];
  a2[1] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
  a2[2] = a1[2];
//             std::cout <<"\n ** k = "<<a2[0] << "; b = " << a2[1] << "; d = " << a2[2] << "; a = " << a1[1] - a2[1]<< "; c = " << a1[2] - a2[2] << ";" << std::endl;
}

template <class Type>
Parabola <Type>  get_parabola_equation( const PointT <Type> p1,const PointT <Type> p2,const PointT <Type> p3) {
    Type x1 = p1.x, x2 = p2.x, x3 = p3.x;
    Type y1 = p1.y, y2 = p2.y, y3 = p3.y;
    Type k,b,d;
    Type det = p1.x * p1.x * (p2.x - p3.x) -p1.x* (p2.x*p2.x - p3.x*p3.x)+ p2.x*p3.x*(p2.x - p3.x);

//     Type det = x1 * x1 * (x2 - x3) -x1* (x2*x2 - x3*x3)+ x2*x3*(x2 - x3) ;
//     Type denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
    if(fabs(det) >0.0000000000000001){
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

    if (fabs(k) < 1.e-16) k = 0 ;
    if (fabs(b) < 1.e-16) b = 0 ;
    if (fabs(d) < 1.e-16) d = 0 ;

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
          Type sum1 = pow(pol1[2], pMax - r)/factorial<Type>(pMax - r) ;
          Type sum = sum1 * (pow(I3[w].second,2 * r + m + 1+n+1-i) - pow(I3[w].first, 2 * r + m + 1+n+1-i))/ (2 * r + m + 1+n+1-i);
          for(int p = 1; p <= pMax - r; p++) {
            sum1 *= pol1[1] * (pMax - r - p + 1) / (pol1[2] * p)   ;
            sum += sum1 * (pow(I3[w].second,2 * r + m+p + 1) - pow(I3[w].first, 2 * r +p+ m + 1))/ (2 * r + m+p + 1) ;
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
void CheckIntersection(int &intersect_number, unsigned int &table_number , std::vector <Type> &intersection, std::vector <Type> &interp_point, const Parabola <Type> &parabola){

  table_number = 9999;
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
//         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
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

/*
cout<< " " << " left " << left << " top "<< top << " right "<< right << " bottom " << bottom  << " table number :"<< table_number << " number of intersection " << intersect_number <<endl;*/

}

template <class Type>
Type find_trig_area_2intersection_formula(const unsigned &m, const unsigned &n, const int &s, const Type &a, Type c,const int &table,  PointT <Type> &p1,  PointT <Type> &p2, const PointT <Type> &p3){   //TODO we have 5 tables. each of them works differently. Work the math first.
    Type area(0);
    Type A1 (0), A2 (0), A3 (0);
    std::vector< std::pair <Type, Type> > I1, I2, I3 ;
    int intersect_number;
    std::vector <Type> pol1(3, 0);
    std::vector <Type> pol2(3, 0);
    std::vector <Type> intersection;
    std::vector <Type> interp_point;    //never used in this function. it was used in interpolation;
    unsigned int table_number = table;
    Parabola <Type> parabola;
    Type det = p1.x * p1.x * (p2.x - p3.x) -p1.x* (p2.x*p2.x - p3.x*p3.x)+ p2.x*p3.x*(p2.x - p3.x) ;// only sort out the points parallel to y axis
    parabola = get_parabola_equation(p1, p2, p3);
    CheckIntersection <Type> (intersect_number, table_number, intersection, interp_point, parabola);

    Type k = parabola.k;
    Type b = parabola.b;
    Type d = parabola.d;
    Type singleintersection;

    bool do_line = 0;

    if (table == 1){ //we only use modified integrals if it is concave down
      if (k>0) {
        do_line = 1;
        Type delta = (b+1.)*(b+1.) - 4*k*d;
//         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
        if (delta >= 0){
              Type sqrtdelta = sqrt(delta);
              int sign = (k > 0) ? 1 : -1;
              for(unsigned i = 0; i < 2; i++) {
                Type x = (- (b+1) - sign * sqrtdelta) / (2 * k);
    //             cout<< "Top x = "<< x<< endl;
                 if(x > 0 && x > p1.x) {
                    p1.x = x;
                }
                sign *= -1;
              }
            }
      }
    }

    else if(table ==2){  // There are six possible cases we have to use modified integrals
     do_line = 1;
     if(k>=0){//concave down (2 possible scenerio)
      if(p1.x < p2.x){ //case (a) take highest p1.x
        Type delta = (b+1.)*(b+1.) - 4*k*d;
        if (delta >= 0){
          Type sqrtdelta = sqrt(delta);
          int sign = (k > 0) ? 1 : -1;    //TODO we can get rid of this if
          for(unsigned i = 0; i < 2; i++) {
            Type x = (- (b+1) - sign * sqrtdelta) / (2 * k);
//          cout<< "Top x = "<< x<< endl;
            if(x > 0 && x > p1.x) {
              p1.x = x;
            }
            sign *= -1;
          }
        }
      }
      else{ //case b and c  take lowest p1.x
        Type delta = (b+1.)*(b+1.) - 4*k*d;
        if (delta >= 0){
              Type sqrtdelta = sqrt(delta);
              int sign = (k > 0) ? -1 : 1;  //this gives us highest p1.x first then lowest.
              for(unsigned i = 0; i < 2; i++) {
                Type x = (- (b+1) - sign * sqrtdelta) / (2 * k);
    //          cout<< "Top x = "<< x<< endl;
                if(x > 0 && x < p1.x) {
                    p1.x = x;
                }
                sign *= -1;
              }
        }
      }

     }
     else{     //concave up k<0 , table 2 case d,e,f
      if(p1.x < p2.x){ //case (d and e) take highest
        Type delta = b*b - 4*k*d;
        if (delta >= 0){
          Type sqrtdelta = sqrt(delta);
          int sign = (k > 0) ? 1 : -1;    //TODO we can get rid of this if
          for(unsigned i = 0; i < 2; i++) {
            Type x = (- b - sign * sqrtdelta) / (2 * k);
//          cout<< "Top x = "<< x<< endl;
            if(x > 0 && x < p2.x) {
              p2.x = x;
            }
            sign *= -1;
          }
        }
      }
      else{  //case f,  x1>x2
        Type delta = (b+1.)*(b+1.) - 4*k*d;
        if (delta >= 0){
          Type sqrtdelta = sqrt(delta);
          int sign = (k > 0) ? -1 : 1;  //this gives us highest x first then lowest.
          for(unsigned i = 0; i < 2; i++) {
            Type x = (- (b+1) - sign * sqrtdelta) / (2 * k);
//          cout<< "Top x = "<< x<< endl;
            if(x > 0 && x > p1.x) {  //highest p2.x
                p1.x = x;
            }
            sign *= -1;
          }
        }
      }
     }
    }
    else if (table == 3){
      if (k<0) {
        do_line = 1;
        Type delta = b*b - 4*k*d;
//         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
        if (delta >= 0){
              Type sqrtdelta = sqrt(delta);
              int sign = (k > 0) ? 1 : -1;
              for(unsigned i = 0; i < 2; i++) {
                Type x = (- b - sign * sqrtdelta) / (2 * k);
    //             cout<< "Top x = "<< x<< endl;
                 if(x < 1 && x < p2.x) {   //lowest p2.x
                    p2.x = x;
                }
                sign *= -1;
              }
            }
      }
    }

    pol2[0] = parabola.k; pol2[1] = parabola.b; pol2[2] = parabola.d;
    if(do_line){ // TODO we donot have to dot them separately. We can do this for each table above

      if (table == 1){
          I1.resize(0);
          I1.resize(1, std::pair<Type, Type>(static_cast<Type>(p1.x), static_cast<Type>(1)));  //not sure if it is taking value. Lets do I1 manually.
          area = trig_integral_A3(m, n, s, a, c, pol2, {{p1.x,static_cast<Type>(1)}}) -  trig_integral_A2(m, n, s, a, c, pol2, {{p1.x,static_cast<Type>(1)}});
        }
      else if (table==2){ //TODO
        if (k>0){
            if (p1.x < p2.x){   //a
              I1.resize(0);
              I1.resize(1, std::pair<Type, Type>(static_cast<Type>(p1.x), static_cast<Type>(p2.x)));
              I3.resize(0);
              I3.resize(1, std::pair<Type, Type>(static_cast<Type>(p2.x), static_cast<Type>(1)));  //not sure if it is taking value. Lets do I1 manually.
              area = trig_integral_A3(m, n, s, a, c, pol2, {{p1.x, p2.x}}) -  trig_integral_A2(m, n, s, a, c, pol2, {{p1.x, p2.x}}) + trig_integral_A3(m, n, s, a, c, pol2, {{p1.x, static_cast<Type>(1)}});
            }
            else{   //b,c
              I1.resize(0);
              I1.resize(1, std::pair<Type, Type>(static_cast<Type>(p2.x), static_cast<Type>(p1.x)));
              I3.resize(0);
              I3.resize(1, std::pair<Type, Type>(static_cast<Type>(0), static_cast<Type>(p2.x)));
              area = trig_integral_A3(m, n, s, a, c, pol2, {{p2.x, p1.x}}) -  trig_integral_A2(m, n, s, a, c, pol2, {{p2.x, p1.x}}) + trig_integral_A3(m, n, s, a, c, pol2, {{ static_cast<Type>(0), p2.x}});
            }
        }
        else {
            if (p1.x<p2.x){   //d,e
              I1.resize(0);
              I1.resize(1, std::pair<Type, Type>(static_cast<Type>(p1.x), static_cast<Type>(p2.x)));
              I3.resize(0);
              I3.resize(1, std::pair<Type, Type>(static_cast<Type>(p2.x), static_cast<Type>(1)));
              area = trig_integral_A3(m, n, s, a, c, pol2, {{p2.x, p1.x}}) -  trig_integral_A2(m, n, s, a, c, pol2, {{p2.x, p1.x}}) + trig_integral_A3(m, n, s, a, c, pol2, {{ static_cast<Type>(0), p2.x}});
            }
            else{  //f
              I1.resize(0);
              I1.resize(1, std::pair<Type, Type>(static_cast<Type>(p2.x), static_cast<Type>(p1.x)));
              I3.resize(0);
              I3.resize(1, std::pair<Type, Type>(static_cast<Type>(0), static_cast<Type>(p2.x)));
              area = trig_integral_A3(m, n, s, a, c, pol2, {{p2.x, p1.x}}) -  trig_integral_A2(m, n, s, a, c, pol2, {{p2.x, p1.x}}) + trig_integral_A3(m, n, s, a, c, pol2, {{ static_cast<Type>(0), p2.x}});
            }
        }

      }

      else if (table == 3){
          I1.resize(0);
          I3.resize(0);
          I1.resize(1, std::pair<Type, Type>(static_cast<Type>(p2.x), static_cast<Type>(1)));
          I3.resize(1, std::pair<Type, Type>(static_cast<Type>(0), static_cast<Type>(p2.x)));
          area = trig_integral_A3(m, n, s, a, c, pol2, I1) -  trig_integral_A2(m, n, s, a, c, pol2, I1) + trig_integral_A3(m, n, s, a, c, pol2, I3);
        }

//       if (table == 1){
//           I1.resize(1, std::pair<Type, Type>(static_cast<Type>(0), static_cast<Type>(p2.x)));
//           I3.resize(1, std::pair<Type, Type>(static_cast<Type>(p2.x), static_cast<Type>(1)));
//           area = trig_integral_A3(m, n, s, a, c, pol2, I1) -  trig_integral_A2(m, n, s, a, c, pol2, I1) + trig_integral_A3(m, n, s, a, c, pol2, I3);
//         }

//       else if (table == 4){
//         I1.resize(1, std::pair<Type, Type>(static_cast<Type>(p1.x), static_cast<Type>(1)));
//         area = trig_integral_A3(m, n, s, a, c, pol2, I1) -  trig_integral_A2(m, n, s, a, c, pol2, I1);
//       }
//
//       else if (table == 6){
//           I1.resize(1, std::pair<Type, Type>(static_cast<Type>(p2.x), static_cast<Type>(1)));
//           I3.resize(1, std::pair<Type, Type>(static_cast<Type>(0), static_cast<Type>(p2.x)));
//           area = trig_integral_A3(m, n, s, a, c, pol2, I1) -  trig_integral_A2(m, n, s, a, c, pol2, I1) + trig_integral_A3(m, n, s, a, c, pol2, I3);
//         }


    }

     else {

        pol1[0] = k+a; pol1[1] = b + c; pol1[2] = d;
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


    return area ;
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
Type trilinier_interpolation(std::vector< std::vector< Type >> & interp_table , const std::vector< Type > &interp_point){

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

    Type x_d = (x-x0)/(x1-x0);
    Type y_d = (y-y0)/(y1-y0);
    Type z_d = (z-z0)/(z1-z0);

  Type c_00 = c_000 * (1-x_d) + c_100 * x_d ;
  Type c_01 = c_001 * (1-x_d) + c_101 * x_d ;
  Type c_10 = c_010 * (1-x_d) + c_110 * x_d ;
  Type c_11 = c_011 * (1-x_d) + c_111 * x_d ;

  Type c_0 = c_00 * (1-y_d) + c_10 * y_d ;
  Type c_1 = c_01 * (1-y_d) + c_11 * y_d ;

  Type cc = c_0 * (1-z_d) + c_1 * z_d ;

  return cc;
}

template <class Type>
void trilinier_interpolation_vector(const std::vector< std::vector< Type >> & interp_table, const std::vector< std::vector< Type >> & interp_table_values  , const std::vector< Type > &interp_point, std::vector< Type > &interp_point_values){

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

  Type x_d = (x-x0)/(x1-x0);
  Type y_d = (y-y0)/(y1-y0);
  Type z_d = (z-z0)/(z1-z0);

  for (unsigned i = 0; i < interp_table_values[0].size(); i++ ){
      Type c_000 = interp_table_values[0][i];
      Type c_001 = interp_table_values[1][i];
      Type c_010 = interp_table_values[2][i];
      Type c_011 = interp_table_values[3][i];
      Type c_100 = interp_table_values[4][i];
      Type c_101 = interp_table_values[5][i];
      Type c_110 = interp_table_values[6][i];
      Type c_111 = interp_table_values[7][i];

      Type c_00 = c_000 * (1-x_d) + c_100 * x_d ;
      Type c_01 = c_001 * (1-x_d) + c_101 * x_d ;
      Type c_10 = c_010 * (1-x_d) + c_110 * x_d ;
      Type c_11 = c_011 * (1-x_d) + c_111 * x_d ;

      Type c_0 = c_00 * (1-y_d) + c_10 * y_d ;
      Type c_1 = c_01 * (1-y_d) + c_11 * y_d ;

      interp_point_values[i] = c_0 * (1-z_d) + c_1 * z_d ;
    }
}


template <class Type>   //TODO change this based on 5 table
void get_p1_p2_p3(const int &table, const std::vector<double> &corner, PointT <Type> &p1, PointT <Type> &p2, PointT <Type> &p3){
    double epsilon = 0.000000000000001;
    Type i1_pm_eps(-1) , i2_pm_eps(-1);

    // std::cout << "Corner " << i << ": (" << corner[0] << ", " << corner[1] << ", " << corner[2] << ") - Print Something\n";

    switch (table) {
        case 0:
            i1_pm_eps = static_cast<Type>(corner[0] - epsilon);
            i2_pm_eps = static_cast<Type>(corner[1] + epsilon);

            p1 = {i1_pm_eps, i1_pm_eps};
            p2 = {i2_pm_eps, i2_pm_eps};
            break;
        case 1:
            i1_pm_eps = static_cast<Type>(corner[0] - epsilon);
            i2_pm_eps = static_cast<Type>(corner[1] + epsilon);
//                    if (i1 == partition ) i1_pm_eps = static_cast<Type>(i1*del_x - epsilon);     //it keeps my i2 in (0,1)
            p1 = {i1_pm_eps, i1_pm_eps};
            p2 = {static_cast<Type>(1), i2_pm_eps};
            break;
        case 2:
            //Do we really need epsilon on this table?
            i1_pm_eps = static_cast<Type>(corner[0] + epsilon);
            i2_pm_eps = static_cast<Type>(corner[1] - epsilon);
            p1 = {i1_pm_eps, i1_pm_eps};
            p2 = {i2_pm_eps, static_cast<Type>(0)};
            break;
        case 3:
            i1_pm_eps = static_cast<Type>(corner[0] + epsilon);
            i2_pm_eps = static_cast<Type>(corner[1] - epsilon);
            p1 = {static_cast<Type>(1), i1_pm_eps};
            p2 = {i2_pm_eps, static_cast<Type>(0)};
            break;
        case 4:
            i1_pm_eps = static_cast<Type>(corner[0] - epsilon);
            i2_pm_eps = static_cast<Type>(corner[1] + epsilon);
            p1 = {i1_pm_eps, static_cast<Type>(0)};
            p2 = {i2_pm_eps, static_cast<Type>(0)};
            break;
            //TODO It is not 3 table it is more like 5 table. fix it.
    }

    p3 = {(p1.x + p2.x)*0.5 , static_cast<Type>(corner[2])};
//       if(fabs(p3.x - p3.y) < epsilon) p3.y = p3.y static_cast<Type>(epsilon);
//     if (fabs(corner[2] - 0) < epsilon ) p3.y = static_cast<Type>(epsilon);
//     if (fabs(corner[2] - 1) < epsilon ) p3.y = static_cast<Type>(1-epsilon);


}


template <class Type>
void find_search_table_trig(const PointT <Type> &p1, const PointT <Type> &p2, const PointT <Type> &p3, unsigned &table_number, Point3D &searchP, PointT <Type> &q1, PointT <Type> &q2, PointT <Type> &q3, bool &vertical){
  double epsilon =0.0000000000001;
  vertical = true;
  bool xSpan = false;
  bool ySpan = false;

  if((p1.x < p3.x && p3.x < p2.x) || ( p2.x > p3.x && p3.x > p1.x)) xSpan = true ;
  if((p1.y < p3.y && p3.y < p2.y) || ( p2.y > p3.y && p3.y > p1.y)) ySpan = true ;
  if(xSpan) {
    if(ySpan) {
      if(fabs(p1.x - p2.x) >= fabs(p1.y - p2.y)) vertical = true;
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

  if(vertical){

      q1 = {(1. - p1.x), p1.y};
      q2 = {(1. - p2.x), p2.y};
      q3 = {(1. - p3.x), p3.y};

      if (fabs(q1.x-q1.y) < epsilon) {
        if (fabs(q2.x - 1.) < epsilon) {
          table_number = 1; searchP.x = static_cast<double>(q1.x); searchP.y = static_cast<double>(q2.y); searchP.z = static_cast<double>(q3.y);}
        else if (fabs(q2.x - q2.y) < epsilon) {
          table_number = 0; searchP.x = static_cast<double>(q1.x); searchP.y = static_cast<double>(q2.x); searchP.z = static_cast<double>(q3.y);}
        else if (fabs(q2.y - 0.) < epsilon) {
          table_number = 2; searchP.x = static_cast<double>(q1.x); searchP.y = static_cast<double>(q2.x); searchP.z = static_cast<double>(q3.y);}
      }
      else if (fabs(q2.x-q2.y) < epsilon) {  // I dont need to do anything for table 0
        if (fabs(q1.x - 1.) < epsilon) {
          table_number = 1; searchP.x = static_cast<double>(q2.x); searchP.y = static_cast<double>(q1.y); searchP.z = static_cast<double>(q3.y);}
        else if (fabs(q1.y - 0.) < epsilon) {
          table_number = 2; searchP.x = static_cast<double>(q2.x); searchP.y = static_cast<double>(q1.x); searchP.z = static_cast<double>(q3.y);}
      }


      else if (fabs(q1.x - 1.) < epsilon){
        if (fabs(q2.y - 0.) < epsilon){
          table_number = 3; searchP.x = static_cast<double>(q1.y); searchP.y = static_cast<double>(q2.x); searchP.z = static_cast<double>(q3.y);}
      }
      else if (fabs(q2.x - 1.) < epsilon){
        if (fabs(q1.y - 0.) < epsilon){
          table_number = 3; searchP.x = static_cast<double>(q2.y); searchP.y = static_cast<double>(q1.x); searchP.z = static_cast<double>(q3.y);}
      }


      else if (fabs(q1.y - 0.) < epsilon){
        if (fabs(q2.y - 0.) < epsilon){table_number = 4; searchP.x = static_cast<double>(q1.x); searchP.y = static_cast<double>(q2.x); searchP.z = static_cast<double>(q3.y
          );}
      }


  }
  else{ //Horizontal

    q1 = {(1. - p1.y), p1.x};
    q2 = {(1. - p2.y), p2.x};
    q3 = {(1. - p3.y), p3.x};

    if (fabs(q1.x-q1.y) < epsilon) {
      if (fabs(q2.y- 0.) < epsilon) {
        table_number = 2; searchP.x = static_cast<double>(q1.x); searchP.y = static_cast<double>(q2.x); searchP.z = static_cast<double>(q3.y);}
      else if (fabs(q2.x- q2.y) < epsilon) {
        table_number = 0; searchP.x = static_cast<double>(q1.x); searchP.y = static_cast<double>(q2.x); searchP.z = static_cast<double>(q3.y);} // No need for table 0 again
      else if (fabs(q2.x- 1.) < epsilon) {
        table_number = 1; searchP.x = static_cast<double>(q1.x); searchP.y = static_cast<double>(q2.y); searchP.z = static_cast<double>(q3.y);}
    }
    else if (fabs(q2.x-q2.y) < epsilon) {
      if (fabs(q1.y- 0.) < epsilon) {
        table_number = 2; searchP.x = static_cast<double>(q2.x); searchP.y = static_cast<double>(q1.x); searchP.z = static_cast<double>(q3.y);}
      else if (fabs(q1.x- 1.) < epsilon) {
        table_number = 1; searchP.x = static_cast<double>(q2.x); searchP.y = static_cast<double>(q1.y); searchP.z = static_cast<double>(q3.y);}
    }


    else if (fabs(q1.x - 1.) < epsilon){
        if (fabs(q2.y - 0.) < epsilon){
          table_number = 3; searchP.x = static_cast<double>(q1.y); searchP.y = static_cast<double>(q2.x); searchP.z = static_cast<double>(q3.y);}
    }
    else if (fabs(q2.x - 1.) < epsilon){
        if (fabs(q1.y - 0.) < epsilon){
          table_number = 3; searchP.x = static_cast<double>(q2.y); searchP.y = static_cast<double>(q1.x); searchP.z = static_cast<double>(q3.y);}
    }


    else if (fabs(q1.x - 0.) < epsilon){
        if (fabs(q2.x - 0.) < epsilon){
          table_number = 4; searchP.x = static_cast<double>(q1.x); searchP.y = static_cast<double>(q2.x); searchP.z = static_cast<double>(q3.y);}
    }
  }
}




double GaussIntegral(const int &xExp, const int &yExp, const double* xg, const double* yg, const std::vector<double> &interp_point_weights, const double* gaussWeight){
  double Integral = 0;
  for(unsigned ig = 0; ig < interp_point_weights.size(); ig++) {
    Integral += pow(xg[ig],xExp) * pow(yg[ig],yExp) * interp_point_weights[ig] * gaussWeight[ig];
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
        if (corners.size() != 8) {
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

        for (size_t i = 0; i < corners.size(); ++i) {
            const auto& corner = corners[i];
            std::vector<double> corner_vec = {corner.x, corner.y, corner.z};
            get_p1_p2_p3(table, corner_vec, p1, p2, p3);

            int count = 0;
            for (unsigned qq = 0; qq <= qM; qq++) {
                for (unsigned jj = 0; jj <= qq; jj++) {
                    unsigned ii = qq - jj;
                    area = find_trig_area_2intersection_formula(ii, jj, s, a, c, table, p1, p2, p3);
                    cornerAreas[i].push_back(static_cast<double>(area));
                    count++;
                }
            }
            (*_Pweights)(s, a, c, table, p1, p2, p3, cornerWeights[i]);
        }
    }

    void subdivideWithRelativeError(int maxDepth, double maxRelativeError, int currentDepth = 0) {
        if (currentDepth >= maxDepth || !isLeaf) {
            getCorners();
            return;
        }

        getCorners();

        // Calculate midpoints for subdivision
        std::vector<Point3D> midpoints(19);
        calculateMidpoints(midpoints);

        std::vector<double> relativeErrors;
        std::vector<double> relativeErrorsOpposite;

        for (const auto& midpoint : midpoints) {
            std::vector<double> interp_point = {midpoint.x, midpoint.y, midpoint.z};

            Type f_area(0);
            Type c = 1;
            PointT<Type> p1, p2, p3;
            get_p1_p2_p3(table, interp_point, p1, p2, p3);

            std::vector<std::vector<double>> interpolation_vector(8);
            int count = 0;
            for (unsigned qq = 0; qq <= qM; qq++) {
                for (unsigned jj = 0; jj <= qq; jj++) {
                    unsigned ii = qq - jj;
                    for (size_t ic = 0; ic < corners.size(); ++ic) {
                        interpolation_vector[ic] = {corners[ic].x, corners[ic].y, corners[ic].z, cornerAreas[ic][count]};
                    }
                    double interp_area = trilinier_interpolation(interpolation_vector, interp_point);
                    f_area = find_trig_area_2intersection_formula(jj, ii, s, a, c, table, p1, p2, p3);
                    double formula_area = static_cast<double>(f_area);
                    double r_error = fabs(formula_area - interp_area) / formula_area;
                    double r_error_opposite = fabs(formula_area - interp_area) / (1.0 / (ii + 1) * (jj + 1) - formula_area);

                    relativeErrors.push_back(r_error);
                    relativeErrorsOpposite.push_back(r_error_opposite);

                    count++;
                }
            }
        }
        relative_error = *std::max_element(relativeErrors.begin(), relativeErrors.end());
        relative_error_opposite = *std::max_element(relativeErrorsOpposite.begin(), relativeErrorsOpposite.end());

        if (depth <= 3 || relative_error > maxRelativeError || relative_error_opposite > maxRelativeError) {
            isLeaf = false;
            children.reserve(children.size() + 8);
            std::vector<std::vector<Point3D>> childCorners = subdivideCorners();
            for (const auto& childCorner : childCorners) {
                children.emplace_back(childCorner, table, depth + 1, qM, _Pweights);
            }

            for (auto& child : children) {
                child.subdivideWithRelativeError(maxDepth, maxRelativeError, currentDepth + 1);
            }
        }
    }

//     OctreeNode* search(const Point3D& point) {
//         if (isLeaf) {
//             return this;
//         }
//
//         for (auto& child : children) {
//             if (child.contains(point)) {
//                 return child.search(point);
//             }
//         }
//
//         return nullptr; // Should not reach here under normal circumstances
//     }



    OctreeNode* search(const Point3D& point) {
      // First check if point is even in this node
      if (contains(point)) {
        std::cout << "\nFound containing node at depth " << depth << ":\n";
        std::cout << "Point: (" << point.x << ", " << point.y << ", " << point.z << ")\n";
        std::cout << "Node corners:\n";
        for (size_t i = 0; i < corners.size(); ++i) {
            std::cout << "Corner " << i << ": ("
                     << corners[i].x << ", "
                     << corners[i].y << ", "
                     << corners[i].z << ")\n";
        }

        // If this is a leaf node, we're done
        if (isLeaf) {
            std::cout << "This is a leaf node - returning\n";
            return this;
        }

        // If not a leaf, check children
        std::cout << "Checking " << children.size() << " children\n";
        for (auto& child : children) {
            OctreeNode* result = child.search(point);
            if (result != nullptr) {
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


//         OctreeNode* search(const Point3D& point) {
//             // First check if point is even in this node
//             if (!contains(point)) {
//                 return nullptr;
//             }
//
//             // If this is a leaf node and contains the point, return this node
//             if (isLeaf) {
//                 return this;
//             }
//
//             // Check each child
//             for (auto& child : children) {
//                 OctreeNode* result = child.search(point);
//                 if (result != nullptr) {
//                     return result;
//                 }
//             }
//
//             // If no child contains the point (shouldn't happen with proper subdivision)
//             // return this node as the best match
//             return this;
//         }


    void saveOctreeToCSV(const std::string& filename) const {
        std::ofstream ofs(filename, std::ios::binary);
        if (!ofs.is_open()) {
            std::cerr << "Error: Unable to open file for writing: " << filename << std::endl;
            return;
        }
        serialize(ofs);
        ofs.close();
    }

    void loadOctreeFromCSV(const std::string& filename) {
        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs.is_open()) {
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

      serializeVector(ofs, corners);
      serializeVector(ofs, cornerAreas);
      serializeVector(ofs, cornerWeights);

      size_t childCount = children.size();
      ofs.write(reinterpret_cast<const char*>(&childCount), sizeof(childCount));
      for(const auto& child : children) {
        child.serialize(ofs);
      }
    }

    void deserialize(std::ifstream& ifs) {
      ifs.read(reinterpret_cast<char*>(&isLeaf), sizeof(isLeaf));
      ifs.read(reinterpret_cast<char*>(&depth), sizeof(depth));

      deserializeVector(ifs, corners);
      deserializeVector(ifs, cornerAreas);
      deserializeVector(ifs, cornerWeights);

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
        for (const auto& point : vec) {
            ofs.write(reinterpret_cast<const char*>(&point.x), sizeof(double));
            ofs.write(reinterpret_cast<const char*>(&point.y), sizeof(double));
            ofs.write(reinterpret_cast<const char*>(&point.z), sizeof(double));
        }
    }

    void deserializeVector(std::ifstream& ifs, std::vector<Point3D>& vec) {
        size_t size;
        ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
        vec.resize(size);
        for (size_t i = 0; i < size; ++i) {
            ifs.read(reinterpret_cast<char*>(&vec[i].x), sizeof(double));
            ifs.read(reinterpret_cast<char*>(&vec[i].y), sizeof(double));
            ifs.read(reinterpret_cast<char*>(&vec[i].z), sizeof(double));
        }
    }



      Point3D middleof(Point3D& point1, Point3D& point2){
        Point3D midpoint{
                (point1.x + point2.x) / 2.,
                (point1.y + point2.y) / 2.,
                (point1.z + point2.z) / 2.
                };
        return midpoint;
      }

    void calculateMidpoints(std::vector<Point3D>& midpoints) {

      midpoints.resize(19);
                                                      // TODO may need to adjust based on your specific requirements

      //midpoint of EDGEs
        midpoints[0] = middleof(corners[0], corners[1]);
        midpoints[1] = middleof(corners[2], corners[3]);
        midpoints[2] = middleof(corners[4], corners[5]);
        midpoints[3] = middleof(corners[6], corners[7]);

        midpoints[4] = middleof(corners[0], corners[2]);
        midpoints[5] = middleof(corners[1], corners[3]);
        midpoints[6] = middleof(corners[4], corners[6]);
        midpoints[7] = middleof(corners[5], corners[7]);

        midpoints[8] = middleof(corners[0], corners[4]);
        midpoints[9] = middleof(corners[1], corners[5]);
        midpoints[10] = middleof(corners[2], corners[6]);
        midpoints[11] = middleof(corners[3], corners[7]);

      //midpoint of faces. I am just considering taking the mid point on one side. I should same if I consider all four corners.

        midpoints[12] = middleof(midpoints[0], midpoints[1]);
        midpoints[13] = middleof(midpoints[0], midpoints[2]);
        midpoints[14] = middleof(midpoints[1], midpoints[3]);
        midpoints[15] = middleof(midpoints[2], midpoints[3]);
        midpoints[16] = middleof(midpoints[4], midpoints[6]);
        midpoints[17] = middleof(midpoints[5], midpoints[7]);

      //inside point
        midpoints[18] = middleof(midpoints[16], midpoints[17]);
    }


    std::vector<std::vector<Point3D>> subdivideCorners() {
        std::vector<std::vector<Point3D>> childCorners(8, std::vector<Point3D>(8));
        std::vector<Point3D> midpoints(19);
        calculateMidpoints(midpoints);

        // Define child corners based on the original corners and midpoints
        // This is a simplified version; you may need to adjust based on your specific subdivision strategy
        childCorners[0] = {corners[0], midpoints[0], midpoints[4], midpoints[12], midpoints[8], midpoints[13], midpoints[16], midpoints[18]};
        childCorners[1] = {midpoints[0], corners[1], midpoints[12], midpoints[5], midpoints[13], midpoints[9], midpoints[18], midpoints[17]};
        childCorners[2] = {midpoints[4], midpoints[12], corners[2], midpoints[1], midpoints[16], midpoints[18], midpoints[10], midpoints[14]};
        childCorners[3] = {midpoints[12], midpoints[5], midpoints[1], corners[3], midpoints[18], midpoints[17], midpoints[14], midpoints[11]};
        childCorners[4] = {midpoints[8], midpoints[13], midpoints[16], midpoints[18], corners[4], midpoints[2], midpoints[6], midpoints[15]};
        childCorners[5] = {midpoints[13], midpoints[9], midpoints[18], midpoints[17], midpoints[2], corners[5], midpoints[15], midpoints[7]};
        childCorners[6] = {midpoints[16], midpoints[18], midpoints[10], midpoints[14], midpoints[6], midpoints[15], corners[6], midpoints[3]};
        childCorners[7] = {midpoints[18], midpoints[17], midpoints[14], midpoints[11], midpoints[15], midpoints[7], midpoints[3], corners[7]};

        return childCorners;
    }



//     bool contains(const Point3D& point) const {// uses tetrahedra decomposition
//     const double EPSILON = 1e-10;
//
//     // Helper function to compute determinant of 3x3 matrix
//     auto det3x3 = [](double a00, double a01, double a02,
//                      double a10, double a11, double a12,
//                      double a20, double a21, double a22) -> double {
//         return a00 * (a11 * a22 - a12 * a21) -
//                a01 * (a10 * a22 - a12 * a20) +
//                a02 * (a10 * a21 - a11 * a20);
//     };
//
//     // Helper function to check if point is inside tetrahedron using barycentric coordinates
//     auto pointInTetrahedron = [EPSILON, &det3x3](const Point3D& p,
//                                                 const Point3D& a,
//                                                 const Point3D& b,
//                                                 const Point3D& c,
//                                                 const Point3D& d) -> bool {
//         double d0 = det3x3(a.x-d.x, b.x-d.x, c.x-d.x,
//                           a.y-d.y, b.y-d.y, c.y-d.y,
//                           a.z-d.z, b.z-d.z, c.z-d.z);
//
//         double d1 = det3x3(p.x-d.x, b.x-d.x, c.x-d.x,
//                           p.y-d.y, b.y-d.y, c.y-d.y,
//                           p.z-d.z, b.z-d.z, c.z-d.z);
//
//         double d2 = det3x3(a.x-d.x, p.x-d.x, c.x-d.x,
//                           a.y-d.y, p.y-d.y, c.y-d.y,
//                           a.z-d.z, p.z-d.z, c.z-d.z);
//
//         double d3 = det3x3(a.x-d.x, b.x-d.x, p.x-d.x,
//                           a.y-d.y, b.y-d.y, p.y-d.y,
//                           a.z-d.z, b.z-d.z, p.z-d.z);
//
//         if (std::abs(d0) < EPSILON) return false;
//
//         double b1 = d1 / d0;
//         double b2 = d2 / d0;
//         double b3 = d3 / d0;
//         double b4 = 1.0 - b1 - b2 - b3;
//
//         return b1 >= -EPSILON && b2 >= -EPSILON && b3 >= -EPSILON && b4 >= -EPSILON;
//     };
//
//     // Decompose hexahedron into five tetrahedra
//     return pointInTetrahedron(point, corners[0], corners[1], corners[2], corners[5]) ||
//            pointInTetrahedron(point, corners[0], corners[2], corners[5], corners[7]) ||
//            pointInTetrahedron(point, corners[2], corners[5], corners[7], corners[6]) ||
//            pointInTetrahedron(point, corners[0], corners[2], corners[3], corners[7]) ||
//            pointInTetrahedron(point, corners[0], corners[4], corners[5], corners[7]);
//     }


        bool contains(const Point3D& point) const {
        const double EPSILON = 1e-10;

        // Helper function to calculate dot product
        auto dot = [](const Point3D& a, const Point3D& b) -> double {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        };

        // Helper function to calculate cross product
        auto cross = [](const Point3D& a, const Point3D& b) -> Point3D {
            return Point3D{
                a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x
            };
        };

        // Helper function to create vector from two points
        auto makeVector = [](const Point3D& from, const Point3D& to) -> Point3D {
            return Point3D{to.x - from.x, to.y - from.y, to.z - from.z};
        };

        // Define the six faces of the hexahedron
        // Each face is defined by four corners in counter-clockwise order
        const std::vector<std::vector<int>> faces = {
            {0, 1, 3, 2},    // front face
            {4, 6, 7, 5},    // back face
            {0, 4, 5, 1},    // bottom face
            {2, 3, 7, 6},    // top face
            {0, 2, 6, 4},    // left face
            {1, 5, 7, 3}     // right face
        };

        // Check if point is on the correct side of all faces
        for (const auto& face : faces) {
            // Get three points from the face to define the plane
            const Point3D& v0 = corners[face[0]];
            const Point3D& v1 = corners[face[1]];
            const Point3D& v2 = corners[face[2]];

            // Calculate face normal using cross product
            Point3D edge1 = makeVector(v0, v1);
            Point3D edge2 = makeVector(v0, v2);
            Point3D normal = cross(edge1, edge2);

            // Calculate signed distance from point to plane
            Point3D toPoint = makeVector(v0, point);
            double signedDist = dot(normal, toPoint);

            // Point must be on the "inside" side of each face
            // For a properly oriented hexahedron, this should be negative
            if (signedDist > EPSILON) {
                return false;
            }
        }

        return true;
    }


//     bool contains(const Point3D& point) const { //uses faces
//         // Check if the point is inside the bounding box defined by the corners
//         Point3D minBound = corners[0];
//         Point3D maxBound = corners[7];
//         for (const auto& corner : corners) {
//             minBound.x = std::min(minBound.x, corner.x);
//             minBound.y = std::min(minBound.y, corner.y);
//             minBound.z = std::min(minBound.z, corner.z);
//             maxBound.x = std::max(maxBound.x, corner.x);
//             maxBound.y = std::max(maxBound.y, corner.y);
//             maxBound.z = std::max(maxBound.z, corner.z);
//         }
//         return point.x >= minBound.x && point.x <= maxBound.x &&
//                point.y >= minBound.y && point.y <= maxBound.y &&
//                point.z >= minBound.z && point.z <= maxBound.z;
//     }

//     bool contains(const Point3D& point) const {
//         We'll use the sign of triple products to determine if the point is on the
//         correct side of each face. A point is inside if it's on the correct side
//         of all faces.
//
//         const double EPSILON = 1e-10;  // Add small epsilon for numerical stability
//
//         Helper function to calculate triple product
//         auto tripleProduct = [](const Point3D& a, const Point3D& b, const Point3D& c) -> double {
//             return (a.x * (b.y * c.z - b.z * c.y) +
//                     a.y * (b.z * c.x - b.x * c.z) +
//                     a.z * (b.x * c.y - b.y * c.x));
//         };
//
//         Helper function to create vector from two points
//         auto makeVector = [](const Point3D& from, const Point3D& to) -> Point3D {
//             return Point3D{to.x - from.x, to.y - from.y, to.z - from.z};
//         };
//
//         Define faces based on your corner ordering
//         Each face is defined by three corners in counter-clockwise order when viewed from outside
//         const std::vector<std::vector<int>> faces = {
//             {0, 1, 3},  // bottom face
//             {0, 2, 1},  // bottom face (second triangle)
//             {4, 7, 5},  // top face
//             {4, 6, 7},  // top face (second triangle)
//             {0, 4, 5},  // front face
//             {0, 5, 1},  // front face (second triangle)
//             {2, 6, 7},  // back face
//             {2, 7, 3},  // back face (second triangle)
//             {0, 3, 7},  // left face
//             {0, 7, 4},  // left face (second triangle)
//             {1, 5, 7},  // right face
//             {1, 7, 3}   // right face (second triangle)
//         };
//
//         for (const auto& face : faces) {
//             Get three points defining the face
//             const Point3D& v0 = corners[face[0]];
//             const Point3D& v1 = corners[face[1]];
//             const Point3D& v2 = corners[face[2]];
//
//             Create vectors
//             Point3D edge1 = makeVector(v0, v1);
//             Point3D edge2 = makeVector(v0, v2);
//             Point3D toPoint = makeVector(v0, point);
//
//             Calculate triple product
//             double tripleP = tripleProduct(edge1, edge2, toPoint);
//
//             If point is on wrong side of any face (with epsilon tolerance), it's outside
//             if (tripleP < -EPSILON) {
//                 return false;
//             }
//         }
//
//         return true;
//     }
};

template <class Type>
void generateAndLoadOctrees(const int &maxDepth, const int &degree, const double &percent, CutFemWeightParabola<double, Type> &Pweights, std::vector<OctreeNode<Type>>& loadedRoots) {


    std::vector<std::vector<Point3D>> initialCorners = {
      {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 1.0, 0.5}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.5}, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
      {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.5}, {0.0, 1.0, 0.0}, {0.0, 1.0, 0.5}, {1.0, 0.0, 0.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
      {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 1.0, 0.5}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.5}, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
      {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.5}, {0.0, 1.0, 0.0}, {0.0, 1.0, 1.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.5}, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
      {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 1.0, 0.5}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.5}, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}},
    };
    for (int ttable = 0; ttable < 5; ++ttable) {
        std::string filename = "save/octree_table_" + std::to_string(ttable) + "_maxdepth_" + std::to_string(maxDepth) + "_per_" + std::to_string(percent) + "_degree_" + std::to_string(degree) + ".csv";

        FILE *fp;
        fp = fopen(filename.c_str(), "r");
        if (fp != NULL) {
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
    for (int ttable = 0; ttable < 5; ++ttable) {

        loadedRoots.emplace_back(initialCorners[ttable], ttable, 0, degree, nullptr);
        loadedRoots.back().loadOctreeFromCSV("save/octree_table_" + std::to_string(ttable) + "_maxdepth_" + std::to_string(maxDepth) + "_per_" + std::to_string(percent) + "_degree_" + std::to_string(degree) + ".csv");
    }
}


// template <class Type>
// void printOctreeStructure(OctreeNode<Type> node, int depth = 0) {
//     for (int i = 0; i < depth; ++i) {
//         std::cout << "  ";
//     }
//     std::cout << "Node Corners: \n";
//     for (int i = 0; i < 8; ++i) {
//         for (int j = 0; j < depth + 1; ++j) {
//             std::cout << "  ";
//         }
//         std::cout << "(" << node->corners[i].x << ", " << node->corners[i].y << ", " << node->corners[i].z << ")\n";
//     }
//
//     if (node->isLeaf) {
//         std::cout << "  Relative error = " << node->relative_error << " depth = " << node->depth << " [Leaf]\n";
//         std::cout << "  Corner Areas and Weights:\n";
//         for (int j = 0; j < 8; ++j) {
//             for (int i = 0; i < depth + 1; ++i) {
//                 std::cout << "  ";
//             }
//             std::cout << "Corner " << j << ":\n";
//             for (int i = 0; i < depth + 2; ++i) {
//                 std::cout << "  ";
//             }
//             std::cout << "Areas = ";
//             for (const auto& area : node->cornerAreas[j]) {
//                 std::cout << area << " ";
//             }
//             std::cout << "\n";
//             for (int i = 0; i < depth + 2; ++i) {
//                 std::cout << "  ";
//             }
//             std::cout << "Weights = ";
//             for (const auto& weight : node->cornerWeights[j]) {
//                 std::cout << weight << " ";
//             }
//             std::cout << "\n";
//         }
//     } else {
//         std::cout << "  Relative error = " << node->relative_error << " depth = " << node->depth << " [Non-Leaf]\n";
//         for (auto& child : children) {
//             printOctreeStructure(child, depth + 1);
//         }
//     }
// }




template <class Type>
void printOctreeStructure(const OctreeNode<Type>& node, int depth = 0) {
    // Print indentation based on depth
    std::string indent(depth * 2, ' ');

    // Print node information
    std::cout << indent << "Node at depth " << node.depth << ": ";
    std::cout << indent << "  Is Leaf: " << (node.isLeaf ? "Yes" : "No") << " \n";
//     std::cout << indent << "  Table: " << node.table << "\n";
//     std::cout << indent << "  Relative Error: " << node.relative_error << "\n";
//     std::cout << indent << "  Relative Error Opposite: " << node.relative_error_opposite << "\n";

    // Print corner information
//     std::cout << indent << "  Corners:\n";
    for (size_t i = 0; i < node.corners.size(); ++i) {
        std::cout << indent << "    Corner " << i << ": ("
                  << node.corners[i].x << ", "
                  << node.corners[i].y << ", "
                  << node.corners[i].z << ")\n";
    }

    // If not a leaf, print children
    if (!node.isLeaf) {
//         std::cout << indent << "  Children:\n";
        for (size_t i = 0; i < node.children.size(); ++i) {
            std::cout << indent << "  Child " << i << ":\n";
            printOctreeStructure(node.children[i], depth + 1);
        }
    }
}

int main() {
  unsigned int m = 0;
  unsigned int n = 0;
  int s = 0;

  std::cout.precision(16);

  typedef cpp_bin_float_oct Type;
  Type k, b, d, a = 0, c = 1, area1, area2, easy_area1, easy_area2,Trig_area1,Trig_area2;
  double Area0 = 0, Area = 0, Ix = 0, Iy = 0,Ixy =0, Ix3 = 0, Ix2y = 0, Ixy2 = 0, Iy3 = 0, Ix2y2 = 0;

  std::vector <Type> pol1(3, 0);
  std::vector <Type> pol2(3, 0);
  clock_t t = clock();

  PointT <Type> p1, p2, p3;  // points in domain D

  //table 0
//   p1 = { 0.7, 0.3 };
//   p2 = { 0.2, 0.8 };
//   p3 = { (p1.x + p2.x) / 2.0, 0.2 };
//   p1 = { static_cast<Type>(0.7), static_cast<Type>(0.3) };
//   p2 = { static_cast<Type>(0.2), static_cast<Type>(0.8) };
//   p3 = { static_cast<Type>((p1.x + p2.x) / 2.0), static_cast<Type>(0.2) };

//
//


//table 1 vertical singleintersection
  p1 = { static_cast<Type>(0), static_cast<Type>(0.3) };
  p2 = { static_cast<Type>(0.2), static_cast<Type>(0.8) };
  p3 = { static_cast<Type>((p1.x + p2.x) / 2.0), static_cast<Type>(0.2) };


// // table 2 vertical multipleintersection. It has 5 different cases.
//     p1 = { static_cast<Type>(.6), static_cast<Type>(0.4) };
//     p2 = { static_cast<Type>(0.2), static_cast<Type>(0.) };
//     p3 = { (p1.x+p2.x)/2., static_cast<Type>(0.55) };



//table 3 vertical singleintersection
//   p1 = { static_cast<Type>(0.), static_cast<Type>(0.6) };
//   p2 = { static_cast<Type>(0.6), static_cast<Type>(0.) };
//   p3 = { static_cast<Type>((p1.x + p2.x) / 2.0), static_cast<Type>(0.4) };



//table 4 vertical singleintersection
//   p1 = { static_cast<Type>(0.3), static_cast<Type>(0.) };
//   p2 = { static_cast<Type>(0.7), static_cast<Type>(0.) };
//   p3 = { static_cast<Type>((p1.x + p2.x) / 2.0), static_cast<Type>(0.4) };



  //test a horizontal table 1 ;
//   p1 = { static_cast<Type>(0.6), static_cast<Type>(0) };
//   p2 = { static_cast<Type>(0.6), static_cast<Type>(0.4) };
//   p3 = { static_cast<Type>(0.55), static_cast<Type>((p1.y + p2.y) / 2.0) };

//     p1 = { static_cast<Type>(0.), static_cast<Type>(0.6) };
//   p2 = { static_cast<Type>(0.4), static_cast<Type>(0.6) };
//   p3 = { static_cast<Type>((p1.x + p2.x) / 2.0), static_cast<Type>(0.55) };


  PointT <Type> q1, q2, q3;   // points in domain D* or D** depending on vertical

  bool vertical;
  unsigned table_number;
  Point3D searchP;

  find_search_table_trig<Type>(p1, p2, p3, table_number, searchP, q1, q2, q3, vertical);


//   q1 = { static_cast<Type>(1.0-0.7), static_cast<Type>(0.3) };
//   q2 = { static_cast<Type>(1.0-0.2), static_cast<Type>(0.8) };
//   q3 = { static_cast<Type>((q1.x + q2.x) / 2.0), static_cast<Type>(0.2) };


//   p1 = { static_cast<Type>(p1.x), static_cast<Type>(p1.y) };
//   p2 = { static_cast<Type>(p2.x), static_cast<Type>(p2.y) };
//   p3 = { static_cast<Type>(p3.x), static_cast<Type>(p3.y) };

  cout<<" table = "<< table_number <<". q = ("<<q1.x<<","<<q1.y<<"), ("<<q2.x<<","<<q2.y<<"), ("<<q3.x<<","<<q3.y<<") " << " searchp = (" << searchP.x << "," << searchP.y <<"," << searchP.z <<") " << "vertical =" << vertical << endl;



  Parabola <Type> parabola = get_parabola_equation(p1,p2,p3);
  std::cout<< "parabola " << parabola.k<<"x^2+"<<parabola.b<<"x+" << parabola.d << " + y = 0 " <<std::endl;

  std::vector<double>weightCF, interp_point_weights, interp_point_integrals;

    CutFemWeightParabola <double, Type> Pweights(TRI, 3, "legendre");
    Pweights(s, a, c, table_number, q1, q2, q3, weightCF);                  // WeightCF is always calculated in Domain D*
//     Pweights(s, a, c, 0, p1, p2, p3, weightCF);

    cout<< " cutfem weight = ";
    for (size_t j = 0; j < weightCF.size(); ++j){
      std::cout << weightCF[j] << ", ";
    }
    cout<<endl;

    const double* gaussWeight =  Pweights.GetGaussWeightPointer();
    const double* xg = Pweights.GetGaussCoordinatePointer(0);
    const double* yg = Pweights.GetGaussCoordinatePointer(1);

    cout<< " gause points = ";
    for (size_t j = 0; j < weightCF.size(); ++j){
      std::cout << "("<<xg[j] << ", " << yg[j] << "), ";
    }
    cout<<endl;
    cout << xg[0] << " " <<yg[0]  ;

      if(vertical){
        Area = GaussIntegral(0, 0, xg, yg, weightCF, gaussWeight);
        Ix   = GaussIntegral(1, 0, xg, yg, weightCF, gaussWeight);
        Iy   = GaussIntegral(0, 1, xg, yg, weightCF, gaussWeight);
        Ixy  = GaussIntegral(1, 1, xg, yg, weightCF, gaussWeight);
        Ix3  = GaussIntegral(3, 0, xg, yg, weightCF, gaussWeight);
        Ix2y = GaussIntegral(2, 1, xg, yg, weightCF, gaussWeight);
        Ixy2 = GaussIntegral(1, 2, xg, yg, weightCF, gaussWeight);
        Iy3  = GaussIntegral(0, 3, xg, yg, weightCF, gaussWeight);
        Ix2y2= GaussIntegral(2, 2, xg, yg, weightCF, gaussWeight);
      }
      else{
        Area = GaussIntegral(0, 0, yg, xg, weightCF, gaussWeight);
        Ix   = GaussIntegral(1, 0, yg, xg, weightCF, gaussWeight);
        Iy   = GaussIntegral(0, 1, yg, xg, weightCF, gaussWeight);
        Ixy  = GaussIntegral(1, 1, yg, xg, weightCF, gaussWeight);
        Ix3  = GaussIntegral(3, 0, yg, xg, weightCF, gaussWeight);
        Ix2y = GaussIntegral(2, 1, yg, xg, weightCF, gaussWeight);
        Ixy2 = GaussIntegral(1, 2, yg, xg, weightCF, gaussWeight);
        Iy3  = GaussIntegral(0, 3, yg, xg, weightCF, gaussWeight);
        Ix2y2= GaussIntegral(2, 2, yg, xg, weightCF, gaussWeight);
      }

    std::cout << "Area0 = " << Area0 << std::endl;
    std::cout << "Area = " << Area << std::endl;
    std::cout << "Ix = " << Ix << std::endl;
    std::cout << "Iy = " << Iy << std::endl;
    std::cout << "Ixy = " << Ixy << std::endl;
    std::cout << "Ix3 = " << Ix3 << std::endl;
    std::cout << "Ix2y = " << Ix2y << std::endl;
    std::cout << "Ixy2 = " << Ixy2 << std::endl;
    std::cout << "Iy3 = " << Iy3 << std::endl;
    std::cout << "Ix2y2 = " << Ix2y2 << std::endl;

    int maxDepth = 5;
    int degree = 3;
    double percent = 0.001;
  //   std::vector<OctreeNode<Type>> roots;
    std::vector<OctreeNode<Type>>loadedRoots;

    generateAndLoadOctrees<Type>(maxDepth, degree, percent, Pweights, loadedRoots);
//     printOctreeStructure(loadedRoots[0]);
//
//     return 1;

//     Point3D originalPoint(0.7,0.2,0.2);
//     Point3D searchP(1. - originalPoint.x , 1. - originalPoint.y, originalPoint.z );
    std::cout << "\nSearch Point: (" << searchP.x << ", " << searchP.y << ", " << searchP.z << ")\n";
    OctreeNode<Type>* result = loadedRoots[table_number].search(searchP);
    cout<< " results =" << result ;
    if(result) {
      std::cout << "Found the smallest sub-cube containing the search point." << std::endl;
      std::cout << "\nSearch Point: (" << searchP.x << ", " << searchP.y << ", " << searchP.z << ")\n";
      std::cout << "Smallest Sub-cube Bounds: ";
      std::cout << "depth : = " << result->depth << " \n";

      for (size_t i = 0; i < result->corners.size(); ++i) {
        std::cout << "    Corner " << i << ": ("
              << result->corners[i].x << ", "
              << result->corners[i].y << ", "
              << result->corners[i].z << ")\n";
      }
//           for (size_t i = 0; i < result->cornerAreas.size(); ++i) {
//             std::cout << "    Corner " << i << " Areas : (" ;
//             for (size_t j = 0; j < result->cornerAreas[i].size(); ++j){
//               std::cout << result->cornerAreas[i][j] << ", ";
//             }
//             std::cout << " )"<<std::endl;
//           }
//
//           for (size_t i = 0; i < result->cornerWeights.size(); ++i) {
//             std::cout << "    Corner " << i << " Weights : (" ;
//             for (size_t j = 0; j < result->cornerWeights[i].size(); ++j){
//               std::cout << result->cornerWeights[i][j] << ", ";
//             }
//             std::cout << " )"<<std::endl;
//           }

      std::vector<double>interp_point = {searchP.x, searchP.y, searchP.z};
      std::vector<std::vector<double>> corners(8, std::vector<double>(3));  // A 2D vector of size 8x3
      for (size_t i = 0; i < result->corners.size(); ++i) {
          corners[i][0] = result->corners[i].x;  // x-coordinate
          corners[i][1] = result->corners[i].y;  // y-coordinate
          corners[i][2] = result->corners[i].z;  // z-coordinate
      }

      std::cout << "\n interp Point: (" << interp_point[0] << ", " << interp_point[1] << ", " << interp_point[2] << ")\n";
      trilinier_interpolation_vector(corners, result->cornerAreas, interp_point, interp_point_weights);
      std::cout << " interpolated integrals = ";
      for (size_t j = 0; j < interp_point_weights.size(); ++j){
          std::cout << interp_point_weights[j] << ", ";
      }
      std::cout << " )"<<std::endl;


      trilinier_interpolation_vector(corners, result->cornerWeights, interp_point, interp_point_weights);
      std::cout << " interpolated weights = ";
      for (size_t j = 0; j < interp_point_weights.size(); ++j){
        std::cout << interp_point_weights[j] << ", ";
      }
      std::cout << " )"<<std::endl;

      if(vertical){
        Area = GaussIntegral(0, 0, xg, yg, interp_point_weights, gaussWeight);
        Ix   = GaussIntegral(1, 0, xg, yg, interp_point_weights, gaussWeight);
        Iy   = GaussIntegral(0, 1, xg, yg, interp_point_weights, gaussWeight);
        Ixy  = GaussIntegral(1, 1, xg, yg, interp_point_weights, gaussWeight);
        Ix3  = GaussIntegral(3, 0, xg, yg, interp_point_weights, gaussWeight);
        Ix2y = GaussIntegral(2, 1, xg, yg, interp_point_weights, gaussWeight);
        Ixy2 = GaussIntegral(1, 2, xg, yg, interp_point_weights, gaussWeight);
        Iy3  = GaussIntegral(0, 3, xg, yg, interp_point_weights, gaussWeight);
        Ix2y2= GaussIntegral(2, 2, xg, yg, interp_point_weights, gaussWeight);
      }
      else{
        Area = GaussIntegral(0, 0, yg, xg, interp_point_weights, gaussWeight);
        Ix   = GaussIntegral(1, 0, yg, xg, interp_point_weights, gaussWeight);
        Iy   = GaussIntegral(0, 1, yg, xg, interp_point_weights, gaussWeight);
        Ixy  = GaussIntegral(1, 1, yg, xg, interp_point_weights, gaussWeight);
        Ix3  = GaussIntegral(3, 0, yg, xg, interp_point_weights, gaussWeight);
        Ix2y = GaussIntegral(2, 1, yg, xg, interp_point_weights, gaussWeight);
        Ixy2 = GaussIntegral(1, 2, yg, xg, interp_point_weights, gaussWeight);
        Iy3  = GaussIntegral(0, 3, yg, xg, interp_point_weights, gaussWeight);
        Ix2y2= GaussIntegral(2, 2, yg, xg, interp_point_weights, gaussWeight);
      }

      std::cout << "Area0 = " << Area0 << std::endl;
      std::cout << "Area = " << Area << std::endl;
      std::cout << "Ix = " << Ix << std::endl;
      std::cout << "Iy = " << Iy << std::endl;
      std::cout << "Ixy = " << Ixy << std::endl;
      std::cout << "Ix3 = " << Ix3 << std::endl;
      std::cout << "Ix2y = " << Ix2y << std::endl;
      std::cout << "Ixy2 = " << Ixy2 << std::endl;
      std::cout << "Iy3 = " << Iy3 << std::endl;
      std::cout << "Ix2y2 = " << Ix2y2 << std::endl;

    }

    std::cout<<" If we use direct formula it will calculate integral of x^m y^n in D* . But we want integral of x^m y^n in D which is equivalent to integral of (1-x)^m y^n in D* . Quadratre rule give us correct result because of the change of basis in polybasisfunction.\n";

    Type direct_area_00 = find_trig_area_2intersection_formula<Type>(0, 0, 0, 0, 1, table_number,  q1,  q2, q3);
    Type direct_area_10 = find_trig_area_2intersection_formula<Type>(1, 0, 0, 0, 1, table_number,  q1,  q2, q3);
    Type direct_area_01 = find_trig_area_2intersection_formula<Type>(0, 1, 0, 0, 1, table_number,  q1,  q2, q3);
    Type direct_area_11 = find_trig_area_2intersection_formula<Type>(1, 1, 0, 0, 1, table_number,  q1,  q2, q3);
    Type direct_area_30 = find_trig_area_2intersection_formula<Type>(3, 0, 0, 0, 1, table_number,  q1,  q2, q3);
    Type direct_area_21 = find_trig_area_2intersection_formula<Type>(2, 1, 0, 0, 1, table_number,  q1,  q2, q3);
    Type direct_area_12 = find_trig_area_2intersection_formula<Type>(1, 2, 0, 0, 1, table_number,  q1,  q2, q3);
    Type direct_area_03 = find_trig_area_2intersection_formula<Type>(0, 3, 0, 0, 1, table_number,  q1,  q2, q3);
    Type direct_area_22 = find_trig_area_2intersection_formula<Type>(2, 2, 0, 0, 1, table_number,  q1,  q2, q3);
    cout << " area = "    << direct_area_00 << endl;
    cout << " area x = "  << direct_area_10 << endl;
    cout << " area y = "  << direct_area_01 << endl;
    cout << " area xy = " << direct_area_11 << endl;
    cout << " area x3 = " << direct_area_30 << endl;
    cout << " area x2y = "<< direct_area_21 << endl;
    cout << " area xy2 = "<< direct_area_12 << endl;
    cout << " area y3 = " << direct_area_03 << endl;
    cout << " area x2y2 ="<< direct_area_22 << endl;

//     Print the octree structure
//     std::cout << "Octree Structure:\n";
//     printOctreeStructure(root);

// printOctreeStructure(loadedRoots[0]);



  t = clock() - t;
  std::cout << "Time taken " << (Type)(t) / CLOCKS_PER_SEC << std::endl;

  return 1;
}



