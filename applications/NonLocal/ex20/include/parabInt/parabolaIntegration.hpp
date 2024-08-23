#ifndef PARABOLAINTEGRATION_HPP_INCLUDED
#define PARABOLAINTEGRATION_HPP_INCLUDED


#include <typeinfo>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <fstream>

using namespace std;

using boost::math::factorial;

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

template <class Type>
void swap(Type& x, Type& y) {
  Type temp = x;
  x = y;
  y = temp;
}

template <class Type>
struct Point {
  Type x;
  Type y;
};

template <class Type>
struct Circle {
  Point<Type> center;
  Type radius;
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

  if(fabs(k) < 1.e-16) k = 0 ;
  if(fabs(b) < 1.e-16) b = 0 ;
  if(fabs(d) < 1.e-16) d = 0 ;

  return {k, b, d};
}

template <class Type>
void CheckIntersection(int &intersect_number, unsigned int &table_number, std::vector <Type> &intersection, std::vector <Type> &interp_point, const Parabola <Type> &parabola) {

  table_number = 9999;
  intersect_number = 0;
  intersection.resize(0);
  interp_point.resize(0);
  Type k = parabola.k;
  Type b = parabola.b;
  Type d = parabola.d;
  Type c = 1;
  int left = 0, top = 0, right = 0, bottom = 0 ;
//   cout<< " parabola I get from solving system of linear equation :  " << parabola.k <<"x^2 + "<< parabola.b <<"x + "<< parabola.d << "+" <<  c << " y=0"  <<endl;

  if(-d >= 0 && -d <= 1) { //LEFT
    intersection.resize(intersection.size() + 2);
    intersection[intersection.size() - 2] = 0;
    intersection[intersection.size() - 1] = -d;
    left = 1 ;
    intersect_number += 1;
    interp_point.resize(interp_point.size() + 1);
    interp_point[interp_point.size() - 1] = -d;
//         cout << " left = " << left ;
  }
  // LEFT-TOP solve kx^2+bx+d-1 =0  ; Table 0
  if(k == 0) {
    Type x = (-1 - d) / b ;
    if(x <= 1 && x >= 0) {
      intersection.resize(intersection.size() + 2);
      intersection[intersection.size() - 2] = x;
      intersection[intersection.size() - 1] = 1;
      interp_point.resize(interp_point.size() + 1);
      interp_point[interp_point.size() - 1] = x;
      top = 1;
      intersect_number += 1;
      if(left == 1) table_number = 0 ;
//           cout << " top = " << top ;
    }
  }
  else {
    Type delta = b * b - 4 * k * (d + 1);
//         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
    if(delta >= 0) {
      Type sqrtdelta = sqrt(delta);
      int sign = (k > 0) ? 1 : -1;

      for(unsigned i = 0; i < 2; i++) {
        Type x = (- b - sign * sqrtdelta) / (2 * k);
        //             cout<< "Top x = "<< x<< endl;
        if(x > 1) break;
        else if(x >= 0) {
          intersection.resize(intersection.size() + 2);
          intersection[intersection.size() - 2] = x;
          intersection[intersection.size() - 1] = 1;
          interp_point.resize(interp_point.size() + 1);
          interp_point[interp_point.size() - 1] = x;
          intersect_number += 1;
          if(top == 1) {
            table_number = 3 ;
          }
          top += 1;
          if(left == 1) {
            table_number = 0 ;
          }
          //               cout << " top = " << top ;
        }
        sign *= -1;
      }
    }
  }
  Type y_1 = -(k + b + d); //LEFT-RIGHT x=1 ; Table 1
  if(y_1 >= 0 && y_1 <= 1) { //TODO check sign when normal change
    intersection.resize(intersection.size() + 2);
    intersection[intersection.size() - 2] = 1;
    intersection[intersection.size() - 1] = y_1;
    interp_point.resize(interp_point.size() + 1);
    interp_point[interp_point.size() - 1] = y_1;
    intersect_number += 1;
    if(left == 1) {
      table_number = 1 ;
    }
    if(top >= 1) {
      table_number = 4 ;
    }
    right = 1 ;
//           cout << " right = " << right ;
  }

  //LEFT-BOTTOM  solve kx^2+bx+d =0 ; Table 2
  if(k == 0) {
    Type x =  -d / b ;
    if(x <= 1 && x >= 0) {
      intersection.resize(intersection.size() + 2);
      intersection[intersection.size() - 2] = x;
      intersection[intersection.size() - 1] = 0;
      interp_point.resize(interp_point.size() + 1);
      interp_point[interp_point.size() - 1] = x;
      intersect_number += 1;
      if(left == 1) {
        table_number = 2 ;
      }
      if(right == 1) {
        table_number = 6 ;
      }
      if(top >= 1) {
        table_number = 5 ;
      }
      bottom = 1;
//             cout << " bottom = " << bottom ;
    }
  }

  else {
    Type delta = b * b - 4 * k * d;
    if(delta >= 0) {
      Type sqrtdelta = sqrt(delta);
      int sign = (k > 0) ? 1 : -1;

      for(unsigned i = 0; i < 2; i++) {
        Type x = (- b - sign * sqrtdelta) / (2 * k);
//               cout << " bottom root = " << x ;
        if(x > 1) break;
        else if(x >= 0) {
          intersection.resize(intersection.size() + 2);
          intersection[intersection.size() - 2] = x;
          intersection[intersection.size() - 1] = 0;
          interp_point.resize(interp_point.size() + 1);
          interp_point[interp_point.size() - 1] = x;
          if(bottom >= 1) {
            table_number = 7 ;
          }
          if(left == 1) {
            table_number = 2 ;
          }
          if(right == 1) {
            table_number = 6 ;
          }
          if(top == 1) {
            table_number = 5 ; // TODO check the table
          }
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
          Type sum1 = pow(pol1[2], pMax - r) / factorial<Type>(pMax - r) ;
          Type sum = sum1 * (pow(I3[w].second, 2 * r + m + 1) - pow(I3[w].first, 2 * r + m + 1)) / (2 * r + m + 1);
          for(int p = 1; p <= pMax - r; p++) {
            sum1 *= pol1[1] * (pMax - r - p + 1) / (pol1[2] * p)   ;
            sum += sum1 * (pow(I3[w].second, 2 * r + m + p + 1) - pow(I3[w].first, 2 * r + p + m + 1)) / (2 * r + m + p + 1) ;
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
void get_p1_p2_p3(const int &table, const std::vector<double> &corner, PointT <Type> &p1, PointT <Type> &p2, PointT <Type> &p3) {
  double epsilon = 0.000000000000001;
  Type i1_pm_eps(-1), i2_pm_eps(-1);

  // std::cout << "Corner " << i << ": (" << corner[0] << ", " << corner[1] << ", " << corner[2] << ") - Print Something\n";

  switch(table) {
  case 0:
//             i1_pm_eps = static_cast<Type>(corner[0]);
//             i2_pm_eps = static_cast<Type>(corner[1]);
//             i2_pm_eps = (corner[1] < epsilon)? static_cast<Type>(epsilon) : static_cast<Type>(corner[1]);

    i1_pm_eps = static_cast<Type>(corner[0] - epsilon);
    i2_pm_eps = static_cast<Type>(corner[1] + epsilon);
//                    if (i1 == partition ) i1_pm_eps = static_cast<Type>(i1*del_x - epsilon);     //it keeps my i2 in (0,1)
    p1 = {static_cast<Type>(0), i1_pm_eps};
    p2 = {i2_pm_eps, static_cast<Type>(1)};
    break;
  case 1:
    //Do we really need epsilon on this table?
    i1_pm_eps = static_cast<Type>(corner[0] + epsilon);
    i2_pm_eps = static_cast<Type>(corner[1] - epsilon);
    p1 = {static_cast<Type>(0), i1_pm_eps};
    p2 = {static_cast<Type>(1), i2_pm_eps};
    break;
  case 2:
    i1_pm_eps = static_cast<Type>(corner[0] + epsilon);
    i2_pm_eps = static_cast<Type>(corner[1] + epsilon);
    p1 = {static_cast<Type>(0), i1_pm_eps};
    p2 = {i2_pm_eps, static_cast<Type>(0)};
    break;
  case 3:
    i1_pm_eps = static_cast<Type>(corner[0] - epsilon);
    i2_pm_eps = static_cast<Type>(corner[1] + epsilon);
    p1 = {i1_pm_eps, static_cast<Type>(1)};
    p2 = {i2_pm_eps, static_cast<Type>(1)};
    break;
  case 4:
    i1_pm_eps = static_cast<Type>(corner[0] - epsilon);
    i2_pm_eps = static_cast<Type>(corner[1] - epsilon);
    p1 = {i1_pm_eps, static_cast<Type>(1)};
    p2 = {static_cast<Type>(1), i2_pm_eps};
    break;
  case 5:
    i1_pm_eps = static_cast<Type>(corner[0] + epsilon);
    i2_pm_eps = static_cast<Type>(corner[1] - epsilon);
    p1 = {i1_pm_eps, static_cast<Type>(1)};
    p2 = {i2_pm_eps, static_cast<Type>(0)};
    break;
  case 6:
    i1_pm_eps = static_cast<Type>(corner[0] + epsilon);
    i2_pm_eps = static_cast<Type>(corner[1] - epsilon);
    p1 = {static_cast<Type>(1), i1_pm_eps};
    p2 = {i2_pm_eps, static_cast<Type>(0)};
    break;
  case 7:
    i1_pm_eps = static_cast<Type>(corner[0] - epsilon);
    i2_pm_eps = static_cast<Type>(corner[1] + epsilon);
    p1 = {i1_pm_eps, static_cast<Type>(0)};
    p2 = {i2_pm_eps, static_cast<Type>(0)};
    break;
  }

  p3 = {(p1.x + p2.x) * 0.5, static_cast<Type>(corner[2])};
  if(fabs(corner[2] - 0) < epsilon) p3.y = static_cast<Type>(epsilon);
  if(fabs(corner[2] - 1) < epsilon) p3.y = static_cast<Type>(1 - epsilon);


}

template <class Type>
std::vector<Point<Type>> findGridIntersectionPoints(Circle<Type> circle, Point<Type> lineStart, Point<Type> lineEnd, std::vector<int> &ltrbNumber, unsigned int &table_number) {


  table_number = 9999;
  Type epsilon = 0.000000000000001; // put it in all the inequality.

  std::vector<Point<Type>> grid_intersections(0);
  // Circle equation: (x - h)^2 + (y - k)^2 = r^2
  Type h = circle.center.x;
  Type k = circle.center.y;
  Type r = circle.radius;


  Type diff_of_sqr, sqrt_diff_of_sqrt, value_p, value_m ;


  //x = i (Left)and circle
  diff_of_sqr = r * r - (lineStart.x - h) * (lineStart.x - h);
//     cout << " diff square = " << diff_of_sqr << endl;
  if(diff_of_sqr > epsilon) {
    sqrt_diff_of_sqrt = sqrt(diff_of_sqr);
    value_p = k + sqrt_diff_of_sqrt;
    value_m = k - sqrt_diff_of_sqrt;
//        cout << value_m << " -- " << lineStart.y << " " << lineEnd.y ;
    if(value_p + epsilon >= lineStart.y && value_p <= lineEnd.y + epsilon) {
      grid_intersections.resize(grid_intersections.size() + 1);
      grid_intersections[grid_intersections.size() - 1] = {lineStart.x, value_p};
      ltrbNumber[0] += 1;// it means the circle intersected the left.
    }
    if(value_m + epsilon >= lineStart.y && value_m <= lineEnd.y + epsilon) {
//         cout << " !# ";
      grid_intersections.resize(grid_intersections.size() + 1);
      grid_intersections[grid_intersections.size() - 1] = {lineStart.x, value_m};
      ltrbNumber[0] += 1;
    }
  }
  else if(fabs(diff_of_sqr - 0) <= epsilon) {
    if(k + epsilon >= lineStart.y && k <= lineEnd.y + epsilon) {
      grid_intersections.resize(grid_intersections.size() + 1);
      grid_intersections[grid_intersections.size() - 1] = {lineStart.x, k};
      ltrbNumber[0] += 1;  // it means the circle intersected the left.
    }
  }

  //y = j+1 (Top) and circle  corner points won't register here because of strickly less than.
  diff_of_sqr = r * r - (lineEnd.y - k) * (lineEnd.y - k);
//     cout << " diff square = " << diff_of_sqr << endl;
  if(diff_of_sqr > epsilon) {
    sqrt_diff_of_sqrt = sqrt(diff_of_sqr);
    value_p = h + sqrt_diff_of_sqrt;
    value_m = h - sqrt_diff_of_sqrt;
    if(value_p > lineStart.x + epsilon && value_p + epsilon < lineEnd.x) {
      grid_intersections.resize(grid_intersections.size() + 1);
      grid_intersections[grid_intersections.size() - 1] = {value_p, lineEnd.y};
      ltrbNumber[1] += 1;
    }
    if(value_m > lineStart.x + epsilon && value_m + epsilon < lineEnd.x) {
      grid_intersections.resize(grid_intersections.size() + 1);
      grid_intersections[grid_intersections.size() - 1] = {value_m, lineEnd.y};
      ltrbNumber[1] += 1;
//         cout << " use negative value in top. value m = " << value_m << "diff "<< value_m - lineEnd.x <<endl;
    }
  }
  else if(fabs(diff_of_sqr - 0) <= epsilon) {
    if(h > lineStart.x + epsilon && h + epsilon < lineEnd.x) {
      grid_intersections.resize(grid_intersections.size() + 1);
      grid_intersections[grid_intersections.size() - 1] = {h, lineEnd.y};
      ltrbNumber[1] += 1;  // it means the circle intersected the left.
    }
  }

  //x = i+1 (Right) and circle
  diff_of_sqr = r * r - (lineEnd.x - h) * (lineEnd.x - h);
//     cout << " diff square = " << diff_of_sqr << endl;
  if(diff_of_sqr > epsilon) {
    sqrt_diff_of_sqrt = sqrt(diff_of_sqr);
    value_p = k + sqrt_diff_of_sqrt;
    value_m = k - sqrt_diff_of_sqrt;
//        cout << value_m << " -- " <<value_p<< " " <<  lineStart.y << " " << lineEnd.y ;
    if(value_p + epsilon >= lineStart.y && value_p <= lineEnd.y + epsilon) {
      grid_intersections.resize(grid_intersections.size() + 1);
      grid_intersections[grid_intersections.size() - 1] = {lineEnd.x, value_p};
      ltrbNumber[2] += 1;
    }
    if(value_m + epsilon >= lineStart.y && value_m <= lineEnd.y + epsilon) {
      grid_intersections.resize(grid_intersections.size() + 1);
      grid_intersections[grid_intersections.size() - 1] = {lineEnd.x, value_m};
      ltrbNumber[2] += 1;
    }
  }
  else if(fabs(diff_of_sqr - 0) <= epsilon) {
    if(k + epsilon >= lineStart.y && k <= lineEnd.y + epsilon) {
      grid_intersections.resize(grid_intersections.size() + 1);
      grid_intersections[grid_intersections.size() - 1] = {lineEnd.x, k};
      ltrbNumber[2] += 1;  // it means the circle intersected the left.
    }
  }

  //y = j (bottom) and circle
  diff_of_sqr = r * r - (lineStart.y - k) * (lineStart.y - k);
//     cout << " diff square = " << diff_of_sqr << endl;
  if(diff_of_sqr > epsilon) {
    sqrt_diff_of_sqrt = sqrt(diff_of_sqr);
    value_p = h + sqrt_diff_of_sqrt;
    value_m = h - sqrt_diff_of_sqrt;
    if(value_p > lineStart.x + epsilon && value_p + epsilon < lineEnd.x) {
      grid_intersections.resize(grid_intersections.size() + 1);
      grid_intersections[grid_intersections.size() - 1] = {value_p, lineStart.y};
      ltrbNumber[3] += 1;
    }
    if(value_m > lineStart.x + epsilon && value_m + epsilon < lineEnd.x) {
      grid_intersections.resize(grid_intersections.size() + 1);
      grid_intersections[grid_intersections.size() - 1] = {value_m, lineStart.y};
      ltrbNumber[3] += 1;
    }
  }
  else if(fabs(diff_of_sqr - 0) <= epsilon) {
    if(h > lineStart.x + epsilon && h + epsilon < lineEnd.x) {
      grid_intersections.resize(grid_intersections.size() + 1);
      grid_intersections[grid_intersections.size() - 1] = {h, lineStart.y};
      ltrbNumber[3] += 1;  // it means the circle intersected the left.
    }
  }



  //TODO check if we have any duplicate points (corner points)
  /*
      for(int ia = 0 ; ia < grid_intersections.size() ; ia++){
        for(int ib = ia+1 ; ib < grid_intersections.size() ; ib++){
          if(grid_intersections[ia].x - grid_intersections[ib].x < 0.0000000001){
            if(grid_intersections[ia].y - grid_intersections[ib].y<0.0000000001)
            grid_intersections.erase(grid_intersections.begin()+ib);
          }
        }
      }*/

  //OR we can only check the corner points and take one value in determinant values in sqrt.


  // TODO if we have horizontal parabola we switch x and y .


  if(grid_intersections.size() == 1) { // it means it touches the corner never enters.
    grid_intersections.resize(0);
    cout << " one intersection deleting. \n " ;
  }

  else if(grid_intersections.size() == 3) { // TODO we have to get rid of one corner.
    cout << "######### get rid of a corner ########## -----------------------------------------------------------------------------------------------" << endl;
  }

  else if(grid_intersections.size() == 2) {

//       cout << "it intersect in 2 points.\n ";
    grid_intersections.resize(3);
    grid_intersections[2] = { (grid_intersections[0].x + grid_intersections[1].x) * 0.5, -1.0};
    diff_of_sqr = r * r - ((grid_intersections[2].x - h) * (grid_intersections[2].x - h));
    if(diff_of_sqr >= 0) {  //it is different for table 5
      sqrt_diff_of_sqrt = sqrt(diff_of_sqr);
      value_p = k + sqrt_diff_of_sqrt;
      value_m = k - sqrt_diff_of_sqrt;
      if(value_p >= lineStart.y && value_p <= lineEnd.y) {
        grid_intersections[2].y = value_p;
      }
      else if(value_m >= lineStart.y && value_m <= lineEnd.y) {
        grid_intersections[2].y = value_m;
      }
//         else {
//           grid_intersections.resize(0);
//           ltrbNumber={0,0,0,0};
//           // meaning we have two intersection but mid_point is outside. it is either 4 intersection or never inters}
//         }
    }

//           cout << "LTRB = " << ltrbNumber[0] << ltrbNumber[1]<<ltrbNumber[2]<<ltrbNumber[3] <<endl;

    if(ltrbNumber[0] == 1) {
      if(ltrbNumber[1] == 1) {
        table_number = 0;
      }
      else if(ltrbNumber[2] == 1) {
        table_number = 1;
      }
      else if(ltrbNumber[3] == 1) {
        table_number = 2;
      }
    }

    else if(ltrbNumber[1] == 2) table_number = 3;
    else if(ltrbNumber[1] == 1) {
      if(ltrbNumber[2] == 1) table_number = 4;
      else if(ltrbNumber[3] == 1) table_number = 5;
    }
    else if(ltrbNumber[2] == 1 && ltrbNumber[3] == 1) table_number = 6;
    else if(ltrbNumber[3] == 2) table_number = 7;

//       swapping x and y when we have Top - Bottom (table-5) into table-1
    if(grid_intersections.size() == 3) {
      cout << " (" << grid_intersections[0].x << "," << grid_intersections[0].y << ") " << " (" << grid_intersections[1].x << "," << grid_intersections[1].y << ") " << " (" << grid_intersections[2].x << "," << grid_intersections[2].y << ") " << endl;
    }

    if(table_number == 5) {
      grid_intersections[2] = { -1.0, (grid_intersections[0].y + grid_intersections[1].y) * 0.5 };
      diff_of_sqr = r * r - ((grid_intersections[2].y - k) * (grid_intersections[2].y - k));
// //         cout << " diff square = " << diff_of_sqr << endl;
      if(diff_of_sqr >= 0) {
        sqrt_diff_of_sqrt = sqrt(diff_of_sqr);
        value_p = h + sqrt_diff_of_sqrt;
        value_m = h - sqrt_diff_of_sqrt;
// //           cout<< "value_m = "<<value_m << " value_p = "<< value_p <<endl;
        if(value_p >= lineStart.x && value_p <= lineEnd.x) {
          grid_intersections[2].x = value_p;
        }
        else if(value_m >= lineStart.x && value_m <= lineEnd.x) {
          grid_intersections[2].x = value_m;
        }
        else {
          grid_intersections.resize(0);
          ltrbNumber = {0, 0, 0, 0};
          // meaning we have two intersection but mid_point is outside. it is either 4 intersection or never inters}
        }
      }
    }

  }

//     cout<< " grid intersection size "<<grid_intersections.size() <<endl;
//     for (int gi = 0 ; gi <grid_intersections.size(); gi++){
//       cout<< " ("<<grid_intersections[gi].x<<"," << grid_intersections[gi].y<<") \n" ;
//     }
//     cout << "table = " << table_number <<endl;



//     cout << "table number aa : " << table_number << " " ;

  //introduce middle point if we have two intersections:
//     if (grid_intersections.size() == 2){
//       grid_intersections.resize(3);
//       grid_intersections[2] = { (grid_intersections[0].x + grid_intersections[1].x)/2.0 , -1.0};
//       diff_of_sqr = r*r - ((grid_intersections[2].x - h)*(grid_intersections[2].x - h));
//       if (diff_of_sqr >=0 ){
//         sqrt_diff_of_sqrt = sqrt(diff_of_sqr);
//         value_p = k + sqrt_diff_of_sqrt;
//         value_m = k - sqrt_diff_of_sqrt;
//         if (value_p +epsilon >=lineStart.y && value_p <= lineEnd.y+epsilon){
//           grid_intersections[2].y = value_m;
//         }
//         else if (value_m >=lineStart.y && value_m <= lineEnd.y ){
//           grid_intersections[2].y = value_m;
//         }
//         else grid_intersections.resize(0);  // meaning we have two intersection but mid_point is outside. it is either 4 intersection or never inters}
//       }
//     }

// TODO checkout the Pushback function. It looks nice and easier.
//     if (discriminant >= 0) {
//         // Calculate the intersection points
//         double x1 = (-b + std::sqrt(discriminant)) / (2 * a);
//         double y1 = m * x1 + c;
//         intersections.push_back({x1, y1});
//
//         double x2 = (-b - std::sqrt(discriminant)) / (2 * a);
//         double y2 = m * x2 + c;
//         intersections.push_back({x2, y2});
//     }

  return grid_intersections;
}

template <class Type>
std::vector<Point<Type>> mapIntersectionPoints(std::vector<Point<Type>>grid_intersections, Point<Type> gridStart, Point<Type> gridEnd, unsigned int &table_number, std::vector< Type > &interp_point, bool &isswap) {

  std::vector<Point<Type>> mapped_intersections(0);
  interp_point.resize(0);
  for(int aa = 0 ; aa < grid_intersections.size(); aa++) {
    mapped_intersections.resize(mapped_intersections.size() + 1);
    mapped_intersections[mapped_intersections.size() - 1] = {((grid_intersections[aa].x - gridStart.x) / (gridEnd.x - gridStart.x)), ((grid_intersections[aa].y - gridStart.y) / (gridEnd.y - gridStart.y))};
    if(mapped_intersections[aa].x < 0) mapped_intersections[aa].x = 0;
    else if(mapped_intersections[aa].x > 1) mapped_intersections[aa].x = 1;
    if(mapped_intersections[aa].y < 0) mapped_intersections[aa].y = 0;
    else if(mapped_intersections[aa].y > 1) mapped_intersections[aa].y = 1;   //TODO check if it should be .y
//       mapped_intersections.push_back({((intersections[aa].x - gridStart.x)/(gridEnd.x -gridStart.x)) , ((intersections[aa].y - gridStart.y)/           //       (gridEnd.y -gridStart.y))});

  }
  for(int mk = 0; mk < mapped_intersections.size(); mk++) {
    cout << "before swapping : (" << mapped_intersections[mk].x << " " << mapped_intersections[mk].y << ") " << endl;
  }


  if(table_number == 5) {
    cout << " swapping values :::::::::" << endl;
    isswap = 1;
    Type swapvalue;
    swapvalue = mapped_intersections[0].x ;
    mapped_intersections[0].x = mapped_intersections[0].y;
    mapped_intersections[0].y = swapvalue;

    swapvalue = mapped_intersections[1].x ;
    mapped_intersections[1].x = mapped_intersections[1].y;
    mapped_intersections[1].y = swapvalue;

    swapvalue = mapped_intersections[2].x ;
    mapped_intersections[2].x = mapped_intersections[2].y;
    mapped_intersections[2].y = swapvalue;

    //TODO fro direct formula we need to exchange p1 and p2

    swapvalue = mapped_intersections[0].x ;
    mapped_intersections[0].x = mapped_intersections[1].x;
    mapped_intersections[1].x = swapvalue;

    swapvalue = mapped_intersections[0].y ;
    mapped_intersections[0].y = mapped_intersections[1].y;
    mapped_intersections[1].y = swapvalue;



    table_number = 1;
    for(int mk = 0; mk < mapped_intersections.size(); mk++) {
      cout << "after swapping : (" << mapped_intersections[mk].x << " " << mapped_intersections[mk].y << ") " << endl;
    }
  }

  if(mapped_intersections.size() == 3) { // Finding interpolation points
    switch(table_number) {
    case 0:
      interp_point.resize(3);
      interp_point = { static_cast<Type>(mapped_intersections[0].y), static_cast<Type>(mapped_intersections[1].x), static_cast<Type>(mapped_intersections[2].y)};
      break;
    case 1:
      interp_point.resize(3);
      interp_point = { static_cast<Type>(mapped_intersections[0].y), static_cast<Type>(mapped_intersections[1].y), static_cast<Type>(mapped_intersections[2].y)};
      break;
    case 2:
      interp_point.resize(3);
      interp_point = { static_cast<Type>(mapped_intersections[0].y), static_cast<Type>(mapped_intersections[1].x), static_cast<Type>(mapped_intersections[2].y)};
      break;
    case 3:
      interp_point.resize(3);
      interp_point = { static_cast<Type>(mapped_intersections[0].x), static_cast<Type>(mapped_intersections[1].x), static_cast<Type>(mapped_intersections[2].y)};
      break;
    case 4:
      interp_point.resize(3);
      interp_point = { static_cast<Type>(mapped_intersections[0].x), static_cast<Type>(mapped_intersections[1].y), static_cast<Type>(mapped_intersections[2].y)};
      break;
    case 5:
      interp_point.resize(3);
      interp_point = { static_cast<Type>(mapped_intersections[0].x), static_cast<Type>(mapped_intersections[1].x), static_cast<Type>(mapped_intersections[2].y)};
      break;
    case 6:
      interp_point.resize(3);
      interp_point = { static_cast<Type>(mapped_intersections[0].y), static_cast<Type>(mapped_intersections[1].x), static_cast<Type>(mapped_intersections[2].y)};
      break;
    case 7:
      interp_point.resize(3);
      interp_point = { static_cast<Type>(mapped_intersections[0].x), static_cast<Type>(mapped_intersections[1].x), static_cast<Type>(mapped_intersections[2].y)};
      break;
    }
  }


  return mapped_intersections;
}

template <class Type>
Type find_area(const unsigned &m, const unsigned &n, const int &s, const Type &a, Type c, const int &table,  PointT <Type> &p1,  PointT <Type> &p2, const PointT <Type> &p3) {
  Type area(0);
  Type A1(0), A2(0), A3(0);
  std::vector< std::pair <Type, Type> > I1, I2, I3 ;
  int intersect_number;
  std::vector <Type> pol1(3, 0);
  std::vector <Type> pol2(3, 0);
  std::vector <Type> intersection;
  std::vector <Type> interp_point;    //never used in this function. it was used in interpolation;
  unsigned int table_number = table;
  Parabola <Type> parabola;
  Type det = p1.x * p1.x * (p2.x - p3.x) - p1.x * (p2.x * p2.x - p3.x * p3.x) + p2.x * p3.x * (p2.x - p3.x) ; // only sort out the points parallel to y axis
  parabola = get_parabola_equation(p1, p2, p3);
  CheckIntersection <Type> (intersect_number, table_number, intersection, interp_point, parabola);


  Type k = parabola.k;
  Type b = parabola.b;
  Type d = parabola.d;
  if(table == 0) {
    if(k > 0) {
      Type delta = b * b - 4 * k * (d + 1);
//         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
      if(delta >= 0) {
        Type sqrtdelta = sqrt(delta);
        int sign = (k > 0) ? 1 : -1;
        for(unsigned i = 0; i < 2; i++) {
          Type x = (- b - sign * sqrtdelta) / (2 * k);
          //             cout<< "Top x = "<< x<< endl;
          if(x > 0 && x < p2.x) {
            p2.x = x;
          }
          sign *= -1;
        }
      }
    }
  }

  if(table == 2) {
    if(k < 0) {
      Type delta = b * b - 4 * k * d;
//         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
      if(delta >= 0) {
        Type sqrtdelta = sqrt(delta);
        int sign = (k > 0) ? 1 : -1;
        for(unsigned i = 0; i < 2; i++) {
          Type x = (- b - sign * sqrtdelta) / (2 * k);
          //             cout<< "Top x = "<< x<< endl;
          if(x > 0 && x < p2.x) {
            p2.x = x;
          }
          sign *= -1;
        }
      }
    }
  }

  if(table == 4) {
    if(k > 0) {
      Type delta = b * b - 4 * k * (d + 1);
//         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
      if(delta >= 0) {
        Type sqrtdelta = sqrt(delta);
        int sign = (k > 0) ? 1 : -1;
        for(unsigned i = 0; i < 2; i++) {
          Type x = (- b - sign * sqrtdelta) / (2 * k);
          //             cout<< "Top x = "<< x<< endl;
          if(x < 1 && x > p1.x) {
            p1.x = x;
          }
          sign *= -1;
        }
      }
    }
  }

  if(table == 6) {
    if(k < 0) {
      Type delta = b * b - 4 * k * d;
//         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
      if(delta >= 0) {
        Type sqrtdelta = sqrt(delta);
        int sign = (k > 0) ? 1 : -1;
        for(unsigned i = 0; i < 2; i++) {
          Type x = (- b - sign * sqrtdelta) / (2 * k);
          //             cout<< "Top x = "<< x<< endl;
          if(x < 1 && x > p2.x) {
            p2.x = x;
          }
          sign *= -1;
        }
      }
    }
  }
//     cout<< "parabola = " << parabola.k <<"x^2 +"<<parabola.b<<"x +" <<parabola.d << "+y=0 "<<endl;  //     cout<< "intersection = " << intersect_number << " table= " << table_number <<endl;
  k = -parabola.k;
  b = -parabola.b;
  d = -parabola.d;
  if(table == 2) {
    area = (1 - d - k * p2.x * p2.x / 3 - b * p2.x / 2) * p2.x - (1 - d - k * p1.x * p1.x / 3 - b * p1.x / 2) * p1.x + 1 - p2.x ;
  }
  else if(table == 6) {
    area = (1 - d - k * p1.x * p1.x / 3 - b * p1.x / 2) * p1.x - (1 - d - k * p2.x * p2.x / 3 - b * p2.x / 2) * p2.x + p2.x ;
  }
  else {
    area = (1 - d - k * p2.x * p2.x / 3 - b * p2.x / 2) * p2.x - (1 - d - k * p1.x * p1.x / 3 - b * p1.x / 2) * p1.x;
  }

  return area ;
}

template <class Type>
Type find_area_2intersection_formula(const unsigned &m, const unsigned &n, const int &s, const Type &a, Type c, const int &table,  PointT <Type> &p1,  PointT <Type> &p2, const PointT <Type> &p3) {
  Type area(0);
  Type A1(0), A2(0), A3(0);
  std::vector< std::pair <Type, Type> > I1, I2, I3 ;
  int intersect_number;
  std::vector <Type> pol1(3, 0);
  std::vector <Type> pol2(3, 0);
  std::vector <Type> intersection;
  std::vector <Type> interp_point;    //never used in this function. it was used in interpolation;
  unsigned int table_number = table;
  Parabola <Type> parabola;
  Type det = p1.x * p1.x * (p2.x - p3.x) - p1.x * (p2.x * p2.x - p3.x * p3.x) + p2.x * p3.x * (p2.x - p3.x); // only sort out the points parallel to y axis
  parabola = get_parabola_equation(p1, p2, p3);
  CheckIntersection <Type> (intersect_number, table_number, intersection, interp_point, parabola);

  Type k = parabola.k;
  Type b = parabola.b;
  Type d = parabola.d;

  bool do_line = 0;

  if(table == 0) {
    if(k > 0) {
      do_line = 1;
      Type delta = b * b - 4 * k * (d + 1);
//         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
      if(delta >= 0) {
        Type sqrtdelta = sqrt(delta);
        int sign = (k > 0) ? 1 : -1;
        for(unsigned i = 0; i < 2; i++) {
          Type x = (- b - sign * sqrtdelta) / (2 * k);
          //             cout<< "Top x = "<< x<< endl;
          if(x > 0 && x < p2.x) {
            p2.x = x;
          }
          sign *= -1;
        }
      }
    }
  }

  if(table == 2) {
    if(k < 0) {
      do_line = 1;
      Type delta = b * b - 4 * k * d;
//         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
      if(delta >= 0) {
        Type sqrtdelta = sqrt(delta);
        int sign = (k > 0) ? 1 : -1;
        for(unsigned i = 0; i < 2; i++) {
          Type x = (- b - sign * sqrtdelta) / (2 * k);
          //             cout<< "Top x = "<< x<< endl;
          if(x > 0 && x < p2.x) {
            p2.x = x;
          }
          sign *= -1;
        }
      }
    }
  }

  if(table == 4) {
    if(k > 0) {
      do_line = 1;
      Type delta = b * b - 4 * k * (d + 1);
//         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
      if(delta >= 0) {
        Type sqrtdelta = sqrt(delta);
        int sign = (k > 0) ? 1 : -1;
        for(unsigned i = 0; i < 2; i++) {
          Type x = (- b - sign * sqrtdelta) / (2 * k);
          //             cout<< "Top x = "<< x<< endl;
          if(x < 1 && x > p1.x) {
            p1.x = x;
          }
          sign *= -1;
        }
      }
    }
  }

  if(table == 6) {
    if(k < 0) {
      do_line = 1;
      Type delta = b * b - 4 * k * d;
//         cout << " k = "<< k << " b = "<< b << " d ="<< d << " delta = " << delta <<endl;
      if(delta >= 0) {
        Type sqrtdelta = sqrt(delta);
        int sign = (k > 0) ? 1 : -1;
        for(unsigned i = 0; i < 2; i++) {
          Type x = (- b - sign * sqrtdelta) / (2 * k);
          //             cout<< "Top x = "<< x<< endl;
          if(x < 1 && x > p2.x) {
            p2.x = x;
          }
          sign *= -1;
        }
      }
    }
  }

  pol2[0] = parabola.k;
  pol2[1] = parabola.b;
  pol2[2] = parabola.d;

//     if (c == -1){
//       pol2[0] *= -1; pol2[1] *= -1; pol2[2] *= -1;
//     }

  if(do_line) {
    if(table == 0) {
      I1.resize(1, std::pair<Type, Type>(static_cast<Type>(0), static_cast<Type>(p2.x)));
      area = easy_integral_A3(m, n, s, a, c, pol2, I1) -  easy_integral_A2(m, n, s, a, c, pol2, I1);
    }

    else if(table == 2) {
      I1.resize(1, std::pair<Type, Type>(static_cast<Type>(0), static_cast<Type>(p2.x)));
      I3.resize(1, std::pair<Type, Type>(static_cast<Type>(p2.x), static_cast<Type>(1)));
      area = easy_integral_A3(m, n, s, a, c, pol2, I1) -  easy_integral_A2(m, n, s, a, c, pol2, I1) + easy_integral_A3(m, n, s, a, c, pol2, I3);
    }

    else if(table == 4) {
      I1.resize(1, std::pair<Type, Type>(static_cast<Type>(p1.x), static_cast<Type>(1)));
      area = easy_integral_A3(m, n, s, a, c, pol2, I1) -  easy_integral_A2(m, n, s, a, c, pol2, I1);
    }

    else if(table == 6) {
      I1.resize(1, std::pair<Type, Type>(static_cast<Type>(p2.x), static_cast<Type>(1)));
      I3.resize(1, std::pair<Type, Type>(static_cast<Type>(0), static_cast<Type>(p2.x)));
      area = easy_integral_A3(m, n, s, a, c, pol2, I1) -  easy_integral_A2(m, n, s, a, c, pol2, I1) + easy_integral_A3(m, n, s, a, c, pol2, I3);
    }
  }

  else {
    pol1[0] = parabola.k ;
    pol1[1] = a + parabola.b;
    pol1[2] = c + parabola.d;
//       if (c == -1) { pol1[0] *= -1; pol1[1] *= -1; pol1[2] *= -1; }

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
  }


  return area ;
}


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
  const Type &x = interp_point[0];
  const Type &y = interp_point[1];
  const Type &z = interp_point[2];

  const Type &x0 = interp_table[0][0];
  const Type &x1 = interp_table[7][0];
  const Type &y0 = interp_table[0][1];
  const Type &y1 = interp_table[7][1];
  const Type &z0 = interp_table[0][2];
  const Type &z1 = interp_table[7][2];

  Type x_d = (x - x0) / (x1 - x0);
  Type oneMx_d = 1. - x_d;
  Type y_d = (y - y0) / (y1 - y0);
  Type oneMy_d = 1. - y_d;
  Type z_d = (z - z0) / (z1 - z0);
  Type oneMz_d = 1. - z_d;

  for(unsigned i = 0; i < interp_table_values[0].size(); i++) {
    const Type & c_000 = interp_table_values[0][i];
    const Type & c_001 = interp_table_values[1][i];
    const Type & c_010 = interp_table_values[2][i];
    const Type & c_011 = interp_table_values[3][i];
    const Type & c_100 = interp_table_values[4][i];
    const Type & c_101 = interp_table_values[5][i];
    const Type & c_110 = interp_table_values[6][i];
    const Type & c_111 = interp_table_values[7][i];

    Type c_00 = c_000 * oneMx_d + c_100 * x_d ;
    Type c_01 = c_001 * oneMx_d + c_101 * x_d ;
    Type c_10 = c_010 * oneMx_d + c_110 * x_d ;
    Type c_11 = c_011 * oneMx_d + c_111 * x_d ;

    Type c_0 = c_00 * oneMy_d + c_10 * y_d ;
    Type c_1 = c_01 * oneMy_d + c_11 * y_d ;

    interp_point_values[i] = c_0 * oneMz_d + c_1 * z_d ;
  }
}


template <class Type>
double findGridArea(const Circle<Type> &circle, const Point<Type> &gridStart, const Point<Type> &gridEnd, std::vector<Point<Type>> grid_intersections, std:: vector <std:: vector <double>> interp_table, const std::vector< Type > &interp_point, const Type &grid_size) {
  double area = 0;
  if(grid_intersections.size() <= 1) {
    // checking wherher the grid is outside or inside the circle . Circle equation: (x - h)^2 + (y - k)^2 = r^2 Find (i+epsilon-h)^2 +(j+epsilon-k)^2 < r^2
    area = ((gridStart.x - circle.center.x) * (gridStart.x - circle.center.x) + (gridStart.y - circle.center.y) * (gridStart.y - circle.center.y) < (circle.radius * circle.radius)) ? 1 : 0 ;
  }
  else if(interp_point.size() == 3) { // here goes our ex_23

    cout << " using interpolation " << endl;
    std::vector< double > double_interp_point = {static_cast<double>(interp_point[0]), static_cast<double>(interp_point[1]), static_cast<double>(interp_point[2])};

    area = trilinier_interpolation<double>(interp_table, double_interp_point) ;


//       cout << "grid area before scaling = "<< area <<endl;
  }

  else {} // direct calculation for 4 intersections. use grid intersections to find the area of for intersection points.

//   cout << "grid size = "<< grid_size <<endl;
  area = area * static_cast<double>(grid_size) * static_cast<double>(grid_size) ;
//   cout << "grid area after scaling = "<< area <<endl;

  return area;


}

struct Point3D {
  double x, y, z;
  Point3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
};





template <class Type>
class OctreeNode {
  public:
    Point3D minBounds, maxBounds;
    bool isLeaf;
    std::vector<OctreeNode*> children;
    std::vector<std::vector<double>> corners;
    std::vector<std::vector<double>> cornerAreas;
    std::vector<std::vector<double>> cornerWeights;
    int table;
    unsigned depth;
    unsigned qM;
    int s = 0;
    Type a = 0;
    double relative_error = -1;
    double relative_error_opposite = -1;
    CutFemWeightParabola <double, Type> *_Pweights;

    OctreeNode(const Point3D& _minBounds, const Point3D& _maxBounds, const int& _table, const int& _depth, const unsigned& _qM, CutFemWeightParabola <double, Type> *Pweights)
      : minBounds(_minBounds), maxBounds(_maxBounds), isLeaf(true), table(_table), depth(_depth), qM(_qM), _Pweights(Pweights) {}

    // Destructor to properly clean up memory
    ~OctreeNode() {
      for(OctreeNode* child : children) {
        delete child;
      }
      children.clear();
    }

// Function to get the eight corners of the node
    void getCorners() {  //TODO initialize it once without using push_back
      corners.resize(8, std::vector<double>(3));
      corners[0] = {minBounds.x, minBounds.y, minBounds.z};
      corners[1] = {minBounds.x, minBounds.y, maxBounds.z};
      corners[2] = {minBounds.x, maxBounds.y, minBounds.z};
      corners[3] = {minBounds.x, maxBounds.y, maxBounds.z};
      corners[4] = {maxBounds.x, minBounds.y, minBounds.z};
      corners[5] = {maxBounds.x, minBounds.y, maxBounds.z};
      corners[6] = {maxBounds.x, maxBounds.y, minBounds.z};
      corners[7] = {maxBounds.x, maxBounds.y, maxBounds.z};

//         for (const auto& corner : corners) {
//             std::cout << "Corner: (" << corner[0] << ", " << corner[1] << ", " << corner[2] << ") - Print Something\n";
//         }
      Type area(0) ;
      Type c(1) ;
      cornerAreas.resize(8);
      cornerWeights.resize(8);
      PointT <Type> p1, p2, p3 ;
      //CutFemWeightParabola <double, Type> Pweight(QUAD, 2, "legendre");

      for(size_t i = 0; i < corners.size(); ++i) {
        const auto& corner = corners[i];
//             std::cout << "Corner " << i << ": (" << corner[0] << ", " << corner[1] << ", " << corner[2] << ") - Print Something\n";searchPoint
        get_p1_p2_p3(table, corner, p1, p2, p3);

        int count = 0;
        for(unsigned qq = 0; qq <= qM; qq++) {
          for(unsigned jj = 0; jj <= qq; jj++) {
            unsigned ii = qq - jj;
            area = find_area_2intersection_formula(ii, jj, s, a, c, table, p1, p2, p3);
            cornerAreas[i].push_back(static_cast<double>(area));
//                     cout<< " f values = " << i << " "<< j << " "<< _f[count] <<endl;
            count++;
          }
        }

        (*_Pweights)(s, a, c, table, p1, p2, p3, cornerWeights[i]);


//             for (size_t ii = 0; ii <= mPn; ++ii){
//               for (size_t jj = 0; ii+jj <= mPn; ++jj){
//                 area = find_area_2intersection_formula(jj,ii,s,a,c,table, p1,p2, p3);
//                 corners[i].push_back(static_cast<double>(area));
//               }
//             }

      }
    }

    void subdivideWithRelativeError(int maxDepth, double maxRelativeError, int currentDepth = 0) {
      if(currentDepth >= maxDepth || !isLeaf) {
        getCorners();
        return;
      }

      getCorners();

      double midX = (minBounds.x + maxBounds.x) / 2.0;
      double midY = (minBounds.y + maxBounds.y) / 2.0;
      double midZ = (minBounds.z + maxBounds.z) / 2.0;
      std::vector<double> relativeErrors;
      std::vector<double> relativeErrorsOpposite;
      std::vector<std::vector<double>>boundary_mid_points(19);

      boundary_mid_points[0] = {midX, minBounds.y, minBounds.z};
      boundary_mid_points[1] = {maxBounds.x, midY, minBounds.z};
      boundary_mid_points[2] = {midX, maxBounds.y, minBounds.z};
      boundary_mid_points[3] = {minBounds.x, midY, minBounds.z};
      boundary_mid_points[4] = {minBounds.x, minBounds.y, midZ};
      boundary_mid_points[5] = {maxBounds.x, minBounds.y, midZ};
      boundary_mid_points[6] = {maxBounds.x, maxBounds.y, midZ};
      boundary_mid_points[7] = {minBounds.x, maxBounds.y, midZ};
      boundary_mid_points[8] = {midX, minBounds.y, maxBounds.z};
      boundary_mid_points[9] = {maxBounds.x, midY, maxBounds.z};
      boundary_mid_points[10] = {midX, maxBounds.y, maxBounds.z};
      boundary_mid_points[11] = {minBounds.x, midY, maxBounds.z};
      boundary_mid_points[12] = {midX, midY, minBounds.z};
      boundary_mid_points[13] = {minBounds.x, midY, midZ};
      boundary_mid_points[14] = {midX, minBounds.y, midZ};
      boundary_mid_points[15] = {midX, midY, maxBounds.z};
      boundary_mid_points[16] = {maxBounds.x, midY, midZ};
      boundary_mid_points[17] = {midX, maxBounds.y, midZ};
      boundary_mid_points[18] = {midX, midY, midZ};

      for(int i = 0; i < boundary_mid_points.size(); ++i) {  //TODO all nineteen points
        // Generate a random 3D point within the bounds
//             double randomX = minBounds.x + static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * (maxBounds.x - minBounds.x);
//             double randomY = minBounds.y + static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * (maxBounds.y - minBounds.y);
//             double randomZ = minBounds.z + static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * (maxBounds.z - minBounds.z);

        //not random using middle boundary points
        double randomX = boundary_mid_points[i][0];
        double randomY = boundary_mid_points[i][1];
        double randomZ = boundary_mid_points[i][2];

        std::vector<double> interp_point = {randomX, randomY, randomZ};


        Type f_area(0);
        Type c = 1 ;
        PointT<Type> p1, p2, p3;
        get_p1_p2_p3(table, interp_point, p1, p2, p3);

        std::vector<std::vector<double>>interpolation_vector(8);
        size_t ii = 0;
        size_t jj = 0;
//             for (size_t ii = 0; ii <= mPn; ++ii){
//               for (size_t jj = 0; ii+jj <= mPn; ++jj){
//
//                 size_t index = 0;
//                 for (size_t mm = 0; mm < ii; ++mm){
//                   index += mPn + 1 - mm ;
//                 }
//
//                 index += jj+3 ;

        int count = 0;
        for(unsigned qq = 0; qq <= qM; qq++) {
          for(unsigned jj = 0; jj <= qq; jj++) {
            unsigned ii = qq - jj;
            for(size_t ic = 0; ic < corners.size(); ++ic) {
              interpolation_vector[ic] = {corners[ic][0], corners[ic][1], corners[ic][2], cornerAreas[ic][count]};
            }
            double interp_area = trilinier_interpolation(interpolation_vector, interp_point);
            f_area = find_area_2intersection_formula(jj, ii, s, a, c, table, p1, p2, p3);
            double formula_area = static_cast<double>(f_area);
            double r_error = fabs(formula_area - interp_area) / formula_area;
            double r_error_opposite = fabs(formula_area - interp_area) / (1.0 / (ii + 1) * (jj + 1) - formula_area);

            relativeErrors.push_back(r_error);
            relativeErrorsOpposite.push_back(r_error_opposite);

            count++;
          }
        }






//                 for (size_t ic = 0; ic < corners.size(); ++ic) {
//                   interpolation_vector[ic] = {corners[ic][0],corners[ic][1],corners[ic][2],cornerAreas[ic][0]};
//                 }
//                 double interp_area = trilinier_interpolation(interpolation_vector, interp_point);
//                 f_area = find_area_2intersection_formula(jj,ii,s,a,c,table, p1,p2, p3);
//                 double formula_area = static_cast<double>(f_area);
//                 double r_error = fabs(formula_area - interp_area) / formula_area;
//                 double r_error_opposite = fabs(formula_area - interp_area) /(1.0/(ii+1)*(jj+1) - formula_area);
//
//                 relativeErrors.push_back(r_error);
//                 relativeErrorsOpposite.push_back(r_error_opposite);
//
//
// //               }
// //             }
//
//
// //                 f_area = find_area_2intersection_formula(m, n, s, a, c, table, p1, p2, p3);
// //                 double formula_area = static_cast<double>(f_area);
//
//     //             double r_error = (m+1)*(n+1) * fabs(formula_area - interp_area) / formula_area;
//     //             double r_error_opposite = fabs(formula_area - interp_area) /(1.0/(m+1)*(n+1) - formula_area);
//
// //                 double r_error = fabs(formula_area - interp_area) / formula_area;
// //                 double r_error_opposite = fabs(formula_area - interp_area) /(1.0/(m+1)*(n+1) - formula_area);
// //
// //                 relativeErrors.push_back(r_error);
// //                 relativeErrorsOpposite.push_back(r_error_opposite);



      }
      relative_error = *std::max_element(relativeErrors.begin(), relativeErrors.end());
      relative_error_opposite = *std::max_element(relativeErrorsOpposite.begin(), relativeErrorsOpposite.end());

      if(depth <= 3 || relative_error > maxRelativeError || relative_error_opposite > maxRelativeError || minBounds.y < 0.003) {
        isLeaf = false;
        children.push_back(new OctreeNode(minBounds, {midX, midY, midZ}, table, depth + 1, qM, _Pweights));
        children.push_back(new OctreeNode({midX, minBounds.y, minBounds.z}, {maxBounds.x, midY, midZ}, table, depth + 1, qM, _Pweights));
        children.push_back(new OctreeNode({minBounds.x, midY, minBounds.z}, {midX, maxBounds.y, midZ}, table, depth + 1, qM, _Pweights));
        children.push_back(new OctreeNode({midX, midY, minBounds.z}, {maxBounds.x, maxBounds.y, midZ}, table, depth + 1, qM, _Pweights));
        children.push_back(new OctreeNode({minBounds.x, minBounds.y, midZ}, {midX, midY, maxBounds.z}, table, depth + 1, qM, _Pweights));
        children.push_back(new OctreeNode({midX, minBounds.y, midZ}, {maxBounds.x, midY, maxBounds.z}, table, depth + 1, qM, _Pweights));
        children.push_back(new OctreeNode({minBounds.x, midY, midZ}, {midX, maxBounds.y, maxBounds.z}, table, depth + 1, qM, _Pweights));
        children.push_back(new OctreeNode({midX, midY, midZ}, maxBounds, table, depth + 1, qM, _Pweights));

        // Recursively check for subdivision in each child node
        for(OctreeNode* child : children) {
          child->subdivideWithRelativeError(maxDepth, maxRelativeError, currentDepth + 1);
        }
      }
    }
    // Recursive search for the smallest sub-cube containing a given point

    OctreeNode* search(const Point3D& point) {
      if(isLeaf) {
        return this;
      }


      if(point.z <= children[0]->maxBounds.z) { // 0->3
        if(point.y <= children[0]->maxBounds.y) { // 0->1
          if(point.x <= children[0]->maxBounds.x)
            return children[0]->search(point);
          else
            return children[1]->search(point);
        }
        else { // 2->3
          if(point.x <= children[2]->maxBounds.x)
            return children[2]->search(point);
          else
            return children[3]->search(point);
        }
      }
      else { // 4->7
        if(point.y <= children[4]->maxBounds.y) { // 4->5
          if(point.x <= children[4]->maxBounds.x)
            return children[4]->search(point);
          else
            return children[5]->search(point);
        }
        else { // 6->7
          if(point.x <= children[6]->maxBounds.x)
            return children[6]->search(point);
          else
            return children[7]->search(point);
        }
      }


      return nullptr; // Should not reach here under normal circumstances
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

    // Function to load the octree structure from a CSV file
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

    // Serialize the OctreeNode
    void serialize(std::ofstream& ofs) const {
      // Serialize the current node
      ofs.write(reinterpret_cast<const char*>(&minBounds), sizeof(minBounds));
      ofs.write(reinterpret_cast<const char*>(&maxBounds), sizeof(maxBounds));
      ofs.write(reinterpret_cast<const char*>(&isLeaf), sizeof(isLeaf));
      ofs.write(reinterpret_cast<const char*>(&depth), sizeof(depth));

//               ofs << minBounds.x << "," << minBounds.y << "," << minBounds.z << ","
//             << maxBounds.x << "," << maxBounds.y << "," << maxBounds.z << ","
//             << isLeaf << "," << depth << std::endl;

      // Serialize the vectors
      serializeVector(ofs, corners);
      serializeVector(ofs, cornerAreas);
      serializeVector(ofs, cornerWeights);

      // Serialize children recursively
      for(OctreeNode* child : children) {
        child->serialize(ofs);
      }
    }

    // Deserialize the OctreeNode
    void deserialize(std::ifstream& ifs) {
      // Deserialize the current node
      ifs.read(reinterpret_cast<char*>(&minBounds), sizeof(minBounds));
      ifs.read(reinterpret_cast<char*>(&maxBounds), sizeof(maxBounds));
      ifs.read(reinterpret_cast<char*>(&isLeaf), sizeof(isLeaf));
      ifs.read(reinterpret_cast<char*>(&depth), sizeof(depth));

      // Deserialize the vectors
      deserializeVector(ifs, corners);
      deserializeVector(ifs, cornerAreas);
      deserializeVector(ifs, cornerWeights);

      // Delete existing children
      for(OctreeNode* child : children) {
        delete child;
      }
      children.clear();

      // Deserialize children recursively
      if(!isLeaf) {
        for(int i = 0; i < 8; ++i) {
          OctreeNode* child = new OctreeNode({0, 0, 0}, {1, 1, 1}, 0, 0, 0, nullptr);
          child->deserialize(ifs);
          children.push_back(child);
        }
      }
    }

    // Helper function to serialize a vector
    void serializeVector(std::ofstream& ofs, const std::vector<std::vector<double>>& vec) const {
      size_t size = vec.size();
      ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));

      for(const auto& innerVec : vec) {
        size_t innerSize = innerVec.size();
        ofs.write(reinterpret_cast<const char*>(&innerSize), sizeof(innerSize));
        ofs.write(reinterpret_cast<const char*>(innerVec.data()), innerSize * sizeof(double));
      }
    }

    // Helper function to deserialize a vector
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
};



template <class Type>
void generateAndLoadOctrees(const int &maxDepth, const int &degree, const double &percent, CutFemWeightParabola <double, Type> &Pweights, std::vector<OctreeNode<Type>>& loadedRoots) {

  for(int ttable = 0; ttable < 8; ++ttable) {
    std::string filename = "save/octree_table_" + std::to_string(ttable) + "_maxdepth_" + std::to_string(maxDepth) + "_per_" + std::to_string(percent) + "_degree_" +  std::to_string(degree) + ".csv";

    FILE *fp;
    fp = fopen(filename.c_str(), "r");
    if(fp != NULL) {
      fclose(fp);
    }
    else {
      cout << "creating the tables" << endl;
      OctreeNode<Type> root({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, ttable, 0, degree, &Pweights);
      if(ttable != 5) {
        root.subdivideWithRelativeError(maxDepth, percent);
      }
      else {
        root.subdivideWithRelativeError(3, 0.1);
      }
      root.saveOctreeToCSV(filename);
      std::cout << "Octree Structure:\n";
    }
  }

  // Load the octree structure and vectors from the CSV file
  loadedRoots.clear();  // Clear any existing data
  loadedRoots.reserve(8);  // Reserve space for 8 octrees
  for(int ttable = 0; ttable < 8; ++ttable) {
    loadedRoots.emplace_back(Point3D{0., 0., 0.}, Point3D{1., 1., 1.}, ttable, 0, degree, nullptr);
    loadedRoots.back().loadOctreeFromCSV("save/octree_table_" + std::to_string(ttable) + "_maxdepth_" + std::to_string(maxDepth) + "_per_" + std::to_string(percent) + "_degree_" +  std::to_string(degree) + ".csv");
  }
}






// // Define an Octree Node
// template <class Type>
// class OctreeNode {
//   public:
//     Point3D minBounds, maxBounds;
//     bool isLeaf;
//     std::vector<OctreeNode*> children;
//     std::vector<std::vector<double>> corners;//(8, std::vector<double>(3));
//     std::vector<std::vector<double>> cornerAreas;
//     std::vector<std::vector<double>> cornerWeights;
//     int table;
//     unsigned depth;
//     unsigned qM;
//     int s = 0;
//     Type a = 0;
//     double relative_error = -1;
//     double relative_error_opposite = -1;
//     CutFemWeightParabola <double, Type> *_Pweights;
//
//     OctreeNode(const Point3D& _minBounds, const Point3D& _maxBounds, const int& _table, const int& _depth, const unsigned& _qM, CutFemWeightParabola <double, Type> *Pweights)
//       : minBounds(_minBounds), maxBounds(_maxBounds), isLeaf(true), table(_table), depth(_depth), qM(_qM), _Pweights(Pweights) {}
//
//     ~OctreeNode() {};
//
//     // Function to get the eight corners of the node
//     void getCorners() {  //TODO initialize it once without using push_back
//       corners.resize(8, std::vector<double>(3));
//       corners[0] = {minBounds.x, minBounds.y, minBounds.z};
//       corners[1] = {minBounds.x, minBounds.y, maxBounds.z};
//       corners[2] = {minBounds.x, maxBounds.y, minBounds.z};
//       corners[3] = {minBounds.x, maxBounds.y, maxBounds.z};
//       corners[4] = {maxBounds.x, minBounds.y, minBounds.z};
//       corners[5] = {maxBounds.x, minBounds.y, maxBounds.z};
//       corners[6] = {maxBounds.x, maxBounds.y, minBounds.z};
//       corners[7] = {maxBounds.x, maxBounds.y, maxBounds.z};
//
// //         for (const auto& corner : corners) {
// //             std::cout << "Corner: (" << corner[0] << ", " << corner[1] << ", " << corner[2] << ") - Print Something\n";
// //         }
//       Type area(0) ;
//       Type c(1) ;
//       cornerAreas.resize(8);
//       cornerWeights.resize(8);
//       PointT <Type> p1, p2, p3 ;
//       //CutFemWeightParabola <double, Type> Pweight(QUAD, 2, "legendre");
//
//       for(size_t i = 0; i < corners.size(); ++i) {
//         const auto& corner = corners[i];
// //             std::cout << "Corner " << i << ": (" << corner[0] << ", " << corner[1] << ", " << corner[2] << ") - Print Something\n";searchPoint
//         get_p1_p2_p3(table, corner, p1, p2, p3);
//
//         int count = 0;
//         for(unsigned qq = 0; qq <= qM; qq++) {
//           for(unsigned jj = 0; jj <= qq; jj++) {
//             unsigned ii = qq - jj;
//             area = find_area_2intersection_formula(ii, jj, s, a, c, table, p1, p2, p3);
//             cornerAreas[i].push_back(static_cast<double>(area));
// //                     cout<< " f values = " << i << " "<< j << " "<< _f[count] <<endl;
//             count++;
//           }
//         }
//
//         (*_Pweights)(s, a, c, table, p1, p2, p3, cornerWeights[i]);
//
//
// //             for (size_t ii = 0; ii <= mPn; ++ii){
// //               for (size_t jj = 0; ii+jj <= mPn; ++jj){
// //                 area = find_area_2intersection_formula(jj,ii,s,a,c,table, p1,p2, p3);
// //                 corners[i].push_back(static_cast<double>(area));
// //               }
// //             }
//
//       }
//     }
//
//     void subdivideWithRelativeError(int maxDepth, double maxRelativeError, int currentDepth = 0) {
//       if(currentDepth >= maxDepth || !isLeaf) {
//         getCorners();
//         return;
//       }
//
//       getCorners();
//
//       double midX = (minBounds.x + maxBounds.x) / 2.0;
//       double midY = (minBounds.y + maxBounds.y) / 2.0;
//       double midZ = (minBounds.z + maxBounds.z) / 2.0;
//       std::vector<double> relativeErrors;
//       std::vector<double> relativeErrorsOpposite;
//       std::vector<std::vector<double>>boundary_mid_points(19);
//
//       boundary_mid_points[0] = {midX, minBounds.y, minBounds.z};
//       boundary_mid_points[1] = {maxBounds.x, midY, minBounds.z};
//       boundary_mid_points[2] = {midX, maxBounds.y, minBounds.z};
//       boundary_mid_points[3] = {minBounds.x, midY, minBounds.z};
//       boundary_mid_points[4] = {minBounds.x, minBounds.y, midZ};
//       boundary_mid_points[5] = {maxBounds.x, minBounds.y, midZ};
//       boundary_mid_points[6] = {maxBounds.x, maxBounds.y, midZ};
//       boundary_mid_points[7] = {minBounds.x, maxBounds.y, midZ};
//       boundary_mid_points[8] = {midX, minBounds.y, maxBounds.z};
//       boundary_mid_points[9] = {maxBounds.x, midY, maxBounds.z};
//       boundary_mid_points[10] = {midX, maxBounds.y, maxBounds.z};
//       boundary_mid_points[11] = {minBounds.x, midY, maxBounds.z};
//       boundary_mid_points[12] = {midX, midY, minBounds.z};
//       boundary_mid_points[13] = {minBounds.x, midY, midZ};
//       boundary_mid_points[14] = {midX, minBounds.y, midZ};
//       boundary_mid_points[15] = {midX, midY, maxBounds.z};
//       boundary_mid_points[16] = {maxBounds.x, midY, midZ};
//       boundary_mid_points[17] = {midX, maxBounds.y, midZ};
//       boundary_mid_points[18] = {midX, midY, midZ};
//
//       for(int i = 0; i < boundary_mid_points.size(); ++i) {  //TODO all nineteen points
//         // Generate a random 3D point within the bounds
// //             double randomX = minBounds.x + static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * (maxBounds.x - minBounds.x);
// //             double randomY = minBounds.y + static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * (maxBounds.y - minBounds.y);
// //             double randomZ = minBounds.z + static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * (maxBounds.z - minBounds.z);
//
//         //not random using middle boundary points
//         double randomX = boundary_mid_points[i][0];
//         double randomY = boundary_mid_points[i][1];
//         double randomZ = boundary_mid_points[i][2];
//
//         std::vector<double> interp_point = {randomX, randomY, randomZ};
//
//
//         Type f_area(0);
//         Type c = 1 ;
//         PointT<Type> p1, p2, p3;
//         get_p1_p2_p3(table, interp_point, p1, p2, p3);
//
//         std::vector<std::vector<double>>interpolation_vector(8);
//         size_t ii = 0;
//         size_t jj = 0;
// //             for (size_t ii = 0; ii <= mPn; ++ii){
// //               for (size_t jj = 0; ii+jj <= mPn; ++jj){
// //
// //                 size_t index = 0;
// //                 for (size_t mm = 0; mm < ii; ++mm){
// //                   index += mPn + 1 - mm ;
// //                 }
// //
// //                 index += jj+3 ;
//
//         int count = 0;
//         for(unsigned qq = 0; qq <= qM; qq++) {
//           for(unsigned jj = 0; jj <= qq; jj++) {
//             unsigned ii = qq - jj;
//             for(size_t ic = 0; ic < corners.size(); ++ic) {
//               interpolation_vector[ic] = {corners[ic][0], corners[ic][1], corners[ic][2], cornerAreas[ic][count]};
//             }
//             double interp_area = trilinier_interpolation(interpolation_vector, interp_point);
//             f_area = find_area_2intersection_formula(jj, ii, s, a, c, table, p1, p2, p3);
//             double formula_area = static_cast<double>(f_area);
//             double r_error = fabs(formula_area - interp_area) / formula_area;
//             double r_error_opposite = fabs(formula_area - interp_area) / (1.0 / (ii + 1) * (jj + 1) - formula_area);
//
//             relativeErrors.push_back(r_error);
//             relativeErrorsOpposite.push_back(r_error_opposite);
//
//             count++;
//           }
//         }
//
//
//
//
//
//
// //                 for (size_t ic = 0; ic < corners.size(); ++ic) {
// //                   interpolation_vector[ic] = {corners[ic][0],corners[ic][1],corners[ic][2],cornerAreas[ic][0]};
// //                 }
// //                 double interp_area = trilinier_interpolation(interpolation_vector, interp_point);
// //                 f_area = find_area_2intersection_formula(jj,ii,s,a,c,table, p1,p2, p3);
// //                 double formula_area = static_cast<double>(f_area);
// //                 double r_error = fabs(formula_area - interp_area) / formula_area;
// //                 double r_error_opposite = fabs(formula_area - interp_area) /(1.0/(ii+1)*(jj+1) - formula_area);
// //
// //                 relativeErrors.push_back(r_error);
// //                 relativeErrorsOpposite.push_back(r_error_opposite);
// //
// //
// // //               }
// // //             }
// //
// //
// // //                 f_area = find_area_2intersection_formula(m, n, s, a, c, table, p1, p2, p3);
// // //                 double formula_area = static_cast<double>(f_area);
// //
// //     //             double r_error = (m+1)*(n+1) * fabs(formula_area - interp_area) / formula_area;
// //     //             double r_error_opposite = fabs(formula_area - interp_area) /(1.0/(m+1)*(n+1) - formula_area);
// //
// // //                 double r_error = fabs(formula_area - interp_area) / formula_area;
// // //                 double r_error_opposite = fabs(formula_area - interp_area) /(1.0/(m+1)*(n+1) - formula_area);
// // //
// // //                 relativeErrors.push_back(r_error);
// // //                 relativeErrorsOpposite.push_back(r_error_opposite);
//
//
//
//       }
//       relative_error = *std::max_element(relativeErrors.begin(), relativeErrors.end());
//       relative_error_opposite = *std::max_element(relativeErrorsOpposite.begin(), relativeErrorsOpposite.end());
//
//       if(depth <= 3 || relative_error > maxRelativeError || relative_error_opposite > maxRelativeError || minBounds.y < 0.003) {
//         isLeaf = false;
//         children.push_back(new OctreeNode(minBounds, {midX, midY, midZ}, table, depth + 1, qM, _Pweights));
//         children.push_back(new OctreeNode({midX, minBounds.y, minBounds.z}, {maxBounds.x, midY, midZ}, table, depth + 1, qM, _Pweights));
//         children.push_back(new OctreeNode({minBounds.x, midY, minBounds.z}, {midX, maxBounds.y, midZ}, table, depth + 1, qM, _Pweights));
//         children.push_back(new OctreeNode({midX, midY, minBounds.z}, {maxBounds.x, maxBounds.y, midZ}, table, depth + 1, qM, _Pweights));
//         children.push_back(new OctreeNode({minBounds.x, minBounds.y, midZ}, {midX, midY, maxBounds.z}, table, depth + 1, qM, _Pweights));
//         children.push_back(new OctreeNode({midX, minBounds.y, midZ}, {maxBounds.x, midY, maxBounds.z}, table, depth + 1, qM, _Pweights));
//         children.push_back(new OctreeNode({minBounds.x, midY, midZ}, {midX, maxBounds.y, maxBounds.z}, table, depth + 1, qM, _Pweights));
//         children.push_back(new OctreeNode({midX, midY, midZ}, maxBounds, table, depth + 1, qM, _Pweights));
//
//         // Recursively check for subdivision in each child node
//         for(OctreeNode* child : children) {
//           child->subdivideWithRelativeError(maxDepth, maxRelativeError, currentDepth + 1);
//         }
//       }
//     }
//     // Recursive search for the smallest sub-cube containing a given point
//
//
//     // OctreeNode* search(const Point3D& point) {
//     //   if(isLeaf) {
//     //     return this;
//     //   }
//     //
//     //   for(OctreeNode* child : children) {
//     //     if(point.x >= child->minBounds.x && point.x <= child->maxBounds.x &&
//     //         point.y >= child->minBounds.y && point.y <= child->maxBounds.y &&
//     //         point.z >= child->minBounds.z && point.z <= child->maxBounds.z) {
//     //       return child->search(point);
//     //     }
//     //   }
//     //
//     //   return nullptr; // Should not reach here under normal circumstances
//     // }
//
//
//     OctreeNode* search(const Point3D& point) {
//       if(isLeaf) {
//         return this;
//       }
//
//
//       if(point.z <= children[0]->maxBounds.z) { // 0->3
//         if(point.y <= children[0]->maxBounds.y) { // 0->1
//           if(point.x <= children[0]->maxBounds.x)
//             return children[0]->search(point);
//           else
//             return children[1]->search(point);
//         }
//         else { // 2->3
//           if(point.x <= children[2]->maxBounds.x)
//             return children[2]->search(point);
//           else
//             return children[3]->search(point);
//         }
//       }
//       else { // 4->7
//         if(point.y <= children[4]->maxBounds.y) { // 4->5
//           if(point.x <= children[4]->maxBounds.x)
//             return children[4]->search(point);
//           else
//             return children[5]->search(point);
//         }
//         else { // 6->7
//           if(point.x <= children[6]->maxBounds.x)
//             return children[6]->search(point);
//           else
//             return children[7]->search(point);
//         }
//       }
//
//
//       return nullptr; // Should not reach here under normal circumstances
//     }
//
//     // Function to save the octree structure to a CSV file
//     void saveOctreeToCSV(const std::string& filename) const {
//       std::ofstream ofs(filename);
//       if(!ofs.is_open()) {
//         std::cerr << "Error: Unable to open file for writing: " << filename << std::endl;
//         return;
//       }
//       // Write the octree structure and vectors to the CSV file recursively
//       serialize(ofs);
//
//       ofs.close();
//     }
//
//     // Function to load the octree structure from a CSV file
//     void loadOctreeFromCSV(const std::string& filename) {
//       std::ifstream ifs(filename);
//       if(!ifs.is_open()) {
//         std::cerr << "Error: Unable to open file for reading: " << filename << std::endl;
//         return;
//       }
//
//       // Read the octree structure and vectors from the CSV file recursively
//       deserialize(ifs);
//
//       ifs.close();
//     }
//
//   private:
//     // Serialize the OctreeNode
//     void serialize(std::ofstream& ofs) const {
//       // Serialize the current node
//       ofs.write(reinterpret_cast<const char*>(&minBounds), sizeof(minBounds));
//       ofs.write(reinterpret_cast<const char*>(&maxBounds), sizeof(maxBounds));
//       ofs.write(reinterpret_cast<const char*>(&isLeaf), sizeof(isLeaf));
//       ofs.write(reinterpret_cast<const char*>(&depth), sizeof(depth));
//
// //               ofs << minBounds.x << "," << minBounds.y << "," << minBounds.z << ","
// //             << maxBounds.x << "," << maxBounds.y << "," << maxBounds.z << ","
// //             << isLeaf << "," << depth << std::endl;
//
//       // Serialize the vectors
//       serializeVector(ofs, corners);
//       serializeVector(ofs, cornerAreas);
//       serializeVector(ofs, cornerWeights);
//
//       // Serialize children recursively
//       for(OctreeNode* child : children) {
//         child->serialize(ofs);
//       }
//     }
//
//     // Deserialize the OctreeNode
//     void deserialize(std::ifstream& ifs) {
//       // Deserialize the current node
//       ifs.read(reinterpret_cast<char*>(&minBounds), sizeof(minBounds));
//       ifs.read(reinterpret_cast<char*>(&maxBounds), sizeof(maxBounds));
//       ifs.read(reinterpret_cast<char*>(&isLeaf), sizeof(isLeaf));
//       ifs.read(reinterpret_cast<char*>(&depth), sizeof(depth));
//
// //               ifs >> minBounds.x >> minBounds.y >> minBounds.z
// //             >> maxBounds.x >> maxBounds.y >> maxBounds.z
// //             >> isLeaf >> depth;
//
//       // Deserialize the vectors
//       deserializeVector(ifs, corners);
//       deserializeVector(ifs, cornerAreas);
//       deserializeVector(ifs, cornerWeights);
//
//       // Deserialize children recursively
//       if(!isLeaf) {
//         children.resize(8);
//         for(int i = 0; i < 8; ++i) {
//           OctreeNode* child = new OctreeNode({0, 0, 0}, {1, 1, 1}, 0, 0, 0, nullptr);
//           child->deserialize(ifs);
//           children[i] = child;
//         }
//       }
//     }
//
//     // Helper function to serialize a vector
//     void serializeVector(std::ofstream& ofs, const std::vector<std::vector<double>>& vec) const {
//       size_t size = vec.size();
//       ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
//
//       for(const auto& innerVec : vec) {
//         size_t innerSize = innerVec.size();
//         ofs.write(reinterpret_cast<const char*>(&innerSize), sizeof(innerSize));
//         ofs.write(reinterpret_cast<const char*>(innerVec.data()), innerSize * sizeof(double));
//       }
//     }
//
//     // Helper function to deserialize a vector
//     void deserializeVector(std::ifstream& ifs, std::vector<std::vector<double>>& vec) {
//       size_t size;
//       ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
//
//       vec.resize(size);
//       for(size_t i = 0; i < size; ++i) {
//         size_t innerSize;
//         ifs.read(reinterpret_cast<char*>(&innerSize), sizeof(innerSize));
//         vec[i].resize(innerSize);
//         ifs.read(reinterpret_cast<char*>(vec[i].data()), innerSize * sizeof(double));
//       }
//     }
// };

template <class Type>
void printOctreeStructure(OctreeNode<Type>* node, int depth = 0) {
  for(int i = 0; i < depth; ++i) {
    std::cout << "  ";
  }

  std::cout << "Node Bounds: ";
  std::cout << "(" << node->minBounds.x << ", " << node->minBounds.y << ", " << node->minBounds.z << ") to ";
  std::cout << "(" << node->maxBounds.x << ", " << node->maxBounds.y << ", " << node->maxBounds.z << ")";

  if(node->isLeaf) {
    std::cout << " relative error = " << node-> relative_error << " depth =" << node->depth << " [Leaf] : ";

//         Print the corner vectors for the leaf node
    std::cout << "  Corners:\n";

    for(int j = 0; j < node->corners.size(); ++j) {
      for(int i = 0; i < depth + 1; ++i) {
        std::cout << "  ";
      }
      const auto& corner = node->corners[j];
      const auto& cornerArea = node->cornerAreas[j];
      const auto& cornerWeight = node->cornerWeights[j];
      for(size_t k = 0; k < corner.size(); ++k) {
        const auto& entry = corner[k];
        std::cout << " " << entry;
      }
      std::cout << " Areas = ";
      for(size_t k = 0; k < cornerArea.size(); ++k) {
        const auto& entry = cornerArea[k];
        std::cout << " " << entry;
      }

      std::cout << " Weights = ";
      for(size_t k = 0; k < cornerWeight.size(); ++k) {
        const auto& entry = cornerWeight[k];
        std::cout << " " << entry;
      }
      std::cout << std::endl;
    }
  }



  /*
          std::cout << "  Corners:\n";
        for (const auto& corner : node->corners) {

            for (const auto& entry : corner) {
                std::cout << " " << entry;
            }

            std::cout << std::endl;
        }
        std::cout << " Areas:\n";
        for (const auto& corner : node->cornerAreas) {

            for (const auto& entry : corner) {
                std::cout << " " << entry;
            }
        }

        std::cout << " weights:\n";
        for (const auto& corner : node->cornerWeights) {

            for (const auto& entry : corner) {
                std::cout << " " << entry;
            }
        }

      }*/




  else {
    std::cout << " relative error = " << node-> relative_error << " " << node->depth << " [Non-Leaf]\n";
    for(OctreeNode<Type>* child : node->children) {
      printOctreeStructure(child, depth + 1);
    }
  }
}


// template <class Type>
// void generateAndLoadOctrees(const int &maxDepth, const int &degree, const double &percent, CutFemWeightParabola <double, Type> &Pweights, /*std::vector<OctreeNode<Type>>& roots,*/ std::vector<OctreeNode<Type>>& loadedRoots) {
//
//   for(int ttable = 0; ttable < 8; ++ttable) {
//     std::string filename = "save/octree_table_" + std::to_string(ttable) + "_maxdepth_" + std::to_string(maxDepth) + "_per_" + std::to_string(percent) + "_degree_" +  std::to_string(degree) + ".csv";
//
//     FILE *fp;
//     fp = fopen(filename.c_str(), "r");
//     if(fp != NULL) {
// //       std::cout << "File " << filename << " already exists. Skipping octree generation." << std::endl;
//       fclose(fp);
//     }
//     else {
//       cout << "creating the tables" << endl;
//       OctreeNode<Type> root({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, ttable, 0, degree, &Pweights);
//       if(ttable != 5) {
//         root.subdivideWithRelativeError(maxDepth, percent);
//       }
//       else {
//         root.subdivideWithRelativeError(3, 0.1);
//       }
//       root.saveOctreeToCSV(filename);
//       std::cout << "Octree Structure:\n";
//       //     roots.push_back(root);
//     }
//   }
//   // Load the octree structure and vectors from the CSV file
//   for(int ttable = 0; ttable < 8; ++ttable) {
//     OctreeNode<Type> loadedRoot({0., 0., 0.}, {1., 1., 1.}, ttable, 0, degree, nullptr);
//     loadedRoot.loadOctreeFromCSV("save/octree_table_" + std::to_string(ttable) + "_maxdepth_" + std::to_string(maxDepth) + "_per_" + std::to_string(percent) + "_degree_" +  std::to_string(degree) + ".csv");
//     loadedRoots.push_back(loadedRoot);
//   }
// }
//
// double GaussIntegral(const int &xExp, const int &yExp, const double* xg, const double* yg, const std::vector<double> &interp_point_weights, const double* gaussWeight) {
//   double Integral = 0;
//   for(unsigned ig = 0; ig < interp_point_weights.size(); ig++) {
//     Integral += pow(xg[ig], xExp) * pow(yg[ig], yExp) * interp_point_weights[ig] * gaussWeight[ig];
//   }
//   return Integral;
// }


int checkVectorRelation(const std::vector<int>& vec1, const std::vector<int>& vec2) {
  // Check if the sizes of the vectors are different
//   if(vec1.size() != vec2.size()) {
//     std::cerr << "Warning: Number of sign do not match" << std::endl;
//     return 0;
//   }

  int equalCount = 0;
  int negativeCount = 0;

  // Iterate through the vectors to count equal and negative elements
  for(size_t i = 0; i < vec1.size(); ++i) {
    if(vec1[i] == vec2[i]) {
      ++equalCount;
    }
    else {
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

void find_search_table(const PointT <double> &q1, const PointT <double> &q2, const PointT <double> &q3, unsigned &table_number, Point3D &searchP) {
  double epsilon = 1.0e-13;
  if(fabs(q1.x - 0) < epsilon) {
    if(fabs(q2.y - 1) < epsilon) {
      table_number = 0;
      searchP.x = q1.y;
      searchP.y = q2.x;
      searchP.z = q3.y;
    }
    else if(fabs(q2.x - 1) < epsilon) {
      table_number = 1;
      searchP.x = q1.y;
      searchP.y = q2.y;
      searchP.z = q3.y;
    }
    else if(fabs(q2.y - 0) < epsilon) {
      table_number = 2;
      searchP.x = q1.y;
      searchP.y = q2.x;
      searchP.z = q3.y;
    }
  }
  else if(fabs(q2.x - 0) < epsilon) {
    if(fabs(q1.y - 1) < epsilon) {
      table_number = 0;
      searchP.x = q2.y;
      searchP.y = q1.x;
      searchP.z = q3.y;
    }
    else if(fabs(q1.x - 1) < epsilon) {
      table_number = 1;
      searchP.x = q2.y;
      searchP.y = q1.y;
      searchP.z = q3.y;
    }
    else if(fabs(q1.y - 0) < epsilon) {
      table_number = 2;
      searchP.x = q2.y;
      searchP.y = q1.x;
      searchP.z = q3.y;
    }
  }
  else if(fabs(q1.x - 1) < epsilon) {
    if(fabs(q2.y - 1) < epsilon) {
      table_number = 4;
      searchP.x = q2.x;
      searchP.y = q1.y;
      searchP.z = q3.y;
    }
    else if(fabs(q2.y - 0) < epsilon) {
      table_number = 6;
      searchP.x = q1.y;
      searchP.y = q2.x;
      searchP.z = q3.y;
    }
  }
  else if(fabs(q2.x - 1) < epsilon) {
    if(fabs(q1.y - 1) < epsilon) {
      table_number = 4;
      searchP.x = q1.x;
      searchP.y = q2.y;
      searchP.z = q3.y;
    }
    else if(fabs(q1.y - 0) < epsilon) {
      table_number = 6;
      searchP.x = q2.y;
      searchP.y = q1.x;
      searchP.z = q3.y;
    }
  }
  else if(fabs(q1.y - 0) < epsilon) {
    if(fabs(q2.y - 1) < epsilon) {
      table_number = 5;
      searchP.x = q2.x;
      searchP.y = q1.x;
      searchP.z = q3.y;
    }
    else if(fabs(q2.y - 0) < epsilon) {
      table_number = 7;
      searchP.x = q1.x;
      searchP.y = q2.x;
      searchP.z = q3.y;
    }
  }
  else if(fabs(q2.y - 0) < epsilon) {
    if(fabs(q1.y - 1) < epsilon) {
      table_number = 5;
      searchP.x = q1.x;
      searchP.y = q2.x;
      searchP.z = q3.y;
    }
  }
  else if(fabs(q1.y - 1) < epsilon) {
    if(fabs(q2.y - 1) < epsilon) {
      table_number = 3;
      searchP.x = q1.x;
      searchP.y = q2.x;
      searchP.z = q3.y;
    }
  }
}







#endif // PARABOLAINTEGRATION_HPP_INCLUDED
