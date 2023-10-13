
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
Parabola <Type>  get_parabola_equation( PointT <Type> p1, PointT <Type> p2, PointT <Type> p3 , Type det) {
    Type x1 = p1.x, x2 = p2.x, x3 = p3.x;
    Type y1 = p1.y, y2 = p2.y, y3 = p3.y;
    Type k,b,d;

    if(det !=0){
       k = -(y1 * (x2 - x3) + y2 * (x3 - x1) + y3 * (x1 - x2)) / det;
       b = -(y1 * (x3*x3 - x2*x2) + y2 * (x1*x1 - x3*x3)+ y3 * ((x2*x2 - x1*x1))) / det;
       d = -(y1 * x2 * x3 * (x2 - x3) + y2 * x3 * x1 * (x3 - x1) + y3 * x1 * x2 * (x1 - x2)) / det;
    }
    else {
       k = -0.;
       b = -1.;
       d = p1.x ;
    }

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

          term *= pterm * r_p1_m2p * (r_p1_m2p + 1) / (p * rMax_mr_pp);
          sum1 += term * pow(pol1[1], r - 2 * p);
        }
        sum1 = sum1 / (factorial<Type>(r) * factorial<Type>(rMax - r));
        sum2 = (r == rMax) ? 0 : sum1 * pow(pol1[0], rMax - r);
        sum1 *= pow(pol1[2], rMax - r);
//         for(unsigned i = 0; i < I2.size(); i++) {

        A2 += sum1 * diff_x_pow[r_pm_p1] +  sum2 * diff_x_pow[rMax_mr_pm_p1] ;

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

        // BEGIN pre evalution of power of U
        std::vector <Type> diff_u_pow(2 * s + 2 * n + m + 3, 0) ;
        Type u1pi = 1. / pow(u1, n);
        Type u2pi = 1. / pow(u2, n);
        for(int pwr = 0; pwr <= n - 1 ; pwr++, u1pi *= u1, u2pi *= u2) {
          int actual_pwr = pwr - n;
          diff_u_pow[pwr] = (u2pi - u1pi) / actual_pwr ;
        }


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

            sum2 += B[q] * diff_u_pow[2 * qMax + p - q];
          }
          A2i += (sum1 + sum2) * pow(-c, m - p) / (factorial<Type>(p) * factorial<Type>(m - p));
        }
        A2 += A2i * ((n % 2 == 0) ? -1 : 1) */* pow(-1, n + 1) **/ factorial<Type>(n) * factorial<Type>(m) / pow(a, m + 1) ;
      }

    }
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


            A[q] *= pow(k[1], q) / (factorial<Type>(q) * factorial<Type>(qMax - q));
            B[q] = A[q] * pow(k[0], qMax - q);
            A[q] *= pow(k[2], qMax - q);

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

// cout<< " " << " left " << left << " top "<< top << " right "<< right << " bottom " << bottom  << " table number :"<< table_number << " number of intersection " << intersect_number <<endl;

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

            parabola = get_parabola_equation(p1, p2, p3, det);

            CheckIntersection <Type> (intersect_number, table_number, intersection, interp_point, parabola);


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

                    pol1[0] = static_cast<Type>(0);
                    pol1[1] = static_cast<Type>(-1);
                    pol1[2] = p1.x;
                    pol2[0] = pol1[0];
                    pol2[1] = pol1[1];
                    pol2[2] = p1.x;
                    c=static_cast<Type>(0);

                }
            }
            else if (intersect_number > 2){
                Type slope = (p3.y-p1.y)/(p3.x-p1.x);
                c=static_cast<Type>(1);
                pol2[0] = static_cast<Type>(0);    //k=0
                pol2[1] = -slope;
                pol2[2] = slope*p1.x - p1.y ;
                pol1[0] = pol2[0];    //k=0
                pol1[1] = pol2[1];
                pol1[2] = pol2[2] + c ;


            }
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
            count ++ ;
//           }         //det ==0 ends
        }
      }
    }
  }
cout << "parabola table created " << parabola_table[4][6][5][7][0][14] << endl;
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

            parabola = get_parabola_equation(p1, p2, p3, det);

            CheckIntersection <Type> (intersect_number, table_number , intersection, interp_point , parabola);



                if (det !=0){
                    pol1[0] = parabola.k;
                    pol1[1] = parabola.b;
                    pol1[2] = parabola.d + c;
                    pol2[0] = parabola.k;
                    pol2[1] = parabola.b;
                    pol2[2] = parabola.d;
                }


//             if (intersect_number < 4) cout << " ----------- determinant = " << det <<" intersection number is something other than 4 check please " << " =/////////////////============///////////////////=======================/////////////////////===============////////////// intersection point = "  << intersect_number<< endl;


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

                 parabola_table_4intersection[table][i1][i2][i3][normal][14] = parabola_table_4intersection[table][i1][i2][i3-1][normal][14] ;
                }
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

//     cout << " original parabola " << k <<"x^2 + "<< b <<"x + "<< d << " + " << c <<"y = 0" << endl;

    bool left=0 , top = 0 , right = 0 , bottom = 0;

//     cout << "normal before = " << normal ;
    if (c<0) normal = (normal+1) % 2 ;
//     cout << " normal after = " << normal <<endl;
    k = k/c; b=b/c; d=d/c; c=static_cast<Type>(1);  //TODO 1 4/24 sign of c is different.



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


        interp_table[0] = {x0,y0,z0,-1 };
        interp_table[1] = {x0,y0,z1,-2 };
        interp_table[2] = {x0,y1,z0,-3 };
        interp_table[3] = {x0,y1,z1,-4 };
        interp_table[4] = {x1,y0,z0,-5 };
        interp_table[5] = {x1,y0,z1,-6 };
        interp_table[6] = {x1,y1,z0,-7 };
        interp_table[7] = {x1,y1,z1,-8 };





      if(intersect_number == 2){
//         cout << " Table used - 2 intersection " << endl;


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

        for(unsigned int i = 0; i <=7 ; i++){

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

//         cout << " modified_normal = "<< modified_normal << endl;
        if (normal == 1){k*=-1;b*=-1;d*=-1;c*=-1;}

        std::vector< Type > grad_given_par{2*k*intersection[4]+b , c}; // will it work like this?



        std::vector< Type > grad_interp_par(2);

        grad_interp_par = {2*parabola_table[table_number][i1_0][i2_0][i3_0][normal][10]*  intersection[4]   +  parabola_table[table_number][i1_0][i2_0][i3_0][normal][11] , parabola_table[table_number][i1_0][i2_0][i3_0][normal][13]};

        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;
        //TODO check what happens if it is 0.

          interp_table[0] = {parabola_table[table_number][i1_0][i2_0][i3_0][modified_normal][1],parabola_table[table_number][i1_0][i2_0][i3_0][modified_normal][2],parabola_table[table_number][i1_0][i2_0][i3_0][modified_normal][3],parabola_table[table_number][i1_0][i2_0][i3_0][modified_normal][14] };

        grad_interp_par = {2*parabola_table[table_number][i1_0][i2_0][i3_1][normal][10]*  intersection[4]   +  parabola_table[table_number][i1_0][i2_0][i3_1][normal][11] , parabola_table[table_number][i1_0][i2_0][i3_1][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;


          interp_table[1] = {parabola_table[table_number][i1_0][i2_0][i3_1][modified_normal][1],parabola_table[table_number][i1_0][i2_0][i3_1][modified_normal][2],parabola_table[table_number][i1_0][i2_0][i3_1][modified_normal][3],parabola_table[table_number][i1_0][i2_0][i3_1][modified_normal][14]};

        grad_interp_par = {2*parabola_table[table_number][i1_0][i2_1][i3_0][normal][10]*  intersection[4]   +  parabola_table[table_number][i1_0][i2_1][i3_0][normal][11] , parabola_table[table_number][i1_0][i2_1][i3_0][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;


          interp_table[2] = {parabola_table[table_number][i1_0][i2_1][i3_0][modified_normal][1],parabola_table[table_number][i1_0][i2_1][i3_0][modified_normal][2],parabola_table[table_number][i1_0][i2_1][i3_0][modified_normal][3],parabola_table[table_number][i1_0][i2_1][i3_0][modified_normal][14] };

        grad_interp_par = {2*parabola_table[table_number][i1_0][i2_1][i3_1][normal][10]* intersection[4]   +  parabola_table[table_number][i1_0][i2_1][i3_1][normal][11] , parabola_table[table_number][i1_0][i2_1][i3_1][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;


          interp_table[3] = {parabola_table[table_number][i1_0][i2_1][i3_1][modified_normal][1],parabola_table[table_number][i1_0][i2_1][i3_1][modified_normal][2],parabola_table[table_number][i1_0][i2_1][i3_1][modified_normal][3],parabola_table[table_number][i1_0][i2_1][i3_1][modified_normal][14] };

        grad_interp_par = {2*parabola_table[table_number][i1_1][i2_0][i3_0][normal][10]* intersection[4]   +  parabola_table[table_number][i1_1][i2_0][i3_0][normal][11] , parabola_table[table_number][i1_1][i2_0][i3_0][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;


          interp_table[4] = {parabola_table[table_number][i1_1][i2_0][i3_0][modified_normal][1],parabola_table[table_number][i1_1][i2_0][i3_0][modified_normal][2],parabola_table[table_number][i1_1][i2_0][i3_0][modified_normal][3],parabola_table[table_number][i1_1][i2_0][i3_0][modified_normal][14] };

        grad_interp_par = {2*parabola_table[table_number][i1_1][i2_0][i3_1][normal][10]* intersection[4]  +  parabola_table[table_number][i1_1][i2_0][i3_1][normal][11] , parabola_table[table_number][i1_1][i2_0][i3_1][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;


          interp_table[5] = {parabola_table[table_number][i1_1][i2_0][i3_1][modified_normal][1],parabola_table[table_number][i1_1][i2_0][i3_1][modified_normal][2],parabola_table[table_number][i1_1][i2_0][i3_1][modified_normal][3],parabola_table[table_number][i1_1][i2_0][i3_1][modified_normal][14] };

        grad_interp_par = {2*parabola_table[table_number][i1_1][i2_1][i3_0][normal][10]*  intersection[4]   +  parabola_table[table_number][i1_1][i2_1][i3_0][normal][11] , parabola_table[table_number][i1_1][i2_1][i3_0][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;


          interp_table[6] = {parabola_table[table_number][i1_1][i2_1][i3_0][modified_normal][1],parabola_table[table_number][i1_1][i2_1][i3_0][modified_normal][2],parabola_table[table_number][i1_1][i2_1][i3_0][modified_normal][3],parabola_table[table_number][i1_1][i2_1][i3_0][modified_normal][14] };

        grad_interp_par = {2*parabola_table[table_number][i1_1][i2_1][i3_1][normal][10]* intersection[4]   +  parabola_table[table_number][i1_1][i2_1][i3_1][normal][11] , parabola_table[table_number][i1_1][i2_1][i3_1][normal][13]};
        modified_normal = (grad_given_par[0]*grad_interp_par[0] + grad_given_par[1]*grad_interp_par[1] < 0)? (normal+1)%2 : normal ;



          interp_table[7] = {parabola_table[table_number][i1_1][i2_1][i3_1][modified_normal][1],parabola_table[table_number][i1_1][i2_1][i3_1][modified_normal][2],parabola_table[table_number][i1_1][i2_1][i3_1][modified_normal][3],parabola_table[table_number][i1_1][i2_1][i3_1][modified_normal][14] };

       }



        }


      else if (intersect_number>2){
        cout << " Table used - 4 intersection " << endl;

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


        interp_table[0] = {parabola_table_4intersection[table_number][i1_0][i2_0][i3_0][normal][1],parabola_table_4intersection[table_number][i1_0][i2_0][i3_0][normal][2],parabola_table_4intersection[table_number][i1_0][i2_0][i3_0][normal][3],parabola_table_4intersection[table_number][i1_0][i2_0][i3_0][normal][14] };
        interp_table[1] = {parabola_table_4intersection[table_number][i1_0][i2_0][i3_1][normal][1],parabola_table_4intersection[table_number][i1_0][i2_0][i3_1][normal][2],parabola_table_4intersection[table_number][i1_0][i2_0][i3_1][normal][3],parabola_table_4intersection[table_number][i1_0][i2_0][i3_1][normal][14]};
        interp_table[2] = {parabola_table_4intersection[table_number][i1_0][i2_1][i3_0][normal][1],parabola_table_4intersection[table_number][i1_0][i2_1][i3_0][normal][2],parabola_table_4intersection[table_number][i1_0][i2_1][i3_0][normal][3],parabola_table_4intersection[table_number][i1_0][i2_1][i3_0][normal][14] };
        interp_table[3] = {parabola_table_4intersection[table_number][i1_0][i2_1][i3_1][normal][1],parabola_table_4intersection[table_number][i1_0][i2_1][i3_1][normal][2],parabola_table_4intersection[table_number][i1_0][i2_1][i3_1][normal][3],parabola_table_4intersection[table_number][i1_0][i2_1][i3_1][normal][14] };
        interp_table[4] = {parabola_table_4intersection[table_number][i1_1][i2_0][i3_0][normal][1],parabola_table_4intersection[table_number][i1_1][i2_0][i3_0][normal][2],parabola_table_4intersection[table_number][i1_1][i2_0][i3_0][normal][3],parabola_table_4intersection[table_number][i1_1][i2_0][i3_0][normal][14] };
        interp_table[5] = {parabola_table_4intersection[table_number][i1_1][i2_0][i3_1][normal][1],parabola_table_4intersection[table_number][i1_1][i2_0][i3_1][normal][2],parabola_table_4intersection[table_number][i1_1][i2_0][i3_1][normal][3],parabola_table_4intersection[table_number][i1_1][i2_0][i3_1][normal][14] };
        interp_table[6] = {parabola_table_4intersection[table_number][i1_1][i2_1][i3_0][normal][1],parabola_table_4intersection[table_number][i1_1][i2_1][i3_0][normal][2],parabola_table_4intersection[table_number][i1_1][i2_1][i3_0][normal][3],parabola_table_4intersection[table_number][i1_1][i2_1][i3_0][normal][14] };
        interp_table[7] = {parabola_table_4intersection[table_number][i1_1][i2_1][i3_1][normal][1],parabola_table_4intersection[table_number][i1_1][i2_1][i3_1][normal][2],parabola_table_4intersection[table_number][i1_1][i2_1][i3_1][normal][3],parabola_table_4intersection[table_number][i1_1][i2_1][i3_1][normal][14] };


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
  //this gives an error error if x0 and ceil giving the same result .

//   check, here x1-x0 should always be 1/partition ? right partitionDistance = x1-x0 = y1-y0 = z1-z0;

//     cout << " x0 x1 " << x0 << " " << x1 <<" " << y0 << " " << y1 << " " << z0 << " " << z1 <<endl;
  Type partitionDistance = static_cast<Type>(1)/partition;

  //Type x_d = (x-x0)/(x1-x0);

  // x_d = log_(x1/x0)(x/x0);

  // y_d = (y1-y0)* x_d + y0 ;

  // y_d =



//   Type y_d = (y-y0)/(y1-y0);
//   Type z_d = (z-z0)/(z1-z0);


  Type x_d = (x-x0)/partitionDistance;
  Type y_d = (y-y0)/partitionDistance;
  Type z_d = (z-z0)/partitionDistance;


  Type c_00 = interp_table[0][3] * (1-x_d) + interp_table[4][3] * x_d ;
  Type c_01 = interp_table[1][3] * (1-x_d) + interp_table[5][3] * x_d ;
  Type c_10 = interp_table[2][3] * (1-x_d) + interp_table[6][3] * x_d ;
  Type c_11 = interp_table[3][3] * (1-x_d) + interp_table[7][3] * x_d ;

  Type c_0 = c_00 * (1-y_d) + c_10 * y_d ;
  Type c_1 = c_01 * (1-y_d) + c_11 * y_d ;

  Type cc = c_0 * (1-z_d) + c_1 * z_d ;

  return cc;
}




struct Point {
    double x;
    double y;
};


struct Circle {
    Point center;
    double radius;
};

// Function to find intersection points of a circle with a line

std::vector<Point> findGridIntersectionPoints(Circle circle, Point lineStart, Point lineEnd , std::vector<int> &ltrbNumber , unsigned int &table_number ) {


    table_number =9999;
    double epsilon = 0.000000000000001; // put it in all the inequality.

    std::vector<Point> grid_intersections(0);
    // Circle equation: (x - h)^2 + (y - k)^2 = r^2
    double h = circle.center.x;
    double k = circle.center.y;
    double r = circle.radius;


    double diff_of_sqr , sqrt_diff_of_sqrt , value_p , value_m ;


    //x = i (Left)and circle
    diff_of_sqr = r*r - (lineStart.x - h)*(lineStart.x - h);
    if (diff_of_sqr > epsilon ){
       sqrt_diff_of_sqrt = sqrt(diff_of_sqr);
       value_p = k + sqrt_diff_of_sqrt;
       value_m = k - sqrt_diff_of_sqrt;
//        cout << value_m << " -- " << lineStart.y << " " << lineEnd.y ;
      if (value_p + epsilon >= lineStart.y && value_p <= lineEnd.y + epsilon){
        grid_intersections.resize(grid_intersections.size()+1);
        grid_intersections[grid_intersections.size()-1] = {lineStart.x,value_p};
        ltrbNumber[0] += 1;// it means the circle intersected the left.
      }
      if (value_m + epsilon >= lineStart.y && value_m <= lineEnd.y+ epsilon){
//         cout << " !# ";
        grid_intersections.resize(grid_intersections.size()+1);
        grid_intersections[grid_intersections.size()-1] = {lineStart.x, value_m};
        ltrbNumber[0] += 1;
      }
    }
    else if (fabs(diff_of_sqr - 0) <= epsilon){
        if (k + epsilon >= lineStart.y && k <= lineEnd.y + epsilon){
          grid_intersections.resize(grid_intersections.size()+1);
          grid_intersections[grid_intersections.size()-1] = {lineStart.x,k};
          ltrbNumber[0] += 1;  // it means the circle intersected the left.
        }
    }

    //y = j+1 (Top) and circle  corner points won't register here because of strickly less than.
    diff_of_sqr = r*r - (lineEnd.y - k)*(lineEnd.y - k);
    if (diff_of_sqr > epsilon ){
       sqrt_diff_of_sqrt = sqrt(diff_of_sqr);
       value_p = h + sqrt_diff_of_sqrt;
       value_m = h - sqrt_diff_of_sqrt;
      if (value_p > lineStart.x + epsilon && value_p + epsilon < lineEnd.x){
        grid_intersections.resize(grid_intersections.size()+1);
        grid_intersections[grid_intersections.size()-1] = {value_p, lineEnd.y};
        ltrbNumber[1] += 1;
      }
      if (value_m > lineStart.x + epsilon && value_m + epsilon < lineEnd.x){
        grid_intersections.resize(grid_intersections.size()+1);
        grid_intersections[grid_intersections.size()-1] = {value_m, lineEnd.y};
        ltrbNumber[1] += 1;
        cout << " use negative value in top. value m = " << value_m << "diff "<< value_m - lineEnd.x <<endl;
      }
    }
      else if (fabs(diff_of_sqr - 0) <= epsilon){
        if (h > lineStart.x + epsilon && h + epsilon < lineEnd.x){
          grid_intersections.resize(grid_intersections.size()+1);
          grid_intersections[grid_intersections.size()-1] = {h,lineEnd.y};
          ltrbNumber[1] += 1;  // it means the circle intersected the left.
        }
    }

    //x = i+1 (Right) and circle
    diff_of_sqr = r*r - (lineEnd.x - h)*(lineEnd.x - h);
    if (diff_of_sqr > epsilon ){
       sqrt_diff_of_sqrt = sqrt(diff_of_sqr);
       value_p = k + sqrt_diff_of_sqrt;
       value_m = k - sqrt_diff_of_sqrt;
//        cout << value_m << " -- " <<value_p<< " " <<  lineStart.y << " " << lineEnd.y ;
      if (value_p + epsilon>=lineStart.y && value_p <= lineEnd.y+ epsilon){
        grid_intersections.resize(grid_intersections.size()+1);
        grid_intersections[grid_intersections.size()-1] = {lineEnd.x, value_p};
        ltrbNumber[2] += 1;
      }
      if (value_m + epsilon>=lineStart.y && value_m <= lineEnd.y+ epsilon){
        grid_intersections.resize(grid_intersections.size()+1);
        grid_intersections[grid_intersections.size()-1] = {lineEnd.x, value_m};
        ltrbNumber[2] += 1;
      }
    }
    else if (fabs(diff_of_sqr - 0) <= epsilon){
        if (k + epsilon >=lineStart.y && k <= lineEnd.y+ epsilon){
          grid_intersections.resize(grid_intersections.size()+1);
          grid_intersections[grid_intersections.size()-1] = {lineEnd.x,k};
          ltrbNumber[2] += 1;  // it means the circle intersected the left.
        }
    }

    //y = j (bottom) and circle
    diff_of_sqr = r*r - (lineStart.y - k)*(lineStart.y - k);
    if (diff_of_sqr > epsilon ){
       sqrt_diff_of_sqrt = sqrt(diff_of_sqr);
       value_p = h + sqrt_diff_of_sqrt;
       value_m = h - sqrt_diff_of_sqrt;
      if (value_p > lineStart.x + epsilon && value_p + epsilon < lineEnd.x){
        grid_intersections.resize(grid_intersections.size()+1);
        grid_intersections[grid_intersections.size()-1] = {value_p, lineStart.y};
        ltrbNumber[3] += 1;
      }
      if (value_m > lineStart.x + epsilon && value_m + epsilon < lineEnd.x){
        grid_intersections.resize(grid_intersections.size()+1);
        grid_intersections[grid_intersections.size()-1] = {value_m, lineStart.y};
        ltrbNumber[3] += 1;
      }
    }
    else if (fabs(diff_of_sqr - 0) <= epsilon){
        if (h > lineStart.x+epsilon && h+ epsilon < lineEnd.x){
          grid_intersections.resize(grid_intersections.size()+1);
          grid_intersections[grid_intersections.size()-1] = {h,lineStart.y};
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


    if (grid_intersections.size() == 1){ // it means it touches the corner never enters.
      grid_intersections.resize(0);
    }

    else if(grid_intersections.size() == 3){ // TODO we have to get rid of one corner.
     cout << "######### get rid of a corner ########## -----------------------------------------------------------------------------------------------" <<endl;
    }

    else if(grid_intersections.size() == 2){
      grid_intersections.resize(3);
      grid_intersections[2] = { (grid_intersections[0].x + grid_intersections[1].x)*0.5 , -1.0};
      diff_of_sqr = r*r - ((grid_intersections[2].x - h)*(grid_intersections[2].x - h));
      if (diff_of_sqr >= 0 ){ //it is different for table 5
        sqrt_diff_of_sqrt = sqrt(diff_of_sqr);
        value_p = k + sqrt_diff_of_sqrt;
        value_m = k - sqrt_diff_of_sqrt;
        if (value_p >=lineStart.y && value_p <= lineEnd.y){
          grid_intersections[2].y = value_p;
        }
        else if (value_m >=lineStart.y && value_m <= lineEnd.y ){
          grid_intersections[2].y = value_m;
        }
        else {
          grid_intersections.resize(0);
          ltrbNumber={0,0,0,0};
          // meaning we have two intersection but mid_point is outside. it is either 4 intersection or never inters}
        }
      }


      if (ltrbNumber[0] == 1){
        if (ltrbNumber[1] == 1) {
          table_number = 0;
        }
        else if(ltrbNumber[2] == 1){
          table_number = 1;
        }
        else if(ltrbNumber[3] == 1){
          table_number = 2;
        }
      }

      else if (ltrbNumber[1] == 2) table_number = 3;
      else if (ltrbNumber[1] == 1){
        if(ltrbNumber[2] == 1) table_number = 4;
        else if(ltrbNumber[3] == 1) table_number = 5;
      }
      else if (ltrbNumber[2] == 1 && ltrbNumber[3] == 1) table_number = 6;
      else if (ltrbNumber[3] == 2) table_number = 7;

//       swapping x and y when we have Top - Bottom (table-5) into table-1

      cout<< " ("<<grid_intersections[0].x<<"," << grid_intersections[0].y<<") " << " ("<<grid_intersections[1].x<<"," << grid_intersections[1].y<<") " << " ("<<grid_intersections[2].x<<"," << grid_intersections[2].y<<") " <<endl;

      if (table_number == 5){
        grid_intersections[2] = { -1.0, (grid_intersections[0].y + grid_intersections[1].y)*0.5 };
        diff_of_sqr = r*r - ((grid_intersections[2].y - h)*(grid_intersections[2].y - h));
        if (diff_of_sqr >= 0 ){
          sqrt_diff_of_sqrt = sqrt(diff_of_sqr);
          value_p = h + sqrt_diff_of_sqrt;
          value_m = h - sqrt_diff_of_sqrt;
          if (value_p >=lineStart.x && value_p <= lineEnd.x){
            grid_intersections[2].x = value_p;
          }
          else if (value_m >=lineStart.x && value_m <= lineEnd.x ){
            grid_intersections[2].x = value_m;
          }
          else {
          grid_intersections.resize(0);
          ltrbNumber={0,0,0,0};
          // meaning we have two intersection but mid_point is outside. it is either 4 intersection or never inters}
          }
        }
      }

    }




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
std::vector<Point> mapIntersectionPoints(std::vector<Point>grid_intersections, Point gridStart, Point gridEnd, unsigned int &table_number , std::vector< Type > &interp_point, bool &isswap){

  std::vector<Point> mapped_intersections(0);
  interp_point.resize(0);
  for(int aa=0 ; aa < grid_intersections.size(); aa++){
        mapped_intersections.resize(mapped_intersections.size()+1);
        mapped_intersections[mapped_intersections.size()-1] = {((grid_intersections[aa].x - gridStart.x)/(gridEnd.x -gridStart.x)), ((grid_intersections[aa].y - gridStart.y)/(gridEnd.y - gridStart.y))};
        if (mapped_intersections[aa].x < 0) mapped_intersections[aa].x = 0;
        else if (mapped_intersections[aa].x > 1) mapped_intersections[aa].x = 1;
        if (mapped_intersections[aa].y < 0) mapped_intersections[aa].y = 0;
        else if (mapped_intersections[aa].y >1) mapped_intersections[aa].x = 1;
//       mapped_intersections.push_back({((intersections[aa].x - gridStart.x)/(gridEnd.x -gridStart.x)) , ((intersections[aa].y - gridStart.y)/           //       (gridEnd.y -gridStart.y))});

  }
for(int mk=0;mk<mapped_intersections.size();mk++){
  cout<< "before swapping : ("<<mapped_intersections[mk].x << " " <<mapped_intersections[mk].y << ") "<<endl;
}


    if(table_number == 5){
      cout << " swapping values :::::::::" <<endl;
      isswap = 1;
      double swapvalue;
      swapvalue = mapped_intersections[0].x ;
      mapped_intersections[0].x = mapped_intersections[0].y;
      mapped_intersections[0].y = swapvalue;

      swapvalue = mapped_intersections[1].x ;
      mapped_intersections[1].x = mapped_intersections[1].y;
      mapped_intersections[1].y = swapvalue;

      swapvalue = mapped_intersections[2].x ;
      mapped_intersections[2].x = mapped_intersections[2].y;
      mapped_intersections[2].y = swapvalue;

      table_number = 1;
      for(int mk=0;mk<mapped_intersections.size();mk++){
        cout<< "after swapping : ("<<mapped_intersections[mk].x << " " <<mapped_intersections[mk].y << ") "<<endl;
      }
    }

  if(mapped_intersections.size() == 3){ // Finding interpolation points
    switch (table_number) {
                case 0:
                  interp_point.resize(3);
                  interp_point = { static_cast<Type>(mapped_intersections[0].y) , static_cast<Type>(mapped_intersections[1].x) , static_cast<Type>(mapped_intersections[2].y)};
                    break;
                case 1:
                  interp_point.resize(3);
                  interp_point = { static_cast<Type>(mapped_intersections[0].y) , static_cast<Type>(mapped_intersections[1].y) , static_cast<Type>(mapped_intersections[2].y)};
                    break;
                case 2:
                  interp_point.resize(3);
                  interp_point = { static_cast<Type>(mapped_intersections[0].y) , static_cast<Type>(mapped_intersections[1].x) , static_cast<Type>(mapped_intersections[2].y)};
                    break;
                case 3:
                  interp_point.resize(3);
                  interp_point = { static_cast<Type>(mapped_intersections[0].x) , static_cast<Type>(mapped_intersections[1].x) , static_cast<Type>(mapped_intersections[2].y)};
                    break;
                case 4:
                  interp_point.resize(3);
                  interp_point = { static_cast<Type>(mapped_intersections[0].x) , static_cast<Type>(mapped_intersections[1].y) , static_cast<Type>(mapped_intersections[2].y)};
                    break;
                case 5:
                  interp_point.resize(3);
                  interp_point = { static_cast<Type>(mapped_intersections[0].x) , static_cast<Type>(mapped_intersections[1].x) , static_cast<Type>(mapped_intersections[2].y)};
                  break;
                case 6:
                  interp_point.resize(3);
                  interp_point = { static_cast<Type>(mapped_intersections[0].y) , static_cast<Type>(mapped_intersections[1].x) , static_cast<Type>(mapped_intersections[2].y)};
                  break;
                case 7:
                  interp_point.resize(3);
                  interp_point = { static_cast<Type>(mapped_intersections[0].x) , static_cast<Type>(mapped_intersections[1].x) , static_cast<Type>(mapped_intersections[2].y)};
                  break;
    }
  }


  return mapped_intersections;
}

template <class Type>   // we only use this if 2 intersection .
void interpolationTable(std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< Type >>>>>> &parabola_table, std::vector< Type > &interp_point , std::vector< std::vector< Type >> & interp_table, const int & partition, const unsigned int &table_number , int &normal){


        int i1_0,i1_1,i2_0,i2_1,i3_0,i3_1 ;
        Type x0,x1,y0,y1,z0,z1;
        x0 = floor(interp_point[0]* partition) / partition;
        x1 = ceil(interp_point[0]* partition) / partition;
        y0 = floor(interp_point[1]* partition) / partition;
        y1 = ceil(interp_point[1]* partition) / partition;
        z0 = floor(interp_point[2]* partition) / partition;
        z1 = ceil(interp_point[2]* partition) / partition;
        cout << "floors 1 " << x0 <<" " << x1<<" "<< y0<<" "<< y1<<" "<< z0<<" " <<z1<< " " << endl;

        interp_table[0] = {x0,y0,z0,-1 };
        interp_table[1] = {x0,y0,z1,-2 };
        interp_table[2] = {x0,y1,z0,-3 };
        interp_table[3] = {x0,y1,z1,-4 };
        interp_table[4] = {x1,y0,z0,-5 };
        interp_table[5] = {x1,y0,z1,-6 };
        interp_table[6] = {x1,y1,z0,-7 };
        interp_table[7] = {x1,y1,z1,-8 };

        cout<< " interpolation points " << interp_point[0] << " " << interp_point[1]<< " " << interp_point[2] << " " <<endl;

        i1_0 = floor(static_cast<double>(interp_point[0] * partition)) ;
        i1_1 = ceil (static_cast<double>(interp_point[0] * partition)) ;
        i2_0 = floor(static_cast<double>(interp_point[1] * partition)) ;
        i2_1 = ceil (static_cast<double>(interp_point[1] * partition)) ;
        i3_0 = floor(static_cast<double>(interp_point[2] * partition)) ;
        i3_1 = ceil (static_cast<double>(interp_point[2] * partition)) ;

                cout << "table number"<< table_number <<" int floor " << i1_0 << " "<< i1_1<<" "<< i2_0 << " "<< i2_1<<" "<< i3_0 << " "<< i3_1<<" " << endl;

//         i1_0 = floor(interp_point[0] * partition) ;
//         i1_1 = ceil (interp_point[0] * partition) ;
//         i2_0 = floor(interp_point[1] * partition) ;
//         i2_1 = ceil (interp_point[1] * partition) ;
//         i3_0 = floor(interp_point[2] * partition) ;
//         i3_1 = ceil (interp_point[2] * partition) ;
/*
        cout<< "type of i1_0 "<< typeid(i1_0).name() << endl;

        cout << "table number"<< table_number <<" int floor " << i1_0 << " "<< i1_1<<" "<< i2_0 << " "<< i2_1<<" "<< i3_0 << " "<< i3_1<<" " << endl;

        cout<< parabola_table[4][6][5][7][0][14] << endl;
        cout<< parabola_table[table_number][i1_0][i2_0][i3_0][normal][14] << "what  happening "<<endl;*/




        interp_table[0][3] = parabola_table[table_number][i1_0][i2_0][i3_0][normal][14] ;
        interp_table[1][3] = parabola_table[table_number][i1_0][i2_0][i3_1][normal][14] ;
        interp_table[2][3] = parabola_table[table_number][i1_0][i2_1][i3_0][normal][14] ;
        interp_table[3][3] = parabola_table[table_number][i1_0][i2_1][i3_1][normal][14] ;
        interp_table[4][3] = parabola_table[table_number][i1_1][i2_0][i3_0][normal][14] ;
        interp_table[5][3] = parabola_table[table_number][i1_1][i2_0][i3_1][normal][14] ;
        interp_table[6][3] = parabola_table[table_number][i1_1][i2_1][i3_0][normal][14] ;
        interp_table[7][3] = parabola_table[table_number][i1_1][i2_1][i3_1][normal][14] ;

        for(int ik=0; ik <= 7;ik++){
         cout<< " ( "<<interp_table[ik][0] << ", " <<interp_table[ik][1] << ", " <<interp_table[ik][2] << ", " <<interp_table[ik][3] << ") " <<endl;
        }
}




template <class Type>
Type findGridArea(const Circle &circle, const Point &gridStart, const Point &gridEnd, std::vector<Point> grid_intersections, std:: vector <std:: vector <Type>> interp_table, const std::vector< Type > &interp_point, const int &partition, const double &grid_size){
  Type area(0);
  if(grid_intersections.size() <= 1){
    // checking wherher the grid is outside or inside the circle . Circle equation: (x - h)^2 + (y - k)^2 = r^2 Find (i+epsilon-h)^2 +(j+epsilon-k)^2 < r^2
    area = ( (gridStart.x - circle.center.x )*(gridStart.x - circle.center.x) + (gridStart.y - circle.center.y)*(gridStart.y - circle.center.y) < (circle.radius*circle.radius)) ? static_cast<Type>(1) : static_cast<Type>(0) ;
  }
  else if (interp_point.size() == 3){ // here goes our ex_23

      cout << " using interpolation " << endl;

      area = trilinier_interpolation<Type>(interp_table , interp_point, partition) ;

//       cout << "grid area before scaling = "<< area <<endl;
  }

  else{} // direct calculation for 4 intersections. use grid intersections to find the area of for intersection points.

//   cout << "grid size = "<< grid_size <<endl;
  area = area * grid_size * grid_size ;
//   cout << "grid area after scaling = "<< area <<endl;

  return area;


}



int main() {

  typedef cpp_bin_float_oct Type;      //     typedef double Type;
//    std::cout.precision(30);
  unsigned int m = 0;
  unsigned int n = 0;
  int s = 0;
  unsigned int partition = 50;
  int normal;
  Type formulaArea(0);
  std::vector<Type> k_values (0);

  std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< Type >>>>>> parabola_table(0) ;
  std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< Type >>>>>> parabola_table_4intersection(0);

  creat_parabola_table(parabola_table, partition, m, n,s);

//   creat_parabola_table<Type>(parabola_table, partition, m,n,s);
  creat_parabola_table_4intersection(parabola_table_4intersection , partition, m, n, s);

    std::vector <Type> given_parabola(4,0);
    std:: vector <Type> intersection(0);
    std:: vector <Type> interp_point(0);
    std:: vector <std:: vector <Type>> interp_table(8, std:: vector <Type> (4));
    std:: vector <vector<Type>> badcases(0);


    unsigned int table_number;
    Type interp_area;
    Type totalerrorID(0);
    int totalSwapped(0);
    int badSwap(0);

//     =======================================================
    int grid_n= 512; // Number of grids
    Type gridArea(0);
    Type area(0);

    // Circle parameters
    Circle circle;
    circle.center = {0.5, 0.5}; // Assuming the circle is centered at (0.5, 0.5)
    circle.radius = 0.379; // Assuming the radius of the circle is 0.4
    double grid_size = 1.0 / grid_n;
    unsigned ccount = 0;
    unsigned formulaUsed =0 ;


//     std::cout << " Type of Point :  " << typeid(circle).name();

    // Loop through each grid and find the intersection points with the circle
    for (int i = 0; i < grid_n; ++i) {
        for (int j = 0; j < grid_n; ++j) {
            gridArea = 0 ;
            bool isswap = 0;
            // Define the corners of the current grid
            Point grid_leftdown = {i * grid_size, j * grid_size};
//             Point grid_lefttop = {i * grid_size, (j+1) * grid_size};
//             Point grid_rightdown = {(i+1) * grid_size, j * grid_size};
            Point grid_righttop = {(i + 1) * grid_size, (j + 1) * grid_size};

            std::vector<int> ltrbNumber{0,0,0,0}; //this vector tells us which axis it intersect and how many times. (left,top,right,bottom)

            //all it does find whether a grid top half or bottom half of the circle


            std::vector<Point> grid_intersections = findGridIntersectionPoints(circle, grid_leftdown, grid_righttop, ltrbNumber , table_number);

            std::vector<Point> mapped_intersections = mapIntersectionPoints(grid_intersections, grid_leftdown, grid_righttop, table_number, interp_point,isswap);
            if(grid_intersections.size() == 3){
              if(isswap){ // 1st and 2nd quadrant change the normal
                normal= ((grid_leftdown.x+grid_righttop.x)*0.5 > circle.center.x )? 1 : 0 ;
                cout << "we used isswap"<<endl;
                totalSwapped++;
              }
              else  normal= ((grid_leftdown.y+grid_righttop.y)*0.5 > circle.center.y )? 1 : 0 ;
//               cout << "normal = "<< normal<<endl;

              interpolationTable(parabola_table, interp_point, interp_table, partition, table_number, normal);
            }

//               This gridArea is the area in referance space (in mapped space) then scaled down.
            gridArea = findGridArea<Type>(circle, grid_leftdown, grid_righttop, grid_intersections,interp_table, interp_point, partition,grid_size);

             // direct integral
            PointT <Type> p1,p2,p3;
            Parabola <Type> parabola ;
            Type c(1) , area_b4(0) , area_aftr(0);

            {
              if (mapped_intersections.size() >= 3){ // find the parabola equation from first three points.
                ccount ++;

                p1 = { static_cast<Type>(mapped_intersections[0].x), static_cast<Type>(mapped_intersections[0].y) };
                p2 = { static_cast<Type>(mapped_intersections[1].x), static_cast<Type>(mapped_intersections[1].y) };
                p3 = { static_cast<Type>(mapped_intersections[2].x), static_cast<Type>(mapped_intersections[2].y) };

                Type det = p1.x * p1.x * (p2.x - p3.x) -p1.x* (p2.x*p2.x - p3.x*p3.x)+ p2.x*p3.x*(p2.x - p3.x) ;

                parabola = get_parabola_equation(p1, p2, p3, det);

                Type k, b, d, a(0);

                std::vector <Type> pol1(3, 0);
                std::vector <Type> pol2(3, 0);
                Type A1 = 0, A2 = 0, A3 = 0;
                Type B1 = 0, B2 = 0, B3 = 0;

                if(normal == 0){
                  c=1;
                  pol1[0] = parabola.k ; pol1[1] = a + parabola.b; pol1[2] = c + parabola.d;  pol2[0] = parabola.k; pol2[1] = parabola.b; pol2[2] = parabola.d;
                }
                else if(normal ==1){
                  c=-1;
                  pol1[0] = -parabola.k ; pol1[1] = -a - parabola.b; pol1[2] = c - parabola.d;  pol2[0] = -parabola.k; pol2[1] = -parabola.b; pol2[2] = -parabola.d;
                }

                k_values.resize(k_values.size()+1);
                k_values[k_values.size()-1] = parabola.k;

//                 cout << " prabola equation from 3 points = " <<pol2[0]<<"x^2 + " << pol2[1]<<" x + "<< pol2[2] << " + " << c << "y "<< endl;


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

                area_b4 = (A1 + A2 + A3);

                area_aftr = (A1 + A2 + A3)*grid_size*grid_size ;

                formulaArea += area_aftr;



                totalerrorID += (gridArea - area_aftr)*grid_n*grid_n;
                if(fabs((gridArea - area_aftr)* grid_n* grid_n) > 0.0005){
                  badcases.resize(badcases.size()+1);
                  badcases[badcases.size()-1].resize(13);
                  badcases[badcases.size()-1] = {grid_size*i, grid_size*j, (gridArea - area_aftr)*grid_n*grid_n,p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, (0.5* fabs(p1.x*(p2.y-p3.y) + p2.x*(p3.y - p1.y) + p3.x*(p1.y-p2.y))),-parabola.k, -parabola.b, -parabola.d};

                  if(isswap){
                    badSwap++;
                  }

//                  if((0.5* fabs(p1.x*(p2.y-p3.y) + p2.x*(p3.y - p1.y) + p3.x*(p1.y-p2.y))) < (0.3/grid_n) /*&& mapped_intersections.size() == 3*/){
//                   gridArea = area1 ;
//                   formulaUsed++;
//                   }

                }

// TODO Using direct formula for bad cases

//                 if((0.5* fabs(p1.x*(p2.y-p3.y) + p2.x*(p3.y - p1.y) + p3.x*(p1.y-p2.y))) < (0.12/grid_n) && mapped_intersections.size() == 3 && partition<grid_n){
//                   gridArea = area1 ;
//                   formulaUsed++;
//                 }

              }
              else formulaArea += gridArea;
            }




            // Output the intersection points of the circle with the current grid



//             for (const Point& intersection : intersections) {
//                 std::cout << "(" << intersection.x << ", " << intersection.y << ") ";
//             }



            area += gridArea;

            if(gridArea!=0 && gridArea/grid_size/grid_size != 1){
              cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<endl;
              cout << "Grid (" << i*grid_size << ", " << j*grid_size << "): ";
              cout <<" table number aa :" << table_number << " grid size =" << grid_size <<endl;
              for (int ap = 0; ap < grid_intersections.size();ap++){
                std::cout << " (" << grid_intersections[ap].x << ", " << grid_intersections[ap].y << ") ";
//                 std::cout << "(" << mapped_intersections[ap].x << "__ " << mapped_intersections[ap].y << ") ";
              }
              cout<< " formula grid area before scaling = " << area_b4 <<endl;
              cout<< " formula grid area after scaling = " << area_aftr <<endl;
              cout<< " formula area total = " << formulaArea <<endl;
              cout << " difference between direct and interpolation : " << (gridArea - area_aftr)*grid_n*grid_n <<endl;
              cout << "LTRB " << ltrbNumber[0] <<ltrbNumber[1]<< ltrbNumber[2]<< ltrbNumber[3]<< endl;
              cout << " 3 points :  " << "(" << p1.x <<" , "<< p1.y<<")" << "(" << p2.x <<" , "<< p2.y<<")"<< "(" << p3.x <<" , "<< p3.y<<")" << endl;
              cout << " prabola equation from 3 points = " <<parabola.k<<"x^2 + " << parabola.b<<" x + "<< parabola.d << " + " << c << "y "<< endl;
              std::cout << " Total area = " << area<<endl;
              cout<< " --------------------------------------------"<<endl;
            }
        }
    }


    cout<< endl;
    for(int xm =0 ; xm < badcases.size(); xm++){
      cout << "Grid=(" << badcases[xm][0] <<"," << badcases[xm][1] << "), area =" << badcases[xm][2] << "; points =("<< badcases[xm][3] << ","<< badcases[xm][4]<< "),(" << badcases[xm][5] << ","<< badcases[xm][6]<< "),("<< badcases[xm][7] << ","<< badcases[xm][8]<< ") triangle area= " << badcases[xm][9] /*<<" parabola= " << badcases[xm][10]<<"x^2+"<<badcases[xm][11]<<"x+"<<badcases[xm][12]*/  <<  endl;
    }
      for(int xm =0 ; xm < badcases.size(); xm++){
      cout << badcases[xm][10]<<"x^2+"<<badcases[xm][11]<<"x+"<<badcases[xm][12] <<endl;
    }

    //find the grids..
/*
      for(int xm =0 ; xm < badcases.size(); xm++){
      cout << "(x-"<<badcases[xm][0]<<")^2 + (y-"<<badcases[xm][1]<<")^2 ="<<"0.004^2" <<endl;
    }*/



    cout<< "The circle went through " << ccount << " grid."  << endl;
    cout<< "Bad cases              =" << badcases.size()  << endl;
    cout<< "Formula used           =" << formulaUsed  << endl;
    cout<< "taotal swapped           =" << totalSwapped  << endl;
    cout<< "bad swapped           =" << badSwap  << endl;

/*    cout <<endl<< " total error from direct formula in referance space: " << totalerrorID <<endl;
        cout <<endl<< " total error from direct formula in physical space: " << totalerrorID*grid_size*grid_size <<endl*/;


    std::cout << "Total area = " << area  << "formulaArea = " <<formulaArea << " original area = " << circle.radius * circle.radius * M_PI << " interpolation error = " << fabs(area - circle.radius * circle.radius * M_PI ) << " formula error = " <<fabs(formulaArea - circle.radius * circle.radius * M_PI ) << endl;

cout<< " =============================================" <<endl;
        cout<< " k values :" <<endl;
        Type summ(0);
// for (int kk=0; kk<k_values.size();kk++){
//
//  cout << k_values[kk] << " " ;
//  summ+=k_values[kk];
//
// }

// cout<< " average k  "
    return 0;
}



























