#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <algorithm>    // std::sort
#include <ctime>
#include <cstdlib>
#include <climits>

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
Parabola <Type>  get_parabola_equation( Point <Type> p1, Point <Type> p2, Point <Type> p3/*, Type &k, Type &b, Type &d*/) {
    Type x1 = p1.x, x2 = p2.x, x3 = p3.x;
    Type y1 = p1.y, y2 = p2.y, y3 = p3.y;
    Type det = x1 * x1 * (x2 - x3) -x1* (x2*x2 - x3*x3)+ x2*x3*(x2 - x3) ;
//     Type denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
    Type k = (y1 * (x2 - x3) + y2 * (x3 - x1) + y3 * (x1 - x2)) / det;
    Type b = (y1 * (x3*x3 - x2*x2) + y2 * (x1*x1 - x3*x3)+ y3 * ((x2*x2 - x1*x1))) / det;
    Type d = (y1 * x2 * x3 * (x2 - x3) + y2 * x3 * x1 * (x3 - x1) + y3 * x1 * x2 * (x1 - x2)) / det;
    return {k, b, d};
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

int main() {

  unsigned int m = 0;
  unsigned int n = 0;
  int s = 0;
  unsigned int count = 0;

  typedef cpp_bin_float_oct Type;
//     typedef double Type;

//  std::cout.precision(20);

  Type c, area ;
  Type a(0);

  std::vector <Type> pol1(3, 0);
  std::vector <Type> pol2(3, 0);
  std::vector< std::vector< Type >> parabola_table(0) ;
  Point <Type> p1, p2, p3 ;
  cout << " Number | (x1,y1) | (x2,y2) | (x3,y3) |  k  |  b  |  d  |  c  | Area |" <<endl;
  clock_t t = clock();

  double partition = 2;
  double del_x = 1/partition;
  for (int table = 1; table <=8; table++){
      cout << "Table " << table << endl;
    for (double i1=0.;i1<=1.;i1+= del_x){
      for (double i2=0.;i2<=1.;i2+=del_x){
        for (double i3=0.;i3<=1.;i3+=del_x){
            switch (table) {
                case 1:
                    p1 = {0, i1};
                    p2 = {i2, 1};
                    break;
                case 2:
                    p1 = {0, i1};
                    p2 = {1,i2};
                    break;
                case 3:
                    p1 = {0, i1};
                    p2 = {i2, 0};
                    break;
                case 4:
                    p1 = {i1,1};
                    p2 = {i2, 1};
                    if(i2 >= i1){p1 = {0, 0};p2 = {0, 0};}
                    break;
                case 5:
                    p1 = {i1,1};
                    p2 = {1, i2};
                    break;
                case 6:
                    p1 = {i1,1};
                    p2 = {i2, 0};
                    break;
                case 7:
                    p1 = {1, i1};
                    p2 = {i2, 0};
                    break;
                case 8:
                    p1 = {i1, 0};
                    p2 = {i2, 0};
                    if(i2 >= i1){p1 = {0, 0};p2 = {0, 0};}
                    break;
            }
//           Point p1 = {0, i1};
//           Point p2 = {i2, 1};
           p3 = {0.5 * (p1.x + p2.x) , i3 };
           c=1 ;

          Type det = p1.x * p1.x * (p2.x - p3.x) -p1.x* (p2.x*p2.x - p3.x*p3.x)+ p2.x*p3.x*(p2.x - p3.x) ;  // only sort out the points parallel to y axis
          cout << det << endl;
          if (det !=0){
            for(int normal=0; normal <=1; normal++ ){
                Type A1 (0), A2 (0), A3 (0);
                Parabola <Type> parabola = get_parabola_equation(p1, p2, p3);
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
                cout << " \n" << count << ". (" << p1.x << ", " << p1.y << ") ," << "(" << p2.x << ", " << p2.y << ") ," << "(" << p3.x << ", " << p3.y << ")  : "  <<  parabola.k << ", " << parabola.b << ", " << parabola.d << ", " << c << ", " << area<<  endl;

                parabola_table.resize(parabola_table.size() + 1);
                parabola_table[count].resize(12);

                parabola_table[count][0] = count;
                parabola_table[count][1] = p1.x;
                parabola_table[count][2] = p1.y;
                parabola_table[count][3] = p2.x;
                parabola_table[count][4] = p2.y;
                parabola_table[count][5] = p3.x;
                parabola_table[count][6] = p3.y;

                parabola_table[count][7] = parabola.k;
                parabola_table[count][8] = parabola.b;
                parabola_table[count][9] = parabola.d;
                parabola_table[count][10] = c;
                parabola_table[count][11] = area;


//                 for (int nm=0; nm<=11; nm++){
//                 cout << parabola_table[count][nm] << " " ;
//                 }
                count ++ ;
            }


          }

        }
      }
    }
  }

  t = clock() - t;
  std::cout << "Time taken " << (Type)(t) / CLOCKS_PER_SEC << std::endl;

    return 0;
}








    //           cout << "The equation of the parabola is: " << parabola.a << "x^2 + " << parabola.b << "x + " << parabola.c << endl;
    //         cout << "The input points are: " << endl;
    //     cout << " (" << p1.x << ", " << p1.y << ") ," << "(" << p2.x << ", " << p2.y << ") ," << "(" << p3.x << ", " << p3.y << ")" << endl;


// void get_parabola_equation(float x1, float y1, float x2, float y2, float x3, float y3, float& a, float& b, float& c) {
//     // Calculate the determinant of the 3x3 matrix
//     float detA = x1 * (y2 - y3) - x2 * (y1 - y3) + x3 * (y1 - y2);
//     if (detA == 0) {
//         cout << "Error: Points are collinear!" << endl;
//         return;
//     }
//
//     // Calculate the coefficients a, b, and c of the parabola equation
//     float A = y1 * (x2 - x3) - y2 * (x1 - x3) + y3 * (x1 - x2);
//     float B = pow(x1, 2) * (y2 - y3) - pow(x2, 2) * (y1 - y3) + pow(x3, 2) * (y1 - y2);
//     float C = pow(x1, 2) * (y3 - y2) + pow(x2, 2) * (y1 - y3) + pow(x3, 2) * (y2 - y1);
//     a = A / detA;
//     b = B / detA;
//     c = C / detA;
// }
//
// int main() {
//     // Example usage of the function
//     float x1 = -0.5, y1 = 3;
//     float x2 = 0, y2 = 4;
//     float x3 = -1, y3 = 3;
//     float a, b, c;
//     get_parabola_equation(x1, y1, x2, y2, x3, y3, a, b, c);
//     cout << "The equation of the parabola passing through (" << x1 << ", " << y1 << "), (" << x2 << ", " << y2 << "), and (" << x3 << ", " << y3 << ") is:" << endl;
//     cout << "y = " << a << "x^2 + " << b << "x + " << c << endl;
//     return 0;
// }
//

/*
int main() {
clock_t t = clock();
std::srand(std::time(NULL));
for(unsigned i=0;i<100000;i++){



}




   t = clock() - t;
   std::cout << "\nTime taken for predetermined cases: " <<(double)(t)/ CLOCKS_PER_SEC << std::endl;
   return 1;
}*/



/*
    Point p1 = {0, 0.5};
    Point p2 = {0.5, 0.75};
    Point p3 = {1, 0.5};*/
//     Point p1 = { -0.5, 3};
//     Point p2 = {0, 4};
//     Point p3 = {-1, 3};
