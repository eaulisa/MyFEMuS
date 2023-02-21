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



// double powr(const double &x, const int &y){
//   double x_p_y = 1.;
//   if (y > 0){
//     for(unsigned i=0; i < y; i++){
//       x_p_y *= x ;
//     }
//   }
//
//   else if(y<0) {
//     for(unsigned i=0; i < -y; i++){
//       x_p_y *= x ;
//     }
//     x_p_y = 1./x_p_y ;
//   }
//
//   else {
//   x_p_y = 1. ;
//   }
//   return x_p_y ;
// }


template <class Type>
void GetIntervalall(const std::vector <Type> &a1, const std::vector <Type> &a2, std::vector< std::pair<Type, Type> > &I1, std::vector< std::pair<Type, Type> > &I2, std::vector<std::pair<Type, Type>> &I3);

template <class Type>
void random_polynomial(std::vector <Type> &a1, std::vector <Type> &a2) {
  a1[0] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
  a1[1] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
  a1[2] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
  a2[0] = a1[0] ;
  a2[1] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
  a2[2] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
//             std::cout <<"\n ** k = "<<a2[0] << "; b = " << a2[1] << "; d = " << a2[2] << "; a = " << a1[1] - a2[1]<< "; c = " << a1[2] - a2[2] << ";" << std::endl;
}

template <class Type>
Type integral_A2(const unsigned &m, const unsigned &n, const int &s, const Type &a, const Type &c, const std::vector <Type> &pol1, const std::vector< std::pair<Type, Type> > &I2) {

  Type A2 = 0;
  if(a == 0) {
    for(unsigned i = 0; i < I2.size(); i++) {
      int pMax = s + n + 1;
// #1
      for(int r = 0; r <= pMax; r++) {
        Type sum = 0;
        Type r_pm_p1 = r + m + 1;
        for(int p = 0; p <= r / 2; p++) {
          sum += (pow(pol1[0], p) * pow(pol1[1], r - 2 * p) * pow(pol1[2], pMax + p - r)) / (factorial<Type>(p) * factorial<Type>(r - 2 * p) * factorial<Type>(pMax + p - r));
        }
        A2 += sum  * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1;
      }
// #2
      for(int r = 0; r < pMax; r++) {
        Type sum = 0;
        Type r_pm_p1 = 2 * pMax - r + m + 1;
        for(int p = 0; p <= r / 2; p++) {
          sum += (pow(pol1[2], p) * pow(pol1[1], r - 2 * p) * pow(pol1[0], pMax + p - r)) / (factorial<Type>(p) * factorial<Type>(r - 2 * p) * factorial<Type>(pMax + p - r));
        }
        A2 += sum  * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1;
      }
    }
    A2 *= pow(-1, n + 1) * factorial<Type>(n) / pow(c, n + 1);

    return A2;
  }
  else {
    std::vector <Type> k(3);
    std::cout.precision(20);
    // std::cout << "AAA "<<pol1[0]<<" "<<pol1[1]<<std::endl;

    k[0] = pol1[0] / (a * a);
    k[1] = pol1[1] / a;
    k[2] = k[0] * c * c - k[1] * c + pol1[2];
    k[1] -= 2 * c * k[0];

    std::vector <Type> A(s + n + 2, 0);
    std::vector <Type> B(s + n + 2, 0);

    unsigned qMax = s + n + 1;

    //   std::cout << " ankor1 " << "k0 = " << k[0] << " k1 = " << k[1] << " k2 = " << k[2] << " qMax = " << qMax << std::endl ;

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
        B[q] = A[q] * (pow(k[1], q) * pow(k[0], s + n + 1 - q)) / (factorial<Type>(q) * factorial<Type>(s + n + 1 - q));
        A[q] *= (pow(k[1], q) * pow(k[2], s + n + 1 - q)) / (factorial<Type>(q) * factorial<Type>(s + n + 1 - q));
        //   std::cout<<"A["<<q<<"] = " << A[q] <<"  B[] ="<< B[q] << std::endl;
//         std::cout << "A[" << q << "] = " << A[q] << "  B[] =" << B[q] << std::endl;
      }
    }

    else {
      Type kterms = (k[0] * k[2]);
      for(int q = 0; q <= qMax; q++) {
        Type term = 1;
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
        B[q] =  A[q] * pow(k[1], q) * pow(k[0], s + n + 1 - q) / (factorial<Type>(q) * factorial<Type>(s + n + 1 - q));
        A[q] *= pow(k[1], q) * pow(k[2], s + n + 1 - q) / (factorial<Type>(q) * factorial<Type>(s + n + 1 - q));
//         std::cout << "A[" << q << "] = " << A[q] << "  B[] =" << B[q] << std::endl;
      }
    }

    //integration starts from here.....
    if(m >= qMax) {
      for(unsigned i = 0; i < I2.size(); i++)  {
        Type u1 = a * I2[i].first + c;
        Type u2 = a * I2[i].second + c;
//         std::cout<< " u1= "<< u1 << std::endl;
//         std::cout<< " u2= "<< u2 << std::endl;

        if(u1 == 0 || u2 == 0) {
          Type c_0 = (a * pol1[1] - pol1[0] * c) / (a * a);
          int pMax = s + n + 1 ;
          // #1
          for(int r = 0; r <= s; r++) {
            Type sum = 0;
            Type r_pm_p1 = r + m + 1;
            for(int p = 0; p <= r; p++) {
              sum += (pow(a, r - 2 * p) * pow(pol1[0], p) * pow(c, s + p - r) * pow(c_0, pMax - p)) / (factorial<Type>(p) * factorial<Type>(r - p) * factorial<Type>(s - r + p) * factorial<Type>(pMax - p));
            }
            A2 += sum  * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1;
          }
          // #2
          for(int r = s + 1; r <= pMax; r++) {
            Type sum = 0;
            Type r_pm_p1 = r + m + 1;
            for(int p = r - s; p <= r; p++) {
              sum += (pow(a, r - 2 * p) * pow(pol1[0], p) * pow(c, s + p - r) * pow(c_0, pMax - p)) / (factorial<Type>(p) * factorial<Type>(r - p) * factorial<Type>(s - r + p) * factorial<Type>(pMax - p));
            }
            A2 += sum  * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1;
          }
          // #3
          for(int r = pMax + 1; r <= pMax + s; r++) {
            Type sum = 0;
            Type r_pm_p1 = r + m + 1;
            for(int p = r - s; p <= pMax; p++) {
              sum += (pow(a, r - 2 * p) * pow(pol1[0], p) * pow(c, s + p - r) * pow(c_0, pMax - p)) / (factorial<Type>(p) * factorial<Type>(r - p) * factorial<Type>(s - r + p) * factorial<Type>(pMax - p));
            }
            A2 += sum  * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1;
          }
          A2 *= pow(-1, n + 1) * factorial<Type>(n) * factorial<Type>(s);
        }
        else {
          //        std::cout<< " I2= "<< I2[i].first << ", "<< I2[i].second << std::endl;
          //        std::cout<< " u1= "<< u1 << std::endl;
          //        std::cout<< " u2= "<< u2 << std::endl;
          // 1
          for(unsigned r = 0; r <= qMax; r++) {
            //         std::cout<< " r= "<< r << std::endl;
            Type sum = 0;
            for(unsigned q = 0; q <= r; q++) {
              //           std::cout<< " q= "<< q << std::endl;
              //           std::cout<< " factorial<Type>(m - r + q)= "<< factorial<Type>(m - r + q) << std::endl;
              //           std::cout<< " factorial<Type>(r - q)= "<< factorial<Type>(r - q) << std::endl;
              sum += A[q] * pow(-c, m - r + q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
//               std::cout<< " sum= "<< sum << std::endl;
            }
            int r_m_n = r - n;
            //         std::cout<< " r_m_n = "<< r_m_n << std::endl;
            A2 += (r != n) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2 / u1));
            //         std::cout<< " pow(u2, r - n)= "<< pow(u2, r_m_n) << std::endl;
            //         std::cout<< " pow(u1, r - n)= "<< pow(u1, r_m_n) << std::endl;
            //         std::cout<< " log(u2)= "<< log(u2) << std::endl;
            //         std::cout<< " log(u1)= "<< log(u1) << std::endl;
            //         std::cout<< " (log(u2/u1)= "<< (log(u2/u1)) << std::endl;
            //
            //         std::cout<< " r-n = "<< r-(n) << std::endl;
          }
//           std::cout << "1. A2= " << A2 << std::endl;

          // 2
          for(unsigned r = qMax + 1; r <= m; r++) {
            Type sum = 0;
            for(unsigned q = 0; q <= qMax; q++) {
              sum += A[q] * pow(-c, m - r + q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = static_cast<int>(r) - static_cast<int>(n);
            A2 += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
          }
//           std::cout << "2. A2= " << A2 << std::endl;

          // 3
          for(unsigned r = m + 1; r <= qMax + m; r++) {
            Type sum = 0;
            for(unsigned q = r - m; q <= qMax; q++) {
              sum += A[q] * pow(-c, m - r + q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = r - n;
            A2 += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
          }
//           std::cout << "3. A2= " << A2 << std::endl;

          // 4

          for(unsigned r = 0; r < qMax; r++) {
            Type sum = 0;
            for(unsigned q = 0; q <= r; q++) {
              sum += B[q] * pow(-c, r - q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            A2 += sum  * (pow(u2, qMax + s + m - r + 1) - pow(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
          }
//           std::cout << "4. A2= " << A2 << std::endl;

          // 5
          for(unsigned r = qMax; r <= m; r++) {
            Type sum = 0;
            for(unsigned q = 0; q < qMax; q++) {
              sum += B[q] * pow(-c, r - q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            A2 += sum  * (pow(u2, qMax + s + m - r + 1) - pow(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
          }
//           std::cout << "5. A2= " << A2 << std::endl;

          // 6
          for(unsigned r = m + 1; r < qMax + m; r++) {
            Type sum = 0;
            for(unsigned q = r - m; q < qMax; q++) {
              sum += B[q] * pow(-c, r - q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            A2 += sum  * (pow(u2, qMax + s + m - r + 1) - pow(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
          }
//           std::cout << "6. A2= " << A2 << std::endl;

          //total
          A2 *= pow(-1, n + 1) * factorial<Type>(n) * factorial<Type>(m) / pow(a, m + 1);
//           std::cout<< "final. A2= "<< A2 << std::endl;
        }
      }
      return A2;
    }

    else {

      for(unsigned i = 0; i < I2.size(); i++)  {
        Type u1 = a * I2[i].first + c;
        Type u2 = a * I2[i].second + c;
        //       std::cout<< " u1= "<< u1 << std::endl;
        //       std::cout<< " u2= "<< u2 << std::endl;
        if(u1 == 0 || u2 == 0) {
          Type c_0 = (a * pol1[1] - pol1[0] * c) / (a * a);
          int pMax = s + n + 1 ;
          // #1
          for(int r = 0; r <= s; r++) {
            Type sum = 0;
            Type r_pm_p1 = r + m + 1;
            for(int p = 0; p <= r; p++) {
              sum += (pow(a, r - 2 * p) * pow(pol1[0], p) * pow(c, s + p - r) * pow(c_0, pMax - p)) / (factorial<Type>(p) * factorial<Type>(r - p) * factorial<Type>(s - r + p) * factorial<Type>(pMax - p));
              //           std::cout << "1sum= " << sum << std::endl;
            }
            A2 += sum  * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1;
            //         std::cout << "11. A2= " << A2 << std::endl;

          }
          // #2
          for(int r = s + 1; r <= pMax; r++) {
            Type sum = 0;
            Type r_pm_p1 = r + m + 1;
            for(int p = r - s; p <= r; p++) {
              sum += (pow(a, r - 2 * p) * pow(pol1[0], p) * pow(c, s + p - r) * pow(c_0, pMax - p)) / (factorial<Type>(p) * factorial<Type>(r - p) * factorial<Type>(s - r + p) * factorial<Type>(pMax - p));
            }
            A2 += sum  * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1;
            //         std::cout << "22. A2= " << A2 << std::endl;
          }
          // #3
          for(int r = pMax + 1; r <= pMax + s; r++) {
            Type sum = 0;
            Type r_pm_p1 = r + m + 1;
            for(int p = r - s; p <= pMax; p++) {
              sum += (pow(a, r - 2 * p) * pow(pol1[0], p) * pow(c, s + p - r) * pow(c_0, pMax - p)) / (factorial<Type>(p) * factorial<Type>(r - p) * factorial<Type>(s - r + p) * factorial<Type>(pMax - p));
            }
            A2 += sum  * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1;
            //         std::cout << "33. A2= " << A2 << std::endl;
          }
          A2 *= pow(-1, n + 1) * factorial<Type>(n) * factorial<Type>(s);
        }
        else {
          // 1
          for(unsigned r = 0; r <= m; r++) {
            Type sum(0);
            for(unsigned q = 0; q <= r; q++) {
              //             std::cout<< " q= "<< q << std::endl;
              //             std::cout<< " factorial<Type>(m - r + q)= "<< factorial<Type>(m - r + q) << std::endl;
              //             std::cout<< " factorial<Type>(r - q)= "<< factorial<Type>(r - q) << std::endl;
              sum += A[q] * pow(-c, m - r + q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
              //             std::cout<< " sum= "<< sum << std::endl;
            }
            int r_m_n = r - n;
            //           std::cout<< " r_m_n = "<< r_m_n << std::endl;
            A2 += (r != n) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2 / u1));
            //           std::cout<< " pow(u2, r - n)= "<< pow(u2, r_m_n) << std::endl;
            //           std::cout<< " pow(u1, r - n)= "<< pow(u1, r_m_n) << std::endl;
            //           std::cout<< " log(u2)= "<< log(u2) << std::endl;
            //           std::cout<< " log(u1)= "<< log(u1) << std::endl;
            //           std::cout<< " (log(u2/u1)= "<< (log(u2/u1)) << std::endl;
            //           std::cout<< " r-n = "<< r-(n) << std::endl;
          }
//           std::cout << "1. A2= " << A2 << std::endl;

          // 2
          for(unsigned r = m + 1; r <= qMax; r++) {
            Type sum(0);
            for(unsigned q = r - m; q <= r; q++) {
              sum += A[q] * pow(-c, m - r + q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = r - n;
            A2 += (r != n) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2 / u1));
          }
//           std::cout << "2. A2= " << A2 << std::endl;

          // 3
          for(unsigned r = qMax + 1; r <= qMax + m; r++) {
            Type sum(0);
            for(unsigned q = r - m; q <= qMax; q++) {
              sum += A[q] * pow(-c, m - r + q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = r - n;
            A2 += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
          }
//           std::cout << "3. A2= " << A2 << std::endl;

          // 4

          for(unsigned r = 0; r <= m; r++) {
            Type sum(0);
            for(unsigned q = 0; q <= r; q++) {
              sum += B[q] * pow(-c, r - q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            A2 += sum  * (pow(u2, qMax + s + m - r + 1) - pow(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
            //std::cout << pow(u2, qMax + s + m - r + 1) <<" "<< pow(u1, qMax + s + m - r + 1) <<std::endl;
          }
//           std::cout << "4. A2= " << A2 << std::endl;

          // 5
          for(unsigned r = m + 1; r < qMax; r++) {
            Type sum(0);
            for(unsigned q = r - m; q <= r; q++) {
              sum += B[q] * pow(-c, r - q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            A2 += sum  * (pow(u2, qMax + s + m - r + 1) - pow(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
          }
//           std::cout << "5. A2= " << A2 << std::endl;

          // 6
          for(unsigned r = qMax; r < qMax + m; r++) {
            Type sum(0);
            for(unsigned q = r - m; q < qMax; q++) {
              sum += B[q] * pow(-c, r - q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            A2 += sum  * (pow(u2, qMax + s + m - r + 1) - pow(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
          }
//           std::cout << "6. A2= " << A2 << std::endl;

          //total
          A2 *= pow(-1, n + 1) * factorial<Type>(n) * factorial<Type>(m) / pow(a, m + 1); // TODO this sign should be checked
// //           std::cout << "final. A2= " << A2 << std::endl;
        }
      }
      return A2;
    }

  }
}

template <class Type>
Type easy_integral_A2(const unsigned &m, const unsigned &n, const int &s, const Type &a, const Type &c, const std::vector <Type> &pol1, const std::vector< std::pair<Type, Type> > &I2) {

  Type A2 = 0;
  if(a == 0) {
    for(unsigned i = 0; i < I2.size(); i++) {
      int pMax = s + n + 1;
// #1
      for(int r = 0; r <= pMax; r++) {
        Type sum = 0;
        Type r_pm_p1 = r + m + 1;
        for(int p = 0; p <= r / 2; p++) {
          sum += (pow(pol1[0], p) * pow(pol1[1], r - 2 * p) * pow(pol1[2], pMax + p - r)) / (factorial<Type>(p) * factorial<Type>(r - 2 * p) * factorial<Type>(pMax + p - r));
        }
        A2 += sum  * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1;
      }
// #2
      for(int r = 0; r < pMax; r++) {
        Type sum = 0;
        Type r_pm_p1 = 2 * pMax - r + m + 1;
        for(int p = 0; p <= r / 2; p++) {
          sum += (pow(pol1[2], p) * pow(pol1[1], r - 2 * p) * pow(pol1[0], pMax + p - r)) / (factorial<Type>(p) * factorial<Type>(r - 2 * p) * factorial<Type>(pMax + p - r));
        }
        A2 += sum  * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1;
      }
    }
    A2 *= pow(-1, n + 1) * factorial<Type>(n) / pow(c, n + 1);

    return A2;
  }
  else {
    std::vector <Type> k(3);
    std::cout.precision(20);
    // std::cout << "AAA "<<pol1[0]<<" "<<pol1[1]<<std::endl;

    k[0] = pol1[0] / (a * a);
    k[1] = pol1[1] / a;
    k[2] = k[0] * c * c - k[1] * c + pol1[2];
    k[1] -= 2 * c * k[0];

    std::vector <Type> A(s + n + 2, 0);
    std::vector <Type> B(s + n + 2, 0);

    unsigned qMax = s + n + 1;

    //   std::cout << " ankor1 " << "k0 = " << k[0] << " k1 = " << k[1] << " k2 = " << k[2] << " qMax = " << qMax << std::endl ;
// pre-evalate A[q] and B[q].
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
        B[q] = A[q] * (pow(k[1], q) * pow(k[0], s + n + 1 - q)) / (factorial<Type>(q) * factorial<Type>(s + n + 1 - q));
        A[q] *= (pow(k[1], q) * pow(k[2], s + n + 1 - q)) / (factorial<Type>(q) * factorial<Type>(s + n + 1 - q));
        //   std::cout<<"A["<<q<<"] = " << A[q] <<"  B[] ="<< B[q] << std::endl;
//         std::cout << "A[" << q << "] = " << A[q] << "  B[] =" << B[q] << std::endl;
      }
    }
    else {
      Type kterms = (k[0] * k[2]);
      for(int q = 0; q <= qMax; q++) {
        Type term = 1;
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
        B[q] =  A[q] * pow(k[1], q) * pow(k[0], s + n + 1 - q) / (factorial<Type>(q) * factorial<Type>(s + n + 1 - q));
        A[q] *= pow(k[1], q) * pow(k[2], s + n + 1 - q) / (factorial<Type>(q) * factorial<Type>(s + n + 1 - q));
//         std::cout << "A[" << q << "] = " << A[q] << "  B[] =" << B[q] << std::endl;
      }
    }

    //integration starts from here.....
      for(unsigned i = 0; i < I2.size(); i++)  {
        Type u1 = a * I2[i].first + c;
        Type u2 = a * I2[i].second + c;
        //       std::cout<< " u1= "<< u1 << std::endl;
        //       std::cout<< " u2= "<< u2 << std::endl;
        if(u1 == 0 || u2 == 0) {
          Type c_0 = (a * pol1[1] - pol1[0] * c) / (a * a);
          int pMax = s + n + 1 ;
          // #1
          for(int r = 0; r <= s; r++) {
            Type sum = 0;
            Type r_pm_p1 = r + m + 1;
            for(int p = 0; p <= r; p++) {
              sum += (pow(a, r - 2 * p) * pow(pol1[0], p) * pow(c, s + p - r) * pow(c_0, pMax - p)) / (factorial<Type>(p) * factorial<Type>(r - p) * factorial<Type>(s - r + p) * factorial<Type>(pMax - p));
              //           std::cout << "1sum= " << sum << std::endl;
            }
            A2 += sum  * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1;
            //         std::cout << "11. A2= " << A2 << std::endl;

          }
          // #2
          for(int r = s + 1; r <= pMax; r++) {
            Type sum = 0;
            Type r_pm_p1 = r + m + 1;
            for(int p = r - s; p <= r; p++) {
              sum += (pow(a, r - 2 * p) * pow(pol1[0], p) * pow(c, s + p - r) * pow(c_0, pMax - p)) / (factorial<Type>(p) * factorial<Type>(r - p) * factorial<Type>(s - r + p) * factorial<Type>(pMax - p));
            }
            A2 += sum  * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1;
            //         std::cout << "22. A2= " << A2 << std::endl;
          }
          // #3
          for(int r = pMax + 1; r <= pMax + s; r++) {
            Type sum = 0;
            Type r_pm_p1 = r + m + 1;
            for(int p = r - s; p <= pMax; p++) {
              sum += (pow(a, r - 2 * p) * pow(pol1[0], p) * pow(c, s + p - r) * pow(c_0, pMax - p)) / (factorial<Type>(p) * factorial<Type>(r - p) * factorial<Type>(s - r + p) * factorial<Type>(pMax - p));
            }
            A2 += sum  * (pow(I2[i].second, r_pm_p1) - pow(I2[i].first, r_pm_p1)) / r_pm_p1;
            //         std::cout << "33. A2= " << A2 << std::endl;
          }
          A2 *= pow(-1, n + 1) * factorial<Type>(n) * factorial<Type>(s);
        }
        else {
          Type C1;
          for(unsigned p = 0; p <= m; p++) {
            Type sum1(0);
            for(unsigned q = 0; q <= qMax; q++) {
              int i = p + q - n;
              sum1 += A[q] * ((i == 0) ? log(u2 / u1) : (pow(u2, i) - pow(u1, i)) / (i));
            }
            Type sum2(0);
            for(unsigned q = 0; q < qMax; q++) {
              int i = 2 * s + n + 2 - q;
              sum2 += B[q] * (pow(u2, i) - pow(u1, i)) / (i);
            }

            C1 = (sum1 + sum2) * pow(-c, m - p) / (factorial<Type>(p) * factorial<Type>(m - p));
          }

//           std::cout << "C1 = " << C1 << std::endl;
          A2 = C1;

          //total
          A2 *= pow(-1, n + 1) * factorial<Type>(n) * factorial<Type>(m) / pow(a, m + 1); // TODO this sign should be checked
// //           std::cout << "final. A2= " << A2 << std::endl;
        }
      }
      return A2;
  }
}

template <class Type>
Type integral_A3(const unsigned &m, const unsigned &n, const int &s, const Type &a, const Type &c, const std::vector <Type> &pol1, const std::vector< std::pair<Type, Type> > &I3) {
  Type A3 = 0.;
  if(a == 0) {
    for(int i = 0; i <= s; i++) {

      for(unsigned w = 0; w < I3.size(); w++) {
        int pMax = s - i;
        // #1
        for(int r = 0; r <= pMax; r++) {
          Type sum = 0;
          Type r_pm_p1 = r + m + 1;
          for(int p = 0; p <= r / 2; p++) {
            sum += (pow(pol1[0], p) * pow(pol1[1], r - 2 * p) * pow(pol1[2], pMax + p - r)) / (factorial<Type>(p) * factorial<Type>(r - 2 * p) * factorial<Type>(pMax + p - r));
          }
          A3 += sum  * (pow(I3[w].second, r_pm_p1) - pow(I3[w].first, r_pm_p1)) / r_pm_p1;
        }
        // #2
        for(int r = 0; r < pMax; r++) {
          Type sum = 0;
          Type r_pm_p1 = 2 * pMax - r + m + 1;
          for(int p = 0; p <= r / 2; p++) {
            sum += (pow(pol1[2], p) * pow(pol1[1], r - 2 * p) * pow(pol1[0], pMax + p - r)) / (factorial<Type>(p) * factorial<Type>(r - 2 * p) * factorial<Type>(pMax + p - r));
          }
          A3 += sum  * (pow(I3[w].second, r_pm_p1) - pow(I3[w].first, r_pm_p1)) / r_pm_p1;
        }
      }
      A3 *= pow(c, i) / ((n + i + 1) * factorial<Type>(i));
    }
    return A3;
  }

  else {
    std::vector <Type> k(3);
    k[0] = pol1[0] / (a * a);
    k[1] = pol1[1] / a;
    k[2] = k[0] * c * c - k[1] * c + pol1[2];
    k[1] -= 2 * c * k[0];

    for(int i = 0; i <= s; i++) {
      std::vector <Type> A(s - i + 1, 0);  // size of all this vector changes.
      std::vector <Type> B(s - i + 1, 0);
      unsigned qMax = s - i;

      //   std::cout<< " ankor1 " << "k0 =" << k[0]<< "k1 =" << k[1]<< "k2 =" << k[2] << "kterms =" << kterms << "qMax =" << qMax << std::endl ;

      if(k[1] != 0) {
        Type kterms = (k[0] * k[2]) / (k[1] * k[1]);
        for(int q = 0; q <= qMax; q++) {
          Type term = 1;
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
          B[q] = A[q] * (pow(k[1], q) * pow(k[0], qMax - q)) / (factorial<Type>(q) * factorial<Type>(qMax - q));
          A[q] *= (pow(k[1], q) * pow(k[2], qMax - q)) / (factorial<Type>(q) * factorial<Type>(qMax - q));
          //         std::cout<<"A["<<q<<"] = " << A[q] <<"  B[] ="<< B[q] << std::endl;
          //         std::cout << "A[" << q << "] = " << A[q] << "  B[] =" << B[q] << std::endl;
        }
      }

      else {
        Type kterms = (k[0] * k[2]);

        for(int q = 0; q <= qMax; q++) {
          Type term = 1;
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

          B[q] = A[q] * pow(k[0], qMax - q) / (factorial<Type>(q) * factorial<Type>(qMax - q));
          A[q] *= pow(k[2], qMax - q) / (factorial<Type>(q) * factorial<Type>(qMax - q));

        }
      }

//       std::cout << " ankor1 " << "k0 =" << k[0] << "k1 =" << k[1] << "k2 =" << k[2]  << "qMax =" << qMax << std::endl ;

      std::vector <Type> A3_part;
      A3_part.resize(i + 1, 0);

      if(m >= qMax) {
        for(unsigned w = 0; w < I3.size(); w++) {
          Type u1 = a * I3[w].first + c;
          Type u2 = a * I3[w].second + c;
//           std::cout << " u1= " << u1 << std::endl;
//           std::cout << " u2= " << u2 << std::endl;
          // 1
          for(unsigned r = 0; r <= qMax; r++) {
            //           std::cout << " r= " << r << std::endl;
            Type sum = 0.;
            for(unsigned q = 0; q <= r; q++) {
//               std::cout << " q= " << q << std::endl;
//               std::cout << " factorial<Type>(m - r + q)= " << factorial<Type>(m - r + q) << std::endl;
//               std::cout << " factorial<Type>(r - q)= " << factorial<Type>(r - q) << std::endl;

              sum += A[q] * pow(-c, m - r + q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));

              //          std::cout<< " sum= "<< sum << std::endl;
            }
            int r_m_n = r + i + 1;      //in A2 this was the power of u after integration. I am keeping the variable.
            //         std::cout<< " r_m_n = "<< r_m_n << std::endl;
            A3_part[i] += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2 / u1));
//                     std::cout<< " pow(u2, r - n)= "<< pow(u2, r_m_n) << std::endl;
//                     std::cout<< " pow(u1, r - n)= "<< pow(u1, r_m_n) << std::endl;
//                     std::cout<< " log(u2)= "<< log(u2) << std::endl;
//                     std::cout<< " log(u1)= "<< log(u1) << std::endl;
//                     std::cout<< " (log(u2/u1)= "<< (log(u2/u1)) << std::endl;
//                     std::cout<< " r-n = "<< r-(n) << std::endl;
          }
//           std::cout << "1. A3_part= " << A3_part[i] << std::endl;

          // 2
          for(unsigned r = qMax + 1; r <= m; r++) {
            Type sum = 0.;
            for(unsigned q = 0; q <= qMax; q++) {
              sum += A[q] * pow(-c, m - r + q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = r + i + 1;
            A3_part[i] += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
          }
//           std::cout << "2. A3_part= " << A3_part[i] << std::endl;

          // 3
          for(unsigned r = m + 1; r <= qMax + m; r++) {
            Type sum = 0.;
            for(unsigned q = r - m; q <= qMax; q++) {
              sum += A[q] * pow(-c, m - r + q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = r + i + 1;
            A3_part[i] += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
          }
//           std::cout << "3. A3_part= " << A3_part[i] << std::endl;

          // 4

          for(unsigned r = 0; r < qMax; r++) {
            Type sum = 0.;
            for(unsigned q = 0; q <= r; q++) {
              sum += B[q] * pow(-c, r - q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = qMax + s + m - r + 1;
            A3_part[i] += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2 / u1));
          }
//           std::cout << "4. A3_part= " << A3_part[i] << std::endl;

          // 5
          for(unsigned r = qMax; r <= m; r++) {
            Type sum = 0.;
            for(unsigned q = 0; q < qMax; q++) {
              sum += B[q] * pow(-c, r - q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = qMax + s + m - r + 1;
            A3_part[i] += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2 / u1));
          }
//           std::cout << "5. A3_part= " << A3_part[i] << std::endl;

          // 6
          for(unsigned r = m + 1; r < qMax + m; r++) {
            Type sum = 0.;
            for(unsigned q = r - m; q < qMax; q++) {
              sum += B[q] * pow(-c, r - q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = qMax + s + m - r + 1;
            A3_part[i] += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2 / u1));
          }
//           std::cout << "6. A3_part= " << A3_part[i] << std::endl;
        }

      }

      else {
        for(unsigned w = 0; w < I3.size(); w++)  {
          Type u1 = a * I3[w].first + c;
          Type u2 = a * I3[w].second + c;
//                 std::cout<< " u1= "<< u1 << std::endl;
//                 std::cout<< " u2= "<< u2 << std::endl;
          // 1
          for(unsigned r = 0; r <= m; r++) {
            //         std::cout<< " r= "<< r << std::endl;
            Type sum = 0.;
            for(unsigned q = 0; q <= r; q++) {
//                         std::cout<< " q= "<< q << std::endl;
//                         std::cout<< " factorial<Type>(m - r + q)= "<< factorial<Type>(m - r + q) << std::endl;
//                         std::cout<< " factorial<Type>(r - q)= "<< factorial<Type>(r - q) << std::endl;
              sum += A[q] * pow(-c, m - r + q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
              //           std::cout<< " sum= "<< sum << std::endl;
            }
            int r_m_n = r + i + 1;      //in A2 this was the power of u after integration. I am keeping the variable.
            //         std::cout<< " r_m_n = "<< r_m_n << std::endl;
            A3_part[i] += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2 / u1));
//                     std::cout<< " pow(u2, r - n)= "<< pow(u2, r_m_n) << std::endl;
//                     std::cout<< " pow(u1, r - n)= "<< pow(u1, r_m_n) << std::endl;
//                     std::cout<< " log(u2)= "<< log(u2) << std::endl;
//                     std::cout<< " log(u1)= "<< log(u1) << std::endl;
//                     std::cout<< " (log(u2/u1)= "<< (log(u2/u1)) << std::endl;
//                     std::cout<< " r-n = "<< r-(n) << std::endl;
          }
//                   std::cout<< "1. A3_part= "<< A3_part[i]<< std::endl;

          // 2
          for(unsigned r = m + 1; r <= qMax; r++) {
            Type sum = 0.;
            for(unsigned q = 0; q <= m; q++) {
              sum += A[q] * pow(-c, m - r + q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = r + i + 1;
            A3_part[i] += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
          }
//                   std::cout<< "2. A3_part= "<< A3_part[i]<< std::endl;

          // 3
          for(unsigned r = qMax + 1; r <= qMax + m; r++) {
            Type sum = 0.;
            for(unsigned q = r - m; q <= qMax; q++) {
              sum += A[q] * pow(-c, m - r + q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = r + i + 1;
            A3_part[i] += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
          }
//                   std::cout<< "3. A3_part= "<< A3_part[i]<< std::endl;

          // 4

          for(unsigned r = 0; r < m; r++) {
            Type sum = 0.;
            for(unsigned q = 0; q <= r; q++) {
              sum += B[q] * pow(-c, r - q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = qMax + s + m - r + 1;
            A3_part[i] += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2 / u1));
          }
//                   std::cout<< "4. A3_part= "<< A3_part[i]<< std::endl;

          // 5
          for(unsigned r = m + 1; r < qMax; r++) {
            Type sum = 0.;
            for(unsigned q = r - m; q <= r; q++) {
              sum += B[q] * pow(-c, r - q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = qMax + s + m - r + 1;
            A3_part[i] += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2 / u1));
          }
//                   std::cout<< "5. A3_part= "<< A3_part[i]<< std::endl;

          // 6
          for(unsigned r = qMax; r < qMax + m; r++) {
            Type sum = 0.;
            for(unsigned q = r - m; q < qMax; q++) {
              sum += B[q] * pow(-c, r - q) / (factorial<Type>(m - r + q) * factorial<Type>(r - q));
            }
            int r_m_n = qMax + s + m - r + 1;
            A3_part[i] += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2 / u1));
          }
//                   std::cout<< "6. A3_part= "<< A3_part[i]<< std::endl;
        }

      }
      //total
      A3_part[i] = A3_part[i] / ((n + i + 1) * factorial<Type>(i));

      A3 += A3_part[i] ;
    }
    A3 *= factorial<Type>(m) / pow(a, m + 1);
//     std::cout<< "final. A3= "<< A3 << std::endl;
    return A3;
  }
}

int main() {
  unsigned int m = 0;
  unsigned int n = 0;
  int s = 0;

  std::cout.precision(20);

 // typedef cpp_bin_float_quad Type;
   typedef double Type;
  Type k, b, d, a, c, area1, area2;
  std::vector <Type> pol1(3, 0);
  std::vector <Type> pol2(3, 0);
  clock_t t = clock();
  std::srand((unsigned)std::time(NULL));
  int count = 0;
  for(unsigned int j = 0; j < 100000; j++) {
    Type A1 = 0, A2 = 0, A3 = 0;
    Type B1 = 0, B2 = 0, B3 = 0;
//     m = (rand() % 6) ;
//     n= (rand() % 6) ;
//     s= (rand() % 3) ;
    random_polynomial(pol1, pol2);
     a = pol1[1] - pol2[1];
     c = pol1[2] - pol2[2];

//      k = 0.85764624125214572459; b = 0.54315816543212092071; d = 1.4313877017383407342; a = -0.0012779049581279622316; c = -1.5795148879194236269;


//     std::vector<std::vector<Type>> sample{{1, -1, 0.25, 0.5, -0.1, 0,0,1}, {0, 0,0,0.5,0,1,0,0}, {1, 0,0,0.5,0,0,0,1},{0.1, 0.2,-1,2,0.1,0.251724386116496,0,0}, {1, -1, 0.25, 0, -0.25, 0,0.3333333333333333,0},{1, -1, 0.25, 0, -0.1, 0,0.2108185106778920348,0.36754446796632406214},{20,-8.4,0.6,-6.5,1.3,0.045806466059993167228,0.136242991647203604,0.5095766326720312378}, {-0.69657011083167508225, -0.69655399150054631008, 1.4832208284564414313, 1, -1.4513087502919645999,0,0.34978109786848710083,0.52319330960346455139}};
//
//     for(unsigned j=0;j < sample.size(); j++){
//       Type A1 = 0, A2 = 0, A3 = 0;
//       Type B1 = 0, B2 = 0, B3 = 0;
//       pol1[0] = sample[j][0];
//       pol1[1] = (sample[j][1]+sample[j][3]);
//       pol1[2] = (sample[j][2]+sample[j][4]);
//       pol2[0] = sample[j][0];
//       pol2[1] = sample[j][1];
//       pol2[2] = sample[j][2];
//       a = sample[j][3];
//       c = sample[j][4];
// k = -1.275057766248964608; b = 1.2482106868402151889; d = 1.7498163020935915135; a = 0.13493826991642743351; c = -1.9456911636263554133;
// k = -1.5517184844015718959; b = 0.74869137571598010084; d = -0.051370344148655711081; a = -4.0164217371740917883e-05; c = 1.5179394621019901557;
//      k = -1.4409885720540716036; b = 0.97870700386292641682; d = 1.6552636733535974756; a = -1.5952086400218350448; c = 0.040815498698882457518;
//     pol1[0] = k; pol1[1] = a + b; pol1[2] = c + d; pol2[0] = k; pol2[1] = b; pol2[2] = d;

        std::vector< std::pair <Type, Type> > I1, I2, I3;
        GetIntervalall(pol1, pol2, I1, I2, I3);

  //       std::cout<< "\nSample " << j+1 << " : " <<std::endl;
  //       std::cout <<"\nm = "<< m << "; n = "<< n << "; s = " << s << "; k = "<<pol2[0] << "; b = " << pol2[1] << "; d = " << pol2[2] << "; a = " << a << "; c = " << c << ";" << std::endl;
  //           for(unsigned i = 0; i < I1.size(); i++) {std::cout << "I1_1 = " << I1[i].first << "; I1_2 = " << I1[i].second << ";" << std::endl;}
  //           for(unsigned i = 0; i < I2.size(); i++) {std::cout << "I2_1 = " << I2[i].first << "; I2_2 = " << I2[i].second << ";" << std::endl;}
  //           for(unsigned i = 0; i < I3.size(); i++) {std::cout << "I3_1 = " << I3[i].first << "; I3_2 = " << I3[i].second << ";" << std::endl;}

          if(I1.size() > 0) {
            A1 = integral_A3(m, n, s, a, c, pol2, I1) -  easy_integral_A2(m, n, s, a, c, pol2, I1);
          }
          if(I2.size() > 0) {
            A2 = easy_integral_A2(m, n, s, a, c, pol2, I2);
          }
          if(I3.size() > 0) {
            A3 = integral_A3(m, n, s, a, c, pol2, I3);
          }

          area1 = A1+A2+A3;
          pol1[0] *=-1; pol1[1] *=-1; pol1[2] *=-1; pol2[0] *=-1; pol2[1] *=-1; pol2[2] *=-1; a *=-1; c *=-1;
          GetIntervalall(pol1, pol2, I1, I2, I3);
          if(I1.size() > 0) {B1 = integral_A3(m, n, s, a, c, pol2, I1) -  easy_integral_A2(m, n, s, a, c, pol2, I1);}
          if(I2.size() > 0) {B2 = easy_integral_A2(m, n, s, a, c, pol2, I2);}
          if(I3.size() > 0) {B3 = integral_A3(m, n, s, a, c, pol2, I3);}
          area2 = B1+B2+B3;

    typedef cpp_bin_float_oct oct;
    //typedef cpp_bin_float_quad oct;
    oct C1 = 0, C2 = 0, C3 = 0;
    oct D1 = 0, D2 = 0, D3 = 0;
    oct area3,area4;

    oct ao = static_cast <oct>(a);
    oct co = static_cast <oct>(c);

    std::vector <oct> pol1o(3);
    pol1o[0] = static_cast <oct>(pol1[0]);
    pol1o[1] = static_cast <oct>(pol1[1]);
    pol1o[2] = static_cast <oct>(pol1[2]);

    std::vector <oct> pol2o(3);
    pol2o[0] = static_cast <oct>(pol2[0]);
    pol2o[1] = static_cast <oct>(pol2[1]);
    pol2o[2] = static_cast <oct>(pol2[2]);

    std::vector< std::pair<oct, oct> > I1o, I2o, I3o;
    GetIntervalall(pol1o, pol2o, I1o, I2o, I3o);

//       std::cout <<"\nm = "<< m << "; n = "<< n << "; s = " << s << "; k = " << pol2o[0] << "; b = " << pol2o[1] << "; d = " << pol2o[2] << "; a = " << ao << "; c = " << co << ";" << std::endl;
//       for(unsigned i = 0; i < I1o.size(); i++) {std::cout << "x1 = " << I1o[i].first << "; x2 = " << I1o[i].second << ";" << std::endl;}
//       for(unsigned i = 0; i < I2o.size(); i++) {std::cout << "y1 = " << I2o[i].first << "; y2 = " << I2o[i].second << ";" << std::endl;}
//       for(unsigned i = 0; i < I3o.size(); i++) {std::cout << "z1 = " << I3o[i].first << "; z2 = " << I3o[i].second << ";" << std::endl;}
//    clock_t t = clock();

      if(I1o.size() > 0){ C1 = integral_A3(m, n, s, ao, co, pol2o, I1o) - easy_integral_A2(m, n, s, ao, co, pol2o, I1o); }
      if(I2o.size() > 0){ C2 = easy_integral_A2(m, n, s, ao, co, pol2o, I2o); }
      if(I3o.size() > 0){ C3 = integral_A3(m, n, s, ao, co, pol2o, I3o); }
//       std::cout << "oct A1= " << B1 << std::endl;
//       std::cout << "oct A2= " << B2 << std::endl;
//       std::cout << "oct A3= " << B3 << std::endl;
//       std::cout << "Area= " << A1 + A2 + A3 << std::endl;
//         t = clock() - t;
//        std::cout << "Time taken for predetermined cases: " << (Type)(t) / CLOCKS_PER_SEC << std::endl;
//        this is for the for loop pass/fail
//        if((abs(B1 - sample[j][5]) > 0.0000000001) || (abs(B2 - sample[j][6]) > 0.0000000001) || (abs(B3 - sample[j][7]) > 0.0000000001)){
//       std::cout << " Failed " << std::endl;
//             std::cout << "\nm = " << m << "; n = " << n << "; s = " << s << "; k = " << pol2[0] << "; b = " << pol2[1] << "; d = " << pol2[2] << "; a = " << a << "; c = " << c << ";" << std::endl;
//       std::cout << "\n diff" << abs(B1 - sample[j][5]) << " " << abs(B2 - sample[j][6]) << " " << abs(B3 - sample[j][6]) << " ---failed " << count + 1 << std::endl;
//       count ++ ;
//        }

      area3 = C1+C2+C3;
      pol1o[0] *=-1; pol1o[1] *=-1; pol1o[2] *=-1; pol2o[0] *=-1; pol2o[1] *=-1; pol2o[2] *=-1; ao *=-1; co *=-1;
      GetIntervalall(pol1o, pol2o, I1o, I2o, I3o);
      if(I1o.size() > 0){ D1 = integral_A3(m, n, s, ao, co, pol2o, I1o) - easy_integral_A2(m, n, s, ao, co, pol2o, I1o); }
      if(I2o.size() > 0){ D2 = easy_integral_A2(m, n, s, ao, co, pol2o, I2o); }
      if(I3o.size() > 0){ D3 = integral_A3(m, n, s, ao, co, pol2o, I3o); }
      area4 = D1+D2+D3;

      Type err = 0.0001;
        if((abs(area1+area2-1) > 0.0000000001) || (abs(area3+area4-1) > 0.00000000000001) || (abs(D1 - A1) > err) || (abs(D2 - A2) > err) || (abs(D3 - A3) > err) ){
            std::cout << "................................ Failed...................................... " << std::endl;
            std::cout << "\nm = " << m << "; n = " << n << "; s = " << s << "; k = " << pol2o[0] << "; b = " << pol2o[1] << "; d = " << pol2o[2] << "; a = " << ao << "; c = " << co << ";" << std::endl;
            for(unsigned i = 0; i < I1o.size(); i++) {std::cout << "I1_1 = " << I1o[i].first << "; I1_2 = " << I1o[i].second << ";" << std::endl;}
            for(unsigned i = 0; i < I2o.size(); i++) {std::cout << "I2_1 = " << I2o[i].first << "; I2_2 = " << I2o[i].second << ";" << std::endl;}
            for(unsigned i = 0; i < I3o.size(); i++) {std::cout << "I3_1 = " << I3o[i].first << "; I3_2 = " << I3o[i].second << ";" << std::endl;}
            std::cout << "double A1= " << A1 << "; oct A1= " << D1 << std::endl;
            std::cout << "double A2= " << A2 << "; oct A2= " << D2 << std::endl;
            std::cout << "double A3= " << A3 << "; oct A3= " << D3 << std::endl;
            std::cout << "\nm = " << m << "; n = " << n << "; s = " << s << "; k = " << pol2[0] << "; b = " << pol2[1] << "; d = " << pol2[2] << "; a = " << a << "; c = " << c << ";" << std::endl;
            for(unsigned i = 0; i < I1.size(); i++) {std::cout << "I1_1 = " << I1[i].first << "; I1_2 = " << I1[i].second << ";" << std::endl;}
            for(unsigned i = 0; i < I2.size(); i++) {std::cout << "I2_1 = " << I2[i].first << "; I2_2 = " << I2[i].second << ";" << std::endl;}
            for(unsigned i = 0; i < I3.size(); i++) {std::cout << "I3_1 = " << I3[i].first << "; I3_2 = " << I3[i].second << ";" << std::endl;}
            std::cout << "double -A1= " << B1 << "; oct -A1= " << C1 << std::endl;
            std::cout << "double -A2= " << B2 << "; oct -A2= " << C2 << std::endl;
            std::cout << "double -A3= " << B3 << "; oct -A3= " << C3 << std::endl;

            std::cout << "\n double area1= " << area1 << " double area2= " << area2 << " double total = " << area1+ area2 << " sum differance "<< abs(area1+area2-1)  << std:: endl;
            std::cout << "oct area1= " << area4 << " oct area2= " << area3 << "; oct total = " << area4 + area3 << " sum differance "<< abs(area3+area4-1) << std:: endl;
            std::cout << "\n differance double vs oct " << abs(D1 - A1) << " " << abs(D2 - A2) << " " << abs(D3 - A3) << " ---failed--- " << count + 1 << std::endl;
            count++;
      }
  }
  t = clock() - t;
  std::cout << "Time taken " << (Type)(t) / CLOCKS_PER_SEC << std::endl;

  return 1;
}


template <class Type>
int Sign(const Type &a) {
  return (a > 0) ? 1 : -1;
}

template <class Type>
void GetIntervalall(const std::vector <Type> &a1, const std::vector <Type> &a2, std::vector< std::pair<Type, Type> > &I1, std::vector< std::pair<Type, Type> > &I2, std::vector<std::pair<Type, Type>> &I3) {
  I1.resize(0);
  I2.resize(0);
  I3.resize(0);
  std::vector <Type> x(6);
  x[0] = 0;
  Type y = 1;
  unsigned cnt = 1;
  Type delta;
  for(unsigned k = 0; k < 2; k++) {
    const std::vector<Type> &a = (k == 0) ? a1 : a2;
//     std::cout << a[0] << " " << a[1] << " " << a[2] << std::endl;
    delta = a[1] * a[1] - 4 * a[0] * a[2];
//     std::cout << "delta " << delta << std::endl;
    if(delta >= 0) {
      Type sign = Sign(a[0]);
      for(unsigned i = 0; i < 2; i++) {

        if(a1[0] != 0) {
          Type y = (-a[1] - sign * sqrt(delta)) / (2 * a[0]);
//           std::cout << "y1 = " << y << std::endl;
          if(y >= 1) break;
          else if(y > 0) {
            x[cnt] = y;
            cnt++;
          }
        }
        else {
          Type y = -a[2] / a[1];
//           std::cout << "y2 = " << y << std::endl;
          if(y < 1 && y > 0) {
            x[cnt] = y;
            cnt++;
            break;
          }
        }
//         std::cout << "y = " << y << std::endl;

        sign *= -1;
      }
    }
  }
  x[cnt] = 1 ;
  cnt++;

  x.resize(cnt);
  std::sort(x.begin(), x.end());
//   for(unsigned i = 0; i < cnt; i++) {
// //     std::cout << "x = " << x[i] << std::endl;
//   }
  for(unsigned i = 0 ; i < cnt - 1 ; i++) {
    Type xm = (x[i] + x[i + 1]) / 2;
    Type f1 = a1[0] * xm * xm + a1[1] * xm + a1[2] ;
    Type f2 = a2[0] * xm * xm + a2[1] * xm + a2[2] ;
    if(f1 > 0) {
      if(f2 > 0) {
        I3.resize(I3.size() + 1, std::pair<Type, Type>(x[i], x[i + 1]));
      }
      else {
        I1.resize(I1.size() + 1, std::pair<Type, Type>(x[i], x[i + 1]));
      }
    }
    else if(f2 > 0) {
      I2.resize(I2.size() + 1, std::pair<Type, Type>(x[i], x[i + 1]));
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
