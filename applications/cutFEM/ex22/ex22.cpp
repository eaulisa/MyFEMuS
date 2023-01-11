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

double integral_A3(const unsigned &m, const unsigned &n, const int &s, const double &a, const double &c, const std::vector <double> &pol1, const std::vector< std::pair<double, double> > &I3) {
  double A3 = 0;
  std::vector <double> k(3);
  k[0] = pol1[0] / (a * a);
  k[1] = pol1[1] / a;
  k[2] = k[0] * c * c - k[1] * c + pol1[2] / a;
  k[1] -= 2 * c * k[0];


  for(int i=0; i <= s; i++){

    std::vector <double> A(s -i + 1, 0);   // size of all this vector changes.
    std::vector <double> B(s -i + 1, 0);


    double kterms = (k[0] * k[2]) / (k[1] * k[1]);
    unsigned qMax = s - i;

//   std::cout<< " ankor1 " << "k0 =" << k[0]<< "k1 =" << k[1]<< "k2 =" << k[2] << "kterms =" << kterms << "qMax =" << qMax << std::endl ;

    for(int q = 0; q <= qMax; q++) {
      double term = 1;
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
      B[q] = A[q] * (pow(k[1], q) * pow(k[0], qMax - q)) / (factorial<double>(q) * factorial<double>(qMax - q));
      A[q] *= (pow(k[1], q) * pow(k[2], qMax - q)) / (factorial<double>(q) * factorial<double>(qMax - q));

    }


    if(m >= qMax) {
      for(unsigned w = 0; w < I3.size(); w++)  {
        double u1 = a * I3[w].first + c;
        double u2 = a * I3[w].second + c;
                        //       std::cout<< " u1= "<< u1 << std::endl;
                        //       std::cout<< " u2= "<< u2 << std::endl;
// 1
        for(unsigned r = 0; r <= qMax; r++) {
                        //         std::cout<< " r= "<< r << std::endl;
          double sum = 0.;
          for(unsigned q = 0; q <= r; q++) {
                        //           std::cout<< " q= "<< q << std::endl;
                        //           std::cout<< " factorial<double>(m - r + q)= "<< factorial<double>(m - r + q) << std::endl;
                        //           std::cout<< " factorial<double>(r - q)= "<< factorial<double>(r - q) << std::endl;
              sum += A[q] * pow(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
                        //           std::cout<< " sum= "<< sum << std::endl;
          }
          int r_m_n = r+i+1;          //in A2 this was the power of u after integration. I am keeping the variable.
                        //         std::cout<< " r_m_n = "<< r_m_n << std::endl;
          A3 += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
                        //         std::cout<< " pow(u2, r - n)= "<< pow(u2, r_m_n) << std::endl;
                        //         std::cout<< " pow(u1, r - n)= "<< pow(u1, r_m_n) << std::endl;
                        //         std::cout<< " log(u2)= "<< log(u2) << std::endl;
                        //         std::cout<< " log(u1)= "<< log(u1) << std::endl;
                        //         std::cout<< " (log(u2/u1)= "<< (log(u2/u1)) << std::endl;
                        //
                        //         std::cout<< " r-n = "<< r-(n) << std::endl;
        }
  //         std::cout<< "1. A3= "<< A3 << std::endl;

  // 2
        for(unsigned r = qMax + 1; r <= m; r++) {
          double sum = 0.;
          for(unsigned q = 0; q <= qMax; q++) {
            sum += A[q] * pow(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
          }
          int r_m_n = r+i+1;
          A3 += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
        }
  //         std::cout<< "2. A3= "<< A3 << std::endl;

  // 3
        for(unsigned r = m + 1; r <= qMax + m; r++) {
          double sum = 0.;
          for(unsigned q = r - m; q <= qMax; q++) {
            sum += A[q] * pow(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
          }
          int r_m_n = r+i+1;
          A3 += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
        }
  //         std::cout<< "3. A3= "<< A3 << std::endl;

  // 4

        for(unsigned r = 0; r < qMax; r++) {
          double sum = 0.;
          for(unsigned q = 0; q <= r; q++) {
            sum += B[q] * pow(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
          }
          int r_m_n = qMax + s + m - r + 1;
          A3 += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
        }
  //         std::cout<< "4. A3= "<< A3 << std::endl;

  // 5
        for(unsigned r = qMax; r <= m; r++) {
          double sum = 0.;
          for(unsigned q = 0; q < qMax; q++) {
            sum += B[q] * pow(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
          }
          int r_m_n = qMax + s + m - r + 1;
          A3 += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
        }
  //         std::cout<< "5. A3= "<< A3 << std::endl;

  // 6
        for(unsigned r = m + 1; r < qMax + m; r++) {
          double sum = 0.;
          for(unsigned q = r - m; q < qMax; q++) {
            sum += B[q] * pow(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
          }
          int r_m_n = qMax + s + m - r + 1;
          A3 += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
        }
  //         std::cout<< "6. A3= "<< A3 << std::endl;

      }
  //total
        A3 *= 1. /((n+i+1)*factorial<double>(i));
  //       std::cout<< "final. A3= "<< A3 << std::endl;
    }



    else {
      for(unsigned w = 0; w < I3.size(); w++)  {
        double u1 = a * I3[w].first + c;
        double u2 = a * I3[w].second + c;
                        //       std::cout<< " u1= "<< u1 << std::endl;
                        //       std::cout<< " u2= "<< u2 << std::endl;
// 1
        for(unsigned r = 0; r <= m; r++) {
                        //         std::cout<< " r= "<< r << std::endl;
          double sum = 0.;
          for(unsigned q = 0; q <= r; q++) {
                        //           std::cout<< " q= "<< q << std::endl;
                        //           std::cout<< " factorial<double>(m - r + q)= "<< factorial<double>(m - r + q) << std::endl;
                        //           std::cout<< " factorial<double>(r - q)= "<< factorial<double>(r - q) << std::endl;
              sum += A[q] * pow(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
                        //           std::cout<< " sum= "<< sum << std::endl;
          }
          int r_m_n = r+i+1;          //in A2 this was the power of u after integration. I am keeping the variable.
                        //         std::cout<< " r_m_n = "<< r_m_n << std::endl;
          A3 += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
                        //         std::cout<< " pow(u2, r - n)= "<< pow(u2, r_m_n) << std::endl;
                        //         std::cout<< " pow(u1, r - n)= "<< pow(u1, r_m_n) << std::endl;
                        //         std::cout<< " log(u2)= "<< log(u2) << std::endl;
                        //         std::cout<< " log(u1)= "<< log(u1) << std::endl;
                        //         std::cout<< " (log(u2/u1)= "<< (log(u2/u1)) << std::endl;
                        //
                        //         std::cout<< " r-n = "<< r-(n) << std::endl;
        }
  //         std::cout<< "1. A3= "<< A3 << std::endl;

  // 2
        for(unsigned r = m + 1; r <= qMax; r++) {
          double sum = 0.;
          for(unsigned q = 0; q <= m; q++) {
            sum += A[q] * pow(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
          }
          int r_m_n = r+i+1;
          A3 += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
        }
  //         std::cout<< "2. A3= "<< A3 << std::endl;

  // 3
        for(unsigned r = qMax + 1; r <= qMax + m; r++) {
          double sum = 0.;
          for(unsigned q = r - m; q <= qMax; q++) {
            sum += A[q] * pow(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
          }
          int r_m_n = r+i+1;
          A3 += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
        }
  //         std::cout<< "3. A3= "<< A3 << std::endl;

  // 4

        for(unsigned r = 0; r < m; r++) {
          double sum = 0.;
          for(unsigned q = 0; q <= r; q++) {
            sum += B[q] * pow(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
          }
          int r_m_n = qMax + s + m - r + 1;
          A3 += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
        }
  //         std::cout<< "4. A3= "<< A3 << std::endl;

  // 5
        for(unsigned r = m + 1; r < qMax; r++) {
          double sum = 0.;
          for(unsigned q = r-m; q <= r; q++) {
            sum += B[q] * pow(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
          }
          int r_m_n = qMax + s + m - r + 1;
          A3 += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
        }
  //         std::cout<< "5. A3= "<< A3 << std::endl;

  // 6
        for(unsigned r = qMax; r < qMax + m; r++) {
          double sum = 0.;
          for(unsigned q = r - m; q < qMax; q++) {
            sum += B[q] * pow(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
          }
          int r_m_n = qMax + s + m - r + 1;
          A3 += (r_m_n != 0) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
        }
  //         std::cout<< "6. A3= "<< A3 << std::endl;

      }
  //total
        A3 *= 1. /((n+i+1)*factorial<double>(i));
  //       std::cout<< "final. A3= "<< A3 << std::endl;
    }

  }
  A3 *= factorial<double>(m) / pow(a , m);
  return A3;

}






double integral_A2(const unsigned &m, const unsigned &n, const int &s, const double &a, const double &c, const std::vector <double> &pol1, const std::vector< std::pair<double, double> > &I2) {
  std::vector <double> k(3);
  k[0] = pol1[0] / (a * a);
  k[1] = pol1[1] / a;
  k[2] = k[0] * c * c - k[1] * c + pol1[2] / a;
  k[1] -= 2 * c * k[0];

  std::vector <double> A(s + n + 2, 0);
  std::vector <double> B(s + n + 2, 0);


  double kterms = (k[0] * k[2]) / (k[1] * k[1]);
  unsigned qMax = s + n + 1;

//   std::cout<< " ankor1 " << "k0 =" << k[0]<< "k1 =" << k[1]<< "k2 =" << k[2] << "kterms =" << kterms << "qMax =" << qMax << std::endl ;

  for(int q = 0; q <= qMax; q++) {
    double term = 1;
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
    B[q] = A[q] * (pow(k[1], q) * pow(k[0], s + n + 1 - q)) / (factorial<double>(q) * factorial<double>(s + n + 1 - q));
    A[q] *= (pow(k[1], q) * pow(k[2], s + n + 1 - q)) / (factorial<double>(q) * factorial<double>(s + n + 1 - q));

  }
  double A2 = 0;
  if(m >= qMax) {
    for(unsigned i = 0; i < I2.size(); i++)  {
      double u1 = a * I2[i].first + c;
      double u2 = a * I2[i].second + c;
//       std::cout<< " u1= "<< u1 << std::endl;
//       std::cout<< " u2= "<< u2 << std::endl;
// 1
      for(unsigned r = 0; r <= qMax; r++) {
//         std::cout<< " r= "<< r << std::endl;
        double sum = 0.;
        for(unsigned q = 0; q <= r; q++) {
//           std::cout<< " q= "<< q << std::endl;
//           std::cout<< " factorial<double>(m - r + q)= "<< factorial<double>(m - r + q) << std::endl;
//           std::cout<< " factorial<double>(r - q)= "<< factorial<double>(r - q) << std::endl;
          sum += A[q] * pow(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
//           std::cout<< " sum= "<< sum << std::endl;
        }
        int r_m_n = r-n;
//         std::cout<< " r_m_n = "<< r_m_n << std::endl;
        A2 += (r != n) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
//         std::cout<< " pow(u2, r - n)= "<< pow(u2, r_m_n) << std::endl;
//         std::cout<< " pow(u1, r - n)= "<< pow(u1, r_m_n) << std::endl;
//         std::cout<< " log(u2)= "<< log(u2) << std::endl;
//         std::cout<< " log(u1)= "<< log(u1) << std::endl;
//         std::cout<< " (log(u2/u1)= "<< (log(u2/u1)) << std::endl;
//
//         std::cout<< " r-n = "<< r-(n) << std::endl;
      }
//         std::cout<< "1. A2= "<< A2 << std::endl;

// 2
      for(unsigned r = qMax + 1; r <= m; r++) {
        double sum = 0.;
        for(unsigned q = 0; q <= qMax; q++) {
          sum += A[q] * pow(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        int r_m_n = r-n;
        A2 += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
      }
//         std::cout<< "2. A2= "<< A2 << std::endl;

// 3
      for(unsigned r = m + 1; r <= qMax + m; r++) {
        double sum = 0.;
        for(unsigned q = r - m; q <= qMax; q++) {
          sum += A[q] * pow(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        int r_m_n = r-n;
        A2 += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
      }
//         std::cout<< "3. A2= "<< A2 << std::endl;

// 4

      for(unsigned r = 0; r < qMax; r++) {
        double sum = 0.;
        for(unsigned q = 0; q <= r; q++) {
          sum += B[q] * pow(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        A2 += sum  * (pow(u2, qMax + s + m - r + 1) - pow(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
      }
//         std::cout<< "4. A2= "<< A2 << std::endl;

// 5
      for(unsigned r = qMax; r <= m; r++) {
        double sum = 0.;
        for(unsigned q = 0; q < qMax; q++) {
          sum += B[q] * pow(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        A2 += sum  * (pow(u2, qMax + s + m - r + 1) - pow(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
      }
//         std::cout<< "5. A2= "<< A2 << std::endl;

// 6
      for(unsigned r = m + 1; r < qMax + m; r++) {
        double sum = 0.;
        for(unsigned q = r - m; q < qMax; q++) {
          sum += B[q] * pow(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        A2 += sum  * (pow(u2, qMax + s + m - r + 1) - pow(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
      }
//         std::cout<< "6. A2= "<< A2 << std::endl;

//total
      A2 *= pow(-1, n) * factorial<double>(n) * factorial<double>(m) / pow(a , m);
//       std::cout<< "final. A2= "<< A2 << std::endl;

    }
    return A2;
  }



  else {

    for(unsigned i = 0; i < I2.size(); i++)  {
      double u1 = a * I2[i].first + c;
      double u2 = a * I2[i].second + c;
//       std::cout<< " u1= "<< u1 << std::endl;
//       std::cout<< " u2= "<< u2 << std::endl;
// 1
      for(unsigned r = 0; r <= m; r++) {
//         std::cout<< " r= "<< r << std::endl;
        double sum = 0.;
        for(unsigned q = 0; q <= r; q++) {
//           std::cout<< " q= "<< q << std::endl;
//           std::cout<< " factorial<double>(m - r + q)= "<< factorial<double>(m - r + q) << std::endl;
//           std::cout<< " factorial<double>(r - q)= "<< factorial<double>(r - q) << std::endl;
          sum += A[q] * pow(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
//           std::cout<< " sum= "<< sum << std::endl;
        }
        int r_m_n = r-n;
//         std::cout<< " r_m_n = "<< r_m_n << std::endl;
        A2 += (r != n) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
//         std::cout<< " pow(u2, r - n)= "<< pow(u2, r_m_n) << std::endl;
//         std::cout<< " pow(u1, r - n)= "<< pow(u1, r_m_n) << std::endl;
//         std::cout<< " log(u2)= "<< log(u2) << std::endl;
//         std::cout<< " log(u1)= "<< log(u1) << std::endl;
//         std::cout<< " (log(u2/u1)= "<< (log(u2/u1)) << std::endl;
//
//         std::cout<< " r-n = "<< r-(n) << std::endl;
      }
//         std::cout<< "1. A2= "<< A2 << std::endl;

// 2
      for(unsigned r = m + 1; r <= qMax; r++) {
        double sum = 0.;
        for(unsigned q = r-m; q <= r; q++) {
          sum += A[q] * pow(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        int r_m_n = r-n;
        A2 += (r != n) ? sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
      }
//         std::cout<< "2. A2= "<< A2 << std::endl;

// 3
      for(unsigned r = qMax + 1; r <= qMax + m; r++) {
        double sum = 0.;
        for(unsigned q = r - m; q <= qMax; q++) {
          sum += A[q] * pow(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        int r_m_n = r-n;
        A2 += sum  * (pow(u2, r_m_n) - pow(u1, r_m_n)) / r_m_n;
      }
//         std::cout<< "3. A2= "<< A2 << std::endl;

// 4

      for(unsigned r = 0; r <= m; r++) {
        double sum = 0.;
        for(unsigned q = 0; q <= r; q++) {
          sum += B[q] * pow(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        A2 += sum  * (pow(u2, qMax + s + m - r + 1) - pow(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
      }
//         std::cout<< "4. A2= "<< A2 << std::endl;

// 5
      for(unsigned r = m+1; r < qMax; r++) {
        double sum = 0.;
        for(unsigned q = r-m; q <= r; q++) {
          sum += B[q] * pow(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        A2 += sum  * (pow(u2, qMax + s + m - r + 1) - pow(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
      }
//         std::cout<< "5. A2= "<< A2 << std::endl;

// 6
      for(unsigned r = qMax; r < qMax + m; r++) {
        double sum = 0.;
        for(unsigned q = r - m; q < qMax; q++) {
          sum += B[q] * pow(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        A2 += sum  * (pow(u2, qMax + s + m - r + 1) - pow(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
      }
//         std::cout<< "6. A2= "<< A2 << std::endl;

//total
      A2 *= pow(-1, n) * factorial<double>(n) * factorial<double>(m) / pow(a , m);
//       std::cout<< "final. A2= "<< A2 << std::endl;

    }
return A2;
  }

}


void GetIntervalall(const std::vector <double> &a1, const std::vector <double> &a2, std::vector< std::pair<double, double> > &I1, std::vector< std::pair<double, double> > &I2, std::vector<std::pair<double, double>> &I3) {

  std::vector <double> x(6);
  x[0] = 0 ;
  unsigned cnt = 1;
  double delta;

  for(unsigned k = 0; k < 2; k++) {
    const std::vector<double> &a = (k == 0) ? a1 : a2;
    delta = a[1] * a[1] - 4 * a[0] * a[2] ;
    if(delta > 0) {
      double sign = 1;
      for(unsigned i = 0; i < 2; i++) {
        double y = (-a[1] - sign * sqrt(delta)) / (2. * a[0]) ;

        if(y >= 1) break;
        else if(y > 0) {
          x[cnt] = y ;
          cnt++;
        }
        sign *= -1.;
      }
    }
  }
  x[cnt] = 1 ;
  cnt++;

  x.resize(cnt);
  std::sort(x.begin(), x.end());
  I1.resize(0);
  I2.resize(0);
  I3.resize(0);
  for(unsigned i = 0 ; i < cnt - 1 ; i++) {
    double xm = (x[i] + x[i + 1]) * 0.5 ;
    double f1 = a1[0] * xm * xm + a1[1] * xm + a1[2] ;
    double f2 = a2[0] * xm * xm + a2[1] * xm + a2[2] ;
    if(f1 > 0) {
      if(f2 > 0) {
        I3.resize(I3.size() + 1, std::pair<double, double>(x[i], x[i + 1]));
      }
      else {
        I1.resize(I1.size() + 1, std::pair<double, double>(x[i], x[i + 1]));
      }
    }
    else if(f2 > 0) {
      I2.resize(I2.size() + 1, std::pair<double, double>(x[i], x[i + 1]));
    }
  }
}


void random_polynomial(std::vector <double> &a1, std::vector <double> &a2) {
  a1[0] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
  a1[1] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
  a1[2] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
  a2[0] = a1[0] ;
  a2[1] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
  a2[2] = ((double(std::rand()) / double(RAND_MAX)) * (4)) - 2;
//         a1[0] = 1;
//     a1[1] = -0.75;
//     a1[2] = .125;
//     a2[0] = a1[0] ;
//     a2[1] = 0.75 ;
//     a2[2] = -0.05;
}



int main() {
  unsigned int m = 9;
  unsigned int n = 5;
  int s = 3;


  std::vector <double> pol1;
  std::vector <double> pol2;

  std::srand(std::time(NULL));
  random_polynomial(pol1, pol2);

  double a = pol2[1] - pol1[1];
  double c = pol2[2] - pol1[2];


//   pol1[0] = -1.53164;
//   pol1[1] = -1.2347;
//   pol1[2] = 1.0401;
//   a = -0.527101;
//   c = -2.00759;

  std::vector< std::pair <double, double> > I1;
  std::vector< std::pair <double, double> > I2;
  std::vector< std::pair <double, double> > I3;
  GetIntervalall(pol1, pol2, I1, I2, I3);

//   I2.resize(1);
//   I2[0].first = 0;
//   I2[0].second = 0.514288;

  std::cout << "k = " << pol1[0] << "; b = " << pol1[1] << "; d = " << pol1[2] << "; a = " << a << "; c = " << c << ";" << std::endl;
  for(unsigned i = 0; i < I2.size(); i++) {
    std::cout << "x1 = " << I2[i].first << "; x2 = " << I2[i].second << ";" << std::endl;
  }
  for(unsigned i = 0; i < I3.size(); i++) {
    std::cout << "y1 = " << I3[i].first << "; y2 = " << I3[i].second << ";" << std::endl;
  }
  {
    clock_t t = clock();
    std::srand(std::time(NULL));
    for(unsigned i = 0; i < 1; i++) {
//      Type A2 = integral_A2(m, n, s, a, c, pol1, I2);
    }
    double A2 = integral_A2(m, n, s, a, c, pol1, I2);
    double A3= integral_A3(m, n, s, a, c, pol1, I3);

    std::cout << "double A2= " << A2 << std::endl;
    std::cout << "double A3= " << A3 << std::endl;
    t = clock() - t;
    std::cout << "\nTime taken for predetermined cases: " << (double)(t) / CLOCKS_PER_SEC << std::endl;
  }

  return 1;
}










