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

double powr(const double &x, const int &y){
  double x_p_y = 1.;
  if (y > 0){
    for(unsigned i=0; i < y; i++){
      x_p_y *= x ;
    }
  }

  else if(y<0) {
    for(unsigned i=0; i < -y; i++){
      x_p_y *= x ;
    }
    x_p_y = 1./x_p_y ;
  }

  else {
  x_p_y = 1. ;
  }
  return x_p_y ;
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
    B[q] = A[q] * (powr(k[1], q) * powr(k[0], s + n + 1 - q)) / (factorial<double>(q) * factorial<double>(s + n + 1 - q));
    A[q] *= (powr(k[1], q) * powr(k[2], s + n + 1 - q)) / (factorial<double>(q) * factorial<double>(s + n + 1 - q));

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
          sum += A[q] * powr(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
//           std::cout<< " sum= "<< sum << std::endl;
        }
        int r_m_n = r-n;
//         std::cout<< " r_m_n = "<< r_m_n << std::endl;
        A2 += (r != n) ? sum  * (powr(u2, r_m_n) - powr(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
//         std::cout<< " powr(u2, r - n)= "<< powr(u2, r_m_n) << std::endl;
//         std::cout<< " powr(u1, r - n)= "<< powr(u1, r_m_n) << std::endl;
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
          sum += A[q] * powr(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        int r_m_n = r-n;
        A2 += sum  * (powr(u2, r_m_n) - powr(u1, r_m_n)) / r_m_n;
      }
//         std::cout<< "2. A2= "<< A2 << std::endl;

// 3
      for(unsigned r = m + 1; r <= qMax + m; r++) {
        double sum = 0.;
        for(unsigned q = r - m; q <= qMax; q++) {
          sum += A[q] * powr(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        int r_m_n = r-n;
        A2 += sum  * (powr(u2, r_m_n) - powr(u1, r_m_n)) / r_m_n;
      }
//         std::cout<< "3. A2= "<< A2 << std::endl;

// 4

      for(unsigned r = 0; r < qMax; r++) {
        double sum = 0.;
        for(unsigned q = 0; q <= r; q++) {
          sum += B[q] * powr(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        A2 += sum  * (powr(u2, qMax + s + m - r + 1) - powr(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
      }
//         std::cout<< "4. A2= "<< A2 << std::endl;

// 5
      for(unsigned r = qMax; r <= m; r++) {
        double sum = 0.;
        for(unsigned q = 0; q < qMax; q++) {
          sum += B[q] * powr(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        A2 += sum  * (powr(u2, qMax + s + m - r + 1) - powr(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
      }
//         std::cout<< "5. A2= "<< A2 << std::endl;

// 6
      for(unsigned r = m + 1; r < qMax + m; r++) {
        double sum = 0.;
        for(unsigned q = r - m; q < qMax; q++) {
          sum += B[q] * powr(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        A2 += sum  * (powr(u2, qMax + s + m - r + 1) - powr(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
      }
//         std::cout<< "6. A2= "<< A2 << std::endl;

//total
      A2 *= powr(-1, n) * factorial<double>(n) * factorial<double>(m) / powr(a , m);
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
          sum += A[q] * powr(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
//           std::cout<< " sum= "<< sum << std::endl;
        }
        int r_m_n = r-n;
//         std::cout<< " r_m_n = "<< r_m_n << std::endl;
        A2 += (r != n) ? sum  * (powr(u2, r_m_n) - powr(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
//         std::cout<< " powr(u2, r - n)= "<< powr(u2, r_m_n) << std::endl;
//         std::cout<< " powr(u1, r - n)= "<< powr(u1, r_m_n) << std::endl;
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
          sum += A[q] * powr(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        int r_m_n = r-n;
        A2 += (r != n) ? sum  * (powr(u2, r_m_n) - powr(u1, r_m_n)) / r_m_n : sum  * (log(u2/u1));
      }
//         std::cout<< "2. A2= "<< A2 << std::endl;

// 3
      for(unsigned r = qMax + 1; r <= qMax + m; r++) {
        double sum = 0.;
        for(unsigned q = r - m; q <= qMax; q++) {
          sum += A[q] * powr(-c, m - r + q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        int r_m_n = r-n;
        A2 += sum  * (powr(u2, r_m_n) - powr(u1, r_m_n)) / r_m_n;
      }
//         std::cout<< "3. A2= "<< A2 << std::endl;

// 4

      for(unsigned r = 0; r <= m; r++) {
        double sum = 0.;
        for(unsigned q = 0; q <= r; q++) {
          sum += B[q] * powr(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        A2 += sum  * (powr(u2, qMax + s + m - r + 1) - powr(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
      }
//         std::cout<< "4. A2= "<< A2 << std::endl;

// 5
      for(unsigned r = m+1; r < qMax; r++) {
        double sum = 0.;
        for(unsigned q = r-m; q <= r; q++) {
          sum += B[q] * powr(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        A2 += sum  * (powr(u2, qMax + s + m - r + 1) - powr(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
      }
//         std::cout<< "5. A2= "<< A2 << std::endl;

// 6
      for(unsigned r = qMax; r < qMax + m; r++) {
        double sum = 0.;
        for(unsigned q = r - m; q < qMax; q++) {
          sum += B[q] * powr(-c, r - q) / (factorial<double>(m - r + q) * factorial<double>(r - q));
        }
        A2 += sum  * (powr(u2, qMax + s + m - r + 1) - powr(u1, qMax + s + m - r + 1)) / (qMax + s + m - r + 1);
      }
//         std::cout<< "6. A2= "<< A2 << std::endl;

//total
      A2 *= powr(-1, n) * factorial<double>(n) * factorial<double>(m) / powr(a , m);
//       std::cout<< "final. A2= "<< A2 << std::endl;

    }
return A2;;
  }

}



int main() {
unsigned int m=4;

unsigned int n = 5;
int s = 3;
double a = 1;
double c = -2;

std::vector <double> pol1{-1,1,1};

std::vector< std::pair<double, double> > I2(1) ;
I2[0].first = 0;
I2[0].second = 1;

clock_t t = clock();
std::srand(std::time(NULL));
for(unsigned i=0;i<100000;i++){

  double A2 = integral_A2(m, n, s, a, c, pol1, I2);

}
   double A2 = integral_A2(m, n, s, a, c, pol1, I2);

   std::cout<< " A2= "<< A2 << std::endl;
   t = clock() - t;
   std::cout << "\nTime taken for predetermined cases: " <<(double)(t)/ CLOCKS_PER_SEC << std::endl;
   return 1;
}





