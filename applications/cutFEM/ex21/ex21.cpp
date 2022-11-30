#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <algorithm>    // std::sort
#include <ctime>
#include <cstdlib>
#include <climits>
using namespace std;

unsigned int factorial(int n){
    long factorial = 1.0;

    if (n < 0)
        cout << "Error! Factorial of a negative number doesn't exist.";
    else {
        for(int i = 1; i <= n; ++i) {
            factorial *= i;
        }
    }
    return factorial ;

}
void integral_A2(const double a, const double c, const std::vector <double> &pol1, std::vector< std::pair<double, double> > &I2){

int m,s,n;
double A2;

std::vector <double> A;
std::vector <double> k(3);
k[0] = pol1[0]/(a*a);
k[1] = pol1[1]/a - 2*pol1[0]*pol1[2]/(a*a);
k[2]=(pol1[0]*c*c - a*pol1[1]*c + a*pol1[2])/(a*a);



for(int q=0; q<= (s+n+1); q++){
 double first_term =(pow(k[1],q)* pow(k[2], s+n+1-q))/(factorial(q) * factorial(s+n+1-q));

  for(int r=1; r <= q/2; r++){
    double next_term = first_term* k[0]*k[2]*(q-2*r+1)*(q-2*r+2)/(r*(s+n+r-q)*k[1]*k[1]);
    A[q] = first_term + next_term ;
    first_term = next_term ;

  }
}
if(m >= s+n+1){


}
else {
}





}










int main (){}


















// double b(const double &p) {
//   return std::max(std::min(p, 1.), 0.);
// }
//
// void GetIntervalall(const std::vector <double> &a1, const std::vector <double> &a2, std::vector< std::pair<double, double> > &I1, std::vector< std::pair<double, double> > &I2, std::vector<std::pair<double, double>> &I3) {
//
//   std::vector <double> x(6);
//   x[0] = 0 ;
//   unsigned cnt = 1;
//   for(unsigned k = 0; k < 2; k++) {
//     const std::vector<double> &a = (k == 0) ? a1 : a2;
//     double delta = a[1] * a[1] - 4 * a[0] * a[2] ;
//     if(delta > 0) {
//       for(unsigned i = 0; i < 2; i++) {
//         double y = (-a[1] - pow(-1, i) * sqrt(delta)) / (2 * a[0]) ;
//         if(y > 0 && y < 1) {
//           x[cnt] = y ;
//           cnt++ ;
//         }
//       }
//     }
//   }
//   x[cnt] = 1 ;
//   cnt++;
//
//   x.resize(cnt);
//   std::sort(x.begin(), x.end());
//   I1.resize(0);
//   I2.resize(0);
//   I3.resize(0);
//   for(unsigned i = 0 ; i < cnt - 1 ; i++) {
//     double xm = (x[i] + x[i + 1]) * 0.5 ;
//     double f1 = a1[0] * xm * xm + a1[1] * xm + a1[2] ;
//     double f2 = a2[0] * xm * xm + a2[1] * xm + a2[2] ;
//     if(f1 > 0) {
//       if(f2 > 0) {
//         I3.resize(I3.size() + 1, std::pair<double, double>(x[i], x[i + 1]));
//       }
//       else {
//         I1.resize(I1.size() + 1, std::pair<double, double>(x[i], x[i + 1]));
//       }
//     }
//     else if(f2 > 0) {
//       I2.resize(I2.size() + 1, std::pair<double, double>(x[i], x[i + 1]));
//     }
//   }
//
// //This part gives segmentation-fault-core-dumped if one of the region is empty. Solved : it had two problem
// // 1. every time it erases the size goes down but i increases.
// //2. for some reason in case of empty I1, it passes the first for loop. (i=0;i<-1;i++) . If we rewrite for loop as
// // for(unsigned i = 1; i < I1.size(); i++) {
// //     if(I1[i-1].second == I1[i].first) {
// //       I1[i-1].second = I1[i].second;
// //       I1.erase(I1.begin() + i );
// //       i--;
// //     }
// //   }
//
//   for(unsigned i = 1; i < I1.size(); i++) {
//     if(I1[i-1].second == I1[i].first) {
//       I1[i-1].second = I1[i].second;
//       I1.erase(I1.begin() + i );
//       i--;
//     }
//   }
//   for(unsigned i = 1; i < I2.size(); i++) {
//     if(I2[i-1].second == I2[i].first) {
//       I2[i-1].second = I2[i].second;
//       I2.erase(I2.begin() + i );
//       i--;
//     }
//   }
//   for(unsigned i = 1; i < I3.size(); i++) {
//     if(I3[i-1].second == I3[i].first) {
//       I3[i-1].second = I3[i].second;
//       I3.erase(I3.begin() + i );
//       i--;
//     }
//   }
// }
//
//
//
// void random_polynomial(std::vector <double> &a1, std::vector <double> &a2){
//     a1[0] = ((double(std::rand()) / double(RAND_MAX)) * (4)) -2;
//     a1[1] = ((double(std::rand()) / double(RAND_MAX)) * (4)) -2;
//     a1[2] = ((double(std::rand()) / double(RAND_MAX)) * (4)) -2;
//     a2[0] = a1[0] ;
//     a2[1] = ((double(std::rand()) / double(RAND_MAX)) * (4)) -2;
//     a2[2] = ((double(std::rand()) / double(RAND_MAX)) * (4)) -2;
// //         a1[0] = 1;
// //     a1[1] = -0.75;
// //     a1[2] = .125;
// //     a2[0] = a1[0] ;
// //     a2[1] = 0.75 ;
// //     a2[2] = -0.05;
//   }
//
// void get_roots(std::vector <double> a, double &delta, double &p, double &q ){
//     delta = a[1] * a[1] - 4 * a[0] * a[2] ;
// //   std::cout <<"delta= "<< delta << std::endl;
//     if(delta > 0) {
//         p = (-a[1] + sqrt(delta)) / (2 * a[0]) ;
//         q = (-a[1] - sqrt(delta)) / (2 * a[0]) ;
//         if(p > q) {
//           std::swap(p,q);
//         }
// //         std::cout <<"roots are  "<< p << " and "<< q << std::endl;
//     }
// }





/*
int main() {
      double p1=0,p2=0,q1=0,q2=0;
      std::vector <double> a2(3), a1(3);
      std::vector< std::pair<double, double> > I1, I2, I3 ;

      clock_t t = clock();
      std::srand(std::time(NULL));
  for(unsigned i=0;i<10;i++){
      random_polynomial(a1,a2);
        std::cout <<"polynomial p(x) = "<< a1[0] << "x^2 + (" << a1[1] << "x) + (" << a1[2] << ") " <<std::endl;
        std::cout <<"polynomial q(x) = "<< a2[0] << "x^2 + (" << a2[1] << "x) + (" << a2[2] << ") " <<std::endl;

      double delta1,delta2;
      get_roots(a1,delta1,p1,p2);
      get_roots(a2,delta2,q1,q2);

        std::cout <<"roots p1 = " << p1 << " & p2 = " << p2 << "d="<<delta1 <<std::endl;
        std::cout <<"roots q1 = " << q1 << " & q2 = " << q2 << "d="<<delta2 <<std::endl;

      GetIntervalall(a1,a2,I1,I2,I3);

      std::cout << "I1= " ;
      for (unsigned i=0 ; i<I1.size(); i++){
         std::cout << "(" << I1[i].first << "," << I1[i].second <<") U " ;
      }
      std::cout << "\nI2= " ;
      for (unsigned i=0 ; i<I2.size(); i++){
         std::cout << "(" << I2[i].first << "," << I2[i].second <<") U " ;
      }

      std::cout << "\nI3= " ;
      for (unsigned i=0 ; i<I3.size(); i++){
         std::cout << "(" << I3[i].first << "," << I3[i].second <<") U " <<std::endl;
      }
      t = clock() - t;
      std::cout << "\nTime taken for generalized algorithm : " << (double)(t)/ CLOCKS_PER_SEC << std::endl;

      if (delta1 > 0){
        if(delta2 > 0){
          GetInterval4Old( p1, p2, q1, q2, a1[0], I1, I2, I3);
        }
        else {
          bool pIsComplex=0;
          GetInterval2(p1, p2, pIsComplex, a1[0], I1, I2, I3);
        }
      }
      else{
        if (delta2 > 0){
          bool pIsComplex=1;
          GetInterval2(q1, q2, pIsComplex, a1[0], I1, I2, I3);
        }
        else {
           GetInterval0(a1[0], I1, I2, I3);
        }
      }
      std::cout << "I1= " ;
      for (unsigned i=0 ; i<I1.size(); i++){
        std::cout << "(" << I1[i].first << "," << I1[i].second <<") U " ;
      }
      std::cout << "\nI2= " ;
      for (unsigned i=0 ; i<I2.size(); i++){
         std::cout << "(" << I2[i].first << "," << I2[i].second <<") U " ;
      }
      std::cout << "\nI3= " ;
      for (unsigned i=0 ; i<I3.size(); i++){
         std::cout << "(" << I3[i].first << "," << I3[i].second <<") U " ;
      }
   t = clock() - t;
   std::cout << "\nTime taken for predetermined cases: " <<(double)(t)/ CLOCKS_PER_SEC << std::endl;
  }
  return 1;
}*/
