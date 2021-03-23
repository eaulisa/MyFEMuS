
#include "FemusInit.hpp"

using namespace femus;

#include <boost/math/special_functions/factorials.hpp>
#include <algorithm>

double LimLi(const int &n, double &x) {
  if(x < 0.) return 0.;
  else if(n == 0 && x == 0.) return -0.5;
  else return -pow(x, n) / boost::math::factorial<double>(n);
}

int main(int argc, char** args) {

  unsigned Nmax = 20;

  std::vector < std::vector< double > > Cl(Nmax);
  std::vector < std::vector< double > > Cr(Nmax);

  std::vector < double > M(Nmax);

  for(unsigned i = 0; i < Nmax; i++) {
    M[i] = (i % 2 == 0) ? 2. / (i + 1.) : 0.;

    Cl[i].resize(i + 1);
    Cr[i].resize(i + 1);

    Cl[i][i] = (i % 2 != 0) ? 1 : -1;
    Cr[i][i] = (i % 2 != 0) ? -1 : 1;
    for(unsigned j = 0; j < i; j++) {
      Cl[i][j] = -Cl[i - 1][j];
      Cr[i][j] = Cr[i - 1][j];
    }
  }

  for(unsigned i = 1; i < Nmax; i++) {
    double f = i;
    Cl[i][1] *= f;
    Cr[i][1] *= f;
    for(unsigned j = 2; j <= i; j++) {
      f *= (i + 1 - j);
      Cl[i][j] *= f;
      Cr[i][j] *= f;
    }
  }


  double a0 = cos(M_PI/4.) * sin(M_PI/3.); 
  double b0 = sin(M_PI/4.) * sin(M_PI/3.);
  double c0 = cos(M_PI/3.); 
  double d0 = 0.;

  std::vector <double> a = {a0, b0, c0};
  double d = d0;
  std::vector <unsigned> m = {0, 1, 0};
  unsigned dim = 3;


  std::vector <double> x(pow(2, dim));
  for(unsigned j = 0; j < pow(2, dim); j++) {
    std::vector < unsigned> ii(dim);
    unsigned jp = j;
    for(unsigned l = 0; l < dim; l++) {
      ii[dim - 1 - l] = jp % 2;
      jp = jp / 2;
    }
    x[j] = d;
    for(unsigned l = 0; l < dim; l++) {
      x[j] += ii[l] * a[l] - (1 - ii[l]) * a[l];
    }
    std::cout << x[j] << " ";
  }
  std::cout << std::endl;


  std::vector <unsigned> i(dim, 0);

  std::vector <double> D(dim, 0);
  D[0] = 1;
  for(unsigned j = 1; j < dim; j++) {
    D[j] = D[j - 1] * fabs(Cl[m[j - 1]][0]) / pow(a[j - 1], 1);
  }


  double HCI = 0.;

  unsigned k = dim - 1;
  for(i[k] = 0; i[0] < m[0] + 1; i[k]++) {
    for(unsigned j = k + 1; j < dim; j++) {
      D[j] = D[j - 1] * fabs(Cl[m[j - 1]][i[j - 1]]) / pow(a[j - 1], i[j - 1] + 1);
    }
    k = dim - 1;
    double DD = D[k] * fabs(Cl[m[k]][i[k]]) / pow(a[k], i[k] + 1);

    std::cout << std::endl << i[0] << " " << i[1] << " " << DD << std::endl << std::flush;

    int s = 0;
    for(unsigned l = 0; l < dim; l++) {
      s += i[l] + 1;
    }
    std::cout << s << std::endl;

    for(unsigned j = 0; j < pow(2, dim); j++) {
      std::vector < unsigned> ii(dim);
      unsigned jp = j;
      for(unsigned l = 0; l < dim; l++) {
        ii[dim - 1 - l] = jp % 2;
        jp = jp / 2;
      }
      unsigned e = dim-1;
      for(unsigned l = 0; l < dim; l++) {
        e += (1-ii[l]) * m[l] + ii[l] * (i[l] + 1);
      }
      double sign = (e % 2 == 0) ? 1 : -1;

      HCI  += DD * sign * LimLi(s,x[j]);

      std::cout<<"A "<< DD <<"\n B "<<sign <<" "<< LimLi(s,x[j]) << " "<< HCI <<std::endl;


    }
//     std::cout << k << std::endl;
//     unsigned p;
//     std::cin >> p;

    while(k > 0 && i[k] == m[k]) {

      i[k--] = 0;
    }

  }

  std::cout << "Integral Value = " << HCI << std::endl;


  return 1;



}
