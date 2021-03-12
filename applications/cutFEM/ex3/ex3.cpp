
#include "FemusInit.hpp"


using namespace femus;
using namespace std;
using std::cout;

#include "../eqPoly/LiSK/lisk.hpp"
#include <algorithm>

int main(int argc, char** args) {

  unsigned Nmax = 20;

  std::vector < std::vector< double > > Cl(Nmax);
  std::vector < std::vector< double > > Cr(Nmax);

  std::vector < double > M(Nmax);

  for(unsigned i = 0; i < Nmax; i++) {
    M[i] = (i % 2 == 0) ? -2. / (i + 1.) : 0.;

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

  double a = 0.34, b = 0.74, c = 0.25, d = 0.5, t = 20., degree = 5;
  unsigned N = degree + 1;
  LiSK::LiSK< complex<double> > *lisk = new LiSK::LiSK< complex<double> > (N + 3);
  {
    // 1D
    if(N > Nmax) {
      std::cout << "Warning degree is too high, for 1D the maximum degree allowed is " << Nmax - 1 << std::endl << std::flush;
      abort();
    }

    double e[2];
    for(unsigned i = 0; i < 2; i++) {
      double sa = (i == 0) ? -a : a;
      e[i] = -exp((sa + d) * t);
    }
    std::vector< double > Li[2];

    for(unsigned i = 0; i < 2; i++) {
      Li[i].resize(N);
      for(unsigned l = 0; l < N; l++) {
        Li[i][l] = (lisk->Li(1 + l, e[i])).real();
      }
    }

    double at2m1 = 1. / (a * t);
    std::vector< double > at2m(N);
    at2m[0] = at2m1;
    for(unsigned i = 1; i < N; i++) {
      at2m[i] = at2m[i - 1] * at2m1;
    }

    std::vector< double > D(N, 0.);
    for(unsigned i = 0; i < N; i++) {
      D[i] = M[i];
      for(unsigned l = 0; l <= i; l++) {
        D[i] +=  -2. * (Cl[i][l] * Li[0][l] + Cr[i][l] * Li[1][l]) * at2m[l];
      }
    }

    std::cout << " 1D \n";
    for(unsigned i = 0; i < N; i++) {
      std::cout << D[i] << std::endl;
    }
  }


  // 2D
  {
    if(N + 1 > Nmax) {
      std::cout << "Warning degree is too high, for 2D the maximum degree allowed is " << Nmax - 2 << std::endl << std::flush;
      abort();
    }

    double e[2][2]; //left or right limits
    for(unsigned i = 0; i < 2; i++) {
      double sa = (i == 0) ? -a : a;
      for(unsigned j = 0; j < 2; j++) {
        double sb = (j == 0) ? -b : b;
        e[i][j] = -exp((sa + sb + d) * t);
      }
    }
    std::vector< double > Li[2][2]; //left or right limits
    for(unsigned i = 0; i < 2; i++) {
      for(unsigned j = 0; j < 2; j++) {
        Li[i][j].resize(N);
        for(unsigned l = 0; l < N; l++) {
          Li[i][j][l] = (lisk->Li(2 + l, e[i][j])).real();
        }
      }
    }


    double at2m1 = 1. / (a * t);
    double bt2m1 = 1. / (b * t);
    std::vector< double > at2m(N);
    std::vector< double > bt2m(N);
    at2m[0] = at2m1;
    bt2m[0] = bt2m1;
    for(unsigned i = 1; i < N; i++) {
      at2m[i] = at2m[i - 1] * at2m1;
      bt2m[i] = bt2m[i - 1] * bt2m1;
    }

    std::vector< std::vector< double > > Dl(N);
    std::vector< std::vector< double > > Dr(N);

    for(unsigned i = 0; i < N; i++) { // int y^i
      Dl[i].assign(N - i, 0.);
      Dr[i].assign(N - i, 0.);
      for(unsigned j = 0; j < N - i; j++) {
        for(unsigned l = 0; l < std::min(N - j, i + 1); l++) { // Li[j + l]
          Dl[i][j] += (Cl[i][l] * Li[0][0][j + l] + Cr[i][l] * Li[0][1][j + l]) * bt2m[l];
          Dr[i][j] += (Cl[i][l] * Li[1][0][j + l] + Cr[i][l] * Li[1][1][j + l]) * bt2m[l];
        }
      }
    }

    std::vector < std::vector< double > > D(N);
    for(unsigned i = 0; i < N; i++) { // int x^i
      D[i].resize(N - i);
      for(unsigned j = 0; j < N - i; j++) { // int f(y) dy
        D[i][j] = M[i] * M[j];
        for(unsigned l = 0; l < std::min(N - j, i + 1); l++) {
          D[i][j] += -2. * (Cl[i][l] * Dl[j][l] + Cr[i][l] * Dr[j][l]) * at2m[l];
        }
      }
    }

    std::cout << " 2D \n";
    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = 0; j < N - i; j++) {
        std::cout << D[i][j] << " ";
      }
      std::cout << std::endl;
    }
  }


  // 3D
  {
    if(N + 2 > Nmax) {
      std::cout << "Warning degree is too high, for 3D the maximum degree allowed is " << Nmax - 3 << std::endl << std::flush;
      abort();
    }


    std::cout << "debug \n \n";

    double e[2][2][2]; //left or right limits
    for(unsigned i = 0; i < 2; i++) {
      double sa = (i == 0) ? -a : a;
      for(unsigned j = 0; j < 2; j++) {
        double sb = (j == 0) ? -b : b;
        for(unsigned k = 0; k < 2; k++) {
          double sc = (k == 0) ? -c : c;
          e[i][j][k] = -exp((sa + sb + sc + d) * t);
        }
      }
    }
    std::vector< double > Li[2][2][2]; //left or right limits
    for(unsigned i = 0; i < 2; i++) {
      for(unsigned j = 0; j < 2; j++) {
        for(unsigned k = 0; k < 2; k++) {
          Li[i][j][k].resize(N);
          for(unsigned l = 0; l < N; l++) {
            Li[i][j][k][l] = (lisk->Li(3 + l, e[i][j][k])).real();
            std::cout <<  Li[i][j][k][l] << "\n";
          }
        }
      }
    }


    double at2m1 = 1. / (a * t);
    double bt2m1 = 1. / (b * t);
    double ct2m1 = 1. / (c * t);
    std::vector< double > at2m(N);
    std::vector< double > bt2m(N);
    std::vector< double > ct2m(N);
    at2m[0] = at2m1;
    bt2m[0] = bt2m1;
    ct2m[0] = ct2m1;
    for(unsigned i = 1; i < N; i++) {
      at2m[i] = at2m[i - 1] * at2m1;
      bt2m[i] = bt2m[i - 1] * bt2m1;
      ct2m[i] = ct2m[i - 1] * ct2m1;
    }
    std::cout << at2m[0] << " " << bt2m[0] << " " << ct2m[0] << std::endl;

    std::vector< std::vector< double > > Dll(N);
    std::vector< std::vector< double > > Dlr(N);
    std::vector< std::vector< double > > Drl(N);
    std::vector< std::vector< double > > Drr(N);

    for(unsigned i = 0; i < N; i++) { //int z^i
      Dll[i].resize(N - i);
      Dlr[i].resize(N - i);
      Drl[i].resize(N - i);
      Drr[i].resize(N - i);
      for(unsigned j = 0; j < N - i; j++) { // Li[j+l]
        for(unsigned l = 0; l < std::min(N - j, i + 1); l++) {
          Dll[i][j] += (Cl[i][l] * Li[0][0][0][j + l] + Cr[i][l] * Li[0][0][1][j + l]) * ct2m[l];
          Dlr[i][j] += (Cl[i][l] * Li[0][1][0][j + l] + Cr[i][l] * Li[0][1][1][j + l]) * ct2m[l];
          Drl[i][j] += (Cl[i][l] * Li[1][0][0][j + l] + Cr[i][l] * Li[1][0][1][j + l]) * ct2m[l];
          Drr[i][j] += (Cl[i][l] * Li[1][1][0][j + l] + Cr[i][l] * Li[1][1][1][j + l]) * ct2m[l];
        }
      }
    }

    std::vector< std::vector< double > > Dl(N);
    std::vector< std::vector< double > > Dr(N);

    for(unsigned i = 0; i < N; i++) { //int y^i
      Dl[i].resize(N - i);
      Dr[i].resize(N - i);
      for(unsigned j = 0; j < N - i; j++) { // int f(z)
        for(unsigned l = 0; l < std::min(N - j, i + 1); l++) {
          Dl[i][j] += (Cl[i][l] * Dll[j][l] + Cr[i][l] * Dlr[j][l]) * bt2m[l];
          Dr[i][j] += (Cl[i][l] * Drl[j][l] + Cr[i][l] * Drr[j][l]) * bt2m[l];
        }
      }
    }


    std::vector< std::vector< std::vector< double > > > D(N);
    for(unsigned i = 0; i < N; i++) {
      D[i].resize(N - i);
      for(unsigned j = 0; j < N - i; j++) {
        D[i][j].resize(N - i - j);
        for(unsigned k = 0; k < N - i - j; k++) {
          D[i][j][k] = M[i] * M[j] * M[k];
          for(unsigned l = 0; l < std::min(N - j - k, i + 1); l++) {
            D[i][j][k] += -2. * (Cl[i][l] * Dl[k+j][l] + Cr[i][l] * Dr[k+j][l]) * at2m[l];
          }
        }
      }
    }


    std::cout << " 3D \n";
    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = 0; j < N - i; j++) {
        for(unsigned k = 0; k < N - i - j; k++) {
          std::cout << D[i][j][k] << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }



  delete lisk;
  return 0;
}


