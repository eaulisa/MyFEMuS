
#include "FemusInit.hpp"

using namespace femus;


#include "../eqPoly/LiSK/lisk.hpp"
#include <algorithm>

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

  //double a = 0.34, b = 0.74, c = 0.25, d = 0.5, t = 20., degree = 4;
  double a = 1., b = 1., c = 1., d = 0., t = 20., degree = 3;
  unsigned N = degree + 1;
  LiSK::LiSK< std::complex<double> > *lisk = new LiSK::LiSK< std::complex<double> > (N + 3);
  {
    // 1D
    if(N > Nmax) {
      std::cout << "Warning degree is too high, for 1D the maximum degree allowed is " << Nmax - 1 << std::endl << std::flush;
      abort();
    }

    double e[2];
    e[0] = -exp((-a + d) * t);
    e[1] = -exp((a + d) * t);

    std::vector< double > Li[2];
    for(unsigned i = 0; i < 2; i++) {
      Li[i].resize(N);
      for(unsigned l = 0; l < N; l++) {
        Li[i][l] = (lisk->Li(1 + l, e[i])).real();
      }
    }

    std::vector< double > at2m(N);
    at2m[0] = 1. / (a * t);
    for(unsigned i = 1; i < N; i++) {
      at2m[i] = at2m[i - 1] * at2m[0];
    }

    std::vector< double > f(N, 0.);
    for(unsigned i = 0; i < N; i++) {
      f[i] = -M[i];
      for(unsigned l = 0; l <= i; l++) {
        f[i] +=  -2. * (Cl[i][l] * Li[0][l] + Cr[i][l] * Li[1][l]) * at2m[l];
      }
    }

    std::cout << " 1D element line, degree " << degree << ", size " << f.size() << std::endl;
    for(unsigned i = 0; i < f.size(); i++) {
      std::cout << "F[" << i << "] = " << f[i] << std::endl;
    }

  }


  // 2D
  {
    if(N + 1 > Nmax) {
      std::cout << "Warning degree is too high, for 2D the maximum degree allowed is " << Nmax - 2 << std::endl << std::flush;
      abort();
    }

    double e[2][2]; //left or right limits
    double sa = -a;
    for(unsigned i = 0; i < 2; i++, sa *= -1.) {
      double sb = -b;
      for(unsigned j = 0; j < 2; j++, sb *= -1.) {
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

    std::vector< double > at2m(N);
    std::vector< double > bt2m(N);
    at2m[0] = 1. / (a * t);
    bt2m[0] = 1. / (b * t);
    for(unsigned i = 1; i < N; i++) {
      at2m[i] = at2m[i - 1] * at2m[0];
      bt2m[i] = bt2m[i - 1] * bt2m[0];
    }

    std::vector < std::vector< double > > D(N);
    for(unsigned i = 0; i < N; i++) {
      D[i].assign(N - i, 0.);
    }

    std::vector< std::vector< double > > Dl = D, Dr = D;

    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = 0; j < N - i; j++) {
        for(unsigned l = 0; l < std::min(N - j, i + 1); l++) { //Transpse[D_s] = Bl.Li_sl + Br.Li_sr
          Dl[j][i] += (Cl[i][l] * Li[0][0][l + j] + Cr[i][l] * Li[0][1][l + j]) * bt2m[l];
          Dr[j][i] += (Cl[i][l] * Li[1][0][l + j] + Cr[i][l] * Li[1][1][l + j]) * bt2m[l];
        }
      }
    }

    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = 0; j < N - i; j++) {
        D[i][j] = -M[i] * M[j];
        for(unsigned l = 0; l < std::min(N - j, i + 1); l++) { //D = Al.D_s + Ar.D_r
          D[i][j] +=  -2. * (Cl[i][l] * Dl[l][j] + Cr[i][l] * Dr[l][j]) * at2m[l];
        }
      }
    }



    unsigned cnt = N * (N + 1) / 2;
    std::vector < std::vector < unsigned > > idx(cnt);
    std::vector < double > f(cnt);
    for(unsigned  i = 0; i < idx.size(); i++) idx[i].resize(2);

    cnt = 0u;
    for(unsigned l = 0 ; l < N; l++) {
      for(int i = l; i >= 0; i--) {

        idx[cnt][0] = i;
        idx[cnt][1] = l - i;

        f[cnt] = D[i][l - i];
        cnt++;
      }
    }

    std::cout << " 2D element quad, degree " << degree << ", size " << f.size() << std::endl;
    for(unsigned  i = 0; i < f.size(); i++) {
      std::cout << "F[" << idx[i][0] << "][" << idx[i][1] << "] = " << f[i] << std::endl;
    }
  }


// 3D
  {
    if(N + 2 > Nmax) {
      std::cout << "Warning degree is too high, for 3D the maximum degree allowed is " << Nmax - 3 << std::endl << std::flush;
      abort();
    }

    double e[2][2][2]; //left or right limits
    double sa = -a;
    for(unsigned i = 0; i < 2; i++, sa *= -1.) {
      double sb = -b;
      for(unsigned j = 0; j < 2; j++, sb *= -1.) {
        double sc = -c;
        for(unsigned k = 0; k < 2; k++, sc *= -1.) {
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
          }
        }
      }
    }

    std::vector< double > at2m(N), bt2m(N), ct2m(N);
    at2m[0] = 1. / (a * t);
    bt2m[0] = 1. / (b * t);
    ct2m[0] = 1. / (c * t);
    for(unsigned i = 1; i < N; i++) {
      at2m[i] = at2m[i - 1] * at2m[0];
      bt2m[i] = bt2m[i - 1] * bt2m[0];
      ct2m[i] = ct2m[i - 1] * ct2m[0];
    }

    std::vector< std::vector< std::vector< double > > > D(N);
    for(unsigned i = 0; i < N; i++) {
      D[i].resize(N - i);
      for(unsigned j = 0; j < N - i; j++) {
        D[i][j].assign(N - i - j, 0.);
      }
    }
    std::vector<std::vector< std::vector< double > > > Dll = D, Dlr = D, Drl = D, Drr = D;

    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = 0; j < N - i; j++) {
        for(unsigned k = 0; k < N - i - j; k++) {
          for(unsigned l = 0; l < std::min(N - j - k, i + 1); l++) { // Transpose[D_st, 1<->2] = Cl.Li_stl + Cr.Li_str
            Dll[j][i][k] += (Cl[i][l] * Li[0][0][0][l + j + k] + Cr[i][l] * Li[0][0][1][l + j + k]) * ct2m[l];
            Dlr[j][i][k] += (Cl[i][l] * Li[0][1][0][l + j + k] + Cr[i][l] * Li[0][1][1][l + j + k]) * ct2m[l];
            Drl[j][i][k] += (Cl[i][l] * Li[1][0][0][l + j + k] + Cr[i][l] * Li[1][0][1][l + j + k]) * ct2m[l];
            Drr[j][i][k] += (Cl[i][l] * Li[1][1][0][l + j + k] + Cr[i][l] * Li[1][1][1][l + j + k]) * ct2m[l];
          }
        }
      }
    }

    std::vector<std::vector< std::vector< double > > > Dl = D, Dr = D;
    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = 0; j < N - i; j++) {
        for(unsigned k = 0; k < N - i - j; k++) {
          for(unsigned l = 0; l < std::min(N - j - k, i + 1); l++) { // Transpose[D_s, 3<->1] = Bl.D_sl + Br.D_sr
            Dl[k][j][i] += (Cl[i][l] * Dll[l][j][k] + Cr[i][l] * Dlr[l][j][k]) * bt2m[l];
            Dr[k][j][i] += (Cl[i][l] * Drl[l][j][k] + Cr[i][l] * Drr[l][j][k]) * bt2m[l];
          }
        }
      }
    }

    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = 0; j < N - i; j++) {
        for(unsigned k = 0; k < N - i - j; k++) {
          D[i][k][j] = -M[i] * M[j] * M[k];
          for(unsigned l = 0; l < std::min(N - j - k, i + 1); l++) { // Transpose[D, 2<->3] = Al.D_l + Ar.D_r
            D[i][k][j] += -2. * (Cl[i][l] * Dl[l][j][k]  + Cr[i][l] * Dr[l][j][k]) * at2m[l];
          }
        }
      }
    }



    unsigned cnt = N * (N + 1) * (N + 2) / 6;
    std::vector < std::vector < unsigned > > idx(cnt);
    std::vector < double > f(cnt);
    for(unsigned  i = 0; i < idx.size(); i++) idx[i].resize(3);

    cnt = 0u;
    for(unsigned l = 0 ; l < N; l++) {
      for(int i = l; i >= 0; i--) {
        for(int j = l - i; j >= 0; j--) {
          idx[cnt][0] = i;
          idx[cnt][1] = j;
          idx[cnt][2] = l - i - j;
          f[cnt] = D[i][j][l - i - j];
          cnt++;
        }
      }
    }

    std::cout << " 3D element hex, degree " << degree << ", size " << f.size() << std::endl;
    for(unsigned  i = 0; i < idx.size(); i++) {
      std::cout << "F[" << idx[i][0] << "][" << idx[i][1] << "][" << idx[i][2] << "] = " << f[i] << std::endl;
    }
  }

  delete lisk;
  return 0;
}




