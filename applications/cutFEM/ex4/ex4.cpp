
#include "FemusInit.hpp"

using namespace femus;

#include <boost/math/special_functions/factorials.hpp>
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

  //double a0 = 0.34, b0 = 0.74, c0 = 0.25, d0 = 0.5, degree = 3;
  double a0 = 1, b0 = 1, c0 = 1., d0 = 0., degree = 3;
  bool surface = true;
  unsigned N = degree + 1;

  {
    // 1D line
    double a = a0;
    double d = d0;
    double scale = sqrt(a * a);
    a /= scale;
    d /= scale;

    if(N > Nmax) {
      std::cout << "Warning degree is too high, for 1D the maximum degree allowed is " << Nmax - 1 << std::endl << std::flush;
      abort();
    }

    double e[2];
    e[0] = std::max(0., -a + d);
    e[1] = std::max(0., a + d);

    std::vector< double > Li[2];
    bool edgeA[2] = {false, false};
    for(unsigned i = 0; i < 2; i++) {
      if(e[i] > 0.) {
        edgeA[i] = true;
        Li[i].resize(N);
        Li[i][0] = (surface) ? -1 : - e[i];
        for(unsigned l = 1; l < N; l++) {
          Li[i][l] = Li[i][l - 1] * e[i] / (1. - surface + l);
        }
      }
    }

    std::vector< double > a2m(N);
    a2m[0] = 1. / a;
    for(unsigned i = 1; i < N; i++) {
      a2m[i] = a2m[i - 1] * a2m[0];
    }

    std::vector< double > f(N, 0.);
    for(unsigned i = 0; i < N; i++) {
      f[i] = 0.;//-M[i];
      for(unsigned l = 0; l <= i; l++) {
        if(edgeA[0]) f[i] +=  -1. * (Cl[i][l] * Li[0][l]) * a2m[l];
        if(edgeA[1]) f[i] +=  -1. * (Cr[i][l] * Li[1][l]) * a2m[l];
      }
    }

    std::cout << " 1D element line, degree " << degree << ", size " << f.size() << std::endl;
    for(unsigned i = 0; i < f.size(); i++) {
      std::cout << "F[" << i << "] = " << f[i] << std::endl;
    }

  }


  // 2D square
  {
    if(N + 1 > Nmax) {
      std::cout << "Warning degree is too high, for 2D the maximum degree allowed is " << Nmax - 2 << std::endl << std::flush;
      abort();
    }

    double a = a0;
    double b = b0;
    double d = d0;
    double scale = sqrt(a * a + b * b);

    a /= scale;
    b /= scale;
    d /= scale;

    double e[2][2]; //left or right limits
    double sa = -a;
    for(unsigned i = 0; i < 2; i++, sa *= -1.) {
      double sb = -b;
      for(unsigned j = 0; j < 2; j++, sb *= -1.) {
        e[i][j] = std::max(0., sa + sb + d);
      }
    }
    std::vector< double > Li[2][2]; //left or right limits
    bool edgeB[2][2] = {{false, false}, {false, false}}; //left or right limits
    for(unsigned i = 0; i < 2; i++) {
      for(unsigned j = 0; j < 2; j++) {
        if(e[i][j] > 0.) {
          edgeB[i][j] = true;
          Li[i][j].resize(N);
          Li[i][j][0] = (surface) ? - e[i][j] : - e[i][j] * e[i][j] / 2.;
          for(unsigned l = 1; l < N; l++) {
            Li[i][j][l] = Li[i][j][l - 1]  * e[i][j] / (2. - surface + l);
          }
        }
      }
    }
    bool edgeA[2];
    edgeA[0] = edgeB[0][0] + edgeB[0][1];
    edgeA[1] = edgeB[1][0] + edgeB[1][1];

    std::vector< double > a2m(N);
    std::vector< double > b2m(N);
    a2m[0] = 1. / a;
    b2m[0] = 1. / b;
    for(unsigned i = 1; i < N; i++) {
      a2m[i] = a2m[i - 1] * a2m[0];
      b2m[i] = b2m[i - 1] * b2m[0];
    }

    std::vector < std::vector< double > > D(N);
    for(unsigned i = 0; i < N; i++) {
      D[i].assign(N - i, 0.);
    }

    std::vector< std::vector< double > > Dl = D, Dr = D;

    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = 0; j < N - i; j++) {
        for(unsigned l = 0; l < std::min(N - j, i + 1); l++) { //Transpse[D_s] = Bl.Li_sl + Br.Li_sr
          if(edgeB[0][0]) Dl[j][i] += (Cl[i][l] * Li[0][0][l + j]) * b2m[l];
          if(edgeB[0][1]) Dl[j][i] += (Cr[i][l] * Li[0][1][l + j]) * b2m[l];
          if(edgeB[1][0]) Dr[j][i] += (Cl[i][l] * Li[1][0][l + j]) * b2m[l];
          if(edgeB[1][1]) Dr[j][i] += (Cr[i][l] * Li[1][1][l + j]) * b2m[l];
        }
      }
    }

    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = 0; j < N - i; j++) {
        D[i][j] = 0.;//-M[i] * M[j];
        for(unsigned l = 0; l < std::min(N - j, i + 1); l++) { //D = Al.D_s + Ar.D_r
          if(edgeA[0]) D[i][j] +=  -1. * (Cl[i][l] * Dl[l][j]) * a2m[l];
          if(edgeA[1]) D[i][j] +=  -1. * (Cr[i][l] * Dr[l][j]) * a2m[l];
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


// 3D cube
  {
    if(N + 2 > Nmax) {
      std::cout << "Warning degree is too high, for 3D the maximum degree allowed is " << Nmax - 3 << std::endl << std::flush;
      abort();
    }


    double a = a0;
    double b = b0;
    double c = c0;
    double d = d0;

    double scale = sqrt(a * a + b * b + c * c);

    a /= scale;
    b /= scale;
    c /= scale;
    d /= scale;

    double e[2][2][2]; //left or right limits
    double sa = -a;
    for(unsigned i = 0; i < 2; i++, sa *= -1.) {
      double sb = -b;
      for(unsigned j = 0; j < 2; j++, sb *= -1.) {
        double sc = -c;
        for(unsigned k = 0; k < 2; k++, sc *= -1.) {
          e[i][j][k] =  std::max(0., sa + sb + sc + d);
        }
      }
    }

    std::vector< double > Li[2][2][2]; //left or right limits
    bool edgeC[2][2][2] = {{{false, false}, {false, false}}, {{false, false}, {false, false}}};
    for(unsigned i = 0; i < 2; i++) {
      for(unsigned j = 0; j < 2; j++) {
        for(unsigned k = 0; k < 2; k++) {
          if(e[i][j][k] > 0.) {
            edgeC[i][j][k] = true;
            Li[i][j][k].resize(N);
            Li[i][j][k][0] = (surface) ? -e[i][j][k] * e[i][j][k] / 2. : -e[i][j][k] * e[i][j][k] * e[i][j][k] / 6.;
            for(unsigned l = 1; l < N; l++) {
              Li[i][j][k][l] = Li[i][j][k][l - 1] * e[i][j][k] / (3. - surface + l);
            }
          }
        }
      }
    }

    bool edgeB[2][2];
    edgeB[0][0] = edgeC[0][0][0] + edgeC[0][0][1];
    edgeB[0][1] = edgeC[0][1][0] + edgeC[0][1][1];
    edgeB[1][0] = edgeC[1][0][0] + edgeC[1][0][1];
    edgeB[1][1] = edgeC[1][1][0] + edgeC[1][1][1];

    bool edgeA[2];
    edgeA[0] = edgeB[0][0] + edgeB[0][1];
    edgeA[1] = edgeB[1][0] + edgeB[1][1];

    std::vector< double > a2m(N), b2m(N), c2m(N);
    a2m[0] = 1. / a;
    b2m[0] = 1. / b;
    c2m[0] = 1. / c;
    for(unsigned i = 1; i < N; i++) {
      a2m[i] = a2m[i - 1] * a2m[0];
      b2m[i] = b2m[i - 1] * b2m[0];
      c2m[i] = c2m[i - 1] * c2m[0];
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
            if(edgeC[0][0][0]) Dll[j][i][k] += (Cl[i][l] * Li[0][0][0][l + j + k]) * c2m[l];
            if(edgeC[0][0][1]) Dll[j][i][k] += (Cr[i][l] * Li[0][0][1][l + j + k]) * c2m[l];
            if(edgeC[0][1][0]) Dlr[j][i][k] += (Cl[i][l] * Li[0][1][0][l + j + k]) * c2m[l];
            if(edgeC[0][1][1]) Dlr[j][i][k] += (Cr[i][l] * Li[0][1][1][l + j + k]) * c2m[l];
            if(edgeC[1][0][0]) Drl[j][i][k] += (Cl[i][l] * Li[1][0][0][l + j + k]) * c2m[l];
            if(edgeC[1][0][1]) Drl[j][i][k] += (Cr[i][l] * Li[1][0][1][l + j + k]) * c2m[l];
            if(edgeC[1][1][0]) Drr[j][i][k] += (Cl[i][l] * Li[1][1][0][l + j + k]) * c2m[l];
            if(edgeC[1][1][1]) Drr[j][i][k] += (Cr[i][l] * Li[1][1][1][l + j + k]) * c2m[l];
          }
        }
      }
    }

    std::vector<std::vector< std::vector< double > > > Dl = D, Dr = D;
    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = 0; j < N - i; j++) {
        for(unsigned k = 0; k < N - i - j; k++) {
          for(unsigned l = 0; l < std::min(N - j - k, i + 1); l++) { // Transpose[D_s, 3<->1] = Bl.D_sl + Br.D_sr
            if(edgeB[0][0]) Dl[k][j][i] += (Cl[i][l] * Dll[l][j][k]) * b2m[l];
            if(edgeB[0][1]) Dl[k][j][i] += (Cr[i][l] * Dlr[l][j][k]) * b2m[l];
            if(edgeB[1][0]) Dr[k][j][i] += (Cl[i][l] * Drl[l][j][k]) * b2m[l];
            if(edgeB[1][1]) Dr[k][j][i] += (Cr[i][l] * Drr[l][j][k]) * b2m[l];
          }
        }
      }
    }

    for(unsigned i = 0; i < N; i++) {
      for(unsigned j = 0; j < N - i; j++) {
        for(unsigned k = 0; k < N - i - j; k++) {
          D[i][k][j] = - 0 * M[i] * M[j] * M[k];
          for(unsigned l = 0; l < std::min(N - j - k, i + 1); l++) { // Transpose[D, 2<->3] = Al.D_l + Ar.D_r
            if(edgeA[0]) D[i][k][j] += -1. * (Cl[i][l] * Dl[l][j][k]) * a2m[l];
            if(edgeA[1]) D[i][k][j] += -1. * (Cr[i][l] * Dr[l][j][k]) * a2m[l];
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

  return 0;
}






