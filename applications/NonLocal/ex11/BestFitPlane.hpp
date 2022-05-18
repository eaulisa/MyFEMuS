
#include "MyEigenFunctions.hpp"

void GetNormalTriBF(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double & R, std::vector<double> &a, double & d,  std::vector<double> &xm, std::vector<double> &b, double & db, unsigned & cut) {

  const unsigned &dim =  xv.size();
  const unsigned &nve =  xv[0].size();

  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];




  std::vector<double> A(2, 0.);
  std::vector<std::vector<double>> xe(4, std::vector<double>(2));
  double D = 0.;
  unsigned intMax = 2;
  unsigned nEdge = 0;
  unsigned cnt = 0;

  for(unsigned i = 0; i < nve; i++) {
    unsigned ip1 = (i + 1) % nve;
    A[0] = xv[1][ip1] - xv[1][i];
    A[1] = - xv[0][ip1] + xv[0][i];
    D = - A[0] * xv[0][i] - A[1] * xv[1][i];


    std::vector<double> inters(intMax, 0.);
    unsigned dir = (fabs(A[0]) > fabs(A[1])) ? 1 : 0 ;
    unsigned dirp1 = (dir + 1) % 2;

    double iMax = std::max(xv[dir][ip1], xv[dir][i]);
    double iMin = std::min(xv[dir][ip1], xv[dir][i]);

    double delta = ((A[0] * A[0] + A[1] * A[1]) * R * R) - (D + A[0] * xg[0] + A[1] * xg[1]) * (D + A[0] * xg[0] + A[1] * xg[1]);
    a.resize(dim);
    if(delta > 0.) {
      inters[0] = (- A[dir] * (D + A[dirp1] * xg[dirp1]) + A[dirp1] * (A[dirp1] * xg[dir] - sqrt(delta))) / (A[0] * A[0] + A[1] * A[1]);
      inters[1] = (- A[dir] * (D + A[dirp1] * xg[dirp1]) + A[dirp1] * (A[dirp1] * xg[dir] + sqrt(delta))) / (A[0] * A[0] + A[1] * A[1]);
      unsigned nInt = 0;
      unsigned jInt = 2;
      for(unsigned j = 0; j < intMax; j++) {
        if(inters[j] < iMax && inters[j] > iMin) {
          nInt++;
          jInt = j;
        }
      }
      if(nInt == 1) {
        xe[cnt][dir] = inters[jInt];
        xe[cnt][dirp1] = (- D - A[dir] * xe[cnt][dir]) / A[dirp1];
        cnt++;
      }
    }
  }
  if(cnt == 0) {
    cut = (R * R - (xv[0][0] - xg[0]) * (xv[0][0] - xg[0]) - (xv[1][0] - xg[1]) * (xv[1][0] - xg[1]) > 0) ? 0 : 2;
    return;
  }
  else if(cnt == 4) {
    cut = 0;
    return;
  }
  else if(cnt == 2) {
    cut = 1;
    std::vector<double> w = {1., 1., 4.};
    std::vector<double> theta(3);
    xe.resize(3);

    theta[0] = atan2(xe[0][1]  - xg[1], xe[0][0] - xg[0]);
    theta[1] = atan2(xe[1][1]  - xg[1], xe[1][0] - xg[0]);

    if(theta[0] > theta[1]) {
      std::swap(theta[0], theta[1]);
    }
    double DT = theta[1] - theta[0];
    if(DT > M_PI) {
      std::swap(theta[0], theta[1]);
      theta[1] += 2. * M_PI;
      DT = theta[1] - theta[0];
    }

    theta[2] = theta[0] + 0.5 * DT;

    xe[2][0] = xg[0] + R * cos(theta[2]);
    xe[2][1] = xg[1] + R * sin(theta[2]);

    std::vector<double> N = {-cos(theta[2]), -sin(theta[2])};

    femus::FindBestFit(xe, w, N, a, d);

    xm.resize(dim);
    xm[0] = -(d * N[0] + a[1] * (N[0] * xg[1] - N[1] * xg[0])) / (a[0] * N[0] + a[1] * N[1]);
    xm[1] = -(d * N[1] + a[0] * (N[1] * xg[0] - N[0] * xg[1])) / (a[0] * N[0] + a[1] * N[1]);

    std::vector<double> xi(dim);

    std::vector < std::vector < double > > J(2, std::vector<double>(2));
    J[0][0] = (-x1 + x2);
    J[0][1] = (-x1 + x3);

    J[1][0] = (-y1 + y2);
    J[1][1] = (-y1 + y3);

    double den = (x3 * y1 - x1 * y3 + x2 * J[1][1] - y2 * J[0][1]);

    xi[0] = (x3 * y1 - x1 * y3 + xm[0] * J[1][1] - xm[1] * J[0][1]) / den;
    xi[1] = (x1 * y2 - x2 * y1 - xm[0] * J[1][0] + xm[1] * J[0][0]) / den;

    b.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        b[k] += J[j][k] * a[j];
      }
    }

    double bNorm = sqrt(b[0] * b[0] + b[1] * b[1]);
    b[0] /= bNorm;
    b[1] /= bNorm;
    db = - b[0] * xi[0] - b[1] * xi[1];

  }
}


void GetNormalTetBF(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double & R, std::vector<double> &a, double & d,  std::vector<double> &xm, std::vector<double> &a2, double & d2, double & volume,  unsigned & cut) {

  const unsigned dim =  3;
  const unsigned nve =  4;

  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& x4 = xv[0][3];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];
  const double& y4 = xv[1][3];
  const double& z1 = xv[2][0];
  const double& z2 = xv[2][1];
  const double& z3 = xv[2][2];
  const double& z4 = xv[2][3];

  double hx = (fabs(x2 - x1) + fabs(x3 - x2) + fabs(x3 - x1) + fabs(x4 - x1) + fabs(x4 - x2) + fabs(x4 - x3)) / 6.;
  double hy = (fabs(y2 - y1) + fabs(y3 - y2) + fabs(y3 - y1) + fabs(y4 - y1) + fabs(y4 - y2) + fabs(y4 - y3)) / 6.;
  double hz = (fabs(z2 - z1) + fabs(z3 - z2) + fabs(z3 - z1) + fabs(z4 - z1) + fabs(z4 - z2) + fabs(z4 - z3)) / 6.;
  double h = sqrt(hx * hx + hy * hy + hz * hz);
  double eps = 1.0e-10 * h;

  std::vector<double> A(dim, 0.);
  std::vector < std::vector <double> > y(6, std::vector<double>(dim));
  std::vector < unsigned > i0(6);
  double D = 0.;
  unsigned intMax = 2;
  unsigned nEdge = 0;
  unsigned cnt = 0;

  for(unsigned i = 0; i < nve - 1; i++) {
    for(unsigned j = i + 1; j < nve; j++) {
      for(unsigned k = 0; k < dim; k++) {
        A[k] = xv[k][j] - xv[k][i];
      }

      std::vector<double> inters(intMax, 0.);
      unsigned dir = (fabs(A[0]) > fabs(A[1])) ? ((fabs(A[0]) > fabs(A[2])) ? 0 : 2) : ((fabs(A[1]) > fabs(A[2])) ? 1 : 2) ;
      unsigned dirp1 = (dir + 1) % dim;
      unsigned dirp2 = (dir + 2) % dim;

      double iMax = std::max(xv[dir][j], xv[dir][i]);
      double iMin = std::min(xv[dir][j], xv[dir][i]);

      double Axdi[3] = {A[1]* (xg[2] - xv[2][i]) - A[2] * (xg[1] - xv[1][i]),
                        A[2]* (xg[0] - xv[0][i]) - A[0] * (xg[2] - xv[2][i]),
                        A[0]* (xg[1] - xv[1][i]) - A[1] * (xg[0] - xv[0][i])
                       } ;

      double Addi = A[0] * (xg[0] - xv[0][i]) +
                    A[1] * (xg[1] - xv[1][i]) +
                    A[2] * (xg[2] - xv[2][i]);


      double den = A[0] * A[0] + A[1] * A[1] + A[2] * A[2];

      double delta = (- (A[dir] * A[dir]) *
                      (Axdi[0] * Axdi[0] + Axdi[1] * Axdi[1] + Axdi[2] * Axdi[2] - R * R * den));

      double var = den * xv[dir][i] + A[dir] * Addi;

      a.resize(dim);
      if(delta > 0.) {
        inters[0] = (var - sqrt(delta)) / den;
        inters[1] = (var + sqrt(delta)) / den;
        unsigned nInt = 0;
        unsigned jInt = 2;
        for(unsigned ii = 0; ii < intMax; ii++) {
          if(inters[ii] < iMax && inters[ii] > iMin) {
            nInt++;
            jInt = ii;
          }
        }
        if(nInt == 1) {
          y[cnt][dir] = inters[jInt];
          y[cnt][dirp1] = xv[dirp1][i] + A[dirp1] * (y[cnt][dir] - xv[dir][i]) / A[dir];
          y[cnt][dirp2] = xv[dirp2][i] + A[dirp2] * (y[cnt][dir] - xv[dir][i]) / A[dir];
          i0[cnt] = (i + j) - (i == 0);
          cnt++;
        }
      }
    }
  }
  if(cnt == 0) {
    cut = (R * R - (xv[0][0] - xg[0]) * (xv[0][0] - xg[0]) - (xv[1][0] - xg[1]) * (xv[1][0] - xg[1])  - (xv[2][0] - xg[2]) * (xv[2][0] - xg[2]) > 0) ? 0 : 2;
    return;
  }
  else if(cnt > 4) { // bad: small ball
    cut = 0;
    return;
  }
  else if(cnt == 4 || cnt == 3) {
    cut = 1;

    std::vector<double> w(cnt + 1, 1.);
    w[cnt] = 2 * cnt;

    y.resize(cnt + 1);

    std::vector <double> yg(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned i = 0; i < cnt; i++) {
        yg[k] += y[i][k];
      }
      yg[k] /= cnt;
    }

    std::vector < double > N(dim);// this is an approximate normal pointing toward the outside of the ball in the physical element
    for(unsigned k = 0; k < dim; k++) {
      N[k] =  yg[k] - xg[k];
    }
    
    double det = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
    for(unsigned k = 0; k < dim; k++) {
      N[k] /=det;
    }
        
    double theta = atan2(N[1], N[0]);
    double phi = acos(N[2]);

    y[cnt][0] = xg[0] + R * sin(phi) * cos(theta);
    y[cnt][1] = xg[1] + R * sin(phi) * sin(theta);
    y[cnt][2] = xg[2] + R * cos(phi);

    for(unsigned k = 0; k < dim; k++) {
      N[k] = N[k];
      std::cout<<" "<<N[k];
    }
    std::cout << std::endl;
    femus::FindBestFit(y, w, N, a, d); // a is an BF normal pointing toward the outside of the ball in the physical element

    det = a[0] * N[0] + a[1] * N[1] + a[2] * N[2];
    xm.resize(dim);
    xm[0] = -(d * N[0] + a[1] * (N[0] * xg[1] - N[1] * xg[0]) + a[2] * (N[0] * xg[2] - N[2] * xg[0])) / det;
    xm[1] = -(d * N[1] + a[2] * (N[1] * xg[2] - N[2] * xg[1]) + a[0] * (N[1] * xg[0] - N[0] * xg[1])) / det;
    xm[2] = -(d * N[2] + a[0] * (N[2] * xg[0] - N[0] * xg[2]) + a[1] * (N[2] * xg[1] - N[1] * xg[2])) / det;

    

    std::vector<double> xi(dim);
    std::vector < std::vector < double > > J(3, std::vector<double>(3));
    
    J[0][0] = (-x1 + x2);
    J[0][1] = (-x1 + x3);
    J[0][2] = (-x1 + x4);

    J[1][0] = (-y1 + y2);
    J[1][1] = (-y1 + y3);
    J[1][2] = (-y1 + y4);

    J[2][0] = (-z1 + z2);
    J[2][1] = (-z1 + z3);
    J[2][2] = (-z1 + z4);

    double den =   J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
                   - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
                   + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

    volume = den / 6.;

    xi[0] = -(x3 * y4 * z1 - x3 * xm[1] * z1 - x1 * y4 * z3 + x1 * xm[1] * z3 - x3 * y1 * z4 + x1 * y3 * z4 - x1 * xm[1] * z4 + x3 * xm[1] * z4 +
              xm[0] * (y3 * z1 - y4 * z1 - y1 * z3 + y4 * z3 + y1 * z4 - y3 * z4) +
              x3 * y1 * xm[2] - x1 * y3 * xm[2] + x1 * y4 * xm[2] - x3 * y4 * xm[2] +
              x4 * (xm[1] * z1 + y1 * z3 - xm[1] * z3 - y1 * xm[2] + y3 * (-z1 + xm[2]))) / den;

    xi[1] = -(-(x2 * y4 * z1) + x2 * xm[1] * z1 + x1 * y4 * z2 - x1 * xm[1] * z2 + x2 * y1 * z4 - x1 * y2 * z4 + x1 * xm[1] * z4 - x2 * xm[1] * z4 +
              xm[0] * (-(y2 * z1) + y4 * z1 + y1 * z2 - y4 * z2 - y1 * z4 + y2 * z4) +
              (-(x2 * y1) + x1 * y2 - x1 * y4 + x2 * y4) * xm[2] +
              x4 * (-(xm[1] * z1) - y1 * z2 + xm[1] * z2 + y2 * (z1 - xm[2]) + y1 * xm[2])) / den;


    xi[2] = -(x2 * y3 * z1 - x2 * xm[1] * z1 - x1 * y3 * z2 + x1 * xm[1] * z2 - x2 * y1 * z3 + x1 * y2 * z3 - x1 * xm[1] * z3 + x2 * xm[1] * z3 +
              xm[0] * (y2 * z1 - y3 * z1 - y1 * z2 + y3 * z2 + y1 * z3 - y2 * z3) +
              x2 * y1 * xm[2] - x1 * y2 * xm[2] + x1 * y3 * xm[2] - x2 * y3 * xm[2] +
              x3 * (xm[1] * z1 + y1 * z2 - xm[1] * z2 - y1 * xm[2] + y2 * (-z1 + xm[2]))) / den;


    a2.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        a2[k] -= J[j][k] * a[j]; // this normal has to point toward the inside of the ball in the parent element (a needs to change sign)
      }
    }
    double bNorm = sqrt(a2[0] * a2[0] + a2[1] * a2[1] + a2[2] * a2[2]);
    a2[0] /= bNorm;
    a2[1] /= bNorm;
    a2[2] /= bNorm;
    d2 = - a2[0] * xi[0] - a2[1] * xi[1] - a2[2] * xi[2];

  }
}

