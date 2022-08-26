
#include "MyEigenFunctions.hpp"

void GetNormalTriBF(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double & R, std::vector<double> &a, double & d,  std::vector<double> &xm, std::vector<double> &b, double & db, unsigned & cut) {

  const unsigned &dim =  xv.size();
  const unsigned nve =  3;

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
  const unsigned ndge =  6;

  std::vector < std::vector <double> > y(ndge, std::vector<double>(dim));
  Eigen::Vector3d A;
  unsigned cnt = 0;

  std::vector < Eigen::Vector3d> x(nve);
  std::vector < double > dist2(nve, R * R);
  for(unsigned i = 0; i < nve; i++) {
    for(unsigned k = 0; k < dim; k++) {
      x[i][k] = xv[k][i] - xg[k];
      dist2[i] -= x[i][k] * x[i][k];
    }
  }

  std::vector<std::vector <unsigned>> vt = {{1, 2, 3}, {2, 3}, {3}};

  for(unsigned i = 0; i < vt.size(); i++) {
    for(unsigned jj = 0; jj < vt[i].size(); jj++) {
      unsigned j = vt[i][jj];
      if(dist2[i] * dist2[j] <= 0) {
        for(unsigned k = 0; k < dim; k++) {
          A[k] = x[j][k] - x[i][k];
        }
        std::vector<double> inters(2, 0.);
        unsigned dir = (fabs(A[0]) > fabs(A[1])) ? ((fabs(A[0]) > fabs(A[2])) ? 0 : 2) : ((fabs(A[1]) > fabs(A[2])) ? 1 : 2) ;
        unsigned dirp1 = (dir + 1) % dim;
        unsigned dirp2 = (dir + 2) % dim;

        double yMax = std::max(x[j][dir], x[i][dir]);
        double yMin = std::min(x[j][dir], x[i][dir]);

        Eigen::Vector3d Axdi = -A.cross(x[i]);
        double Addi = -A.dot(x[i]);
        double den = A.dot(A);

        double delta = - (A[dir] * A[dir]) * (Axdi.dot(Axdi) - R * R * den);
        double var = den * x[i][dir] + A[dir] * Addi;

        a.resize(dim);
        if(delta > 0.) {
          bool test = false;
          for(unsigned ii = 0; ii < 2; ii++) {
            double yi = (var + ((ii == 0) ? -1 : +1) * sqrt(delta)) / den;
            if(yMin < yi && yi < yMax) { // this should exclude intersections in the vertices
              if(test == false) {
                y[cnt][dir] = yi;
                y[cnt][dirp1] = x[i][dirp1] + A[dirp1] * (y[cnt][dir] - x[i][dir]) / A[dir];
                y[cnt][dirp2] = x[i][dirp2] + A[dirp2] * (y[cnt][dir] - x[i][dir]) / A[dir];
                cnt++;
                test = true;
              }
              else { // this is to remove double intersection on the same edge
                cnt--;
              }
            }
          }
        }
      }
    }
  }
  if(cnt == 0) {
    cut = (dist2[0] > 0) ? 0 : 2;
    return;
  }
  else if(cnt > 4) { // bad: small ball
    cut = 0;
    return;
  }
  else if(cnt == 4 || cnt == 3) {// these I am not sure are the right limits
    cut = 1;

    std::vector<double> w(cnt + 1, 1.);
    w[cnt] = 2. * cnt;
    y.resize(cnt + 1, {0., 0., 0.});

    std::vector <double> N(dim, 0); // this is an approximate normal pointing toward the outside of the ball in the physical element
    double det = 0.;
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned i = 0; i < cnt; i++) {
        N[k] += y[i][k];
      }
      det += N[k] * N[k];
    }
    det = sqrt(det);
    for(unsigned k = 0; k < dim; k++) {
      N[k] /=  det;
    }

    double theta = atan2(N[1], N[0]);
    double phi = acos(N[2]);

    y[cnt][0] = R * sin(phi) * cos(theta);
    y[cnt][1] = R * sin(phi) * sin(theta);
    y[cnt][2] = R * cos(phi);

    femus::FindBestFit(y, w, N, a, d); // a is the unit BF normal pointing toward the outside of the ball in the physical element

    xm.resize(dim);
    xm[0] = -(d * a[0]);
    xm[1] = -(d * a[1]);
    xm[2] = -(d * a[2]);

    const double& x1 = x[0][0];
    const double& x2 = x[1][0];
    const double& x3 = x[2][0];
    const double& x4 = x[3][0];
    const double& y1 = x[0][1];
    const double& y2 = x[1][1];
    const double& y3 = x[2][1];
    const double& y4 = x[3][1];
    const double& z1 = x[0][2];
    const double& z2 = x[1][2];
    const double& z3 = x[2][2];
    const double& z4 = x[3][2];

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


void GetNormalHexBF(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double & R, std::vector<double> &a, double & d,  std::vector<double> &xm, std::vector<double> &a2, double & d2, double & volume,  unsigned & cut, const elem_type *fem) {

  const unsigned dim =  3;
  const unsigned nve =  8;
  const unsigned ndge =  12;

  std::vector < std::vector <double> > y(ndge, std::vector<double>(dim));
  Eigen::Vector3d A;
  unsigned cnt = 0;

  std::vector < Eigen::Vector3d> x(nve);
  std::vector < double > dist2(nve, R * R);
  for(unsigned i = 0; i < nve; i++) {
    for(unsigned k = 0; k < dim; k++) {
      x[i][k] = xv[k][i] - xg[k];
      dist2[i] -= x[i][k] * x[i][k];
    }
  }

  std::vector<std::vector <unsigned>> vt = {{1, 3, 4}, {2, 5}, {3, 6}, {7}, {5, 7}, {6}, {7}};

  for(unsigned i = 0; i < vt.size(); i++) {
    for(unsigned jj = 0; jj < vt[i].size(); jj++) {
      unsigned j = vt[i][jj];
      if(dist2[i] * dist2[j] <= 0) {
        for(unsigned k = 0; k < dim; k++) {
          A[k] = x[j][k] - x[i][k];
        }
        std::vector<double> inters(2, 0.);
        unsigned dir = (fabs(A[0]) > fabs(A[1])) ? ((fabs(A[0]) > fabs(A[2])) ? 0 : 2) : ((fabs(A[1]) > fabs(A[2])) ? 1 : 2) ;
        unsigned dirp1 = (dir + 1) % dim;
        unsigned dirp2 = (dir + 2) % dim;

        double yMax = std::max(x[j][dir], x[i][dir]);
        double yMin = std::min(x[j][dir], x[i][dir]);

        Eigen::Vector3d Axdi = -A.cross(x[i]);
        double Addi = -A.dot(x[i]);
        double den = A.dot(A);

        double delta = - (A[dir] * A[dir]) * (Axdi.dot(Axdi) - R * R * den);
        double var = den * x[i][dir] + A[dir] * Addi;

        a.resize(dim);
        if(delta > 0.) {
          bool test = false;
          for(unsigned ii = 0; ii < 2; ii++) {
            double yi = (var + ((ii == 0) ? -1 : +1) * sqrt(delta)) / den;
            if(yMin < yi && yi < yMax) { // this should exclude intersections in the vertices
              if(test == false) {
                y[cnt][dir] = yi;
                y[cnt][dirp1] = x[i][dirp1] + A[dirp1] * (y[cnt][dir] - x[i][dir]) / A[dir];
                y[cnt][dirp2] = x[i][dirp2] + A[dirp2] * (y[cnt][dir] - x[i][dir]) / A[dir];
                cnt++;
                test = true;
              }
              else { // this is to remove double intersection on the same edge
                cnt--;
              }
            }
          }
        }
      }
    }
  }
  if(cnt == 0) {
    cut = (dist2[0] > 0) ? 0 : 2;
    return;
  }
  else if(cnt > 6) { // bad: small ball
    cut = 0;
    return;
  }
  else if(cnt >= 3 && cnt <= 6) { // these I am not sure are the right limits
    cut = 1;

    std::vector<double> w(cnt + 1, 1.);
    w[cnt] = 2. * cnt;
    y.resize(cnt + 1, {0., 0., 0.});

    std::vector <double> N(dim, 0); // this is an approximate normal pointing toward the outside of the ball in the physical element
    double det = 0.;
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned i = 0; i < cnt; i++) {
        N[k] += y[i][k];
      }
      det += N[k] * N[k];
    }
    det = sqrt(det);
    for(unsigned k = 0; k < dim; k++) {
      N[k] /=  det;
    }

    double theta = atan2(N[1], N[0]);
    double phi = acos(N[2]);

    y[cnt][0] = R * sin(phi) * cos(theta);
    y[cnt][1] = R * sin(phi) * sin(theta);
    y[cnt][2] = R * cos(phi);

    femus::FindBestFit(y, w, N, a, d); // a is the unit BF normal pointing toward the outside of the ball in the physical element

    xm.resize(dim);
    xm[0] = -(d * a[0]) + xg[0];
    xm[1] = -(d * a[1]) + xg[1];
    xm[2] = -(d * a[2]) + xg[2];

    /* Differently from the triangle and the tetrahedron, a plane in the physical element does not map into a plane in the parent element.
       This below is mapping of the BF physical plane to an approximate plane in the parent element */
    
    short unsigned hex = 0;
    unsigned linear = 0;

    std::vector < std::vector < std::vector <double > > > aP(1);
    ProjectNodalToPolynomialCoefficients(aP[linear], xv, hex, linear) ;
    std::vector <double> xi;
    GetClosestPointInReferenceElement(xv, xm, hex, xi);
    bool inverseMapping = GetInverseMapping(linear, hex, aP, xm, xi, 100);

    std::vector < std::vector <double> > J, JI;
    fem->GetJacobianMatrix(xv, xi, det, J, JI);

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




//     std::vector<double> w(cnt, 1.);
//     unsigned cntp = (cnt * (cnt - 1)) / 2;
//     w.resize(cnt + cntp, (2. * cnt) / cntp);
//     y.resize(cnt + cntp);
//
//
//     std::cout << cntp<<" ";
//     unsigned cnt0 = cnt;
//     std::vector <double> N(dim, 0);
//     for(unsigned i = 0; i < cnt0 - 1; i++) {
//       for(unsigned j = i + 1; j < cnt0; j++) {
//         std::cout<< cnt <<" ";
//         std::vector <double> Nij(dim, 0); // this is an approximate normal pointing toward the outside of the ball in the physical element
//         double det = 0.;
//         for(unsigned k = 0; k < dim; k++) {
//           Nij[k] += 0.5 * (y[i][k] + y[j][k]);
//           det += Nij[k] * Nij[k];
//         }
//         det = sqrt(det);
//         for(unsigned k = 0; k < dim; k++) {
//           Nij[k] /=  det;
//           N[k] += Nij[k];
//         }
//         double theta = atan2(Nij[1], Nij[0]);
//         double phi = acos(Nij[2]);
//         y[cnt][0] = R * sin(phi) * cos(theta);
//         y[cnt][1] = R * sin(phi) * sin(theta);
//         y[cnt][2] = R * cos(phi);
//         cnt++;
//       }
//     }

