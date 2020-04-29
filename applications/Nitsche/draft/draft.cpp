#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <ctime>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include </usr/include/eigen3/Eigen/src/Core/util/DisableStupidWarnings.h>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <fstream>
#include<cmath>

// valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all ./executable
//



double a0, a1, a3, a5, a7, a9;
void SetConstants(const double &eps);
double GetDistance(const Eigen::VectorXd &x);
void GetGaussPointsWeights(unsigned &N, Eigen::VectorXd &xg, Eigen::VectorXd &wg);

void Cheb(const unsigned & m, Eigen::VectorXd &xg, Eigen::MatrixXd &C);
void GetParticlesOnBox(const double & a, const double & b, const unsigned & n1, const unsigned& dim, Eigen::MatrixXd &x, Eigen::MatrixXd &xL);
void GetParticleOnDisk(const double & a, const double & b, const unsigned & n1, const unsigned& dim, Eigen::MatrixXd &x, Eigen::MatrixXd &xL);
void AssembleMatEigen(double& a, double& b, const unsigned& m, const unsigned& dim, const unsigned& np, Eigen::Tensor<double, 3, Eigen::RowMajor>  &PmX, Eigen::MatrixXd &Pg, Eigen::VectorXd &wg, Eigen::MatrixXd &A, Eigen::VectorXd &F);
void AssembleMatEigenOnDisk(double& a, double& b, const unsigned int& m, const unsigned int& dim, const unsigned int& np,
                            Eigen::Tensor<double, 3, Eigen::RowMajor>& PmX, Eigen::MatrixXd& Pg, Eigen::VectorXd &xg, Eigen::VectorXd &wg,
                            Eigen::MatrixXd& A, Eigen::VectorXd& F);
void SolWeightEigen(Eigen::MatrixXd &A, Eigen::VectorXd &F, Eigen::VectorXd &wp, Eigen::VectorXd &w_new);
void GetChebXInfo(const unsigned& m, const unsigned& dim, const unsigned& np, Eigen::MatrixXd &xL, Eigen::Tensor<double, 3, Eigen::RowMajor>& PmX);
void Testing(double& a, double& b, const unsigned& m, const unsigned& dim, Eigen::MatrixXd &x, Eigen::VectorXd &w_new);

int main(int argc, char** args) {

  unsigned dim = 2;
  double a = 0;
  double b = 1;
  unsigned n1 = 10;// particles in one direction
  unsigned np = pow(n1, dim);
  unsigned NG = 30; // Gauss points
  unsigned m = 4; // Polynomial degree
  Eigen::MatrixXd x;
  Eigen::MatrixXd xL;
  //GetParticlesOnBox(a, b, n1, dim, x, xL);

  GetParticleOnDisk(a, b, n1, dim, x, xL);
  unsigned nq = x.row(0).size();

  Eigen::Tensor<double, 3, Eigen::RowMajor> PmX;
  //GetChebXInfo(m, dim, np, xL, PmX);
  GetChebXInfo(m, dim, nq, xL, PmX);

  Eigen::VectorXd xg;
  Eigen::VectorXd wg;
  GetGaussPointsWeights(NG, xg, wg);

  Eigen::MatrixXd Pg;
  Cheb(m, xg, Pg);

  Eigen::MatrixXd A;
  Eigen::VectorXd F;
  //AssembleMatEigen(a, b, m, dim, np, PmX, Pg, wg, A, F);
  AssembleMatEigenOnDisk(a, b, m, dim, nq, PmX, Pg, xg, wg, A, F);


  //Eigen::VectorXd wp(np);
  Eigen::VectorXd wp(nq);
  wp.fill(pow(b - a, dim) / np);

  Eigen::VectorXd w_new;
  SolWeightEigen(A, F, wp, w_new);

  Testing(a, b, m, dim, x, w_new);


  return 0;
}


void GetChebXInfo(const unsigned& m, const unsigned& dim, const unsigned& np, Eigen::MatrixXd &xL, Eigen::Tensor<double, 3, Eigen::RowMajor>& PmX) {

  PmX.resize(dim, m + 1, np);
  Eigen::MatrixXd Ptemp;
  Eigen::VectorXd xtemp;
  for(unsigned k = 0; k < dim; k++) {
    xtemp = xL.row(k);
    Cheb(m, xtemp, Ptemp);
    for(unsigned i = 0; i < m + 1; i++) {
      for(unsigned j = 0; j < np; j++) {
        PmX(k, i, j) = Ptemp(i, j);
      }
    }
  }
}


void AssembleMatEigen(double& a, double& b, const unsigned& m, const unsigned& dim, const unsigned& np, Eigen::Tensor<double, 3, Eigen::RowMajor>  &PmX, Eigen::MatrixXd &Pg, Eigen::VectorXd &wg, Eigen::MatrixXd &A, Eigen::VectorXd &F) {

  A.resize(pow(m + 1, dim), np);
  F.resize(pow(m + 1, dim));
  Eigen::VectorXi I(dim);
  Eigen::VectorXi N(dim);


  for(unsigned k = 0; k < dim ; k++) {
    N(k) = pow(m + 1, dim - k - 1);
  }

  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    I(0) = t / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = t % N(k - 1);
      I(k) = pk / N(k); // dimensional index over on the space of polynomaials
    }
    for(unsigned j = 0; j < np; j++) {
      double r = 1;
      for(unsigned k = 0; k < dim; k++) {
        r = r * PmX(k, I[k], j);
      }
      A(t, j) = r ;
    }

  }
  unsigned ng = Pg.row(0).size();
  Eigen::VectorXi J(dim);
  Eigen::VectorXi NG(dim);
  Eigen::VectorXd y(dim);


  for(unsigned k = 0; k < dim ; k++) {
    NG(k) = pow(ng, dim - k - 1);
  }

  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    I(0) = t / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = t % N(k - 1);
      I(k) = pk / N(k); // dimensional index over on the space of polynomaials
    }
    F(t) = 0.;
    for(unsigned g = 0; g < pow(ng, dim) ; g++) { // multidimensional index on the space of polynomaials
      J(0) = g / NG(0);
      for(unsigned k = 1; k < dim ; k++) {
        unsigned pk = g % NG(k - 1);
        J(k) = pk / NG(k); // dimensional index over on the space of polynomaials
      }
      double value = 1.;
      for(unsigned k = 0; k < dim ; k++) {
        value *= (b - a) / 2. * Pg(I(k), J(k)) * wg(J(k)) ;
      }
      F(t) += value;
    }
  }

}

void AssembleMatEigenOnDisk(double& a, double& b, const unsigned int& m, const unsigned int& dim, const unsigned int& np,
                            Eigen::Tensor<double, 3, Eigen::RowMajor>& PmX, Eigen::MatrixXd& Pg, Eigen::VectorXd &xg, Eigen::VectorXd &wg,
                            Eigen::MatrixXd& A, Eigen::VectorXd& F) {

  A.resize(pow(m + 1, dim), np);
  F.resize(pow(m + 1, dim));
  Eigen::VectorXi I(dim);
  Eigen::VectorXi N(dim);


  for(unsigned k = 0; k < dim ; k++) {
    N(k) = pow(m + 1, dim - k - 1);
  }

  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    I(0) = t / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = t % N(k - 1);
      I(k) = pk / N(k); // dimensional index over on the space of polynomaials
    }
    for(unsigned j = 0; j < np; j++) {
      double r = 1;
      for(unsigned k = 0; k < dim; k++) {
        r = r * PmX(k, I[k], j);
      }
      A(t, j) = r ;
    }

  }

  double a0;
  double a1;
  double a3;
  double a5;
  double a7;
  double a9;
  double eps = (b - a) / 10.;

  a0 = 0.5; // 128./256.;
  a1 = pow(eps, -1.) * 1.23046875; // 315/256.;
  a3 = -pow(eps, -3.) * 1.640625; //420./256.;
  a5 = pow(eps, -5.) * 1.4765625; // 378./256.;
  a7 = -pow(eps, -7.) * 0.703125; // 180./256.;
  a9 = pow(eps, -9.) * 0.13671875; // 35./256.;

  unsigned ng = Pg.row(0).size();
  Eigen::VectorXi J(dim);
  Eigen::VectorXi NG(dim);
  Eigen::VectorXd y(dim);


  for(unsigned k = 0; k < dim ; k++) {
    NG(k) = pow(ng, dim - k - 1);
  }

  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    I(0) = t / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = t % N(k - 1);
      I(k) = pk / N(k); // dimensional index over on the space of polynomaials
    }
    F(t) = 0.;
    for(unsigned g = 0; g < pow(ng, dim) ; g++) { // multidimensional index on the space of polynomaials
      J(0) = g / NG(0);
      for(unsigned k = 1; k < dim ; k++) {
        unsigned pk = g % NG(k - 1);
        J(k) = pk / NG(k); // dimensional index over on the space of polynomaials
      }
      double value = 1.;
      for(unsigned k = 0; k < dim ; k++) {
        value *= (b - a) / 2. * Pg(I(k), J(k)) * wg(J(k)) ;
        y[k] = a + ((xg(J[k]) + 1.) / 2.) * (b - a) ;
      }
      double dg1 = GetDistance(y);
      double dg2 = dg1 * dg1;
      double xi;
      if(dg1 < -eps)
        xi = 0.;
      else if(dg1 > eps) {
        xi = 1.;
      }
      else {
        xi = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
      }
      F(t) += xi * value;
    }
  }

}


void SolWeightEigen(Eigen::MatrixXd &A, Eigen::VectorXd &F, Eigen::VectorXd &wp, Eigen::VectorXd &w_new) {

  Eigen::VectorXd y_temp = (A * A.transpose()).partialPivLu().solve(F - A * wp);
  Eigen::VectorXd w_new_temp = A.transpose() * y_temp;
  w_new = w_new_temp + wp;

}


void Testing(double &a, double &b, const unsigned &m, const unsigned &dim, Eigen::MatrixXd &x, Eigen::VectorXd &w_new) {

  double s = 0;
  for(unsigned i = 0; i < w_new.size(); i++) {
    double r = 1;
    for(unsigned k = 0; k < dim; k++) {
      r = r * x(k, i);
    }
    s = s  /*pow(r, m)*/ + w_new(i);
  }
  //double err = (s - pow((pow(b, m + 1) - pow(a, m + 1)) / (m + 1), dim)) / s;
  std::cout << s << " " << M_PI * 0.75 * 0.75 / 4. << std::endl;
  double err = (s - M_PI * 0.75 * 0.75 / 4.) / s;
  std::cout << "Error is: " << err << std::endl;

}


void Cheb(const unsigned &m, Eigen::VectorXd &xg, Eigen::MatrixXd &C) {

  C.resize(xg.size(), m + 1);
  for(unsigned i = 0; i < xg.size(); i++) {
    C(i, 0) = 1;
    C(i, 1) = xg(i);
    for(unsigned j = 2; j <= m; j++) {
      C(i, j) =  2 * xg(i) * C(i, j - 1) - C(i, j - 2);
    }
  }
  C.transposeInPlace();

}

void  GetParticlesOnBox(const double &a, const double &b, const unsigned &n1, const unsigned &dim, Eigen::MatrixXd &x, Eigen::MatrixXd &xL) {
  double h = (b - a) / n1;
  x.resize(dim, pow(n1, dim));
  Eigen::VectorXi I(dim);
  Eigen::VectorXi N(dim);

  for(unsigned k = 0; k < dim ; k++) {
    N(k) = pow(n1, dim - k - 1);
  }

  for(unsigned p = 0; p < pow(n1, dim) ; p++) {
    I(0) = 1 + p / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = p % N(k - 1);
      I(k) = 1 + pk / N(k);
    }
    //std::cout << I(0) << " " << I(1) << std::endl;

    for(unsigned k = 0; k < dim ; k++) {
      std::srand(std::time(0));
      double r = 2 * ((double) rand() / (RAND_MAX)) - 1;
      x(k, p) = a + h / 2 + (I(k) - 1) * h; // + 0.1 * r;
    }
  }

  xL.resize(dim, pow(n1, dim));
  Eigen::MatrixXd ID;
  ID.resize(dim, pow(n1, dim));
  ID.fill(1.);
  xL = (2. / (b - a)) * x - ((b + a) / (b - a)) * ID;


}


void GetGaussPointsWeights(unsigned &N, Eigen::VectorXd &xg, Eigen::VectorXd &wg) {
  unsigned N1 = N ;
  unsigned N2 = N + 1;
  Eigen::VectorXd xu;
  xu.setLinSpaced(N1, -1, 1);
  xg.resize(N1);
  for(unsigned i = 0; i <= N - 1; i++) {
    xg(i) = cos((2 * i + 1) * M_PI / (2 * (N - 1) + 2)) + (0.27 / N1) * sin(M_PI * xu(i) * (N - 1) / N2) ;
  }
  Eigen::MatrixXd L(N1, N2);
  L.fill(0.);
  Eigen::VectorXd Lp(N1);
  Lp.fill(0.);
  Eigen::VectorXd y0(xg.size());
  y0.fill(2);
  double eps = 1e-15;
  Eigen::VectorXd d = xg - y0;

  double max = d.cwiseAbs().maxCoeff();

  while(max > eps) {

    L.col(0).fill(1.);
    L.col(1) = xg;

    for(unsigned k = 2; k < N2; k++) {
      for(unsigned i = 0; i < N1; i++) {
        L(i, k) = ((2 * k - 1) * xg(i) * L(i, k - 1) - (k - 1) * L(i, k - 2)) / k;
      }
    }


    for(unsigned i = 0; i < N1; i++) {
      Lp(i) = N2 * (L(i, N1 - 1) - xg(i) * L(i, N2 - 1)) / (1 - xg(i) * xg(i));
    }

    y0 = xg;

    for(unsigned i = 0; i < N1; i++) {
      xg(i) =  y0(i) - L(i, N2 - 1) / Lp(i);
    }

    d = xg - y0;

    max = d.cwiseAbs().maxCoeff();
  }

  wg.resize(N1);

  for(unsigned i = 0; i < N1; i++) {
    double r = double(N2) / double(N1);
    wg(i) = (2) / ((1 - xg(i) * xg(i)) * Lp(i) * Lp(i)) * r * r;
  }
}




double GetDistance(const Eigen::VectorXd &x) {

  double radius = 0.75;
  Eigen::VectorXd xc(x.size());
  xc.fill(0.);

  double rx = 0;
  for(unsigned i = 0; i < x.size(); i++) {
    rx += (x[i] - xc[i]) * (x[i] - xc[i]);
  }
  return radius - sqrt(rx);

}

void  GetParticleOnDisk(const double &a, const double &b, const unsigned &n1, const unsigned &dim, Eigen::MatrixXd &x, Eigen::MatrixXd &xL) {
  double h = (b - a) / n1;
  x.resize(dim, pow(n1, dim));
  Eigen::VectorXi I(dim);
  Eigen::VectorXi N(dim);

  for(unsigned k = 0; k < dim ; k++) {
    N(k) = pow(n1, dim - k - 1);
  }
  double q = 0;

  Eigen::VectorXd y(dim);

  for(unsigned p = 0; p < pow(n1, dim) ; p++) {
    I(0) = 1 + p / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = p % N(k - 1);
      I(k) = 1 + pk / N(k);
    }

    for(unsigned k = 0; k < dim ; k++) {

      //std::srand(std::time(0));
      //double r = 2 * ((double) rand() / (RAND_MAX)) - 1;

      y[k] = a + h / 2 + (I(k) - 1) * h; // + 0.1 * r;

    }


    if(GetDistance(y) > 0) {
      x.col(q) = y;
//       for(unsigned k = 0; k < dim ; k++) {
//           //x(k,q) = y[k];
//       }
      q++;
    }


  }

  x.conservativeResize(dim, q);

  xL.resize(dim, pow(n1, dim));
  Eigen::MatrixXd ID;
  ID.resize(dim, pow(n1, dim));
  ID.fill(1.);
  xL = (2. / (b - a)) * x - ((b + a) / (b - a)) * ID;

}


void SetConstants(const double & eps) {
  a0 = 0.5; // 128./256.;
  a1 = pow(eps, -1.) * 1.23046875; // 315/256.;
  a3 = -pow(eps, -3.) * 1.640625; //420./256.;
  a5 = pow(eps, -5.) * 1.4765625; // 378./256.;
  a7 = -pow(eps, -7.) * 0.703125; // 180./256.;
  a9 = pow(eps, -9.) * 0.13671875; // 35./256.;
}






















