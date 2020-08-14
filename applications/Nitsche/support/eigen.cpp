

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <ctime>
#include <fstream>
#include <cmath>
#include "Marker.hpp"
#include "Line.hpp"
#include <algorithm>

#include "eigen.hpp"
#include "sharedFunctions.hpp"
 
void Testing(double & a, double & b, const unsigned & m, const unsigned & dim, Eigen::MatrixXd & x,
             Eigen::VectorXd & w_new, std::vector<double> &dist, const double & eps, double & QuadSum, double & IntSum) {

  
  double deps = (b - a) * eps; // eps1

  double a0 = 0.5; // 128./256.;
  double a1 = pow(deps, -1.) * 1.23046875; // 315/256.;
  double a3 = -pow(deps, -3.) * 1.640625; //420./256.;
  double a5 = pow(deps, -5.) * 1.4765625; // 378./256.;
  double a7 = -pow(deps, -7.) * 0.703125; // 180./256.;
  double a9 = pow(deps, -9.) * 0.13671875; // 35./256.;

  QuadSum = 0.;
  IntSum = 0.;

  for(unsigned i = 0; i < w_new.size(); i++) {

    double dg1 = dist[i];
    double dg2 = dg1 * dg1;
    double xi;
    if(dg1 < -deps)
      xi = 0.;
    else if(dg1 > deps) {
      xi = 1.;
    }
    else {
      xi = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
    }

    double r = 0.75 * 0.75 ;
    for(unsigned k = 0; k < dim; k++) {
      r +=  - x(k, i) * x(k, i);
    }
    IntSum += xi * r * w_new(i);

    r = 1.;
    for(unsigned k = 0; k < dim; k++) {
      r *=  x(k, i);
    }
    QuadSum += pow(r, m) * w_new(i);
  }

}


void Cheb(const unsigned & m, Eigen::VectorXd & xg, Eigen::MatrixXd & C) {

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




void GetChebGaussF(const unsigned & dim, const unsigned & m, std::vector<double> &VxL, std::vector<double> &VxU, Eigen::MatrixXd & Pg,  Eigen::VectorXd & wg, Eigen::VectorXd & F) {

  F.resize(pow(m + 1, dim));
  F.setZero();
  Eigen::VectorXi I(dim);
  Eigen::VectorXi N(dim);
  Eigen::VectorXi J(dim);
  Eigen::VectorXi NG(dim);

  unsigned ng = Pg.row(0).size();

  for(unsigned k = 0; k < dim ; k++) {
    N(k) = pow(m + 1, dim - k - 1);
  }

  for(unsigned k = 0; k < dim ; k++) {
    NG(k) = pow(ng, dim - k - 1);
  }

  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    I(0) = t / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = t % N(k - 1);
      I(k) = pk / N(k);
    }
    F(t) = 0.;
    for(unsigned g = 0; g < pow(ng, dim) ; g++) { // gauss loop
      J(0) = g / NG(0);
      for(unsigned k = 1; k < dim ; k++) {
        unsigned pk = g % NG(k - 1);
        J(k) = pk / NG(k);
      }
      double value = 1.;
      unsigned ig = 0;
      double jac = 1.;

      for(unsigned k = 0; k < dim ; k++) {
        value *= 0.5 * (VxU[k] - VxL[k]) * Pg(I(k), J(k)) * wg(J(k));
      }
      F(t) += value;
    }
  }

}


void GetChebGaussF(const unsigned & dim, const unsigned & m, const std::vector<double> &jac, Eigen::MatrixXd & Pg,  Eigen::VectorXd & wg, Eigen::VectorXd & F) {

  F.resize(pow(m + 1, dim));
  F.setZero();
  Eigen::VectorXi I(dim);
  Eigen::VectorXi N(dim);
  Eigen::VectorXi J(dim);
  Eigen::VectorXi NG(dim);

  unsigned ng = Pg.row(0).size();

  for(unsigned k = 0; k < dim ; k++) {
    N(k) = pow(m + 1, dim - k - 1);
  }

  for(unsigned k = 0; k < dim ; k++) {
    NG(k) = pow(ng, dim - k - 1);
  }

  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    I(0) = t / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = t % N(k - 1);
      I(k) = pk / N(k);
    }
    F(t) = 0.;
    for(unsigned g = 0; g < pow(ng, dim) ; g++) { // gauss loop
      J(0) = g / NG(0);
      for(unsigned k = 1; k < dim ; k++) {
        unsigned pk = g % NG(k - 1);
        J(k) = pk / NG(k);
      }

      double value = jac[g];
      for(unsigned k = 0; k < dim ; k++) {
        value *= Pg(I(k), J(k)) * wg(J(k));
      }

      F(t) += value;
    }
  }

}



void GetChebXInfo(const unsigned & m, const unsigned & dim, const unsigned & np, Eigen::MatrixXd & xL, Eigen::Tensor<double, 3, Eigen::RowMajor>& PmX) {
  // xL is taken in reference coordinate system
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



void GetMultiDimChebMatrix(const unsigned & dim, const unsigned & m, const unsigned & np, Eigen::Tensor<double, 3, Eigen::RowMajor>  &PmX, Eigen::MatrixXd & A) {


  A.resize(pow(m + 1, dim), np);
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
        r *= PmX(k, I[k], j);
      }
      A(t, j) = r ;
    }
  }


}




void SolWeightEigen(Eigen::MatrixXd & A, Eigen::VectorXd & F, Eigen::VectorXd & wP, Eigen::VectorXd & w_new) {

  w_new.resize(wP.size());

  Eigen::VectorXd y = A.transpose() * (A * A.transpose()).partialPivLu().solve(F - A * wP);
  w_new = y + wP;


}


void PrintMarkers(const unsigned & dim, const Eigen::MatrixXd & xP, const std::vector <double> &dist,
                  const Eigen::VectorXd wP, const Eigen::VectorXd & w_new, const unsigned & l, const unsigned & t) {

  std::ofstream fout;

  char filename[100];
  sprintf(filename, "marker%d.txt", l);

  if(t == 0) {
    fout.open(filename);
  }
  else {
    fout.open(filename, std::ios::app);
  }

  for(unsigned i = 0; i < w_new.size(); i++) {

    for(unsigned k = 0; k < dim; k++) {
      fout << xP(k, i) << " ";
    }
    fout <<  dist[i] << " " << wP[i] << " " << w_new[i] << std::endl;
  }

  fout.close();

}



// N-point gauss quadrature points and weights by finding the roots of Legendre polynomaial
void GetGaussPointsWeights(unsigned & N, Eigen::VectorXd & xg, Eigen::VectorXd & wg) {
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

  //Newton step for finding the roots
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

  // compute the weights from the roots
  wg.resize(N1);
  for(unsigned i = 0; i < N1; i++) {
    double r = double(N2) / double(N1);
    wg(i) = (2) / ((1 - xg(i) * xg(i)) * Lp(i) * Lp(i)) * r * r;
  }
}


double GetDistance(const Eigen::VectorXd & x) {

  double radius = 0.75;
  Eigen::VectorXd xc(x.size());
  xc.fill(0.);

  double rx = 0;
  for(unsigned i = 0; i < x.size(); i++) {
    rx += (x[i] - xc[i]) * (x[i] - xc[i]);
  }
  return radius - sqrt(rx);

}


double get_g(const double & r, const double & T, const unsigned & n) {
  double rn = pow(r, n);
  return (-1. + r) * (-1 + r * rn + T - r * T) / (1. + (-1. + n * (-1 + r)) * rn);
}

double get_r(const double & T, const unsigned & n) {
  double r0 = 2.;
  double r = 0;
  while(fabs(r - r0) > 1.0e-10) {
    r0 = r;
    r = r0 - get_g(r0, T, n);
  }
  return r;
}



void GetInterfaceElementEigenvalues(MultiLevelSolution& mlSol, Line* line3, Line* lineI, const double &deps) {

  double a0 = 0.5; // 128./256.;
  double a1 = pow(deps, -1.) * 1.23046875; // 315/256.;
  double a3 = -pow(deps, -3.) * 1.640625; //420./256.;
  double a5 = pow(deps, -5.) * 1.4765625; // 378./256.;
  double a7 = -pow(deps, -7.) * 0.703125; // 180./256.;
  double a9 = pow(deps, -9.) * 0.13671875; // 35./256.;
  
    
    
  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel(level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level);
  unsigned iproc  = msh->processor_id();

  unsigned eflagIndex = mlSol.GetIndex("eflag");

  unsigned CMIndex[2];
  CMIndex[0] = mlSol.GetIndex("CM1");
  CMIndex[1] = mlSol.GetIndex("CM2");

  unsigned CLIndex[2];
  CLIndex[0] = mlSol.GetIndex("CL1");
  CLIndex[1] = mlSol.GetIndex("CL2");

  unsigned solIndex = mlSol.GetIndex("VX1");
  unsigned soluType = mlSol.GetSolutionType(solIndex);

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  std::vector < std::vector<double> >  x(dim);

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector < std::vector < std::vector <double > > > aP(3);

  std::vector<Marker*> particle3 = line3->GetParticles();
  std::vector<unsigned> markerOffset3 = line3->GetMarkerOffset();
  unsigned imarker3 = markerOffset3[iproc];


  std::vector<Marker*> particleI = lineI->GetParticles();
  std::vector<unsigned> markerOffsetI = lineI->GetMarkerOffset();
  unsigned imarkerI = markerOffsetI[iproc];

  sol->_Sol[CMIndex[0]]->zero();
  sol->_Sol[CMIndex[1]]->zero();

  sol->_Sol[CLIndex[0]]->zero();
  sol->_Sol[CLIndex[1]]->zero();


  Eigen::MatrixXd AM;
  Eigen::MatrixXd AL;
  Eigen::MatrixXd BM[2];
  Eigen::MatrixXd BL[2];


  clock_t eigenTime = 0;

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));
    if(eFlag == 1) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);  // number of solution element dofs

      unsigned sizeAll = dim * nDofu;

      AM.resize(sizeAll, sizeAll);
      AM.setZero();

      AL.resize(sizeAll, sizeAll);
      AL.setZero();

      for(unsigned k = 0; k < 2; k++) {
        BM[k].resize(sizeAll, sizeAll);
        BM[k].setZero();

        BL[k].resize(sizeAll, sizeAll);
        BL[k].setZero();
      }

//       aM.assign(sizeAll * sizeAll, 0.);
//       bM[0].assign(sizeAll * sizeAll, 0.);
//       bM[1].assign(sizeAll * sizeAll, 0.);
//
//       aL.assign(sizeAll * sizeAll, 0.);
//       bL[0].assign(sizeAll * sizeAll, 0.);
//       bL[1].assign(sizeAll * sizeAll, 0.);

      for(int k = 0; k < dim; k++) {
        x[k].resize(nDofu);
      }

      for(unsigned i = 0; i < nDofu; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof
        for(unsigned k = 0; k < dim; k++) {
          x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
        }
      }

      while(imarker3 < markerOffset3[iproc + 1] && iel > particle3[imarker3]->GetMarkerElement()) {
        imarker3++;
      }
      while(imarker3 < markerOffset3[iproc + 1] && iel == particle3[imarker3]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle3[imarker3]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle3[imarker3]->GetMarkerMass();

        double dg1 = particle3[imarker3]->GetMarkerDistance();

        double dg2 = dg1 * dg1;
        double chi;
        if(dg1 < -deps)
          chi = 0.;
        else if(dg1 > deps) {
          chi = 1.;
        }
        else {
          chi = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
        }

        for(unsigned k = 0; k < dim; k++) {
          for(unsigned i = 0; i < nDofu; i++) {
            for(unsigned l = 0; l < dim; l++) {
              for(unsigned j = 0; j < nDofu; j++) {
                BM[0](nDofu * k + i, k * nDofu + j) += (1. - chi) * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + l] * weight;
                BM[0](nDofu * k + i, l * nDofu + j) += (1. - chi) * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + k] * weight;

                BL[0](nDofu * k + i, l * nDofu + j) += (1. - chi) * phi_x[i * dim + k] * phi_x[j * dim + l] * weight;

                BM[1](nDofu * k + i, k * nDofu + j) += chi * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + l] * weight;
                BM[1](nDofu * k + i, l * nDofu + j) += chi * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + k] * weight;

                BL[1](nDofu * k + i, l * nDofu + j) += chi * phi_x[i * dim + k] * phi_x[j * dim + l] * weight;

              }
            }
          }
        }
        imarker3++;
      }



      // interface
      while(imarkerI < markerOffsetI[iproc + 1] && iel > particleI[imarkerI]->GetMarkerElement()) {
        imarkerI++;
      }
      while(imarkerI < markerOffsetI[iproc + 1] && iel == particleI[imarkerI]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particleI[imarkerI]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);

        std::vector <std::vector < double > > T;
        particleI[imarkerI]->GetMarkerTangent(T);

        double weight;
        std::vector < double > N(dim);
        if(dim == 2) {
          N[0] =  T[0][1];
          N[1] = -T[0][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1]);
          N[0] /= weight;
          N[1] /= weight;
        }
        else {
          N[0] = T[0][1] * T[1][2] - T[0][2] * T[1][1];
          N[1] = T[0][2] * T[1][0] - T[0][0] * T[1][2];
          N[2] = T[0][0] * T[1][1] - T[0][1] * T[1][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
          N[0] /= weight;
          N[1] /= weight;
          N[2] /= weight;
        }

        // *** phi_i loop ***

        for(unsigned k = 0; k < dim; k++) {
          for(int i = 0; i < nDofu; i++) {

            double gradPhiiDotN = 0.;
            for(unsigned l = 0; l < dim; l++) {
              gradPhiiDotN += phi_x[i * dim + l] * N[l];
            }
            for(int j = 0; j < nDofu; j++) {
              for(unsigned l = 0; l < dim; l++) {

                AM(nDofu * k + i, k * nDofu + j) += 0.5 * gradPhiiDotN * 0.5 * N[l] *  phi_x[j * dim + l]  * weight;
                AM(nDofu * k + i, l * nDofu + j) += 0.5 * gradPhiiDotN * 0.5 * N[l] *  phi_x[j * dim + k]  * weight;

                AL(nDofu * k + i, l * nDofu + j) += phi_x[i * dim + k] * phi_x[j * dim + l] * weight;

              }
              for(unsigned l1 = 0; l1 < dim; l1++) {
                for(unsigned l2 = 0; l2 < dim; l2++) {
                  AM(nDofu * k + i, l1 * nDofu + j) += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  phi_x[j * dim + l2]  * weight;
                  AM(nDofu * k + i, l2 * nDofu + j) += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  phi_x[j * dim + l1]  * weight;
                }
              }

            }
          } // end phi_i loop
        }
        imarkerI++;
      }

      double perturbation = 1e-10;

      std::cout << "======================EIGEN===================================" << std::endl;

      clock_t start = clock();

      Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
      double inf = 1e+10;

      for(unsigned k = 0; k < 2; k++) {
        double BM0Lk = BM[k].norm();

        for(unsigned i = 0; i < sizeAll; i++) {
          BM[k](i, i) += perturbation * BM0Lk;
        }

        ges.compute(AM, BM[k], false);
        std::complex < double > temp;
        Eigen::VectorXcd eig;

        eig = ges.eigenvalues();
        double emax0 = 0.;
        for(unsigned i = 0; i < sizeAll; i++) {
          temp = eig(i);
          if(fabs(real(temp)) > emax0 && fabs(real(temp)) < inf) {
            emax0 = fabs(real(temp));
          }
        }
        std::cout << iel << " " << emax0 << std::endl;
        sol->_Sol[CMIndex[k]]->set(iel, emax0);
      }

      for(unsigned k = 0; k < 2; k++) {
        double norm = BL[k].norm();

        for(unsigned i = 0; i < sizeAll; i++) {
          BL[k](i, i) += perturbation * norm;
        }

        ges.compute(AL, BL[k], false);
        std::complex < double > temp;
        Eigen::VectorXcd eig;

        eig = ges.eigenvalues();
        double emax0 = 0.;
        for(unsigned i = 0; i < sizeAll; i++) {
          temp = eig(i);
          if(fabs(real(temp)) > emax0 && fabs(real(temp)) < inf) {
            emax0 = fabs(real(temp));
          }
        }
        std::cout << iel << " " << emax0 << std::endl;
        sol->_Sol[CLIndex[k]]->set(iel, emax0);
      }
      eigenTime += (clock() - start);
    } // end of eflag loop
  } //end of element loop

  sol->_Sol[CMIndex[0]]->close();
  sol->_Sol[CMIndex[1]]->close();

  sol->_Sol[CLIndex[0]]->close();
  sol->_Sol[CLIndex[1]]->close();

  //std::cout << std::endl << "petsc TIME:\t" << static_cast<double>(petscTime) / CLOCKS_PER_SEC << std::endl;
  std::cout << std::endl << "Eigen TIME:\t" << static_cast<double>(eigenTime) / CLOCKS_PER_SEC << std::endl;

}



void GetInterfaceElementEigenvaluesAD(MultiLevelSolution& mlSol, Line* line3, Line* lineI, const double &deps) {

  double a0 = 0.5; // 128./256.;
  double a1 = pow(deps, -1.) * 1.23046875; // 315/256.;
  double a3 = -pow(deps, -3.) * 1.640625; //420./256.;
  double a5 = pow(deps, -5.) * 1.4765625; // 378./256.;
  double a7 = -pow(deps, -7.) * 0.703125; // 180./256.;
  double a9 = pow(deps, -9.) * 0.13671875; // 35./256.;

  adept::Stack& s = FemusInit::_adeptStack;

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel(level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level);
  unsigned iproc  = msh->processor_id();

  unsigned eflagIndex = mlSol.GetIndex("eflag");

  unsigned CMIndex[2];
  CMIndex[0] = mlSol.GetIndex("CM1");
  CMIndex[1] = mlSol.GetIndex("CM2");

  unsigned CLIndex[2];
  CLIndex[0] = mlSol.GetIndex("CL1");
  CLIndex[1] = mlSol.GetIndex("CL2");

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  std::vector < unsigned > solIndex(dim);
  solIndex[0] = mlSol.GetIndex("VX1");
  solIndex[1] = mlSol.GetIndex("VX2");
  if(dim == 3) {
    solIndex[2] = mlSol.GetIndex("VX3");
  }
  unsigned soluType = mlSol.GetSolutionType(solIndex[0]);


  std::vector < std::vector<adept::adouble> >  solV(dim);
  std::vector < std::vector<double> >  x(dim);

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector < std::vector < std::vector <double > > > aP(3);

  std::vector<Marker*> particle3 = line3->GetParticles();
  std::vector<unsigned> markerOffset3 = line3->GetMarkerOffset();
  unsigned imarker3 = markerOffset3[iproc];


  std::vector<Marker*> particleI = lineI->GetParticles();
  std::vector<unsigned> markerOffsetI = lineI->GetMarkerOffset();
  unsigned imarkerI = markerOffsetI[iproc];

  sol->_Sol[CMIndex[0]]->zero();
  sol->_Sol[CMIndex[1]]->zero();

  sol->_Sol[CLIndex[0]]->zero();
  sol->_Sol[CLIndex[1]]->zero();


  Eigen::MatrixXd AM;
  Eigen::MatrixXd AL;
  Eigen::MatrixXd BM[2];
  Eigen::MatrixXd BL[2];


  std::vector< adept::adouble > resAM;
  std::vector< adept::adouble > resAL;;
  std::vector< adept::adouble > resBM[2];
  std::vector< adept::adouble > resBL[2];

  std::vector< double > JacAM;
  std::vector< double > JacAL;;
  std::vector< double > JacBM[2];
  std::vector< double > JacBL[2];


  clock_t eigenTime = 0;

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));
    if(eFlag == 1) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);  // number of solution element dofs

      const unsigned sizeAll = dim * nDofu;

      AM.resize(sizeAll, sizeAll);
      AM.setZero();

      AL.resize(sizeAll, sizeAll);
      AL.setZero();

      for(unsigned k = 0; k < 2; k++) {
        BM[k].resize(sizeAll, sizeAll);
        BM[k].setZero();

        BL[k].resize(sizeAll, sizeAll);
        BL[k].setZero();
      }

      for(int k = 0; k < dim; k++) {
        x[k].resize(nDofu);
        solV[k].resize(nDofu);
      }

      resAM.assign(sizeAll, 0.);
      resAL.assign(sizeAll, 0.);
      resBM[0].assign(sizeAll, 0.);
      resBM[1].assign(sizeAll, 0.);
      resBL[0].assign(sizeAll, 0.);
      resBL[1].assign(sizeAll, 0.);

      for(unsigned i = 0; i < nDofu; i++) {

        unsigned uDof  = msh->GetSolutionDof(i, iel, soluType);
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
        for(unsigned k = 0; k < dim; k++) {
          solV[k][i] = (*sol->_Sol[solIndex[k]])(uDof);
          x[k][i] = (*msh->_topology->_Sol[k])(xDof);
        }
      }

      s.new_recording();

      while(imarker3 < markerOffset3[iproc + 1] && iel > particle3[imarker3]->GetMarkerElement()) {
        imarker3++;
      }
      while(imarker3 < markerOffset3[iproc + 1] && iel == particle3[imarker3]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle3[imarker3]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle3[imarker3]->GetMarkerMass();

        std::vector<adept::adouble> solVg(dim, 0.);
        std::vector < std::vector<adept::adouble> > gradSolVg(dim);
        for(unsigned k = 0; k < dim; k++) {
          gradSolVg[k].assign(dim, 0.);
          for(unsigned i = 0; i < nDofu; i++) {
            solVg[k] += phi[i] * solV[k][i];
            for(unsigned l = 0; l < dim; l++) {
              gradSolVg[k][l] += phi_x[i * dim + l] * solV[k][i];
            }
          }
        }

        double dg1 = particle3[imarker3]->GetMarkerDistance();

        double dg2 = dg1 * dg1;
        double chi;
        if(dg1 < -deps)
          chi = 0.;
        else if(dg1 > deps) {
          chi = 1.;
        }
        else {
          chi = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
        }

        for(unsigned k = 0; k < dim; k++) {
          for(unsigned i = 0; i < nDofu; i++) {
            for(unsigned l = 0; l < dim; l++) {

              resBM[0][nDofu * k + i] += (1. - chi) * 0.5 * phi_x[i * dim + l] * gradSolVg[k][l] * weight;
              resBM[0][nDofu * k + i] += (1. - chi) * 0.5 * phi_x[i * dim + l] * gradSolVg[l][k] * weight;
              resBL[0][nDofu * k + i] += (1. - chi) * phi_x[i * dim + k] * gradSolVg[l][l] * weight;

              resBM[1][nDofu * k + i] += chi * 0.5 * phi_x[i * dim + l] * gradSolVg[k][l]  * weight;
              resBM[1][nDofu * k + i] += chi * 0.5 * phi_x[i * dim + l]  * gradSolVg[l][k] * weight;
              resBL[1][nDofu * k + i] += chi * phi_x[i * dim + k] * gradSolVg[l][l] * weight;

//               for(unsigned j = 0; j < nDofu; j++) {
//                 BM[0](nDofu * k + i, k * nDofu + j) += (1. - chi) * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + l] * weight;
//                 BM[0](nDofu * k + i, l * nDofu + j) += (1. - chi) * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + k] * weight;
//
//                 BL[0](nDofu * k + i, l * nDofu + j) += (1. - chi) * phi_x[i * dim + k] * phi_x[j * dim + l] * weight;
//
//                 BM[1](nDofu * k + i, k * nDofu + j) += chi * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + l] * weight;
//                 BM[1](nDofu * k + i, l * nDofu + j) += chi * 0.5 * phi_x[i * dim + l] * phi_x[j * dim + k] * weight;
//
//                 BL[1](nDofu * k + i, l * nDofu + j) += chi * phi_x[i * dim + k] * phi_x[j * dim + l] * weight;
//
//               }
            }
          }
        }
        imarker3++;
      }



      // interface
      while(imarkerI < markerOffsetI[iproc + 1] && iel > particleI[imarkerI]->GetMarkerElement()) {
        imarkerI++;
      }
      while(imarkerI < markerOffsetI[iproc + 1] && iel == particleI[imarkerI]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particleI[imarkerI]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);

        std::vector<adept::adouble> solVg(dim, 0.);
        std::vector < std::vector<adept::adouble> > gradSolVg(dim);
        for(unsigned k = 0; k < dim; k++) {
          gradSolVg[k].assign(dim, 0.);
          for(unsigned i = 0; i < nDofu; i++) {
            solVg[k] += phi[i] * solV[k][i];
            for(unsigned l = 0; l < dim; l++) {
              gradSolVg[k][l] += phi_x[i * dim + l] * solV[k][i];
            }
          }
        }

        std::vector <std::vector < double > > T;
        particleI[imarkerI]->GetMarkerTangent(T);

        double weight;
        std::vector < double > N(dim);
        if(dim == 2) {
          N[0] =  T[0][1];
          N[1] = -T[0][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1]);
          N[0] /= weight;
          N[1] /= weight;
        }
        else {
          N[0] = T[0][1] * T[1][2] - T[0][2] * T[1][1];
          N[1] = T[0][2] * T[1][0] - T[0][0] * T[1][2];
          N[2] = T[0][0] * T[1][1] - T[0][1] * T[1][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
          N[0] /= weight;
          N[1] /= weight;
          N[2] /= weight;
        }

        // *** phi_i loop ***

        for(unsigned k = 0; k < dim; k++) {
          for(int i = 0; i < nDofu; i++) {

            double gradPhiiDotN = 0.;
            for(unsigned l = 0; l < dim; l++) {
              gradPhiiDotN += phi_x[i * dim + l] * N[l];
            }

            for(unsigned l = 0; l < dim; l++) {

              resAM[nDofu * k + i] += 0.5 * gradPhiiDotN * 0.5 * N[l] *  gradSolVg[k][l]  * weight;
              resAM[nDofu * k + i] += 0.5 * gradPhiiDotN * 0.5 * N[l] *  gradSolVg[l][k]  * weight;
              resAL[nDofu * k + i] += phi_x[i * dim + k] * gradSolVg[l][l] * weight;
            }
            for(unsigned l1 = 0; l1 < dim; l1++) {
              for(unsigned l2 = 0; l2 < dim; l2++) {
                resAM[nDofu * k + i] += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  gradSolVg[l1][l2] * weight;
                resAM[nDofu * k + i] += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  gradSolVg[l2][l1] * weight;
              }
            }

//             for(int j = 0; j < nDofu; j++) {
//               for(unsigned l = 0; l < dim; l++) {
//
//                 AM(nDofu * k + i, k * nDofu + j) += 0.5 * gradPhiiDotN * 0.5 * N[l] *  phi_x[j * dim + l]  * weight;
//                 AM(nDofu * k + i, l * nDofu + j) += 0.5 * gradPhiiDotN * 0.5 * N[l] *  phi_x[j * dim + k]  * weight;
//
//                 AL(nDofu * k + i, l * nDofu + j) += phi_x[i * dim + k] * phi_x[j * dim + l] * weight;
//
//               }
//               for(unsigned l1 = 0; l1 < dim; l1++) {
//                 for(unsigned l2 = 0; l2 < dim; l2++) {
//                   AM(nDofu * k + i, l1 * nDofu + j) += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  phi_x[j * dim + l2]  * weight;
//                   AM(nDofu * k + i, l2 * nDofu + j) += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  phi_x[j * dim + l1]  * weight;
//                 }
//               }
//             }
          } // end phi_i loop
        }
        imarkerI++;
      }



      JacBM[0].resize(sizeAll * sizeAll);
      JacBM[1].resize(sizeAll * sizeAll);
      JacBL[0].resize(sizeAll * sizeAll);
      JacBL[1].resize(sizeAll * sizeAll);


      for(unsigned  k = 0; k < dim; k++) {
        s.independent(&solV[k][0], nDofu);
      }
      s.dependent(&resAM[0], sizeAll);
      JacAM.resize(sizeAll * sizeAll);
      s.jacobian(&JacAM[0], true);
      s.clear_dependents();

      s.dependent(&resAL[0], sizeAll);
      JacAL.resize(sizeAll * sizeAll);
      s.jacobian(&JacAL[0], true);
      s.clear_dependents();

      s.dependent(&resBM[0][0], sizeAll);
      JacBM[0].resize(sizeAll * sizeAll);
      s.jacobian(&JacBM[0][0], true);
      s.clear_dependents();

      s.dependent(&resBM[1][0], sizeAll);
      JacBM[1].resize(sizeAll * sizeAll);
      s.jacobian(&JacBM[1][0], true);
      s.clear_dependents();

      s.dependent(&resBL[0][0], sizeAll);
      JacBL[0].resize(sizeAll * sizeAll);
      s.jacobian(&JacBL[0][0], true);
      s.clear_dependents();

      s.dependent(&resBL[1][0], sizeAll);
      JacBL[1].resize(sizeAll * sizeAll);
      s.jacobian(&JacBL[1][0], true);
      s.clear_dependents();

      s.clear_independents();

      for(unsigned i = 0; i < sizeAll; i++) {
        for(unsigned j = 0; j < sizeAll; j++) {
          AM(i, j) = JacAM[i * sizeAll + j];
          AL(i, j) = JacAL[i * sizeAll + j];
          BM[0](i, j) = JacBM[0][i * sizeAll + j];
          BM[1](i, j) = JacBM[1][i * sizeAll + j];
          BL[0](i, j) = JacBL[0][i * sizeAll + j];
          BL[1](i, j) = JacBL[1][i * sizeAll + j];
        }
      }

      double perturbation = 1e-10;

      std::cout << "======================EIGEN===================================" << std::endl;

      clock_t start = clock();

      Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
      double inf = 1e+10;

      for(unsigned k = 0; k < 2; k++) {
        double BM0Lk = BM[k].norm();

        for(unsigned i = 0; i < sizeAll; i++) {
          BM[k](i, i) += perturbation * BM0Lk;
        }

        ges.compute(AM, BM[k], false);
        std::complex < double > temp;
        Eigen::VectorXcd eig;

        eig = ges.eigenvalues();
        double emax0 = 0.;
        for(unsigned i = 0; i < sizeAll; i++) {
          temp = eig(i);
          if(fabs(real(temp)) > emax0 && fabs(real(temp)) < inf) {
            emax0 = fabs(real(temp));
          }
        }
        std::cout << iel << " " << emax0 << std::endl;
        sol->_Sol[CMIndex[k]]->set(iel, emax0);
      }

      for(unsigned k = 0; k < 2; k++) {
        double norm = BL[k].norm();

        for(unsigned i = 0; i < sizeAll; i++) {
          BL[k](i, i) += perturbation * norm;
        }

        ges.compute(AL, BL[k], false);
        std::complex < double > temp;
        Eigen::VectorXcd eig;

        eig = ges.eigenvalues();
        double emax0 = 0.;
        for(unsigned i = 0; i < sizeAll; i++) {
          temp = eig(i);
          if(fabs(real(temp)) > emax0 && fabs(real(temp)) < inf) {
            emax0 = fabs(real(temp));
          }
        }
        std::cout << iel << " " << emax0 << std::endl;
        sol->_Sol[CLIndex[k]]->set(iel, emax0);
      }
      eigenTime += (clock() - start);
    } // end of eflag loop
  } //end of element loop

  sol->_Sol[CMIndex[0]]->close();
  sol->_Sol[CMIndex[1]]->close();

  sol->_Sol[CLIndex[0]]->close();
  sol->_Sol[CLIndex[1]]->close();

  //std::cout << std::endl << "petsc TIME:\t" << static_cast<double>(petscTime) / CLOCKS_PER_SEC << std::endl;
  std::cout << std::endl << "Eigen TIME:\t" << static_cast<double>(eigenTime) / CLOCKS_PER_SEC << std::endl;
}



void GetParticleWeights(MultiLevelSolution & mlSol, Line* line3, Line* lineI) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel(level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level);
  unsigned iproc  = msh->processor_id();

  unsigned eflagIndex = mlSol.GetIndex("eflag");

  unsigned soluType = 0; //Linear Lagrange

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  std::vector < std::vector<double> >  x(dim);

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector<Marker*> particle3 = line3->GetParticles();
  std::vector<unsigned> markerOffset3 = line3->GetMarkerOffset();

  unsigned imarker3 = markerOffset3[iproc];

  unsigned m = 3;  // Chebyshev degree

  unsigned i0, i1;
  if(dim == 3) {
    i0 = 0;
    i1 = 3;
  }
  else if(dim == 2) {
    i0 = 3;
    i1 = 5;
  }
  else {
    i0 = 5;
    i1 = 6;
  }

//grab the gauss points with elemtype and degree
  std::string name[6] = {"hex", "tet", "wedge", "quad", "tri", "line"};

  unsigned ng[6];
  Eigen::MatrixXd xg[6];
  Eigen::VectorXd wg[6];
  std::vector < double > jac[6];
  Eigen::MatrixXd Pg[6];

  for(unsigned i = i0; i < i1; i++) {

    const Gauss *gauss = new  Gauss(name[i].c_str(), "fourth");
    ng[i] = gauss->GetGaussPointsNumber();

    std::vector< const double * > Xg(dim);

    for(unsigned k = 0; k < dim; k++) {
      Xg[k] = gauss->GetGaussCoordinatePointer(k);
    }

    xg[i].resize(dim, ng[i]);
    for(unsigned k = 0; k < dim ; k++) {
      for(unsigned j = 0; j < ng[i]; j++) {
        xg[i](k, j) = Xg[k][j];
      }
    }

    const double *Wg = gauss->GetGaussWeightsPointer();
    wg[i].resize(ng[i]);
    for(unsigned j = 0; j < ng[i]; j++) {
      wg[i](j) = Wg[j];
    }

    jac[i].resize(ng[i]);

    Eigen::Tensor<double, 3, Eigen::RowMajor> PmG;
    GetChebXInfo(m, dim, ng[i], xg[i], PmG);

    GetMultiDimChebMatrix(dim, m, ng[i], PmG, Pg[i]);

    delete gauss;
  }

  Eigen::VectorXd F; // holds P_{n}(x_g) * wg * J(xg), n = 0,1,..m, g = 1,2,..ng  multidimensional
  F.resize(pow(m + 1, dim));
  Eigen::MatrixXd A; // // holds P_{n}(x_p) * wp , n = 0,1,..m, p = 1,2,..np multidimensional

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));

    if(eFlag == 1) {
      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);  // number of solution element dofs

      for(int k = 0; k < dim; k++) {
        x[k].resize(nDofu);
      }

      for(unsigned i = 0; i < nDofu; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof
        for(unsigned k = 0; k < dim; k++) {
          x[k][i] = (*msh->_topology->_Sol[k])(xDof);     // global extraction and local storage for the element coordinates
        }
      }


      for(unsigned ig = 0; ig < ng[ielGeom] ; ig++) { // gauss loop to get Jacobians
        std::vector <double> xi(dim);
        for(unsigned k = 0; k < dim; k++) {
          xi[k] = xg[ielGeom](k, ig);
        }
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, jac[ielGeom][ig], phi, phi_x);
      }

      //Assemble F
      F.setZero();
      for(unsigned i = 0; i < pow(m + 1, dim); i++) {
        for(unsigned j = 0; j < ng[ielGeom] ; j++) {
          F(i) += Pg[ielGeom](i, j) * jac[ielGeom][j] * wg[ielGeom](j);
        }
      }

      // identify the first particle inside iel
      while(imarker3 < markerOffset3[iproc + 1] && iel > particle3[imarker3]->GetMarkerElement()) {
        imarker3++;
      }
      unsigned imarker0 = imarker3;
      // loop on all particles inside iel to find how many particles are in iel
      while(imarker3 < markerOffset3[iproc + 1] && iel == particle3[imarker3]->GetMarkerElement()) {
        imarker3++;
      }
      unsigned nmarker = imarker3 - imarker0;
      imarker3 = imarker0;

      unsigned cnt = 0;
      Eigen::MatrixXd xP(dim, nmarker);
      Eigen::VectorXd wP(nmarker);


      // loop on all particles inside iel
      while(imarker3 < markerOffset3[iproc + 1] && iel == particle3[imarker3]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle3[imarker3]->GetMarkerLocalCoordinates();
        //msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle3[imarker3]->GetMarkerMass();

        for(unsigned k = 0; k < dim; k++) {
          xP(k, cnt) = xi[k];
        }
        wP[cnt] = weight;

        cnt++;
        imarker3++;
      }

      Eigen::Tensor<double, 3, Eigen::RowMajor> PmX;
      GetChebXInfo(m, dim, nmarker, xP, PmX);

      GetMultiDimChebMatrix(dim, m, nmarker, PmX, A); //multidimensional Chebyshev polynomial evaluation in particle points up to m

      Eigen::VectorXd w_new;
      SolWeightEigen(A, F, wP, w_new); // New weights for iel are avaliable at this point

//       for(unsigned j = 0; j < nmarker; j++) {
//         std::cout << xP(0, j) << " " << xP(1, j) << " " << wP(j) << " " << w_new(j) << std::endl;
//       }

      // loop on all particles inside iel to attach the optimized weights to iel particles.
      imarker3 = imarker0;
      cnt = 0;
      while(imarker3 < markerOffset3[iproc + 1] && iel == particle3[imarker3]->GetMarkerElement()) {

        particle3[imarker3]->SetMarkerMass(w_new[cnt]);

        imarker3++;
        cnt++;
      }

    } // end of interface loop
  } // end of iel loop

}

