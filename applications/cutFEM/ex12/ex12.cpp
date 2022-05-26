#include </usr/include/eigen3/Eigen/Core>
#include </usr/include/eigen3/Eigen/SVD>
#include <vector>
#include <iostream>
#include "MyEigenFunctions.hpp"

//using namespace femus;
using namespace std;
using namespace Eigen;


void FindBestFit(const std::vector < double > &xp, const std::vector < double > &w, const std::vector < double > &N, std::vector < double > &a, double &d);

int main() {


  unsigned dim = 2;
  unsigned np = 3;
  std::vector <double> a;
  double d;

  std::vector < std::vector<double> > xp = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},{ 0.5, 0.5, 0.1}};
  std::vector<double> w = {1, 1, 1, 1, 8};
  femus::FindBestFit(xp, w, {0., 0., 1}, a, d);

  std::cout << a[0] << " " << a[1] <<" " << a[2]<<" "<< d<<std::endl;
  
//   xp = {{0, 0}, {1, 1}, {2, 0}};
//   w = {2, 3, 4};
//   femus::FindBestFit(xp, w, {-0.1, -1}, a, d);
// 
//   std::cout << a[0] << " " << a[1] <<" "<< d<<std::endl;
// 
//   
//   xp = {{0, 0}, {1, 1}, {1, 1}, {2, 0},{ 2, 0}};
//   w = {1, 1, 1, 1, 1};
//   femus::FindBestFit(xp, w, {0.1, 1}, a, d);
// 
//   std::cout << a[0] << " " << a[1] <<" "<< d<<std::endl;
//   
//   femus::FindBestFit(xp, boost::none, {0.1, 1}, a, d);
// 
//   std::cout << a[0] << " " << a[1] <<" "<< d<<std::endl;


  
/*
  std::vector<double> yp = {0, 0, 1, 1, 1, 1, 2, 0, 2, 0};
  w = {2, 1, 2, 3, 1};
  FindBestFit(yp, w, {0.1, 1}, a, d);*/

  
  return 1;
}


void FindBestFit(const std::vector < double > &xp, const std::vector < double > &w, const std::vector < double > &N, std::vector < double > &a, double &d) {

  const unsigned& dim = N.size();
  a.resize(dim);

  unsigned np = xp.size() / dim;
  MatrixXd m(np, dim);

  std::vector < double > xg(dim, 0.);

  //Calculate centroid
  unsigned cnt = 0;
  double wSum = 0;
  for(unsigned i = 0; i < np; i++) {
    wSum += w[i];
    for(unsigned j = 0; j < dim; j++, cnt++) {
      xg[j] += w[i] * xp[cnt];
    }
  }
  for(unsigned j = 0; j < dim; j++) {
    xg[j] /= wSum;
  }

  //Fill matrix to be passed to JacobiSVD
  cnt = 0;
  for(unsigned i = 0; i < np; i++) {
    for(unsigned j = 0; j < dim; j++, cnt++) {
      m(i, j) = sqrt(w[i]) * (xp[cnt] - xg[j]);
    }
  }

  JacobiSVD<MatrixXd> svd(m, ComputeThinU | ComputeThinV);
  MatrixXd v = svd.matrixV();

  // use singular vector associated with min singular vector
  double aDotN = 0.;
  for(unsigned i = 0; i < dim; i++) {
    a[i] = v(i, dim - 1);
    aDotN += a[i] * N[i];
  }

  //Rotate normal by pi if Normal dot coefficents is less than zero
  if(aDotN < 0) {
    for(unsigned i = 0; i < dim; i++) {
      a[i] *= -1.;
    }
  }


  //Calculate constant d in ax+by+d=0 or ax+by+cz+d=0
  d = 0.;
  for(unsigned i = 0; i < dim; i++) {
    d -= a[i] * xg[i];
    std::cout << a[i] << " ";
  }
  std::cout << d << std::endl;


}









