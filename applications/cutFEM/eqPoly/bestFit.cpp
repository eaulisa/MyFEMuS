#include <iostream>
#include <iomanip>

#include "./LiSK/lisk.hpp"
#include <cmath>       /* exp */
#include "bestFit.hpp"
#include <Eigen/Core>
#include <Eigen/SVD>

//using namespace femus;
using namespace std;
using namespace Eigen;


std::vector < double > BestFit::FindBestFit(const std::vector < double > &pts, const std::vector < double > &Npts, const unsigned &dim) {


  std::vector < double > bestfit(dim + 1, 0.);
  unsigned cnt = 0;
  double normaldotcoefficients = 0.;
  double d = 0;
  unsigned numberofpoints = pts.size() / dim;
  MatrixXd m(numberofpoints, dim);
  std::vector < double > N(dim, 0.);
  std::vector < double > centroid(dim, 0.);

//Calculate average Normal and centroid from points
  for(unsigned i = 0; i < numberofpoints; i++) {

    for(unsigned j = 0; j < dim; j++, cnt++) {
      centroid[j] += pts[cnt];
      N[j] += Npts[cnt];
    }

  }
  for(unsigned j = 0; j < dim; j++) {
    N[j] /= numberofpoints;
    centroid[j] /= numberofpoints;
  }

  cnt = 0;

  //Fill matrix to be passed to JacobiSVD
  for(unsigned i = 0; i < numberofpoints; i++) {

    for(unsigned j = 0; j < dim; j++, cnt++) {
      m(i, j) = pts[cnt] - centroid[j];
    }

  }

  JacobiSVD<MatrixXd> svd(m, ComputeThinU | ComputeThinV);
  MatrixXd v = svd.matrixV();


//If dim = 2 and line of best fit is desired, use singualr vector associated with max singular vector

  if(dim <= 2) {


    for(unsigned i = 0; i < dim; i++) {

      bestfit[i] = v(i, 0);
      normaldotcoefficients += bestfit[i] * N[i];
    }

  }

  //If dim = 3 and plane of best fit is desired, use singular vector associated with min singular vector
  if(dim == 3) {

    for(unsigned i = 0; i < dim; i++) {

      bestfit[i] = v(i, dim - 1);
      normaldotcoefficients += bestfit[i] * N[i];
    }

  }

//   for(unsigned i = 0; i < dim; i++) {
//
//     std::cout << " coefficent before dot product "  << bestfit[i] << endl;
//   }

  //Rotate normal by pi if Normal dot coefficents is less than zero
  if(normaldotcoefficients < 0) {

    for(unsigned i = 0; i < dim; i++) {
      bestfit[i] *= -1.;

    }

  }


//
//Calculate constant d in ax+by+d=0 or ax+by+cz+d=0
  for(unsigned i = 0; i < dim; i++) {

    d -= bestfit[i] * centroid[i];
  }

  bestfit[dim] = d;

//   for(unsigned i = 0; i < dim + 1; i++) {
//
//     std::cout << " coefficent "  << bestfit[i] << endl;
//   }




  //std::cout << bestfit[0] * bestfit[0] + bestfit[2] * bestfit[2] + bestfit[1] * bestfit[1]  << "   = norm squared" << endl;
  //std::cout << v << "   v matrix" << endl;
  //std::cout << bestfit[0] * 1. + bestfit[1] * 2.  + bestfit[2] << "  check equation" << std::endl;

  return bestfit;


}









