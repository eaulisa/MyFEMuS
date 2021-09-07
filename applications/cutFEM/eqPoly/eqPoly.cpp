#include "eqPoly.hpp"
#include <boost/math/special_functions/factorials.hpp>
#include <algorithm>
#include <iostream>
#include <stdlib.h>

using namespace boost;

using boost::math::factorial;



//TODO hexahedron = 0, tet = 1, wedge =2, quad = 3, tri = 4, line = 5, point = 6


double EquivalentPolynomial::IntegrationValue(const int &s, const unsigned &dim, const unsigned &degree,  std::vector < double > &normal, const double &d, const unsigned &element) {

  
  double value = 0.;
  unsigned int n = degree + 1;
  unsigned int floorj = 1;
  
  


  if(dim > 1) {

    n = (dim == 3) ? (((degree + 1) * (degree + 2) * (degree + 3)) / 6) : (((degree + 1) * (degree + 2)) / 2);

  }

  _polycoeffs.resize(n);


  std::vector < unsigned > powers(dim);
  std::vector < std::vector< double > > legendreCoefficients{{1.}, {1.}, { 1.5, -0.5}, { 2.5, -1.5}};
  std::vector < std::vector< unsigned > > position{{0}, {1}, {0, 2}, {1, 3}};
  std::vector < unsigned > size = {1, 1, 2, 2};

  // My thought is to calculate the integrals involving the monomials and Li_s, then construct the coefficeints as a linear combination of the evaluated integrals, since the moment matrix will be the identity, after Gram-Schmidt. Then we should be able to evaluate the integral of the the equivalent polynomial and the shape function via Gauss Quadrature.

  if(element == 5 || element == 3 || element == 0) {


    if(dim == 1) {
      std::vector < double > monomials1D(n);
      _coefficients.resize(n);
      powers = {degree};
      

      for(int j = 0; j < n; j++) {
          std::cout <<  "  s = " << s << "  powers = " << powers[0] << "  normal = " << normal[0] << "  d = " << d << std::endl;

        //monomials1D[j] = EquivalentPolynomial::HyperCubeA(s, powers, normal, d);


        for(int l = 0; l < size[n - 1]; l++) {



          _coefficients[j] += legendreCoefficients[degree][l] * monomials1D[position[n - 1][l]];



        }
        
        std::cout << " " << _coefficients[j] << std::endl;


      }

    }

    if(dim == 2) {

      std::vector < std::vector< double > > monomials2D(n);
      std::vector < std::vector< double > > c_star(n);


      for(unsigned i = 0; i < degree; i++) {
        for(unsigned j = 0; j < degree; j++) {
          powers = {i, j};
          //monomials2D[i][j] = EquivalentPolynomial::HyperCubeA(s, powers, normal, d);
          for(int g = 0; i < size[i]; g++) {
            for(int h = 0; j < size[j]; h++) {

              c_star[i][j] += legendreCoefficients[i][g] * legendreCoefficients[j][h] * monomials2D[position[i][g]][position[j][h]];

            }
          }
        }
      }
    }



  }


  return value;
}


//template<typename Type>
double EquivalentPolynomial::LimLi(const int &n, const double & x) {
  if(x < 0.) return 0.;
  else if(n != 0) return -pow(x, n) / factorial<double>(n);
  else if(x > 0.) return -1.;
  else return -0.5;
}

//template<typename Type>
double EquivalentPolynomial::HyperCubeB(const unsigned & n, const int &s, std::vector<unsigned> &m,
           const std::vector <double> &a, const std::vector <double> &ma,
           const double & d, const double & md) {

  double HCI = 0.;
  double aI = 1. / a[n];

  int sl = (m[n] % 2 == 1) ? 1 : -1; // this is (-1)^(m-1)
  int sr = 1; // this is (-1)^(j-1) for j = 1, ..., m + 1
  double c = aI; // this is m!/(m+1-j)! 1/a^j for j = 1,...,m + 1

  for(int j = 1; j <= m[n] + 1;  c *= aI * (m[n] + 1 - j), sr = -sr, j++) {
    HCI += c * (sl * HyperCubeA(n - 1, s + j, m, a, ma, -a[n] + d, a[n] - d) + sr * HyperCubeA(n - 1, s + j, m, a, ma, a[n] + d, -a[n] - d));
  }
  return HCI;
}

//template<typename Type>
double EquivalentPolynomial::HyperCubeA(const unsigned & n, const int &s, std::vector<unsigned> &m,
           const std::vector <double> &a, const std::vector <double> &ma,
           const double & d, const double & md) {
  //typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  //TODO if(n == 0)  return LSI(s, m[0], a[0], d);

  switch(s) {
    case -1: // interface integral
      if(d <= 0) {
        return HyperCubeB(n, s, m, a, ma, d, md);
      }
      else {
        return HyperCubeB(n, s, m, ma, a, md, d);
      }
      break;
    case 0: // step function integral
      if(d <= 0) {
        return HyperCubeB(n, s, m, a, ma, d, md);
      }
      else {
        double fullIntegral = 1.;
        for(unsigned i = 0; i <= n; i++) {
          // TODO fullIntegral *= (m[i] % 2 == 0) ? 2 / Type(m[i] + 1.) : 0;
        }
        return (fullIntegral - HyperCubeB(n, s, m, ma, a, md, d));
      }
      break;
    default: // all other cases, whatever they mean
      if(d < fabs(a[n])) {
        std::cout << "!";
        return HyperCubeB(n, s, m, a, ma, d, md);
      }
//       else {
//         std::cout << "*";
//         return HyperCubeC(n, s, m, a, ma, d, md);
//       }
  }
}

// LSI(const int &s, const unsigned &m, double a, double d) {
//   //typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;
// 
//   switch(s) {
//     case -1:
//       return LSIm1(m, a, d);
//       break;
// 
//     case 0:
//       return LSI0(m, a, d);
//       break;
// 
//     default:
// 
//       if(a == 0) {
//         return -LimLi(s, d) * (1 + pow(-1, m)) / Type(m + 1);
//       }
// 
//       Type INT(0);
// 
//       Type x(a + d);
//       Type y(-a + d);
// 
//       if(d < fabs(a)) { // in all these cases no-significant digit cancellation occurs
//         if(x > 0.) {
//           Type xDa = x / a;
//           Type c1 = s + 1;
//           Type c2 = m + 1;
//           Type fx =  pow(x, c1) / (a * factorial<Type>(s + 1));
//           for(unsigned i = 0; i < m + 1; i++) {
//             INT += fx;
//             fx *= - xDa * (--c2) / (++c1);
//           }
//         }
//         else if(y > 0.) {
//           Type yDa = y / a;
//           Type c1 = s + 1;
//           Type c2 = m + 1;
//           Type fy =  - pow(-1., m) * pow(y, c1) / (a * factorial<Type>(s + 1));
//           for(unsigned i = 0; i < m + 1; i++) {
//             INT += fy;
//             fy *= yDa * (--c2) / (++c1);
//           }
//         }
//       }
//       else { //alternative formula to avoid significant digit cancellations when s>1, and (a+d) and (-a+d) are non-negative and d >> a
// 
//         Type px = 1.;
//         Type py = pow(-1, m + s);
// 
//         Type c1 = m + 1 + s;
//         Type c2 = 1;
//         Type f1 = pow(-a, s) / factorial<Type>(m + 1 + s);
// 
//         for(int i = 0; i < s + 1; i++) {
//           INT += f1 * (px + py);
//           f1 *= -(c1--) / (a * c2++);
//           px *= x;
//           py *= -y;
//         }
//         INT *= factorial<Type>(m);
//       }
//       return INT;
//       break;
//   }
// }


//************************************************************************************************************************


/*void EquivalentPolynomial::MatrixVectorMultiply(const std::vector<std::vector <double>> &A, const std::vector < complex < double > > &bv, std::vector < complex < double > > &xv) {

  xv.resize(bv.size());
  for(unsigned i = 0; i < bv.size(); i++) {
    xv[i] = 0;
    for(unsigned j = 0; j < bv.size(); j++) {
      xv[i] += A[i][j] * bv[j].real();
      //std::cout << xv[i] <<  "  made it here in MatrixVectorMultiply" << endl;
    }
    //std::cout << "  MatrixVectorMultiply " << bv[i] <<  std::endl;
  }

  std::cout <<  "  made it here OUT OG MatrixVectorMultiply" <<  _dim << endl;
}

//This function takes a vector of points as inputs and calculates the best fit plane for those points

// void EquivalentPolynomial::FindBestFit(const std::vector < double > &pts, const std::vector < double > &Npts, const unsigned &dim) {
//
//
//   unsigned cnt = 0;
//   double normaldotcoefficients = 0.;
//   double d = 0;
//   unsigned numberofpoints = pts.size() / dim;
//   MatrixXd m(numberofpoints, dim);
//   _bestfit.resize(dim + 1);
//   std::vector < double > N(dim, 0.);
//   std::vector < double > centroid(dim, 0.);
//
// //Calculate average Normal and centroid from points
//   for(unsigned i = 0; i < numberofpoints; i++) {
//
//     for(unsigned j = 0; j < dim; j++, cnt++) {
//       centroid[j] += pts[cnt];
//       N[j] += Npts[cnt];
//     }
//
//   }
//   for(unsigned j = 0; j < dim; j++) {
//     N[j] /= numberofpoints;
//     centroid[j] /= numberofpoints;
//   }
//
//   cnt = 0;
//
//   //Fill matrix to be passed to JacobiSVD
//   for(unsigned i = 0; i < numberofpoints; i++) {
//
//     for(unsigned j = 0; j < dim; j++, cnt++) {
//       m(i, j) = pts[cnt] - centroid[j];
//     }
//
//   }
//
//   JacobiSVD<MatrixXd> svd(m, ComputeThinU | ComputeThinV);
//   MatrixXd v = svd.matrixV();
//
//
// //If dim = 2 and line of best fit is desired, use singualr vector associated with max singular vector
//
//   if(dim <= 2) {
//
//
//     for(unsigned i = 0; i < dim; i++) {
//
//       _bestfit[i] = v(i, 0);
//       normaldotcoefficients += _bestfit[i] * N[i];
//     }
//
//   }
//
//   //If dim = 3 and plane of best fit is desired, use singualr vector associated with min singular vector
//   if(dim == 3) {
//
//     for(unsigned i = 0; i < dim; i++) {
//
//       _bestfit[i] = v(i, dim - 1);
//       normaldotcoefficients += _bestfit[i] * N[i];
//     }
//
//   }
//
// //   for(unsigned i = 0; i < dim; i++) {
// //
// //     std::cout << " coefficent before dot product "  << _bestfit[i] << endl;
// //   }
//
//   //Rotate normal by pi if Normal dot coefficents is less than zero
//   if(normaldotcoefficients < 0) {
//
//     for(unsigned i = 0; i < dim; i++) {
//       _bestfit[i] *= -1.;
//
//     }
//
//   }
//
// //   double numerator = 0;
// //   double denominator = 0;
// //
// //   for(unsigned i = 0; i < numberofpoints - 1; i++) {
// //
// //     numerator += (pts[2*i] - centroid[0]) * (pts[2*i + 1] - centroid[1]);
// //     denominator += (pts[2*i + 1] - centroid[1]) * (pts[2*i + 1] - centroid[1]) - ((pts[2*i] - centroid[0]) * (pts[2*i] - centroid[0]));
// //
// //   }
// //
// //   if(denominator != 0.) {
// //
// //     _bestfit[0] = cos(0.5 * atan(2*numerator / denominator));
// //     _bestfit[1] = sin(0.5 * atan(2*numerator / denominator));
// //
// //   }
// //
// //Calculate constant d in ax+by+d=0 or ax+by+cz+d=0
//   for(unsigned i = 0; i < dim; i++) {
//
//     d -= _bestfit[i] * centroid[i];
//   }
//
//   _bestfit[dim] = d;
//
// //   for(unsigned i = 0; i < dim + 1; i++) {
// //
// //     std::cout << " coefficent "  << _bestfit[i] << endl;
// //   }
//
//
//
//
//   //std::cout << _bestfit[0] * _bestfit[0] + _bestfit[2] * _bestfit[2] + _bestfit[1] * _bestfit[1]  << "   = norm squared" << endl;
//   //std::cout << v << "   v matrix" << endl;
//   //std::cout << _bestfit[0] * 1. + _bestfit[1] * 2.  + _bestfit[2] << "  check equation" << std::endl;
//
//
//
// }


double EquivalentPolynomial::GetValue(std::vector <double> &x, unsigned & element) {

  unsigned numvalues = x.size();
  unsigned numcoefficients = _coefficients.size();
  //std::vector < double > values(numvalues);
  //std::fill(values.begin(), values.end(), 0.);
  double value = 0.;

  if(_dim == 1) {

    for(unsigned k = 0; k < numvalues; k++) {

      for(unsigned i = 0; i < numcoefficients; i++) {

        value += _coefficients[i].real() * pow(x[k], i);

      }

    }

  }

  if(_dim == 2) {

    if(element == 3) {
      if(_degree == 2) {

        value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] +
                _coefficients[3].real() * x[0] * x[0] + _coefficients[4].real() * x[1] * x[1] + _coefficients[5].real() * x[0] * x[1];

      }

      if(_degree == 3) {

        value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] +
                _coefficients[3].real() * x[0] * x[0] + _coefficients[4].real() * x[1] * x[1] + _coefficients[5].real() * x[0] * x[1]  +
                _coefficients[6].real() * x[0] * x[0] * x[0] + _coefficients[7].real() * x[1] * x[1] * x[1] +
                _coefficients[8].real() * x[0] * x[0] * x[1] + _coefficients[9].real() * x[1] * x[1] * x[0];

      }

      if(_degree == 4) {

        value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] +
                _coefficients[3].real() * x[0] * x[0] + _coefficients[4].real() * x[1] * x[1] + _coefficients[5].real() * x[0] * x[1] +
                _coefficients[6].real() * x[0] * x[0] * x[0] + _coefficients[7].real() * x[1] * x[1] * x[1] + _coefficients[8].real() * x[0] * x[0] * x[1] + _coefficients[9].real() * x[1] * x[1] * x[0]  + _coefficients[10].real() * x[0] * x[0] * x[0] * x[0] +
                _coefficients[11].real() * x[1] * x[1] * x[1] * x[1] + _coefficients[12].real() * x[0] * x[0] * x[1] * x[1] +
                _coefficients[13].real() * x[0] * x[0] * x[0] * x[1] + _coefficients[14].real() * x[1] * x[1] * x[1] * x[0];

      }

      if(_degree == 5) {

        value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] +
                _coefficients[3].real() * x[0] * x[0] + _coefficients[4].real() * x[1] * x[1] + _coefficients[5].real() * x[0] * x[0] * x[0] +
                _coefficients[6].real() * x[1] * x[1] * x[1] + _coefficients[7].real() * x[0] * x[0] * x[0] * x[0] + _coefficients[8].real() * x[1] * x[1] * x[1] * x[1] +
                _coefficients[9].real() * pow(x[0], 5) + _coefficients[10].real() * pow(x[1], 5) +
                _coefficients[11].real() * x[0] * x[0] * x[1] + _coefficients[12].real() * x[0] * x[0] * x[0] * x[1] +
                _coefficients[13].real() * x[0] * x[0] * x[0] * x[0] * x[1] + _coefficients[14].real() * x[1] * x[1] * x[0] +
                _coefficients[15].real() * x[1] * x[1] * x[1] * x[0] + _coefficients[16].real() * x[1] * x[1] * x[1] * x[1] * x[0] +
                _coefficients[17].real() * x[0] * x[0] * x[0] * x[1] * x[1] + _coefficients[18].real() * x[1] * x[1] * x[1] * x[0] * x[0] +
                _coefficients[19].real() * x[0] * x[1] + _coefficients[20].real() * x[0] * x[0] * x[1] * x[1];

      }

    }

    if(element == 4) {

      value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] +
              _coefficients[3].real() * x[0] * x[0] + _coefficients[4].real() * x[1] * x[1] + _coefficients[5].real() * x[0] * x[1];



    }
  }

  if(_dim == 3) {

    if(element == 0) {

      value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] +
              _coefficients[3].real() * x[2] + _coefficients[4].real() * x[0] * x[1] + _coefficients[5].real() * x[0] * x[2] +
              _coefficients[6].real() * x[1] * x[2] + _coefficients[7].real() * x[0] * x[0] +  _coefficients[8].real() * x[1] * x[1] +
              _coefficients[9].real() * x[2] * x[2];
    }
  }

  return value;

}









*/



/*
#include "eqPoly.hpp"
#include <time.h>
//#include "FemusInit.hpp"
//#include "Elem.hpp"
#include <stdlib.h>



//using namespace femus;
using namespace std;
//using namespace Eigen;
using std::cout;

int main(int argc, char** args) {

//   //FemusInit mpinit(argc, args, MPI_COMM_WORLD);
// 
//   //SetCoefficients takes ( dim, degree of equivalent polynomial, rho, vector < a, b, c > for the dicontinuity ( point, line, or plane ), and element (3=triangle/tet, 4=square/cube
//   // a*x + b*y + c*z + d = 0, and d) as inputs
//   EquivalentPolynomial eqP;
// 
//   //
//    eqP.SetCoefficients(3, 2, 2, std::vector<double> {1., 2., 1.}, 0., 3);
//    eqP.PrintCoefficients();
// //   std::cout << eqP.GetValue(std::vector<double> {0.5, 0.5}) << " " << std::endl;
// 
// 
//   std::vector < double >points {1.,2.,4.,5.,7.,8.,-5.,23.,12.,15.,14.,14.};
//   //std::vector < double >points(2000);
//   std::vector < double >normal {1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0.};
//   unsigned dim = 3;
//   eqP.FindBestFit(points, normal, dim);
//   std::vector < double > onepoint {1.,2.,3.};
//   unsigned element = 0;
//   std::cout << eqP.GetValue(onepoint, element) << "  value" << endl;
//   
  
  
        
  EquivalentPolynomial eqP;
  
  std::vector < double > N {1.};
  
  
  eqP.IntegrationValue(0, 1, 5, N, 0., 5 );
  
  
  
  
  return 0;
}*/


