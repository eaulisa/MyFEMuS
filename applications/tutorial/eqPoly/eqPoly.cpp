//
//  main.cpp
//  example
//
//  Created by Sebastian on 18/05/16.
//  Copyright © 2016 Sebastian Kirchner. All rights reserved.
//
//  This file is part of LiSK.
//
//  LiSK is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  LiSK is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with LiSK.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <iomanip>
#include "./LiSK/lisk.hpp"



#include <iostream>
#include <stdio.h>      /* printf */
#include <math.h>       /* exp */
using namespace std;




// //Fake function just to be a place holder for the Lisk Li function
// int Li(int i, double j)
// {
//   return i * j;
// }

// This function accepts: the diminsion, dim; polynomial degree, degree; generalized Heavyside funcion paprameter, p;
// discontinuity (point, line, or plane) ax + by + cz = d, where b and c are set to zero by default;
// The function returns an array of the equivalent polynomial coefficients, coefficients
complex<double> * getCoefficients (const int &dim, const int &degree, const double &p, const double & a, const double & d, const double b = 0, const double c = 0)
{
  LiSK::LiSK< complex<double> > lisk(5);

  if(dim == 1) {

    //TODO

    if(degree == 3) {

      double x1r = -exp(p - d * p) ;
      const complex<double> x1(x1r, 0.);

      double x2r = -exp((-1 - d) * p);
      const complex<double> x2(x2r, 0.);
      
      double x3r = -exp((-1 - d) * p);
      const complex<double> x3(x2r, 0.);


      static complex < double > coefficients[4];

      coefficients[0] = (1./(4. * pow(p,3))) * (14. * pow(p,3) + 9. * pow(p,2) * log(1. + exp((-1. + d) * p)) + 

                        15. * pow(p,2) * log(1. + exp(-(1. + d) * p)) - 9. * pow(p,2) * log(1. + exp((1. + d) * p)) - 

                        15. * pow(p,2) * log(1. + exp((p - d* p))) - 30. * p * lisk.Li(2, x2) - 

                        30. * p * lisk.Li(2, x1) - (30. * lisk.Li(3, x2)) + 30. * lisk.Li(3, x1));

                        

      coefficients[1] = (-15.*(pow(p,3)*log(1. + exp((-1. - d)*p)) + 

                    pow(p,3)*log(1. + exp(p - d*p)) - 

                    8.*pow(p,2)*lisk.Li(2,x2) + 

                    8.*pow(p,2)*lisk.Li(2,x1) - 

                    21.*p*lisk.Li(3,x2) - 

                    21.*p*lisk.Li(3,x1) - 

                    21.*lisk.Li(4,x2) + 

                    21.*lisk.Li(4,x1)))/(2.*pow(p,4));


      coefficients[2] = -7.5 - (45. * log(1. + exp((-1. - d) * p)))/(4. * p) - 

                        (15. * log(1. + exp((-1. + d) * p)))/(4. * p) + 

                        (15. * log(1. + exp((1. + d) * p)))/(4. * p) + 

                        (45. * log(1. + exp(p - d * p)))/(4. * p) + 

                        (45. * lisk.Li(2,-exp((-1. - d) * p)))/(2. * pow(p,2)) + 

                        (45. * lisk.Li(2,-exp(p - d * p)))/(2. * pow(p,2)) + 

                        (45. * lisk.Li(3,-exp((-1 - d) * p)))/(2. * pow(p,3)) - 

                        (45. * lisk.Li(3,-exp(p - d * p)))/(2. * pow(p,3));


      coefficients[3] = (35. * log(1. + exp((-1. - d) * p)))/(2. * p) + 

                        (35. * log(1. + exp(p - d * p)))/(2. * p) - 

                        (105. * lisk.Li(2,-exp((-1. - d) * p)))/pow(p,2) + 

                        (105. * lisk.Li(2,-exp(p - d * p)))/pow(p,2) - 

                        (525. * lisk.Li(3,-exp((-1. - d) * p)))/(2. * pow(p,3)) - 

                        (525. * lisk.Li(3,-exp(p - d * p)))/(2. * pow(p,3)) - 

                        (525. * lisk.Li(4,-exp((-1. - d) * p)))/(2. * pow(p,4)) + 

                        (525. * lisk.Li(4,-exp(p - d * p)))/(2. * pow(p,4));

      return coefficients;
    }

  }

  else if(dim == 2) {
    double coefficients[((degree + 1) * (degree + 2)) / 2];

    //TODO

  }

  else if (dim == 3) {
    double coefficients[((degree + 1) * (degree + 2) * (degree + 3)) / 6];

    //TODO


  }


}

int main ()
{
  complex < double > * ptr = getCoefficients(1, 3, 50, 1, 0.7);
  std::cout<< ptr[0] << " "<< ptr[1] << " "<<  ptr[2]<< " "<<  ptr[3] << std::endl;
  std::cout<< ptr[0].real() << " "<< ptr[1].real() << " "<<  ptr[2].real()<< " "<<  ptr[3].real() << std::endl;
  return 0;
}

/*
using namespace cln;*/

int main1(int argc, const char * argv[]) {

  using namespace std;

  complex<double> x(1 / 4., 1 / 4.);
  double xr = 1.;
  double xi = 1.;
  x = complex <double> (xr, xi);

  x = complex <double> (real(x), imag(x));

  x.real(2.);

  std::cout << x << std::endl;


  // Double precision example for Li_n and Li_22
  try {
    //Create a LiSK object. Constants are prepared for
    //computations of Li_n up to n=20
    LiSK::LiSK<complex<double>> lisk(20);

    complex<double> x(1 / 4., 1 / 4.);
    double xr = 1.;
    double xi = 1.;
    x = complex <double> (xr, xi);
    const vector<int> weights = {1, 2, 3, 4, 5, 6, 10, 15, 20};

    // Compute Li_n(x) for n=1,2,3,4,5,6,10,15,20 at x=1/4+I*1/4
    for (auto n : weights) {
      cout << "Li(" << n << ",x) = " << setprecision(17) << lisk.Li(n, x) << endl;
    }

    // Compute Li_22(x,y) at x=y=1/4+I*1/4
    cout << "\nLi_22(x,y) = " << setprecision(17) << lisk.Li22(x, x) << endl << endl;


  }
  catch (runtime_error& e) {
    cout << e.what() << endl;
  }


  // Arbitrary precision example for Li_n and Li_22
  try {
    //Create a LiSK object. Constants are prepared for
    //computations of Li_n up to n=20 with 34 digit precision
    LiSK::LiSK<cln::cl_N> lisk(20, 34);

    const cln::cl_N x = cln::complex("1/4", "1/4");
    const vector<int> weights = {1, 2, 3, 4, 5, 6, 10, 15, 20};

    // Compute Li_n(x) for n=1,2,3,4,5,6,10,15,20 at x=1/4+I*1/4
    for (auto n : weights) {
      cout << "Li(" << n << ",x) = " << setprecision(17) << lisk.Li(n, x) << endl;
    }

    // Compute Li_22(x,y) at x=y=1/4+I*1/4
    cout << "\nLi_22(x,y) = " << setprecision(17) << lisk.Li22(x, x) << endl << endl;


  }
  catch (runtime_error& e) {
    cout << e.what() << endl;
  }

  return 0;
}


