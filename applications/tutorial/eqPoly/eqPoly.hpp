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
#include <cmath>       /* exp */

using namespace std;

class EquivalentPolynomial {
  public:
    EquivalentPolynomial() {
      _lisk = new LiSK::LiSK< complex<double> > (5);
    }
    ~EquivalentPolynomial() {
      delete _lisk;
    }

    void SetCoefficients(const unsigned &dim, const unsigned &degree, const double &p, const std::vector < double > &c, const double & d);
    
    const std::vector < complex < double > > &GetCoefficients() {
      return _coefficients;
    };

    void PrintCoefficients() {
      for(unsigned i = 0; i < _coefficients.size(); i++) {
        std::cout << _coefficients[i].real() << " ";
      }
      std::cout << std::endl;
    }

    double GetValue(std::vector <double> x) {
      double value;
      if(_dim == 1) {
        //value =
      }


      return value;
    }

  private:
    std::vector < complex < double > > _coefficients;
    LiSK::LiSK< complex<double> > *_lisk;
    unsigned _dim, _degree;
};

void EquivalentPolynomial::SetCoefficients(const unsigned &dim, const unsigned &degree, const double &p, const std::vector < double > &c, const double & d) {

  _dim = dim;
  _degree = degree;
  if(dim == 1) {
    if(degree <= 3) {
      _coefficients.resize(4);

      double x1 = -exp(p - d * p) ;
      double x2 = -exp((-1 - d) * p);
      double x3 = -exp((-1 - d) * p);

      _coefficients[0] = (1. / (4. * pow(p, 3))) * (14. * pow(p, 3) + 9. * pow(p, 2) * log(1. + exp((-1. + d) * p)) +
                                                    15. * pow(p, 2) * log(1. + exp(-(1. + d) * p)) - 9. * pow(p, 2) * log(1. + exp((1. + d) * p)) -
                                                    15. * pow(p, 2) * log(1. + exp((p - d * p))) - 30. * p * _lisk->Li(2, x2) -
                                                    30. * p * _lisk->Li(2, x1) - (30. * _lisk->Li(3, x2)) + 30. * _lisk->Li(3, x1));

      _coefficients[1] = (-15.*(pow(p, 3) * log(1. + exp((-1. - d) * p)) +
                                pow(p, 3) * log(1. + exp(p - d * p)) -
                                8.*pow(p, 2) * _lisk->Li(2, x2) +
                                8.*pow(p, 2) * _lisk->Li(2, x1) -
                                21.* p * _lisk->Li(3, x2) -
                                21.* p * _lisk->Li(3, x1) -
                                21.* _lisk->Li(4, x2) +
                                21.* _lisk->Li(4, x1))) / (2.*pow(p, 4));

      _coefficients[2] = -7.5 - (45. * log(1. + exp((-1. - d) * p))) / (4. * p) -
                         (15. * log(1. + exp((-1. + d) * p))) / (4. * p) +
                         (15. * log(1. + exp((1. + d) * p))) / (4. * p) +
                         (45. * log(1. + exp(p - d * p))) / (4. * p) +
                         (45. * _lisk->Li(2, -exp((-1. - d) * p))) / (2. * pow(p, 2)) +
                         (45. * _lisk->Li(2, -exp(p - d * p))) / (2. * pow(p, 2)) +
                         (45. * _lisk->Li(3, -exp((-1 - d) * p))) / (2. * pow(p, 3)) -
                         (45. * _lisk->Li(3, -exp(p - d * p))) / (2. * pow(p, 3));

      _coefficients[3] = (35. * log(1. + exp((-1. - d) * p))) / (2. * p) +
                         (35. * log(1. + exp(p - d * p))) / (2. * p) -
                         (105. * _lisk->Li(2, -exp((-1. - d) * p))) / pow(p, 2) +
                         (105. * _lisk->Li(2, -exp(p - d * p))) / pow(p, 2) -
                         (525. * _lisk->Li(3, -exp((-1. - d) * p))) / (2. * pow(p, 3)) -
                         (525. * _lisk->Li(3, -exp(p - d * p))) / (2. * pow(p, 3)) -
                         (525. * _lisk->Li(4, -exp((-1. - d) * p))) / (2. * pow(p, 4)) +
                         (525. * _lisk->Li(4, -exp(p - d * p))) / (2. * pow(p, 4));
    }
    else {
      std::cout << "Degree " << degree << " has not been implemented for dimension " << dim << std::endl;
      abort();
    }
  }
  else if(dim == 2) {
    if(degree <= 2) {
      unsigned maxDegree = 2u;
      unsigned n = ((maxDegree + 1u) * (maxDegree + 2u)) / 2u;
      _coefficients.resize(n);
      _lisk = new LiSK::LiSK< complex<double> > (5);

      double a = c[0];
      double b = c[1];

      _coefficients[0] = (-10. *  pow(a, 3) *  pow(b, 3) *  pow(p, 4) -  pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, - exp((-b - d) * p)) +
                          pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, - exp((a - b - d) * p)) +  pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, - exp((b - d) * p)) -
                          pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, - exp((a + b - d) * p)) +  pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, - exp((-b + d) * p)) +
                          15. *  pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, - exp((-a - b + d) * p)) -
                          16. *  pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, - exp((a - b + d) * p)) -
                          pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, - exp((b + d) * p)) -
                          15. *  pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, - exp((-a + b + d) * p)) +
                          16. *  pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, - exp((a + b + d) * p)) + 30. *  pow(a, 2) * b * p * _lisk->Li(3, - exp((-a - b + d) * p)) +
                          30. * a *  pow(b, 2) * p * _lisk->Li(3, - exp((-a - b + d) * p)) - 30. *  pow(a, 2) * b * p * _lisk->Li(3, - exp((a - b + d) * p)) +
                          30. * a *  pow(b, 2) * p * _lisk->Li(3, - exp((a - b + d) * p)) + 30. *  pow(a, 2) * b * p * _lisk->Li(3, - exp((-a + b + d) * p)) -
                          30. * a *  pow(b, 2) * p * _lisk->Li(3, - exp((-a + b + d) * p)) - 30. *  pow(a, 2) * b * p * _lisk->Li(3, - exp((a + b + d) * p)) -
                          30. * a *  pow(b, 2) * p * _lisk->Li(3, - exp((a + b + d) * p)) + 30. *  pow(a, 2) * _lisk->Li(4, - exp((-a - b + d) * p)) +
                          30. *  pow(b, 2) * _lisk->Li(4, - exp((-a - b + d) * p)) - 30. *  pow(a, 2) * _lisk->Li(4, - exp((a - b + d) * p)) -
                          30. *  pow(b, 2) * _lisk->Li(4, - exp((a - b + d) * p)) - 30. *  pow(a, 2) * _lisk->Li(4, - exp((-a + b + d) * p)) -
                          30. *  pow(b, 2) * _lisk->Li(4, - exp((-a + b + d) * p)) + 30. *  pow(a, 2) * _lisk->Li(4, - exp((a + b + d) * p)) +
                          30. *  pow(b, 2) * _lisk->Li(4, - exp((a + b + d) * p))) / (8. *  pow(a, 3) *  pow(b, 3) *  pow(p, 4));

      _coefficients[1] = (3. * (a * p * _lisk->Li(2, -exp((-a - b + d) * p)) + a * p * _lisk->Li(2, -exp((a - b + d) * p)) - a * p * _lisk->Li(2, -exp((-a + b + d) * p)) -
                                a * p * _lisk->Li(2, -exp((a + b + d) * p)) + _lisk->Li(3, -exp((-a - b + d) * p)) - _lisk->Li(3, -exp((a - b + d) * p)) -
                                _lisk->Li(3, -exp((-a + b + d) * p)) + _lisk->Li(3, -exp((a + b + d) * p)))) / (2. * pow(a, 2) * b * pow(p, 3));


      _coefficients[2] = (-3. * (a * b * pow(p, 2) * log(1. + exp((-a - b + d) * p)) + a * b * pow(p, 2) * log(1. + exp((a - b + d) * p)) +
                                 a * b * pow(p, 2) * log(1. + exp((-a + b + d) * p)) - a * b * pow(p, 2) *
                                 log(exp((-a - b + d) * p) * (1. + exp((a + b - d) * p)) * (1. + exp((-a + b + d) * p))) +
                                 a * b * pow(p, 2) * log(1. + exp((a + b + d) * p)) - a * b * pow(p, 2) * log((1. + exp((a - b + d) * p)) * (1. + exp((a + b + d) * p))) -
                                 b * p * _lisk->Li(2, -exp((-a - b + d) * p)) + b * p * _lisk->Li(2, -exp((a - b + d) * p)) - b * p * _lisk->Li(2, -exp((-a + b + d) * p)) +
                                 b * p * _lisk->Li(2, -exp((a + b + d) * p)) - _lisk->Li(3, -exp((-a - b + d) * p)) + _lisk->Li(3, -exp((a - b + d) * p)) +
                                 _lisk->Li(3, -exp((-a + b + d) * p)) - _lisk->Li(3, -exp((a + b + d) * p)))) / (2. * a * pow(b, 2) * pow(p, 3));


      _coefficients[3] = (-15. * (2. * pow(a, 3) * b * pow(p, 4) + pow(a, 2) * pow(p, 2) * _lisk->Li(2, -exp((-b - d) * p)) -
                                  pow(a, 2) * pow(p, 2) * _lisk->Li(2, -exp((a - b - d) * p)) - pow(a, 2) * pow(p, 2) * _lisk->Li(2, -exp((b - d) * p)) +
                                  pow(a, 2) * pow(p, 2) * _lisk->Li(2, -exp((a + b - d) * p)) - pow(a, 2) * pow(p, 2) * _lisk->Li(2, -exp((-b + d) * p)) +
                                  3. * pow(a, 2) * pow(p, 2) * _lisk->Li(2, -exp((-a - b + d) * p)) - 2. * pow(a, 2) * pow(p, 2) * _lisk->Li(2, -exp((a - b + d) * p)) +
                                  pow(a, 2) * pow(p, 2) * _lisk->Li(2, -exp((b + d) * p)) - 3. * pow(a, 2) * pow(p, 2) * _lisk->Li(2, -exp((-a + b + d) * p)) +
                                  2. * pow(a, 2) * pow(p, 2) * _lisk->Li(2, -exp((a + b + d) * p)) + 6. * a * p * _lisk->Li(3, -exp((-a - b + d) * p)) +
                                  6. * a * p * _lisk->Li(3, -exp((a - b + d) * p)) - 6. * a * p * _lisk->Li(3, -exp((-a + b + d) * p)) - 6. * a * p * _lisk->Li(3, -exp((a + b + d) * p)) +
                                  6. * _lisk->Li(4, -exp((-a - b + d) * p)) - 6. * _lisk->Li(4, -exp((a - b + d) * p)) - 6. * _lisk->Li(4, -exp((-a + b + d) * p)) +
                                  6. * _lisk->Li(4, -exp((a + b + d) * p)))) / (8. * pow(a, 3) * b * pow(p, 4));


      _coefficients[4] = (15. * (2. * a * pow(b, 3) * pow(p, 4) + pow(b, 2) * pow(p, 2) * _lisk->Li(2, -exp((-b - d) * p)) -
                                 pow(b, 2) * pow(p, 2) * _lisk->Li(2, -exp((a - b - d) * p)) - pow(b, 2) * pow(p, 2) * _lisk->Li(2, -exp((b - d) * p)) +
                                 pow(b, 2) * pow(p, 2) * _lisk->Li(2, -exp((a + b - d) * p)) - pow(b, 2) * pow(p, 2) * _lisk->Li(2, -exp((-b + d) * p)) +
                                 pow(b, 2) * pow(p, 2) * _lisk->Li(2, -exp((a - b + d) * p)) + pow(b, 2) * pow(p, 2) * _lisk->Li(2, -exp((b + d) * p)) -
                                 pow(b, 2) * pow(p, 2) * _lisk->Li(2, -exp((a + b + d) * p)) - 3. * b * p * _lisk->Li(3, -exp((-a - b + d) * p)) +
                                 3. * b * p * _lisk->Li(3, -exp((a - b + d) * p)) - 3. * b * p * _lisk->Li(3, -exp((-a + b + d) * p)) + 3. * b * p * _lisk->Li(3, -exp((a + b + d) * p)) -
                                 3. * _lisk->Li(4, -exp((-a - b + d) * p)) + 3. * _lisk->Li(4, -exp((a - b + d) * p)) + 3. * _lisk->Li(4, -exp((-a + b + d) * p)) -
                                 3. * _lisk->Li(4, -exp((a + b + d) * p)))) / (4. * a * pow(b, 3) * pow(p, 4));


      _coefficients[5] = (-3. * (2. * pow(a, 3) * b * pow(p, 4) - 3. * pow(a, 2) * b * pow(p, 3) * log(1. + exp((a - b - d) * p)) -
                                 3. * pow(a, 2) * b * pow(p, 3) * log(1. + exp((a + b - d) * p)) + 3. * pow(a, 2) * b * pow(p, 3) * log(1. + exp((a - b + d) * p)) +
                                 3. * pow(a, 2) * b * pow(p, 3) * log(exp((-a - b + d) * p) * (1. + exp((a + b - d) * p)) * (1. + exp((-a + b + d) * p))) +
                                 3. * pow(a, 2) * b * pow(p, 3) * log(1. + exp((a + b + d) * p)) -
                                 3. * pow(a, 2) * b * pow(p, 3) * log((1. + exp((a - b + d) * p)) * (1. + exp((a + b + d) * p))) -
                                 6. * a * b * pow(p, 2) * _lisk->Li(2, -exp((a - b - d) * p)) - 6. * a * b * pow(p, 2) * _lisk->Li(2, -exp((a + b - d) * p)) +
                                 6. * a * b * pow(p, 2) * _lisk->Li(2, -exp((a - b + d) * p)) + 6. * a * b * pow(p, 2) * _lisk->Li(2, -exp((a + b + d) * p)) -
                                 6. * b * p * _lisk->Li(3, -exp((-b - d) * p)) + 6. * b * p * _lisk->Li(3, -exp((a - b - d) * p)) - 6. * b * p * _lisk->Li(3, -exp((b - d) * p)) +
                                 6. * b * p * _lisk->Li(3, -exp((a + b - d) * p)) + 6. * b * p * _lisk->Li(3, -exp((-b + d) * p)) + 6. * a * p * _lisk->Li(3, -exp((-a - b + d) * p)) +
                                 6. * a * p * _lisk->Li(3, -exp((a - b + d) * p)) - 6. * b * p * _lisk->Li(3, -exp((a - b + d) * p)) + 6. * b * p * _lisk->Li(3, -exp((b + d) * p)) -
                                 6. * a * p * _lisk->Li(3, -exp((-a + b + d) * p)) - 6. * a * p * _lisk->Li(3, -exp((a + b + d) * p)) - 6. * b * p * _lisk->Li(3, -exp((a + b + d) * p)) +
                                 6. * _lisk->Li(4, -exp((-a - b + d) * p)) - 6. * _lisk->Li(4, -exp((a - b + d) * p)) - 6. * _lisk->Li(4, -exp((-a + b + d) * p)) +
                                 6. * _lisk->Li(4, -exp((a + b + d) * p)))) / (4. * pow(a, 2) * pow(b, 2) * pow(p, 4));

    }
    else {
      std::cout << "Degree " << degree << " has not been implemented for dimension " << dim << std::endl;
      abort();
    }
  }
  else if(dim == 3) {
    {
      std::cout << "Degree " << degree << " has not been implemented for dimension " << dim << std::endl;
      abort();
    }
  }
  else {
    std::cout << "Wrong Dimension " << dim << std::endl;
    abort();
  }
}

