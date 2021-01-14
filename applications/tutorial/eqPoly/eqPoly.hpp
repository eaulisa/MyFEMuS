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
      
        unsigned numvalues = x.size();
        unsigned numcoefficients = _coefficients.size();
        //std::vector < double > values(numvalues);
        //std::fill(values.begin(), values.end(), 0.);
        double value = 0.;
        
        if(_dim == 1) {
            
            for(unsigned k = 0; k < numvalues; k++){
                
                for(unsigned i = 0; i < numcoefficients; i++) {
            
                    value += _coefficients[i].real() * pow(x[k], i);
                
                }
            
            }   
          
        }
        
        if(_dim == 2){
            
            if(_degree == 2){
            
                value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] + 
                        _coefficients[3].real() * pow(x[0], 2) + _coefficients[4].real() * pow(x[1], 2) + _coefficients[5].real() * x[0] * x[1];
                
            }
            
            if(_degree == 3){
            
                value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] + 
                        _coefficients[3].real() * pow(x[0], 2) + _coefficients[4].real() * pow(x[1], 2) + _coefficients[5].real() * pow(x[0], 3) + _coefficients[6].real() * pow(x[1], 3) + _coefficients[7].real() * pow(x[0], 2) * x[1] + _coefficients[8].real() * pow(x[1], 2) * x[0] + _coefficients[9].real() * x[0] * x[1];
                
            }
            
            if(_degree == 4){
            
                value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] + 
                        _coefficients[3].real() * pow(x[0], 2) + _coefficients[4].real() * pow(x[1], 2) + _coefficients[5].real() * pow(x[0], 3) + _coefficients[6].real() * pow(x[1], 3) + _coefficients[7].real() * pow(x[0], 4) + _coefficients[8].real() * pow(x[1], 4) + _coefficients[9].real() * pow(x[0], 2) * x[1] + _coefficients[10].real() * pow(x[0], 3) * x[1] + _coefficients[11].real() * pow(x[1], 2) * x[0] + _coefficients[12].real() * pow(x[1], 3) * x[0] +  _coefficients[13].real() * x[0] * x[1] + _coefficients[14].real() * pow(x[0], 2) * pow(x[1], 2);
                
            }
            
            if(_degree == 5){
            
                value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] + 
                        _coefficients[3].real() * pow(x[0], 2) + _coefficients[4].real() * pow(x[1], 2) + _coefficients[5].real() * pow(x[0], 3) + _coefficients[6].real() * pow(x[1], 3) + _coefficients[7].real() * pow(x[0], 4) + _coefficients[8].real() * pow(x[1], 4) + _coefficients[9].real() * pow(x[0], 5) + _coefficients[10].real() * pow(x[1], 5) + _coefficients[11].real() * pow(x[0], 2) * x[1] + _coefficients[12].real() * pow(x[0], 3) * x[1] + _coefficients[13].real() * pow(x[0], 4) * x[1] + _coefficients[14].real() * pow(x[1], 2) * x[0] + _coefficients[15].real() * pow(x[1], 3) * x[0] + _coefficients[16].real() * pow(x[1], 4) * x[0] + _coefficients[17].real() * pow(x[0], 3) * pow(x[1], 2) + _coefficients[18].real() * pow(x[1], 3) * pow(x[0], 2) + _coefficients[19].real() * x[0] * x[1] + _coefficients[20].real() * pow(x[0], 2) * pow(x[1], 2);
                
            }
            
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
    if(degree == 3) {
      _coefficients.resize(4);

      double x1 = -exp(p - d * p) ;
      double x2 = -exp((-1 - d) * p);
      double x3 = -exp((-1 + d) * p);
      double x4 = -exp((1 + d) * p);

      _coefficients[0] = (1. / (4. * pow(p, 3))) * (14. * pow(p, 3) + 
                            9. * pow(p, 2) * log(1. - x3) +
                            15. * pow(p, 2) * log(1. - x2) - 
                            9. * pow(p, 2) * log(1. - x4) -
                            15. * pow(p, 2) * log(1. - x1) - 
                            30. * p * _lisk->Li(2, x2) -
                            30. * p * _lisk->Li(2, x1) - (30. * _lisk->Li(3, x2)) + 
                            30. * _lisk->Li(3, x1));

      _coefficients[1] = (-15.*(pow(p, 3) * log(1. - x2) +
                            pow(p, 3) * log(1. - x1) -
                            8.*pow(p, 2) * _lisk->Li(2, x2) +
                            8.*pow(p, 2) * _lisk->Li(2, x1) -
                            21.* p * _lisk->Li(3, x2) -
                            21.* p * _lisk->Li(3, x1) -
                            21.* _lisk->Li(4, x2) +
                            21.* _lisk->Li(4, x1))) / (2.*pow(p, 4));

      _coefficients[2] = -7.5 - (45. * log(1. - x2)) / (4. * p) -
                            (15. * log(1. - x3)) / (4. * p) +
                            (15. * log(1. - x4)) / (4. * p) +
                            (45. * log(1. - x1)) / (4. * p) +
                            (45. * _lisk->Li(2, x2)) / (2. * pow(p, 2)) +
                            (45. * _lisk->Li(2, x1)) / (2. * pow(p, 2)) +
                            (45. * _lisk->Li(3, x2)) / (2. * pow(p, 3)) -
                            (45. * _lisk->Li(3, x1)) / (2. * pow(p, 3));

      _coefficients[3] = (35. * log(1. - x2)) / (2. * p) +
                            (35. * log(1. - x1)) / (2. * p) -
                            (105. * _lisk->Li(2, x2)) / pow(p, 2) +
                            (105. * _lisk->Li(2, x1)) / pow(p, 2) -
                            (525. * _lisk->Li(3, x2)) / (2. * pow(p, 3)) -
                            (525. * _lisk->Li(3, x1)) / (2. * pow(p, 3)) -
                            (525. * _lisk->Li(4, x2)) / (2. * pow(p, 4)) +
                            (525. * _lisk->Li(4, x1)) / (2. * pow(p, 4));
    }
    else {
      std::cout << "Degree " << degree << " has not been implemented for dimension " << dim << std::endl;
      abort();
    }
  }
  else if(dim == 2) {
    unsigned n = ((degree + 1u) * (degree + 2u)) / 2u;
    _coefficients.resize(n);
    _lisk = new LiSK::LiSK< complex<double> > (7);
    double a = c[0];
    double b = c[1];
    double x1 = - exp((-b - d) * p);
    double x2 = - exp((a - b - d) * p);
    double x3 = - exp((b - d) * p);
    double x4 = - exp((a + b - d) * p);
    double x5 = - exp((-b + d) * p);
    double x6 = - exp((-a - b + d) * p);
    double x7 = - exp((a - b + d) * p);
    double x8 = - exp((b + d) * p);
    double x9 = - exp((-a + b + d) * p);
    double x10 = - exp((a + b + d) * p);
    double x11 = - exp((b + d) * p);
    double x12 = - exp((-a - b - d) * p);
    double x13 = - exp((a + b) * p);
    double x14 = - exp((-a + b - d));
    
   
      
    if(degree <= 2) {
      

      
      
      _coefficients[0] = (-10. *  pow(a, 3) *  pow(b, 3) *  pow(p, 4) -  
                            pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, x1) +
                            pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, x2) + 
                            pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, x3) -
                            pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, x4) +  
                            pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, x5) +
                            15. *  pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, x6) -
                            16. *  pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, x7) -
                            pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, x8) -
                            15. *  pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, x9) +
                            16. *  pow(a, 2) *  pow(b, 2) *  pow(p, 2) * _lisk->Li(2, x10) + 
                            30. *  pow(a, 2) * b * p * _lisk->Li(3, x6) +
                            30. * a *  pow(b, 2) * p * _lisk->Li(3, x6) - 
                            30. *  pow(a, 2) * b * p * _lisk->Li(3, x7) +
                            30. * a *  pow(b, 2) * p * _lisk->Li(3, x7) + 
                            30. *  pow(a, 2) * b * p * _lisk->Li(3, x9) -
                            30. * a *  pow(b, 2) * p * _lisk->Li(3, x9) - 
                            30. *  pow(a, 2) * b * p * _lisk->Li(3, x10) -
                            30. * a *  pow(b, 2) * p * _lisk->Li(3, x10) + 
                            30. *  pow(a, 2) * _lisk->Li(4, x6) +
                            30. *  pow(b, 2) * _lisk->Li(4, x6) - 
                            30. *  pow(a, 2) * _lisk->Li(4, x7) -
                            30. *  pow(b, 2) * _lisk->Li(4, x7) - 
                            30. *  pow(a, 2) * _lisk->Li(4, x9) -
                            30. *  pow(b, 2) * _lisk->Li(4, x9) + 
                            30. *  pow(a, 2) * _lisk->Li(4, x10) +
                            30. *  pow(b, 2) * _lisk->Li(4, x10)) / 
                            (8. *  pow(a, 3) *  pow(b, 3) *  pow(p, 4));
      
      _coefficients[1] = (3. * (a * p * _lisk->Li(2, x6) + 
                            a * p * _lisk->Li(2, x7) - a * p * 
                            _lisk->Li(2, x9) -
                            a * p * _lisk->Li(2, x10) + 
                            _lisk->Li(3, x6) - 
                            _lisk->Li(3, x7) -
                            _lisk->Li(3, x9) + 
                            _lisk->Li(3, x10))) / (2. * pow(a, 2) * b * pow(p, 3));
      
      
      _coefficients[2] = (-3. * (a * b * pow(p, 2) * log(1. - x6) + 
                            a * b * pow(p, 2) * log(1. - x7) +
                            a * b * pow(p, 2) * log(1. - x9) - 
                            a * b * pow(p, 2) *
                            log(-x6 * (1. - x4) * 
                            (1. - x9)) +
                            a * b * pow(p, 2) * log(1. - x10) - a * b * pow(p, 2) * 
                            log((1. - x7) * (1. - x10)) -
                            b * p * _lisk->Li(2, x6) + 
                            b * p * _lisk->Li(2, x7) - 
                            b * p * _lisk->Li(2, x9) +
                            b * p * _lisk->Li(2, x10) - 
                            _lisk->Li(3, x6) + 
                            _lisk->Li(3, x7) +
                            _lisk->Li(3, x9) - _lisk->Li(3, x10))) /
                            (2. * a * pow(b, 2) * pow(p, 3));
      
      
      _coefficients[3] = (-15. * (2. * pow(a, 3) * b * pow(p, 4) + 
                            pow(a, 2) * pow(p, 2) * _lisk->Li(2, x1) -
                            pow(a, 2) * pow(p, 2) * _lisk->Li(2, x2) - 
                            pow(a, 2) * pow(p, 2) * _lisk->Li(2, x3) +
                            pow(a, 2) * pow(p, 2) * _lisk->Li(2, x4) - 
                            pow(a, 2) * pow(p, 2) * _lisk->Li(2, x5) +
                            3. * pow(a, 2) * pow(p, 2) * _lisk->Li(2, x6) - 
                            2. * pow(a, 2) * pow(p, 2) * _lisk->Li(2, x7) +
                            pow(a, 2) * pow(p, 2) * _lisk->Li(2, x11) - 
                            3. * pow(a, 2) * pow(p, 2) * _lisk->Li(2, x9) +
                            2. * pow(a, 2) * pow(p, 2) * _lisk->Li(2, x10) + 
                            6. * a * p * _lisk->Li(3, x6) +
                            6. * a * p * _lisk->Li(3, x7) - 
                            6. * a * p * _lisk->Li(3, x9) - 
                            6. * a * p * _lisk->Li(3, x10) +
                            6. * _lisk->Li(4, x6) - 
                            6. * _lisk->Li(4, x7) - 
                            6. * _lisk->Li(4, x9) +
                            6. * _lisk->Li(4, x10))) / 
                            (8. * pow(a, 3) * b * pow(p, 4));
        
      
      _coefficients[4] = (15. * (2. * a * pow(b, 3) * pow(p, 4) + pow(b, 2) * pow(p, 2) * _lisk->Li(2, x1) -
                            pow(b, 2) * pow(p, 2) * _lisk->Li(2, x2) - 
                            pow(b, 2) * pow(p, 2) * _lisk->Li(2, x3) +
                            pow(b, 2) * pow(p, 2) * _lisk->Li(2, x4) - 
                            pow(b, 2) * pow(p, 2) * _lisk->Li(2, x5) +
                            pow(b, 2) * pow(p, 2) * _lisk->Li(2, x7) + 
                            pow(b, 2) * pow(p, 2) * _lisk->Li(2, x11) -
                            pow(b, 2) * pow(p, 2) * _lisk->Li(2, x10) - 
                            3. * b * p * _lisk->Li(3, x6) +
                            3. * b * p * _lisk->Li(3, x7) - 
                            3. * b * p * _lisk->Li(3, x9) + 
                            3. * b * p * _lisk->Li(3, x10) -
                            3. * _lisk->Li(4, x6) + 
                            3. * _lisk->Li(4, x7) + 
                            3. * _lisk->Li(4, x9) -
                            3. * _lisk->Li(4, x10))) / 
                            (4. * a * pow(b, 3) * pow(p, 4));
      
      
      _coefficients[5] = (-3. * (2. * pow(a, 3) * b * pow(p, 4) - 
                            3. * pow(a, 2) * b * pow(p, 3) * log(1. - x2) -
                            3. * pow(a, 2) * b * pow(p, 3) * log(1. - x4) + 
                            3. * pow(a, 2) * b * pow(p, 3) * log(1. - x7) +
                            3. * pow(a, 2) * b * pow(p, 3) * log(-x6 * 
                            (1. - x4) * (1. - x9)) +
                            3. * pow(a, 2) * b * pow(p, 3) * log(1. - x10) -
                            3. * pow(a, 2) * b * pow(p, 3) * log((1. - x7) * 
                            (1. - x10)) -
                            6. * a * b * pow(p, 2) * _lisk->Li(2, x2) - 
                            6. * a * b * pow(p, 2) * _lisk->Li(2, x4) +
                            6. * a * b * pow(p, 2) * _lisk->Li(2, x7) + 
                            6. * a * b * pow(p, 2) * _lisk->Li(2, x10) -
                            6. * b * p * _lisk->Li(3, x1) + 
                            6. * b * p * _lisk->Li(3, x2) - 
                            6. * b * p * _lisk->Li(3, x3) +
                            6. * b * p * _lisk->Li(3, x4) + 
                            6. * b * p * _lisk->Li(3, x5) + 
                            6. * a * p * _lisk->Li(3, x6) +
                            6. * a * p * _lisk->Li(3, x7) - 
                            6. * b * p * _lisk->Li(3, x7) + 
                            6. * b * p * _lisk->Li(3, x11) -
                            6. * a * p * _lisk->Li(3, x9) - 
                            6. * a * p * _lisk->Li(3, x10) - 
                            6. * b * p * _lisk->Li(3, x10) +
                            6. * _lisk->Li(4, x6) - 
                            6. * _lisk->Li(4, x7) - 6. * _lisk->Li(4, x9) +
                            6. * _lisk->Li(4, x10))) / 
                            (4. * pow(a, 2) * pow(b, 2) * pow(p, 4));

    }
    
    else if (degree == 3) {
        
       _coefficients[0] =  (-10. * pow(a,3) * pow(b,3) * pow(p,4) - pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(2,x1) + pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(2,x2) + 
                            pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(2,x3) - pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(2,x4) + 
                            pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(2,x5) + 15. * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(2,x6) - 
                            16. * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(2,x7) - pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(2,x8) - 
                            15. * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(2,x9) + 16. * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(2,x10) + 
                            30. * pow(a,2) * b * p * _lisk->Li(3,x6) + 30. * a * pow(b,2) * p * _lisk->Li(3,x6) - 30. * pow(a,2) * b * p * _lisk->Li(3,x7) + 
                            30. * a * pow(b,2) * p * _lisk->Li(3,x7) + 30. * pow(a,2) * b * p * _lisk->Li(3,x9) - 30. * a * pow(b,2) * p * _lisk->Li(3,x9) - 
                            30. * pow(a,2) * b * p * _lisk->Li(3,x10) - 30. * a * pow(b,2) * p * _lisk->Li(3,x10) + 30. * pow(a,2) * _lisk->Li(4,x6) + 
                            30. * pow(b,2) * _lisk->Li(4,x6) - 30. * pow(a,2) * _lisk->Li(4,x7) - 30. * pow(b,2) * _lisk->Li(4,x7) - 
                            30. * pow(a,2) * _lisk->Li(4,x9) - 30. * pow(b,2) * _lisk->Li(4,x9) + 30. * pow(a,2) * _lisk->Li(4,x10) + 
                            30. * pow(b,2) * _lisk->Li(4,x10))/(8. * pow(a,3) * pow(b,3) * pow(p,4));
                            
                            
       _coefficients[1] = (-15. * (2 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(2,x6) + 2 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(2,x7) - 
                            2. * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(2,x9) - 2 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(2,x10) + 
                            3. * pow(a,3) * b * pow(p,2) * _lisk->Li(3,x6) + 9. * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x6) + 
                            3. * pow(a,3) * b * pow(p,2) * _lisk->Li(3,x7) - 9. * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x7) + 
                            3. * pow(a,3) * b * pow(p,2) * _lisk->Li(3,x9) - 9. * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x9) + 
                            3. * pow(a,3) * b * pow(p,2) * _lisk->Li(3,x10) + 9. * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x10) + 
                            3. * pow(a,3) * p * _lisk->Li(4,x6) + 3. * pow(a,2) * b * p * _lisk->Li(4,x6) + 21 * a * pow(b,2) * p * _lisk->Li(4,x6) + 
                            3. * pow(a,3) * p * _lisk->Li(4,x7) - 3. * pow(a,2) * b * p * _lisk->Li(4,x7) + 21 * a * pow(b,2) * p * _lisk->Li(4,x7) - 
                            3.* pow(a,3) * p * _lisk->Li(4,x9) + 3. * pow(a,2) * b * p * _lisk->Li(4,x9) - 21 * a * pow(b,2) * p * _lisk->Li(4,x9) - 
                            3. * pow(a,3) * p * _lisk->Li(4,x10) - 3. * pow(a,2) * b * p * _lisk->Li(4,x10) - 21 * a * pow(b,2) * p * _lisk->Li(4,x10) + 
                            3. * pow(a,2) * _lisk->Li(5,x6) + 21. * pow(b,2) * _lisk->Li(5,x6) - 3. * pow(a,2) * _lisk->Li(5,x7) - 
                            21. * pow(b,2) * _lisk->Li(5,x7) - 3. * pow(a,2) * _lisk->Li(5,x9) - 21. * pow(b,2) * _lisk->Li(5,x9) + 
                            3. * pow(a,2) * _lisk->Li(5,x10) + 21. * pow(b,2) * _lisk->Li(5,x10)))/(4. * pow(a,4) * pow(b,3) * pow(p,5)); 
                            
       _coefficients[2] = (-15. * (pow(a,4) * pow(b,3) * pow(p,5) - 2 * pow(a,3) * pow(b,3) * pow(p,4) * log(1 - x2) - 2 * pow(a,3) * pow(b,3) * pow(p,4) * log(1 - x4) + 
                            12. * pow(a,3) * pow(b,3) * pow(p,4) * log(1 - x6) + 10 * pow(a,3) * pow(b,3) * pow(p,4) * log(1 - x7) + 
                            12. * pow(a,3) * pow(b,3) * pow(p,4) * log(1 - x9) - 
                            10. * pow(a,3) * pow(b,3) * pow(p,4) * log(-x6 * (1 - x4) * (1 - x9)) + 
                            10. * pow(a,3) * pow(b,3) * pow(p,4) * log(1 - x10) - 10 * pow(a,3) * pow(b,3) * pow(p,4) * log((1 - x7) * (1 - x10)) - 
                            6. * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(2,x2) - 6 * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(2,x4) + 
                            2. * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(2,x6) - 8 * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(2,x7) + 
                            2. * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(2,x9) - 8 * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(2,x10) + 
                            12. * a * pow(b,3) * pow(p,2) * _lisk->Li(3,x2) + 12 * a * pow(b,3) * pow(p,2) * _lisk->Li(3,x4) + 
                            36. * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x6) - 36 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x7) + 
                            12. * a * pow(b,3) * pow(p,2) * _lisk->Li(3,x7) - 36 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x9) + 
                            36. * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x10) + 12 * a * pow(b,3) * pow(p,2) * _lisk->Li(3,x10) + 12 * pow(b,3) * p * _lisk->Li(4,x1) - 
                            12. * pow(b,3) * p * _lisk->Li(4,x2) + 12. * pow(b,3) * p * _lisk->Li(4,x3) - 12 * pow(b,3) * p * _lisk->Li(4,x4) + 
                            12. * pow(b,3) * p * _lisk->Li(4,x5) + 84. * pow(a,2) * b * p * _lisk->Li(4,x6) + 12 * a * pow(b,2) * p * _lisk->Li(4,x6) - 
                            84. * pow(a,2) * b * p * _lisk->Li(4,x7) + 12. * a * pow(b,2) * p * _lisk->Li(4,x7) - 12 * pow(b,3) * p * _lisk->Li(4,x7) + 
                            12. * pow(b,3) * p * _lisk->Li(4,x8) + 84. * pow(a,2) * b * p * _lisk->Li(4,x9) - 12 * a * pow(b,2) * p * _lisk->Li(4,x9) - 
                            84. * pow(a,2) * b * p * _lisk->Li(4,x10) - 12 * a * pow(b,2) * p * _lisk->Li(4,x10) - 12 * pow(b,3) * p * _lisk->Li(4,x10) + 
                            84. * pow(a,2) * _lisk->Li(5,x6) + 12. * pow(b,2) * _lisk->Li(5,x6) - 84 * pow(a,2) * _lisk->Li(5,x7) - 
                            12. * pow(b,2) * _lisk->Li(5,x7) - 84. * pow(a,2) * _lisk->Li(5,x9) - 12 * pow(b,2) * _lisk->Li(5,x9) + 
                            84. * pow(a,2) * _lisk->Li(5,x10) + 12. * pow(b,2) * _lisk->Li(5,x10)))/(16. * pow(a,3) * pow(b,4) * pow(p,5));                     
      
       _coefficients[3] = (-15. * (2 * pow(a,3) * b * pow(p,4) + pow(a,2) * pow(p,2) * _lisk->Li(2,x1) - pow(a,2) * pow(p,2) * _lisk->Li(2,x2) - 
                            pow(a,2) * pow(p,2) * _lisk->Li(2,x3) + pow(a,2) * pow(p,2) * _lisk->Li(2,x4) - pow(a,2) * pow(p,2) * _lisk->Li(2,x5) + 
                            3 * pow(a,2) * pow(p,2) * _lisk->Li(2,x6) - 2 * pow(a,2) * pow(p,2) * _lisk->Li(2,x7) + pow(a,2) * pow(p,2) * _lisk->Li(2,x8) - 
                            3 * pow(a,2) * pow(p,2) * _lisk->Li(2,x9) + 2 * pow(a,2) * pow(p,2) * _lisk->Li(2,x10) + 6 * a * p * _lisk->Li(3,x6) + 
                            6 * a * p * _lisk->Li(3,x7) - 6 * a * p * _lisk->Li(3,x9) - 6 * a * p * _lisk->Li(3,x10) + 6. * _lisk->Li(4,x6) - 
                            6. * _lisk->Li(4,x7) - 6. * _lisk->Li(4,x9) + 6. * _lisk->Li(4,x10)))/(8. * pow(a,3) * b * pow(p,4));  
                            
                            
       _coefficients[4] = (15. * (2 * a * pow(b,3) * pow(p,4) + pow(b,2) * pow(p,2) * _lisk->Li(2,x1) - pow(b,2) * pow(p,2) * _lisk->Li(2,x2) - 
                            pow(b,2) * pow(p,2) * _lisk->Li(2,x3) + pow(b,2) * pow(p,2) * _lisk->Li(2,x4) - pow(b,2) * pow(p,2) * _lisk->Li(2,x5) + 
                            pow(b,2) * pow(p,2) * _lisk->Li(2,x7) + pow(b,2) * pow(p,2) * _lisk->Li(2,x8) - pow(b,2) * pow(p,2) * _lisk->Li(2,x10) - 
                            3 * b * p * _lisk->Li(3,x6) + 3 * b * p * _lisk->Li(3,x7) - 3 * b * p * _lisk->Li(3,x9) + 3 * b * p * _lisk->Li(3,x10) - 
                            3. * _lisk->Li(4,x6) + 3. * _lisk->Li(4,x7) + 3. * _lisk->Li(4,x9) - 3. * _lisk->Li(4,x10)))/(4. * a * pow(b,3) * pow(p,4));          
                            
                            
       _coefficients[5] = (35. * (pow(a,3) * pow(p,3) * _lisk->Li(2,x6) + pow(a,3) * pow(p,3) * _lisk->Li(2,x7) - pow(a,3) * pow(p,3) * _lisk->Li(2,x9) - 
                            pow(a,3) * pow(p,3) * _lisk->Li(2,x10) + 6 * pow(a,2) * pow(p,2) * _lisk->Li(3,x6) - 6 * pow(a,2) * pow(p,2) * _lisk->Li(3,x7) - 
                            6 * pow(a,2) * pow(p,2) * _lisk->Li(3,x9) + 6 * pow(a,2) * pow(p,2) * _lisk->Li(3,x10) + 15 * a * p * _lisk->Li(4,x6) + 
                            15 * a * p * _lisk->Li(4,x7) - 15 * a * p * _lisk->Li(4,x9) - 15 * a * p * _lisk->Li(4,x10) + 15. * _lisk->Li(5,x6) - 
                            15. * _lisk->Li(5,x7) - 15. * _lisk->Li(5,x9) + 15. * _lisk->Li(5,x10)))/(4. * pow(a,4) * b * pow(p,5));      
                            
                            
       _coefficients[6] = (35. * (3 * a * pow(b,3) * pow(p,4) * log(1 - x6) + 3 * a * pow(b,3) * pow(p,4) * log(1 - x7) + 3 * a * pow(b,3) * pow(p,4) * log(1 - x9) - 
                            3 * a * pow(b,3) * pow(p,4) * log(-x6 * (1 - x4) * (1 - x9)) + 3 * a * pow(b,3) * pow(p,4) * log(1 - x10) - 
                            3 * a * pow(b,3) * pow(p,4) * log((1 - x7) * (1 - x10)) + 2 * pow(b,3) * pow(p,3) * _lisk->Li(2,x6) - 
                            2 * pow(b,3) * pow(p,3) * _lisk->Li(2,x7) + 2 * pow(b,3) * pow(p,3) * _lisk->Li(2,x9) - 2 * pow(b,3) * pow(p,3) * _lisk->Li(2,x10) + 
                            12 * pow(b,2) * pow(p,2) * _lisk->Li(3,x6) - 12 * pow(b,2) * pow(p,2) * _lisk->Li(3,x7) - 
                            12 * pow(b,2) * pow(p,2) * _lisk->Li(3,x9) + 12 * pow(b,2) * pow(p,2) * _lisk->Li(3,x10) + 30 * b * p * _lisk->Li(4,x6) - 
                            30 * b * p * _lisk->Li(4,x7) + 30 * b * p * _lisk->Li(4,x9) - 30 * b * p * _lisk->Li(4,x10) + 30. * _lisk->Li(5,x6) - 
                            30. * _lisk->Li(5,x7) - 30. * _lisk->Li(5,x9) + 30. * _lisk->Li(5,x10)))/(8. * a * pow(b,4) * pow(p,5));
                            
                            
       _coefficients[7] = (45. * (pow(a,4) * b * pow(p,5) - 2 * pow(a,3) * b * pow(p,4) * log(1 - x2) - 2 * pow(a,3) * b * pow(p,4) * log(1 - x4) + 
                            2 * pow(a,3) * b * pow(p,4) * log(1 - x6) + 2 * pow(a,3) * b * pow(p,4) * log(1 - x9) - 6 * pow(a,2) * b * pow(p,3) * _lisk->Li(2,x2) - 
                            6 * pow(a,2) * b * pow(p,3) * _lisk->Li(2,x4) - 2 * pow(a,2) * b * pow(p,3) * _lisk->Li(2,x6) - 
                            4 * pow(a,2) * b * pow(p,3) * _lisk->Li(2,x7) - 2 * pow(a,2) * b * pow(p,3) * _lisk->Li(2,x9) - 
                            4 * pow(a,2) * b * pow(p,3) * _lisk->Li(2,x10) + 12 * a * b * pow(p,2) * _lisk->Li(3,x2) + 12 * a * b * pow(p,2) * _lisk->Li(3,x4) + 
                            4 * pow(a,2) * pow(p,2) * _lisk->Li(3,x6) - 4 * pow(a,2) * pow(p,2) * _lisk->Li(3,x7) + 12 * a * b * pow(p,2) * _lisk->Li(3,x7) - 
                            4 * pow(a,2) * pow(p,2) * _lisk->Li(3,x9) + 4 * pow(a,2) * pow(p,2) * _lisk->Li(3,x10) + 12 * a * b * pow(p,2) * _lisk->Li(3,x10) + 
                            12 * b * p * _lisk->Li(4,x1) - 12 * b * p * _lisk->Li(4,x2) + 12 * b * p * _lisk->Li(4,x3) - 12 * b * p * _lisk->Li(4,x4) + 
                            12 * b * p * _lisk->Li(4,x5) + 12 * a * p * _lisk->Li(4,x6) + 12 * a * p * _lisk->Li(4,x7) - 12 * b * p * _lisk->Li(4,x7) + 
                            12 * b * p * _lisk->Li(4,x8) - 12 * a * p * _lisk->Li(4,x9) - 12 * a * p * _lisk->Li(4,x10) - 12 * b * p * _lisk->Li(4,x10) + 
                            12. * _lisk->Li(5,x6) - 12. * _lisk->Li(5,x7) - 12. * _lisk->Li(5,x9) + 12. * _lisk->Li(5,x10)))/(16. * pow(a,3) * pow(b,2) * pow(p,5));            
        
        
       _coefficients[8] = (45. * (a * pow(b,2) * pow(p,3) * _lisk->Li(2,x6) + a * pow(b,2) * pow(p,3) * _lisk->Li(2,x7) - 
                            a * pow(b,2) * pow(p,3) * _lisk->Li(2,x9) - a * pow(b,2) * pow(p,3) * _lisk->Li(2,x10) + 3 * a * b * pow(p,2) * _lisk->Li(3,x6) + 
                            pow(b,2) * pow(p,2) * _lisk->Li(3,x6) + 3 * a * b * pow(p,2) * _lisk->Li(3,x7) - pow(b,2) * pow(p,2) * _lisk->Li(3,x7) + 
                            3 * a * b * pow(p,2) * _lisk->Li(3,x9) - pow(b,2) * pow(p,2) * _lisk->Li(3,x9) + 3 * a * b * pow(p,2) * _lisk->Li(3,x10) + 
                            pow(b,2) * pow(p,2) * _lisk->Li(3,x10) + 3 * a * p * _lisk->Li(4,x6) + 3 * b * p * _lisk->Li(4,x6) + 
                            3 * a * p * _lisk->Li(4,x7) - 3 * b * p * _lisk->Li(4,x7) - 3 * a * p * _lisk->Li(4,x9) + 3 * b * p * _lisk->Li(4,x9) - 
                            3 * a * p * _lisk->Li(4,x10) - 3 * b * p * _lisk->Li(4,x10) + 3. * _lisk->Li(5,x6) - 3. * _lisk->Li(5,x7) - 
                            3. * _lisk->Li(5,x9) + 3. * _lisk->Li(5,x10)))/(4. * pow(a,2) * pow(b,3) * pow(p,5)); 
                            
                            
                            
       _coefficients[9] = (-3. * (2 * pow(a,3) * b * pow(p,4) - 3 * pow(a,2) * b * pow(p,3) * log(1 - x2) - 3 * pow(a,2) * b * pow(p,3) * log(1 - x4) + 
                            3 * pow(a,2) * b * pow(p,3) * log(1 - x7) + 3 * pow(a,2) * b * pow(p,3) * log(-x6 * (1 - x4) * (1 - x9)) + 
                            3 * pow(a,2) * b * pow(p,3) * log(1 - x10) - 3 * pow(a,2) * b * pow(p,3) * log((1 - x7) * (1 - x10)) - 
                            6 * a * b * pow(p,2) * _lisk->Li(2,x2) - 6 * a * b * pow(p,2) * _lisk->Li(2,x4) + 6 * a * b * pow(p,2) * _lisk->Li(2,x7) + 
                            6 * a * b * pow(p,2) * _lisk->Li(2,x10) - 6 * b * p * _lisk->Li(3,x1) + 6 * b * p * _lisk->Li(3,x2) - 6 * b * p * _lisk->Li(3,x3) + 
                            6 * b * p * _lisk->Li(3,x4) + 6 * b * p * _lisk->Li(3,x5) + 6 * a * p * _lisk->Li(3,x6) + 6 * a * p * _lisk->Li(3,x7) - 
                            6 * b * p * _lisk->Li(3,x7) + 6 * b * p * _lisk->Li(3,x8) - 6 * a * p * _lisk->Li(3,x9) - 6 * a * p * _lisk->Li(3,x10) - 
                            6 * b * p * _lisk->Li(3,x10) + 6. * _lisk->Li(4,x6) - 6. * _lisk->Li(4,x7) - 6. * _lisk->Li(4,x9) + 
                            6. * _lisk->Li(4,x10)))/(4. * pow(a,2) * pow(b,2) * pow(p,4)); 
        
    }
    
    
    else if (_degree == 4) {
     
       _coefficients[0] = (-75. * pow(a,6) * pow(b,4) * pow(p,6) - 1856 * pow(a,5) * pow(b,5) * pow(p,6) - 300 * pow(a,5) * pow(b,4) * pow(p,5) * log(1 - x12) - 
                            300 * pow(a,5) * pow(b,4) * pow(p,5) * log(1 - x2) + 300 * pow(a,5) * pow(b,4) * pow(p,5) * log(1 - x4) + 
                            300 * pow(a,5) * pow(b,4) * pow(p,5) * log(1 - x7) + 300 * pow(a,5) * pow(b,4) * pow(p,5) * log((1 - x10)/(1 - x7)) + 
                            300 * pow(a,5) * pow(b,4) * pow(p,5) * log((-x13 + exp((2 * b + d) * p))/(-x13 + exp(d * p))) - 
                            864 * pow(a,4) * pow(b,4) * pow(p,4) * _lisk->Li(2,x1) + 900 * pow(a,4) * pow(b,4) * pow(p,4) * _lisk->Li(2,x12) - 
                            36 * pow(a,4) * pow(b,4) * pow(p,4) * _lisk->Li(2,x2) + 864 * pow(a,4) * pow(b,4) * pow(p,4) * _lisk->Li(2,x3) + 
                            36 * pow(a,4) * pow(b,4) * pow(p,4) * _lisk->Li(2,x4) + 864 * pow(a,4) * pow(b,4) * pow(p,4) * _lisk->Li(2,x5) - 
                            540 * pow(a,4) * pow(b,4) * pow(p,4) * _lisk->Li(2,x6) + 576 * pow(a,4) * pow(b,4) * pow(p,4) * _lisk->Li(2,x7) - 
                            864 * pow(a,4) * pow(b,4) * pow(p,4) * _lisk->Li(2,x8) + 540 * pow(a,4) * pow(b,4) * pow(p,4) * _lisk->Li(2,x9) + 
                            324 * pow(a,4) * pow(b,4) * pow(p,4) * _lisk->Li(2,x10) + 1800 * pow(a,3) * pow(b,4) * pow(p,3) * _lisk->Li(3,x12) + 
                            1800 * pow(a,3) * pow(b,4) * pow(p,3) * _lisk->Li(3,x2) - 1800 * pow(a,3) * pow(b,4) * pow(p,3) * _lisk->Li(3,x4) - 
                            2880 * pow(a,4) * pow(b,3) * pow(p,3) * _lisk->Li(3,x6) - 1080 * pow(a,3) * pow(b,4) * pow(p,3) * _lisk->Li(3,x6) + 
                            2880 * pow(a,4) * pow(b,3) * pow(p,3) * _lisk->Li(3,x7) - 2880 * pow(a,3) * pow(b,4) * pow(p,3) * _lisk->Li(3,x7) - 
                            2880 * pow(a,4) * pow(b,3) * pow(p,3) * _lisk->Li(3,x9) + 1080 * pow(a,3) * pow(b,4) * pow(p,3) * _lisk->Li(3,x9) + 
                            2880 * pow(a,4) * pow(b,3) * pow(p,3) * _lisk->Li(3,x10) + 1080 * pow(a,3) * pow(b,4) * pow(p,3) * _lisk->Li(3,x10) + 
                            1800 * pow(a,2) * pow(b,4) * pow(p,2) * _lisk->Li(4,x12) - 1800 * pow(a,2) * pow(b,4) * pow(p,2) * _lisk->Li(4,x2) - 
                            1800 * pow(a,2) * pow(b,4) * pow(p,2) * _lisk->Li(4,x3) + 1800 * pow(a,2) * pow(b,4) * pow(p,2) * _lisk->Li(4,x4) - 
                            1800 * pow(a,2) * pow(b,4) * pow(p,2) * _lisk->Li(4,x5) - 10440 * pow(a,4) * pow(b,2) * pow(p,2) * _lisk->Li(4,x6) - 
                            3600 * pow(a,3) * pow(b,3) * pow(p,2) * _lisk->Li(4,x6) - 8640 * pow(a,2) * pow(b,4) * pow(p,2) * _lisk->Li(4,x6) + 
                            10440 * pow(a,4) * pow(b,2) * pow(p,2) * _lisk->Li(4,x7) - 3600 * pow(a,3) * pow(b,3) * pow(p,2) * _lisk->Li(4,x7) + 
                            10440 * pow(a,2) * pow(b,4) * pow(p,2) * _lisk->Li(4,x7) + 10440 * pow(a,4) * pow(b,2) * pow(p,2) * _lisk->Li(4,x9) - 
                            3600 * pow(a,3) * pow(b,3) * pow(p,2) * _lisk->Li(4,x9) + 8640 * pow(a,2) * pow(b,4) * pow(p,2) * _lisk->Li(4,x9) - 
                            10440 * pow(a,4) * pow(b,2) * pow(p,2) * _lisk->Li(4,x10) - 3600 * pow(a,3) * pow(b,3) * pow(p,2) * _lisk->Li(4,x10) - 
                            8640 * pow(a,2) * pow(b,4) * pow(p,2) * _lisk->Li(4,x10) - 22680 * pow(a,4) * b * p * _lisk->Li(5,x6) - 
                            3600 * pow(a,3) * pow(b,2) * p * _lisk->Li(5,x6) - 3600 * pow(a,2) * pow(b,3) * p * _lisk->Li(5,x6) - 
                            22680 * a * pow(b,4) * p * _lisk->Li(5,x6) + 22680 * pow(a,4) * b * p * _lisk->Li(5,x7) - 3600 * pow(a,3) * pow(b,2) * p * _lisk->Li(5,x7) + 
                            3600 * pow(a,2) * pow(b,3) * p * _lisk->Li(5,x7) - 22680 * a * pow(b,4) * p * _lisk->Li(5,x7) - 22680 * pow(a,4) * b * p * _lisk->Li(5,x9) + 
                            3600 * pow(a,3) * pow(b,2) * p * _lisk->Li(5,x9) - 3600 * pow(a,2) * pow(b,3) * p * _lisk->Li(5,x9) + 
                            22680 * a * pow(b,4) * p * _lisk->Li(5,x9) + 22680 * pow(a,4) * b * p * _lisk->Li(5,x10) + 3600 * pow(a,3) * pow(b,2) * p * _lisk->Li(5,x10) + 
                            3600 * pow(a,2) * pow(b,3) * p * _lisk->Li(5,x10) + 22680 * a * pow(b,4) * p * _lisk->Li(5,x10) - 22680 * pow(a,4) * _lisk->Li(6,x6) - 
                            3600 * pow(a,2) * pow(b,2) * _lisk->Li(6,x6) - 22680 * pow(b,4) * _lisk->Li(6,x6) + 22680 * pow(a,4) * _lisk->Li(6,x7) + 
                            3600 * pow(a,2) * pow(b,2) * _lisk->Li(6,x7) + 22680 * pow(b,4) * _lisk->Li(6,x7) + 22680 * pow(a,4) * _lisk->Li(6,x9) + 
                            3600 * pow(a,2) * pow(b,2) * _lisk->Li(6,x9) + 22680 * pow(b,4) * _lisk->Li(6,x9) - 22680 * pow(a,4) * _lisk->Li(6,x10) - 
                            3600 * pow(a,2) * pow(b,2) * _lisk->Li(6,x10) - 22680 * pow(b,4) * _lisk->Li(6,x10))/(128. * pow(a,5) * pow(b,5) * pow(p,6));
       
       _coefficients[1] = (-15. * (2 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(2,x6) + 2 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(2,x7) - 
                            2 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(2,x9) - 2 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(2,x10) + 
                            3 * pow(a,3) * b * pow(p,2) * _lisk->Li(3,x6) + 9 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x6) + 
                            3 * pow(a,3) * b * pow(p,2) * _lisk->Li(3,x7) - 9 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x7) + 
                            3 * pow(a,3) * b * pow(p,2) * _lisk->Li(3,x9) - 9 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x9) + 
                            3 * pow(a,3) * b * pow(p,2) * _lisk->Li(3,x10) + 9 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x10) + 
                            3 * pow(a,3) * p * _lisk->Li(4,x6) + 3 * pow(a,2) * b * p * _lisk->Li(4,x6) + 21 * a * pow(b,2) * p * _lisk->Li(4,x6) + 
                            3 * pow(a,3) * p * _lisk->Li(4,x7) - 3 * pow(a,2) * b * p * _lisk->Li(4,x7) + 21 * a * pow(b,2) * p * _lisk->Li(4,x7) - 
                            3 * pow(a,3) * p * _lisk->Li(4,x9) + 3 * pow(a,2) * b * p * _lisk->Li(4,x9) - 21 * a * pow(b,2) * p * _lisk->Li(4,x9) - 
                            3 * pow(a,3) * p * _lisk->Li(4,x10) - 3 * pow(a,2) * b * p * _lisk->Li(4,x10) - 21 * a * pow(b,2) * p * _lisk->Li(4,x10) + 
                            3 * pow(a,2) * _lisk->Li(5,x6) + 21 * pow(b,2) * _lisk->Li(5,x6) - 3 * pow(a,2) * _lisk->Li(5,x7) - 
                            21 * pow(b,2) * _lisk->Li(5,x7) - 3 * pow(a,2) * _lisk->Li(5,x9) - 21 * pow(b,2) * _lisk->Li(5,x9) + 
                            3 * pow(a,2) * _lisk->Li(5,x10) + 21 * pow(b,2) * _lisk->Li(5,x10)))/(4. * pow(a,4) * pow(b,3) * pow(p,5));
       
       _coefficients[2] = (-15. * (pow(a,4) * pow(b,3) * pow(p,5) - 2 * pow(a,3) * pow(b,3) * pow(p,4) * log(1 - x2) - 2 * pow(a,3) * pow(b,3) * pow(p,4) * log(1 - x4) + 
                            12 * pow(a,3) * pow(b,3) * pow(p,4) * log(1 - x6) + 10 * pow(a,3) * pow(b,3) * pow(p,4) * log(1 - x7) + 
                            12 * pow(a,3) * pow(b,3) * pow(p,4) * log(1 - x9) - 
                            10 * pow(a,3) * pow(b,3) * pow(p,4) * log(-x6 * (1 - x4) * (1 - x9)) + 
                            10 * pow(a,3) * pow(b,3) * pow(p,4) * log(1 - x10) - 10 * pow(a,3) * pow(b,3) * pow(p,4) * log((1 - x7) * (1 - x10)) - 
                            6 * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(2,x2) - 6 * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(2,x4) + 
                            2 * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(2,x6) - 8 * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(2,x7) + 
                            2 * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(2,x9) - 8 * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(2,x10) + 
                            12 * a * pow(b,3) * pow(p,2) * _lisk->Li(3,x2) + 12 * a * pow(b,3) * pow(p,2) * _lisk->Li(3,x4) + 
                            36 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x6) - 36 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x7) + 
                            12 * a * pow(b,3) * pow(p,2) * _lisk->Li(3,x7) - 36 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x9) + 
                            36 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(3,x10) + 12 * a * pow(b,3) * pow(p,2) * _lisk->Li(3,x10) + 12 * pow(b,3) * p * _lisk->Li(4,x1) - 
                            12 * pow(b,3) * p * _lisk->Li(4,x2) + 12 * pow(b,3) * p * _lisk->Li(4,x3) - 12 * pow(b,3) * p * _lisk->Li(4,x4) + 
                            12 * pow(b,3) * p * _lisk->Li(4,x5) + 84 * pow(a,2) * b * p * _lisk->Li(4,x6) + 12 * a * pow(b,2) * p * _lisk->Li(4,x6) - 
                            84 * pow(a,2) * b * p * _lisk->Li(4,x7) + 12 * a * pow(b,2) * p * _lisk->Li(4,x7) - 12 * pow(b,3) * p * _lisk->Li(4,x7) + 
                            12 * pow(b,3) * p * _lisk->Li(4,x8) + 84 * pow(a,2) * b * p * _lisk->Li(4,x9) - 12 * a * pow(b,2) * p * _lisk->Li(4,x9) - 
                            84 * pow(a,2) * b * p * _lisk->Li(4,x10) - 12 * a * pow(b,2) * p * _lisk->Li(4,x10) - 12 * pow(b,3) * p * _lisk->Li(4,x10) + 
                            84 * pow(a,2) * _lisk->Li(5,x6) + 12 * pow(b,2) * _lisk->Li(5,x6) - 84 * pow(a,2) * _lisk->Li(5,x7) - 
                            12 * pow(b,2) * _lisk->Li(5,x7) - 84 * pow(a,2) * _lisk->Li(5,x9) - 12 * pow(b,2) * _lisk->Li(5,x9) + 
                            84 * pow(a,2) * _lisk->Li(5,x10) + 12 * pow(b,2) * _lisk->Li(5,x10)))/(16. * pow(a,3) * pow(b,4) * pow(p,5));
       
       _coefficients[3] = (45. * (5 * pow(a,6) * pow(b,2) * pow(p,6) - 20 * pow(a,5) * pow(b,3) * pow(p,6) + 20 * pow(a,5) * pow(b,2) * pow(p,5) * log(1 - x12) + 
                            20 * pow(a,5) * pow(b,2) * pow(p,5) * log(1 - x2) - 20 * pow(a,5) * pow(b,2) * pow(p,5) * log(1 - x4) - 
                            20 * pow(a,5) * pow(b,2) * pow(p,5) * log(1 - x7) - 20 * pow(a,5) * pow(b,2) * pow(p,5) * log((1 - x10)/(1 - x7)) - 
                            20 * pow(a,5) * pow(b,2) * pow(p,5) * log((-x13 + exp((2 * b + d) * p))/(-x13 + exp(d * p))) - 
                            10 * pow(a,4) * pow(b,2) * pow(p,4) * _lisk->Li(2,x1) - 60 * pow(a,4) * pow(b,2) * pow(p,4) * _lisk->Li(2,x12) + 
                            70 * pow(a,4) * pow(b,2) * pow(p,4) * _lisk->Li(2,x2) + 10 * pow(a,4) * pow(b,2) * pow(p,4) * _lisk->Li(2,x3) - 
                            70 * pow(a,4) * pow(b,2) * pow(p,4) * _lisk->Li(2,x4) + 10 * pow(a,4) * pow(b,2) * pow(p,4) * _lisk->Li(2,x5) - 
                            6 * pow(a,4) * pow(b,2) * pow(p,4) * _lisk->Li(2,x6) - 64 * pow(a,4) * pow(b,2) * pow(p,4) * _lisk->Li(2,x7) - 
                            10 * pow(a,4) * pow(b,2) * pow(p,4) * _lisk->Li(2,x8) + 6 * pow(a,4) * pow(b,2) * pow(p,4) * _lisk->Li(2,x9) + 
                            4 * pow(a,4) * pow(b,2) * pow(p,4) * _lisk->Li(2,x10) - 120 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(3,x12) - 
                            120 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(3,x2) + 120 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(3,x4) + 
                            80 * pow(a,4) * b * pow(p,3) * _lisk->Li(3,x6) + 408 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(3,x6) - 
                            80 * pow(a,4) * b * pow(p,3) * _lisk->Li(3,x7) + 528 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(3,x7) + 
                            80 * pow(a,4) * b * pow(p,3) * _lisk->Li(3,x9) - 408 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(3,x9) - 
                            80 * pow(a,4) * b * pow(p,3) * _lisk->Li(3,x10) - 408 * pow(a,3) * pow(b,2) * pow(p,3) * _lisk->Li(3,x10) - 
                            120 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(4,x12) + 120 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(4,x2) + 
                            120 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(4,x3) - 120 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(4,x4) + 
                            120 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(4,x5) + 80 * pow(a,4) * pow(p,2) * _lisk->Li(4,x6) + 
                            240 * pow(a,3) * b * pow(p,2) * _lisk->Li(4,x6) + 2088 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(4,x6) - 
                            80 * pow(a,4) * pow(p,2) * _lisk->Li(4,x7) + 240 * pow(a,3) * b * pow(p,2) * _lisk->Li(4,x7) - 
                            2208 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(4,x7) - 80 * pow(a,4) * pow(p,2) * _lisk->Li(4,x9) + 
                            240 * pow(a,3) * b * pow(p,2) * _lisk->Li(4,x9) - 2088 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(4,x9) + 
                            80 * pow(a,4) * pow(p,2) * _lisk->Li(4,x10) + 240 * pow(a,3) * b * pow(p,2) * _lisk->Li(4,x10) + 
                            2088 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(4,x10) + 240 * pow(a,3) * p * _lisk->Li(5,x6) + 240 * pow(a,2) * b * p * _lisk->Li(5,x6) + 
                            5040 * a * pow(b,2) * p * _lisk->Li(5,x6) + 240 * pow(a,3) * p * _lisk->Li(5,x7) - 240 * pow(a,2) * b * p * _lisk->Li(5,x7) + 
                            5040 * a * pow(b,2) * p * _lisk->Li(5,x7) - 240 * pow(a,3) * p * _lisk->Li(5,x9) + 240 * pow(a,2) * b * p * _lisk->Li(5,x9) - 
                            5040 * a * pow(b,2) * p * _lisk->Li(5,x9) - 240 * pow(a,3) * p * _lisk->Li(5,x10) - 240 * pow(a,2) * b * p * _lisk->Li(5,x10) - 
                            5040 * a * pow(b,2) * p * _lisk->Li(5,x10) + 240 * pow(a,2) * _lisk->Li(6,x6) + 5040 * pow(b,2) * _lisk->Li(6,x6) - 
                            240 * pow(a,2) * _lisk->Li(6,x7) - 5040 * pow(b,2) * _lisk->Li(6,x7) - 240 * pow(a,2) * _lisk->Li(6,x9) - 
                            5040 * pow(b,2) * _lisk->Li(6,x9) + 240 * pow(a,2) * _lisk->Li(6,x10) + 5040 * pow(b,2) * _lisk->Li(6,x10)))/
                        (128. * pow(a,5) * pow(b,3) * pow(p,6));
       
       _coefficients[4] = (45. * (5 * pow(a,4) * pow(b,4) * pow(p,6) + 372 * pow(a,3) * pow(b,5) * pow(p,6) + 
                            20 * pow(a,3) * pow(b,4) * pow(p,5) * log(1 - x12) + 20 * pow(a,3) * pow(b,4) * pow(p,5) * log(1 - x2) - 
                            20 * pow(a,3) * pow(b,4) * pow(p,5) * log(1 - x4) - 20 * pow(a,3) * pow(b,4) * pow(p,5) * log(1 - x7) - 
                            20 * pow(a,3) * pow(b,4) * pow(p,5) * log((1 - x10)/(1 - x7)) - 
                            20 * pow(a,3) * pow(b,4) * pow(p,5) * log((-x13 + exp((2 * b + d) * p))/(-x13 + exp(d * p))) + 
                            186 * pow(a,2) * pow(b,4) * pow(p,4) * _lisk->Li(2,x1) - 60 * pow(a,2) * pow(b,4) * pow(p,4) * _lisk->Li(2,x12) - 
                            126 * pow(a,2) * pow(b,4) * pow(p,4) * _lisk->Li(2,x2) - 186 * pow(a,2) * pow(b,4) * pow(p,4) * _lisk->Li(2,x3) + 
                            126 * pow(a,2) * pow(b,4) * pow(p,4) * _lisk->Li(2,x4) - 186 * pow(a,2) * pow(b,4) * pow(p,4) * _lisk->Li(2,x5) + 
                            190 * pow(a,2) * pow(b,4) * pow(p,4) * _lisk->Li(2,x6) - 64 * pow(a,2) * pow(b,4) * pow(p,4) * _lisk->Li(2,x7) + 
                            186 * pow(a,2) * pow(b,4) * pow(p,4) * _lisk->Li(2,x8) - 190 * pow(a,2) * pow(b,4) * pow(p,4) * _lisk->Li(2,x9) + 
                            4 * pow(a,2) * pow(b,4) * pow(p,4) * _lisk->Li(2,x10) - 120 * a * pow(b,4) * pow(p,3) * _lisk->Li(3,x12) - 
                            120 * a * pow(b,4) * pow(p,3) * _lisk->Li(3,x2) + 120 * a * pow(b,4) * pow(p,3) * _lisk->Li(3,x4) + 
                            528 * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(3,x6) - 40 * a * pow(b,4) * pow(p,3) * _lisk->Li(3,x6) - 
                            528 * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(3,x7) + 80 * a * pow(b,4) * pow(p,3) * _lisk->Li(3,x7) + 
                            528 * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(3,x9) + 40 * a * pow(b,4) * pow(p,3) * _lisk->Li(3,x9) - 
                            528 * pow(a,2) * pow(b,3) * pow(p,3) * _lisk->Li(3,x10) + 40 * a * pow(b,4) * pow(p,3) * _lisk->Li(3,x10) - 
                            120 * pow(b,4) * pow(p,2) * _lisk->Li(4,x12) + 120 * pow(b,4) * pow(p,2) * _lisk->Li(4,x2) + 120 * pow(b,4) * pow(p,2) * _lisk->Li(4,x3) - 
                            120 * pow(b,4) * pow(p,2) * _lisk->Li(4,x4) + 120 * pow(b,4) * pow(p,2) * _lisk->Li(4,x5) + 
                            2208 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(4,x6) + 240 * a * pow(b,3) * pow(p,2) * _lisk->Li(4,x6) - 
                            40 * pow(b,4) * pow(p,2) * _lisk->Li(4,x6) - 2208 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(4,x7) + 
                            240 * a * pow(b,3) * pow(p,2) * _lisk->Li(4,x7) - 80 * pow(b,4) * pow(p,2) * _lisk->Li(4,x7) - 
                            2208 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(4,x9) + 240 * a * pow(b,3) * pow(p,2) * _lisk->Li(4,x9) + 
                            40 * pow(b,4) * pow(p,2) * _lisk->Li(4,x9) + 2208 * pow(a,2) * pow(b,2) * pow(p,2) * _lisk->Li(4,x10) + 
                            240 * a * pow(b,3) * pow(p,2) * _lisk->Li(4,x10) - 40 * pow(b,4) * pow(p,2) * _lisk->Li(4,x10) + 5040 * pow(a,2) * b * p * _lisk->Li(5,x6) + 
                            240 * a * pow(b,2) * p * _lisk->Li(5,x6) + 240 * pow(b,3) * p * _lisk->Li(5,x6) - 5040 * pow(a,2) * b * p * _lisk->Li(5,x7) + 
                            240 * a * pow(b,2) * p * _lisk->Li(5,x7) - 240 * pow(b,3) * p * _lisk->Li(5,x7) + 5040 * pow(a,2) * b * p * _lisk->Li(5,x9) - 
                            240 * a * pow(b,2) * p * _lisk->Li(5,x9) + 240 * pow(b,3) * p * _lisk->Li(5,x9) - 5040 * pow(a,2) * b * p * _lisk->Li(5,x10) - 
                            240 * a * pow(b,2) * p * _lisk->Li(5,x10) - 240 * pow(b,3) * p * _lisk->Li(5,x10) + 5040 * pow(a,2) * _lisk->Li(6,x6) + 
                            240 * pow(b,2) * _lisk->Li(6,x6) - 5040 * pow(a,2) * _lisk->Li(6,x7) - 240 * pow(b,2) * _lisk->Li(6,x7) - 
                            5040 * pow(a,2) * _lisk->Li(6,x9) - 240 * pow(b,2) * _lisk->Li(6,x9) + 5040 * pow(a,2) * _lisk->Li(6,x10) + 
                            240 * pow(b,2) * _lisk->Li(6,x10)))/(128. * pow(a,3) * pow(b,5) * pow(p,6));
       
       _coefficients[5] = (35. * (pow(a,3) * pow(p,3) * _lisk->Li(2,x6) + pow(a,3) * pow(p,3) * _lisk->Li(2,x7) - pow(a,3) * pow(p,3) * _lisk->Li(2,x9) - 
                            pow(a,3) * pow(p,3) * _lisk->Li(2,x10) + 6 * pow(a,2) * pow(p,2) * _lisk->Li(3,x6) - 6 * pow(a,2) * pow(p,2) * _lisk->Li(3,x7) - 
                            6 * pow(a,2) * pow(p,2) * _lisk->Li(3,x9) + 6 * pow(a,2) * pow(p,2) * _lisk->Li(3,x10) + 15 * a * p * _lisk->Li(4,x6) + 
                            15 * a * p * _lisk->Li(4,x7) - 15 * a * p * _lisk->Li(4,x9) - 15 * a * p * _lisk->Li(4,x10) + 15. * _lisk->Li(5,x6) - 
                            15. * _lisk->Li(5,x7) - 15. * _lisk->Li(5,x9) + 15. * _lisk->Li(5,x10)))/(4. * pow(a,4) * b * pow(p,5));
       
       _coefficients[6] = (35. * (3 * a * pow(b,3) * pow(p,4) * log(1 - x6) + 3 * a * pow(b,3) * pow(p,4) * log(1 - x7) + 3 * a * pow(b,3) * pow(p,4) * log(1 - x9) - 
                            3 * a * pow(b,3) * pow(p,4) * log(-x6 * (1 - x4) * (1 - x9)) + 3 * a * pow(b,3) * pow(p,4) * log(1 - x10) - 
                            3 * a * pow(b,3) * pow(p,4) * log((1 - x7) * (1 - x10)) + 2 * pow(b,3) * pow(p,3) * _lisk->Li(2,x6) - 
                            2 * pow(b,3) * pow(p,3) * _lisk->Li(2,x7) + 2 * pow(b,3) * pow(p,3) * _lisk->Li(2,x9) - 2 * pow(b,3) * pow(p,3) * _lisk->Li(2,x10) + 
                            12 * pow(b,2) * pow(p,2) * _lisk->Li(3,x6) - 12 * pow(b,2) * pow(p,2) * _lisk->Li(3,x7) - 
                            12 * pow(b,2) * pow(p,2) * _lisk->Li(3,x9) + 12 * pow(b,2) * pow(p,2) * _lisk->Li(3,x10) + 30 * b * p * _lisk->Li(4,x6) - 
                            30 * b * p * _lisk->Li(4,x7) + 30 * b * p * _lisk->Li(4,x9) - 30 * b * p * _lisk->Li(4,x10) + 30. * _lisk->Li(5,x6) - 
                            30. * _lisk->Li(5,x7) - 30. * _lisk->Li(5,x9) + 30. * _lisk->Li(5,x10)))/(8. * a * pow(b,4) * pow(p,5));
       
       _coefficients[7] = (315. * (6 * pow(a,5) * b * pow(p,6) + 3 * pow(a,4) * pow(p,4) * _lisk->Li(2,x1) - 3 * pow(a,4) * pow(p,4) * _lisk->Li(2,x2) - 
                            3 * pow(a,4) * pow(p,4) * _lisk->Li(2,x3) + 3 * pow(a,4) * pow(p,4) * _lisk->Li(2,x4) - 3 * pow(a,4) * pow(p,4) * _lisk->Li(2,x5) - 
                            5 * pow(a,4) * pow(p,4) * _lisk->Li(2,x6) + 8 * pow(a,4) * pow(p,4) * _lisk->Li(2,x7) + 3 * pow(a,4) * pow(p,4) * _lisk->Li(2,x8) + 
                            5 * pow(a,4) * pow(p,4) * _lisk->Li(2,x9) - 8 * pow(a,4) * pow(p,4) * _lisk->Li(2,x10) - 80 * pow(a,3) * pow(p,3) * _lisk->Li(3,x6) - 
                            80 * pow(a,3) * pow(p,3) * _lisk->Li(3,x7) + 80 * pow(a,3) * pow(p,3) * _lisk->Li(3,x9) + 80 * pow(a,3) * pow(p,3) * _lisk->Li(3,x10) - 
                            360 * pow(a,2) * pow(p,2) * _lisk->Li(4,x6) + 360 * pow(a,2) * pow(p,2) * _lisk->Li(4,x7) + 
                            360 * pow(a,2) * pow(p,2) * _lisk->Li(4,x9) - 360 * pow(a,2) * pow(p,2) * _lisk->Li(4,x10) - 840 * a * p * _lisk->Li(5,x6) - 
                            840 * a * p * _lisk->Li(5,x7) + 840 * a * p * _lisk->Li(5,x9) + 840 * a * p * _lisk->Li(5,x10) - 840. * _lisk->Li(6,x6) + 
                            840. * _lisk->Li(6,x7) + 840. * _lisk->Li(6,x9) - 840. * _lisk->Li(6,x10)))/(128. * pow(a,5) * b * pow(p,6));
       
       _coefficients[8] = (-315. * (54 * a * pow(b,5) * pow(p,6) + 27. * pow(b,4) * pow(p,4) * _lisk->Li(2,x1) - 27. * pow(b,4) * pow(p,4) * _lisk->Li(2,x2) - 
                            27. * pow(b,4) * pow(p,4) * _lisk->Li(2,x3) + 27 * pow(b,4) * pow(p,4) * _lisk->Li(2,x4) - 27 * pow(b,4) * pow(p,4) * _lisk->Li(2,x5) + 
                            35 * pow(b,4) * pow(p,4) * _lisk->Li(2,x6) - 8 * pow(b,4) * pow(p,4) * _lisk->Li(2,x7) + 27 * pow(b,4) * pow(p,4) * _lisk->Li(2,x8) - 
                            35 * pow(b,4) * pow(p,4) * _lisk->Li(2,x9) + 8 * pow(b,4) * pow(p,4) * _lisk->Li(2,x10) + 80 * pow(b,3) * pow(p,3) * _lisk->Li(3,x6) - 
                            80 * pow(b,3) * pow(p,3) * _lisk->Li(3,x7) + 80 * pow(b,3) * pow(p,3) * _lisk->Li(3,x9) - 80 * pow(b,3) * pow(p,3) * _lisk->Li(3,x10) + 
                            360 * pow(b,2) * pow(p,2) * _lisk->Li(4,x6) - 360 * pow(b,2) * pow(p,2) * _lisk->Li(4,x7) - 
                            360 * pow(b,2) * pow(p,2) * _lisk->Li(4,x9) + 360 * pow(b,2) * pow(p,2) * _lisk->Li(4,x10) + 840 * b * p * _lisk->Li(5,x6) - 
                            840 * b * p * _lisk->Li(5,x7) + 840 * b * p * _lisk->Li(5,x9) - 840 * b * p * _lisk->Li(5,x10) + 840. * _lisk->Li(6,x6) - 
                            840. * _lisk->Li(6,x7) - 840. * _lisk->Li(6,x9) + 840. * _lisk->Li(6,x10)))/(128. * a * pow(b,5) * pow(p,6));
       
       _coefficients[9] = (45. * (pow(a,4) * b * pow(p,5) - 2 * pow(a,3) * b * pow(p,4) * log(1 - x2) - 2 * pow(a,3) * b * pow(p,4) * log(1 - x4) + 
                            2 * pow(a,3) * b * pow(p,4) * log(1 - x6) + 2 * pow(a,3) * b * pow(p,4) * log(1 - x9) - 6 * pow(a,2) * b * pow(p,3) * _lisk->Li(2,x2) - 
                            6 * pow(a,2) * b * pow(p,3) * _lisk->Li(2,x4) - 2 * pow(a,2) * b * pow(p,3) * _lisk->Li(2,x6) - 
                            4 * pow(a,2) * b * pow(p,3) * _lisk->Li(2,x7) - 2 * pow(a,2) * b * pow(p,3) * _lisk->Li(2,x9) - 
                            4 * pow(a,2) * b * pow(p,3) * _lisk->Li(2,x10) + 12 * a * b * pow(p,2) * _lisk->Li(3,x2) + 12 * a * b * pow(p,2) * _lisk->Li(3,x4) + 
                            4 * pow(a,2) * pow(p,2) * _lisk->Li(3,x6) - 4 * pow(a,2) * pow(p,2) * _lisk->Li(3,x7) + 12 * a * b * pow(p,2) * _lisk->Li(3,x7) - 
                            4 * pow(a,2) * pow(p,2) * _lisk->Li(3,x9) + 4 * pow(a,2) * pow(p,2) * _lisk->Li(3,x10) + 12 * a * b * pow(p,2) * _lisk->Li(3,x10) + 
                            12 * b * p * _lisk->Li(4,x1) - 12 * b * p * _lisk->Li(4,x2) + 12 * b * p * _lisk->Li(4,x3) - 12 * b * p * _lisk->Li(4,x4) + 
                            12 * b * p * _lisk->Li(4,x5) + 12 * a * p * _lisk->Li(4,x6) + 12 * a * p * _lisk->Li(4,x7) - 12 * b * p * _lisk->Li(4,x7) + 
                            12 * b * p * _lisk->Li(4,x8) - 12 * a * p * _lisk->Li(4,x9) - 12 * a * p * _lisk->Li(4,x10) - 12 * b * p * _lisk->Li(4,x10) + 
                            12. * _lisk->Li(5,x6) - 12. * _lisk->Li(5,x7) - 12. * _lisk->Li(5,x9) + 12. * _lisk->Li(5,x10)))/(16. * pow(a,3) * pow(b,2) * pow(p,5));
       
       _coefficients[10] = (45. * (a * pow(b,2) * pow(p,3) * _lisk->Li(2,x6) + a * pow(b,2) * pow(p,3) * _lisk->Li(2,x7) - 
                            a * pow(b,2) * pow(p,3) * _lisk->Li(2,x9) - a * pow(b,2) * pow(p,3) * _lisk->Li(2,x10) + 3 * a * b * pow(p,2) * _lisk->Li(3,x6) + 
                            pow(b,2) * pow(p,2) * _lisk->Li(3,x6) + 3 * a * b * pow(p,2) * _lisk->Li(3,x7) - pow(b,2) * pow(p,2) * _lisk->Li(3,x7) + 
                            3 * a * b * pow(p,2) * _lisk->Li(3,x9) - pow(b,2) * pow(p,2) * _lisk->Li(3,x9) + 3 * a * b * pow(p,2) * _lisk->Li(3,x10) + 
                            pow(b,2) * pow(p,2) * _lisk->Li(3,x10) + 3 * a * p * _lisk->Li(4,x6) + 3 * b * p * _lisk->Li(4,x6) + 
                            3 * a * p * _lisk->Li(4,x7) - 3 * b * p * _lisk->Li(4,x7) - 3 * a * p * _lisk->Li(4,x9) + 3 * b * p * _lisk->Li(4,x9) - 
                            3 * a * p * _lisk->Li(4,x10) - 3 * b * p * _lisk->Li(4,x10) + 3. * _lisk->Li(5,x6) - 3. * _lisk->Li(5,x7) - 
                            3. * _lisk->Li(5,x9) + 3. * _lisk->Li(5,x10)))/(4. * pow(a,2) * pow(b,3) * pow(p,5));
       
       _coefficients[11] = (105. * (2 * pow(a,5) * b * pow(p,6) - pow(a,4) * b * pow(p,5) * log(1 - x2) - pow(a,4) * b * pow(p,5) * log(1 - x4) + 
                            pow(a,4) * b * pow(p,5) * log(1 - x7) + pow(a,4) * b * pow(p,5) * log(-x6 * (1 - x4) * (1 - x9)) + 
                            pow(a,4) * b * pow(p,5) * log(1 - x10) - pow(a,4) * b * pow(p,5) * log((1 - x7) * (1 - x10)) + 
                            8 * pow(a,3) * b * pow(p,4) * _lisk->Li(2,x2) + 8 * pow(a,3) * b * pow(p,4) * _lisk->Li(2,x4) - 
                            8 * pow(a,3) * b * pow(p,4) * _lisk->Li(2,x7) - 8 * pow(a,3) * b * pow(p,4) * _lisk->Li(2,x10) - 12 * pow(a,2) * b * pow(p,3) * _lisk->Li(3,x1) - 
                            48 * pow(a,2) * b * pow(p,3) * _lisk->Li(3,x2) - 12 * pow(a,2) * b * pow(p,3) * _lisk->Li(3,x3) - 
                            48 * pow(a,2) * b * pow(p,3) * _lisk->Li(3,x4) + 12 * pow(a,2) * b * pow(p,3) * _lisk->Li(3,x5) - 8 * pow(a,3) * pow(p,3) * _lisk->Li(3,x6) - 
                            8 * pow(a,3) * pow(p,3) * _lisk->Li(3,x7) + 48 * pow(a,2) * b * pow(p,3) * _lisk->Li(3,x7) + 12 * pow(a,2) * b * pow(p,3) * _lisk->Li(3,x8) + 
                            8 * pow(a,3) * pow(p,3) * _lisk->Li(3,x9) + 8 * pow(a,3) * pow(p,3) * _lisk->Li(3,x10) + 48 * pow(a,2) * b * pow(p,3) * _lisk->Li(3,x10) + 
                            120 * a * b * pow(p,2) * _lisk->Li(4,x2) + 120 * a * b * pow(p,2) * _lisk->Li(4,x4) - 48 * pow(a,2) * pow(p,2) * _lisk->Li(4,x6) + 
                            48 * pow(a,2) * pow(p,2) * _lisk->Li(4,x7) - 120 * a * b * pow(p,2) * _lisk->Li(4,x7) + 48 * pow(a,2) * pow(p,2) * _lisk->Li(4,x9) - 
                            48 * pow(a,2) * pow(p,2) * _lisk->Li(4,x10) - 120 * a * b * pow(p,2) * _lisk->Li(4,x10) + 120 * b * p * _lisk->Li(5,x1) - 
                            120 * b * p * _lisk->Li(5,x2) + 120 * b * p * _lisk->Li(5,x3) - 120 * b * p * _lisk->Li(5,x4) - 120 * b * p * _lisk->Li(5,x5) - 
                            120 * a * p * _lisk->Li(5,x6) - 120 * a * p * _lisk->Li(5,x7) + 120 * b * p * _lisk->Li(5,x7) - 120 * b * p * _lisk->Li(5,x8) + 
                            120 * a * p * _lisk->Li(5,x9) + 120 * a * p * _lisk->Li(5,x10) + 120 * b * p * _lisk->Li(5,x10) - 120. * _lisk->Li(6,x6) + 
                            120. * _lisk->Li(6,x7) + 120. * _lisk->Li(6,x9) - 120.* _lisk->Li(6,x10)))/(32. * pow(a,4) * pow(b,2) * pow(p,6));
       
       _coefficients[12] = (-105. * (3 * pow(a,3) * pow(b,3) * pow(p,6) - 2 * pow(a,2) * pow(b,3) * pow(p,5) * log(1 + exp((a - b - d) * p)) + 5 * pow(a,2) * pow(b,3) * pow(p,5) * log(1 + exp((-a + b - d) * p)) - 2 * pow(a,2) * pow(b,3) * pow(p,5) * log(1 + exp((a + b - d) * p)) - 3 * pow(a,2) * pow(b,3) * pow(p,5) * log(1 + exp((a - b + d) * p)) + 2 * pow(a,2) * pow(b,3) * pow(p,5) * log(exp((-a - b + d) * p) * (1 + exp((a + b - d) * p)) * (1 + exp((-a + b + d) * p))) + 2 * pow(a,2) * pow(b,3) * pow(p,5) * log(1 + exp((a + b + d) * p)) - 2 * pow(a,2) * pow(b,3) * pow(p,5) * log((1 + exp((a - b + d) * p)) * (1 + exp((a + b + d) * p))) - 4 * a * pow(b,3) * pow(p,4) * _lisk->Li(2,-exp((a - b - d) * p)) - 10 * a * pow(b,3) * pow(p,4) * _lisk->Li(2,-exp((-a + b - d) * p)) - 4 * a * pow(b,3) * pow(p,4) * _lisk->Li(2,-exp((a + b - d) * p)) - 6 * a * pow(b,3) * pow(p,4) * _lisk->Li(2,-exp((a - b + d) * p)) + 4 * a * pow(b,3) * pow(p,4) * _lisk->Li(2,-exp((a + b + d) * p)) - 4 * pow(b,3) * pow(p,3) * _lisk->Li(3,-exp((-b - d) * p)) + 4 * pow(b,3) * pow(p,3) * _lisk->Li(3,-exp((a - b - d) * p)) + 6 * pow(b,3) * pow(p,3) * _lisk->Li(3,-exp((b - d) * p)) - 10 * pow(b,3) * pow(p,3) * _lisk->Li(3,-exp((-a + b - d) * p)) + 4 * pow(b,3) * pow(p,3) * _lisk->Li(3,-exp((a + b - d) * p)) - 
        6 * pow(b,3) * pow(p,3) * _lisk->Li(3,-exp((-b + d) * p)) + 24 * a * pow(b,2) * pow(p,3) * _lisk->Li(3,-exp((-a - b + d) * p)) + 24 * a * pow(b,2) * pow(p,3) * _lisk->Li(3,-exp((a - b + d) * p)) + 6 * pow(b,3) * pow(p,3) * _lisk->Li(3,-exp((a - b + d) * p)) + 4 * pow(b,3) * pow(p,3) * _lisk->Li(3,-exp((b + d) * p)) - 24 * a * pow(b,2) * pow(p,3) * _lisk->Li(3,-exp((-a + b + d) * p)) - 24 * a * pow(b,2) * pow(p,3) * _lisk->Li(3,-exp((a + b + d) * p)) - 4 * pow(b,3) * pow(p,3) * _lisk->Li(3,-exp((a + b + d) * p)) + 60 * a * b * pow(p,2) * _lisk->Li(4,-exp((-a - b + d) * p)) + 24 * pow(b,2) * pow(p,2) * _lisk->Li(4,-exp((-a - b + d) * p)) + 60 * a * b * pow(p,2) * _lisk->Li(4,-exp((a - b + d) * p)) - 24 * pow(b,2) * pow(p,2) * _lisk->Li(4,-exp((a - b + d) * p)) + 60 * a * b * pow(p,2) * _lisk->Li(4,-exp((-a + b + d) * p)) - 24 * pow(b,2) * pow(p,2) * _lisk->Li(4,-exp((-a + b + d) * p)) + 60 * a * b * pow(p,2) * _lisk->Li(4,-exp((a + b + d) * p)) + 24 * pow(b,2) * pow(p,2) * _lisk->Li(4,-exp((a + b + d) * p)) + 60 * a * p * _lisk->Li(5,-exp((-a - b + d) * p)) + 60 * b * p * _lisk->Li(5,-exp((-a - b + d) * p)) +                   60 * a * p * _lisk->Li(5,-exp((a - b + d) * p)) - 60 * b * p * _lisk->Li(5,-exp((a - b + d) * p)) - 60 * a * p * _lisk->Li(5,-exp((-a + b + d) * p)) + 60 * b * p * _lisk->Li(5,-exp((-a + b + d) * p)) - 60 * a * p * _lisk->Li(5,-exp((a + b + d) * p)) - 60 * b * p * _lisk->Li(5,-exp((a + b + d) * p)) + 60. * _lisk->Li(6,-exp((-a - b + d) * p)) - 60. * _lisk->Li(6,-exp((a - b + d) * p)) - 60. * _lisk->Li(6,-exp((-a + b + d) * p)) + 60. * _lisk->Li(6,-exp((a + b + d) * p))))/(16. * pow(a,2) * pow(b,4) * pow(p,6));
    
       _coefficients[13] = (51. * a)/(8. * b) - (117 * log(1 + exp((a - b - d) * p)))/(32. * b * p) + (315 * log(1 + exp((-a + b - d) * p)))/(16. * b * p) - (117 * log(1 + exp((a + b - d) * p)))/(32. * b * p) - 
    (513 * log(1 + exp((a - b + d) * p)))/(32. * b * p) + (117 * log(exp((-a - b + d) * p) * (1 + exp((a + b - d) * p)) * (1 + exp((-a + b + d) * p))))/(32. * b * p) + 
    (117 * log(1 + exp((a + b + d) * p)))/(32. * b * p) - (117 * log((1 + exp((a - b + d) * p)) * (1 + exp((a + b + d) * p))))/(32. * b * p) - (27. * _lisk->Li(2,-exp((a - b - d) * p)))/(a * b * pow(p,2)) - 
    (315. * _lisk->Li(2,-exp((-a + b - d) * p)))/(8. * a * b * pow(p,2)) - (27. * _lisk->Li(2,-exp((a + b - d) * p)))/(a * b * pow(p,2)) - (99. * _lisk->Li(2,-exp((a - b + d) * p)))/(8. * a * b * pow(p,2)) + 
    (27. * _lisk->Li(2,-exp((a + b + d) * p)))/(a * b * pow(p,2)) + (99. * _lisk->Li(3,-exp((-b - d) * p)))/(8. * pow(a,2) * b * pow(p,3)) + 
    (423. * _lisk->Li(3,-exp((a - b - d) * p)))/(4. * pow(a,2) * b * pow(p,3)) + (207. * _lisk->Li(3,-exp((b - d) * p)))/(4. * pow(a,2) * b * pow(p,3)) - 
    (315. * _lisk->Li(3,-exp((-a + b - d) * p)))/(8. * pow(a,2) * b * pow(p,3)) + (423. * _lisk->Li(3,-exp((a + b - d) * p)))/(4. * pow(a,2) * b * pow(p,3)) - 
    (207. * _lisk->Li(3,-exp((-b + d) * p)))/(4. * pow(a,2) * b * pow(p,3)) + (423. * _lisk->Li(3,-exp((-a - b + d) * p)))/(4. * a * pow(b,2) * pow(p,3)) + 
    (423. * _lisk->Li(3,-exp((a - b + d) * p)))/(4. * a * pow(b,2) * pow(p,3)) - (531. * _lisk->Li(3,-exp((a - b + d) * p)))/(8. * pow(a,2) * b * pow(p,3)) - 
    (99. * _lisk->Li(3,-exp((b + d) * p)))/(8. * pow(a,2) * b * pow(p,3)) - (423. * _lisk->Li(3,-exp((-a + b + d) * p)))/(4. * a * pow(b,2) * pow(p,3)) - 
    (423. * _lisk->Li(3,-exp((a + b + d) * p)))/(4. * a * pow(b,2) * pow(p,3)) - (423. * _lisk->Li(3,-exp((a + b + d) * p)))/(4. * pow(a,2) * b * pow(p,3)) - 
    (945. * _lisk->Li(4,-exp((a - b - d) * p)))/(4. * pow(a,3) * b * pow(p,4)) - (945. * _lisk->Li(4,-exp((a + b - d) * p)))/(4. * pow(a,3) * b * pow(p,4)) + 
    (945. * _lisk->Li(4,-exp((-a - b + d) * p)))/(4. * a * pow(b,3) * pow(p,4)) + (369. * _lisk->Li(4,-exp((-a - b + d) * p)))/(2. * pow(a,2) * pow(b,2) * pow(p,4)) + 
    (945. * _lisk->Li(4,-exp((a - b + d) * p)))/(4. * a * pow(b,3) * pow(p,4)) - (369. * _lisk->Li(4,-exp((a - b + d) * p)))/(2. * pow(a,2) * pow(b,2) * pow(p,4)) + 
    (945. * _lisk->Li(4,-exp((a - b + d) * p)))/(4. * pow(a,3) * b * pow(p,4)) + (945. * _lisk->Li(4,-exp((-a + b + d) * p)))/(4. * a * pow(b,3) * pow(p,4)) - 
    (369. * _lisk->Li(4,-exp((-a + b + d) * p)))/(2. * pow(a,2) * pow(b,2) * pow(p,4)) + (945. * _lisk->Li(4,-exp((a + b + d) * p)))/(4. * a * pow(b,3) * pow(p,4)) + 
    (369. * _lisk->Li(4,-exp((a + b + d) * p)))/(2. * pow(a,2) * pow(b,2) * pow(p,4)) + (945. * _lisk->Li(4,-exp((a + b + d) * p)))/(4. * pow(a,3) * b * pow(p,4)) - 
    (945. * _lisk->Li(5,-exp((-b - d) * p)))/(4. * pow(a,4) * b * pow(p,5)) + (945. * _lisk->Li(5,-exp((a - b - d) * p)))/(4. * pow(a,4) * b * pow(p,5)) - 
    (945. * _lisk->Li(5,-exp((b - d) * p)))/(4. * pow(a,4) * b * pow(p,5)) + (945. * _lisk->Li(5,-exp((a + b - d) * p)))/(4. * pow(a,4) * b * pow(p,5)) + 
    (945. * _lisk->Li(5,-exp((-b + d) * p)))/(4. * pow(a,4) * b * pow(p,5)) + (945. * _lisk->Li(5,-exp((-a - b + d) * p)))/(4. * a * pow(b,4) * pow(p,5)) + 
    (945. * _lisk->Li(5,-exp((-a - b + d) * p)))/(4. * pow(a,2) * pow(b,3) * pow(p,5)) + (945. * _lisk->Li(5,-exp((-a - b + d) * p)))/(4. * pow(a,3) * pow(b,2) * pow(p,5)) + 
    (945. * _lisk->Li(5,-exp((a - b + d) * p)))/(4. * a * pow(b,4) * pow(p,5)) - (945. * _lisk->Li(5,-exp((a - b + d) * p)))/(4. * pow(a,2) * pow(b,3) * pow(p,5)) + 
    (945. * _lisk->Li(5,-exp((a - b + d) * p)))/(4. * pow(a,3) * pow(b,2) * pow(p,5)) - (945. * _lisk->Li(5,-exp((a - b + d) * p)))/(4. * pow(a,4) * b * pow(p,5)) + 
    (945. * _lisk->Li(5,-exp((b + d) * p)))/(4. * pow(a,4) * b * pow(p,5)) - (945. * _lisk->Li(5,-exp((-a + b + d) * p)))/(4. * a * pow(b,4) * pow(p,5)) + 
    (945. * _lisk->Li(5,-exp((-a + b + d) * p)))/(4. * pow(a,2) * pow(b,3) * pow(p,5)) - (945. * _lisk->Li(5,-exp((-a + b + d) * p)))/(4. * pow(a,3) * pow(b,2) * pow(p,5)) - 
    (945. * _lisk->Li(5,-exp((a + b + d) * p)))/(4. * a * pow(b,4) * pow(p,5)) - (945. * _lisk->Li(5,-exp((a + b + d) * p)))/(4. * pow(a,2) * pow(b,3) * pow(p,5)) - 
    (945. * _lisk->Li(5,-exp((a + b + d) * p)))/(4. * pow(a,3) * pow(b,2) * pow(p,5)) - (945. * _lisk->Li(5,-exp((a + b + d) * p)))/(4. * pow(a,4) * b * pow(p,5)) + 
    (945. * _lisk->Li(6,-exp((-a - b + d) * p)))/(4. * pow(a,2) * pow(b,4) * pow(p,6)) + (945. * _lisk->Li(6,-exp((-a - b + d) * p)))/(4. * pow(a,4) * pow(b,2) * pow(p,6)) - 
    (945. * _lisk->Li(6,-exp((a - b + d) * p)))/(4. * pow(a,2) * pow(b,4) * pow(p,6)) - (945. * _lisk->Li(6,-exp((a - b + d) * p)))/(4. * pow(a,4) * pow(b,2) * pow(p,6)) - 
    (945. * _lisk->Li(6,-exp((-a + b + d) * p)))/(4. * pow(a,2) * pow(b,4) * pow(p,6)) - (945. * _lisk->Li(6,-exp((-a + b + d) * p)))/(4. * pow(a,4) * pow(b,2) * pow(p,6)) + 
    (945. * _lisk->Li(6,-exp((a + b + d) * p)))/(4. * pow(a,2) * pow(b,4) * pow(p,6)) + (945. * _lisk->Li(6,-exp((a + b + d) * p)))/(4. * pow(a,4) * pow(b,2) * pow(p,6));
       
       _coefficients[14] = (-225. * (3 * pow(a,4) * pow(b,2) * pow(p,6) + 16 * pow(a,3) * pow(b,3) * pow(p,6) + 12 * pow(a,3) * pow(b,2) * pow(p,5) * log(1 - x12) + 
                            12 * pow(a,3) * pow(b,2) * pow(p,5) * log(1 - x2) - 12 * pow(a,3) * pow(b,2) * pow(p,5) * log(1 - x4) - 
                            12 * pow(a,3) * pow(b,2) * pow(p,5) * log(1 - x7) - 12 * pow(a,3) * pow(b,2) * pow(p,5) * log((1 - x10)/(1 - x7)) - 
                            12 * pow(a,3) * pow(b,2) * pow(p,5) * log((-x13 + exp((2 * b + d) * p))/(-x13 + exp(d * p))) + 
                            8 * pow(a,2) * pow(b,2) * pow(p,4) * _lisk->Li(2,x1) - 36 * pow(a,2) * pow(b,2) * pow(p,4) * _lisk->Li(2,x12) + 
                            28 * pow(a,2) * pow(b,2) * pow(p,4) * _lisk->Li(2,x2) - 8 * pow(a,2) * pow(b,2) * pow(p,4) * _lisk->Li(2,x3) - 
                            28 * pow(a,2) * pow(b,2) * pow(p,4) * _lisk->Li(2,x4) - 8 * pow(a,2) * pow(b,2) * pow(p,4) * _lisk->Li(2,x5) - 
                            12 * pow(a,2) * pow(b,2) * pow(p,4) * _lisk->Li(2,x6) - 16 * pow(a,2) * pow(b,2) * pow(p,4) * _lisk->Li(2,x7) + 
                            8 * pow(a,2) * pow(b,2) * pow(p,4) * _lisk->Li(2,x8) + 12 * pow(a,2) * pow(b,2) * pow(p,4) * _lisk->Li(2,x9) - 
                            20 * pow(a,2) * pow(b,2) * pow(p,4) * _lisk->Li(2,x10) - 72 * a * pow(b,2) * pow(p,3) * _lisk->Li(3,x12) - 
                            72 * a * pow(b,2) * pow(p,3) * _lisk->Li(3,x2) + 72 * a * pow(b,2) * pow(p,3) * _lisk->Li(3,x4) + 
                            48 * pow(a,2) * b * pow(p,3) * _lisk->Li(3,x6) - 24 * a * pow(b,2) * pow(p,3) * _lisk->Li(3,x6) - 
                            48 * pow(a,2) * b * pow(p,3) * _lisk->Li(3,x7) + 48 * a * pow(b,2) * pow(p,3) * _lisk->Li(3,x7) + 
                            48 * pow(a,2) * b * pow(p,3) * _lisk->Li(3,x9) + 24 * a * pow(b,2) * pow(p,3) * _lisk->Li(3,x9) - 
                            48 * pow(a,2) * b * pow(p,3) * _lisk->Li(3,x10) + 24 * a * pow(b,2) * pow(p,3) * _lisk->Li(3,x10) - 
                            72 * pow(b,2) * pow(p,2) * _lisk->Li(4,x12) + 72 * pow(b,2) * pow(p,2) * _lisk->Li(4,x2) + 72 * pow(b,2) * pow(p,2) * _lisk->Li(4,x3) - 
                            72 * pow(b,2) * pow(p,2) * _lisk->Li(4,x4) + 72 * pow(b,2) * pow(p,2) * _lisk->Li(4,x5) + 48 * pow(a,2) * pow(p,2) * _lisk->Li(4,x6) + 
                            144 * a * b * pow(p,2) * _lisk->Li(4,x6) - 24 * pow(b,2) * pow(p,2) * _lisk->Li(4,x6) - 48 * pow(a,2) * pow(p,2) * _lisk->Li(4,x7) + 
                            144 * a * b * pow(p,2) * _lisk->Li(4,x7) - 48 * pow(b,2) * pow(p,2) * _lisk->Li(4,x7) - 48 * pow(a,2) * pow(p,2) * _lisk->Li(4,x9) + 
                            144 * a * b * pow(p,2) * _lisk->Li(4,x9) + 24 * pow(b,2) * pow(p,2) * _lisk->Li(4,x9) + 48 * pow(a,2) * pow(p,2) * _lisk->Li(4,x10) + 
                            144 * a * b * pow(p,2) * _lisk->Li(4,x10) - 24 * pow(b,2) * pow(p,2) * _lisk->Li(4,x10) + 144 * a * p * _lisk->Li(5,x6) + 
                            144 * b * p * _lisk->Li(5,x6) + 144 * a * p * _lisk->Li(5,x7) - 144 * b * p * _lisk->Li(5,x7) - 144 * a * p * _lisk->Li(5,x9) + 
                            144 * b * p * _lisk->Li(5,x9) - 144 * a * p * _lisk->Li(5,x10) - 144 * b * p * _lisk->Li(5,x10) + 144. * _lisk->Li(6,x6) - 
                            144. * _lisk->Li(6,x7) - 144. * _lisk->Li(6,x9) + 144. * _lisk->Li(6,x10)))/(128. * pow(a,3) * pow(b,3) * pow(p,6));
       
    
       
        
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

