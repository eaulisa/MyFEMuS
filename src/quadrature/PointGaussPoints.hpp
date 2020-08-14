/*=========================================================================

 Program: FEMUS
 Module: FemusInit
 Authors: Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_quadrature_PointGaussPoints_hpp__
#define __femus_quadrature_PointGaussPoints_hpp__

#include <vector>
#include <string>

namespace femus {
  
  class point_gauss {
  public:
    static const unsigned GaussPoints[5];
    static const double *Gauss[5];  
    static const double Gauss0[2][1];
    static const double Gauss1[2][1];
    static const double Gauss2[2][1];
    static const double Gauss3[2][1];
    static const double Gauss4[2][1];
  };  
     
} //end namespace femus     


#endif
