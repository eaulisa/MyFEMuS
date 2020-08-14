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

#ifndef __femus_quadrature_TriGaussPoints_hpp__
#define __femus_quadrature_TriGaussPoints_hpp__

#include <vector>
#include <string>

namespace femus {

  class tri_gauss {
  public:
    static const unsigned GaussPoints[7];
    static const double *Gauss[7];  
    static const double Gauss0[3][1];
    static const double Gauss1[3][4];
    static const double Gauss2[3][7];
    static const double Gauss3[3][13];
    static const double Gauss4[3][19];
    static const double Gauss5[3][28];
    static const double Gauss6[3][37];
  };
     
} //end namespace femus     


#endif
