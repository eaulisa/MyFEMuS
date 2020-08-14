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

#ifndef __femus_quadrature_LineGaussPoints_hpp__
#define __femus_quadrature_LineGaussPoints_hpp__

#include <vector>
#include <string>

namespace femus {
  
  class line_gauss {
  public:
    static const unsigned GaussPoints[6];
    static const double *Gauss[6];  
    static const double Gauss0[2][1];
    static const double Gauss1[2][2];
    static const double Gauss2[2][3];
    static const double Gauss3[2][4];
    static const double Gauss4[2][5];
    static const double Gauss5[2][6];
    static const double Gauss6[2][7];
    static const double Gauss7[2][8];
    static const double Gauss8[2][9];
    static const double Gauss9[2][10];
    static const double Gauss10[2][11];
    static const double Gauss11[2][12];
    static const double Gauss12[2][13];
    static const double Gauss13[2][14];
    static const double Gauss14[2][15];
    static const double Gauss15[2][16];
    static const double Gauss16[2][17];
    static const double Gauss17[2][18];
    static const double Gauss18[2][19];
    static const double Gauss19[2][20];
        
  };  

     
} //end namespace femus     


#endif
