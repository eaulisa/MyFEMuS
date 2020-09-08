/*=========================================================================

Program: FEMUS
Module: Gauss
Authors: Eugenio Aulisa, Giorgio Bornia, Erdi Kara

Copyright (c) FEMTTU
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_quadrature_LineGaussLobattoPoints_hpp__
#define __femus_quadrature_LineGaussLobattoPoints_hpp__

namespace femus {

 class line_gaussLobatto {
 public:
 static const unsigned GaussPoints[38];
 static const double *Gauss[38];
 static const double Gauss0[2][2];
 static const double Gauss1[2][2];
 static const double Gauss2[2][3];
 static const double Gauss3[2][3];
 static const double Gauss4[2][4];
 static const double Gauss5[2][4];
 static const double Gauss6[2][5];
 static const double Gauss7[2][5];
 static const double Gauss8[2][6];
 static const double Gauss9[2][6];
 static const double Gauss10[2][7];
 static const double Gauss11[2][7];
 static const double Gauss12[2][8];
 static const double Gauss13[2][8];
 static const double Gauss14[2][9];
 static const double Gauss15[2][9];
 static const double Gauss16[2][10];
 static const double Gauss17[2][10];
 static const double Gauss18[2][11];
 static const double Gauss19[2][11];
 static const double Gauss20[2][12];
 static const double Gauss21[2][12];
 static const double Gauss22[2][13];
 static const double Gauss23[2][13];
 static const double Gauss24[2][14];
 static const double Gauss25[2][14];
 static const double Gauss26[2][15];
 static const double Gauss27[2][15];
 static const double Gauss28[2][16];
 static const double Gauss29[2][16];
 static const double Gauss30[2][17];
 static const double Gauss31[2][17];
 static const double Gauss32[2][18];
 static const double Gauss33[2][18];
 static const double Gauss34[2][19];
 static const double Gauss35[2][19];
 static const double Gauss36[2][20];
 static const double Gauss37[2][20];
 };

 };

#endif