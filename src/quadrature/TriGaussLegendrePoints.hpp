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

#ifndef __femus_quadrature_TriGaussLegendrePoints_hpp__
#define __femus_quadrature_TriGaussLegendrePoints_hpp__

namespace femus {

 class tri_gaussLegendre {
 public:
 static const unsigned GaussPoints[38];
 static const double *Gauss[38];
 static const double Gauss0[3][1];
 static const double Gauss1[3][1];
 static const double Gauss2[3][4];
 static const double Gauss3[3][4];
 static const double Gauss4[3][7];
 static const double Gauss5[3][7];
 static const double Gauss6[3][13];
 static const double Gauss7[3][13];
 static const double Gauss8[3][19];
 static const double Gauss9[3][19];
 static const double Gauss10[3][28];
 static const double Gauss11[3][28];
 static const double Gauss12[3][37];
 static const double Gauss13[3][37];
 static const double Gauss14[3][64];
 static const double Gauss15[3][81];
 static const double Gauss16[3][81];
 static const double Gauss17[3][100];
 static const double Gauss18[3][100];
 static const double Gauss19[3][121];
 static const double Gauss20[3][121];
 static const double Gauss21[3][144];
 static const double Gauss22[3][144];
 static const double Gauss23[3][169];
 static const double Gauss24[3][169];
 static const double Gauss25[3][196];
 static const double Gauss26[3][196];
 static const double Gauss27[3][225];
 static const double Gauss28[3][225];
 static const double Gauss29[3][256];
 static const double Gauss30[3][256];
 static const double Gauss31[3][289];
 static const double Gauss32[3][289];
 static const double Gauss33[3][324];
 static const double Gauss34[3][324];
 static const double Gauss35[3][361];
 static const double Gauss36[3][361];
 static const double Gauss37[3][400];
 };

 };

#endif