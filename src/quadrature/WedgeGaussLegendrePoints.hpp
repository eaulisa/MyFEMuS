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

#ifndef __femus_quadrature_WedgeGaussLegendrePoints_hpp__
#define __femus_quadrature_WedgeGaussLegendrePoints_hpp__

namespace femus {

 class wedge_gaussLegendre {
 public:
 static const unsigned GaussPoints[38];
 static const double *Gauss[38];
 static const double Gauss0[4][1];
 static const double Gauss1[4][1];
 static const double Gauss2[4][8];
 static const double Gauss3[4][8];
 static const double Gauss4[4][21];
 static const double Gauss5[4][21];
 static const double Gauss6[4][52];
 static const double Gauss7[4][52];
 static const double Gauss8[4][95];
 static const double Gauss9[4][95];
 static const double Gauss10[4][168];
 static const double Gauss11[4][168];
 static const double Gauss12[4][259];
 static const double Gauss13[4][259];
 static const double Gauss14[4][512];
 static const double Gauss15[4][648];
 static const double Gauss16[4][729];
 static const double Gauss17[4][900];
 static const double Gauss18[4][1000];
 static const double Gauss19[4][1210];
 static const double Gauss20[4][1331];
 static const double Gauss21[4][1584];
 static const double Gauss22[4][1728];
 static const double Gauss23[4][2028];
 static const double Gauss24[4][2197];
 static const double Gauss25[4][2548];
 static const double Gauss26[4][2744];
 static const double Gauss27[4][3150];
 static const double Gauss28[4][3375];
 static const double Gauss29[4][3840];
 static const double Gauss30[4][4096];
 static const double Gauss31[4][4624];
 static const double Gauss32[4][4913];
 static const double Gauss33[4][5508];
 static const double Gauss34[4][5832];
 static const double Gauss35[4][6498];
 static const double Gauss36[4][6859];
 static const double Gauss37[4][7600];
 };

 };

#endif