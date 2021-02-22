/*=========================================================================

  Program: FEMUS
  Module: Gauss
  Authors: Eugenio Aulisa, Giorgio Bornia

  Copyright (c) FEMTTU
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/

#include "PointGaussPoints.hpp"
#include <iostream>
#include <stdlib.h>
#include <string.h>

namespace femus {
// ************** POINT ***************
  const unsigned point_gauss::GaussPoints[38] = {1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19};
  const double * point_gauss::Gauss[38] = {Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], 
                                          Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], 
                                          Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], 
                                          Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0], Gauss0[0]};
  
  
  //first row-weights, second row: x-coordinates
  const double point_gauss::Gauss0[2][1] = {{1}, {0}};


} //end namespace femus
