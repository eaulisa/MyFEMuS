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

#ifndef __femus_quadrature_GaussPoints_hpp__
#define __femus_quadrature_GaussPoints_hpp__

#include <vector>
#include <string>

#include "PointGaussPoints.hpp"

#include "HexGaussLegendrePoints.hpp"
#include "WedgeGaussLegendrePoints.hpp"
#include "TetGaussLegendrePoints.hpp"
#include "QuadGaussLegendrePoints.hpp"
#include "TriGaussLegendrePoints.hpp"
#include "LineGaussLegendrePoints.hpp"

#include "HexGaussLobattoPoints.hpp"
#include "WedgeGaussLobattoPoints.hpp"
#include "TetGaussLobattoPoints.hpp"
#include "QuadGaussLobattoPoints.hpp"
#include "TriGaussLobattoPoints.hpp"
#include "LineGaussLobattoPoints.hpp"

#include "GeomElTypeEnum.hpp"

namespace femus {
  
  class Gauss {
     
  public:

    Gauss(const char *geom_elem, const char *order_gauss, const char *gauss_type = "legendre" );
    Gauss(const GeomElType &GeomElemType, const unsigned &gauss_order, const char *gauss_type = "legendre" );
    
  inline const double *  GetGaussWeightsPointer() const {
    return _GaussWeight;
  };
  
  inline const double *  GetGaussCoordinatePointer (const unsigned &k) const {
    return _GaussWeight + (k+1) * _GaussPoints;
  };
  
  
  inline const double  GetGaussWeight(const unsigned ig) const {
    return _GaussWeight[ig];
  };
  
  inline const unsigned GetGaussPointsNumber() const {
      return _GaussPoints;
  };     

  inline const std::string  GetGaussOrderString() const {
    return _order;
  };
  
  inline unsigned GetGaussOrderIdx() const {
    return _gauss_order;
  };
  
  protected:
    
    unsigned _gauss_order;
    std::string _order;
    unsigned _GaussPoints;  
    const double *_GaussWeight;
   
  };
     

     
} //end namespace femus     


#endif
