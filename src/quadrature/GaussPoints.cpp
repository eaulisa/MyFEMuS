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

#include "GaussPoints.hpp"
#include <iostream>
#include <stdlib.h>
#include <string.h>

namespace femus {

  unsigned GetGaussOrder(const char* order_gauss) {
    for(unsigned i = 0; i < numberName.size(); i++) {
      if(!strcmp(order_gauss, numberName[i].c_str())) {
        return i;
      }
    }
    std::cout << order_gauss << " is not a valid option for the Gauss points\n";
    abort();
  }


  Gauss::Gauss(const GeomElType &GeomElemType, const unsigned &gauss_order, const char *gauss_type): _gauss_order(gauss_order) {
    _order = "unknown";

    if(_gauss_order >= 37) {
      std::cout << "WARNING! **** Gauss order rule greater than 37 is not available. Set it to 37 **** WARNING!" << std::endl;
      _gauss_order = 37;
    }

    if(!strcmp(gauss_type, "legendre")) {
      if(GeomElemType == HEX) {
        _GaussWeight = hex_gaussLegendre::Gauss[_gauss_order];
        _GaussPoints = hex_gaussLegendre::GaussPoints[_gauss_order];
      }
      else if(GeomElemType == WEDGE) {
        _GaussWeight = wedge_gaussLegendre::Gauss[_gauss_order];
        _GaussPoints = wedge_gaussLegendre::GaussPoints[_gauss_order];
      }
      else if(GeomElemType == TET) {
        _GaussWeight = tet_gaussLegendre::Gauss[_gauss_order];
        _GaussPoints = tet_gaussLegendre::GaussPoints[_gauss_order];
      }
      else if(GeomElemType == QUAD) {
        _GaussWeight = quad_gaussLegendre::Gauss[_gauss_order];
        _GaussPoints = quad_gaussLegendre::GaussPoints[_gauss_order];
      }
      else if(GeomElemType == TRI) {
        _GaussWeight = tri_gaussLegendre::Gauss[_gauss_order];
        _GaussPoints = tri_gaussLegendre::GaussPoints[_gauss_order];
      }
      else if(GeomElemType == LINE) {
        _GaussWeight = line_gaussLegendre::Gauss[_gauss_order];
        _GaussPoints = line_gaussLegendre::GaussPoints[_gauss_order];
      }
      else if(GeomElemType == POINT) {
        _GaussWeight = point_gauss::Gauss[_gauss_order];
        _GaussPoints = point_gauss::GaussPoints[_gauss_order];
      }
      else {
        std::cout << GeomElemType << " is not a valid option" << std::endl;
        abort();
      }
    }
    else if(!strcmp(gauss_type, "lobatto")) {
      if(GeomElemType == HEX) {
        _GaussWeight = hex_gaussLobatto::Gauss[_gauss_order];
        _GaussPoints = hex_gaussLobatto::GaussPoints[_gauss_order];
      }
      else if(GeomElemType == WEDGE) {
        _GaussWeight = wedge_gaussLobatto::Gauss[_gauss_order];
        _GaussPoints = wedge_gaussLobatto::GaussPoints[_gauss_order];
      }
      else if(GeomElemType == TET) {
        _GaussWeight = tet_gaussLobatto::Gauss[_gauss_order];
        _GaussPoints = tet_gaussLobatto::GaussPoints[_gauss_order];
      }
      else if(GeomElemType == QUAD) {
        _GaussWeight = quad_gaussLobatto::Gauss[_gauss_order];
        _GaussPoints = quad_gaussLobatto::GaussPoints[_gauss_order];
      }
      else if(GeomElemType == TRI) {
        _GaussWeight = tri_gaussLobatto::Gauss[_gauss_order];
        _GaussPoints = tri_gaussLobatto::GaussPoints[_gauss_order];
      }
      else if(GeomElemType == LINE) {
        _GaussWeight = line_gaussLobatto::Gauss[_gauss_order];
        _GaussPoints = line_gaussLobatto::GaussPoints[_gauss_order];
      }
      else if(GeomElemType == POINT) {
        _GaussWeight = point_gauss::Gauss[_gauss_order];
        _GaussPoints = point_gauss::GaussPoints[_gauss_order];
      }
      else {
        std::cout << GeomElemType << " is not a valid option" << std::endl;
        abort();
      }
    }
    else {
      std::cout << gauss_type << " is not a valid option" << std::endl;
      abort();
    }
  }




  Gauss::Gauss(const char *geom_elem, const char *order_gauss, const char *gauss_type) : _order(order_gauss) {
//     if(!strcmp(order_gauss, "zero")) {
//       _gauss_order = 0;
//     }
//     else if(!strcmp(order_gauss, "first")) {
//       _gauss_order = 1;
//     }
//     else if(!strcmp(order_gauss, "second")) {
//       _gauss_order = 2;
//     }
//     else if(!strcmp(order_gauss, "third")) {
//       _gauss_order = 3;
//     }
//     else if(!strcmp(order_gauss, "fourth")) {
//       _gauss_order = 4;
//     }
//     else if(!strcmp(order_gauss, "fifth")) {
//       _gauss_order = 5;
//     }
//     else if(!strcmp(order_gauss, "sixth")) {
//       _gauss_order = 6;
//     }
//     else if(!strcmp(order_gauss, "seventh")) {
//       _gauss_order = 7;
//     }
//     else if(!strcmp(order_gauss, "eighth")) {
//       _gauss_order = 8;
//     }
//     else if(!strcmp(order_gauss, "ninth")) {
//       _gauss_order = 9;
//     }
//     else if(!strcmp(order_gauss, "tenth")) {
//       _gauss_order = 10;
//     }
//     else if(!strcmp(order_gauss, "eleventh")) {
//       _gauss_order = 11;
//     }
//     else if(!strcmp(order_gauss, "twelfth")) {
//       _gauss_order = 12;
//     }
//     else if(!strcmp(order_gauss, "thirteenth")) {
//       _gauss_order = 13;
//     }
//     else if(!strcmp(order_gauss, "fourteenth")) {
//       _gauss_order = 14;
//     }
//     else if(!strcmp(order_gauss, "fifteenth")) {
//       _gauss_order = 15;
//     }
//     else if(!strcmp(order_gauss, "sixteenth")) {
//       _gauss_order = 16;
//     }
//     else if(!strcmp(order_gauss, "seventeenth")) {
//       _gauss_order = 17;
//     }
//     else if(!strcmp(order_gauss, "eighteenth")) {
//       _gauss_order = 18;
//     }
//     else if(!strcmp(order_gauss, "nineteenth")) {
//       _gauss_order = 19;
//     }
//     else if(!strcmp(order_gauss, "twentieth")) {
//       _gauss_order = 20;
//     }
//     else if(!strcmp(order_gauss, "twenty first")) {
//       _gauss_order = 21;
//     }
//     else if(!strcmp(order_gauss, "twenty second")) {
//       _gauss_order = 22;
//     }
//     else if(!strcmp(order_gauss, "twenty third")) {
//       _gauss_order = 23;
//     }
//     else if(!strcmp(order_gauss, "twenty fourth")) {
//       _gauss_order = 24;
//     }
//     else if(!strcmp(order_gauss, "twenty fifth")) {
//       _gauss_order = 25;
//     }
//     else if(!strcmp(order_gauss, "twenty sixth")) {
//       _gauss_order = 26;
//     }
//     else if(!strcmp(order_gauss, "twenty seventh")) {
//       _gauss_order = 27;
//     }
//     else if(!strcmp(order_gauss, "twenty eighth")) {
//       _gauss_order = 28;
//     }
//     else if(!strcmp(order_gauss, "twenty ninth")) {
//       _gauss_order = 29;
//     }
//     else if(!strcmp(order_gauss, "thirtieth")) {
//       _gauss_order = 30;
//     }
//     else if(!strcmp(order_gauss, "thirty first")) {
//       _gauss_order = 31;
//     }
//     else if(!strcmp(order_gauss, "thirty second")) {
//       _gauss_order = 32;
//     }
//     else if(!strcmp(order_gauss, "thirty third")) {
//       _gauss_order = 33;
//     }
//     else if(!strcmp(order_gauss, "thirty fourth")) {
//       _gauss_order = 34;
//     }
//     else if(!strcmp(order_gauss, "thirty fifth")) {
//       _gauss_order = 35;
//     }
//     else if(!strcmp(order_gauss, "thirty sixth")) {
//       _gauss_order = 36;
//     }
//     else if(!strcmp(order_gauss, "thirty seventh")) {
//       _gauss_order = 37;
//     }
//     else {
//       std::cout << order_gauss << " is not a valid option for the Gauss points of " << geom_elem << std::endl;
//       abort();
//     }
    
    _gauss_order = GetGaussOrder(order_gauss);

    if(!strcmp(gauss_type, "legendre")) {
      if(!strcmp(geom_elem, "hex")) {
        _GaussWeight = hex_gaussLegendre::Gauss[_gauss_order];
        _GaussPoints = hex_gaussLegendre::GaussPoints[_gauss_order];
      }
      else if(!strcmp(geom_elem, "wedge")) {
        _GaussWeight = wedge_gaussLegendre::Gauss[_gauss_order];
        _GaussPoints = wedge_gaussLegendre::GaussPoints[_gauss_order];
      }
      else if(!strcmp(geom_elem, "tet")) {
        _GaussWeight = tet_gaussLegendre::Gauss[_gauss_order];
        _GaussPoints = tet_gaussLegendre::GaussPoints[_gauss_order];
      }
      else if(!strcmp(geom_elem, "quad")) {
        _GaussWeight = quad_gaussLegendre::Gauss[_gauss_order];
        _GaussPoints = quad_gaussLegendre::GaussPoints[_gauss_order];
      }
      else if(!strcmp(geom_elem, "tri")) {
        _GaussWeight = tri_gaussLegendre::Gauss[_gauss_order];
        _GaussPoints = tri_gaussLegendre::GaussPoints[_gauss_order];
      }
      else if(!strcmp(geom_elem, "line")) {
        _GaussWeight = line_gaussLegendre::Gauss[_gauss_order];
        _GaussPoints = line_gaussLegendre::GaussPoints[_gauss_order];
      }
      else if(!strcmp(geom_elem, "point")) {
        _GaussWeight = point_gauss::Gauss[_gauss_order];
        _GaussPoints = point_gauss::GaussPoints[_gauss_order];
      }
      else {
        std::cout << geom_elem << " is not a valid option" << std::endl;
        abort();
      }
    }
    else if(!strcmp(gauss_type, "lobatto")) {
      if(!strcmp(geom_elem, "hex")) {
        _GaussWeight = hex_gaussLobatto::Gauss[_gauss_order];
        _GaussPoints = hex_gaussLobatto::GaussPoints[_gauss_order];
      }
      else if(!strcmp(geom_elem, "wedge")) {
        _GaussWeight = wedge_gaussLobatto::Gauss[_gauss_order];
        _GaussPoints = wedge_gaussLobatto::GaussPoints[_gauss_order];
      }
      else if(!strcmp(geom_elem, "tet")) {
        _GaussWeight = tet_gaussLobatto::Gauss[_gauss_order];
        _GaussPoints = tet_gaussLobatto::GaussPoints[_gauss_order];
      }
      else if(!strcmp(geom_elem, "quad")) {
        _GaussWeight = quad_gaussLobatto::Gauss[_gauss_order];
        _GaussPoints = quad_gaussLobatto::GaussPoints[_gauss_order];
      }
      else if(!strcmp(geom_elem, "tri")) {
        _GaussWeight = tri_gaussLobatto::Gauss[_gauss_order];
        _GaussPoints = tri_gaussLobatto::GaussPoints[_gauss_order];
      }
      else if(!strcmp(geom_elem, "line")) {
        _GaussWeight = line_gaussLobatto::Gauss[_gauss_order];
        _GaussPoints = line_gaussLobatto::GaussPoints[_gauss_order];
      }
      else if(!strcmp(geom_elem, "point")) {
        _GaussWeight = point_gauss::Gauss[_gauss_order];
        _GaussPoints = point_gauss::GaussPoints[_gauss_order];
      }
      else {
        std::cout << geom_elem << " is not a valid option" << std::endl;
        abort();
      }
    }
    else {
      std::cout << gauss_type << " is not a valid option" << std::endl;
      abort();
    }
  }
} //end namespace femus
