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

  Gauss::Gauss(const char *geom_elem, const char *order_gauss, const char *gauss_type) : _order(order_gauss) {
    if(!strcmp(order_gauss, "zero")) {
      gauss_order = 0;
    }
    else if(!strcmp(order_gauss, "first")) {
      gauss_order = 1;
    }
    else if(!strcmp(order_gauss, "second")) {
      gauss_order = 2;
    }
    else if(!strcmp(order_gauss, "third")) {
      gauss_order = 3;
    }
    else if(!strcmp(order_gauss, "fourth")) {
      gauss_order = 4;
    }
    else if(!strcmp(order_gauss, "fifth")) {
      gauss_order = 5;
    }
    else if(!strcmp(order_gauss, "sixth")) {
      gauss_order = 6;
    }
    else if(!strcmp(order_gauss, "seventh")) {
      gauss_order = 7;
    }
    else if(!strcmp(order_gauss, "eighth")) {
      gauss_order = 8;
    }
    else if(!strcmp(order_gauss, "ninth")) {
      gauss_order = 9;
    }
    else if(!strcmp(order_gauss, "tenth")) {
      gauss_order = 10;
    }
    else if(!strcmp(order_gauss, "eleventh")) {
      gauss_order = 11;
    }
    else if(!strcmp(order_gauss, "twelfth")) {
      gauss_order = 12;
    }
    else if(!strcmp(order_gauss, "thirteenth")) {
      gauss_order = 13;
    }
    else if(!strcmp(order_gauss, "fourteenth")) {
      gauss_order = 14;
    }
    else if(!strcmp(order_gauss, "fifteenth")) {
      gauss_order = 15;
    }
    else if(!strcmp(order_gauss, "sixteenth")) {
      gauss_order = 16;
    }
    else if(!strcmp(order_gauss, "seventeenth")) {
      gauss_order = 17;
    }
    else if(!strcmp(order_gauss, "eighteenth")) {
      gauss_order = 18;
    }
    else if(!strcmp(order_gauss, "nineteenth")) {
      gauss_order = 19;
    }
    else if(!strcmp(order_gauss, "twentieth")) {
      gauss_order = 20;
    }
    else if(!strcmp(order_gauss, "twenty first")) {
      gauss_order = 21;
    }
    else if(!strcmp(order_gauss, "twenty second")) {
      gauss_order = 22;
    }
    else if(!strcmp(order_gauss, "twenty third")) {
      gauss_order = 23;
    }
    else if(!strcmp(order_gauss, "twenty fourth")) {
      gauss_order = 24;
    }
    else if(!strcmp(order_gauss, "twenty fifth")) {
      gauss_order = 25;
    }
    else if(!strcmp(order_gauss, "twenty sixth")) {
      gauss_order = 26;
    }
    else if(!strcmp(order_gauss, "twenty seventh")) {
      gauss_order = 27;
    }
    else if(!strcmp(order_gauss, "twenty eighth")) {
      gauss_order = 28;
    }
    else if(!strcmp(order_gauss, "twenty ninth")) {
      gauss_order = 29;
    }
    else if(!strcmp(order_gauss, "thirtieth")) {
      gauss_order = 30;
    }
    else if(!strcmp(order_gauss, "thirty first")) {
      gauss_order = 31;
    }
    else if(!strcmp(order_gauss, "thirty second")) {
      gauss_order = 32;
    }
    else if(!strcmp(order_gauss, "thirty third")) {
      gauss_order = 33;
    }
    else if(!strcmp(order_gauss, "thirty fourth")) {
      gauss_order = 34;
    }
    else if(!strcmp(order_gauss, "thirty fifth")) {
      gauss_order = 35;
    }
    else if(!strcmp(order_gauss, "thirty sixth")) {
      gauss_order = 36;
    }
    else if(!strcmp(order_gauss, "thirty seventh")) {
      gauss_order = 37;
    }
    else {
      std::cout << order_gauss << "is not a valid option for the Gauss points of" << geom_elem << std::endl;
      abort();
    }
    
    if(!strcmp(gauss_type, "legendre")) {
      if(!strcmp(geom_elem, "hex")) {
        GaussWeight = hex_gaussLegendre::Gauss[gauss_order];
        GaussPoints = hex_gaussLegendre::GaussPoints[gauss_order];
      }
      else if(!strcmp(geom_elem, "wedge")) {
        GaussWeight = wedge_gaussLegendre::Gauss[gauss_order];
        GaussPoints = wedge_gaussLegendre::GaussPoints[gauss_order];
      }
      else if(!strcmp(geom_elem, "tet")) {
        GaussWeight = tet_gaussLegendre::Gauss[gauss_order];
        GaussPoints = tet_gaussLegendre::GaussPoints[gauss_order];
      }
      else if(!strcmp(geom_elem, "quad")) {
        GaussWeight = quad_gaussLegendre::Gauss[gauss_order];
        GaussPoints = quad_gaussLegendre::GaussPoints[gauss_order];
      }
      else if(!strcmp(geom_elem, "tri")) {
        GaussWeight = tri_gaussLegendre::Gauss[gauss_order];
        GaussPoints = tri_gaussLegendre::GaussPoints[gauss_order];
      }
      else if(!strcmp(geom_elem, "line")) {
        GaussWeight = line_gaussLegendre::Gauss[gauss_order];
        GaussPoints = line_gaussLegendre::GaussPoints[gauss_order];
      }
      else if(!strcmp(geom_elem, "point")) {
        GaussWeight = point_gauss::Gauss[gauss_order];
        GaussPoints = point_gauss::GaussPoints[gauss_order];
      }
      else {
        std::cout << geom_elem << " is not a valid option" << std::endl;
        abort();
      }
    }
    else if(!strcmp(gauss_type, "lobatto")) {
      if(!strcmp(geom_elem, "hex")) {
        GaussWeight = hex_gaussLobatto::Gauss[gauss_order];
        GaussPoints = hex_gaussLobatto::GaussPoints[gauss_order];
      }
      else if(!strcmp(geom_elem, "wedge")) {
        GaussWeight = wedge_gaussLobatto::Gauss[gauss_order];
        GaussPoints = wedge_gaussLobatto::GaussPoints[gauss_order];
      }
      else if(!strcmp(geom_elem, "tet")) {
        GaussWeight = tet_gaussLobatto::Gauss[gauss_order];
        GaussPoints = tet_gaussLobatto::GaussPoints[gauss_order];
      }
      else if(!strcmp(geom_elem, "quad")) {
        GaussWeight = quad_gaussLobatto::Gauss[gauss_order];
        GaussPoints = quad_gaussLobatto::GaussPoints[gauss_order];
      }
      else if(!strcmp(geom_elem, "tri")) {
        GaussWeight = tri_gaussLobatto::Gauss[gauss_order];
        GaussPoints = tri_gaussLobatto::GaussPoints[gauss_order];
      }
      else if(!strcmp(geom_elem, "line")) {
        GaussWeight = line_gaussLobatto::Gauss[gauss_order];
        GaussPoints = line_gaussLobatto::GaussPoints[gauss_order];
      }
      else if(!strcmp(geom_elem, "point")) {
        GaussWeight = point_gauss::Gauss[gauss_order];
        GaussPoints = point_gauss::GaussPoints[gauss_order];
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
