#ifndef __femus_fem_hpp__
#define __femus_fem_hpp__

#include "ElemType.hpp"


namespace femus {

  class Fem {
    public:
      Fem(const unsigned &gaussOrder, const unsigned &dim) {
        _dim = dim;
        _gaussOrder = gaussOrder;
        build();
      }

      void resize(const unsigned &gaussOrder, const unsigned &dim) {
        if(gaussOrder != _gaussOrder || _dim != dim) {
          clear();
          _dim = dim;
          _gaussOrder = gaussOrder;
          build();
        }
      }

      void build() {
        if(_dim > 3) {
          std::cout << "wrong dimension " << _dim << " in Fem class\n";
          abort();
        }

        if(_gaussOrder > 37) {
          std::cout << "wrong Gauss Order " << _gaussOrder << " in Fem class\n";
          abort();
        }


        if(_dim > 2) {
          _finiteElement[0][0] = new const elem_type_3D("hex", "linear", numberName[_gaussOrder].c_str());
          _finiteElement[0][1] = new const elem_type_3D("hex", "quadratic", numberName[_gaussOrder].c_str());
          _finiteElement[0][2] = new const elem_type_3D("hex", "biquadratic", numberName[_gaussOrder].c_str());
          _finiteElement[0][3] = new const elem_type_3D("hex", "constant", numberName[_gaussOrder].c_str());
          _finiteElement[0][4] = new const elem_type_3D("hex", "disc_linear", numberName[_gaussOrder].c_str());


          _finiteElement[1][0] = new const elem_type_3D("tet", "linear", numberName[_gaussOrder].c_str());
          _finiteElement[1][1] = new const elem_type_3D("tet", "quadratic", numberName[_gaussOrder].c_str());
          _finiteElement[1][2] = new const elem_type_3D("tet", "biquadratic", numberName[_gaussOrder].c_str());
          _finiteElement[1][3] = new const elem_type_3D("tet", "constant", numberName[_gaussOrder].c_str());
          _finiteElement[1][4] = new const elem_type_3D("tet", "disc_linear", numberName[_gaussOrder].c_str());


          _finiteElement[2][0] = new const elem_type_3D("wedge", "linear", numberName[_gaussOrder].c_str());
          _finiteElement[2][1] = new const elem_type_3D("wedge", "quadratic", numberName[_gaussOrder].c_str());
          _finiteElement[2][2] = new const elem_type_3D("wedge", "biquadratic", numberName[_gaussOrder].c_str());
          _finiteElement[2][3] = new const elem_type_3D("wedge", "constant", numberName[_gaussOrder].c_str());
          _finiteElement[2][4] = new const elem_type_3D("wedge", "disc_linear", numberName[_gaussOrder].c_str());
        }
        if(_dim > 1) {
          _finiteElement[3][0] = new const elem_type_2D("quad", "linear", numberName[_gaussOrder].c_str());
          _finiteElement[3][1] = new const elem_type_2D("quad", "quadratic", numberName[_gaussOrder].c_str());
          _finiteElement[3][2] = new const elem_type_2D("quad", "biquadratic", numberName[_gaussOrder].c_str());
          _finiteElement[3][3] = new const elem_type_2D("quad", "constant", numberName[_gaussOrder].c_str());
          _finiteElement[3][4] = new const elem_type_2D("quad", "disc_linear", numberName[_gaussOrder].c_str());

          _finiteElement[4][0] = new const elem_type_2D("tri", "linear", numberName[_gaussOrder].c_str());
          _finiteElement[4][1] = new const elem_type_2D("tri", "quadratic", numberName[_gaussOrder].c_str());
          _finiteElement[4][2] = new const elem_type_2D("tri", "biquadratic", numberName[_gaussOrder].c_str());
          _finiteElement[4][3] = new const elem_type_2D("tri", "constant", numberName[_gaussOrder].c_str());
          _finiteElement[4][4] = new const elem_type_2D("tri", "disc_linear", numberName[_gaussOrder].c_str());
        }

        _finiteElement[5][0] = new const elem_type_1D("line", "linear", numberName[_gaussOrder].c_str());
        _finiteElement[5][1] = new const elem_type_1D("line", "quadratic", numberName[_gaussOrder].c_str());
        _finiteElement[5][2] = new const elem_type_1D("line", "biquadratic", numberName[_gaussOrder].c_str());
        _finiteElement[5][3] = new const elem_type_1D("line", "constant", numberName[_gaussOrder].c_str());
        _finiteElement[5][4] = new const elem_type_1D("line", "disc_linear", numberName[_gaussOrder].c_str());
      }
      ~Fem() {
        clear();
      }
      void clear() {
        if(_dim > 2) {
          for(unsigned i = 0; i < 3; i++) {
            for(unsigned j = 0; j < 5; j++) {
              delete _finiteElement[i][j];
            }
          }
        }

        if(_dim > 1) {
          for(unsigned i = 3; i < 5; i++) {
            for(unsigned j = 0; j < 5; j++) {
              delete _finiteElement[i][j];
            }
          }
        }

        for(unsigned j = 0; j < 5; j++) {
          delete _finiteElement[5][j];
        }
      }
      
      unsigned GetGaussOrder(){ return _gaussOrder;}
      const elem_type * GetFiniteElement(const unsigned &geom, const unsigned &type ){ return _finiteElement[geom][type];}

    private:
      unsigned _gaussOrder;
      unsigned _dim;
      const elem_type *_finiteElement[6][5];

  };
}

#endif
