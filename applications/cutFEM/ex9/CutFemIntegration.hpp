#include "GeomElTypeEnum.hpp"

#include "Line.hpp"
#include "Square.hpp"
#include "Cube.hpp"
#include "HyperCube.hpp"
#include "Triangle.hpp"
#include "Tetrahedron.hpp"
#include "Prism.hpp"

#include "GramSchmidt.hpp"

#include "GaussPoints.hpp"

#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>

namespace boost {
  namespace multiprecision {
    typedef number < backends::cpp_bin_float < 24, backends::digit_base_2, void, boost::int16_t, -126, 127 >, et_off >         cpp_bin_float_single;
    typedef number < backends::cpp_bin_float < 53, backends::digit_base_2, void, boost::int16_t, -1022, 1023 >, et_off >       cpp_bin_float_double;
    typedef number < backends::cpp_bin_float < 64, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >     cpp_bin_float_double_extended;
    typedef number < backends::cpp_bin_float < 113, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >    cpp_bin_float_quad;
    typedef number < backends::cpp_bin_float < 237, backends::digit_base_2, void, boost::int32_t, -262142, 262143 >, et_off >  cpp_bin_float_oct;
  }
} // namespaces

using boost::multiprecision::cpp_bin_float_oct;
using boost::multiprecision::cpp_bin_float_quad;

using namespace femus;




template <class Type>
class CutFEM_Integral {
  public:
    CutFEM_Integral(const char* geom_elem, const unsigned &qM, const char* gauss_type) {
      if ( !strcmp(geom_elem, "hex") )         _GeomElemType = HEX;
      else if ( !strcmp(geom_elem, "tet") )    _GeomElemType = TET;
      else if ( !strcmp(geom_elem, "wedge") )  _GeomElemType = WEDGE;
      else if ( !strcmp(geom_elem, "quad") )   _GeomElemType = QUAD;
      else if ( !strcmp(geom_elem, "tri") )    _GeomElemType = TRI;
      else if ( !strcmp(geom_elem, "line") )   _GeomElemType = LINE;
      else {
        std::cout << " No " << geom_elem << " implemented" << std::endl;
        abort();
      }
      const char* order_gauss = 2 * qM;
      _gauss = new Gauss(geom_elem, order_gauss, gauss_type);
    }

    ~CutFEM_Integral() {
      clear();
    };

    void clear() {

//         TODO
    };

    Type operator()(const GeomElType &geom, const unsigned &qM, const char* gauss_type, const int &s, const std::vector<Type> &a, const Type &d);


  protected:
//     void foComputation( std::vector< std::vector<Type> > &f, const unsigned &qM );

  private:
    GeomElType _GeomElemType;
    const Gauss _gauss;
    std::vector<std::vector<Type>> _ATA;
};

std::vector<std::string> geomName = {"hex", "tet", "wedge", "quad", "tri", "line", "point"};

template <class Type>
Type CutFEM_Integral<Type>::operator()(const GeomElType &geom, const unsigned &qM, const char* gauss_type, const int &s, const std::vector <Type> &a, const Type & d) {

  std::vector< std::vector<Type> > f(1, std::vector<Type>((qM + 1) * (qM + 2) * (qM + 3) / 6));

  const unsigned dim = 3; //TODO put here the function to have the space dim from geom_el
//     TODO very bad programming here, put something more sofisticated.
  if ( geom == HEX || geom == QUAD || geom = LINE ) {
    HCImap <Type, Type> obj(dim, qM, 0);
//       TODO is it possible to put HCImap, TTImap, ... as variables? In this way we can pass it to foComputation.
  }
  else if ( geom == TET ) {
    TTImap <Type, Type> obj(dim, qM, 0);
  }
  else if ( geom == TRI ) {
    TRImap <Type, Type> obj(dim, qM, 0);
  }
  else if( geom == WEDGE ) {
  }
  else {
    std::cout << " No " << geom << " implemented" << std::endl;
    abort();
  }

//   TODO this should be replaced with the function foComputation if possible
  {

    unsigned count = 0;

    if(dim == 3) {
      for(unsigned q = 0; q <= qM; q++) {
        for(int ii = q; ii >= 0; ii--) {
          for(int jj = q - ii; jj >= 0; jj--) {
            unsigned i = static_cast<unsigned>(ii);
            unsigned j = static_cast<unsigned>(jj);
            unsigned k = q - i - j;
            if( geom == WEDGE ) Prism<Type, Type>(s, {i, j, k}, a, d);
            else f[0][count] = obj(s, {i, j, k}, a, d);
            count++;
          }
        }
      }

    }
    else if (dim == 2) {
      for(unsigned q = 0; q <= qM; q++) {
        for(unsigned j = 0; j <= q; j++) {
          unsigned i = q - j;
          f[0][count] = obj(s, {i, j}, a, d);
          count++;
        }
      }
    }
    else if (dim == 1) {
        for(unsigned i = 0; i <= qM; i++) {
      f[0][i] = obj(s, {i}, a, d);
    }
    }
    else {
      std::cout << " Dimension =" << dim << " not admissible" << std::endl;
      abort();
    }
  }

  if(_ATA.dimension != (qM + 1) * (qM + 1)) {
    Get_GS_ATA_Matrix(_GeomElemType, qM, _ATA, false);
  }

  std::vector< std::vector<Type> > Co = MatrixMatrixMultiply(f, _ATA);


};

// template <class Type>
// void CutFEM_Integral<Type>::foComputation( std::vector< std::vector<Type> > &f, const unsigned &qM ) {
// 
//   unsigned count = 0;
// 
//   for(unsigned q = 0; q <= qM; q++) {
//     for(int ii = q; ii >= 0; ii--) {
//       for(int jj = q - ii; jj >= 0; jj--) {
//         unsigned i = static_cast<unsigned>(ii);
//         unsigned j = static_cast<unsigned>(jj);
//         unsigned k = q - i - j;
// //         f[0][count] = /*tet3*/(-1, {i, j, k}, {0/sqrt(2), 1/sqrt(2), -1/sqrt(2)}, 0); TODO read above
//         count++;
//       }
//     }
//   }
// };

