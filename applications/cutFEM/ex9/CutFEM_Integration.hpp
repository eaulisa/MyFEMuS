#include "GeomElTypeEnum.hpp"

#include "Line.hpp"
#include "Square.hpp"
#include "Cube.hpp"
#include "HyperCube.hpp"
#include "Triangle.hpp"
#include "Tetrahedron.hpp"
#include "Prism.hpp" 

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
     CutFEM_Integral(const char* geom_elem, const unsigned &qM, const char* gauss_type){
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
    
    Type operator()(const char* geom_elem, const unsigned &qM, const char* gauss_type, const int &s, const std::vector<unsigned> &m, const std::vector<Type> &a, const Type &d);
    
    
protected:
    std::vector< std::vector<Type> > foComputation();  
    
private: 
    GeomElType _GeomElemType;
    const Gauss _gauss;
//     Gauss* _gauss;
};

template <class Type>
Type CutFEM_Integral<Type>::operator()(const char* geom_elem, const unsigned &qM, const char* gauss_type, const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type & d) {
    
    
    
    std::vector< std::vector<Type> > f(1, std::vector<Type>((qM + 1) * (qM + 2) * (qM + 3) / 6));
};

template <class Type>
std::vector< std::vector<Type> > CutFEM_Integral<Type>::foComputation() {
    
};
