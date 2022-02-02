#ifndef __femus_cut_fem_int_hpp__
#define __femus_cut_fem_int_hpp__

#include "GeomElTypeEnum.hpp"

#include "Line.hpp"
#include "Square.hpp"
#include "Cube.hpp"
#include "HyperCube.hpp"
#include "Triangle.hpp"
#include "Tetrahedron.hpp"
#include "Prism.hpp"

#include "cutFem.hpp"

// #include "GramSchmidt.hpp"

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
class CutFemIntegral {
  public:
    CutFemIntegral(const char* geom_elem, const unsigned &qM, const char* gauss_type) {
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


      if ( _GeomElemType == HEX || _GeomElemType == _GeomElemType || _GeomElemType = TET) _dim = 3;
      if ( _GeomElemType == QUAD || _GeomElemType == TRI ) _dim = 2;
      if ( _GeomElemType == LINE) _dim = 1;

      if ( _GeomElemType == HEX || _GeomElemType == QUAD || _GeomElemType = LINE ) {
        *_obj = new HCImap <Type, Type> (_dim, qM, 0);
      }
      else if ( _GeomElemType == TET ) {
        *_obj = new TTImap <Type, Type> (_dim, qM, 0);
      }
      else if ( _GeomElemType == TRI ) {
        *_obj = new TRImap <Type, Type> (_dim, qM, 0);
      }
      else if( _GeomElemType == WEDGE ) {
      }
      else {
        std::cout << " No " << _GeomElemType << " implemented" << std::endl;
        abort();
      }

      _orderGauss = 2 * qM;
      _gauss = new Gauss(geom_elem, _orderGauss, gauss_type);

      _gn = _gauss.GetGaussPointsNumber();
      for(unsigned k = 0; k < _dim; k++) {
        _xgp[k] = _gauss.GetGaussCoordinatePointer(k);
      }

      _L = 1;
      for( unsigned idim = 0; idim < _dim; idim++) _L *= ( (qM + 1 + idim) / (idim + 1) );

    }

    ~CutFemIntegral() {
      clear();
    };

    void clear() {

//         TODO
    };

    Type operator()(const GeomElType &geom, const unsigned &qM, const char* gauss_type, const int &s, const std::vector<Type> &a, const Type &d);


  protected:
//     void foComputation( std::vector< std::vector<Type> > &f, const unsigned &qM );
    void polyBasis( std::vector<Type> &bo, const unsigned &qM, std::vector<Type> GaussCoords );

  private:
    unsigned _dim;
    GeomElType _GeomElemType;
    const Gauss _gauss;
    std::vector<std::vector<Type>> _ATA;
    unsigned _L;

    const char* _orderGauss;
    unsigned _gn;
    double *_xgp;

    cutFEMmap <Type, Type> *_obj;

};

std::vector<std::string> geomName = {"hex", "tet", "wedge", "quad", "tri", "line", "point"};

template <class Type>
Type CutFemIntegral<Type>::operator()(const GeomElType &geom, const unsigned &qM, const char* gauss_type, const int &s, const std::vector <Type> &a, const Type & d) {





  std::vector< std::vector<Type> > f(1, std::vector<Type>(_L));

  unsigned count = 0;

  if(_dim == 3) {
    for(unsigned q = 0; q <= qM; q++) {
      for(int ii = q; ii >= 0; ii--) {
        for(int jj = q - ii; jj >= 0; jj--) {
          unsigned i = static_cast<unsigned>(ii);
          unsigned j = static_cast<unsigned>(jj);
          unsigned k = q - i - j;
          if( geom == WEDGE ) Prism<Type, Type>(s, {i, j, k}, a, d);
          else f[0][count] = (*_obj)(s, {i, j, k}, a, d);
          count++;
        }
      }
    }

  }
  else if (_dim == 2) {
    for(unsigned q = 0; q <= qM; q++) {
      for(unsigned j = 0; j <= q; j++) {
        unsigned i = q - j;
        f[0][count] = (*_obj)(s, {i, j}, a, d);
        count++;
      }
    }
  }
  else if (_dim == 1) {
    for(unsigned i = 0; i <= qM; i++) {
      f[0][i] = (*_obj)(s, {i}, a, d);
    }
  }
  else {
    std::cout << " Dimension =" << _dim << " not admissible" << std::endl;
    abort();
  }

  if(_ATA.dimension != _L) { // TODO verify L
    Get_GS_ATA_Matrix(_GeomElemType, qM, _ATA, false);
  }

  std::vector< std::vector<Type> > Co = MatrixMatrixMultiply(f, _ATA);

  std::vector<Type> bo;
  std::vector<Type> GaussCoords(_dim); //TODO

  std::vector <Type> weightCF;
  weightCF.assign(_gn, 0);

  for( unsigned ig = 0; ig < _gn; ig++) {
    for ( unsigned idim = 0; idim < _dim; idim++ ) GaussCoords[idim] = *(_xgp[idim] + ig);  //TODO verify if it is correct
    polyBasis( bo, qM, _dim, GaussCoords );
    for( unsigned i = 0; i <= _L; i++ ) {
      weightCF[ig] += Co[0][i] * bo[i]; //TODO
    }
  }


};


template <class Type>
void CutFemIntegral<Type>::polyBasis( std::vector<Type> &bo, const unsigned &qM, std::vector<Type> GaussCoords ) {
  unsigned count = 0;

  if(_dim == 3) {
    for(unsigned q = 0; q <= qM; q++) {
      for(int ii = q; ii >= 0; ii--) {
        for(int jj = q - ii; jj >= 0; jj--) {
          unsigned i = static_cast<unsigned>(ii);
          unsigned j = static_cast<unsigned>(jj);
          unsigned k = q - i - j;
          bo[count] = pow(GaussCoords[0], i) * pow(GaussCoords[1], j) * pow(GaussCoords[2], k); //TODO pay attention here, to verify
          count++;
        }
      }
    }

  }
  else if (_dim == 2) {
    for(unsigned q = 0; q <= qM; q++) {
      for(unsigned j = 0; j <= q; j++) {
        unsigned i = q - j;
        bo[count] = pow(GaussCoords[0], i) * pow(GaussCoords[1], j);
        count++;
      }
    }
  }
  else if (_dim == 1) {
    for(unsigned i = 0; i <= qM; i++) {
      bo[i] = pow(GaussCoords[0], i);
    }
  }
  else {
    std::cout << " Dimension =" << _dim << " not admissible" << std::endl;
    abort();
  }

}

// template <class Type>
// void CutFemIntegral<Type>::foComputation( std::vector< std::vector<Type> > &f, const unsigned &qM ) {
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


#endif

