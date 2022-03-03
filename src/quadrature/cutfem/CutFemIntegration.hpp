#ifndef __femus_cut_fem_int_hpp__
#define __femus_cut_fem_int_hpp__

#include "GeomElTypeEnum.hpp"

#include "CutLine.hpp"
#include "CutHyperCube.hpp"
#include "CutTriangle.hpp"
#include "CutTetrahedron.hpp"
#include "CutPrism.hpp"

#include "CutFem.hpp"

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

std::vector<std::string> geomName = {"hex", "tet", "wedge", "quad", "tri", "line", "point"};

template <class TypeIO, class TypeA>
class CutFemIntegral {
  public:
    CutFemIntegral(const GeomElType &geomElemType, const unsigned &qM, const std::string &gaussType) {
      _geomElemType = geomElemType;
      _qM = qM;
      _gaussType = gaussType;

      if(_geomElemType == HEX || _geomElemType == WEDGE || _geomElemType == TET) _dim = 3;
      if(_geomElemType == QUAD || _geomElemType == TRI) _dim = 2;
      if(_geomElemType == LINE) _dim = 1;
      this->build();
    }

    void build() {
      _gaussOrder = 2 * _qM;
      if(_geomElemType == HEX || _geomElemType == QUAD || _geomElemType == LINE) {
        _obj = new HCImap <TypeA, TypeA> (_dim, _qM, 0);
      }
      else if(_geomElemType == TET) {
        _obj = new TTImap <TypeA, TypeA> (_dim, _qM, 0);
      }
      else if(_geomElemType == TRI) {
        _obj = new TRImap <TypeA, TypeA> (_dim, _qM, 0);
      }
      else if(_geomElemType == WEDGE) {
      }
      else {
        std::cout << " No " << _geomElemType << " implemented" << std::endl;
        abort();
      }
      _gaussOrder = 2 * _qM;
      _gauss = new Gauss(_geomElemType, _gaussOrder, _gaussType.c_str());

      _gn = _gauss->GetGaussPointsNumber();
      //_xgp.resize(_dim);
      for(unsigned k = 0; k < _dim; k++) {
        _xgp[k] = _gauss->GetGaussCoordinatePointer(k);
      }

      _L = 1;
      for(unsigned idim = 0; idim < _dim; idim++) {
        _L *= (_qM + 1 + idim);
        _L /= (idim + 1);
      }
      Get_GS_ATA_Matrix(_geomElemType, _qM, _ATA, false);

      _f.resize(1, std::vector<cpp_bin_float_oct>(_L));
    }

    ~CutFemIntegral() {
      clear();
      delete _gauss;
      if(_geomElemType != WEDGE) delete _obj;
    };

    void clear() {
      if(_geomElemType != WEDGE) _obj->clear();
    };

    void operator()(const unsigned &qM, const int &s, const std::vector <TypeIO> &a, const TypeIO & d, std::vector <TypeIO> &weightCF, bool &wMap = 0);

    const double* GetGaussWeightPointer() {
      return _gauss->GetGaussWeightsPointer();
    };
    const double* GetGaussCoordinatePointer(const unsigned &k) {
      return _gauss->GetGaussCoordinatePointer(k);
    };


  protected:
    void polyBasis(const unsigned &qM, const std::vector<TypeIO> &x, std::vector<TypeIO> &bo);

  private:
    unsigned _dim;
    GeomElType _geomElemType;
    Gauss *_gauss;
    std::vector<std::vector<cpp_bin_float_oct>> _ATA;
    unsigned _L;

    unsigned _qM;
    std::string _gaussType;
    unsigned _gaussOrder;
    unsigned _gn;
    const double * _xgp[3];

    CutFEMmap <TypeA> *_obj;

    std::vector< std::vector<cpp_bin_float_oct> > _f;

    std::map < std::pair<std::vector<TypeIO>, TypeIO>, std::vector<TypeIO> > _WeightMap;
    typename std::map < std::pair<std::vector<TypeIO>, TypeIO>, std::vector<TypeIO> >::iterator _it;
    std::pair<std::vector<TypeIO>, TypeIO> _key;

};





template <class TypeIO, class TypeA>
void CutFemIntegral<TypeIO, TypeA>::operator()(const unsigned &qM, const int &s, const std::vector <TypeIO> &a, const TypeIO & d,  std::vector <TypeIO> &weightCF, bool &wMap) {

  if(_qM != qM) {
    clear();
    delete _gauss;
    _WeightMap.clear();
    _qM = qM;
    build();
  }

  if(wMap) {
    _key = std::make_pair(a, d);
    _it = _WeightMap.find(_key);
    if(_it != _WeightMap.end()) {
      weightCF.clear();
      weightCF = _WeightMap[_key];
      return;
    }
  }

  //std::vector< std::vector<TypeIO> > f(1, std::vector<TypeIO>(_L));
//   std::vector< std::vector<cpp_bin_float_oct> > f(1, std::vector<cpp_bin_float_oct>(_L));

  std::vector< TypeA > aA(_dim);
  for(unsigned i = 0; i < _dim; i++) aA[i] = static_cast<TypeA>(a[i]);
  TypeA dA = static_cast<TypeA>(d);


  unsigned count = 0;

  if(_dim == 3) {
    for(unsigned q = 0; q <= _qM; q++) {
      for(int ii = q; ii >= 0; ii--) {
        for(int jj = q - ii; jj >= 0; jj--) {
          unsigned i = static_cast<unsigned>(ii);
          unsigned j = static_cast<unsigned>(jj);
          unsigned k = q - i - j;
          _f[0][count] = static_cast<cpp_bin_float_oct>((_geomElemType == WEDGE) ?
                                                        Prism<TypeA, TypeA>(s, {i, j, k}, aA, dA) :
                                                        (*_obj)(s, {i, j, k}, aA, dA));
          count++;
        }
      }
    }

  }
  else if(_dim == 2) {
    for(unsigned q = 0; q <= _qM; q++) {
      for(unsigned j = 0; j <= q; j++) {
        unsigned i = q - j;
        _f[0][count] = static_cast<cpp_bin_float_oct>((*_obj)(s, {i, j}, aA, dA));
        count++;
      }
    }
  }
  else if(_dim == 1) {
    for(unsigned i = 0; i <= _qM; i++) {
      _f[0][i] = static_cast<cpp_bin_float_oct>((*_obj)(s, {i}, aA, dA));
    }
  }
  else {
    std::cout << " Dimension =" << _dim << " not admissible" << std::endl;
    abort();
  }

  std::vector< std::vector<TypeIO> > Co = MatrixMatrixMultiply<TypeIO, cpp_bin_float_oct>(_f, _ATA);

  std::vector<TypeIO> bo(_L);
  std::vector<TypeIO> x(_dim); //TODO

  unsigned TEType = 2;
  if(_geomElemType == TET) {
    TEType = (std::max(fabs(a[1] + a[0]), fabs(a[2] - a[1]))  >= fabs(a[0] - a[2])) ? 0 : 1;
  }



  weightCF.assign(_gn, 0);
  for(unsigned ig = 0; ig < _gn; ig++) {



    if(_geomElemType == LINE || _geomElemType == QUAD || _geomElemType == HEX) {
      for(unsigned k = 0; k < _dim; k++)  x[k] = 0.5 * (1. + _xgp[k][ig]);
    }
    else if(_geomElemType == TRI) {
      x[0] = static_cast<TypeIO>(1. - _xgp[0][ig]);
      x[1] = static_cast<TypeIO>(_xgp[1][ig]);
    }
    else if(_geomElemType == WEDGE) {
      x[0] = static_cast<TypeIO>(1. - _xgp[0][ig]);
      x[1] = static_cast<TypeIO>(_xgp[1][ig]);
      x[2] = static_cast<TypeIO>(0.5 * (1. + _xgp[2][ig]));
    }
    else {
      x[0] = static_cast<TypeIO>(_xgp[0][ig] + _xgp[1][ig] + _xgp[2][ig]);
      if(TEType == 0) {
        x[1] = static_cast<TypeIO>(_xgp[1][ig] + _xgp[2][ig]);
        x[2] = static_cast<TypeIO>(_xgp[2][ig]);
      }
      else {
        x[1] = static_cast<TypeIO>(_xgp[2][ig] + _xgp[0][ig]);
        x[2] = static_cast<TypeIO>(_xgp[0][ig]);
      }
    }

    polyBasis(_qM, x, bo);

    TypeA weight(0);
    for(unsigned i = 0; i < _L; i++) {
      weight += Co[0][i] * bo[i];
    }
    weightCF[ig] = static_cast<TypeIO>(weight);
  }

  if(wMap) {
    _WeightMap[_key] = weightCF;
  }
};


template <class TypeIO, class TypeA>
void CutFemIntegral<TypeIO, TypeA>::polyBasis(const unsigned & qM, const std::vector<TypeIO> &x, std::vector<TypeIO> &bo) {
  unsigned count = 0;

  if(_dim == 3) {

    std::vector<std::vector<std::vector<TypeIO>>> bOld(qM + 1);
    std::vector<std::vector<std::vector<TypeIO>>> bNew(qM + 1);
    for(unsigned q = 0; q <= qM; q++) {
      bOld[q].assign(qM + 1, std::vector<TypeIO>(qM + 1));
      bNew[q].assign(qM + 1, std::vector<TypeIO>(qM + 1));
    }

    bo[0] = bNew[0][0][0] = 1.;
    for(unsigned q = 1; q <= qM; q++) {
      bOld.swap(bNew);
      for(int ii = q; ii >= 0; ii--) {
        for(int jj = q - ii; jj >= 0; jj--) {
          count++;
          unsigned i = static_cast<unsigned>(ii);
          unsigned j = static_cast<unsigned>(jj);
          unsigned k = q - i - j;
          if(i > j) {
            if(i > k) {
              bo[count] = bNew[i][j][k] = bOld[i - 1][j][k] * x[0];
            }
            else {
              bo[count] = bNew[i][j][k] = bOld[i][j][k - 1] * x[2];
            }
          }
          else {
            if(j > k) {
              bo[count] = bNew[i][j][k] = bOld[i][j - 1][k] * x[1];
            }
            else {
              bo[count] = bNew[i][j][k] = bOld[i][j][k - 1] * x[2];
            }
          }

          //bo[count] = pow(x[0], i) * pow(x[1], j) * pow(x[2], k); //TO BE OPTIMIZED
        }
      }
    }
  }
  else if(_dim == 2) {
    std::vector<std::vector<TypeIO>> bOld(qM + 1, std::vector<TypeIO>(qM + 1));
    std::vector<std::vector<TypeIO>> bNew(qM + 1, std::vector<TypeIO>(qM + 1));
    bo[0] = bNew[0][0] = 1.;
    for(unsigned q = 1; q <= qM; q++) {
      bOld.swap(bNew);
      for(unsigned j = 0; j <= q; j++) {
        count++;
        unsigned i = q - j;
        if(i > j) bo[count] = bNew[i][j] = bOld[i - 1][j] * x[0];
        else bo[count] = bNew[i][j] = bOld[i][j - 1] * x[1];

      }
    }
  }
  else if(_dim == 1) {
    bo[0] = 1;
    for(unsigned i = 1; i <= qM; i++) {
      bo[i] = bo[i - 1] * x[0];
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

