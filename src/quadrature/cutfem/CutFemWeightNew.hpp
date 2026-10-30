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

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

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

static const std::vector<std::string> geomName = {"hex", "tet", "wedge", "quad", "tri", "line", "point"};

template <class TypeIO, class TypeA>
class CutFemWeight {
  public:
    CutFemWeight(const GeomElType &geomElemType, const unsigned &qM, const std::string &gaussType)
      : _geomElemType(geomElemType),
        _gauss(NULL),
        _qM(qM),
        _gaussType(gaussType),
        _cntCall(0),
        _obj(NULL),
        _tetBaseType(0) {

      if(_geomElemType == HEX || _geomElemType == WEDGE || _geomElemType == TET) _dim = 3;
      else if(_geomElemType == QUAD || _geomElemType == TRI) _dim = 2;
      else if(_geomElemType == LINE) _dim = 1;
      else {
        std::cout << " No " << _geomElemType << " implemented" << std::endl;
        abort();
      }

      build();
    }

    void build() {

      _cntCall = 0;

      if(_geomElemType == HEX || _geomElemType == QUAD || _geomElemType == LINE) {
        _obj = new HCImap<TypeA, TypeA>(_dim, _qM, 0);
      }
      else if(_geomElemType == TET) {
        _obj = new TTImap<TypeA, TypeA>(_dim, _qM, 0);
      }
      else if(_geomElemType == TRI) {
        _obj = new TRImap<TypeA, TypeA>(_dim, _qM, 0);
      }
      else if(_geomElemType == WEDGE) {
        _obj = NULL;
      }
      else {
        std::cout << " No " << _geomElemType << " implemented" << std::endl;
        abort();
      }

      _gaussOrder = 2 * _qM;
      _gauss = new Gauss(_geomElemType, _gaussOrder, _gaussType.c_str());

      _gn = _gauss->GetGaussPointsNumber();

      for(unsigned k = 0; k < _dim; k++) {
        _xgp[k] = _gauss->GetGaussCoordinatePointer(k);
      }

      _L = 1;
      for(unsigned idim = 0; idim < _dim; idim++) {
        _L *= (_qM + 1 + idim);
        _L /= (idim + 1);
      }

      std::vector<std::vector<cpp_bin_float_oct> > ATAoct;
      Get_GS_ATA_Matrix(_geomElemType, _qM, ATAoct, false);

      _ATA.resize(ATAoct.size(), ATAoct.size());
      for(unsigned i = 0; i < ATAoct.size(); i++) {
        for(unsigned j = 0; j < ATAoct[i].size(); j++) {
          _ATA(i, j) = static_cast<TypeA>(ATAoct[i][j]);
        }
      }

      _f.resize(_L);
      _Co.resize(_L);

      // B0 depends only on the element, quadrature rule, polynomial order,
      // and tetrahedron orientation.  Tetrahedra need both orientations;
      // every other element needs one matrix.
      const unsigned numberOfB0Matrices = (_geomElemType == TET) ? 2u : 1u;
      _B0.resize(numberOfB0Matrices);

      std::vector<double> bo(_L);
      std::vector<double> x(_dim);

      for(unsigned iMatrix = 0; iMatrix < numberOfB0Matrices; iMatrix++) {
        const unsigned TETtype = iMatrix + 1u;

        _B0[iMatrix].resize(_gn, _L);

        for(unsigned ig = 0; ig < _gn; ig++) {

          if(_geomElemType == LINE || _geomElemType == QUAD || _geomElemType == HEX) {
            for(unsigned k = 0; k < _dim; k++) {
              x[k] = 0.5 * (1. + _xgp[k][ig]);
            }
          }
          else if(_geomElemType == TRI) {
            x[0] = 1. - _xgp[0][ig];
            x[1] = _xgp[1][ig];
          }
          else if(_geomElemType == WEDGE) {
            x[0] = 1. - _xgp[0][ig];
            x[1] = _xgp[1][ig];
            x[2] = 0.5 * (1. + _xgp[2][ig]);
          }
          else { // TET
            x[0] = _xgp[0][ig] + _xgp[1][ig] + _xgp[2][ig];
            if(TETtype == 1) {
              x[1] = _xgp[1][ig] + _xgp[2][ig];
              x[2] = _xgp[2][ig];
            }
            else { // TETtype == 2
              x[1] = _xgp[2][ig] + _xgp[0][ig];
              x[2] = _xgp[0][ig];
            }
          }

          PolyBasis(x, bo);

          for(unsigned i = 0; i < _L; i++) {
            _B0[iMatrix](ig, i) = static_cast<TypeIO>(bo[i]);
          }
        }
      }
    }

    ~CutFemWeight() {
      clear();
      delete _gauss;
      delete _obj;
    }

    void clear() {
      if(_obj != NULL) _obj->clear();
    }

    void operator()(const int &s, const std::vector<TypeIO> &a, const TypeIO &d,
                    std::vector<TypeIO> &weightCF);

    void GetWeightWithMap(const int &s, const std::vector<TypeIO> &a, const TypeIO &d,
                          std::vector<TypeIO> &weightCF);

    void UpdateQuadratureRule(const unsigned &qM) {
      clear();
      delete _gauss;
      delete _obj;
      _gauss = NULL;
      _obj = NULL;
      ClearMap();

      _qM = qM;
      build();
    }

    unsigned GetCutFEMQuadratureOrder() {
      return _qM;
    }

    unsigned GetGaussQuadratureOrder() {
      return _gaussOrder;
    }

    unsigned GetDimension() {
      return _dim;
    }

    unsigned GetGaussQuadraturePointNumber() {
      return _gn;
    }

    const double* GetGaussWeightPointer() {
      return _gauss->GetGaussWeightsPointer();
    }

    const double* GetGaussCoordinatePointer(const unsigned &k) {
      return _gauss->GetGaussCoordinatePointer(k);
    }

    const unsigned& GetCounter() {
      return _cntCall;
    }

    void ClearMap() {
      for(unsigned i = 0; i < _WeightMap.size(); i++) _WeightMap[i].clear();
      _WeightMap.clear();
    }

    void SetTetBaseType(const unsigned &value) {
      _tetBaseType = value;
    }

  protected:
    void PolyBasis(const std::vector<double> &x, std::vector<double> &bo);

  private:
    unsigned _dim;
    GeomElType _geomElemType;
    Gauss *_gauss;

    typedef Eigen::Matrix<TypeIO, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> B0Matrix;

    Eigen::Matrix<TypeA, Eigen::Dynamic, Eigen::Dynamic> _ATA;
    std::vector<B0Matrix> _B0;
    Eigen::Matrix<TypeA, Eigen::Dynamic, 1> _f;
    Eigen::Matrix<TypeA, Eigen::Dynamic, 1> _Co;
    unsigned _L;

    unsigned _qM;
    std::string _gaussType;
    unsigned _gaussOrder;
    unsigned _gn;
    const double *_xgp[3];

    unsigned _cntCall;

    CutFEMmap<TypeA> *_obj;

    unsigned _tetBaseType;

    std::vector<std::map<std::pair<std::vector<float>, float>, std::vector<TypeIO> > > _WeightMap;
    typename std::map<std::pair<std::vector<float>, float>, std::vector<TypeIO> >::iterator _it;
    std::pair<std::vector<float>, float> _key;

    std::vector<std::vector<TypeIO> > _bOld2d;
    std::vector<std::vector<TypeIO> > _bNew2d;

    std::vector<std::vector<std::vector<TypeIO> > > _bOld3d;
    std::vector<std::vector<std::vector<TypeIO> > > _bNew3d;
};

template <class TypeIO, class TypeA>
void CutFemWeight<TypeIO, TypeA>::GetWeightWithMap(const int &s,
    const std::vector<TypeIO> &a,
    const TypeIO &d,
    std::vector<TypeIO> &weightCF) {
  _WeightMap.resize(s + 2);
  std::vector<float> af(a.begin(), a.end());
  const float df = static_cast<float>(d);
  _key = std::make_pair(af, df);
  _it = _WeightMap[s + 1].find(_key);

  if(_it != _WeightMap[s + 1].end()) {
    weightCF = _it->second;
    return;
  }

  operator()(s, a, d, weightCF);
  _WeightMap[s + 1][_key] = weightCF;
}

template <class TypeIO, class TypeA>
void CutFemWeight<TypeIO, TypeA>::operator()(const int &s,
    const std::vector<TypeIO> &a,
    const TypeIO &d,
    std::vector<TypeIO> &weightCF) {

  _cntCall++;

  std::vector<TypeA> aA(_dim);
  for(unsigned i = 0; i < _dim; i++) {
    aA[i] = static_cast<TypeA>(a[i]);
  }
  const TypeA dA = static_cast<TypeA>(d);

  if(_obj != NULL) _obj->SetBaseType(_tetBaseType);

  unsigned count = 0;

  if(_dim == 3) {
    for(unsigned q = 0; q <= _qM; q++) {
      for(int ii = q; ii >= 0; ii--) {
        for(int jj = q - ii; jj >= 0; jj--) {
          const unsigned i = static_cast<unsigned>(ii);
          const unsigned j = static_cast<unsigned>(jj);
          const unsigned k = q - i - j;
          _f[count] = (_geomElemType == WEDGE)
                      ? Prism<TypeA, TypeA>(s, {i, j, k}, aA, dA)
                      : (*_obj)(s, {i, j, k}, aA, dA);
          count++;
        }
      }
    }
  }
  else if(_dim == 2) {
    for(unsigned q = 0; q <= _qM; q++) {
      for(unsigned j = 0; j <= q; j++) {
        const unsigned i = q - j;
        _f[count] = (*_obj)(s, {i, j}, aA, dA);
        count++;
      }
    }
  }
  else if(_dim == 1) {
    for(unsigned i = 0; i <= _qM; i++) {
      _f[i] = (*_obj)(s, {i}, aA, dA);
    }
  }
  else {
    std::cout << " Dimension =" << _dim << " not admissible" << std::endl;
    abort();
  }

  _Co.noalias() = _ATA * _f;

  unsigned B0index = 0;

  if(_geomElemType == TET) {
    unsigned TETtype = _tetBaseType;

    if(TETtype == 0) {
      TETtype = (std::max(std::fabs(a[1] - a[0]), std::fabs(a[2] - a[1]))
                 >= std::fabs(a[0] - a[2])) ? 1u : 2u;
    }

    if(TETtype != 1 && TETtype != 2) {
      std::cout << " TET base type = " << TETtype << " not admissible" << std::endl;
      abort();
    }

    B0index = TETtype - 1u;
  }

  weightCF.resize(_gn);
  for(unsigned ig = 0; ig < _gn; ig++) {
    const TypeIO *b0Row = _B0[B0index].data() + ig * _L;
    TypeA weight = 0.;
    for(unsigned i = 0; i < _L; i++) {
      weight += b0Row[i] * _Co[i];
    }
    weightCF[ig] = static_cast<TypeIO>(weight);
  }
}

template <class TypeIO, class TypeA>
void CutFemWeight<TypeIO, TypeA>::PolyBasis(const std::vector<double> &x,
    std::vector<double> &bo) {
  unsigned count = 0;

  if(_dim == 3) {
    _bOld3d.resize(_qM + 1);
    _bNew3d.resize(_qM + 1);
    for(unsigned q = 0; q <= _qM; q++) {
      _bOld3d[q].assign(_qM + 1, std::vector<TypeIO>(_qM + 1));
      _bNew3d[q].assign(_qM + 1, std::vector<TypeIO>(_qM + 1));
    }

    bo[0] = _bNew3d[0][0][0] = 1.;
    for(unsigned q = 1; q <= _qM; q++) {
      _bOld3d.swap(_bNew3d);
      for(int ii = q; ii >= 0; ii--) {
        for(int jj = q - ii; jj >= 0; jj--) {
          count++;
          const unsigned i = static_cast<unsigned>(ii);
          const unsigned j = static_cast<unsigned>(jj);
          const unsigned k = q - i - j;

          if(i > j) {
            if(i > k) {
              bo[count] = _bNew3d[i][j][k] = _bOld3d[i - 1][j][k] * x[0];
            }
            else {
              bo[count] = _bNew3d[i][j][k] = _bOld3d[i][j][k - 1] * x[2];
            }
          }
          else {
            if(j > k) {
              bo[count] = _bNew3d[i][j][k] = _bOld3d[i][j - 1][k] * x[1];
            }
            else {
              bo[count] = _bNew3d[i][j][k] = _bOld3d[i][j][k - 1] * x[2];
            }
          }
        }
      }
    }
  }
  else if(_dim == 2) {
    _bOld2d.assign(_qM + 1, std::vector<TypeIO>(_qM + 1));
    _bNew2d.assign(_qM + 1, std::vector<TypeIO>(_qM + 1));
    bo[0] = _bNew2d[0][0] = 1.;

    for(unsigned q = 1; q <= _qM; q++) {
      _bOld2d.swap(_bNew2d);
      for(unsigned j = 0; j <= q; j++) {
        count++;
        const unsigned i = q - j;
        if(i > j) bo[count] = _bNew2d[i][j] = _bOld2d[i - 1][j] * x[0];
        else bo[count] = _bNew2d[i][j] = _bOld2d[i][j - 1] * x[1];
      }
    }
  }
  else if(_dim == 1) {
    bo[0] = 1.;
    for(unsigned i = 1; i <= _qM; i++) {
      bo[i] = bo[i - 1] * x[0];
    }
  }
  else {
    std::cout << " Dimension =" << _dim << " not admissible" << std::endl;
    abort();
  }
}

#endif

// #ifndef __femus_cut_fem_int_hpp__
// #define __femus_cut_fem_int_hpp__
//
// #include "GeomElTypeEnum.hpp"
//
// #include "CutLine.hpp"
// #include "CutHyperCube.hpp"
// #include "CutTriangle.hpp"
// #include "CutTetrahedron.hpp"
// #include "CutPrism.hpp"
//
// #include "CutFem.hpp"
//
// #include "GramSchmidt.hpp"
//
// #include "GaussPoints.hpp"
//
// #include <Eigen/Dense>
//
// #include <algorithm>
// #include <cmath>
// #include <cstdlib>
// #include <iostream>
// #include <map>
// #include <string>
// #include <utility>
// #include <vector>
//
// #include <boost/multiprecision/cpp_bin_float.hpp>
//
// namespace boost {
//   namespace multiprecision {
//     typedef number < backends::cpp_bin_float < 24, backends::digit_base_2, void, boost::int16_t, -126, 127 >, et_off >         cpp_bin_float_single;
//     typedef number < backends::cpp_bin_float < 53, backends::digit_base_2, void, boost::int16_t, -1022, 1023 >, et_off >       cpp_bin_float_double;
//     typedef number < backends::cpp_bin_float < 64, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >     cpp_bin_float_double_extended;
//     typedef number < backends::cpp_bin_float < 113, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >    cpp_bin_float_quad;
//     typedef number < backends::cpp_bin_float < 237, backends::digit_base_2, void, boost::int32_t, -262142, 262143 >, et_off >  cpp_bin_float_oct;
//   }
// } // namespaces
//
// using boost::multiprecision::cpp_bin_float_oct;
// using boost::multiprecision::cpp_bin_float_quad;
//
// using namespace femus;
//
// static const std::vector<std::string> geomName = {"hex", "tet", "wedge", "quad", "tri", "line", "point"};
//
// template <class TypeIO, class TypeA>
// class CutFemWeight {
//   public:
//     CutFemWeight(const GeomElType &geomElemType, const unsigned &qM, const std::string &gaussType)
//       : _geomElemType(geomElemType),
//         _gauss(NULL),
//         _qM(qM),
//         _gaussType(gaussType),
//         _cntCall(0),
//         _obj(NULL),
//         _tetBaseType(0) {
//
//       if (_geomElemType == HEX || _geomElemType == WEDGE || _geomElemType == TET) _dim = 3;
//       else if (_geomElemType == QUAD || _geomElemType == TRI) _dim = 2;
//       else if (_geomElemType == LINE) _dim = 1;
//       else {
//         std::cout << " No " << _geomElemType << " implemented" << std::endl;
//         abort();
//       }
//
//       build();
//     }
//
//     void build() {
//
//       _cntCall = 0;
//
//       if (_geomElemType == HEX || _geomElemType == QUAD || _geomElemType == LINE) {
//         _obj = new HCImap<TypeA, TypeA>(_dim, _qM, 0);
//       }
//       else if (_geomElemType == TET) {
//         _obj = new TTImap<TypeA, TypeA>(_dim, _qM, 0);
//       }
//       else if (_geomElemType == TRI) {
//         _obj = new TRImap<TypeA, TypeA>(_dim, _qM, 0);
//       }
//       else if (_geomElemType == WEDGE) {
//         _obj = NULL;
//       }
//       else {
//         std::cout << " No " << _geomElemType << " implemented" << std::endl;
//         abort();
//       }
//
//       _gaussOrder = 2 * _qM;
//       _gauss = new Gauss(_geomElemType, _gaussOrder, _gaussType.c_str());
//
//       _gn = _gauss->GetGaussPointsNumber();
//
//       for (unsigned k = 0; k < _dim; k++) {
//         _xgp[k] = _gauss->GetGaussCoordinatePointer(k);
//       }
//
//       _L = 1;
//       for (unsigned idim = 0; idim < _dim; idim++) {
//         _L *= (_qM + 1 + idim);
//         _L /= (idim + 1);
//       }
//
//       std::vector<std::vector<cpp_bin_float_oct> > ATAoct;
//       Get_GS_ATA_Matrix(_geomElemType, _qM, ATAoct, false);
//
//       Eigen::Matrix<TypeA, Eigen::Dynamic, Eigen::Dynamic> ATA(ATAoct.size(), ATAoct.size());
//       for (unsigned i = 0; i < ATAoct.size(); i++) {
//         for (unsigned j = 0; j < ATAoct[i].size(); j++) {
//           ATA(i, j) = static_cast<TypeA>(ATAoct[i][j]);
//         }
//       }
//
//       _f.resize(_L);
//       _weight.resize(_gn);
//
//       // B0 depends only on the element, quadrature rule, polynomial order,
//       // and tetrahedron orientation.  Tetrahedra need both orientations;
//       // every other element needs one matrix.
//       const unsigned numberOfB0ATAMatrices = (_geomElemType == TET) ? 2u : 1u;
//       _B0ATA.resize(numberOfB0ATAMatrices);
//
//       std::vector<double> bo(_L);
//       std::vector<double> x(_dim);
//       Eigen::Matrix<TypeA, Eigen::Dynamic, Eigen::Dynamic> B0(_gn, _L);
//
//       for (unsigned iMatrix = 0; iMatrix < numberOfB0ATAMatrices; iMatrix++) {
//         const unsigned TETtype = iMatrix + 1u;
//
//         for (unsigned ig = 0; ig < _gn; ig++) {
//
//           if (_geomElemType == LINE || _geomElemType == QUAD || _geomElemType == HEX) {
//             for (unsigned k = 0; k < _dim; k++) {
//               x[k] = 0.5 * (1. + _xgp[k][ig]);
//             }
//           }
//           else if (_geomElemType == TRI) {
//             x[0] = 1. - _xgp[0][ig];
//             x[1] = _xgp[1][ig];
//           }
//           else if (_geomElemType == WEDGE) {
//             x[0] = 1. - _xgp[0][ig];
//             x[1] = _xgp[1][ig];
//             x[2] = 0.5 * (1. + _xgp[2][ig]);
//           }
//           else { // TET
//             x[0] = _xgp[0][ig] + _xgp[1][ig] + _xgp[2][ig];
//             if (TETtype == 1) {
//               x[1] = _xgp[1][ig] + _xgp[2][ig];
//               x[2] = _xgp[2][ig];
//             }
//             else { // TETtype == 2
//               x[1] = _xgp[2][ig] + _xgp[0][ig];
//               x[2] = _xgp[0][ig];
//             }
//           }
//
//           PolyBasis(x, bo);
//
//           for (unsigned i = 0; i < _L; i++) {
//             B0(ig, i) = static_cast<TypeA>(bo[i]);
//           }
//         }
//
//         _B0ATA[iMatrix].resize(_gn, _L);
//         _B0ATA[iMatrix].noalias() = B0 * ATA;
//       }
//     }
//
//     ~CutFemWeight() {
//       clear();
//       delete _gauss;
//       delete _obj;
//     }
//
//     void clear() {
//       if (_obj != NULL) _obj->clear();
//     }
//
//     void operator()(const int &s, const std::vector<TypeIO> &a, const TypeIO &d,
//                     std::vector<TypeIO> &weightCF);
//
//     void GetWeightWithMap(const int &s, const std::vector<TypeIO> &a, const TypeIO &d,
//                           std::vector<TypeIO> &weightCF);
//
//     void UpdateQuadratureRule(const unsigned &qM) {
//       clear();
//       delete _gauss;
//       delete _obj;
//       _gauss = NULL;
//       _obj = NULL;
//       ClearMap();
//
//       _qM = qM;
//       build();
//     }
//
//     unsigned GetCutFEMQuadratureOrder() {
//       return _qM;
//     }
//
//     unsigned GetGaussQuadratureOrder() {
//       return _gaussOrder;
//     }
//
//     unsigned GetDimension() {
//       return _dim;
//     }
//
//     unsigned GetGaussQuadraturePointNumber() {
//       return _gn;
//     }
//
//     const double* GetGaussWeightPointer() {
//       return _gauss->GetGaussWeightsPointer();
//     }
//
//     const double* GetGaussCoordinatePointer(const unsigned &k) {
//       return _gauss->GetGaussCoordinatePointer(k);
//     }
//
//     const unsigned& GetCounter() {
//       return _cntCall;
//     }
//
//     void ClearMap() {
//       for (unsigned i = 0; i < _WeightMap.size(); i++) _WeightMap[i].clear();
//       _WeightMap.clear();
//     }
//
//     void SetTetBaseType(const unsigned &value) {
//       _tetBaseType = value;
//     }
//
//   protected:
//     void PolyBasis(const std::vector<double> &x, std::vector<double> &bo);
//
//   private:
//     unsigned _dim;
//     GeomElType _geomElemType;
//     Gauss *_gauss;
//
//     std::vector<Eigen::Matrix<TypeA, Eigen::Dynamic, Eigen::Dynamic> > _B0ATA;
//     Eigen::Matrix<TypeA, Eigen::Dynamic, 1> _f;
//     Eigen::Matrix<TypeA, Eigen::Dynamic, 1> _weight;
//     unsigned _L;
//
//     unsigned _qM;
//     std::string _gaussType;
//     unsigned _gaussOrder;
//     unsigned _gn;
//     const double *_xgp[3];
//
//     unsigned _cntCall;
//
//     CutFEMmap<TypeA> *_obj;
//
//     unsigned _tetBaseType;
//
//     std::vector<std::map<std::pair<std::vector<float>, float>, std::vector<TypeIO> > > _WeightMap;
//     typename std::map<std::pair<std::vector<float>, float>, std::vector<TypeIO> >::iterator _it;
//     std::pair<std::vector<float>, float> _key;
//
//     std::vector<std::vector<TypeIO> > _bOld2d;
//     std::vector<std::vector<TypeIO> > _bNew2d;
//
//     std::vector<std::vector<std::vector<TypeIO> > > _bOld3d;
//     std::vector<std::vector<std::vector<TypeIO> > > _bNew3d;
// };
//
// template <class TypeIO, class TypeA>
// void CutFemWeight<TypeIO, TypeA>::GetWeightWithMap(const int &s,
//                                                    const std::vector<TypeIO> &a,
//                                                    const TypeIO &d,
//                                                    std::vector<TypeIO> &weightCF) {
//   _WeightMap.resize(s + 2);
//   std::vector<float> af(a.begin(), a.end());
//   const float df = static_cast<float>(d);
//   _key = std::make_pair(af, df);
//   _it = _WeightMap[s + 1].find(_key);
//
//   if (_it != _WeightMap[s + 1].end()) {
//     weightCF = _it->second;
//     return;
//   }
//
//   operator()(s, a, d, weightCF);
//   _WeightMap[s + 1][_key] = weightCF;
// }
//
// template <class TypeIO, class TypeA>
// void CutFemWeight<TypeIO, TypeA>::operator()(const int &s,
//                                              const std::vector<TypeIO> &a,
//                                              const TypeIO &d,
//                                              std::vector<TypeIO> &weightCF) {
//
//   _cntCall++;
//
//   std::vector<TypeA> aA(_dim);
//   for (unsigned i = 0; i < _dim; i++) {
//     aA[i] = static_cast<TypeA>(a[i]);
//   }
//   const TypeA dA = static_cast<TypeA>(d);
//
//   if (_obj != NULL) _obj->SetBaseType(_tetBaseType);
//
//   unsigned count = 0;
//
//   if (_dim == 3) {
//     for (unsigned q = 0; q <= _qM; q++) {
//       for (int ii = q; ii >= 0; ii--) {
//         for (int jj = q - ii; jj >= 0; jj--) {
//           const unsigned i = static_cast<unsigned>(ii);
//           const unsigned j = static_cast<unsigned>(jj);
//           const unsigned k = q - i - j;
//           _f[count] = (_geomElemType == WEDGE)
//                       ? Prism<TypeA, TypeA>(s, {i, j, k}, aA, dA)
//                       : (*_obj)(s, {i, j, k}, aA, dA);
//           count++;
//         }
//       }
//     }
//   }
//   else if (_dim == 2) {
//     for (unsigned q = 0; q <= _qM; q++) {
//       for (unsigned j = 0; j <= q; j++) {
//         const unsigned i = q - j;
//         _f[count] = (*_obj)(s, {i, j}, aA, dA);
//         count++;
//       }
//     }
//   }
//   else if (_dim == 1) {
//     for (unsigned i = 0; i <= _qM; i++) {
//       _f[i] = (*_obj)(s, {i}, aA, dA);
//     }
//   }
//   else {
//     std::cout << " Dimension =" << _dim << " not admissible" << std::endl;
//     abort();
//   }
//
//   unsigned B0ATAindex = 0;
//
//   if (_geomElemType == TET) {
//     unsigned TETtype = _tetBaseType;
//
//     if (TETtype == 0) {
//       TETtype = (std::max(std::fabs(a[1] - a[0]), std::fabs(a[2] - a[1]))
//                  >= std::fabs(a[0] - a[2])) ? 1u : 2u;
//     }
//
//     if (TETtype != 1 && TETtype != 2) {
//       std::cout << " TET base type = " << TETtype << " not admissible" << std::endl;
//       abort();
//     }
//
//     B0ATAindex = TETtype - 1u;
//   }
//
//   _weight.noalias() = _B0ATA[B0ATAindex] * _f;
//
//   weightCF.resize(_gn);
//   for (unsigned ig = 0; ig < _gn; ig++) {
//     weightCF[ig] = static_cast<TypeIO>(_weight[ig]);
//   }
// }
//
// template <class TypeIO, class TypeA>
// void CutFemWeight<TypeIO, TypeA>::PolyBasis(const std::vector<double> &x,
//                                             std::vector<double> &bo) {
//   unsigned count = 0;
//
//   if (_dim == 3) {
//     _bOld3d.resize(_qM + 1);
//     _bNew3d.resize(_qM + 1);
//     for (unsigned q = 0; q <= _qM; q++) {
//       _bOld3d[q].assign(_qM + 1, std::vector<TypeIO>(_qM + 1));
//       _bNew3d[q].assign(_qM + 1, std::vector<TypeIO>(_qM + 1));
//     }
//
//     bo[0] = _bNew3d[0][0][0] = 1.;
//     for (unsigned q = 1; q <= _qM; q++) {
//       _bOld3d.swap(_bNew3d);
//       for (int ii = q; ii >= 0; ii--) {
//         for (int jj = q - ii; jj >= 0; jj--) {
//           count++;
//           const unsigned i = static_cast<unsigned>(ii);
//           const unsigned j = static_cast<unsigned>(jj);
//           const unsigned k = q - i - j;
//
//           if (i > j) {
//             if (i > k) {
//               bo[count] = _bNew3d[i][j][k] = _bOld3d[i - 1][j][k] * x[0];
//             }
//             else {
//               bo[count] = _bNew3d[i][j][k] = _bOld3d[i][j][k - 1] * x[2];
//             }
//           }
//           else {
//             if (j > k) {
//               bo[count] = _bNew3d[i][j][k] = _bOld3d[i][j - 1][k] * x[1];
//             }
//             else {
//               bo[count] = _bNew3d[i][j][k] = _bOld3d[i][j][k - 1] * x[2];
//             }
//           }
//         }
//       }
//     }
//   }
//   else if (_dim == 2) {
//     _bOld2d.assign(_qM + 1, std::vector<TypeIO>(_qM + 1));
//     _bNew2d.assign(_qM + 1, std::vector<TypeIO>(_qM + 1));
//     bo[0] = _bNew2d[0][0] = 1.;
//
//     for (unsigned q = 1; q <= _qM; q++) {
//       _bOld2d.swap(_bNew2d);
//       for (unsigned j = 0; j <= q; j++) {
//         count++;
//         const unsigned i = q - j;
//         if (i > j) bo[count] = _bNew2d[i][j] = _bOld2d[i - 1][j] * x[0];
//         else bo[count] = _bNew2d[i][j] = _bOld2d[i][j - 1] * x[1];
//       }
//     }
//   }
//   else if (_dim == 1) {
//     bo[0] = 1.;
//     for (unsigned i = 1; i <= _qM; i++) {
//       bo[i] = bo[i - 1] * x[0];
//     }
//   }
//   else {
//     std::cout << " Dimension =" << _dim << " not admissible" << std::endl;
//     abort();
//   }
// }
//
// #endif
