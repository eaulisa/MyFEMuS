#ifndef __femus_cut_fem_int_hpp__
#define __femus_cut_fem_int_hpp__

#include "GeomElTypeEnum.hpp"

#include "CutLine.hpp"
#include "CutHyperCube.hpp"
#include "CutTriangle.hpp"
#include "CutTetrahedron.hpp"
#include "CutPrism.hpp"

#include "CutFem.hpp"
#include "FELevelSetMoments.hpp"
#include "Hex27LevelSetMoments.hpp"
#include "AdaptiveLineSliceMoments.hpp"
#include "FemusAlgoim.hpp"

#include "GramSchmidt.hpp"

#include "GaussPoints.hpp"


#include "SquareParabolaMoment.hpp"

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
    CutFemWeight(const GeomElType &geomElemType, const unsigned &qM, const std::string &gaussType) : _geomElemType(geomElemType), _gauss(NULL), _qM(qM), _gaussType(gaussType), _cntCall(0), _obj(NULL), _tetBaseType(0) {

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
      _nodalWeightCacheValid = false;

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

      Eigen::Matrix<TypeA,Eigen::Dynamic,Eigen::Dynamic> gramInverse(ATAoct.size(),ATAoct.size());
      for(unsigned i = 0; i < ATAoct.size(); i++) {
        for(unsigned j = 0; j < ATAoct[i].size(); j++) {
          gramInverse(i,j) = static_cast<TypeA>(ATAoct[i][j]);
        }
      }

      _f.resize(_L);
      _fPair.resize(_L, 2);

      // Cache the combined moment-to-Gauss map in TypeA for each orientation.
      const unsigned numberOfB0Matrices = (_geomElemType == TET) ? 2u : 1u;
      _momentToWeight.resize(numberOfB0Matrices);

      std::vector<TypeA> x(_dim);
      std::vector<std::vector<TypeA>> powers(_dim,std::vector<TypeA>(_qM+1U,TypeA(1)));

      for(unsigned iMatrix = 0; iMatrix < numberOfB0Matrices; iMatrix++) {
        const unsigned TETtype = iMatrix + 1u;

        Eigen::Matrix<TypeA,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> basis(_gn,_L);

        for(unsigned ig = 0; ig < _gn; ig++) {

          if(_geomElemType == LINE || _geomElemType == QUAD || _geomElemType == HEX) {
            for(unsigned k = 0; k < _dim; k++) {
              x[k] = (TypeA(1)+TypeA(_xgp[k][ig]))/TypeA(2);
            }
          }
          else if(_geomElemType == TRI) {
            x[0] = TypeA(1)-TypeA(_xgp[0][ig]);
            x[1] = TypeA(_xgp[1][ig]);
          }
          else if(_geomElemType == WEDGE) {
            x[0] = TypeA(1)-TypeA(_xgp[0][ig]);
            x[1] = TypeA(_xgp[1][ig]);
            x[2] = (TypeA(1)+TypeA(_xgp[2][ig]))/TypeA(2);
          }
          else { // TET
            x[0] = TypeA(_xgp[0][ig])+TypeA(_xgp[1][ig])+TypeA(_xgp[2][ig]);
            if(TETtype == 1) {
              x[1] = TypeA(_xgp[1][ig])+TypeA(_xgp[2][ig]);
              x[2] = TypeA(_xgp[2][ig]);
            }
            else { // TETtype == 2
              x[1] = TypeA(_xgp[2][ig])+TypeA(_xgp[0][ig]);
              x[2] = TypeA(_xgp[0][ig]);
            }
          }

          for(unsigned d = 0; d < _dim; ++d) {
            for(unsigned i = 1; i <= _qM; ++i) powers[d][i] = powers[d][i-1U]*x[d];
          }
          unsigned index = 0;
          for(unsigned q = 0; q <= _qM; ++q) {
            if(_dim == 1) basis(ig,index++) = powers[0][q];
            else {
              for(int ii = q; ii >= 0; --ii) {
                unsigned i = static_cast<unsigned>(ii);
                if(_dim == 2) basis(ig,index++) = powers[0][i]*powers[1][q-i];
                else {
                  for(int jj = q-i; jj >= 0; --jj) {
                    unsigned j = static_cast<unsigned>(jj);
                    basis(ig,index++) = powers[0][i]*powers[1][j]*powers[2][q-i-j];
                  }
                }
              }
            }
          }
        }
        _momentToWeight[iMatrix].noalias() = basis*gramInverse;
      }
    }

    ~CutFemWeight() {
      clear();
      delete _gauss;
      delete _obj;
    }

    void clear() {
      if(_obj != NULL) _obj->clear();
      _nodalMomentCacheValid = false;
      _nodalWeightCacheValid = false;
    }

    void operator()(const int &s, const std::vector<TypeIO> &a, const TypeIO &d, std::vector<TypeIO> &weightCF);


    void operator()(const int &s, const std::vector<TypeA> &a, const TypeA &d, const bool &xSquared, std::vector<TypeIO> &weightCF);
    void operator()(const int &s, const std::vector<TypeA> &aX, const TypeA &dX, const std::vector<TypeA> &aY, const TypeA &dY, const bool &xSquared, std::vector<TypeIO> &weightCF);

    void GetWeightWithMap(const int &s, const std::vector<TypeIO> &a, const TypeIO &d, std::vector<TypeIO> &weightCF);

    void GetWeightFromMoments(const std::vector<std::vector<std::vector<double>>>& moments, std::vector<TypeIO>& weightCF) {
      if(_geomElemType != HEX || moments.size() < _qM+1U) throw std::invalid_argument("Expected unit-cube HEX moments.");
      unsigned count = 0;
      for(unsigned q = 0; q <= _qM; ++q) {
        for(int ii = q; ii >= 0; --ii) {
          for(int jj = q-static_cast<unsigned>(ii); jj >= 0; --jj) {
            unsigned i = static_cast<unsigned>(ii), j = static_cast<unsigned>(jj), k = q-i-j;
            if(moments[i].size() <= j || moments[i][j].size() <= k || !std::isfinite(moments[i][j][k])) throw std::invalid_argument("Invalid HEX moment table.");
            _f[count++] = static_cast<TypeA>(moments[i][j][k]);
          }
        }
      }
      ++_cntCall;
      EvaluateMomentWeights(weightCF);
    }

    void GetWeightFromMoments(const std::vector<std::vector<double>>& moments, std::vector<TypeIO>& weightCF) {
      if (_dim != 2 || (_geomElemType != QUAD && _geomElemType != TRI) || moments.size() < _qM + 1) {
        throw std::invalid_argument("Expected canonical QUAD or TRI moments.");
      }
      unsigned count = 0;
      for (unsigned q = 0; q <= _qM; ++q) {
        for (unsigned j = 0; j <= q; ++j) {
          unsigned i = q - j;
          if (moments[i].size() <= j || !std::isfinite(moments[i][j])) {
            throw std::invalid_argument("Invalid canonical moment table.");
          }
          _f[count++] = static_cast<TypeA>(moments[i][j]);
        }
      }
      ++_cntCall;
      EvaluateMomentWeights(weightCF);
    }

    void GetWeightFromMoments(const std::vector<double>& moments, std::vector<TypeIO>& weightCF) {
      if (_geomElemType != LINE || moments.size() < _qM + 1) {
        throw std::invalid_argument("Expected canonical LINE moments.");
      }
      for (unsigned i = 0; i <= _qM; ++i) {
        if (!std::isfinite(moments[i])) {
          throw std::invalid_argument("Invalid canonical LINE moment.");
        }
        _f[i] = static_cast<TypeA>(moments[i]);
      }
      ++_cntCall;
      EvaluateMomentWeights(weightCF);
    }

    void GetWeightFromFE(const int& s, const std::vector<std::vector<double>>& coefficient, std::vector<TypeIO>& weightCF, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      if (_geomElemType != QUAD || (s != 0 && s != -1)) {
        throw std::invalid_argument("QUAD9 FE weights require QUAD and s=0 or s=-1.");
      }
      if (coefficient.size() != 3) {
        throw std::invalid_argument("Expected a 3 by 3 QUAD9 coefficient matrix.");
      }
      bool planar = true;
      for (unsigned i = 0; i < 3; ++i) {
        if (coefficient[i].size() != 3) {
          throw std::invalid_argument("Expected a 3 by 3 QUAD9 coefficient matrix.");
        }
        for (unsigned j = 0; j < 3; ++j) {
          if (!std::isfinite(coefficient[i][j])) {
            throw std::invalid_argument("A QUAD9 coefficient is not finite.");
          }
          if (i + j > 1 && coefficient[i][j] != 0.0) {
            planar = false;
          }
        }
      }
      if (planar && (coefficient[1][0] != 0.0 || coefficient[0][1] != 0.0)) {
        std::vector<TypeIO> normal = {static_cast<TypeIO>(coefficient[1][0]), static_cast<TypeIO>(coefficient[0][1])};
        operator()(s, normal, static_cast<TypeIO>(coefficient[0][0]), weightCF);
        return;
      }
      std::vector<std::vector<double>> positive;
      std::vector<std::vector<double>> delta;
      FELevelSetMoments::IntegrateQuad9Canonical(coefficient, _qM, positive, delta, absoluteTolerance, relativeTolerance, -1.0, 1.0, -1.0, 1.0, s == -1);
      GetWeightFromMoments((s == 0) ? positive : delta, weightCF);
    }

    void GetWeightFromFE(const int& s, const std::vector<double>& coefficient, std::vector<TypeIO>& weightCF, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      if (s != 0 && s != -1) {
        throw std::invalid_argument("FE weights require s=0 or s=-1.");
      }
      if (_geomElemType == TRI) {
        if (coefficient.size() != 7) {
          throw std::invalid_argument("Expected seven TRI7 coefficients.");
        }
        for (unsigned k = 0; k < 7; ++k) {
          if (!std::isfinite(coefficient[k])) {
            throw std::invalid_argument("A TRI7 coefficient is not finite.");
          }
        }
        if (coefficient[3] == 0.0 && coefficient[4] == 0.0 && coefficient[5] == 0.0 && coefficient[6] == 0.0 && (coefficient[1] != 0.0 || coefficient[2] != 0.0)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(coefficient[1]), static_cast<TypeIO>(coefficient[2])};
          operator()(s, normal, static_cast<TypeIO>(coefficient[0]), weightCF);
          return;
        }
        std::vector<std::vector<double>> positive;
        std::vector<std::vector<double>> delta;
        FELevelSetMoments::IntegrateTri7Canonical(coefficient, _qM, positive, delta, absoluteTolerance, relativeTolerance, s == -1);
        GetWeightFromMoments((s == 0) ? positive : delta, weightCF);
      }
      else if (_geomElemType == LINE && coefficient.size() == 3) {
        for (unsigned k = 0; k < 3; ++k) {
          if (!std::isfinite(coefficient[k])) {
            throw std::invalid_argument("A LINE3 coefficient is not finite.");
          }
        }
        if (coefficient[0] == 0.0 && coefficient[1] != 0.0) {
          std::vector<TypeIO> normal(1, static_cast<TypeIO>(coefficient[1]));
          operator()(s, normal, static_cast<TypeIO>(coefficient[2]), weightCF);
          return;
        }
        std::vector<double> positive;
        std::vector<double> delta;
        FELevelSetMoments::Integrate1DCanonical(coefficient[0], coefficient[1], coefficient[2], _qM, positive, delta, -1.0, 1.0, s == -1);
        GetWeightFromMoments((s == 0) ? positive : delta, weightCF);
      }
      else {
        throw std::invalid_argument("Expected TRI7 coefficients or LINE3 coefficients {a,b,c}.");
      }
    }

    void GetWeightFromCanonicalFE(const int& s, const std::vector<std::vector<double>>& coefficient, std::vector<TypeIO>& weightCF, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      if (_geomElemType != QUAD || (s != 0 && s != -1)) {
        throw std::invalid_argument("Canonical QUAD9 weights require QUAD and s=0 or s=-1.");
      }
      if (coefficient.size() != 3) {
        throw std::invalid_argument("Expected a 3 by 3 unit-square QUAD9 coefficient matrix.");
      }
      bool planar = true;
      for (unsigned i = 0; i < 3; ++i) {
        if (coefficient[i].size() != 3) {
          throw std::invalid_argument("Expected a 3 by 3 unit-square QUAD9 coefficient matrix.");
        }
        for (unsigned j = 0; j < 3; ++j) {
          if (!std::isfinite(coefficient[i][j])) {
            throw std::invalid_argument("A QUAD9 coefficient is not finite.");
          }
          if (i + j > 1 && coefficient[i][j] != 0.0) {
            planar = false;
          }
        }
      }
      if (planar && (coefficient[1][0] != 0.0 || coefficient[0][1] != 0.0)) {
        std::vector<TypeIO> normal = {static_cast<TypeIO>(0.5L * coefficient[1][0]), static_cast<TypeIO>(0.5L * coefficient[0][1])};
        operator()(s, normal, static_cast<TypeIO>(static_cast<long double>(coefficient[0][0]) + 0.5L * coefficient[1][0] + 0.5L * coefficient[0][1]), weightCF);
        return;
      }
      std::vector<std::vector<double>> positive;
      std::vector<std::vector<double>> delta;
      FELevelSetMoments::IntegrateQuad9Unit(coefficient, _qM, positive, delta, absoluteTolerance, relativeTolerance, s == -1);
      GetWeightFromMoments((s == 0) ? positive : delta, weightCF);
    }

    void GetWeightFromCanonicalFE(const int& s, const std::vector<double>& coefficient, std::vector<TypeIO>& weightCF, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      if (s != 0 && s != -1) {
        throw std::invalid_argument("Canonical FE weights require s=0 or s=-1.");
      }
      if (_geomElemType == TRI) {
        if (coefficient.size() != 7) {
          throw std::invalid_argument("Expected reflected TRI7 coefficients of 1,x,y,x*x,x*y,y*y,x*y*(x-y).");
        }
        for (unsigned k = 0; k < 7; ++k) {
          if (!std::isfinite(coefficient[k])) {
            throw std::invalid_argument("A TRI7 coefficient is not finite.");
          }
        }
        if (coefficient[3] == 0.0 && coefficient[4] == 0.0 && coefficient[5] == 0.0 && coefficient[6] == 0.0 && (coefficient[1] != 0.0 || coefficient[2] != 0.0)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(-coefficient[1]), static_cast<TypeIO>(coefficient[2])};
          operator()(s, normal, static_cast<TypeIO>(static_cast<long double>(coefficient[0]) + coefficient[1]), weightCF);
          return;
        }
        std::vector<std::vector<double>> positive;
        std::vector<std::vector<double>> delta;
        FELevelSetMoments::IntegrateTri7Reflected(coefficient, _qM, positive, delta, absoluteTolerance, relativeTolerance, s == -1);
        GetWeightFromMoments((s == 0) ? positive : delta, weightCF);
      }
      else if (_geomElemType == LINE && coefficient.size() == 3) {
        for (unsigned k = 0; k < 3; ++k) {
          if (!std::isfinite(coefficient[k])) {
            throw std::invalid_argument("A LINE3 coefficient is not finite.");
          }
        }
        if (coefficient[0] == 0.0 && coefficient[1] != 0.0) {
          std::vector<TypeIO> normal(1, static_cast<TypeIO>(0.5L * coefficient[1]));
          operator()(s, normal, static_cast<TypeIO>(static_cast<long double>(coefficient[2]) + 0.5L * coefficient[1]), weightCF);
          return;
        }
        std::vector<double> positive;
        std::vector<double> delta;
        FELevelSetMoments::Integrate1D(coefficient[0], coefficient[1], coefficient[2], _qM, positive, delta, 0.0, 1.0, s == -1);
        GetWeightFromMoments((s == 0) ? positive : delta, weightCF);
      }
      else {
        throw std::invalid_argument("Expected reflected TRI7 coefficients or unit-interval LINE3 coefficients {a,b,c}.");
      }
    }

    void GetWeightFromCanonicalFE(const int& s, const std::vector<std::vector<std::vector<double>>>& coefficient, std::vector<TypeIO>& weightCF, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      if(_geomElemType != HEX || (s != 0 && s != -1)) throw std::invalid_argument("HEX27 weights require HEX and s=0 or s=-1.");
      std::vector<long double> checked;
      long double scale;
      FELevelSetMoments::Hex27Normalize(coefficient,checked,scale);
      bool planar = true;
      for(unsigned i = 0; i < 3; ++i) {
        for(unsigned j = 0; j < 3; ++j) {
          for(unsigned k = 0; k < 3; ++k) if(i+j+k > 1 && coefficient[i][j][k] != 0.) planar = false;
        }
      }
      if(planar && (coefficient[1][0][0] != 0. || coefficient[0][1][0] != 0. || coefficient[0][0][1] != 0.)) {
        std::vector<TypeIO> normal = {static_cast<TypeIO>(0.5L*coefficient[1][0][0]),static_cast<TypeIO>(0.5L*coefficient[0][1][0]),static_cast<TypeIO>(0.5L*coefficient[0][0][1])};
        operator()(s,normal,static_cast<TypeIO>(static_cast<long double>(coefficient[0][0][0])+0.5L*coefficient[1][0][0]+0.5L*coefficient[0][1][0]+0.5L*coefficient[0][0][1]),weightCF);
        return;
      }
      std::vector<std::vector<std::vector<double>>> positive, delta;
      FELevelSetMoments::IntegrateHex27Unit(coefficient,_qM,positive,delta,absoluteTolerance,relativeTolerance,s == -1);
      GetWeightFromMoments(s == 0 ? positive : delta,weightCF);
    }

    void GetWeightFromNodalValues(const int& s, const std::vector<double>& nodalValues, std::vector<TypeIO>& weightCF, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      if(_geomElemType == HEX) {
        std::vector<std::vector<std::vector<double>>> coefficient;
        FELevelSetMoments::Hex27UnitCoefficientsFromNodalValues(nodalValues,coefficient);
        GetWeightFromCanonicalFE(s,coefficient,weightCF,absoluteTolerance,relativeTolerance);
      }
      else if (_geomElemType == QUAD) {
        std::vector<std::vector<double>> coefficient;
        FELevelSetMoments::Quad9UnitCoefficientsFromNodalValues(nodalValues, coefficient);
        GetWeightFromCanonicalFE(s, coefficient, weightCF, absoluteTolerance, relativeTolerance);
      }
      else if (_geomElemType == TRI) {
        std::vector<double> coefficient;
        FELevelSetMoments::Tri7ReflectedCoefficientsFromNodalValues(nodalValues, coefficient);
        GetWeightFromCanonicalFE(s, coefficient, weightCF, absoluteTolerance, relativeTolerance);
      }
      else if (_geomElemType == LINE) {
        double a = 0.0;
        double b = 0.0;
        double c = 0.0;
        FELevelSetMoments::Line3CoefficientsFromNodalValues(nodalValues, a, b, c);
        std::vector<double> coefficient = {a, b, c};
        GetWeightFromCanonicalFE(s, coefficient, weightCF, absoluteTolerance, relativeTolerance);
      }
      else {
        throw std::invalid_argument("FE nodal weights support LINE3, QUAD9, TRI7, and HEX27.");
      }
    }

    void GetWeightsFromNodalValues(const std::vector<double>& nodalValues, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      GetNodalWeightPair(nodalValues,positiveWeight,deltaWeight,absoluteTolerance,relativeTolerance,false);
    }

    void GetWeightsFromNodalValuesByLineSlices(const std::vector<double>& nodalValues, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      GetNodalWeightPair(nodalValues,positiveWeight,deltaWeight,absoluteTolerance,relativeTolerance,true);
    }

    void GetWeightsFromNodalValuesByAlgoim(const std::vector<double>& nodalValues, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight, unsigned quadratureOrder = 20, unsigned subdivisions = 0) {
      if(quadratureOrder == 0 || quadratureOrder > 100 || subdivisions > 32) {
        throw std::invalid_argument("Invalid Algoim quadrature order or subdivisions.");
      }
      GetNodalWeightPair(nodalValues,positiveWeight,deltaWeight,1.0e-12,1.0e-12,false,quadratureOrder,subdivisions,&FemusAlgoim::IntegrateNodalMoments);
    }

    void UpdateQuadratureRule(const unsigned &qM) {
      bool retainNodalMoments = _nodalMomentCacheValid;
      clear();
      delete _gauss;
      delete _obj;
      _gauss = NULL;
      _obj = NULL;
      ClearMap();

      _qM = qM;
      build();
      _nodalMomentCacheValid = retainNodalMoments;
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
      _nodalMomentCacheValid = false;
      _nodalWeightCacheValid = false;
    }

    void SetTetBaseType(const unsigned &value) {
      _tetBaseType = value;
    }

  protected:
    void PolyBasis(const std::vector<double> &x, std::vector<double> &bo);

  private:
    void GetNodalWeightPair(const std::vector<double>& nodalValues, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight, double absoluteTolerance, double relativeTolerance, bool useLineSlices, unsigned algoimOrder = 0, unsigned algoimSubdivisions = 1, void (*algoimIntegrator)(const std::vector<double>&, unsigned, std::vector<double>&, std::vector<double>&, unsigned, unsigned) = nullptr) {
      if(&positiveWeight == &deltaWeight) {
        throw std::invalid_argument("Heaviside and delta weights require distinct output vectors.");
      }
      if(!(absoluteTolerance > 0.) || !(relativeTolerance >= 0.) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance)) {
        throw std::invalid_argument("Invalid nodal integration tolerance.");
      }
      useLineSlices = useLineSlices && (_geomElemType == LINE || _geomElemType == QUAD || _geomElemType == TRI);
      bool sameField = algoimOrder == _nodalMomentAlgoimOrder && algoimSubdivisions == _nodalMomentAlgoimSubdivisions && useLineSlices == _nodalMomentUsesLineSlices && _nodalMomentCacheValid && nodalValues == _nodalValues && absoluteTolerance == _nodalAbsoluteTolerance && relativeTolerance == _nodalRelativeTolerance;
      if(sameField && _nodalWeightCacheValid) {
        positiveWeight = _nodalPositiveWeight;
        deltaWeight = _nodalDeltaWeight;
        return;
      }
      _nodalWeightCacheValid = false;
      if(sameField && _qM <= _nodalMomentDegree) {
        _fPair = _nodalMomentCache.topRows(_L);
      }
      else {
        IntegrateNodalMomentPair(nodalValues, absoluteTolerance, relativeTolerance, useLineSlices, algoimOrder, algoimSubdivisions, algoimIntegrator);
        _nodalMomentCache = _fPair;
        _nodalValues = nodalValues;
        _nodalMomentDegree = _qM;
        _nodalAbsoluteTolerance = absoluteTolerance;
        _nodalRelativeTolerance = relativeTolerance;
        _nodalMomentUsesLineSlices = useLineSlices;
        _nodalMomentAlgoimOrder = algoimOrder;
        _nodalMomentAlgoimSubdivisions = algoimSubdivisions;
        _nodalMomentCacheValid = true;
      }
      ++_cntCall;
      EvaluateMomentWeightPair(positiveWeight, deltaWeight);
      _nodalPositiveWeight = positiveWeight;
      _nodalDeltaWeight = deltaWeight;
      _nodalWeightCacheValid = true;
    }

    void StoreAlgoimMomentPair(const std::vector<double>& positive, const std::vector<double>& delta) {
      if(positive.size() != _L || delta.size() != _L) {
        throw std::runtime_error("Algoim returned an inconsistent moment table.");
      }
      for(unsigned i = 0; i < _L; ++i) {
        _fPair(i,0) = static_cast<TypeA>(positive[i]);
        _fPair(i,1) = static_cast<TypeA>(delta[i]);
      }
    }

    void StoreNodalMomentPair(const std::vector<std::vector<std::vector<double>>>& positive, const std::vector<std::vector<std::vector<double>>>& delta) {
      unsigned count = 0;
      for(unsigned q = 0; q <= _qM; ++q) {
        for(int ii = q; ii >= 0; --ii) {
          for(int jj = q-static_cast<unsigned>(ii); jj >= 0; --jj) {
            unsigned i = static_cast<unsigned>(ii), j = static_cast<unsigned>(jj), k = q-i-j;
            _fPair(count,0) = static_cast<TypeA>(positive[i][j][k]);
            _fPair(count,1) = static_cast<TypeA>(delta[i][j][k]);
            ++count;
          }
        }
      }
    }

    void StoreNodalMomentPair(const std::vector<std::vector<double>>& positive, const std::vector<std::vector<double>>& delta) {
      unsigned count = 0;
      for(unsigned q = 0; q <= _qM; ++q) {
        for(unsigned j = 0; j <= q; ++j) {
          unsigned i = q - j;
          _fPair(count, 0) = static_cast<TypeA>(positive[i][j]);
          _fPair(count, 1) = static_cast<TypeA>(delta[i][j]);
          ++count;
        }
      }
    }

    void IntegratePlanarMomentPair(const std::vector<TypeIO>& normal, const TypeIO& offset) {
      std::vector<TypeA> a(_dim);
      for(unsigned d = 0; d < _dim; ++d) a[d] = static_cast<TypeA>(normal[d]);
      TypeA constant = static_cast<TypeA>(offset);
      if(_dim == 3) {
        unsigned count = 0;
        for(unsigned q = 0; q <= _qM; ++q) {
          for(int ii = q; ii >= 0; --ii) {
            for(int jj = q-static_cast<unsigned>(ii); jj >= 0; --jj) {
              unsigned i = static_cast<unsigned>(ii), j = static_cast<unsigned>(jj);
              std::vector<unsigned> powers = {i,j,q-i-j};
              _fPair(count,0) = (*_obj)(0,powers,a,constant);
              _fPair(count,1) = (*_obj)(-1,powers,a,constant);
              ++count;
            }
          }
        }
      }
      else if(_dim == 2) {
        unsigned count = 0;
        for(unsigned q = 0; q <= _qM; ++q) {
          for(unsigned j = 0; j <= q; ++j) {
            std::vector<unsigned> powers = {q - j, j};
            _fPair(count, 0) = (*_obj)(0, powers, a, constant);
            _fPair(count, 1) = (*_obj)(-1, powers, a, constant);
            ++count;
          }
        }
      }
      else {
        for(unsigned i = 0; i <= _qM; ++i) {
          std::vector<unsigned> powers(1, i);
          _fPair(i, 0) = (*_obj)(0, powers, a, constant);
          _fPair(i, 1) = (*_obj)(-1, powers, a, constant);
        }
      }
    }

    void IntegrateNodalMomentPair(const std::vector<double>& nodalValues, double absoluteTolerance, double relativeTolerance, bool useLineSlices, unsigned algoimOrder = 0, unsigned algoimSubdivisions = 1, void (*algoimIntegrator)(const std::vector<double>&, unsigned, std::vector<double>&, std::vector<double>&, unsigned, unsigned) = nullptr) {
      if(_geomElemType == HEX) {
        std::vector<std::vector<std::vector<double>>> coefficient;
        FELevelSetMoments::Hex27UnitCoefficientsFromNodalValues(nodalValues,coefficient);
        bool planar = true;
        for(unsigned i = 0; i < 3; ++i) {
          for(unsigned j = 0; j < 3; ++j) {
            for(unsigned k = 0; k < 3; ++k) if(i+j+k > 1 && coefficient[i][j][k] != 0.) planar = false;
          }
        }
        if(planar && (coefficient[1][0][0] != 0. || coefficient[0][1][0] != 0. || coefficient[0][0][1] != 0.)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(0.5L*coefficient[1][0][0]),static_cast<TypeIO>(0.5L*coefficient[0][1][0]),static_cast<TypeIO>(0.5L*coefficient[0][0][1])};
          IntegratePlanarMomentPair(normal,static_cast<TypeIO>(static_cast<long double>(coefficient[0][0][0])+0.5L*coefficient[1][0][0]+0.5L*coefficient[0][1][0]+0.5L*coefficient[0][0][1]));
          return;
        }
        if(algoimOrder != 0) {
          std::vector<double> positiveMoments;
          std::vector<double> deltaMoments;
          algoimIntegrator(nodalValues,_qM,positiveMoments,deltaMoments,algoimOrder,algoimSubdivisions);
          StoreAlgoimMomentPair(positiveMoments,deltaMoments);
          return;
        }
        std::vector<std::vector<std::vector<double>>> positive, delta;
        FELevelSetMoments::IntegrateHex27Unit(coefficient,_qM,positive,delta,absoluteTolerance,relativeTolerance);
        StoreNodalMomentPair(positive,delta);
      }
      else if(_geomElemType == QUAD) {
        std::vector<std::vector<double>> coefficient;
        FELevelSetMoments::Quad9UnitCoefficientsFromNodalValues(nodalValues, coefficient);
        bool planar = true;
        for(unsigned i = 0; i < 3; ++i) {
          for(unsigned j = 0; j < 3; ++j) {
            if(i + j > 1 && coefficient[i][j] != 0.) planar = false;
          }
        }
        if(planar && (coefficient[1][0] != 0. || coefficient[0][1] != 0.)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(0.5L * coefficient[1][0]), static_cast<TypeIO>(0.5L * coefficient[0][1])};
          IntegratePlanarMomentPair(normal, static_cast<TypeIO>(static_cast<long double>(coefficient[0][0]) + 0.5L * coefficient[1][0] + 0.5L * coefficient[0][1]));
          return;
        }
        if(algoimOrder != 0) {
          std::vector<double> positiveMoments;
          std::vector<double> deltaMoments;
          algoimIntegrator(nodalValues,_qM,positiveMoments,deltaMoments,algoimOrder,algoimSubdivisions);
          StoreAlgoimMomentPair(positiveMoments,deltaMoments);
          return;
        }
        std::vector<std::vector<double>> positive;
        std::vector<std::vector<double>> delta;
        if(useLineSlices) FELevelSetMoments::IntegrateQuad9NodalByLineSlices(nodalValues,_qM,positive,delta,absoluteTolerance,relativeTolerance);
        else FELevelSetMoments::IntegrateQuad9Unit(coefficient, _qM, positive, delta, absoluteTolerance, relativeTolerance);
        StoreNodalMomentPair(positive, delta);
      }
      else if(_geomElemType == TRI) {
        std::vector<double> coefficient;
        FELevelSetMoments::Tri7ReflectedCoefficientsFromNodalValues(nodalValues, coefficient);
        if(coefficient[3] == 0. && coefficient[4] == 0. && coefficient[5] == 0. && coefficient[6] == 0. && (coefficient[1] != 0. || coefficient[2] != 0.)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(-coefficient[1]), static_cast<TypeIO>(coefficient[2])};
          IntegratePlanarMomentPair(normal, static_cast<TypeIO>(static_cast<long double>(coefficient[0]) + coefficient[1]));
          return;
        }
        if(algoimOrder != 0) {
          std::vector<double> positiveMoments;
          std::vector<double> deltaMoments;
          algoimIntegrator(nodalValues,_qM,positiveMoments,deltaMoments,algoimOrder,algoimSubdivisions);
          StoreAlgoimMomentPair(positiveMoments,deltaMoments);
          return;
        }
        std::vector<std::vector<double>> positive;
        std::vector<std::vector<double>> delta;
        if(useLineSlices) FELevelSetMoments::IntegrateTri7NodalByLineSlices(nodalValues,_qM,positive,delta,absoluteTolerance,relativeTolerance);
        else FELevelSetMoments::IntegrateTri7Reflected(coefficient, _qM, positive, delta, absoluteTolerance, relativeTolerance);
        StoreNodalMomentPair(positive, delta);
      }
      else if(_geomElemType == LINE) {
        double a = 0.;
        double b = 0.;
        double c = 0.;
        FELevelSetMoments::Line3CoefficientsFromNodalValues(nodalValues, a, b, c);
        if(a == 0. && b != 0.) {
          std::vector<TypeIO> normal(1, static_cast<TypeIO>(0.5L * b));
          IntegratePlanarMomentPair(normal, static_cast<TypeIO>(static_cast<long double>(c) + 0.5L * b));
          return;
        }
        std::vector<double> positive;
        std::vector<double> delta;
        if(useLineSlices) FELevelSetMoments::IntegrateLine3NodalByLineSlices(nodalValues,_qM,positive,delta,absoluteTolerance,relativeTolerance);
        else FELevelSetMoments::Integrate1D(a, b, c, _qM, positive, delta);
        for(unsigned i = 0; i <= _qM; ++i) {
          _fPair(i, 0) = static_cast<TypeA>(positive[i]);
          _fPair(i, 1) = static_cast<TypeA>(delta[i]);
        }
      }
      else {
        throw std::invalid_argument("Paired nodal weights support LINE3, QUAD9, TRI7, and HEX27.");
      }
    }

    void EvaluateMomentWeightPair(std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight) {
      positiveWeight.resize(_gn);
      deltaWeight.resize(_gn);
      for(unsigned ig = 0; ig < _gn; ++ig) {
        const TypeA* basisRow = _momentToWeight[0].data()+ig*_L;
        TypeA positive = 0.;
        TypeA delta = 0.;
        for(unsigned i = 0; i < _L; ++i) {
          positive += basisRow[i]*_fPair(i,0);
          delta += basisRow[i]*_fPair(i,1);
        }
        positiveWeight[ig] = static_cast<TypeIO>(positive);
        deltaWeight[ig] = static_cast<TypeIO>(delta);
      }
    }

    void EvaluateMomentWeights(std::vector<TypeIO>& weightCF, unsigned basisIndex = 0) {
      weightCF.resize(_gn);
      for(unsigned ig = 0; ig < _gn; ++ig) {
        const TypeA* basisRow = _momentToWeight[basisIndex].data()+ig*_L;
        TypeA value = 0.;
        for(unsigned i = 0; i < _L; ++i) {
          value += basisRow[i]*_f[i];
        }
        weightCF[ig] = static_cast<TypeIO>(value);
      }
    }

    unsigned _dim;
    GeomElType _geomElemType;
    Gauss *_gauss;

    std::vector<Eigen::Matrix<TypeA, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> _momentToWeight;
    Eigen::Matrix<TypeA, Eigen::Dynamic, 1> _f;
    Eigen::Matrix<TypeA, Eigen::Dynamic, 2> _fPair;
    Eigen::Matrix<TypeA, Eigen::Dynamic, 2> _nodalMomentCache;
    bool _nodalMomentUsesLineSlices = false;
    unsigned _nodalMomentAlgoimOrder = 0;
    unsigned _nodalMomentAlgoimSubdivisions = 1;
    bool _nodalMomentCacheValid = false;
    bool _nodalWeightCacheValid = false;
    unsigned _nodalMomentDegree = 0;
    double _nodalAbsoluteTolerance = 0.;
    double _nodalRelativeTolerance = 0.;
    std::vector<double> _nodalValues;
    std::vector<TypeIO> _nodalPositiveWeight;
    std::vector<TypeIO> _nodalDeltaWeight;
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
void CutFemWeight<TypeIO, TypeA>::GetWeightWithMap(const int &s, const std::vector<TypeIO> &a, const TypeIO &d, std::vector<TypeIO> &weightCF) {
  if(_WeightMap.size() < static_cast<unsigned>(s+2)) _WeightMap.resize(s+2);
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
void CutFemWeight<TypeIO, TypeA>::operator()(const int &s, const std::vector<TypeIO> &a, const TypeIO &d, std::vector<TypeIO> &weightCF) {

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
          _f[count] = (_geomElemType == WEDGE) ? Prism<TypeA, TypeA>(s, {i, j, k}, aA, dA) : (*_obj)(s, {i, j, k}, aA, dA);
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

  unsigned B0index = 0;

  if(_geomElemType == TET) {
    unsigned TETtype = _tetBaseType;

    if(TETtype == 0) {
      TETtype = (std::max(std::fabs(a[1] - a[0]), std::fabs(a[2] - a[1])) >= std::fabs(a[0] - a[2])) ? 1u : 2u;
    }

    if(TETtype != 1 && TETtype != 2) {
      std::cout << " TET base type = " << TETtype << " not admissible" << std::endl;
      abort();
    }

    B0index = TETtype - 1u;
  }

  EvaluateMomentWeights(weightCF, B0index);
}



template <class TypeIO, class TypeA>
void CutFemWeight<TypeIO, TypeA>::operator()(const int &s, const std::vector<TypeA> &a, const TypeA &d, const bool &xSquared, std::vector<TypeIO> &weightCF) {
  if(_geomElemType != QUAD) throw std::invalid_argument("Quadratic CutFemWeight requires QUAD");
  if(a.size() != 3) throw std::invalid_argument("Quadratic CutFemWeight requires a={k,b,c}");
  if(s < -1) throw std::invalid_argument("Quadratic CutFemWeight requires s >= -1");

  _cntCall++;

  unsigned count = 0;
  for(unsigned q = 0; q <= _qM; q++) {
    for(unsigned j = 0; j <= q; j++) {
      unsigned i = q - j;
      _f[count] = SquareParabolaMoment<TypeA>(s, {i, j}, a, d, xSquared);
      count++;
    }
  }

  EvaluateMomentWeights(weightCF);
}

template <class TypeIO, class TypeA>
void CutFemWeight<TypeIO, TypeA>::operator()(const int &s, const std::vector<TypeA> &aX, const TypeA &dX, const std::vector<TypeA> &aY, const TypeA &dY, const bool &xSquared, std::vector<TypeIO> &weightCF) {
  if(xSquared) operator()(s, aX, dX, true, weightCF);
  else operator()(s, aY, dY, false, weightCF);
}

template <class TypeIO, class TypeA>
void CutFemWeight<TypeIO, TypeA>::PolyBasis(const std::vector<double> &x, std::vector<double> &bo) {
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
