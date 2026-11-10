#ifndef FEMUS_CUT_FEM_WEIGHT_LS_HPP
#define FEMUS_CUT_FEM_WEIGHT_LS_HPP

#include "GeomElTypeEnum.hpp"
#include "GaussPoints.hpp"
#include "gram_smithLS.hpp"
#include "FELevelSetMoments.hpp"
#include "Hex27LevelSetMoments.hpp"
#include "Hex27DirectionSelection.hpp"
#include "AdaptiveLineSliceMoments.hpp"
#include "Quad9DirectionSelection.hpp"
#include "Tri7DirectionSelection.hpp"
#include "Tet15LevelSet.hpp"
#include "Tet15LevelSetMoments.hpp"
#include "UnitElementMoments.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
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

template <class TypeIO, class TypeA>
class CutFemWeight_LS {
  public:
    CutFemWeight_LS(const GeomElType &geomElemType, const unsigned &qM, const std::string &gaussType) : _geomElemType(ValidateGeometry(geomElemType)), _qM(qM), _dim(_geomElemType == LINE ? 1U : (_geomElemType == QUAD || _geomElemType == TRI ? 2U : 3U)), _gauss(_geomElemType, 2U * _qM, gaussType.c_str()), _gn(_gauss.GetGaussPointsNumber()), _L(1) {
      for(unsigned d = 0; d < _dim; ++d) {
        _L *= _qM + 1U + d;
        _L /= d + 1U;
      }
      _fPair.resize(_L, 2);
      BuildMomentToWeight();
    }

    void SetIntegrationDirection(const unsigned direction) {
      const unsigned maximum = _geomElemType == LINE ? 1U : (_geomElemType == QUAD ? 2U : (_geomElemType == TRI ? 3U : (_geomElemType == TET ? 4U : 6U)));
      if(direction > maximum) throw std::invalid_argument("Invalid level-set integration direction.");
      if(direction == _integrationDirection) return;
      _integrationDirection = direction;
      _nodalWeightCacheValid = false;
    }

    void SetDirectionValidation(const unsigned mode) {
      if(mode > 2) throw std::invalid_argument("Level-set direction validation must be 0, 1, or 2.");
      if(mode == _directionValidation) return;
      _directionValidation = mode;
      _nodalWeightCacheValid = false;
    }

    void ClearCache() {
      _nodalWeightCacheValid = false;
      _nodalValues.clear();
      _nodalGeometryCoordinates.clear();
      _nodalPositiveWeight.clear();
      _nodalDeltaWeight.clear();
    }

    void GetWeightsFromNodalValues(const std::vector<double>& nodalValues, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool referenceSurface = false) {
      GetNodalWeightPair(nodalValues, positiveWeight, deltaWeight, absoluteTolerance, relativeTolerance, referenceSurface);
    }

    // Surface moments include |cof(J) grad G| and use the reference Gauss weight during assembly.
    void GetPhysicalSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodalValues, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& surfaceWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      GetNodalWeightPair(nodalValues, positiveWeight, surfaceWeight, absoluteTolerance, relativeTolerance, true, &coordinates);
    }

    unsigned GetCutFEMQuadratureOrder() const {
      return _qM;
    }

    unsigned GetGaussQuadratureOrder() const {
      return 2U * _qM;
    }

    unsigned GetDimension() const {
      return _dim;
    }

    unsigned GetGaussQuadraturePointNumber() const {
      return _gn;
    }

    const double* GetGaussWeightPointer() const {
      return _gauss.GetGaussWeightsPointer();
    }

    const double* GetGaussCoordinatePointer(const unsigned &k) const {
      return _gauss.GetGaussCoordinatePointer(k);
    }

  private:
    static GeomElType ValidateGeometry(const GeomElType geomElemType) {
      if(geomElemType != LINE && geomElemType != QUAD && geomElemType != TRI && geomElemType != HEX && geomElemType != TET) throw std::invalid_argument("Unsupported level-set geometry.");
      return geomElemType;
    }

    // Evaluates the unit-monomial Gram inverse through its shifted Legendre factorization.
    void BuildTensorMomentToWeight() {
      std::vector<std::vector<boost::multiprecision::cpp_bin_float_oct>> coefficient(_qM+1U,std::vector<boost::multiprecision::cpp_bin_float_oct>(_qM+1U,0));
      coefficient[0][0] = 1;
      if(_qM != 0) { coefficient[1][0] = -1; coefficient[1][1] = 2; }
      for(unsigned n = 2; n <= _qM; ++n) for(unsigned k = 0; k <= n; ++k) coefficient[n][k] = ((2U*n-1U)*(2*(k == 0 ? boost::multiprecision::cpp_bin_float_oct(0) : coefficient[n-1U][k-1U])-coefficient[n-1U][k])-(n-1U)*coefficient[n-2U][k])/n;
      _momentToWeight.resize(_gn,_L);
      for(unsigned ig = 0; ig < _gn; ++ig) {
        std::vector<std::vector<boost::multiprecision::cpp_bin_float_oct>> polynomials(_dim,std::vector<boost::multiprecision::cpp_bin_float_oct>(_qM+1U,1));
        for(unsigned d = 0; d < _dim; ++d) {
          const boost::multiprecision::cpp_bin_float_oct point(_gauss.GetGaussCoordinatePointer(d)[ig]);
          if(_qM != 0) polynomials[d][1] = point;
          for(unsigned n = 2; n <= _qM; ++n) polynomials[d][n] = ((2U*n-1U)*point*polynomials[d][n-1U]-(n-1U)*polynomials[d][n-2U])/n;
        }
        std::vector<boost::multiprecision::cpp_bin_float_oct> row(_L,0);
        for(unsigned i = 0; i <= _qM; ++i) for(unsigned j = 0; j <= (_dim > 1 ? _qM-i : 0U); ++j) for(unsigned k = 0; k <= (_dim > 2 ? _qM-i-j : 0U); ++k) {
          boost::multiprecision::cpp_bin_float_oct value = (2U*i+1U)*polynomials[0][i];
          if(_dim > 1) value *= (2U*j+1U)*polynomials[1][j];
          if(_dim > 2) value *= (2U*k+1U)*polynomials[2][k];
          for(unsigned a = 0; a <= i; ++a) for(unsigned b = 0; b <= j; ++b) for(unsigned c = 0; c <= k; ++c) {
            boost::multiprecision::cpp_bin_float_oct term = value*coefficient[i][a];
            if(_dim > 1) term *= coefficient[j][b];
            if(_dim > 2) term *= coefficient[k][c];
            row[FELevelSetMoments::UnitMomentIndex(_dim,a,b,c)] += term;
          }
        }
        for(unsigned n = 0; n < _L; ++n) _momentToWeight(ig,n) = static_cast<TypeA>(row[n]);
      }
    }

    void BuildMomentToWeight() {
      if(_geomElemType == LINE || _geomElemType == QUAD || _geomElemType == HEX) { BuildTensorMomentToWeight(); return; }
      std::vector<std::vector<boost::multiprecision::cpp_bin_float_oct>> ATAoct;
      CutFemWeightLSBasis::Get_GS_ATA_Matrix(_geomElemType, _qM, ATAoct);
      std::vector<const double*> xgp(_dim);
      for(unsigned d = 0; d < _dim; ++d) xgp[d] = _gauss.GetGaussCoordinatePointer(d);
      std::vector<boost::multiprecision::cpp_bin_float_oct> x(_dim);
      std::vector<std::vector<boost::multiprecision::cpp_bin_float_oct>> powers(_dim, std::vector<boost::multiprecision::cpp_bin_float_oct>(_qM + 1U, 1));
      _momentToWeight.resize(_gn, _L);
      for(unsigned ig = 0; ig < _gn; ++ig) {
        for(unsigned d = 0; d < _dim; ++d) x[d] = boost::multiprecision::cpp_bin_float_oct(xgp[d][ig]);
        for(unsigned d = 0; d < _dim; ++d) {
          for(unsigned i = 1; i <= _qM; ++i) powers[d][i] = powers[d][i - 1U] * x[d];
        }
        unsigned index = 0;
        std::vector<boost::multiprecision::cpp_bin_float_oct> basis(_L, 0);
        for(unsigned q = 0; q <= _qM; ++q) {
          if(_dim == 1) basis[index++] = powers[0][q];
          else {
            for(int ii = q; ii >= 0; --ii) {
              unsigned i = static_cast<unsigned>(ii);
              if(_dim == 2) basis[index++] = powers[0][i] * powers[1][q - i];
              else {
                for(int jj = q - i; jj >= 0; --jj) {
                  unsigned j = static_cast<unsigned>(jj);
                  basis[index++] = powers[0][i] * powers[1][j] * powers[2][q - i - j];
                }
              }
            }
          }
        }
        for(unsigned column = 0; column < _L; ++column) {
          boost::multiprecision::cpp_bin_float_oct value = 0;
          for(unsigned row = 0; row < _L; ++row) value += basis[row]*ATAoct[row][column];
          _momentToWeight(ig, column) = static_cast<TypeA>(value);
        }
      }
    }

    void GetNodalWeightPair(const std::vector<double>& nodalValues, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight, double absoluteTolerance, double relativeTolerance, bool referenceSurface, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
      if(&positiveWeight == &deltaWeight) throw std::invalid_argument("Heaviside and delta weights require distinct output vectors.");
      if(!(absoluteTolerance > 0.) || !(relativeTolerance >= 0.) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance)) throw std::invalid_argument("Invalid nodal integration tolerance.");
      if(geometryCoordinates != nullptr) FELevelSetMoments::ValidateSurfaceGeometryCoordinates(*geometryCoordinates, _dim, _geomElemType == LINE ? 3U : (_geomElemType == TRI ? 7U : (_geomElemType == QUAD ? 9U : (_geomElemType == TET ? 15U : 27U))));
      const bool sameGeometry = geometryCoordinates == nullptr ? _nodalGeometryCoordinates.empty() : _nodalGeometryCoordinates == *geometryCoordinates;
      if(_nodalWeightCacheValid && sameGeometry && referenceSurface == _nodalUsesReferenceSurface && nodalValues == _nodalValues && absoluteTolerance == _nodalAbsoluteTolerance && relativeTolerance == _nodalRelativeTolerance) {
        positiveWeight = _nodalPositiveWeight;
        deltaWeight = _nodalDeltaWeight;
        return;
      }
      _nodalWeightCacheValid = false;
      IntegrateNodalMomentPair(nodalValues, absoluteTolerance, relativeTolerance, referenceSurface, geometryCoordinates);
      _nodalValues = nodalValues;
      if(geometryCoordinates == nullptr) _nodalGeometryCoordinates.clear();
      else _nodalGeometryCoordinates = *geometryCoordinates;
      _nodalAbsoluteTolerance = absoluteTolerance;
      _nodalRelativeTolerance = relativeTolerance;
      _nodalUsesReferenceSurface = referenceSurface;
      EvaluateMomentWeightPair(positiveWeight, deltaWeight);
      _nodalPositiveWeight = positiveWeight;
      _nodalDeltaWeight = deltaWeight;
      _nodalWeightCacheValid = true;
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

    void IntegrateNodalMomentPair(const std::vector<double>& nodalValues, double absoluteTolerance, double relativeTolerance, bool referenceSurface, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
      if(_geomElemType == LINE || _geomElemType == QUAD || _geomElemType == HEX) {
        std::vector<long double> values;
        std::vector<long double> errors;
        FELevelSetMoments::IntegrateUnitTensorNodalMoments(nodalValues,_dim,_qM,values,errors,absoluteTolerance,relativeTolerance,referenceSurface,geometryCoordinates,_integrationDirection,_directionValidation);
        const long double surfaceFactor = geometryCoordinates != nullptr ? 1.0L/(1U<<_dim) : (referenceSurface ? 0.5L : 1.0L);
        for(unsigned n = 0; n < _L; ++n) { _fPair(n,0) = static_cast<TypeA>(values[n]); _fPair(n,1) = static_cast<TypeA>(surfaceFactor*values[_L+n]); }
      }
      else if(_geomElemType == TET) {
        std::vector<std::vector<std::vector<double>>> positive, delta;
        FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodalValues,_qM,positive,delta,absoluteTolerance,relativeTolerance,true,24,_integrationDirection,200000,nullptr,_directionValidation,nullptr,referenceSurface,nullptr,geometryCoordinates);
        StoreNodalMomentPair(positive,delta);
      }
      else if(_geomElemType == TRI) {
        std::vector<std::vector<double>> positive, delta;
        FELevelSetMoments::IntegrateTri7NodalSelected(nodalValues,_qM,positive,delta,absoluteTolerance,relativeTolerance,true,false,_directionValidation,_integrationDirection,nullptr,true,referenceSurface,nullptr,geometryCoordinates);
        StoreNodalMomentPair(positive,delta);
      }
    }

    void EvaluateMomentWeightPair(std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight) {
      std::vector<TypeIO> positiveResult(_gn);
      std::vector<TypeIO> deltaResult(_gn);
      for(unsigned ig = 0; ig < _gn; ++ig) {
        const TypeA* basisRow = _momentToWeight.data()+ig*_L;
        TypeA positive = 0.;
        TypeA delta = 0.;
        TypeA positiveCorrection = 0.;
        TypeA deltaCorrection = 0.;
        for(unsigned i = 0; i < _L; ++i) {
          const TypeA positiveTerm = basisRow[i]*_fPair(i,0)-positiveCorrection;
          const TypeA deltaTerm = basisRow[i]*_fPair(i,1)-deltaCorrection;
          const TypeA positiveNext = positive+positiveTerm;
          const TypeA deltaNext = delta+deltaTerm;
          positiveCorrection = (positiveNext-positive)-positiveTerm;
          deltaCorrection = (deltaNext-delta)-deltaTerm;
          positive = positiveNext;
          delta = deltaNext;
        }
        positiveResult[ig] = static_cast<TypeIO>(positive);
        deltaResult[ig] = static_cast<TypeIO>(delta);
        if(!Eigen::numext::isfinite(positiveResult[ig]) || !Eigen::numext::isfinite(deltaResult[ig])) throw std::overflow_error("A level-set quadrature correction exceeds the output range.");
      }
      positiveWeight.swap(positiveResult);
      deltaWeight.swap(deltaResult);
    }

    const GeomElType _geomElemType;
    const unsigned _qM;
    const unsigned _dim;
    femus::Gauss _gauss;
    const unsigned _gn;
    unsigned _L;
    Eigen::Matrix<TypeA, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _momentToWeight;
    Eigen::Matrix<TypeA, Eigen::Dynamic, 2> _fPair;
    unsigned _directionValidation = 1;
    unsigned _integrationDirection = 0;
    bool _nodalUsesReferenceSurface = false;
    bool _nodalWeightCacheValid = false;
    double _nodalAbsoluteTolerance = 0.;
    double _nodalRelativeTolerance = 0.;
    std::vector<double> _nodalValues;
    std::vector<std::vector<double>> _nodalGeometryCoordinates;
    std::vector<TypeIO> _nodalPositiveWeight;
    std::vector<TypeIO> _nodalDeltaWeight;
};

#endif
