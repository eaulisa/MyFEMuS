#ifndef FEMUS_CUT_FEM_WEIGHT_LS_HPP
#define FEMUS_CUT_FEM_WEIGHT_LS_HPP

#include "GeomElTypeEnum.hpp"
#include "GaussPoints.hpp"
#include "gram_smithLS.hpp"
#include "PlanarMoments_LS.hpp"
#include "FELevelSetMoments.hpp"
#include "Hex27LevelSetMoments.hpp"
#include "Hex27DirectionSelection.hpp"
#include "AdaptiveLineSliceMoments.hpp"
#include "Quad9DirectionSelection.hpp"
#include "Tri7DirectionSelection.hpp"
#include "Tet15LevelSet.hpp"
#include "Tet15LevelSetMoments.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <memory>
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
      if(_geomElemType == TET) {
        _obj.reset(new CutFemWeightLSDetail::TTImap<TypeA, TypeA>());
        _obj->SetBaseType(1);
      }
      else if(_geomElemType == TRI) _obj.reset(new CutFemWeightLSDetail::TRImap<TypeA, TypeA>());
      else _obj.reset(new CutFemWeightLSDetail::HCImap<TypeA, TypeA>());
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
      _obj->clear();
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

    void BuildMomentToWeight() {
      std::vector<std::vector<boost::multiprecision::cpp_bin_float_oct>> ATAoct;
      CutFemWeightLSBasis::Get_GS_ATA_Matrix(_geomElemType, _qM, ATAoct);
      Eigen::Matrix<TypeA, Eigen::Dynamic, Eigen::Dynamic> gramInverse(ATAoct.size(), ATAoct.size());
      for(unsigned i = 0; i < ATAoct.size(); ++i) {
        for(unsigned j = 0; j < ATAoct[i].size(); ++j) gramInverse(i, j) = static_cast<TypeA>(ATAoct[i][j]);
      }
      std::vector<const double*> xgp(_dim);
      for(unsigned d = 0; d < _dim; ++d) xgp[d] = _gauss.GetGaussCoordinatePointer(d);
      std::vector<TypeA> x(_dim);
      std::vector<std::vector<TypeA>> powers(_dim, std::vector<TypeA>(_qM + 1U, TypeA(1)));
      Eigen::Matrix<TypeA, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> basis(_gn, _L);
      for(unsigned ig = 0; ig < _gn; ++ig) {
        if(_geomElemType == LINE || _geomElemType == QUAD || _geomElemType == HEX) {
          for(unsigned d = 0; d < _dim; ++d) x[d] = (TypeA(1) + TypeA(xgp[d][ig])) / TypeA(2);
        }
        else if(_geomElemType == TRI) {
          x[0] = TypeA(1) - TypeA(xgp[0][ig]);
          x[1] = TypeA(xgp[1][ig]);
        }
        else {
          x[0] = TypeA(xgp[0][ig]) + TypeA(xgp[1][ig]) + TypeA(xgp[2][ig]);
          x[1] = TypeA(xgp[1][ig]) + TypeA(xgp[2][ig]);
          x[2] = TypeA(xgp[2][ig]);
        }
        for(unsigned d = 0; d < _dim; ++d) {
          for(unsigned i = 1; i <= _qM; ++i) powers[d][i] = powers[d][i - 1U] * x[d];
        }
        unsigned index = 0;
        for(unsigned q = 0; q <= _qM; ++q) {
          if(_dim == 1) basis(ig, index++) = powers[0][q];
          else {
            for(int ii = q; ii >= 0; --ii) {
              unsigned i = static_cast<unsigned>(ii);
              if(_dim == 2) basis(ig, index++) = powers[0][i] * powers[1][q - i];
              else {
                for(int jj = q - i; jj >= 0; --jj) {
                  unsigned j = static_cast<unsigned>(jj);
                  basis(ig, index++) = powers[0][i] * powers[1][j] * powers[2][q - i - j];
                }
              }
            }
          }
        }
      }
      _momentToWeight.noalias() = basis * gramInverse;
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

    void ScalePlanarReferenceSurfaceMoments(const std::vector<TypeIO>& normal) {
      TypeA squaredNorm = 0;
      for(unsigned d = 0; d < _dim; ++d) squaredNorm += static_cast<TypeA>(normal[d])*static_cast<TypeA>(normal[d]);
      TypeA norm = Eigen::numext::sqrt(squaredNorm);
      for(unsigned i = 0; i < _L; ++i) _fPair(i,1) *= norm;
    }

    void IntegrateLineReferenceSurfaceMomentPair(double a, double b, double c, bool zeroLeft, bool zeroRight) {
      long double scale = std::max(std::abs(static_cast<long double>(a)),std::max(std::abs(static_cast<long double>(b)),std::abs(static_cast<long double>(c))));
      if(scale == 0.0L) throw std::domain_error("The delta of an identically zero field is undefined.");
      long double quadratic = static_cast<long double>(a)/scale;
      long double linear = static_cast<long double>(b)/scale;
      long double constant = static_cast<long double>(c)/scale;
      std::vector<long double> positive;
      std::vector<long double> delta;
      FELevelSetMoments::Integrate1DInternal(quadratic,linear,constant,_qM,0.0L,1.0L,false,false,positive,delta);
      std::vector<long double> roots;
      if(zeroLeft && zeroRight) roots = {0.0L,1.0L};
      else if(zeroLeft) {
        roots.push_back(0.0L);
        if(quadratic != 0.0L) roots.push_back(-linear/quadratic);
      }
      else if(zeroRight) {
        roots.push_back(1.0L);
        if(quadratic != 0.0L) roots.push_back(constant/quadratic);
      }
      else roots = FELevelSetMoments::QuadraticRoots(quadratic,linear,constant);
      for(unsigned rootIndex = 0; rootIndex < roots.size(); ++rootIndex) {
        long double point = roots[rootIndex];
        if(point < 0.0L || point > 1.0L) continue;
        long double derivative = zeroRight && point == 1.0L ? quadratic-constant : 2.0L*quadratic*point+linear;
        if(derivative == 0.0L) throw std::domain_error("A repeated 1D zero has no regular surface moment.");
        long double weight = point == 0.0L || point == 1.0L ? 0.25L : 0.5L;
        long double power = 1.0L;
        for(unsigned i = 0; i <= _qM; ++i) {
          delta[i] += weight*power;
          power *= point;
        }
      }
      for(unsigned i = 0; i <= _qM; ++i) {
        _fPair(i,0) = static_cast<TypeA>(static_cast<double>(positive[i]));
        _fPair(i,1) = static_cast<TypeA>(delta[i]);
      }
    }

    void IntegrateNodalMomentPair(const std::vector<double>& nodalValues, double absoluteTolerance, double relativeTolerance, bool referenceSurface, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
      if(_geomElemType == TET) {
        std::vector<std::vector<std::vector<double>>> coefficient;
        FELevelSetMoments::Tet15ReferenceCoefficientsFromNodalValues(nodalValues,coefficient);
        bool planar = true;
        for(unsigned i = 0; i < 5; ++i) {
          for(unsigned j = 0; j < 5; ++j) {
            for(unsigned k = 0; k < 5; ++k) if(i+j+k > 1 && coefficient[i][j][k] != 0.0) planar = false;
          }
        }
        if(geometryCoordinates == nullptr && planar && (coefficient[1][0][0] != 0.0 || coefficient[0][1][0] != 0.0 || coefficient[0][0][1] != 0.0)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(coefficient[1][0][0]),static_cast<TypeIO>(coefficient[0][1][0]),static_cast<TypeIO>(coefficient[0][0][1])};
          IntegratePlanarMomentPair(normal,static_cast<TypeIO>(coefficient[0][0][0]));
          if(referenceSurface) ScalePlanarReferenceSurfaceMoments(normal);
        }
        else {
          std::vector<std::vector<std::vector<double>>> positive;
          std::vector<std::vector<std::vector<double>>> delta;
          FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodalValues,_qM,positive,delta,absoluteTolerance,relativeTolerance,true,24,_integrationDirection,200000,nullptr,_directionValidation,nullptr,referenceSurface,nullptr,geometryCoordinates);
          StoreNodalMomentPair(positive,delta);
        }
      }
      else if(_geomElemType == HEX) {
        std::vector<std::vector<std::vector<double>>> coefficient;
        FELevelSetMoments::Hex27UnitCoefficientsFromNodalValues(nodalValues,coefficient);
        bool planar = true;
        for(unsigned i = 0; i < 3; ++i) {
          for(unsigned j = 0; j < 3; ++j) {
            for(unsigned k = 0; k < 3; ++k) if(i+j+k > 1 && coefficient[i][j][k] != 0.) planar = false;
          }
        }
        if(geometryCoordinates == nullptr && planar && (coefficient[1][0][0] != 0. || coefficient[0][1][0] != 0. || coefficient[0][0][1] != 0.)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(0.5L*coefficient[1][0][0]),static_cast<TypeIO>(0.5L*coefficient[0][1][0]),static_cast<TypeIO>(0.5L*coefficient[0][0][1])};
          IntegratePlanarMomentPair(normal,static_cast<TypeIO>(static_cast<long double>(coefficient[0][0][0])+0.5L*coefficient[1][0][0]+0.5L*coefficient[0][1][0]+0.5L*coefficient[0][0][1]));
          if(referenceSurface) ScalePlanarReferenceSurfaceMoments(normal);
          return;
        }
        std::vector<std::vector<std::vector<double>>> positive, delta;
        FELevelSetMoments::IntegrateHex27NodalSelected(nodalValues,_qM,positive,delta,absoluteTolerance,relativeTolerance,true,_directionValidation,_integrationDirection,nullptr,true,nullptr,referenceSurface,nullptr,geometryCoordinates);
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
        if(geometryCoordinates == nullptr && planar && (coefficient[1][0] != 0. || coefficient[0][1] != 0.)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(0.5L * coefficient[1][0]), static_cast<TypeIO>(0.5L * coefficient[0][1])};
          IntegratePlanarMomentPair(normal, static_cast<TypeIO>(static_cast<long double>(coefficient[0][0]) + 0.5L * coefficient[1][0] + 0.5L * coefficient[0][1]));
          if(referenceSurface) ScalePlanarReferenceSurfaceMoments(normal);
          return;
        }
        std::vector<std::vector<double>> positive;
        std::vector<std::vector<double>> delta;
        FELevelSetMoments::IntegrateQuad9NodalSelected(nodalValues, _qM, positive, delta, absoluteTolerance, relativeTolerance, true, false, _directionValidation, _integrationDirection, nullptr, true, referenceSurface, nullptr, geometryCoordinates);
        StoreNodalMomentPair(positive, delta);
      }
      else if(_geomElemType == TRI) {
        std::vector<double> coefficient;
        FELevelSetMoments::Tri7ReflectedCoefficientsFromNodalValues(nodalValues, coefficient);
        if(geometryCoordinates == nullptr && coefficient[3] == 0. && coefficient[4] == 0. && coefficient[5] == 0. && coefficient[6] == 0. && (coefficient[1] != 0. || coefficient[2] != 0.)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(-coefficient[1]), static_cast<TypeIO>(coefficient[2])};
          IntegratePlanarMomentPair(normal, static_cast<TypeIO>(static_cast<long double>(coefficient[0]) + coefficient[1]));
          if(referenceSurface) ScalePlanarReferenceSurfaceMoments(normal);
          return;
        }
        std::vector<std::vector<double>> positive;
        std::vector<std::vector<double>> delta;
        FELevelSetMoments::IntegrateTri7NodalSelected(nodalValues, _qM, positive, delta, absoluteTolerance, relativeTolerance, true, false, _directionValidation, _integrationDirection, nullptr, true, referenceSurface, nullptr, geometryCoordinates);
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
          if(referenceSurface) ScalePlanarReferenceSurfaceMoments(normal);
          return;
        }
        if(referenceSurface) {
          IntegrateLineReferenceSurfaceMomentPair(a,b,c,nodalValues[0] == 0.0,nodalValues[1] == 0.0);
          return;
        }
        std::vector<double> positive;
        std::vector<double> delta;
        FELevelSetMoments::Integrate1D(a, b, c, _qM, positive, delta);
        for(unsigned i = 0; i <= _qM; ++i) {
          _fPair(i, 0) = static_cast<TypeA>(positive[i]);
          _fPair(i, 1) = static_cast<TypeA>(delta[i]);
        }
      }
      else {
        throw std::invalid_argument("Paired nodal weights support LINE3, QUAD9, TRI7, HEX27, and TET15.");
      }
    }

    void EvaluateMomentWeightPair(std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight) {
      positiveWeight.resize(_gn);
      deltaWeight.resize(_gn);
      for(unsigned ig = 0; ig < _gn; ++ig) {
        const TypeA* basisRow = _momentToWeight.data()+ig*_L;
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

    const GeomElType _geomElemType;
    const unsigned _qM;
    const unsigned _dim;
    femus::Gauss _gauss;
    const unsigned _gn;
    unsigned _L;
    std::unique_ptr<CutFemWeightLSDetail::CutFEMmap<TypeA>> _obj;
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
