#ifndef FEMUS_CUT_FEM_WEIGHT_LS_HPP
#define FEMUS_CUT_FEM_WEIGHT_LS_HPP

#include "GeomElTypeEnum.hpp"
#include "GaussPoints.hpp"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "gram_smithLS.hpp"
#include <Eigen/Dense>
#include "CutFemLevelSetData.hpp"
#include "CutFemPreparedMoments.hpp"


namespace boost {
  namespace multiprecision {
    typedef number < backends::cpp_bin_float < 24, backends::digit_base_2, void, boost::int16_t, -126, 127 >, et_off >         cpp_bin_float_single;
    typedef number < backends::cpp_bin_float < 53, backends::digit_base_2, void, boost::int16_t, -1022, 1023 >, et_off >       cpp_bin_float_double;
    typedef number < backends::cpp_bin_float < 64, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >     cpp_bin_float_double_extended;
    typedef number < backends::cpp_bin_float < 113, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >    cpp_bin_float_quad;
    typedef number < backends::cpp_bin_float < 237, backends::digit_base_2, void, boost::int32_t, -262142, 262143 >, et_off >  cpp_bin_float_oct;
  }
} // namespaces

template <class TypeIO, class TypeA>
class CutFemWeight_LS {
  public:
    CutFemWeight_LS(const GeomElType &geomElemType, const unsigned &qM, const std::string &gaussType) : _geomElemType(ValidateGeometry(geomElemType)), _qM(qM), _dim(_geomElemType == LINE ? 1U : (_geomElemType == QUAD || _geomElemType == TRI ? 2U : 3U)), _gauss(_geomElemType, 2U * _qM, gaussType.c_str()), _gn(_gauss.GetGaussPointsNumber()), _L(1) {
      for(unsigned d = 0; d < _dim; ++d) {
        _L *= _qM + 1U + d;
        _L /= d + 1U;
      }
      _fPair.resize(_L, 2);
      FELevelSetMoments::LineIntegration::FullStripRules();
      FELevelSetMoments::GaussKronrodRules::All();
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

    void SetOuterQuadratureRules(const FELevelSetMoments::OuterQuadratureRules& rules) {
      rules.Validate();
      if(_outerQuadratureRules.SameAs(rules)) return;
      _outerQuadratureRules = rules;
      _nodalWeightCacheValid = false;
    }

    const FELevelSetMoments::OuterQuadratureRules& GetOuterQuadratureRules() const {
      return _outerQuadratureRules;
    }

    void ClearCache() {
      _nodalWeightCacheValid = false;
      _nodalValues.clear();
      _nodalGeometryCoordinates.clear();
      _nodalPositiveWeight.clear();
      for(unsigned mode = 0; mode < _nodalDeltaWeightCache.size(); ++mode) _nodalDeltaWeightCache[mode].clear();
    }

    void GetWeightsFromNodalValues(const std::vector<double>& nodalValues, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool referenceSurface = false) {
      GetNodalWeightPair(nodalValues, positiveWeight, deltaWeight, absoluteTolerance, relativeTolerance, referenceSurface, _outerQuadratureRules);
    }

    // Surface moments include |cof(J) grad G| and use the reference Gauss weight during assembly.
    void GetPhysicalSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodalValues, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& surfaceWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      GetNodalWeightPair(nodalValues, positiveWeight, surfaceWeight, absoluteTolerance, relativeTolerance, true, _outerQuadratureRules, &coordinates);
    }

    void GetWeightsFromNodalValues(const std::vector<double>& nodalValues, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight, double absoluteTolerance, double relativeTolerance, bool referenceSurface, const FELevelSetMoments::OuterQuadratureRules& rules) {
      GetNodalWeightPair(nodalValues,positiveWeight,deltaWeight,absoluteTolerance,relativeTolerance,referenceSurface,rules);
    }

    void GetPhysicalSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodalValues, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& surfaceWeight, double absoluteTolerance, double relativeTolerance, const FELevelSetMoments::OuterQuadratureRules& rules) {
      GetNodalWeightPair(nodalValues,positiveWeight,surfaceWeight,absoluteTolerance,relativeTolerance,true,rules,&coordinates);
    }

    // Fits delta(G)/|cof(J) grad_ref G| and uses only reference Gauss weights during assembly.
    void GetTangentialSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodalValues, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& tangentialWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      GetNodalWeightPair(nodalValues,positiveWeight,tangentialWeight,absoluteTolerance,relativeTolerance,true,_outerQuadratureRules,&coordinates,true);
    }

    void GetTangentialSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodalValues, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& tangentialWeight, double absoluteTolerance, double relativeTolerance, const FELevelSetMoments::OuterQuadratureRules& rules) {
      GetNodalWeightPair(nodalValues,positiveWeight,tangentialWeight,absoluteTolerance,relativeTolerance,true,rules,&coordinates,true);
    }


    void GetWeightsFromExternalMoments(const std::vector<long double>& positiveMoments, const std::vector<long double>& surfaceMoments, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& surfaceWeight) {

      if(positiveMoments.size() != _L || surfaceMoments.size() != _L) throw std::invalid_argument("Invalid external moment size.");

      if(&positiveWeight == &surfaceWeight) throw std::invalid_argument("Weight outputs must be distinct.");

      for(unsigned n = 0; n < _L; ++n) if(!std::isfinite(positiveMoments[n]) || !std::isfinite(surfaceMoments[n])) throw std::invalid_argument("Nonfinite external moment.");
      ClearCache();
      _haveIntegratedMoments = false;
      for(unsigned n = 0; n < _L; ++n) { _fPair(n,0) = static_cast<TypeA>(positiveMoments[n]); _fPair(n,1) = static_cast<TypeA>(surfaceMoments[n]); }

      EvaluateMomentWeightPair(positiveWeight,surfaceWeight);
    }

    // Returns -1 or +1 only for a certified phase; zero retains every possible cut.
    int GetLevelSetPhaseFromNodalValues(const std::vector<double>& nodalValues) const {
      _preparedData.Reset(_geomElemType,nodalValues,_maximumCertificationDepth);
      return _preparedData.GetPhase();
    }

    // Fits the three signed measures delta(G) (cof(J) grad_ref G)_a / |cof(J) grad_ref G|.
    void GetNormalSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodalValues, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positiveWeight, std::vector<std::vector<TypeIO>>& normalWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      GetNormalSurfaceWeightsFromNodalGeometry(nodalValues,coordinates,positiveWeight,normalWeight,absoluteTolerance,relativeTolerance,_outerQuadratureRules);
    }

    void GetNormalSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodalValues, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positiveWeight, std::vector<std::vector<TypeIO>>& normalWeight, double absoluteTolerance, double relativeTolerance, const FELevelSetMoments::OuterQuadratureRules& rules) {
      if(_dim != 3) throw std::invalid_argument("Three normal surface weights require a three-dimensional element.");
      for(unsigned component = 0; component < normalWeight.size(); ++component) if(&positiveWeight == &normalWeight[component]) throw std::invalid_argument("Volume and normal weights require distinct output vectors.");
      std::vector<TypeIO> positive;
      std::vector<std::vector<TypeIO>> result(3);
      for(unsigned component = 0; component < 3; ++component) GetNodalWeightPair(nodalValues,positive,result[component],absoluteTolerance,relativeTolerance,true,rules,&coordinates,false,static_cast<int>(component));
      positiveWeight.swap(positive);
      normalWeight.swap(result);
    }

    void GetWeightsFromLevelSetData(const CutFemLevelSetData& data, std::vector<TypeIO>& positive, std::vector<TypeIO>& delta, double absolute = 1e-12, double relative = 1e-12, bool referenceSurface = false) { GetPreparedWeightPair(data,positive,delta,absolute,relative,referenceSurface,_outerQuadratureRules); }
    void GetPhysicalSurfaceWeightsFromLevelSetData(const CutFemLevelSetData& data, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positive, std::vector<TypeIO>& surface, double absolute = 1e-12, double relative = 1e-12) { GetPreparedWeightPair(data,positive,surface,absolute,relative,true,_outerQuadratureRules,&coordinates); }
    void GetTangentialSurfaceWeightsFromLevelSetData(const CutFemLevelSetData& data, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positive, std::vector<TypeIO>& surface, double absolute = 1e-12, double relative = 1e-12) { GetPreparedWeightPair(data,positive,surface,absolute,relative,true,_outerQuadratureRules,&coordinates,true); }
    void GetNormalSurfaceWeightsFromLevelSetData(const CutFemLevelSetData& data, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positive, std::vector<std::vector<TypeIO>>& surface, double absolute = 1e-12, double relative = 1e-12) {
      if(_dim != 3) throw std::invalid_argument("Three normal surface weights require a three-dimensional element.");
      for(unsigned component = 0; component < surface.size(); ++component) if(&positive == &surface[component]) throw std::invalid_argument("Volume and normal weights require distinct output vectors.");
      std::vector<TypeIO> resultPositive;
      std::vector<std::vector<TypeIO>> result(3);
      for(unsigned component = 0; component < 3; ++component) GetPreparedWeightPair(data,resultPositive,result[component],absolute,relative,true,_outerQuadratureRules,&coordinates,false,static_cast<int>(component));
      positive.swap(resultPositive);
      surface.swap(result);
    }

    // Returns the last actually integrated pair; cached weight returns do not replace it.
    void GetLastIntegratedMoments(std::vector<long double>& positive, std::vector<long double>& surface) const {
      if(!_haveIntegratedMoments || &positive == &surface) throw std::logic_error("Distinct outputs and a completed integration are required.");
      positive.resize(_L);
      surface.resize(_L);
      for(unsigned n = 0; n < _L; ++n) { positive[n] = static_cast<long double>(_fPair(n,0)); surface[n] = static_cast<long double>(_fPair(n,1)); }
    }

    const FELevelSetMoments::PreparedStatistics& GetPreparedStatistics() const { return _preparedWorkspace.statistics; }
    const CutFemLevelSetData& GetPreparedLevelSetData() const { if(!_preparedData.IsInitialized()) throw std::logic_error("No level-set record has been prepared."); return _preparedData; }
    void SetMaximumCertificationDepth(unsigned depth) { if(depth > 32) throw std::invalid_argument("Invalid certification depth."); _maximumCertificationDepth = depth; if(_preparedData.IsInitialized()) _preparedData.SetMaximumDepth(depth); }

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
    void GetPreparedWeightPair(const CutFemLevelSetData& data, std::vector<TypeIO>& positive, std::vector<TypeIO>& surface, double absolute, double relative, bool referenceSurface, const FELevelSetMoments::OuterQuadratureRules& rules, const std::vector<std::vector<double>>* coordinates = nullptr, bool inverseSurface = false, int normalComponent = -1) {
      if(data.GetGeometry() != _geomElemType || data.GetDimension() != _dim || data.GetRootDimension() != _dim || data.GetNodalValues().empty()) throw std::invalid_argument("Prepared weights require the matching native nodal element record.");
      const CutFemLevelSetData* previous = _activePreparedData;
      _activePreparedData = &data;
      try { GetNodalWeightPair(data.GetNodalValues(),positive,surface,absolute,relative,referenceSurface,rules,coordinates,inverseSurface,normalComponent); }
      catch(...) { _activePreparedData = previous; throw; }
      _activePreparedData = previous;
    }

    static GeomElType ValidateGeometry(const GeomElType geomElemType) {
      if(geomElemType != LINE && geomElemType != QUAD && geomElemType != TRI && geomElemType != HEX && geomElemType != TET && geomElemType != WEDGE) throw std::invalid_argument("Unsupported level-set geometry.");
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
        for(unsigned d = 0; d < _dim; ++d) x[d] = _geomElemType == WEDGE && d == 2 ? (boost::multiprecision::cpp_bin_float_oct(xgp[d][ig])+1)/2 : boost::multiprecision::cpp_bin_float_oct(xgp[d][ig]);
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

    void GetNodalWeightPair(const std::vector<double>& nodalValues, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight, double absoluteTolerance, double relativeTolerance, bool referenceSurface, const FELevelSetMoments::OuterQuadratureRules& rules, const std::vector<std::vector<double>>* geometryCoordinates = nullptr, bool inverseSurface = false, int normalComponent = -1) {
      const char* geometryNames[6] = {"HEX27","TET15","WEDGE21","QUAD9","TRI7","LINE3"};
      FELevelSetMoments::QuadratureWarningScope warningScope(std::string(geometryNames[static_cast<unsigned>(_geomElemType)])+" q="+std::to_string(_qM)+(inverseSurface ? " tangential surface" : (referenceSurface ? " surface" : " delta")));
      rules.Validate();
      if(&positiveWeight == &deltaWeight) throw std::invalid_argument("Heaviside and delta weights require distinct output vectors.");
      if(!(absoluteTolerance > 0.) || !(relativeTolerance >= 0.) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance)) throw std::invalid_argument("Invalid nodal integration tolerance.");
      if(geometryCoordinates != nullptr) FELevelSetMoments::ValidateSurfaceGeometryCoordinates(*geometryCoordinates, _dim, _geomElemType == LINE ? 3U : (_geomElemType == TRI ? 7U : (_geomElemType == QUAD ? 9U : (_geomElemType == TET ? 15U : (_geomElemType == WEDGE ? 21U : 27U)))));
      const bool sameGeometry = geometryCoordinates == nullptr ? _nodalGeometryCoordinates.empty() : _nodalGeometryCoordinates == *geometryCoordinates;
      const bool sameInput = _nodalWeightCacheValid && _cachedOuterQuadratureRules.SameAs(rules) && sameGeometry && nodalValues == _nodalValues && absoluteTolerance == _nodalAbsoluteTolerance && relativeTolerance == _nodalRelativeTolerance;
      const unsigned mode = normalComponent >= 0 ? 3U+static_cast<unsigned>(normalComponent) : (inverseSurface ? 2U : (referenceSurface ? 1U : 0U));
      if(sameInput && !_nodalDeltaWeightCache[mode].empty()) {
        positiveWeight = _nodalPositiveWeight;
        deltaWeight = _nodalDeltaWeightCache[mode];
        return;
      }
      const bool reusePositive = sameInput && _dim == 3;
      if(!sameInput) for(unsigned cachedMode = 0; cachedMode < _nodalDeltaWeightCache.size(); ++cachedMode) _nodalDeltaWeightCache[cachedMode].clear();
      _nodalWeightCacheValid = false;
      IntegrateNodalMomentPair(nodalValues,absoluteTolerance,relativeTolerance,referenceSurface,rules,geometryCoordinates,inverseSurface,normalComponent,!reusePositive);
      _cachedOuterQuadratureRules = rules;
      _nodalValues = nodalValues;
      if(geometryCoordinates == nullptr) _nodalGeometryCoordinates.clear();
      else _nodalGeometryCoordinates = *geometryCoordinates;
      _nodalAbsoluteTolerance = absoluteTolerance;
      _nodalRelativeTolerance = relativeTolerance;
      EvaluateMomentWeightPair(positiveWeight,deltaWeight,reusePositive,absoluteTolerance,relativeTolerance);
      _nodalPositiveWeight = positiveWeight;
      _nodalDeltaWeightCache[mode] = deltaWeight;
      _nodalWeightCacheValid = true;
    }

    void IntegrateNodalMomentPair(const std::vector<double>& nodalValues, double absoluteTolerance, double relativeTolerance, bool referenceSurface, const FELevelSetMoments::OuterQuadratureRules& rules, const std::vector<std::vector<double>>* geometryCoordinates = nullptr, bool inverseSurface = false, int normalComponent = -1, bool computePositive = true) {
      if(_activePreparedData == nullptr) _preparedData.Reset(_geomElemType,nodalValues,_maximumCertificationDepth);
      const CutFemLevelSetData& data = _activePreparedData == nullptr ? _preparedData : *_activePreparedData;
      const bool tensor = _geomElemType == LINE || _geomElemType == QUAD || _geomElemType == HEX;
      const long double measureScale = data.GetUnitToReferenceMeasureScale();
      const long double surfaceFactor = normalComponent >= 0 ? 1.0L : (geometryCoordinates != nullptr ? (inverseSurface ? measureScale : 1/measureScale) : (referenceSurface && tensor ? 0.5L : 1.0L));
      if(inverseSurface) absoluteTolerance /= measureScale*measureScale;
      if(_geomElemType == LINE) {
        std::vector<long double> values, errors;
        FELevelSetMoments::IntegrateLine(data,_qM,referenceSurface,inverseSurface,values,errors,computePositive);
        FELevelSetMoments::CheckQuadratureTolerance("Prepared LINE moments","moment error exceeds requested tolerance",values,errors,absoluteTolerance,absoluteTolerance,relativeTolerance);
        for(unsigned n = 0; n < _L; ++n) { _fPair(n,0) = static_cast<TypeA>(values[n]); _fPair(n,1) = static_cast<TypeA>(surfaceFactor*values[_L+n]); }
        _haveIntegratedMoments = true;
        return;
      }
      if(_dim == 2) {
        std::vector<long double> values, errors;
        FELevelSetMoments::IntegratePlanar(data,_qM,absoluteTolerance,relativeTolerance,referenceSurface,geometryCoordinates,_integrationDirection,_directionValidation,rules,inverseSurface,normalComponent,computePositive,_preparedWorkspace,values,errors);
        for(unsigned n = 0; n < _L; ++n) { if(computePositive) _fPair(n,0) = static_cast<TypeA>(values[n]); _fPair(n,1) = static_cast<TypeA>(surfaceFactor*values[_L+n]); }
        _haveIntegratedMoments = true;
        return;
      }
      std::vector<long double> values, errors;
      FELevelSetMoments::IntegrateSpatial(data,_qM,absoluteTolerance,relativeTolerance,referenceSurface,geometryCoordinates,_integrationDirection,_directionValidation,rules,inverseSurface,normalComponent,computePositive,_preparedWorkspace,values,errors);
      for(unsigned n = 0; n < _L; ++n) { if(computePositive) _fPair(n,0) = static_cast<TypeA>(values[n]); _fPair(n,1) = static_cast<TypeA>(surfaceFactor*values[_L+n]); }
      _haveIntegratedMoments = true;
    }

    void EvaluateMomentWeightPair(std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight, bool reusePositive = false, double absoluteTolerance = 0, double relativeTolerance = 0) {
      std::vector<TypeIO> positiveResult = reusePositive ? _nodalPositiveWeight : std::vector<TypeIO>(_gn);
      std::vector<TypeIO> deltaResult(_gn);
      const bool checkRounding = absoluteTolerance > 0 && absoluteTolerance+relativeTolerance < 8*static_cast<long double>(std::numeric_limits<TypeIO>::epsilon());
      long double positiveRounding = 0, deltaRounding = 0;
      const long double measureScale = _geomElemType == LINE || _geomElemType == QUAD || _geomElemType == HEX ? static_cast<long double>(1U<<_dim) : (_geomElemType == WEDGE ? 2 : 1);
      for(unsigned ig = 0; ig < _gn; ++ig) {
        const TypeA* basisRow = _momentToWeight.data()+ig*_L;
        TypeA positive = 0.;
        TypeA delta = 0.;
        TypeA positiveCorrection = 0.;
        TypeA deltaCorrection = 0.;
        for(unsigned i = 0; i < _L; ++i) {
          const TypeA deltaTerm = basisRow[i]*_fPair(i,1)-deltaCorrection;
          const TypeA deltaNext = delta+deltaTerm;
          deltaCorrection = (deltaNext-delta)-deltaTerm;
          delta = deltaNext;
          if(!reusePositive) {
            const TypeA positiveTerm = basisRow[i]*_fPair(i,0)-positiveCorrection;
            const TypeA positiveNext = positive+positiveTerm;
            positiveCorrection = (positiveNext-positive)-positiveTerm;
            positive = positiveNext;
          }
        }
        if(!reusePositive) positiveResult[ig] = static_cast<TypeIO>(positive);
        deltaResult[ig] = static_cast<TypeIO>(delta);
        if(!Eigen::numext::isfinite(positiveResult[ig]) || !Eigen::numext::isfinite(deltaResult[ig])) throw std::overflow_error("A level-set quadrature correction exceeds the output range.");
        if(checkRounding) {
          const long double weight = std::abs(GetGaussWeightPointer()[ig])/measureScale;
          deltaRounding += weight*std::abs(static_cast<long double>(delta-static_cast<TypeA>(deltaResult[ig])));
          if(!reusePositive) positiveRounding += weight*std::abs(static_cast<long double>(positive-static_cast<TypeA>(positiveResult[ig])));
        }
      }
      if(checkRounding) {
        long double positiveTolerance = std::numeric_limits<long double>::infinity(), deltaTolerance = positiveTolerance;
        for(unsigned moment = 0; moment < _L; ++moment) { positiveTolerance = std::min(positiveTolerance,absoluteTolerance+relativeTolerance*std::abs(static_cast<long double>(_fPair(moment,0)))); deltaTolerance = std::min(deltaTolerance,absoluteTolerance+relativeTolerance*std::abs(static_cast<long double>(_fPair(moment,1)))); }
        if(positiveRounding > positiveTolerance) FELevelSetMoments::WarnQuadrature("FEM output","rounded positive weights exceed the requested moment error bound",0,positiveRounding,positiveTolerance);
        if(deltaRounding > deltaTolerance) FELevelSetMoments::WarnQuadrature("FEM output","rounded surface weights exceed the requested moment error bound",_L,deltaRounding,deltaTolerance);
      }
      positiveWeight.swap(positiveResult);
      deltaWeight.swap(deltaResult);
    }

    mutable CutFemLevelSetData _preparedData;
    const CutFemLevelSetData* _activePreparedData = nullptr;
    bool _haveIntegratedMoments = false;
    FELevelSetMoments::PreparedMomentWorkspace _preparedWorkspace;
    unsigned _maximumCertificationDepth = 1;

    const GeomElType _geomElemType;
    const unsigned _qM;
    const unsigned _dim;
    femus::Gauss _gauss;
    const unsigned _gn;
    unsigned _L;
    Eigen::Matrix<TypeA, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _momentToWeight;
    Eigen::Matrix<TypeA, Eigen::Dynamic, 2> _fPair;
    FELevelSetMoments::OuterQuadratureRules _outerQuadratureRules;
    FELevelSetMoments::OuterQuadratureRules _cachedOuterQuadratureRules;
    unsigned _directionValidation = 1;
    unsigned _integrationDirection = 0;
    bool _nodalWeightCacheValid = false;
    double _nodalAbsoluteTolerance = 0.;
    double _nodalRelativeTolerance = 0.;
    std::vector<double> _nodalValues;
    std::vector<std::vector<double>> _nodalGeometryCoordinates;
    std::vector<TypeIO> _nodalPositiveWeight;
    std::array<std::vector<TypeIO>,6> _nodalDeltaWeightCache;
};




#endif
