#ifndef FEMUS_CUT_FEM_LEVEL_SET_DATA_HPP
#define FEMUS_CUT_FEM_LEVEL_SET_DATA_HPP

#include "GeomElTypeEnum.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class CutFemLevelSetData {
  public:
    struct LineRoot {
      long double coordinate = 0;
      long double derivative = 0;
      unsigned multiplicity = 0;
      int boundary = 0;
    };

    struct LineAnalysis {
      unsigned degree = 0;
      unsigned rootCount = 0;
      unsigned positiveIntervalCount = 0;
      bool identicallyZero = false;
      bool compensated = false;
      std::array<LineRoot,2> roots;
      std::array<std::array<long double,2>,2> positiveIntervals = {};
    };

    struct GraphCertificate {
      std::array<long double,2> bounds = {0,0};
      int sign = 0;
      bool strict = false;
      int zeroFiber = 0;
    };

    struct FiberProjection {
      unsigned dimension = 0;
      unsigned degree = 4;
      std::array<std::vector<long double>,3> coefficients;
      std::vector<long double> discriminant;
      std::array<std::array<long double,2>,5> compensatedDiscriminant = {};
      bool compensatedReady = false;
    };

    // FEMuS nodal ordering is preserved; all polynomial coordinates are unit coordinates.
    CutFemLevelSetData(GeomElType geometry, const std::vector<double>& nodalValues, unsigned maximumDepth = 1) {
      Reset(geometry, nodalValues, maximumDepth);
    }

    // An empty instance is reusable output storage for Restrict().
    CutFemLevelSetData() = default;

    void Reset(GeomElType geometry, const std::vector<double>& nodalValues, unsigned maximumDepth = 1) {
      ValidateGeometry(geometry);
      ValidateDepth(maximumDepth);
      if(nodalValues.size() != NodeCount(geometry)) throw std::invalid_argument("Incorrect FEMuS level-set node count.");
      for(double value : nodalValues) if(!std::isfinite(value)) throw std::invalid_argument("Level-set nodal values must be finite.");
      if(_initialized && _root && geometry == _geometry && nodalValues == _nodalValues) {
        SetMaximumDepth(maximumDepth);
        return;
      }
      CutFemLevelSetData replacement;
      replacement._geometry = geometry;
      replacement._dimension = Dimension(geometry);
      replacement._rootDimension = replacement._dimension;
      replacement._parentDimension = replacement._dimension;
      replacement._maximumDepth = maximumDepth;
      replacement._nodalValues = nodalValues;
      replacement._initialized = true;
      replacement._root = true;
      replacement._nativeSpace = true;
      replacement._meshLevel = _meshLevel;
      replacement._elementId = _elementId;
      replacement._nodalVersion = _nodalVersion + 1;
      replacement._revision = _revision + 1;
      replacement._parentOrigin.assign(replacement._dimension, 0);
      replacement._rootOrigin.assign(replacement._dimension, 0);
      replacement._parentMap.assign(replacement._dimension, std::vector<long double>(replacement._dimension, 0));
      replacement._rootMap = replacement._parentMap;
      for(unsigned axis = 0; axis < replacement._dimension; ++axis) {
        replacement._parentMap[axis][axis] = 1;
        replacement._rootMap[axis][axis] = 1;
      }
      replacement.SetBernsteinOrders();
      replacement.Classify();
      for(unsigned i = 0; i < _workCounts.size(); ++i) replacement._workCounts[i] += _workCounts[i];
      *this = std::move(replacement);
    }

    // Imports unit-coordinate coefficients directly for polynomial-level integration and testing.
    void ResetPowerCoefficients(GeomElType geometry, const std::vector<long double>& coefficient, unsigned maximumDepth = 1) {
      ValidateGeometry(geometry);
      ValidateDepth(maximumDepth);
      const unsigned dimension = Dimension(geometry);
      if(coefficient.size() != PolynomialSize(dimension,2)) throw std::invalid_argument("Incorrect unit power-coefficient count.");
      long double maximum = 0;
      for(unsigned n = 0; n < coefficient.size(); ++n) { if(!std::isfinite(coefficient[n])) throw std::invalid_argument("Power coefficients must be finite."); const std::array<unsigned,3> power = PowerExponents(n,dimension,2); if(coefficient[n] != 0 && ((geometry == TET && power[0]+power[1]+power[2] > 4) || (geometry == WEDGE && power[0]+power[1] > 3))) throw std::invalid_argument("Polynomial exceeds the supported reference layout."); maximum = std::max(maximum,std::abs(coefficient[n])); }
      CutFemLevelSetData replacement;
      replacement._initialized = true;
      replacement._geometry = geometry;
      replacement._dimension = replacement._rootDimension = replacement._parentDimension = dimension;
      replacement._maximumDepth = maximumDepth;
      replacement._power = coefficient;
      replacement._powerReady = true;
      replacement._nodalVersion = _nodalVersion+1;
      replacement._revision = _revision+1;
      replacement._meshLevel = _meshLevel;
      replacement._elementId = _elementId;
      replacement._parentOrigin.assign(dimension,0);
      replacement._rootOrigin.assign(dimension,0);
      replacement._parentMap.assign(dimension,std::vector<long double>(dimension,0));
      for(unsigned axis = 0; axis < dimension; ++axis) { replacement._parentMap[axis][axis] = 1; Differentiate(coefficient,dimension,2,axis,replacement._rootDerivatives[axis]); }
      replacement._rootMap = replacement._parentMap;
      replacement.SetBernsteinOrders();
      if(geometry == TRI && coefficient[8] != 0) replacement._bernsteinOrders = {4,4,0};
      if(maximum > 0) { int exponent = 0; std::frexp(maximum,&exponent); replacement._normalization = std::ldexp(1.0L,exponent-1); }
      replacement._workCounts = _workCounts;
      *this = std::move(replacement);
    }

    void SetMaximumDepth(unsigned maximumDepth) {
      RequireInitialized();
      ValidateDepth(maximumDepth);
      if(maximumDepth == _maximumDepth && _classified) return;
      _maximumDepth = maximumDepth;
      _classified = false;
      Classify();
    }

    void SetContext(unsigned meshLevel, std::uint64_t elementId, std::uint64_t nodalVersion) {
      RequireInitialized();
      _meshLevel = meshLevel;
      _elementId = elementId;
      _nodalVersion = nodalVersion;
    }
    bool IsInitialized() const {
      return _initialized;
    }
    GeomElType GetGeometry() const {
      RequireInitialized();
      return _geometry;
    }
    unsigned GetDimension() const {
      RequireInitialized();
      return _dimension;
    }
    unsigned GetRootDimension() const {
      RequireInitialized();
      return _rootDimension;
    }
    unsigned GetMaximumDepth() const {
      RequireInitialized();
      return _maximumDepth;
    }
    unsigned GetMeshLevel() const {
      RequireInitialized();
      return _meshLevel;
    }
    std::uint64_t GetElementId() const {
      RequireInitialized();
      return _elementId;
    }
    std::uint64_t GetNodalVersion() const {
      RequireInitialized();
      return _nodalVersion;
    }
    std::uint64_t GetRevision() const {
      RequireInitialized();
      return _revision;
    }
    int GetPhase() const {
      Classify();
      return _phase;
    }
    bool IsCutCandidate() const {
      return GetPhase() == 0;
    }
    bool IsDegenerate() const {
      RequireInitialized();
      return _degenerate;
    }
    std::string GetClassification() const {
      Classify();
      return _classification;
    }
    // True for a nodal polynomial in the active chart, rather than an ancestor-only source.
    bool HasNodalPolynomialSource() const { RequireInitialized(); return _root || (_source && _restrictionSteps == 0 && _source->nodes != 0); }
    bool HasPowerCoefficients() const {
      return _powerReady;
    }
    bool HasBernsteinCoefficients() const {
      return _bernsteinReady;
    }
    const std::vector<double>& GetNodalValues() const {
      RequireInitialized();
      return _nodalValues;
    }

    // Power coefficient index is i*3^(dimension-1)+j*3^(dimension-2)+k for the active coordinates.
    const std::vector<long double>& GetPowerCoefficients() const {
      EnsurePower();
      return _power;
    }
    const std::vector<long double>& GetNormalizedPowerCoefficients() const {
      EnsurePower();
      if(!_normalizedReady) {
        _normalized.resize(_power.size());
        for(unsigned i = 0; i < _power.size(); ++i) _normalized[i] = _power[i] / _normalization;
        _normalizedReady = true;
      }
      return _normalized;
    }
    long double GetNormalizationScale() const {
      EnsurePower();
      return _normalization;
    }
    const std::vector<long double>& GetBernsteinCoefficients() const {
      EnsureBernstein();
      return _bernstein;
    }
    std::array<unsigned, 3> GetBernsteinOrders() const {
      RequireInitialized();
      return _bernsteinOrders;
    }
    const std::vector<std::array<unsigned, 3>>& GetBernsteinIndices() const {
      EnsureBernstein();
      return _bernsteinIndices;
    }

    // Counts are nodal expansions, power-to-Bernstein conversions, local derivatives, and restrictions performed.
    std::array<unsigned long long, 4> GetWorkCounts() const {
      return _workCounts;
    }

    long double Evaluate(const std::vector<long double>& point) const {
      EnsurePower();
      ValidatePoint(point);
      return EvaluateUnchecked(_power.data(), _dimension, 2, point.data());
    }

    long double EvaluateBernstein(const std::vector<long double>& point) const {
      EnsureBernstein();
      ValidatePoint(point);
      std::array<long double,3> unit = {0,0,0};
      for(unsigned axis = 0; axis < _dimension; ++axis) unit[axis] = point[axis];
      const bool simplex = _geometry == TRI || _geometry == TET || _geometry == WEDGE;
      const unsigned baseDimension = _geometry == TET ? 3U : 2U;
      static const std::array<std::array<unsigned,35>,3> multipliers = BuildBernsteinMultipliers();
      const unsigned multiplierLayout = _geometry == TET ? 1U : (_bernsteinOrders[0] == 4 ? 2U : 0U);
      std::array<std::array<long double,5>,3> factors;
      std::array<long double,5> remainingPowers;
      remainingPowers[0] = 1;
      if(simplex) {
        long double remaining = 1;
        for(unsigned axis = 0; axis < baseDimension; ++axis) remaining -= unit[axis];
        for(unsigned power = 1; power <= _bernsteinOrders[0]; ++power) remainingPowers[power] = remainingPowers[power-1]*remaining;
      }
      for(unsigned axis = 0; axis < _dimension; ++axis) {
        if(simplex && axis < baseDimension) { factors[axis][0] = 1; for(unsigned power = 1; power <= _bernsteinOrders[axis]; ++power) factors[axis][power] = factors[axis][power-1]*unit[axis]; }
        else for(unsigned power = 0; power <= _bernsteinOrders[axis]; ++power) factors[axis][power] = Bernstein(_bernsteinOrders[axis],power,unit[axis]);
      }
      long double result = 0, correction = 0;
      for(unsigned n = 0; n < _bernstein.size(); ++n) {
        const std::array<unsigned,3>& index = _bernsteinIndices[n];
        long double value = _bernstein[n];
        if(simplex) { unsigned sum = 0; for(unsigned axis = 0; axis < baseDimension; ++axis) sum += index[axis]; const unsigned baseIndex = _geometry == WEDGE ? n/(_bernsteinOrders[2]+1) : n; value *= multipliers[multiplierLayout][baseIndex]*remainingPowers[_bernsteinOrders[0]-sum]; }
        for(unsigned axis = 0; axis < _dimension; ++axis) value *= factors[axis][index[axis]];
        Accumulate(value,result,correction);
      }
      return result;
    }

    std::array<long double, 2> GetBounds() const {
      EnsureBernstein();
      return {_bernsteinExtrema[0]-_bernsteinPadding,_bernsteinExtrema[1]+_bernsteinPadding};
    }

    const std::vector<long double>& GetDerivativeCoefficients(unsigned axis) const {
      EnsurePower();
      if(axis >= _dimension) throw std::out_of_range("Derivative axis exceeds the active dimension.");
      if(!_derivativeReady[axis]) {
        Differentiate(_power, _dimension, 2, axis, _derivatives[axis]);
        _derivativeReady[axis] = true;
        ++_workCounts[2];
      }
      return _derivatives[axis];
    }

    long double EvaluateDerivative(unsigned axis, const std::vector<long double>& point) const {
      ValidatePoint(point);
      return EvaluateUnchecked(GetDerivativeCoefficients(axis).data(), _dimension, 2, point.data());
    }
    long double EvaluateSecondDerivative(unsigned first, unsigned second, const std::vector<long double>& point) const {
      ValidatePoint(point);
      if(first >= _dimension || second >= _dimension) throw std::out_of_range("Derivative axis exceeds the active dimension.");
      if(first > second) std::swap(first, second);
      const unsigned index = first * (7 - first) / 2 + second - first;
      if(!_secondDerivativeReady[index]) {
        Differentiate(GetDerivativeCoefficients(first), _dimension, 2, second, _secondDerivatives[index]);
        _secondDerivativeReady[index] = true;
        ++_workCounts[2];
      }
      return EvaluateUnchecked(_secondDerivatives[index].data(), _dimension, 2, point.data());
    }

    std::vector<long double> GetGradient(const std::vector<long double>& point) const {
      std::array<long double,3> result;
      GetGradient(point,result);
      return std::vector<long double>(result.begin(),result.begin()+_dimension);
    }

    void GetGradient(const std::vector<long double>& point, std::array<long double,3>& result) const {
      ValidatePoint(point);
      result.fill(0);
      for(unsigned axis = 0; axis < _dimension; ++axis) result[axis] = EvaluateUnchecked(GetDerivativeCoefficients(axis).data(),_dimension,2,point.data());
    }

    // Root gradients retain every original reference-coordinate component after dimension reduction.
    std::vector<long double> GetRootGradient(const std::vector<long double>& point) const {
      std::array<long double,3> result;
      GetRootGradient(point,result);
      return std::vector<long double>(result.begin(),result.begin()+_rootDimension);
    }

    void GetRootGradient(const std::vector<long double>& point, std::array<long double,3>& result) const {
      ValidatePoint(point);
      result.fill(0);
      for(unsigned axis = 0; axis < _rootDimension; ++axis) result[axis] = EvaluateUnchecked(GetRootDerivative(axis).data(),_dimension,2,point.data());
    }
    const std::vector<long double>& GetRootDerivativeCoefficients(unsigned axis) const { return GetRootDerivative(axis); }

    std::array<long double, 2> GetDerivativeBounds(unsigned axis) const {
      RequireInitialized();
      if(axis >= _dimension) throw std::out_of_range("Derivative axis exceeds the active dimension.");
      if(!_derivativeBoundsReady[axis]) {
        EnsureBernstein();
        const bool simplexAxis = _geometry == TRI || _geometry == TET || (_geometry == WEDGE && axis < 2);
        const unsigned degree = _bernsteinOrders[axis];
        const long double factor = degree;
        long double minimum = std::numeric_limits<long double>::infinity(), maximum = -minimum, norm = 0;
        long double cornerMinimum = minimum, cornerMaximum = maximum;
        for(unsigned n = 0; n < _bernsteinIndices.size(); ++n) {
          std::array<unsigned,3> index = _bernsteinIndices[n];
          const unsigned total = _geometry == TET ? index[0]+index[1]+index[2] : index[0]+index[1];
          if(simplexAxis ? total >= degree : index[axis] >= degree) continue;
          ++index[axis];
          const long double value = factor*(_bernstein[BernsteinIndex(_geometry,_bernsteinOrders,index)]-_bernstein[n]);
          if(_geometry == QUAD && (index[1-axis] == 0 || index[1-axis] == _bernsteinOrders[1-axis]) && (index[axis] == 1 || index[axis] == degree)) { cornerMinimum = std::min(cornerMinimum,value); cornerMaximum = std::max(cornerMaximum,value); }
          minimum = std::min(minimum,value);
          maximum = std::max(maximum,value);
          norm += std::abs(value);
        }
        const long double padding = 8192*std::numeric_limits<long double>::epsilon()*norm+2*degree*_bernsteinPadding;
        if(_geometry == HEX && minimum != maximum && minimum <= padding && maximum >= -padding) RefineHexDerivativeBounds(axis,padding,minimum,maximum);
        _derivativeBounds[axis] = {minimum-padding,maximum+padding};
        _derivativeCornerSigns[axis] = static_cast<unsigned char>((cornerMinimum < -padding ? 1U : 0U)|(cornerMaximum > padding ? 2U : 0U));
        _derivativeBoundsReady[axis] = true;
      }
      return _derivativeBounds[axis];
    }

    int GetGraphDirection() const {
      RequireInitialized();
      if(_graphDirectionReady) return _graphDirection;
      int result = -1;
      long double best = -1;
      for(unsigned axis = 0; axis < _dimension; ++axis) {
        std::array<long double,2> bounds = GetDerivativeBounds(axis);
        if(_geometry == QUAD && bounds[0] <= 0 && bounds[1] >= 0 && _derivativeCornerSigns[axis] != 3) bounds = PrepareGraphCertificate(axis,false).bounds;
        const long double lower = bounds[0] > 0 ? bounds[0] : (bounds[1] < 0 ? -bounds[1] : 0);
        const long double upper = std::max(std::abs(bounds[0]), std::abs(bounds[1]));
        const long double score = upper > 0 ? lower / upper : 0;
        if(lower > 0 && score > best) {
          best = score;
          result = static_cast<int>(axis);
        }
      }
      _graphDirection = result;
      _graphDirectionReady = true;
      return _graphDirection;
    }

    const GraphCertificate& GetGraphCertificate(unsigned axis) const;
    int GetZeroFiberStatus(unsigned axis) const;
    int GetGraphDirection(bool allowWeak) const;
    bool IsIdenticallyZero() const;
    bool DependsOn(unsigned axis) const;
    int GetReflectionParity(unsigned axis) const;
    unsigned GetBoundaryCount() const;
    bool IsZeroOnBoundary(unsigned boundary) const;
    unsigned GetZeroBoundaryMask() const;
    void RestrictBoundary(unsigned boundary, CutFemLevelSetData& child) const;
    void PermuteAxes(const std::array<unsigned,3>& axes, CutFemLevelSetData& chart) const;
    void PermuteSimplexVertices(const std::array<unsigned,4>& vertices, CutFemLevelSetData& chart) const;
    void TransformParentComponents(unsigned degree, const std::vector<std::vector<long double>>& components, unsigned& outputDegree, std::vector<std::vector<long double>>& output) const;
    const FiberProjection& GetFiberProjection(unsigned axis) const;
    bool GetCompensatedDiscriminantRoots(unsigned axis, long double left, long double right, std::vector<long double>& roots) const;
    void AnalyzeLocalFiber(unsigned axis, long double coordinate, LineAnalysis& analysis) const;
    long double EvaluateLocalFiberDiscriminant(unsigned axis, long double coordinate) const;
    void GetLocalShiftedFiberDiscriminant(unsigned axis, long double endpoint, std::array<long double,5>& coefficient) const;
    const std::vector<long double>& GetUnitPowerCoefficients() const;
    static void ExpandNodalPolynomial(GeomElType geometry, const std::vector<long double>& nodal, std::vector<long double>& coefficient) { ValidateGeometry(geometry); if(nodal.size() != NodeCount(geometry)) throw std::invalid_argument("Incorrect nodal polynomial size."); for(long double value : nodal) if(!std::isfinite(value)) throw std::invalid_argument("Nodal polynomial values must be finite."); ExpandUnit(geometry,nodal,coefficient); }
    long double GetUnitToReferenceMeasureScale() const;
    long double EvaluateFiberDiscriminant(unsigned axis, const std::vector<long double>& point) const;
    long double EvaluateFiberDiscriminant(unsigned axis, long double point) const;
    void GetFiberBoundaryPolynomial(unsigned axis, long double constant, long double slope, std::vector<long double>& coefficient) const;
    void GetShiftedFiberDiscriminant(unsigned axis, long double endpoint, std::array<long double,5>& coefficient) const;

    const std::vector<long double>& GetParentOrigin() const {
      RequireInitialized();
      return _parentOrigin;
    }
    const std::vector<std::vector<long double>>& GetRootMap() const { RequireInitialized(); return _rootMap; }
    const std::vector<std::vector<long double>>& GetParentMap() const {
      RequireInitialized();
      return _parentMap;
    }
    long double GetParentMeasureScale() const {
      RequireInitialized();
      return _parentMeasure;
    }
    std::vector<long double> MapToParent(const std::vector<long double>& point) const {
      ValidatePoint(point);
      return ApplyMap(_parentOrigin, _parentMap, point);
    }
    std::vector<long double> MapToRoot(const std::vector<long double>& point) const {
      ValidatePoint(point);
      return ApplyMap(_rootOrigin, _rootMap, point);
    }

    // Restriction fixes an active reference coordinate and fills an existing lower-dimensional data object.
    void Restrict(unsigned axis, long double coordinate, CutFemLevelSetData& child) const { Restrict(axis,coordinate,child,_rootGradientAvailable); }
    void Restrict(unsigned axis, long double coordinate, CutFemLevelSetData& child, bool retainRootGradient) const;
    bool HasRootGradientContext() const { RequireInitialized(); return _rootGradientAvailable; }

    // External polynomial components keep their component count and their full supplied degree under restriction.
    void RestrictComponents(unsigned axis, long double coordinate, unsigned degree, const std::vector<std::vector<long double>>& components, std::vector<std::vector<long double>>& output) const {
      std::array<long double,3> origin, scales;
      GeomElType childGeometry;
      RestrictionMap(axis, coordinate, childGeometry, origin, scales);
      if(&components == &output) throw std::invalid_argument("Restriction input and output components must be distinct.");
      output.resize(components.size());
      for(unsigned component = 0; component < components.size(); ++component) {
        ValidatePolynomial(components[component], _dimension, degree);
        RestrictPower(components[component], _dimension, degree, axis, coordinate, origin, scales, output[component]);
      }
    }

    std::vector<long double> GetLineRoots() const;
    unsigned GetLineRoots(std::array<long double,2>& roots) const;
    const LineAnalysis& GetLineAnalysis() const;

    static long double GetRegularRootDeltaFactor(const LineRoot& root) {
      if(!std::isfinite(root.coordinate) || !std::isfinite(root.derivative)) throw std::invalid_argument("Root data must be finite.");
      if(root.multiplicity != 1 || root.derivative == 0) throw std::domain_error("A singular line root requires endpoint handling.");
      return 1/std::abs(root.derivative);
    }

    long double GetLineDeltaFactor(const LineRoot& root) const {
      RequireInitialized();
      if(_dimension != 1) throw std::logic_error("A root density requires line data.");
      return GetRegularRootDeltaFactor(root);
    }

    long double EvaluateLine(long double coordinate) const { EnsurePower(); ValidateLineCoordinate(coordinate); return EvaluateQuadratic(_power.data(),1,coordinate); }

    long double GetLineDeltaFactor(long double coordinate) const {
      EnsurePower();
      ValidateLineCoordinate(coordinate);
      if(_lineAnalysisReady) for(unsigned n = 0; n < _lineAnalysis->rootCount; ++n) if(coordinate == _lineAnalysis->roots[n].coordinate) return GetLineDeltaFactor(_lineAnalysis->roots[n]);
      const long double derivative = 2*_power[2]*coordinate+_power[1];
      if(derivative == 0) throw std::domain_error("A singular line root requires endpoint handling.");
      return 1/std::abs(derivative);
    }

    long double GetLineSurfaceData(long double coordinate, std::array<long double,3>& rootGradient) const {
      if(_lineAnalysisReady) for(unsigned n = 0; n < _lineAnalysis->rootCount; ++n) if(coordinate == _lineAnalysis->roots[n].coordinate) return GetLineSurfaceData(_lineAnalysis->roots[n],rootGradient);
      const long double inverseDerivative = GetLineDeltaFactor(coordinate);
      rootGradient.fill(0);
      for(unsigned axis = 0; axis < _rootDimension; ++axis) rootGradient[axis] = EvaluateQuadratic(GetRootDerivative(axis).data(),1,coordinate);
      return inverseDerivative;
    }

    long double GetLineReferenceSurfaceFactor(long double coordinate) const {
      std::array<long double,3> gradient;
      const long double inverseDerivative = GetLineSurfaceData(coordinate,gradient);
      return inverseDerivative*SurfaceMetric(gradient);
    }

    long double GetLinePhysicalSurfaceFactor(long double coordinate, const std::array<long double,3>& cofactorGradient) const { const long double inverseDerivative = GetLineDeltaFactor(coordinate); return inverseDerivative*SurfaceMetric(cofactorGradient); }
    long double GetLineTangentialSurfaceFactor(long double coordinate, const std::array<long double,3>& cofactorGradient) const { const long double inverseDerivative = GetLineDeltaFactor(coordinate); return inverseDerivative/SurfaceMetric(cofactorGradient); }

    long double GetLineNormalSurfaceFactor(long double coordinate, const std::array<long double,3>& cofactorGradient, unsigned component) const {
      ValidateSurfaceComponent(component);
      const long double inverseDerivative = GetLineDeltaFactor(coordinate);
      const long double norm = SurfaceMetric(cofactorGradient);
      return inverseDerivative*(cofactorGradient[component]/norm);
    }

    void GetLineNormalSurfaceFactors(long double coordinate, const std::array<long double,3>& cofactorGradient, std::array<long double,3>& normalFactors) const {
      const long double inverseDerivative = GetLineDeltaFactor(coordinate);
      const long double norm = SurfaceMetric(cofactorGradient);
      for(unsigned axis = 0; axis < normalFactors.size(); ++axis) normalFactors[axis] = axis < _rootDimension ? inverseDerivative*(cofactorGradient[axis]/norm) : 0;
    }

    long double GetLinePhysicalSurfaceFactorFromCofactor(long double coordinate, const std::array<std::array<long double,3>,3>& cofactor) const {
      std::array<long double,3> transformed;
      const long double inverseDerivative = GetLinePhysicalSurfaceData(coordinate,cofactor,transformed);
      return inverseDerivative*SurfaceMetric(transformed);
    }

    long double GetLineTangentialSurfaceFactorFromCofactor(long double coordinate, const std::array<std::array<long double,3>,3>& cofactor) const {
      std::array<long double,3> transformed;
      const long double inverseDerivative = GetLinePhysicalSurfaceData(coordinate,cofactor,transformed);
      return inverseDerivative/SurfaceMetric(transformed);
    }

    long double GetLineNormalSurfaceFactorFromCofactor(long double coordinate, const std::array<std::array<long double,3>,3>& cofactor, unsigned component) const {
      ValidateSurfaceComponent(component);
      std::array<long double,3> transformed;
      const long double inverseDerivative = GetLinePhysicalSurfaceData(coordinate,cofactor,transformed);
      const long double norm = SurfaceMetric(transformed);
      return inverseDerivative*(transformed[component]/norm);
    }

    void GetLineNormalSurfaceFactorsFromCofactor(long double coordinate, const std::array<std::array<long double,3>,3>& cofactor, std::array<long double,3>& normalFactors) const {
      std::array<long double,3> transformed;
      const long double inverseDerivative = GetLinePhysicalSurfaceData(coordinate,cofactor,transformed);
      const long double norm = SurfaceMetric(transformed);
      for(unsigned axis = 0; axis < normalFactors.size(); ++axis) normalFactors[axis] = axis < _rootDimension ? inverseDerivative*(transformed[axis]/norm) : 0;
    }

    // Scalars are delta, reference surface, physical surface, and inverse metric; normals are the signed delta measures.
    void GetLineAllSurfaceFactors(long double coordinate, const std::array<long double,3>& cofactorGradient, std::array<long double,4>& scalarFactors, std::array<long double,3>& normalFactors) const {
      std::array<long double,3> rootGradient;
      const long double inverseDerivative = GetLineSurfaceData(coordinate,rootGradient);
      CompleteSurfaceFactors(inverseDerivative,rootGradient,cofactorGradient,scalarFactors,normalFactors);
    }

    void GetLineAllSurfaceFactorsFromCofactor(long double coordinate, const std::array<std::array<long double,3>,3>& cofactor, std::array<long double,4>& scalarFactors, std::array<long double,3>& normalFactors) const {
      std::array<long double,3> rootGradient, transformed = {0,0,0};
      const long double inverseDerivative = GetLineSurfaceData(coordinate,rootGradient);
      TransformSurfaceGradient(rootGradient,cofactor,transformed);
      CompleteSurfaceFactors(inverseDerivative,rootGradient,transformed,scalarFactors,normalFactors);
    }

    long double GetLineSurfaceData(const LineRoot& root, std::array<long double,3>& rootGradient) const {
      const long double inverseDerivative = GetLineDeltaFactor(root);
      rootGradient.fill(0);
      for(unsigned axis = 0; axis < _rootDimension; ++axis) rootGradient[axis] = EvaluateQuadratic(GetRootDerivative(axis).data(),1,root.coordinate);
      unsigned selected = 0;
      for(unsigned axis = 1; axis < _rootDimension; ++axis) if(std::abs(_rootMap[axis][0]) > std::abs(_rootMap[selected][0])) selected = axis;
      if(_rootMap[selected][0] == 0) throw std::domain_error("A collapsed line chart has no surface density.");
      long double remainder = root.derivative;
      for(unsigned axis = 0; axis < _rootDimension; ++axis) if(axis != selected) remainder -= _rootMap[axis][0]*rootGradient[axis];
      rootGradient[selected] = remainder/_rootMap[selected][0];
      return inverseDerivative;
    }

    long double GetLineReferenceSurfaceFactor(const LineRoot& root) const {
      std::array<long double,3> gradient;
      const long double inverseDerivative = GetLineSurfaceData(root,gradient);
      return inverseDerivative*SurfaceMetric(gradient);
    }

    long double GetLinePhysicalSurfaceFactor(const LineRoot& root, const std::array<long double,3>& cofactorGradient) const { const long double inverseDerivative = GetLineDeltaFactor(root); return inverseDerivative*SurfaceMetric(cofactorGradient); }
    long double GetLineTangentialSurfaceFactor(const LineRoot& root, const std::array<long double,3>& cofactorGradient) const { const long double inverseDerivative = GetLineDeltaFactor(root); return inverseDerivative/SurfaceMetric(cofactorGradient); }

    long double GetLineNormalSurfaceFactor(const LineRoot& root, const std::array<long double,3>& cofactorGradient, unsigned component) const {
      ValidateSurfaceComponent(component);
      const long double inverseDerivative = GetLineDeltaFactor(root);
      const long double norm = SurfaceMetric(cofactorGradient);
      return inverseDerivative*(cofactorGradient[component]/norm);
    }

    void GetLineNormalSurfaceFactors(const LineRoot& root, const std::array<long double,3>& cofactorGradient, std::array<long double,3>& normalFactors) const {
      const long double inverseDerivative = GetLineDeltaFactor(root);
      const long double norm = SurfaceMetric(cofactorGradient);
      for(unsigned axis = 0; axis < normalFactors.size(); ++axis) normalFactors[axis] = axis < _rootDimension ? inverseDerivative*(cofactorGradient[axis]/norm) : 0;
    }

    long double GetLinePhysicalSurfaceFactorFromCofactor(const LineRoot& root, const std::array<std::array<long double,3>,3>& cofactor) const {
      std::array<long double,3> transformed;
      const long double inverseDerivative = GetLinePhysicalSurfaceData(root,cofactor,transformed);
      return inverseDerivative*SurfaceMetric(transformed);
    }

    long double GetLineTangentialSurfaceFactorFromCofactor(const LineRoot& root, const std::array<std::array<long double,3>,3>& cofactor) const {
      std::array<long double,3> transformed;
      const long double inverseDerivative = GetLinePhysicalSurfaceData(root,cofactor,transformed);
      return inverseDerivative/SurfaceMetric(transformed);
    }

    long double GetLineNormalSurfaceFactorFromCofactor(const LineRoot& root, const std::array<std::array<long double,3>,3>& cofactor, unsigned component) const {
      ValidateSurfaceComponent(component);
      std::array<long double,3> transformed;
      const long double inverseDerivative = GetLinePhysicalSurfaceData(root,cofactor,transformed);
      const long double norm = SurfaceMetric(transformed);
      return inverseDerivative*(transformed[component]/norm);
    }

    void GetLineNormalSurfaceFactorsFromCofactor(const LineRoot& root, const std::array<std::array<long double,3>,3>& cofactor, std::array<long double,3>& normalFactors) const {
      std::array<long double,3> transformed;
      const long double inverseDerivative = GetLinePhysicalSurfaceData(root,cofactor,transformed);
      const long double norm = SurfaceMetric(transformed);
      for(unsigned axis = 0; axis < normalFactors.size(); ++axis) normalFactors[axis] = axis < _rootDimension ? inverseDerivative*(transformed[axis]/norm) : 0;
    }

    // Scalars are delta, reference surface, physical surface, and inverse metric; normals are the signed delta measures.
    void GetLineAllSurfaceFactors(const LineRoot& root, const std::array<long double,3>& cofactorGradient, std::array<long double,4>& scalarFactors, std::array<long double,3>& normalFactors) const {
      std::array<long double,3> rootGradient;
      const long double inverseDerivative = GetLineSurfaceData(root,rootGradient);
      CompleteSurfaceFactors(inverseDerivative,rootGradient,cofactorGradient,scalarFactors,normalFactors);
    }

    void GetLineAllSurfaceFactorsFromCofactor(const LineRoot& root, const std::array<std::array<long double,3>,3>& cofactor, std::array<long double,4>& scalarFactors, std::array<long double,3>& normalFactors) const {
      std::array<long double,3> rootGradient, transformed = {0,0,0};
      const long double inverseDerivative = GetLineSurfaceData(root,rootGradient);
      TransformSurfaceGradient(rootGradient,cofactor,transformed);
      CompleteSurfaceFactors(inverseDerivative,rootGradient,transformed,scalarFactors,normalFactors);
    }

    static long double EvaluatePolynomial(const std::vector<long double>& coefficient, unsigned dimension, unsigned degree, const std::vector<long double>& point) {
      if(dimension < 1 || dimension > 3 || point.size() != dimension || coefficient.size() != PolynomialSize(dimension, degree)) throw std::invalid_argument("Invalid polynomial evaluation layout.");
      for(long double coordinate : point) if(!std::isfinite(coordinate)) throw std::invalid_argument("Polynomial evaluation coordinates must be finite.");
      return EvaluateUnchecked(coefficient.data(),dimension,degree,point.data());
    }

  private:
    void ValidateSurfaceComponent(unsigned component) const { RequireInitialized(); if(component >= _rootDimension) throw std::out_of_range("Surface component exceeds the original dimension."); }

    long double SurfaceMetric(const std::array<long double,3>& gradient) const {
      long double norm = 0;
      for(unsigned axis = 0; axis < _rootDimension; ++axis) { if(!std::isfinite(gradient[axis])) throw std::invalid_argument("Surface gradient components must be finite."); norm = std::hypot(norm,gradient[axis]); }
      if(norm == 0) throw std::domain_error("The surface metric is zero.");
      return norm;
    }

    void TransformSurfaceGradient(const std::array<long double,3>& gradient, const std::array<std::array<long double,3>,3>& cofactor, std::array<long double,3>& transformed) const {
      transformed.fill(0);
      for(unsigned physical = 0; physical < _rootDimension; ++physical) for(unsigned reference = 0; reference < _rootDimension; ++reference) { if(!std::isfinite(cofactor[physical][reference])) throw std::invalid_argument("Surface cofactor entries must be finite."); transformed[physical] += cofactor[physical][reference]*gradient[reference]; }
    }

    long double GetLinePhysicalSurfaceData(long double coordinate, const std::array<std::array<long double,3>,3>& cofactor, std::array<long double,3>& transformed) const {
      std::array<long double,3> gradient;
      const long double inverseDerivative = GetLineSurfaceData(coordinate,gradient);
      TransformSurfaceGradient(gradient,cofactor,transformed);
      return inverseDerivative;
    }

    long double GetLinePhysicalSurfaceData(const LineRoot& root, const std::array<std::array<long double,3>,3>& cofactor, std::array<long double,3>& transformed) const {
      std::array<long double,3> gradient;
      const long double inverseDerivative = GetLineSurfaceData(root,gradient);
      TransformSurfaceGradient(gradient,cofactor,transformed);
      return inverseDerivative;
    }

    void CompleteSurfaceFactors(long double inverseDerivative, const std::array<long double,3>& rootGradient, const std::array<long double,3>& cofactorGradient, std::array<long double,4>& scalarFactors, std::array<long double,3>& normalFactors) const {
      long double referenceNorm = 0, physicalNorm = 0;
      for(unsigned axis = 0; axis < _rootDimension; ++axis) { if(!std::isfinite(cofactorGradient[axis])) throw std::invalid_argument("Surface cofactor components must be finite."); referenceNorm = std::hypot(referenceNorm,rootGradient[axis]); physicalNorm = std::hypot(physicalNorm,cofactorGradient[axis]); }
      if(physicalNorm == 0) throw std::domain_error("The physical surface metric is zero.");
      scalarFactors = {inverseDerivative,referenceNorm*inverseDerivative,physicalNorm*inverseDerivative,inverseDerivative/physicalNorm};
      for(unsigned axis = 0; axis < normalFactors.size(); ++axis) normalFactors[axis] = axis < _rootDimension ? inverseDerivative*(cofactorGradient[axis]/physicalNorm) : 0;
    }

    struct BernsteinRestrictionPlan {
      bool simplex = false;
      unsigned degree = 0;
      std::array<unsigned,3> orders = {0,0,0};
      std::vector<std::array<unsigned,3>> indices;
      std::vector<std::array<unsigned,3>> terms;
      std::vector<unsigned> offsets;
    };

    static std::array<BernsteinRestrictionPlan,19> BuildRestrictionPlans();
    const BernsteinRestrictionPlan& GetRestrictionPlan(unsigned axis) const {
      static const std::array<BernsteinRestrictionPlan,19> plans = BuildRestrictionPlans();
      const unsigned base = _geometry == HEX ? 0U : (_geometry == TET ? 3U : (_geometry == WEDGE ? 6U : (_geometry == QUAD ? (_bernsteinOrders[0] == 3 ? 15U : (_bernsteinOrders[1] == 3 ? 17U : 9U)) : (_bernsteinOrders[0] == 4 ? 13U : 11U))));
      return plans[base+axis];
    }

    static long double EvaluateQuadratic(const long double* coefficient, unsigned stride, long double coordinate) { return (coefficient[2*stride]*coordinate+coefficient[stride])*coordinate+coefficient[0]; }

    static long double EvaluateUnchecked(const long double* coefficient, unsigned dimension, unsigned degree, const long double* point) {
      if(degree == 2) {
        if(dimension == 1) return EvaluateQuadratic(coefficient,1,point[0]);
        long double rows[3];
        for(unsigned i = 0; i < 3; ++i) {
          if(dimension == 2) rows[i] = EvaluateQuadratic(coefficient+3*i,1,point[1]);
          else { long double lines[3]; for(unsigned j = 0; j < 3; ++j) lines[j] = EvaluateQuadratic(coefficient+9*i+3*j,1,point[2]); rows[i] = EvaluateQuadratic(lines,1,point[1]); }
        }
        return EvaluateQuadratic(rows,1,point[0]);
      }
      const unsigned side = degree + 1;
      long double result = 0;
      if(dimension == 1) {
        result = coefficient[degree];
        for(unsigned i = degree; i > 0; --i) result = result * point[0] + coefficient[i - 1];
        return result;
      }
      for(unsigned i = side; i > 0; --i) {
        long double row = 0;
        for(unsigned j = side; j > 0; --j) {
          long double line;
          if(dimension == 2) line = coefficient[(i - 1) * side + j - 1];
          else {
            line = coefficient[((i - 1) * side + j - 1) * side + degree];
            for(unsigned k = degree; k > 0; --k) line = line * point[2] + coefficient[((i - 1) * side + j - 1) * side + k - 1];
          }
          row = j == side ? line : row * point[1] + line;
        }
        result = i == side ? row : result * point[0] + row;
      }
      return result;
    }

    static unsigned Dimension(GeomElType geometry) {
      return geometry == LINE ? 1U : (geometry == QUAD || geometry == TRI ? 2U : 3U);
    }
    static unsigned NodeCount(GeomElType geometry) {
      return geometry == LINE ? 3U : (geometry == QUAD ? 9U : (geometry == TRI ? 7U : (geometry == HEX ? 27U : (geometry == TET ? 15U : 21U))));
    }
    static void ValidateGeometry(GeomElType geometry) {
      if(geometry != LINE && geometry != QUAD && geometry != TRI && geometry != HEX && geometry != TET && geometry != WEDGE) throw std::invalid_argument("Unsupported level-set element.");
    }
    static void ValidateDepth(unsigned depth) {
      if(depth > 32) throw std::invalid_argument("Maximum detection depth must lie between 0 and 32.");
    }
    void RequireInitialized() const {
      if(!_initialized) throw std::logic_error("The level-set data object is not initialized.");
    }
    static unsigned PolynomialSize(unsigned dimension, unsigned degree) {
      if(dimension < 1 || dimension > 3 || degree > 32) throw std::invalid_argument("Invalid polynomial dimensions or degree.");
      unsigned result = 1;
      for(unsigned axis = 0; axis < dimension; ++axis) result *= degree + 1;
      return result;
    }
    static void ValidatePolynomial(const std::vector<long double>& coefficient, unsigned dimension, unsigned degree) {
      if(coefficient.size() != PolynomialSize(dimension, degree)) throw std::invalid_argument("Incorrect polynomial coefficient count.");
      for(long double value : coefficient) if(!std::isfinite(value)) throw std::invalid_argument("Polynomial coefficients must be finite.");
    }
    static long double Power(long double value, unsigned exponent) {
      long double result = 1;
      for(unsigned i = 0; i < exponent; ++i) result *= value;
      return result;
    }
    static long double Binomial(unsigned n, unsigned k) {
      if(k > n) return 0;
      k = std::min(k,n-k);
      if(k == 0) return 1;
      if(k == 1) return n;
      if(k == 2) return n*(n-1)/2U;
      if(k == 3) return n*(n-1)*(n-2)/6U;
      long double result = 1;
      for(unsigned i = 1; i <= k; ++i) result *= static_cast<long double>(n-i+1)/i;
      return result;
    }
    static long double Falling(unsigned n, unsigned k) {
      if(k > n) return 0;
      long double result = 1;
      for(unsigned i = 0; i < k; ++i) result *= n - i;
      return result;
    }
    static long double Bernstein(unsigned degree, unsigned index, long double coordinate) {
      return Binomial(degree, index) * Power(coordinate, index) * Power(1 - coordinate, degree - index);
    }
    static long double Multinomial(unsigned degree, const std::array<unsigned, 3>& index, unsigned dimension) {
      long double result = 1;
      for(unsigned axis = 0; axis < dimension; ++axis) {
        result *= Binomial(degree, index[axis]);
        degree -= index[axis];
      }
      return result;
    }
    static std::array<std::array<unsigned,35>,3> BuildBernsteinMultipliers() {
      const GeomElType geometries[3] = {TRI,TET,TRI};
      const std::array<std::array<unsigned,3>,3> orders = {{{3,3,0},{4,4,4},{4,4,0}}};
      std::array<std::array<unsigned,35>,3> result = {};
      for(unsigned layout = 0; layout < 3; ++layout) { std::vector<std::array<unsigned,3>> indices; BuildBernsteinIndices(geometries[layout],orders[layout],indices); for(unsigned n = 0; n < indices.size(); ++n) result[layout][n] = static_cast<unsigned>(Multinomial(orders[layout][0],indices[n],Dimension(geometries[layout]))); }
      return result;
    }
    static void Accumulate(long double term, long double& sum, long double& correction) {
      const long double adjusted = term - correction;
      const long double next = sum + adjusted;
      correction = (next - sum) - adjusted;
      sum = next;
    }
    static long double Padding(const std::vector<long double>& coefficient) {
      long double norm = 0;
      for(long double value : coefficient) norm += std::abs(value);
      return 8192 * std::numeric_limits<long double>::epsilon() * norm;
    }
    static unsigned PowerIndex(const std::array<unsigned, 3>& index, unsigned dimension, unsigned degree) {
      unsigned result = index[0];
      for(unsigned axis = 1; axis < dimension; ++axis) result = result * (degree + 1) + index[axis];
      return result;
    }
    static std::array<unsigned, 3> PowerExponents(unsigned index, unsigned dimension, unsigned degree) {
      std::array<unsigned, 3> result = {0, 0, 0};
      for(unsigned axis = dimension; axis > 0; --axis) {
        result[axis - 1] = index % (degree + 1);
        index /= degree + 1;
      }
      return result;
    }
    static std::vector<long double> ApplyMap(const std::vector<long double>& origin, const std::vector<std::vector<long double>>& map, const std::vector<long double>& point) {
      std::vector<long double> result = origin;
      for(unsigned row = 0; row < map.size(); ++row) for(unsigned column = 0; column < point.size(); ++column) result[row] += map[row][column] * point[column];
      return result;
    }

    void ValidatePoint(const std::vector<long double>& point) const {
      RequireInitialized();
      if(point.size() != _dimension) throw std::invalid_argument("Point dimension does not match the level-set data.");
      for(long double coordinate : point) if(!std::isfinite(coordinate)) throw std::invalid_argument("Point coordinates must be finite.");
      long double sum = 0;
      for(unsigned axis = 0; axis < _dimension; ++axis) {
        const bool simplexAxis = _geometry == TRI || _geometry == TET || (_geometry == WEDGE && axis < 2);
        const long double lower = 0;
        const long double padding = 64 * std::numeric_limits<long double>::epsilon();
        if(point[axis] < lower - padding || point[axis] > 1 + padding) throw std::out_of_range("Point lies outside the reference element.");
        if(simplexAxis) sum += point[axis];
      }
      if((_geometry == TRI || _geometry == TET || _geometry == WEDGE) && sum > 1 + 128 * std::numeric_limits<long double>::epsilon()) throw std::out_of_range("Point lies outside the reference simplex.");
    }

    void ValidateLineCoordinate(long double coordinate) const {
      RequireInitialized();
      if(_dimension != 1) throw std::logic_error("A line operation requires one-dimensional data.");
      if(!std::isfinite(coordinate)) throw std::invalid_argument("Line coordinates must be finite.");
      const long double padding = 64*std::numeric_limits<long double>::epsilon();
      if(coordinate < -padding || coordinate > 1+padding) throw std::out_of_range("Point lies outside the reference line.");
    }

    static void Differentiate(const std::vector<long double>& coefficient, unsigned dimension, unsigned degree, unsigned axis, std::vector<long double>& result) {
      result.assign(coefficient.size(), 0);
      for(unsigned n = 0; n < coefficient.size(); ++n) {
        std::array<unsigned, 3> index = PowerExponents(n, dimension, degree);
        if(index[axis] == 0) continue;
        const long double value = index[axis] * coefficient[n];
        --index[axis];
        result[PowerIndex(index, dimension, degree)] = value;
      }
    }

    void SetBernsteinOrders() {
      _bernsteinOrders = {0, 0, 0};
      for(unsigned axis = 0; axis < _dimension; ++axis) _bernsteinOrders[axis] = _geometry == TET ? 4U : ((_geometry == TRI || (_geometry == WEDGE && axis < 2)) ? 3U : 2U);
    }
    void EnsurePower() const;
    void EnsureBernstein() const;
    void RefineHexDerivativeBounds(unsigned axis, long double padding, long double& minimum, long double& maximum) const {
      const unsigned stride[3] = {9,3,1};
      const unsigned first = stride[(axis+1)%3], second = stride[(axis+2)%3];
      // The derivative is affine between these faces, and the shared midpoint control of each quadratic split is the mean of its neighbors.
      static constexpr std::array<unsigned,2> pairs[] = {{0,0},{2,2},{0,1},{1,2}};
      long double lower = std::numeric_limits<long double>::infinity(), upper = -lower;
      for(unsigned side = 0; side < 2; ++side) {
        long double face[3][3];
        for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) {
          const unsigned index = side*stride[axis]+i*first+j*second;
          face[i][j] = 2*(_bernstein[index+stride[axis]]-_bernstein[index]);
        }
        for(const std::array<unsigned,2>& row : pairs) for(const std::array<unsigned,2>& column : pairs) {
          const long double control = 0.25L*face[row[0]][column[0]]+0.25L*face[row[0]][column[1]]+0.25L*face[row[1]][column[0]]+0.25L*face[row[1]][column[1]];
          lower = std::min(lower,control);
          upper = std::max(upper,control);
          if(lower <= padding && upper >= -padding) return;
        }
      }
      minimum = lower;
      maximum = upper;
    }
    void UpdateBernsteinBounds(long double inheritedPadding) const {
      long double minimum = _bernstein[0], maximum = minimum, norm = 0;
      for(long double value : _bernstein) { minimum = std::min(minimum,value); maximum = std::max(maximum,value); norm += std::abs(value); }
      _bernsteinExtrema = {minimum,maximum};
      _bernsteinPadding = inheritedPadding+8192*std::numeric_limits<long double>::epsilon()*norm;
    }
    void Classify() const;
    const std::vector<long double>& GetRootDerivative(unsigned axis) const {
      if(axis >= _rootDimension) throw std::out_of_range("Root derivative axis exceeds the original dimension.");
      if(!_rootGradientAvailable) throw std::logic_error("This slice omitted original-gradient context; restrict again from a prepared ancestor.");
      if(_root) return GetDerivativeCoefficients(axis);
      return _rootDerivatives[axis];
    }
    template<class Real> static void ExpandTensor(GeomElType geometry, const std::vector<Real>& nodal, std::vector<long double>& coefficient);
    static void ExpandTriangle(const std::vector<long double>& nodal, std::vector<long double>& coefficient);
    template<class Real> static void ExpandTetrahedron(const std::vector<Real>& nodal, std::vector<long double>& coefficient);
    static void AddBarycentricProduct(const std::array<std::array<long double, 4>, 4>& barycentric, const std::vector<unsigned>& factors, long double multiplier, std::vector<long double>& coefficient);
    static void ConvertToBernstein(const std::vector<long double>& power, GeomElType geometry, const std::array<unsigned, 3>& orders, std::vector<long double>& controls);
    static void BuildBernsteinIndices(GeomElType geometry, const std::array<unsigned, 3>& orders, std::vector<std::array<unsigned, 3>>& indices);
    static std::array<std::vector<std::array<unsigned,3>>,11> BuildCachedBernsteinIndices() {
      const GeomElType geometries[11] = {HEX,TET,WEDGE,QUAD,TRI,TRI,QUAD,LINE,LINE,LINE,QUAD};
      const std::array<std::array<unsigned,3>,11> orders = {{{2,2,2},{4,4,4},{3,3,2},{2,2,0},{3,3,0},{4,4,0},{3,2,0},{2,0,0},{3,0,0},{4,0,0},{2,3,0}}};
      std::array<std::vector<std::array<unsigned,3>>,11> result;
      for(unsigned layout = 0; layout < result.size(); ++layout) BuildBernsteinIndices(geometries[layout],orders[layout],result[layout]);
      return result;
    }
    static const std::vector<std::array<unsigned,3>>& CachedBernsteinIndices(GeomElType geometry, const std::array<unsigned,3>& orders) {
      static const std::array<std::vector<std::array<unsigned,3>>,11> indices = BuildCachedBernsteinIndices();
      const unsigned layout = geometry == HEX ? 0U : (geometry == TET ? 1U : (geometry == WEDGE ? 2U : (geometry == QUAD ? (orders[0] == 3 ? 6U : (orders[1] == 3 ? 10U : 3U)) : (geometry == TRI ? (orders[0] == 4 ? 5U : 4U) : orders[0]+5U))));
      return indices[layout];
    }
    static unsigned BernsteinIndex(GeomElType geometry, const std::array<unsigned, 3>& orders, const std::array<unsigned, 3>& index);
    static void AffinePower(const std::vector<long double>& input, unsigned dimension, unsigned degree, unsigned axis, long double origin, long double scale, std::vector<long double>& output);
    static bool CertifySign(GeomElType geometry, const std::array<unsigned,3>& orders, const long double* controls, unsigned count, int sign, long double padding, unsigned depth);
    static bool SplitBernstein(GeomElType geometry, const std::array<unsigned,3>& orders, const long double* controls, unsigned count, std::array<long double,35>& left, std::array<long double,35>& right);
    void RestrictionMap(unsigned axis, long double coordinate, GeomElType& childGeometry, std::array<long double,3>& origins, std::array<long double,3>& scales) const;
    static void RestrictPower(const std::vector<long double>& input, unsigned dimension, unsigned degree, unsigned axis, long double coordinate, const std::array<long double,3>& origins, const std::array<long double,3>& scales, std::vector<long double>& output);
    void RestrictBernstein(unsigned axis, long double coordinate, CutFemLevelSetData& child) const;

    struct PolynomialSource {
      GeomElType geometry = LINE;
      unsigned dimension = 0;
      long double norm = 0;
      mutable std::array<long double,27> power = {};
      mutable std::array<long double,27> tail = {};
      std::array<long double,27> nodal = {};
      std::vector<long double> nodalTail;
      std::array<long double,2> GetNodal(unsigned node) const { return {nodal[node],nodalTail.empty() ? 0 : nodalTail[node]}; }
      unsigned nodes = 0;
      mutable bool compensatedReady = false;
    };

    void PreserveChartSource(CutFemLevelSetData& chart) const;
    void InvalidateAnalysis();
    void EnsureSource() const;
    void EnsureCompensatedSource() const;
    void EnsureCompensatedProjection(unsigned axis) const;
    bool GetLocalFiberCoefficients(unsigned axis, long double coordinate, std::array<std::array<long double,2>,3>& coefficient) const;
    void MapCompensatedSourcePoint(const std::array<long double,3>& coordinate, std::array<std::array<long double,2>,3>& point, std::array<std::array<long double,2>,3>* direction = nullptr, unsigned axis = 0) const;
    void MapCompensatedSourceCoordinates(const std::array<std::array<long double,2>,3>& coordinate, std::array<std::array<long double,2>,3>& point, std::array<std::array<long double,2>,3>* direction, unsigned axis) const;
    void GetCompensatedFiberCoefficients(unsigned axis, const std::array<long double,3>& coordinate, std::array<std::array<long double,2>,3>& coefficient) const;
    void GetCompensatedLineCoefficients(std::array<std::array<long double,2>,3>& coefficient) const;
    std::array<long double,2> CompensatedSourceValue(const std::array<long double,3>& point) const;
    std::array<long double,2> CompensatedSourceValueAt(const std::array<std::array<long double,2>,3>& coordinate) const;
    static std::array<long double,2> CompensatedAdd(const std::array<long double,2>& first, const std::array<long double,2>& second);
    static std::array<long double,2> CompensatedMultiply(const std::array<long double,2>& first, const std::array<long double,2>& second);
    static std::array<long double,2> CompensatedEvaluate(const long double* coefficient, unsigned dimension, unsigned degree, const std::array<std::array<long double,2>,3>& point);
    static void CompensatedPolynomialRoots(const std::array<std::array<long double,2>,5>& coefficient, unsigned degree, long double left, long double right, std::vector<long double>& roots);
    static long double QuadraticDiscriminant(const std::array<std::array<long double,2>,3>& coefficient);
    static unsigned SolveQuadratic(long double a, long double b, long double c, long double discriminant, std::array<long double,2>& roots);
    static const std::vector<std::array<long double,3>>& ReferenceNodes(GeomElType geometry);
    static void ExpandUnit(GeomElType geometry, const std::vector<long double>& nodal, std::vector<long double>& coefficient);
    static std::array<std::array<std::array<long double,27>,27>,6> BuildInterpolationTables();
    static void ExpandChartNodal(GeomElType geometry, const std::vector<long double>& nodal, std::vector<long double>& coefficient);
    void InitializeChart(CutFemLevelSetData& chart) const;
    void ComposeChartMaps(CutFemLevelSetData& chart) const;
    void EnsureBoundaryData() const;
    void PropagateBoundaryData(unsigned axis, long double coordinate, CutFemLevelSetData& child) const;
    int ZeroFiberStatus(unsigned axis) const;
    const GraphCertificate& PrepareGraphCertificate(unsigned axis, bool includeWeak) const;

    mutable std::shared_ptr<const PolynomialSource> _source;
    std::shared_ptr<PolynomialSource> _chartSourceStorage;
    unsigned _restrictionSteps = 0;
    std::array<unsigned,2> _restrictionAxes = {0,0};
    std::array<long double,2> _restrictionCoordinates = {0,0};
    mutable bool _lineAnalysisReady = false;
    mutable std::shared_ptr<LineAnalysis> _lineAnalysis;
    mutable bool _zeroReady = false;
    mutable bool _identicallyZero = false;
    mutable bool _boundaryReady = false;
    mutable unsigned _zeroBoundaryMask = 0;
    mutable bool _dependenceReady = false;
    mutable std::array<bool,3> _depends = {false,false,false};
    mutable std::array<int,3> _parity = {0,0,0};
    mutable std::array<bool,3> _graphCertificateReady = {false,false,false};
    mutable std::array<bool,3> _graphWeakReady = {false,false,false};
    mutable std::shared_ptr<std::array<GraphCertificate,3>> _graphCertificates;
    mutable std::array<bool,3> _zeroFiberReady = {false,false,false};
    mutable std::array<int,3> _zeroFiberStatus = {0,0,0};
    mutable std::array<bool,3> _projectionReady = {false,false,false};
    mutable std::array<std::shared_ptr<FiberProjection>,3> _projections;
    std::vector<long double> _chartNodal;
    mutable std::vector<long double> _componentTerm;
    mutable std::vector<long double> _componentNext;

    bool _initialized = false;
    bool _nativeSpace = false;
    bool _root = false;
    bool _rootGradientAvailable = true;
    mutable bool _classified = false;
    bool _degenerate = false;
    GeomElType _geometry = LINE;
    unsigned _dimension = 0;
    unsigned _rootDimension = 0;
    unsigned _parentDimension = 0;
    unsigned _maximumDepth = 1;
    unsigned _meshLevel = 0;
    std::uint64_t _elementId = 0;
    std::uint64_t _nodalVersion = 0;
    std::uint64_t _revision = 0;
    mutable int _phase = 0;
    mutable std::string _classification = "uninitialized";
    std::vector<double> _nodalValues;
    mutable bool _powerReady = false;
    mutable bool _bernsteinReady = false;
    mutable bool _normalizedReady = false;
    mutable long double _normalization = 1;
    mutable std::vector<long double> _power;
    mutable std::vector<long double> _normalized;
    mutable std::vector<long double> _bernstein;
    mutable long double _bernsteinPadding = 0;
    mutable std::array<long double,2> _bernsteinExtrema = {0,0};
    mutable std::vector<std::array<unsigned, 3>> _bernsteinIndices;
    std::array<unsigned, 3> _bernsteinOrders = {0, 0, 0};
    mutable std::array<bool, 3> _derivativeReady = {false, false, false};
    mutable std::array<bool, 3> _derivativeBoundsReady = {false, false, false};
    mutable std::array<unsigned char,3> _derivativeCornerSigns = {0,0,0};
    mutable bool _graphDirectionReady = false;
    mutable int _graphDirection = -1;
    mutable std::array<std::vector<long double>, 3> _derivatives;
    mutable std::array<bool, 6> _secondDerivativeReady = {false, false, false, false, false, false};
    mutable std::array<std::vector<long double>, 6> _secondDerivatives;
    mutable std::array<std::array<long double, 2>, 3> _derivativeBounds;
    std::array<std::vector<long double>, 3> _rootDerivatives;
    std::vector<long double> _parentOrigin;
    std::vector<std::vector<long double>> _parentMap;
    std::vector<long double> _rootOrigin;
    std::vector<std::vector<long double>> _rootMap;
    long double _parentMeasure = 1;
    mutable std::array<unsigned long long, 4> _workCounts = {0, 0, 0, 0};
};

template<class Real> inline void CutFemLevelSetData::ExpandTensor(GeomElType geometry, const std::vector<Real>& nodal, std::vector<long double>& coefficient) {
  const unsigned dimension = Dimension(geometry);
  const unsigned line[3][3] = {{0, 0, 0}, {2, 0, 0}, {1, 0, 0}};
  const unsigned quad[9][3] = {{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {1, 0, 0}, {2, 1, 0}, {1, 2, 0}, {0, 1, 0}, {1, 1, 0}};
  const unsigned hex[27][3] = {{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, {2, 2, 2}, {0, 2, 2}, {1, 0, 0}, {2, 1, 0}, {1, 2, 0}, {0, 1, 0}, {1, 0, 2}, {2, 1, 2}, {1, 2, 2}, {0, 1, 2}, {0, 0, 1}, {2, 0, 1}, {2, 2, 1}, {0, 2, 1}, {1, 0, 1}, {2, 1, 1}, {1, 2, 1}, {0, 1, 1}, {1, 1, 0}, {1, 1, 2}, {1, 1, 1}};
  const long double interpolation[3][3] = {{1, 0, 0}, {-3, 4, -1}, {2, -4, 2}};
  coefficient.assign(PolynomialSize(dimension, 2), 0);
  for(unsigned n = 0; n < nodal.size(); ++n) {
    const unsigned* index = dimension == 1 ? line[n] : (dimension == 2 ? quad[n] : hex[n]);
    coefficient[PowerIndex({index[0], index[1], index[2]}, dimension, 2)] = nodal[n];
  }
  std::vector<long double> result(coefficient.size());
  for(unsigned axis = 0; axis < dimension; ++axis) {
    for(unsigned n = 0; n < coefficient.size(); ++n) {
      std::array<unsigned, 3> index = PowerExponents(n, dimension, 2);
      const unsigned power = index[axis];
      long double sum = 0;
      long double correction = 0;
      for(unsigned sample = 0; sample < 3; ++sample) {
        index[axis] = sample;
        Accumulate(interpolation[power][sample]*coefficient[PowerIndex(index, dimension, 2)], sum, correction);
      }
      result[n] = sum;
    }
    coefficient.swap(result);
  }
}

inline void CutFemLevelSetData::ExpandTriangle(const std::vector<long double>& nodal, std::vector<long double>& coefficient) {
  const long double bubble = 3 * (nodal[0] + nodal[1] + nodal[2]) - 12 * (nodal[3] + nodal[4] + nodal[5]) + 27 * nodal[6];
  coefficient.assign(9, 0);
  coefficient[0] = nodal[0];
  coefficient[3] = -3 * nodal[0] - nodal[1] + 4 * nodal[3];
  coefficient[1] = -3 * nodal[0] - nodal[2] + 4 * nodal[5];
  coefficient[6] = 2 * nodal[0] + 2 * nodal[1] - 4 * nodal[3];
  coefficient[2] = 2 * nodal[0] + 2 * nodal[2] - 4 * nodal[5];
  coefficient[4] = 4 * nodal[0] - 4 * nodal[3] + 4 * nodal[4] - 4 * nodal[5] + bubble;
  coefficient[7] = coefficient[5] = -bubble;
}

inline void CutFemLevelSetData::AddBarycentricProduct(const std::array<std::array<long double, 4>, 4>& barycentric, const std::vector<unsigned>& factors, long double multiplier, std::vector<long double>& coefficient) {
  if(multiplier == 0) return;
  std::vector<long double> product(27, 0), next(27, 0);
  product[0] = multiplier;
  for(unsigned factor : factors) {
    std::fill(next.begin(), next.end(), 0);
    for(unsigned n = 0; n < product.size(); ++n) {
      if(product[n] == 0) continue;
      const std::array<unsigned, 3> index = PowerExponents(n, 3, 2);
      next[n] += product[n] * barycentric[factor][0];
      for(unsigned axis = 0; axis < 3; ++axis) {
        if(barycentric[factor][axis + 1] == 0) continue;
        if(index[axis] == 2) throw std::logic_error("TET15 basis expansion exceeded its coordinate degree.");
        std::array<unsigned, 3> target = index;
        ++target[axis];
        next[PowerIndex(target, 3, 2)] += product[n] * barycentric[factor][axis + 1];
      }
    }
    product.swap(next);
  }
  for(unsigned n = 0; n < product.size(); ++n) coefficient[n] += product[n];
}

template<class Real> inline void CutFemLevelSetData::ExpandTetrahedron(const std::vector<Real>& nodal, std::vector<long double>& coefficient) {
  const std::array<std::array<long double, 4>, 4> barycentric = {{{1, -1, -1, -1}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}};
  const unsigned edges[6][2] = {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}};
  const unsigned faces[4][3] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
  coefficient.assign(27, 0);
  long double interior = 256 * static_cast<long double>(nodal[14]);
  for(unsigned i = 0; i < 4; ++i) {
    AddBarycentricProduct(barycentric, {i, i}, 2 * static_cast<long double>(nodal[i]), coefficient);
    AddBarycentricProduct(barycentric, {i}, -static_cast<long double>(nodal[i]), coefficient);
    interior -= 4 * static_cast<long double>(nodal[i]);
  }
  for(unsigned edge = 0; edge < 6; ++edge) {
    AddBarycentricProduct(barycentric, {edges[edge][0], edges[edge][1]}, 4 * static_cast<long double>(nodal[4 + edge]), coefficient);
    interior += 32 * static_cast<long double>(nodal[4 + edge]);
  }
  for(unsigned face = 0; face < 4; ++face) {
    long double value = 27 * static_cast<long double>(nodal[10 + face]);
    for(unsigned vertex = 0; vertex < 3; ++vertex) value += 3 * static_cast<long double>(nodal[faces[face][vertex]]);
    for(unsigned edge = 0; edge < 6; ++edge) {
      unsigned contained = 0;
      for(unsigned vertex = 0; vertex < 3; ++vertex) contained += edges[edge][0] == faces[face][vertex] || edges[edge][1] == faces[face][vertex];
      if(contained == 2) value -= 12 * static_cast<long double>(nodal[4 + edge]);
    }
    AddBarycentricProduct(barycentric, {faces[face][0], faces[face][1], faces[face][2]}, value, coefficient);
    interior -= 108 * static_cast<long double>(nodal[10 + face]);
  }
  AddBarycentricProduct(barycentric, {0, 1, 2, 3}, interior, coefficient);
}

inline void CutFemLevelSetData::EnsurePower() const {
  RequireInitialized();
  if(_powerReady) return;
  if(!_root) throw std::logic_error("A slice must receive its coefficients directly from its parent.");
  if(_geometry == LINE || _geometry == QUAD || _geometry == HEX) ExpandTensor(_geometry, _nodalValues, _power);
  else if(_geometry == TRI) ExpandTriangle(std::vector<long double>(_nodalValues.begin(), _nodalValues.end()), _power);
  else if(_geometry == TET) ExpandTetrahedron(_nodalValues, _power);
  else {
    const unsigned nodes[3][7] = {{0, 1, 2, 6, 7, 8, 18}, {12, 13, 14, 15, 16, 17, 20}, {3, 4, 5, 9, 10, 11, 19}};
    const long double interpolation[3][3] = {{1, 0, 0}, {-3, 4, -1}, {2, -4, 2}};
    _power.assign(27, 0);
    for(unsigned plane = 0; plane < 3; ++plane) {
      std::vector<long double> nodal(7), triangle;
      for(unsigned n = 0; n < 7; ++n) nodal[n] = _nodalValues[nodes[plane][n]];
      ExpandTriangle(nodal, triangle);
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) _power[9 * i + 3 * j + k] += interpolation[k][plane] * triangle[3 * i + j];
    }
  }
  long double maximum = 0;
  for(long double value : _power) {
    if(!std::isfinite(value)) throw std::overflow_error("A level-set power coefficient is not finite.");
    maximum = std::max(maximum, std::abs(value));
  }
  if(maximum > 0) {
    int exponent = 0;
    std::frexp(maximum, &exponent);
    _normalization = std::ldexp(1.0L, exponent - 1);
  }
  else _normalization = 1;
  _powerReady = true;
  ++_workCounts[0];
}

inline void CutFemLevelSetData::AffinePower(const std::vector<long double>& input, unsigned dimension, unsigned degree, unsigned axis, long double origin, long double scale, std::vector<long double>& output) {
  if(&input != &output) output = input;
  if(origin == 0 && scale == 1) return;
  const unsigned side = degree+1;
  unsigned stride = 1;
  for(unsigned next = axis+1; next < dimension; ++next) stride *= side;
  std::array<long double,33> scalePowers, originPowers;
  scalePowers[0] = originPowers[0] = 1;
  for(unsigned n = 1; n <= degree; ++n) { scalePowers[n] = scalePowers[n-1]*scale; originPowers[n] = originPowers[n-1]*origin; }
  if(origin == 0) { for(unsigned n = 0; n < input.size(); ++n) output[n] = input[n]*scalePowers[(n/stride)%side]; return; }
  for(unsigned block = 0; block < input.size(); block += stride*side) for(unsigned offset = 0; offset < stride; ++offset) {
    const unsigned base = block+offset;
    std::array<long double,33> source, transformed;
    for(unsigned power = 0; power <= degree; ++power) source[power] = output[base+power*stride];
    for(unsigned power = 0; power <= degree; ++power) { long double sum = 0, correction = 0; for(unsigned old = power; old <= degree; ++old) Accumulate(source[old]*Binomial(old,power)*originPowers[old-power]*scalePowers[power],sum,correction); transformed[power] = sum; }
    for(unsigned power = 0; power <= degree; ++power) output[base+power*stride] = transformed[power];
  }
}

inline void CutFemLevelSetData::BuildBernsteinIndices(GeomElType geometry, const std::array<unsigned, 3>& orders, std::vector<std::array<unsigned, 3>>& indices) {
  indices.clear();
  const unsigned dimension = Dimension(geometry);
  for(unsigned i = 0; i <= orders[0]; ++i) for(unsigned j = 0; j <= (dimension == 1 ? 0U : ((geometry == TRI || geometry == TET || geometry == WEDGE) ? orders[0] - i : orders[1])); ++j) for(unsigned k = 0; k <= (dimension < 3 ? 0U : (geometry == TET ? orders[0] - i - j : orders[2])); ++k) indices.push_back({i, j, k});
}

inline unsigned CutFemLevelSetData::BernsteinIndex(GeomElType geometry, const std::array<unsigned, 3>& orders, const std::array<unsigned, 3>& index) {
  if(geometry == LINE) return index[0];
  if(geometry == TRI || geometry == WEDGE) {
    const unsigned base = index[0] * (2 * orders[0] + 3 - index[0]) / 2 + index[1];
    return geometry == TRI ? base : base * (orders[2] + 1) + index[2];
  }
  if(geometry == TET) return static_cast<unsigned>(Binomial(orders[0] + 3, 3) - Binomial(orders[0] - index[0] + 3, 3)) + index[1] * (2 * (orders[0] - index[0]) + 3 - index[1]) / 2 + index[2];
  return geometry == QUAD ? index[0] * (orders[1] + 1) + index[1] : (index[0] * (orders[1] + 1) + index[1]) * (orders[2] + 1) + index[2];
}

inline void CutFemLevelSetData::ConvertToBernstein(const std::vector<long double>& power, GeomElType geometry, const std::array<unsigned, 3>& orders, std::vector<long double>& controls) {
  const unsigned dimension = Dimension(geometry);
  const std::vector<long double>& unit = power;
  const std::vector<std::array<unsigned,3>>& indices = CachedBernsteinIndices(geometry,orders);
  controls.assign(indices.size(), 0);
  for(unsigned n = 0; n < indices.size(); ++n) {
    long double correction = 0;
    for(unsigned m = 0; m < unit.size(); ++m) {
      if(unit[m] == 0) continue;
      const std::array<unsigned, 3> exponent = PowerExponents(m, dimension, 2);
      long double factor = 1;
      if(geometry == TRI || geometry == TET || geometry == WEDGE) {
        const unsigned baseDimension = geometry == TET ? 3U : 2U;
        unsigned total = 0;
        for(unsigned axis = 0; axis < baseDimension; ++axis) {
          total += exponent[axis];
          factor *= Falling(indices[n][axis], exponent[axis]);
        }
        if(total > orders[0]) throw std::logic_error("Polynomial exceeds the simplex Bernstein order.");
        factor /= Falling(orders[0], total);
        if(geometry == WEDGE) factor *= Binomial(indices[n][2], exponent[2]) / Binomial(orders[2], exponent[2]);
      }
      else for(unsigned axis = 0; axis < dimension; ++axis) factor *= Binomial(indices[n][axis], exponent[axis]) / Binomial(orders[axis], exponent[axis]);
      Accumulate(unit[m]*factor, controls[n], correction);
    }
  }
}

inline void CutFemLevelSetData::EnsureBernstein() const {
  EnsurePower();
  if(_bernsteinReady) return;
  ConvertToBernstein(_power, _geometry, _bernsteinOrders, _bernstein);
  _bernsteinIndices = CachedBernsteinIndices(_geometry,_bernsteinOrders);
  UpdateBernsteinBounds(0);
  _bernsteinReady = true;
  ++_workCounts[1];
}

inline void CutFemLevelSetData::Classify() const {
  RequireInitialized();
  if(_classified) return;
  bool positive = false;
  bool negative = false;
  bool zero = false;
  for(double value : _nodalValues) {
    positive = positive || value > 0;
    negative = negative || value < 0;
    zero = zero || value == 0;
  }
  _phase = 0;
  _classification = "uncertified";
  if(positive && negative) {
    _classification = "crossing";
    _classified = true;
    return;
  }
  if(zero) {
    _classification = "contact";
    _classified = true;
    return;
  }
  EnsureBernstein();
  const long double padding = _bernsteinPadding;
  const long double minimum = _bernsteinExtrema[0];
  const long double maximum = _bernsteinExtrema[1];
  if(minimum > padding) {
    _phase = 1;
    _classification = "positive";
    _classified = true;
    return;
  }
  if(maximum < -padding) {
    _phase = -1;
    _classification = "negative";
    _classified = true;
    return;
  }
  if(_nodalValues.empty()) {
    for(unsigned n = 0; n < _bernsteinIndices.size(); ++n) {
      const std::array<unsigned, 3>& index = _bernsteinIndices[n];
      bool vertex = false;
      if(_geometry == TRI || _geometry == TET || _geometry == WEDGE) {
        const unsigned baseDimension = _geometry == TET ? 3U : 2U;
        unsigned sum = 0;
        unsigned maximumIndex = 0;
        for(unsigned axis = 0; axis < baseDimension; ++axis) {
          sum += index[axis];
          maximumIndex = std::max(maximumIndex, index[axis]);
        }
        vertex = sum == 0 || maximumIndex == _bernsteinOrders[0];
        if(_geometry == WEDGE) vertex = vertex && (index[2] == 0 || index[2] == _bernsteinOrders[2]);
      }
      else {
        vertex = true;
        for(unsigned axis = 0; axis < _dimension; ++axis) vertex = vertex && (index[axis] == 0 || index[axis] == _bernsteinOrders[axis]);
      }
      if(vertex) {
        positive = positive || _bernstein[n] > 0;
        negative = negative || _bernstein[n] < 0;
        zero = zero || _bernstein[n] == 0;
      }
    }
    if(positive && negative) {
      _classification = "crossing";
      _classified = true;
      return;
    }
    if(zero) {
      _classification = "contact";
      _classified = true;
      return;
    }
  }
  const int sign = positive ? 1 : -1;
  if(CertifySign(_geometry, _bernsteinOrders, _bernstein.data(), static_cast<unsigned>(_bernstein.size()), sign, padding, _maximumDepth)) {
    _phase = sign;
    _classification = sign > 0 ? "positive" : "negative";
  }
  _classified = true;
}

inline bool CutFemLevelSetData::CertifySign(GeomElType geometry, const std::array<unsigned,3>& orders, const long double* controls, unsigned count, int sign, long double padding, unsigned depth) {
  bool certified = true;
  for(unsigned n = 0; n < count; ++n) certified = certified && sign*controls[n] > padding;
  if(certified) return true;
  if(depth == 0) return false;
  if(geometry == TRI || geometry == TET || geometry == WEDGE) {
    const unsigned baseDimension = geometry == TET ? 3U : 2U;
    for(unsigned vertex = 0; vertex <= baseDimension; ++vertex) for(unsigned height = 0; height < (geometry == WEDGE ? 2U : 1U); ++height) {
        std::array<unsigned, 3> index = {0, 0, 0};
        if(vertex > 0) index[vertex - 1] = orders[0];
        if(geometry == WEDGE) index[2] = height * orders[2];
        if(sign * controls[BernsteinIndex(geometry, orders, index)] <= padding) return false;
      }
  }
  else {
    const unsigned dimension = Dimension(geometry);
    for(unsigned vertex = 0; vertex < (1U << dimension); ++vertex) {
      std::array<unsigned, 3> index = {0, 0, 0};
      for(unsigned axis = 0; axis < dimension; ++axis) index[axis] = ((vertex >> axis) & 1U) * orders[axis];
      if(sign * controls[BernsteinIndex(geometry, orders, index)] <= padding) return false;
    }
  }
  std::array<long double,35> left, right;
  if(!SplitBernstein(geometry,orders,controls,count,left,right)) return false;
  return CertifySign(geometry,orders,left.data(),count,sign,padding,depth-1) && CertifySign(geometry,orders,right.data(),count,sign,padding,depth-1);
}

inline bool CutFemLevelSetData::SplitBernstein(GeomElType geometry, const std::array<unsigned,3>& orders, const long double* controls, unsigned count, std::array<long double,35>& left, std::array<long double,35>& right) {
  const unsigned dimension = Dimension(geometry);
  const bool simplex = geometry == TRI || geometry == TET || geometry == WEDGE;
  const unsigned baseDimension = geometry == TET ? 3U : 2U;
  static constexpr long double dyad[5] = {1,0.5L,0.25L,0.125L,0.0625L};
  if(!simplex) {
    std::array<unsigned,3> strides = {1,1,1};
    for(unsigned axis = dimension-1; axis > 0; --axis) strides[axis-1] = strides[axis]*(orders[axis]+1);
    long double variation = 0;
    unsigned selectedAxis = 0;
    for(unsigned n = 0; n < count; ++n) for(unsigned axis = 0; axis < dimension; ++axis) { if((n/strides[axis])%(orders[axis]+1) == orders[axis]) continue; const long double difference = std::abs(controls[n]-controls[n+strides[axis]]); if(difference > variation) { variation = difference; selectedAxis = axis; } }
    if(variation == 0) return false;
    const unsigned stride = strides[selectedAxis], degree = orders[selectedAxis];
    for(unsigned n = 0; n < count; ++n) {
      const unsigned index = (n/stride)%(degree+1), base = n-index*stride;
      long double first = 0, second = 0, firstCorrection = 0, secondCorrection = 0;
      for(unsigned k = 0; k <= index; ++k) Accumulate(Binomial(index,k)*dyad[index]*controls[base+k*stride],first,firstCorrection);
      for(unsigned k = 0; k <= degree-index; ++k) Accumulate(Binomial(degree-index,k)*dyad[degree-index]*controls[base+(index+k)*stride],second,secondCorrection);
      left[n] = first;
      right[n] = second;
    }
    return true;
  }
  const std::vector<std::array<unsigned,3>>& indices = CachedBernsteinIndices(geometry,orders);
  long double variation = 0;
  unsigned selectedAxis = 0;
  unsigned selectedFirst = 0;
  unsigned selectedSecond = 1;
  bool triangleSplit = simplex;
  for(unsigned n = 0; n < indices.size(); ++n) {
    const std::array<unsigned, 3>& index = indices[n];
    for(unsigned axis = 0; axis < dimension; ++axis) {
      if(simplex && !(geometry == WEDGE && axis == 2)) continue;
      if(index[axis] >= orders[axis]) continue;
      std::array<unsigned, 3> next = index;
      ++next[axis];
      const long double difference = std::abs(controls[n] - controls[BernsteinIndex(geometry, orders, next)]);
      if(difference > variation) {
        variation = difference;
        selectedAxis = axis;
        triangleSplit = false;
      }
    }
    if(!simplex) continue;
    std::array<unsigned, 4> barycentric = {orders[0], 0, 0, 0};
    for(unsigned axis = 0; axis < baseDimension; ++axis) {
      barycentric[axis + 1] = index[axis];
      barycentric[0] -= index[axis];
    }
    for(unsigned first = 0; first <= baseDimension; ++first) for(unsigned second = first + 1; second <= baseDimension; ++second) {
        if(barycentric[second] == 0) continue;
        std::array<unsigned, 4> shifted = barycentric;
        ++shifted[first];
        --shifted[second];
        std::array<unsigned, 3> next = index;
        for(unsigned axis = 0; axis < baseDimension; ++axis) next[axis] = shifted[axis + 1];
        const long double difference = std::abs(controls[n] - controls[BernsteinIndex(geometry, orders, next)]);
        if(difference > variation) {
          variation = difference;
          selectedFirst = first;
          selectedSecond = second;
          triangleSplit = true;
        }
      }
  }
  if(variation == 0) return false;
  for(unsigned n = 0; n < indices.size(); ++n) {
    long double firstValue = 0, secondValue = 0, firstCorrection = 0, secondCorrection = 0;
    if(!triangleSplit) {
      std::array<unsigned, 3> source = indices[n];
      const unsigned index = source[selectedAxis];
      const unsigned degree = orders[selectedAxis];
      for(unsigned k = 0; k <= index; ++k) {
        source[selectedAxis] = k;
        Accumulate(Binomial(index, k)*dyad[index]*controls[BernsteinIndex(geometry, orders, source)], firstValue, firstCorrection);
      }
      for(unsigned k = 0; k <= degree - index; ++k) {
        source[selectedAxis] = index + k;
        Accumulate(Binomial(degree - index, k)*dyad[degree-index]*controls[BernsteinIndex(geometry, orders, source)], secondValue, secondCorrection);
      }
    }
    else {
      std::array<unsigned, 4> barycentric = {orders[0], 0, 0, 0};
      for(unsigned axis = 0; axis < baseDimension; ++axis) {
        barycentric[axis + 1] = indices[n][axis];
        barycentric[0] -= indices[n][axis];
      }
      const unsigned firstCount = barycentric[selectedFirst];
      const unsigned secondCount = barycentric[selectedSecond];
      for(unsigned k = 0; k <= secondCount; ++k) {
        std::array<unsigned, 4> shifted = barycentric;
        shifted[selectedFirst] += k;
        shifted[selectedSecond] -= k;
        std::array<unsigned, 3> source = indices[n];
        for(unsigned axis = 0; axis < baseDimension; ++axis) source[axis] = shifted[axis + 1];
        Accumulate(Binomial(secondCount, k)*dyad[secondCount]*controls[BernsteinIndex(geometry, orders, source)], firstValue, firstCorrection);
      }
      for(unsigned k = 0; k <= firstCount; ++k) {
        std::array<unsigned, 4> shifted = barycentric;
        shifted[selectedFirst] -= k;
        shifted[selectedSecond] += k;
        std::array<unsigned, 3> source = indices[n];
        for(unsigned axis = 0; axis < baseDimension; ++axis) source[axis] = shifted[axis + 1];
        Accumulate(Binomial(firstCount, k)*dyad[firstCount]*controls[BernsteinIndex(geometry, orders, source)], secondValue, secondCorrection);
      }
    }
    left[n] = firstValue;
    right[n] = secondValue;
  }
  return true;
}

inline void CutFemLevelSetData::RestrictionMap(unsigned axis, long double coordinate, GeomElType& childGeometry, std::array<long double,3>& origins, std::array<long double,3>& scales) const {
  RequireInitialized();
  if(_dimension < 2) throw std::logic_error("A line has no positive-dimensional slice.");
  if(axis >= _dimension) throw std::out_of_range("Slice axis exceeds the active dimension.");
  if(!std::isfinite(coordinate)) throw std::invalid_argument("Slice coordinate must be finite.");
  const bool simplexAxis = _geometry == TRI || _geometry == TET || (_geometry == WEDGE && axis < 2);
  if(coordinate < 0 || coordinate > 1) throw std::out_of_range("Slice coordinate lies outside the reference element.");
  childGeometry = _geometry == HEX ? QUAD : (_geometry == TET ? TRI : (_geometry == WEDGE ? (axis == 2 ? TRI : QUAD) : LINE));
  origins.fill(0);
  scales.fill(1);
  const long double width = 1 - coordinate;
  if(_geometry == TET) std::fill(scales.begin(), scales.end(), width);
  if(_geometry == TRI || (_geometry == WEDGE && axis < 2)) scales[0] = width;
}

inline void CutFemLevelSetData::RestrictPower(const std::vector<long double>& input, unsigned dimension, unsigned degree, unsigned axis, long double coordinate, const std::array<long double,3>& origins, const std::array<long double,3>& scales, std::vector<long double>& output) {
  const unsigned side = degree+1;
  unsigned stride = 1;
  for(unsigned next = axis+1; next < dimension; ++next) stride *= side;
  output.resize(PolynomialSize(dimension - 1, degree));
  for(unsigned n = 0; n < output.size(); ++n) {
    const unsigned base = (n/stride)*stride*side+n%stride;
    if(degree == 2) { output[n] = EvaluateQuadratic(input.data()+base,stride,coordinate); continue; }
    long double value = input[base+degree*stride];
    for(unsigned power = degree; power > 0; --power) {
      value = std::fma(value,coordinate,input[base+(power-1)*stride]);
    }
    output[n] = value;
  }
  for(unsigned childAxis = 0; childAxis < dimension - 1; ++childAxis) AffinePower(output, dimension - 1, degree, childAxis, origins[childAxis], scales[childAxis], output);
}

inline std::array<CutFemLevelSetData::BernsteinRestrictionPlan,19> CutFemLevelSetData::BuildRestrictionPlans() {
  const GeomElType geometries[8] = {HEX,TET,WEDGE,QUAD,TRI,TRI,QUAD,QUAD};
  const std::array<std::array<unsigned,3>,8> orders = {{{2,2,2},{4,4,4},{3,3,2},{2,2,0},{3,3,0},{4,4,0},{3,2,0},{2,3,0}}};
  std::array<BernsteinRestrictionPlan,19> result;
  unsigned nextPlan = 0;
  for(unsigned layout = 0; layout < 8; ++layout) {
    const GeomElType geometry = geometries[layout];
    const unsigned dimension = Dimension(geometry);
    for(unsigned axis = 0; axis < dimension; ++axis) {
      BernsteinRestrictionPlan& plan = result[nextPlan++];
      plan.simplex = geometry == TRI || geometry == TET || (geometry == WEDGE && axis < 2);
      plan.degree = orders[layout][axis];
      const GeomElType childGeometry = geometry == HEX ? QUAD : (geometry == TET ? TRI : (geometry == WEDGE ? (axis == 2 ? TRI : QUAD) : LINE));
      if(plan.simplex) {
        if(geometry == TET) plan.orders = {plan.degree,plan.degree,0};
        else if(geometry == WEDGE) plan.orders = {plan.degree,orders[layout][2],0};
        else plan.orders = {plan.degree,0,0};
      }
      else { unsigned next = 0; for(unsigned parentAxis = 0; parentAxis < dimension; ++parentAxis) if(parentAxis != axis) plan.orders[next++] = orders[layout][parentAxis]; }
      BuildBernsteinIndices(childGeometry,plan.orders,plan.indices);
      plan.offsets.push_back(0);
      for(const std::array<unsigned,3>& index : plan.indices) {
        if(!plan.simplex) {
          std::array<unsigned,3> source = {0,0,0};
          unsigned next = 0;
          for(unsigned parentAxis = 0; parentAxis < dimension; ++parentAxis) if(parentAxis != axis) source[parentAxis] = index[next++];
          for(unsigned k = 0; k <= plan.degree; ++k) { source[axis] = k; plan.terms.push_back({BernsteinIndex(geometry,orders[layout],source),k,1}); }
        }
        else {
          const unsigned baseDimension = geometry == TET ? 3U : 2U;
          const unsigned first = index[0], second = geometry == TET ? index[1] : 0U, zeroth = plan.degree-first-second;
          for(unsigned a = 0; a <= zeroth; ++a) for(unsigned b = 0; b <= first; ++b) for(unsigned c = 0; c <= second; ++c) {
            const unsigned total = a+b+c;
            std::array<unsigned,3> source = {0,0,0};
            source[axis] = plan.degree-total;
            unsigned next = 0;
            for(unsigned parentAxis = 0; parentAxis < baseDimension; ++parentAxis) if(parentAxis != axis) source[parentAxis] = next++ == 0 ? b : c;
            if(geometry == WEDGE) source[2] = index[1];
            const unsigned multiplier = static_cast<unsigned>(Binomial(zeroth,a)*Binomial(first,b)*Binomial(second,c));
            plan.terms.push_back({BernsteinIndex(geometry,orders[layout],source),plan.degree-total,multiplier});
          }
        }
        plan.offsets.push_back(static_cast<unsigned>(plan.terms.size()));
      }
    }
  }
  return result;
}

inline void CutFemLevelSetData::RestrictBernstein(unsigned axis, long double coordinate, CutFemLevelSetData& child) const {
  const BernsteinRestrictionPlan& plan = GetRestrictionPlan(axis);
  if(child._bernsteinOrders != plan.orders || child._bernsteinIndices.size() != plan.indices.size()) { child._bernsteinOrders = plan.orders; child._bernsteinIndices = plan.indices; }
  child._bernstein.resize(plan.indices.size());
  std::array<long double,5> weights, fixedPowers, widthPowers;
  if(plan.simplex) {
    fixedPowers[0] = widthPowers[0] = 1;
    for(unsigned k = 1; k <= plan.degree; ++k) { fixedPowers[k] = fixedPowers[k-1]*coordinate; widthPowers[k] = widthPowers[k-1]*(1-coordinate); }
  }
  else { const long double unit = coordinate; for(unsigned k = 0; k <= plan.degree; ++k) weights[k] = Bernstein(plan.degree,k,unit); }
  for(unsigned n = 0; n < child._bernstein.size(); ++n) {
    long double value = 0, correction = 0;
    for(unsigned term = plan.offsets[n]; term < plan.offsets[n+1]; ++term) {
      const std::array<unsigned,3>& entry = plan.terms[term];
      const long double factor = plan.simplex ? entry[2]*fixedPowers[entry[1]]*widthPowers[plan.degree-entry[1]] : weights[entry[1]];
      Accumulate(factor*_bernstein[entry[0]],value,correction);
    }
    child._bernstein[n] = value;
  }
  child.UpdateBernsteinBounds(_bernsteinPadding);
  child._bernsteinReady = true;
}

inline void CutFemLevelSetData::Restrict(unsigned axis, long double coordinate, CutFemLevelSetData& child, bool retainRootGradient) const {
  if(&child == this) throw std::invalid_argument("Restriction requires distinct parent and child data.");
  GeomElType childGeometry;
  std::array<long double,3> origins, scales;
  RestrictionMap(axis, coordinate, childGeometry, origins, scales);
  EnsurePower();
  EnsureBernstein();
  if(retainRootGradient) for(unsigned rootAxis = 0; rootAxis < _rootDimension; ++rootAxis) GetRootDerivative(rootAxis);
  child._initialized = true;
  child._root = false;
  child._nativeSpace = _nativeSpace;
  child._rootGradientAvailable = retainRootGradient;
  child._geometry = childGeometry;
  child._dimension = _dimension - 1;
  child._rootDimension = _rootDimension;
  child._parentDimension = _dimension;
  child._maximumDepth = _maximumDepth;
  child._meshLevel = _meshLevel;
  child._elementId = _elementId;
  child._nodalVersion = _nodalVersion;
  ++child._revision;
  child._nodalValues.clear();
  child._classified = false;
  child._normalizedReady = false;
  child._normalized.clear();
  child._normalization = _normalization;
  child._derivativeReady.fill(false);
  child._secondDerivativeReady.fill(false);
  child._derivativeBoundsReady.fill(false);
  child._graphDirectionReady = false;
  child.InvalidateAnalysis();
  EnsureSource();
  child._source = _source;
  child._restrictionSteps = _restrictionSteps+1;
  child._restrictionAxes = _restrictionAxes;
  child._restrictionCoordinates = _restrictionCoordinates;
  child._restrictionAxes[_restrictionSteps] = axis;
  child._restrictionCoordinates[_restrictionSteps] = coordinate;
  RestrictPower(_power, _dimension, 2, axis, coordinate, origins, scales, child._power);
  child._powerReady = true;
  if(retainRootGradient) for(unsigned rootAxis = 0; rootAxis < _rootDimension; ++rootAxis) RestrictPower(GetRootDerivative(rootAxis), _dimension, 2, axis, coordinate, origins, scales, child._rootDerivatives[rootAxis]);
  else for(std::vector<long double>& derivative : child._rootDerivatives) derivative.clear();
  child._parentOrigin.assign(_dimension, 0);
  child._parentOrigin[axis] = coordinate;
  child._parentMap.resize(_dimension);
  for(std::vector<long double>& row : child._parentMap) row.assign(child._dimension, 0);
  child._parentMeasure = 1;
  unsigned next = 0;
  for(unsigned parentAxis = 0; parentAxis < _dimension; ++parentAxis) if(parentAxis != axis) {
      child._parentOrigin[parentAxis] = origins[next];
      child._parentMap[parentAxis][next] = scales[next];
      child._parentMeasure *= std::abs(scales[next]);
      ++next;
    }
  child._degenerate = _degenerate || child._parentMeasure == 0;
  child._rootOrigin.resize(_rootDimension);
  for(unsigned rootAxis = 0; rootAxis < _rootDimension; ++rootAxis) { child._rootOrigin[rootAxis] = _rootOrigin[rootAxis]; for(unsigned parentAxis = 0; parentAxis < _dimension; ++parentAxis) child._rootOrigin[rootAxis] += _rootMap[rootAxis][parentAxis]*child._parentOrigin[parentAxis]; }
  child._rootMap.resize(_rootDimension);
  for(unsigned rootAxis = 0; rootAxis < _rootDimension; ++rootAxis) {
    child._rootMap[rootAxis].assign(child._dimension, 0);
    unsigned local = 0;
    for(unsigned parentAxis = 0; parentAxis < _dimension; ++parentAxis) if(parentAxis != axis) { child._rootMap[rootAxis][local] = _rootMap[rootAxis][parentAxis]*scales[local]; ++local; }
  }
  RestrictBernstein(axis, coordinate, child);
  PropagateBoundaryData(axis,coordinate,child);
  ++_workCounts[3];
}

inline std::vector<long double> CutFemLevelSetData::GetLineRoots() const {
  std::array<long double,2> roots;
  const unsigned count = GetLineRoots(roots);
  return std::vector<long double>(roots.begin(),roots.begin()+count);
}

inline unsigned CutFemLevelSetData::GetLineRoots(std::array<long double,2>& roots) const {
  const LineAnalysis& analysis = GetLineAnalysis();
  if(analysis.identicallyZero) throw std::domain_error("An identically zero line has no isolated root set.");
  roots.fill(0);
  for(unsigned n = 0; n < analysis.rootCount; ++n) roots[n] = analysis.roots[n].coordinate;
  return analysis.rootCount;
}

#include "CutFemLevelSetDataAnalysis.hpp"
#include "CutFemLevelSetDataCharts.hpp"

#endif
