#ifndef FEMUS_EX47_CUT_FEM_LEVEL_SET_DATA_HPP
#define FEMUS_EX47_CUT_FEM_LEVEL_SET_DATA_HPP

#include "GeomElTypeEnum.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class CutFemLevelSetData {
public:
  // Nodal inputs use FEMuS ordering and reference coordinates; the default allows one sign-certification bisection.
  CutFemLevelSetData(GeomElType geometry, const std::vector<double>& nodalValues, unsigned maximumDepth = 1) { Reset(geometry,nodalValues,maximumDepth); }

  // An empty instance is reusable output storage for Restrict().
  CutFemLevelSetData() = default;

  void Reset(GeomElType geometry, const std::vector<double>& nodalValues, unsigned maximumDepth = 1) {
    ValidateGeometry(geometry);
    ValidateDepth(maximumDepth);
    if(nodalValues.size() != NodeCount(geometry)) throw std::invalid_argument("Incorrect FEMuS level-set node count.");
    for(double value : nodalValues) if(!std::isfinite(value)) throw std::invalid_argument("Level-set nodal values must be finite.");
    if(_initialized && _root && geometry == _geometry && nodalValues == _nodalValues) { SetMaximumDepth(maximumDepth); return; }
    CutFemLevelSetData replacement;
    replacement._geometry = geometry;
    replacement._dimension = Dimension(geometry);
    replacement._rootDimension = replacement._dimension;
    replacement._parentDimension = replacement._dimension;
    replacement._maximumDepth = maximumDepth;
    replacement._nodalValues = nodalValues;
    replacement._initialized = true;
    replacement._root = true;
    replacement._meshLevel = _meshLevel;
    replacement._elementId = _elementId;
    replacement._nodalVersion = _nodalVersion+1;
    replacement._revision = _revision+1;
    replacement._parentOrigin.assign(replacement._dimension,0);
    replacement._rootOrigin.assign(replacement._dimension,0);
    replacement._parentMap.assign(replacement._dimension,std::vector<long double>(replacement._dimension,0));
    replacement._rootMap = replacement._parentMap;
    for(unsigned axis = 0; axis < replacement._dimension; ++axis) { replacement._parentMap[axis][axis] = 1; replacement._rootMap[axis][axis] = 1; }
    replacement.SetBernsteinOrders();
    replacement.Classify();
    for(unsigned i = 0; i < _workCounts.size(); ++i) replacement._workCounts[i] += _workCounts[i];
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

  void SetContext(unsigned meshLevel, std::uint64_t elementId, std::uint64_t nodalVersion) { RequireInitialized(); _meshLevel = meshLevel; _elementId = elementId; _nodalVersion = nodalVersion; }
  bool IsInitialized() const { return _initialized; }
  GeomElType GetGeometry() const { RequireInitialized(); return _geometry; }
  unsigned GetDimension() const { RequireInitialized(); return _dimension; }
  unsigned GetRootDimension() const { RequireInitialized(); return _rootDimension; }
  unsigned GetMaximumDepth() const { RequireInitialized(); return _maximumDepth; }
  unsigned GetMeshLevel() const { RequireInitialized(); return _meshLevel; }
  std::uint64_t GetElementId() const { RequireInitialized(); return _elementId; }
  std::uint64_t GetNodalVersion() const { RequireInitialized(); return _nodalVersion; }
  std::uint64_t GetRevision() const { RequireInitialized(); return _revision; }
  int GetPhase() const { RequireInitialized(); return _phase; }
  bool IsCutCandidate() const { return GetPhase() == 0; }
  bool IsDegenerate() const { RequireInitialized(); return _degenerate; }
  std::string GetClassification() const { RequireInitialized(); return _classification; }
  bool HasPowerCoefficients() const { return _powerReady; }
  bool HasBernsteinCoefficients() const { return _bernsteinReady; }
  const std::vector<double>& GetNodalValues() const { RequireInitialized(); return _nodalValues; }

  // Power coefficient index is i*3^(dimension-1)+j*3^(dimension-2)+k for the active coordinates.
  const std::vector<long double>& GetPowerCoefficients() const { EnsurePower(); return _power; }
  const std::vector<long double>& GetNormalizedPowerCoefficients() const { EnsurePower(); if(!_normalizedReady) { _normalized.resize(_power.size()); for(unsigned i = 0; i < _power.size(); ++i) _normalized[i] = _power[i]/_normalization; _normalizedReady = true; } return _normalized; }
  long double GetNormalizationScale() const { EnsurePower(); return _normalization; }
  const std::vector<long double>& GetBernsteinCoefficients() const { EnsureBernstein(); return _bernstein; }
  std::array<unsigned,3> GetBernsteinOrders() const { RequireInitialized(); return _bernsteinOrders; }
  const std::vector<std::array<unsigned,3>>& GetBernsteinIndices() const { EnsureBernstein(); return _bernsteinIndices; }

  // Counts are nodal expansions, power-to-Bernstein conversions, local derivatives, and restrictions performed.
  std::array<unsigned long long,4> GetWorkCounts() const { return _workCounts; }

  long double Evaluate(const std::vector<long double>& point) const { EnsurePower(); ValidatePoint(point); return EvaluatePolynomial(_power,_dimension,2,point); }

  long double EvaluateBernstein(const std::vector<long double>& point) const {
    EnsureBernstein();
    ValidatePoint(point);
    std::vector<long double> unit = point;
    if(_geometry == LINE || _geometry == QUAD || _geometry == HEX) for(long double& coordinate : unit) coordinate = (coordinate+1)/2;
    if(_geometry == WEDGE) unit[2] = (unit[2]+1)/2;
    long double result = 0;
    long double correction = 0;
    for(unsigned n = 0; n < _bernstein.size(); ++n) {
      const std::array<unsigned,3>& index = _bernsteinIndices[n];
      long double value = _bernstein[n];
      if(_geometry == TRI || _geometry == TET || _geometry == WEDGE) {
        const unsigned baseDimension = _geometry == TET ? 3U : 2U;
        const unsigned degree = _bernsteinOrders[0];
        unsigned sum = 0;
        long double remaining = 1;
        for(unsigned axis = 0; axis < baseDimension; ++axis) { sum += index[axis]; remaining -= unit[axis]; }
        value *= Multinomial(degree,index,baseDimension)*Power(remaining,degree-sum);
        for(unsigned axis = 0; axis < baseDimension; ++axis) value *= Power(unit[axis],index[axis]);
        if(_geometry == WEDGE) value *= Bernstein(_bernsteinOrders[2],index[2],unit[2]);
      }
      else for(unsigned axis = 0; axis < _dimension; ++axis) value *= Bernstein(_bernsteinOrders[axis],index[axis],unit[axis]);
      Accumulate(value,result,correction);
    }
    return result;
  }

  std::array<long double,2> GetBounds() const { EnsureBernstein(); return {*std::min_element(_bernstein.begin(),_bernstein.end())-_bernsteinPadding,*std::max_element(_bernstein.begin(),_bernstein.end())+_bernsteinPadding}; }

  const std::vector<long double>& GetDerivativeCoefficients(unsigned axis) const {
    EnsurePower();
    if(axis >= _dimension) throw std::out_of_range("Derivative axis exceeds the active dimension.");
    if(!_derivativeReady[axis]) { Differentiate(_power,_dimension,2,axis,_derivatives[axis]); _derivativeReady[axis] = true; ++_workCounts[2]; }
    return _derivatives[axis];
  }

  long double EvaluateDerivative(unsigned axis, const std::vector<long double>& point) const { ValidatePoint(point); return EvaluatePolynomial(GetDerivativeCoefficients(axis),_dimension,2,point); }
  long double EvaluateSecondDerivative(unsigned first, unsigned second, const std::vector<long double>& point) const {
    ValidatePoint(point);
    if(first >= _dimension || second >= _dimension) throw std::out_of_range("Derivative axis exceeds the active dimension.");
    if(first > second) std::swap(first,second);
    const unsigned index = first*(7-first)/2+second-first;
    if(!_secondDerivativeReady[index]) { Differentiate(GetDerivativeCoefficients(first),_dimension,2,second,_secondDerivatives[index]); _secondDerivativeReady[index] = true; ++_workCounts[2]; }
    return EvaluatePolynomial(_secondDerivatives[index],_dimension,2,point);
  }

  std::vector<long double> GetGradient(const std::vector<long double>& point) const { ValidatePoint(point); std::vector<long double> result(_dimension); for(unsigned axis = 0; axis < _dimension; ++axis) result[axis] = EvaluateDerivative(axis,point); return result; }

  // Root gradients retain every original reference-coordinate component after dimension reduction.
  std::vector<long double> GetRootGradient(const std::vector<long double>& point) const { ValidatePoint(point); std::vector<long double> result(_rootDimension); for(unsigned axis = 0; axis < _rootDimension; ++axis) result[axis] = EvaluatePolynomial(GetRootDerivative(axis),_dimension,2,point); return result; }

  std::array<long double,2> GetDerivativeBounds(unsigned axis) const {
    const std::vector<long double>& derivative = GetDerivativeCoefficients(axis);
    if(!_derivativeBoundsReady[axis]) { EnsureBernstein(); std::vector<long double> controls; ConvertToBernstein(derivative,_geometry,_bernsteinOrders,controls); const long double padding = Padding(controls)+2*_bernsteinOrders[axis]*_bernsteinPadding; _derivativeBounds[axis] = {*std::min_element(controls.begin(),controls.end())-padding,*std::max_element(controls.begin(),controls.end())+padding}; _derivativeBoundsReady[axis] = true; }
    return _derivativeBounds[axis];
  }

  int GetGraphDirection() const {
    RequireInitialized();
    int result = -1;
    long double best = -1;
    for(unsigned axis = 0; axis < _dimension; ++axis) { const std::array<long double,2> bounds = GetDerivativeBounds(axis); const long double lower = bounds[0] > 0 ? bounds[0] : (bounds[1] < 0 ? -bounds[1] : 0); const long double upper = std::max(std::abs(bounds[0]),std::abs(bounds[1])); const long double score = upper > 0 ? lower/upper : 0; if(lower > 0 && score > best) { best = score; result = static_cast<int>(axis); } }
    return result;
  }

  const std::vector<long double>& GetParentOrigin() const { RequireInitialized(); return _parentOrigin; }
  const std::vector<std::vector<long double>>& GetParentMap() const { RequireInitialized(); return _parentMap; }
  long double GetParentMeasureScale() const { RequireInitialized(); return _parentMeasure; }
  std::vector<long double> MapToParent(const std::vector<long double>& point) const { ValidatePoint(point); return ApplyMap(_parentOrigin,_parentMap,point); }
  std::vector<long double> MapToRoot(const std::vector<long double>& point) const { ValidatePoint(point); return ApplyMap(_rootOrigin,_rootMap,point); }

  // Restriction fixes an active reference coordinate and fills an existing lower-dimensional data object.
  void Restrict(unsigned axis, long double coordinate, CutFemLevelSetData& child) const;

  // External polynomial components keep their component count and their full supplied degree under restriction.
  void RestrictComponents(unsigned axis, long double coordinate, unsigned degree, const std::vector<std::vector<long double>>& components, std::vector<std::vector<long double>>& output) const {
    std::vector<long double> origin;
    std::vector<long double> scales;
    GeomElType childGeometry;
    RestrictionMap(axis,coordinate,childGeometry,origin,scales);
    if(&components == &output) throw std::invalid_argument("Restriction input and output components must be distinct.");
    output.resize(components.size());
    for(unsigned component = 0; component < components.size(); ++component) { ValidatePolynomial(components[component],_dimension,degree); RestrictPower(components[component],_dimension,degree,axis,coordinate,origin,scales,output[component]); }
  }

  std::vector<long double> GetLineRoots() const;

  static long double EvaluatePolynomial(const std::vector<long double>& coefficient, unsigned dimension, unsigned degree, const std::vector<long double>& point) {
    if(dimension < 1 || dimension > 3 || point.size() != dimension || coefficient.size() != PolynomialSize(dimension,degree)) throw std::invalid_argument("Invalid polynomial evaluation layout.");
    for(long double coordinate : point) if(!std::isfinite(coordinate)) throw std::invalid_argument("Polynomial evaluation coordinates must be finite.");
    const unsigned side = degree+1;
    long double result = 0;
    if(dimension == 1) { result = coefficient[degree]; for(unsigned i = degree; i > 0; --i) result = result*point[0]+coefficient[i-1]; return result; }
    for(unsigned i = side; i > 0; --i) {
      long double row = 0;
      for(unsigned j = side; j > 0; --j) {
        long double line;
        if(dimension == 2) line = coefficient[(i-1)*side+j-1];
        else { line = coefficient[((i-1)*side+j-1)*side+degree]; for(unsigned k = degree; k > 0; --k) line = line*point[2]+coefficient[((i-1)*side+j-1)*side+k-1]; }
        row = j == side ? line : row*point[1]+line;
      }
      result = i == side ? row : result*point[0]+row;
    }
    return result;
  }

private:
  static unsigned Dimension(GeomElType geometry) { return geometry == LINE ? 1U : (geometry == QUAD || geometry == TRI ? 2U : 3U); }
  static unsigned NodeCount(GeomElType geometry) { return geometry == LINE ? 3U : (geometry == QUAD ? 9U : (geometry == TRI ? 7U : (geometry == HEX ? 27U : (geometry == TET ? 15U : 21U)))); }
  static void ValidateGeometry(GeomElType geometry) { if(geometry != LINE && geometry != QUAD && geometry != TRI && geometry != HEX && geometry != TET && geometry != WEDGE) throw std::invalid_argument("Unsupported level-set element."); }
  static void ValidateDepth(unsigned depth) { if(depth > 32) throw std::invalid_argument("Maximum detection depth must lie between 0 and 32."); }
  void RequireInitialized() const { if(!_initialized) throw std::logic_error("The level-set data object is not initialized."); }
  static unsigned PolynomialSize(unsigned dimension, unsigned degree) { if(dimension < 1 || dimension > 3 || degree > 32) throw std::invalid_argument("Invalid polynomial dimensions or degree."); unsigned result = 1; for(unsigned axis = 0; axis < dimension; ++axis) result *= degree+1; return result; }
  static void ValidatePolynomial(const std::vector<long double>& coefficient, unsigned dimension, unsigned degree) { if(coefficient.size() != PolynomialSize(dimension,degree)) throw std::invalid_argument("Incorrect polynomial coefficient count."); for(long double value : coefficient) if(!std::isfinite(value)) throw std::invalid_argument("Polynomial coefficients must be finite."); }
  static long double Power(long double value, unsigned exponent) { long double result = 1; for(unsigned i = 0; i < exponent; ++i) result *= value; return result; }
  static long double Binomial(unsigned n, unsigned k) { if(k > n) return 0; long double result = 1; for(unsigned i = 1; i <= k; ++i) result *= static_cast<long double>(n-i+1)/i; return result; }
  static long double Falling(unsigned n, unsigned k) { if(k > n) return 0; long double result = 1; for(unsigned i = 0; i < k; ++i) result *= n-i; return result; }
  static long double Bernstein(unsigned degree, unsigned index, long double coordinate) { return Binomial(degree,index)*Power(coordinate,index)*Power(1-coordinate,degree-index); }
  static long double Multinomial(unsigned degree, const std::array<unsigned,3>& index, unsigned dimension) { long double result = 1; for(unsigned axis = 0; axis < dimension; ++axis) { result *= Binomial(degree,index[axis]); degree -= index[axis]; } return result; }
  static void Accumulate(long double term, long double& sum, long double& correction) { const long double adjusted = term-correction; const long double next = sum+adjusted; correction = (next-sum)-adjusted; sum = next; }
  static long double Padding(const std::vector<long double>& coefficient) { long double norm = 0; for(long double value : coefficient) norm += std::abs(value); return 8192*std::numeric_limits<long double>::epsilon()*norm; }
  static unsigned PowerIndex(const std::array<unsigned,3>& index, unsigned dimension, unsigned degree) { unsigned result = index[0]; for(unsigned axis = 1; axis < dimension; ++axis) result = result*(degree+1)+index[axis]; return result; }
  static std::array<unsigned,3> PowerExponents(unsigned index, unsigned dimension, unsigned degree) { std::array<unsigned,3> result = {0,0,0}; for(unsigned axis = dimension; axis > 0; --axis) { result[axis-1] = index%(degree+1); index /= degree+1; } return result; }
  static std::vector<long double> ApplyMap(const std::vector<long double>& origin, const std::vector<std::vector<long double>>& map, const std::vector<long double>& point) { std::vector<long double> result = origin; for(unsigned row = 0; row < map.size(); ++row) for(unsigned column = 0; column < point.size(); ++column) result[row] += map[row][column]*point[column]; return result; }

  void ValidatePoint(const std::vector<long double>& point) const {
    RequireInitialized();
    if(point.size() != _dimension) throw std::invalid_argument("Point dimension does not match the level-set data.");
    for(long double coordinate : point) if(!std::isfinite(coordinate)) throw std::invalid_argument("Point coordinates must be finite.");
    long double sum = 0;
    for(unsigned axis = 0; axis < _dimension; ++axis) {
      const bool simplexAxis = _geometry == TRI || _geometry == TET || (_geometry == WEDGE && axis < 2);
      const long double lower = simplexAxis ? 0 : -1;
      const long double padding = 64*std::numeric_limits<long double>::epsilon();
      if(point[axis] < lower-padding || point[axis] > 1+padding) throw std::out_of_range("Point lies outside the reference element.");
      if(simplexAxis) sum += point[axis];
    }
    if((_geometry == TRI || _geometry == TET || _geometry == WEDGE) && sum > 1+128*std::numeric_limits<long double>::epsilon()) throw std::out_of_range("Point lies outside the reference simplex.");
  }

  static void Differentiate(const std::vector<long double>& coefficient, unsigned dimension, unsigned degree, unsigned axis, std::vector<long double>& result) {
    result.assign(coefficient.size(),0);
    for(unsigned n = 0; n < coefficient.size(); ++n) { std::array<unsigned,3> index = PowerExponents(n,dimension,degree); if(index[axis] == 0) continue; const long double value = index[axis]*coefficient[n]; --index[axis]; result[PowerIndex(index,dimension,degree)] = value; }
  }

  void SetBernsteinOrders() { _bernsteinOrders = {0,0,0}; for(unsigned axis = 0; axis < _dimension; ++axis) _bernsteinOrders[axis] = _geometry == TET ? 4U : ((_geometry == TRI || (_geometry == WEDGE && axis < 2)) ? 3U : 2U); }
  void EnsurePower() const;
  void EnsureBernstein() const;
  void Classify();
  const std::vector<long double>& GetRootDerivative(unsigned axis) const { if(axis >= _rootDimension) throw std::out_of_range("Root derivative axis exceeds the original dimension."); if(_root) return GetDerivativeCoefficients(axis); return _rootDerivatives[axis]; }
  static void ExpandTensor(GeomElType geometry, const std::vector<double>& nodal, std::vector<long double>& coefficient);
  static void ExpandTriangle(const std::vector<long double>& nodal, std::vector<long double>& coefficient);
  static void ExpandTetrahedron(const std::vector<double>& nodal, std::vector<long double>& coefficient);
  static void AddBarycentricProduct(const std::array<std::array<long double,4>,4>& barycentric, const std::vector<unsigned>& factors, long double multiplier, std::vector<long double>& coefficient);
  static void ConvertToBernstein(const std::vector<long double>& power, GeomElType geometry, const std::array<unsigned,3>& orders, std::vector<long double>& controls);
  static void BuildBernsteinIndices(GeomElType geometry, const std::array<unsigned,3>& orders, std::vector<std::array<unsigned,3>>& indices);
  static unsigned BernsteinIndex(GeomElType geometry, const std::array<unsigned,3>& orders, const std::array<unsigned,3>& index);
  static void AffinePower(const std::vector<long double>& input, unsigned dimension, unsigned degree, unsigned axis, long double origin, long double scale, std::vector<long double>& output);
  static bool CertifySign(GeomElType geometry, const std::array<unsigned,3>& orders, const std::vector<long double>& controls, int sign, long double padding, unsigned depth);
  static void SplitBernstein(GeomElType geometry, const std::array<unsigned,3>& orders, const std::vector<long double>& controls, std::vector<long double>& left, std::vector<long double>& right);
  void RestrictionMap(unsigned axis, long double coordinate, GeomElType& childGeometry, std::vector<long double>& origins, std::vector<long double>& scales) const;
  static void RestrictPower(const std::vector<long double>& input, unsigned dimension, unsigned degree, unsigned axis, long double coordinate, const std::vector<long double>& origins, const std::vector<long double>& scales, std::vector<long double>& output);
  void RestrictBernstein(unsigned axis, long double coordinate, CutFemLevelSetData& child) const;

  bool _initialized = false;
  bool _root = false;
  bool _classified = false;
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
  int _phase = 0;
  std::string _classification = "uninitialized";
  std::vector<double> _nodalValues;
  mutable bool _powerReady = false;
  mutable bool _bernsteinReady = false;
  mutable bool _normalizedReady = false;
  mutable long double _normalization = 1;
  mutable std::vector<long double> _power;
  mutable std::vector<long double> _normalized;
  mutable std::vector<long double> _bernstein;
  mutable long double _bernsteinPadding = 0;
  mutable std::vector<std::array<unsigned,3>> _bernsteinIndices;
  std::array<unsigned,3> _bernsteinOrders = {0,0,0};
  mutable std::array<bool,3> _derivativeReady = {false,false,false};
  mutable std::array<bool,3> _derivativeBoundsReady = {false,false,false};
  mutable std::array<std::vector<long double>,3> _derivatives;
  mutable std::array<bool,6> _secondDerivativeReady = {false,false,false,false,false,false};
  mutable std::array<std::vector<long double>,6> _secondDerivatives;
  mutable std::array<std::array<long double,2>,3> _derivativeBounds;
  std::array<std::vector<long double>,3> _rootDerivatives;
  std::vector<long double> _parentOrigin;
  std::vector<std::vector<long double>> _parentMap;
  std::vector<long double> _rootOrigin;
  std::vector<std::vector<long double>> _rootMap;
  long double _parentMeasure = 1;
  mutable std::array<unsigned long long,4> _workCounts = {0,0,0,0};
};

inline void CutFemLevelSetData::ExpandTensor(GeomElType geometry, const std::vector<double>& nodal, std::vector<long double>& coefficient) {
  const unsigned dimension = Dimension(geometry);
  const unsigned line[3][3] = {{0,0,0},{2,0,0},{1,0,0}};
  const unsigned quad[9][3] = {{0,0,0},{2,0,0},{2,2,0},{0,2,0},{1,0,0},{2,1,0},{1,2,0},{0,1,0},{1,1,0}};
  const unsigned hex[27][3] = {{0,0,0},{2,0,0},{2,2,0},{0,2,0},{0,0,2},{2,0,2},{2,2,2},{0,2,2},{1,0,0},{2,1,0},{1,2,0},{0,1,0},{1,0,2},{2,1,2},{1,2,2},{0,1,2},{0,0,1},{2,0,1},{2,2,1},{0,2,1},{1,0,1},{2,1,1},{1,2,1},{0,1,1},{1,1,0},{1,1,2},{1,1,1}};
  const long double interpolation[3][3] = {{0,1,0},{-0.5L,0,0.5L},{0.5L,-1,0.5L}};
  coefficient.assign(PolynomialSize(dimension,2),0);
  for(unsigned n = 0; n < nodal.size(); ++n) { const unsigned* index = dimension == 1 ? line[n] : (dimension == 2 ? quad[n] : hex[n]); coefficient[PowerIndex({index[0],index[1],index[2]},dimension,2)] = nodal[n]; }
  std::vector<long double> result(coefficient.size());
  for(unsigned axis = 0; axis < dimension; ++axis) {
    for(unsigned n = 0; n < coefficient.size(); ++n) { std::array<unsigned,3> index = PowerExponents(n,dimension,2); const unsigned power = index[axis]; long double sum = 0; long double correction = 0; for(unsigned sample = 0; sample < 3; ++sample) { index[axis] = sample; Accumulate(interpolation[power][sample]*coefficient[PowerIndex(index,dimension,2)],sum,correction); } result[n] = sum; }
    coefficient.swap(result);
  }
}

inline void CutFemLevelSetData::ExpandTriangle(const std::vector<long double>& nodal, std::vector<long double>& coefficient) {
  const long double bubble = 3*(nodal[0]+nodal[1]+nodal[2])-12*(nodal[3]+nodal[4]+nodal[5])+27*nodal[6];
  coefficient.assign(9,0);
  coefficient[0] = nodal[0];
  coefficient[3] = -3*nodal[0]-nodal[1]+4*nodal[3];
  coefficient[1] = -3*nodal[0]-nodal[2]+4*nodal[5];
  coefficient[6] = 2*nodal[0]+2*nodal[1]-4*nodal[3];
  coefficient[2] = 2*nodal[0]+2*nodal[2]-4*nodal[5];
  coefficient[4] = 4*nodal[0]-4*nodal[3]+4*nodal[4]-4*nodal[5]+bubble;
  coefficient[7] = coefficient[5] = -bubble;
}

inline void CutFemLevelSetData::AddBarycentricProduct(const std::array<std::array<long double,4>,4>& barycentric, const std::vector<unsigned>& factors, long double multiplier, std::vector<long double>& coefficient) {
  if(multiplier == 0) return;
  std::vector<long double> product(27,0), next(27,0);
  product[0] = multiplier;
  for(unsigned factor : factors) {
    std::fill(next.begin(),next.end(),0);
    for(unsigned n = 0; n < product.size(); ++n) {
      if(product[n] == 0) continue;
      const std::array<unsigned,3> index = PowerExponents(n,3,2);
      next[n] += product[n]*barycentric[factor][0];
      for(unsigned axis = 0; axis < 3; ++axis) {
        if(barycentric[factor][axis+1] == 0) continue;
        if(index[axis] == 2) throw std::logic_error("TET15 basis expansion exceeded its coordinate degree.");
        std::array<unsigned,3> target = index;
        ++target[axis];
        next[PowerIndex(target,3,2)] += product[n]*barycentric[factor][axis+1];
      }
    }
    product.swap(next);
  }
  for(unsigned n = 0; n < product.size(); ++n) coefficient[n] += product[n];
}

inline void CutFemLevelSetData::ExpandTetrahedron(const std::vector<double>& nodal, std::vector<long double>& coefficient) {
  const std::array<std::array<long double,4>,4> barycentric = {{{1,-1,-1,-1},{0,1,0,0},{0,0,1,0},{0,0,0,1}}};
  const unsigned edges[6][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
  const unsigned faces[4][3] = {{0,1,2},{0,1,3},{1,2,3},{0,2,3}};
  coefficient.assign(27,0);
  long double interior = 256*static_cast<long double>(nodal[14]);
  for(unsigned i = 0; i < 4; ++i) { AddBarycentricProduct(barycentric,{i,i},2*static_cast<long double>(nodal[i]),coefficient); AddBarycentricProduct(barycentric,{i},-static_cast<long double>(nodal[i]),coefficient); interior -= 4*static_cast<long double>(nodal[i]); }
  for(unsigned edge = 0; edge < 6; ++edge) { AddBarycentricProduct(barycentric,{edges[edge][0],edges[edge][1]},4*static_cast<long double>(nodal[4+edge]),coefficient); interior += 32*static_cast<long double>(nodal[4+edge]); }
  for(unsigned face = 0; face < 4; ++face) {
    long double value = 27*static_cast<long double>(nodal[10+face]);
    for(unsigned vertex = 0; vertex < 3; ++vertex) value += 3*static_cast<long double>(nodal[faces[face][vertex]]);
    for(unsigned edge = 0; edge < 6; ++edge) { unsigned contained = 0; for(unsigned vertex = 0; vertex < 3; ++vertex) contained += edges[edge][0] == faces[face][vertex] || edges[edge][1] == faces[face][vertex]; if(contained == 2) value -= 12*static_cast<long double>(nodal[4+edge]); }
    AddBarycentricProduct(barycentric,{faces[face][0],faces[face][1],faces[face][2]},value,coefficient);
    interior -= 108*static_cast<long double>(nodal[10+face]);
  }
  AddBarycentricProduct(barycentric,{0,1,2,3},interior,coefficient);
}

inline void CutFemLevelSetData::EnsurePower() const {
  RequireInitialized();
  if(_powerReady) return;
  if(!_root) throw std::logic_error("A slice must receive its coefficients directly from its parent.");
  if(_geometry == LINE || _geometry == QUAD || _geometry == HEX) ExpandTensor(_geometry,_nodalValues,_power);
  else if(_geometry == TRI) ExpandTriangle(std::vector<long double>(_nodalValues.begin(),_nodalValues.end()),_power);
  else if(_geometry == TET) ExpandTetrahedron(_nodalValues,_power);
  else {
    const unsigned nodes[3][7] = {{0,1,2,6,7,8,18},{12,13,14,15,16,17,20},{3,4,5,9,10,11,19}};
    const long double interpolation[3][3] = {{0,1,0},{-0.5L,0,0.5L},{0.5L,-1,0.5L}};
    _power.assign(27,0);
    for(unsigned plane = 0; plane < 3; ++plane) { std::vector<long double> nodal(7), triangle; for(unsigned n = 0; n < 7; ++n) nodal[n] = _nodalValues[nodes[plane][n]]; ExpandTriangle(nodal,triangle); for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) _power[9*i+3*j+k] += interpolation[k][plane]*triangle[3*i+j]; }
  }
  long double maximum = 0;
  for(long double value : _power) { if(!std::isfinite(value)) throw std::overflow_error("A level-set power coefficient is not finite."); maximum = std::max(maximum,std::abs(value)); }
  if(maximum > 0) { int exponent = 0; std::frexp(maximum,&exponent); _normalization = std::ldexp(1.0L,exponent-1); }
  else _normalization = 1;
  _powerReady = true;
  ++_workCounts[0];
}

inline void CutFemLevelSetData::AffinePower(const std::vector<long double>& input, unsigned dimension, unsigned degree, unsigned axis, long double origin, long double scale, std::vector<long double>& output) {
  if(&input != &output) output = input;
  if(origin == 0 && scale == 1) return;
  std::array<long double,33> scalePowers, originPowers;
  scalePowers[0] = originPowers[0] = 1;
  for(unsigned n = 1; n <= degree; ++n) { scalePowers[n] = scalePowers[n-1]*scale; originPowers[n] = originPowers[n-1]*origin; }
  if(origin == 0) { for(unsigned n = 0; n < input.size(); ++n) output[n] = input[n]*scalePowers[PowerExponents(n,dimension,degree)[axis]]; return; }
  for(unsigned n = 0; n < input.size(); ++n) {
    std::array<unsigned,3> index = PowerExponents(n,dimension,degree);
    if(index[axis] != 0) continue;
    std::array<long double,33> source, transformed;
    for(unsigned power = 0; power <= degree; ++power) { index[axis] = power; source[power] = output[PowerIndex(index,dimension,degree)]; }
    for(unsigned power = 0; power <= degree; ++power) { long double sum = 0; long double correction = 0; for(unsigned old = power; old <= degree; ++old) Accumulate(source[old]*Binomial(old,power)*originPowers[old-power]*scalePowers[power],sum,correction); transformed[power] = sum; }
    for(unsigned power = 0; power <= degree; ++power) { index[axis] = power; output[PowerIndex(index,dimension,degree)] = transformed[power]; }
  }
}

inline void CutFemLevelSetData::BuildBernsteinIndices(GeomElType geometry, const std::array<unsigned,3>& orders, std::vector<std::array<unsigned,3>>& indices) {
  indices.clear();
  const unsigned dimension = Dimension(geometry);
  for(unsigned i = 0; i <= orders[0]; ++i) for(unsigned j = 0; j <= (dimension == 1 ? 0U : ((geometry == TRI || geometry == TET || geometry == WEDGE) ? orders[0]-i : orders[1])); ++j) for(unsigned k = 0; k <= (dimension < 3 ? 0U : (geometry == TET ? orders[0]-i-j : orders[2])); ++k) indices.push_back({i,j,k});
}

inline unsigned CutFemLevelSetData::BernsteinIndex(GeomElType geometry, const std::array<unsigned,3>& orders, const std::array<unsigned,3>& index) {
  if(geometry == LINE) return index[0];
  if(geometry == TRI || geometry == WEDGE) { const unsigned base = index[0]*(2*orders[0]+3-index[0])/2+index[1]; return geometry == TRI ? base : base*(orders[2]+1)+index[2]; }
  if(geometry == TET) return static_cast<unsigned>(Binomial(orders[0]+3,3)-Binomial(orders[0]-index[0]+3,3))+index[1]*(2*(orders[0]-index[0])+3-index[1])/2+index[2];
  return geometry == QUAD ? index[0]*(orders[1]+1)+index[1] : (index[0]*(orders[1]+1)+index[1])*(orders[2]+1)+index[2];
}

inline void CutFemLevelSetData::ConvertToBernstein(const std::vector<long double>& power, GeomElType geometry, const std::array<unsigned,3>& orders, std::vector<long double>& controls) {
  const unsigned dimension = Dimension(geometry);
  std::vector<long double> unit = power, scratch;
  for(unsigned axis = 0; axis < dimension; ++axis) if(geometry == LINE || geometry == QUAD || geometry == HEX || (geometry == WEDGE && axis == 2)) { AffinePower(unit,dimension,2,axis,-1,2,scratch); unit.swap(scratch); }
  std::vector<std::array<unsigned,3>> indices;
  BuildBernsteinIndices(geometry,orders,indices);
  controls.assign(indices.size(),0);
  for(unsigned n = 0; n < indices.size(); ++n) {
    long double correction = 0;
    for(unsigned m = 0; m < unit.size(); ++m) {
      if(unit[m] == 0) continue;
      const std::array<unsigned,3> exponent = PowerExponents(m,dimension,2);
      long double factor = 1;
      if(geometry == TRI || geometry == TET || geometry == WEDGE) {
        const unsigned baseDimension = geometry == TET ? 3U : 2U;
        unsigned total = 0;
        for(unsigned axis = 0; axis < baseDimension; ++axis) { total += exponent[axis]; factor *= Falling(indices[n][axis],exponent[axis]); }
        if(total > orders[0]) throw std::logic_error("Polynomial exceeds the simplex Bernstein order.");
        factor /= Falling(orders[0],total);
        if(geometry == WEDGE) factor *= Binomial(indices[n][2],exponent[2])/Binomial(orders[2],exponent[2]);
      }
      else for(unsigned axis = 0; axis < dimension; ++axis) factor *= Binomial(indices[n][axis],exponent[axis])/Binomial(orders[axis],exponent[axis]);
      Accumulate(unit[m]*factor,controls[n],correction);
    }
  }
}

inline void CutFemLevelSetData::EnsureBernstein() const {
  EnsurePower();
  if(_bernsteinReady) return;
  ConvertToBernstein(_power,_geometry,_bernsteinOrders,_bernstein);
  BuildBernsteinIndices(_geometry,_bernsteinOrders,_bernsteinIndices);
  _bernsteinPadding = Padding(_bernstein);
  _bernsteinReady = true;
  ++_workCounts[1];
}

inline void CutFemLevelSetData::Classify() {
  RequireInitialized();
  if(_classified) return;
  bool positive = false;
  bool negative = false;
  bool zero = false;
  for(double value : _nodalValues) { positive = positive || value > 0; negative = negative || value < 0; zero = zero || value == 0; }
  _phase = 0;
  _classification = "uncertified";
  if(positive && negative) { _classification = "crossing"; _classified = true; return; }
  if(zero) { _classification = "contact"; _classified = true; return; }
  EnsureBernstein();
  const long double padding = _bernsteinPadding;
  const long double minimum = *std::min_element(_bernstein.begin(),_bernstein.end());
  const long double maximum = *std::max_element(_bernstein.begin(),_bernstein.end());
  if(minimum > padding) { _phase = 1; _classification = "positive"; _classified = true; return; }
  if(maximum < -padding) { _phase = -1; _classification = "negative"; _classified = true; return; }
  if(_nodalValues.empty()) {
    for(unsigned n = 0; n < _bernsteinIndices.size(); ++n) {
      const std::array<unsigned,3>& index = _bernsteinIndices[n];
      bool vertex = false;
      if(_geometry == TRI || _geometry == TET || _geometry == WEDGE) {
        const unsigned baseDimension = _geometry == TET ? 3U : 2U;
        unsigned sum = 0;
        unsigned maximumIndex = 0;
        for(unsigned axis = 0; axis < baseDimension; ++axis) { sum += index[axis]; maximumIndex = std::max(maximumIndex,index[axis]); }
        vertex = sum == 0 || maximumIndex == _bernsteinOrders[0];
        if(_geometry == WEDGE) vertex = vertex && (index[2] == 0 || index[2] == _bernsteinOrders[2]);
      }
      else { vertex = true; for(unsigned axis = 0; axis < _dimension; ++axis) vertex = vertex && (index[axis] == 0 || index[axis] == _bernsteinOrders[axis]); }
      if(vertex) { positive = positive || _bernstein[n] > 0; negative = negative || _bernstein[n] < 0; zero = zero || _bernstein[n] == 0; }
    }
    if(positive && negative) { _classification = "crossing"; _classified = true; return; }
    if(zero) { _classification = "contact"; _classified = true; return; }
  }
  const int sign = positive ? 1 : -1;
  if(CertifySign(_geometry,_bernsteinOrders,_bernstein,sign,padding,_maximumDepth)) { _phase = sign; _classification = sign > 0 ? "positive" : "negative"; }
  _classified = true;
}

inline bool CutFemLevelSetData::CertifySign(GeomElType geometry, const std::array<unsigned,3>& orders, const std::vector<long double>& controls, int sign, long double padding, unsigned depth) {
  bool certified = true;
  for(long double value : controls) certified = certified && sign*value > padding;
  if(certified) return true;
  if(depth == 0) return false;
  if(geometry == TRI || geometry == TET || geometry == WEDGE) {
    const unsigned baseDimension = geometry == TET ? 3U : 2U;
    for(unsigned vertex = 0; vertex <= baseDimension; ++vertex) for(unsigned height = 0; height < (geometry == WEDGE ? 2U : 1U); ++height) { std::array<unsigned,3> index = {0,0,0}; if(vertex > 0) index[vertex-1] = orders[0]; if(geometry == WEDGE) index[2] = height*orders[2]; if(sign*controls[BernsteinIndex(geometry,orders,index)] <= padding) return false; }
  }
  else {
    const unsigned dimension = Dimension(geometry);
    for(unsigned vertex = 0; vertex < (1U<<dimension); ++vertex) { std::array<unsigned,3> index = {0,0,0}; for(unsigned axis = 0; axis < dimension; ++axis) index[axis] = ((vertex>>axis)&1U)*orders[axis]; if(sign*controls[BernsteinIndex(geometry,orders,index)] <= padding) return false; }
  }
  std::vector<long double> left, right;
  SplitBernstein(geometry,orders,controls,left,right);
  if(left.empty()) return false;
  return CertifySign(geometry,orders,left,sign,padding,depth-1) && CertifySign(geometry,orders,right,sign,padding,depth-1);
}

inline void CutFemLevelSetData::SplitBernstein(GeomElType geometry, const std::array<unsigned,3>& orders, const std::vector<long double>& controls, std::vector<long double>& left, std::vector<long double>& right) {
  const unsigned dimension = Dimension(geometry);
  const bool simplex = geometry == TRI || geometry == TET || geometry == WEDGE;
  const unsigned baseDimension = geometry == TET ? 3U : 2U;
  std::vector<std::array<unsigned,3>> indices;
  BuildBernsteinIndices(geometry,orders,indices);
  long double variation = 0;
  unsigned selectedAxis = 0;
  unsigned selectedFirst = 0;
  unsigned selectedSecond = 1;
  bool triangleSplit = simplex;
  for(unsigned n = 0; n < indices.size(); ++n) {
    const std::array<unsigned,3>& index = indices[n];
    for(unsigned axis = 0; axis < dimension; ++axis) {
      if(simplex && !(geometry == WEDGE && axis == 2)) continue;
      if(index[axis] >= orders[axis]) continue;
      std::array<unsigned,3> next = index;
      ++next[axis];
      const long double difference = std::abs(controls[n]-controls[BernsteinIndex(geometry,orders,next)]);
      if(difference > variation) { variation = difference; selectedAxis = axis; triangleSplit = false; }
    }
    if(!simplex) continue;
    std::array<unsigned,4> barycentric = {orders[0],0,0,0};
    for(unsigned axis = 0; axis < baseDimension; ++axis) { barycentric[axis+1] = index[axis]; barycentric[0] -= index[axis]; }
    for(unsigned first = 0; first <= baseDimension; ++first) for(unsigned second = first+1; second <= baseDimension; ++second) {
      if(barycentric[second] == 0) continue;
      std::array<unsigned,4> shifted = barycentric;
      ++shifted[first];
      --shifted[second];
      std::array<unsigned,3> next = index;
      for(unsigned axis = 0; axis < baseDimension; ++axis) next[axis] = shifted[axis+1];
      const long double difference = std::abs(controls[n]-controls[BernsteinIndex(geometry,orders,next)]);
      if(difference > variation) { variation = difference; selectedFirst = first; selectedSecond = second; triangleSplit = true; }
    }
  }
  if(variation == 0) return;
  left.resize(controls.size());
  right.resize(controls.size());
  for(unsigned n = 0; n < indices.size(); ++n) {
    long double firstValue = 0, secondValue = 0, firstCorrection = 0, secondCorrection = 0;
    if(!triangleSplit) {
      std::array<unsigned,3> source = indices[n];
      const unsigned index = source[selectedAxis];
      const unsigned degree = orders[selectedAxis];
      for(unsigned k = 0; k <= index; ++k) { source[selectedAxis] = k; Accumulate(Binomial(index,k)*std::ldexp(1.0L,-static_cast<int>(index))*controls[BernsteinIndex(geometry,orders,source)],firstValue,firstCorrection); }
      for(unsigned k = 0; k <= degree-index; ++k) { source[selectedAxis] = index+k; Accumulate(Binomial(degree-index,k)*std::ldexp(1.0L,-static_cast<int>(degree-index))*controls[BernsteinIndex(geometry,orders,source)],secondValue,secondCorrection); }
    }
    else {
      std::array<unsigned,4> barycentric = {orders[0],0,0,0};
      for(unsigned axis = 0; axis < baseDimension; ++axis) { barycentric[axis+1] = indices[n][axis]; barycentric[0] -= indices[n][axis]; }
      const unsigned firstCount = barycentric[selectedFirst];
      const unsigned secondCount = barycentric[selectedSecond];
      for(unsigned k = 0; k <= secondCount; ++k) { std::array<unsigned,4> shifted = barycentric; shifted[selectedFirst] += k; shifted[selectedSecond] -= k; std::array<unsigned,3> source = indices[n]; for(unsigned axis = 0; axis < baseDimension; ++axis) source[axis] = shifted[axis+1]; Accumulate(Binomial(secondCount,k)*std::ldexp(1.0L,-static_cast<int>(secondCount))*controls[BernsteinIndex(geometry,orders,source)],firstValue,firstCorrection); }
      for(unsigned k = 0; k <= firstCount; ++k) { std::array<unsigned,4> shifted = barycentric; shifted[selectedFirst] -= k; shifted[selectedSecond] += k; std::array<unsigned,3> source = indices[n]; for(unsigned axis = 0; axis < baseDimension; ++axis) source[axis] = shifted[axis+1]; Accumulate(Binomial(firstCount,k)*std::ldexp(1.0L,-static_cast<int>(firstCount))*controls[BernsteinIndex(geometry,orders,source)],secondValue,secondCorrection); }
    }
    left[n] = firstValue;
    right[n] = secondValue;
  }
}

inline void CutFemLevelSetData::RestrictionMap(unsigned axis, long double coordinate, GeomElType& childGeometry, std::vector<long double>& origins, std::vector<long double>& scales) const {
  RequireInitialized();
  if(_dimension < 2) throw std::logic_error("A line has no positive-dimensional slice.");
  if(axis >= _dimension) throw std::out_of_range("Slice axis exceeds the active dimension.");
  if(!std::isfinite(coordinate)) throw std::invalid_argument("Slice coordinate must be finite.");
  const bool simplexAxis = _geometry == TRI || _geometry == TET || (_geometry == WEDGE && axis < 2);
  if(coordinate < (simplexAxis ? 0 : -1) || coordinate > 1) throw std::out_of_range("Slice coordinate lies outside the reference element.");
  childGeometry = _geometry == HEX ? QUAD : (_geometry == TET ? TRI : (_geometry == WEDGE ? (axis == 2 ? TRI : QUAD) : LINE));
  origins.assign(_dimension-1,0);
  scales.assign(_dimension-1,1);
  const long double width = 1-coordinate;
  if(_geometry == TET) std::fill(scales.begin(),scales.end(),width);
  if(_geometry == TRI || (_geometry == WEDGE && axis < 2)) origins[0] = scales[0] = width/2;
}

inline void CutFemLevelSetData::RestrictPower(const std::vector<long double>& input, unsigned dimension, unsigned degree, unsigned axis, long double coordinate, const std::vector<long double>& origins, const std::vector<long double>& scales, std::vector<long double>& output) {
  output.resize(PolynomialSize(dimension-1,degree));
  for(unsigned n = 0; n < output.size(); ++n) {
    const std::array<unsigned,3> childIndex = PowerExponents(n,dimension-1,degree);
    std::array<unsigned,3> parentIndex = {0,0,0};
    unsigned childAxis = 0;
    for(unsigned parentAxis = 0; parentAxis < dimension; ++parentAxis) if(parentAxis != axis) parentIndex[parentAxis] = childIndex[childAxis++];
    parentIndex[axis] = degree;
    long double value = input[PowerIndex(parentIndex,dimension,degree)];
    for(unsigned power = degree; power > 0; --power) { parentIndex[axis] = power-1; value = value*coordinate+input[PowerIndex(parentIndex,dimension,degree)]; }
    output[n] = value;
  }
  for(unsigned childAxis = 0; childAxis < dimension-1; ++childAxis) AffinePower(output,dimension-1,degree,childAxis,origins[childAxis],scales[childAxis],output);
}

inline void CutFemLevelSetData::RestrictBernstein(unsigned axis, long double coordinate, CutFemLevelSetData& child) const {
  const bool simplexSlice = _geometry == TRI || _geometry == TET || (_geometry == WEDGE && axis < 2);
  child._bernsteinOrders = {0,0,0};
  if(simplexSlice) {
    if(_geometry == TET) child._bernsteinOrders = {_bernsteinOrders[0],_bernsteinOrders[0],0};
    else if(_geometry == WEDGE) child._bernsteinOrders = {_bernsteinOrders[0],_bernsteinOrders[2],0};
    else child._bernsteinOrders = {_bernsteinOrders[0],0,0};
  }
  else { unsigned next = 0; for(unsigned parentAxis = 0; parentAxis < _dimension; ++parentAxis) if(parentAxis != axis) child._bernsteinOrders[next++] = _bernsteinOrders[parentAxis]; }
  BuildBernsteinIndices(child._geometry,child._bernsteinOrders,child._bernsteinIndices);
  child._bernstein.resize(child._bernsteinIndices.size());
  if(!simplexSlice) {
    const long double unit = (coordinate+1)/2;
    const unsigned degree = _bernsteinOrders[axis];
    std::array<long double,33> weights;
    for(unsigned k = 0; k <= degree; ++k) weights[k] = Bernstein(degree,k,unit);
    for(unsigned n = 0; n < child._bernstein.size(); ++n) {
      std::array<unsigned,3> source = {0,0,0};
      unsigned next = 0;
      for(unsigned parentAxis = 0; parentAxis < _dimension; ++parentAxis) if(parentAxis != axis) source[parentAxis] = child._bernsteinIndices[n][next++];
      long double value = 0, correction = 0;
      for(unsigned k = 0; k <= degree; ++k) { source[axis] = k; Accumulate(weights[k]*_bernstein[BernsteinIndex(_geometry,_bernsteinOrders,source)],value,correction); }
      child._bernstein[n] = value;
    }
  }
  else {
    const unsigned degree = _bernsteinOrders[0];
    const unsigned baseDimension = _geometry == TET ? 3U : 2U;
    std::array<long double,33> fixedPowers, widthPowers;
    fixedPowers[0] = widthPowers[0] = 1;
    for(unsigned k = 1; k <= degree; ++k) { fixedPowers[k] = fixedPowers[k-1]*coordinate; widthPowers[k] = widthPowers[k-1]*(1-coordinate); }
    for(unsigned n = 0; n < child._bernstein.size(); ++n) {
      const std::array<unsigned,3>& childIndex = child._bernsteinIndices[n];
      const unsigned first = childIndex[0];
      const unsigned second = _geometry == TET ? childIndex[1] : 0;
      const unsigned zeroth = degree-first-second;
      long double value = 0, correction = 0;
      for(unsigned a = 0; a <= zeroth; ++a) for(unsigned b = 0; b <= first; ++b) for(unsigned c = 0; c <= second; ++c) {
        const unsigned total = a+b+c;
        std::array<unsigned,3> source = {0,0,0};
        source[axis] = degree-total;
        unsigned next = 0;
        for(unsigned parentAxis = 0; parentAxis < baseDimension; ++parentAxis) if(parentAxis != axis) source[parentAxis] = next++ == 0 ? b : c;
        if(_geometry == WEDGE) source[2] = childIndex[1];
        const long double factor = Binomial(zeroth,a)*Binomial(first,b)*Binomial(second,c)*fixedPowers[degree-total]*widthPowers[total];
        Accumulate(factor*_bernstein[BernsteinIndex(_geometry,_bernsteinOrders,source)],value,correction);
      }
      child._bernstein[n] = value;
    }
  }
  child._bernsteinPadding = _bernsteinPadding+Padding(child._bernstein);
  child._bernsteinReady = true;
}

inline void CutFemLevelSetData::Restrict(unsigned axis, long double coordinate, CutFemLevelSetData& child) const {
  if(&child == this) throw std::invalid_argument("Restriction requires distinct parent and child data.");
  GeomElType childGeometry;
  std::vector<long double> origins, scales;
  RestrictionMap(axis,coordinate,childGeometry,origins,scales);
  EnsurePower();
  EnsureBernstein();
  for(unsigned rootAxis = 0; rootAxis < _rootDimension; ++rootAxis) GetRootDerivative(rootAxis);
  child._initialized = true;
  child._root = false;
  child._geometry = childGeometry;
  child._dimension = _dimension-1;
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
  RestrictPower(_power,_dimension,2,axis,coordinate,origins,scales,child._power);
  child._powerReady = true;
  for(unsigned rootAxis = 0; rootAxis < _rootDimension; ++rootAxis) RestrictPower(GetRootDerivative(rootAxis),_dimension,2,axis,coordinate,origins,scales,child._rootDerivatives[rootAxis]);
  child._parentOrigin.assign(_dimension,0);
  child._parentOrigin[axis] = coordinate;
  child._parentMap.resize(_dimension);
  for(std::vector<long double>& row : child._parentMap) row.assign(child._dimension,0);
  child._parentMeasure = 1;
  unsigned next = 0;
  for(unsigned parentAxis = 0; parentAxis < _dimension; ++parentAxis) if(parentAxis != axis) { child._parentOrigin[parentAxis] = origins[next]; child._parentMap[parentAxis][next] = scales[next]; child._parentMeasure *= std::abs(scales[next]); ++next; }
  child._degenerate = _degenerate || child._parentMeasure == 0;
  child._rootOrigin = ApplyMap(_rootOrigin,_rootMap,child._parentOrigin);
  child._rootMap.resize(_rootDimension);
  for(unsigned rootAxis = 0; rootAxis < _rootDimension; ++rootAxis) { child._rootMap[rootAxis].assign(child._dimension,0); for(unsigned local = 0; local < child._dimension; ++local) for(unsigned parentAxis = 0; parentAxis < _dimension; ++parentAxis) child._rootMap[rootAxis][local] += _rootMap[rootAxis][parentAxis]*child._parentMap[parentAxis][local]; }
  RestrictBernstein(axis,coordinate,child);
  child.Classify();
  ++_workCounts[3];
}

inline std::vector<long double> CutFemLevelSetData::GetLineRoots() const {
  EnsurePower();
  if(_dimension != 1) throw std::logic_error("Line roots require one-dimensional data.");
  const long double a = _power[2]/_normalization;
  const long double b = _power[1]/_normalization;
  const long double c = _power[0]/_normalization;
  std::vector<long double> candidates, result;
  if(a == 0) { if(b == 0) { if(c == 0) throw std::domain_error("An identically zero line has no isolated root set."); return result; } candidates.push_back(-c/b); }
  else {
    const long double square = b*b;
    const long double product = 4*a*c;
    const long double discriminant = (square-product)+(std::fma(b,b,-square)-std::fma(4*a,c,-product));
    if(discriminant < 0) return result;
    if(discriminant == 0) candidates.push_back(-b/(2*a));
    else { const long double q = -0.5L*(b+std::copysign(std::sqrt(discriminant),b)); candidates.push_back(q/a); candidates.push_back(c/q); }
  }
  std::sort(candidates.begin(),candidates.end());
  const long double padding = 64*std::numeric_limits<long double>::epsilon();
  for(long double root : candidates) if(std::isfinite(root) && root >= -1-padding && root <= 1+padding) { root = std::max(-1.0L,std::min(1.0L,root)); if(result.empty() || root != result.back()) result.push_back(root); }
  return result;
}

#endif
