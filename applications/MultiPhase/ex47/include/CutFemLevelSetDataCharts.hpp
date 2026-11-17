#ifndef FEMUS_EX47_CUT_FEM_LEVEL_SET_DATA_CHARTS_HPP
#define FEMUS_EX47_CUT_FEM_LEVEL_SET_DATA_CHARTS_HPP

inline const std::vector<std::array<long double,3>>& CutFemLevelSetData::ReferenceNodes(GeomElType geometry) {
  static const std::array<std::vector<std::array<long double,3>>,6> nodes = {{{{{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1},{0,-1,-1},{1,0,-1},{0,1,-1},{-1,0,-1},{0,-1,1},{1,0,1},{0,1,1},{-1,0,1},{-1,-1,0},{1,-1,0},{1,1,0},{-1,1,0},{0,-1,0},{1,0,0},{0,1,0},{-1,0,0},{0,0,-1},{0,0,1},{0,0,0}}},{{{0,0,0},{1,0,0},{0,1,0},{0,0,1},{0.5L,0,0},{0.5L,0.5L,0},{0,0.5L,0},{0,0,0.5L},{0.5L,0,0.5L},{0,0.5L,0.5L},{1.0L/3,1.0L/3,0},{1.0L/3,0,1.0L/3},{1.0L/3,1.0L/3,1.0L/3},{0,1.0L/3,1.0L/3},{0.25L,0.25L,0.25L}}},{{{0,0,-1},{1,0,-1},{0,1,-1},{0,0,1},{1,0,1},{0,1,1},{0.5L,0,-1},{0.5L,0.5L,-1},{0,0.5L,-1},{0.5L,0,1},{0.5L,0.5L,1},{0,0.5L,1},{0,0,0},{1,0,0},{0,1,0},{0.5L,0,0},{0.5L,0.5L,0},{0,0.5L,0},{1.0L/3,1.0L/3,-1},{1.0L/3,1.0L/3,1},{1.0L/3,1.0L/3,0}}},{{{-1,-1,0},{1,-1,0},{1,1,0},{-1,1,0},{0,-1,0},{1,0,0},{0,1,0},{-1,0,0},{0,0,0}}},{{{0,0,0},{1,0,0},{0,1,0},{0.5L,0,0},{0.5L,0.5L,0},{0,0.5L,0},{1.0L/3,1.0L/3,0}}},{{{-1,0,0},{1,0,0},{0,0,0}}}}};
  return nodes[static_cast<unsigned>(geometry)];
}

inline void CutFemLevelSetData::ExpandNative(GeomElType geometry, const std::vector<long double>& nodal, std::vector<long double>& coefficient) {
  if(geometry == LINE || geometry == QUAD || geometry == HEX) ExpandTensor(geometry,nodal,coefficient);
  else if(geometry == TRI) ExpandTriangle(nodal,coefficient);
  else if(geometry == TET) ExpandTetrahedron(nodal,coefficient);
  else {
    const unsigned nodes[3][7] = {{0,1,2,6,7,8,18},{12,13,14,15,16,17,20},{3,4,5,9,10,11,19}};
    const long double interpolation[3][3] = {{0,1,0},{-0.5L,0,0.5L},{0.5L,-1,0.5L}};
    coefficient.assign(27,0);
    for(unsigned plane = 0; plane < 3; ++plane) {
      std::vector<long double> values(7), triangle;
      for(unsigned n = 0; n < 7; ++n) values[n] = nodal[nodes[plane][n]];
      ExpandTriangle(values,triangle);
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) coefficient[9*i+3*j+k] += interpolation[k][plane]*triangle[3*i+j];
    }
  }
}

inline unsigned CutFemLevelSetData::GetBoundaryCount() const {
  RequireInitialized();
  return _geometry == LINE ? 2U : (_geometry == TRI ? 3U : (_geometry == QUAD || _geometry == TET ? 4U : (_geometry == WEDGE ? 5U : 6U)));
}

inline std::array<std::array<std::array<long double,27>,27>,6> CutFemLevelSetData::BuildInterpolationTables() {
  std::array<std::array<std::array<long double,27>,27>,6> tables = {};
  for(GeomElType geometry : {HEX,TET,WEDGE,QUAD,TRI,LINE}) {
    std::vector<long double> nodal(NodeCount(geometry),0), coefficient;
    for(unsigned node = 0; node < nodal.size(); ++node) { nodal[node] = 1; ExpandNative(geometry,nodal,coefficient); std::copy(coefficient.begin(),coefficient.end(),tables[static_cast<unsigned>(geometry)][node].begin()); nodal[node] = 0; }
  }
  return tables;
}

inline void CutFemLevelSetData::ExpandChartNodal(GeomElType geometry, const std::vector<long double>& nodal, std::vector<long double>& coefficient) {
  static const std::array<std::array<std::array<long double,27>,27>,6> tables = BuildInterpolationTables();
  coefficient.assign(PolynomialSize(Dimension(geometry),2),0);
  std::array<long double,27> corrections = {};
  for(unsigned node = 0; node < nodal.size(); ++node) for(unsigned power = 0; power < coefficient.size(); ++power) if(tables[static_cast<unsigned>(geometry)][node][power] != 0) Accumulate(nodal[node]*tables[static_cast<unsigned>(geometry)][node][power],coefficient[power],corrections[power]);
}

inline void CutFemLevelSetData::EnsureBoundaryData() const {
  RequireInitialized();
  if(_boundaryReady) return;
  if(IsIdenticallyZero()) { _zeroBoundaryMask = (1U<<GetBoundaryCount())-1; _boundaryReady = true; return; }
  if(_root) {
    const std::vector<std::array<long double,3>>& nodes = ReferenceNodes(_geometry);
    for(unsigned boundary = 0; boundary < GetBoundaryCount(); ++boundary) {
      bool zero = true;
      for(unsigned n = 0; n < nodes.size(); ++n) {
        long double distance;
        if(_geometry == TRI || _geometry == TET || (_geometry == WEDGE && boundary < 3)) {
          distance = boundary == 0 ? 1-nodes[n][0]-nodes[n][1]-(_geometry == TET ? nodes[n][2] : 0) : nodes[n][boundary-1];
        }
        else { const unsigned axis = _geometry == WEDGE ? 2U : boundary/2; const unsigned side = _geometry == WEDGE ? boundary-3 : boundary%2; distance = nodes[n][axis]-(side == 0 ? -1 : 1); }
        if(std::abs(distance) <= 64*std::numeric_limits<long double>::epsilon() && _nodalValues[n] != 0) { zero = false; break; }
      }
      if(zero) _zeroBoundaryMask |= 1U<<boundary;
    }
  }
  else {
    EnsureBernstein();
    for(unsigned boundary = 0; boundary < GetBoundaryCount(); ++boundary) {
      bool zero = true;
      for(unsigned n = 0; n < _bernstein.size(); ++n) {
        const std::array<unsigned,3>& index = _bernsteinIndices[n];
        bool onBoundary;
        if(_geometry == TRI || _geometry == TET || (_geometry == WEDGE && boundary < 3)) onBoundary = boundary == 0 ? index[0]+index[1]+(_geometry == TET ? index[2] : 0) == _bernsteinOrders[0] : index[boundary-1] == 0;
        else { const unsigned axis = _geometry == WEDGE ? 2U : boundary/2; const unsigned side = _geometry == WEDGE ? boundary-3 : boundary%2; onBoundary = index[axis] == (side == 0 ? 0 : _bernsteinOrders[axis]); }
        if(onBoundary && _bernstein[n] != 0) { zero = false; break; }
      }
      if(zero && _source && (_zeroBoundaryMask&(1U<<boundary)) == 0) for(const std::array<long double,3>& point : ReferenceNodes(_geometry)) {
        long double distance;
        if(_geometry == TRI || _geometry == TET || (_geometry == WEDGE && boundary < 3)) distance = boundary == 0 ? 1-point[0]-point[1]-(_geometry == TET ? point[2] : 0) : point[boundary-1];
        else { const unsigned axis = _geometry == WEDGE ? 2U : boundary/2, side = _geometry == WEDGE ? boundary-3 : boundary%2; distance = point[axis]-(side == 0 ? -1 : 1); }
        if(std::abs(distance) <= 64*std::numeric_limits<long double>::epsilon()) { const std::array<long double,2> value = CompensatedSourceValue(point); if(value[0] != 0 || value[1] != 0) { zero = false; break; } }
      }
      if(zero) _zeroBoundaryMask |= 1U<<boundary;
    }
  }
  _boundaryReady = true;
}

inline unsigned CutFemLevelSetData::GetZeroBoundaryMask() const { EnsureBoundaryData(); return _zeroBoundaryMask; }

inline bool CutFemLevelSetData::IsZeroOnBoundary(unsigned boundary) const {
  if(boundary >= GetBoundaryCount()) throw std::out_of_range("Boundary index exceeds the element boundary count.");
  EnsureBoundaryData();
  return (_zeroBoundaryMask&(1U<<boundary)) != 0;
}

inline void CutFemLevelSetData::PropagateBoundaryData(unsigned axis, long double coordinate, CutFemLevelSetData& child) const {
  const bool simplex = _geometry == TRI || _geometry == TET || (_geometry == WEDGE && axis < 2);
  if(simplex ? coordinate == 0 : (coordinate == -1 || coordinate == 1)) {
    const unsigned boundary = simplex ? axis+1 : (_geometry == WEDGE ? 3U+(coordinate == 1) : 2*axis+(coordinate == 1));
    if(IsZeroOnBoundary(boundary)) {
      child._zeroReady = child._identicallyZero = true;
      child._boundaryReady = true;
      child._zeroBoundaryMask = (1U<<child.GetBoundaryCount())-1;
      std::fill(child._power.begin(),child._power.end(),0);
      std::fill(child._bernstein.begin(),child._bernstein.end(),0);
      child.UpdateBernsteinBounds(_bernsteinPadding);
      return;
    }
  }
  if(!_boundaryReady) return;
  if(_geometry == TET || (_geometry == WEDGE && axis == 2)) {
    if((_zeroBoundaryMask&1U) != 0) child._zeroBoundaryMask |= 1U;
    unsigned next = 1;
    for(unsigned direction = 0; direction < (_geometry == TET ? 3U : 2U); ++direction) if(direction != axis) { if((_zeroBoundaryMask&(1U<<(direction+1))) != 0) child._zeroBoundaryMask |= 1U<<next; ++next; }
  }
  else if(_geometry == TRI || (_geometry == WEDGE && axis < 2)) {
    const unsigned other = 1-axis;
    if((_zeroBoundaryMask&(1U<<(other+1))) != 0) child._zeroBoundaryMask |= 1U;
    if((_zeroBoundaryMask&1U) != 0) child._zeroBoundaryMask |= 2U;
    if(_geometry == WEDGE) { if((_zeroBoundaryMask&8U) != 0) child._zeroBoundaryMask |= 4U; if((_zeroBoundaryMask&16U) != 0) child._zeroBoundaryMask |= 8U; }
  }
  else {
    unsigned next = 0;
    for(unsigned direction = 0; direction < _dimension; ++direction) if(direction != axis) { if((_zeroBoundaryMask&(1U<<(2*direction))) != 0) child._zeroBoundaryMask |= 1U<<(2*next); if((_zeroBoundaryMask&(1U<<(2*direction+1))) != 0) child._zeroBoundaryMask |= 1U<<(2*next+1); ++next; }
  }
}

inline void CutFemLevelSetData::InitializeChart(CutFemLevelSetData& chart) const {
  if(&chart == this) throw std::invalid_argument("Chart input and output must be distinct.");
  EnsurePower();
  EnsureBernstein();
  chart._initialized = true;
  chart._root = false;
  chart._rootGradientAvailable = _rootGradientAvailable;
  chart._geometry = _geometry;
  chart._dimension = _dimension;
  chart._rootDimension = _rootDimension;
  chart._parentDimension = _dimension;
  chart._maximumDepth = _maximumDepth;
  chart._meshLevel = _meshLevel;
  chart._elementId = _elementId;
  chart._nodalVersion = _nodalVersion;
  ++chart._revision;
  chart._nodalValues.clear();
  chart._classified = false;
  chart._normalization = _normalization;
  chart._normalizedReady = false;
  chart._powerReady = true;
  chart._bernsteinReady = true;
  chart._derivativeReady.fill(false);
  chart._secondDerivativeReady.fill(false);
  chart._derivativeBoundsReady.fill(false);
  chart._graphDirectionReady = false;
  chart._degenerate = _degenerate;
  chart._parentMeasure = 1;
  chart._source.reset();
  chart._restrictionSteps = 0;
  chart.InvalidateAnalysis();
  chart._parentOrigin.assign(_dimension,0);
  chart._parentMap.resize(_dimension);
  for(std::vector<long double>& row : chart._parentMap) row.assign(_dimension,0);
}

inline void CutFemLevelSetData::ComposeChartMaps(CutFemLevelSetData& chart) const {
  chart._rootOrigin = _rootOrigin;
  chart._rootMap.resize(_rootDimension);
  for(unsigned rootAxis = 0; rootAxis < _rootDimension; ++rootAxis) {
    chart._rootMap[rootAxis].assign(chart._dimension,0);
    for(unsigned parentAxis = 0; parentAxis < _dimension; ++parentAxis) {
      chart._rootOrigin[rootAxis] += _rootMap[rootAxis][parentAxis]*chart._parentOrigin[parentAxis];
      for(unsigned local = 0; local < chart._dimension; ++local) chart._rootMap[rootAxis][local] += _rootMap[rootAxis][parentAxis]*chart._parentMap[parentAxis][local];
    }
  }
}

inline void CutFemLevelSetData::PermuteAxes(const std::array<unsigned,3>& axes, CutFemLevelSetData& chart) const {
  RequireInitialized();
  unsigned mask = 0;
  for(unsigned axis = 0; axis < _dimension; ++axis) { if(axes[axis] >= _dimension || (mask&(1U<<axes[axis])) != 0) throw std::invalid_argument("Invalid active-axis permutation."); mask |= 1U<<axes[axis]; }
  if(_geometry == WEDGE && axes[2] != 2) throw std::invalid_argument("A wedge axis permutation must preserve its height axis.");
  InitializeChart(chart);
  chart._power.resize(_power.size());
  for(unsigned n = 0; n < _power.size(); ++n) {
    const std::array<unsigned,3> source = PowerExponents(n,_dimension,2);
    std::array<unsigned,3> target = {0,0,0};
    for(unsigned axis = 0; axis < _dimension; ++axis) target[axis] = source[axes[axis]];
    const unsigned index = PowerIndex(target,_dimension,2);
    chart._power[index] = _power[n];
  }
  chart._bernsteinOrders = {0,0,0};
  for(unsigned axis = 0; axis < _dimension; ++axis) { chart._bernsteinOrders[axis] = _bernsteinOrders[axes[axis]]; chart._parentMap[axes[axis]][axis] = 1; }
  chart._bernsteinIndices = CachedBernsteinIndices(_geometry,chart._bernsteinOrders);
  chart._bernstein.resize(_bernstein.size());
  for(unsigned n = 0; n < chart._bernstein.size(); ++n) {
    std::array<unsigned,3> source = {0,0,0};
    for(unsigned axis = 0; axis < _dimension; ++axis) source[axes[axis]] = chart._bernsteinIndices[n][axis];
    chart._bernstein[n] = _bernstein[BernsteinIndex(_geometry,_bernsteinOrders,source)];
  }
  if(_rootGradientAvailable) for(unsigned component = 0; component < _rootDimension; ++component) {
    const std::vector<long double>& derivative = GetRootDerivative(component);
    chart._rootDerivatives[component].resize(_power.size());
    for(unsigned n = 0; n < _power.size(); ++n) { const std::array<unsigned,3> source = PowerExponents(n,_dimension,2); std::array<unsigned,3> target = {0,0,0}; for(unsigned axis = 0; axis < _dimension; ++axis) target[axis] = source[axes[axis]]; chart._rootDerivatives[component][PowerIndex(target,_dimension,2)] = derivative[n]; }
  }
  else for(std::vector<long double>& derivative : chart._rootDerivatives) derivative.clear();
  if(_boundaryReady) {
    if(_geometry == TRI || _geometry == TET || _geometry == WEDGE) { chart._zeroBoundaryMask = _zeroBoundaryMask&1U; for(unsigned axis = 0; axis < (_geometry == WEDGE ? 2U : _dimension); ++axis) if((_zeroBoundaryMask&(1U<<(axes[axis]+1))) != 0) chart._zeroBoundaryMask |= 1U<<(axis+1); if(_geometry == WEDGE) chart._zeroBoundaryMask |= _zeroBoundaryMask&24U; }
    else for(unsigned axis = 0; axis < _dimension; ++axis) for(unsigned side = 0; side < 2; ++side) if((_zeroBoundaryMask&(1U<<(2*axes[axis]+side))) != 0) chart._zeroBoundaryMask |= 1U<<(2*axis+side);
  }
  ComposeChartMaps(chart);
  chart.UpdateBernsteinBounds(_bernsteinPadding);
  chart.Classify();
}

inline void CutFemLevelSetData::PermuteSimplexVertices(const std::array<unsigned,4>& vertices, CutFemLevelSetData& chart) const {
  RequireInitialized();
  if(_geometry != TRI && _geometry != TET && _geometry != WEDGE) throw std::logic_error("A barycentric chart requires a simplex or wedge.");
  const unsigned baseDimension = _geometry == TET ? 3U : 2U;
  unsigned mask = 0;
  for(unsigned vertex = 0; vertex <= baseDimension; ++vertex) { if(vertices[vertex] > baseDimension || (mask&(1U<<vertices[vertex])) != 0) throw std::invalid_argument("Invalid barycentric vertex permutation."); mask |= 1U<<vertices[vertex]; }
  InitializeChart(chart);
  for(unsigned row = 0; row < baseDimension; ++row) {
    chart._parentOrigin[row] = vertices[0] == row+1 ? 1 : 0;
    for(unsigned column = 0; column < baseDimension; ++column) chart._parentMap[row][column] = (vertices[column+1] == row+1 ? 1 : 0)-chart._parentOrigin[row];
  }
  if(_geometry == WEDGE) chart._parentMap[2][2] = 1;
  const std::vector<std::array<long double,3>>& nodes = ReferenceNodes(_geometry);
  chart._chartNodal.resize(nodes.size());
  std::vector<long double>& nodal = chart._chartNodal;
  for(unsigned n = 0; n < nodes.size(); ++n) {
    std::array<long double,3> point = {0,0,0};
    for(unsigned row = 0; row < _dimension; ++row) { point[row] = chart._parentOrigin[row]; for(unsigned column = 0; column < _dimension; ++column) point[row] += chart._parentMap[row][column]*nodes[n][column]; }
    if(_root) {
      unsigned match = nodes.size();
      for(unsigned candidate = 0; candidate < nodes.size(); ++candidate) { bool equal = true; for(unsigned axis = 0; axis < _dimension; ++axis) equal = equal && std::abs(nodes[candidate][axis]-point[axis]) <= 64*std::numeric_limits<long double>::epsilon(); if(equal) { match = candidate; break; } }
      if(match == nodes.size()) throw std::logic_error("A supported vertex permutation must preserve the FEMuS nodes.");
      nodal[n] = _nodalValues[match];
    }
    else nodal[n] = EvaluateUnchecked(_power.data(),_dimension,2,point.data());
  }
  ExpandChartNodal(_geometry,nodal,chart._power);
  if(_rootGradientAvailable) for(unsigned component = 0; component < _rootDimension; ++component) {
    const std::vector<long double>& derivative = GetRootDerivative(component);
    for(unsigned n = 0; n < nodes.size(); ++n) { std::array<long double,3> point = {0,0,0}; for(unsigned row = 0; row < _dimension; ++row) { point[row] = chart._parentOrigin[row]; for(unsigned column = 0; column < _dimension; ++column) point[row] += chart._parentMap[row][column]*nodes[n][column]; } nodal[n] = EvaluateUnchecked(derivative.data(),_dimension,2,point.data()); }
    ExpandChartNodal(_geometry,nodal,chart._rootDerivatives[component]);
  }
  else for(std::vector<long double>& derivative : chart._rootDerivatives) derivative.clear();
  chart._bernsteinOrders = _bernsteinOrders;
  chart._bernsteinIndices = _bernsteinIndices;
  chart._bernstein.resize(_bernstein.size());
  for(unsigned n = 0; n < chart._bernstein.size(); ++n) {
    std::array<unsigned,4> exponent = {_bernsteinOrders[0],0,0,0}, sourceExponent = {0,0,0,0};
    for(unsigned axis = 0; axis < baseDimension; ++axis) { exponent[axis+1] = _bernsteinIndices[n][axis]; exponent[0] -= exponent[axis+1]; }
    for(unsigned vertex = 0; vertex <= baseDimension; ++vertex) sourceExponent[vertices[vertex]] = exponent[vertex];
    std::array<unsigned,3> source = {sourceExponent[1],sourceExponent[2],sourceExponent[3]};
    if(_geometry == WEDGE) source[2] = _bernsteinIndices[n][2];
    chart._bernstein[n] = _bernstein[BernsteinIndex(_geometry,_bernsteinOrders,source)];
  }
  EnsureBoundaryData();
  for(unsigned vertex = 0; vertex <= baseDimension; ++vertex) if((_zeroBoundaryMask&(1U<<vertices[vertex])) != 0) chart._zeroBoundaryMask |= 1U<<vertex;
  if(_geometry == WEDGE) chart._zeroBoundaryMask |= _zeroBoundaryMask&24U;
  ComposeChartMaps(chart);
  chart.UpdateBernsteinBounds(_bernsteinPadding+Padding(chart._power));
  chart.Classify();
}

inline void CutFemLevelSetData::RestrictBoundary(unsigned boundary, CutFemLevelSetData& child) const {
  if(&child == this) throw std::invalid_argument("Boundary input and output must be distinct.");
  if(boundary >= GetBoundaryCount()) throw std::out_of_range("Boundary index exceeds the element boundary count.");
  if(_dimension == 1) throw std::logic_error("A line boundary is a point, not a positive-dimensional data object.");
  EnsureBoundaryData();
  if(_geometry == TRI || _geometry == TET || (_geometry == WEDGE && boundary < 3)) {
    if(boundary != 0) { Restrict(boundary-1,0,child); return; }
    CutFemLevelSetData chart;
    PermuteSimplexVertices({1,0,2,3},chart);
    chart.Restrict(0,0,child);
    const std::vector<long double> origin = child._parentOrigin;
    const std::vector<std::vector<long double>> map = child._parentMap;
    for(unsigned row = 0; row < _dimension; ++row) {
      child._parentOrigin[row] = chart._parentOrigin[row];
      for(unsigned parentAxis = 0; parentAxis < _dimension; ++parentAxis) child._parentOrigin[row] += chart._parentMap[row][parentAxis]*origin[parentAxis];
      for(unsigned local = 0; local < child._dimension; ++local) { child._parentMap[row][local] = 0; for(unsigned parentAxis = 0; parentAxis < _dimension; ++parentAxis) child._parentMap[row][local] += chart._parentMap[row][parentAxis]*map[parentAxis][local]; }
    }
    if(child._dimension == 1) { long double square = 0; for(unsigned row = 0; row < _dimension; ++row) square += child._parentMap[row][0]*child._parentMap[row][0]; child._parentMeasure = std::sqrt(square); }
    else { std::array<long double,3> normal; for(unsigned row = 0; row < 3; ++row) normal[row] = child._parentMap[(row+1)%3][0]*child._parentMap[(row+2)%3][1]-child._parentMap[(row+2)%3][0]*child._parentMap[(row+1)%3][1]; child._parentMeasure = std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]); }
  }
  else { const unsigned axis = _geometry == WEDGE ? 2U : boundary/2; const unsigned side = _geometry == WEDGE ? boundary-3 : boundary%2; Restrict(axis,side == 0 ? -1 : 1,child); }
}

inline const std::vector<long double>& CutFemLevelSetData::GetUnitPowerCoefficients() const {
  EnsurePower();
  if(_unitPowerReady) return _unitPower;
  _unitPower = _power;
  for(unsigned axis = 0; axis < _dimension; ++axis) if(_geometry == LINE || _geometry == QUAD || _geometry == HEX || (_geometry == WEDGE && axis == 2)) AffinePower(_unitPower,_dimension,2,axis,-1,2,_unitPower);
  _unitPowerReady = true;
  return _unitPower;
}

inline void CutFemLevelSetData::TransformParentComponents(unsigned degree, const std::vector<std::vector<long double>>& components, unsigned& outputDegree, std::vector<std::vector<long double>>& output) const {
  RequireInitialized();
  if(&components == &output) throw std::invalid_argument("Chart component input and output must be distinct.");
  unsigned transformedDegree = degree;
  for(const std::vector<long double>& component : components) {
    ValidatePolynomial(component,_parentDimension,degree);
    for(unsigned n = 0; n < component.size(); ++n) if(component[n] != 0) {
      const std::array<unsigned,3> index = PowerExponents(n,_parentDimension,degree);
      for(unsigned column = 0; column < _dimension; ++column) { unsigned power = 0; for(unsigned row = 0; row < _parentDimension; ++row) if(_parentMap[row][column] != 0) power += index[row]; transformedDegree = std::max(transformedDegree,power); }
    }
  }
  const unsigned size = PolynomialSize(_dimension,transformedDegree);
  outputDegree = transformedDegree;
  output.resize(components.size());
  bool permutation = _dimension == _parentDimension;
  bool barycentric = _dimension == _parentDimension;
  int mixedRow = -1;
  std::array<unsigned,3> mixedAxes = {0,0,0};
  unsigned mixedCount = 0;
  std::array<unsigned,3> permutationAxes = {0,0,0};
  for(unsigned row = 0; row < _parentDimension; ++row) {
    unsigned nonzero = 0;
    for(unsigned column = 0; column < _dimension; ++column) if(_parentMap[row][column] != 0) { ++nonzero; permutationAxes[row] = column; permutation = permutation && std::abs(_parentMap[row][column]) == 1; }
    permutation = permutation && nonzero == 1 && _parentOrigin[row] == 0;
    if(_parentOrigin[row] == 1 && mixedRow < 0) {
      mixedRow = static_cast<int>(row);
      for(unsigned column = 0; column < _dimension; ++column) if(_parentMap[row][column] != 0) { barycentric = barycentric && _parentMap[row][column] == -1; mixedAxes[mixedCount++] = column; }
      barycentric = barycentric && mixedCount >= 2;
    }
    else barycentric = barycentric && nonzero == 1 && _parentOrigin[row] == 0 && _parentMap[row][permutationAxes[row]] == 1;
  }
  barycentric = barycentric && mixedRow >= 0;
  if(!permutation && !barycentric) { _componentTerm.resize(size); _componentNext.resize(size); }
  for(unsigned component = 0; component < components.size(); ++component) {
    output[component].assign(size,0);
    for(unsigned n = 0; n < components[component].size(); ++n) if(components[component][n] != 0) {
      const std::array<unsigned,3> index = PowerExponents(n,_parentDimension,degree);
      if(permutation) {
        std::array<unsigned,3> target = {0,0,0};
        long double value = components[component][n];
        for(unsigned row = 0; row < _parentDimension; ++row) { target[permutationAxes[row]] = index[row]; if(_parentMap[row][permutationAxes[row]] < 0 && index[row]%2U != 0) value = -value; }
        output[component][PowerIndex(target,_dimension,transformedDegree)] = value;
        continue;
      }
      if(barycentric) {
        std::array<unsigned,3> base = {0,0,0};
        for(unsigned row = 0; row < _parentDimension; ++row) if(static_cast<int>(row) != mixedRow) base[permutationAxes[row]] += index[row];
        const unsigned power = index[static_cast<unsigned>(mixedRow)];
        for(unsigned first = 0; first <= power; ++first) for(unsigned second = 0; second <= power-first; ++second) for(unsigned third = 0; third <= (mixedCount == 3 ? power-first-second : 0); ++third) {
          std::array<unsigned,3> target = base;
          target[mixedAxes[0]] += first;
          target[mixedAxes[1]] += second;
          if(mixedCount == 3) target[mixedAxes[2]] += third;
          const long double multiplier = Binomial(power,first)*Binomial(power-first,second)*Binomial(power-first-second,third)*((first+second+third)%2U == 0 ? 1 : -1);
          output[component][PowerIndex(target,_dimension,transformedDegree)] += components[component][n]*multiplier;
        }
        continue;
      }
      std::fill(_componentTerm.begin(),_componentTerm.end(),0);
      _componentTerm[0] = components[component][n];
      for(unsigned row = 0; row < _parentDimension; ++row) for(unsigned factor = 0; factor < index[row]; ++factor) {
        std::fill(_componentNext.begin(),_componentNext.end(),0);
        for(unsigned old = 0; old < size; ++old) if(_componentTerm[old] != 0) {
          const std::array<unsigned,3> exponent = PowerExponents(old,_dimension,transformedDegree);
          _componentNext[old] += _componentTerm[old]*_parentOrigin[row];
          for(unsigned column = 0; column < _dimension; ++column) if(_parentMap[row][column] != 0) { std::array<unsigned,3> target = exponent; ++target[column]; if(target[column] > transformedDegree) throw std::logic_error("Chart polynomial degree calculation failed."); _componentNext[PowerIndex(target,_dimension,transformedDegree)] += _componentTerm[old]*_parentMap[row][column]; }
        }
        _componentTerm.swap(_componentNext);
      }
      for(unsigned target = 0; target < size; ++target) output[component][target] += _componentTerm[target];
    }
  }
}

inline long double CutFemLevelSetData::GetUnitToReferenceMeasureScale() const {
  RequireInitialized();
  return _geometry == LINE || _geometry == QUAD || _geometry == HEX ? static_cast<long double>(1U<<_dimension) : (_geometry == WEDGE ? 2 : 1);
}

#endif
