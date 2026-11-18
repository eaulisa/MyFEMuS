#ifndef FEMUS_EX47_CUT_FEM_LEVEL_SET_DATA_ANALYSIS_HPP
#define FEMUS_EX47_CUT_FEM_LEVEL_SET_DATA_ANALYSIS_HPP

inline void CutFemLevelSetData::InvalidateAnalysis() {
  _lineAnalysisReady = false;
  _zeroReady = false;
  _boundaryReady = false;
  _zeroBoundaryMask = 0;
  _dependenceReady = false;
  _graphCertificateReady.fill(false);
  _graphWeakReady.fill(false);
  _zeroFiberReady.fill(false);
  _projectionReady.fill(false);
  _unitPowerReady = false;
}

inline std::array<long double,2> CutFemLevelSetData::CompensatedAdd(const std::array<long double,2>& first, const std::array<long double,2>& second) {
  const long double sum = first[0]+second[0];
  const long double part = sum-first[0];
  const long double error = (first[0]-(sum-part))+(second[0]-part)+first[1]+second[1];
  const long double value = sum+error;
  const long double tail = value-sum;
  return {value,(sum-(value-tail))+(error-tail)};
}

inline std::array<long double,2> CutFemLevelSetData::CompensatedMultiply(const std::array<long double,2>& first, const std::array<long double,2>& second) {
  const long double product = first[0]*second[0];
  long double error = std::fma(first[0],second[0],-product);
  if(first[1] != 0) error += first[1]*second[0];
  if(second[1] != 0) error += first[0]*second[1]+first[1]*second[1];
  return CompensatedAdd({product,0},{error,0});
}

inline std::array<long double,2> CutFemLevelSetData::CompensatedEvaluate(const long double* coefficient, unsigned dimension, unsigned degree, const std::array<std::array<long double,2>,3>& point) {
  unsigned stride = 1;
  for(unsigned axis = 1; axis < dimension; ++axis) stride *= degree+1;
  std::array<long double,2> value = {0,0};
  std::array<std::array<long double,2>,3> next = point;
  next[0] = point[1];
  next[1] = point[2];
  for(unsigned power = degree+1; power > 0; --power) {
    const std::array<long double,2> term = dimension == 1 ? std::array<long double,2>{{coefficient[power-1],0}} : CompensatedEvaluate(coefficient+(power-1)*stride,dimension-1,degree,next);
    value = CompensatedAdd(CompensatedMultiply(value,point[0]),term);
  }
  return value;
}

inline long double CutFemLevelSetData::QuadraticDiscriminant(const std::array<std::array<long double,2>,3>& coefficient) {
  const std::array<long double,2> value = CompensatedAdd(CompensatedMultiply(coefficient[1],coefficient[1]),CompensatedMultiply({-4,0},CompensatedMultiply(coefficient[2],coefficient[0])));
  return value[0]+value[1];
}

inline unsigned CutFemLevelSetData::SolveQuadratic(long double a, long double b, long double c, long double discriminant, std::array<long double,2>& roots) {
  if(a == 0) { if(b == 0) return 0; roots[0] = -c/b; return 1; }
  if(discriminant < 0) return 0;
  if(discriminant == 0) { roots[0] = -b/(2*a); return 1; }
  const long double q = -0.5L*(b+std::copysign(std::sqrt(discriminant),b));
  roots = {q/a,c/q};
  if(roots[1] < roots[0]) std::swap(roots[0],roots[1]);
  return 2;
}

inline void CutFemLevelSetData::EnsureSource() const {
  if(_source) return;
  EnsurePower();
  std::shared_ptr<PolynomialSource> source = std::make_shared<PolynomialSource>();
  source->geometry = _geometry;
  source->dimension = _dimension;
  std::copy(_power.begin(),_power.end(),source->power.begin());
  for(long double value : _power) source->norm += std::abs(value);
  _source = source;
}

inline void CutFemLevelSetData::MapCompensatedSourcePoint(const std::array<long double,3>& coordinate, std::array<std::array<long double,2>,3>& point, std::array<std::array<long double,2>,3>* direction, unsigned axis) const {
  point = {{{coordinate[0],0},{coordinate[1],0},{coordinate[2],0}}};
  if(direction != nullptr) { *direction = {}; (*direction)[axis] = {1,0}; }
  if(!_source || _restrictionSteps == 0) return;
  std::array<GeomElType,3> geometries = {_source->geometry,LINE,LINE};
  for(unsigned step = 0; step < _restrictionSteps; ++step) geometries[step+1] = geometries[step] == HEX ? QUAD : (geometries[step] == TET ? TRI : (geometries[step] == WEDGE ? (_restrictionAxes[step] == 2 ? TRI : QUAD) : LINE));
  for(unsigned step = _restrictionSteps; step > 0; --step) {
    const unsigned fixedAxis = _restrictionAxes[step-1];
    const long double fixed = _restrictionCoordinates[step-1];
    const GeomElType geometry = geometries[step-1];
    const std::array<long double,2> width = CompensatedAdd({1,0},{-fixed,0});
    std::array<std::array<long double,2>,3> parent = {}, parentDirection = {};
    parent[fixedAxis] = {fixed,0};
    unsigned next = 0;
    for(unsigned parentAxis = 0; parentAxis < Dimension(geometry); ++parentAxis) if(parentAxis != fixedAxis) {
      parent[parentAxis] = point[next];
      if(direction != nullptr) parentDirection[parentAxis] = (*direction)[next];
      if(geometry == TET) { parent[parentAxis] = CompensatedMultiply(width,point[next]); if(direction != nullptr) parentDirection[parentAxis] = CompensatedMultiply(width,(*direction)[next]); }
      else if((geometry == TRI || (geometry == WEDGE && fixedAxis < 2)) && next == 0) { const std::array<long double,2> scale = CompensatedMultiply(width,{0.5L,0}); parent[parentAxis] = CompensatedMultiply(scale,CompensatedAdd({1,0},point[next])); if(direction != nullptr) parentDirection[parentAxis] = CompensatedMultiply(scale,(*direction)[next]); }
      ++next;
    }
    point = parent;
    if(direction != nullptr) *direction = parentDirection;
  }
}

inline std::array<long double,2> CutFemLevelSetData::CompensatedSourceValue(const std::array<long double,3>& coordinate) const {
  std::array<std::array<long double,2>,3> point;
  MapCompensatedSourcePoint(coordinate,point);
  return _source && _restrictionSteps != 0 ? CompensatedEvaluate(_source->power.data(),_source->dimension,2,point) : CompensatedEvaluate(_power.data(),_dimension,2,point);
}

inline void CutFemLevelSetData::GetCompensatedFiberCoefficients(unsigned axis, const std::array<long double,3>& coordinate, std::array<std::array<long double,2>,3>& coefficient) const {
  std::array<long double,3> origin = coordinate;
  origin[axis] = 0;
  std::array<std::array<long double,2>,3> point, direction;
  MapCompensatedSourcePoint(origin,point,&direction,axis);
  const bool ancestor = _source && _restrictionSteps != 0;
  const unsigned dimension = ancestor ? _source->dimension : _dimension;
  const long double* power = ancestor ? _source->power.data() : _power.data();
  unsigned freeAxis = dimension;
  for(unsigned sourceAxis = 0; sourceAxis < dimension; ++sourceAxis) if(direction[sourceAxis][0] != 0 || direction[sourceAxis][1] != 0) { if(freeAxis != dimension) throw std::logic_error("A compensated fiber requires one free source coordinate."); freeAxis = sourceAxis; }
  coefficient = {};
  if(freeAxis == dimension) { coefficient[0] = CompensatedSourceValue(origin); coefficient[0][0] /= _normalization; coefficient[0][1] /= _normalization; return; }
  std::array<std::array<std::array<long double,2>,3>,3> powers;
  for(unsigned sourceAxis = 0; sourceAxis < dimension; ++sourceAxis) if(sourceAxis != freeAxis) { powers[sourceAxis][0] = {1,0}; powers[sourceAxis][1] = point[sourceAxis]; powers[sourceAxis][2] = CompensatedMultiply(point[sourceAxis],point[sourceAxis]); }
  const unsigned count = dimension == 1 ? 3U : (dimension == 2 ? 9U : 27U);
  for(unsigned n = 0; n < count; ++n) if(power[n] != 0) {
    const std::array<unsigned,3> exponent = PowerExponents(n,dimension,2);
    std::array<long double,2> term = {power[n],0};
    for(unsigned sourceAxis = 0; sourceAxis < dimension; ++sourceAxis) if(sourceAxis != freeAxis && exponent[sourceAxis] != 0) term = CompensatedMultiply(term,powers[sourceAxis][exponent[sourceAxis]]);
    coefficient[exponent[freeAxis]] = CompensatedAdd(coefficient[exponent[freeAxis]],term);
  }
  const std::array<long double,2> a = coefficient[2], b = coefficient[1], c = coefficient[0], scale = direction[freeAxis];
  coefficient[0] = CompensatedAdd(CompensatedMultiply(CompensatedAdd(CompensatedMultiply(a,point[freeAxis]),b),point[freeAxis]),c);
  coefficient[1] = CompensatedMultiply(scale,CompensatedAdd(CompensatedMultiply(CompensatedMultiply({2,0},a),point[freeAxis]),b));
  coefficient[2] = CompensatedMultiply(a,CompensatedMultiply(scale,scale));
  for(std::array<long double,2>& value : coefficient) { value[0] /= _normalization; value[1] /= _normalization; }
}

inline void CutFemLevelSetData::GetCompensatedLineCoefficients(std::array<std::array<long double,2>,3>& coefficient) const {
  if(!_source || _restrictionSteps == 0) { for(unsigned n = 0; n < 3; ++n) coefficient[n] = {_power[n]/_normalization,0}; return; }
  GetCompensatedFiberCoefficients(0,{0,0,0},coefficient);
}

inline bool CutFemLevelSetData::IsIdenticallyZero() const {
  EnsurePower();
  if(_zeroReady) return _identicallyZero;
  _identicallyZero = true;
  for(long double value : _power) if(value != 0) { _identicallyZero = false; break; }
  if(_identicallyZero && _source && _restrictionSteps != 0) for(const std::array<long double,3>& point : ReferenceNodes(_geometry)) { const std::array<long double,2> value = CompensatedSourceValue(point); if(value[0] != 0 || value[1] != 0) { _identicallyZero = false; break; } }
  _zeroReady = true;
  return _identicallyZero;
}

inline const CutFemLevelSetData::LineAnalysis& CutFemLevelSetData::GetLineAnalysis() const {
  EnsurePower();
  if(_dimension != 1) throw std::logic_error("Line analysis requires one-dimensional data.");
  if(_lineAnalysisReady) return *_lineAnalysis;
  EnsureBoundaryData();
  if(!_lineAnalysis || !_lineAnalysis.unique()) _lineAnalysis = std::make_shared<LineAnalysis>();
  *_lineAnalysis = LineAnalysis();
  std::array<std::array<long double,2>,3> coefficient = {{{_power[0]/_normalization,0},{_power[1]/_normalization,0},{_power[2]/_normalization,0}}};
  long double a = coefficient[2][0], b = coefficient[1][0], c = coefficient[0][0];
  const long double square = b*b, product = 4*a*c;
  long double discriminant = square-product;
  const long double coefficientNorm = std::abs(a)+std::abs(b)+std::abs(c), sourceNorm = _source ? _source->norm/_normalization : 0;
  const long double threshold = std::max(2048*std::numeric_limits<long double>::epsilon()*(std::abs(square)+std::abs(product)),0.001L*coefficientNorm*sourceNorm);
  if((a != 0 && std::abs(discriminant) <= threshold) || (_source && coefficientNorm <= 4096*std::numeric_limits<long double>::epsilon()*sourceNorm)) {
    GetCompensatedLineCoefficients(coefficient);
    a = coefficient[2][0]+coefficient[2][1];
    b = coefficient[1][0]+coefficient[1][1];
    c = coefficient[0][0]+coefficient[0][1];
    discriminant = QuadraticDiscriminant(coefficient);
    _lineAnalysis->compensated = true;
  }
  else if(product > 0 && std::abs(discriminant) < 0.125L*(std::abs(square)+std::abs(product))) discriminant = QuadraticDiscriminant(coefficient);
  if(_zeroBoundaryMask == 3) { b = 0; c = -a; coefficient[1] = {b,0}; coefficient[0] = {c,0}; discriminant = 4*a*a; }
  else if(_zeroBoundaryMask != 0) {
    const long double endpoint = (_zeroBoundaryMask&1U) != 0 ? -1 : 1;
    coefficient[0] = CompensatedAdd({-coefficient[2][0],-coefficient[2][1]},CompensatedMultiply({-endpoint,0},coefficient[1]));
    c = coefficient[0][0]+coefficient[0][1];
    discriminant = QuadraticDiscriminant(coefficient);
  }
  if(_zeroReady && _identicallyZero) a = b = c = 0;
  _lineAnalysis->degree = a != 0 ? 2U : (b != 0 ? 1U : 0U);
  _lineAnalysis->identicallyZero = a == 0 && b == 0 && c == 0;
  if(!_lineAnalysis->identicallyZero) {
    std::array<long double,2> candidates;
    const unsigned count = SolveQuadratic(a,b,c,discriminant,candidates);
    const long double slope = a == 0 ? b : (discriminant > 0 ? std::sqrt(discriminant) : 0);
    const long double padding = 64*std::numeric_limits<long double>::epsilon();
    for(unsigned n = 0; n < count; ++n) if(std::isfinite(candidates[n]) && candidates[n] >= -1-padding && candidates[n] <= 1+padding) {
      LineRoot& root = _lineAnalysis->roots[_lineAnalysis->rootCount++];
      root.coordinate = std::max(-1.0L,std::min(1.0L,candidates[n]));
      root.multiplicity = a != 0 && discriminant == 0 ? 2U : 1U;
      root.derivative = _normalization*(a == 0 ? slope : std::copysign(slope,(n == 0 ? -a : a)));
      root.boundary = root.coordinate == -1 ? -1 : (root.coordinate == 1 ? 1 : 0);
    }
    long double left = -1;
    for(unsigned n = 0; n <= _lineAnalysis->rootCount; ++n) {
      const long double right = n == _lineAnalysis->rootCount ? 1 : _lineAnalysis->roots[n].coordinate;
      if(right > left) {
        const bool positive = a == 0 ? (b == 0 ? c > 0 : (b > 0 ? left >= candidates[0] : right <= candidates[0])) : (discriminant <= 0 || right <= candidates[0] || left >= candidates[1] ? a > 0 : a < 0);
        if(positive) _lineAnalysis->positiveIntervals[_lineAnalysis->positiveIntervalCount++] = {left,right};
      }
      left = right;
    }
  }
  _lineAnalysisReady = true;
  return *_lineAnalysis;
}

inline bool CutFemLevelSetData::DependsOn(unsigned axis) const {
  EnsurePower();
  if(axis >= _dimension) throw std::out_of_range("Dependence axis exceeds the active dimension.");
  if(!_dependenceReady) {
    _depends.fill(false);
    _parity.fill(1);
    std::array<bool,3> even = {true,true,true}, odd = {true,true,true};
    for(unsigned n = 0; n < _power.size(); ++n) if(_power[n] != 0) {
      const std::array<unsigned,3> index = PowerExponents(n,_dimension,2);
      for(unsigned direction = 0; direction < _dimension; ++direction) { _depends[direction] = _depends[direction] || index[direction] != 0; if(index[direction]%2U == 0) odd[direction] = false; else even[direction] = false; }
    }
    for(unsigned direction = 0; direction < _dimension; ++direction) _parity[direction] = even[direction] ? 1 : (odd[direction] ? -1 : 0);
    _dependenceReady = true;
  }
  return _depends[axis];
}

inline int CutFemLevelSetData::GetReflectionParity(unsigned axis) const {
  DependsOn(axis);
  if(_geometry == TRI || _geometry == TET || (_geometry == WEDGE && axis < 2)) return 0;
  return _parity[axis];
}

inline const CutFemLevelSetData::FiberProjection& CutFemLevelSetData::GetFiberProjection(unsigned axis) const {
  EnsurePower();
  if(_dimension < 2 || axis >= _dimension) throw std::out_of_range("Projection requires a valid axis in dimension two or three.");
  if(_projectionReady[axis]) return *_projections[axis];
  if(!_projections[axis] || !_projections[axis].unique()) _projections[axis] = std::make_shared<FiberProjection>();
  FiberProjection& projection = *_projections[axis];
  projection.dimension = _dimension-1;
  const unsigned size = PolynomialSize(projection.dimension,2);
  for(std::vector<long double>& coefficient : projection.coefficients) coefficient.assign(size,0);
  for(unsigned n = 0; n < _power.size(); ++n) {
    const std::array<unsigned,3> source = PowerExponents(n,_dimension,2);
    std::array<unsigned,3> target = {0,0,0};
    unsigned next = 0;
    for(unsigned direction = 0; direction < _dimension; ++direction) if(direction != axis) target[next++] = source[direction];
    projection.coefficients[source[axis]][PowerIndex(target,projection.dimension,2)] = _power[n]/_normalization;
  }
  projection.discriminant.assign(PolynomialSize(projection.dimension,4),0);
  for(unsigned i = 0; i < size; ++i) for(unsigned j = 0; j < size; ++j) {
    std::array<unsigned,3> index = PowerExponents(i,projection.dimension,2);
    const std::array<unsigned,3> second = PowerExponents(j,projection.dimension,2);
    for(unsigned direction = 0; direction < projection.dimension; ++direction) index[direction] += second[direction];
    projection.discriminant[PowerIndex(index,projection.dimension,4)] += projection.coefficients[1][i]*projection.coefficients[1][j]-4*projection.coefficients[2][i]*projection.coefficients[0][j];
  }
  _projectionReady[axis] = true;
  return projection;
}

inline long double CutFemLevelSetData::EvaluateFiberDiscriminant(unsigned axis, const std::vector<long double>& point) const {
  if(_dimension == 2 && point.size() == 1) return EvaluateFiberDiscriminant(axis,point[0]);
  const FiberProjection& projection = GetFiberProjection(axis);
  if(point.size() != projection.dimension) throw std::invalid_argument("Discriminant point has an incorrect dimension.");
  std::array<std::array<long double,2>,3> coordinates = {};
  for(unsigned direction = 0; direction < projection.dimension; ++direction) { if(!std::isfinite(point[direction])) throw std::invalid_argument("Discriminant coordinates must be finite."); coordinates[direction] = {point[direction],0}; }
  std::array<std::array<long double,2>,3> coefficient;
  for(unsigned power = 0; power < 3; ++power) coefficient[power] = CompensatedEvaluate(projection.coefficients[power].data(),projection.dimension,2,coordinates);
  return QuadraticDiscriminant(coefficient);
}

inline long double CutFemLevelSetData::EvaluateFiberDiscriminant(unsigned axis, long double point) const {
  if(_dimension != 2 || !std::isfinite(point)) throw std::invalid_argument("A scalar discriminant requires a finite planar coordinate.");
  if(axis >= 2) throw std::out_of_range("Discriminant axis exceeds the planar dimension.");
  if(_zeroReady && _identicallyZero) return 0;
  if(_source && _restrictionSteps != 0) {
    std::array<long double,3> coordinate = {0,0,0};
    coordinate[1-axis] = point;
    std::array<std::array<long double,2>,3> coefficient;
    GetCompensatedFiberCoefficients(axis,coordinate,coefficient);
    return QuadraticDiscriminant(coefficient);
  }
  const FiberProjection& projection = GetFiberProjection(axis);
  const std::array<std::array<long double,2>,3> coordinates = {{{point,0},{0,0},{0,0}}};
  std::array<std::array<long double,2>,3> coefficient;
  for(unsigned power = 0; power < 3; ++power) coefficient[power] = CompensatedEvaluate(projection.coefficients[power].data(),1,2,coordinates);
  return QuadraticDiscriminant(coefficient);
}

inline void CutFemLevelSetData::GetFiberBoundaryPolynomial(unsigned axis, long double constant, long double slope, std::vector<long double>& coefficient) const {
  if(_dimension != 2) throw std::logic_error("An affine boundary polynomial requires planar data.");
  if(!std::isfinite(constant) || !std::isfinite(slope)) throw std::invalid_argument("Boundary chart values must be finite.");
  const FiberProjection& projection = GetFiberProjection(axis);
  coefficient.assign(5,0);
  for(unsigned power = 0; power < 3; ++power) for(unsigned outer = 0; outer < 3; ++outer) for(unsigned exponent = 0; exponent <= power; ++exponent) coefficient[outer+exponent] += projection.coefficients[power][outer]*Binomial(power,exponent)*Power(constant,power-exponent)*Power(slope,exponent);
}

inline void CutFemLevelSetData::GetShiftedFiberDiscriminant(unsigned axis, long double endpoint, std::array<long double,5>& coefficient) const {
  if(_dimension != 2 || !std::isfinite(endpoint)) throw std::invalid_argument("A shifted discriminant requires a finite planar endpoint.");
  const FiberProjection& projection = GetFiberProjection(axis);
  coefficient.fill(0);
  for(unsigned power = 0; power < 5; ++power) for(unsigned old = power; old < 5; ++old) coefficient[power] += projection.discriminant[old]*Binomial(old,power)*Power(endpoint,old-power);
  coefficient[0] = EvaluateFiberDiscriminant(axis,endpoint);
}

inline int CutFemLevelSetData::ZeroFiberStatus(unsigned axis) const {
  if(_dimension == 1) return IsIdenticallyZero() ? 1 : -1;
  const FiberProjection& projection = GetFiberProjection(axis);
  if(_dimension != 2) {
    for(const std::vector<long double>& coefficient : projection.coefficients) { bool constant = coefficient[0] != 0; for(unsigned n = 1; n < coefficient.size(); ++n) constant = constant && coefficient[n] == 0; if(constant) return -1; }
    return IsIdenticallyZero() ? 1 : 0;
  }
  unsigned selected = 0;
  long double largest = 0;
  for(unsigned power = 0; power < 3; ++power) {
    const std::vector<long double>& coefficient = projection.coefficients[power];
    if(coefficient[0] != 0 && coefficient[1] == 0 && coefficient[2] == 0) return -1;
    long double norm = 0;
    for(long double value : coefficient) norm += std::abs(value);
    if(norm > largest) { selected = power; largest = norm; }
  }
  if(largest == 0) return 1;
  const std::vector<long double>& selectedCoefficient = projection.coefficients[selected];
  const std::array<std::array<long double,2>,3> polynomial = {{{selectedCoefficient[0],0},{selectedCoefficient[1],0},{selectedCoefficient[2],0}}};
  std::array<long double,2> roots;
  const unsigned count = SolveQuadratic(selectedCoefficient[2],selectedCoefficient[1],selectedCoefficient[0],QuadraticDiscriminant(polynomial),roots);
  const long double lower = _geometry == TRI ? 0 : -1;
  bool uncertain = false;
  const long double padding = 64*std::numeric_limits<long double>::epsilon();
  for(unsigned n = 0; n < count; ++n) if(roots[n] >= lower-padding && roots[n] <= 1+padding) {
    bool excluded = false, exact = true;
    const std::array<std::array<long double,2>,3> point = {{{roots[n],0},{0,0},{0,0}}};
    for(const std::vector<long double>& coefficient : projection.coefficients) {
      const std::array<long double,2> value = CompensatedEvaluate(coefficient.data(),1,2,point);
      const long double magnitude = std::abs(value[0]+value[1]);
      excluded = excluded || magnitude > 256*std::numeric_limits<long double>::epsilon()*(std::abs(coefficient[0])+std::abs(coefficient[1])+std::abs(coefficient[2]));
      exact = exact && value[0] == 0 && value[1] == 0;
    }
    if(!excluded && exact && roots[n] >= lower && roots[n] <= 1) return 1;
    if(!excluded) uncertain = true;
  }
  return uncertain ? 0 : -1;
}

inline int CutFemLevelSetData::GetZeroFiberStatus(unsigned axis) const {
  RequireInitialized();
  if(axis >= _dimension) throw std::out_of_range("Fiber axis exceeds the active dimension.");
  if(!_zeroFiberReady[axis]) { _zeroFiberStatus[axis] = ZeroFiberStatus(axis); _zeroFiberReady[axis] = true; }
  return _zeroFiberStatus[axis];
}

inline const CutFemLevelSetData::GraphCertificate& CutFemLevelSetData::GetGraphCertificate(unsigned axis) const { return PrepareGraphCertificate(axis,true); }

inline const CutFemLevelSetData::GraphCertificate& CutFemLevelSetData::PrepareGraphCertificate(unsigned axis, bool includeWeak) const {
  if(axis >= _dimension) throw std::out_of_range("Graph axis exceeds the active dimension.");
  if(_graphCertificateReady[axis] && (!includeWeak || _graphWeakReady[axis] || (*_graphCertificates)[axis].strict)) return (*_graphCertificates)[axis];
  if(!_graphCertificates) _graphCertificates = std::make_shared<std::array<GraphCertificate,3>>();
  else if(!_graphCertificates.unique()) _graphCertificates = std::make_shared<std::array<GraphCertificate,3>>(*_graphCertificates);
  GraphCertificate& certificate = (*_graphCertificates)[axis];
  certificate = GraphCertificate();
  certificate.bounds = GetDerivativeBounds(axis);
  if(certificate.bounds[0] > 0 || certificate.bounds[1] < 0) { certificate.sign = certificate.bounds[0] > 0 ? 1 : -1; certificate.strict = true; certificate.zeroFiber = -1; }
  else if(_geometry == QUAD) {
    EnsurePower();
    long double minimum = std::numeric_limits<long double>::infinity(), maximum = -minimum;
    for(unsigned side = 0; side < 2; ++side) {
      const long double fixed = side == 0 ? -1 : 1;
      std::array<long double,3> coefficient;
      for(unsigned power = 0; power < 3; ++power) coefficient[power] = axis == 0 ? _power[3+power]+2*fixed*_power[6+power] : _power[3*power+1]+2*fixed*_power[3*power+2];
      const long double margin = 512*std::numeric_limits<long double>::epsilon()*(std::abs(coefficient[0])+std::abs(coefficient[1])+std::abs(coefficient[2]));
      for(unsigned endpoint = 0; endpoint < 2; ++endpoint) { const long double point = endpoint == 0 ? -1 : 1; long double value = EvaluateQuadratic(coefficient.data(),1,point); if(margin > 0 && std::abs(value) <= margin) { const std::array<long double,2> exact = CompensatedEvaluate(coefficient.data(),1,2,{{{point,0},{0,0},{0,0}}}); value = exact[0]+exact[1]; } minimum = std::min(minimum,value); maximum = std::max(maximum,value); }
      if(coefficient[2] != 0) { const long double stationary = -coefficient[1]/(2*coefficient[2]); if(stationary > -1 && stationary < 1) { long double value = EvaluateQuadratic(coefficient.data(),1,stationary); if(margin > 0 && std::abs(value) <= margin) value = -QuadraticDiscriminant({{{coefficient[0],0},{coefficient[1],0},{coefficient[2],0}}})/(4*coefficient[2]); minimum = std::min(minimum,value); maximum = std::max(maximum,value); } }
    }
    const long double padding = 4*Padding(_power)+_bernsteinPadding;
    certificate.bounds = {minimum-padding,maximum+padding};
    certificate.sign = certificate.bounds[0] > 0 ? 1 : (certificate.bounds[1] < 0 ? -1 : 0);
    certificate.strict = certificate.sign != 0;
    if(includeWeak && !certificate.strict && ((minimum >= 0 && maximum > padding) || (maximum <= 0 && minimum < -padding))) certificate.sign = maximum > padding ? 1 : -1;
    if(certificate.sign != 0) certificate.zeroFiber = certificate.strict ? -1 : GetZeroFiberStatus(axis);
  }
  else if(includeWeak) {
    EnsureBernstein();
    const unsigned degree = _bernsteinOrders[axis];
    const bool simplex = _geometry == TRI || _geometry == TET || (_geometry == WEDGE && axis < 2);
    long double minimum = std::numeric_limits<long double>::infinity(), maximum = -minimum;
    for(const std::array<unsigned,3>& index : _bernsteinIndices) {
      if(simplex ? index[0]+index[1]+(_geometry == TET ? index[2] : 0) >= degree : index[axis] >= degree) continue;
      std::array<unsigned,3> next = index;
      ++next[axis];
      const long double value = _bernstein[BernsteinIndex(_geometry,_bernsteinOrders,next)]-_bernstein[BernsteinIndex(_geometry,_bernsteinOrders,index)];
      minimum = std::min(minimum,value);
      maximum = std::max(maximum,value);
    }
    if((minimum >= 0 && maximum > 0) || (maximum <= 0 && minimum < 0)) { certificate.sign = maximum > 0 ? 1 : -1; certificate.zeroFiber = GetZeroFiberStatus(axis); }
  }
  _graphCertificateReady[axis] = true;
  _graphWeakReady[axis] = includeWeak;
  return certificate;
}

inline int CutFemLevelSetData::GetGraphDirection(bool allowWeak) const {
  const int strict = GetGraphDirection();
  if(strict >= 0 || !allowWeak) return strict;
  for(unsigned axis = 0; axis < _dimension; ++axis) { const GraphCertificate& certificate = GetGraphCertificate(axis); if(certificate.sign != 0 && certificate.zeroFiber == -1) return static_cast<int>(axis); }
  return -1;
}

#endif
