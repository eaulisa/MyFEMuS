#include "NativeCompatibility.hpp"

namespace Ex47NativeTests {

struct OperationTimer { const char* method; WorkCounters before; std::chrono::steady_clock::time_point start; OperationTimer(const char* name) : method(name), before(Counters()), start(std::chrono::steady_clock::now()) {} ~OperationTimer() { Record(method,start,true,&before); } };

void LineIntegration::FullStripMoments(unsigned degree, const FELevelSetMoments::LineIntegrationStrip& strip, long double left, long double right, std::vector<long double>& values, std::vector<long double>& errors, const std::array<long double,5>* coordinates, const std::vector<unsigned char>* requested) { const FELevelSetMoments::Ex47Prepared::LineIntegrationStrip local = {strip.left,strip.right,strip.lowerConstant,strip.lowerSlope,strip.upperConstant,strip.upperSlope,strip.fullBoundaryDelta,strip.certifiedLowerZero,strip.certifiedUpperZero}; const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now(); FELevelSetMoments::Ex47Prepared::LineIntegration::FullStripMoments(degree,local,left,right,values,errors,coordinates,requested); Record("full_strip",start,true); }

struct CompatibilityWork {
  CutFemLevelSetData data;
  CutFemLevelSetData embedding;
  FELevelSetMoments::Ex47Prepared::PreparedMomentWorkspace workspace;
  std::vector<std::vector<long double>> flatMetric;
  unsigned metricDegree = 2;
  unsigned compared = 0;
};

CompatibilityWork& Work(GeomElType geometry) { static thread_local std::array<CompatibilityWork,6> work; return work[static_cast<unsigned>(geometry)]; }

unsigned Dimension(GeomElType geometry) { return geometry == LINE ? 1U : (geometry == QUAD || geometry == TRI ? 2U : 3U); }

unsigned MomentCount(unsigned dimension, unsigned degree) { return dimension == 1 ? degree+1 : (dimension == 2 ? (degree+1)*(degree+2)/2 : (degree+1)*(degree+2)*(degree+3)/6); }

template<class Real> std::vector<long double> Flatten(const std::vector<std::vector<Real>>& input) { if(input.size() != 3) throw std::invalid_argument("A planar polynomial requires three rows."); std::vector<long double> output(9); for(unsigned i = 0; i < 3; ++i) { if(input[i].size() != 3) throw std::invalid_argument("A planar polynomial requires three columns."); for(unsigned j = 0; j < 3; ++j) { output[3*i+j] = input[i][j]; if(!std::isfinite(output[3*i+j])) throw std::invalid_argument("Nonfinite polynomial coefficient."); } } return output; }

template<class Real> std::vector<long double> Flatten(const std::vector<std::vector<std::vector<Real>>>& input) { if(input.size() != 3) throw std::invalid_argument("A spatial polynomial requires three planes."); std::vector<long double> output(27); for(unsigned i = 0; i < 3; ++i) { if(input[i].size() != 3) throw std::invalid_argument("Invalid polynomial rows."); for(unsigned j = 0; j < 3; ++j) { if(input[i][j].size() != 3) throw std::invalid_argument("Invalid polynomial columns."); for(unsigned k = 0; k < 3; ++k) { output[9*i+3*j+k] = input[i][j][k]; if(!std::isfinite(output[9*i+3*j+k])) throw std::invalid_argument("Nonfinite coefficient."); } } } return output; }

std::vector<std::vector<long double>> Matrix(const std::vector<long double>& input, unsigned side = 3) { std::vector<std::vector<long double>> output(side,std::vector<long double>(side)); for(unsigned i = 0; i < side; ++i) for(unsigned j = 0; j < side; ++j) output[i][j] = input[side*i+j]; return output; }

void Validate(unsigned degree, long double absolute, long double relative, unsigned direction, unsigned maximumDirection, unsigned depth, unsigned budget) { if(degree > 64 || !(absolute > 0) || !(relative >= 0) || !std::isfinite(absolute) || !std::isfinite(relative) || direction > maximumDirection || depth > 60 || budget < 3) throw std::invalid_argument("Invalid prepared integration request."); }

void ResetPolynomial(CompatibilityWork& work, GeomElType geometry, const std::vector<long double>& coefficient) { if(!work.data.IsInitialized() || !work.data.GetNodalValues().empty() || work.data.GetGeometry() != geometry || work.data.GetPowerCoefficients() != coefficient) work.data.ResetPowerCoefficients(geometry,coefficient); }

std::vector<long double> AffinePolynomial(const std::vector<long double>& input, unsigned dimension, unsigned degree, unsigned axes, long double origin, long double scale) {
  std::vector<long double> output = input;
  const unsigned side = degree+1;
  for(unsigned axis = 0; axis < dimension; ++axis) if((axes&(1U<<axis)) != 0) {
    unsigned stride = 1;
    for(unsigned next = axis+1; next < dimension; ++next) stride *= side;
    std::vector<long double> transformed(output.size(),0);
    for(unsigned n = 0; n < output.size(); ++n) if(output[n] != 0) { const unsigned exponent = (n/stride)%side, base = n-exponent*stride; for(unsigned power = 0; power <= exponent; ++power) transformed[base+power*stride] += output[n]*FELevelSetMoments::ReferenceBinomial(exponent,power)*FELevelSetMoments::ReferencePower(origin,exponent-power)*FELevelSetMoments::ReferencePower(scale,power); }
    output.swap(transformed);
  }
  return output;
}

void ReferenceOutput(unsigned dimension, unsigned degree, const std::vector<long double>& input, const std::vector<long double>& inputErrors, std::vector<long double>& output, std::vector<long double>& errors) {
  const unsigned count = input.size()/2;
  output.assign(input.size(),0);
  errors.assign(input.size(),0);
  for(unsigned pair = 0; pair < 2; ++pair) for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= (dimension > 1 ? degree-i : 0U); ++j) for(unsigned k = 0; k <= (dimension > 2 ? degree-i-j : 0U); ++k) {
    const unsigned target = pair*count+(dimension == 2 ? FELevelSetMoments::ReferenceMomentIndex2D(degree,i,j) : FELevelSetMoments::UnitMomentIndex(dimension,i,j,k));
    for(unsigned a = 0; a <= i; ++a) for(unsigned b = 0; b <= j; ++b) for(unsigned c = 0; c <= k; ++c) {
      const unsigned source = pair*count+(dimension == 2 ? FELevelSetMoments::ReferenceMomentIndex2D(degree,a,b) : FELevelSetMoments::UnitMomentIndex(dimension,a,b,c));
      const long double factor = (1U<<dimension)*FELevelSetMoments::ReferenceBinomial(i,a)*FELevelSetMoments::ReferenceBinomial(j,b)*FELevelSetMoments::ReferenceBinomial(k,c)*FELevelSetMoments::ReferencePower(2,a+b+c)*((i+j+k-a-b-c)%2U == 0 ? 1 : -1);
      output[target] += factor*input[source];
      errors[target] += std::abs(factor)*inputErrors[source];
    }
  }
}

void Metric(const CutFemLevelSetData& data, bool surface, const std::vector<std::vector<double>>* inverse, const std::vector<std::vector<double>>* geometry, CompatibilityWork& work, std::vector<std::vector<long double>>& metric, unsigned& degree) {
  degree = 2;
  metric.clear();
  if(!surface && (inverse != nullptr || geometry != nullptr)) throw std::invalid_argument("A physical metric requires surface integration.");
  if(inverse != nullptr && geometry != nullptr) throw std::invalid_argument("Provide either nodal geometry or an inverse Jacobian.");
  if(geometry != nullptr) { FELevelSetMoments::ValidateSurfaceGeometryCoordinates(*geometry,data.GetDimension(),data.GetGeometry() == LINE ? 3U : (data.GetGeometry() == TRI ? 7U : (data.GetGeometry() == QUAD ? 9U : (data.GetGeometry() == TET ? 15U : (data.GetGeometry() == WEDGE ? 21U : 27U))))); FELevelSetMoments::Ex47Prepared::PrepareMetric(data,*geometry,work.workspace); metric = work.workspace.metric; degree = work.workspace.metricDegree; }
  if(inverse != nullptr) { FELevelSetMoments::ValidateSurfaceInverseJacobian(*inverse,data.GetDimension()); metric.assign(data.GetDimension(),std::vector<long double>(data.GetPowerCoefficients().size(),0)); for(unsigned i = 0; i < data.GetDimension(); ++i) for(unsigned j = 0; j < data.GetDimension(); ++j) { const std::vector<long double>& gradient = data.GetRootDerivativeCoefficients(j); for(unsigned n = 0; n < gradient.size(); ++n) metric[i][n] += (*inverse)[j][i]*gradient[n]; } }
}

void FlattenMetric(const std::vector<std::vector<std::vector<long double>>>& input, std::vector<std::vector<long double>>& output, unsigned& degree) { degree = 0; for(const std::vector<std::vector<long double>>& component : input) { degree = std::max(degree,static_cast<unsigned>(component.size()-1)); for(const std::vector<long double>& row : component) degree = std::max(degree,static_cast<unsigned>(row.size()-1)); } output.assign(input.size(),std::vector<long double>((degree+1)*(degree+1),0)); for(unsigned c = 0; c < input.size(); ++c) for(unsigned i = 0; i < input[c].size(); ++i) for(unsigned j = 0; j < input[c][i].size(); ++j) output[c][(degree+1)*i+j] = input[c][i][j]; }

void PlanarData(const CutFemLevelSetData& source, unsigned degree, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned depth, unsigned budget, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned& selected, const std::vector<std::vector<long double>>* metric, unsigned metricDegree, const std::array<long double,5>* parent, const std::vector<unsigned char>* requested, unsigned direction, bool fullBoundary, const FELevelSetMoments::OuterQuadratureRules& rules, bool inverseSurface, int normalComponent, CompatibilityWork& work, unsigned validation = 0) {
  rules.Validate();
  Validate(degree,absoluteH,relative,direction,2,depth,budget);
  if(!(absoluteD > 0) || !std::isfinite(absoluteD)) throw std::invalid_argument("Invalid delta tolerance.");
  const CutFemLevelSetData* field = &source;
  if(metric != nullptr && metric->size() > source.GetRootDimension()) { std::vector<long double> embedded(27,0); const std::vector<long double>& power = source.GetPowerCoefficients(); for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) embedded[9*i+3*j] = power[3*i+j]; work.embedding.ResetPowerCoefficients(source.GetGeometry() == TRI ? TET : HEX,embedded); work.embedding.Restrict(2,0,work.workspace.chart,false); field = &work.workspace.chart; }
  const long double normalization = field->GetNormalizationScale();
  const bool surface = metric != nullptr;
  work.workspace.plane.Bind(*field,normalization,surface,inverseSurface,normalComponent,metric,metricDegree,&work.workspace.statistics);
  if(surface) FELevelSetMoments::Ex47Prepared::PlanarPolicy(*field,metric,metricDegree,work.workspace.policy);
  const std::vector<unsigned> boundary = work.workspace.plane.Boundaries();
  const std::vector<std::vector<long double>> normalized = FELevelSetMoments::Ex47Prepared::PlanarCoefficients(*field,normalization);
  const WorkCounters before = Counters();
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  try { FELevelSetMoments::Ex47Prepared::PlanarEndpointIntegration::Integrate(normalized,degree,source.GetGeometry() == TRI,computeDelta,absoluteH,absoluteD*normalization,relative,depth,budget,values,errors,evaluations,accepted,selected,surface ? &work.workspace.policy : nullptr,parent,requested,direction,fullBoundary,&boundary,rules,inverseSurface,normalComponent,&work.workspace.plane); } catch(...) { Record("planar_kernel",start,false,&before); throw; }
  work.compared = 0;
  if(validation == 2) work.compared = FELevelSetMoments::Ex47Prepared::PlanarEndpointIntegration::Compare(normalized,degree,source.GetGeometry() == TRI,computeDelta,absoluteH,absoluteD*normalization,relative,depth,budget,values,errors,evaluations,accepted,selected,surface ? &work.workspace.policy : nullptr,parent,fullBoundary,&boundary,rules,inverseSurface,normalComponent,&work.workspace.plane);
  Record("planar_kernel",start,true,&before);
  const unsigned count = values.size()/2;
  for(unsigned n = count; n < values.size(); ++n) { values[n] /= normalization; errors[n] /= normalization; }
}

void SpatialData(const CutFemLevelSetData& source, unsigned degree, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned depth, unsigned budget, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<long double>>* metric, unsigned metricDegree, unsigned inner, unsigned outer, bool computePositive, const FELevelSetMoments::OuterQuadratureRules& rules, bool inverseSurface, int normalComponent, CompatibilityWork& work, unsigned validation = 0) {
  Validate(degree,absoluteH,relative,0,0,depth,budget);
  FELevelSetMoments::Ex47Prepared::SurfaceEndpointData& data = work.workspace.surface;
  const long double normalization = source.GetNormalizationScale();
  data.prepared = &source;
  data.normalization = normalization;
  data.coefficient = source.GetNormalizedPowerCoefficients();
  data.referenceSurface = metric != nullptr;
  data.physicalMetric = metric;
  data.metricDegree = metricDegree;
  data.gradient = metric == nullptr ? std::vector<std::vector<long double>>() : *metric;
  data.statistics = &work.workspace.statistics;
  data.degree = degree;
  data.maximumDepth = depth;
  data.maximumEvaluations = budget;
  data.evaluations = data.accepted = 0;
  data.graphAxis = 3;
  data.preferredInner = inner;
  data.preferredOuter = outer;
  data.innerCounts.fill(0);
  data.outerCounts.fill(0);
  data.retainUncertified = false;
  data.wedge = source.GetGeometry() == WEDGE;
  data.computePositive = computePositive;
  data.computeDelta = computeDelta;
  data.inverseSurface = inverseSurface;
  data.normalComponent = normalComponent;
  data.absoluteH = absoluteH;
  data.absoluteD = absoluteD*normalization;
  data.relative = relative;
  data.quadratureRules = rules;
  data.zeroFaces.resize(source.GetBoundaryCount());
  for(unsigned n = 0; n < data.zeroFaces.size(); ++n) data.zeroFaces[n] = source.IsZeroOnBoundary(n);
  const WorkCounters before = Counters();
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  try { FELevelSetMoments::Ex47Prepared::SurfaceEndpointIntegration::Integrate(data,source.GetGeometry() == TET,values,errors); } catch(...) { Record("spatial_kernel",start,false,&before); throw; }
  if(validation == 2) FELevelSetMoments::Ex47Prepared::SurfaceEndpointIntegration::Compare(data,source.GetGeometry() == TET,values,errors);
  Record("spatial_kernel",start,true,&before);
  const unsigned count = values.size()/2;
  for(unsigned n = count; n < values.size(); ++n) { values[n] /= normalization; errors[n] /= normalization; }
}

void OutputPlanar(unsigned degree, const std::vector<long double>& values, const std::vector<long double>& errors, long double absolute, long double relative, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, std::vector<double>* estimate, unsigned evaluations, unsigned accepted, unsigned selected) {
  if(&positive == &delta) throw std::invalid_argument("Moment outputs must be distinct.");
  const unsigned count = (degree+1)*(degree+2)/2;
  std::vector<std::vector<double>> h(degree+1), d(degree+1);
  long double maximumH = 0, maximumD = 0;
  for(unsigned i = 0; i <= degree; ++i) { h[i].resize(degree-i+1); d[i].resize(degree-i+1); for(unsigned j = 0; j <= degree-i; ++j) { const unsigned n = FELevelSetMoments::ReferenceMomentIndex2D(degree,i,j); h[i][j] = values[n]; d[i][j] = values[count+n]; if(!std::isfinite(h[i][j]) || !std::isfinite(d[i][j])) throw std::overflow_error("Moment output exceeds double range."); const long double eh = errors[n]+std::abs(values[n]-static_cast<long double>(h[i][j])), ed = errors[count+n]+std::abs(values[count+n]-static_cast<long double>(d[i][j])); maximumH = std::max(maximumH,eh); maximumD = std::max(maximumD,ed); if(eh > absolute+relative*std::abs(values[n])) FELevelSetMoments::WarnQuadrature("Planar output","rounded volume moment exceeds tolerance",n,eh,absolute+relative*std::abs(values[n])); if(ed > absolute+relative*std::abs(values[count+n])) FELevelSetMoments::WarnQuadrature("Planar output","rounded interface moment exceeds tolerance",count+n,ed,absolute+relative*std::abs(values[count+n])); } }
  positive.swap(h);
  delta.swap(d);
  if(estimate != nullptr) *estimate = {static_cast<double>(maximumH),static_cast<double>(maximumD),static_cast<double>(accepted),static_cast<double>(evaluations),static_cast<double>(selected),0,1};
}

void OutputSpatial(unsigned degree, bool computeDelta, const std::vector<long double>& values, const std::vector<long double>& errors, long double absolute, long double relative, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, std::vector<double>* estimate, const CompatibilityWork& work) { FELevelSetMoments::SurfaceEndpointData output; output.degree = degree; output.computeDelta = computeDelta; output.absoluteH = output.absoluteD = absolute; output.relative = relative; output.powers = work.workspace.surface.powers; output.evaluations = work.workspace.surface.evaluations; output.accepted = work.workspace.surface.accepted; output.outerCounts = work.workspace.surface.outerCounts; FELevelSetMoments::SurfaceEndpointIntegration::Output(output,1,values,errors,positive,delta,estimate); }

void SelectPlanar(const CutFemLevelSetData& source, unsigned degree, bool computeDelta, bool fullBoundary, unsigned direction, unsigned validation, unsigned depth, unsigned budget, long double absolute, long double relative, bool surface, const std::vector<std::vector<double>>* inverse, const std::vector<std::vector<double>>* geometry, const std::vector<std::vector<std::vector<long double>>>* gradient, const FELevelSetMoments::OuterQuadratureRules& rules, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, std::vector<double>* estimate, CompatibilityWork& work, bool inverseSurface = false, int normalComponent = -1) {
  FELevelSetMoments::QuadratureWarningScope scope("Prepared planar compatibility");
  if(&positive == &delta || direction > (source.GetGeometry() == TRI ? 3U : 2U) || validation > 2) throw std::invalid_argument("Invalid selected planar request.");
  std::vector<std::vector<long double>> metric;
  unsigned metricDegree = 2;
  Metric(source,surface,inverse,geometry,work,metric,metricDegree);
  if(gradient != nullptr) { FlattenMetric(*gradient,metric,metricDegree); if(source.GetGeometry() == QUAD) for(std::vector<long double>& component : metric) component = AffinePolynomial(component,2,metricDegree,3,-1,2); }
  else if(surface && metric.empty()) { metric.resize(source.GetRootDimension()); for(unsigned c = 0; c < metric.size(); ++c) metric[c] = source.GetRootDerivativeCoefficients(c); }
  if(source.GetGeometry() == QUAD && gradient == nullptr) for(std::vector<long double>& component : metric) for(long double& value : component) value /= geometry != nullptr ? 4 : 2;
  const CutFemLevelSetData* field = &source;
  const bool barycentric = source.GetGeometry() == TRI && direction == 3;
  if(barycentric) { source.PermuteSimplexVertices({1,0,2,3},work.workspace.chart); field = &work.workspace.chart; if(!metric.empty()) { unsigned transformedDegree = 0; field->TransformParentComponents(metricDegree,metric,transformedDegree,work.flatMetric); metric = work.flatMetric; metricDegree = transformedDegree; } }
  const long double amplification = barycentric ? FELevelSetMoments::ReferencePower(3,degree) : 1;
  std::vector<long double> values, errors;
  unsigned evaluations = 0, accepted = 0, selected = 0;
  const std::array<long double,5> parent = source.GetGeometry() == QUAD ? std::array<long double,5>{{-1,2,-1,2,4}} : std::array<long double,5>{{0,1,0,1,1}};
  PlanarData(*field,degree,computeDelta,absolute/amplification,absolute/amplification,relative/amplification,depth,budget,values,errors,evaluations,accepted,selected,metric.empty() ? nullptr : &metric,metricDegree,&parent,nullptr,barycentric ? 1U : direction,fullBoundary,rules,inverseSurface,normalComponent,work,validation);
  if(barycentric) { std::vector<long double> original, originalErrors; FELevelSetMoments::EndpointSimplexChart::RestoreTriangle(degree,values,errors,original,originalErrors); values.swap(original); errors.swap(originalErrors); }
  OutputPlanar(degree,values,errors,absolute,relative,positive,delta,estimate,evaluations,accepted,selected);
  if(estimate != nullptr) (*estimate)[5] = work.compared;
}

void SelectSpatial(const CutFemLevelSetData& source, unsigned degree, bool computeDelta, unsigned direction, unsigned validation, unsigned depth, unsigned budget, long double absolute, long double relative, bool surface, const std::vector<std::vector<double>>* inverse, const std::vector<std::vector<double>>* geometry, const FELevelSetMoments::OuterQuadratureRules& rules, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, std::vector<double>* estimate, CompatibilityWork& work, bool inverseSurface = false, int normalComponent = -1, bool computePositive = true) {
  FELevelSetMoments::QuadratureWarningScope scope("Prepared spatial compatibility");
  if(&positive == &delta || direction > (source.GetGeometry() == TET ? 4U : 6U) || validation > 2) throw std::invalid_argument("Invalid selected spatial request.");
  std::vector<std::vector<long double>> metric;
  unsigned metricDegree = 2;
  Metric(source,surface,inverse,geometry,work,metric,metricDegree);
  if(surface && metric.empty()) { metric.resize(3); for(unsigned c = 0; c < 3; ++c) metric[c] = source.GetRootDerivativeCoefficients(c); }
  if(source.GetGeometry() == HEX) for(std::vector<long double>& component : metric) for(long double& value : component) value /= geometry != nullptr ? 8 : 2;
  const CutFemLevelSetData* field = &source;
  const bool barycentric = source.GetGeometry() == TET && direction == 1;
  if(barycentric) { source.PermuteSimplexVertices({1,0,2,3},work.workspace.chart); field = &work.workspace.chart; if(!metric.empty()) { unsigned transformedDegree = 0; field->TransformParentComponents(metricDegree,metric,transformedDegree,work.flatMetric); metric = work.flatMetric; metricDegree = transformedDegree; } }
  unsigned inner = 3, outer = 3;
  if(direction != 0) { if(source.GetGeometry() == TET) outer = direction == 1 ? 0U : direction-2; else { unsigned middle; FELevelSetMoments::Hex27SelectedAxes(direction,inner,middle,outer); } }
  const long double amplification = barycentric ? FELevelSetMoments::ReferencePower(4,degree) : 1;
  std::vector<long double> values, errors;
  SpatialData(*field,degree,computeDelta,absolute/amplification,absolute/amplification,relative/amplification,depth,budget,values,errors,metric.empty() ? nullptr : &metric,metricDegree,inner,outer,computePositive,rules,inverseSurface,normalComponent,work,validation);
  if(barycentric) { std::vector<long double> original, originalErrors; FELevelSetMoments::EndpointSimplexChart::RestoreTetrahedron(degree,values,errors,original,originalErrors); values.swap(original); errors.swap(originalErrors); }
  if(source.GetGeometry() == HEX) { std::vector<long double> native, nativeErrors; ReferenceOutput(3,degree,values,errors,native,nativeErrors); values.swap(native); errors.swap(nativeErrors); }
  OutputSpatial(degree,computeDelta,values,errors,absolute,relative,positive,delta,estimate,work);
}

void IntegrateQuad9NodalSelected(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance, double relativeTolerance, bool computeDelta, bool fullBoundaryDelta, unsigned validationMode, unsigned integrationOrder, std::vector<double>* estimatedError, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<double>>* geometryCoordinates, const FELevelSetMoments::OuterQuadratureRules& rules) {
  CompatibilityWork& work = Work(QUAD); work.data.Reset(QUAD,nodalValues); SelectPlanar(work.data,degree,computeDelta,fullBoundaryDelta,integrationOrder,validationMode,28,200000,absoluteTolerance,relativeTolerance,referenceSurface,inverseJacobian,geometryCoordinates,nullptr,rules,positive,delta,estimatedError,work);
}

void IntegrateTri7NodalSelected(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance, double relativeTolerance, bool computeDelta, bool fullBoundaryDelta, unsigned validationMode, unsigned integrationDirection, std::vector<double>* estimatedError, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<double>>* geometryCoordinates, const FELevelSetMoments::OuterQuadratureRules& rules, bool inverseSurface, int normalComponent) {
  CompatibilityWork& work = Work(TRI); work.data.Reset(TRI,nodalValues); SelectPlanar(work.data,degree,computeDelta,fullBoundaryDelta,integrationDirection,validationMode,28,200000,absoluteTolerance,relativeTolerance,referenceSurface,inverseJacobian,geometryCoordinates,nullptr,rules,positive,delta,estimatedError,work,inverseSurface,normalComponent);
}

void IntegrateQuad9ReferenceSelected(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance, double relativeTolerance, bool computeDelta, bool fullBoundaryDelta, unsigned maximumDepth, unsigned integrationOrder, unsigned maximumEvaluations, std::vector<double>* estimatedError, unsigned validationMode, const std::vector<unsigned>* zeroBoundary, bool referenceSurface, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<double>>* geometryCoordinates, const FELevelSetMoments::OuterQuadratureRules& rules) {
  CompatibilityWork& work = Work(QUAD); ResetPolynomial(work,QUAD,AffinePolynomial(Flatten(coefficient),2,2,3,-1,2)); SelectPlanar(work.data,degree,computeDelta,fullBoundaryDelta,integrationOrder,validationMode,maximumDepth,maximumEvaluations,absoluteTolerance,relativeTolerance,referenceSurface,inverseJacobian,geometryCoordinates,surfaceGradient,rules,positive,delta,estimatedError,work);
}

void IntegrateTri7ReferenceSelected(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance, double relativeTolerance, bool computeDelta, bool fullBoundaryDelta, unsigned maximumDepth, unsigned integrationDirection, unsigned maximumEvaluations, std::vector<double>* estimatedError, unsigned validationMode, const std::vector<unsigned>* zeroBoundary, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<double>>* geometryCoordinates, const FELevelSetMoments::OuterQuadratureRules& rules, bool inverseSurface, int normalComponent) {
  const std::vector<long double> flattened = Flatten(coefficient); if(flattened[8] != 0 || flattened[7] != flattened[5]) throw std::invalid_argument("TRI7 requires a quadratic polynomial with a barycentric bubble."); if(zeroBoundary != nullptr && zeroBoundary->size() != 3) throw std::invalid_argument("TRI7 requires three boundary flags.");
  CompatibilityWork& work = Work(TRI); ResetPolynomial(work,TRI,Flatten(coefficient)); SelectPlanar(work.data,degree,computeDelta,fullBoundaryDelta,integrationDirection,validationMode,maximumDepth,maximumEvaluations,absoluteTolerance,relativeTolerance,referenceSurface,inverseJacobian,geometryCoordinates,nullptr,rules,positive,delta,estimatedError,work,inverseSurface,normalComponent);
}

void IntegrateHex27NodalSelected(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, double absoluteTolerance, double relativeTolerance, bool computeDelta, unsigned validationMode, unsigned integrationOrder, std::vector<double>* estimatedError, std::vector<double>* directionInfo, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<double>>* geometryCoordinates, const FELevelSetMoments::OuterQuadratureRules& rules) {
  CompatibilityWork& work = Work(HEX); work.data.Reset(HEX,nodalValues); SelectSpatial(work.data,degree,computeDelta,integrationOrder,validationMode,24,200000,absoluteTolerance,relativeTolerance,referenceSurface,inverseJacobian,geometryCoordinates,rules,positive,delta,estimatedError,work); if(directionInfo != nullptr) *directionInfo = {1,1,0};
}

void IntegrateHex27ReferenceSelected(const std::vector<std::vector<std::vector<double>>>& coefficient, unsigned degree, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, double absoluteTolerance, double relativeTolerance, bool computeDelta, unsigned maximumDepth, unsigned integrationOrder, unsigned maximumEvaluations, std::vector<double>* estimatedError, unsigned validationMode, const std::vector<unsigned>* zeroFaces, std::vector<double>* directionInfo, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<double>>* geometryCoordinates, const FELevelSetMoments::OuterQuadratureRules& rules) {
  CompatibilityWork& work = Work(HEX); ResetPolynomial(work,HEX,AffinePolynomial(Flatten(coefficient),3,2,7,-1,2)); SelectSpatial(work.data,degree,computeDelta,integrationOrder,validationMode,maximumDepth,maximumEvaluations,absoluteTolerance,relativeTolerance,referenceSurface,inverseJacobian,geometryCoordinates,rules,positive,delta,estimatedError,work); if(directionInfo != nullptr) *directionInfo = {1,1,0};
}

void IntegrateTet15NodalByTriangleSlices(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, double absoluteTolerance, double relativeTolerance, bool computeDelta, unsigned maximumDepth, unsigned integrationAxis, unsigned maximumEvaluations, std::vector<double>* estimatedError, unsigned directionValidation, std::vector<double>* directionInfo, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian, const std::vector<std::vector<double>>* geometryCoordinates, const FELevelSetMoments::OuterQuadratureRules& rules, bool inverseSurface, int normalComponent, bool computePositive, std::vector<std::vector<long double>>* gradientCache) {
  CompatibilityWork& work = Work(TET); work.data.Reset(TET,nodalValues); SelectSpatial(work.data,degree,computeDelta,integrationAxis,directionValidation,maximumDepth,maximumEvaluations,absoluteTolerance,relativeTolerance,referenceSurface,inverseJacobian,geometryCoordinates,rules,positive,delta,estimatedError,work,inverseSurface,normalComponent,computePositive); if(directionInfo != nullptr) *directionInfo = {1,1,0};
}

void IntegrateUnitTensorNodalMoments(const std::vector<double>& nodal, unsigned dimension, unsigned degree, std::vector<long double>& values, std::vector<long double>& errors, long double absolute, long double relative, bool surface, const std::vector<std::vector<double>>* geometry, unsigned direction, unsigned validation, const FELevelSetMoments::OuterQuadratureRules& rules, bool inverseSurface, int normalComponent, bool computePositive, std::vector<std::vector<long double>>* gradientCache) {
  OperationTimer timer("unit_kernel");
  if(dimension < 1 || dimension > 3) throw std::invalid_argument("Invalid tensor dimension.");
  const GeomElType geometryType = dimension == 1 ? LINE : (dimension == 2 ? QUAD : HEX);
  CompatibilityWork& work = Work(geometryType);
  work.data.Reset(geometryType,nodal);
  Validate(degree,absolute,relative,direction,dimension == 1 ? 1U : (dimension == 2 ? 2U : 6U),28,200000);
  if(dimension == 1) FELevelSetMoments::Ex47Prepared::IntegrateLine(work.data,degree,surface,inverseSurface,values,errors,computePositive);
  else if(dimension == 2) FELevelSetMoments::Ex47Prepared::IntegratePlanar(work.data,degree,absolute,relative,surface,geometry,direction,validation,rules,inverseSurface,normalComponent,computePositive,work.workspace,values,errors);
  else FELevelSetMoments::Ex47Prepared::IntegrateSpatial(work.data,degree,absolute,relative,surface,geometry,direction,validation,rules,inverseSurface,normalComponent,computePositive,work.workspace,values,errors);

}

void IntegrateWedge21NodalMoments(const std::vector<double>& nodal, unsigned degree, std::vector<long double>& values, std::vector<long double>& errors, long double absolute, long double relative, bool referenceSurface, const std::vector<std::vector<double>>* geometry, unsigned direction, unsigned validation, const FELevelSetMoments::OuterQuadratureRules& rules, bool inverseSurface, int normalComponent, bool computePositive, std::vector<std::vector<long double>>* gradientCache) {
  OperationTimer timer("wedge_kernel");
  CompatibilityWork& work = Work(WEDGE); work.data.Reset(WEDGE,nodal); FELevelSetMoments::Ex47Prepared::IntegrateSpatial(work.data,degree,absolute,relative,referenceSurface,geometry,direction,validation,rules,inverseSurface,normalComponent,computePositive,work.workspace,values,errors);
}

int NodalLevelSetPhase(GeomElType geometry, const std::vector<double>& nodal) {
  OperationTimer timer("phase");
  const char* configured = std::getenv("EX47_NATIVE_DETECTION_DEPTH"); CutFemLevelSetData data(geometry,nodal,configured == nullptr ? 1U : static_cast<unsigned>(std::strtoul(configured,nullptr,10))); return data.GetPhase();
}

void NodalLevelSetBernsteinControls(GeomElType geometry, const std::vector<double>& nodal, std::vector<long double>& controls, long double& padding) {
  OperationTimer timer("bernstein");
  CutFemLevelSetData data(geometry,nodal); controls = data.GetBernsteinCoefficients(); const std::array<long double,2> bounds = data.GetBounds(); padding = std::max(*std::min_element(controls.begin(),controls.end())-bounds[0],bounds[1]-*std::max_element(controls.begin(),controls.end()));
}

void Line3ReferenceCoefficientsFromNodalValues(const std::vector<double>& nodalValues, double& a, double& b, double& c) {
  CutFemLevelSetData data(LINE,nodalValues); const std::vector<long double> power = AffinePolynomial(data.GetPowerCoefficients(),1,2,1,0.5L,0.5L); a = power[2]; b = power[1]; c = power[0];
}

void Quad9ReferenceCoefficientsFromNodalValues(const std::vector<double>& nodalValues, std::vector<std::vector<double>>& coefficient) {
  CutFemLevelSetData data(QUAD,nodalValues); const std::vector<long double> power = AffinePolynomial(data.GetPowerCoefficients(),2,2,3,0.5L,0.5L); coefficient.assign(3,std::vector<double>(3)); for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) coefficient[i][j] = power[3*i+j];
}

void Tri7ReferenceCoefficientsFromNodalValues(const std::vector<double>& nodalValues, std::vector<std::vector<double>>& coefficient) {
  CutFemLevelSetData data(TRI,nodalValues); const std::vector<long double>& power = data.GetPowerCoefficients(); coefficient.assign(3,std::vector<double>(3)); for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) coefficient[i][j] = power[3*i+j];
}

void Hex27ReferenceCoefficientsFromNodalValues(const std::vector<double>& nodalValues, std::vector<std::vector<std::vector<double>>>& coefficient) {
  CutFemLevelSetData data(HEX,nodalValues); const std::vector<long double> power = AffinePolynomial(data.GetPowerCoefficients(),3,2,7,0.5L,0.5L); coefficient.assign(3,std::vector<std::vector<double>>(3,std::vector<double>(3))); for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) coefficient[i][j][k] = power[9*i+3*j+k];
}

void Wedge21UnitPolynomialFromNodalValues(const std::vector<double>& nodal, std::vector<long double>& coefficient) {
  CutFemLevelSetData data(WEDGE,nodal); coefficient = data.GetUnitPowerCoefficients();
}

void UnitTensorPolynomialFromNodalValues(const std::vector<double>& nodal, unsigned dimension, std::vector<long double>& coefficient) {
  if(dimension < 1 || dimension > 3) throw std::invalid_argument("Invalid tensor dimension."); CutFemLevelSetData data(dimension == 1 ? LINE : (dimension == 2 ? QUAD : HEX),nodal); const std::vector<long double>& power = data.GetUnitPowerCoefficients(); coefficient.assign(27,0); for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < (dimension > 1 ? 3U : 1U); ++j) for(unsigned k = 0; k < (dimension > 2 ? 3U : 1U); ++k) coefficient[9*i+3*j+k] = power[dimension == 1 ? i : (dimension == 2 ? 3*i+j : 9*i+3*j+k)];
}

void Tet15PolynomialFromNodalValues(const std::vector<double>& nodalValues, bool duffy, std::vector<long double>& coefficient) {
  CutFemLevelSetData data(TET,nodalValues); const std::vector<long double>& power = data.GetPowerCoefficients(); coefficient.assign(125,0); for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) { if(!duffy) coefficient[25*i+5*j+k] = power[9*i+3*j+k]; else if(power[9*i+3*j+k] != 0) for(unsigned a = 0; a <= i; ++a) for(unsigned b = 0; b <= j; ++b) coefficient[25*(i+j+k)+5*(a+j+k)+b+k] += power[9*i+3*j+k]*FELevelSetMoments::ReferenceBinomial(i,a)*FELevelSetMoments::ReferenceBinomial(j,b)*((a+b)%2 == 0 ? 1 : -1); }
}

void Integrate1DInternal(long double a, long double b, long double c, unsigned degree, long double left, long double right, bool computeDelta, bool fullBoundaryDelta, std::vector<long double>& positive, std::vector<long double>& delta) {
  if(&positive == &delta || degree > 64 || !(right > left)) throw std::invalid_argument("Invalid LINE integration request.");
  const long double origin = left, scale = right-left;
  CompatibilityWork& work = Work(LINE);
  ResetPolynomial(work,LINE,{std::fma(std::fma(a,origin,b),origin,c),scale*std::fma(2*a,origin,b),a*scale*scale});
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  const CutFemLevelSetData::LineAnalysis& analysis = work.data.GetLineAnalysis();
  if(analysis.identicallyZero && computeDelta) throw std::domain_error("An identically zero line has no regular delta measure.");
  std::vector<long double> h(degree+1,0), d(degree+1,0);
  if(analysis.identicallyZero) for(unsigned p = 0; p <= degree; ++p) h[p] = 0.5L*FELevelSetMoments::SurfaceEndpointIntegration::IntegralPower(left,right,p);
  for(unsigned interval = 0; interval < analysis.positiveIntervalCount; ++interval) { const long double lower = origin+scale*analysis.positiveIntervals[interval][0], upper = origin+scale*analysis.positiveIntervals[interval][1]; for(unsigned p = 0; p <= degree; ++p) h[p] += FELevelSetMoments::SurfaceEndpointIntegration::IntegralPower(lower,upper,p); }
  if(computeDelta) for(unsigned n = 0; n < analysis.rootCount; ++n) { const CutFemLevelSetData::LineRoot& root = analysis.roots[n]; const long double coordinate = origin+scale*root.coordinate, weight = scale*(root.boundary != 0 && !fullBoundaryDelta ? 0.5L : 1.0L)*work.data.GetLineDeltaFactor(root); long double power = 1; for(unsigned p = 0; p <= degree; ++p) { d[p] += weight*power; power *= coordinate; } }
  positive.swap(h);
  delta.swap(d);
  Record("line_kernel",start,true);
}

struct StripWork {
  CompatibilityWork work;
  std::vector<std::vector<long double>> original;
  std::vector<std::vector<std::vector<long double>>> originalMetric;
  std::array<long double,4> bounds = {0,0,0,0};
  bool inverse = false;
  int component = -1;
  long double normalization = 1;
  std::vector<std::vector<long double>> normalized;
  FELevelSetMoments::Ex47Prepared::LineIntegrationStrip strip;
};

StripWork& PrepareStrip(const std::vector<std::vector<long double>>& coefficient, const FELevelSetMoments::LineIntegrationStrip& strip, const std::vector<std::vector<std::vector<long double>>>* gradient, bool inverse, int component) {
  static thread_local StripWork cache;
  const long double lower = std::min(strip.lowerConstant+strip.lowerSlope*strip.left,strip.lowerConstant+strip.lowerSlope*strip.right), upper = std::max(strip.upperConstant+strip.upperSlope*strip.left,strip.upperConstant+strip.upperSlope*strip.right);
  const std::array<long double,4> bounds = {strip.left,strip.right-strip.left,lower,upper-lower};
  if(!(bounds[1] > 0) || !(bounds[3] > 0)) throw std::invalid_argument("Invalid strip bounds.");
  const std::vector<std::vector<std::vector<long double>>> empty;
  if(!cache.work.data.IsInitialized() || cache.original != coefficient || cache.bounds != bounds || cache.originalMetric != (gradient == nullptr ? empty : *gradient) || cache.inverse != inverse || cache.component != component) {
    cache.original = coefficient;
    cache.bounds = bounds;
    cache.originalMetric = gradient == nullptr ? empty : *gradient;
    cache.inverse = inverse;
    cache.component = component;
    std::vector<std::vector<long double>> chart;
    FELevelSetMoments::ReferencePolynomialAffine2D(coefficient,bounds[0],bounds[1],bounds[2],bounds[3],chart);
    cache.work.data.ResetPowerCoefficients(QUAD,Flatten(chart));
    if(gradient != nullptr) { std::vector<std::vector<std::vector<long double>>> transformed(gradient->size()); for(unsigned c = 0; c < gradient->size(); ++c) FELevelSetMoments::ReferencePolynomialAffine2D((*gradient)[c],bounds[0],bounds[1],bounds[2],bounds[3],transformed[c]); FlattenMetric(transformed,cache.work.flatMetric,cache.work.metricDegree); }
    cache.normalization = cache.work.data.GetNormalizationScale();
    cache.normalized = FELevelSetMoments::Ex47Prepared::PlanarCoefficients(cache.work.data,cache.normalization);
    cache.work.workspace.plane.Bind(cache.work.data,cache.normalization,gradient != nullptr,inverse,component,gradient == nullptr ? nullptr : &cache.work.flatMetric,cache.work.metricDegree,&cache.work.workspace.statistics);
  }
  cache.strip.left = 0;
  cache.strip.right = 1;
  cache.strip.lowerConstant = (strip.lowerConstant+strip.lowerSlope*bounds[0]-bounds[2])/bounds[3];
  cache.strip.lowerSlope = strip.lowerSlope*bounds[1]/bounds[3];
  cache.strip.upperConstant = (strip.upperConstant+strip.upperSlope*bounds[0]-bounds[2])/bounds[3];
  cache.strip.upperSlope = strip.upperSlope*bounds[1]/bounds[3];
  cache.strip.fullBoundaryDelta = strip.fullBoundaryDelta;
  cache.strip.certifiedLowerZero = strip.certifiedLowerZero || FELevelSetMoments::LineIntegration::BoundaryIsZero(coefficient,strip.lowerConstant,strip.lowerSlope);
  cache.strip.certifiedUpperZero = strip.certifiedUpperZero || FELevelSetMoments::LineIntegration::BoundaryIsZero(coefficient,strip.upperConstant,strip.upperSlope);
  return cache;
}

std::array<long double,5> StripParent(const std::array<long double,4>& bounds, const std::array<long double,5>* parent, bool integrateOuter) { const std::array<long double,5> input = parent == nullptr ? std::array<long double,5>{0,1,0,1,1} : *parent; return {input[0]+input[1]*bounds[0],input[1]*bounds[1],input[2]+input[3]*bounds[2],input[3]*bounds[3],input[4]*bounds[3]*(integrateOuter ? bounds[1] : 1)}; }

void LineIntegration::FastGraphSample(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const FELevelSetMoments::LineIntegrationStrip& strip, long double x, bool computeDelta, bool lowerZero, bool upperZero, std::vector<long double>& values, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>& active, long double mappedDerivative, const std::array<long double,5>* coordinates, int derivativeSign, bool inverseSurface, int normalComponent) {
  StripWork& cache = PrepareStrip(coefficient,strip,surfaceGradient,inverseSurface,normalComponent); const std::array<long double,5> parent = StripParent(cache.bounds,coordinates,false); const long double point = (x-cache.bounds[0])/cache.bounds[1]; const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now(); FELevelSetMoments::Ex47Prepared::LineIntegration::FastGraphSample(cache.normalized,degree,cache.strip,point,computeDelta,lowerZero,upperZero,values,nullptr,active,mappedDerivative < 0 ? mappedDerivative : mappedDerivative*cache.bounds[3]/cache.normalization,&parent,derivativeSign,inverseSurface,normalComponent,&cache.work.workspace.plane,1); const unsigned count = values.size()/2; for(unsigned n = count; n < values.size(); ++n) values[n] /= cache.normalization; Record("graph_sample",start,true);
}

void LineIntegration::IntegrateStrip(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const FELevelSetMoments::LineIntegrationStrip& strip, bool graph, bool mapped, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::array<long double,5>* coordinates, const std::vector<unsigned char>* requested, const FELevelSetMoments::OuterQuadratureRules& rules, bool inverseSurface, int normalComponent) {
  StripWork& cache = PrepareStrip(coefficient,strip,surfaceGradient,inverseSurface,normalComponent); const std::array<long double,5> parent = StripParent(cache.bounds,coordinates,true); std::vector<long double> result(values.size(),0), resultErrors(errors.size(),0); const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now(); FELevelSetMoments::Ex47Prepared::LineIntegration::IntegrateStrip(cache.normalized,degree,cache.strip,graph,mapped,computeDelta,absoluteH,absoluteD*cache.normalization,relative,maximumDepth,maximumEvaluations,evaluations,accepted,result,resultErrors,nullptr,&parent,requested,rules,inverseSurface,normalComponent,&cache.work.workspace.plane,1); const unsigned count = result.size()/2; for(unsigned n = 0; n < result.size(); ++n) { const long double factor = n >= count ? 1/cache.normalization : 1; values[n] += factor*result[n]; errors[n] += factor*resultErrors[n]; } Record("strip_kernel",start,true);
}

void TriangleIntegration::Integrate(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned& selected, const std::vector<std::vector<std::vector<long double>>>* gradient, const std::array<long double,5>* coordinates, const std::vector<unsigned char>* requested, unsigned direction, bool fullBoundaryDelta, const std::vector<unsigned>* zeroBoundary, const FELevelSetMoments::OuterQuadratureRules& rules, bool inverseSurface, int normalComponent) {
  PlanarEndpointIntegration::Integrate(coefficient,degree,true,computeDelta,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,selected,gradient,coordinates,requested,direction,fullBoundaryDelta,zeroBoundary,rules,inverseSurface,normalComponent);
}

void PlanarEndpointIntegration::Integrate(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool triangle, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned& selected, const std::vector<std::vector<std::vector<long double>>>* gradient, const std::array<long double,5>* coordinates, const std::vector<unsigned char>* requested, unsigned direction, bool fullBoundaryDelta, const std::vector<unsigned>* zeroBoundary, const FELevelSetMoments::OuterQuadratureRules& rules, bool inverseSurface, int normalComponent) {
  CompatibilityWork& work = Work(triangle ? TRI : QUAD);
  std::vector<std::vector<long double>> chart = coefficient;
  std::array<long double,5> parent = coordinates == nullptr ? std::array<long double,5>{0,1,0,1,1} : *coordinates;
  std::vector<std::vector<std::vector<long double>>> chartGradient;

  ResetPolynomial(work,triangle ? TRI : QUAD,Flatten(chart));
  std::vector<std::vector<long double>> metric;
  unsigned metricDegree = 2;
  if(gradient != nullptr) FlattenMetric(*gradient,metric,metricDegree);
  PlanarData(work.data,degree,computeDelta,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,selected,gradient == nullptr ? nullptr : &metric,metricDegree,&parent,requested,direction,fullBoundaryDelta,rules,inverseSurface,normalComponent,work);

}

unsigned PlanarEndpointIntegration::Compare(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool triangle, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, const std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned selected, const std::vector<std::vector<std::vector<long double>>>* gradient, const std::array<long double,5>* coordinates, bool fullBoundaryDelta, const std::vector<unsigned>* zeroBoundary, const FELevelSetMoments::OuterQuadratureRules& rules, bool inverseSurface, int normalComponent) {
    bool depends[2] = {false,false}, constantMetric = true;
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) if(coefficient[i][j] != 0) { depends[0] = depends[0] || i != 0; depends[1] = depends[1] || j != 0; }
    if(gradient != nullptr) for(const std::vector<std::vector<long double>>& component : *gradient) for(unsigned i = 0; i < component.size(); ++i) for(unsigned j = 0; j < component[i].size(); ++j) if(i+j != 0 && component[i][j] != 0) constantMetric = false;
    if(!triangle && depends[0] != depends[1] && constantMetric) return ::FELevelSetMoments::PlanarEndpointIntegration::Compare(coefficient,degree,triangle,computeDelta,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,selected,gradient,coordinates,fullBoundaryDelta,zeroBoundary,rules,inverseSurface,normalComponent);

  std::vector<long double> second, secondErrors; unsigned independent = 0; try { Integrate(coefficient,degree,triangle,computeDelta,0.1L*absoluteH,0.1L*absoluteD,0.1L*relative,maximumDepth,maximumEvaluations,second,secondErrors,evaluations,accepted,independent,gradient,coordinates,nullptr,selected == 0 ? 0 : 3-selected,fullBoundaryDelta,zeroBoundary,rules,inverseSurface,normalComponent); } catch(const std::runtime_error& error) { FELevelSetMoments::WarnQuadrature("Independent planar validation",std::string("validation unavailable: ")+error.what(),0,std::numeric_limits<long double>::infinity(),absoluteD); return 0; } const unsigned count = values.size()/2; for(unsigned n = 0; n < values.size(); ++n) { errors[n] += secondErrors[n]+std::abs(values[n]-second[n]); if(errors[n] > (n < count ? absoluteH : absoluteD)+relative*std::abs(values[n])) FELevelSetMoments::WarnQuadrature("Independent planar validation","endpoint projections disagree",n,errors[n],(n < count ? absoluteH : absoluteD)+relative*std::abs(values[n])); } return independent != selected;
}

bool IntegrateCanonicalGraph(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned order, bool triangle, bool computeDelta, long double derivativeBound, long double absoluteH, long double absoluteDelta, long double relativeTolerance, unsigned maximumDepth, std::vector<long double>& total, std::vector<long double>& totalError, unsigned& accepted, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient) {
  const unsigned count = (degree+1)*(degree+2)/2; std::vector<long double> values(2*count,0), errors(2*count,0); FELevelSetMoments::LineIntegrationStrip strip = {0,1,0,0,1,triangle ? -1.0L : 0.0L}; std::vector<std::vector<long double>> field = coefficient; std::vector<std::vector<std::vector<long double>>> transposed; if(order == 2) { for(unsigned i = 0; i < 3; ++i) for(unsigned j = i+1; j < 3; ++j) std::swap(field[i][j],field[j][i]); if(surfaceGradient != nullptr) { FELevelSetMoments::SurfaceGradientTranspose(*surfaceGradient,transposed); surfaceGradient = &transposed; } } unsigned evaluations = 0; LineIntegration::IntegrateStrip(field,degree,strip,true,false,computeDelta,absoluteH,absoluteDelta,relativeTolerance,maximumDepth,200000,evaluations,accepted,values,errors,surfaceGradient); if(order == 2) for(unsigned i = 0; i <= degree; ++i) for(unsigned j = i+1; j <= degree-i; ++j) for(unsigned pair = 0; pair < 2; ++pair) { std::swap(values[pair*count+FELevelSetMoments::ReferenceMomentIndex2D(degree,i,j)],values[pair*count+FELevelSetMoments::ReferenceMomentIndex2D(degree,j,i)]); std::swap(errors[pair*count+FELevelSetMoments::ReferenceMomentIndex2D(degree,i,j)],errors[pair*count+FELevelSetMoments::ReferenceMomentIndex2D(degree,j,i)]); } total.swap(values); totalError.swap(errors); return true;
}

void IntegrateCanonicalGraphInterval(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned order, long double left, long double right, bool computeDelta, long double absoluteH, long double absoluteDelta, long double relativeTolerance, unsigned depth, unsigned maximumDepth, std::vector<long double>& total, std::vector<long double>& totalError, unsigned& accepted, bool triangle, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>* active) {
  const unsigned count = (degree+1)*(degree+2)/2; FELevelSetMoments::LineIntegrationStrip strip = {left,right,0,0,1,triangle ? -1.0L : 0.0L}; std::vector<std::vector<long double>> field = coefficient; std::vector<std::vector<std::vector<long double>>> transposed; if(order == 2) { for(unsigned i = 0; i < 3; ++i) for(unsigned j = i+1; j < 3; ++j) std::swap(field[i][j],field[j][i]); if(surfaceGradient != nullptr) { FELevelSetMoments::SurfaceGradientTranspose(*surfaceGradient,transposed); surfaceGradient = &transposed; } } unsigned evaluations = 0; std::vector<long double> values(2*count,0), errors(2*count,0); LineIntegration::IntegrateStrip(field,degree,strip,true,false,computeDelta,absoluteH,absoluteDelta,relativeTolerance,maximumDepth,200000,evaluations,accepted,values,errors,surfaceGradient,nullptr,active); if(order == 2) for(unsigned i = 0; i <= degree; ++i) for(unsigned j = i+1; j <= degree-i; ++j) for(unsigned pair = 0; pair < 2; ++pair) { std::swap(values[pair*count+FELevelSetMoments::ReferenceMomentIndex2D(degree,i,j)],values[pair*count+FELevelSetMoments::ReferenceMomentIndex2D(degree,j,i)]); std::swap(errors[pair*count+FELevelSetMoments::ReferenceMomentIndex2D(degree,i,j)],errors[pair*count+FELevelSetMoments::ReferenceMomentIndex2D(degree,j,i)]); } for(unsigned n = 0; n < total.size(); ++n) { total[n] += values[n]; totalError[n] += errors[n]; }
}

void SurfaceEndpointIntegration::Integrate(FELevelSetMoments::SurfaceEndpointData& data, bool tetrahedron, std::vector<long double>& values, std::vector<long double>& errors) {
  const GeomElType geometry = tetrahedron ? TET : (data.wedge ? WEDGE : HEX);
  CompatibilityWork& work = Work(geometry);
  const unsigned axes = !tetrahedron && !data.wedge && !data.unitCube ? 7U : 0U;
  const long double measure = axes == 7 ? 8 : 1;
  ResetPolynomial(work,geometry,AffinePolynomial(data.coefficient,3,2,axes,-1,2));
  std::vector<std::vector<long double>> metric = data.gradient;
  unsigned metricDegree = metric.empty() ? 2U : FELevelSetMoments::SurfacePolynomialCubeSize(metric[0].size())-1;
  for(std::vector<long double>& component : metric) component = AffinePolynomial(component,3,metricDegree,axes,-1,2);
  SpatialData(work.data,data.degree,data.computeDelta,data.absoluteH/measure,data.absoluteD/measure,data.relative,data.maximumDepth,data.maximumEvaluations,values,errors,metric.empty() ? nullptr : &metric,metricDegree,data.preferredInner,data.preferredOuter,data.computePositive,data.quadratureRules,data.inverseSurface,data.normalComponent,work);
  if(axes != 0) { std::vector<long double> native, nativeErrors; ReferenceOutput(3,data.degree,values,errors,native,nativeErrors); values.swap(native); errors.swap(nativeErrors); }
  data.powers = work.workspace.surface.powers;
  data.evaluations += work.workspace.surface.evaluations;
  data.accepted += work.workspace.surface.accepted;
  data.innerCounts = work.workspace.surface.innerCounts;
  data.outerCounts = work.workspace.surface.outerCounts;
  data.graphAxis = work.workspace.surface.graphAxis;
}

void SurfaceEndpointIntegration::Adaptive(FELevelSetMoments::SurfaceEndpointData& data, const FELevelSetMoments::SurfaceEndpointDomain& cell, unsigned outer, long double left, long double right, int mapping, long double endpoint, long double delta, unsigned exponent, long double budget, unsigned depth, const std::vector<unsigned char>& requested, std::vector<long double>& total, std::vector<long double>& totalErrors) {
  const GeomElType geometry = cell.tetrahedron ? TET : (cell.wedge ? WEDGE : HEX);
  const unsigned axes = !cell.tetrahedron && !cell.wedge && !data.unitCube ? 7U : 0U;
  const long double measure = axes == 7 ? 8 : 1;
  const FELevelSetMoments::OuterInterval legacyInterval = {mapping == 0 ? left : (mapping > 0 ? endpoint : endpoint-delta),mapping == 0 ? right : (mapping > 0 ? endpoint+delta : endpoint),endpoint,delta,mapping,exponent,true};
  if(mapping != 0 && !FELevelSetMoments::OuterIntervalQuadrature::Resolved(FELevelSetMoments::GaussKronrodRules::Get(3),legacyInterval,left,right)) throw FELevelSetMoments::OuterQuadratureResolutionError("The requested coordinate chart has no resolved endpoint rule.");
  CompatibilityWork& work = Work(geometry);
  ResetPolynomial(work,geometry,AffinePolynomial(data.coefficient,3,2,axes,-1,2));
  std::vector<std::vector<long double>> metric = data.gradient;
  const unsigned metricDegree = metric.empty() ? 2U : FELevelSetMoments::SurfacePolynomialCubeSize(metric[0].size())-1;
  for(std::vector<long double>& component : metric) component = AffinePolynomial(component,3,metricDegree,axes,-1,2);
  FELevelSetMoments::Ex47Prepared::SurfaceEndpointData& prepared = work.workspace.surface;
  prepared.prepared = &work.data;
  prepared.normalization = work.data.GetNormalizationScale();
  prepared.coefficient = work.data.GetNormalizedPowerCoefficients();
  prepared.referenceSurface = !metric.empty();
  prepared.physicalMetric = metric.empty() ? nullptr : &metric;
  prepared.metricDegree = metricDegree;
  prepared.gradient = metric;
  prepared.statistics = &work.workspace.statistics;
  prepared.degree = data.degree;
  prepared.maximumDepth = data.maximumDepth;
  prepared.maximumEvaluations = data.maximumEvaluations;
  prepared.evaluations = data.evaluations;
  prepared.accepted = data.accepted;
  prepared.graphAxis = data.graphAxis;
  prepared.preferredInner = data.preferredInner;
  prepared.preferredOuter = data.preferredOuter;
  prepared.wedge = cell.wedge;
  prepared.computePositive = data.computePositive;
  prepared.computeDelta = data.computeDelta;
  prepared.inverseSurface = data.inverseSurface;
  prepared.normalComponent = data.normalComponent;
  prepared.absoluteH = data.absoluteH/measure;
  prepared.absoluteD = data.absoluteD*prepared.normalization/measure;
  prepared.relative = data.relative;
  prepared.quadratureRules = data.quadratureRules;
  prepared.powers = data.powers;
  prepared.domainVolume = data.domainVolume/measure;
  FELevelSetMoments::Ex47Prepared::SurfaceEndpointDomain domain;
  domain.tetrahedron = cell.tetrahedron;
  domain.wedge = cell.wedge;
  domain.vertices = cell.vertices;
  for(std::array<long double,3>& vertex : domain.vertices) for(unsigned axis = 0; axis < 3; ++axis) if((axes&(1U<<axis)) != 0) vertex[axis] = (vertex[axis]+1)/2;
  std::vector<unsigned char> active = requested;
  const unsigned count = data.powers.size();
  if(axes != 0) for(unsigned n = 0; n < count; ++n) for(unsigned pair = 0; pair < 2; ++pair) if(requested[pair*count+n]) { const std::array<unsigned,3>& power = data.powers[n]; for(unsigned i = (axes&1U) != 0 ? 0U : power[0]; i <= power[0]; ++i) for(unsigned j = (axes&2U) != 0 ? 0U : power[1]; j <= power[1]; ++j) for(unsigned k = (axes&4U) != 0 ? 0U : power[2]; k <= power[2]; ++k) active[pair*count+FELevelSetMoments::HexMomentIndex(i,j,k)] = 1; }
  std::vector<long double> values(total.size(),0), errors(totalErrors.size(),0);
  const bool transformedOuter = (axes&(1U<<outer)) != 0;
  FELevelSetMoments::Ex47Prepared::SurfaceEndpointIntegration::Adaptive(prepared,domain,outer,mapping == 0 && transformedOuter ? (left+1)/2 : left,mapping == 0 && transformedOuter ? (right+1)/2 : right,mapping,transformedOuter ? (endpoint+1)/2 : endpoint,transformedOuter ? delta/2 : delta,exponent,budget,depth,active,values,errors);
  for(unsigned n = count; n < values.size(); ++n) { values[n] /= prepared.normalization; errors[n] /= prepared.normalization; }
  if(axes != 0) { std::vector<long double> native, nativeErrors; ReferenceOutput(3,data.degree,values,errors,native,nativeErrors); values.swap(native); errors.swap(nativeErrors); }
  for(unsigned n = 0; n < requested.size(); ++n) if(requested[n]) { total[n] += values[n]; totalErrors[n] += errors[n]; }
  data.evaluations = prepared.evaluations;
  data.accepted = prepared.accepted;
}

}
