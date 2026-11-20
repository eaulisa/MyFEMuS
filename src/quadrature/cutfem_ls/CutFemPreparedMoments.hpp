#ifndef FEMUS_PREPARED_MOMENTS_HPP
#define FEMUS_PREPARED_MOMENTS_HPP

#include "CutFemLevelSetData.hpp"
#include "UnitElementMoments.hpp"
#include "Hex27DirectionSelection.hpp"
#include "EndpointSimplexChart.hpp"
#include "SurfaceEndpointIntegration.hpp"

namespace FELevelSetMoments {

inline void IntegrateLine(const CutFemLevelSetData& data, unsigned degree, bool referenceSurface, bool inverseSurface, std::vector<long double>& values, std::vector<long double>& errors, bool computePositive = true) {
  const CutFemLevelSetData::LineAnalysis& analysis = data.GetLineAnalysis();
  if(analysis.identicallyZero) throw std::domain_error("The delta of an identically zero LINE field is undefined.");
  const unsigned count = degree+1;
  values.assign(2*count,0);
  errors.assign(2*count,0);
  if(computePositive) for(unsigned interval = 0; interval < analysis.positiveIntervalCount; ++interval) {
    const long double left = analysis.positiveIntervals[interval][0], right = analysis.positiveIntervals[interval][1];
    long double rightPower = 1, geometric = 0;
    for(unsigned power = 0; power <= degree; ++power) { geometric = rightPower+left*geometric; values[power] += (right-left)*geometric/(power+1); rightPower *= right; }
  }
  const std::array<std::array<long double,3>,3> identity = {{{1,0,0},{0,1,0},{0,0,1}}};
  for(unsigned n = 0; n < analysis.rootCount; ++n) {
    const CutFemLevelSetData::LineRoot& root = analysis.roots[n];
    const long double factor = root.boundary == 0 ? 1.0L : 0.5L;
    const long double weight = factor*(inverseSurface ? data.GetLineTangentialSurfaceFactorFromCofactor(root,identity) : (referenceSurface ? data.GetLineReferenceSurfaceFactor(root) : data.GetLineDeltaFactor(root)));
    const long double unit = root.coordinate;
    long double power = 1;
    for(unsigned k = 0; k <= degree; ++k) { values[count+k] += weight*power; power *= unit; }
  }
  for(unsigned n = 0; n < values.size(); ++n) errors[n] = 128*std::numeric_limits<long double>::epsilon()*std::abs(values[n]);
}

struct PreparedMomentWorkspace {
  PreparedPlane plane;
  SurfaceEndpointData surface;
  PreparedStatistics statistics;
  CutFemLevelSetData chart;
  std::vector<double> chartNodes;
  std::vector<double> metricNodes;
  std::vector<std::vector<double>> metricGeometry;
  std::vector<std::vector<long double>> metric;
  std::vector<std::vector<long double>> chartMetric;
  std::vector<std::vector<long double>> referenceMetric;
  std::vector<std::vector<std::vector<long double>>> policy;
  unsigned metricDegree = 0;
  unsigned chartMetricDegree = 0;
  unsigned long long metricRevision = 0;
  unsigned long long chartMetricRevision = 0;
};

inline unsigned PolynomialIndex(unsigned dimension, unsigned degree, unsigned i, unsigned j = 0, unsigned k = 0) { return dimension == 1 ? i : (dimension == 2 ? (degree+1)*i+j : ((degree+1)*i+j)*(degree+1)+k); }

inline void TrimMetric(unsigned dimension, std::vector<std::vector<long double>>& metric, unsigned& degree) {
  unsigned maximum = 0;
  for(const std::vector<long double>& component : metric) for(unsigned n = 0; n < component.size(); ++n) if(component[n] != 0) { unsigned index = n; for(unsigned axis = 0; axis < dimension; ++axis) { maximum = std::max(maximum,index%(degree+1)); index /= degree+1; } }
  if(maximum == degree) return;
  unsigned size = 1;
  for(unsigned axis = 0; axis < dimension; ++axis) size *= maximum+1;
  for(std::vector<long double>& component : metric) {
    std::vector<long double> trimmed(size);
    for(unsigned i = 0; i <= maximum; ++i) for(unsigned j = 0; j <= (dimension > 1 ? maximum : 0U); ++j) for(unsigned k = 0; k <= (dimension > 2 ? maximum : 0U); ++k) trimmed[PolynomialIndex(dimension,maximum,i,j,k)] = component[PolynomialIndex(dimension,degree,i,j,k)];
    component.swap(trimmed);
  }
  degree = maximum;
}

inline void PrepareMetric(const CutFemLevelSetData& data, const std::vector<std::vector<double>>& geometry, PreparedMomentWorkspace& workspace) {
  if(!data.GetNodalValues().empty() && workspace.metricNodes == data.GetNodalValues() && workspace.metricGeometry == geometry && !workspace.metric.empty()) return;
  const unsigned dimension = data.GetDimension();
  std::vector<std::vector<long double>> coordinates(dimension), gradient(dimension);
  for(unsigned component = 0; component < dimension; ++component) {
    std::vector<long double> centered(geometry[component].size());
    for(unsigned node = 0; node < centered.size(); ++node) centered[node] = static_cast<long double>(geometry[component][node])-static_cast<long double>(geometry[component][0]);
    CutFemLevelSetData::ExpandNodalPolynomial(data.GetGeometry(),centered,coordinates[component]);
    gradient[component] = data.GetDerivativeCoefficients(component);
  }
  if(dimension == 2) {
    std::vector<std::vector<std::vector<long double>>> g(2,std::vector<std::vector<long double>>(3,std::vector<long double>(3))), x = g;
    for(unsigned component = 0; component < 2; ++component) for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) { g[component][i][j] = gradient[component][3*i+j]; x[component][i][j] = coordinates[component][3*i+j]; }
    SurfaceGradientCofactor2D(g,x,1,1);
    workspace.metricDegree = g[0].size()-1;
    workspace.metric.assign(2,std::vector<long double>((workspace.metricDegree+1)*(workspace.metricDegree+1)));
    for(unsigned component = 0; component < 2; ++component) for(unsigned i = 0; i <= workspace.metricDegree; ++i) for(unsigned j = 0; j <= workspace.metricDegree; ++j) workspace.metric[component][(workspace.metricDegree+1)*i+j] = g[component][i][j];
  }
  else {
    SurfaceGradientCofactor3D(gradient,coordinates,1);
    workspace.metric = gradient;
    workspace.metricDegree = SurfacePolynomialCubeSize(gradient[0].size())-1;
  }
  TrimMetric(dimension,workspace.metric,workspace.metricDegree);
  workspace.metricNodes = data.GetNodalValues();
  workspace.metricGeometry = geometry;
  ++workspace.metricRevision;
  ++workspace.statistics.geometryPreparations;
}

inline void PlanarPolicy(const CutFemLevelSetData& data, const std::vector<std::vector<long double>>* metric, unsigned degree, std::vector<std::vector<std::vector<long double>>>& output) {
  const unsigned components = metric == nullptr ? data.GetRootDimension() : metric->size();
  if(metric == nullptr) degree = 2;
  output.resize(components);
  for(unsigned component = 0; component < components; ++component) {
    const std::vector<long double>& source = metric == nullptr ? data.GetRootDerivativeCoefficients(component) : (*metric)[component];
    output[component].resize(degree+1);
    for(unsigned i = 0; i <= degree; ++i) { output[component][i].resize(degree+1); for(unsigned j = 0; j <= degree; ++j) output[component][i][j] = source[(degree+1)*i+j]; }
  }
}

inline void IntegratePlanar(const CutFemLevelSetData& source, unsigned degree, double absolute, double relative, bool referenceSurface, const std::vector<std::vector<double>>* geometry, unsigned direction, unsigned validation, const OuterQuadratureRules& rules, bool inverseSurface, int normalComponent, bool computePositive, PreparedMomentWorkspace& workspace, std::vector<long double>& values, std::vector<long double>& errors) {
  const CutFemLevelSetData* field = &source;
  const bool barycentric = source.GetGeometry() == TRI && direction == 3;
  if(barycentric) { if(!workspace.chart.IsInitialized() || source.GetNodalValues().empty() || workspace.chartNodes != source.GetNodalValues()) { source.PermuteSimplexVertices({1,0,2,3},workspace.chart); workspace.chartNodes = source.GetNodalValues(); workspace.chartMetricRevision = 0; } field = &workspace.chart; }
  const long double normalization = source.GetNormalizationScale();
  const std::vector<std::vector<long double>>* metric = nullptr;
  unsigned metricDegree = 0;
  if(geometry != nullptr) {
    PrepareMetric(source,*geometry,workspace);
    metric = &workspace.metric;
    metricDegree = workspace.metricDegree;
    if(barycentric) { if(workspace.chartMetricRevision != workspace.metricRevision) { field->TransformParentComponents(metricDegree,*metric,workspace.chartMetricDegree,workspace.chartMetric); TrimMetric(field->GetDimension(),workspace.chartMetric,workspace.chartMetricDegree); workspace.chartMetricRevision = workspace.metricRevision; } metric = &workspace.chartMetric; metricDegree = workspace.chartMetricDegree; }
  }
  workspace.plane.Bind(*field,normalization,referenceSurface,inverseSurface,normalComponent,metric,metricDegree,&workspace.statistics);
  if(referenceSurface) PlanarPolicy(*field,metric,metricDegree,workspace.policy);
  const unsigned count = (degree+1)*(degree+2)/2;
  std::vector<unsigned char> requested(2*count,1);
  if(!computePositive) std::fill(requested.begin(),requested.begin()+count,0);
  const std::array<long double,5> parent = {0,1,0,1,1};
  const std::vector<unsigned> boundary = workspace.plane.Boundaries();
  const long double volume = 1, amplification = barycentric ? ReferencePower(3,degree) : 1;
  std::vector<long double> moments, momentErrors;
  unsigned evaluations = 0, accepted = 0, selected = 0;
  const std::vector<std::vector<long double>> coefficient = PlanarCoefficients(*field,normalization);
  PlanarEndpointIntegration::Integrate(coefficient,degree,source.GetGeometry() == TRI,true,absolute*volume/amplification,absolute*normalization*volume/amplification,relative/amplification,28,200000,moments,momentErrors,evaluations,accepted,selected,referenceSurface ? &workspace.policy : nullptr,&parent,&requested,barycentric ? 1U : direction,false,&boundary,rules,inverseSurface,normalComponent,&workspace.plane);
  if(validation == 2) PlanarEndpointIntegration::Compare(coefficient,degree,source.GetGeometry() == TRI,true,absolute*volume/amplification,absolute*normalization*volume/amplification,relative/amplification,28,200000,moments,momentErrors,evaluations,accepted,selected,referenceSurface ? &workspace.policy : nullptr,&parent,false,&boundary,rules,inverseSurface,normalComponent,&workspace.plane);
  if(barycentric) { std::vector<long double> original, originalErrors; EndpointSimplexChart::RestoreTriangle(degree,moments,momentErrors,original,originalErrors); moments.swap(original); momentErrors.swap(originalErrors); }
  std::vector<long double> h(count), d(count);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) { const unsigned original = ReferenceMomentIndex2D(degree,i,j), target = UnitMomentIndex(2,i,j); h[target] = moments[original]; d[target] = moments[count+original]/normalization; }
  values = h;
  values.insert(values.end(),d.begin(),d.end());
  std::vector<long double> eh(count), ed(count);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) { const unsigned original = ReferenceMomentIndex2D(degree,i,j), target = UnitMomentIndex(2,i,j); eh[target] = momentErrors[original]; ed[target] = momentErrors[count+original]/normalization; }
  errors = eh;
  errors.insert(errors.end(),ed.begin(),ed.end());
}

inline void IntegrateSpatial(const CutFemLevelSetData& source, unsigned degree, double absolute, double relative, bool referenceSurface, const std::vector<std::vector<double>>* geometry, unsigned direction, unsigned validation, const OuterQuadratureRules& rules, bool inverseSurface, int normalComponent, bool computePositive, PreparedMomentWorkspace& workspace, std::vector<long double>& values, std::vector<long double>& errors) {
  const CutFemLevelSetData* field = &source;
  const bool barycentric = source.GetGeometry() == TET && direction == 1;
  if(barycentric) { if(!workspace.chart.IsInitialized() || source.GetNodalValues().empty() || workspace.chartNodes != source.GetNodalValues()) { source.PermuteSimplexVertices({1,0,2,3},workspace.chart); workspace.chartNodes = source.GetNodalValues(); workspace.chartMetricRevision = 0; } field = &workspace.chart; }
  const long double normalization = field->GetNormalizationScale();
  const std::vector<std::vector<long double>>* metric = nullptr;
  unsigned metricDegree = 0;
  if(geometry != nullptr) {
    PrepareMetric(source,*geometry,workspace);
    metric = &workspace.metric;
    metricDegree = workspace.metricDegree;
    if(barycentric) { if(workspace.chartMetricRevision != workspace.metricRevision) { field->TransformParentComponents(metricDegree,*metric,workspace.chartMetricDegree,workspace.chartMetric); TrimMetric(field->GetDimension(),workspace.chartMetric,workspace.chartMetricDegree); workspace.chartMetricRevision = workspace.metricRevision; } metric = &workspace.chartMetric; metricDegree = workspace.chartMetricDegree; }
  }
  if(referenceSurface && geometry == nullptr && source.GetGeometry() == WEDGE) {
    workspace.referenceMetric.resize(3);
    for(unsigned component = 0; component < 3; ++component) workspace.referenceMetric[component] = field->GetRootDerivativeCoefficients(component);
    for(long double& value : workspace.referenceMetric[2]) value *= 0.5L;
    metric = &workspace.referenceMetric;
    metricDegree = 2;
  }
  SurfaceEndpointData& data = workspace.surface;
  data.prepared = field;
  data.physicalMetric = metric;
  data.metricDegree = metricDegree;
  data.referenceSurface = referenceSurface;
  data.normalization = normalization;
  data.statistics = &workspace.statistics;
  data.coefficient = field->GetNormalizedPowerCoefficients();
  data.gradient.clear();
  if(referenceSurface) { if(metric != nullptr) data.gradient = *metric; else { data.gradient.resize(3); for(unsigned component = 0; component < 3; ++component) data.gradient[component] = field->GetRootDerivativeCoefficients(component); } }
  data.degree = degree;
  data.evaluations = data.accepted = 0;
  data.graphAxis = data.preferredInner = data.preferredOuter = 3;
  data.innerCounts.fill(0);
  data.outerCounts.fill(0);
  data.retainUncertified = false;
  data.wedge = source.GetGeometry() == WEDGE;
  data.computePositive = computePositive;
  data.computeDelta = true;
  data.inverseSurface = inverseSurface;
  data.normalComponent = normalComponent;
  data.quadratureRules = rules;
  data.zeroFaces.resize(field->GetBoundaryCount());
  for(unsigned side = 0; side < data.zeroFaces.size(); ++side) data.zeroFaces[side] = field->IsZeroOnBoundary(side);
  if(direction != 0) { if(source.GetGeometry() == TET) data.preferredOuter = direction == 1 ? 0U : direction-2; else { unsigned middle; Hex27SelectedAxes(direction,data.preferredInner,middle,data.preferredOuter); } }
  const long double measure = 1, amplification = barycentric ? ReferencePower(4,degree) : 1;
  const long double positiveScale = normalComponent >= 0 ? source.GetUnitToReferenceMeasureScale() : 1;
  data.absoluteH = absolute*measure/(amplification*positiveScale*positiveScale);
  data.absoluteD = absolute*normalization*measure/amplification;
  data.relative = relative/amplification;
  std::vector<long double> moments, momentErrors;
  SurfaceEndpointIntegration::Integrate(data,source.GetGeometry() == TET,moments,momentErrors);
  if(validation == 2) SurfaceEndpointIntegration::Compare(data,source.GetGeometry() == TET,moments,momentErrors);
  if(barycentric) { std::vector<long double> original, originalErrors; EndpointSimplexChart::RestoreTetrahedron(degree,moments,momentErrors,original,originalErrors); moments.swap(original); momentErrors.swap(originalErrors); }
  const unsigned count = moments.size()/2;
  std::vector<long double> h(moments.begin(),moments.begin()+count), d(moments.begin()+count,moments.end());
  for(long double& value : d) value /= normalization;
  values = h;
  values.insert(values.end(),d.begin(),d.end());
  std::vector<long double> eh(momentErrors.begin(),momentErrors.begin()+count), ed(momentErrors.begin()+count,momentErrors.end());
  for(long double& value : ed) value /= normalization;
  errors = eh;
  errors.insert(errors.end(),ed.begin(),ed.end());
}

}


#endif
