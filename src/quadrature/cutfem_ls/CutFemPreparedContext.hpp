#ifndef FEMUS_PREPARED_CONTEXT_HPP
#define FEMUS_PREPARED_CONTEXT_HPP

#include "CutFemLevelSetData.hpp"

namespace FELevelSetMoments {

struct PreparedStatistics {
  unsigned long long lineRestrictions = 0;
  unsigned long long surfaceRestrictions = 0;
  unsigned long long metricRestrictions = 0;
  unsigned long long rootAnalyses = 0;
  unsigned long long endpointSamples = 0;
  unsigned long long geometryPreparations = 0;
};

struct PhysicalSurfaceMetric {
  std::array<long double,3> gradient = {};
  long double norm = 0;

  long double Factor(bool inverse, int component) const { return component >= 0 ? gradient[component]/norm : (inverse ? 1/norm : norm); }
};

struct SharedGraphStrips;
struct EndpointStripGeometry;

struct PreparedPlane {
  SharedGraphStrips* graphs = nullptr;
  EndpointStripGeometry* endpoints = nullptr;
  unsigned surfaceBlocks = 0;
  const CutFemLevelSetData* data = nullptr;
  const std::vector<std::vector<long double>>* metric = nullptr;
  long double normalization = 1;
  unsigned metricDegree = 0;
  bool referenceSurface = false;
  bool inverseSurface = false;
  int normalComponent = -1;
  PreparedStatistics* statistics = nullptr;
  mutable std::vector<long double> point;
  mutable std::array<CutFemLevelSetData::LineAnalysis,2> localAnalysis;
  mutable std::array<long double,2> analysisCoordinate = {0,0};
  mutable std::array<bool,2> analysisReady = {false,false};

  void Bind(const CutFemLevelSetData& field, long double scale, bool surface, bool inverse, int component, const std::vector<std::vector<long double>>* physicalMetric, unsigned physicalDegree, PreparedStatistics* counts) {
    graphs = nullptr;
    endpoints = nullptr;
    surfaceBlocks = 0;
    data = &field;
    metric = physicalMetric;
    normalization = scale;
    referenceSurface = surface;
    inverseSurface = inverse;
    normalComponent = component;
    metricDegree = physicalDegree;
    statistics = counts;
    if(surface && physicalMetric == nullptr) point.resize(2);
    analysisReady.fill(false);
  }

  static unsigned InnerAxis(unsigned order) { if(order < 1 || order > 2) throw std::invalid_argument("Invalid planar direction."); return order == 1 ? 1U : 0U; }

  const CutFemLevelSetData::LineAnalysis& Analyze(unsigned order, long double coordinate) const {
    const unsigned index = order-1;
    if(!analysisReady[index] || analysisCoordinate[index] != coordinate) { data->AnalyzeLocalFiber(InnerAxis(order),coordinate,localAnalysis[index]); analysisCoordinate[index] = coordinate; analysisReady[index] = true; if(statistics != nullptr) ++statistics->rootAnalyses; }
    return localAnalysis[index];
  }

  long double RootCoordinate(unsigned order, const CutFemLevelSetData::LineRoot& root) const { return root.coordinate; }

  long double WorkingDerivative(unsigned order, const CutFemLevelSetData::LineRoot& root) const { return root.derivative/normalization; }

  static long double EvaluateMetricRow(const long double* row, unsigned degree, long double coordinate) {
    long double value = row[degree];
    for(unsigned power = degree; power > 0; --power) value = value*coordinate+row[power-1];
    return value;
  }

  PhysicalSurfaceMetric EvaluateMetric(unsigned order, long double x, long double y) const {
    const long double first = InnerAxis(order) == 0 ? y : x, second = InnerAxis(order) == 0 ? x : y;
    const unsigned size = metricDegree+1;
    PhysicalSurfaceMetric result;
    for(unsigned component = 0; component < metric->size(); ++component) {
      const long double* polynomial = (*metric)[component].data();
      long double value = EvaluateMetricRow(polynomial+metricDegree*size,metricDegree,second);
      for(unsigned row = metricDegree; row > 0; --row) value = value*first+EvaluateMetricRow(polynomial+(row-1)*size,metricDegree,second);
      result.gradient[component] = value;
      result.norm = component == 0 ? std::abs(value) : std::hypot(result.norm,value);
    }
    if(!(result.norm > 0) || !std::isfinite(result.norm)) throw std::domain_error("A local fiber has no regular surface metric.");
    return result;
  }

  long double Density(unsigned order, long double x, long double y, long double signedDerivative) const {
    if(!referenceSurface) { CutFemLevelSetData::LineRoot root; root.coordinate = y; root.derivative = signedDerivative; root.multiplicity = 1; return CutFemLevelSetData::GetRegularRootDeltaFactor(root); }
    if(metric != nullptr) return EvaluateMetric(order,x,y).Factor(inverseSurface,normalComponent)/std::abs(signedDerivative);
    const unsigned inner = InnerAxis(order);
    point[inner] = y;
    point[1-inner] = x;
    std::array<long double,3> gradient;
    data->GetRootGradient(point,gradient);
    const std::vector<std::vector<long double>>& map = data->GetRootMap();
    unsigned selected = 0;
    for(unsigned component = 1; component < data->GetRootDimension(); ++component) if(std::abs(map[component][inner]) > std::abs(map[selected][inner])) selected = component;
    if(map[selected][inner] == 0) throw std::domain_error("A collapsed local fiber has no surface density.");
    long double remainder = signedDerivative*normalization;
    for(unsigned component = 0; component < data->GetRootDimension(); ++component) if(component != selected) remainder -= map[component][inner]*gradient[component];
    gradient[selected] = remainder/map[selected][inner];
    long double norm = 0;
    for(unsigned component = 0; component < data->GetRootDimension(); ++component) norm = std::hypot(norm,gradient[component]);
    return norm/std::abs(signedDerivative);
  }

  long double Discriminant(unsigned order, long double x) const { const long double factor = data->GetNormalizationScale()/normalization; return factor*factor*data->EvaluateLocalFiberDiscriminant(InnerAxis(order),x); }

  long double DiscriminantDerivative(unsigned order, long double x) const { std::array<long double,5> shifted; ShiftedDiscriminant(order,x,shifted); return shifted[1]; }

  std::vector<long double> BoundaryPolynomial(unsigned order, long double constant, long double slope) const { std::vector<long double> result; data->GetFiberBoundaryPolynomial(InnerAxis(order),constant,slope,result); const long double factor = data->GetNormalizationScale()/normalization; for(long double& value : result) value *= factor; return result; }

  std::vector<long double> DiscriminantPolynomial(unsigned order) const { std::vector<long double> result = data->GetFiberProjection(InnerAxis(order)).discriminant; const long double factor = data->GetNormalizationScale()/normalization; for(long double& value : result) value *= factor*factor; return result; }

  void ShiftedDiscriminant(unsigned order, long double endpoint, std::array<long double,5>& output) const { data->GetLocalShiftedFiberDiscriminant(InnerAxis(order),endpoint,output); const long double factor = data->GetNormalizationScale()/normalization; for(long double& value : output) value *= factor*factor; }

  bool ZeroBoundary(unsigned order, bool upper) const { const unsigned inner = InnerAxis(order); return data->IsZeroOnBoundary(data->GetGeometry() == TRI ? (upper ? 0U : inner+1) : 2*inner+static_cast<unsigned>(upper)); }

  bool ZeroBoundary(unsigned order, long double constant, long double slope) const {
    const bool triangle = data->GetGeometry() == TRI;
    if(slope == 0 && constant == 0.0L) return ZeroBoundary(order,false);
    if(constant == 1 && slope == (triangle ? -1.0L : 0.0L)) return ZeroBoundary(order,true);
    const std::vector<long double> coefficient = BoundaryPolynomial(order,constant,slope);
    for(long double value : coefficient) if(value != 0) return false;
    return true;
  }

  std::vector<unsigned> Boundaries() const {
    if(data->GetGeometry() == TRI) return {static_cast<unsigned>(data->IsZeroOnBoundary(1)),static_cast<unsigned>(data->IsZeroOnBoundary(2)),static_cast<unsigned>(data->IsZeroOnBoundary(0))};
    return {static_cast<unsigned>(data->IsZeroOnBoundary(0)),static_cast<unsigned>(data->IsZeroOnBoundary(1)),static_cast<unsigned>(data->IsZeroOnBoundary(2)),static_cast<unsigned>(data->IsZeroOnBoundary(3))};
  }
};

inline std::vector<std::vector<long double>> PlanarCoefficients(const CutFemLevelSetData& data, long double normalization) {
  const std::vector<long double>& power = data.GetPowerCoefficients();
  std::vector<std::vector<long double>> result(3,std::vector<long double>(3));
  for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) result[i][j] = power[3*i+j]/normalization;
  return result;
}

}


#endif
