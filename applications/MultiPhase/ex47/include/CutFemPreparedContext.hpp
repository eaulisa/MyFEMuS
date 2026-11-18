#ifndef FEMUS_EX47_PREPARED_CONTEXT_HPP
#define FEMUS_EX47_PREPARED_CONTEXT_HPP

#include "CutFemLevelSetData.hpp"

namespace FELevelSetMoments {
namespace Ex47Prepared {

struct PreparedStatistics {
  unsigned long long lineRestrictions = 0;
  unsigned long long surfaceRestrictions = 0;
  unsigned long long metricRestrictions = 0;
  unsigned long long rootAnalyses = 0;
  unsigned long long endpointSamples = 0;
  unsigned long long geometryPreparations = 0;
};

struct PreparedPlane {
  const CutFemLevelSetData* data = nullptr;
  const std::vector<std::vector<long double>>* metric = nullptr;
  long double normalization = 1;
  unsigned metricDegree = 0;
  bool referenceSurface = false;
  bool inverseSurface = false;
  int normalComponent = -1;
  PreparedStatistics* statistics = nullptr;
  mutable std::array<CutFemLevelSetData,2> lines;
  mutable std::array<std::vector<std::vector<long double>>,2> lineMetric;
  mutable std::array<bool,2> lineReady = {false,false};
  mutable std::array<long double,2> lineCoordinate = {0,0};

  void Bind(const CutFemLevelSetData& field, long double scale, bool surface, bool inverse, int component, const std::vector<std::vector<long double>>* physicalMetric, unsigned physicalDegree, PreparedStatistics* counts) {
    data = &field;
    metric = physicalMetric;
    normalization = scale;
    referenceSurface = surface;
    inverseSurface = inverse;
    normalComponent = component;
    metricDegree = physicalDegree;
    statistics = counts;
    lineReady.fill(false);
  }

  static unsigned InnerAxis(unsigned order) { if(order < 1 || order > 2) throw std::invalid_argument("Invalid planar direction."); return order == 1 ? 1U : 0U; }

  const CutFemLevelSetData& Fiber(unsigned order, long double coordinate) const {
    const unsigned index = order-1;
    if(!lineReady[index] || lineCoordinate[index] != coordinate) {
      const unsigned fixed = 1-InnerAxis(order);
      data->Restrict(fixed,coordinate,lines[index],referenceSurface && metric == nullptr);
      if(statistics != nullptr) ++statistics->lineRestrictions;
      if(metric != nullptr) { data->RestrictComponents(fixed,coordinate,metricDegree,*metric,lineMetric[index]); if(statistics != nullptr) ++statistics->metricRestrictions; }
      lineReady[index] = true;
      lineCoordinate[index] = coordinate;
    }
    return lines[index];
  }

  long double RootCoordinate(unsigned order, const CutFemLevelSetData::LineRoot& root) const { const CutFemLevelSetData& line = lines[order-1]; return line.GetParentOrigin()[InnerAxis(order)]+line.GetParentMap()[InnerAxis(order)][0]*root.coordinate; }

  long double WorkingDerivative(unsigned order, const CutFemLevelSetData::LineRoot& root) const { return root.derivative/(normalization*lines[order-1].GetParentMap()[InnerAxis(order)][0]); }

  long double Density(unsigned order, long double x, long double y, long double signedDerivative) const {
    if(!referenceSurface) { CutFemLevelSetData::LineRoot root; root.coordinate = y; root.derivative = signedDerivative; root.multiplicity = 1; return CutFemLevelSetData::GetRegularRootDeltaFactor(root); }
    const CutFemLevelSetData& line = Fiber(order,x);
    const long double origin = line.GetParentOrigin()[InnerAxis(order)], scale = line.GetParentMap()[InnerAxis(order)][0];
    CutFemLevelSetData::LineRoot root;
    root.coordinate = scale == 0 ? 0 : (y-origin)/scale;
    root.derivative = signedDerivative*normalization*(scale == 0 ? 1 : scale);
    root.multiplicity = 1;
    long double density;
    if(!referenceSurface) density = line.GetLineDeltaFactor(root);
    else if(metric == nullptr) density = line.GetLineReferenceSurfaceFactor(root);
    else {
      std::array<long double,3> h = {0,0,0};
      const std::vector<std::vector<long double>>& polynomial = lineMetric[order-1];
      for(unsigned component = 0; component < polynomial.size(); ++component) { long double value = polynomial[component][metricDegree]; for(unsigned power = metricDegree; power > 0; --power) value = value*root.coordinate+polynomial[component][power-1]; h[component] = value; }
      if(normalComponent >= 0) density = line.GetLineNormalSurfaceFactor(root,h,static_cast<unsigned>(normalComponent));
      else if(inverseSurface) density = line.GetLineTangentialSurfaceFactor(root,h);
      else density = line.GetLinePhysicalSurfaceFactor(root,h);
    }
    return normalization*(scale == 0 ? 1 : scale)*density;
  }

  long double Discriminant(unsigned order, long double x) const { const long double factor = data->GetNormalizationScale()/normalization; return factor*factor*data->EvaluateFiberDiscriminant(InnerAxis(order),x); }

  long double DiscriminantDerivative(unsigned order, long double x) const { const std::vector<long double>& coefficient = data->GetFiberProjection(InnerAxis(order)).discriminant; long double value = 0; for(unsigned power = coefficient.size()-1; power > 0; --power) value = value*x+power*coefficient[power]; const long double factor = data->GetNormalizationScale()/normalization; return factor*factor*value; }

  std::vector<long double> BoundaryPolynomial(unsigned order, long double constant, long double slope) const { std::vector<long double> result; data->GetFiberBoundaryPolynomial(InnerAxis(order),constant,slope,result); const long double factor = data->GetNormalizationScale()/normalization; for(long double& value : result) value *= factor; return result; }

  std::vector<long double> DiscriminantPolynomial(unsigned order) const { std::vector<long double> result = data->GetFiberProjection(InnerAxis(order)).discriminant; const long double factor = data->GetNormalizationScale()/normalization; for(long double& value : result) value *= factor*factor; return result; }

  void ShiftedDiscriminant(unsigned order, long double endpoint, std::array<long double,5>& output) const { data->GetShiftedFiberDiscriminant(InnerAxis(order),endpoint,output); const long double factor = data->GetNormalizationScale()/normalization; for(long double& value : output) value *= factor*factor; }

  bool ZeroBoundary(unsigned order, bool upper) const { const unsigned inner = InnerAxis(order); return data->IsZeroOnBoundary(data->GetGeometry() == TRI ? (upper ? 0U : inner+1) : 2*inner+static_cast<unsigned>(upper)); }

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
}

#endif
