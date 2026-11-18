#include "UnitElementMoments.hpp"
#include "Tri7DirectionSelection.hpp"
#include "Tet15LevelSetMoments.hpp"
#include "Basis.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#if __has_include("OuterIntervalQuadrature.hpp")
#define FEMUS_UNIFIED_OUTER 1
#endif

long double OuterBenchmarkField(const std::string& field, const std::vector<long double>& point) {
  const unsigned dimension = point.size();
  if(field == "graph" || field == "parabola") { long double value = point.back()-0.1L-0.65L*point[0]*point[0]; if(dimension == 3) value += 0.45L*point[0]*point[0]-0.1L*point[1]*point[1]; return value; }
  if(field == "tiny_plane") return point[0]-(1.0L-1.0L/1048576.0L);
  if(field == "full") return 1;
  if(field == "empty") return -1;
  if(field == "plane") return point.back()-0.375L;
  if(field == "boundary_low") return point.back();
  if(field == "boundary_high") return point.back()-1;
  if(field == "oblique") { long double value = -0.43L; for(unsigned d = 0; d < dimension; ++d) value += point[d]; return value; }
  if(field == "slab") return (point[0]-0.21L)*(point[0]-0.68L);
  if(field == "near_plane") { long double value = point.back()-0.25L; for(unsigned d = 0; d+1U < dimension; ++d) value += 1e-8L*point[d]*point[d]; return value; }
  if(field == "bubble") { long double value = 1; long double remaining = 1; for(unsigned d = 0; d < dimension; ++d) { value *= point[d]; remaining -= point[d]; } return value*remaining-(dimension == 2 ? 1.0L/64 : 1.0L/512); }
  const long double radius = field == "tiny_ball" ? 1.0L/128 : (dimension == 2 ? 0.15L : 0.10L);
  const long double centers[3] = {field == "touch_ball" ? radius : (field == "near_tangent" ? radius+1e-12L : 0.29L),0.23L,0.19L};
  const long double axes[3] = {0.18L,0.07L,0.05L};
  const bool ellipse = field == "ellipse" || field == "rotated_ellipse";
  std::vector<long double> shifted(dimension);
  for(unsigned d = 0; d < dimension; ++d) shifted[d] = point[d]-centers[d];
  if(field == "rotated_ellipse") { const long double first = shifted[0]; shifted[0] = 0.8L*first+0.6L*shifted[1]; shifted[1] = -0.6L*first+0.8L*shifted[1]; }
  long double value = ellipse ? 1 : radius*radius;
  for(unsigned d = 0; d < dimension; ++d) value -= ellipse ? shifted[d]*shifted[d]/(axes[d]*axes[d]) : shifted[d]*shifted[d];
  return value;
}

void OuterBenchmarkGeometry(const std::vector<long double>& point, std::vector<long double>& physical) {
  physical = point;
  physical[0] += 0.1L*point[0]*point[0];
  physical[1] += 0.12L*point[0]*point[0];
  if(point.size() == 3) physical[2] += 0.08L*point[0]*point[1];
}

void OuterBenchmarkString(const std::string& value) {
  std::cout << '"';
  for(unsigned i = 0; i < value.size(); ++i) { if(value[i] == '\\' || value[i] == '"') std::cout << '\\'; if(value[i] == '\n') std::cout << "\\n"; else std::cout << value[i]; }
  std::cout << '"';
}

void OuterBenchmarkMoments(unsigned element, unsigned degree, const std::vector<double>& nodal, bool surface, const std::vector<std::vector<double>>* geometry, double tolerance, unsigned validation, unsigned direction, std::vector<long double>& values) {
  const unsigned dimension = element < 2 ? 2U : 3U;
  if(element%2U == 0) { std::vector<long double> errors; FELevelSetMoments::IntegrateUnitTensorNodalMoments(nodal,dimension,degree,values,errors,tolerance,tolerance,surface,geometry,direction,validation); return; }
  if(element == 1) {
    std::vector<std::vector<double>> h;
    std::vector<std::vector<double>> d;
    FELevelSetMoments::IntegrateTri7NodalSelected(nodal,degree,h,d,tolerance,tolerance,true,false,validation,direction,nullptr,surface,nullptr,geometry);
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    values.assign(2U*count,0);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) { const unsigned n = FELevelSetMoments::UnitMomentIndex(2,i,j); values[n] = h[i][j]; values[count+n] = d[i][j]; }
  }
  else {
    std::vector<std::vector<std::vector<double>>> h;
    std::vector<std::vector<std::vector<double>>> d;
    FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,degree,h,d,tolerance,tolerance,true,24,direction,200000,nullptr,validation,nullptr,surface,nullptr,geometry);
    const unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
    values.assign(2U*count,0);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) { const unsigned n = FELevelSetMoments::UnitMomentIndex(3,i,j,k); values[n] = h[i][j][k]; values[count+n] = d[i][j][k]; }
  }
}

int main(int argc, char** argv) {
  if(argc != 3) return 64;
  const unsigned element = std::strtoul(argv[1],nullptr,10);
  const unsigned degree = std::strtoul(argv[2],nullptr,10);
  if(element > 3 || degree > 10) return 64;
  const unsigned dimension = element < 2 ? 2U : 3U;
  const bool simplex = element%2U != 0;
  femus::QuadBiquadratic quad;
  femus::TriBiquadratic tri;
  femus::HexBiquadratic hex;
  femus::TetBiquadratic tet;
  femus::basis* basis = element == 0 ? static_cast<femus::basis*>(&quad) : (element == 1 ? static_cast<femus::basis*>(&tri) : (element == 2 ? static_cast<femus::basis*>(&hex) : static_cast<femus::basis*>(&tet)));
  const unsigned nodes[4] = {9,7,27,15};
  std::cout << std::setprecision(21);
  std::string field;
  unsigned metric;
  unsigned direction;
  unsigned repetitions;
  double tolerance;
  while(std::cin >> field >> tolerance >> metric >> direction >> repetitions) {
    std::vector<double> nodal(nodes[element]);
    std::vector<std::vector<double>> coordinates(dimension,std::vector<double>(nodes[element]));
    for(unsigned node = 0; node < nodes[element]; ++node) {
      std::vector<long double> point(dimension);
      for(unsigned axis = 0; axis < dimension; ++axis) point[axis] = simplex ? basis->GetXcoarse(node)[axis] : (basis->GetXcoarse(node)[axis]+1)/2;
      nodal[node] = static_cast<double>(OuterBenchmarkField(field,point));
      std::vector<long double> physical;
      OuterBenchmarkGeometry(point,physical);
      for(unsigned axis = 0; axis < dimension; ++axis) coordinates[axis][node] = static_cast<double>(physical[axis]);
    }
    std::vector<double> seconds;
    std::vector<long double> values;
    std::string failure;
    for(unsigned repeat = 0; repeat < repetitions; ++repeat) {
      FELevelSetMoments::LineIntegration::Statistics() = FELevelSetMoments::LineIntegrationStatistics();
      FELevelSetMoments::EndpointStripQuadrature::Statistics() = FELevelSetMoments::EndpointStripStatistics();
      FELevelSetMoments::SurfaceEndpointIntegration::Statistics() = FELevelSetMoments::SurfaceEndpointStatistics();
#ifdef FEMUS_UNIFIED_OUTER
      for(unsigned dim = 2; dim <= 3; ++dim) FELevelSetMoments::OuterIntervalQuadrature::Statistics(dim) = FELevelSetMoments::OuterIntervalStatistics();
#endif
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      try { OuterBenchmarkMoments(element,degree,nodal,metric != 0,metric == 2 ? &coordinates : nullptr,tolerance,1,direction,values); }
      catch(const std::exception& error) { failure = error.what(); }
      seconds.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
      if(!failure.empty()) break;
    }
    std::sort(seconds.begin(),seconds.end());
    const FELevelSetMoments::LineIntegrationStatistics& line = FELevelSetMoments::LineIntegration::Statistics();
    const FELevelSetMoments::EndpointStripStatistics& curve = FELevelSetMoments::EndpointStripQuadrature::Statistics();
    const FELevelSetMoments::SurfaceEndpointStatistics& surface = FELevelSetMoments::SurfaceEndpointIntegration::Statistics();
    unsigned long long rootSamples = line.fastSamples+curve.stripSamples+curve.bulkSamples;
#ifndef FEMUS_UNIFIED_OUTER
    rootSamples += line.fullSamples;
#endif
    std::cout << "{\"counter_schema\":2,\"field\":"; OuterBenchmarkString(field);
    std::cout << ",\"element\":" << element << ",\"degree\":" << degree << ",\"tolerance\":" << tolerance << ",\"metric\":" << metric << ",\"direction\":" << direction << ",\"seconds\":" << seconds[seconds.size()/2U] << ",\"root_samples\":" << rootSamples << ",\"slices\":" << surface.quadCalls+surface.triCalls << ",\"outer_samples\":" << surface.outerSamples << ",\"statistics\":[";
#ifdef FEMUS_UNIFIED_OUTER
    for(unsigned dim = 2; dim <= 3; ++dim) {
      if(dim != 2) std::cout << ',';
      const FELevelSetMoments::OuterIntervalStatistics& stat = FELevelSetMoments::OuterIntervalQuadrature::Statistics(dim);
      std::cout << '[' << stat.regularPanels << ',' << stat.transformedPanels << ',' << stat.emergency7 << ',' << stat.emergency3 << ',' << stat.subdivisions << ',' << stat.coordinateCollapses << ',' << stat.resolutionLimits << ',' << stat.samples << ',' << stat.precisionSamples << ',' << stat.acceptedMoments << ',' << stat.skippedMoments << ',' << stat.regularIntervals << ',' << stat.oneSingularIntervals << ',' << stat.twoSingularIntervals << ',' << stat.higherOrderMaps << ']';
    }
#elif defined(FEMUS_BASELINE_OUTER_COUNTERS)
    std::cout << '[' << line.benchmarkRegular15+curve.benchmarkRegular15 << ',' << line.benchmarkMapped15+curve.benchmarkMapped15 << ",0,0," << line.benchmarkSubdivisions+curve.stripBisections+curve.bulkBisections << ",0,0," << rootSamples << ",0,null,null,null,null,null,null],[" << surface.benchmarkRegular15 << ',' << surface.benchmarkMapped15 << ',' << surface.benchmarkEmergency7 << ',' << surface.benchmarkEmergency3 << ',' << surface.benchmarkSubdivisions << ',' << surface.benchmarkCollapses << ",0," << surface.outerSamples << ",0,null,null,null,null,null,null]";
#else
    std::cout << "[0," << line.regularGraphPanels+line.mappedPanels+curve.rule15 << ",0,0," << curve.stripBisections+curve.bulkBisections << ",0,0," << rootSamples << ",0,null,null,null,null,null,null]";
#endif
    std::cout << "],\"legacy\": [";
#ifdef FEMUS_UNIFIED_OUTER
    std::cout << "0,0,0";
#else
    std::cout << curve.rule31 << ',' << line.pairPanels << ',';
#ifdef FEMUS_BASELINE_OUTER_COUNTERS
    std::cout << surface.benchmarkRegular31+surface.benchmarkMapped31;
#else
    std::cout << "null";
#endif
#endif
    std::cout << "],\"values\":[";
    if(failure.empty()) for(unsigned n = 0; n < values.size(); ++n) { if(n != 0) std::cout << ','; std::cout << values[n]; }
    std::cout << "],\"error\":";
    if(failure.empty()) std::cout << "null"; else OuterBenchmarkString(failure);
    std::cout << "}\n" << std::flush;
  }
}
