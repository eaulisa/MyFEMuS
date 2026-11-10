#include "Quad9DirectionSelection.hpp"
#include "Tri7DirectionSelection.hpp"
#include "Hex27DirectionSelection.hpp"
#include "Tet15LevelSetMoments.hpp"
#include "Basis.hpp"
#ifdef FEMUS_UNIT_COORDINATE_MOMENTS
#include "UnitElementMoments.hpp"
#endif
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>

struct DimensionGraph {
  const char* name;
  double cubeOffset;
  double simplexOffset;
  double firstLinear;
  double firstQuadratic;
  double secondLinear;
  double secondQuadratic;
  double cross;
  double mixed;
  double dependentQuadratic;
};

const DimensionGraph dimensionGraphs[10] = {
  {"axis_plane",0.11,0.18,0,0,0,0,0,0,0},
  {"tilted_plane",-0.04,0.18,0.20,0,-0.12,0,0,0,0},
  {"near_plane",0.05,0.22,0.15,1.0e-6,0.10,2.0e-6,0,0,0},
  {"gentle_convex",-0.10,0.16,0.05,0.10,0.07,0.06,0,0,0},
  {"gentle_concave",0.12,0.24,-0.10,-0.15,0.04,-0.10,0,0,0},
  {"convex",-0.15,0.12,0.16,0.30,-0.08,0.23,0,0,0},
  {"saddle",0.07,0.20,0.12,0.25,0.09,-0.22,0,0,0},
  {"mixed_curvature",-0.08,0.18,-0.12,0.10,0.15,0.07,0.20,0,0},
  {"implicit_quadratic",0.03,0.21,0.16,0.11,-0.13,0.08,-0.06,0,0.15},
  {"variable_slope",-0.05,0.19,0.08,0.12,0.14,0.07,0.10,0.20,-0.05}
};

volatile double dimensionSink = 0.0;

double DimensionField(unsigned element, unsigned index, const double* point) {
  const DimensionGraph& graph = dimensionGraphs[index];
  const bool three = element >= 2;
  const bool simplex = element%2 != 0;
  const double x = point[0];
  const double y = three ? point[1] : 0.0;
  const double dependent = point[three ? 2 : 1];
  return dependent-(simplex ? graph.simplexOffset : graph.cubeOffset)-graph.firstLinear*x-graph.firstQuadratic*x*x-graph.secondLinear*y-graph.secondQuadratic*y*y-graph.cross*x*y+graph.mixed*x*dependent+graph.dependentQuadratic*dependent*dependent;
}

std::vector<double> DimensionNodal(unsigned element, unsigned index) {
  femus::QuadBiquadratic quad;
  femus::TriBiquadratic tri;
  femus::HexBiquadratic hex;
  femus::TetBiquadratic tet;
  const unsigned counts[4] = {9,7,27,15};
  std::vector<double> nodal(counts[element]);
  for(unsigned node = 0; node < nodal.size(); ++node) {
    const double* point = element == 0 ? quad.GetXcoarse(node) : (element == 1 ? tri.GetXcoarse(node) : (element == 2 ? hex.GetXcoarse(node) : tet.GetXcoarse(node)));
    nodal[node] = DimensionField(element,index,point);
  }
  return nodal;
}

double DimensionIntegrate(unsigned element, const std::vector<double>& nodal, unsigned degree) {
#ifdef FEMUS_UNIT_COORDINATE_MOMENTS
  if(element == 0 || element == 2) {
    std::vector<long double> values;
    std::vector<long double> errors;
    FELevelSetMoments::IntegrateUnitTensorNodalMoments(nodal,element == 0 ? 2U : 3U,degree,values,errors,1.0e-10L,1.0e-10L,true);
    const unsigned count = values.size()/2U;
    if(!(values[0] > 0) || !(values[count] > 0)) throw std::runtime_error("A unit graph benchmark returned an empty cut.");
    long double checksum = 0;
    for(unsigned n = 0; n < values.size(); ++n) checksum += (n+1U)*values[n];
    if(!std::isfinite(checksum)) throw std::runtime_error("A unit graph benchmark returned nonfinite moments.");
    return static_cast<double>(checksum);
  }
#endif
  std::vector<std::vector<double>> positive2;
  std::vector<std::vector<double>> surface2;
  std::vector<std::vector<std::vector<double>>> positive3;
  std::vector<std::vector<std::vector<double>>> surface3;
  if(element == 0) FELevelSetMoments::IntegrateQuad9NodalSelected(nodal,degree,positive2,surface2,1.0e-10,1.0e-10,true,false,1,0,nullptr,true);
  else if(element == 1) FELevelSetMoments::IntegrateTri7NodalSelected(nodal,degree,positive2,surface2,1.0e-10,1.0e-10,true,false,1,0,nullptr,true);
  else if(element == 2) FELevelSetMoments::IntegrateHex27NodalSelected(nodal,degree,positive3,surface3,1.0e-10,1.0e-10,true,1,0,nullptr,nullptr,true);
  else FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,degree,positive3,surface3,1.0e-10,1.0e-10,true,24,0,200000,nullptr,1,nullptr,true);
  double checksum = 0.0;
  if(element < 2) {
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) checksum += (1U+i+3U*j)*positive2[i][j]+(2U+2U*i+j)*surface2[i][j];
    if(!(positive2[0][0] > 0.0) || !(surface2[0][0] > 0.0)) throw std::runtime_error("A graph benchmark returned an empty cut.");
  }
  else {
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) checksum += (1U+i+3U*j+5U*k)*positive3[i][j][k]+(2U+2U*i+j+3U*k)*surface3[i][j][k];
    if(!(positive3[0][0][0] > 0.0) || !(surface3[0][0][0] > 0.0)) throw std::runtime_error("A graph benchmark returned an empty cut.");
  }
  if(!std::isfinite(checksum)) throw std::runtime_error("A graph benchmark returned nonfinite moments.");
  return checksum;
}

void DimensionMeasure(unsigned element, unsigned index, unsigned degree) {
  const char* names[4] = {"QUAD","TRI","HEX","TET"};
  const std::vector<double> nodal = DimensionNodal(element,index);
  FELevelSetMoments::SurfaceEndpointIntegration::Statistics() = FELevelSetMoments::SurfaceEndpointStatistics();
  dimensionSink += 1.0e-300*DimensionIntegrate(element,nodal,degree);
  const FELevelSetMoments::SurfaceEndpointStatistics counts = FELevelSetMoments::SurfaceEndpointIntegration::Statistics();
  unsigned calls = 1;
  while(calls < 8192) {
    const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
    for(unsigned repeat = 0; repeat < calls; ++repeat) dimensionSink += 1.0e-300*DimensionIntegrate(element,nodal,degree);
    if(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count() >= 0.008) break;
    calls *= 2;
  }
  std::vector<double> seconds;
  for(unsigned round = 0; round < 5; ++round) {
    const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
    for(unsigned repeat = 0; repeat < calls; ++repeat) dimensionSink += 1.0e-300*DimensionIntegrate(element,nodal,degree);
    seconds.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()/calls);
  }
  std::sort(seconds.begin(),seconds.end());
  std::cout << names[element] << ',' << degree << ',' << index << ',' << dimensionGraphs[index].name << ',' << seconds[2] << ',' << calls << ',' << counts.outerSamples << ',' << counts.quadCalls+counts.triCalls << ',' << counts.endpointStrips << ',' << counts.graphCells << '\n' << std::flush;
}

int main(int argc, char** argv) {
  const unsigned degree = argc > 1 ? std::strtoul(argv[1],nullptr,10) : 2U;
  if(degree > 10) return 64;
  std::cout << std::setprecision(17);
  std::cout << "element,degree,case,name,seconds,calls,slices,inner_2d_calls,endpoint_strips,graph_regions\n";
  try {
    for(unsigned index = 0; index < 10; ++index) for(unsigned element = 0; element < 4; ++element) DimensionMeasure(element,index,degree);
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
  return dimensionSink == 123456789.0 ? 1 : 0;
}
