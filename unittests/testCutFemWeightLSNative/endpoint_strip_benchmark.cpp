#include "UnitElementMoments.hpp"
#include "Tri7DirectionSelection.hpp"
#include "Basis.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

const char* const tangentNames[23] = {"axis_graph","tilted_graph","curved_graph","quadratic_graph","rational_graph","enriched_graph","near_plane_graph","swapped_graph","diagonal_graph","weak_graph","circle","clipped_circle","edge_circle","rotated_ellipse","two_sheets","hyperbola","bubble","perturbed_circle","tensor_lobes","hyperbolic_sheets","near_parallel_sheets","corner_circle","tiny_circle"};
volatile long double tangentSink = 0;

std::array<double,9> TangentCircle(double x, double y, double radius) { return {radius*radius-x*x-y*y,2*y,-1,2*x,0,0,-1,0,0}; }

std::array<double,9> TangentCoefficient(unsigned index) {
  std::array<double,9> c = {};
  if(index < 10) {
    c[0] = -0.3; c[1] = 1;
    if(index == 1) c[3] = 0.2;
    if(index == 2) c[6] = 0.2;
    if(index == 3) { c[2] = 0.15; c[3] = 0.1; }
    if(index == 4) { c[4] = 0.2; c[6] = 0.1; }
    if(index == 5) { c[4] = 0.1; c[5] = c[7] = -0.1; }
    if(index == 6) c[6] = 1e-8;
    if(index == 7) { c[1] = 0; c[3] = 1; c[2] = 0.1; }
    if(index == 8) c[3] = 1;
    if(index == 9) { c[0] = -0.1; c[1] = 0; c[2] = 1; }
    return c;
  }
  if(index == 10) return TangentCircle(0.29,0.23,0.15);
  if(index == 11) return TangentCircle(0.70,0.40,0.45);
  if(index == 12) return TangentCircle(0.13,0.74,0.22);
  if(index == 13) {
    const double cosine = std::cos(0.45), sine = std::sin(0.45), x = 0.4, y = 0.35;
    const double a = cosine*cosine/(0.28*0.28)+sine*sine/(0.12*0.12), b = 2*cosine*sine*(1/(0.28*0.28)-1/(0.12*0.12)), d = sine*sine/(0.28*0.28)+cosine*cosine/(0.12*0.12);
    return {1-a*x*x-b*x*y-d*y*y,2*d*y+b*x,-d,2*a*x+b*y,-b,0,-a,0,0};
  }
  if(index == 14) { c[0] = 0.21*0.68; c[3] = -0.89; c[6] = 1; }
  if(index == 15) { c[0] = 0.225; c[1] = -0.4; c[3] = -0.6; c[4] = 1; }
  if(index == 16) { c[0] = -1.0/64; c[4] = 1; c[5] = c[7] = -1; }
  if(index == 17) { c = TangentCircle(0.3,0.3,0.18); c[4] = 0.3; c[5] = c[7] = -0.3; }
  if(index == 18) { const double x[3] = {0.16,-1,1}, y[3] = {0.1875,-1,1}; for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) c[3*i+j] = x[i]*y[j]; c[0] -= 0.001; }
  if(index == 19) { c[0] = 0.1375; c[3] = -0.8; c[6] = 1; c[1] = 0.05; c[2] = -0.05; }
  if(index == 20) { c[0] = 0.21*0.68+0.25e-7; c[3] = -0.89; c[6] = 1; c[1] = -1e-7; c[2] = 1e-7; }
  if(index == 21) return TangentCircle(0.25,0.25,std::sqrt(0.125));
  if(index == 22) return TangentCircle(0.29,0.23,1.0/128);
  return c;
}

void TangentIntegrate(unsigned element, unsigned degree, const std::vector<double>& nodal, const std::vector<std::vector<double>>& geometry, const std::string& mode, double tolerance, std::vector<long double>& values) {
  if(element == 0) { std::vector<long double> errors; FELevelSetMoments::IntegrateUnitTensorNodalMoments(nodal,2,degree,values,errors,tolerance,tolerance,mode != "raw",mode == "physical" ? &geometry : nullptr); return; }
  std::vector<std::vector<double>> h;
  std::vector<std::vector<double>> d;
  FELevelSetMoments::IntegrateTri7NodalSelected(nodal,degree,h,d,tolerance,tolerance,true,false,1,0,nullptr,mode != "raw",nullptr,mode == "physical" ? &geometry : nullptr);
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  values.assign(2U*count,0);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) { const unsigned n = FELevelSetMoments::UnitMomentIndex(2,i,j); values[n] = h[i][j]; values[count+n] = d[i][j]; }
}

void TangentIntegrateQuadReference(unsigned degree, const std::vector<double>& nodal, const std::vector<std::vector<double>>& geometry, const std::string& mode, double tolerance, std::vector<long double>& values) {
  std::vector<std::vector<double>> h;
  std::vector<std::vector<double>> d;
  FELevelSetMoments::IntegrateQuad9NodalSelected(nodal,degree,h,d,tolerance,tolerance,true,false,1,0,nullptr,mode != "raw",nullptr,mode == "physical" ? &geometry : nullptr);
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  values.assign(2U*count,0);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) { const unsigned n = FELevelSetMoments::UnitMomentIndex(2,i,j); values[n] = h[i][j]; values[count+n] = d[i][j]; }
}

int main(int argc, char** argv) {
  const unsigned degree = argc > 1 ? std::strtoul(argv[1],nullptr,10) : 2U;
  const std::string mode = argc > 2 ? argv[2] : "surface";
  const unsigned rounds = argc > 3 ? std::strtoul(argv[3],nullptr,10) : 3U;
  const double tolerance = argc > 4 ? std::strtod(argv[4],nullptr) : 1e-10;
  const bool measure = argc > 5 && std::strtoul(argv[5],nullptr,10) != 0;
  const int onlyElement = argc > 6 ? std::strtol(argv[6],nullptr,10) : -1;
  const int onlyField = argc > 7 ? std::strtol(argv[7],nullptr,10) : -1;
  if(degree > 10 || rounds == 0 || (mode != "raw" && mode != "surface" && mode != "physical")) return 64;
  femus::QuadBiquadratic quad;
  femus::TriBiquadratic tri;
  std::cout << std::setprecision(21);
  unsigned failures = 0;
  for(unsigned round = 0; round < rounds; ++round) for(unsigned element = 0; element < 2; ++element) for(unsigned field = 0; field < 23; ++field) {
    if(element == 1 && field == 18) continue;
    if((onlyElement >= 0 && element != static_cast<unsigned>(onlyElement)) || (onlyField >= 0 && field != static_cast<unsigned>(onlyField))) continue;
    const std::array<double,9> coefficient = TangentCoefficient(field);
    const unsigned nodes = element == 0 ? 9U : 7U;
    std::vector<double> nodal(nodes);
    std::vector<std::vector<double>> geometry(2,std::vector<double>(nodes));
    for(unsigned node = 0; node < nodes; ++node) {
      const double* point = element == 0 ? quad.GetXcoarse(node) : tri.GetXcoarse(node);
      const long double x = element == 0 ? (point[0]+1)/2 : point[0], y = element == 0 ? (point[1]+1)/2 : point[1];
      long double value = 0;
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) value += coefficient[3*i+j]*FELevelSetMoments::ReferencePower(x,i)*FELevelSetMoments::ReferencePower(y,j);
      nodal[node] = value;
      geometry[0][node] = x+0.1L*x*x; geometry[1][node] = y+0.12L*x*x;
    }
    FELevelSetMoments::LineIntegration::Statistics() = FELevelSetMoments::LineIntegrationStatistics();
#ifdef FEMUS_ENDPOINT_STRIP_QUADRATURE_HPP
    FELevelSetMoments::EndpointStripQuadrature::Statistics() = FELevelSetMoments::EndpointStripStatistics();
#endif
    std::vector<long double> values;
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    std::string failure;
    try { TangentIntegrate(element,degree,nodal,geometry,mode,tolerance,values); }
    catch(const std::exception& error) { failure = error.what(); ++failures; }
    double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
#ifdef FEMUS_ENDPOINT_STRIP_QUADRATURE_HPP
    const FELevelSetMoments::EndpointStripStatistics statistics = FELevelSetMoments::EndpointStripQuadrature::Statistics();
#endif
    unsigned calls = 1;
    if(measure && failure.empty()) {
      std::vector<double> shifted = nodal;
      for(unsigned n = 0; n < shifted.size(); ++n) shifted[n] += 0.00000000002910383045673370361328125;
      unsigned sequence = 1;
      std::vector<long double> scratch;
      while(calls < 8192) {
        const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        for(unsigned n = 0; n < calls; ++n) { TangentIntegrate(element,degree,(sequence++%2U == 0 ? nodal : shifted),geometry,mode,tolerance,scratch); tangentSink += 1e-300L*scratch[0]; }
        if(std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count() >= 0.008) break;
        calls *= 2;
      }
      std::vector<double> samples;
      for(unsigned sample = 0; sample < 5; ++sample) {
        const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        for(unsigned n = 0; n < calls; ++n) { TangentIntegrate(element,degree,(sequence++%2U == 0 ? nodal : shifted),geometry,mode,tolerance,scratch); tangentSink += 1e-300L*scratch[0]; }
        samples.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()/calls);
      }
      std::sort(samples.begin(),samples.end());
      seconds = samples[2];
    }
    std::cout << "T," << element << ',' << degree << ',' << mode << ',' << field << ',' << tangentNames[field] << ',' << round << ',' << seconds << ',' << (failure.empty() ? "PASS" : "FAIL");
#ifdef FEMUS_ENDPOINT_STRIP_QUADRATURE_HPP
    std::cout << ',' << statistics.attempts << ',' << statistics.completed << ',' << statistics.rejected << ',' << statistics.stripSamples << ',' << statistics.stripBisections << ',' << statistics.bulkSamples << ',' << statistics.geometryBuilds << ',' << statistics.geometryHits << ',' << statistics.minimumDelta << ',' << statistics.maximumDelta;
#else
    std::cout << ",0,0,0,0,0,0,0,0,0,0";
#endif
#ifdef FEMUS_ENDPOINT_STRIP_QUADRATURE_HPP
    std::cout << ',' << statistics.bulkBisections << ',' << statistics.maximumStripDepth << ',' << statistics.maximumBulkDepth << ',' << statistics.rule15 << ',' << statistics.rule31;
#else
    std::cout << ",0,0,0,0,0";
#endif
    std::cout << ',' << calls << '\n';
    if(!failure.empty()) { std::cerr << element << ',' << field << ',' << failure << '\n'; continue; }
    if(round == 0) {
      std::cout << "C," << element << ',' << field;
      for(unsigned n = 0; n < coefficient.size(); ++n) std::cout << ',' << coefficient[n];
      std::cout << '\n';
      const unsigned count = values.size()/2U;
      for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) { const unsigned n = FELevelSetMoments::UnitMomentIndex(2,i,j); std::cout << "M," << element << ',' << degree << ',' << mode << ',' << field << ',' << i << ',' << j << ',' << values[n] << ',' << values[count+n] << '\n'; }
      if(element == 0) {
        try {
          std::vector<long double> referenceValues;
          TangentIntegrateQuadReference(degree,nodal,geometry,mode,tolerance,referenceValues);
          for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) { const unsigned n = FELevelSetMoments::UnitMomentIndex(2,i,j); std::cout << "F," << element << ',' << degree << ',' << mode << ',' << field << ',' << i << ',' << j << ',' << referenceValues[n] << ',' << referenceValues[count+n] << '\n'; }
        }
        catch(const std::exception& error) { ++failures; std::cerr << element << ',' << field << ',' << error.what() << '\n'; }
      }
    }
    std::cout << std::flush;
  }
  return failures == 0 ? 0 : 1;
}
