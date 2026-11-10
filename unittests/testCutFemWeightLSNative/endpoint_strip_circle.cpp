#include "UnitElementMoments.hpp"
#include "Tri7DirectionSelection.hpp"
#include "EndpointStripQuadrature.hpp"
#include "Basis.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

volatile long double endpointCircleSink = 0;

void EndpointCircleIntegrate(const std::vector<double>& nodal, bool triangle, unsigned degree, const std::string& method, long double width, std::vector<long double>& values) {
  if(method != "mapped" && method != "baseline") throw std::invalid_argument("Invalid endpoint quadrature test method.");
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  if(method == "baseline") {
    if(!triangle) { std::vector<long double> errors; FELevelSetMoments::IntegrateUnitTensorNodalMoments(nodal,2,degree,values,errors,1e-10L,1e-10L,true); return; }
    std::vector<std::vector<double>> positive;
    std::vector<std::vector<double>> surface;
    FELevelSetMoments::IntegrateTri7NodalSelected(nodal,degree,positive,surface,1e-10,1e-10,true,false,1,0,nullptr,true);
    values.assign(2U*count,0);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) { const unsigned n = FELevelSetMoments::UnitMomentIndex(2,i,j); values[n] = positive[i][j]; values[count+n] = surface[i][j]; }
    return;
  }
  std::vector<std::vector<long double>> coefficient(3,std::vector<long double>(3,0));
  if(triangle) { std::vector<std::vector<double>> field; FELevelSetMoments::Tri7ReferenceCoefficientsFromNodalValues(nodal,field); for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) coefficient[i][j] = field[i][j]; }
  else {
    std::vector<long double> tensor;
    FELevelSetMoments::UnitTensorPolynomialFromNodalValues(nodal,2,tensor);
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) coefficient[i][j] = tensor[FELevelSetMoments::Hex27Index(i,j,0)];
  }
  std::vector<std::vector<std::vector<long double>>> gradient(2,std::vector<std::vector<long double>>(3,std::vector<long double>(3,0)));
  for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) { if(i > 0) gradient[0][i-1U][j] = i*coefficient[i][j]; if(j > 0) gradient[1][i][j-1U] = j*coefficient[i][j]; }
  FELevelSetMoments::EndpointStripQuadrature::StripWidth() = width;
  FELevelSetMoments::EndpointStripQuadrature::PreferredDirection() = 2;
  std::vector<long double> native;
  std::vector<long double> errors;
  unsigned evaluations = 0;
  unsigned accepted = 0;
  unsigned selected = 0;
  const bool complete = FELevelSetMoments::EndpointStripQuadrature::Try(coefficient,degree,triangle,true,1e-12L,1e-12L,1e-12L,32,200000,native,errors,evaluations,accepted,selected,&gradient);
  if(!complete) throw std::runtime_error(FELevelSetMoments::EndpointStripQuadrature::Statistics().failure);
  values.assign(2U*count,0);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) { const unsigned n = FELevelSetMoments::UnitMomentIndex(2,i,j); const unsigned source = LineSliceUtilities::MomentIndex(degree,i,j); values[n] = native[source]; values[count+n] = native[count+source]; }
}

int main(int argc, char** argv) {
  const std::string method = argc > 1 ? argv[1] : "mapped";
  const long double width = argc > 2 ? std::strtold(argv[2],nullptr) : 0.04393398282201787134L;
  const unsigned degree = argc > 3 ? std::strtoul(argv[3],nullptr,10) : 10U;
  const bool triangle = argc > 4 && std::strtoul(argv[4],nullptr,10) != 0;
  const unsigned repetitions = argc > 5 ? std::strtoul(argv[5],nullptr,10) : 16U;
  if(degree > 10 || repetitions == 0 || width <= 0 || (method != "mapped" && method != "baseline" && method != "switched") || (method == "switched" && triangle)) return 64;
  femus::QuadBiquadratic quad;
  femus::TriBiquadratic tri;
  std::vector<double> nodal(triangle ? 7U : 9U);
  for(unsigned n = 0; n < nodal.size(); ++n) {
    const double* point = triangle ? tri.GetXcoarse(n) : quad.GetXcoarse(n);
    const long double x = triangle ? point[0] : (point[0]+1)/2;
    const long double y = triangle ? point[1] : (point[1]+1)/2;
    nodal[n] = 0.15L*0.15L-(x-0.29L)*(x-0.29L)-(y-0.23L)*(y-0.23L);
  }
  FELevelSetMoments::EndpointStripQuadrature::Statistics() = FELevelSetMoments::EndpointStripStatistics();
  FELevelSetMoments::LineIntegration::Statistics() = FELevelSetMoments::LineIntegrationStatistics();
  std::vector<long double> values;
  try {
    const std::chrono::steady_clock::time_point firstStart = std::chrono::steady_clock::now();
    EndpointCircleIntegrate(nodal,triangle,degree,method,width,values);
    const double firstSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-firstStart).count();
    const FELevelSetMoments::EndpointStripStatistics statistics = FELevelSetMoments::EndpointStripQuadrature::Statistics();
    const FELevelSetMoments::LineIntegrationStatistics baselineStatistics = FELevelSetMoments::LineIntegration::Statistics();
    std::vector<double> shifted = nodal;
    for(unsigned n = 0; n < shifted.size(); ++n) shifted[n] += 0.00000000002910383045673370361328125;
    std::vector<double> samples;
    std::vector<long double> scratch;
    for(unsigned batch = 0; batch < 11; ++batch) {
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      for(unsigned n = 0; n < repetitions; ++n) { EndpointCircleIntegrate(n%2U == 0 ? shifted : nodal,triangle,degree,method,width,scratch); endpointCircleSink += scratch[0]*1e-300L; }
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/repetitions;
      if(batch > 1) samples.push_back(seconds);
    }
    std::sort(samples.begin(),samples.end());
    std::cout << std::setprecision(21);
    std::cout << "T," << method << ',' << triangle << ',' << degree << ',' << width << ',' << samples[4] << ',' << firstSeconds << ',' << statistics.stripSamples << ',' << statistics.bulkSamples << ',' << statistics.stripBisections << ',' << statistics.bulkBisections << ',' << statistics.maximumStripDepth << ',' << statistics.maximumBulkDepth << ',' << statistics.minimumBulkWidth << ',' << statistics.rule15 << ',' << statistics.rule31 << ',' << statistics.geometryBuilds << ',' << statistics.rejected << ',' << baselineStatistics.fastSamples+baselineStatistics.fullSamples << '\n';
    const unsigned count = values.size()/2U;
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) { const unsigned n = FELevelSetMoments::UnitMomentIndex(2,i,j); std::cout << "M," << i << ',' << j << ',' << values[n] << ',' << values[count+n] << '\n'; }
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
  return 0;
}
