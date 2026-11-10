#include "LineIntegration.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

unsigned rootChecks = 0;
long double rootMaximumError = 0.0L;
volatile long double rootSink = 0.0L;

void RootCheck(long double value, long double reference) {
  const long double error = std::abs(value-reference)/(1.0L+std::abs(reference));
  rootMaximumError = std::max(rootMaximumError,error);
  if(!std::isfinite(value) || error > 1.0e-13L) throw std::runtime_error("A certified single-root moment differs from the general solver.");
  ++rootChecks;
}

void RootProbe(long double a, long double b, long double c, long double lower, long double upper, int sign, bool lowerZero = false, bool upperZero = false) {
  const unsigned degree = 10;
  const FELevelSetMoments::LineIntegrationStrip strip = {0,1,lower,0,upper,0};
  std::vector<std::vector<long double>> coefficient(3,std::vector<long double>(3,0.0L));
  coefficient[0][0] = c;
  coefficient[0][1] = b;
  coefficient[0][2] = a;
  const std::vector<unsigned char> active((degree+1U)*(degree+2U),1);
  const std::array<long double,5> coordinates = {-1,2,-0.75L,1.5L,3};
  std::vector<long double> previous;
  std::vector<long double> single;
  FELevelSetMoments::LineIntegration::FastGraphSample(coefficient,degree,strip,0.37L,true,lowerZero,upperZero,previous,nullptr,active,-1,&coordinates);
  FELevelSetMoments::LineIntegration::FastGraphSample(coefficient,degree,strip,0.37L,true,lowerZero,upperZero,single,nullptr,active,-1,&coordinates,sign);
  if(previous.size() != single.size()) throw std::runtime_error("A certified root changed the requested moment count.");
  for(unsigned i = 0; i < single.size(); ++i) RootCheck(single[i],previous[i]);
}

double RootMeasure(unsigned degree, bool single) {
  const FELevelSetMoments::LineIntegrationStrip strip = {0,1,0,0,1,0};
  const std::vector<std::vector<long double>> coefficient = {{-0.3L,1.2L,0.3L},{0.1L,0.02L,0.05L},{0.1L,0,0}};
  std::vector<std::vector<std::vector<long double>>> gradient;
  std::vector<std::vector<double>> input(3,std::vector<double>(3,0.0));
  for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) input[i][j] = coefficient[i][j];
  FELevelSetMoments::SurfaceGradientFromPolynomial(input,1,1,gradient);
  const std::vector<unsigned char> active((degree+1U)*(degree+2U),1);
  const std::array<long double,5> coordinates = {-1,2,-1,2,4};
  std::vector<long double> values;
  unsigned calls = 16;
  std::vector<double> timings;
  for(unsigned batch = 0; batch < 6; ++batch) {
    const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
    for(unsigned repeat = 0; repeat < calls; ++repeat) {
      const long double x = ((repeat%997U)+0.5L)/997.0L;
      FELevelSetMoments::LineIntegration::FastGraphSample(coefficient,degree,strip,x,true,false,false,values,&gradient,active,-1,&coordinates,single ? 1 : 0);
      for(unsigned index = 0; index < values.size(); ++index) rootSink += 1.0e-300L*values[index];
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    if(batch == 0 && seconds < 0.015 && calls < 131072) { calls *= 2; --batch; continue; }
    if(batch != 0) timings.push_back(seconds/calls);
  }
  std::sort(timings.begin(),timings.end());
  return timings[timings.size()/2];
}

void RootDispatchProbe() {
  const unsigned degree = 10;
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  const FELevelSetMoments::LineIntegrationStrip strip = {0,1,0,0,1,0};
  const std::vector<std::vector<long double>> coefficient = {{-0.4L,1.8L,1},{0,0,0},{-0.6L,-0.3L,0}};
  const unsigned long long before = FELevelSetMoments::LineIntegration::Statistics().singleRootSamples;
  std::vector<long double> previous(2U*count,0.0L);
  std::vector<long double> previousErrors(2U*count,0.0L);
  std::vector<long double> single(2U*count,0.0L);
  std::vector<long double> errors(2U*count,0.0L);
  unsigned evaluations = 0;
  unsigned accepted = 0;
  if(!FELevelSetMoments::LineIntegration::IntegrateFastGraph(coefficient,degree,1,false,true,1.5L,1.0e-13L,1.0e-13L,1.0e-13L,24,previous,previousErrors,accepted)) throw std::runtime_error("The independent graph reference did not complete.");
  FELevelSetMoments::LineIntegration::IntegrateStrip(coefficient,degree,strip,true,false,true,1.0e-13L,1.0e-13L,1.0e-13L,24,200000,evaluations,accepted,single,errors);
  for(unsigned index = 0; index < single.size(); ++index) RootCheck(single[index],previous[index]);
  const FELevelSetMoments::LineIntegrationStatistics statistics = FELevelSetMoments::LineIntegration::Statistics();
  if(statistics.singleRootSamples-before != statistics.fastSamples || statistics.pairPanels == 0) throw std::runtime_error("The cached root branch was not passed through both graph quadrature rules.");
}

int main(int argc, char** argv) {
  try {
    RootProbe(0,1,-0.3L,0,1,1);
    RootProbe(0,-1,0.3L,0,1,-1);
    RootProbe(1,2,-0.25L,0,1,1);
    RootProbe(-1,-2,0.25L,0,1,-1);
    RootProbe(1,-3,1,0,1,-1);
    RootProbe(-1,3,-1,0,1,1);
    RootProbe(1,0,-2,1,2,1);
    RootProbe(1,0,-2,-2,-1,-1);
    RootProbe(1.0e-18L,1,-0.4L,0,1,1);
    RootProbe(-1.0e-18L,-1,0.4L,0,1,-1);
    RootProbe(0.2L,1,0,0,1,1,true,false);
    RootProbe(1,1,-2,0,1,1,false,true);
    RootProbe(1,-1,0.16L,0,0.4L,-1);
    RootProbe(1,-1,0.16L,0.6L,1,1);
    RootProbe(1,2,1,0,1,1);
    RootProbe(1,2,2,0,1,1);
    RootProbe(1,-3,3,0,1,-1);
    RootProbe(1.0e-12L,2.0e-12L,-0.25e-12L,0,1,1);
    if(FELevelSetMoments::LineIntegration::Statistics().singleRootSamples != 18) throw std::runtime_error("The certified single-root path was not exercised.");
    RootDispatchProbe();
    std::cout << std::setprecision(17) << rootChecks << " certified single-root checks passed; maximum normalized difference " << rootMaximumError << '\n';
    if(argc > 1) for(unsigned index = 0; index < 3; ++index) {
      const unsigned degree = index == 0 ? 2U : (index == 1 ? 5U : 10U);
      const double previous = RootMeasure(degree,false);
      const double single = RootMeasure(degree,true);
      std::cout << "ROOT_TIMING," << degree << ',' << previous << ',' << single << ',' << previous/single << '\n';
    }
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
