#include "OuterIntervalQuadrature.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

struct IntervalContractContext {
  unsigned evaluations = 0;
  unsigned mode = 0;
  long double delta = 1;
  long double endpoint = 0;
  std::vector<unsigned> visits = std::vector<unsigned>(22,0);
};

void IntervalContractRequire(bool condition, const char* message) {
  if(!condition) throw std::runtime_error(message);
}

void IntervalContractSample(IntervalContractContext& context, const FELevelSetMoments::OuterIntervalPoint& point, const std::vector<unsigned char>& active, std::vector<long double>& values, std::vector<long double>& errors, long double absoluteH, long double absoluteD, long double relative) {
  ++context.evaluations;
  values.assign(active.size(),0);
  errors.assign(active.size(),0);
  if(context.mode >= 4) IntervalContractRequire(point.coordinate != context.endpoint,"A collapsed endpoint reached the sampler.");
  for(unsigned n = 0; n < active.size(); ++n) if(active[n]) {
    ++context.visits[n];
    const unsigned degree = n%11U;
    long double value = FELevelSetMoments::ReferencePower(point.coordinate,degree)*point.jacobian;
    if(context.mode == 1) value /= std::sqrt(point.offset);
    if(context.mode == 2) value /= std::sqrt(point.coordinate*(1-point.coordinate));
    if(context.mode == 3) value /= std::pow(point.offset,0.75L);
    if(context.mode == 6) value /= std::pow(point.offset,5.0L/6.0L);
    if(context.mode == 4) value = point.jacobian;
    if(context.mode == 7) value = point.jacobian/std::sqrt(point.coordinate-context.endpoint);
    if(context.mode == 5) value = FELevelSetMoments::ReferencePower(point.offset/context.delta,degree)*point.jacobian;
    values[n] = value;
  }
}

void IntervalContractRun(unsigned leftExponent, unsigned rightExponent, unsigned mode, unsigned dimension) {
  FELevelSetMoments::OuterIntervalQuadrature::Statistics(dimension) = FELevelSetMoments::OuterIntervalStatistics();
  const std::vector<FELevelSetMoments::OuterInterval> plan = FELevelSetMoments::OuterIntervalQuadrature::Plan(0,1,leftExponent,rightExponent,dimension);
  IntervalContractRequire(plan.size() == (leftExponent > 1 && rightExponent > 1 ? 2U : 1U),"Incorrect map count.");
  for(unsigned part = 0; part < plan.size(); ++part) IntervalContractRequire(plan[part].mapping == 0 || plan[part].delta == 1.0L/plan.size(),"Incorrect full-interval mapping width.");
  IntervalContractContext context;
  context.mode = mode;
  std::vector<unsigned char> active(22,1);
  std::vector<long double> values(22,0);
  std::vector<long double> errors(22,0);
  unsigned accepted = 0;
  for(unsigned part = 0; part < plan.size(); ++part) FELevelSetMoments::OuterIntervalQuadrature::Adaptive(context,IntervalContractSample,plan[part],0,1,1e-12L/plan.size(),1e-12L/plan.size(),1e-12L,0,40,200000,context.evaluations,accepted,dimension,active,values,errors);
  long double maximum = 0;
  for(unsigned n = 0; n < values.size(); ++n) {
    const unsigned degree = n%11U;
    long double exact = 1.0L/(degree+1U);
    if(mode == 1) exact = 1.0L/(degree+0.5L);
    if(mode == 2) exact = std::tgamma(degree+0.5L)*std::tgamma(0.5L)/std::tgamma(degree+1.0L);
    if(mode == 3) exact = 1.0L/(degree+0.25L);
    if(mode == 6) exact = 1.0L/(degree+1.0L/6.0L);
    maximum = std::max(maximum,std::abs(values[n]-exact));
    IntervalContractRequire(std::abs(values[n]-exact) <= 1e-11L,"An analytic interval moment failed.");
  }
  const FELevelSetMoments::OuterIntervalStatistics& statistics = FELevelSetMoments::OuterIntervalQuadrature::Statistics(dimension);
  if(leftExponent == 1 && rightExponent == 1) IntervalContractRequire(statistics.transformedPanels == 0 && statistics.regularPanels > 0,"A regular interval was transformed.");
  else IntervalContractRequire(statistics.regularPanels == 0 && statistics.transformedPanels > 0,"A singular child lost its mapping.");
  IntervalContractRequire(statistics.emergency7 == 0 && statistics.emergency3 == 0,"An emergency rule ran without collapse.");
  if(mode == 0 && (leftExponent > 1 || rightExponent > 1)) IntervalContractRequire(context.visits[0] < context.visits[10],"Completed moments were resampled in children.");
  std::cout << "analytic," << dimension << ',' << leftExponent << ',' << rightExponent << ',' << mode << ',' << maximum << ',' << statistics.regularPanels << ',' << statistics.transformedPanels << ',' << statistics.subdivisions << ',' << context.evaluations << '\n';
}

void IntervalContractEmergency(unsigned dimension, long double ulps, unsigned mode, bool failure, unsigned expectedRule) {
  FELevelSetMoments::OuterIntervalQuadrature::Statistics(dimension) = FELevelSetMoments::OuterIntervalStatistics();
  const long double delta = ulps*std::numeric_limits<long double>::epsilon();
  const FELevelSetMoments::OuterInterval interval = {1,1+delta,1,delta,1,2,true};
  IntervalContractContext context;
  context.mode = mode;
  context.delta = delta;
  context.endpoint = 1;
  std::vector<unsigned char> active(22,1);
  std::vector<long double> values(22,0);
  std::vector<long double> errors(22,0);
  unsigned accepted = 0;
  bool rejected = false;
  try { FELevelSetMoments::OuterIntervalQuadrature::Adaptive(context,IntervalContractSample,interval,0,1,1e-12L*delta,1e-12L*delta,1e-12L,0,30,10000,context.evaluations,accepted,dimension,active,values,errors); }
  catch(const std::runtime_error&) { rejected = true; }
  IntervalContractRequire(rejected == failure,"Incorrect emergency-rule acceptance.");
  const FELevelSetMoments::OuterIntervalStatistics& statistics = FELevelSetMoments::OuterIntervalQuadrature::Statistics(dimension);
  IntervalContractRequire(statistics.coordinateCollapses == 1 && statistics.subdivisions == 0,"A collapsed interval was subdivided.");
  if(expectedRule == 7) IntervalContractRequire(statistics.emergency7 == 1 && statistics.emergency3 == 0,"The first resolvable emergency rule was skipped.");
  if(expectedRule == 3) IntervalContractRequire(statistics.emergency7 == 0 && statistics.emergency3 == 1,"The last resolvable emergency rule was skipped.");
  if(!failure) IntervalContractRequire(values[0] > 0 && std::abs(values[0]-delta) <= 1e-12L*delta,"An emergency cut was lost.");
  std::cout << "emergency," << dimension << ',' << ulps << ',' << mode << ',' << rejected << ',' << statistics.emergency7 << ',' << statistics.emergency3 << ',' << statistics.resolutionLimits << '\n';
}

void IntervalContractBudgets(unsigned dimension) {
  std::vector<long double> values(2,1);
  std::vector<long double> errors(2,0.2L);
  std::vector<unsigned char> active(2,1);
  std::vector<long double> total(2,0);
  std::vector<long double> totalErrors(2,0);
  IntervalContractRequire(FELevelSetMoments::OuterIntervalQuadrature::Accept(values,errors,0.25L,0.1L,0,active,total,totalErrors,dimension),"Independent surface tolerance was lost.");
  IntervalContractRequire(active[0] == 0 && active[1] == 1 && total[0] == 1 && total[1] == 0,"Volume and surface acceptance were coupled.");
  errors[1] = 0.05L;
  IntervalContractRequire(!FELevelSetMoments::OuterIntervalQuadrature::Accept(values,errors,0.25L,0.1L,0,active,total,totalErrors,dimension) && total[0] == 1 && total[1] == 1,"A completed volume moment was accepted twice.");
}

int main() {
  std::cout << std::setprecision(18);
  try {
    for(unsigned dimension = 2; dimension <= 3; ++dimension) {
      IntervalContractBudgets(dimension);
      IntervalContractRun(1,1,0,dimension);
      IntervalContractRun(2,1,1,dimension);
      IntervalContractRun(1,2,0,dimension);
      IntervalContractRun(2,2,2,dimension);
      IntervalContractRun(4,1,3,dimension);
      IntervalContractRun(6,1,6,dimension);
      IntervalContractEmergency(dimension,4096,4,false,7);
      IntervalContractEmergency(dimension,128,4,false,3);
      IntervalContractEmergency(dimension,1,4,true,0);
      IntervalContractEmergency(dimension,4096,5,true,0);
      IntervalContractEmergency(dimension,4096,7,true,0);
      IntervalContractEmergency(dimension,64,4,false,3);
      IntervalContractEmergency(dimension,64,7,true,0);
    }
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
