#include "LineIntegration.hpp"
#define private public
#include "TriangleIntegration.hpp"
#undef private
#include "Quad9DirectionSelection.hpp"
#include "Tri7DirectionSelection.hpp"
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

unsigned reviewChecks = 0;
long double maximumBudgetFraction = 0.0L;

void ReviewCheck(bool condition, const std::string& name) {
  ++reviewChecks;
  if(!condition) throw std::runtime_error(name);
}

long double ReviewPower(long double value, unsigned degree) {
  long double answer = 1.0L;
  for(unsigned k = 0; k < degree; ++k) answer *= value;
  return answer;
}

long double AngularMoment(unsigned first, unsigned second) {
  if((first&1U) || (second&1U)) return 0.0L;
  return 2.0L*std::tgamma(0.5L*(first+1U))*std::tgamma(0.5L*(second+1U))/std::tgamma(0.5L*(first+second)+1.0L);
}

long double CircleMoment(unsigned i, unsigned j, long double cx, long double cy, long double radius, bool delta) {
  long double value = 0.0L;
  for(unsigned a = 0; a <= i; ++a) for(unsigned b = 0; b <= j; ++b) {
    const long double multiplier = FELevelSetMoments::ReferenceBinomial(i,a)*FELevelSetMoments::ReferenceBinomial(j,b)*ReviewPower(cx,i-a)*ReviewPower(cy,j-b)*AngularMoment(a,b);
    value += multiplier*ReviewPower(radius,a+b+(delta ? 0U : 2U))*(delta ? 0.5L : 1.0L/(a+b+2U));
  }
  return value;
}

void CheckMoment(long double computed, long double exact, long double tolerance, const std::string& label) {
  const long double error = std::abs(computed-exact);
  const long double budget = tolerance*(1.0L+std::abs(exact));
  maximumBudgetFraction = std::max(maximumBudgetFraction,error/budget);
  if(!std::isfinite(computed) || error > budget) {
    std::cerr << std::setprecision(21) << label << " value=" << computed << " exact=" << exact << " error=" << error << " budget=" << budget << '\n';
    throw std::runtime_error("Momentwise accuracy failure: "+label);
  }
  ++reviewChecks;
}

void CheckTriangleGeometry() {
  const std::vector<FELevelSetMoments::TriangleDomainVertex> triangle = {{1.0L,1.0L},{0.0L,0.0L},{1.0L,0.0L}};
  const unsigned geometryDegrees[2] = {8,10};
  for(unsigned degreeIndex = 0; degreeIndex < 2; ++degreeIndex) {
    const unsigned degree = geometryDegrees[degreeIndex];
    long double estimate = 0.0L;
    const long double value = FELevelSetMoments::TriangleIntegration::TriangleMoment(triangle,0,degree,&estimate);
    const long double exact = 1.0L/((degree+1.0L)*(degree+2.0L));
    ReviewCheck(value > 0.0L && std::abs(value-exact) <= estimate,"Stable high-degree full-triangle moment and estimator");
  }
  const long double s = std::ldexp(1.0L,-30);
  const std::vector<FELevelSetMoments::TriangleDomainVertex> small = {{0.125L,0.125L},{0.125L+s,0.125L},{0.125L+s,0.125L+s}};
  std::vector<std::vector<FELevelSetMoments::TriangleDomainVertex>> retained;
  FELevelSetMoments::TriangleIntegration::TriangulatePolygon(small,retained);
  ReviewCheck(retained.size() == 1,"Tiny positive-area triangle was discarded");
  std::vector<std::vector<long double>> field(3,std::vector<long double>(3,0.0L));
  field[1][0] = 1.0L;
  field[0][1] = -1.0L;
  field[0][0] = -0.5L*s;
  std::vector<long double> values;
  std::vector<long double> errors;
  unsigned evaluations = 0;
  unsigned accepted = 0;
  unsigned direction = 0;
  FELevelSetMoments::TriangleIntegration::Integrate(field,0,small,true,1.0e-15L,1.0e-15L,1.0e-12L,28,200000,0,10,values,errors,evaluations,accepted,direction);
  CheckMoment(values[1],0.5L*s,1.0e-15L,"Tiny-triangle delta");
  field.assign(3,std::vector<long double>(3,0.0L));
  field[0][1] = 1.0L;
  field[0][0] = -0.5L;
  const std::vector<std::vector<FELevelSetMoments::TriangleDomainVertex>> children = {{{0.0L,0.0L},{1.0L,0.0L},{1.0L,0.5L}},{{0.0L,0.0L},{1.0L,0.5L},{0.5L,0.5L}},{{0.5L,0.5L},{1.0L,0.5L},{1.0L,1.0L}}};
  const unsigned degree = 6;
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  std::vector<long double> sum(2U*count,0.0L);
  for(unsigned child = 0; child < children.size(); ++child) {
    FELevelSetMoments::TriangleIntegration::Integrate(field,degree,children[child],true,1.0e-12L,1.0e-12L,1.0e-12L,28,200000,0,10,values,errors,evaluations,accepted,direction);
    for(unsigned k = 0; k < sum.size(); ++k) sum[k] += values[k];
  }
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
    const unsigned k = FELevelSetMoments::MomentIndex(degree,i,j);
    const long double exact = ReviewPower(0.5L,j)*(1.0L-ReviewPower(0.5L,i+1U))/(i+1U);
    CheckMoment(sum[count+k],exact,1.0e-12L,"Artificial-edge delta ownership");
  }
  std::cout << "REVIEW partition, tiny-region, degree-8/10 moment regressions passed\n";
}

void CheckFallbackState() {
  const unsigned degree = 2;
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  std::vector<std::vector<long double>> coefficient(3,std::vector<long double>(3,0.0L));
  coefficient[0][0] = -0.25L;
  coefficient[0][1] = 1.0L;
  coefficient[2][1] = 1.0L;
  const FELevelSetMoments::LineIntegrationStrip strip = {0.0L,1.0L,0.0L,0.0L,1.0L,0.0L};
  for(unsigned mode = 0; mode < 2; ++mode) {
    std::vector<long double> values(2U*count,0.0L);
    std::vector<long double> errors(2U*count,0.0L);
    std::vector<unsigned char> active(2U*count,0);
    values[0] = 123.0L;
    active[count] = 1;
    unsigned evaluations = 0;
    unsigned accepted = 0;
    FELevelSetMoments::LineIntegration::Statistics() = FELevelSetMoments::LineIntegrationStatistics();
    FELevelSetMoments::LineIntegration::Adaptive(coefficient,degree,strip,0.0L,1.0L,0.0L,1.0L,false,true,mode == 0,true,false,false,1.0e-13L,1.0e-13L,1.0e-13L,8,24,200000,evaluations,accepted,values,errors,nullptr,active);
    ReviewCheck(values[0] == 123.0L,"Accepted moment was recomputed during fallback");
    CheckMoment(values[count],0.25L*std::acos(-1.0L),1.0e-13L,"Selective fallback integral");
    if(mode == 0) ReviewCheck(FELevelSetMoments::LineIntegration::Statistics().taylorPanels > 0,"Taylor fallback was not exercised");
    else ReviewCheck(FELevelSetMoments::LineIntegration::Statistics().fullSamples > 0,"Full GK fallback was not exercised");
  }
  std::cout << "REVIEW selective cross-stage fallback regressions passed\n";
}

void CheckExplicitValidationMode() {
  const unsigned degree = 10;
  std::vector<std::vector<double>> field(3,std::vector<double>(3,0.0));
  std::vector<std::vector<double>> positive;
  std::vector<std::vector<double>> delta;
  std::vector<double> estimate;
  field[0][0] = -0.1;
  field[1][0] = 1.0;
  field[0][1] = 0.5;
  FELevelSetMoments::IntegrateQuad9ReferenceSelected(field,degree,positive,delta,1.0e-12,1.0e-12,true,false,28,0,200000,&estimate,2,nullptr,true);
  ReviewCheck(estimate.size() >= 7 && estimate[5] != 0.0,"QUAD validation mode 2 did not perform an independent comparison");
  field.assign(3,std::vector<double>(3,0.0));
  field[0][0] = -0.25;
  field[1][0] = 1.0;
  field[0][1] = 1.0;
  estimate.clear();
  FELevelSetMoments::IntegrateTri7ReferenceSelected(field,degree,positive,delta,1.0e-12,1.0e-12,true,false,28,0,200000,&estimate,2,nullptr,true);
  ReviewCheck(estimate.size() >= 7 && estimate[5] != 0.0,"TRI validation mode 2 did not perform an independent comparison");
  std::cout << "REVIEW explicit heavy validation mode passed on regular graphs\n";
}

void CheckCircleConvergence() {
  const unsigned degrees[4] = {0,2,6,10};
  const double tolerances[3] = {1.0e-8,1.0e-10,1.0e-12};
  unsigned cases = 0;
  unsigned failures = 0;
  for(unsigned shape = 0; shape < 2; ++shape) for(unsigned geometry = 0; geometry < 3; ++geometry) for(unsigned d = 0; d < 4; ++d) for(unsigned t = 0; t < 3; ++t) for(unsigned sign = 0; sign < 2; ++sign) {
    const unsigned degree = degrees[d];
    const long double cx = shape == 0 ? 0.125L : 0.25L;
    const long double cy = shape == 0 ? -0.125L : 0.25L;
    const long double radius = geometry == 0 ? 0.125L : (geometry == 1 ? 0.03125L : 0.00390625L);
    const long double scale = sign == 0 ? 1.0L : -8.0L;
    std::vector<std::vector<double>> field(3,std::vector<double>(3,0.0));
    field[0][0] = static_cast<double>(scale*(cx*cx+cy*cy-radius*radius));
    field[1][0] = static_cast<double>(-2.0L*scale*cx);
    field[0][1] = static_cast<double>(-2.0L*scale*cy);
    field[2][0] = static_cast<double>(scale);
    field[0][2] = static_cast<double>(scale);
    const std::vector<std::vector<double>> original = field;
    std::vector<std::vector<double>> positive;
    std::vector<std::vector<double>> delta;
    std::vector<double> estimate;
    const std::string label = std::string(shape == 0 ? "QUAD" : "TRI")+" circle geometry="+std::to_string(geometry)+" degree="+std::to_string(degree)+" tolerance="+std::to_string(t)+" sign="+std::to_string(sign);
    try {
      try {
      if(shape == 0) FELevelSetMoments::IntegrateQuad9ReferenceSelected(field,degree,positive,delta,tolerances[t],tolerances[t],true,false,28,0,200000,&estimate,1,nullptr,true);
      else FELevelSetMoments::IntegrateTri7ReferenceSelected(field,degree,positive,delta,tolerances[t],tolerances[t],true,false,28,0,200000,&estimate,1,nullptr,true);
    }
      catch(const std::exception& exception) { throw std::runtime_error(label+": "+exception.what()); }
    ReviewCheck(field == original,"Integration modified the parent coefficient tensor");
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
      const long double disk = CircleMoment(i,j,cx,cy,radius,false);
      const long double full = shape == 0 ? (i%2U == 0 && j%2U == 0 ? 4.0L/((i+1.0L)*(j+1.0L)) : 0.0L) : std::tgamma(i+1.0L)*std::tgamma(j+1.0L)/std::tgamma(i+j+3.0L);
      CheckMoment(positive[i][j],sign == 0 ? full-disk : disk,tolerances[t],label+" H");
      CheckMoment(delta[i][j],CircleMoment(i,j,cx,cy,radius,true)/std::abs(scale),tolerances[t],label+" D");
    }
    ReviewCheck(estimate.size() >= 7,"Difficult cut did not report integration diagnostics");
    ++cases;
    }
    catch(const std::exception& exception) {
      std::cerr << "REVIEW CIRCLE CASE FAIL " << label << ": " << exception.what() << '\n';
      ++failures;
    }
  }
  std::cout << "REVIEW circle convergence cases=" << cases << " failures=" << failures << " tolerances=1e-8,1e-10,1e-12\n";
  if(failures != 0) throw std::runtime_error("Circle convergence failures="+std::to_string(failures));
}

void CheckAnalyticReferenceFile(const char* path) {
  std::ifstream input(path);
  if(!input) throw std::runtime_error("Cannot open independent analytic reference table");
  unsigned shape;
  unsigned i;
  unsigned j;
  double beta;
  long double h;
  long double d;
  unsigned cases = 0;
  while(input >> shape >> beta) {
    std::vector<std::vector<double>> field(3,std::vector<double>(3,0.0));
    if(shape == 0) {
      field[0][0] = -0.25;
      field[2][0] = 0.125;
      field[0][1] = 1.0;
      field[2][1] = beta;
    }
    else {
      field[0][0] = -0.25;
      field[0][1] = 1.0;
      field[1][1] = beta;
      field[2][1] = -beta;
      field[1][2] = -beta;
    }
    std::vector<std::vector<double>> positive;
    std::vector<std::vector<double>> delta;
    if(shape == 0) FELevelSetMoments::IntegrateQuad9ReferenceSelected(field,6,positive,delta,1.0e-12,1.0e-12,true,false,28,0,200000,nullptr,1,nullptr,true);
    else FELevelSetMoments::IntegrateTri7ReferenceSelected(field,6,positive,delta,1.0e-12,1.0e-12,true,false,28,0,200000,nullptr,1,nullptr,true);
    for(unsigned k = 0; k < 28; ++k) {
      if(!(input >> i >> j >> h >> d)) throw std::runtime_error("Malformed reference table");
      CheckMoment(positive[i][j],h,1.0e-12L,"Independent high-precision graph H");
      CheckMoment(delta[i][j],d,1.0e-12L,"Independent high-precision graph D");
    }
    ++cases;
  }
  ReviewCheck(cases == 4,"Independent reference cases missing");
  std::cout << "REVIEW independent 60-digit graph and TRI bubble references passed\n";
}

int main(int argc, char** argv) {
  try {
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    CheckTriangleGeometry();
    CheckFallbackState();
    CheckExplicitValidationMode();
    CheckCircleConvergence();
    if(argc > 1) CheckAnalyticReferenceFile(argv[1]);
    std::cout << std::setprecision(16) << "REVIEW PASS checks=" << reviewChecks << " maximum_actual_error_fraction=" << maximumBudgetFraction << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count() << '\n';
    return 0;
  }
  catch(const std::exception& exception) {
    std::cerr << "REVIEW FAIL " << exception.what() << '\n';
    return 1;
  }
}
