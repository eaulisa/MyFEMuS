#include "LineIntegration.hpp"
#include "TriangleIntegration.hpp"
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
  for(unsigned degree = 8; degree <= 10; degree += 2) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    std::vector<long double> values(2U*count,0);
    std::vector<long double> errors(2U*count,0);
    const FELevelSetMoments::LineIntegrationStrip strip = {0,1,0,0,0,1};
    FELevelSetMoments::LineIntegration::FullStripMoments(degree,strip,0,1,values,errors);
    const unsigned n = FELevelSetMoments::ReferenceMomentIndex2D(degree,0,degree);
    ReviewCheck(values[n] > 0 && std::abs(values[n]-1.0L/((degree+1U)*(degree+2U))) <= errors[n],"Full-triangle high-degree moment and estimator");
  }
  const long double size = std::ldexp(1.0L,-30);
  std::vector<std::vector<long double>> field(3,std::vector<long double>(3,0));
  field[0][0] = 0.5L*size;
  field[1][0] = field[0][1] = -size;
  const std::array<long double,5> small = {0.125L+size,-size,0.125L,size,size*size};
  std::vector<long double> values;
  std::vector<long double> errors;
  unsigned evaluations = 0;
  unsigned accepted = 0;
  unsigned direction = 0;
  FELevelSetMoments::TriangleIntegration::Integrate(field,0,true,1e-15L,1e-15L,1e-12L,28,200000,values,errors,evaluations,accepted,direction,nullptr,&small);
  CheckMoment(values[1],0.5L*size,1e-15L,"Tiny-triangle delta");
  ReviewCheck(values[0] > 0,"Tiny positive triangle contribution was discarded");
  field[0][0] = -0.5L;
  field[0][1] = 1;
  field[1][0] = 0;
  const std::array<long double,5> parent = {1,-1,0,1,1};
  const unsigned degree = 6;
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  FELevelSetMoments::TriangleIntegration::Integrate(field,degree,true,1e-12L,1e-12L,1e-12L,28,200000,values,errors,evaluations,accepted,direction,nullptr,&parent);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
    const unsigned n = FELevelSetMoments::ReferenceMomentIndex2D(degree,i,j);
    CheckMoment(values[count+n],ReviewPower(0.5L,j)*(1-ReviewPower(0.5L,i+1U))/(i+1U),1e-12L,"Full-triangle trace ownership");
  }
  std::cout << "REVIEW full triangle, tiny region, and degree-8/10 moments passed\n";
}

void CheckSelectiveRefinement() {
  const unsigned degree = 2;
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  std::vector<std::vector<long double>> coefficient(3,std::vector<long double>(3,0.0L));
  coefficient[0][0] = -0.25L;
  coefficient[0][1] = 1.0L;
  coefficient[2][1] = 100.0L;
  const FELevelSetMoments::LineIntegrationStrip strip = {0.0L,1.0L,0.0L,0.0L,1.0L,0.0L};
  for(unsigned repeat = 0; repeat < 1; ++repeat) {
    std::vector<long double> values(2U*count,0.0L);
    std::vector<long double> errors(2U*count,0.0L);
    std::vector<unsigned char> active(2U*count,0);
    values[0] = 123.0L;
    active[count] = 1;
    unsigned evaluations = 0;
    unsigned accepted = 0;
    FELevelSetMoments::LineIntegration::Statistics() = FELevelSetMoments::LineIntegrationStatistics();
    FELevelSetMoments::LineIntegration::Adaptive(coefficient,degree,strip,0.0L,1.0L,0.0L,1.0L,false,true,true,false,false,1.0e-13L,1.0e-13L,1.0e-13L,8,24,200000,evaluations,accepted,values,errors,nullptr,active);
    ReviewCheck(values[0] == 123.0L,"Accepted moment was recomputed during selective refinement");
    CheckMoment(values[count],std::atan(10.0L)/10.0L,1.0e-13L,"Selective graph integral");
  }
  std::cout << "REVIEW selective refinement regressions passed\n";
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
  FELevelSetMoments::IntegrateQuad9ReferenceSelected(field,degree,positive,delta,1.0e-12,1.0e-12,true,false,28,0,200000,&estimate,2,nullptr);
  ReviewCheck(estimate.size() >= 7 && estimate[5] != 0.0,"QUAD validation mode 2 did not perform an independent comparison");
  field.assign(3,std::vector<double>(3,0.0));
  field[0][0] = -0.25;
  field[1][0] = 1.0;
  field[0][1] = 1.0;
  estimate.clear();
  FELevelSetMoments::IntegrateTri7ReferenceSelected(field,degree,positive,delta,1.0e-12,1.0e-12,true,false,28,0,200000,&estimate,2,nullptr);
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
      if(shape == 0) FELevelSetMoments::IntegrateQuad9ReferenceSelected(field,degree,positive,delta,tolerances[t],tolerances[t],true,false,28,0,200000,&estimate,1,nullptr);
      else FELevelSetMoments::IntegrateTri7ReferenceSelected(field,degree,positive,delta,tolerances[t],tolerances[t],true,false,28,0,200000,&estimate,1,nullptr);
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

void CheckMappedDiscriminant() {
  std::vector<std::vector<long double>> field(3,std::vector<long double>(3,0.0L));
  field[1][0] = -1.0L;
  field[2][0] = 1.0L;
  field[0][2] = 1.0L;
  field[2][2] = 1.0L;
  const FELevelSetMoments::LineIntegrationDiscriminant discriminant = FELevelSetMoments::LineIntegration::MappedDiscriminant(field,0.0L,1.0L);
  ReviewCheck(discriminant.leftZero && discriminant.rightZero,"Mapped quartic tangencies were not retained");
  const long double parameters[7] = {std::ldexp(1.0L,-20),0.001L,0.1L,0.5L,0.9L,0.999L,1.0L-std::ldexp(1.0L,-20)};
  for(unsigned k = 0; k < 7; ++k) {
    const long double sine = std::sin(0.5L*std::acos(-1.0L)*parameters[k]);
    const long double cosine = std::cos(0.5L*std::acos(-1.0L)*parameters[k]);
    const long double point = sine*sine;
    const long double exact = 2.0L*sine*cosine*std::sqrt(1.0L+point*point);
    CheckMoment(FELevelSetMoments::LineIntegration::MappedDerivativeMagnitude(discriminant,parameters[k]),exact,1.0e-17L,"Factored mapped quartic discriminant");
  }
  std::cout << "REVIEW mapped quartic endpoint arithmetic passed\n";
}

void CheckTranslatedEllipses() {
  const unsigned degrees[2] = {6U,10U};
  unsigned cases = 0;
  for(unsigned shape = 0; shape < 2; ++shape) for(unsigned geometry = 0; geometry < 6; ++geometry) for(unsigned d = 0; d < 2; ++d) for(unsigned sign = 0; sign < 2; ++sign) {
    const unsigned degree = degrees[d];
    const long double centerX = 0.213L+0.017L*geometry;
    const long double centerY = 0.187L+0.011L*geometry;
    const long double radius = geometry%2U == 0 ? 0.037L : 0.0037L;
    const long double a = 0.7L+0.13L*geometry;
    const long double b = 1.4L+0.17L*geometry;
    const std::string label = std::string(shape == 0 ? "QUAD" : "TRI")+" ellipse geometry="+std::to_string(geometry)+" degree="+std::to_string(degree)+" sign="+std::to_string(sign);
    const long double multiplier = sign == 0 ? 1.0L : -8.0L;
    std::vector<std::vector<double>> field(3,std::vector<double>(3,0.0));
    field[0][0] = static_cast<double>(multiplier*(a*centerX*centerX+b*centerY*centerY-radius*radius));
    field[1][0] = static_cast<double>(-2.0L*multiplier*a*centerX);
    field[0][1] = static_cast<double>(-2.0L*multiplier*b*centerY);
    field[2][0] = static_cast<double>(multiplier*a);
    field[0][2] = static_cast<double>(multiplier*b);
    const long double aa = field[2][0]/multiplier;
    const long double bb = field[0][2]/multiplier;
    const long double cx = -field[1][0]/(2.0L*multiplier*aa);
    const long double cy = -field[0][1]/(2.0L*multiplier*bb);
    const long double squaredRadius = std::fma(aa*cx,cx,std::fma(bb*cy,cy,-field[0][0]/multiplier));
    const long double rx = std::sqrt(squaredRadius/aa);
    const long double ry = std::sqrt(squaredRadius/bb);
    std::vector<std::vector<double>> positive;
    std::vector<std::vector<double>> delta;
    if(shape == 0) FELevelSetMoments::IntegrateQuad9ReferenceSelected(field,degree,positive,delta,1.0e-12,1.0e-12,true,false,28,0,200000,nullptr,1,nullptr);
    else FELevelSetMoments::IntegrateTri7ReferenceSelected(field,degree,positive,delta,1.0e-12,1.0e-12,true,false,28,0,200000,nullptr,1,nullptr);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
      long double disk = 0.0L;
      long double interface = 0.0L;
      for(unsigned u = 0; u <= i; ++u) for(unsigned v = 0; v <= j; ++v) {
        const long double term = FELevelSetMoments::ReferenceBinomial(i,u)*FELevelSetMoments::ReferenceBinomial(j,v)*ReviewPower(cx,i-u)*ReviewPower(cy,j-v)*ReviewPower(rx,u)*ReviewPower(ry,v)*AngularMoment(u,v);
        disk += term*rx*ry/(u+v+2U);
        interface += term/(2.0L*std::sqrt(aa*bb)*std::abs(multiplier));
      }
      const long double full = shape == 0 ? (i%2U == 0 && j%2U == 0 ? 4.0L/((i+1.0L)*(j+1.0L)) : 0.0L) : std::tgamma(i+1.0L)*std::tgamma(j+1.0L)/std::tgamma(i+j+3.0L);
      CheckMoment(positive[i][j],sign == 0 ? full-disk : disk,1.0e-12L,label+" H");
      CheckMoment(delta[i][j],interface,1.0e-12L,label+" D");
    }
    ++cases;
  }
  std::cout << "REVIEW translated ellipse cases=" << cases << " failures=0\n";
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
    if(shape == 0) FELevelSetMoments::IntegrateQuad9ReferenceSelected(field,6,positive,delta,1.0e-12,1.0e-12,true,false,28,0,200000,nullptr,1,nullptr);
    else FELevelSetMoments::IntegrateTri7ReferenceSelected(field,6,positive,delta,1.0e-12,1.0e-12,true,false,28,0,200000,nullptr,1,nullptr);
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
    CheckSelectiveRefinement();
    CheckExplicitValidationMode();
    CheckCircleConvergence();
    CheckMappedDiscriminant();
    CheckTranslatedEllipses();
    if(argc > 1) CheckAnalyticReferenceFile(argv[1]);
    std::cout << std::setprecision(16) << "REVIEW PASS checks=" << reviewChecks << " maximum_actual_error_fraction=" << maximumBudgetFraction << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count() << '\n';
    return 0;
  }
  catch(const std::exception& exception) {
    std::cerr << "REVIEW FAIL " << exception.what() << '\n';
    return 1;
  }
}
