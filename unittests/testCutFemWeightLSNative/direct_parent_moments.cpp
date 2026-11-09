#include "LineIntegration.hpp"
#include <iostream>

unsigned directChecks = 0;
long double directMaximumError = 0.0L;

long double DirectInterval(unsigned degree, long double origin, long double scale, long double left, long double right) {
  return (FELevelSetMoments::ReferencePower(origin+scale*right,degree+1U)-FELevelSetMoments::ReferencePower(origin+scale*left,degree+1U))/(scale*(degree+1U));
}

void DirectCheck(long double actual, long double expected) {
  const long double error = std::abs(actual-expected);
  directMaximumError = std::max(directMaximumError,error);
  if(!std::isfinite(actual) || error > 2.0e-12L) throw std::runtime_error("A direct parent-coordinate moment is incorrect.");
  ++directChecks;
}

void DirectProbe(const std::array<long double,5>& coordinates, bool sparse, bool taylor) {
  const unsigned degree = taylor ? 2U : 10U;
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  const long double cut = 0.6L;
  const FELevelSetMoments::LineIntegrationStrip strip = {0,1,0,0,1,0};
  std::vector<std::vector<long double>> coefficient(3,std::vector<long double>(3,0.0L));
  coefficient[0][0] = -cut;
  coefficient[0][1] = 1.0L;
  std::vector<unsigned char> requested(2U*count,sparse ? 0 : 1);
  if(sparse) {
    requested[FELevelSetMoments::ReferenceMomentIndex2D(degree,3,2)] = 1;
    requested[count+FELevelSetMoments::ReferenceMomentIndex2D(degree,2,4)] = 1;
  }
  std::vector<long double> values(2U*count,0.0L);
  std::vector<long double> errors(2U*count,0.0L);
  unsigned evaluations = 0;
  unsigned accepted = 0;
  if(taylor) FELevelSetMoments::LineIntegration::TaylorPanel(coefficient,degree,strip,0,1,true,values,errors,nullptr,&requested,&coordinates);
  else FELevelSetMoments::LineIntegration::IntegrateStrip(coefficient,degree,strip,true,false,true,true,1.0e-12L,1.0e-12L,1.0e-12L,24,200000,evaluations,accepted,values,errors,nullptr,&coordinates,&requested);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
    const unsigned index = FELevelSetMoments::ReferenceMomentIndex2D(degree,i,j);
    const long double first = DirectInterval(i,coordinates[0],coordinates[1],0,1);
    const long double h = coordinates[4]*first*DirectInterval(j,coordinates[2],coordinates[3],cut,1);
    const long double d = coordinates[4]*first*FELevelSetMoments::ReferencePower(coordinates[2]+coordinates[3]*cut,j);
    DirectCheck(values[index],requested[index] ? h : 0.0L);
    DirectCheck(values[count+index],requested[count+index] ? d : 0.0L);
  }
}

int main() {
  try {
    DirectProbe({-1,2,-1,2,4},false,false);
    DirectProbe({1,-2,1,-2,4},false,false);
    DirectProbe({-0.75L,2,-1,2,4},true,false);
    DirectProbe({-1,2,-1,2,4},false,true);
    std::cout << directChecks << " direct parent-coordinate checks passed; maximum absolute error " << directMaximumError << '\n';
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
