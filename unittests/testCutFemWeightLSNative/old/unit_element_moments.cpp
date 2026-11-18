#include "CutFemWeight_LS.hpp"
#include "Basis.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

unsigned unitChecks = 0;
long double unitMomentError = 0;
long double unitWeightError = 0;

void UnitCheck(long double actual, long double expected, long double tolerance) {
  ++unitChecks;
  if(!std::isfinite(actual) || std::abs(actual-expected) > tolerance) { std::cerr << std::setprecision(18) << "actual=" << actual << " expected=" << expected << " tolerance=" << tolerance << '\n'; throw std::runtime_error("A unit-element check failed."); }
}

void UnitProbe(unsigned dimension, unsigned degree) {
  femus::LineBiquadratic line;
  femus::QuadBiquadratic quad;
  femus::HexBiquadratic hex;
  femus::basis* basis = dimension == 1 ? static_cast<femus::basis*>(&line) : (dimension == 2 ? static_cast<femus::basis*>(&quad) : static_cast<femus::basis*>(&hex));
  const unsigned nodes = dimension == 1 ? 3U : (dimension == 2 ? 9U : 27U);
  std::vector<double> nodal(nodes);
  for(unsigned n = 0; n < nodes; ++n) nodal[n] = (basis->GetXcoarse(n)[dimension-1U]+1)/2-0.375;
  const unsigned count = dimension == 1 ? degree+1U : (dimension == 2 ? (degree+1U)*(degree+2U)/2U : (degree+1U)*(degree+2U)*(degree+3U)/6U);
  std::vector<long double> expected(2U*count,1);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= (dimension > 1 ? degree-i : 0U); ++j) for(unsigned k = 0; k <= (dimension > 2 ? degree-i-j : 0U); ++k) {
    const unsigned powers[3] = {i,j,k};
    const unsigned index = FELevelSetMoments::UnitMomentIndex(dimension,i,j,k);
    for(unsigned d = 0; d < dimension; ++d) {
      expected[index] *= d == dimension-1U ? (1-FELevelSetMoments::ReferencePower(0.375L,powers[d]+1U))/(powers[d]+1U) : 1.0L/(powers[d]+1U);
      expected[count+index] *= d == dimension-1U ? FELevelSetMoments::ReferencePower(0.375L,powers[d]) : 1.0L/(powers[d]+1U);
    }
  }
  const double scales[3] = {1,0.00000095367431640625,1048576};
  for(unsigned scale = 0; scale < 3; ++scale) {
    std::vector<double> input = nodal;
    for(unsigned n = 0; n < nodes; ++n) input[n] *= scales[scale];
    std::vector<long double> values;
    std::vector<long double> errors;
    FELevelSetMoments::IntegrateUnitTensorNodalMoments(input,dimension,degree,values,errors,1e-12L,1e-12L,true);
    for(unsigned n = 0; n < values.size(); ++n) { unitMomentError = std::max(unitMomentError,std::abs(values[n]-expected[n])); UnitCheck(values[n],expected[n],2e-12L); }
  }
  const GeomElType element = dimension == 1 ? LINE : (dimension == 2 ? QUAD : HEX);
  const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
  CutFemWeight_LS<long double,long double> rule(element,degree,"legendre");
  std::vector<long double> volume;
  std::vector<long double> surface;
  rule.GetWeightsFromNodalValues(nodal,volume,surface,1e-12,1e-12,true);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= (dimension > 1 ? degree-i : 0U); ++j) for(unsigned k = 0; k <= (dimension > 2 ? degree-i-j : 0U); ++k) {
    const unsigned powers[3] = {i,j,k};
    const unsigned index = FELevelSetMoments::UnitMomentIndex(dimension,i,j,k);
    long double h = 0;
    long double s = 0;
    for(unsigned point = 0; point < volume.size(); ++point) {
      long double monomial = 1;
      for(unsigned d = 0; d < dimension; ++d) monomial *= FELevelSetMoments::ReferencePower((static_cast<long double>(rule.GetGaussCoordinatePointer(d)[point])+1)/2,powers[d]);
      h += rule.GetGaussWeightPointer()[point]*monomial*volume[point];
      s += rule.GetGaussWeightPointer()[point]*monomial*surface[point];
    }
    const long double exactH = (1U<<dimension)*expected[index];
    const long double exactS = (1U<<(dimension-1U))*expected[count+index];
    unitWeightError = std::max(unitWeightError,std::max(std::abs(h-exactH),std::abs(s-exactS)));
    UnitCheck(h,exactH,1e-10L);
    UnitCheck(s,exactS,1e-10L);
  }
  std::cout << "unit fit dimension=" << dimension << " degree=" << degree << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count() << '\n';
  for(unsigned n = 0; n < nodes; ++n) nodal[n] = 0.03*(static_cast<int>((n*17U+3U)%19U)-9);
  std::vector<long double> coefficient;
  FELevelSetMoments::UnitTensorPolynomialFromNodalValues(nodal,dimension,coefficient);
  for(unsigned n = 0; n < nodes; ++n) {
    std::array<long double,3> point = {0,0,0};
    for(unsigned d = 0; d < dimension; ++d) point[d] = (basis->GetXcoarse(n)[d]+1)/2;
    UnitCheck(FELevelSetMoments::SurfaceEndpointIntegration::Value(coefficient,point),nodal[n],2e-15L);
  }
  for(unsigned sample = 0; sample < 13; ++sample) {
    std::vector<double> native(dimension);
    std::array<long double,3> point = {0,0,0};
    for(unsigned d = 0; d < dimension; ++d) { point[d] = (1.0L+(sample*(d+3U)+d*5U)%29U)/31; native[d] = 2*point[d]-1; }
    long double value = 0;
    std::vector<long double> gradient(dimension,0);
    for(unsigned n = 0; n < nodes; ++n) {
      value += nodal[n]*basis->eval_phi(n,native);
      for(unsigned d = 0; d < dimension; ++d) gradient[d] += 2*nodal[n]*basis->eval_dphidxyz(d,basis->GetIND(n),native.data());
    }
    UnitCheck(FELevelSetMoments::SurfaceEndpointIntegration::Value(coefficient,point),value,3e-14L);
    for(unsigned d = 0; d < dimension; ++d) UnitCheck(FELevelSetMoments::SurfaceEndpointIntegration::Value(FELevelSetMoments::SurfaceEndpointIntegration::Derivative(coefficient,d),point),gradient[d],6e-14L);
  }
}

void UnitBoundaryProbe(unsigned dimension) {
  femus::LineBiquadratic line;
  femus::QuadBiquadratic quad;
  femus::HexBiquadratic hex;
  femus::basis* basis = dimension == 1 ? static_cast<femus::basis*>(&line) : (dimension == 2 ? static_cast<femus::basis*>(&quad) : static_cast<femus::basis*>(&hex));
  const unsigned nodes = dimension == 1 ? 3U : (dimension == 2 ? 9U : 27U);
  for(unsigned axis = 0; axis < dimension; ++axis) for(unsigned side = 0; side < 2; ++side) {
    std::vector<double> nodal(nodes);
    for(unsigned n = 0; n < nodes; ++n) nodal[n] = (basis->GetXcoarse(n)[axis]+1)/2-side;
    std::vector<long double> values;
    std::vector<long double> errors;
    FELevelSetMoments::IntegrateUnitTensorNodalMoments(nodal,dimension,2,values,errors,1e-10L,1e-10L,true,nullptr,0,2);
    UnitCheck(values[0],side == 0 ? 1 : 0,1e-10L);
    UnitCheck(values[values.size()/2U],0.5L,1e-10L);
  }
}

int main() {
  try {
    for(unsigned dimension = 1; dimension <= 3; ++dimension) UnitBoundaryProbe(dimension);
    for(unsigned dimension = 1; dimension <= 3; ++dimension) for(unsigned degree = 2; degree <= 10; degree += degree == 2 ? 3U : 5U) UnitProbe(dimension,degree);
    std::cout << std::setprecision(16) << unitChecks << " unit-element checks passed; moment_error=" << unitMomentError << " fitted_error=" << unitWeightError << '\n';
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
