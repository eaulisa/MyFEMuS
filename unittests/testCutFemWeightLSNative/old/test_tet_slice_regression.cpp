#include "Tet15LevelSetMoments.hpp"
#include "Basis.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <vector>

unsigned sliceChecks = 0;
unsigned sliceFailures = 0;
unsigned sliceRejectedCharts = 0;
long double sliceMaximumError = 0.0L;

void CheckSlice(bool condition, const char* description) {
  ++sliceChecks;
  if(condition) return;
  ++sliceFailures;
  if(sliceFailures <= 6) std::cerr << description << '\n';
}

long double SlicePower(long double value, unsigned degree) {
  long double result = 1.0L;
  for(unsigned k = 0; k < degree; ++k) result *= value;
  return result;
}

void CheckSliceValue(long double actual, long double expected, long double tolerance) {
  const long double error = std::abs(actual-expected);
  sliceMaximumError = std::max(sliceMaximumError,error);
  CheckSlice(std::isfinite(actual) && error <= tolerance,"A first-slice moment differs from its geometric reference.");
}

void CheckTetOutputRounding() {
  femus::TetBiquadratic basis;
  std::vector<double> nodal(basis._nc);
  for(unsigned node = 0; node < nodal.size(); ++node) nodal[node] = basis.GetXcoarse(node)[0]+1.0;
  std::vector<std::vector<std::vector<double>>> positive(1,std::vector<std::vector<double>>(1,std::vector<double>(1,17.0)));
  std::vector<std::vector<std::vector<double>>> delta(1,std::vector<std::vector<double>>(1,std::vector<double>(1,19.0)));
  std::vector<double> estimated;
  std::ostringstream warnings;
  std::streambuf* previous = std::cerr.rdbuf(warnings.rdbuf());
  bool rejected = false;
  try { FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,0,positive,delta,8.0e-18,0.0,true,24,0,200000,&estimated,0); }
  catch(const std::exception&) { rejected = true; }
  std::cerr.rdbuf(previous);
  CheckSlice(!rejected,"An unattainable output tolerance aborted instead of returning its quadrature.");
  CheckSliceValue(positive[0][0][0],1.0L/6.0L,2.0e-17L);
  CheckSlice(delta[0][0][0] == 0,"Uncertified full-cell quadrature acquired an interface mass.");
  CheckSlice(estimated.size() >= 2 && estimated[0] > 8.0e-18,"Unmet output accuracy was hidden in the error estimate.");
  CheckSlice(warnings.str().find("rounded volume moment exceeds tolerance") != std::string::npos,"An unattainable output tolerance omitted its rounding warning.");
  bool completed = false;
  try {
    FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,0,positive,delta,2.0e-17,0.0,true,24,0,200000,&estimated,0);
    completed = true;
  }
  catch(const std::exception&) { CheckSlice(false,"An attainable final double-output tolerance was rejected."); }
  if(!completed) return;
  CheckSliceValue(positive[0][0][0],1.0L/6.0L,2.0e-17L);
  CheckSlice(delta[0][0][0] == 0.0,"A full positive TET acquired an interface mass.");
  CheckSlice(estimated.size() >= 2,"A completed TET rule omitted its error estimates.");
  if(estimated.size() < 2) return;
  const long double outputError = std::abs(static_cast<long double>(positive[0][0][0])-1.0L/6.0L);
  CheckSlice(std::isfinite(estimated[0]) && estimated[0] >= outputError && estimated[0] <= 2.0e-17,"The reported TET error omits final double rounding.");
  CheckSlice(std::isfinite(estimated[1]) && estimated[1] >= 0.0 && estimated[1] <= 2.0e-17,"A full positive TET has an invalid delta error estimate.");
}

void CheckThinCapSlice() {
  femus::TetBiquadratic basis;
  std::vector<double> nodal(15);
  const long double cut = 1.0L-1.0L/256;
  for(unsigned node = 0; node < nodal.size(); ++node) nodal[node] = basis.GetXcoarse(node)[0]-static_cast<double>(cut);
  std::vector<long double> original;
  FELevelSetMoments::Tet15PolynomialFromNodalValues(nodal,false,original);
  std::vector<long double> coefficient(27);
  bool enriched = false;
  for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) { coefficient[9*i+3*j+k] = original[25*i+5*j+k]; enriched = enriched || (i+j+k > 2 && coefficient[9*i+3*j+k] != 0); }
  CheckSlice(enriched,"The thin-cap fixture lost its FE enrichment");
  const long double parameter = 0.25L*std::acos(-1.0L)*(1-0.991455371120812639206854697526329L);
  const long double height = SlicePower(std::sin(parameter),2)/256;
  const long double width = 1-height;
  std::vector<std::vector<long double>> chart;
  FELevelSetMoments::ReferencePolynomialAffine2D(FELevelSetMoments::SurfaceEndpointIntegration::SlicePolynomial(coefficient,1,height),0,width,0,width,chart);
  const unsigned degree = 5;
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  for(unsigned requested = 0; requested <= 2; requested += 2) {
    std::vector<long double> values;
    std::vector<long double> errors;
    unsigned evaluations = 0;
    unsigned accepted = 0;
    unsigned selected = 0;
    FELevelSetMoments::TriangleIntegration::Integrate(chart,degree,true,1e-13L,1e-13L,1e-13L,32,200000,values,errors,evaluations,accepted,selected,nullptr,nullptr,nullptr,requested);
    const long double root = cut/width;
    const long double thickness = 1-root;
    CheckSlice(values[count] > 0,"The thin-cap slice lost its interface contribution");
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
      long double h = 0;
      for(unsigned k = 0; k <= i; ++k) h += FELevelSetMoments::ReferenceBinomial(i,k)*(k%2U == 0 ? 1 : -1)*SlicePower(thickness,j+k+2U)/((j+1U)*(j+k+2U));
      const long double d = SlicePower(root,i)*SlicePower(thickness,j+1U)/((j+1U)*width);
      const unsigned n = FELevelSetMoments::ReferenceMomentIndex2D(degree,i,j);
      CheckSliceValue(values[n],h,1e-13L);
      CheckSliceValue(values[count+n],d,1e-13L);
      CheckSlice(errors[n] <= 1e-13L && errors[count+n] <= 1e-13L,"A thin-cap slice exceeded its error budget");
    }
  }
}

int main() {
  try { CheckThinCapSlice(); CheckTetOutputRounding(); }
  catch(const std::exception& error) { CheckSlice(false,error.what()); }
  std::cout << "Thin-cap and output-rounding checks=" << sliceChecks << " failures=" << sliceFailures << " maximum_error=" << sliceMaximumError << '\n';
  return sliceFailures == 0 ? 0 : 1;
}
