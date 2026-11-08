#include "Tet15LevelSetMoments.hpp"
#include "Basis.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
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

void CheckCorrectedIntervals() {
  const std::vector<std::vector<long double>> field = {{0.0L},{0.0L},{1.0L}};
  const long double left = 0.00390607476535280978686L;
  const long double right = 0.00390607476535280979448L;
  const long double leftOffset = 7.55356123471419197954e-21L;
  const long double rightOffset = -7.3894370506072447955e-21L;
  CheckSlice((right-left)+(rightOffset-leftOffset) < 0.0L,"The recorded reversed interval was not reproduced.");
  bool rejected = false;
  std::vector<long double> values;
  try { FELevelSetMoments::LineSliceSample(field,{}, {},0,0,left,right,leftOffset,rightOffset,0.5L,true,0,1,false,false,false,false,values); }
  catch(const std::exception&) { rejected = true; }
  CheckSlice(rejected,"A reversed corrected interval was accepted silently.");
  rejected = false;
  try { FELevelSetMoments::LineSliceSample(field,{}, {},0,0,left,left,0.0L,0.0L,0.5L,true,0,1,false,false,false,false,values); }
  catch(const std::exception&) { rejected = true; }
  CheckSlice(rejected,"A zero-width corrected interval was accepted silently.");
  const long double thinRight = std::nextafter(0.5L,1.0L);
  const long double width = thinRight-0.5L;
  try {
    FELevelSetMoments::LineSliceSample(field,{}, {},0,0,0.5L,thinRight,0.0L,0.0L,0.5L,true,0,1,false,false,false,false,values);
    CheckSlice(values.size() == 2 && std::isfinite(values[0]) && values[0] > 0.0L,"A valid thin interval lost its positive measure.");
    if(values.size() == 2) CheckSliceValue(values[0],width,8.0L*std::numeric_limits<long double>::epsilon()*width);
  }
  catch(const std::exception&) { CheckSlice(false,"A valid thin interval was rejected."); }
}

void BuildFirstSlice(std::vector<std::vector<long double>>& normalized, long double& scale, long double& height) {
  femus::TetBiquadratic basis;
  std::vector<double> nodal(basis._nc);
  for(unsigned node = 0; node < nodal.size(); ++node) nodal[node] = basis.GetXcoarse(node)[0]-(1.0-1.0/256.0);
  std::vector<std::vector<long double>> field;
  std::vector<std::vector<long double>> unit;
  FELevelSetMoments::Tet15SliceField(nodal,0,0.5L,field,unit);
  const long double halfWidth = 0.25L*std::acos(-1.0L);
  const long double parameter = halfWidth-halfWidth*0.991455371120812639206854697526329L;
  height = SlicePower(std::sin(parameter),2)/256.0L;
  std::vector<long double> coefficient(7);
  scale = 0.0L;
  for(unsigned term = 0; term < 7; ++term) {
    coefficient[term] = FELevelSetMoments::PolynomialValue(unit[term],height);
    scale = std::max(scale,std::abs(coefficient[term]));
  }
  int exponent = 0;
  std::frexp(scale,&exponent);
  scale = std::ldexp(1.0L,exponent-1);
  normalized.assign(3,std::vector<long double>(3,0.0L));
  normalized[0][0] = coefficient[0]/scale;
  normalized[1][0] = coefficient[1]/scale;
  normalized[0][1] = coefficient[2]/scale;
  normalized[2][0] = coefficient[3]/scale;
  normalized[1][1] = coefficient[4]/scale;
  normalized[0][2] = coefficient[5]/scale;
  normalized[2][1] = coefficient[6]/scale;
  normalized[1][2] = -coefficient[6]/scale;
  CheckSlice(coefficient[6] != 0.0L,"The first-slice fixture lost its small enriched coefficient.");
}

void CheckFirstSliceMoments(const std::vector<long double>& values, const std::vector<long double>& errors, long double scale, long double height) {
  const unsigned degree = 5;
  const unsigned count = 21;
  CheckSlice(values.size() == 2U*count && errors.size() == values.size(),"A first-slice moment vector has an invalid size.");
  if(values.size() != 2U*count || errors.size() != values.size()) return;
  const long double thickness = (1.0L/256.0L-height)/(1.0L-height);
  CheckSlice(values[count] > 0.0L,"A first-slice interface returned a false zero delta mass.");
  for(unsigned i = 0; i <= degree; ++i) {
    for(unsigned j = 0; j <= degree-i; ++j) {
      const unsigned index = LineSliceUtilities::MomentIndex(degree,i,j);
      const long double expectedH = SlicePower(thickness,i+j+2U)/((j+1U)*(i+j+2U));
      const long double expectedD = SlicePower(thickness,i+j+1U)/(2.0L*(1.0L-height)*(j+1U));
      const long double toleranceH = 1.0e-13L*(1.0L+std::abs(expectedH));
      const long double toleranceD = 1.0e-13L*(1.0L+std::abs(expectedD));
      CheckSliceValue(values[index],expectedH,toleranceH);
      CheckSliceValue(values[count+index]/scale,expectedD,toleranceD);
      CheckSlice(std::isfinite(errors[index]) && errors[index] >= 0.0L && errors[index] <= toleranceH,"A first-slice volume error estimate exceeds the original final tolerance.");
      CheckSlice(std::isfinite(errors[count+index]) && errors[count+index] >= 0.0L && errors[count+index]/scale <= toleranceD,"A first-slice delta error estimate exceeds the original final tolerance.");
    }
  }
}

void ProbeFirstSliceChart(const std::vector<std::vector<long double>>& normalized, long double scale, long double height, unsigned chart) {
  std::vector<std::vector<long double>> coefficient;
  FELevelSetMoments::LineSliceTriangleChart(normalized,chart,coefficient);
  long double lower = 0.0L;
  long double upper = 0.0L;
  FELevelSetMoments::LineSliceTriangleDerivativeBounds(coefficient,lower,upper);
  const long double amplification = chart == 0 ? 1.0L : SlicePower(chart == 1 ? 3.0L : 2.0L,5);
  const long double absoluteTolerance = chart == 0 ? static_cast<double>(0.002L*0.5L*1.0e-13L/SlicePower(3.0L,5)) : 1.0e-13L;
  const long double relativeTolerance = chart == 0 ? static_cast<double>(0.02L*1.0e-13L) : 1.0e-13L;
  const std::vector<unsigned> boundary(3,0);
  std::vector<long double> raw;
  std::vector<long double> rawError;
  unsigned evaluations = 0;
  unsigned accepted = 0;
  bool completed = false;
  try {
    FELevelSetMoments::LineSliceDirection(coefficient,5,true,1,true,false,absoluteTolerance/amplification,absoluteTolerance*scale/amplification,relativeTolerance/amplification,32,200000,raw,rawError,evaluations,accepted,lower > 0.0L ? 1U : 0U,&boundary);
    completed = true;
  }
  catch(const std::exception&) {
    if(chart == 0) ++sliceRejectedCharts;
    else CheckSlice(false,"A regular first-slice chart failed at the original final tolerance.");
  }
  if(!completed) return;
  std::vector<long double> values;
  std::vector<long double> errors;
  FELevelSetMoments::LineSliceTriangleRestoreMoments(5,chart,raw,rawError,values,errors);
  CheckFirstSliceMoments(values,errors,scale,height);
}

void ProbeFirstSliceSelection(const std::vector<std::vector<long double>>& normalized, long double scale, long double height, unsigned preferred) {
  const std::vector<unsigned> boundary(3,0);
  std::vector<long double> values;
  std::vector<long double> errors;
  unsigned evaluations = 0;
  unsigned accepted = 0;
  unsigned selected = 0;
  bool completed = false;
  try {
    FELevelSetMoments::LineSliceTriangleDirections(normalized,scale,5,true,false,1.0e-13L,1.0e-13L,32,preferred,200000,&boundary,values,errors,evaluations,accepted,selected);
    completed = true;
  }
  catch(const std::exception&) { CheckSlice(false,"First-slice chart selection failed at the original final tolerance."); }
  if(completed) CheckFirstSliceMoments(values,errors,scale,height);
}

void CheckTetOutputRounding() {
  femus::TetBiquadratic basis;
  std::vector<double> nodal(basis._nc);
  for(unsigned node = 0; node < nodal.size(); ++node) nodal[node] = basis.GetXcoarse(node)[0]+1.0;
  std::vector<std::vector<std::vector<double>>> positive(1,std::vector<std::vector<double>>(1,std::vector<double>(1,17.0)));
  std::vector<std::vector<std::vector<double>>> delta(1,std::vector<std::vector<double>>(1,std::vector<double>(1,19.0)));
  std::vector<double> estimated;
  bool rejected = false;
  try { FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,0,positive,delta,8.0e-18,0.0,true,24,0,200000,&estimated,0); }
  catch(const std::exception&) { rejected = true; }
  CheckSlice(rejected,"An unattainable final double-output tolerance was accepted.");
  CheckSlice(positive == std::vector<std::vector<std::vector<double>>>(1,std::vector<std::vector<double>>(1,std::vector<double>(1,17.0))),"Rejected TET output rounding changed positive moments.");
  CheckSlice(delta == std::vector<std::vector<std::vector<double>>>(1,std::vector<std::vector<double>>(1,std::vector<double>(1,19.0))),"Rejected TET output rounding changed delta moments.");
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

int main() {
  std::cout << std::setprecision(18);
  CheckCorrectedIntervals();
  std::vector<std::vector<long double>> normalized;
  long double scale = 0.0L;
  long double height = 0.0L;
  BuildFirstSlice(normalized,scale,height);
  for(unsigned chart = 0; chart < 3; ++chart) ProbeFirstSliceChart(normalized,scale,height,chart);
  for(unsigned preferred = 0; preferred <= 3; ++preferred) ProbeFirstSliceSelection(normalized,scale,height,preferred);
  CheckTetOutputRounding();
  std::cout << "First-slice regressions: checks=" << sliceChecks << " failures=" << sliceFailures << " explicit_chart_rejections=" << sliceRejectedCharts << " maximum_absolute_error=" << sliceMaximumError << '\n';
  return sliceFailures == 0 ? 0 : 1;
}
