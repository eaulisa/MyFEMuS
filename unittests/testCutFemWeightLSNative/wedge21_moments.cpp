#include "CutFemWeight_LS.hpp"
#include "Basis.hpp"
#include <chrono>
#include <iostream>
#include <memory>

unsigned checks = 0;
long double maximumError = 0;

void CheckWedge(long double actual, long double exact, long double tolerance) {
  const long double error = std::abs(actual-exact);
  maximumError = std::max(maximumError,error);
  ++checks;
  if(!std::isfinite(actual) || error > tolerance) { std::cerr << "WEDGE21 check " << checks << ": " << actual << " expected " << exact << " error " << error << '\n'; throw std::runtime_error("WEDGE21 regression failed."); }
}

double WedgeField(unsigned field, double x, double y, double t) {
  if(field == 0) return 1;
  if(field == 1) return -1;
  if(field == 2) return t-0.375;
  if(field == 3) return x-0.25;
  if(field == 4) return y-0.25;
  if(field == 5) return x+y-0.5;
  if(field == 6) return t-(x+y)*(x+y);
  if(field == 7) return t-0.25-0.5*x*y*(1-x-y);
  if(field == 8) return 0.015625-(x-0.25)*(x-0.25)-(y-0.25)*(y-0.25)-(t-0.5)*(t-0.5);
  if(field == 9) return x;
  if(field == 10) return t;
  if(field == 11) return 1-x-y;
  if(field == 12) return 1-t;
  return 0.0625-(x-0.25)*(x-0.25)-(y-0.25)*(y-0.25);
}

std::vector<double> WedgeNodal(unsigned field) {
  femus::WedgeBiquadratic basis;
  std::vector<double> nodal(21);
  for(unsigned n = 0; n < 21; ++n) { const double* x = basis.GetXcoarse(n); nodal[n] = WedgeField(field,x[0],x[1],(x[2]+1)*0.5); }
  return nodal;
}

long double WedgeExact(unsigned field, unsigned i, unsigned j, unsigned k, bool delta) {
  const long double triangle = FELevelSetMoments::ReferenceTriangleMoment(i,j);
  if(field < 2) return delta || field == 1 ? 0 : triangle/(k+1U);
  if(field == 2) return triangle*(delta ? FELevelSetMoments::ReferencePower(0.375L,k) : (1-FELevelSetMoments::ReferencePower(0.375L,k+1U))/(k+1U));
  if(field == 3 || field == 4) {
    if(field == 4) std::swap(i,j);
    if(delta) return FELevelSetMoments::ReferencePower(0.25L,i)*FELevelSetMoments::ReferencePower(0.75L,j+1U)/((j+1U)*(k+1U));
    long double result = 0;
    for(unsigned n = 0; n <= j+1U; ++n) result += (n%2U == 0 ? 1 : -1)*FELevelSetMoments::ReferenceBinomial(j+1U,n)*(1-FELevelSetMoments::ReferencePower(0.25L,i+n+1U))/(i+n+1U);
    return result/((j+1U)*(k+1U));
  }
  if(field == 5) return triangle*(delta ? (i+j+2U)*FELevelSetMoments::ReferencePower(0.5L,i+j+1U) : 1-FELevelSetMoments::ReferencePower(0.5L,i+j+2U))/(k+1U);
  if(field == 6) return triangle*(delta ? static_cast<long double>(i+j+2U)/(i+j+2U*k+2U) : (1-static_cast<long double>(i+j+2U)/(i+j+2U*k+4U))/(k+1U));
  throw std::logic_error("No exact monomial formula.");
}

void WedgeCoefficients() {
  femus::WedgeBiquadratic basis;
  for(unsigned field = 0; field < 22; ++field) {
    std::vector<double> nodal(21,0);
    if(field < 21) nodal[field] = 1;
    else for(unsigned n = 0; n < 21; ++n) nodal[n] = std::sin(n+0.125);
    std::vector<long double> coefficient;
    FELevelSetMoments::Wedge21UnitPolynomialFromNodalValues(nodal,coefficient);
    for(unsigned sample = 0; sample < 20; ++sample) {
      std::vector<double> point = {((sample*13U+3U)%101U)/100.0,0,((sample*31U+7U)%101U)/50.0-1};
      point[1] = (1-point[0])*((sample*17U+1U)%101U)/100.0;
      const std::array<long double,3> unit = {point[0],point[1],(point[2]+1)*0.5};
      long double expected = 0;
      for(unsigned n = 0; n < 21; ++n) expected += nodal[n]*basis.eval_phi(basis.GetIND(n),point.data());
      CheckWedge(FELevelSetMoments::SurfaceEndpointIntegration::Value(coefficient,unit),expected,3e-14L);
      for(unsigned axis = 0; axis < 3; ++axis) {
        expected = 0;
        for(unsigned n = 0; n < 21; ++n) expected += nodal[n]*basis.eval_dphidxyz(axis,basis.GetIND(n),point.data());
        const std::vector<long double> derivative = FELevelSetMoments::SurfaceEndpointIntegration::Derivative(coefficient,axis);
        CheckWedge(FELevelSetMoments::SurfaceEndpointIntegration::Value(derivative,unit)*(axis == 2 ? 0.5L : 1),expected,8e-14L);
      }
    }
  }
}

void WedgeMoments() {
  const unsigned degree = 10;
  const unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  for(unsigned field = 0; field <= 6; ++field) {
    std::vector<long double> values, errors;
    FELevelSetMoments::IntegrateWedge21NodalMoments(WedgeNodal(field),degree,values,errors,1e-10L,1e-10L);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) for(unsigned pair = 0; pair < 2; ++pair) CheckWedge(values[pair*count+FELevelSetMoments::HexMomentIndex(i,j,k)],WedgeExact(field,i,j,k,pair != 0),3e-9L);
  }
  for(unsigned field = 2; field <= 6; ++field) for(unsigned direction = 1; direction <= 6; ++direction) {
    std::vector<long double> values, errors;
    FELevelSetMoments::IntegrateWedge21NodalMoments(WedgeNodal(field),2,values,errors,1e-9L,1e-9L,false,nullptr,direction);
    CheckWedge(values[0],field == 7 ? 0.375L-0.5L/120 : WedgeExact(field,0,0,0,false),3e-9L);
    CheckWedge(values[10],field == 7 ? 0.5L : WedgeExact(field,0,0,0,true),3e-9L);
  }
  std::vector<long double> enriched, errors;
  FELevelSetMoments::IntegrateWedge21NodalMoments(WedgeNodal(7),2,enriched,errors,1e-9L,1e-9L);
  CheckWedge(enriched[0],0.375L-0.5L/120,3e-9L);
  CheckWedge(enriched[10],0.5L,3e-9L);
}

void WedgeWeights() {
  CutFemWeight_LS<double,long double> rule(WEDGE,5,"legendre");
  for(unsigned field = 0; field <= 5; ++field) for(unsigned mode = 0; mode < 3; ++mode) {
    std::vector<double> positive, delta;
    const std::vector<double> nodal = WedgeNodal(field);
    femus::WedgeBiquadratic basis;
    std::vector<std::vector<double>> geometry(3,std::vector<double>(21));
    for(unsigned n = 0; n < 21; ++n) for(unsigned d = 0; d < 3; ++d) geometry[d][n] = basis.GetXcoarse(n)[d]*(d+2);
    if(mode == 2) rule.GetPhysicalSurfaceWeightsFromNodalGeometry(nodal,geometry,positive,delta,1e-10,1e-10);
    else rule.GetWeightsFromNodalValues(nodal,positive,delta,1e-10,1e-10,mode == 1);
    const long double metric = field == 2 ? (mode == 2 ? 3 : 0.5L) : (field == 3 ? (mode == 2 ? 12 : 1) : (field == 4 ? (mode == 2 ? 8 : 1) : (mode == 2 ? std::sqrt(208.0L) : std::sqrt(2.0L))));
    for(unsigned i = 0; i <= 5; ++i) for(unsigned j = 0; j <= 5-i; ++j) for(unsigned k = 0; k <= 5-i-j; ++k) {
      long double h = 0, d = 0;
      for(unsigned n = 0; n < positive.size(); ++n) {
        const long double monomial = FELevelSetMoments::ReferencePower(rule.GetGaussCoordinatePointer(0)[n],i)*FELevelSetMoments::ReferencePower(rule.GetGaussCoordinatePointer(1)[n],j)*FELevelSetMoments::ReferencePower((rule.GetGaussCoordinatePointer(2)[n]+1)*0.5L,k)*rule.GetGaussWeightPointer()[n];
        h += monomial*positive[n]; d += monomial*delta[n];
      }
      CheckWedge(h,2*WedgeExact(field,i,j,k,false),2e-9L);
      CheckWedge(d,2*WedgeExact(field,i,j,k,true)*(mode == 0 ? 1 : metric),2e-9L);
    }
  }
}

void WedgeSurfaces() {
  const long double pi = std::acos(-1.0L);
  for(unsigned field = 8; field <= 13; ++field) {
    std::vector<long double> values, errors;
    FELevelSetMoments::IntegrateWedge21NodalMoments(WedgeNodal(field),0,values,errors,1e-8L,1e-8L);
    const long double expectedH = field == 8 ? pi/384 : (field == 13 ? pi/16 : 0.5L);
    const long double expectedD = field == 8 ? pi/4 : (field == 13 ? pi : (field == 9 || field == 11 ? 0.5L : 0.25L));
    CheckWedge(values[0],expectedH,4e-8L);
    CheckWedge(values[1],expectedD,4e-8L);
  }
  std::vector<long double> values, errors;
  FELevelSetMoments::IntegrateWedge21NodalMoments(WedgeNodal(2),0,values,errors,1e-10L,1e-10L,false,nullptr,0,2);
  CheckWedge(values[1],0.5L,1e-9L);
}

void WedgeGeometryAndScaling() {
  femus::WedgeBiquadratic basis;
  CutFemWeight_LS<double,long double> rule(WEDGE,3,"legendre");
  std::vector<std::vector<double>> geometry(3,std::vector<double>(21));
  std::vector<double> nodal = WedgeNodal(3);
  for(unsigned node = 0; node < 21; ++node) {
    const double* point = basis.GetXcoarse(node);
    const double t = (point[2]+1)*0.5;
    geometry[0][node] = point[0]+t*t;
    geometry[1][node] = point[1];
    geometry[2][node] = point[2];
  }
  std::vector<double> positive, surface;
  rule.GetPhysicalSurfaceWeightsFromNodalGeometry(nodal,geometry,positive,surface,1e-10,1e-10);
  long double area = 0;
  for(unsigned n = 0; n < surface.size(); ++n) area += rule.GetGaussWeightPointer()[n]*surface[n];
  CheckWedge(area,0.75L*(std::sqrt(2.0L)+std::asinh(1.0L)),3e-9L);
  nodal = WedgeNodal(2);
  for(unsigned test = 0; test < 3; ++test) {
    const double multiplier = test == 0 ? 1e-12 : (test == 1 ? 1e12 : -2);
    std::vector<double> scaled = nodal;
    for(unsigned n = 0; n < scaled.size(); ++n) scaled[n] *= multiplier;
    std::vector<long double> values, errors;
    FELevelSetMoments::IntegrateWedge21NodalMoments(scaled,2,values,errors,1e-9L,1e-9L);
    CheckWedge(values[0],multiplier > 0 ? 0.3125L : 0.1875L,3e-9L);
    CheckWedge(values[10]*std::abs(multiplier),0.5L,3e-9L);
    FELevelSetMoments::IntegrateWedge21NodalMoments(scaled,2,values,errors,1e-9L,1e-9L,true);
    CheckWedge(values[10],0.25L,3e-9L);
  }
  for(unsigned n = 0; n < nodal.size(); ++n) nodal[n] = (basis.GetXcoarse(n)[2]+1)*0.5-(1-1e-8);
  std::vector<long double> values, errors;
  FELevelSetMoments::IntegrateWedge21NodalMoments(nodal,2,values,errors,1e-12L,1e-10L);
  CheckWedge(values[0],0.5e-8L,2e-15L);
  CheckWedge(values[10],0.5L,3e-9L);
  FELevelSetMoments::OuterQuadratureRules rules;
  rules.regular2D = 7; rules.tangent2D = 31; rules.regular3D = 7; rules.tangent3D = 31;
  FELevelSetMoments::IntegrateWedge21NodalMoments(WedgeNodal(6),2,values,errors,1e-9L,1e-9L,false,nullptr,0,1,rules);
  CheckWedge(values[0],0.25L,3e-9L);
  CheckWedge(values[10],0.5L,3e-9L);
}

int main() {
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  WedgeCoefficients(); WedgeMoments(); WedgeWeights(); WedgeSurfaces(); WedgeGeometryAndScaling();
  std::cout << "WEDGE21: " << checks << " checks; maximum absolute error=" << maximumError << "; seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count() << '\n';
}
