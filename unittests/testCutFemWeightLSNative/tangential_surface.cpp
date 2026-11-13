#include "CutFemWeight_LS.hpp"
#include "TangentialSurfaceIntegration.hpp"
#include "Basis.hpp"
#include <chrono>
#include <iostream>
#include <memory>

unsigned checks = 0;
long double maximumError = 0;

void CheckTangential(long double actual, long double expected, long double tolerance) {
  ++checks;
  const long double error = std::abs(actual-expected)/(1+std::abs(expected));
  maximumError = std::max(maximumError,error);
  if(!std::isfinite(actual) || error > tolerance) { std::cerr << "Tangential check " << checks << ": " << actual << " expected " << expected << " scaled error " << error << '\n'; throw std::runtime_error("Tangential surface regression failed."); }
}

std::unique_ptr<femus::basis> TangentialBasis(GeomElType geometry) {
  if(geometry == LINE) return std::unique_ptr<femus::basis>(new femus::LineBiquadratic);
  if(geometry == QUAD) return std::unique_ptr<femus::basis>(new femus::QuadBiquadratic);
  if(geometry == TRI) return std::unique_ptr<femus::basis>(new femus::TriBiquadratic);
  if(geometry == HEX) return std::unique_ptr<femus::basis>(new femus::HexBiquadratic);
  if(geometry == TET) return std::unique_ptr<femus::basis>(new femus::TetBiquadratic);
  return std::unique_ptr<femus::basis>(new femus::WedgeBiquadratic);
}

unsigned TangentialDimension(GeomElType geometry) {
  return geometry == LINE ? 1U : (geometry == QUAD || geometry == TRI ? 2U : 3U);
}

long double TangentialIntegral(const CutFemWeight_LS<double,long double>& rule, const std::vector<double>& correction, unsigned i = 0, unsigned j = 0, unsigned k = 0) {
  long double result = 0;
  for(unsigned n = 0; n < correction.size(); ++n) {
    long double value = rule.GetGaussWeightPointer()[n]*correction[n]*FELevelSetMoments::ReferencePower(rule.GetGaussCoordinatePointer(0)[n],i);
    if(rule.GetDimension() > 1) value *= FELevelSetMoments::ReferencePower(rule.GetGaussCoordinatePointer(1)[n],j);
    if(rule.GetDimension() > 2) value *= FELevelSetMoments::ReferencePower(rule.GetGaussCoordinatePointer(2)[n],k);
    result += value;
  }
  return result;
}

void TangentialPlanes(GeomElType geometry, unsigned degree) {
  const unsigned dimension = TangentialDimension(geometry);
  std::unique_ptr<femus::basis> basis = TangentialBasis(geometry);
  CutFemWeight_LS<double,long double> rule(geometry,degree,"legendre");
  std::vector<double> nodal(basis->_nc), positive, correction, surface, saved;
  std::vector<std::vector<double>> coordinates(dimension,std::vector<double>(basis->_nc));
  for(unsigned n = 0; n < nodal.size(); ++n) {
    const double* point = basis->GetXcoarse(n);
    nodal[n] = 2*(point[0]-0.25);
    for(unsigned d = 0; d < dimension; ++d) coordinates[d][n] = (d+1)*point[d];
  }
  rule.SetDirectionValidation(1);
  rule.GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,correction,1e-10,1e-10);
  saved = correction;
  long double transverseScale = 1;
  for(unsigned d = 1; d < dimension; ++d) transverseScale *= d+1U;
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= (dimension > 1 ? degree-i : 0U); ++j) for(unsigned k = 0; k <= (dimension > 2 ? degree-i-j : 0U); ++k) {
    long double exact = FELevelSetMoments::ReferencePower(0.25L,i)/(4*transverseScale);
    if(geometry == QUAD || geometry == HEX) exact *= j%2U == 0 ? 2.0L/(j+1U) : 0;
    if(geometry == HEX || geometry == WEDGE) exact *= k%2U == 0 ? 2.0L/(k+1U) : 0;
    if(geometry == TRI || geometry == WEDGE) exact *= FELevelSetMoments::ReferencePower(0.75L,j+1U)/(j+1U);
    if(geometry == TET) exact *= FELevelSetMoments::ReferencePower(0.75L,j+k+2U)*FELevelSetMoments::ReferenceTriangleMoment(j,k);
    CheckTangential(TangentialIntegral(rule,correction,i,j,k),exact,2e-8L);
  }
  rule.GetPhysicalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,surface,1e-10,1e-10);
  CheckTangential(TangentialIntegral(rule,surface),TangentialIntegral(rule,saved)*4*transverseScale*transverseScale,2e-9L);
  rule.GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,correction,1e-10,1e-10);
  for(unsigned n = 0; n < saved.size(); ++n) CheckTangential(correction[n],saved[n],0);
  FELevelSetMoments::OuterQuadratureRules rules;
  rules.regular2D = rules.regular3D = 7;
  rules.tangent2D = rules.tangent3D = 31;
  rule.GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,correction,1e-10,1e-10,rules);
  CheckTangential(TangentialIntegral(rule,correction),TangentialIntegral(rule,saved),2e-9L);
  const double multipliers[3] = {0.125,-4,1024};
  for(unsigned test = 0; test < 3; ++test) {
    std::vector<double> scaled = nodal;
    for(unsigned n = 0; n < scaled.size(); ++n) scaled[n] *= multipliers[test];
    rule.GetTangentialSurfaceWeightsFromNodalGeometry(scaled,coordinates,positive,correction,1e-10,1e-10);
    CheckTangential(TangentialIntegral(rule,correction)*multipliers[test]*multipliers[test],TangentialIntegral(rule,saved),2e-9L);
  }
  for(unsigned n = 0; n < nodal.size(); ++n) nodal[n] = 1;
  rule.GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,correction,1e-10,1e-10);
  CheckTangential(TangentialIntegral(rule,correction),0,0);
}

void TangentialMap(const std::vector<double>& point, bool curved, std::vector<double>& physical, std::vector<std::vector<double>>& jacobian) {
  const unsigned dimension = point.size();
  const double x = point[0], y = point[1], z = dimension == 3 ? point[2] : 0;
  const double bend = curved ? 1 : 0;
  physical = {1.3*x+0.2*y+0.05*bend*y*z,1.1*y+0.1*z+0.04*bend*x*x};
  jacobian = {{1.3,0.2+0.05*bend*z},{0.08*bend*x,1.1}};
  if(dimension == 3) {
    physical.push_back(0.9*z+0.03*bend*x*y);
    jacobian[0].push_back(0.05*bend*y);
    jacobian[1].push_back(0.1);
    jacobian.push_back({0.03*bend*y,0.03*bend*x,0.9});
  }
}

double TangentialGraph(GeomElType geometry, const std::vector<double>& point, std::vector<double>& gradient) {
  const unsigned dimension = point.size();
  const double x = point[0], y = point[1];
  gradient.assign(dimension,0.);
  gradient[dimension-1U] = 1;
  if(geometry == TRI) { gradient[0] = 0.25-0.1*(1-2*x); return point[1]-(0.25*(1-x)+0.1*x*(1-x)); }
  if(geometry == QUAD) { gradient[0] = -0.3*x; return point[1]-(0.2+0.15*x*x); }
  if(geometry == TET) { gradient[0] = gradient[1] = 0.25; return point[2]-0.25*(1-x-y); }
  gradient[0] = -0.15*y;
  gradient[1] = -0.15*x;
  return point[2]-(0.2+0.15*x*y);
}

void TangentialShapes(const std::vector<double>& point, std::vector<double>& gradient) {
  const unsigned dimension = point.size();
  gradient.assign(4U*dimension,0.);
  gradient[0] = 1;
  gradient[dimension+1U] = 1;
  gradient[2U*dimension] = 2*point[0];
  gradient[3U*dimension] = point[1];
  gradient[3U*dimension+1U] = point[0];
  if(dimension == 3) { gradient[3U*dimension] *= point[2]; gradient[3U*dimension+1U] *= point[2]; gradient[3U*dimension+2U] = point[0]*point[1]; }
}

double TangentialPhysicalGradients(const std::vector<std::vector<double>>& jacobian, const std::vector<double>& gradient, const std::vector<double>& shapeGradient, std::vector<double>& projected) {
  const unsigned dimension = gradient.size();
  std::vector<std::vector<double>> augmented(dimension,std::vector<double>(2U*dimension,0.));
  for(unsigned i = 0; i < dimension; ++i) for(unsigned j = 0; j < dimension; ++j) { augmented[i][j] = jacobian[i][j]; augmented[i][dimension+j] = i == j ? 1 : 0; }
  double determinant = 1.;
  for(unsigned j = 0; j < dimension; ++j) {
    unsigned pivot = j;
    for(unsigned i = j+1U; i < dimension; ++i) if(std::abs(augmented[i][j]) > std::abs(augmented[pivot][j])) pivot = i;
    if(pivot != j) { augmented[pivot].swap(augmented[j]); determinant *= -1; }
    const double diagonal = augmented[j][j];
    determinant *= diagonal;
    for(unsigned k = 0; k < 2U*dimension; ++k) augmented[j][k] /= diagonal;
    for(unsigned i = 0; i < dimension; ++i) if(i != j) {
      const double multiplier = augmented[i][j];
      for(unsigned k = 0; k < 2U*dimension; ++k) augmented[i][k] -= multiplier*augmented[j][k];
    }
  }
  std::vector<double> physicalGradient(dimension,0.);
  for(unsigned a = 0; a < dimension; ++a) for(unsigned b = 0; b < dimension; ++b) physicalGradient[a] += augmented[b][dimension+a]*gradient[b];
  double squared = 0.;
  for(unsigned a = 0; a < dimension; ++a) squared += physicalGradient[a]*physicalGradient[a];
  projected.assign(shapeGradient.size(),0.);
  for(unsigned i = 0; i < shapeGradient.size()/dimension; ++i) {
    double dot = 0.;
    for(unsigned a = 0; a < dimension; ++a) { for(unsigned b = 0; b < dimension; ++b) projected[i*dimension+a] += augmented[b][dimension+a]*shapeGradient[i*dimension+b]; dot += projected[i*dimension+a]*physicalGradient[a]; }
    for(unsigned a = 0; a < dimension; ++a) projected[i*dimension+a] -= physicalGradient[a]*dot/squared;
  }
  return std::abs(determinant)*std::sqrt(squared);
}

void TangentialGraphCase(GeomElType geometry, bool curved) {
  const unsigned dimension = TangentialDimension(geometry);
  const unsigned tangentSize = dimension == 2 ? 1U : 3U;
  std::unique_ptr<femus::basis> basis = TangentialBasis(geometry);
  CutFemWeight_LS<double,long double> rule(geometry,8,"legendre");
  std::vector<double> nodal(basis->_nc), positive, correction;
  std::vector<std::vector<double>> coordinates(dimension,std::vector<double>(basis->_nc));
  for(unsigned n = 0; n < nodal.size(); ++n) {
    const double* node = basis->GetXcoarse(n);
    const std::vector<double> point(node,node+dimension);
    std::vector<double> physical, gradient;
    std::vector<std::vector<double>> jacobian;
    TangentialMap(point,curved,physical,jacobian);
    nodal[n] = TangentialGraph(geometry,point,gradient);
    for(unsigned d = 0; d < dimension; ++d) coordinates[d][n] = physical[d];
  }
  rule.SetDirectionValidation(0);
  rule.GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,correction,1e-10,1e-10);
  std::vector<long double> fitted(4U*dimension+17U,0), reference(fitted.size(),0), residualPlus(16,0), residualMinus(16,0);
  const double step = 1e-6, cnew = 0.25, cold = 0.75, sigma = 0.7, dt = 0.03;
  for(unsigned n = 0; n < correction.size(); ++n) {
    std::vector<double> point(dimension), physical, gradient, shapeGradient, force, tangent;
    std::vector<std::vector<double>> jacobian;
    for(unsigned d = 0; d < dimension; ++d) point[d] = rule.GetGaussCoordinatePointer(d)[n];
    TangentialMap(point,curved,physical,jacobian);
    TangentialGraph(geometry,point,gradient);
    TangentialShapes(point,shapeGradient);
    const double squared = FELevelSetMoments::TangentialSurfaceFactors(jacobian,gradient,shapeGradient,force,tangent);
    const long double weight = rule.GetGaussWeightPointer()[n]*correction[n];
    for(unsigned k = 0; k < force.size(); ++k) fitted[k] += weight*force[k];
    for(unsigned i = 0; i < 4; ++i) for(unsigned j = 0; j < 4; ++j) for(unsigned a = 0; a < tangentSize; ++a) fitted[4U*dimension+4U*i+j] += weight*tangent[i*tangentSize+a]*tangent[j*tangentSize+a];
    fitted.back() += weight*squared;
    for(unsigned j = 0; j < 4; ++j) for(unsigned a = 0; a < tangentSize; ++a) {
      double plus = 0., minus = 0.;
      for(unsigned k = 0; k < 4; ++k) {
        const double current = 0.2+0.1*k, previous = -0.3+0.05*k;
        plus += (cnew*(current+(k == j ? step : 0))-cold*previous)*tangent[k*tangentSize+a];
        minus += (cnew*(current-(k == j ? step : 0))-cold*previous)*tangent[k*tangentSize+a];
      }
      for(unsigned i = 0; i < 4; ++i) { residualPlus[4U*i+j] -= cnew*sigma*dt*weight*plus*tangent[i*tangentSize+a]; residualMinus[4U*i+j] -= cnew*sigma*dt*weight*minus*tangent[i*tangentSize+a]; }
    }
  }
  femus::Gauss line(LINE,36,"legendre");
  const bool simplex = geometry == TRI || geometry == TET || geometry == WEDGE;
  for(unsigned ix = 0; ix < line.GetGaussPointsNumber(); ++ix) for(unsigned iy = 0; iy < (dimension == 3 ? line.GetGaussPointsNumber() : 1U); ++iy) {
    std::vector<double> point(dimension,0), physical, gradient, shapeGradient;
    std::vector<std::vector<double>> jacobian;
    point[0] = simplex ? 0.5*(line.GetGaussCoordinatePointer(0)[ix]+1) : line.GetGaussCoordinatePointer(0)[ix];
    long double weight = line.GetGaussWeightsPointer()[ix]*(simplex ? 0.5 : 1);
    if(dimension == 3) { point[1] = simplex ? 0.5*(1-point[0])*(line.GetGaussCoordinatePointer(0)[iy]+1) : line.GetGaussCoordinatePointer(0)[iy]; weight *= line.GetGaussWeightsPointer()[iy]*(simplex ? 0.5*(1-point[0]) : 1); }
    point[dimension-1U] = -TangentialGraph(geometry,point,gradient);
    TangentialMap(point,curved,physical,jacobian);
    TangentialShapes(point,shapeGradient);
    std::vector<double> projected;
    const double area = TangentialPhysicalGradients(jacobian,gradient,shapeGradient,projected);
    for(unsigned i = 0; i < 4; ++i) for(unsigned a = 0; a < dimension; ++a) reference[i*dimension+a] += weight*area*projected[i*dimension+a];
    for(unsigned i = 0; i < 4; ++i) for(unsigned j = 0; j < 4; ++j) for(unsigned a = 0; a < dimension; ++a) reference[4U*dimension+4U*i+j] += weight*area*projected[i*dimension+a]*projected[j*dimension+a];
    reference.back() += weight*area;
  }
  for(unsigned n = 0; n < fitted.size(); ++n) CheckTangential(fitted[n],reference[n],2e-8L);
  for(unsigned i = 0; i < 4; ++i) for(unsigned j = 0; j < 4; ++j) CheckTangential(-(residualPlus[4U*i+j]-residualMinus[4U*i+j])/(2*step),cnew*cnew*sigma*dt*fitted[4U*dimension+4U*i+j],1e-10L);
}

void TangentialClosed(GeomElType geometry) {
  const unsigned dimension = TangentialDimension(geometry);
  std::unique_ptr<femus::basis> basis = TangentialBasis(geometry);
  CutFemWeight_LS<double,long double> rule(geometry,3,"legendre");
  std::vector<double> nodal(basis->_nc), positive, correction;
  std::vector<std::vector<double>> coordinates(dimension,std::vector<double>(basis->_nc));
  const double radius = 0.125, center = 0.25;
  for(unsigned n = 0; n < nodal.size(); ++n) {
    nodal[n] = radius*radius;
    const double* point = basis->GetXcoarse(n);
    for(unsigned d = 0; d < dimension; ++d) { nodal[n] -= (point[d]-center)*(point[d]-center); coordinates[d][n] = point[d]; }
  }
  rule.SetDirectionValidation(0);
  rule.GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,correction,1e-9,1e-9);
  const long double mass = dimension == 2 ? std::acos(-1.0L)/(2*radius) : std::acos(-1.0L);
  CheckTangential(TangentialIntegral(rule,correction),mass,2e-8L);
  CheckTangential(TangentialIntegral(rule,correction,1),center*mass,2e-8L);
  CheckTangential(TangentialIntegral(rule,correction,2),(center*center+radius*radius/dimension)*mass,2e-8L);
}

int main() {
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  const GeomElType geometries[6] = {LINE,QUAD,TRI,HEX,TET,WEDGE};
  for(unsigned n = 0; n < 6; ++n) {
    for(unsigned degree = 2; degree <= 10; degree = degree == 2 ? 5U : (degree == 5 ? 10U : 11U)) TangentialPlanes(geometries[n],degree);
    if(n != 0) { TangentialGraphCase(geometries[n],false); TangentialGraphCase(geometries[n],true); TangentialClosed(geometries[n]); }
    std::cout << "Tangential element " << geometries[n] << " passed" << std::endl;
  }
  std::cout << "Tangential surface: " << checks << " checks; maximum scaled error=" << maximumError << "; seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count() << '\n';
}
