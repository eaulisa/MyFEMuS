#include "CutFemWeight_LS.hpp"
#include "TangentialSurfaceIntegration.hpp"
#include "Basis.hpp"
#include <chrono>
#include <iostream>
#include <memory>

unsigned normalChecks = 0;
long double normalMaximumError = 0;

void CheckNormal(long double actual, long double expected, long double tolerance) {
  ++normalChecks;
  const long double error = std::abs(actual-expected)/(1+std::abs(expected));
  normalMaximumError = std::max(normalMaximumError,error);
  if(!std::isfinite(actual) || error > tolerance) throw std::runtime_error("Normal surface check "+std::to_string(normalChecks)+": scaled error "+std::to_string(static_cast<double>(error)));
}

long double NormalIntegral(const CutFemWeight_LS<double,long double>& rule, const std::vector<double>& weight, unsigned i, unsigned j, unsigned k) {
  long double value = 0;
  for(unsigned n = 0; n < weight.size(); ++n) value += rule.GetGaussWeightPointer()[n]*weight[n]*FELevelSetMoments::ReferencePower(rule.GetGaussCoordinatePointer(0)[n],i)*FELevelSetMoments::ReferencePower(rule.GetGaussCoordinatePointer(1)[n],j)*FELevelSetMoments::ReferencePower(rule.GetGaussCoordinatePointer(2)[n],k);
  return value;
}

void NormalPlane(GeomElType geometry) {
  std::unique_ptr<femus::basis> basis;
  if(geometry == HEX) basis.reset(new femus::HexBiquadratic);
  else if(geometry == TET) basis.reset(new femus::TetBiquadratic);
  else basis.reset(new femus::WedgeBiquadratic);
  CutFemWeight_LS<double,long double> rule(geometry,4,"legendre");
  const double matrix[3][3] = {{1.1,0.2,0.05},{0.1,0.9,0.1},{-0.05,0.03,1.2}};
  const double normal[3] = {0.048,-0.043,0.97};
  const double norm = std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  std::vector<double> nodal(basis->_nc), positive, scalar;
  std::vector<std::vector<double>> coordinates(3,std::vector<double>(basis->_nc)), weights;
  for(unsigned n = 0; n < nodal.size(); ++n) {
    const double* point = basis->GetXcoarse(n);
    nodal[n] = 2*(point[2]-0.2);
    for(unsigned d = 0; d < 3; ++d) for(unsigned k = 0; k < 3; ++k) coordinates[d][n] += matrix[d][k]*point[k];
  }
  rule.GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,weights,1e-12,1e-12);
  for(unsigned i = 0; i <= 4; ++i) for(unsigned j = 0; j <= 4-i; ++j) for(unsigned k = 0; k <= 4-i-j; ++k) {
    long double area;
    if(geometry == HEX) area = FELevelSetMoments::ReferenceIntervalMoment(i)*FELevelSetMoments::ReferenceIntervalMoment(j);
    else area = FELevelSetMoments::ReferenceTriangleMoment(i,j)*FELevelSetMoments::ReferencePower(geometry == TET ? 0.8L : 1.0L,i+j+2U);
    for(unsigned d = 0; d < 3; ++d) CheckNormal(NormalIntegral(rule,weights[d],i,j,k),0.5L*area*FELevelSetMoments::ReferencePower(0.2L,k)*normal[d]/norm,2e-10L);
  }
  const std::vector<std::vector<double>> saved = weights;
  rule.GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,scalar,1e-12,1e-12);
  rule.GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,weights,1e-12,1e-12);
  for(unsigned d = 0; d < 3; ++d) for(unsigned n = 0; n < weights[d].size(); ++n) CheckNormal(weights[d][n],saved[d][n],0);
  const double multipliers[3] = {-4,0.125,1024};
  for(unsigned test = 0; test < 3; ++test) {
    std::vector<double> scaled = nodal;
    for(unsigned n = 0; n < scaled.size(); ++n) scaled[n] *= multipliers[test];
    rule.GetNormalSurfaceWeightsFromNodalGeometry(scaled,coordinates,positive,weights,1e-12,1e-12);
    for(unsigned d = 0; d < 3; ++d) CheckNormal(NormalIntegral(rule,weights[d],0,0,0)*multipliers[test],NormalIntegral(rule,saved[d],0,0,0),2e-10L);
  }
  if(geometry == TET) {
    rule.SetIntegrationDirection(1);
    rule.GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,weights,1e-11,1e-11);
    for(unsigned d = 0; d < 3; ++d) CheckNormal(NormalIntegral(rule,weights[d],0,0,0),NormalIntegral(rule,saved[d],0,0,0),2e-9L);
  }
}

long double NormalPolynomial(const std::vector<long double>& coefficient, const double* point, int derivative = -1) {
  long double value = 0;
  for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) {
    const unsigned exponent[3] = {i,j,k};
    if(derivative >= 0 && exponent[derivative] == 0) continue;
    long double term = coefficient[FELevelSetMoments::Hex27Index(i,j,k)]*(derivative < 0 ? 1U : exponent[derivative]);
    for(unsigned d = 0; d < 3; ++d) term *= FELevelSetMoments::ReferencePower(point[d],exponent[d]-(derivative == static_cast<int>(d)));
    value += term;
  }
  return value;
}

void NormalField(const std::vector<std::vector<long double>>& map, const std::vector<long double>& field, const double* point, std::vector<std::vector<double>>& jacobian, std::vector<double>& gradient, std::vector<double>& shapeGradient) {
  femus::HexBiquadratic basis;
  jacobian.assign(3,std::vector<double>(3));
  gradient.resize(3);
  shapeGradient.resize(81);
  for(unsigned a = 0; a < 3; ++a) {
    gradient[a] = NormalPolynomial(field,point,a);
    for(unsigned d = 0; d < 3; ++d) jacobian[d][a] = NormalPolynomial(map[d],point,a);
    for(unsigned i = 0; i < 27; ++i) shapeGradient[3U*i+a] = basis.eval_dphidxyz(a,basis.GetIND(i),point);
  }
}

void NormalCurvedHex() {
  femus::HexBiquadratic basis;
  std::vector<std::vector<long double>> map(3,std::vector<long double>(27,0));
  map[0][FELevelSetMoments::Hex27Index(1,0,0)] = map[1][FELevelSetMoments::Hex27Index(0,1,0)] = map[2][FELevelSetMoments::Hex27Index(0,0,1)] = 1;
  map[0][FELevelSetMoments::Hex27Index(2,2,2)] = 0.008L;
  map[0][FELevelSetMoments::Hex27Index(2,1,2)] = 0.004L;
  map[1][FELevelSetMoments::Hex27Index(2,2,2)] = 0.006L;
  map[1][FELevelSetMoments::Hex27Index(1,2,2)] = -0.003L;
  map[2][FELevelSetMoments::Hex27Index(2,2,2)] = -0.005L;
  map[2][FELevelSetMoments::Hex27Index(2,2,1)] = 0.002L;
  std::vector<long double> field(27,0);
  field[0] = -0.15L;
  field[FELevelSetMoments::Hex27Index(0,0,1)] = 1;
  field[FELevelSetMoments::Hex27Index(2,0,0)] = -0.07L;
  field[FELevelSetMoments::Hex27Index(1,1,0)] = -0.05L;
  field[FELevelSetMoments::Hex27Index(2,2,2)] = 0.015L;
  field[FELevelSetMoments::Hex27Index(2,1,2)] = 0.01L;
  field[FELevelSetMoments::Hex27Index(1,2,2)] = 0.008L;
  field[FELevelSetMoments::Hex27Index(2,2,1)] = 0.006L;
  std::vector<double> nodal(27), positive;
  std::vector<std::vector<double>> coordinates(3,std::vector<double>(27)), weights;
  for(unsigned n = 0; n < 27; ++n) {
    const double* point = basis.GetXcoarse(n);
    nodal[n] = NormalPolynomial(field,point);
    for(unsigned d = 0; d < 3; ++d) coordinates[d][n] = NormalPolynomial(map[d],point);
  }
  CutFemWeight_LS<double,long double> rule(HEX,13,"legendre");
  rule.SetDirectionValidation(0);
  rule.GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,weights,1e-13,1e-13);
  std::vector<long double> fitted(81,0), reference(81,0);
  for(unsigned n = 0; n < rule.GetGaussQuadraturePointNumber(); ++n) {
    const double point[3] = {rule.GetGaussCoordinatePointer(0)[n],rule.GetGaussCoordinatePointer(1)[n],rule.GetGaussCoordinatePointer(2)[n]};
    std::vector<std::vector<double>> jacobian;
    std::vector<double> gradient, shapes, force, tangent;
    NormalField(map,field,point,jacobian,gradient,shapes);
    FELevelSetMoments::TangentialSurfaceFactors(jacobian,gradient,shapes,force,tangent);
    for(unsigned i = 0; i < 27; ++i) for(unsigned d = 0; d < 3; ++d) fitted[3U*i+d] += rule.GetGaussWeightPointer()[n]*(weights[(d+1U)%3U][n]*tangent[3U*i+(d+2U)%3U]-weights[(d+2U)%3U][n]*tangent[3U*i+(d+1U)%3U]);
  }
  std::vector<long double> points, quadrature;
  FELevelSetMoments::LineIntegration::GaussRule(40,points,quadrature);
  for(unsigned first = 0; first < points.size(); ++first) for(unsigned second = 0; second < points.size(); ++second) {
    double point[3] = {static_cast<double>(2*points[first]-1),static_cast<double>(2*points[second]-1),0.2};
    for(unsigned iteration = 0; iteration < 12; ++iteration) point[2] -= NormalPolynomial(field,point)/NormalPolynomial(field,point,2);
    std::vector<std::vector<double>> jacobian;
    std::vector<double> gradient, shapes;
    NormalField(map,field,point,jacobian,gradient,shapes);
    double cofactor[3][3];
    double normal[3] = {0,0,0};
    for(unsigned a = 0; a < 3; ++a) for(unsigned b = 0; b < 3; ++b) { cofactor[a][b] = jacobian[(a+1U)%3U][(b+1U)%3U]*jacobian[(a+2U)%3U][(b+2U)%3U]-jacobian[(a+1U)%3U][(b+2U)%3U]*jacobian[(a+2U)%3U][(b+1U)%3U]; normal[a] += cofactor[a][b]*gradient[b]; }
    const double determinant = jacobian[0][0]*cofactor[0][0]+jacobian[0][1]*cofactor[0][1]+jacobian[0][2]*cofactor[0][2];
    const double norm = std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    for(unsigned a = 0; a < 3; ++a) normal[a] /= norm;
    for(unsigned i = 0; i < 27; ++i) {
      double physical[3] = {0,0,0};
      double dot = 0;
      for(unsigned a = 0; a < 3; ++a) { for(unsigned b = 0; b < 3; ++b) physical[a] += cofactor[a][b]*shapes[3U*i+b]/determinant; dot += physical[a]*normal[a]; }
      for(unsigned a = 0; a < 3; ++a) reference[3U*i+a] += 4*quadrature[first]*quadrature[second]*norm/std::abs(gradient[2])*(physical[a]-normal[a]*dot);
    }
  }
  for(unsigned n = 0; n < fitted.size(); ++n) CheckNormal(fitted[n],reference[n],2e-8L);
}

void NormalSmallSphere() {
  femus::HexBiquadratic basis;
  const double center[3] = {0.173,0.219,0.137};
  const double radius = 0.05;
  std::vector<double> nodal(27), positive;
  std::vector<std::vector<double>> coordinates(3,std::vector<double>(27)), weights;
  for(unsigned n = 0; n < 27; ++n) {
    const double* point = basis.GetXcoarse(n);
    nodal[n] = -radius*radius;
    for(unsigned d = 0; d < 3; ++d) { coordinates[d][n] = point[d]; nodal[n] += (point[d]-center[d])*(point[d]-center[d]); }
    CheckNormal(nodal[n] > 0,1,0);
  }
  CutFemWeight_LS<double,long double> rule(HEX,4,"legendre");
  CheckNormal(rule.GetLevelSetPhaseFromNodalValues(nodal),0,0);
  rule.SetDirectionValidation(0);
  rule.GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,weights,1e-11,1e-11);
  CheckNormal(NormalIntegral(rule,positive,0,0,0),8-4*std::acos(-1.0L)*radius*radius*radius/3,2e-10L);
  for(unsigned a = 0; a < 3; ++a) {
    CheckNormal(NormalIntegral(rule,weights[a],0,0,0),0,2e-10L);
    for(unsigned d = 0; d < 3; ++d) CheckNormal(NormalIntegral(rule,weights[a],d == 0,d == 1,d == 2),a == d ? 2*std::acos(-1.0L)*radius*radius/3 : 0,2e-10L);
  }
  nodal.assign(27,1);
  rule.GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,weights,1e-11,1e-11);
  CheckNormal(NormalIntegral(rule,positive,0,0,0),8,2e-10L);
  for(unsigned a = 0; a < 3; ++a) for(unsigned n = 0; n < weights[a].size(); ++n) CheckNormal(weights[a][n],0,0);
}

int main() {
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  try {
    NormalPlane(HEX);
    NormalPlane(TET);
    NormalPlane(WEDGE);
    NormalSmallSphere();
    NormalCurvedHex();
    std::cout << "Normal surface: " << normalChecks << " checks; maximum scaled error=" << normalMaximumError << "; seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count() << '\n';
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
