#include "CutFemWeightRules_LS.hpp"
#include "Basis.hpp"
#include <iostream>

unsigned reuseChecks = 0;

void CheckSurfaceReuse(bool condition) {
  ++reuseChecks;
  if(!condition) throw std::runtime_error("Surface rule reuse check "+std::to_string(reuseChecks)+" failed.");
}

void SurfaceReuseWeights(CutFemWeight_LS<double,long double>& rule, const std::vector<double>& nodal, const std::vector<std::vector<double>>& coordinates, std::vector<std::vector<double>>& weights) {
  weights.resize(6);
  std::vector<std::vector<double>> normal;
  rule.GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,weights[0],weights[1],1e-10,1e-10);
  rule.GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,weights[0],normal,1e-10,1e-10);
  for(unsigned component = 0; component < 3; ++component) weights[component+2U] = normal[component];
  rule.GetPhysicalSurfaceWeightsFromNodalGeometry(nodal,coordinates,weights[0],weights[5],1e-10,1e-10);
}

void CheckSurfaceReuseWeights(const std::vector<std::vector<double>>& actual, const std::vector<std::vector<double>>& expected) {
  CheckSurfaceReuse(actual.size() == expected.size());
  for(unsigned mode = 0; mode < actual.size(); ++mode) {
    CheckSurfaceReuse(actual[mode].size() == expected[mode].size());
    for(unsigned point = 0; point < actual[mode].size(); ++point) CheckSurfaceReuse(std::abs(actual[mode][point]-expected[mode][point]) <= 2e-10*(1+std::abs(expected[mode][point])));
  }
}

void SurfaceReuseElement(CutFemWeightRules_LS<double,long double>& rules, GeomElType geometry) {
  std::unique_ptr<femus::basis> basis;
  if(geometry == HEX) basis.reset(new femus::HexBiquadratic);
  else if(geometry == TET) basis.reset(new femus::TetBiquadratic);
  else basis.reset(new femus::WedgeBiquadratic);
  CutFemWeight_LS<double,long double>* rule = rules.GetRule(geometry);
  const double* gauss = rule->GetGaussWeightPointer();
  const double* points[3] = {rule->GetGaussCoordinatePointer(0),rule->GetGaussCoordinatePointer(1),rule->GetGaussCoordinatePointer(2)};
  std::vector<double> nodal(basis->_nc);
  std::vector<std::vector<double>> coordinates(3,std::vector<double>(basis->_nc)), initial, repeated, changed;
  for(unsigned n = 0; n < nodal.size(); ++n) {
    const double* point = basis->GetXcoarse(n);
    nodal[n] = point[2]-0.2-0.02*point[0];
    coordinates[0][n] = 1.1*point[0]+0.03*point[1]*point[2];
    coordinates[1][n] = point[1]+0.02*point[0]*point[2];
    coordinates[2][n] = 0.9*point[2]+0.01*point[0]*point[1];
  }
  SurfaceReuseWeights(*rule,nodal,coordinates,initial);
  for(unsigned repeat = 0; repeat < 3; ++repeat) {
    FELevelSetMoments::SurfaceEndpointIntegration::Statistics() = FELevelSetMoments::SurfaceEndpointStatistics();
    SurfaceReuseWeights(*rules.GetRule(geometry),nodal,coordinates,repeated);
    CheckSurfaceReuseWeights(repeated,initial);
    CheckSurfaceReuse(FELevelSetMoments::SurfaceEndpointIntegration::Statistics().integrationPasses == 0);
    CheckSurfaceReuse(rule == rules.GetRule(geometry) && gauss == rule->GetGaussWeightPointer());
    for(unsigned axis = 0; axis < 3; ++axis) CheckSurfaceReuse(points[axis] == rule->GetGaussCoordinatePointer(axis));
  }
  rule->ClearCache();
  SurfaceReuseWeights(*rule,nodal,coordinates,repeated);
  CheckSurfaceReuseWeights(repeated,initial);
  for(unsigned n = 0; n < nodal.size(); ++n) nodal[n] -= 0.03;
  SurfaceReuseWeights(*rule,nodal,coordinates,changed);
  CheckSurfaceReuse(changed[0] != initial[0]);
  rule->ClearCache();
  SurfaceReuseWeights(*rule,nodal,coordinates,repeated);
  CheckSurfaceReuseWeights(repeated,changed);
  for(unsigned n = 0; n < coordinates[0].size(); ++n) coordinates[0][n] *= 1.2;
  SurfaceReuseWeights(*rule,nodal,coordinates,changed);
  rule->ClearCache();
  SurfaceReuseWeights(*rule,nodal,coordinates,repeated);
  CheckSurfaceReuseWeights(repeated,changed);
  CheckSurfaceReuse(rule->GetCutFEMQuadratureOrder() == rules.GetOrder());
  CheckSurfaceReuse(rule->GetOuterQuadratureRules().SameAs(rules.GetOuterRules()));
}

int main() {
  try {
    FELevelSetMoments::OuterQuadratureRules outerRules;
    outerRules.regular2D = outerRules.regular3D = 7;
    outerRules.tangent2D = outerRules.tangent3D = 31;
    CutFemWeightRules_LS<double,long double> rules(4,{HEX,TET,WEDGE,HEX},outerRules);
    SurfaceReuseElement(rules,HEX);
    SurfaceReuseElement(rules,TET);
    SurfaceReuseElement(rules,WEDGE);
    bool rejected = false;
    try { rules.GetRule(QUAD); }
    catch(const std::logic_error&) { rejected = true; }
    CheckSurfaceReuse(rejected);
    std::cout << "Surface rule reuse: " << reuseChecks << " checks passed.\n";
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
