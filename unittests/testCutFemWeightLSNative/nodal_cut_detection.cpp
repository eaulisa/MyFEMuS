#include "NodalCutDetection.hpp"
#include "Basis.hpp"
#include <iostream>
#include <memory>

unsigned detectionChecks = 0;

void CheckDetection(bool condition) {
  ++detectionChecks;
  if(!condition) throw std::runtime_error("Nodal cut-detection check "+std::to_string(detectionChecks)+" failed.");
}

std::unique_ptr<femus::basis> DetectionBasis(GeomElType geometry) {
  if(geometry == LINE) return std::unique_ptr<femus::basis>(new femus::LineBiquadratic);
  if(geometry == QUAD) return std::unique_ptr<femus::basis>(new femus::QuadBiquadratic);
  if(geometry == TRI) return std::unique_ptr<femus::basis>(new femus::TriBiquadratic);
  if(geometry == HEX) return std::unique_ptr<femus::basis>(new femus::HexBiquadratic);
  if(geometry == TET) return std::unique_ptr<femus::basis>(new femus::TetBiquadratic);
  return std::unique_ptr<femus::basis>(new femus::WedgeBiquadratic);
}

void DetectionExamples() {
  std::vector<long double> controls;
  long double padding;
  const std::vector<double> zeroControl = {1,1,0.5};
  FELevelSetMoments::NodalLevelSetBernsteinControls(LINE,zeroControl,controls,padding);
  CheckDetection(*std::min_element(controls.begin(),controls.end()) == 0);
  CheckDetection(FELevelSetMoments::NodalLevelSetPhase(LINE,zeroControl) == 1);
  const std::vector<double> mixedControls = {1,1,0.45};
  FELevelSetMoments::NodalLevelSetBernsteinControls(LINE,mixedControls,controls,padding);
  CheckDetection(*std::min_element(controls.begin(),controls.end()) < 0);
  CheckDetection(FELevelSetMoments::NodalLevelSetPhase(LINE,mixedControls) == 1);
  const std::vector<double> crossing = {-0.25,0.75,0.25};
  FELevelSetMoments::NodalLevelSetBernsteinControls(LINE,crossing,controls,padding);
  CheckDetection(std::find(controls.begin(),controls.end(),0) == controls.end());
  CheckDetection(FELevelSetMoments::NodalLevelSetPhase(LINE,crossing) == 0);
}

void DetectionElement(GeomElType geometry) {
  std::unique_ptr<femus::basis> basis = DetectionBasis(geometry);
  const unsigned dimension = geometry == LINE ? 1U : (geometry == QUAD || geometry == TRI ? 2U : 3U);
  const double center[3] = {0.173,0.219,0.137};
  const double radii[3] = {1e-2,1e-4,1e-6};
  std::vector<double> nodal(basis->_nc);
  for(unsigned test = 0; test < 3; ++test) {
    for(unsigned n = 0; n < nodal.size(); ++n) {
      const double* point = basis->GetXcoarse(n);
      nodal[n] = -radii[test]*radii[test];
      for(unsigned d = 0; d < dimension; ++d) nodal[n] += (point[d]-center[d])*(point[d]-center[d]);
      CheckDetection(nodal[n] > 0);
    }
    CheckDetection(FELevelSetMoments::NodalLevelSetPhase(geometry,nodal) == 0);
    double value = 0;
    for(unsigned n = 0; n < nodal.size(); ++n) value += nodal[n]*basis->eval_phi(basis->GetIND(n),center);
    CheckDetection(value < 0);
    for(unsigned n = 0; n < nodal.size(); ++n) nodal[n] *= -1e-120;
    CheckDetection(FELevelSetMoments::NodalLevelSetPhase(geometry,nodal) == 0);
  }
  for(unsigned n = 0; n < nodal.size(); ++n) {
    const double* point = basis->GetXcoarse(n);
    nodal[n] = 0.02;
    for(unsigned d = 0; d < dimension; ++d) nodal[n] += (point[d]-center[d])*(point[d]-center[d]);
  }
  CheckDetection(FELevelSetMoments::NodalLevelSetPhase(geometry,nodal) == 1);
  for(unsigned n = 0; n < nodal.size(); ++n) nodal[n] *= -1;
  CheckDetection(FELevelSetMoments::NodalLevelSetPhase(geometry,nodal) == -1);
  nodal.assign(nodal.size(),1e-200);
  CheckDetection(FELevelSetMoments::NodalLevelSetPhase(geometry,nodal) == 1);
  nodal.assign(nodal.size(),-1e200);
  CheckDetection(FELevelSetMoments::NodalLevelSetPhase(geometry,nodal) == -1);
  nodal.assign(nodal.size(),0);
  CheckDetection(FELevelSetMoments::NodalLevelSetPhase(geometry,nodal) == 0);
}

int main() {
  try {
    DetectionExamples();
    const GeomElType geometries[6] = {LINE,QUAD,TRI,HEX,TET,WEDGE};
    for(unsigned n = 0; n < 6; ++n) DetectionElement(geometries[n]);
    std::cout << "Nodal cut detection: " << detectionChecks << " checks passed.\n";
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
