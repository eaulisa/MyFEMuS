#include "CutFemWeight_LS.hpp"
#include "Basis.hpp"
#include <iostream>
#include <memory>
#include <sstream>

unsigned warningChecks = 0;

void WarningRequire(bool condition, const char* problem) {
  ++warningChecks;
  if(!condition) throw std::runtime_error(problem);
}

void WarningClose(long double value, long double expected, long double tolerance) {
  WarningRequire(std::isfinite(value) && std::abs(value-expected) <= tolerance*(1+std::abs(expected)),"Best available quadrature lost or duplicated an integration contribution.");
}

unsigned WarningLines(const std::string& text) {
  unsigned count = 0;
  for(unsigned n = 0; n < text.size(); ++n) count += text[n] == '\n';
  return count;
}

void WarningSample(std::vector<unsigned>& context, const FELevelSetMoments::OuterIntervalPoint& point, const std::vector<unsigned char>& active, std::vector<long double>& values, std::vector<long double>& errors, long double absoluteH, long double absoluteD, long double relative) {
  ++context[1];
  values.assign(active.size(),0);
  errors.assign(active.size(),0);
  for(unsigned n = 0; n < active.size(); ++n) if(active[n]) values[n] = point.jacobian*(context[0] == 0 ? (n == 0 ? std::exp(10*point.coordinate) : 1) : (n+1U)*FELevelSetMoments::ReferencePower(point.coordinate,4));
}

void WarningPanels() {
  const FELevelSetMoments::OuterInterval interval = {0,1,0,0,0,1,true};
  const std::vector<unsigned char> active(2,1);
  std::vector<unsigned> context = {0,0};
  std::vector<long double> values(2,0), errors(2,0);
  unsigned accepted = 0;
  FELevelSetMoments::OuterIntervalQuadrature::Adaptive(context,WarningSample,interval,0,1,1e-15L,1e-15L,0,0,0,1000,context[1],accepted,2,active,values,errors);
  WarningClose(values[0],(std::exp(10.0L)-1)/10,2e-12L);
  WarningClose(values[1],1,1e-15L);
  WarningRequire(errors[0] > 1e-15L,"The unmet panel error was hidden.");
  context = {1,0}; values.assign(2,0); errors.assign(2,0); accepted = 0;
  for(unsigned part = 0; part < 2; ++part) FELevelSetMoments::OuterIntervalQuadrature::Adaptive(context,WarningSample,interval,0.5L*part,0.5L*(part+1U),1e-30L,1e-30L,0,0,30,3,context[1],accepted,2,active,values,errors);
  WarningClose(values[0],0.2L,1e-15L);
  WarningClose(values[1],0.4L,1e-15L);
  WarningRequire(context[1] == 6,"An exhausted budget triggered further adaptive refinement.");
}

void WarningPlanar() {
  std::vector<std::vector<long double>> field(3,std::vector<long double>(3,0));
  field[0][1] = 1; field[2][0] = -1;
  std::vector<long double> values, errors;
  unsigned evaluations = 0, accepted = 0, selected = 0;
  FELevelSetMoments::PlanarEndpointIntegration::Integrate(field,1,false,true,1e-30L,1e-30L,0,0,3,values,errors,evaluations,accepted,selected);
  WarningClose(values[0],2.0L/3,1e-15L);
  WarningClose(values[3],1,1e-15L);
  WarningRequire(errors[0] > 1e-30L,"The failed planar tolerance was marked satisfied.");
  std::vector<unsigned char> mask(6,0); mask[0] = mask[3] = 1;
  evaluations = accepted = selected = 0;
  FELevelSetMoments::PlanarEndpointIntegration::Integrate(field,1,false,true,1e-30L,1e-30L,0,0,3,values,errors,evaluations,accepted,selected,nullptr,nullptr,&mask);
  WarningClose(values[0],2.0L/3,1e-15L);
  WarningClose(values[3],1,1e-15L);
  for(unsigned n = 0; n < values.size(); ++n) if(!mask[n]) WarningClose(values[n],0,0);
  field.assign(3,std::vector<long double>(3,0));
  field[0][0] = 0.015625L-0.5L; field[1][0] = field[0][1] = 1; field[2][0] = field[0][2] = -1;
  evaluations = accepted = selected = 0;
  FELevelSetMoments::PlanarEndpointIntegration::Integrate(field,1,false,true,1e-30L,1e-30L,0,0,3,values,errors,evaluations,accepted,selected);
  WarningClose(values[0],std::acos(-1.0L)*0.015625L,0.002L);
  WarningClose(values[3],std::acos(-1.0L),0.002L);
  field.assign(3,std::vector<long double>(3,0)); field[0][0] = -0.25L; field[0][1] = 1;
  evaluations = accepted = selected = 0;
  FELevelSetMoments::PlanarEndpointIntegration::Integrate(field,1,false,true,1e-10L,1e-10L,1e-10L,24,1000,values,errors,evaluations,accepted,selected);
  values[0] += 0.01L;
  const long double retained = values[0];
  FELevelSetMoments::PlanarEndpointIntegration::Compare(field,1,false,true,1e-12L,1e-12L,0,24,1000,values,errors,evaluations,accepted,selected,nullptr,nullptr);
  WarningClose(values[0],retained,0);
}

void WarningSurface(bool tetrahedron) {
  FELevelSetMoments::SurfaceEndpointData data;
  data.coefficient.assign(27,0);
  data.coefficient[FELevelSetMoments::Hex27Index(0,0,1)] = 1;
  if(tetrahedron) { data.coefficient[0] = -0.25L; data.coefficient[FELevelSetMoments::Hex27Index(1,0,0)] = data.coefficient[FELevelSetMoments::Hex27Index(0,1,0)] = 0.25L; }
  else data.coefficient[FELevelSetMoments::Hex27Index(2,0,0)] = -1;
  data.degree = 1;
  data.unitCube = !tetrahedron;
  data.maximumDepth = 0;
  data.maximumEvaluations = 3;
  data.absoluteH = data.absoluteD = 1e-30L;
  data.relative = 0;
  std::vector<long double> values, errors;
  FELevelSetMoments::SurfaceEndpointIntegration::Integrate(data,tetrahedron,values,errors);
  WarningClose(values[0],tetrahedron ? 0.125L : 2.0L/3,1e-14L);
  WarningClose(values[4],tetrahedron ? 0.5L : 1,1e-14L);
  WarningRequire(data.absoluteH == 1e-30L && data.absoluteD == 1e-30L && data.relative == 0 && data.requested.empty(),"A tolerance warning corrupted the caller's integration state.");
}

std::unique_ptr<femus::basis> WarningBasis(GeomElType geometry) {
  if(geometry == LINE) return std::unique_ptr<femus::basis>(new femus::LineBiquadratic);
  if(geometry == QUAD) return std::unique_ptr<femus::basis>(new femus::QuadBiquadratic);
  if(geometry == TRI) return std::unique_ptr<femus::basis>(new femus::TriBiquadratic);
  if(geometry == HEX) return std::unique_ptr<femus::basis>(new femus::HexBiquadratic);
  if(geometry == TET) return std::unique_ptr<femus::basis>(new femus::TetBiquadratic);
  return std::unique_ptr<femus::basis>(new femus::WedgeBiquadratic);
}

void WarningWeights(GeomElType geometry, std::ostringstream& diagnostics) {
  const unsigned dimension = geometry == LINE ? 1U : (geometry == QUAD || geometry == TRI ? 2U : 3U);
  std::unique_ptr<femus::basis> basis = WarningBasis(geometry);
  CutFemWeight_LS<double,long double> rule(geometry,2,"legendre");
  rule.SetDirectionValidation(0);
  std::vector<double> nodal(basis->_nc), positive, correction;
  std::vector<std::vector<double>> coordinates(dimension,std::vector<double>(basis->_nc));
  for(unsigned n = 0; n < nodal.size(); ++n) { const double* point = basis->GetXcoarse(n); nodal[n] = 2*(point[0]-0.25); for(unsigned d = 0; d < dimension; ++d) coordinates[d][n] = point[d]; }
  long double transverse = geometry == LINE ? 1 : (geometry == QUAD ? 2 : (geometry == TRI ? 0.75L : (geometry == HEX ? 4 : (geometry == TET ? 0.28125L : 1.5L))));
  unsigned previousLines = WarningLines(diagnostics.str());
  rule.GetWeightsFromNodalValues(nodal,positive,correction,1e-9,1e-9);
  WarningRequire(WarningLines(diagnostics.str()) == previousLines,"A converged quadrature emitted a warning.");
  for(unsigned mode = 0; mode < 3; ++mode) {
    if(mode == 0) rule.GetWeightsFromNodalValues(nodal,positive,correction,1e-30,0);
    if(mode == 1) rule.GetPhysicalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,correction,1e-30,0);
    if(mode == 2) rule.GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,correction,1e-30,0);
    long double mass = 0;
    for(unsigned n = 0; n < correction.size(); ++n) { WarningRequire(std::isfinite(correction[n]) && std::isfinite(positive[n]),"A warning returned nonfinite FEM weights."); mass += rule.GetGaussWeightPointer()[n]*correction[n]; }
    WarningClose(mass,transverse*(mode == 0 ? 0.5L : (mode == 1 ? 1 : 0.25L)),2e-11L);
    const unsigned lines = WarningLines(diagnostics.str());
    if(geometry != LINE || mode != 2) WarningRequire(lines == previousLines+1U,"Tolerance failures were not aggregated into one warning per FEM rule.");
    previousLines = lines;
  }
}

void WarningCanonical() {
  std::vector<std::vector<long double>> field(3,std::vector<long double>(3,0));
  field[0][1] = 1; field[2][0] = -1;
  std::vector<long double> values(6,0), errors(6,0);
  unsigned accepted = 0;
  FELevelSetMoments::IntegrateCanonicalGraphInterval(field,1,1,0,1,true,1e-30L,1e-30L,0,0,0,values,errors,accepted);
  WarningClose(values[0],2.0L/3,1e-15L);
  WarningClose(values[3],1,1e-15L);
}

int main() {
  std::ostringstream diagnostics;
  std::streambuf* previous = std::cerr.rdbuf(diagnostics.rdbuf());
  try {
    WarningPanels(); WarningPlanar(); WarningSurface(false); WarningSurface(true); WarningCanonical();
    const GeomElType geometries[6] = {LINE,QUAD,TRI,HEX,TET,WEDGE};
    for(unsigned n = 0; n < 6; ++n) WarningWeights(geometries[n],diagnostics);
    const std::string text = diagnostics.str();
    WarningRequire(text.find("maximum depth") != std::string::npos && text.find("evaluation budget") != std::string::npos,"A resource-limit warning did not describe the failure.");
    WarningRequire(text.find("estimated error=") != std::string::npos && text.find("requested tolerance=") != std::string::npos && text.find("returning the best available quadrature") != std::string::npos,"Tolerance diagnostics omitted the error, target, or recovery action.");
    WarningRequire(text.find("rounded") != std::string::npos && text.find("disagreement") != std::string::npos,"Rounding or independent validation failures were not reported.");
  }
  catch(const std::exception& error) { std::cerr.rdbuf(previous); std::cerr << error.what() << '\n' << diagnostics.str(); return 1; }
  std::cerr.rdbuf(previous);
  std::cout << "Tolerance warnings: " << warningChecks << " checks passed; forced depth, budget, precision, rounding, validation, and all native weight modes retained usable quadrature.\n";
}
