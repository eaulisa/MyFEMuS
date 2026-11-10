#include "SurfaceEndpointIntegration.hpp"
#include <iostream>
#include <iomanip>

std::vector<long double> NestedCoefficient(bool tet, unsigned multiplicity) {
  std::vector<long double> coefficient(27,0);
  coefficient[FELevelSetMoments::Hex27Index(1,0,0)] = 1;
  coefficient[FELevelSetMoments::Hex27Index(0,0,2)] = 1;
  coefficient[FELevelSetMoments::Hex27Index(0,2,1)] = -1;
  if(tet) {
    coefficient[FELevelSetMoments::Hex27Index(1,1,1)] = -1;
    coefficient[FELevelSetMoments::Hex27Index(0,1,2)] = -1;
    if(multiplicity == 3) { coefficient[FELevelSetMoments::Hex27Index(0,1,1)] = -1; coefficient[FELevelSetMoments::Hex27Index(0,2,0)] = 0.25L; }
  }
  else if(multiplicity == 3) { coefficient[FELevelSetMoments::Hex27Index(0,1,1)] = -2; coefficient[FELevelSetMoments::Hex27Index(0,2,0)] = 1; }
  return coefficient;
}

int main(int argc, char** argv) {
  const unsigned degree = argc > 1 ? std::strtoul(argv[1],nullptr,10) : 2;
  const long double tolerance = argc > 2 ? std::strtold(argv[2],nullptr) : 1e-9L;
  if(degree > 10 || !(tolerance > 0)) return 64;
  unsigned failures = 0;
  std::cout << std::setprecision(17);
  for(unsigned tet = 0; tet < 2; ++tet) for(unsigned multiplicity = 3; multiplicity <= 4; ++multiplicity) for(unsigned surface = 0; surface < 2; ++surface) {
    FELevelSetMoments::SurfaceEndpointData graph;
    graph.coefficient = NestedCoefficient(tet != 0,multiplicity);
    graph.unitCube = true;
    graph.degree = degree;
    graph.absoluteH = graph.absoluteD = graph.relative = 0.01L*tolerance;
    graph.maximumEvaluations = 2000000;
    if(surface != 0) for(unsigned axis = 0; axis < 3; ++axis) graph.gradient.push_back(FELevelSetMoments::SurfaceEndpointIntegration::Derivative(graph.coefficient,axis));
    std::vector<long double> reference;
    std::vector<long double> referenceErrors;
    std::vector<long double> values;
    std::vector<long double> errors;
    try {
      FELevelSetMoments::SurfaceEndpointIntegration::Integrate(graph,tet != 0,reference,referenceErrors);
      FELevelSetMoments::SurfaceEndpointData endpoint = graph;
      endpoint.evaluations = endpoint.accepted = 0;
      endpoint.absoluteH = endpoint.absoluteD = endpoint.relative = tolerance;
      endpoint.preferredOuter = 0;
      endpoint.preferredInner = 2;
      FELevelSetMoments::SurfaceEndpointIntegration::Integrate(endpoint,tet != 0,values,errors);
      long double volumeError = 0;
      long double surfaceError = 0;
      for(unsigned n = 0; n < values.size(); ++n) {
        if(n < values.size()/2U) volumeError = std::max(volumeError,std::abs(values[n]-reference[n]));
        else surfaceError = std::max(surfaceError,std::abs(values[n]-reference[n]));
      }
      if(volumeError > 3*tolerance || surfaceError > 3*tolerance) throw std::runtime_error("Nested contact projections disagree.");
      std::cout << "PASS," << tet << ',' << multiplicity << ',' << surface << ',' << volumeError << ',' << surfaceError << ',' << endpoint.evaluations << '\n';
    }
    catch(const std::exception& error) { ++failures; std::cout << "UNSUPPORTED," << tet << ',' << multiplicity << ',' << surface << ',' << error.what() << '\n'; }
  }
  return failures == 0 ? 0 : 1;
}
