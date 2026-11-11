#ifndef FEMUS_TRIANGLE_INTEGRATION_HPP
#define FEMUS_TRIANGLE_INTEGRATION_HPP

#include "PlanarEndpointIntegration.hpp"

namespace FELevelSetMoments {

class TriangleIntegration {
public:
  // Integrates the full standard triangle and its exact parent-coordinate moments.
  static void Integrate(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned& selected, const std::vector<std::vector<std::vector<long double>>>* gradient = nullptr, const std::array<long double,5>* coordinates = nullptr, const std::vector<unsigned char>* requested = nullptr, unsigned direction = 0, bool fullBoundaryDelta = false, const std::vector<unsigned>* zeroBoundary = nullptr, const OuterQuadratureRules& rules = OuterQuadratureRules()) {
    PlanarEndpointIntegration::Integrate(coefficient,degree,true,computeDelta,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,selected,gradient,coordinates,requested,direction,fullBoundaryDelta,zeroBoundary,rules);
  }
};

}

#endif
