#ifndef FEMUS_LINE_INTEGRATION_HPP
#define FEMUS_LINE_INTEGRATION_HPP

#include "LineSliceUtilities.hpp"
#include "FELevelSetMoments.hpp"
#include "LineDomainIntegration.hpp"

namespace FELevelSetMoments {

class LineIntegration : public LineDomainIntegration {
public:
  // Centralizes the one-dimensional integration kernels reused by higher-dimensional cut elements.
  static void IntegratePolynomial(long double a, long double b, long double c, unsigned degree, long double left, long double right, bool computeDelta, bool fullBoundaryDelta, std::vector<long double>& positive, std::vector<long double>& delta) {
    Integrate1DInternal(a,b,c,degree,left,right,computeDelta,fullBoundaryDelta,positive,delta);
  }

  static bool IntegrateFastGraph(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned order, bool triangle, bool computeDelta, long double derivativeBound, long double absoluteH, long double absoluteDelta, long double relativeTolerance, unsigned maximumDepth, std::vector<long double>& values, std::vector<long double>& errors, unsigned& accepted, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
    return IntegrateCanonicalGraph(coefficient,degree,order,triangle,computeDelta,derivativeBound,absoluteH,absoluteDelta,relativeTolerance,maximumDepth,values,errors,accepted,surfaceGradient);
  }




};

}

#endif
