#ifndef FEMUS_LINE_INTEGRATION_HPP
#define FEMUS_LINE_INTEGRATION_HPP

#include "AdaptiveLineSliceMoments.hpp"
#include "FELevelSetMoments.hpp"

namespace FELevelSetMoments {

class LineIntegration {
public:
  // Centralizes the one-dimensional integration kernels reused by higher-dimensional cut elements.
  static void IntegratePolynomial(long double a, long double b, long double c, unsigned degree, long double left, long double right, bool computeDelta, bool fullBoundaryDelta, std::vector<long double>& positive, std::vector<long double>& delta) {
    Integrate1DInternal(a,b,c,degree,left,right,computeDelta,fullBoundaryDelta,positive,delta);
  }

  static bool IntegrateFastGraph(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned order, bool triangle, bool computeDelta, long double derivativeBound, long double absoluteH, long double absoluteDelta, long double relativeTolerance, unsigned maximumDepth, std::vector<long double>& values, std::vector<long double>& errors, unsigned& accepted, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
    return IntegrateCanonicalGraph(coefficient,degree,order,triangle,computeDelta,derivativeBound,absoluteH,absoluteDelta,relativeTolerance,maximumDepth,values,errors,accepted,surfaceGradient);
  }

  static void IntegrateTaylor(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned triangle, const std::vector<long double>& domain, bool computePositive, bool computeDelta, bool fullBoundaryDelta, long double absolutePositive, long double absoluteDelta, long double relativeTolerance, unsigned seriesOrder, unsigned maximumDepth, unsigned integrationOrder, std::vector<long double>& values, std::vector<long double>& errors, unsigned& accepted, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
    IntegrateNormalizedCell(coefficient,degree,triangle,domain,computePositive,computeDelta,fullBoundaryDelta,absolutePositive,absoluteDelta,relativeTolerance,seriesOrder,maximumDepth,integrationOrder,values,errors,accepted,surfaceGradient);
  }

  static void IntegrateFullGaussKronrod(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool triangle, unsigned order, bool computeDelta, bool fullBoundaryDelta, long double absoluteH, long double absoluteDelta, long double relativeTolerance, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned guaranteedMonotoneOrder = 0, const std::vector<unsigned>* zeroBoundary = nullptr, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
    LineSliceDirection(coefficient,degree,triangle,order,computeDelta,fullBoundaryDelta,absoluteH,absoluteDelta,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,guaranteedMonotoneOrder,zeroBoundary,surfaceGradient);
  }
};

}

#endif
