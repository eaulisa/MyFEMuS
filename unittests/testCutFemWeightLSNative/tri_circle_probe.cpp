#include "Tri7DirectionSelection.hpp"
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

int main() {
  unsigned failures = 0;
  const unsigned degrees[2] = {6U,10U};
  const long double radii[2] = {0.03125L,0.00390625L};
  const double tolerances[2] = {1.0e-10,1.0e-12};
  for(unsigned g = 0; g < 2; ++g) for(unsigned d = 0; d < 2; ++d) for(unsigned t = 0; t < 2; ++t) for(unsigned sign = 0; sign < 2; ++sign) {
    const long double cx = 0.25L;
    const long double cy = 0.25L;
    const long double radius = radii[g];
    const long double scale = sign == 0 ? 1.0L : -8.0L;
    std::vector<std::vector<double>> field(3,std::vector<double>(3,0.0));
    field[0][0] = static_cast<double>(scale*(cx*cx+cy*cy-radius*radius));
    field[1][0] = static_cast<double>(-2.0L*scale*cx);
    field[0][1] = static_cast<double>(-2.0L*scale*cy);
    field[2][0] = static_cast<double>(scale);
    field[0][2] = static_cast<double>(scale);
    std::vector<std::vector<double>> positive;
    std::vector<std::vector<double>> delta;
    std::vector<double> estimate;
#ifdef FEMUS_LINE_DOMAIN_INTEGRATION_HPP
    FELevelSetMoments::LineIntegration::Statistics() = FELevelSetMoments::LineIntegrationStatistics();
#endif
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    try {
      FELevelSetMoments::IntegrateTri7ReferenceSelected(field,degrees[d],positive,delta,tolerances[t],tolerances[t],true,false,28,0,200000,&estimate,1,nullptr,true);
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
      const long double area = std::acos(-1.0L)*radius*radius;
      const long double exactH = sign == 0 ? 0.5L-area : area;
      const long double exactD = std::acos(-1.0L)/std::abs(scale);
      if(!std::isfinite(positive[0][0]) || !std::isfinite(delta[0][0]) || std::abs(positive[0][0]-exactH) > tolerances[t]*(1.0L+std::abs(exactH)) || std::abs(delta[0][0]-exactD) > tolerances[t]*(1.0L+std::abs(exactD))) throw std::runtime_error("Analytic circle area or delta mass exceeds tolerance.");
      std::cout << "CASE g=" << g << " degree=" << degrees[d] << " tol=" << tolerances[t] << " sign=" << sign << " OK time=" << seconds;
      if(estimate.size() >= 7) std::cout << " accepted=" << estimate[2] << " evaluations=" << estimate[3] << " selected=" << estimate[4] << " compared=" << estimate[5] << " attempts=" << estimate[6];
#ifdef FEMUS_LINE_DOMAIN_INTEGRATION_HPP
      const FELevelSetMoments::LineIntegrationStatistics& statistics = FELevelSetMoments::LineIntegration::Statistics();
      std::cout << " fast=" << statistics.fastSamples << " full=" << statistics.fullSamples << " taylor=" << statistics.taylorPanels << " roots=" << statistics.rootSolves << " fallback=" << statistics.fallbackComponents << " mapped=" << statistics.mappedPanels << " regular=" << statistics.regularGraphPanels << " pairs=" << statistics.pairPanels << " firstFallback=" << statistics.firstFallbackIndex << " worstPairIndex=" << statistics.worstPairIndex << " worstPairRatio=" << static_cast<double>(statistics.worstPairRatio) << " pairValue=" << static_cast<double>(statistics.worstPairValue) << " pairError=" << static_cast<double>(statistics.worstPairError) << " pairTolerance=" << static_cast<double>(statistics.worstPairTolerance) << " pairLeft=" << static_cast<double>(statistics.worstPairLeft) << " pairRight=" << static_cast<double>(statistics.worstPairRight);
#endif
      std::cout << "\n";
    }
    catch(const std::exception& e) {
      ++failures;
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
      std::cout << "CASE g=" << g << " degree=" << degrees[d] << " tol=" << tolerances[t] << " sign=" << sign << " FAIL time=" << seconds << " error=" << e.what();
#ifdef FEMUS_LINE_DOMAIN_INTEGRATION_HPP
      const FELevelSetMoments::LineIntegrationStatistics& statistics = FELevelSetMoments::LineIntegration::Statistics();
      std::cout << " fast=" << statistics.fastSamples << " full=" << statistics.fullSamples << " taylor=" << statistics.taylorPanels << " roots=" << statistics.rootSolves << " fallback=" << statistics.fallbackComponents << " mapped=" << statistics.mappedPanels << " regular=" << statistics.regularGraphPanels << " pairs=" << statistics.pairPanels << " firstFallback=" << statistics.firstFallbackIndex << " worstPairIndex=" << statistics.worstPairIndex << " worstPairRatio=" << static_cast<double>(statistics.worstPairRatio) << " pairValue=" << static_cast<double>(statistics.worstPairValue) << " pairError=" << static_cast<double>(statistics.worstPairError) << " pairTolerance=" << static_cast<double>(statistics.worstPairTolerance) << " pairLeft=" << static_cast<double>(statistics.worstPairLeft) << " pairRight=" << static_cast<double>(statistics.worstPairRight);
#endif
      std::cout << "\n";
    }
  }
  std::cout << "TRI circle probe cases=16 failures=" << failures << '\n';
  return failures == 0 ? 0 : 1;
}
