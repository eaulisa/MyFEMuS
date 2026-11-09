#include "Quad9DirectionSelection.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

unsigned quadrantChecks = 0;
long double quadrantError = 0.0L;

long double QuadrantPower(long double value, unsigned degree) {
  long double result = 1.0L;
  for(unsigned k = 0; k < degree; ++k) result *= value;
  return result;
}

long double QuadrantIntegral(unsigned degree, long double lower, long double upper) {
  return (QuadrantPower(upper,degree+1U)-QuadrantPower(lower,degree+1U))/(degree+1U);
}

long double QuadrantBinomial(unsigned degree, unsigned index) {
  long double result = 1.0L;
  for(unsigned k = 1; k <= index; ++k) result *= static_cast<long double>(degree-index+k)/k;
  return result;
}

void CheckQuadrantValue(double actual, long double expected, unsigned field, unsigned i, unsigned j) {
  const long double error = std::abs(static_cast<long double>(actual)-expected);
  quadrantError = std::max(quadrantError,error);
  if(!std::isfinite(actual) || error > 2.0e-12L) {
    std::cerr << "field=" << field << " i=" << i << " j=" << j << " actual=" << actual << " expected=" << expected << " error=" << error << '\n';
    throw std::runtime_error("A native quadrant moment is incorrect.");
  }
  ++quadrantChecks;
}

void CheckQuadrants(unsigned field, unsigned direction, bool taylor, bool fullBoundary) {
  std::vector<std::vector<double>> coefficient(3,std::vector<double>(3,0.0));
  coefficient[1][0] = 1.0;
  if(field == 0) coefficient[0][0] = -0.25;
  if(field == 2) coefficient[0][0] = 1.0;
  if(field == 3) coefficient[0][0] = -1.0;
  if(field == 4) { coefficient[1][0] = -1.0; coefficient[2][0] = 1.0; }
  if(field == 5 || field == 6) { coefficient[1][0] = 0.0; coefficient[2][0] = 1.0; coefficient[0][0] = field == 5 ? -1.0 : -0.25; }
  if(field == 7) { coefficient[0][0] = -0.25; coefficient[0][2] = -0.25; }
  std::vector<std::vector<double>> positive;
  std::vector<std::vector<double>> delta;
  FELevelSetMoments::IntegrateQuad9ReferenceSelected(coefficient,10,positive,delta,1.0e-12,1.0e-12,true,fullBoundary,28,direction,200000,nullptr,0,nullptr,taylor);
  const long double boundary = fullBoundary ? 1.0L : 0.5L;
  for(unsigned i = 0; i <= 10; ++i) {
    for(unsigned j = 0; j <= 10-i; ++j) {
      const long double transverse = QuadrantIntegral(j,-1.0L,1.0L);
      long double expectedH = 0.0L;
      long double expectedD = 0.0L;
      if(field <= 3) {
        const long double cut = field == 0 ? 0.25L : (field == 1 ? 0.0L : (field == 2 ? -1.0L : 1.0L));
        expectedH = QuadrantIntegral(i,cut,1.0L)*transverse;
        expectedD = QuadrantPower(cut,i)*transverse*(field >= 2 ? boundary : 1.0L);
      }
      if(field == 4) { expectedH = QuadrantIntegral(i,-1.0L,0.0L)*transverse; expectedD = (QuadrantPower(0.0L,i)+boundary)*transverse; }
      if(field == 5) expectedD = 0.5L*boundary*(QuadrantPower(-1.0L,i)+1.0L)*transverse;
      if(field == 6) { expectedH = (QuadrantIntegral(i,-1.0L,-0.5L)+QuadrantIntegral(i,0.5L,1.0L))*transverse; expectedD = (QuadrantPower(-0.5L,i)+QuadrantPower(0.5L,i))*transverse; }
      if(field == 7) {
        long double graphH = 0.0L;
        long double graphD = 0.0L;
        for(unsigned k = 0; k <= i+1U; ++k) graphH += QuadrantBinomial(i+1U,k)*QuadrantIntegral(j+2U*k,-1.0L,1.0L);
        for(unsigned k = 0; k <= i; ++k) graphD += QuadrantBinomial(i,k)*QuadrantIntegral(j+2U*k,-1.0L,1.0L);
        expectedH = (transverse-QuadrantPower(0.25L,i+1U)*graphH)/(i+1U);
        expectedD = QuadrantPower(0.25L,i)*graphD;
      }
      CheckQuadrantValue(positive[i][j],expectedH,field,i,j);
      CheckQuadrantValue(delta[i][j],expectedD,field,i,j);
    }
  }
}

void CheckLocalSplitDirections() {
  std::vector<std::vector<double>> coefficient(3,std::vector<double>(3,0.0));
  coefficient[2][0] = 1.0;
  coefficient[1][0] = 0.2;
  coefficient[0][2] = 1.0;
  coefficient[0][1] = -0.6;
  coefficient[0][0] = -0.3;
  unsigned directions[4] = {};
  unsigned split = 0;
  for(unsigned firstSide = 0; firstSide < 2; ++firstSide) {
    const long double firstSign = firstSide == 0 ? -1.0L : 1.0L;
    for(unsigned secondSide = 0; secondSide < 2; ++secondSide) {
      const long double secondSign = secondSide == 0 ? -1.0L : 1.0L;
      std::vector<std::vector<long double>> local;
      FELevelSetMoments::ReferencePolynomialAffine2D(coefficient,0.5L,0.5L*firstSign,0.5L,0.5L*secondSign,local);
      std::vector<std::vector<long double>> normalized;
      long double scale = 0.0L;
      LineSliceUtilities::NormalizeCanonicalPolynomial(local,normalized,scale);
      long double lower[2] = {};
      long double upper[2] = {};
      directions[split++] = FELevelSetMoments::Quad9SelectedOrder(normalized,lower,upper);
    }
  }
  bool first = false;
  bool second = false;
  for(unsigned k = 0; k < 4; ++k) {
    first = first || directions[k] == 1;
    second = second || directions[k] == 2;
  }
  if(!first || !second) throw std::runtime_error("QUAD split direction selection did not vary locally.");
  ++quadrantChecks;
}

int main() {
  CheckLocalSplitDirections();
  for(unsigned direction = 1; direction <= 2; ++direction) {
    for(unsigned taylor = 0; taylor < 2; ++taylor) {
      for(unsigned full = 0; full < 2; ++full) {
        for(unsigned field = 0; field < 8; ++field) CheckQuadrants(field,direction,taylor != 0,full != 0);
      }
    }
  }
  std::cout << quadrantChecks << " degree-10 quadrant checks passed; maximum absolute error " << quadrantError << '\n';
}
