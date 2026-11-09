#include "LineIntegration.hpp"
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

long double LineTestPower(long double x, unsigned degree) {
  long double value = 1.0L;
  for(unsigned k = 0; k < degree; ++k) value *= x;
  return value;
}

long double LineTestIntegral(long double left, long double right, unsigned degree) {
  return (LineTestPower(right,degree+1U)-LineTestPower(left,degree+1U))/(degree+1U);
}

int main() {
  const unsigned degree = 10;
  std::vector<long double> positive;
  std::vector<long double> delta;
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  for(unsigned repeat = 0; repeat < 10000; ++repeat) FELevelSetMoments::LineIntegration::IntegratePolynomial(1.0L,0.0L,-0.25L,degree,-1.0L,1.0L,true,false,positive,delta);
  const std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
  long double maximumError = 0.0L;
  for(unsigned k = 0; k <= degree; ++k) {
    const long double exactPositive = LineTestIntegral(-1.0L,-0.5L,k)+LineTestIntegral(0.5L,1.0L,k);
    const long double exactDelta = LineTestPower(-0.5L,k)+LineTestPower(0.5L,k);
    maximumError = std::max(maximumError,std::abs(positive[k]-exactPositive));
    maximumError = std::max(maximumError,std::abs(delta[k]-exactDelta));
  }
  if(maximumError > 5.0e-18L) throw std::runtime_error("LINE quadratic moment regression failed.");
  const double seconds = std::chrono::duration<double>(finish-start).count();
  std::cout << "LINE degree-10 exact moments passed; maximum absolute error " << maximumError << "; 10000 integrations " << seconds << " s\n";
}
