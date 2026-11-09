#include "TriangleIntegration.hpp"
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

void CheckValue(long double actual, long double expected, long double tolerance, const char* name) {
  const long double error = std::abs(actual-expected);
  if(!std::isfinite(actual) || error > tolerance) {
    std::cerr << name << " actual=" << actual << " expected=" << expected << " error=" << error << '\n';
    throw std::runtime_error("TRI parent-coordinate splitting regression failed.");
  }
}

int main() {
  const long double cx = 0.65L;
  const long double cy = 0.35L;
  const long double radius = 0.20L;
  const long double pi = std::acos(-1.0L);
  std::vector<std::vector<long double>> coefficient(3,std::vector<long double>(3,0.0L));
  coefficient[0][0] = cx*cx+cy*cy-radius*radius;
  coefficient[1][0] = -2.0L*cx;
  coefficient[0][1] = -2.0L*cy;
  coefficient[2][0] = 1.0L;
  coefficient[0][2] = 1.0L;
  const std::vector<std::vector<long double>> original = coefficient;
  const std::vector<FELevelSetMoments::TriangleDomainVertex> vertices = {{0.0L,0.0L},{1.0L,0.0L},{1.0L,1.0L}};
  std::vector<long double> values;
  std::vector<long double> errors;
  unsigned evaluations = 0;
  unsigned accepted = 0;
  unsigned direction = 0;
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  FELevelSetMoments::TriangleIntegration::Integrate(coefficient,1,vertices,true,1.0e-11L,1.0e-11L,1.0e-11L,28,200000,0,10,values,errors,evaluations,accepted,direction);
  const std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
  if(coefficient != original) throw std::runtime_error("TRI splitting modified the parent level-set polynomial.");
  const long double circleArea = pi*radius*radius;
  const long double expectedArea = 0.5L-circleArea;
  const long double expectedX = 1.0L/3.0L-circleArea*cx;
  const long double expectedY = 1.0L/6.0L-circleArea*cy;
  const long double expectedDelta = pi;
  const long double expectedDeltaX = pi*cx;
  const long double expectedDeltaY = pi*cy;
  const unsigned count = 3;
  CheckValue(values[0],expectedArea,3.0e-10L,"H00");
  CheckValue(values[1],expectedY,3.0e-10L,"H01");
  CheckValue(values[2],expectedX,3.0e-10L,"H10");
  CheckValue(values[count],expectedDelta,3.0e-10L,"D00");
  CheckValue(values[count+1U],expectedDeltaY,3.0e-10L,"D01");
  CheckValue(values[count+2U],expectedDeltaX,3.0e-10L,"D10");
  const double seconds = std::chrono::duration<double>(finish-start).count();
  std::cout << "TRI parent-coordinate circle split passed; evaluations=" << evaluations << " accepted=" << accepted << " direction=" << direction << " time=" << seconds << " s\n";
}
