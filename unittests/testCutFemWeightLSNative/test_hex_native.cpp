#include "Hex27DirectionSelection.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstdlib>

long double Interval(unsigned power) { return power%2U == 0 ? 2.0L/(power+1U) : 0.0L; }

void CheckPlane(unsigned degree, double tolerance, double offset = 0.25) {
  const int node[27][3] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1},{0,-1,-1},{1,0,-1},{0,1,-1},{-1,0,-1},{0,-1,1},{1,0,1},{0,1,1},{-1,0,1},{-1,-1,0},{1,-1,0},{1,1,0},{-1,1,0},{0,-1,0},{1,0,0},{0,1,0},{-1,0,0},{0,0,-1},{0,0,1},{0,0,0}};
  std::vector<double> nodal(27);
  std::vector<std::vector<double>> geometry(3,std::vector<double>(27));
  for(unsigned n = 0; n < 27; ++n) {
    nodal[n] = node[n][0]-offset;
    for(unsigned axis = 0; axis < 3; ++axis) geometry[axis][n] = (axis+2.0)*node[n][axis];
  }
  std::vector<std::vector<std::vector<double>>> positive;
  std::vector<std::vector<std::vector<double>>> delta;
  std::vector<double> error;
  for(unsigned mode = 0; mode < 3; ++mode) {
    FELevelSetMoments::IntegrateHex27NodalSelected(nodal,degree,positive,delta,tolerance,tolerance,true,1,0,&error,true,nullptr,mode != 0,nullptr,mode == 2 ? &geometry : nullptr);
    long double mismatch = 0.0L;
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        for(unsigned k = 0; k <= degree-i-j; ++k) {
          const long double h = (1.0L-std::pow(static_cast<long double>(offset),i+1U))/(i+1U)*Interval(j)*Interval(k);
          const long double d = std::pow(static_cast<long double>(offset),i)*Interval(j)*Interval(k)*(mode == 2 ? 12.0L : 1.0L)*(std::abs(offset) == 1.0 ? 0.5L : 1.0L);
          mismatch = std::max(mismatch,std::abs(h-positive[i][j][k]));
          mismatch = std::max(mismatch,std::abs(d-delta[i][j][k]));
        }
      }
    }
    std::cout << "plane offset=" << offset << " mode=" << mode << " maxerr=" << mismatch << " estimated=" << error[0] << ',' << error[1] << " samples=" << error[3] << std::endl;
    if(mismatch > 8.0L*tolerance) throw std::runtime_error("Plane moments mismatch.");
  }
}

void CheckCurvedGeometry(unsigned degree, double tolerance) {
  const int node[27][3] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1},{0,-1,-1},{1,0,-1},{0,1,-1},{-1,0,-1},{0,-1,1},{1,0,1},{0,1,1},{-1,0,1},{-1,-1,0},{1,-1,0},{1,1,0},{-1,1,0},{0,-1,0},{1,0,0},{0,1,0},{-1,0,0},{0,0,-1},{0,0,1},{0,0,0}};
  std::vector<double> nodal(27);
  std::vector<std::vector<double>> geometry(3,std::vector<double>(27));
  for(unsigned n = 0; n < 27; ++n) {
    nodal[n] = node[n][2];
    for(unsigned axis = 0; axis < 3; ++axis) geometry[axis][n] = node[n][axis];
    geometry[2][n] += 0.3*node[n][0]*node[n][0];
  }
  std::vector<std::vector<std::vector<double>>> positive;
  std::vector<std::vector<std::vector<double>>> delta;
  std::vector<double> error;
  FELevelSetMoments::IntegrateHex27NodalSelected(nodal,degree,positive,delta,tolerance,tolerance,true,1,0,&error,true,nullptr,true,nullptr,&geometry);
  const long double expected = 2.0L*(std::sqrt(1.36L)+std::asinh(0.6L)/0.6L);
  const long double mismatch = std::abs(expected-delta[0][0][0]);
  std::cout << "curved_geometry area=" << delta[0][0][0] << " expected=" << expected << " maxerr=" << mismatch << " estimated=" << error[1] << " samples=" << error[3] << std::endl;
  if(mismatch > 8.0L*tolerance || std::abs(positive[0][0][0]-4.0) > tolerance) throw std::runtime_error("Curved physical area mismatch.");
}

void CheckParaboloid(unsigned degree, double tolerance) {
  std::vector<std::vector<std::vector<double>>> coefficient(3,std::vector<std::vector<double>>(3,std::vector<double>(3,0.0)));
  coefficient[0][0][0] = -0.2;
  coefficient[0][0][1] = 1.0;
  coefficient[2][0][0] = -0.1;
  coefficient[0][2][0] = -0.15;
  std::vector<std::vector<std::vector<double>>> positive;
  std::vector<std::vector<std::vector<double>>> delta;
  std::vector<double> error;
  FELevelSetMoments::IntegrateHex27ReferenceSelected(coefficient,degree,positive,delta,tolerance,tolerance,true,24,0,200000,&error,1,true);
  const long double h = 4.0L-0.2L*4.0L-(0.1L+0.15L)*4.0L/3.0L;
  const long double d = 4.0L;
  std::cout << "paraboloid H=" << positive[0][0][0] << " D=" << delta[0][0][0] << " expected=" << h << ',' << d << " estimated=" << error[0] << ',' << error[1] << " samples=" << error[3] << std::endl;
  if(std::abs(h-positive[0][0][0]) > 8.0L*tolerance || std::abs(d-delta[0][0][0]) > 8.0L*tolerance) throw std::runtime_error("Paraboloid mismatch.");
}

int main(int argc, char** argv) {
  unsigned degree = argc > 1 ? std::strtoul(argv[1],nullptr,10) : 5;
  double tolerance = argc > 2 ? std::strtod(argv[2],nullptr) : 1.0e-12;
  std::cout << std::setprecision(17);
  try {
    CheckPlane(degree,tolerance);
    CheckPlane(degree,tolerance,0.0);
    CheckPlane(degree,tolerance,-1.0);
    CheckPlane(degree,tolerance,1.0);
    CheckCurvedGeometry(degree,tolerance);
    CheckParaboloid(degree,tolerance);
  }
  catch(const std::exception& error) { std::cerr << error.what() << std::endl; return 1; }
  return 0;
}
