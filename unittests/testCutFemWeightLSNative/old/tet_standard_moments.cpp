#include "Tet15LevelSetMoments.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

const double tetBarycentric[15][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},{0.5,0.5,0,0},{0,0.5,0.5,0},{0.5,0,0.5,0},{0.5,0,0,0.5},{0,0.5,0,0.5},{0,0,0.5,0.5},{1.0/3.0,1.0/3.0,1.0/3.0,0},{1.0/3.0,1.0/3.0,0,1.0/3.0},{0,1.0/3.0,1.0/3.0,1.0/3.0},{1.0/3.0,0,1.0/3.0,1.0/3.0},{0.25,0.25,0.25,0.25}};
unsigned tetChecks = 0;
unsigned tetCases = 0;
long double tetMaximumError = 0.0L;

long double TetFactorial(unsigned n) {
  long double value = 1.0L;
  for(unsigned i = 2; i <= n; ++i) value *= i;
  return value;
}

long double TetPower(long double value, unsigned exponent) {
  long double result = 1.0L;
  for(unsigned i = 0; i < exponent; ++i) result *= value;
  return result;
}

long double TetFullMoment(unsigned i, unsigned j, unsigned k) {
  return TetFactorial(i)*TetFactorial(j)*TetFactorial(k)/TetFactorial(i+j+k+3U);
}

long double TetPositiveMoment(unsigned vertex, long double cut, unsigned i, unsigned j, unsigned k) {
  const long double width = 1.0L-cut;
  if(vertex == 0) return TetPower(width,i+j+k+3U)*TetFullMoment(i,j,k);
  const unsigned powers[3] = {i,j,k};
  const unsigned fixed = powers[vertex-1U];
  const unsigned first = powers[vertex == 1 ? 1U : 0U];
  const unsigned second = powers[vertex == 3 ? 1U : 2U];
  long double value = 0.0L;
  for(unsigned r = 0; r <= fixed; ++r) value += TetFactorial(fixed)/TetFactorial(fixed-r)*TetPower(cut,fixed-r)*TetPower(width,r+first+second+3U)*TetFactorial(first)*TetFactorial(second)/TetFactorial(r+first+second+3U);
  return value;
}

long double TetPlaneMoment(unsigned vertex, long double cut, unsigned i, unsigned j, unsigned k) {
  const long double width = 1.0L-cut;
  if(vertex == 0) return TetPower(width,i+j+k+2U)*TetFactorial(i)*TetFactorial(j)*TetFactorial(k)/TetFactorial(i+j+k+2U);
  const unsigned powers[3] = {i,j,k};
  const unsigned fixed = powers[vertex-1U];
  const unsigned first = powers[vertex == 1 ? 1U : 0U];
  const unsigned second = powers[vertex == 3 ? 1U : 2U];
  return TetPower(cut,fixed)*TetPower(width,first+second+2U)*TetFactorial(first)*TetFactorial(second)/TetFactorial(first+second+2U);
}

void CheckTetValue(double actual, long double expected, const std::string& label) {
  const long double error = std::abs(static_cast<long double>(actual)-expected);
  tetMaximumError = std::max(tetMaximumError,error);
  if(!std::isfinite(actual) || error > 1.0e-10L) {
    std::cerr << label << " actual=" << actual << " expected=" << expected << " error=" << error << '\n';
    throw std::runtime_error("A standard tetrahedral moment is incorrect.");
  }
  ++tetChecks;
}

void CheckTetFull(unsigned direction) {
  std::vector<double> nodal(15,1.0);
  std::vector<std::vector<std::vector<double>>> positive;
  std::vector<std::vector<std::vector<double>>> delta;
  FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,3,positive,delta,1.0e-11,1.0e-11,true,24,direction);
  for(unsigned i = 0; i <= 3; ++i) {
    for(unsigned j = 0; j <= 3-i; ++j) {
      for(unsigned k = 0; k <= 3-i-j; ++k) {
        CheckTetValue(positive[i][j][k],TetFullMoment(i,j,k),"full H");
        CheckTetValue(delta[i][j][k],0.0L,"full delta");
      }
    }
  }
  ++tetCases;
}

void CheckTetCut(unsigned vertex, double cut, unsigned direction, bool surface, bool quadratic, bool physical) {
  std::vector<double> nodal(15);
  std::vector<std::vector<double>> coordinates(3,std::vector<double>(15));
  for(unsigned node = 0; node < 15; ++node) {
    const double value = 3.0*tetBarycentric[node][vertex];
    nodal[node] = quadratic ? value*value-9.0*cut*cut : value-3.0*cut;
    for(unsigned d = 0; d < 3; ++d) coordinates[d][node] = (d+2.0)*tetBarycentric[node][d+1U];
  }
  std::vector<std::vector<std::vector<double>>> positive;
  std::vector<std::vector<std::vector<double>>> delta;
  FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,2,positive,delta,1.0e-11,1.0e-11,true,24,direction,200000,nullptr,1,nullptr,surface,nullptr,physical ? &coordinates : nullptr);
  long double multiplier = quadratic ? 1.0L/(18.0L*cut) : 1.0L/3.0L;
  if(surface) multiplier = vertex == 0 ? std::sqrt(3.0L) : 1.0L;
  if(physical) multiplier = vertex == 0 ? std::sqrt(244.0L) : (vertex == 1 ? 12.0L : (vertex == 2 ? 8.0L : 6.0L));
  if(cut == 0.0) multiplier *= 0.5L;
  for(unsigned i = 0; i <= 2; ++i) {
    for(unsigned j = 0; j <= 2-i; ++j) {
      for(unsigned k = 0; k <= 2-i-j; ++k) {
        CheckTetValue(positive[i][j][k],TetPositiveMoment(vertex,cut,i,j,k),"cut H");
        CheckTetValue(delta[i][j][k],multiplier*TetPlaneMoment(vertex,cut,i,j,k),"cut delta/surface");
      }
    }
  }
  ++tetCases;
}

int main() {
  for(unsigned direction = 1; direction <= 4; ++direction) CheckTetFull(direction);
  for(unsigned vertex = 0; vertex < 4; ++vertex) {
    CheckTetCut(vertex,0.5,vertex+1U,false,false,false);
    CheckTetCut(vertex,0.5,vertex+1U,true,false,false);
  }
  for(unsigned direction = 1; direction <= 4; ++direction) {
    CheckTetCut(1,0.5,direction,false,false,false);
    CheckTetCut(1,0.5,direction,true,false,false);
  }
  for(unsigned vertex = 0; vertex < 2; ++vertex) {
    CheckTetCut(vertex,0.0,vertex+1U,false,false,false);
    CheckTetCut(vertex,0.0,vertex+1U,true,false,false);
    CheckTetCut(vertex,0.5,vertex+1U,true,false,true);
  }
  CheckTetCut(1,0.5,2,false,true,false);
  CheckTetCut(1,0.5,2,true,true,false);
  std::cout << tetChecks << " standard TET moment checks passed in " << tetCases << " cases; maximum absolute error " << tetMaximumError << ".\n";
}
