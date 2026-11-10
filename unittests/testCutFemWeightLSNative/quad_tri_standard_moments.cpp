#include "Quad9DirectionSelection.hpp"
#include "Tri7DirectionSelection.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

const double quadNodes[9][2] = {{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0},{0,0}};
const double triNodes[7][2] = {{0,0},{1,0},{0,1},{0.5,0},{0.5,0.5},{0,0.5},{1.0/3.0,1.0/3.0}};
unsigned momentChecks = 0;
unsigned momentCases = 0;
long double largestError = 0.0L;

long double ProbePower(long double value, unsigned exponent) {
  long double result = 1.0L;
  for(unsigned k = 0; k < exponent; ++k) result *= value;
  return result;
}

long double ProbeFactorial(unsigned index) {
  long double result = 1.0L;
  for(unsigned k = 2; k <= index; ++k) result *= k;
  return result;
}

long double ProbeInterval(unsigned exponent, long double lower, long double upper) {
  return (ProbePower(upper,exponent+1U)-ProbePower(lower,exponent+1U))/(exponent+1U);
}

long double ProbeTriangle(unsigned i, unsigned j) {
  return ProbeFactorial(i)*ProbeFactorial(j)/ProbeFactorial(i+j+2U);
}

long double ProbeTriangleCut(unsigned i, unsigned j, long double cut) {
  long double result = 0.0L;
  for(unsigned k = 0; k <= i; ++k) result += ProbeFactorial(i)/ProbeFactorial(i-k)*ProbePower(cut,i-k)*ProbePower(1.0L-cut,k+j+2U)*ProbeFactorial(j)/ProbeFactorial(k+j+2U);
  return result;
}

void ProbeCheck(double actual, long double expected, const std::string& label) {
  const long double error = std::abs(static_cast<long double>(actual)-expected);
  largestError = std::max(largestError,error);
  if(!std::isfinite(actual) || error > 5.0e-11L) {
    std::cerr << label << " actual=" << actual << " expected=" << expected << " error=" << error << '\n';
    throw std::runtime_error("A native two-dimensional moment is incorrect.");
  }
  ++momentChecks;
}

void ProbeCut(bool triangle, unsigned axis, double cut, bool quadratic, unsigned direction, unsigned metric) {
  const unsigned nodes = triangle ? 7U : 9U;
  std::vector<double> nodal(nodes);
  std::vector<std::vector<double>> coordinates(2,std::vector<double>(nodes));
  const std::vector<std::vector<double>> inverse = {{0.5,0},{0,1.0/3.0}};
  for(unsigned node = 0; node < nodes; ++node) {
    const double x = triangle ? triNodes[node][0] : quadNodes[node][0];
    const double y = triangle ? triNodes[node][1] : quadNodes[node][1];
    const double value = 3.0*(axis == 0 ? x : y);
    nodal[node] = quadratic ? value*value-9.0*cut*cut : value-3.0*cut;
    coordinates[0][node] = 2.0*x+0.5*x*y;
    coordinates[1][node] = 3.0*y+0.25*x*y;
  }
  std::vector<std::vector<double>> positive;
  std::vector<std::vector<double>> delta;
  std::vector<double> estimated;
  if(triangle) FELevelSetMoments::IntegrateTri7NodalSelected(nodal,3,positive,delta,1.0e-11,1.0e-11,true,false,0,direction,&estimated,metric != 0,metric == 3 ? &inverse : nullptr,metric == 2 ? &coordinates : nullptr);
  else FELevelSetMoments::IntegrateQuad9NodalSelected(nodal,3,positive,delta,1.0e-11,1.0e-11,true,false,0,direction,&estimated,metric != 0,metric == 3 ? &inverse : nullptr,metric == 2 ? &coordinates : nullptr);
  for(unsigned i = 0; i <= 3; ++i) {
    for(unsigned j = 0; j <= 3-i; ++j) {
      const unsigned fixed = axis == 0 ? i : j;
      const unsigned free = axis == 0 ? j : i;
      long double expectedH = triangle ? ProbeTriangleCut(fixed,free,cut) : ProbeInterval(fixed,cut,1.0L)*ProbeInterval(free,-1.0L,1.0L);
      long double expectedD = triangle ? ProbePower(cut,fixed)*ProbeInterval(free,0.0L,1.0L-cut) : ProbePower(cut,fixed)*ProbeInterval(free,-1.0L,1.0L);
      if(quadratic && !triangle) {
        expectedH += ProbeInterval(fixed,-1.0L,-cut)*ProbeInterval(free,-1.0L,1.0L);
        expectedD += ProbePower(-cut,fixed)*ProbeInterval(free,-1.0L,1.0L);
      }
      long double multiplier = quadratic ? 1.0L/(18.0L*cut) : 1.0L/3.0L;
      if(metric == 1) multiplier = 1.0L;
      if(metric == 2) multiplier = axis == 0 ? std::sqrt(ProbePower(0.5L*cut,2)+ProbePower(3.0L+0.25L*cut,2)) : std::sqrt(ProbePower(2.0L+0.5L*cut,2)+ProbePower(0.25L*cut,2));
      if(metric == 3) multiplier = axis == 0 ? 0.5L : 1.0L/3.0L;
      if(cut == (triangle ? 0.0 : -1.0)) multiplier *= 0.5L;
      ProbeCheck(positive[i][j],expectedH,"cut H");
      ProbeCheck(delta[i][j],multiplier*expectedD,"cut delta/surface");
    }
  }
  ++momentCases;
}

void ProbeUniform(bool triangle, double constant) {
  std::vector<double> nodal(triangle ? 7U : 9U,constant);
  std::vector<std::vector<double>> positive;
  std::vector<std::vector<double>> delta;
  if(triangle) FELevelSetMoments::IntegrateTri7NodalSelected(nodal,5,positive,delta,1.0e-12,1.0e-12,constant != 0.0);
  else FELevelSetMoments::IntegrateQuad9NodalSelected(nodal,5,positive,delta,1.0e-12,1.0e-12,constant != 0.0);
  const long double factor = constant > 0.0 ? 1.0L : (constant == 0.0 ? 0.5L : 0.0L);
  for(unsigned i = 0; i <= 5; ++i) {
    for(unsigned j = 0; j <= 5-i; ++j) {
      ProbeCheck(positive[i][j],factor*(triangle ? ProbeTriangle(i,j) : ProbeInterval(i,-1,1)*ProbeInterval(j,-1,1)),"uniform H");
      ProbeCheck(delta[i][j],0.0L,"uniform delta");
    }
  }
  ++momentCases;
}

void ProbeTolerance(bool triangle) {
  std::vector<double> nodal(triangle ? 7U : 9U,1.0);
  std::vector<std::vector<double>> positive(1,std::vector<double>(1,17.0));
  std::vector<std::vector<double>> delta(1,std::vector<double>(1,19.0));
  bool rejected = false;
  try {
    if(triangle) FELevelSetMoments::IntegrateTri7NodalSelected(nodal,3,positive,delta,1.0e-22,0.0);
    else FELevelSetMoments::IntegrateQuad9NodalSelected(nodal,3,positive,delta,1.0e-22,0.0);
  }
  catch(const std::runtime_error&) { rejected = true; }
  if(!rejected || positive.size() != 1 || positive[0][0] != 17.0 || delta.size() != 1 || delta[0][0] != 19.0) throw std::runtime_error("An unattainable tolerance was not rejected atomically.");
  ++momentChecks;
  ++momentCases;
}

int main() {
  for(unsigned shape = 0; shape < 2; ++shape) {
    ProbeUniform(shape != 0,1.0);
    ProbeUniform(shape != 0,0.0);
    ProbeUniform(shape != 0,-1.0);
    ProbeTolerance(shape != 0);
    for(unsigned direction = 1; direction <= (shape != 0 ? 3U : 2U); ++direction) {
      for(unsigned repeat = 0; repeat < 1; ++repeat) {
        for(unsigned axis = 0; axis < 2; ++axis) {
          for(unsigned metric = 0; metric < 4; ++metric) ProbeCut(shape != 0,axis,0.25,false,direction,metric);
          ProbeCut(shape != 0,axis,shape != 0 ? 0.0 : -1.0,false,direction,1);
        }
        ProbeCut(shape != 0,0,0.5,true,direction,0);
        ProbeCut(shape != 0,0,0.5,true,direction,1);
      }
    }
  }
  std::cout << momentChecks << " native QUAD/TRI checks passed in " << momentCases << " cases; maximum absolute error " << largestError << '\n';
}
