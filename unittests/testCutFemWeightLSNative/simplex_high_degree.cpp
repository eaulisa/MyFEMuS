#include "Tri7DirectionSelection.hpp"
#include "Tet15LevelSetMoments.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

const double simplexTriangleX[7] = {0,1,0,0.5,0.5,0,1.0/3.0};
const double simplexTetrahedronX[15] = {0,1,0,0,0.5,0.5,0,0,0.5,0,1.0/3.0,1.0/3.0,1.0/3.0,0,0.25};
unsigned simplexChecks = 0;
unsigned simplexFailures = 0;
long double simplexMaximumError = 0.0L;

long double SimplexFactorial(unsigned degree) {
  long double value = 1.0L;
  for(unsigned k = 2; k <= degree; ++k) value *= k;
  return value;
}

long double SimplexPower(long double value, unsigned degree) {
  long double result = 1.0L;
  for(unsigned k = 0; k < degree; ++k) result *= value;
  return result;
}

long double SimplexBinomial(unsigned degree, unsigned index) {
  return SimplexFactorial(degree)/(SimplexFactorial(index)*SimplexFactorial(degree-index));
}

long double SimplexExact(unsigned dimension, unsigned field, bool delta, unsigned i, unsigned j, unsigned k) {
  if(field == 0 && delta) return 0.0L;
  const long double cut = field == 1 ? 0.5L : 0.0L;
  if(delta) return (field == 2 ? 0.5L : 1.0L)*SimplexPower(cut,i)*SimplexPower(1.0L-cut,j+k+dimension-1U)*SimplexFactorial(j)*SimplexFactorial(k)/(3.0L*SimplexFactorial(j+k+dimension-1U));
  long double value = 0.0L;
  for(unsigned a = 0; a <= i; ++a) value += SimplexFactorial(i)/SimplexFactorial(i-a)*SimplexPower(cut,i-a)*SimplexPower(1.0L-cut,a+j+k+dimension)*SimplexFactorial(j)*SimplexFactorial(k)/SimplexFactorial(a+j+k+dimension);
  return value;
}

long double SimplexExpected(unsigned dimension, unsigned field, bool delta, unsigned i, unsigned j, unsigned k) {
#ifdef BASELINE_CHART
  long double result = 0.0L;
  if(dimension == 2) {
    for(unsigned a = 0; a <= i; ++a) result += SimplexBinomial(i,a)*SimplexPower(-1.0L,a)*SimplexExact(dimension,field,delta,a,j,0);
  }
  else {
    for(unsigned a = 0; a <= i; ++a) {
      for(unsigned b = 0; b <= i-a; ++b) {
        for(unsigned c = 0; c <= j; ++c) result += SimplexFactorial(i)/(SimplexFactorial(a)*SimplexFactorial(b)*SimplexFactorial(i-a-b))*SimplexBinomial(j,c)*SimplexExact(dimension,field,delta,a,b+c,i-a-b+j-c+k);
      }
    }
  }
  return result;
#else
  return SimplexExact(dimension,field,delta,i,j,k);
#endif
}

void CheckSimplexValue(double actual, long double expected) {
  const long double error = std::abs(static_cast<long double>(actual)-expected);
  simplexMaximumError = std::max(simplexMaximumError,error);
  if(!std::isfinite(actual) || error > 2.0e-12L) throw std::runtime_error("A high-degree simplex moment is incorrect.");
  ++simplexChecks;
}

void CheckSimplex(unsigned dimension, unsigned degree, unsigned field) {
  const unsigned nodes = dimension == 2 ? 7U : 15U;
  std::vector<double> nodal(nodes,1.0);
  for(unsigned node = 0; node < nodes; ++node) if(field != 0) nodal[node] = 3.0*(dimension == 2 ? simplexTriangleX[node] : simplexTetrahedronX[node])-(field == 1 ? 1.5 : 0.0);
  std::vector<double> error;
  try {
    if(dimension == 2) {
      std::vector<std::vector<double>> positive;
      std::vector<std::vector<double>> delta;
      FELevelSetMoments::IntegrateTri7NodalSelected(nodal,degree,positive,delta,1.0e-12,1.0e-12,true,false,1,0,&error);
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) {
          CheckSimplexValue(positive[i][j],SimplexExpected(dimension,field,false,i,j,0));
          CheckSimplexValue(delta[i][j],SimplexExpected(dimension,field,true,i,j,0));
        }
      }
    }
    else {
      std::vector<std::vector<std::vector<double>>> positive;
      std::vector<std::vector<std::vector<double>>> delta;
      FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,degree,positive,delta,1.0e-12,1.0e-12,true,24,0,200000,&error);
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) {
          for(unsigned k = 0; k <= degree-i-j; ++k) {
            CheckSimplexValue(positive[i][j][k],SimplexExpected(dimension,field,false,i,j,k));
            CheckSimplexValue(delta[i][j][k],SimplexExpected(dimension,field,true,i,j,k));
          }
        }
      }
    }
    std::cout << "PASS dimension=" << dimension << " degree=" << degree << " field=" << field << " estimatedH=" << error[0] << " estimatedD=" << error[1] << '\n';
  }
  catch(const std::exception& exception) {
    ++simplexFailures;
    std::cout << "FAIL dimension=" << dimension << " degree=" << degree << " field=" << field << " reason=" << exception.what() << '\n';
  }
}

int main() {
  for(unsigned dimension = 2; dimension <= 3; ++dimension) {
    for(unsigned degree = 8; degree <= 10; degree += 2) {
      for(unsigned field = 0; field < 3; ++field) CheckSimplex(dimension,degree,field);
    }
  }
  std::cout << simplexChecks << " high-degree simplex checks completed; failures=" << simplexFailures << "; maximum absolute error=" << simplexMaximumError << '\n';
  return simplexFailures == 0 ? 0 : 1;
}
