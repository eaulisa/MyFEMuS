#include "FELevelSetMoments.hpp"
#include "Hex27LevelSetMoments.hpp"
#include "Tet15LevelSet.hpp"
#include "Basis.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

femus::basis* CoefficientBasis(unsigned element) {
  if(element == 0) return new femus::LineBiquadratic();
  if(element == 1) return new femus::QuadBiquadratic();
  if(element == 2) return new femus::TriBiquadratic();
  if(element == 3) return new femus::HexBiquadratic();
  return new femus::TetBiquadratic();
}

long double CoefficientPower(long double value, unsigned exponent) {
  long double result = 1.0L;
  for(unsigned i = 0; i < exponent; ++i) result *= value;
  return result;
}

int main() {
  unsigned checks = 0;
  long double maximum = 0.0L;
  for(unsigned element = 0; element < 5; ++element) {
    std::unique_ptr<femus::basis> basis(CoefficientBasis(element));
    const unsigned dimension = element == 0 ? 1U : (element < 3 ? 2U : 3U);
    const bool simplex = element == 2 || element == 4;
    for(unsigned field = 0; field < 12; ++field) {
      std::vector<double> nodal(basis->_nc);
      for(unsigned node = 0; node < nodal.size(); ++node) nodal[node] = field == 0 ? 1.0 : std::sin((node+1.0)*(field+0.25))+0.0625*(node%3U);
      unsigned side = element == 4 ? 5U : 3U;
      std::vector<long double> coefficient(side*side*side,0.0L);
      if(element == 0) {
        double a,b,c;
        FELevelSetMoments::Line3ReferenceCoefficientsFromNodalValues(nodal,a,b,c);
        coefficient[0] = c;
        coefficient[side*side] = b;
        coefficient[2U*side*side] = a;
      }
      else if(dimension == 2) {
        std::vector<std::vector<double>> polynomial;
        if(element == 1) FELevelSetMoments::Quad9ReferenceCoefficientsFromNodalValues(nodal,polynomial);
        else FELevelSetMoments::Tri7ReferenceCoefficientsFromNodalValues(nodal,polynomial);
        for(unsigned i = 0; i < side; ++i) {
          for(unsigned j = 0; j < side; ++j) coefficient[side*side*i+side*j] = polynomial[i][j];
        }
      }
      else if(element == 3) {
        std::vector<std::vector<std::vector<double>>> polynomial;
        FELevelSetMoments::Hex27ReferenceCoefficientsFromNodalValues(nodal,polynomial);
        for(unsigned i = 0; i < side; ++i) {
          for(unsigned j = 0; j < side; ++j) {
            for(unsigned k = 0; k < side; ++k) coefficient[side*side*i+side*j+k] = polynomial[i][j][k];
          }
        }
      }
      else FELevelSetMoments::Tet15PolynomialFromNodalValues(nodal,false,coefficient);
      for(unsigned sample = 0; sample < 80; ++sample) {
        std::vector<double> point(dimension);
        double width = 1.0;
        for(unsigned d = 0; d < dimension; ++d) {
          const double fraction = static_cast<double>((37U*sample+23U*d+11U)%101U)/100.0;
          point[d] = simplex ? width*fraction : 2.0*fraction-1.0;
          if(simplex) width -= point[d];
        }
        long double value = 0.0L;
        std::vector<long double> gradient(dimension,0.0L);
        for(unsigned i = 0; i < side; ++i) {
          for(unsigned j = 0; j < side; ++j) {
            for(unsigned k = 0; k < side; ++k) {
              const long double c = coefficient[side*side*i+side*j+k];
              if(c == 0.0L) continue;
              const unsigned exponent[3] = {i,j,k};
              long double term = c;
              for(unsigned d = 0; d < dimension; ++d) term *= CoefficientPower(point[d],exponent[d]);
              value += term;
              for(unsigned axis = 0; axis < dimension; ++axis) {
                if(exponent[axis] == 0) continue;
                long double derivative = c*exponent[axis];
                for(unsigned d = 0; d < dimension; ++d) derivative *= CoefficientPower(point[d],exponent[d]-(d == axis ? 1U : 0U));
                gradient[axis] += derivative;
              }
            }
          }
        }
        long double feValue = 0.0L;
        std::vector<long double> feGradient(dimension,0.0L);
        for(unsigned node = 0; node < nodal.size(); ++node) {
          feValue += static_cast<long double>(nodal[node])*basis->eval_phi(node,point);
          for(unsigned axis = 0; axis < dimension; ++axis) feGradient[axis] += static_cast<long double>(nodal[node])*basis->eval_dphidxyz(axis,basis->GetIND(node),point.data());
        }
        maximum = std::max(maximum,std::abs(value-feValue));
        if(std::abs(value-feValue) > 2.0e-12L) throw std::runtime_error("Native coefficient value differs from FEMuS shape evaluation.");
        ++checks;
        for(unsigned axis = 0; axis < dimension; ++axis) {
          maximum = std::max(maximum,std::abs(gradient[axis]-feGradient[axis]));
          if(std::abs(gradient[axis]-feGradient[axis]) > 5.0e-12L) throw std::runtime_error("Native coefficient gradient differs from FEMuS shape evaluation.");
          ++checks;
        }
      }
    }
  }
  std::cout << std::setprecision(17) << "Native coefficients: " << checks << " independent FEMuS value/gradient checks; maximum absolute error=" << maximum << '\n';
}
