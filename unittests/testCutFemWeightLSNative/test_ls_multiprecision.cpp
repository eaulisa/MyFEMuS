#include "CutFemWeight_LS.hpp"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
  CutFemWeight_LS<double,boost::multiprecision::cpp_bin_float_oct> rule(QUAD,2,"legendre");
  std::vector<double> positive;
  std::vector<double> surface;
  rule.GetWeightsFromNodalValues({-1.,1.,1.,-1.,0.,1.,0.,-1.,0.},positive,surface,1.e-8,1.e-8,true);
  double area = 0.;
  double length = 0.;
  for(unsigned i = 0; i < positive.size(); ++i) {
    area += positive[i] * rule.GetGaussWeightPointer()[i];
    length += surface[i] * rule.GetGaussWeightPointer()[i];
  }
  if(std::abs(area - 2.) > 1.e-12 || std::abs(length - 2.) > 1.e-12) throw std::runtime_error("Multiprecision level-set integral failed.");
  std::cout << "Standalone LS cpp_bin_float_oct: half-square area=" << area << ", reference interface length=" << length << ".\n";
  return 0;
}
