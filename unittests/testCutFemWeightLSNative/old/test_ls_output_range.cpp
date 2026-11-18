#include "CutFemWeight_LS.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
  CutFemWeight_LS<double,long double> rule(LINE,0,"legendre");
  std::vector<double> positive(2,17.0), delta(2,19.0);
  bool rejected = false;
  try { rule.GetWeightsFromNodalValues({-1.0e-320,1.0e-320,0.0},positive,delta); }
  catch(const std::overflow_error&) { rejected = true; }
  if(!rejected || positive != std::vector<double>(2,17.0) || delta != std::vector<double>(2,19.0)) throw std::runtime_error("An overflowing correction changed output vectors.");
  rule.GetWeightsFromNodalValues({-1.0e-320,1.0e-320,0.0},positive,delta,1.0e-12,1.0e-12,true);
  if(std::abs(positive[0]-0.5) > 1.0e-14 || std::abs(delta[0]-0.5) > 1.0e-14) throw std::runtime_error("Surface recovery after rejected raw delta failed.");
  std::vector<double> cachedPositive, cachedDelta;
  rule.GetWeightsFromNodalValues({-1.0e-320,1.0e-320,0.0},cachedPositive,cachedDelta,1.0e-12,1.0e-12,true);
  if(cachedPositive != positive || cachedDelta != delta) throw std::runtime_error("Cache changed after range recovery.");
  CutFemWeight_LS<boost::multiprecision::cpp_bin_float_oct,long double> wide(LINE,0,"legendre");
  std::vector<boost::multiprecision::cpp_bin_float_oct> widePositive, wideDelta;
  wide.GetWeightsFromNodalValues({-1.0e-320,1.0e-320,0.0},widePositive,wideDelta);
  if(!Eigen::numext::isfinite(wideDelta[0]) || !(wideDelta[0] > boost::multiprecision::cpp_bin_float_oct("1e319"))) throw std::runtime_error("A supported wider output range was rejected.");
  std::cout << "Output-range checks passed: overflow rejection, atomic outputs, surface recovery, cache, and multiprecision output.\n";
}
