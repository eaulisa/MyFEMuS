#include "CutFemWeight.hpp"
#include "CutFemWeight_LS.hpp"
#include <type_traits>
#include <vector>

static_assert(!std::is_base_of<CutFemWeight<double,long double>,CutFemWeight_LS<double,long double>>::value,"Level-set class inherits from the planar class.");
static_assert(!std::is_convertible<CutFemWeight_LS<double,long double>*,CutFemWeight<double,long double>*>::value,"Level-set class exposes planar conversion.");

int main() {
  CutFemWeight<double,boost::multiprecision::cpp_bin_float_oct> planar(QUAD,1,"legendre");
  CutFemWeight_LS<double,boost::multiprecision::cpp_bin_float_oct> levelSet(QUAD,1,"legendre");
  std::vector<double> positive;
  std::vector<double> surface;
  const std::vector<double> nodal = {-1.,1.,1.,-1.,0.,1.,0.,-1.,0.};
  levelSet.SetIntegrationDirection(0);
  levelSet.SetDirectionValidation(0);
  levelSet.GetWeightsFromNodalValues(nodal,positive,surface,1.e-8,1.e-8,true);
  planar(0,{1.,0.},0.,positive);
  return 0;
}
