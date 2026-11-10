#include "EndpointProjection.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

long double EndpointModel(unsigned multiplicity, unsigned exponent, long double parameter) {
  return exponent*std::pow(parameter,exponent*(0.5L+1.0L/multiplicity)-1);
}

int main() {
  unsigned checks = 0;
  for(unsigned multiplicity = 2; multiplicity <= 4; ++multiplicity) {
    std::vector<std::vector<long double>> constraint(multiplicity+1U,std::vector<long double>(1,0));
    constraint[0] = {0,1};
    constraint[multiplicity][0] = -1;
    const unsigned exponent = FELevelSetMoments::EndpointProjection::ContactExponent(constraint,0,false);
    if(exponent != (multiplicity == 3 ? 6U : multiplicity)) throw std::runtime_error("Projected contact exponent is incorrect.");
    long double integral = 0;
    for(unsigned q = 0; q < boost::math::quadrature::gauss_kronrod<long double,31>::abscissa().size(); ++q) for(unsigned side = 0; side < (q == 0 ? 1U : 2U); ++side) {
      const long double point = 0.5L+(side == 0 ? -0.5L : 0.5L)*boost::math::quadrature::gauss_kronrod<long double,31>::abscissa()[q];
      integral += 0.5L*boost::math::quadrature::gauss_kronrod<long double,31>::weights()[q]*EndpointModel(multiplicity,exponent,point);
    }
    if(std::abs(integral-1.0L/(0.5L+1.0L/multiplicity)) > 1e-17L) throw std::runtime_error("Transformed analytical endpoint quadrature is incorrect.");
    if(!std::isfinite(EndpointModel(multiplicity,2,1e-15L))) throw std::runtime_error("The quadratic endpoint map did not regularize the model.");
    checks += 3;
  }
  for(unsigned contact = 3; contact <= 4; ++contact) {
    std::vector<long double> coefficient(27,0);
    coefficient[FELevelSetMoments::Hex27Index(1,0,0)] = 1;
    coefficient[FELevelSetMoments::Hex27Index(0,0,2)] = 1;
    coefficient[FELevelSetMoments::Hex27Index(0,2,1)] = 1;
    if(contact == 3) { coefficient[FELevelSetMoments::Hex27Index(0,1,1)] = 2; coefficient[FELevelSetMoments::Hex27Index(0,2,0)] = 1; }
    const std::vector<std::vector<long double>> discriminant = FELevelSetMoments::EndpointProjection::Discriminant(FELevelSetMoments::EndpointProjection::Field(coefficient,0,2));
    if(discriminant.size() != 5 || FELevelSetMoments::EndpointProjection::ContactExponent(discriminant,0,false) != (contact == 3 ? 6U : 4U)) throw std::runtime_error("The full triquadratic discriminant was truncated.");
    ++checks;
  }
  for(unsigned contact = 3; contact <= 4; ++contact) {
    std::vector<long double> coefficient(27,0);
    coefficient[FELevelSetMoments::Hex27Index(1,0,0)] = 1;
    coefficient[FELevelSetMoments::Hex27Index(0,0,2)] = 1;
    coefficient[FELevelSetMoments::Hex27Index(1,1,1)] = -1;
    coefficient[FELevelSetMoments::Hex27Index(0,2,1)] = -1;
    coefficient[FELevelSetMoments::Hex27Index(0,1,2)] = -1;
    if(contact == 3) { coefficient[FELevelSetMoments::Hex27Index(0,1,1)] = 1; coefficient[FELevelSetMoments::Hex27Index(0,2,0)] = 0.25L; }
    const std::vector<std::vector<long double>> discriminant = FELevelSetMoments::EndpointProjection::Discriminant(FELevelSetMoments::EndpointProjection::Field(coefficient,0,2));
    if(discriminant.size() != 5 || FELevelSetMoments::EndpointProjection::ContactExponent(discriminant,0,true) != (contact == 3 ? 6U : 4U)) throw std::runtime_error("An enriched TET contact exponent is incorrect.");
    ++checks;
  }
  bool rejected = false;
  if(FELevelSetMoments::EndpointProjection::ContactExponent({{0.25L+1e-17L,1},{-1},{1},{1e-17L}},0,false,true) != 2) throw std::runtime_error("An outer-event rounding residual concealed a quadratic contact.");
  ++checks;
  try { FELevelSetMoments::EndpointProjection::ContactExponent({{0,0,1},{0},{-1}},0,false); }
  catch(const std::runtime_error&) { rejected = true; }
  if(!rejected) throw std::runtime_error("A non-transverse projection was reported as regular.");
  std::cout << checks+1U << " endpoint exponent and exact-discriminant checks passed\n";
}
