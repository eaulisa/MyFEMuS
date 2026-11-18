#include "CutFemWeight.hpp"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

void Vector(const std::vector<double>& values) {
  std::cout << '[';
  for(unsigned i = 0; i < values.size(); ++i) { if(i != 0) std::cout << ','; if(std::isfinite(values[i])) std::cout << values[i]; else std::cout << "null"; }
  std::cout << ']';
}

int main(int argc, char** argv) {
  if(argc != 6) return 64;
  const unsigned element = std::strtoul(argv[1],nullptr,10);
  const unsigned degree = std::strtoul(argv[2],nullptr,10);
  const std::string field = argv[3];
  const double scale = std::strtod(argv[4],nullptr);
  const unsigned orientation = std::strtoul(argv[5],nullptr,10);
  const GeomElType elements[5] = {LINE,QUAD,TRI,HEX,TET};
  const unsigned dimension = element == 0 ? 1U : (element < 3 ? 2U : 3U);
  const bool simplex = element == 2 || element == 4;
  std::vector<double> normal(dimension,0.);
  double location = simplex ? 0.25 : 0.125;
  if(field == "boundary_low") location = simplex ? 0. : -1.;
  if(field == "boundary_high") location = 1.;
  if(field == "tiny") location = 1.-1./256.;
  normal[field == "plane_y" && dimension > 1 ? 1U : 0U] = 1.;
  if(field == "diagonal") { normal.assign(dimension,3.); location = 1.5; }
  for(unsigned d = 0; d < dimension; ++d) normal[d] *= scale;
  const double offset = -scale * location;
  CutFemWeight<double,long double> rule(elements[element],degree,"legendre");
  rule.SetTetBaseType(orientation);
  std::vector<double> h,delta;
  rule(0,normal,offset,h);
  rule(-1,normal,offset,delta);
  std::cout << std::setprecision(21) << "{\"h\":";
  Vector(h);
  std::cout << ",\"delta\":";
  Vector(delta);
  std::cout << ",\"gauss\":[";
  for(unsigned i = 0; i < rule.GetGaussQuadraturePointNumber(); ++i) { if(i != 0) std::cout << ','; std::cout << rule.GetGaussWeightPointer()[i]; }
  std::cout << "],\"points\":[";
  for(unsigned i = 0; i < rule.GetGaussQuadraturePointNumber(); ++i) {
    if(i != 0) std::cout << ',';
    std::cout << '[';
    for(unsigned d = 0; d < dimension; ++d) { if(d != 0) std::cout << ','; std::cout << rule.GetGaussCoordinatePointer(d)[i]; }
    std::cout << ']';
  }
  std::cout << "]}\n";
  return 0;
}
