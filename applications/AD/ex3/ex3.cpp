#include "FemusInit.hpp"
#include "adept.h"
#include <vector>

using namespace femus;

int main(int argc, char** args)
{
  adept::Stack& s = FemusInit::_adeptStack;
  
  adept::adouble x, y, z;
  adept::adouble f, g, h;
  
  x = 2;
  y = 3;
  z = 5;
  
  s.new_recording();
  
  f = x * (1. - y);
  g = x * y * (1. - z);
  h = x * y * z;
  
  s.dependent(&f, 1);
  s.dependent(&g, 1);
  s.dependent(&h, 1);
  s.independent(&x, 1);
  s.independent(&y, 1);
  s.independent(&z, 1);
  
  double jac[3][3];
  
  s.jacobian(jac[0], true);
  
  for (unsigned i = 0; i < 3; i++) {
    for(unsigned j = 0; j < 3; j++) {
      std::cout << jac[ i][ j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  s.clear_independents();
  s.clear_dependents();

  return 0;


}