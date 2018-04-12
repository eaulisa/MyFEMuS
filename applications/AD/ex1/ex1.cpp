
#include "FemusInit.hpp"
#include "adept.h"

#include <vector>

using namespace femus;

void example1();
void example2();


int main(int argc, char** args)
{

  example2();

  return 0;
}


void example1()
{

  adept::Stack& s = FemusInit::_adeptStack;

  adept::adouble x, y, z;
  adept::adouble f, g, h;

  x = 2;
  y = 3;
  z = 1;

  s.new_recording();

  f = x * y * z;
  h = x * x;

  g = sin(f + h);

  s.dependent(&f, 1);
  s.dependent(&h, 1);
  s.dependent(&g, 1);
  s.independent(&x, 1);
  s.independent(&y, 1);
  s.independent(&z, 1);

  double jac[9];

  s.jacobian(jac, true);

  for (unsigned i = 0; i < 3; i++) {
    for(unsigned j = 0; j < 3; j++) {
      std::cout << jac[ i * 3 + j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  s.clear_independents();
  s.clear_dependents();
}


void example2()
{

  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned n1 = 3;
  std::vector < adept::adouble > f(n1);

  const unsigned n2 = 10;
  std::vector < adept::adouble > x(n2);
  

  for(unsigned i = 0; i < n2; i++) {
    x[i] = i * i;
  }

  s.new_recording();

  f[0] = 0.;
  for(unsigned j = 0; j < n2; j++) {
    f[0] += j * x[j];
  }
  
  for(unsigned i = 1; i < n1; i++) {
    f[i] = f[i-1] * f[i-1];
    for(unsigned j = 0; j < n1; j++) {
      f[i] += j * i * x[j] ;
    }
  }

  s.dependent(&f[0], n1);
  s.independent(&x[0], n2);

  double jac[n1][n2];

  s.jacobian(jac[0], true);

  for (unsigned i = 0; i < n1; i++) {
    for(unsigned j = 0; j < n2; j++) {
      std::cout << jac[i][ j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  s.clear_independents();
  s.clear_dependents();
}
