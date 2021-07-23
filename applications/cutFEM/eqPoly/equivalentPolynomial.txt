#include <iostream>
#include <stdio.h>      /* printf */
#include <math.h>       /* exp */
using namespace std;




//Fake function just to be a place holder for the Lisk Li function
int Li(int i, double j)
{
    return i * j;
}

// This function accepts: the diminsion, dim; polynomial degree, degree; generalized Heavyside funcion paprameter, p;
//discontinuity (point, line, or plane) ax + by + cz = d, where b and c are set to zero by default;
//The function returns an array of the equivalent polynomial coefficients, coefficients
double * getCoefficients (int dim, int degree, double p, double a, double d, double b = 0, double c = 0)
{
  if(dim == 1) {
      
       //TODO
      
      if(degree == 3) {
          
          static double coefficients[4];
       
          coefficients[0] = 2*(-1 + (2*p + log(1 + exp(-1 + d)*p)) - log(1 + exp((1 + d)*p)))/p;

          coefficients[1] = (2*(p*log(1 + exp((-1 - d)*p)) + p*log(1 + exp(p - d*p)) - Li(2,-exp((-1 - d)*p)) + Li(2,-exp(p - d*p))))/pow(p,2);

          coefficients[2] =-0.6666666666666666 + (2*(-(pow(p,2)*log(1 + exp((-1 - d)*p))) + pow(p,2)*log(1 + exp(p - d*p)) + 
            2*p*Li(2,-exp((-1 - d)*p)) + 2*p*Li(2,-exp(p - d*p)) + 2*Li(3,-exp((-1 - d)*p)) - 
            2*Li(3,-exp(p - d*p))))/pow(p,3);

          coefficients[3] = (2*(pow(p,3)*log(1 + exp((-1 - d)*p)) + pow(p,3)*log(1 + exp(p - d*p)) - 3*pow(p,2)*Li(2,-exp((-1 - d)*p)) + 
            3*pow(p,2)*Li(2,-exp(p - d*p)) - 6*p*Li(3,-exp((-1 - d)*p)) - 6*p*Li(3,-exp(p - d*p)) - 
            6*Li(4,-exp((-1 - d)*p)) + 6*Li(4,-exp(p - d*p))))/pow(p,4);
          
          return coefficients;  
      }
      
  }
  
  else if(dim == 2) {
      double coefficients[((degree + 1) * (degree + 2)) / 2];
  
      //TODO
      
  }
  
  else if (dim == 3) {
      double coefficients[((degree + 1) * (degree + 2) * (degree + 3)) / 6];
      
       //TODO
  
  
  }
  
  
}

int main ()
{
  double* ptr = getCoefficients(1, 3, 100, 1, 0.5);
  printf("%d %d %d %d", ptr[0], ptr[1],ptr[2], ptr[3]); 
  return 0;
}
