#include "LineIntegration.hpp"
#include <iomanip>
#include <iostream>

long double LineMeshPower(long double value, unsigned exponent) {
  long double result = 1;
  for(unsigned n = 0; n < exponent; ++n) result *= value;
  return result;
}

int main() {
  const long double radius = 0.613L;
  const long double center = 0.071L;
  std::cout << std::setprecision(18) << "divisions,volume_error,surface_error,status\n";
  for(unsigned divisions = 4; divisions <= 64; divisions *= 2) {
    std::vector<long double> positive(11,0);
    std::vector<long double> surface(11,0);
    const long double scale = 1.0L/divisions;
    for(unsigned cell = 0; cell < divisions; ++cell) {
      const long double origin = -1+(2*cell+1)*scale;
      std::vector<double> nodal(3);
      const int node[3] = {-1,1,0};
      for(unsigned n = 0; n < 3; ++n) nodal[n] = radius*radius-LineMeshPower(origin+scale*node[n]-center,2);
      double a;
      double b;
      double c;
      FELevelSetMoments::Line3ReferenceCoefficientsFromNodalValues(nodal,a,b,c);
      std::vector<long double> h;
      std::vector<long double> d;
      FELevelSetMoments::LineIntegration::IntegratePolynomial(a,b,c,10,-1,1,true,false,h,d);
      for(unsigned n = 0; n <= 10; ++n) for(unsigned k = 0; k <= n; ++k) {
        const long double factor = scale*LineSliceUtilities::Binomial(n,k)*LineMeshPower(origin,n-k)*LineMeshPower(scale,k);
        positive[n] += factor*h[k];
        surface[n] += 2*radius*factor*d[k];
      }
    }
    long double volumeError = 0;
    long double surfaceError = 0;
    for(unsigned n = 0; n <= 10; ++n) {
      const long double h = (LineMeshPower(center+radius,n+1)-LineMeshPower(center-radius,n+1))/(n+1);
      const long double d = LineMeshPower(center+radius,n)+LineMeshPower(center-radius,n);
      volumeError = std::max(volumeError,std::abs(positive[n]-h));
      surfaceError = std::max(surfaceError,std::abs(surface[n]-d));
    }
    if(volumeError > 1e-12L || surfaceError > 1e-12L) throw std::runtime_error("A LINE mesh moment failed its analytical reference.");
    std::cout << divisions << ',' << volumeError << ',' << surfaceError << ",PASS\n";
  }
}
