#include "FELevelSetMoments.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

void EndpointLineMeasure(unsigned graph, unsigned repetitions, double& seconds, long double& checksum) {
  std::vector<long double> h;
  std::vector<long double> d;
  std::vector<long double> values;
  std::vector<long double> errors;
  std::vector<std::vector<long double>> coefficient(3,std::vector<long double>(3,0));
  coefficient[0][0] = -0.3L;
  coefficient[0][1] = 1;
  coefficient[2][0] = -0.1L;
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  for(unsigned n = 0; n < repetitions; ++n) {
    if(graph) {
      unsigned accepted = 0;
      values.assign(132,0);
      errors.assign(132,0);
      if(!FELevelSetMoments::IntegrateCanonicalGraph(coefficient,10,1,false,true,1,1e-10L,1e-10L,1e-10L,24,values,errors,accepted)) throw std::runtime_error("The retained fast graph kernel rejected its graph.");
    }
    else FELevelSetMoments::Integrate1DInternal(1,0,-0.25L,10,-1,1,true,false,h,d);
  }
  seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/repetitions;
  checksum = graph ? values[0] : h[0]+d[0];
}

int main(int argc, char** argv) {
  if(argc != 3) return 64;
  const unsigned graph = std::strtoul(argv[1],nullptr,10);
  const unsigned repetitions = std::strtoul(argv[2],nullptr,10);
  if(graph > 1 || repetitions == 0) return 64;
  double seconds;
  long double checksum;
  EndpointLineMeasure(graph,1,seconds,checksum);
  EndpointLineMeasure(graph,repetitions,seconds,checksum);
  const long double exact = graph ? 2.0L/3 : 3;
  if(std::abs(checksum-exact) > 1e-15L) throw std::runtime_error("A retained LINE kernel failed its analytical volume reference.");
  std::cout << std::setprecision(21) << "{\"seconds\":" << seconds << ",\"checksum\":" << checksum << "}\n";
}
