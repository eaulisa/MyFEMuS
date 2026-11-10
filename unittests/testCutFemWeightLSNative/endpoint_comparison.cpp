#include "Quad9DirectionSelection.hpp"
#include "Tri7DirectionSelection.hpp"
#include "Hex27DirectionSelection.hpp"
#include "Tet15LevelSetMoments.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

long double ComparisonField(unsigned element, unsigned field, const std::vector<double>& point) {
  const bool simplex = element == 2 || element == 4;
  const unsigned dimension = element == 0 ? 1U : (element < 3 ? 2U : 3U);
  const long double last = point[dimension-1U];
  if(field == 0 || field == 1 || field == 4) {
    long double value = last-(field == 4 ? 1.0L-1.0L/256 : (simplex ? 0.2L : 0.15L));
    if(field == 1) for(unsigned axis = 0; axis+1U < dimension; ++axis) value += 1e-8L*point[axis]*point[axis];
    return value;
  }
  if(field == 3) return (last-(simplex ? 0.2L : -0.4L))*(last-(simplex ? 0.6L : 0.4L));
  const long double center = simplex ? 0.25L : 0;
  const long double radius = field == 5 ? 1.0L/128 : (simplex ? 0.125L : 0.6L);
  long double value = radius*radius;
  for(unsigned axis = 0; axis < dimension; ++axis) value -= (point[axis]-center)*(point[axis]-center);
  return value;
}

std::vector<double> ComparisonNodal(unsigned element, unsigned field) {
  const unsigned hex[27][3] = {{0,0,0},{2,0,0},{2,2,0},{0,2,0},{0,0,2},{2,0,2},{2,2,2},{0,2,2},{1,0,0},{2,1,0},{1,2,0},{0,1,0},{1,0,2},{2,1,2},{1,2,2},{0,1,2},{0,0,1},{2,0,1},{2,2,1},{0,2,1},{1,0,1},{2,1,1},{1,2,1},{0,1,1},{1,1,0},{1,1,2},{1,1,1}};
  const double tet[15][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},{0.5,0,0},{0.5,0.5,0},{0,0.5,0},{0,0,0.5},{0.5,0,0.5},{0,0.5,0.5},{1.0/3,1.0/3,0},{1.0/3,0,1.0/3},{1.0/3,1.0/3,1.0/3},{0,1.0/3,1.0/3},{0.25,0.25,0.25}};
  const double tri[7][2] = {{0,0},{1,0},{0,1},{0.5,0},{0.5,0.5},{0,0.5},{1.0/3,1.0/3}};
  const unsigned quad[9] = {0,1,2,3,8,9,10,11,24};
  const unsigned counts[5] = {3,9,7,27,15};
  std::vector<double> result(counts[element]);
  for(unsigned node = 0; node < result.size(); ++node) {
    std::vector<double> point(3,0);
    if(element == 0) point[0] = node == 0 ? -1 : (node == 1 ? 1 : 0);
    else if(element == 1 || element == 3) for(unsigned axis = 0; axis < 3; ++axis) point[axis] = static_cast<int>(hex[element == 1 ? quad[node] : node][axis])-1;
    else if(element == 2) for(unsigned axis = 0; axis < 2; ++axis) point[axis] = tri[node][axis];
    else for(unsigned axis = 0; axis < 3; ++axis) point[axis] = tet[node][axis];
    result[node] = static_cast<double>(ComparisonField(element,field,point));
  }
  return result;
}

void ComparisonIntegrate(unsigned element, const std::vector<double>& nodal, unsigned degree, double tolerance, std::vector<long double>& values) {
  std::vector<std::vector<double>> h2;
  std::vector<std::vector<double>> d2;
  std::vector<std::vector<std::vector<double>>> h3;
  std::vector<std::vector<std::vector<double>>> d3;
  if(element == 0) {
    double a;
    double b;
    double c;
    FELevelSetMoments::Line3ReferenceCoefficientsFromNodalValues(nodal,a,b,c);
    std::vector<long double> positive;
    std::vector<long double> delta;
    FELevelSetMoments::Integrate1DInternal(a,b,c,degree,-1,1,true,false,positive,delta);
    values = positive;
    const long double derivative = a == 0 ? std::abs(b) : std::sqrt(std::fma(-4.0L*a,c,static_cast<long double>(b)*b));
    for(unsigned n = 0; n < delta.size(); ++n) delta[n] *= derivative;
    values.insert(values.end(),delta.begin(),delta.end());
    return;
  }
#ifdef FEMUS_ENDPOINT_ONLY
  if(element == 1) FELevelSetMoments::IntegrateQuad9NodalSelected(nodal,degree,h2,d2,tolerance,tolerance,true,false,1,0,nullptr,true);
  if(element == 2) FELevelSetMoments::IntegrateTri7NodalSelected(nodal,degree,h2,d2,tolerance,tolerance,true,false,1,0,nullptr,true);
  if(element == 3) FELevelSetMoments::IntegrateHex27NodalSelected(nodal,degree,h3,d3,tolerance,tolerance,true,1,0,nullptr,nullptr,true);
#else
  if(element == 1) FELevelSetMoments::IntegrateQuad9NodalSelected(nodal,degree,h2,d2,tolerance,tolerance,true,false,1,0,nullptr,true,true);
  if(element == 2) FELevelSetMoments::IntegrateTri7NodalSelected(nodal,degree,h2,d2,tolerance,tolerance,true,false,1,0,nullptr,true,true);
  if(element == 3) FELevelSetMoments::IntegrateHex27NodalSelected(nodal,degree,h3,d3,tolerance,tolerance,true,1,0,nullptr,true,nullptr,true);
#endif
  if(element == 4) FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,degree,h3,d3,tolerance,tolerance,true,24,0,200000,nullptr,1,nullptr,true);
  values.clear();
  for(unsigned pair = 0; pair < 2; ++pair) for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
    if(element < 3) values.push_back(pair == 0 ? h2[i][j] : d2[i][j]);
    else for(unsigned k = 0; k <= degree-i-j; ++k) values.push_back(pair == 0 ? h3[i][j][k] : d3[i][j][k]);
  }
}

int main(int argc, char** argv) {
  if(argc != 6) return 64;
  const unsigned element = std::strtoul(argv[1],nullptr,10);
  const unsigned degree = std::strtoul(argv[2],nullptr,10);
  const unsigned field = std::strtoul(argv[3],nullptr,10);
  const double tolerance = std::strtod(argv[4],nullptr);
  const unsigned repetitions = std::strtoul(argv[5],nullptr,10);
  if(element > 4 || degree > 10 || field > 5 || repetitions == 0) return 64;
  try {
    const std::vector<double> nodal = ComparisonNodal(element,field);
    std::vector<long double> values;
    ComparisonIntegrate(element,nodal,degree,tolerance,values);
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    for(unsigned repeat = 0; repeat < repetitions; ++repeat) ComparisonIntegrate(element,nodal,degree,tolerance,values);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/repetitions;
    std::cout << std::setprecision(21) << "{\"seconds\":" << seconds << ",\"moments\":[";
    for(unsigned n = 0; n < values.size(); ++n) std::cout << (n == 0 ? "" : ",") << values[n];
    std::cout << "]}\n";
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
