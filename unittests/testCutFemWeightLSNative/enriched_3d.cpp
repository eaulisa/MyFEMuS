#include "Hex27DirectionSelection.hpp"
#include "Tet15LevelSetMoments.hpp"
#include <cstdlib>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
  if(argc != 5) return 64;
  const bool tetrahedron = std::string(argv[1]) == "TET";
  const unsigned degree = std::strtoul(argv[2],nullptr,10);
  const std::string mode = argv[3];
  const bool surface = mode.find("surface") != std::string::npos;
  const bool negative = mode.find("negative") != std::string::npos;
  const double tolerance = std::strtod(argv[4],nullptr);
  if(degree > 10 || !(tolerance > 0.0)) return 64;
  try {
    std::vector<std::vector<std::vector<double>>> positive;
    std::vector<std::vector<std::vector<double>>> delta;
    if(tetrahedron) {
      std::vector<double> nodal(15,-1.0/512.0);
      nodal[14] = 1.0/512.0;
      if(negative) for(unsigned node = 0; node < nodal.size(); ++node) nodal[node] = -nodal[node];
      FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,degree,positive,delta,tolerance,tolerance,true,24,0,200000,nullptr,1,nullptr,surface);
    }
    else {
      std::vector<std::vector<std::vector<double>>> coefficient(3,std::vector<std::vector<double>>(3,std::vector<double>(3,0.0)));
      coefficient[0][0][0] = 0.35;
      coefficient[0][0][2] = 0.2;
      coefficient[2][0][0] = coefficient[0][2][0] = -1.0;
      coefficient[2][0][2] = coefficient[0][2][2] = -0.8;
      coefficient[2][2][0] = 0.2;
      if(negative) for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) coefficient[i][j][k] = -coefficient[i][j][k];
      FELevelSetMoments::IntegrateHex27ReferenceSelected(coefficient,degree,positive,delta,tolerance,tolerance,true,24,0,200000,nullptr,1,true,nullptr,nullptr,surface);
    }
    std::cout << std::setprecision(21);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) std::cout << i << ' ' << j << ' ' << k << ' ' << positive[i][j][k] << ' ' << delta[i][j][k] << '\n';
  }
  catch(const std::exception& exception) { std::cerr << exception.what() << '\n'; return 1; }
}
