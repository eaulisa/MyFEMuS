#include "CutFemWeightRules_LS.hpp"
#include "Basis.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

double SurfaceBenchmarkElapsed(const std::chrono::steady_clock::time_point& start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}

void SurfaceBenchmarkNodal(unsigned cell, std::vector<double>& nodal, std::vector<std::vector<double>>& coordinates) {
  femus::HexBiquadratic basis;
  nodal.resize(27);
  coordinates.assign(3,std::vector<double>(27));
  for(unsigned n = 0; n < 27; ++n) {
    const double* point = basis.GetXcoarse(n);
    const double x = point[0], y = point[1], z = point[2];
    nodal[n] = z-0.2-0.003*cell-0.05*x*x-0.025*x*y+0.01*x*x*y*y*z*z;
    coordinates[0][n] = x+0.008*x*x*y*y*z*z;
    coordinates[1][n] = y+0.006*x*x*y*z*z;
    coordinates[2][n] = z+0.004*x*y*y*z*z;
  }
}

int main() {
  const std::chrono::steady_clock::time_point startup = std::chrono::steady_clock::now();
  CutFemWeightRules_LS<double,long double> rules(13,{HEX},FELevelSetMoments::OuterQuadratureRules());
  CutFemWeight_LS<double,long double>& rule = *rules.GetRule(HEX);
  const double setup = SurfaceBenchmarkElapsed(startup);
  const double* gauss = rule.GetGaussWeightPointer();
  double timings[2] = {0,0};
  long double checksum = 0;
  for(unsigned mode = 0; mode < 2; ++mode) {
    rule.ClearCache();
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    for(unsigned cell = 0; cell < 8; ++cell) {
      std::vector<double> nodal, positive, tangential;
      std::vector<std::vector<double>> coordinates, normal;
      SurfaceBenchmarkNodal(mode == 0 ? 0 : cell,nodal,coordinates);
      rule.GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,tangential,1e-8,1e-8);
      rule.GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,normal,1e-8,1e-8);
      if(gauss != rule.GetGaussWeightPointer() || rule.GetCutFEMQuadratureOrder() != 13) throw std::runtime_error("The fixed quadrature changed during cell evaluation.");
      for(unsigned n = 0; n < positive.size(); ++n) checksum += gauss[n]*(positive[n]+tangential[n]+normal[0][n]+normal[1][n]+normal[2][n]);
    }
    timings[mode] = SurfaceBenchmarkElapsed(start);
  }
  std::cout << std::setprecision(16) << "{\"setup_seconds\":" << setup << ",\"same_cell_seconds\":" << timings[0] << ",\"different_cells_seconds\":" << timings[1] << ",\"checksum\":" << checksum << "}\n";
}
