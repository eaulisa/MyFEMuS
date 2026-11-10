#include "CutFemWeight_LS.hpp"
#include "Basis.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

long double WideMeshField(const std::vector<long double>& point) {
  const long double center[3] = {0.071L,-0.043L,0.029L};
  long double value = 0.613L*0.613L;
  for(unsigned d = 0; d < point.size(); ++d) value -= (point[d]-center[d])*(point[d]-center[d]);
  return value;
}

int main(int argc, char** argv) {
  if(argc != 4) return 64;
  const unsigned element = std::strtoul(argv[1],nullptr,10);
  const unsigned degree = std::strtoul(argv[2],nullptr,10);
  const unsigned maximum = std::strtoul(argv[3],nullptr,10);
  if(element > 3 || degree < 2 || degree > 10 || maximum > 64) return 64;
  const unsigned dimension = element < 2 ? 2U : 3U;
  const bool simplex = element%2U != 0;
  const GeomElType types[4] = {QUAD,TRI,HEX,TET};
  const unsigned nodes[4] = {9,7,27,15};
  const unsigned permutations[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
  femus::QuadBiquadratic quad;
  femus::TriBiquadratic tri;
  femus::HexBiquadratic hex;
  femus::TetBiquadratic tet;
  femus::basis* basis = element == 0 ? static_cast<femus::basis*>(&quad) : (element == 1 ? static_cast<femus::basis*>(&tri) : (element == 2 ? static_cast<femus::basis*>(&hex) : static_cast<femus::basis*>(&tet)));
  CutFemWeight_LS<long double,long double> rule(types[element],degree,"legendre");
  std::cout << std::setprecision(18) << "element,degree,divisions,cells,cut_cells,volume_error,surface_error,volume_x2_error,surface_x2_error,quad_splits,triangle_splits,cell_splits,seconds,status\n";
  unsigned failures = 0;
  for(unsigned divisions = 2; divisions <= maximum; divisions *= 2) {
    const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
    const long double h = 2.0L/divisions;
    const long double center[3] = {0.071L,-0.043L,0.029L};
    const long double r2 = 0.613L*0.613L;
    long double sums[4] = {};
    unsigned cells = 0;
    unsigned cutCells = 0;
    FELevelSetMoments::LineIntegration::Statistics() = FELevelSetMoments::LineIntegrationStatistics();
    FELevelSetMoments::GraphCellIntegration::Statistics() = FELevelSetMoments::GraphCellStatistics();
    try {
      for(unsigned i = 0; i < divisions; ++i) for(unsigned j = 0; j < divisions; ++j) for(unsigned k = 0; k < (dimension == 3 ? divisions : 1U); ++k) {
        const unsigned grid[3] = {i,j,k};
        long double minimum = 0;
        long double maximumDistance = 0;
        for(unsigned d = 0; d < dimension; ++d) {
          const long double low = -1+h*grid[d];
          const long double high = low+h;
          const long double distance = center[d] < low ? low-center[d] : (center[d] > high ? center[d]-high : 0);
          minimum += distance*distance;
          maximumDistance += std::max((low-center[d])*(low-center[d]),(high-center[d])*(high-center[d]));
        }
        for(unsigned side = 0; side < (simplex ? (dimension == 2 ? 2U : 6U) : 1U); ++side) {
          ++cells;
          if(minimum >= r2) continue;
          long double vertex[4][3] = {};
          if(dimension == 2) {
            if(side == 0) { vertex[1][0] = 1; vertex[2][1] = 1; }
            else { vertex[0][0] = vertex[0][1] = 1; vertex[1][1] = 1; vertex[2][0] = 1; }
          }
          else { vertex[1][permutations[side][0]] = 1; vertex[2][permutations[side][0]] = vertex[2][permutations[side][1]] = 1; vertex[3][0] = vertex[3][1] = vertex[3][2] = 1; }
          const bool cut = maximumDistance > r2;
          std::vector<long double> volume(rule.GetGaussQuadraturePointNumber(),1);
          std::vector<long double> surface(volume.size(),0);
          if(cut) {
            ++cutCells;
            std::vector<double> nodal(nodes[element]);
            std::vector<std::vector<double>> geometry(dimension,std::vector<double>(nodes[element]));
            for(unsigned n = 0; n < nodes[element]; ++n) {
              std::vector<long double> point(dimension);
              for(unsigned d = 0; d < dimension; ++d) {
                long double coordinate = simplex ? vertex[0][d] : (basis->GetXcoarse(n)[d]+1)/2;
                if(simplex) for(unsigned axis = 0; axis < dimension; ++axis) coordinate += (vertex[axis+1U][d]-vertex[0][d])*basis->GetXcoarse(n)[axis];
                point[d] = -1+h*grid[d]+h*coordinate;
                geometry[d][n] = point[d];
              }
              nodal[n] = WideMeshField(point);
            }
            const std::chrono::steady_clock::time_point cellStart = std::chrono::steady_clock::now();
            if(std::getenv("FEMUS_WIDE_TRACE")) std::cerr << "begin cell " << divisions << ' ' << cells << std::endl;
            rule.GetPhysicalSurfaceWeightsFromNodalGeometry(nodal,geometry,volume,surface,1e-11,1e-11);
            if(std::getenv("FEMUS_WIDE_TRACE")) std::cerr << "end cell " << divisions << ' ' << cells << " seconds " << std::chrono::duration<double>(std::chrono::steady_clock::now()-cellStart).count() << std::endl;
          }
          long double determinant = 1;
          for(unsigned d = 0; d < dimension; ++d) determinant *= simplex ? h : h/2;
          for(unsigned point = 0; point < volume.size(); ++point) {
            long double x = simplex ? vertex[0][0] : (rule.GetGaussCoordinatePointer(0)[point]+1)/2;
            if(simplex) for(unsigned axis = 0; axis < dimension; ++axis) x += (vertex[axis+1U][0]-vertex[0][0])*rule.GetGaussCoordinatePointer(axis)[point];
            x = -1+h*i+h*x;
            const long double weight = rule.GetGaussWeightPointer()[point];
            sums[0] += determinant*weight*volume[point]; sums[1] += weight*surface[point]; sums[2] += determinant*weight*volume[point]*x*x; sums[3] += weight*surface[point]*x*x;
          }
        }
      }
      const long double exactVolume = dimension == 2 ? std::acos(-1.0L)*r2 : 4*std::acos(-1.0L)*r2*0.613L/3;
      const long double exactSurface = dimension*exactVolume/0.613L;
      const long double exact[4] = {exactVolume,exactSurface,exactVolume*(center[0]*center[0]+r2/(dimension+2)),exactSurface*(center[0]*center[0]+r2/dimension)};
      bool success = true;
      for(unsigned n = 0; n < 4; ++n) success = success && std::isfinite(sums[n]) && std::abs(sums[n]-exact[n]) < 2e-9L;
      if(!success) ++failures;
      std::cout << element << ',' << degree << ',' << divisions << ',' << cells << ',' << cutCells;
      for(unsigned n = 0; n < 4; ++n) std::cout << ',' << std::abs(sums[n]-exact[n]);
      std::cout << ',' << FELevelSetMoments::LineIntegration::Statistics().splitCells << ',' << FELevelSetMoments::LineIntegration::Statistics().triangleSplitCells << ',' << FELevelSetMoments::GraphCellIntegration::Statistics().splitCells << ',' << std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count() << ',' << (success ? "PASS" : "FAIL") << std::endl;
    }
    catch(const std::exception& error) { ++failures; std::cerr << "element=" << element << " q=" << degree << " n=" << divisions << ": " << error.what() << '\n'; }
  }
  return failures == 0 ? 0 : 1;
}
