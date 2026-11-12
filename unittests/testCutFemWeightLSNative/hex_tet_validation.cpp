#include "Hex27DirectionSelection.hpp"
#include "Tet15LevelSetMoments.hpp"
#include "Basis.hpp"
#include <array>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

struct ValidationCell3D {
  std::vector<double> nodal;
  long double determinant;
  std::array<long double,3> x;
  long double origin;
};

long double ValidationField3D(const std::string& field, const std::array<long double,3>& point, bool tet) {
  const long double x = point[0];
  const long double y = point[1];
  const long double z = point[2];
  if(field == "plane") return z-0.25L;
  if(field == "tiny") return x-(1.0L-1.0L/256.0L);
  if(field == "near_plane") return z-0.25L+1.0e-8L*(x*x+y*y);
  if(field == "graph") return tet ? z-0.2L-0.08L*x*x+0.02L*x*y+0.768L*x*y*z*(1.0L-x-y-z) : z+0.05L-0.12L*x*x-0.07L*y*y+0.03L*x*x*z*z;
  const long double radius = field == "small_sphere" ? (tet ? 1.0L/128.0L : 1.0L/64.0L) : (tet ? 0.125L : 0.6L);
  const long double center = tet ? 0.25L : 0.0L;
  return radius*radius-(x-center)*(x-center)-(y-center)*(y-center)-(z-center)*(z-center);
}

void ValidationIntegrate3D(bool tet, const std::vector<double>& nodal, unsigned degree, double tolerance, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta) {
  if(tet) FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,degree,positive,delta,tolerance,tolerance,true,24,0,200000,nullptr,1);
  else FELevelSetMoments::IntegrateHex27NodalSelected(nodal,degree,positive,delta,tolerance,tolerance,true,1,0,nullptr);
}

void ValidationCase3D(bool tet, unsigned degree, const std::string& field, double tolerance) {
  femus::TetBiquadratic tetBasis;
  femus::HexBiquadratic hexBasis;
  const unsigned nodes = tet ? 15U : 27U;
  std::vector<double> nodal(nodes);
  for(unsigned node = 0; node < nodes; ++node) {
    const double* point = tet ? tetBasis.GetXcoarse(node) : hexBasis.GetXcoarse(node);
    nodal[node] = ValidationField3D(field,{point[0],point[1],point[2]},tet);
  }
  std::vector<double> seconds;
  double h = 0.0;
  double d = 0.0;
  for(unsigned repeat = 0; repeat < 3; ++repeat) {
    std::vector<std::vector<std::vector<double>>> positive;
    std::vector<std::vector<std::vector<double>>> delta;
    const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
    ValidationIntegrate3D(tet,nodal,degree,tolerance,positive,delta);
    seconds.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count());
    h = positive[0][0][0];
    d = delta[0][0][0];
  }
  std::sort(seconds.begin(),seconds.end());
  if(field == "sphere" || field == "small_sphere") {
    const long double radius = field == "small_sphere" ? (tet ? 1.0L/128.0L : 1.0L/64.0L) : (tet ? 0.125L : 0.6L);
    if(std::abs(h-4.0L*std::acos(-1.0L)*radius*radius*radius/3.0L) > 8.0L*tolerance || std::abs(d-2.0L*std::acos(-1.0L)*radius) > 8.0L*tolerance) throw std::runtime_error("The sphere benchmark failed its independent reference.");
  }
  std::cout << "PASS," << (tet ? "TET" : "HEX") << ',' << field << ',' << degree << ',' << tolerance << ',' << h << ',' << d << ',' << seconds[1] << '\n';
}

long double ValidationMeshField3D(const std::array<long double,3>& point) {
  const long double radius = 0.613L;
  return radius*radius-(point[0]-0.071L)*(point[0]-0.071L)-(point[1]+0.043L)*(point[1]+0.043L)-(point[2]-0.029L)*(point[2]-0.029L);
}

std::vector<ValidationCell3D> ValidationMesh3D(bool tet, unsigned divisions) {
  femus::TetBiquadratic tetBasis;
  femus::HexBiquadratic hexBasis;
  const unsigned permutation[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
  const long double h = 2.0L/divisions;
  std::vector<ValidationCell3D> cells;
  for(unsigned i = 0; i < divisions; ++i) for(unsigned j = 0; j < divisions; ++j) for(unsigned k = 0; k < divisions; ++k) for(unsigned side = 0; side < (tet ? 6U : 1U); ++side) {
    ValidationCell3D cell;
    cell.nodal.resize(tet ? 15U : 27U);
    cell.determinant = tet ? h*h*h : h*h*h/8.0L;
    std::array<long double,3> origin = {-1.0L+h*i,-1.0L+h*j,-1.0L+h*k};
    long double vertex[4][3] = {};
    vertex[1][permutation[side][0]] = 1.0L;
    vertex[2][permutation[side][0]] = vertex[2][permutation[side][1]] = 1.0L;
    vertex[3][0] = vertex[3][1] = vertex[3][2] = 1.0L;
    cell.origin = tet ? origin[0] : origin[0]+h/2.0L;
    cell.x = tet ? std::array<long double,3>{h*vertex[1][0],h*vertex[2][0],h*vertex[3][0]} : std::array<long double,3>{h/2.0L,0.0L,0.0L};
    for(unsigned node = 0; node < cell.nodal.size(); ++node) {
      const double* local = tet ? tetBasis.GetXcoarse(node) : hexBasis.GetXcoarse(node);
      std::array<long double,3> point;
      for(unsigned axis = 0; axis < 3; ++axis) point[axis] = origin[axis]+h*(tet ? vertex[1][axis]*local[0]+vertex[2][axis]*local[1]+vertex[3][axis]*local[2] : (1.0L+local[axis])/2.0L);
      cell.nodal[node] = ValidationMeshField3D(point);
    }
    cells.push_back(cell);
  }
  return cells;
}

long double ValidationPhysicalX2(const ValidationCell3D& cell, const std::vector<std::vector<std::vector<double>>>& moment) {
  long double value = cell.origin*cell.origin*moment[0][0][0];
  for(unsigned first = 0; first < 3; ++first) {
    unsigned powers[3] = {0,0,0};
    ++powers[first];
    value += 2.0L*cell.origin*cell.x[first]*moment[powers[0]][powers[1]][powers[2]];
    for(unsigned second = 0; second < 3; ++second) {
      ++powers[second];
      value += cell.x[first]*cell.x[second]*moment[powers[0]][powers[1]][powers[2]];
      --powers[second];
    }
  }
  return value;
}

void ValidationRefinement3D(bool tet, unsigned degree, unsigned divisions, double tolerance) {
  const std::vector<ValidationCell3D> cells = ValidationMesh3D(tet,divisions);
  long double volume = 0.0L;
  long double surface = 0.0L;
  long double volumeX2 = 0.0L;
  long double surfaceX2 = 0.0L;
  const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
  for(unsigned cell = 0; cell < cells.size(); ++cell) {
    std::vector<std::vector<std::vector<double>>> positive;
    std::vector<std::vector<std::vector<double>>> delta;
    try { ValidationIntegrate3D(tet,cells[cell].nodal,degree,tolerance,positive,delta); }
    catch(const std::exception& error) { throw std::runtime_error("Mesh cell "+std::to_string(cell)+": "+error.what()); }
    volume += cells[cell].determinant*positive[0][0][0];
    surface += 2.0L*0.613L*cells[cell].determinant*delta[0][0][0];
    volumeX2 += cells[cell].determinant*ValidationPhysicalX2(cells[cell],positive);
    surfaceX2 += 2.0L*0.613L*cells[cell].determinant*ValidationPhysicalX2(cells[cell],delta);
  }
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
  const long double exactVolume = 4.0L*std::acos(-1.0L)*0.613L*0.613L*0.613L/3.0L;
  const long double exactSurface = 4.0L*std::acos(-1.0L)*0.613L*0.613L;
  const long double errors[4] = {std::abs(volume-exactVolume),std::abs(surface-exactSurface),std::abs(volumeX2-exactVolume*(0.071L*0.071L+0.613L*0.613L/5.0L)),std::abs(surfaceX2-exactSurface*(0.071L*0.071L+0.613L*0.613L/3.0L))};
  for(unsigned n = 0; n < 4; ++n) if(!std::isfinite(errors[n]) || errors[n] > 2.0e-9L) throw std::runtime_error("A sphere mesh-refinement reference failed.");
  std::cout << "PASS," << (tet ? "TET" : "HEX") << ",mesh," << degree << ',' << divisions << ',' << cells.size();
  for(unsigned n = 0; n < 4; ++n) std::cout << ',' << errors[n];
  std::cout << ',' << seconds << '\n';
}

int main(int argc, char** argv) {
  if(argc != 6) return 64;
  const std::string mode = argv[1];
  const bool tet = std::string(argv[2]) == "TET";
  const unsigned degree = std::strtoul(argv[3],nullptr,10);
  const double tolerance = std::strtod(argv[5],nullptr);
  if(degree > 10 || degree < 2 || !(tolerance > 0.0)) return 64;
  std::cout << std::setprecision(18);
  try {
    if(mode == "case") ValidationCase3D(tet,degree,argv[4],tolerance);
    else if(mode == "mesh") ValidationRefinement3D(tet,degree,std::strtoul(argv[4],nullptr,10),tolerance);
    else return 64;
  }
  catch(const std::exception& exception) { std::cerr << exception.what() << '\n'; return 1; }
}
