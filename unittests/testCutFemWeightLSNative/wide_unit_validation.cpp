#include "CutFemWeight_LS.hpp"
#include "Basis.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>

long double WideField(const std::string& field, const std::vector<long double>& point) {
  const unsigned dimension = point.size();
  if(field == "full") return 1;
  if(field == "empty") return -1;
  if(field == "plane") return point.back()-0.375L;
  if(field == "boundary_low") return point.back();
  if(field == "boundary_high") return point.back()-1;
  if(field == "oblique") { long double value = -0.43L; for(unsigned d = 0; d < dimension; ++d) value += point[d]; return value; }
  if(field == "slab") return (point[0]-0.21L)*(point[0]-0.68L);
  if(field == "near_plane") { long double value = point.back()-0.25L; for(unsigned d = 0; d+1U < dimension; ++d) value += 1e-8L*point[d]*point[d]; return value; }
  if(field == "bubble") { long double value = 1; long double remaining = 1; for(unsigned d = 0; d < dimension; ++d) { value *= point[d]; remaining -= point[d]; } return value*remaining-(dimension == 2 ? 1.0L/64 : 1.0L/512); }
  const long double radius = field == "tiny_ball" ? 1.0L/128 : (dimension == 2 ? 0.15L : 0.10L);
  const long double centers[3] = {field == "touch_ball" ? radius : 0.29L,0.23L,0.19L};
  const long double axes[3] = {0.18L,0.07L,0.05L};
  long double value = field == "ellipse" ? 1 : radius*radius;
  for(unsigned d = 0; d < dimension; ++d) { const long double x = point[d]-centers[d]; value -= field == "ellipse" ? x*x/(axes[d]*axes[d]) : x*x; }
  return value;
}

void WideGeometry(const std::vector<long double>& point, std::vector<long double>& physical) {
  physical = point;
  physical[0] += 0.1L*point[0]*point[0];
  physical[1] += 0.12L*point[0]*point[0];
  if(point.size() == 3) physical[2] += 0.08L*point[0]*point[1];
}

void WideString(const std::string& value) {
  std::cout << '"';
  for(unsigned i = 0; i < value.size(); ++i) { if(value[i] == '\\' || value[i] == '"') std::cout << '\\'; if(value[i] == '\n') std::cout << "\\n"; else std::cout << value[i]; }
  std::cout << '"';
}

void WideMoments(unsigned element, unsigned degree, const std::vector<double>& nodal, bool surface, const std::vector<std::vector<double>>* geometry, double tolerance, unsigned validation, unsigned direction, std::vector<long double>& values) {
  const unsigned dimension = element < 2 ? 2U : 3U;
  if(element%2U == 0) { std::vector<long double> errors; FELevelSetMoments::IntegrateUnitTensorNodalMoments(nodal,dimension,degree,values,errors,tolerance,tolerance,surface,geometry,direction,validation); return; }
  if(element == 1) {
    std::vector<std::vector<double>> h;
    std::vector<std::vector<double>> d;
    FELevelSetMoments::IntegrateTri7NodalSelected(nodal,degree,h,d,tolerance,tolerance,true,false,validation,direction,nullptr,surface,nullptr,geometry);
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    values.assign(2U*count,0);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) { const unsigned n = FELevelSetMoments::UnitMomentIndex(2,i,j); values[n] = h[i][j]; values[count+n] = d[i][j]; }
  }
  else {
    std::vector<std::vector<std::vector<double>>> h;
    std::vector<std::vector<std::vector<double>>> d;
    FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,degree,h,d,tolerance,tolerance,true,24,direction,200000,nullptr,validation,nullptr,surface,nullptr,geometry);
    const unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
    values.assign(2U*count,0);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) { const unsigned n = FELevelSetMoments::UnitMomentIndex(3,i,j,k); values[n] = h[i][j][k]; values[count+n] = d[i][j][k]; }
  }
}

void WideBasis(CutFemWeight_LS<long double,long double>& rule, bool simplex, std::vector<std::vector<long double>>& values, std::vector<long double>& first) {
  const unsigned dimension = rule.GetDimension();
  const unsigned degree = rule.GetCutFEMQuadratureOrder();
  const unsigned count = dimension == 2 ? (degree+1U)*(degree+2U)/2U : (degree+1U)*(degree+2U)*(degree+3U)/6U;
  const unsigned points = rule.GetGaussQuadraturePointNumber();
  values.assign(count,std::vector<long double>(points));
  first.resize(points);
  for(unsigned point = 0; point < points; ++point) {
    std::vector<std::vector<long double>> power(dimension,std::vector<long double>(degree+1U,1));
    for(unsigned d = 0; d < dimension; ++d) { const long double x = rule.GetGaussCoordinatePointer(d)[point]; const long double u = simplex ? x : (x+1)/2; if(d == 0) first[point] = u; for(unsigned n = 1; n <= degree; ++n) power[d][n] = power[d][n-1U]*u; }
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= (dimension == 3 ? degree-i-j : 0U); ++k) values[FELevelSetMoments::UnitMomentIndex(dimension,i,j,k)][point] = power[0][i]*power[1][j]*(dimension == 3 ? power[2][k] : 1);
  }
}

int main(int argc, char** argv) {
  if(argc != 3) return 64;
  const unsigned element = std::strtoul(argv[1],nullptr,10);
  const unsigned degree = std::strtoul(argv[2],nullptr,10);
  if(element > 3 || degree > 10) return 64;
  const unsigned dimension = element < 2 ? 2U : 3U;
  const bool simplex = element%2U != 0;
  const GeomElType types[4] = {QUAD,TRI,HEX,TET};
  femus::QuadBiquadratic quad;
  femus::TriBiquadratic tri;
  femus::HexBiquadratic hex;
  femus::TetBiquadratic tet;
  femus::basis* basis = element == 0 ? static_cast<femus::basis*>(&quad) : (element == 1 ? static_cast<femus::basis*>(&tri) : (element == 2 ? static_cast<femus::basis*>(&hex) : static_cast<femus::basis*>(&tet)));
  const unsigned nodes[4] = {9,7,27,15};
  std::unique_ptr<CutFemWeight_LS<long double,long double>> rule;
  std::vector<std::vector<long double>> monomials;
  std::vector<long double> first;
  std::cout << std::setprecision(21);
  std::string field;
  double scale;
  double tolerance;
  unsigned geometry;
  unsigned fitting;
  unsigned validation;
  unsigned direction;
  while(std::cin >> field >> scale >> tolerance >> geometry >> fitting >> validation >> direction) {
    const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
    double setup = 0;
    std::string stage = "reconstruction";
    std::cout << "{\"field\":"; WideString(field);
    std::cout << ",\"scale\":" << scale << ",\"tolerance\":" << tolerance << ",\"geometry\":" << geometry << ",\"fitting\":" << fitting << ",\"validation\":" << validation << ",\"direction\":" << direction << std::flush;
    try {
      std::vector<double> nodal(nodes[element]);
      std::vector<std::vector<double>> coordinates(dimension,std::vector<double>(nodes[element]));
      for(unsigned node = 0; node < nodes[element]; ++node) {
        std::vector<long double> point(dimension);
        for(unsigned d = 0; d < dimension; ++d) point[d] = simplex ? basis->GetXcoarse(node)[d] : (basis->GetXcoarse(node)[d]+1)/2;
        nodal[node] = scale*WideField(field,point);
        std::vector<long double> physical;
        WideGeometry(point,physical);
        for(unsigned d = 0; d < dimension; ++d) coordinates[d][node] = physical[d];
      }
      FELevelSetMoments::LineIntegration::Statistics() = FELevelSetMoments::LineIntegrationStatistics();
      FELevelSetMoments::OuterIntervalQuadrature::Statistics(2) = FELevelSetMoments::OuterIntervalStatistics();
      FELevelSetMoments::OuterIntervalQuadrature::Statistics(3) = FELevelSetMoments::OuterIntervalStatistics();
      FELevelSetMoments::SurfaceEndpointIntegration::Statistics() = FELevelSetMoments::SurfaceEndpointStatistics();
      std::vector<long double> raw;
      std::vector<long double> surface;
      std::vector<long double> physical;
      stage = "raw";
      WideMoments(element,degree,nodal,false,nullptr,tolerance,validation,direction,raw);
      stage = "surface";
      WideMoments(element,degree,nodal,true,nullptr,tolerance,validation,direction,surface);
      stage = "physical";
      if(geometry != 0) WideMoments(element,degree,nodal,true,&coordinates,tolerance,validation,direction,physical);
      const unsigned count = raw.size()/2U;
      const FELevelSetMoments::SurfaceEndpointStatistics graphStats = FELevelSetMoments::SurfaceEndpointIntegration::Statistics();
      std::vector<long double> h;
      std::vector<long double> d;
      std::vector<long double> s;
      std::vector<long double> p;
      if(fitting != 0) {
        stage = "fitting";
        if(!rule) { const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now(); rule.reset(new CutFemWeight_LS<long double,long double>(types[element],degree,"legendre")); WideBasis(*rule,simplex,monomials,first); setup = std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count(); }
        rule->SetDirectionValidation(validation);
        rule->SetIntegrationDirection(direction);
        rule->GetWeightsFromNodalValues(nodal,h,d,tolerance,tolerance,false);
        std::vector<long double> unused;
        rule->GetWeightsFromNodalValues(nodal,unused,s,tolerance,tolerance,true);
        if(geometry != 0) rule->GetPhysicalSurfaceWeightsFromNodalGeometry(nodal,coordinates,unused,p,tolerance,tolerance);
      }
      std::cout << ",\"endpoint_strips\":" << graphStats.endpointStrips << ",\"general_surfaces\":" << graphStats.generalSurfaces << ",\"quad_calls\":" << graphStats.quadCalls << ",\"tri_calls\":" << graphStats.triCalls << ",\"mapped_panels\":" << FELevelSetMoments::OuterIntervalQuadrature::Statistics(2).transformedPanels << ",\"moments\":[";
      for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= (dimension == 3 ? degree-i-j : 0U); ++k) {
        const unsigned n = FELevelSetMoments::UnitMomentIndex(dimension,i,j,k);
        if(n != 0) std::cout << ',';
        std::cout << '[' << i << ',' << j << ',' << k << ',' << raw[n] << ',' << raw[count+n] << ',' << surface[count+n] << ',';
        if(geometry != 0) std::cout << physical[count+n]; else std::cout << "null";
        long double sums[5] = {};
        if(fitting != 0) for(unsigned point = 0; point < h.size(); ++point) {
          const long double weight = rule->GetGaussWeightPointer()[point]*monomials[n][point];
          sums[0] += weight*h[point]; sums[1] += weight*d[point]; sums[2] += weight*s[point];
          if(geometry != 0) { sums[3] += weight*p[point]; sums[4] += weight*h[point]*(1+0.2L*first[point]); }
        }
        const long double volumeMeasure = simplex ? 1 : (1U<<dimension);
        const long double surfaceMeasure = simplex ? 1 : (1U<<(dimension-1U));
        for(unsigned mode = 0; mode < 5; ++mode) { std::cout << ','; if(fitting == 0 || (mode >= 3 && geometry == 0)) std::cout << "null"; else std::cout << sums[mode]/(mode == 2 ? surfaceMeasure : (mode == 3 ? 1 : volumeMeasure)); }
        std::cout << ']';
      }
      std::cout << "],\"setup_seconds\":" << setup << ",\"seconds\":" << std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count() << ",\"error\":null}\n";
    }
    catch(const std::exception& error) { std::cout << ",\"seconds\":" << std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count() << ",\"error\":"; WideString(error.what()); std::cout << ",\"stage\":"; WideString(stage); std::cout << "}\n"; }
    std::cout << std::flush;
  }
}
