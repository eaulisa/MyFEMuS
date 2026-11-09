#include "GraphCellIntegration.hpp"
#include "Hex27DirectionSelection.hpp"
#include "Tet15LevelSetMoments.hpp"
#include "Basis.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

unsigned graphChecks = 0;
unsigned graphCases = 0;
long double graphMaximumError = 0.0L;

long double GraphFactorial(unsigned degree) {
  long double value = 1.0L;
  for(unsigned n = 2; n <= degree; ++n) value *= n;
  return value;
}

void GraphCheck(long double actual, long double expected, long double tolerance, const char* label) {
  const long double error = std::abs(actual-expected);
  graphMaximumError = std::max(graphMaximumError,error);
  ++graphChecks;
  if(!std::isfinite(actual) || error > tolerance) {
    std::cerr << label << " actual=" << actual << " expected=" << expected << " error=" << error << " tolerance=" << tolerance << '\n';
    throw std::runtime_error("A 3D graph integration check failed.");
  }
}

void GraphGeometry() {
  FELevelSetMoments::GraphCellRegion cell;
  cell.tetrahedron = true;
  cell.vertices = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
  std::vector<FELevelSetMoments::GraphCellRegion> regions = {cell};
  for(unsigned step = 0; step < 5; ++step) {
    std::vector<FELevelSetMoments::GraphCellRegion> children;
    const long double cut = 0.17L+0.13L*step;
    for(unsigned n = 0; n < regions.size(); ++n) {
      FELevelSetMoments::GraphCellIntegration::Clip(regions[n],step%3U,cut,true,children);
      FELevelSetMoments::GraphCellIntegration::Clip(regions[n],step%3U,cut,false,children);
    }
    regions.swap(children);
    for(unsigned i = 0; i <= 5; ++i) for(unsigned j = 0; j <= 5-i; ++j) for(unsigned k = 0; k <= 5-i-j; ++k) {
      long double actual = 0.0L;
      for(unsigned n = 0; n < regions.size(); ++n) actual += FELevelSetMoments::GraphCellIntegration::TetMoment(regions[n],i,j,k);
      GraphCheck(actual,GraphFactorial(i)*GraphFactorial(j)*GraphFactorial(k)/GraphFactorial(i+j+k+3U),2.0e-17L,"clipped tet moment");
    }
  }
  std::vector<long double> parent(27,0.0L);
  parent[0] = 0.13L;
  parent[FELevelSetMoments::Hex27Index(2,1,0)] = -1.7L;
  parent[FELevelSetMoments::Hex27Index(1,1,2)] = 2.1L;
  parent[FELevelSetMoments::Hex27Index(1,2,0)] = 0.9L;
  for(unsigned n = 0; n < regions.size(); ++n) {
    const std::vector<long double> local = FELevelSetMoments::GraphCellIntegration::SimplexPolynomial(parent,regions[n]);
    long double minimum;
    long double maximum;
    long double padding;
    FELevelSetMoments::GraphCellIntegration::Range(parent,regions[n],4,minimum,maximum,padding);
    for(unsigned sample = 0; sample < 7; ++sample) {
      const long double u = (sample+1.0L)/32.0L;
      const long double v = (7.0L-sample)/32.0L;
      const long double w = 0.125L;
      std::array<long double,3> point;
      for(unsigned d = 0; d < 3; ++d) point[d] = regions[n].vertices[0][d]+u*(regions[n].vertices[1][d]-regions[n].vertices[0][d])+v*(regions[n].vertices[2][d]-regions[n].vertices[0][d])+w*(regions[n].vertices[3][d]-regions[n].vertices[0][d]);
      const long double original = FELevelSetMoments::GraphCellIntegration::Value(parent,point);
      long double mapped = 0.0L;
      for(unsigned i = 0; i <= 4; ++i) for(unsigned j = 0; j <= 4-i; ++j) for(unsigned k = 0; k <= 4-i-j; ++k) mapped += local[FELevelSetMoments::GraphCellIntegration::Index(i,j,k)]*FELevelSetMoments::ReferencePower(u,i)*FELevelSetMoments::ReferencePower(v,j)*FELevelSetMoments::ReferencePower(w,k);
      GraphCheck(original,mapped,2.0e-17L,"parent polynomial restriction");
      GraphCheck(original,std::max(minimum,std::min(maximum,original)),padding+2.0e-17L,"Bernstein enclosure");
    }
  }
  ++graphCases;
}

void GraphIntegrate(bool tet, const std::vector<long double>& coefficient, unsigned degree, double tolerance, bool surface, double scaling, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, std::vector<double>& error) {
  if(tet) {
    femus::TetBiquadratic basis;
    std::vector<double> nodal(15);
    for(unsigned node = 0; node < 15; ++node) {
      const double* point = basis.GetXcoarse(node);
      nodal[node] = scaling*static_cast<double>(FELevelSetMoments::GraphCellIntegration::Value(coefficient,{point[0],point[1],point[2]}));
    }
    FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,degree,positive,delta,tolerance,tolerance,true,24,0,200000,&error,1,nullptr,surface);
  }
  else {
    std::vector<std::vector<std::vector<double>>> field(3,std::vector<std::vector<double>>(3,std::vector<double>(3)));
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) field[i][j][k] = scaling*coefficient[FELevelSetMoments::Hex27Index(i,j,k)];
    FELevelSetMoments::IntegrateHex27ReferenceSelected(field,degree,positive,delta,tolerance,tolerance,true,24,0,200000,&error,1,true,nullptr,nullptr,surface);
  }
}

long double GraphPlaneH(long double cut, unsigned fixed, unsigned first, unsigned second) {
  long double value = 0.0L;
  for(unsigned r = 0; r <= fixed; ++r) value += GraphFactorial(fixed)/GraphFactorial(fixed-r)*FELevelSetMoments::ReferencePower(cut,fixed-r)*FELevelSetMoments::ReferencePower(1.0L-cut,r+first+second+3U)*GraphFactorial(first)*GraphFactorial(second)/GraphFactorial(r+first+second+3U);
  return value;
}

void GraphPlanes(bool tet, unsigned degree) {
  const std::vector<long double> cuts = tet ? std::vector<long double>{0.0L,0.25L,1.0L-1.0L/256.0L} : std::vector<long double>{-1.0L,0.0L,0.25L,1.0L-1.0L/256.0L,1.0L};
  for(unsigned axis = 0; axis < 3; ++axis) for(unsigned c = 0; c < cuts.size(); ++c) {
    std::vector<long double> coefficient(27,0.0L);
    coefficient[0] = -cuts[c];
    unsigned power[3] = {0,0,0};
    power[axis] = 1;
    coefficient[FELevelSetMoments::Hex27Index(power[0],power[1],power[2])] = 1.0L;
    std::vector<std::vector<std::vector<double>>> h;
    std::vector<std::vector<std::vector<double>>> d;
    std::vector<double> error;
    GraphIntegrate(tet,coefficient,degree,1.0e-11,false,1.0,h,d,error);
    unsigned first;
    unsigned second;
    FELevelSetMoments::Hex27OtherAxes(axis,first,second);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) {
      const unsigned powers[3] = {i,j,k};
      const long double expectedH = tet ? GraphPlaneH(cuts[c],powers[axis],powers[first],powers[second]) : FELevelSetMoments::GraphCellIntegration::IntegralPower(cuts[c],1.0L,powers[axis])*FELevelSetMoments::GraphCellIntegration::IntegralPower(-1.0L,1.0L,powers[first])*FELevelSetMoments::GraphCellIntegration::IntegralPower(-1.0L,1.0L,powers[second]);
      long double expectedD = FELevelSetMoments::ReferencePower(cuts[c],powers[axis]);
      if(tet) expectedD *= FELevelSetMoments::ReferencePower(1.0L-cuts[c],powers[first]+powers[second]+2U)*GraphFactorial(powers[first])*GraphFactorial(powers[second])/GraphFactorial(powers[first]+powers[second]+2U);
      else expectedD *= FELevelSetMoments::GraphCellIntegration::IntegralPower(-1.0L,1.0L,powers[first])*FELevelSetMoments::GraphCellIntegration::IntegralPower(-1.0L,1.0L,powers[second]);
      if(cuts[c] == (tet ? 0.0L : -1.0L) || cuts[c] == 1.0L) expectedD *= 0.5L;
      GraphCheck(h[i][j][k],expectedH,4.0e-11L,"plane H");
      GraphCheck(d[i][j][k],expectedD,4.0e-11L,"plane delta");
    }
    ++graphCases;
  }
}

long double GraphCenteredBallMoment(unsigned i, unsigned j, unsigned k, long double radius) {
  if((i|j|k)%2U != 0 || i%2U != 0 || j%2U != 0 || k%2U != 0) return 0.0L;
  return 2.0L*FELevelSetMoments::ReferencePower(radius,i+j+k+3U)*std::tgamma(0.5L*(i+1U))*std::tgamma(0.5L*(j+1U))*std::tgamma(0.5L*(k+1U))/((i+j+k+3U)*std::tgamma(0.5L*(i+j+k+3U)));
}

long double GraphBallMoment(unsigned i, unsigned j, unsigned k, const std::array<long double,3>& center, long double radius, bool surface, bool delta) {
  long double value = 0.0L;
  for(unsigned a = 0; a <= i; a += 2) for(unsigned b = 0; b <= j; b += 2) for(unsigned c = 0; c <= k; c += 2) {
    long double term = FELevelSetMoments::ReferenceBinomial(i,a)*FELevelSetMoments::ReferenceBinomial(j,b)*FELevelSetMoments::ReferenceBinomial(k,c)*FELevelSetMoments::ReferencePower(center[0],i-a)*FELevelSetMoments::ReferencePower(center[1],j-b)*FELevelSetMoments::ReferencePower(center[2],k-c)*GraphCenteredBallMoment(a,b,c,radius);
    if(delta) term *= (a+b+c+3U)/(2.0L*radius*radius);
    if(delta && surface) term *= 2.0L*radius;
    value += term;
  }
  return value;
}

void GraphBall(bool tet, unsigned degree, const std::array<long double,3>& center, long double radius, double tolerance, bool surface, double scale) {
  std::vector<long double> coefficient(27,0.0L);
  coefficient[0] = radius*radius;
  for(unsigned axis = 0; axis < 3; ++axis) {
    coefficient[0] -= center[axis]*center[axis];
    unsigned power[3] = {0,0,0};
    power[axis] = 1;
    coefficient[FELevelSetMoments::Hex27Index(power[0],power[1],power[2])] = 2.0L*center[axis];
    power[axis] = 2;
    coefficient[FELevelSetMoments::Hex27Index(power[0],power[1],power[2])] = -1.0L;
  }
  std::vector<std::vector<std::vector<double>>> h;
  std::vector<std::vector<std::vector<double>>> d;
  std::vector<double> error;
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  GraphIntegrate(tet,coefficient,degree,tolerance,surface,scale,h,d,error);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) {
    GraphCheck(h[i][j][k],GraphBallMoment(i,j,k,center,radius,false,false),8.0L*tolerance,"ball H");
    GraphCheck(d[i][j][k],GraphBallMoment(i,j,k,center,radius,surface,true)/(surface ? 1.0L : std::abs(scale)),8.0L*tolerance,"ball surface/delta");
  }
  const double milliseconds = 1000.0*std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  std::cout << (tet ? "TET" : "HEX") << " ball degree=" << degree << " radius=" << radius << " surface=" << surface << " scale=" << scale << " tolerance=" << tolerance << " milliseconds=" << milliseconds << " samples=" << error[3] << " estimated=" << error[0] << ',' << error[1] << std::endl;
  ++graphCases;
}

void GraphBoundaryPolynomial() {
  std::vector<double> nodal(15);
  for(unsigned node = 0; node < 15; ++node) nodal[node] = std::sin(0.7*(node+1U));
  std::vector<long double> original;
  FELevelSetMoments::Tet15PolynomialFromNodalValues(nodal,false,original);
  std::vector<long double> coefficient(27);
  for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) coefficient[FELevelSetMoments::Hex27Index(i,j,k)] = original[FELevelSetMoments::Tet15Index(i,j,k)];
  for(unsigned outer = 0; outer < 3; ++outer) for(unsigned inner = 0; inner < 3; ++inner) if(inner != outer) {
    const std::vector<std::vector<long double>> face = FELevelSetMoments::GraphCellIntegration::TetrahedronBoundary(coefficient,outer,inner);
    for(unsigned i = 0; i <= 8; ++i) for(unsigned j = 0; j <= 8-i; ++j) {
      std::array<long double,3> point;
      point[outer] = i/8.0L;
      point[inner] = j/8.0L;
      point[3U-outer-inner] = 1.0L-point[outer]-point[inner];
      long double value = 0.0L;
      for(unsigned first = 0; first < 3; ++first) for(unsigned second = 0; second < 3; ++second) value += face[first][second]*FELevelSetMoments::ReferencePower(point[outer],first)*FELevelSetMoments::ReferencePower(point[inner],second);
      GraphCheck(value,FELevelSetMoments::GraphCellIntegration::Value(coefficient,point),2.0e-13L,"TET15 original face restriction");
    }
  }
  ++graphCases;
}

void GraphObliquePlane(unsigned degree) {
  for(unsigned boundary = 0; boundary < 2; ++boundary) {
    const long double cut = boundary == 0 ? 0.5L : 1.0L;
    std::vector<long double> coefficient(27,0.0L);
    coefficient[0] = cut;
    coefficient[FELevelSetMoments::Hex27Index(1,0,0)] = -1.0L;
    coefficient[FELevelSetMoments::Hex27Index(0,1,0)] = -1.0L;
    coefficient[FELevelSetMoments::Hex27Index(0,0,1)] = -1.0L;
    std::vector<std::vector<std::vector<double>>> h;
    std::vector<std::vector<std::vector<double>>> d;
    std::vector<double> error;
    femus::TetBiquadratic basis;
    std::vector<double> nodal(15);
    for(unsigned node = 0; node < 15; ++node) {
      const double* point = basis.GetXcoarse(node);
      nodal[node] = node == 12 ? static_cast<double>(cut-1.0L) : static_cast<double>(FELevelSetMoments::GraphCellIntegration::Value(coefficient,{point[0],point[1],point[2]}));
    }
    FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,degree,h,d,1.0e-11,1.0e-11,true,24,0,200000,&error);
    std::cout << "TET oblique cut=" << cut << " degree=" << degree << " passed" << std::endl;
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) {
      const long double numerator = GraphFactorial(i)*GraphFactorial(j)*GraphFactorial(k);
      const long double expectedH = FELevelSetMoments::ReferencePower(cut,i+j+k+3U)*numerator/GraphFactorial(i+j+k+3U);
      const long double expectedD = (boundary == 0 ? 1.0L : 0.5L)*FELevelSetMoments::ReferencePower(cut,i+j+k+2U)*numerator/GraphFactorial(i+j+k+2U);
      GraphCheck(h[i][j][k],expectedH,4.0e-11L,"oblique plane H");
      GraphCheck(d[i][j][k],expectedD,4.0e-11L,"oblique plane delta");
    }
    ++graphCases;
  }
}

void GraphScales(bool tet, unsigned degree) {
  std::vector<long double> coefficient(27,0.0L);
  coefficient[0] = -0.2L;
  coefficient[FELevelSetMoments::Hex27Index(0,0,1)] = 1.0L;
  coefficient[FELevelSetMoments::Hex27Index(2,0,0)] = -0.08L;
  if(tet) {
    coefficient[FELevelSetMoments::Hex27Index(1,1,1)] = 0.768L;
    coefficient[FELevelSetMoments::Hex27Index(2,1,1)] = -0.768L;
    coefficient[FELevelSetMoments::Hex27Index(1,2,1)] = -0.768L;
    coefficient[FELevelSetMoments::Hex27Index(1,1,2)] = -0.768L;
  }
  else coefficient[FELevelSetMoments::Hex27Index(2,2,2)] = 0.07L;
  const double scales[3] = {1.0/1048576.0,1048576.0,-2.0};
  for(unsigned surface = 0; surface < 2; ++surface) {
    std::vector<std::vector<std::vector<double>>> baselineH;
    std::vector<std::vector<std::vector<double>>> baselineD;
    std::vector<double> error;
    GraphIntegrate(tet,coefficient,degree,1.0e-10,surface,1.0,baselineH,baselineD,error);
    for(unsigned scaling = 0; scaling < 3; ++scaling) {
      std::vector<std::vector<std::vector<double>>> h;
      std::vector<std::vector<std::vector<double>>> d;
      GraphIntegrate(tet,coefficient,degree,1.0e-10,surface,scales[scaling],h,d,error);
      for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) {
        const long double full = tet ? GraphFactorial(i)*GraphFactorial(j)*GraphFactorial(k)/GraphFactorial(i+j+k+3U) : FELevelSetMoments::GraphCellIntegration::IntegralPower(-1.0L,1.0L,i)*FELevelSetMoments::GraphCellIntegration::IntegralPower(-1.0L,1.0L,j)*FELevelSetMoments::GraphCellIntegration::IntegralPower(-1.0L,1.0L,k);
        const long double expectedH = scales[scaling] > 0.0 ? baselineH[i][j][k] : full-baselineH[i][j][k];
        const long double expectedD = baselineD[i][j][k]/(surface ? 1.0L : std::abs(scales[scaling]));
        GraphCheck(h[i][j][k],expectedH,8.0e-10L,"scaled H and complement");
        GraphCheck(d[i][j][k],expectedD,8.0e-10L*std::max(1.0L,std::abs(expectedD)),"scaled surface/delta");
      }
      ++graphCases;
    }
  }
}

int main(int argc, char** argv) {
  const unsigned degree = argc > 1 ? std::strtoul(argv[1],nullptr,10) : 6;
  if(degree > 10) return 2;
  std::cout << std::setprecision(16);
  try {
    GraphGeometry();
    GraphBoundaryPolynomial();
    GraphObliquePlane(degree);
    GraphPlanes(false,degree);
    GraphPlanes(true,degree);
    GraphScales(false,degree);
    GraphScales(true,degree);
    GraphBall(false,degree,{0,0,0},0.5L,1.0e-10,false,1.0);
    GraphBall(false,degree,{0.125L,-0.0625L,0.25L},0.375L,1.0e-10,true,1.0);
    GraphBall(true,degree,{0.25L,0.25L,0.25L},0.125L,1.0e-10,false,1.0);
    GraphBall(true,degree,{0.25L,0.25L,0.25L},0.125L,1.0e-10,true,1.0);
    const FELevelSetMoments::GraphCellStatistics& statistics = FELevelSetMoments::GraphCellIntegration::Statistics();
    if(statistics.quadCalls == 0 || statistics.triCalls == 0 || statistics.splitCells == 0 || statistics.fallbackCells != 0) throw std::runtime_error("The 3D graph hierarchy was not exercised successfully.");
    const FELevelSetMoments::LineIntegrationStatistics& line = FELevelSetMoments::LineDomainIntegration::Statistics();
    std::cout << "line_fast=" << line.fastSamples << " line_full=" << line.fullSamples << " taylor_panels=" << line.taylorPanels << " bounded_regions=" << statistics.boundedCells << " relaxed_slices=" << statistics.relaxedSlices << std::endl;
    std::cout << "checks=" << graphChecks << " cases=" << graphCases << " max_error=" << graphMaximumError << " splits=" << statistics.splitCells << " graphs=" << statistics.graphCells << " weak=" << statistics.weakGraphs << " quad=" << statistics.quadCalls << " tri=" << statistics.triCalls << " accepted_components=" << statistics.acceptedComponents << " skipped_components=" << statistics.skippedComponents << '\n';
  }
  catch(const std::exception& exception) { std::cerr << exception.what() << '\n'; return 1; }
}
