#include "Quad9DirectionSelection.hpp"
#include "Tri7DirectionSelection.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

struct RefinementCell {
  std::vector<double> values;
  long double determinant;
};

long double RefinementLevelSet(long double x, long double y) {
  return 0.613L*0.613L-(x-0.071L)*(x-0.071L)-(y+0.043L)*(y+0.043L);
}

std::vector<RefinementCell> RefinementCells(unsigned divisions, bool triangle) {
  const long double h = 2.0L/divisions;
  const long double quad[9][2] = {{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0},{0,0}};
  const long double tri[7][2] = {{0,0},{1,0},{0,1},{0.5L,0},{0.5L,0.5L},{0,0.5L},{1.0L/3.0L,1.0L/3.0L}};
  std::vector<RefinementCell> cells;
  cells.reserve((triangle ? 2U : 1U)*divisions*divisions);
  for(unsigned i = 0; i < divisions; ++i) for(unsigned j = 0; j < divisions; ++j) for(unsigned side = 0; side < (triangle ? 2U : 1U); ++side) {
    RefinementCell cell;
    cell.values.resize(triangle ? 7U : 9U);
    cell.determinant = triangle ? h*h : 0.25L*h*h;
    for(unsigned node = 0; node < cell.values.size(); ++node) {
      const long double x = -1.0L+h*i+(triangle ? h*(side == 0 ? tri[node][0] : 1.0L-tri[node][0]) : 0.5L*h*(1.0L+quad[node][0]));
      const long double y = -1.0L+h*j+(triangle ? h*(side == 0 ? tri[node][1] : 1.0L-tri[node][1]) : 0.5L*h*(1.0L+quad[node][1]));
      cell.values[node] = static_cast<double>(RefinementLevelSet(x,y));
    }
    cells.push_back(cell);
  }
  return cells;
}

int main() {
  unsigned failures = 0;
  const unsigned degrees[3] = {2U,6U,10U};
  const long double exactArea = std::acos(-1.0L)*0.613L*0.613L;
  const long double exactLength = 2.0L*std::acos(-1.0L)*0.613L;
  std::cout << std::setprecision(18) << "element,divisions,cells,degree,area_error,length_error,total_seconds,status\n";
  for(unsigned shape = 0; shape < 2; ++shape) for(unsigned divisions = 4; divisions <= 64; divisions *= 2U) {
    const bool triangle = shape != 0;
    const std::vector<RefinementCell> cells = RefinementCells(divisions,triangle);
    for(unsigned degreeIndex = 0; degreeIndex < 3; ++degreeIndex) {
      const unsigned degree = degrees[degreeIndex];
      try {
        long double area = 0.0L;
        long double length = 0.0L;
        const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
        for(unsigned cell = 0; cell < cells.size(); ++cell) {
          std::vector<std::vector<double>> positive;
          std::vector<std::vector<double>> delta;
          if(triangle) FELevelSetMoments::IntegrateTri7NodalSelected(cells[cell].values,degree,positive,delta,1.0e-12,1.0e-12,true,false,1,0,nullptr,true);
          else FELevelSetMoments::IntegrateQuad9NodalSelected(cells[cell].values,degree,positive,delta,1.0e-12,1.0e-12,true,false,1,0,nullptr,true);
          area += cells[cell].determinant*positive[0][0];
          length += 2.0L*0.613L*cells[cell].determinant*delta[0][0];
        }
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        const long double areaError = std::abs(area-exactArea);
        const long double lengthError = std::abs(length-exactLength);
        const bool success = std::isfinite(area) && std::isfinite(length) && areaError <= 1.0e-11L && lengthError <= 1.0e-11L;
        if(!success) ++failures;
        std::cout << (triangle ? "TRI" : "QUAD") << ',' << divisions << ',' << cells.size() << ',' << degree << ',' << areaError << ',' << lengthError << ',' << seconds << ',' << (success ? "PASS" : "FAIL") << '\n';
      }
      catch(const std::exception& exception) {
        ++failures;
        std::cerr << (triangle ? "TRI" : "QUAD") << " divisions=" << divisions << " degree=" << degree << ": " << exception.what() << '\n';
      }
    }
  }
  return failures == 0 ? 0 : 1;
}
