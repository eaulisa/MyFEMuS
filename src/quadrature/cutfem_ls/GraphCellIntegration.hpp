#ifndef FEMUS_GRAPH_CELL_INTEGRATION_HPP
#define FEMUS_GRAPH_CELL_INTEGRATION_HPP

#include "Quad9DirectionSelection.hpp"
#include "TriangleIntegration.hpp"
#include "Hex27LevelSetMoments.hpp"
#include <array>
#include <map>
#include <set>

namespace FELevelSetMoments {

struct GraphCellStatistics {
  unsigned long long splitCells = 0;
  unsigned long long graphCells = 0;
  unsigned long long weakGraphs = 0;
  unsigned long long fallbackCells = 0;
  unsigned long long quadCalls = 0;
  unsigned long long triCalls = 0;
  unsigned long long outerSamples = 0;
  unsigned long long eventBuilds = 0;
  unsigned long long acceptedComponents = 0;
  unsigned long long skippedComponents = 0;
  unsigned long long relaxedSlices = 0;
  unsigned long long boundedCells = 0;
  unsigned directions[3] = {0,0,0};
};

struct GraphCellRegion {
  bool tetrahedron = false;
  std::vector<std::array<long double,3>> vertices;
  std::vector<std::vector<unsigned>> faces;
  unsigned depth = 0;
};

struct GraphCellData {
  std::vector<long double> coefficient;
  std::vector<std::vector<long double>> gradient;
  std::vector<std::array<unsigned,3>> powers;
  unsigned degree = 0;
  unsigned maximumDepth = 24;
  unsigned maximumEvaluations = 200000;
  unsigned evaluations = 0;
  unsigned accepted = 0;
  unsigned graphAxis = 0;
  unsigned preferredInner = 3;
  unsigned preferredOuter = 3;
  bool unitCube = false;
  long double domainVolume = 1.0L;
  std::vector<unsigned char> requested;
  std::array<unsigned,3> innerCounts = {0,0,0};
  std::array<unsigned,3> outerCounts = {0,0,0};
  bool computeDelta = true;
  long double absoluteH = 0.0L;
  long double absoluteD = 0.0L;
  long double relative = 0.0L;
};

class GraphCellIntegration {
public:
  // Partitions the parent polynomial geometrically and reuses the native QUAD/TRI integrators on every slice.
  static GraphCellStatistics& Statistics() {
    static thread_local GraphCellStatistics data;
    return data;
  }

  static long double Determinant(const std::array<long double,3>& first, const std::array<long double,3>& second, const std::array<long double,3>& third) {
    return std::fma(first[0],std::fma(second[1],third[2],-second[2]*third[1]),std::fma(-first[1],std::fma(second[0],third[2],-second[2]*third[0]),first[2]*std::fma(second[0],third[1],-second[1]*third[0])));
  }

  static std::array<long double,3> Center(const GraphCellRegion& cell) {
    std::array<long double,3> center = {0,0,0};
    for(unsigned vertex = 0; vertex < cell.vertices.size(); ++vertex) for(unsigned axis = 0; axis < 3; ++axis) center[axis] += cell.vertices[vertex][axis]/cell.vertices.size();
    return center;
  }

  static long double Volume(const GraphCellRegion& cell) {
    if(!cell.tetrahedron) return (cell.vertices[1][0]-cell.vertices[0][0])*(cell.vertices[1][1]-cell.vertices[0][1])*(cell.vertices[1][2]-cell.vertices[0][2]);
    if(!cell.faces.empty()) {
      const std::array<long double,3> center = Center(cell);
      long double volume = 0.0L;
      for(unsigned face = 0; face < cell.faces.size(); ++face) for(unsigned triangle = 1; triangle+1U < cell.faces[face].size(); ++triangle) {
        std::array<long double,3> edges[3];
        const unsigned vertices[3] = {cell.faces[face][0],cell.faces[face][triangle],cell.faces[face][triangle+1U]};
        for(unsigned vertex = 0; vertex < 3; ++vertex) for(unsigned axis = 0; axis < 3; ++axis) edges[vertex][axis] = cell.vertices[vertices[vertex]][axis]-center[axis];
        volume += std::abs(Determinant(edges[0],edges[1],edges[2]))/6.0L;
      }
      return volume;
    }
    std::array<long double,3> edges[3];
    for(unsigned i = 0; i < 3; ++i) for(unsigned d = 0; d < 3; ++d) edges[i][d] = cell.vertices[i+1U][d]-cell.vertices[0][d];
    return std::abs(Determinant(edges[0],edges[1],edges[2]))/6.0L;
  }

  static void Bounds(const GraphCellRegion& cell, std::array<long double,3>& lower, std::array<long double,3>& upper) {
    lower = upper = cell.vertices[0];
    for(unsigned i = 1; i < cell.vertices.size(); ++i) for(unsigned d = 0; d < 3; ++d) {
      lower[d] = std::min(lower[d],cell.vertices[i][d]);
      upper[d] = std::max(upper[d],cell.vertices[i][d]);
    }
  }

  static long double Value(const std::vector<long double>& coefficient, const std::array<long double,3>& point) {
    long double value = 0.0L;
    for(unsigned i = 3; i > 0; --i) {
      long double row = 0.0L;
      for(unsigned j = 3; j > 0; --j) {
        long double line = 0.0L;
        for(unsigned k = 3; k > 0; --k) line = std::fma(line,point[2],coefficient[Hex27Index(i-1U,j-1U,k-1U)]);
        row = std::fma(row,point[1],line);
      }
      value = std::fma(value,point[0],row);
    }
    return value;
  }

  static std::vector<long double> Derivative(const std::vector<long double>& coefficient, unsigned axis) {
    std::vector<long double> result(27,0.0L);
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) {
      unsigned powers[3] = {i,j,k};
      if(powers[axis] == 0) continue;
      const long double value = powers[axis]*coefficient[Hex27Index(i,j,k)];
      --powers[axis];
      result[Hex27Index(powers[0],powers[1],powers[2])] = value;
    }
    return result;
  }

  static std::vector<long double> BoxPolynomial(const std::vector<long double>& coefficient, const std::array<long double,3>& lower, const std::array<long double,3>& upper) {
    std::vector<long double> result = coefficient;
    for(unsigned axis = 0; axis < 3; ++axis) {
      std::vector<long double> next(27,0.0L);
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) {
        const unsigned powers[3] = {i,j,k};
        unsigned target[3] = {i,j,k};
        for(unsigned a = 0; a <= powers[axis]; ++a) {
          target[axis] = a;
          const long double multiplier = ReferenceBinomial(powers[axis],a)*ReferencePower(lower[axis],powers[axis]-a)*ReferencePower(upper[axis]-lower[axis],a);
          next[Hex27Index(target[0],target[1],target[2])] = std::fma(multiplier,result[Hex27Index(i,j,k)],next[Hex27Index(target[0],target[1],target[2])]);
        }
      }
      result.swap(next);
    }
    return result;
  }

  static unsigned Index(unsigned i, unsigned j, unsigned k) {
    return 25U*i+5U*j+k;
  }

  static std::vector<long double> SimplexPolynomial(const std::vector<long double>& coefficient, const GraphCellRegion& cell) {
    std::vector<long double> result(125,0.0L);
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) {
      const long double multiplier = coefficient[Hex27Index(i,j,k)];
      if(multiplier == 0.0L) continue;
      if(i+j+k > 4) throw std::logic_error("The tetrahedral parent polynomial exceeds total degree four.");
      const unsigned powers[3] = {i,j,k};
      std::vector<long double> product(125,0.0L);
      product[0] = multiplier;
      unsigned degree = 0;
      for(unsigned axis = 0; axis < 3; ++axis) for(unsigned repetition = 0; repetition < powers[axis]; ++repetition) {
        std::vector<long double> next(125,0.0L);
        for(unsigned a = 0; a <= degree; ++a) for(unsigned b = 0; b <= degree-a; ++b) for(unsigned c = 0; c <= degree-a-b; ++c) {
          const long double value = product[Index(a,b,c)];
          next[Index(a,b,c)] = std::fma(cell.vertices[0][axis],value,next[Index(a,b,c)]);
          next[Index(a+1U,b,c)] = std::fma(cell.vertices[1][axis]-cell.vertices[0][axis],value,next[Index(a+1U,b,c)]);
          next[Index(a,b+1U,c)] = std::fma(cell.vertices[2][axis]-cell.vertices[0][axis],value,next[Index(a,b+1U,c)]);
          next[Index(a,b,c+1U)] = std::fma(cell.vertices[3][axis]-cell.vertices[0][axis],value,next[Index(a,b,c+1U)]);
        }
        product.swap(next);
        ++degree;
      }
      for(unsigned n = 0; n < 125; ++n) result[n] += product[n];
    }
    return result;
  }

  static void Range(const std::vector<long double>& coefficient, const GraphCellRegion& cell, unsigned degree, long double& minimum, long double& maximum, long double& padding) {
    std::vector<long double> controls;
    long double norm = 0.0L;
    if(cell.tetrahedron && cell.faces.empty()) {
      const std::vector<long double> power = SimplexPolynomial(coefficient,cell);
      for(unsigned n = 0; n < power.size(); ++n) norm += std::abs(power[n]);
      for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) {
        long double value = 0.0L;
        for(unsigned a = 0; a <= i; ++a) for(unsigned b = 0; b <= j; ++b) for(unsigned c = 0; c <= k; ++c) value += power[Index(a,b,c)]*ReferenceBinomial(i,a)*ReferenceBinomial(j,b)*ReferenceBinomial(k,c)/(ReferenceBinomial(degree,a+b+c)*ReferenceBinomial(a+b+c,a)*ReferenceBinomial(b+c,b));
        controls.push_back(value);
      }
    }
    else {
      std::array<long double,3> lower;
      std::array<long double,3> upper;
      Bounds(cell,lower,upper);
      const std::vector<long double> power = BoxPolynomial(coefficient,lower,upper);
      for(unsigned n = 0; n < power.size(); ++n) norm += std::abs(power[n]);
      controls = Hex27Bernstein(power);
    }
    minimum = *std::min_element(controls.begin(),controls.end());
    maximum = *std::max_element(controls.begin(),controls.end());
    long double parentNorm = 0.0L;
    for(unsigned n = 0; n < coefficient.size(); ++n) parentNorm += std::abs(coefficient[n]);
    padding = 512.0L*std::numeric_limits<long double>::epsilon()*std::max(norm,parentNorm);
  }

  static int Phase(const GraphCellData& data, const GraphCellRegion& cell) {
    long double minimum;
    long double maximum;
    long double padding;
    Range(data.coefficient,cell,4,minimum,maximum,padding);
    if(minimum > padding || (!data.computeDelta && minimum >= 0.0L && maximum > 0.0L)) return 1;
    if(maximum < -padding || (!data.computeDelta && maximum <= 0.0L && minimum < 0.0L)) return -1;
    return 0;
  }

  static std::array<long double,3> Intersect(const std::array<long double,3>& first, const std::array<long double,3>& second, unsigned axis, long double coordinate) {
    if(coordinate == first[axis]) return first;
    if(coordinate == second[axis]) return second;
    const long double ratio = (coordinate-first[axis])/(second[axis]-first[axis]);
    std::array<long double,3> point;
    for(unsigned d = 0; d < 3; ++d) point[d] = std::fma(ratio,second[d]-first[d],first[d]);
    point[axis] = coordinate;
    return point;
  }

  static void Clip(const GraphCellRegion& cell, unsigned axis, long double coordinate, bool lowerSide, std::vector<GraphCellRegion>& output) {
    if(!cell.tetrahedron) {
      GraphCellRegion child = cell;
      if(lowerSide) child.vertices[1][axis] = std::min(child.vertices[1][axis],coordinate);
      else child.vertices[0][axis] = std::max(child.vertices[0][axis],coordinate);
      ++child.depth;
      if(Volume(child) > 0.0L) output.push_back(child);
      return;
    }
    const std::vector<std::vector<unsigned>> faces = cell.faces.empty() ? std::vector<std::vector<unsigned>>{{0,1,2},{0,1,3},{0,2,3},{1,2,3}} : cell.faces;
    unsigned inside = 0;
    for(unsigned vertex = 0; vertex < cell.vertices.size(); ++vertex) if(lowerSide ? cell.vertices[vertex][axis] <= coordinate : cell.vertices[vertex][axis] >= coordinate) ++inside;
    if(inside == 0) return;
    if(inside == cell.vertices.size()) {
      output.push_back(cell);
      return;
    }
    GraphCellRegion child;
    child.tetrahedron = true;
    child.depth = cell.depth+1U;
    std::map<std::pair<unsigned,unsigned>,std::array<long double,3>> intersections;
    std::vector<unsigned> cut;
    for(unsigned face = 0; face < faces.size(); ++face) {
      std::vector<unsigned> polygon;
      for(unsigned edge = 0; edge < faces[face].size(); ++edge) {
        const unsigned first = faces[face][edge];
        const unsigned second = faces[face][(edge+1U)%faces[face].size()];
        const bool firstInside = lowerSide ? cell.vertices[first][axis] <= coordinate : cell.vertices[first][axis] >= coordinate;
        const bool secondInside = lowerSide ? cell.vertices[second][axis] <= coordinate : cell.vertices[second][axis] >= coordinate;
        if(firstInside) polygon.push_back(VertexIndex(child,cell.vertices[first]));
        if(firstInside != secondInside) {
          const std::pair<unsigned,unsigned> key = std::minmax(first,second);
          if(intersections.find(key) == intersections.end()) intersections[key] = Intersect(cell.vertices[key.first],cell.vertices[key.second],axis,coordinate);
          const unsigned vertex = VertexIndex(child,intersections[key]);
          polygon.push_back(vertex);
          cut.push_back(vertex);
        }
        if(firstInside && cell.vertices[first][axis] == coordinate) cut.push_back(VertexIndex(child,cell.vertices[first]));
      }
      polygon.erase(std::unique(polygon.begin(),polygon.end()),polygon.end());
      if(polygon.size() > 1 && polygon.front() == polygon.back()) polygon.pop_back();
      if(polygon.size() >= 3) child.faces.push_back(polygon);
    }
    std::sort(cut.begin(),cut.end());
    cut.erase(std::unique(cut.begin(),cut.end()),cut.end());
    if(cut.size() >= 3) {
      unsigned first;
      unsigned second;
      Hex27OtherAxes(axis,first,second);
      long double cx = 0.0L;
      long double cy = 0.0L;
      for(unsigned n = 0; n < cut.size(); ++n) {
        cx += child.vertices[cut[n]][first]/cut.size();
        cy += child.vertices[cut[n]][second]/cut.size();
      }
      for(unsigned i = 0; i < cut.size(); ++i) for(unsigned j = i+1U; j < cut.size(); ++j) if(std::atan2(child.vertices[cut[j]][second]-cy,child.vertices[cut[j]][first]-cx) < std::atan2(child.vertices[cut[i]][second]-cy,child.vertices[cut[i]][first]-cx)) std::swap(cut[i],cut[j]);
      child.faces.push_back(cut);
    }
    if(child.vertices.size() >= 4 && Volume(child) > 0.0L) output.push_back(child);
  }

  static unsigned VertexIndex(GraphCellRegion& cell, const std::array<long double,3>& point) {
    for(unsigned vertex = 0; vertex < cell.vertices.size(); ++vertex) if(cell.vertices[vertex] == point) return vertex;
    cell.vertices.push_back(point);
    return cell.vertices.size()-1U;
  }

  static std::vector<GraphCellRegion> Slab(const GraphCellRegion& cell, unsigned axis, long double left, long double right) {
    std::vector<GraphCellRegion> first;
    std::vector<GraphCellRegion> result;
    Clip(cell,axis,left,false,first);
    for(unsigned i = 0; i < first.size(); ++i) Clip(first[i],axis,right,true,result);
    return result;
  }

  static long double IntegralPower(long double left, long double right, unsigned degree) {
    long double result = 0.0L;
    long double power = 1.0L;
    for(unsigned i = 0; i <= degree; ++i) {
      result = power+left*result;
      power *= right;
    }
    return (right-left)*result/(degree+1U);
  }

  static long double TetMoment(const GraphCellRegion& cell, unsigned first, unsigned second, unsigned third) {
    if(!cell.faces.empty()) {
      const std::array<long double,3> center = Center(cell);
      long double value = 0.0L;
      for(unsigned face = 0; face < cell.faces.size(); ++face) for(unsigned triangle = 1; triangle+1U < cell.faces[face].size(); ++triangle) {
        GraphCellRegion tetrahedron;
        tetrahedron.tetrahedron = true;
        tetrahedron.vertices = {center,cell.vertices[cell.faces[face][0]],cell.vertices[cell.faces[face][triangle]],cell.vertices[cell.faces[face][triangle+1U]]};
        value += TetMoment(tetrahedron,first,second,third);
      }
      return value;
    }
    bool standard = true;
    for(unsigned vertex = 0; vertex < 4; ++vertex) for(unsigned axis = 0; axis < 3; ++axis) standard = standard && cell.vertices[vertex][axis] == (vertex == axis+1U ? 1.0L : 0.0L);
    if(standard) {
      long double value = 1.0L;
      for(unsigned n = 2; n <= first; ++n) value *= n;
      for(unsigned n = 2; n <= second; ++n) value *= n;
      for(unsigned n = 2; n <= third; ++n) value *= n;
      for(unsigned n = 2; n <= first+second+third+3U; ++n) value /= n;
      return value;
    }
    std::vector<long double> coefficients(1,1.0L);
    const unsigned powers[3] = {first,second,third};
    unsigned degree = 0;
    for(unsigned axis = 0; axis < 3; ++axis) for(unsigned repetition = 0; repetition < powers[axis]; ++repetition) {
      ++degree;
      std::vector<long double> next((degree+1U)*(degree+2U)*(degree+3U)/6U,0.0L);
      for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= degree-i-j; ++k) {
        const unsigned zero = degree-i-j-k;
        long double value = 0.0L;
        if(zero != 0) value += zero*cell.vertices[0][axis]*coefficients[HexMomentIndex(i,j,k)];
        if(i != 0) value += i*cell.vertices[1][axis]*coefficients[HexMomentIndex(i-1U,j,k)];
        if(j != 0) value += j*cell.vertices[2][axis]*coefficients[HexMomentIndex(i,j-1U,k)];
        if(k != 0) value += k*cell.vertices[3][axis]*coefficients[HexMomentIndex(i,j,k-1U)];
        next[HexMomentIndex(i,j,k)] = value/degree;
      }
      coefficients.swap(next);
    }
    long double sum = 0.0L;
    for(unsigned i = 0; i < coefficients.size(); ++i) sum += coefficients[i];
    return Volume(cell)*sum/coefficients.size();
  }

  static void Full(const GraphCellData& data, const GraphCellRegion& cell, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<unsigned char>& active, long double multiplier = 1.0L) {
    const unsigned count = data.powers.size();
    for(unsigned n = 0; n < count; ++n) {
      if(!active[n]) continue;
      long double value = cell.tetrahedron ? TetMoment(cell,data.powers[n][0],data.powers[n][1],data.powers[n][2]) : IntegralPower(cell.vertices[0][0],cell.vertices[1][0],data.powers[n][0])*IntegralPower(cell.vertices[0][1],cell.vertices[1][1],data.powers[n][1])*IntegralPower(cell.vertices[0][2],cell.vertices[1][2],data.powers[n][2]);
      value *= multiplier;
      values[n] += value;
      errors[n] += (64.0L+16.0L*data.degree)*std::numeric_limits<long double>::epsilon()*std::abs(value);
    }
  }

  static std::vector<std::array<long double,3>> Corners(const GraphCellRegion& cell) {
    if(cell.tetrahedron) return cell.vertices;
    std::vector<std::array<long double,3>> result(8);
    for(unsigned n = 0; n < 8; ++n) for(unsigned d = 0; d < 3; ++d) result[n][d] = cell.vertices[(n>>d)&1U][d];
    return result;
  }

  static std::vector<std::vector<unsigned>> Edges(const GraphCellRegion& cell) {
    std::vector<std::vector<unsigned>> result;
    if(!cell.faces.empty()) {
      std::set<std::pair<unsigned,unsigned>> edges;
      for(unsigned face = 0; face < cell.faces.size(); ++face) for(unsigned edge = 0; edge < cell.faces[face].size(); ++edge) edges.insert(std::minmax(cell.faces[face][edge],cell.faces[face][(edge+1U)%cell.faces[face].size()]));
      for(std::set<std::pair<unsigned,unsigned>>::const_iterator edge = edges.begin(); edge != edges.end(); ++edge) result.push_back({edge->first,edge->second});
    }
    else if(cell.tetrahedron) {
      for(unsigned i = 0; i < 4; ++i) for(unsigned j = i+1U; j < 4; ++j) result.push_back({i,j});
    }
    else for(unsigned n = 0; n < 8; ++n) for(unsigned axis = 0; axis < 3; ++axis) if(((n>>axis)&1U) == 0) result.push_back({n,n|(1U<<axis)});
    return result;
  }

  static bool PointLess(const std::array<long double,3>& first, const std::array<long double,3>& second) {
    if(first[0] != second[0]) return first[0] < second[0];
    if(first[1] != second[1]) return first[1] < second[1];
    return first[2] < second[2];
  }

  static void SlicePolygon(const GraphCellRegion& cell, unsigned axis, long double coordinate, std::vector<TriangleDomainVertex>& polygon) {
    const std::vector<std::array<long double,3>> vertices = Corners(cell);
    const std::vector<std::vector<unsigned>> edges = Edges(cell);
    std::vector<std::array<long double,3>> intersections;
    for(unsigned n = 0; n < edges.size(); ++n) {
      const std::array<long double,3> first = vertices[edges[n][0]];
      const std::array<long double,3> second = vertices[edges[n][1]];
      if(first[axis] == coordinate) intersections.push_back(first);
      if(second[axis] == coordinate) intersections.push_back(second);
      if((first[axis] < coordinate && second[axis] > coordinate) || (first[axis] > coordinate && second[axis] < coordinate)) intersections.push_back(Intersect(first,second,axis,coordinate));
    }
    std::sort(intersections.begin(),intersections.end(),PointLess);
    intersections.erase(std::unique(intersections.begin(),intersections.end()),intersections.end());
    polygon.clear();
    unsigned first;
    unsigned second;
    Hex27OtherAxes(axis,first,second);
    for(unsigned n = 0; n < intersections.size(); ++n) polygon.push_back({intersections[n][first],intersections[n][second]});
    if(polygon.size() < 3) return;
    long double cx = 0.0L;
    long double cy = 0.0L;
    for(unsigned n = 0; n < polygon.size(); ++n) {
      cx += polygon[n].x/polygon.size();
      cy += polygon[n].y/polygon.size();
    }
    for(unsigned i = 0; i < polygon.size(); ++i) for(unsigned j = i+1U; j < polygon.size(); ++j) if(std::atan2(polygon[j].y-cy,polygon[j].x-cx) < std::atan2(polygon[i].y-cy,polygon[i].x-cx)) std::swap(polygon[i],polygon[j]);
  }

  static std::vector<std::vector<long double>> SlicePolynomial(const std::vector<long double>& coefficient, unsigned axis, long double coordinate) {
    const unsigned size = SurfacePolynomialCubeSize(coefficient.size());
    unsigned first;
    unsigned second;
    Hex27OtherAxes(axis,first,second);
    std::vector<std::vector<long double>> result(size,std::vector<long double>(size,0.0L));
    for(unsigned i = 0; i < size; ++i) for(unsigned j = 0; j < size; ++j) {
      unsigned powers[3] = {0,0,0};
      powers[first] = i;
      powers[second] = j;
      for(unsigned k = size; k > 0; --k) {
        powers[axis] = k-1U;
        result[i][j] = std::fma(result[i][j],coordinate,coefficient[size*size*powers[0]+size*powers[1]+powers[2]]);
      }
    }
    return result;
  }

  static void Sample(GraphCellData& data, const GraphCellRegion& cell, unsigned axis, long double coordinate, const std::vector<unsigned char>& active, std::vector<long double>& values, std::vector<long double>& errors, bool& cut, long double allocatedArea = -1.0L) {
    if(data.evaluations >= data.maximumEvaluations) throw std::runtime_error("The 3D graph slice budget was exhausted.");
    ++data.evaluations;
    ++Statistics().outerSamples;
    unsigned first;
    unsigned second;
    Hex27OtherAxes(axis,first,second);
    const unsigned count = data.powers.size();
    unsigned degree = 0;
    for(unsigned n = 0; n < count; ++n) if(active[n] || active[count+n]) degree = std::max(degree,data.powers[n][first]+data.powers[n][second]);
    const unsigned sliceCount = (degree+1U)*(degree+2U)/2U;
    std::vector<unsigned char> required(2U*sliceCount,0);
    for(unsigned n = 0; n < count; ++n) for(unsigned pair = 0; pair < 2; ++pair) if(active[pair*count+n]) required[pair*sliceCount+ReferenceMomentIndex2D(degree,data.powers[n][first],data.powers[n][second])] = 1;

    const std::vector<std::vector<long double>> coefficient = SlicePolynomial(data.coefficient,axis,coordinate);
    std::vector<std::vector<std::vector<long double>>> gradient;
    for(unsigned n = 0; n < data.gradient.size(); ++n) gradient.push_back(SlicePolynomial(data.gradient[n],axis,coordinate));
    required[0] = 1;
    if(data.computeDelta) required[sliceCount] = 1;
    std::vector<long double> moments(2U*sliceCount,0.0L);
    std::vector<long double> momentErrors(2U*sliceCount,0.0L);
    long double area = 0.0L;
    if(cell.tetrahedron) {
      std::vector<TriangleDomainVertex> polygon;
      SlicePolygon(cell,axis,coordinate,polygon);
      std::vector<std::vector<TriangleDomainVertex>> triangles;
      for(unsigned n = 1; n+1U < polygon.size(); ++n) {
        std::vector<TriangleDomainVertex> triangle = {polygon[0],polygon[n],polygon[n+1U]};
        const bool transpose = data.graphAxis != first;
        if(transpose) for(unsigned vertex = 0; vertex < 3; ++vertex) std::swap(triangle[vertex].x,triangle[vertex].y);
        for(unsigned i = 0; i < 3; ++i) for(unsigned j = i+1U; j < 3; ++j) if(triangle[j].y < triangle[i].y) std::swap(triangle[i],triangle[j]);
        std::vector<std::vector<TriangleDomainVertex>> pieces;
        if(triangle[1].y > triangle[0].y && triangle[1].y < triangle[2].y) {
          const long double ratio = (triangle[1].y-triangle[0].y)/(triangle[2].y-triangle[0].y);
          const TriangleDomainVertex point = {std::fma(ratio,triangle[2].x-triangle[0].x,triangle[0].x),triangle[1].y};
          pieces = {{triangle[0],triangle[1],point},{triangle[1],triangle[2],point}};
        }
        else pieces.push_back(triangle);
        for(unsigned piece = 0; piece < pieces.size(); ++piece) {
          if(transpose) for(unsigned vertex = 0; vertex < 3; ++vertex) std::swap(pieces[piece][vertex].x,pieces[piece][vertex].y);
          triangles.push_back(pieces[piece]);
        }
      }
      std::array<long double,3> lower;
      std::array<long double,3> upper;
      Bounds(cell,lower,upper);
      const long double averageArea = Volume(cell)/(upper[axis]-lower[axis]);
      for(unsigned n = 0; n < triangles.size(); ++n) {
        std::vector<TriangleDomainVertex> triangle = triangles[n];
        const long double measure = 0.5L*std::abs(std::fma(triangle[1].x-triangle[0].x,triangle[2].y-triangle[0].y,-(triangle[2].x-triangle[0].x)*(triangle[1].y-triangle[0].y)));
        if(measure == 0.0L) continue;
        area += measure;
        long double firstMinimum = triangle[0].x;
        long double firstMaximum = firstMinimum;
        long double secondMinimum = triangle[0].y;
        long double secondMaximum = secondMinimum;
        for(unsigned vertex = 1; vertex < 3; ++vertex) {
          firstMinimum = std::min(firstMinimum,triangle[vertex].x);
          firstMaximum = std::max(firstMaximum,triangle[vertex].x);
          secondMinimum = std::min(secondMinimum,triangle[vertex].y);
          secondMaximum = std::max(secondMaximum,triangle[vertex].y);
        }
        const long double firstScale = firstMaximum-firstMinimum;
        const long double secondScale = secondMaximum-secondMinimum;
        const long double jacobian = firstScale*secondScale;
        const long double allocation = 0.02L*(allocatedArea > 0.0L ? allocatedArea : averageArea)/triangles.size();
        std::vector<std::vector<long double>> chart;
        ReferencePolynomialAffine2D(coefficient,firstMinimum,firstScale,secondMinimum,secondScale,chart);
        std::vector<std::vector<std::vector<long double>>> chartGradient(gradient.size());
        for(unsigned component = 0; component < gradient.size(); ++component) ReferencePolynomialAffine2D(gradient[component],firstMinimum,firstScale,secondMinimum,secondScale,chartGradient[component]);
        for(unsigned vertex = 0; vertex < 3; ++vertex) {
          triangle[vertex].x = (triangle[vertex].x-firstMinimum)/firstScale;
          triangle[vertex].y = (triangle[vertex].y-secondMinimum)/secondScale;
        }
        std::vector<long double> native;
        std::vector<long double> nativeErrors;
        unsigned evaluations = 0;
        unsigned accepted = 0;
        unsigned selected = 0;
        ++Statistics().triCalls;
        const std::array<long double,5> parent = {firstMinimum,firstScale,secondMinimum,secondScale,jacobian};
        long double relaxation = 1.0L;
        for(unsigned retry = 0; retry < 5; ++retry) {
          try {
            TriangleIntegration::Integrate(chart,degree,triangle,data.computeDelta,relaxation*allocation*data.absoluteH,relaxation*allocation*data.absoluteD,0.02L*data.relative,std::min(data.maximumDepth,8U),data.maximumEvaluations,0,10,native,nativeErrors,evaluations,accepted,selected,gradient.empty() ? nullptr : &chartGradient,&parent,&required,data.graphAxis == first ? 2U : 1U);
            break;
          }
          catch(const std::runtime_error&) {
            if(retry == 4) throw;
            relaxation *= 4096.0L;
            evaluations = accepted = 0;
            ++Statistics().relaxedSlices;
          }
        }
        std::vector<long double> triangleMoments;
        std::vector<long double> triangleErrors;
        triangleMoments.swap(native);
        triangleErrors.swap(nativeErrors);
        for(unsigned k = 0; k < moments.size(); ++k) {
          moments[k] += triangleMoments[k];
          momentErrors[k] += triangleErrors[k]+8.0L*std::numeric_limits<long double>::epsilon()*std::abs(triangleMoments[k]);
        }
      }
    }
    else {
      const long double firstScale = 0.5L*(cell.vertices[1][first]-cell.vertices[0][first]);
      const long double secondScale = 0.5L*(cell.vertices[1][second]-cell.vertices[0][second]);
      const long double firstOrigin = cell.vertices[0][first]+firstScale;
      const long double secondOrigin = cell.vertices[0][second]+secondScale;
      area = 4.0L*firstScale*secondScale;
      std::vector<std::vector<long double>> chart;
      ReferencePolynomialAffine2D(coefficient,cell.vertices[0][first],2.0L*firstScale,cell.vertices[0][second],2.0L*secondScale,chart);
      std::vector<std::vector<std::vector<long double>>> chartGradient(gradient.size());
      for(unsigned n = 0; n < gradient.size(); ++n) ReferencePolynomialAffine2D(gradient[n],cell.vertices[0][first],2.0L*firstScale,cell.vertices[0][second],2.0L*secondScale,chartGradient[n]);
      long double lower[2];
      long double upper[2];
      const unsigned selected = Quad9SelectedOrder(chart,lower,upper);
      std::vector<long double> native;
      std::vector<long double> nativeErrors;
      unsigned evaluations = 0;
      unsigned accepted = 0;
      ++Statistics().quadCalls;
      const long double allocation = allocatedArea > 0.0L ? allocatedArea : area;
      const std::array<long double,5> parent = {cell.vertices[0][first],2.0L*firstScale,cell.vertices[0][second],2.0L*secondScale,area};
      Quad9SelectedReferenceOneDirection(chart,data.absoluteD/data.absoluteH,degree,selected,lower[2U-selected],data.computeDelta,false,0.02L*data.absoluteH*allocation,0.02L*data.relative,data.maximumDepth,data.maximumEvaluations,nullptr,native,nativeErrors,evaluations,accepted,true,gradient.empty() ? nullptr : &chartGradient,true,&parent,&required,data.graphAxis == first ? 2U : 1U);
      moments.swap(native);
      momentErrors.swap(nativeErrors);
    }
    cut = data.computeDelta ? moments[sliceCount] > 0.0L : (moments[0] > 0.0L && moments[0] < area);
    values.assign(2U*count,0.0L);
    errors.assign(2U*count,0.0L);
    for(unsigned n = 0; n < count; ++n) {
      if(!active[n] && !active[count+n]) {
        ++Statistics().skippedComponents;
        continue;
      }
      const unsigned index = ReferenceMomentIndex2D(degree,data.powers[n][first],data.powers[n][second]);
      const long double power = ReferencePower(coordinate,data.powers[n][axis]);
      for(unsigned pair = 0; pair < 2; ++pair) if(active[pair*count+n]) {
        values[pair*count+n] = power*moments[pair*sliceCount+index];
        errors[pair*count+n] = std::abs(power)*momentErrors[pair*sliceCount+index]+8.0L*std::numeric_limits<long double>::epsilon()*std::abs(values[pair*count+n]);
      }
    }
    if(!active[0]) values[0] = moments[0];
  }

  static std::vector<long double> EdgePolynomial(const std::vector<long double>& coefficient, const std::array<long double,3>& first, const std::array<long double,3>& second) {
    std::vector<long double> result(7,0.0L);
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) {
      const long double value = coefficient[Hex27Index(i,j,k)];
      if(value == 0.0L) continue;
      const unsigned powers[3] = {i,j,k};
      std::vector<long double> product(1,value);
      for(unsigned d = 0; d < 3; ++d) for(unsigned n = 0; n < powers[d]; ++n) product = LineSliceUtilities::PolynomialProduct(product,{first[d],second[d]-first[d]});
      for(unsigned n = 0; n < product.size(); ++n) result[n] += product[n];
    }
    LineSliceUtilities::Trim(result);
    return result;
  }

  static std::vector<std::vector<long double>> FacePolynomial(const std::vector<long double>& coefficient, const std::array<long double,3>& origin, const std::array<long double,3>& first, const std::array<long double,3>& second) {
    std::vector<std::vector<long double>> result(7,std::vector<long double>(7,0.0L));
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) {
      const long double value = coefficient[Hex27Index(i,j,k)];
      if(value == 0.0L) continue;
      const unsigned powers[3] = {i,j,k};
      std::vector<std::vector<long double>> product(7,std::vector<long double>(7,0.0L));
      product[0][0] = value;
      unsigned degree = 0;
      for(unsigned d = 0; d < 3; ++d) for(unsigned n = 0; n < powers[d]; ++n) {
        std::vector<std::vector<long double>> next(7,std::vector<long double>(7,0.0L));
        for(unsigned a = 0; a <= degree; ++a) for(unsigned b = 0; b <= degree-a; ++b) {
          next[a][b] = std::fma(origin[d],product[a][b],next[a][b]);
          next[a+1U][b] = std::fma(first[d],product[a][b],next[a+1U][b]);
          next[a][b+1U] = std::fma(second[d],product[a][b],next[a][b+1U]);
        }
        product.swap(next);
        ++degree;
      }
      for(unsigned a = 0; a <= degree; ++a) for(unsigned b = 0; b <= degree-a; ++b) result[a][b] += product[a][b];
    }
    Hex27TrimBivariate(result);
    return result;
  }

  static std::vector<std::vector<long double>> TetrahedronBoundary(const std::vector<long double>& coefficient, unsigned outer, unsigned inner) {
    const long double firstBubble = -coefficient[Hex27Index(0,2,1)];
    const long double secondBubble = -coefficient[Hex27Index(2,0,1)];
    const long double thirdBubble = -coefficient[Hex27Index(2,1,0)];
    const long double interiorBubble = -coefficient[Hex27Index(2,1,1)];
    std::vector<long double> face(27,0.0L);
    for(unsigned i = 0; i <= 2; ++i) for(unsigned j = 0; j <= 2-i; ++j) for(unsigned k = 0; k <= 2-i-j; ++k) face[Hex27Index(i,j,k)] = coefficient[Hex27Index(i,j,k)];
    face[Hex27Index(0,1,1)] -= firstBubble;
    face[Hex27Index(1,0,1)] -= secondBubble;
    face[Hex27Index(1,1,0)] -= thirdBubble;
    face[Hex27Index(1,1,1)] = coefficient[Hex27Index(1,1,1)]+firstBubble+secondBubble+thirdBubble-interiorBubble;
    const unsigned dependent = 3U-outer-inner;
    std::array<long double,3> origin = {0,0,0};
    std::array<long double,3> first = {0,0,0};
    std::array<long double,3> second = {0,0,0};
    origin[dependent] = 1.0L;
    first[inner] = second[outer] = 1.0L;
    first[dependent] = second[dependent] = -1.0L;
    const std::vector<std::vector<long double>> polynomial = FacePolynomial(face,origin,first,second);
    std::vector<std::vector<long double>> result(3,std::vector<long double>(3,0.0L));
    for(unsigned i = 0; i < polynomial.size(); ++i) for(unsigned j = 0; j < polynomial[i].size(); ++j) if(i < 3 && j < 3) result[j][i] = polynomial[i][j];
    return result;
  }

  static std::vector<long double> Events(const GraphCellData& data, const GraphCellRegion& cell, unsigned axis, bool& complete, bool& mapped) {
    ++Statistics().eventBuilds;
    std::array<long double,3> lower;
    std::array<long double,3> upper;
    Bounds(cell,lower,upper);
    std::vector<long double> events = {lower[axis],upper[axis]};
    const std::vector<std::array<long double,3>> vertices = Corners(cell);
    const std::vector<std::vector<unsigned>> edges = Edges(cell);
    for(unsigned n = 0; n < vertices.size(); ++n) events.push_back(vertices[n][axis]);
    for(unsigned n = 0; n < edges.size(); ++n) {
      const std::array<long double,3> first = vertices[edges[n][0]];
      const std::array<long double,3> second = vertices[edges[n][1]];
      const std::vector<long double> polynomial = EdgePolynomial(data.coefficient,first,second);
      std::vector<long double> roots;
      LineSliceUtilities::PolynomialRoots(polynomial,0.0L,1.0L,roots);
      for(unsigned k = 0; k < roots.size(); ++k) events.push_back(std::fma(roots[k],second[axis]-first[axis],first[axis]));
    }
    complete = true;
    mapped = false;
    if(!cell.tetrahedron) {
      for(unsigned fixed = 0; fixed < 3; ++fixed) if(fixed != axis) for(unsigned side = 0; side < 2; ++side) {
        std::vector<std::vector<long double>> face = SlicePolynomial(data.coefficient,fixed,cell.vertices[side][fixed]);
        unsigned first;
        unsigned second;
        Hex27OtherAxes(fixed,first,second);
        if(second == axis) for(unsigned i = 0; i < 3; ++i) for(unsigned j = i+1U; j < 3; ++j) std::swap(face[i][j],face[j][i]);
        const unsigned inner = first == axis ? second : first;
        const LineIntegrationStrip strip = {lower[axis],upper[axis],lower[inner],0.0L,upper[inner],0.0L};
        std::vector<long double> faceEvents;
        std::vector<long double> tangencies;
        LineIntegration::CollectEvents(face,strip,faceEvents,&tangencies);
        mapped = mapped || !tangencies.empty();
        events.insert(events.end(),faceEvents.begin(),faceEvents.end());
      }
    }
    else {
      const std::vector<std::vector<unsigned>> faces = cell.faces.empty() ? std::vector<std::vector<unsigned>>{{0,1,2},{0,1,3},{0,2,3},{1,2,3}} : cell.faces;
      for(unsigned face = 0; face < faces.size(); ++face) {
        unsigned fixed = 3;
        for(unsigned d = 0; d < 3; ++d) {
          bool constant = true;
          for(unsigned vertex = 1; vertex < faces[face].size(); ++vertex) constant = constant && cell.vertices[faces[face][vertex]][d] == cell.vertices[faces[face][0]][d];
          if(constant) fixed = d;
        }
        if(fixed == axis) continue;
        unsigned first;
        unsigned second;
        Hex27OtherAxes(fixed < 3 ? fixed : axis,first,second);
        unsigned inner = first;
        std::vector<std::vector<long double>> polynomial;
        if(fixed < 3) {
          polynomial = SlicePolynomial(data.coefficient,fixed,cell.vertices[faces[face][0]][fixed]);
          if(second == axis) for(unsigned i = 0; i < 3; ++i) for(unsigned j = i+1U; j < 3; ++j) std::swap(polynomial[i][j],polynomial[j][i]);
          inner = first == axis ? second : first;
        }
        else polynomial = TetrahedronBoundary(data.coefficient,axis,inner);
        const LineIntegrationStrip strip = {lower[axis],upper[axis],lower[inner],0.0L,upper[inner],0.0L};
        std::vector<long double> faceEvents;
        std::vector<long double> tangencies;
        LineIntegration::CollectEvents(polynomial,strip,faceEvents,&tangencies);
        mapped = mapped || !tangencies.empty();
        events.insert(events.end(),faceEvents.begin(),faceEvents.end());
      }
    }
    LineSliceUtilities::SortDistinct(events);
    return events;
  }

  static void Adaptive(GraphCellData& data, const GraphCellRegion& cell, unsigned axis, long double left, long double right, long double parameterLeft, long double parameterRight, bool mapped, bool projectionComplete, long double budgetVolume, unsigned depth, std::vector<long double>& total, std::vector<long double>& totalErrors, const std::vector<unsigned char>& requested) {
    std::vector<unsigned char> active = requested;
    bool needed = false;
    for(unsigned n = 0; n < active.size(); ++n) needed = needed || active[n] != 0;
    if(!needed) return;
    const long double pi = std::acos(-1.0L);
    const long double heightLeft = LineIntegration::MappedPoint(left,right,parameterLeft,mapped);
    const long double heightRight = LineIntegration::MappedPoint(left,right,parameterRight,mapped);
    if(!(heightRight > heightLeft)) throw std::runtime_error("The 3D graph interval reached floating point resolution.");
    const std::vector<GraphCellRegion> slabs = Slab(cell,axis,heightLeft,heightRight);
    int phase = 2;
    for(unsigned n = 0; n < slabs.size(); ++n) {
      const int childPhase = Phase(data,slabs[n]);
      if(phase == 2) phase = childPhase;
      else if(childPhase != phase) phase = 0;
    }
    if(phase != 0 && phase != 2) {
      if(phase > 0) for(unsigned n = 0; n < slabs.size(); ++n) Full(data,slabs[n],total,totalErrors,active);
      ++data.accepted;
      return;
    }
    static const long double nodes[8] = {0.991455371120812639206854697526329L,0.949107912342758524526189684047851L,0.864864423359769072789712788640926L,0.741531185599394439863864773280788L,0.586087235467691130294144838258730L,0.405845151377397166906606412076961L,0.207784955007898467600689403773245L,0.0L};
    static const long double highWeights[8] = {0.022935322010529224963732008058970L,0.063092092629978553290700663189204L,0.104790010322250183839876322541518L,0.140653259715525918745189590510238L,0.169004726639267902826583426598550L,0.190350578064785409913256402421014L,0.204432940075298892414161999234649L,0.209482141084727828012999174891714L};
    static const long double lowWeights[4] = {0.129484966168869693270611432679082L,0.279705391489276667901467771423780L,0.381830050505118944950369775488975L,0.417959183673469387755102040816327L};
    const unsigned count = data.powers.size();
    std::vector<long double> high(2U*count,0.0L);
    std::vector<long double> low(2U*count,0.0L);
    std::vector<long double> inner(2U*count,0.0L);
    std::vector<long double> absolute(2U*count,0.0L);
    std::vector<std::vector<long double>> samples;
    std::vector<long double> weights;
    const long double middle = 0.5L*(parameterLeft+parameterRight);
    const long double halfWidth = 0.5L*(parameterRight-parameterLeft);
    bool anyCut = false;
    int consensus = 2;
    for(unsigned node = 0; node < 8; ++node) for(unsigned side = 0; side < (node == 7 ? 1U : 2U); ++side) {
      const long double parameter = middle+(side == 0 ? -1.0L : 1.0L)*halfWidth*nodes[node];
      const long double point = LineIntegration::MappedPoint(left,right,parameter,mapped);
      const long double jacobian = mapped ? 0.5L*pi*(right-left)*std::sin(pi*parameter) : right-left;
      std::vector<long double> values;
      std::vector<long double> errors;
      bool cut = false;
      Sample(data,cell,axis,point,active,values,errors,cut,budgetVolume/(heightRight-heightLeft));
      anyCut = anyCut || cut;
      if(!cut) {
        const int samplePhase = values[0] == 0.0L ? -1 : 1;
        if(consensus == 2) consensus = samplePhase;
        else if(consensus != samplePhase) consensus = 0;
      }
      for(unsigned n = 0; n < active.size(); ++n) if(active[n]) {
        values[n] *= jacobian;
        errors[n] *= std::abs(jacobian);
        high[n] += highWeights[node]*values[n];
        inner[n] += highWeights[node]*errors[n];
        absolute[n] += highWeights[node]*std::abs(values[n]);
        if(node%2U != 0) low[n] += lowWeights[node/2U]*values[n];
      }
      samples.push_back(values);
      weights.push_back(highWeights[node]);
    }
    bool remaining = false;
    for(unsigned n = 0; n < active.size(); ++n) {
      if(!active[n]) continue;
      long double variation = 0.0L;
      for(unsigned q = 0; q < samples.size(); ++q) variation += weights[q]*std::abs(samples[q][n]-0.5L*high[n]);
      long double error = std::abs(high[n]-low[n]);
      if(error > 0.0L && variation > 0.0L) error = variation*std::min(1.0L,std::pow(200.0L*error/variation,1.5L));
      error = halfWidth*(std::max(error,128.0L*std::numeric_limits<long double>::epsilon()*absolute[n])+inner[n]);
      high[n] *= halfWidth;
      const long double tolerance = (n < count ? data.absoluteH : data.absoluteD)*budgetVolume+data.relative*std::abs(high[n]);
      if((anyCut || (projectionComplete && consensus != 0)) && std::isfinite(high[n]) && std::isfinite(error) && error <= tolerance) {
        total[n] += high[n];
        totalErrors[n] += error;
        active[n] = 0;
        ++Statistics().acceptedComponents;
      }
      else remaining = true;
    }
    if(!remaining) {
      ++data.accepted;
      return;
    }
    if(depth >= data.maximumDepth || middle == parameterLeft || middle == parameterRight) throw std::runtime_error("The 3D graph moments did not reach their requested tolerance.");
    Adaptive(data,cell,axis,left,right,parameterLeft,middle,mapped,projectionComplete,0.5L*budgetVolume,depth+1U,total,totalErrors,active);
    Adaptive(data,cell,axis,left,right,middle,parameterRight,mapped,projectionComplete,0.5L*budgetVolume,depth+1U,total,totalErrors,active);
  }

  static unsigned Direction(const GraphCellData& data, const GraphCellRegion& cell, bool& weak, long double* energy) {
    unsigned selected = 3;
    long double best = -1.0L;
    weak = false;
    for(unsigned axis = 0; axis < 3; ++axis) {
      long double minimum;
      long double maximum;
      long double padding;
      Range(Derivative(data.coefficient,axis),cell,3,minimum,maximum,padding);
      energy[axis] = std::max(std::abs(minimum),std::abs(maximum));
      const long double lower = minimum > padding ? minimum-padding : (maximum < -padding ? -maximum-padding : 0.0L);
      if(lower == 0.0L) continue;
      const bool monotone = (minimum >= -padding && maximum > padding) || (maximum <= padding && minimum < -padding);
      if(!monotone) continue;
      const long double score = axis == data.preferredInner ? std::numeric_limits<long double>::infinity() : (lower > 0.0L ? 2.0L+lower/std::max(energy[axis],std::numeric_limits<long double>::min()) : energy[axis]/(1.0L+energy[axis]));
      if(score > best) {
        best = score;
        selected = axis;
        weak = lower == 0.0L;
      }
    }
    return selected;
  }

  static bool CriticalSplit(const GraphCellData& data, const GraphCellRegion& cell, unsigned axis, long double& coordinate) {
    std::array<long double,3> lower;
    std::array<long double,3> upper;
    Bounds(cell,lower,upper);
    std::array<long double,3> point;
    for(unsigned d = 0; d < 3; ++d) point[d] = 0.5L*(lower[d]+upper[d]);
    std::vector<std::vector<long double>> gradient(3);
    for(unsigned d = 0; d < 3; ++d) gradient[d] = Derivative(data.coefficient,d);
    for(unsigned iteration = 0; iteration < 8; ++iteration) {
      std::array<long double,3> values;
      std::array<long double,3> columns[3];
      for(unsigned i = 0; i < 3; ++i) {
        values[i] = Value(gradient[i],point);
        for(unsigned j = 0; j < 3; ++j) columns[j][i] = Value(Derivative(gradient[i],j),point);
      }
      const long double determinant = Determinant(columns[0],columns[1],columns[2]);
      if(determinant == 0.0L || !std::isfinite(determinant)) return false;
      const long double step[3] = {Determinant(values,columns[1],columns[2])/determinant,Determinant(columns[0],values,columns[2])/determinant,Determinant(columns[0],columns[1],values)/determinant};
      long double largest = 0.0L;
      for(unsigned d = 0; d < 3; ++d) {
        point[d] -= step[d];
        largest = std::max(largest,std::abs(step[d]));
        if(!std::isfinite(point[d]) || point[d] < lower[d] || point[d] > upper[d]) return false;
      }
      if(largest <= 64.0L*std::numeric_limits<long double>::epsilon()) break;
    }
    std::vector<long double> polynomial(3,0.0L);
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) {
      const unsigned powers[3] = {i,j,k};
      long double value = data.coefficient[Hex27Index(i,j,k)];
      for(unsigned d = 0; d < 3; ++d) if(d != axis) value *= ReferencePower(point[d],powers[d]);
      polynomial[powers[axis]] += value;
    }
    const std::vector<long double> roots = LineSliceUtilities::QuadraticRoots(polynomial[2],polynomial[1],polynomial[0]);
    long double best = std::numeric_limits<long double>::infinity();
    bool found = false;
    for(unsigned root = 0; root < roots.size(); ++root) {
      const long double candidate = 0.5L*(point[axis]+roots[root]);
      if(candidate <= lower[axis]+0.05L*(upper[axis]-lower[axis]) || candidate >= upper[axis]-0.05L*(upper[axis]-lower[axis]) || roots[root] == point[axis]) continue;
      const long double distance = std::abs(candidate-0.5L*(lower[axis]+upper[axis]));
      if(distance < best) {
        best = distance;
        coordinate = candidate;
        found = true;
      }
    }
    return found;
  }

  static void Split(const GraphCellData& data, const GraphCellRegion& cell, unsigned preferred, std::vector<GraphCellRegion>& children) {
    std::array<long double,3> lower;
    std::array<long double,3> upper;
    Bounds(cell,lower,upper);
    std::array<long double,3> center;
    for(unsigned d = 0; d < 3; ++d) center[d] = 0.5L*(lower[d]+upper[d]);
    if(preferred > 2) {
      long double energy[3];
      bool weak;
      Direction(data,cell,weak,energy);
      preferred = energy[1] > energy[0] ? 1U : 0U;
      if(energy[2] > energy[preferred]) preferred = 2;
    }
    const std::vector<long double> derivative = Derivative(data.coefficient,preferred);
    unsigned axis = 0;
    long double largest = 0.0L;
    for(unsigned d = 0; d < 3; ++d) {
      const long double variation = std::abs(Value(Derivative(derivative,d),center))*(upper[d]-lower[d]);
      if(variation > largest) {
        largest = variation;
        axis = d;
      }
    }
    if(largest == 0.0L) {
      for(unsigned d = 1; d < 3; ++d) if(upper[d]-lower[d] > upper[axis]-lower[axis]) axis = d;
    }
    long double coordinate = center[axis];
    const long double slope = Value(Derivative(derivative,axis),center);
    if(slope != 0.0L) {
      const long double root = center[axis]-Value(derivative,center)/slope;
      if(root > lower[axis]+0.05L*(upper[axis]-lower[axis]) && root < upper[axis]-0.05L*(upper[axis]-lower[axis])) coordinate = root;
    }
    CriticalSplit(data,cell,axis,coordinate);
    const long double resolution = 256.0L*std::numeric_limits<long double>::epsilon()*std::max(1.0L,std::max(std::abs(lower[axis]),std::abs(upper[axis])));
    for(unsigned vertex = 0; vertex < cell.vertices.size(); ++vertex) if(cell.vertices[vertex][axis] > lower[axis]+0.05L*(upper[axis]-lower[axis]) && cell.vertices[vertex][axis] < upper[axis]-0.05L*(upper[axis]-lower[axis]) && std::abs(coordinate-cell.vertices[vertex][axis]) <= resolution) coordinate = cell.vertices[vertex][axis];
    Clip(cell,axis,coordinate,true,children);
    Clip(cell,axis,coordinate,false,children);
    if(children.size() < 2) throw std::runtime_error("The 3D graph partition did not produce two positive-volume regions.");
  }

  static bool Independent(const GraphCellData& data, unsigned axis) {
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) {
      const unsigned powers[3] = {i,j,k};
      if(powers[axis] != 0 && data.coefficient[Hex27Index(i,j,k)] != 0.0L) return false;
    }
    for(unsigned n = 0; n < data.gradient.size(); ++n) {
      const unsigned size = SurfacePolynomialCubeSize(data.gradient[n].size());
      for(unsigned i = 0; i < size; ++i) for(unsigned j = 0; j < size; ++j) for(unsigned k = 0; k < size; ++k) {
        const unsigned powers[3] = {i,j,k};
        if(powers[axis] != 0 && data.gradient[n][size*size*i+size*j+k] != 0.0L) return false;
      }
    }
    return true;
  }

  static bool BoundSmallRegion(const GraphCellData& data, const GraphCellRegion& cell, unsigned inner, std::vector<long double>& errors, const std::vector<unsigned char>& active) {
    std::array<long double,3> regionLower;
    std::array<long double,3> regionUpper;
    Bounds(cell,regionLower,regionUpper);
    const long double boxVolume = (regionUpper[0]-regionLower[0])*(regionUpper[1]-regionLower[1])*(regionUpper[2]-regionLower[2]);
    const long double boundH = Volume(cell)+64.0L*std::numeric_limits<long double>::epsilon()*boxVolume;
    if(boundH > 0.0001L*data.absoluteH*data.domainVolume || (data.computeDelta && inner > 2)) return false;
    long double boundD = 0.0L;
    if(data.computeDelta) {
      long double minimum;
      long double maximum;
      long double padding;
      Range(Derivative(data.coefficient,inner),cell,3,minimum,maximum,padding);
      const long double derivative = minimum > padding ? minimum-padding : -maximum-padding;
      if(!(derivative > 0.0L)) return false;
      unsigned first;
      unsigned second;
      Hex27OtherAxes(inner,first,second);
      std::array<long double,3> lower;
      std::array<long double,3> upper;
      Bounds(cell,lower,upper);
      long double projection = (upper[first]-lower[first])*(upper[second]-lower[second]);
      if(cell.tetrahedron && cell.faces.empty()) {
        const long double boxArea = projection;
        projection = 0.0L;
        for(unsigned omitted = 0; omitted < 4; ++omitted) {
          std::vector<unsigned> face;
          for(unsigned vertex = 0; vertex < 4; ++vertex) if(vertex != omitted) face.push_back(vertex);
          const long double a = cell.vertices[face[1]][first]-cell.vertices[face[0]][first];
          const long double b = cell.vertices[face[1]][second]-cell.vertices[face[0]][second];
          const long double c = cell.vertices[face[2]][first]-cell.vertices[face[0]][first];
          const long double d = cell.vertices[face[2]][second]-cell.vertices[face[0]][second];
          projection += 0.25L*std::abs(std::fma(a,d,-b*c));
        }
        projection += 64.0L*std::numeric_limits<long double>::epsilon()*boxArea;
      }
      long double metric = data.gradient.empty() ? 1.0L : 0.0L;
      for(unsigned component = 0; component < data.gradient.size(); ++component) for(unsigned n = 0; n < data.gradient[component].size(); ++n) metric += std::abs(data.gradient[component][n]);
      boundD = projection*metric/derivative;
      if(boundD > 0.0001L*data.absoluteD*data.domainVolume) return false;
    }
    const unsigned count = data.powers.size();
    for(unsigned n = 0; n < count; ++n) {
      if(active[n]) errors[n] += boundH;
      if(active[count+n]) errors[count+n] += boundD;
    }
    ++Statistics().boundedCells;
    return true;
  }

  static int RefinedPhase(const GraphCellData& data, const GraphCellRegion& cell) {
    std::array<long double,3> lower;
    std::array<long double,3> upper;
    Bounds(cell,lower,upper);
    const std::vector<long double> power = BoxPolynomial(data.coefficient,lower,upper);
    std::vector<long double> controls = Hex27Bernstein(power);
    unsigned first = 0;
    for(unsigned axis = 1; axis < 3; ++axis) if(upper[axis]-lower[axis] > upper[first]-lower[first]) first = axis;
    unsigned second = first == 0 ? 1U : 0U;
    for(unsigned axis = 0; axis < 3; ++axis) if(axis != first && upper[axis]-lower[axis] > upper[second]-lower[second]) second = axis;
    unsigned budget = 64;
    const int phase = Hex27BernsteinPhase(controls,first,second,0,8,budget,data.computeDelta);
    if(phase == 0) return 0;
    long double parentNorm = 0.0L;
    long double localNorm = 0.0L;
    for(unsigned n = 0; n < data.coefficient.size(); ++n) parentNorm += std::abs(data.coefficient[n]);
    for(unsigned n = 0; n < power.size(); ++n) localNorm += std::abs(power[n]);
    const long double margin = 512.0L*std::numeric_limits<long double>::epsilon()*std::max(parentNorm,localNorm);
    for(unsigned n = 0; n < controls.size(); ++n) controls[n] += phase > 0 ? -margin : margin;
    budget = 64;
    return Hex27BernsteinPhase(controls,first,second,0,8,budget,data.computeDelta) == phase ? phase : 0;
  }

  static void Partition(GraphCellData& data, const GraphCellRegion& cell, std::vector<long double>& total, std::vector<long double>& totalErrors) {
    const unsigned count = data.powers.size();
    std::vector<unsigned char> active = data.requested.empty() ? std::vector<unsigned char>(2U*count,1) : data.requested;
    if(!data.computeDelta) std::fill(active.begin()+count,active.end(),0);
    const int phase = Phase(data,cell);
    if(phase != 0) {
      if(phase > 0) Full(data,cell,total,totalErrors,active);
      ++data.accepted;
      return;
    }
    bool weak;
    long double energy[3];
    const unsigned inner = Direction(data,cell,weak,energy);
    if(inner > 2) {
      const int refinedPhase = RefinedPhase(data,cell);
      if(refinedPhase != 0) {
        if(refinedPhase > 0) Full(data,cell,total,totalErrors,active);
        ++data.accepted;
        return;
      }
    }
    if(BoundSmallRegion(data,cell,inner,totalErrors,active)) return;
    std::runtime_error failure("The 3D region has no certified graph direction.");
    if(inner < 3) {
      unsigned outer = inner == 0 ? 1U : 0U;
      for(unsigned axis = 0; axis < 3; ++axis) if(axis != inner && energy[axis] < energy[outer]) outer = axis;
      if(data.preferredOuter < 3 && data.preferredOuter != inner) outer = data.preferredOuter;
      data.graphAxis = inner;
      ++Statistics().graphCells;
      if(weak) ++Statistics().weakGraphs;
      ++Statistics().directions[inner];
      ++data.innerCounts[inner];
      ++data.outerCounts[outer];
      for(unsigned attempt = 0; attempt < 2; ++attempt) {
        if(attempt != 0) outer = 3U-inner-outer;
        std::vector<long double> values(2U*count,0.0L);
        std::vector<long double> errors(2U*count,0.0L);
        try {
          if(!cell.tetrahedron && Independent(data,outer)) {
            bool cut = false;
            Sample(data,cell,outer,0.5L*(cell.vertices[0][outer]+cell.vertices[1][outer]),active,values,errors,cut,data.domainVolume/(cell.vertices[1][outer]-cell.vertices[0][outer]));
            const std::vector<long double> sourceValues = values;
            const std::vector<long double> sourceErrors = errors;
            for(unsigned n = 0; n < count; ++n) {
              const long double factor = IntegralPower(cell.vertices[0][outer],cell.vertices[1][outer],data.powers[n][outer]);
              unsigned powers[3] = {data.powers[n][0],data.powers[n][1],data.powers[n][2]};
              powers[outer] = 0;
              const unsigned source = HexMomentIndex(powers[0],powers[1],powers[2]);
              for(unsigned pair = 0; pair < 2; ++pair) {
                values[pair*count+n] = factor*sourceValues[pair*count+source];
                errors[pair*count+n] = std::abs(factor)*sourceErrors[pair*count+source];
              }
            }
            ++data.accepted;
          }
          else {
            bool complete;
            bool mapped;
            const std::vector<long double> events = Events(data,cell,outer,complete,mapped);
            const long double intervalBudget = data.domainVolume/(events.size()-1U);
            for(unsigned n = 1; n < events.size(); ++n) if(events[n] > events[n-1U]) Adaptive(data,cell,outer,events[n-1U],events[n],0.0L,1.0L,mapped,complete,intervalBudget,0,values,errors,active);
          }
          for(unsigned n = 0; n < total.size(); ++n) {
            total[n] += values[n];
            totalErrors[n] += errors[n]+8.0L*std::numeric_limits<long double>::epsilon()*std::abs(values[n]);
          }
          return;
        }
        catch(const std::runtime_error& exception) {
          failure = std::runtime_error(exception.what());
          if(data.evaluations >= data.maximumEvaluations) throw;
        }
      }
    }
    if(cell.depth >= data.maximumDepth) {
      ++Statistics().fallbackCells;
      throw failure;
    }
    std::vector<GraphCellRegion> children;
    Split(data,cell,inner,children);
    ++Statistics().splitCells;
    for(unsigned n = 0; n < children.size(); ++n) Partition(data,children[n],total,totalErrors);
  }

  static void Symmetry(const GraphCellData& data, std::vector<long double>& values, std::vector<long double>& errors) {
    const unsigned count = data.powers.size();
    for(unsigned axis = 0; axis < 3; ++axis) {
      bool even = true;
      bool odd = true;
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) {
        const unsigned powers[3] = {i,j,k};
        if(data.coefficient[Hex27Index(i,j,k)] == 0.0L) continue;
        if(powers[axis]%2U != 0) even = false;
        else odd = false;
      }
      bool metric = true;
      for(unsigned n = 0; n < data.gradient.size(); ++n) {
        const unsigned size = SurfacePolynomialCubeSize(data.gradient[n].size());
        bool componentEven = true;
        bool componentOdd = true;
        for(unsigned i = 0; i < size; ++i) for(unsigned j = 0; j < size; ++j) for(unsigned k = 0; k < size; ++k) {
          const unsigned powers[3] = {i,j,k};
          if(data.gradient[n][size*size*i+size*j+k] == 0.0L) continue;
          if(powers[axis]%2U != 0) componentEven = false;
          else componentOdd = false;
        }
        metric = metric && (componentEven || componentOdd);
      }
      for(unsigned n = 0; n < count; ++n) if(data.powers[n][axis]%2U != 0) {
        if(even) values[n] = errors[n] = 0.0L;
        if((even || odd) && metric) values[count+n] = errors[count+n] = 0.0L;
      }
    }
  }

  static void Integrate(GraphCellData& data, bool tetrahedron, std::vector<long double>& values, std::vector<long double>& errors) {
    const unsigned count = (data.degree+1U)*(data.degree+2U)*(data.degree+3U)/6U;
    data.powers.assign(count,{0,0,0});
    for(unsigned i = 0; i <= data.degree; ++i) for(unsigned j = 0; j <= data.degree-i; ++j) for(unsigned k = 0; k <= data.degree-i-j; ++k) data.powers[HexMomentIndex(i,j,k)] = {i,j,k};
    GraphCellRegion cell;
    cell.tetrahedron = tetrahedron;
    cell.vertices = tetrahedron ? std::vector<std::array<long double,3>>{{0,0,0},{1,0,0},{0,1,0},{0,0,1}} : (data.unitCube ? std::vector<std::array<long double,3>>{{0,0,0},{1,1,1}} : std::vector<std::array<long double,3>>{{-1,-1,-1},{1,1,1}});
    values.assign(2U*count,0.0L);
    errors.assign(2U*count,0.0L);
    long double norm = 0.0L;
    for(unsigned n = 0; n < data.coefficient.size(); ++n) norm += std::abs(data.coefficient[n]);
    if(norm == 0.0L) {
      if(data.computeDelta) throw std::domain_error("The delta of an identically zero 3D field is undefined.");
      const std::vector<unsigned char> active(2U*count,1);
      Full(data,cell,values,errors,active,0.5L);
      return;
    }
    const long double absoluteH = data.absoluteH;
    const long double absoluteD = data.absoluteD;
    const long double relative = data.relative;
    data.domainVolume = Volume(cell);
    data.requested.assign(2U*count,1);
    if(!data.computeDelta) std::fill(data.requested.begin()+count,data.requested.end(),0);
    long double budgetFactor = 1.0L;
    for(unsigned pass = 0; pass < 8; ++pass) {
      data.absoluteH = budgetFactor*absoluteH/data.domainVolume;
      data.absoluteD = budgetFactor*absoluteD/data.domainVolume;
      data.relative = budgetFactor*relative;
      std::vector<long double> passValues(2U*count,0.0L);
      std::vector<long double> passErrors(2U*count,0.0L);
      Partition(data,cell,passValues,passErrors);
      if(!tetrahedron && !data.unitCube) Symmetry(data,passValues,passErrors);
      long double worst = 0.0L;
      bool remaining = false;
      for(unsigned n = 0; n < values.size(); ++n) if(data.requested[n]) {
        if(!std::isfinite(passValues[n]) || !std::isfinite(passErrors[n]) || passErrors[n] < 0.0L) throw std::runtime_error("A 3D graph moment is not finite.");
        const long double tolerance = (n < count ? absoluteH : absoluteD)+relative*std::abs(passValues[n]);
        if(passErrors[n] <= 0.9L*tolerance) {
          values[n] = passValues[n];
          errors[n] = passErrors[n];
          data.requested[n] = 0;
        }
        else {
          remaining = true;
          worst = std::max(worst,passErrors[n]/tolerance);
        }
      }
      if(!remaining) {
        data.absoluteH = absoluteH;
        data.absoluteD = absoluteD;
        data.relative = relative;
        data.requested.clear();
        return;
      }
      budgetFactor /= std::max(2.0L,2.0L*worst);
    }
    throw std::runtime_error("The accumulated 3D graph moment error exceeds tolerance.");
  }

  static void Output(const GraphCellData& data, long double scale, const std::vector<long double>& values, const std::vector<long double>& errors, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, std::vector<double>* estimate) {
    const unsigned count = data.powers.size();
    std::vector<std::vector<std::vector<double>>> outputH(data.degree+1U);
    std::vector<std::vector<std::vector<double>>> outputD(data.degree+1U);
    long double errorH = 0.0L;
    long double errorD = 0.0L;
    for(unsigned i = 0; i <= data.degree; ++i) {
      outputH[i].resize(data.degree-i+1U);
      outputD[i].resize(data.degree-i+1U);
      for(unsigned j = 0; j <= data.degree-i; ++j) {
        outputH[i][j].resize(data.degree-i-j+1U);
        outputD[i][j].resize(data.degree-i-j+1U);
        for(unsigned k = 0; k <= data.degree-i-j; ++k) {
          const unsigned n = HexMomentIndex(i,j,k);
          const long double h = values[n];
          const long double d = data.computeDelta ? values[count+n]/scale : 0.0L;
          outputH[i][j][k] = static_cast<double>(h);
          outputD[i][j][k] = static_cast<double>(d);
          if(!std::isfinite(outputH[i][j][k]) || !std::isfinite(outputD[i][j][k])) throw std::overflow_error("A 3D graph moment exceeds double precision range.");
          const long double hError = errors[n]+std::abs(h-static_cast<long double>(outputH[i][j][k]));
          const long double dError = data.computeDelta ? errors[count+n]/scale+std::abs(d-static_cast<long double>(outputD[i][j][k])) : 0.0L;
          if(hError > data.absoluteH+data.relative*std::abs(h) || dError > data.absoluteD/scale+data.relative*std::abs(d)) throw std::runtime_error("A rounded 3D graph moment exceeds its error budget.");
          errorH = std::max(errorH,hError);
          errorD = std::max(errorD,dError);
        }
      }
    }
    positive.swap(outputH);
    delta.swap(outputD);
    if(estimate != nullptr) {
      unsigned selected = 0;
      for(unsigned axis = 1; axis < 3; ++axis) if(data.outerCounts[axis] > data.outerCounts[selected]) selected = axis;
      *estimate = {static_cast<double>(errorH),static_cast<double>(errorD),static_cast<double>(data.accepted),static_cast<double>(data.evaluations),data.evaluations == 0 ? 0.0 : static_cast<double>(selected+1U)};
    }
  }
};

}

#endif
