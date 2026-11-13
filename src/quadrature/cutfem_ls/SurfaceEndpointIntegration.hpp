#ifndef FEMUS_SURFACE_ENDPOINT_INTEGRATION_HPP
#define FEMUS_SURFACE_ENDPOINT_INTEGRATION_HPP

#include "Quad9DirectionSelection.hpp"
#include "TriangleIntegration.hpp"
#include "Hex27LevelSetMoments.hpp"
#include <array>
#include <map>
#include <set>
#include "EndpointProjection.hpp"

namespace FELevelSetMoments {

struct SurfaceEndpointStatistics {
  unsigned long long endpointStrips = 0;
  unsigned long long graphCells = 0;
  unsigned long long weakGraphs = 0;
  unsigned long long generalSurfaces = 0;
  unsigned long long quadCalls = 0;
  unsigned long long triCalls = 0;
  unsigned long long outerSamples = 0;
  unsigned long long eventBuilds = 0;
  unsigned long long integrationPasses = 0;
  unsigned long long acceptedComponents = 0;
  unsigned long long skippedComponents = 0;


  unsigned directions[3] = {0,0,0};
};

struct SurfaceEndpointDomain {
  bool tetrahedron = false;
  bool wedge = false;
  std::vector<std::array<long double,3>> vertices;

};

struct SurfaceEndpointData {
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
  bool wedge = false;
  std::vector<unsigned> zeroFaces;
  long double domainVolume = 1.0L;
  std::vector<unsigned char> requested;
  std::array<unsigned,3> innerCounts = {0,0,0};
  std::array<unsigned,3> outerCounts = {0,0,0};
  bool computeDelta = true;
  bool inverseSurface = false;
  long double absoluteH = 0.0L;
  long double absoluteD = 0.0L;
  long double relative = 0.0L;
  OuterQuadratureRules quadratureRules;
};

struct SurfaceOuterContext {
  SurfaceEndpointData* data = nullptr;
  const SurfaceEndpointDomain* cell = nullptr;
  unsigned outer = 0;

};

class SurfaceEndpointIntegration {
public:
  // Integrates one parent element through full QUAD/TRI slices and endpoint quadrature.
  static SurfaceEndpointStatistics& Statistics() {
    static thread_local SurfaceEndpointStatistics data;
    return data;
  }

  static long double Determinant(const std::array<long double,3>& first, const std::array<long double,3>& second, const std::array<long double,3>& third) {
    return std::fma(first[0],std::fma(second[1],third[2],-second[2]*third[1]),std::fma(-first[1],std::fma(second[0],third[2],-second[2]*third[0]),first[2]*std::fma(second[0],third[1],-second[1]*third[0])));
  }



  static long double Volume(const SurfaceEndpointDomain& cell) {
    if(cell.wedge) return 0.5L;
    if(!cell.tetrahedron) return (cell.vertices[1][0]-cell.vertices[0][0])*(cell.vertices[1][1]-cell.vertices[0][1])*(cell.vertices[1][2]-cell.vertices[0][2]);

    std::array<long double,3> edges[3];
    for(unsigned i = 0; i < 3; ++i) for(unsigned d = 0; d < 3; ++d) edges[i][d] = cell.vertices[i+1U][d]-cell.vertices[0][d];
    return std::abs(Determinant(edges[0],edges[1],edges[2]))/6.0L;
  }

  static void Bounds(const SurfaceEndpointDomain& cell, std::array<long double,3>& lower, std::array<long double,3>& upper) {
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

  static std::vector<long double> SimplexPolynomial(const std::vector<long double>& coefficient, const SurfaceEndpointDomain& cell) {
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

  static void Range(const std::vector<long double>& coefficient, const SurfaceEndpointDomain& cell, unsigned degree, long double& minimum, long double& maximum, long double& padding) {
    std::vector<long double> controls;
    long double norm = 0.0L;
    if(cell.tetrahedron) {
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

  static int Phase(const SurfaceEndpointData& data, const SurfaceEndpointDomain& cell) {
    if(data.computeDelta) for(unsigned n = 0; n < data.zeroFaces.size(); ++n) if(data.zeroFaces[n] != 0) return 0;
    long double minimum;
    long double maximum;
    long double padding;
    Range(data.coefficient,cell,4,minimum,maximum,padding);
    if(minimum > padding || (!data.computeDelta && minimum >= 0.0L && maximum > 0.0L)) return 1;
    if(maximum < -padding || (!data.computeDelta && maximum <= 0.0L && minimum < 0.0L)) return -1;
    return 0;
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

  static long double TetMoment(const SurfaceEndpointDomain& cell, unsigned first, unsigned second, unsigned third) {

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

  static void Full(const SurfaceEndpointData& data, const SurfaceEndpointDomain& cell, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<unsigned char>& active, long double multiplier = 1.0L) {
    const unsigned count = data.powers.size();
    for(unsigned n = 0; n < count; ++n) {
      if(!active[n]) continue;
      long double value = cell.wedge ? ReferenceTriangleMoment(data.powers[n][0],data.powers[n][1])/(data.powers[n][2]+1U) : cell.tetrahedron ? TetMoment(cell,data.powers[n][0],data.powers[n][1],data.powers[n][2]) : IntegralPower(cell.vertices[0][0],cell.vertices[1][0],data.powers[n][0])*IntegralPower(cell.vertices[0][1],cell.vertices[1][1],data.powers[n][1])*IntegralPower(cell.vertices[0][2],cell.vertices[1][2],data.powers[n][2]);
      value *= multiplier;
      values[n] += value;
      errors[n] += (64.0L+16.0L*data.degree)*std::numeric_limits<long double>::epsilon()*std::abs(value);
    }
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

  static unsigned Direction(const SurfaceEndpointData& data, const SurfaceEndpointDomain& cell, bool& weak, long double* energy) {
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
      if(lower == 0.0L || axis == data.preferredOuter) continue;
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

  static bool Independent(const SurfaceEndpointData& data, unsigned axis) {
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

  static void Symmetry(const SurfaceEndpointData& data, std::vector<long double>& values, std::vector<long double>& errors) {
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


  static void Sample(SurfaceEndpointData& data, const SurfaceEndpointDomain& cell, unsigned outer, long double coordinate, const std::vector<unsigned char>& active, std::vector<long double>& values, std::vector<long double>& errors, long double allocation) {
    if(data.evaluations >= data.maximumEvaluations) throw std::runtime_error("3D endpoint evaluation budget exhausted.");
    ++data.evaluations;
    ++Statistics().outerSamples;
    unsigned first;
    unsigned second;
    Hex27OtherAxes(outer,first,second);
    unsigned degree = 0;
    const unsigned count = data.powers.size();
    for(unsigned n = 0; n < count; ++n) if(active[n] || active[count+n]) degree = std::max(degree,data.powers[n][first]+data.powers[n][second]);
    const unsigned sliceCount = (degree+1U)*(degree+2U)/2U;
    const unsigned preferred = data.graphAxis < 3 ? data.graphAxis : data.preferredInner;
    const unsigned requestedDirection = preferred < 3 && preferred != outer ? (preferred == first ? 2U : 1U) : 0;
    std::vector<unsigned char> required(2U*sliceCount,0);
    for(unsigned n = 0; n < count; ++n) for(unsigned pair = 0; pair < 2; ++pair) if(active[pair*count+n]) required[pair*sliceCount+ReferenceMomentIndex2D(degree,data.powers[n][first],data.powers[n][second])] = 1;
    const long double origin = cell.tetrahedron || data.unitCube || cell.wedge ? 0 : -1;
    const bool triangleSlice = cell.tetrahedron || (cell.wedge && outer == 2);
    const long double width = cell.tetrahedron || (cell.wedge && outer < 2) ? 1-coordinate : (data.unitCube || cell.wedge ? 1 : 2);
    const long double secondWidth = cell.wedge ? 1 : width;
    if(!(width > 0)) throw std::runtime_error("A tetrahedral endpoint slice reached coordinate resolution.");
    const std::vector<std::vector<long double>> slice = SlicePolynomial(data.coefficient,outer,coordinate);
    std::vector<std::vector<long double>> chart;
    ReferencePolynomialAffine2D(slice,origin,width,origin,secondWidth,chart);
    std::vector<std::vector<std::vector<long double>>> gradient(data.gradient.size());
    for(unsigned component = 0; component < data.gradient.size(); ++component) ReferencePolynomialAffine2D(SlicePolynomial(data.gradient[component],outer,coordinate),origin,width,origin,secondWidth,gradient[component]);
    std::vector<unsigned> boundary;
    if(!data.zeroFaces.empty()) {
      if(cell.wedge) boundary = outer == 2 ? std::vector<unsigned>{data.zeroFaces[0],data.zeroFaces[1],data.zeroFaces[2]} : std::vector<unsigned>{data.zeroFaces[first],data.zeroFaces[2],data.zeroFaces[3],data.zeroFaces[4]};
      else if(cell.tetrahedron) boundary = {data.zeroFaces[first+1U],data.zeroFaces[second+1U],data.zeroFaces[0]};
      else boundary = {data.zeroFaces[2U*first],data.zeroFaces[2U*first+1U],data.zeroFaces[2U*second],data.zeroFaces[2U*second+1U]};
    }
    const std::array<long double,5> parent = {origin,width,origin,secondWidth,width*secondWidth};
    std::vector<long double> moments;
    std::vector<long double> momentErrors;
    unsigned evaluations = 0;
    unsigned accepted = 0;
    unsigned selected = 0;
    try {
      if(triangleSlice) {
        ++Statistics().triCalls;
        TriangleIntegration::Integrate(chart,degree,data.computeDelta,0.01L*allocation*data.absoluteH,0.01L*allocation*data.absoluteD,0.01L*data.relative,data.maximumDepth,data.maximumEvaluations-data.evaluations,moments,momentErrors,evaluations,accepted,selected,gradient.empty() ? nullptr : &gradient,&parent,&required,requestedDirection,false,boundary.empty() ? nullptr : &boundary,data.quadratureRules,data.inverseSurface);
      }
      else {
        ++Statistics().quadCalls;
        Quad9SelectedDirect(chart,degree,requestedDirection,data.computeDelta,0.01L*allocation*data.absoluteH,0.01L*allocation*data.absoluteD,0.01L*data.relative,data.maximumDepth,data.maximumEvaluations-data.evaluations,moments,momentErrors,evaluations,accepted,gradient.empty() ? nullptr : &gradient,parent,&required,boundary.empty() ? nullptr : &boundary,data.quadratureRules,data.inverseSurface);
      }
    }
    catch(...) { data.evaluations += evaluations; throw; }
    data.evaluations += evaluations;
    values.assign(2U*count,0);
    errors.assign(2U*count,0);
    for(unsigned n = 0; n < count; ++n) {
      if(!active[n] && !active[count+n]) { ++Statistics().skippedComponents; continue; }
      const unsigned source = ReferenceMomentIndex2D(degree,data.powers[n][first],data.powers[n][second]);
      const long double power = ReferencePower(coordinate,data.powers[n][outer]);
      for(unsigned pair = 0; pair < 2; ++pair) if(active[pair*count+n]) {
        values[pair*count+n] = power*moments[pair*sliceCount+source];
        errors[pair*count+n] = std::abs(power)*momentErrors[pair*sliceCount+source]+8*std::numeric_limits<long double>::epsilon()*std::abs(values[pair*count+n]);
      }
    }
  }

  static void OuterSample(SurfaceOuterContext& context, const OuterIntervalPoint& point, const std::vector<unsigned char>& active, std::vector<long double>& values, std::vector<long double>& errors, long double absoluteH, long double absoluteD, long double relative) {
    Sample(*context.data,*context.cell,context.outer,point.coordinate,active,values,errors,absoluteH/context.data->absoluteH);
    for(unsigned n = 0; n < active.size(); ++n) if(active[n]) { values[n] *= point.jacobian; errors[n] *= point.jacobian; }
  }

  static void Adaptive(SurfaceEndpointData& data, const SurfaceEndpointDomain& cell, unsigned outer, long double left, long double right, int mapping, long double endpoint, long double delta, unsigned exponent, long double budget, unsigned depth, const std::vector<unsigned char>& requested, std::vector<long double>& total, std::vector<long double>& totalErrors) {
    SurfaceOuterContext context;
    context.data = &data;
    context.cell = &cell;
    context.outer = outer;
    const OuterInterval interval = {mapping == 0 ? left : (mapping > 0 ? endpoint : endpoint-delta),mapping == 0 ? right : (mapping > 0 ? endpoint+delta : endpoint),endpoint,delta,mapping,exponent,true};
    OuterIntervalQuadrature::Adaptive(context,OuterSample,interval,left,right,budget*data.absoluteH,budget*data.absoluteD,data.relative,depth,data.maximumDepth,data.maximumEvaluations,data.evaluations,data.accepted,3,requested,total,totalErrors,data.quadratureRules);
  }

  static void IntegrateParent(SurfaceEndpointData& data, const SurfaceEndpointDomain& cell, std::vector<long double>& total, std::vector<long double>& totalErrors) {
    const std::vector<unsigned char> active = data.requested;
    const int phase = Phase(data,cell);
    if(phase != 0) { if(phase > 0) Full(data,cell,total,totalErrors,active); ++data.accepted; return; }
    if(data.preferredOuter < 3) {
      long double axial = 0;
      long double transverse = 0;
      for(unsigned axis = 0; axis < 3; ++axis) {
        const std::vector<long double> derivative = Derivative(data.coefficient,axis);
        long double norm = 0;
        for(unsigned n = 0; n < derivative.size(); ++n) norm += std::abs(derivative[n]);
        if(axis == data.preferredOuter) axial = norm; else transverse += norm;
      }
      if(axial > 0 && transverse <= 1e-10L*axial) data.preferredOuter = 3;
    }
    bool weak = false;
    long double energy[3];
    const unsigned inner = Direction(data,cell,weak,energy);
    data.graphAxis = inner;
    if(inner < 3) ++Statistics().graphCells; else ++Statistics().generalSurfaces;
    std::vector<unsigned> axes;
    for(unsigned axis = 0; axis < 3; ++axis) if(axis != inner) axes.push_back(axis);
    for(unsigned i = 0; i < axes.size(); ++i) for(unsigned j = i+1U; j < axes.size(); ++j) if((axes[j] == data.preferredOuter) || (axes[i] != data.preferredOuter && energy[axes[j]] < energy[axes[i]])) std::swap(axes[i],axes[j]);
    std::runtime_error failure("No supported 3D endpoint projection.");
    for(unsigned attempt = 0; attempt < axes.size(); ++attempt) {
      const unsigned outer = axes[attempt];
      if(data.preferredOuter < 3 && outer != data.preferredOuter) continue;
      std::vector<long double> values(total.size(),0);
      std::vector<long double> errors(total.size(),0);
      try {
        ++data.outerCounts[outer];
        if(inner < 3) ++data.innerCounts[inner];
        if(!cell.tetrahedron && (!cell.wedge || outer == 2) && Independent(data,outer)) {
          std::vector<unsigned char> dependencies = active;
          for(unsigned n = 0; n < data.powers.size(); ++n) {
            unsigned powers[3] = {data.powers[n][0],data.powers[n][1],data.powers[n][2]};
            powers[outer] = 0;
            for(unsigned pair = 0; pair < 2; ++pair) if(active[pair*data.powers.size()+n]) dependencies[pair*data.powers.size()+HexMomentIndex(powers[0],powers[1],powers[2])] = 1;
          }
          Sample(data,cell,outer,0,dependencies,values,errors,data.domainVolume/(data.unitCube || cell.wedge ? 1 : 2));
          const std::vector<long double> source = values;
          const std::vector<long double> sourceError = errors;
          for(unsigned n = 0; n < data.powers.size(); ++n) {
            unsigned powers[3] = {data.powers[n][0],data.powers[n][1],data.powers[n][2]};
            const long double factor = IntegralPower(data.unitCube || cell.wedge ? 0 : -1,1,powers[outer]);
            powers[outer] = 0;
            const unsigned index = HexMomentIndex(powers[0],powers[1],powers[2]);
            for(unsigned pair = 0; pair < 2; ++pair) {
              values[pair*data.powers.size()+n] = factor*source[pair*data.powers.size()+index];
              errors[pair*data.powers.size()+n] = std::abs(factor)*sourceError[pair*data.powers.size()+index];
            }
          }
          ++data.accepted;
        }
        else {
          std::vector<long double> chart = data.coefficient;
          if(!cell.tetrahedron && !data.unitCube && !cell.wedge) ReferenceCubePolynomialToChart(data.coefficient,chart);
          ++Statistics().eventBuilds;
          const std::vector<EndpointProjectionEvent> events = EndpointProjection::Build(chart,cell.tetrahedron,outer,inner < 3 ? inner : data.preferredInner,inner < 3,cell.wedge);
          const long double origin = cell.tetrahedron || data.unitCube || cell.wedge ? 0 : -1;
          const long double scale = cell.tetrahedron || data.unitCube || cell.wedge ? 1 : 2;
          const long double budget = data.domainVolume/(events.size()-1U);
          for(unsigned n = 1; n < events.size(); ++n) {
            const long double left = std::fma(scale,events[n-1U].point,origin);
            const long double right = std::fma(scale,events[n].point,origin);
            if(!(right > left)) continue;
            const std::vector<OuterInterval> plan = OuterIntervalQuadrature::Plan(left,right,events[n-1U].exponent,events[n].exponent,3);
            for(unsigned part = 0; part < plan.size(); ++part) {
              Adaptive(data,cell,outer,plan[part].mapping == 0 ? plan[part].left : 0,plan[part].mapping == 0 ? plan[part].right : 1,plan[part].mapping,plan[part].endpoint,plan[part].delta,plan[part].exponent,budget/plan.size(),0,active,values,errors);
              if(plan[part].mapping != 0) ++Statistics().endpointStrips;
            }
          }
        }
        for(unsigned n = 0; n < total.size(); ++n) { total[n] += values[n]; totalErrors[n] += errors[n]; }
        return;
      }
      catch(const std::runtime_error& error) { failure = std::runtime_error(error.what()); if(data.evaluations >= data.maximumEvaluations) throw; }
    }
    throw failure;
  }

  static void Integrate(SurfaceEndpointData& data, bool tetrahedron, std::vector<long double>& values, std::vector<long double>& errors) {
    data.quadratureRules.Validate();
    const unsigned count = (data.degree+1U)*(data.degree+2U)*(data.degree+3U)/6U;
    data.powers.assign(count,{0,0,0});
    for(unsigned i = 0; i <= data.degree; ++i) for(unsigned j = 0; j <= data.degree-i; ++j) for(unsigned k = 0; k <= data.degree-i-j; ++k) data.powers[HexMomentIndex(i,j,k)] = {i,j,k};
    SurfaceEndpointDomain cell;
    cell.tetrahedron = tetrahedron;
    cell.wedge = data.wedge;
    cell.vertices = tetrahedron ? std::vector<std::array<long double,3>>{{0,0,0},{1,0,0},{0,1,0},{0,0,1}} : (data.unitCube || data.wedge ? std::vector<std::array<long double,3>>{{0,0,0},{1,1,1}} : std::vector<std::array<long double,3>>{{-1,-1,-1},{1,1,1}});
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
    if(!tetrahedron && !data.unitCube && !data.wedge) {
      std::vector<long double> parityValues(2U*count,0);
      std::vector<long double> parityErrors(2U*count,1);
      Symmetry(data,parityValues,parityErrors);
      for(unsigned n = 0; n < data.requested.size(); ++n) if(parityErrors[n] == 0) data.requested[n] = 0;
    }
    long double budgetFactor = 1.0L;
    for(unsigned pass = 0; pass < 8; ++pass) {
      data.absoluteH = budgetFactor*absoluteH/data.domainVolume;
      data.absoluteD = budgetFactor*absoluteD/data.domainVolume;
      data.relative = budgetFactor*relative;
      std::vector<long double> passValues(2U*count,0.0L);
      std::vector<long double> passErrors(2U*count,0.0L);
      ++Statistics().integrationPasses;
      IntegrateParent(data,cell,passValues,passErrors);
      if(!tetrahedron && !data.unitCube && !data.wedge) Symmetry(data,passValues,passErrors);
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

  static void Compare(SurfaceEndpointData& data, bool tetrahedron, const std::vector<long double>& values, std::vector<long double>& errors) {
    SurfaceEndpointData independent = data;
    independent.evaluations = independent.accepted = 0;
    unsigned outer = 0;
    for(unsigned axis = 1; axis < 3; ++axis) if(data.outerCounts[axis] > data.outerCounts[outer]) outer = axis;
    independent.preferredOuter = (outer+1U)%3U;
    if(independent.preferredOuter == data.graphAxis) independent.preferredOuter = (outer+2U)%3U;
    std::vector<long double> second;
    std::vector<long double> secondErrors;
    Integrate(independent,tetrahedron,second,secondErrors);
    for(unsigned n = 0; n < values.size(); ++n) {
      errors[n] += secondErrors[n]+std::abs(values[n]-second[n]);
      if(errors[n] > (n < values.size()/2U ? data.absoluteH : data.absoluteD)+data.relative*std::abs(values[n])) throw std::runtime_error("Independent 3D endpoint projections do not agree within tolerance.");
    }
    data.evaluations += independent.evaluations;
  }

  static void Output(const SurfaceEndpointData& data, long double scale, const std::vector<long double>& values, const std::vector<long double>& errors, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, std::vector<double>* estimate) {
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
