#ifndef FEMUS_PREPARED_SURFACE_ENDPOINT_INTEGRATION_HPP
#define FEMUS_PREPARED_SURFACE_ENDPOINT_INTEGRATION_HPP


#include "PlanarEndpointIntegration.hpp"
#include "Hex27LevelSetMoments.hpp"
#include <array>
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
  const CutFemLevelSetData* prepared = nullptr;
  const std::vector<std::vector<long double>>* physicalMetric = nullptr;
  unsigned metricDegree = 0;
  bool referenceSurface = false;
  long double normalization = 1;
  CutFemLevelSetData slice;
  PreparedPlane plane;
  PreparedStatistics* statistics = nullptr;
  std::vector<std::vector<long double>> sliceMetric;
  std::vector<std::vector<std::vector<long double>>> slicePolicy;
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
  bool wedge = false;
  std::vector<unsigned> zeroFaces;
  long double domainVolume = 1.0L;
  std::vector<unsigned char> requested;
  std::array<unsigned,3> innerCounts = {0,0,0};
  std::array<unsigned,3> outerCounts = {0,0,0};
  bool computeDelta = true;
  bool computePositive = true;
  bool inverseSurface = false;
  int normalComponent = -1;
  bool retainUncertified = false;
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

  static int Phase(const SurfaceEndpointData& data, const SurfaceEndpointDomain& cell) {
    if(data.computeDelta) for(unsigned n = 0; n < data.zeroFaces.size(); ++n) if(data.zeroFaces[n] != 0) return 0;
    return data.prepared->GetPhase();
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
      long double value = cell.wedge ? ReferenceTriangleMoment(data.powers[n][0],data.powers[n][1])*IntegralPower(0,1,data.powers[n][2]) : cell.tetrahedron ? TetMoment(cell,data.powers[n][0],data.powers[n][1],data.powers[n][2]) : IntegralPower(cell.vertices[0][0],cell.vertices[1][0],data.powers[n][0])*IntegralPower(cell.vertices[0][1],cell.vertices[1][1],data.powers[n][1])*IntegralPower(cell.vertices[0][2],cell.vertices[1][2],data.powers[n][2]);
      value *= multiplier;
      values[n] += value;
      errors[n] += (64.0L+16.0L*data.degree)*std::numeric_limits<long double>::epsilon()*std::abs(value);
    }
  }

  static unsigned Direction(const SurfaceEndpointData& data, const SurfaceEndpointDomain& cell, bool& weak, long double* energy) {
    unsigned selected = 3;
    long double best = -1.0L;
    weak = false;
    for(unsigned axis = 0; axis < 3; ++axis) {
      const CutFemLevelSetData::GraphCertificate& certificate = data.prepared->GetGraphCertificate(axis);
      const long double minimum = certificate.bounds[0], maximum = certificate.bounds[1];
      energy[axis] = std::max(std::abs(minimum),std::abs(maximum));
      const long double lower = minimum > 0 ? minimum : (maximum < 0 ? -maximum : 0);
      if(!certificate.strict || lower == 0 || axis == data.preferredOuter) continue;
      const long double score = axis == data.preferredInner ? std::numeric_limits<long double>::infinity() : (lower > 0 ? 2+lower/std::max(energy[axis],std::numeric_limits<long double>::min()) : energy[axis]/(1+energy[axis]));
      if(score > best) {
        best = score;
        selected = axis;
        weak = lower == 0.0L;
      }
    }
    return selected;
  }

  static std::vector<long double> RestrictHexHalf(const std::vector<long double>& coefficients, unsigned degree, unsigned axis, unsigned side) {
    const unsigned size = degree+1, stride = axis == 0 ? size*size : (axis == 1 ? size : 1);
    std::vector<long double> result(coefficients.size());
    for(unsigned i = 0; i < coefficients.size(); ++i) {
      const unsigned power = i/stride%size;
      result[i] = std::ldexp(coefficients[i],-static_cast<int>(power));
      if(side == 1) for(unsigned k = power+1; k <= degree; ++k) result[i] += ReferenceBinomial(k,power)*std::ldexp(coefficients[i+(k-power)*stride],-static_cast<int>(k));
    }
    return result;
  }

  // Definite transverse Hessians admit extrema but exclude saddle points of the slicing height.
  static bool DefiniteHexSlices(const CutFemLevelSetData& field, long double padding) {
    std::array<long double,2> hessian[3][3];
    for(unsigned i = 0; i < 3; ++i) {
      CutFemLevelSetData derivative;
      derivative.ResetPowerCoefficients(HEX,field.GetDerivativeCoefficients(i));
      for(unsigned j = i; j < 3; ++j) hessian[i][j] = hessian[j][i] = derivative.GetDerivativeBounds(j);
    }
    for(unsigned outer = 0; outer < 3; ++outer) {
      const unsigned i = (outer+1)%3, j = (outer+2)%3;
      const long double mixed = std::max(std::abs(hessian[i][j][0]),std::abs(hessian[i][j][1]))+padding;
      if(std::min(hessian[i][i][0],hessian[j][j][0]) <= mixed && std::max(hessian[i][i][1],hessian[j][j][1]) >= -mixed) return false;
    }
    return true;
  }

  // Store the split axes in preorder, with 3 marking an uncut box or a box with a signed derivative.
  static bool PlanHexGraphs(const CutFemLevelSetData& field, unsigned remainingDepth, long double padding, std::vector<unsigned>& plan) {
    const std::array<long double,2> value = field.GetBounds();
    if(value[0] > padding || value[1] < -padding) { plan.push_back(3); return true; }
    unsigned axis = 0;
    long double variation = -1;
    for(unsigned d = 0; d < 3; ++d) {
      const std::array<long double,2> bounds = field.GetDerivativeBounds(d);
      if(bounds[0] > padding || bounds[1] < -padding) { plan.push_back(3); return true; }
      const long double spread = bounds[1]-bounds[0];
      if(spread > variation) { variation = spread; axis = d; }
    }
    if(remainingDepth == 0) return false;
    plan.push_back(axis);
    for(unsigned side = 0; side < 2; ++side) {
      CutFemLevelSetData child;
      child.ResetPowerCoefficients(HEX,RestrictHexHalf(field.GetPowerCoefficients(),2,axis,side));
      if(!PlanHexGraphs(child,remainingDepth-1,padding,plan)) return false;
    }
    return true;
  }

  static bool Independent(const SurfaceEndpointData& data, unsigned axis) {
    if(data.prepared->DependsOn(axis)) return false;
    for(unsigned n = 0; n < data.gradient.size(); ++n) {
      const unsigned size = SurfacePolynomialCubeSize(data.gradient[n].size());
      for(unsigned i = 0; i < size; ++i) for(unsigned j = 0; j < size; ++j) for(unsigned k = 0; k < size; ++k) {
        const unsigned powers[3] = {i,j,k};
        if(powers[axis] != 0 && data.gradient[n][size*size*i+size*j+k] != 0.0L) return false;
      }
    }
    return true;
  }

  static void Sample(SurfaceEndpointData& data, const SurfaceEndpointDomain& cell, unsigned outer, long double coordinate, const std::vector<unsigned char>& active, std::vector<long double>& values, std::vector<long double>& errors, long double allocation) {
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
    data.prepared->Restrict(outer,coordinate,data.slice,data.referenceSurface && data.physicalMetric == nullptr);
    if(data.statistics != nullptr) ++data.statistics->surfaceRestrictions;
    const bool triangleSlice = data.slice.GetGeometry() == TRI;
    if(!(data.slice.GetParentMeasureScale() > 0)) throw OuterQuadratureResolutionError("A prepared 3D slice reached coordinate resolution.");
    if(data.physicalMetric != nullptr) { data.prepared->RestrictComponents(outer,coordinate,data.metricDegree,*data.physicalMetric,data.sliceMetric); if(data.statistics != nullptr) ++data.statistics->metricRestrictions; }
    data.plane.Bind(data.slice,data.normalization,data.referenceSurface,data.inverseSurface,data.normalComponent,data.physicalMetric == nullptr ? nullptr : &data.sliceMetric,data.metricDegree,data.statistics);
    const std::vector<std::vector<long double>> chart = PlanarCoefficients(data.slice,data.normalization);
    if(data.referenceSurface) {
      const unsigned components = data.physicalMetric == nullptr ? data.slice.GetRootDimension() : data.sliceMetric.size(), size = data.physicalMetric == nullptr ? 3U : data.metricDegree+1;
      data.slicePolicy.resize(components);
      for(unsigned component = 0; component < components; ++component) {
        const std::vector<long double>& field = data.physicalMetric == nullptr ? data.slice.GetRootDerivativeCoefficients(component) : data.sliceMetric[component];
        data.slicePolicy[component].resize(size);
        for(unsigned i = 0; i < size; ++i) { data.slicePolicy[component][i].resize(size); for(unsigned j = 0; j < size; ++j) data.slicePolicy[component][i][j] = field[size*i+j]; }
      }
    }
    const std::vector<unsigned> boundary = data.plane.Boundaries();
    const std::vector<long double>& origin = data.slice.GetParentOrigin();
    const std::vector<std::vector<long double>>& map = data.slice.GetParentMap();
    std::array<long double,5> parent = {origin[first],map[first][0],origin[second],map[second][1],data.slice.GetParentMeasureScale()};
    long double parentOuter = coordinate;
    if(!cell.tetrahedron && !cell.wedge && data.domainVolume != 1) {
      const std::array<long double,3>& lower = cell.vertices[0];
      const std::array<long double,3>& upper = cell.vertices[1];
      parent = {lower[first],upper[first]-lower[first],lower[second],upper[second]-lower[second],data.domainVolume};
      parentOuter = std::fma(upper[outer]-lower[outer],coordinate,lower[outer]);
    }
    std::vector<long double> moments;
    std::vector<long double> momentErrors;
    unsigned evaluations = 0;
    unsigned accepted = 0;
    unsigned selected = 0;
    try {
      if(triangleSlice) ++Statistics().triCalls; else ++Statistics().quadCalls;
      PlanarEndpointIntegration::Integrate(chart,degree,triangleSlice,data.computeDelta,0.01L*allocation*data.absoluteH,0.01L*allocation*data.absoluteD,0.01L*data.relative,data.maximumDepth,std::max(3U,data.maximumEvaluations > data.evaluations ? data.maximumEvaluations-data.evaluations : 0U),moments,momentErrors,evaluations,accepted,selected,data.referenceSurface ? &data.slicePolicy : nullptr,&parent,&required,requestedDirection,false,&boundary,data.quadratureRules,data.inverseSurface,data.normalComponent,&data.plane);
    }
    catch(...) { data.evaluations += evaluations; throw; }
    data.evaluations += evaluations;
    values.assign(2U*count,0);
    errors.assign(2U*count,0);
    for(unsigned n = 0; n < count; ++n) {
      if(!active[n] && !active[count+n]) { ++Statistics().skippedComponents; continue; }
      const unsigned source = ReferenceMomentIndex2D(degree,data.powers[n][first],data.powers[n][second]);
      const long double power = ReferencePower(parentOuter,data.powers[n][outer]);
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
        const std::vector<long double>& derivative = data.prepared->GetDerivativeCoefficients(axis);
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
    std::vector<long double> best(total.size(),std::numeric_limits<long double>::quiet_NaN());
    std::vector<long double> bestErrors(total.size(),std::numeric_limits<long double>::infinity());
    bool available = false;
    bool accuracyFailure = false;
    unsigned accuracyOuter = 3;
    for(unsigned attempt = 0; attempt < axes.size(); ++attempt) {
      const unsigned outer = axes[attempt];
      if(data.preferredOuter < 3 && outer != data.preferredOuter) continue;
      std::vector<long double> values(total.size(),0);
      std::vector<long double> errors(total.size(),0);
      const unsigned long long previousWarnings = QuadratureWarningCount();
      try {
        QuadratureAccuracyRetryScope accuracyScope(!data.retainUncertified);
        ++data.outerCounts[outer];
        if(inner < 3) ++data.innerCounts[inner];
        if(!cell.tetrahedron && (!cell.wedge || outer == 2) && Independent(data,outer)) {
          std::vector<unsigned char> dependencies = active;
          for(unsigned n = 0; n < data.powers.size(); ++n) {
            unsigned powers[3] = {data.powers[n][0],data.powers[n][1],data.powers[n][2]};
            powers[outer] = 0;
            for(unsigned pair = 0; pair < 2; ++pair) if(active[pair*data.powers.size()+n]) dependencies[pair*data.powers.size()+HexMomentIndex(powers[0],powers[1],powers[2])] = 1;
          }
          Sample(data,cell,outer,0,dependencies,values,errors,data.domainVolume/2);
          const std::vector<long double> source = values;
          const std::vector<long double> sourceError = errors;
          for(unsigned n = 0; n < data.powers.size(); ++n) {
            unsigned powers[3] = {data.powers[n][0],data.powers[n][1],data.powers[n][2]};
            const long double factor = IntegralPower(cell.vertices[0][outer],cell.vertices[1][outer],powers[outer])/(cell.vertices[1][outer]-cell.vertices[0][outer]);
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
          ++Statistics().eventBuilds;
          const std::vector<EndpointProjectionEvent> events = EndpointProjection::Build(data.coefficient,cell.tetrahedron,outer,inner < 3 ? inner : data.preferredInner,inner < 3,cell.wedge,data.prepared);
          const long double budget = data.domainVolume/(events.size()-1U);
          for(unsigned n = 1; n < events.size(); ++n) {
            const long double left = events[n-1U].point;
            const long double right = events[n].point;
            if(!(right > left)) continue;
            const std::vector<OuterInterval> plan = OuterIntervalQuadrature::Plan(left,right,events[n-1U].exponent,events[n].exponent,3);
            for(unsigned part = 0; part < plan.size(); ++part) {
              Adaptive(data,cell,outer,plan[part].mapping == 0 ? plan[part].left : 0,plan[part].mapping == 0 ? plan[part].right : 1,plan[part].mapping,plan[part].endpoint,plan[part].delta,plan[part].exponent,budget/plan.size(),0,active,values,errors);
              if(plan[part].mapping != 0) ++Statistics().endpointStrips;
            }
          }
        }
        KeepBestQuadrature(values,errors,active,best,bestErrors);
        available = true;
        if(QuadratureWarningCount() != previousWarnings && data.evaluations < data.maximumEvaluations && data.preferredOuter >= 3 && attempt+1U < axes.size()) continue;
        const bool certified = QuadratureWarningCount() == previousWarnings;
        for(unsigned n = 0; n < total.size(); ++n) if(active[n]) { total[n] += certified ? values[n] : best[n]; totalErrors[n] += certified ? errors[n] : bestErrors[n]; }
        return;
      }
      catch(const QuadratureAccuracyRetry& error) { failure = std::runtime_error(error.what()); accuracyFailure = true; if(accuracyOuter >= 3) accuracyOuter = outer; if(data.evaluations >= data.maximumEvaluations) break; }
      catch(const OuterQuadratureResolutionError& error) { failure = std::runtime_error(error.what()); accuracyFailure = true; if(accuracyOuter >= 3) accuracyOuter = outer; if(data.evaluations >= data.maximumEvaluations) break; }
      catch(const std::runtime_error& error) { failure = std::runtime_error(error.what()); if(data.evaluations >= data.maximumEvaluations && !available) throw; }
    }
    if(accuracyFailure && !data.retainUncertified) {
      const unsigned previousOuter = data.preferredOuter;
      data.retainUncertified = true;
      data.preferredOuter = accuracyOuter;
      try { IntegrateParent(data,cell,total,totalErrors); }
      catch(...) { data.retainUncertified = false; data.preferredOuter = previousOuter; throw; }
      data.retainUncertified = false;
      data.preferredOuter = previousOuter;
      return;
    }
    if(available) {
      std::vector<unsigned char> remaining = active;
      RetainUncertifiedQuadrature("3D projection",std::string("alternative projections failed; ")+failure.what(),best,bestErrors,data.absoluteH*data.domainVolume,data.absoluteD*data.domainVolume,data.relative,remaining,total,totalErrors);
      return;
    }
    throw failure;
  }

  static void IntegrateHexGraphs(SurfaceEndpointData& data, const SurfaceEndpointDomain& cell, const std::vector<unsigned>& plan, std::size_t& index, std::vector<long double>& values, std::vector<long double>& errors) {
    const unsigned axis = plan[index++];
    if(axis == 3) { IntegrateParent(data,cell,values,errors); return; }
    const long double middle = (cell.vertices[0][axis]+cell.vertices[1][axis])/2;
    for(unsigned side = 0; side < 2; ++side) {
      SurfaceEndpointDomain box = cell;
      box.vertices[1-side][axis] = middle;
      CutFemLevelSetData field;
      field.ResetPowerCoefficients(HEX,RestrictHexHalf(data.prepared->GetPowerCoefficients(),2,axis,side));
      SurfaceEndpointData child = data;
      child.prepared = &field;
      child.domainVolume = data.domainVolume/2;
      child.coefficient = field.GetNormalizedPowerCoefficients();
      child.metricDegree = data.physicalMetric == nullptr ? 2U : data.metricDegree;
      for(unsigned n = 0; n < child.gradient.size(); ++n) child.gradient[n] = RestrictHexHalf(data.gradient[n],child.metricDegree,axis,side);
      child.physicalMetric = child.referenceSurface ? &child.gradient : nullptr;
      for(unsigned n = 0; n < child.zeroFaces.size(); ++n) child.zeroFaces[n] = field.IsZeroOnBoundary(n);
      IntegrateHexGraphs(child,box,plan,index,values,errors);
      data.evaluations = child.evaluations;
      data.accepted = child.accepted;
      data.innerCounts = child.innerCounts;
      data.outerCounts = child.outerCounts;
    }
  }

  static void Integrate(SurfaceEndpointData& data, bool tetrahedron, std::vector<long double>& values, std::vector<long double>& errors) {
    QuadratureWarningScope warningScope("3D endpoint integration");
    data.quadratureRules.Validate();
    const unsigned count = (data.degree+1U)*(data.degree+2U)*(data.degree+3U)/6U;
    data.powers.assign(count,{0,0,0});
    for(unsigned i = 0; i <= data.degree; ++i) for(unsigned j = 0; j <= data.degree-i; ++j) for(unsigned k = 0; k <= data.degree-i-j; ++k) data.powers[HexMomentIndex(i,j,k)] = {i,j,k};
    SurfaceEndpointDomain cell;
    cell.tetrahedron = tetrahedron;
    cell.wedge = data.wedge;
    cell.vertices = tetrahedron ? std::vector<std::array<long double,3>>{{0,0,0},{1,0,0},{0,1,0},{0,0,1}} : (data.wedge ? std::vector<std::array<long double,3>>{{0,0,0},{1,1,1}} : std::vector<std::array<long double,3>>{{0,0,0},{1,1,1}});
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
    if(!data.computePositive) std::fill(data.requested.begin(),data.requested.begin()+count,0);
    if(!data.computeDelta) std::fill(data.requested.begin()+count,data.requested.end(),0);
    std::vector<long double> best(2U*count,std::numeric_limits<long double>::quiet_NaN());
    std::vector<long double> bestErrors(2U*count,std::numeric_limits<long double>::infinity());
    std::vector<unsigned> plan;
    if(!tetrahedron && !cell.wedge && data.preferredOuter == 3 && Phase(data,cell) == 0 && data.prepared->GetGraphDirection() < 0) {
      const long double padding = 8192*std::numeric_limits<long double>::epsilon()*norm*data.prepared->GetNormalizationScale();
      if(!DefiniteHexSlices(*data.prepared,padding) && !PlanHexGraphs(*data.prepared,data.maximumDepth,padding,plan)) plan.clear();
    }
    long double budgetFactor = 1.0L;
    for(unsigned pass = 0; pass < 8; ++pass) {
      data.absoluteH = budgetFactor*absoluteH/data.domainVolume;
      data.absoluteD = budgetFactor*absoluteD/data.domainVolume;
      data.relative = budgetFactor*relative;
      std::vector<long double> passValues(2U*count,0.0L);
      std::vector<long double> passErrors(2U*count,0.0L);
      ++Statistics().integrationPasses;
      std::size_t index = 0;
      if(plan.empty()) IntegrateParent(data,cell,passValues,passErrors);
      else IntegrateHexGraphs(data,cell,plan,index,passValues,passErrors);
      KeepBestQuadrature(passValues,passErrors,data.requested,best,bestErrors);
      long double worst = 0.0L;
      bool remaining = false;
      for(unsigned n = 0; n < values.size(); ++n) if(data.requested[n]) {
        if(!std::isfinite(passValues[n])) throw std::runtime_error("A 3D graph moment is not finite.");
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
      if(pass == 7 || !std::isfinite(worst) || data.evaluations >= data.maximumEvaluations) {
        RetainUncertifiedQuadrature("3D moments",data.evaluations >= data.maximumEvaluations ? "evaluation budget exhausted before accumulated tolerance" : "accumulated error exceeds tolerance after selective retries",best,bestErrors,absoluteH,absoluteD,relative,data.requested,values,errors);
        data.absoluteH = absoluteH;
        data.absoluteD = absoluteD;
        data.relative = relative;
        data.requested.clear();
        return;
      }
      budgetFactor /= std::max(2.0L,2.0L*worst);
    }
  }

  static void Compare(SurfaceEndpointData& data, bool tetrahedron, const std::vector<long double>& values, std::vector<long double>& errors) {
    QuadratureWarningScope warningScope("Independent 3D validation");
    SurfaceEndpointData independent = data;
    independent.evaluations = independent.accepted = 0;
    unsigned outer = 0;
    for(unsigned axis = 1; axis < 3; ++axis) if(data.outerCounts[axis] > data.outerCounts[outer]) outer = axis;
    independent.preferredOuter = (outer+1U)%3U;
    if(independent.preferredOuter == data.graphAxis) independent.preferredOuter = (outer+2U)%3U;
    std::vector<long double> second;
    std::vector<long double> secondErrors;
    try { Integrate(independent,tetrahedron,second,secondErrors); }
    catch(const std::runtime_error& error) { WarnQuadrature("Independent 3D validation",std::string("validation unavailable: ")+error.what(),0,std::numeric_limits<long double>::infinity(),data.absoluteD); data.evaluations += independent.evaluations; return; }
    for(unsigned n = 0; n < values.size(); ++n) {
      errors[n] += secondErrors[n]+std::abs(values[n]-second[n]);
      if(errors[n] > (n < values.size()/2U ? data.absoluteH : data.absoluteD)+data.relative*std::abs(values[n])) WarnQuadrature("Independent 3D validation","endpoint projections disagree",n,errors[n],(n < values.size()/2U ? data.absoluteH : data.absoluteD)+data.relative*std::abs(values[n]));
    }
    data.evaluations += independent.evaluations;
  }


};

}



#endif
