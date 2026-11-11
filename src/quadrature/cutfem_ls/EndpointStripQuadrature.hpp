#ifndef FEMUS_ENDPOINT_STRIP_QUADRATURE_HPP
#define FEMUS_ENDPOINT_STRIP_QUADRATURE_HPP

#include "LineIntegration.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace FELevelSetMoments {

struct EndpointStripStatistics {
  unsigned long long stripSamples = 0;
  unsigned long long attempts = 0;
  unsigned long long completed = 0;
  unsigned long long rejected = 0;
  unsigned long long classifications = 0;
  unsigned long long bulkSamples = 0;
  unsigned long long inactiveComponents = 0;
  unsigned long long geometryBuilds = 0;
  unsigned long long geometryHits = 0;
  unsigned directions[2] = {0,0};
  std::string failure;
};

struct EndpointStripInterval {
  long double left = 0;
  long double right = 1;
  bool leftTangent = false;
  bool rightTangent = false;
  int mapping = 0;
  long double endpoint = 0;
  long double delta = 0;
  long double discriminantConstant = 0;
  std::vector<long double> reducedDiscriminant;
  std::vector<int> signs;
  std::vector<int> boundary;
  std::vector<long double> positive;
  long double phase = 0;
};

struct EndpointStripDirection {
  std::vector<std::vector<long double>> coefficient;
  std::vector<std::vector<long double>> field;
  std::vector<std::vector<std::vector<long double>>> gradient;
  std::vector<std::vector<std::vector<long double>>> coordinates;
  std::vector<std::vector<double>> inverseJacobian;
  long double metricScale = 1;
  std::vector<EndpointStripInterval> intervals;
  bool lowerZero = false;
  bool upperZero = false;
  bool fullBoundaryDelta = false;
  unsigned triangle = 0;
  std::array<long double,5> parent = {0,1,0,1,1};
  std::vector<unsigned char> requested;
  unsigned order = 1;
  long double score = std::numeric_limits<long double>::infinity();
};

struct EndpointCurveContext {
  const EndpointStripDirection* data = nullptr;
  EndpointStripInterval interval;
  unsigned degree = 0;
  unsigned* evaluations = nullptr;
};

class EndpointStripQuadrature {
public:
  // Integrates the exact curve with selective quadrature after quadratic endpoint maps.
  static EndpointStripStatistics& Statistics() { static thread_local EndpointStripStatistics statistics; return statistics; }

  static unsigned& PreferredDirection() { static thread_local unsigned direction = 0; return direction; }

  static bool IsGraph(const std::vector<std::vector<long double>>& coefficient, bool triangle) {
    for(unsigned direction = 0; direction < (triangle ? 3U : 2U); ++direction) {
      std::vector<std::vector<long double>> derivative(3,std::vector<long double>(3,0));
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) { if(i != 0 && direction != 1) derivative[i-1U][j] += i*coefficient[i][j]; if(j != 0 && direction != 0) derivative[i][j-1U] += j*coefficient[i][j]; }
      long double minimum = std::numeric_limits<long double>::infinity();
      long double maximum = -minimum;
      if(triangle) {
        unsigned degree = 0;
        for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) if(derivative[i][j] != 0) degree = std::max(degree,i+j);
        for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
          long double value = 0;
          for(unsigned a = 0; a <= i && a < 3; ++a) for(unsigned b = 0; b <= j && b < 3; ++b) { long double multiplier = 1; for(unsigned n = 0; n < a; ++n) multiplier *= static_cast<long double>(i-n)/(degree-n); for(unsigned n = 0; n < b; ++n) multiplier *= static_cast<long double>(j-n)/(degree-a-n); value += multiplier*derivative[a][b]; }
          minimum = std::min(minimum,value); maximum = std::max(maximum,value);
        }
      }
      else {
        long double bernstein[3][3];
        LineSliceUtilities::CanonicalBernstein(derivative,bernstein);
        for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) { minimum = std::min(minimum,bernstein[i][j]); maximum = std::max(maximum,bernstein[i][j]); }
      }
      const long double padding = 4096*std::numeric_limits<long double>::epsilon()*std::max(1.0L,std::max(std::abs(minimum),std::abs(maximum)));
      if((minimum >= -padding && maximum > padding) || (maximum <= padding && minimum < -padding)) return true;
    }
    return false;
  }

  static bool Try(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool triangle, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned& selected, const std::vector<std::vector<std::vector<long double>>>* gradient = nullptr, long double metricScale = 1, const std::vector<std::vector<std::vector<long double>>>* coordinates = nullptr, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::array<long double,5>* parent = nullptr, const std::vector<unsigned char>* requested = nullptr, unsigned requestedDirection = 0, bool fullBoundaryDelta = false, const OuterQuadratureRules& rules = OuterQuadratureRules()) {
    rules.Validate();
    if(requestedDirection > 2) throw std::invalid_argument("Invalid endpoint direction.");
    ++Statistics().attempts;
    static thread_local std::vector<std::vector<long double>> cachedCoefficient;
    static thread_local bool cachedTriangle = false;
    static thread_local std::array<EndpointStripDirection,2> cachedDirections;
    if(cachedCoefficient != coefficient || cachedTriangle != triangle) {
      for(unsigned d = 0; d < 2; ++d) {
        cachedDirections[d] = EndpointStripDirection();
        try { BuildDirection(coefficient,triangle,d+1U,nullptr,cachedDirections[d]); }
        catch(const std::runtime_error& error) { cachedDirections[d].score = std::numeric_limits<long double>::infinity(); Statistics().failure = error.what(); }
      }
      cachedCoefficient = coefficient;
      cachedTriangle = triangle;
      ++Statistics().geometryBuilds;
    }
    else ++Statistics().geometryHits;
    std::array<EndpointStripDirection,2> directions = cachedDirections;
    if(gradient != nullptr) { directions[0].gradient = *gradient; SurfaceGradientTranspose(*gradient,directions[1].gradient); }
    for(unsigned d = 0; d < 2; ++d) { directions[d].metricScale = metricScale; if(inverseJacobian != nullptr) directions[d].inverseJacobian = *inverseJacobian; }
    if(coordinates != nullptr) { directions[0].coordinates = *coordinates; SurfaceGradientTranspose(*coordinates,directions[1].coordinates); }
    for(unsigned d = 0; d < 2; ++d) { directions[d].fullBoundaryDelta = fullBoundaryDelta; if(parent != nullptr) directions[d].parent = LineIntegration::MomentCoordinates(parent,d+1U); directions[d].requested = LineIntegration::MomentRequest(degree,requested,d+1U,computeDelta); }
    const unsigned first = requestedDirection != 0 ? requestedDirection-1U : PreferredDirection() == 0 ? (directions[1].score < directions[0].score ? 1U : 0U) : PreferredDirection()-1U;
    std::string resolutionFailure;
    for(unsigned attempt = 0; attempt < (requestedDirection == 0 ? 2U : 1U); ++attempt) {
      const unsigned axis = attempt == 0 ? first : 1U-first;
      if(!std::isfinite(directions[axis].score)) continue;
      std::vector<long double> result;
      std::vector<long double> resultError;
      unsigned localEvaluations = 0;
      unsigned localAccepted = 0;
      try {
        IntegrateDirection(directions[axis],degree,computeDelta,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations-evaluations,result,resultError,localEvaluations,localAccepted,rules);
        if(axis == 1) {
          const unsigned count = result.size()/2U;
          for(unsigned i = 0; i <= degree; ++i) for(unsigned j = i+1U; j <= degree-i; ++j) for(unsigned p = 0; p < 2; ++p) { std::swap(result[p*count+LineSliceUtilities::MomentIndex(degree,i,j)],result[p*count+LineSliceUtilities::MomentIndex(degree,j,i)]); std::swap(resultError[p*count+LineSliceUtilities::MomentIndex(degree,i,j)],resultError[p*count+LineSliceUtilities::MomentIndex(degree,j,i)]); }
        }
        values.swap(result);
        errors.swap(resultError);
        evaluations += localEvaluations;
        accepted += localAccepted;
        selected = axis+1U;
        ++Statistics().completed;
        ++Statistics().directions[axis];
        return true;
      }
      catch(const OuterQuadratureResolutionError& error) { evaluations += localEvaluations; Statistics().failure = error.what(); resolutionFailure = error.what(); if(evaluations >= maximumEvaluations) break; }
      catch(const std::runtime_error& error) { evaluations += localEvaluations; Statistics().failure = error.what(); if(evaluations >= maximumEvaluations) break; }
    }
    ++Statistics().rejected;
    if(!resolutionFailure.empty() && evaluations < maximumEvaluations) throw OuterQuadratureResolutionError(resolutionFailure);
    return false;
  }

private:
  static long double Derivative(const std::vector<std::vector<long double>>& coefficient, long double x, long double y, unsigned dx, unsigned dy) {
    long double value = 0;
    for(unsigned i = dx; i < 3; ++i) for(unsigned j = dy; j < 3; ++j) {
      long double term = coefficient[i][j]*ReferencePower(x,i-dx)*ReferencePower(y,j-dy);
      for(unsigned k = 0; k < dx; ++k) term *= i-k;
      for(unsigned k = 0; k < dy; ++k) term *= j-k;
      value += term;
    }
    return value;
  }

  static std::vector<std::vector<long double>> Field(const std::vector<std::vector<long double>>& coefficient) {
    std::vector<std::vector<long double>> field(3,std::vector<long double>(3));
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) field[2U-j][i] = coefficient[i][j];
    for(unsigned n = 0; n < 3; ++n) LineSliceUtilities::Trim(field[n]);
    return field;
  }

  static void Classify(const EndpointStripDirection& data, long double x, EndpointStripInterval& interval) {
    const std::vector<long double> domain = {0,1,0,1};
    IntervalBranches(data.field,x,data.triangle,domain,data.lowerZero,data.upperZero,true,interval.signs,interval.boundary,interval.positive,interval.phase);
    ++Statistics().classifications;
  }

  static bool Tangent(const EndpointStripDirection& data, long double x) {
    const long double a = PolynomialValue(data.field[0],x);
    const long double b = PolynomialValue(data.field[1],x);
    if(a == 0) return false;
    const long double y = -b/(2*a);
    const long double padding = 4096*std::numeric_limits<long double>::epsilon();
    if(y < -padding || y > (data.triangle ? 1-x : 1)+padding) return false;
    const long double derivative = Derivative(data.coefficient,x,y,1,0);
    const long double norm = LineSliceUtilities::PolynomialNorm(data.field[0])+LineSliceUtilities::PolynomialNorm(data.field[1])+LineSliceUtilities::PolynomialNorm(data.field[2]);
    if(std::abs(derivative) <= padding*norm) throw std::runtime_error("A tangency strip requires a regular simple fold.");
    return true;
  }

  static void BuildDirection(const std::vector<std::vector<long double>>& coefficient, bool triangle, unsigned order, const std::vector<std::vector<std::vector<long double>>>* gradient, EndpointStripDirection& data) {
    data.coefficient = coefficient;
    data.triangle = triangle ? 1U : 0U;
    data.order = order;
    if(order == 2) for(unsigned i = 0; i < 3; ++i) for(unsigned j = i+1U; j < 3; ++j) std::swap(data.coefficient[i][j],data.coefficient[j][i]);
    if(LineIntegration::HasVerticalFiber(data.coefficient,1)) throw std::runtime_error("A projection contains an entire zero fiber.");
    data.field = Field(data.coefficient);
    const long double fieldNorm = std::max(LineSliceUtilities::PolynomialNorm(data.field[0]),std::max(LineSliceUtilities::PolynomialNorm(data.field[1]),LineSliceUtilities::PolynomialNorm(data.field[2])));
    if(LineSliceUtilities::PolynomialNorm(data.field[0])+LineSliceUtilities::PolynomialNorm(data.field[1]) <= 64*std::sqrt(std::numeric_limits<long double>::epsilon())*fieldNorm) throw std::runtime_error("A nearly parallel projection retains the transverse direction.");
    if(gradient != nullptr) { if(order == 1) data.gradient = *gradient; else SurfaceGradientTranspose(*gradient,data.gradient); }
    const LineIntegrationStrip strip = {0,1,0,0,1,triangle ? -1.0L : 0.0L};
    data.lowerZero = LineIntegration::BoundaryIsZero(data.coefficient,0,0);
    data.upperZero = LineIntegration::BoundaryIsZero(data.coefficient,1,strip.upperSlope);
    std::vector<long double> events;
    std::vector<long double> candidates;
    LineIntegration::CollectEvents(data.coefficient,strip,events,&candidates,false);
    LineSliceUtilities::PolynomialRoots(data.field[0],0,1,events);
    LineSliceUtilities::SortDistinct(events);
    std::vector<long double> tangencies;
    for(unsigned n = 0; n < candidates.size(); ++n) if(Tangent(data,candidates[n])) tangencies.push_back(candidates[n]);
    data.score = 0;
    for(unsigned n = 1; n < events.size(); ++n) {
      if(!(events[n] > events[n-1U])) continue;
      EndpointStripInterval interval;
      interval.left = events[n-1U];
      interval.right = events[n];
      interval.leftTangent = std::find(tangencies.begin(),tangencies.end(),interval.left) != tangencies.end();
      interval.rightTangent = std::find(tangencies.begin(),tangencies.end(),interval.right) != tangencies.end();
      const long double middle = interval.left+0.5L*(interval.right-interval.left);
      Classify(data,middle,interval);
      if(!interval.signs.empty()) {
        data.score += 1+2*(interval.leftTangent+interval.rightTangent);
        const std::vector<long double> roots = QuadraticRoots(PolynomialValue(data.field[0],middle),PolynomialValue(data.field[1],middle),PolynomialValue(data.field[2],middle));
        for(unsigned r = 0; r < roots.size(); ++r) if(roots[r] > 0 && roots[r] < (triangle ? 1-middle : 1)) {
          const long double dx = Derivative(data.coefficient,middle,roots[r],1,0);
          const long double dy = Derivative(data.coefficient,middle,roots[r],0,1);
          data.score += 0.05L*std::min(10.0L,std::abs(dx)/std::max(std::abs(dy),std::numeric_limits<long double>::min()));
        }
      }
      data.intervals.push_back(interval);
    }
  }

  static void Sample(const EndpointStripDirection& data, const EndpointStripInterval& interval, unsigned degree, const OuterIntervalPoint& point, const std::vector<unsigned char>& active, std::vector<long double>& values) {
    const unsigned count = active.size()/2U;
    const long double parameter = point.parameter;
    const long double offset = point.offset;
    const long double x = point.coordinate;
    const long double jacobian = point.jacobian;
    const long double a = PolynomialValue(data.field[0],x);
    const long double b = PolynomialValue(data.field[1],x);
    const long double c = PolynomialValue(data.field[2],x);
    const long double discriminant = std::fma(-4*a,c,b*b);
    const long double reduced = interval.mapping == 0 ? 0 : interval.discriminantConstant/(parameter*parameter)+interval.mapping*interval.delta*PolynomialValue(interval.reducedDiscriminant,offset);
    if(interval.mapping != 0 && !(reduced > 0)) throw std::runtime_error("An endpoint map requires a positive simple-fold discriminant.");
    const long double regularDerivative = interval.mapping == 0 ? 0 : std::sqrt(reduced);
    const long double derivative = interval.mapping == 0 ? (a == 0 ? std::abs(b) : std::sqrt(std::max(0.0L,discriminant))) : parameter*regularDerivative;
    const long double surfaceRatio = interval.mapping == 0 ? 1.0L/derivative : (2*interval.delta/regularDerivative)*(jacobian/(2*interval.delta*parameter));
    const long double upper = data.triangle ? 1-x : 1;
    const long double parentX = std::fma(data.parent[1],x,data.parent[0]);
    std::vector<long double> roots;
    std::vector<long double> cuts(1,0);
    for(unsigned r = 0; r < interval.signs.size(); ++r) {
      long double root = interval.boundary[r] < 0 ? 0 : upper;
      if(interval.boundary[r] == 0) {
        const long double denominator = b+interval.signs[r]*derivative;
        root = a == 0 ? -c/b : (interval.signs[r]*b >= 0 && denominator != 0 ? -2*c/denominator : (-b+interval.signs[r]*derivative)/(2*a));
      }
      if(!std::isfinite(root) || root < -1e-12L || root > upper+1e-12L) throw std::runtime_error("Cached root topology is inconsistent with a quadrature node.");
      roots.push_back(root);
      if(interval.boundary[r] == 0) cuts.push_back(root);
    }
    cuts.push_back(upper);
    std::vector<long double> h(degree+1U,0);
    std::vector<long double> d(degree+1U,0);
    for(unsigned r = 0; r < interval.positive.size(); ++r) if(interval.positive[r] != 0) {
      const long double parentLeft = std::fma(data.parent[3],cuts[r],data.parent[2]);
      const long double parentRight = std::fma(data.parent[3],cuts[r+1U],data.parent[2]);
      long double power = 1;
      long double geometric = 0;
      for(unsigned j = 0; j <= degree; ++j) { geometric = power+parentLeft*geometric; h[j] += (cuts[r+1U]-cuts[r])*geometric/(j+1U); power *= parentRight; }
    }
    for(unsigned r = 0; r < roots.size(); ++r) {
      if(!(derivative > 0)) throw std::runtime_error("A bulk node reached a singular fiber.");
      long double weight = (interval.boundary[r] == 0 || data.fullBoundaryDelta ? 1.0L : 0.5L)*surfaceRatio;
      if(!data.gradient.empty()) weight *= SurfaceGradientNorm(data.gradient,x,roots[r]);
      const long double parentRoot = std::fma(data.parent[3],roots[r],data.parent[2]);
      long double power = 1;
      for(unsigned j = 0; j <= degree; ++j) { d[j] += weight*power; power *= parentRoot; }
    }
    values.assign(active.size(),0);
    long double power = 1;
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned n = LineSliceUtilities::MomentIndex(degree,i,j);
        if(active[n]) values[n] = data.parent[4]*jacobian*power*h[j]; else ++Statistics().inactiveComponents;
        if(active[count+n]) values[count+n] = data.parent[4]*power*d[j]; else ++Statistics().inactiveComponents;
      }
      power *= parentX;
    }
    if(interval.mapping == 0) ++Statistics().bulkSamples; else ++Statistics().stripSamples;
  }

  static void CurveSample(EndpointCurveContext& context, const OuterIntervalPoint& point, const std::vector<unsigned char>& active, std::vector<long double>& values, std::vector<long double>& errors, long double absoluteH, long double absoluteD, long double relative) {
    ++*context.evaluations;
    Sample(*context.data,context.interval,context.degree,point,active,values);
    errors.assign(values.size(),0);
  }

  static void IntegrateDirection(const EndpointStripDirection& data, unsigned degree, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, const OuterQuadratureRules& rules = OuterQuadratureRules()) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    values.assign(2U*count,0);
    errors.assign(values.size(),0);
    const long double intervalBudget = 1.0L/std::max(1U,static_cast<unsigned>(data.intervals.size()));
    const LineIntegrationStrip strip = {0,1,0,0,1,data.triangle ? -1.0L : 0.0L};
    for(unsigned k = 0; k < data.intervals.size(); ++k) {
      const EndpointStripInterval& interval = data.intervals[k];
      if(interval.signs.empty()) { if(interval.phase > 0) LineIntegration::FullStripMoments(degree,strip,interval.left,interval.right,values,errors,&data.parent,&data.requested); ++accepted; continue; }
      std::vector<unsigned char> active = data.requested;
      if(!computeDelta) std::fill(active.begin()+count,active.end(),0);
      const unsigned leftExponent = interval.leftTangent || LineIntegration::NearEndpoint(data.coefficient,interval.left,1,interval.right-interval.left) ? 2U : 1U;
      const unsigned rightExponent = interval.rightTangent || LineIntegration::NearEndpoint(data.coefficient,interval.right,-1,interval.right-interval.left) ? 2U : 1U;
      const std::vector<OuterInterval> plan = OuterIntervalQuadrature::Plan(interval.left,interval.right,leftExponent,rightExponent,2,interval.leftTangent,interval.rightTangent);
      for(unsigned part = 0; part < plan.size(); ++part) {
        EndpointCurveContext context;
        context.data = &data;
        context.interval = interval;
        context.degree = degree;
        context.evaluations = &evaluations;
        context.interval.mapping = plan[part].mapping;
        context.interval.endpoint = plan[part].endpoint;
        context.interval.delta = plan[part].delta;
        const LineIntegrationDiscriminant discriminant = LineIntegration::EndpointDiscriminant(data.coefficient,plan[part]);
        context.interval.discriminantConstant = discriminant.constant;
        context.interval.reducedDiscriminant = discriminant.coefficient;
        OuterIntervalQuadrature::Adaptive(context,CurveSample,plan[part],plan[part].mapping == 0 ? plan[part].left : 0,plan[part].mapping == 0 ? plan[part].right : 1,intervalBudget*absoluteH/plan.size(),intervalBudget*absoluteD/plan.size(),relative,0,maximumDepth,maximumEvaluations,evaluations,accepted,2,active,values,errors,rules);
      }
    }
    for(unsigned n = 0; n < values.size(); ++n) if(data.requested[n] && (!std::isfinite(values[n]) || !std::isfinite(errors[n]) || errors[n] > (n < count ? absoluteH : absoluteD)+relative*std::abs(values[n]))) throw std::runtime_error("Accumulated non-graph moment error exceeds tolerance.");
  }
};

}

#endif
