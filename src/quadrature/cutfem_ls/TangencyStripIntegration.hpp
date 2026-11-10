#ifndef FEMUS_TANGENCY_STRIP_INTEGRATION_HPP
#define FEMUS_TANGENCY_STRIP_INTEGRATION_HPP

#include "LineIntegration.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace FELevelSetMoments {

struct TangencyStripStatistics {
  unsigned long long attempts = 0;
  unsigned long long completed = 0;
  unsigned long long rejected = 0;
  unsigned long long classifications = 0;
  unsigned long long analyticPanels = 0;
  unsigned long long modelRefinements = 0;
  unsigned long long bulkSamples = 0;
  unsigned long long inactiveComponents = 0;
  unsigned long long rule15 = 0;
  unsigned long long rule31 = 0;
  unsigned long long tangencies = 0;
  unsigned long long geometryBuilds = 0;
  unsigned long long geometryHits = 0;
  unsigned directions[2] = {0,0};
  std::string failure;
  long double minimumDelta = 1;
  long double maximumDelta = 0;
};

struct TangencyStripInterval {
  long double left = 0;
  long double right = 1;
  bool leftTangent = false;
  bool rightTangent = false;
  std::vector<int> signs;
  std::vector<int> boundary;
  std::vector<long double> positive;
  long double phase = 0;
};

struct TangencyStripDirection {
  std::vector<std::vector<long double>> coefficient;
  std::vector<std::vector<long double>> field;
  std::vector<std::vector<std::vector<long double>>> gradient;
  std::vector<std::vector<std::vector<long double>>> coordinates;
  std::vector<std::vector<double>> inverseJacobian;
  long double metricScale = 1;
  std::vector<TangencyStripInterval> intervals;
  unsigned triangle = 0;
  unsigned order = 1;
  long double score = std::numeric_limits<long double>::infinity();
};

class TangencyStripIntegration {
public:
  // Non-graph unit cells use cached root topology, quadratic endpoint models, and selective GK in the regular intervals.
  static TangencyStripStatistics& Statistics() { static thread_local TangencyStripStatistics statistics; return statistics; }

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

  static bool Try(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool triangle, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned& selected, const std::vector<std::vector<std::vector<long double>>>* gradient = nullptr, long double metricScale = 1, const std::vector<std::vector<std::vector<long double>>>* coordinates = nullptr, const std::vector<std::vector<double>>* inverseJacobian = nullptr) {
    if(degree > 10) return false;
    ++Statistics().attempts;
    static thread_local std::vector<std::vector<long double>> cachedCoefficient;
    static thread_local bool cachedTriangle = false;
    static thread_local std::array<TangencyStripDirection,2> cachedDirections;
    if(cachedCoefficient != coefficient || cachedTriangle != triangle) {
      for(unsigned d = 0; d < 2; ++d) {
        cachedDirections[d] = TangencyStripDirection();
        try { BuildDirection(coefficient,triangle,d+1U,nullptr,cachedDirections[d]); }
        catch(const std::runtime_error& error) { cachedDirections[d].score = std::numeric_limits<long double>::infinity(); Statistics().failure = error.what(); }
      }
      cachedCoefficient = coefficient;
      cachedTriangle = triangle;
      ++Statistics().geometryBuilds;
    }
    else ++Statistics().geometryHits;
    std::array<TangencyStripDirection,2> directions = cachedDirections;
    if(gradient != nullptr) { directions[0].gradient = *gradient; SurfaceGradientTranspose(*gradient,directions[1].gradient); }
    for(unsigned d = 0; d < 2; ++d) { directions[d].metricScale = metricScale; if(inverseJacobian != nullptr) directions[d].inverseJacobian = *inverseJacobian; }
    if(coordinates != nullptr) { directions[0].coordinates = *coordinates; SurfaceGradientTranspose(*coordinates,directions[1].coordinates); }
    const unsigned first = directions[1].score < directions[0].score ? 1U : 0U;
    for(unsigned attempt = 0; attempt < 2; ++attempt) {
      const unsigned axis = attempt == 0 ? first : 1U-first;
      if(!std::isfinite(directions[axis].score)) continue;
      std::vector<long double> result;
      std::vector<long double> resultError;
      unsigned localEvaluations = 0;
      unsigned localAccepted = 0;
      try {
        IntegrateDirection(directions[axis],degree,computeDelta,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,result,resultError,localEvaluations,localAccepted);
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
      catch(const std::runtime_error& error) { Statistics().failure = error.what(); }
    }
    ++Statistics().rejected;
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

  static void Classify(const TangencyStripDirection& data, long double x, TangencyStripInterval& interval) {
    const std::vector<long double> domain = {0,1,0,1};
    IntervalBranches(data.field,x,data.triangle,domain,false,false,true,interval.signs,interval.boundary,interval.positive,interval.phase);
    ++Statistics().classifications;
  }

  static bool Tangent(const TangencyStripDirection& data, long double x) {
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

  static void BuildDirection(const std::vector<std::vector<long double>>& coefficient, bool triangle, unsigned order, const std::vector<std::vector<std::vector<long double>>>* gradient, TangencyStripDirection& data) {
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
    if(LineIntegration::BoundaryIsZero(data.coefficient,0,0) || LineIntegration::BoundaryIsZero(data.coefficient,1,strip.upperSlope)) throw std::runtime_error("A zero boundary retains the existing trace integrator.");
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
      TangencyStripInterval interval;
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

  static std::vector<std::vector<long double>> QuadraticModel(const TangencyStripDirection& data, long double endpoint) {
    const long double x = endpoint;
    const long double y = -PolynomialValue(data.field[1],x)/(2*PolynomialValue(data.field[0],x));
    const long double dx = Derivative(data.coefficient,x,y,1,0);
    const long double yy = 0.5L*Derivative(data.coefficient,x,y,0,2);
    std::vector<std::vector<long double>> model(3,std::vector<long double>(3,0));
    model[0][2] = yy;
    model[1][0] = dx;
    model[0][1] = -2*yy*y;
    model[0][0] = -dx*x+yy*y*y;
    return model;
  }

  static void ModelMetric(const TangencyStripDirection& data, const std::vector<std::vector<long double>>& coefficient, std::vector<std::vector<std::vector<long double>>>& gradient) {
    if(data.gradient.empty()) return;
    gradient.assign(2,std::vector<std::vector<long double>>(3,std::vector<long double>(3,0)));
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) { if(i != 0) gradient[0][i-1U][j] = data.metricScale*i*coefficient[i][j]; if(j != 0) gradient[1][i][j-1U] = data.metricScale*j*coefficient[i][j]; }
    if(!data.coordinates.empty()) SurfaceGradientCofactor2D(gradient,data.coordinates,1,1);
    else if(!data.inverseJacobian.empty()) { if(data.order == 2) gradient[0].swap(gradient[1]); SurfaceGradientPhysical(gradient,data.inverseJacobian); }
  }

  static bool ProductIntegral(const std::vector<long double>& first, const std::vector<long double>& second, long double& value, long double& error) {
    const unsigned order = first.size()-1U;
    long double norm = 0;
    long double tail = 0;
    long double previous = 0;
    value = 0;
    for(unsigned n = 0; n <= order; ++n) {
      long double coefficient = 0;
      for(unsigned k = 0; k <= n; ++k) coefficient += first[k]*second[n-k];
      const long double term = coefficient/(n+1U);
      value += term;
      norm += std::abs(term);
      if(n+8U > order) tail += std::abs(term);
      else if(n+16U > order) previous += std::abs(term);
    }
    const long double rounding = 256*std::numeric_limits<long double>::epsilon()*norm;
    error = 16*tail+rounding;
    return std::isfinite(value) && std::isfinite(error) && (tail <= 0.5L*previous || tail <= rounding);
  }

  static bool AnalyticStrip(const TangencyStripDirection& data, const TangencyStripInterval& interval, const std::vector<std::vector<long double>>& coefficient, unsigned degree, long double endpoint, long double delta, bool left, const std::vector<unsigned char>& active, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<std::vector<long double>>>* modelMetric = nullptr) {
    const unsigned order = std::max(24U,2U*degree+16U);
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    const std::vector<std::vector<long double>> field = Field(coefficient);
    std::vector<long double> x = SeriesConstant(endpoint,order);
    x[2] = left ? delta : -delta;
    std::vector<long double> discriminant = DiscriminantPolynomial(field);
    const long double criticalResidual = std::abs(PolynomialValue(discriminant,endpoint));
    const long double criticalSlope = std::abs(PolynomialValue(PolynomialDerivative(discriminant),endpoint));
    if(!(criticalSlope > 0)) return false;
    long double uncertainty = criticalResidual/criticalSlope;
    std::vector<long double> reduced = DivideLinear(discriminant,endpoint);
    if(!left) for(unsigned n = 0; n < reduced.size(); ++n) reduced[n] = -reduced[n];
    std::vector<long double> regular = SeriesSquareRoot(SeriesPolynomial(reduced,x));
    std::vector<long double> squareRoot(order+1U,0);
    for(unsigned n = 1; n <= order; ++n) squareRoot[n] = std::sqrt(delta)*regular[n-1U];
    std::vector<long double> a = SeriesPolynomial(field[0],x);
    std::vector<long double> b = SeriesPolynomial(field[1],x);
    if(coefficient != data.coefficient) {
      const long double y = -PolynomialValue(data.field[1],endpoint)/(2*PolynomialValue(data.field[0],endpoint));
      const long double dx = Derivative(data.coefficient,endpoint,y,1,0);
      const long double xx = coefficient[2][0];
      const long double xy = coefficient[1][1];
      const long double yy = coefficient[0][2];
      std::vector<long double> reducedModel = SeriesConstant((left ? 1.0L : -1.0L)*(-4*yy*dx),order);
      reducedModel[2] = delta*(xy*xy-4*yy*xx);
      regular = SeriesSquareRoot(reducedModel);
      for(unsigned n = 1; n <= order; ++n) squareRoot[n] = std::sqrt(delta)*regular[n-1U];
      a = SeriesConstant(yy,order);
      b = SeriesConstant(-2*yy*y,order);
      b[2] = (left ? delta : -delta)*xy;
      uncertainty = 0;
    }
    if(a[0] == 0) return false;
    for(unsigned n = 0; n <= order; ++n) a[n] *= 2;
    TangencyStripDirection model = data;
    model.coefficient = coefficient;
    model.field = field;
    TangencyStripInterval topology = interval;
    if(coefficient != data.coefficient) Classify(model,endpoint+(left ? 0.5L*delta : -0.5L*delta),topology);
    std::vector<long double> endpointEvents;
    const LineIntegrationStrip strip = {left ? endpoint : endpoint-delta,left ? endpoint+delta : endpoint,0,0,1,data.triangle ? -1.0L : 0.0L};
    if(coefficient != data.coefficient) LineIntegration::CollectEvents(coefficient,strip,endpointEvents,nullptr,false);
    const long double padding = 4096*std::numeric_limits<long double>::epsilon()*std::max(1.0L,std::abs(endpoint));
    for(unsigned n = 0; n < endpointEvents.size(); ++n) if(endpointEvents[n] > strip.left+padding && endpointEvents[n] < strip.right-padding) return false;
    std::vector<long double> lower = SeriesConstant(0,order);
    std::vector<long double> upper = SeriesConstant(1,order);
    if(data.triangle) for(unsigned n = 0; n <= order; ++n) upper[n] -= x[n];
    std::vector<std::vector<long double>> roots;
    std::vector<std::vector<long double>> cuts(1,lower);
    for(unsigned r = 0; r < topology.signs.size(); ++r) {
      std::vector<long double> numerator(b);
      for(unsigned n = 0; n <= order; ++n) numerator[n] = -numerator[n]+topology.signs[r]*squareRoot[n];
      roots.push_back(topology.boundary[r] < 0 ? lower : (topology.boundary[r] > 0 ? upper : SeriesQuotient(numerator,a)));
      if(topology.boundary[r] == 0) cuts.push_back(roots.back());
    }
    cuts.push_back(upper);
    if(cuts.size() != topology.positive.size()+1U) return false;
    std::vector<long double> jacobian = SeriesConstant(0,order);
    jacobian[1] = 2*delta;
    std::vector<long double> common = SeriesQuotient(SeriesConstant(2*std::sqrt(delta),order),regular);
    std::vector<std::vector<long double>> weights(roots.size(),common);
    const std::vector<std::vector<std::vector<long double>>>& metric = modelMetric == nullptr ? data.gradient : *modelMetric;
    for(unsigned r = 0; r < roots.size(); ++r) {
      if(topology.boundary[r] != 0) for(unsigned n = 0; n <= order; ++n) weights[r][n] *= 0.5L;
      if(!metric.empty()) weights[r] = SeriesProduct(weights[r],SurfaceGradientSeriesNorm(metric,x,roots[r]));
    }
    std::vector<std::vector<long double>> cutPowers(cuts.size(),SeriesConstant(1,order));
    std::vector<std::vector<long double>> rootPowers(roots.size(),SeriesConstant(1,order));
    std::vector<std::vector<long double>> positive(degree+1U,SeriesConstant(0,order));
    std::vector<std::vector<long double>> deltaSlice = positive;
    for(unsigned j = 0; j <= degree; ++j) {
      for(unsigned r = 0; r < cuts.size(); ++r) cutPowers[r] = SeriesProduct(cutPowers[r],cuts[r]);
      for(unsigned r = 0; r < topology.positive.size(); ++r) if(topology.positive[r] != 0) for(unsigned n = 0; n <= order; ++n) positive[j][n] += (cutPowers[r+1U][n]-cutPowers[r][n])/(j+1U);
      positive[j] = SeriesProduct(positive[j],jacobian);
      for(unsigned r = 0; r < roots.size(); ++r) {
        const std::vector<long double> term = SeriesProduct(rootPowers[r],weights[r]);
        for(unsigned n = 0; n <= order; ++n) deltaSlice[j][n] += term[n];
        rootPowers[r] = SeriesProduct(rootPowers[r],roots[r]);
      }
    }
    values.assign(2U*count,0);
    errors.assign(values.size(),0);
    std::vector<long double> power = SeriesConstant(1,order);
    bool converged = true;
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) for(unsigned p = 0; p < 2; ++p) {
        const unsigned n = p*count+LineSliceUtilities::MomentIndex(degree,i,j);
        if(!active[n]) continue;
        const bool bounded = ProductIntegral(power,p == 0 ? positive[j] : deltaSlice[j],values[n],errors[n]);
        errors[n] += 64*(degree+1U)*uncertainty*std::abs(values[n])/delta;
        if(!bounded) errors[n] = std::numeric_limits<long double>::infinity();
      }
      power = SeriesProduct(power,x);
    }
    ++Statistics().analyticPanels;
    return converged;
  }

  static void Endpoint(const TangencyStripDirection& data, const TangencyStripInterval& interval, unsigned degree, bool left, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, std::vector<long double>& extents, std::vector<long double>& values, std::vector<long double>& errors) {
    const unsigned count = values.size()/2U;
    const long double endpoint = left ? interval.left : interval.right;
    const std::vector<std::vector<long double>> model = QuadraticModel(data,endpoint);
    std::vector<std::vector<std::vector<long double>>> modelMetric;
    ModelMetric(data,model,modelMetric);
    const bool exact = model == data.coefficient;
    const long double modelReserve = data.triangle ? 1.0L : 4*ReferencePower(3.0L,degree);
    std::vector<unsigned char> active(values.size(),1);
    if(!computeDelta) std::fill(active.begin()+count,active.end(),0);
    long double delta = 0.05L*(interval.right-interval.left);
    for(unsigned level = 0; level < 32; ++level) {
      std::vector<long double> reference;
      std::vector<long double> referenceError;
      std::vector<long double> approximate;
      std::vector<long double> approximateError;
      const bool referenceBounded = AnalyticStrip(data,interval,data.coefficient,degree,endpoint,delta,left,active,reference,referenceError);
      const bool modelBounded = exact ? referenceBounded : AnalyticStrip(data,interval,model,degree,endpoint,delta,left,active,approximate,approximateError,&modelMetric);
      if(exact) { approximate = reference; approximateError = referenceError; }
      bool any = false;
      for(unsigned n = 0; n < active.size(); ++n) if(active[n]) {
        const long double error = referenceBounded && modelBounded ? std::abs(approximate[n]-reference[n])+approximateError[n]+referenceError[n] : std::numeric_limits<long double>::infinity();
        const long double tolerance = ((n < count ? absoluteH : absoluteD)+relative*std::abs(approximate.empty() ? 0.0L : approximate[n]))/modelReserve;
        if(error <= tolerance) {
          values[n] += approximate[n]; errors[n] += error; extents[n] = delta; active[n] = 0;
          Statistics().minimumDelta = std::min(Statistics().minimumDelta,delta); Statistics().maximumDelta = std::max(Statistics().maximumDelta,delta);
        }
        else any = true;
      }
      if(!any) return;
      ++Statistics().modelRefinements;
      delta *= 0.5L;
      if(delta < 4096*std::numeric_limits<long double>::epsilon()*std::max(1.0L,std::abs(endpoint))) break;
    }
    throw std::runtime_error("A quadratic endpoint model did not reach the moment tolerance.");
  }

  static void Sample(const TangencyStripDirection& data, const TangencyStripInterval& interval, unsigned degree, long double x, const std::vector<unsigned char>& active, std::vector<long double>& values) {
    const unsigned count = active.size()/2U;
    const long double a = PolynomialValue(data.field[0],x);
    const long double b = PolynomialValue(data.field[1],x);
    const long double c = PolynomialValue(data.field[2],x);
    const long double discriminant = std::fma(-4*a,c,b*b);
    const long double derivative = a == 0 ? std::abs(b) : std::sqrt(std::max(0.0L,discriminant));
    const long double upper = data.triangle ? 1-x : 1;
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
      long double power = 1;
      long double geometric = 0;
      for(unsigned j = 0; j <= degree; ++j) { geometric = power+cuts[r]*geometric; h[j] += (cuts[r+1U]-cuts[r])*geometric/(j+1U); power *= cuts[r+1U]; }
    }
    for(unsigned r = 0; r < roots.size(); ++r) {
      if(!(derivative > 0)) throw std::runtime_error("A bulk node reached a singular fiber.");
      long double weight = (interval.boundary[r] == 0 ? 1.0L : 0.5L)/derivative;
      if(!data.gradient.empty()) weight *= SurfaceGradientNorm(data.gradient,x,roots[r]);
      long double power = 1;
      for(unsigned j = 0; j <= degree; ++j) { d[j] += weight*power; power *= roots[r]; }
    }
    values.assign(active.size(),0);
    long double power = 1;
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned n = LineSliceUtilities::MomentIndex(degree,i,j);
        if(active[n]) values[n] = power*h[j]; else ++Statistics().inactiveComponents;
        if(active[count+n]) values[count+n] = power*d[j]; else ++Statistics().inactiveComponents;
      }
      power *= x;
    }
    ++Statistics().bulkSamples;
  }

  template<unsigned Rule> static void Panel(const TangencyStripDirection& data, const TangencyStripInterval& interval, unsigned degree, long double left, long double right, const std::vector<unsigned char>& active, unsigned maximumEvaluations, unsigned& evaluations, std::vector<long double>& values, std::vector<long double>& errors) {
    const long double middle = left+0.5L*(right-left);
    const long double half = 0.5L*(right-left);
    values.assign(active.size(),0);
    errors.assign(active.size(),0);
    std::vector<long double> low(active.size(),0);
    std::vector<long double> norm(active.size(),0);
    for(unsigned node = 0; node < boost::math::quadrature::gauss_kronrod<long double,Rule>::abscissa().size(); ++node) for(unsigned side = 0; side < (node == 0 ? 1U : 2U); ++side) {
      if(evaluations >= maximumEvaluations) throw std::runtime_error("Non-graph GK evaluation budget exhausted.");
      ++evaluations;
      const long double point = middle+(side == 0 ? -half : half)*boost::math::quadrature::gauss_kronrod<long double,Rule>::abscissa()[node];
      const long double weight = half*boost::math::quadrature::gauss_kronrod<long double,Rule>::weights()[node];
      const long double gaussWeight = node%2U == 0 ? half*boost::math::quadrature::gauss<long double,(Rule-1U)/2U>::weights()[node/2U] : 0;
      std::vector<long double> sample;
      Sample(data,interval,degree,point,active,sample);
      for(unsigned n = 0; n < active.size(); ++n) if(active[n]) { values[n] += weight*sample[n]; low[n] += gaussWeight*sample[n]; norm[n] += weight*std::abs(sample[n]); }
    }
    for(unsigned n = 0; n < active.size(); ++n) if(active[n]) errors[n] = std::max(std::abs(values[n]-low[n]),256*std::numeric_limits<long double>::epsilon()*norm[n]);
    if(Rule == 15) ++Statistics().rule15; else ++Statistics().rule31;
  }

  static void Adaptive(const TangencyStripDirection& data, const TangencyStripInterval& interval, unsigned degree, long double left, long double right, long double absoluteH, long double absoluteD, long double relative, unsigned depth, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, const std::vector<unsigned char>& requested, std::vector<long double>& total, std::vector<long double>& totalError) {
    std::vector<unsigned char> active = requested;
    std::vector<long double> values;
    std::vector<long double> errors;
    Panel<15>(data,interval,degree,left,right,active,maximumEvaluations,evaluations,values,errors);
    if(!LineIntegration::Accept(values,errors,absoluteH,absoluteD,relative,active,total,totalError)) { ++accepted; return; }
    Panel<31>(data,interval,degree,left,right,active,maximumEvaluations,evaluations,values,errors);
    if(!LineIntegration::Accept(values,errors,absoluteH,absoluteD,relative,active,total,totalError)) { ++accepted; return; }
    const long double middle = left+0.5L*(right-left);
    if(depth >= maximumDepth || middle == left || middle == right) throw std::runtime_error("Non-graph GK did not reach tolerance.");
    Adaptive(data,interval,degree,left,middle,0.5L*absoluteH,0.5L*absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,active,total,totalError);
    Adaptive(data,interval,degree,middle,right,0.5L*absoluteH,0.5L*absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,active,total,totalError);
  }

  static void IntegrateDirection(const TangencyStripDirection& data, unsigned degree, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    values.assign(2U*count,0);
    errors.assign(values.size(),0);
    const long double intervalBudget = 1.0L/std::max(1U,static_cast<unsigned>(data.intervals.size()));
    const LineIntegrationStrip strip = {0,1,0,0,1,data.triangle ? -1.0L : 0.0L};
    for(unsigned k = 0; k < data.intervals.size(); ++k) {
      const TangencyStripInterval& interval = data.intervals[k];
      if(interval.signs.empty()) { if(interval.phase > 0) LineIntegration::FullStripMoments(degree,strip,interval.left,interval.right,values,errors); ++accepted; continue; }
      std::vector<long double> leftExtents(values.size(),0);
      std::vector<long double> rightExtents(values.size(),0);
      if(interval.leftTangent) { ++Statistics().tangencies; Endpoint(data,interval,degree,true,computeDelta,0.25L*intervalBudget*absoluteH,0.25L*intervalBudget*absoluteD,relative,leftExtents,values,errors); }
      if(interval.rightTangent) { ++Statistics().tangencies; Endpoint(data,interval,degree,false,computeDelta,0.25L*intervalBudget*absoluteH,0.25L*intervalBudget*absoluteD,relative,rightExtents,values,errors); }
      std::vector<long double> events = {interval.left,interval.right};
      for(unsigned n = 0; n < values.size(); ++n) { if(leftExtents[n] > 0) events.push_back(interval.left+leftExtents[n]); if(rightExtents[n] > 0) events.push_back(interval.right-rightExtents[n]); }
      LineSliceUtilities::SortDistinct(events);
      const long double budget = 0.5L*intervalBudget/std::max(1U,static_cast<unsigned>(events.size())-1U);
      for(unsigned j = 1; j < events.size(); ++j) {
        if(!(events[j] > events[j-1U])) continue;
        const long double middle = events[j-1U]+0.5L*(events[j]-events[j-1U]);
        std::vector<unsigned char> active(values.size(),0);
        bool any = false;
        for(unsigned n = 0; n < values.size(); ++n) if((n < count || computeDelta) && middle >= interval.left+leftExtents[n] && middle <= interval.right-rightExtents[n]) { active[n] = 1; any = true; }
        if(any) Adaptive(data,interval,degree,events[j-1U],events[j],budget*absoluteH,budget*absoluteD,relative,0,maximumDepth,maximumEvaluations,evaluations,accepted,active,values,errors);
      }
    }
    for(unsigned n = 0; n < values.size(); ++n) if(!std::isfinite(values[n]) || !std::isfinite(errors[n]) || errors[n] > (n < count ? absoluteH : absoluteD)+relative*std::abs(values[n])) throw std::runtime_error("Accumulated non-graph moment error exceeds tolerance.");
  }
};

}

#endif
