#ifndef FEMUS_ENDPOINT_STRIP_QUADRATURE_HPP
#define FEMUS_ENDPOINT_STRIP_QUADRATURE_HPP

#include "LineIntegration.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace FELevelSetMoments {

struct EndpointStripStatistics {
  unsigned long long stripSamples = 0;
  unsigned long long stripBisections = 0;
  unsigned long long bulkBisections = 0;
  unsigned maximumStripDepth = 0;
  unsigned maximumBulkDepth = 0;
  long double minimumBulkWidth = 1;
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

struct EndpointStripInterval {
  long double left = 0;
  long double right = 1;
  bool leftTangent = false;
  bool rightTangent = false;
  int mapping = 0;
  long double endpoint = 0;
  long double delta = 0;
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
  unsigned triangle = 0;
  unsigned order = 1;
  long double score = std::numeric_limits<long double>::infinity();
};

class EndpointStripQuadrature {
public:
  // Experiments with selective GK on the original curve after quadratic endpoint coordinate maps.
  static EndpointStripStatistics& Statistics() { static thread_local EndpointStripStatistics statistics; return statistics; }

  static long double& StripWidth() { static thread_local long double width = 0.03L; return width; }
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

  static bool Try(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool triangle, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned& selected, const std::vector<std::vector<std::vector<long double>>>* gradient = nullptr, long double metricScale = 1, const std::vector<std::vector<std::vector<long double>>>* coordinates = nullptr, const std::vector<std::vector<double>>* inverseJacobian = nullptr) {
    if(degree > 10) return false;
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
    const unsigned first = PreferredDirection() == 0 ? (directions[1].score < directions[0].score ? 1U : 0U) : PreferredDirection()-1U;
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

  static void Classify(const EndpointStripDirection& data, long double x, EndpointStripInterval& interval) {
    const std::vector<long double> domain = {0,1,0,1};
    IntervalBranches(data.field,x,data.triangle,domain,false,false,true,interval.signs,interval.boundary,interval.positive,interval.phase);
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

  static EndpointStripInterval MappedStrip(const EndpointStripDirection& data, const EndpointStripInterval& interval, bool left, long double delta) {
    EndpointStripInterval strip = interval;
    strip.mapping = left ? 1 : -1;
    strip.endpoint = left ? interval.left : interval.right;
    strip.delta = delta;
    const std::vector<long double> discriminant = DiscriminantPolynomial(data.field);
    strip.reducedDiscriminant.assign(discriminant.size()-1U,0);
    for(unsigned i = 1; i < discriminant.size(); ++i) for(unsigned k = i; k < discriminant.size(); ++k) {
      long double binomial = 1;
      for(unsigned n = 0; n < i; ++n) binomial *= static_cast<long double>(k-n)/(i-n);
      strip.reducedDiscriminant[i-1U] += binomial*discriminant[k]*ReferencePower(strip.endpoint,k-i);
    }
    Statistics().minimumDelta = std::min(Statistics().minimumDelta,delta);
    Statistics().maximumDelta = std::max(Statistics().maximumDelta,delta);
    ++Statistics().tangencies;
    return strip;
  }

  static void Sample(const EndpointStripDirection& data, const EndpointStripInterval& interval, unsigned degree, long double x, const std::vector<unsigned char>& active, std::vector<long double>& values) {
    const unsigned count = active.size()/2U;
    const long double parameter = x;
    const long double offset = interval.mapping*interval.delta*parameter*parameter;
    if(interval.mapping != 0) x = interval.endpoint+offset;
    const long double jacobian = interval.mapping == 0 ? 1.0L : 2*interval.delta*parameter;
    const long double a = PolynomialValue(data.field[0],x);
    const long double b = PolynomialValue(data.field[1],x);
    const long double c = PolynomialValue(data.field[2],x);
    const long double discriminant = std::fma(-4*a,c,b*b);
    const long double reduced = interval.mapping == 0 ? 0 : interval.mapping*interval.delta*PolynomialValue(interval.reducedDiscriminant,offset);
    if(interval.mapping != 0 && !(reduced > 0)) throw std::runtime_error("An endpoint map requires a positive simple-fold discriminant.");
    const long double regularDerivative = interval.mapping == 0 ? 0 : std::sqrt(reduced);
    const long double derivative = interval.mapping == 0 ? (a == 0 ? std::abs(b) : std::sqrt(std::max(0.0L,discriminant))) : parameter*regularDerivative;
    const long double surfaceRatio = interval.mapping == 0 ? 1.0L/derivative : 2*interval.delta/regularDerivative;
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
      long double weight = (interval.boundary[r] == 0 ? 1.0L : 0.5L)*surfaceRatio;
      if(!data.gradient.empty()) weight *= SurfaceGradientNorm(data.gradient,x,roots[r]);
      long double power = 1;
      for(unsigned j = 0; j <= degree; ++j) { d[j] += weight*power; power *= roots[r]; }
    }
    values.assign(active.size(),0);
    long double power = 1;
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned n = LineSliceUtilities::MomentIndex(degree,i,j);
        if(active[n]) values[n] = jacobian*power*h[j]; else ++Statistics().inactiveComponents;
        if(active[count+n]) values[count+n] = power*d[j]; else ++Statistics().inactiveComponents;
      }
      power *= x;
    }
    if(interval.mapping == 0) ++Statistics().bulkSamples; else ++Statistics().stripSamples;
  }

  template<unsigned Rule> static void Panel(const EndpointStripDirection& data, const EndpointStripInterval& interval, unsigned degree, long double left, long double right, const std::vector<unsigned char>& active, unsigned maximumEvaluations, unsigned& evaluations, std::vector<long double>& values, std::vector<long double>& errors) {
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

  static void Adaptive(const EndpointStripDirection& data, const EndpointStripInterval& interval, unsigned degree, long double left, long double right, long double absoluteH, long double absoluteD, long double relative, unsigned depth, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, const std::vector<unsigned char>& requested, std::vector<long double>& total, std::vector<long double>& totalError) {
    std::vector<unsigned char> active = requested;
    std::vector<long double> values;
    std::vector<long double> errors;
    Panel<15>(data,interval,degree,left,right,active,maximumEvaluations,evaluations,values,errors);
    if(!LineIntegration::Accept(values,errors,absoluteH,absoluteD,relative,active,total,totalError)) { ++accepted; if(interval.mapping == 0) Statistics().minimumBulkWidth = std::min(Statistics().minimumBulkWidth,right-left); return; }
    Panel<31>(data,interval,degree,left,right,active,maximumEvaluations,evaluations,values,errors);
    if(!LineIntegration::Accept(values,errors,absoluteH,absoluteD,relative,active,total,totalError)) { ++accepted; if(interval.mapping == 0) Statistics().minimumBulkWidth = std::min(Statistics().minimumBulkWidth,right-left); return; }
    if(interval.mapping == 0) { ++Statistics().bulkBisections; Statistics().maximumBulkDepth = std::max(Statistics().maximumBulkDepth,depth+1U); }
    else { ++Statistics().stripBisections; Statistics().maximumStripDepth = std::max(Statistics().maximumStripDepth,depth+1U); }
    const long double middle = left+0.5L*(right-left);
    if(depth >= maximumDepth || middle == left || middle == right) throw std::runtime_error("Non-graph GK did not reach tolerance.");
    Adaptive(data,interval,degree,left,middle,0.5L*absoluteH,0.5L*absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,active,total,totalError);
    Adaptive(data,interval,degree,middle,right,0.5L*absoluteH,0.5L*absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,active,total,totalError);
  }

  static void IntegrateDirection(const EndpointStripDirection& data, unsigned degree, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    values.assign(2U*count,0);
    errors.assign(values.size(),0);
    const long double intervalBudget = 1.0L/std::max(1U,static_cast<unsigned>(data.intervals.size()));
    const LineIntegrationStrip strip = {0,1,0,0,1,data.triangle ? -1.0L : 0.0L};
    for(unsigned k = 0; k < data.intervals.size(); ++k) {
      const EndpointStripInterval& interval = data.intervals[k];
      if(interval.signs.empty()) { if(interval.phase > 0) LineIntegration::FullStripMoments(degree,strip,interval.left,interval.right,values,errors); ++accepted; continue; }
      std::vector<unsigned char> active(values.size(),1);
      if(!computeDelta) std::fill(active.begin()+count,active.end(),0);
      const long double width = std::min(StripWidth(),0.45L*(interval.right-interval.left));
      long double left = interval.left;
      long double right = interval.right;
      if(interval.leftTangent) {
        const EndpointStripInterval mapped = MappedStrip(data,interval,true,width);
        Adaptive(data,mapped,degree,0,1,0.25L*intervalBudget*absoluteH,0.25L*intervalBudget*absoluteD,relative,0,maximumDepth,maximumEvaluations,evaluations,accepted,active,values,errors);
        left += width;
      }
      if(interval.rightTangent) {
        const EndpointStripInterval mapped = MappedStrip(data,interval,false,width);
        Adaptive(data,mapped,degree,0,1,0.25L*intervalBudget*absoluteH,0.25L*intervalBudget*absoluteD,relative,0,maximumDepth,maximumEvaluations,evaluations,accepted,active,values,errors);
        right -= width;
      }
      if(right > left) Adaptive(data,interval,degree,left,right,0.5L*intervalBudget*absoluteH,0.5L*intervalBudget*absoluteD,relative,0,maximumDepth,maximumEvaluations,evaluations,accepted,active,values,errors);
    }
    for(unsigned n = 0; n < values.size(); ++n) if(!std::isfinite(values[n]) || !std::isfinite(errors[n]) || errors[n] > (n < count ? absoluteH : absoluteD)+relative*std::abs(values[n])) throw std::runtime_error("Accumulated non-graph moment error exceeds tolerance.");
  }
};

}

#endif
