#ifndef FEMUS_LINE_INTEGRATION_HPP
#define FEMUS_LINE_INTEGRATION_HPP

#include "LineSliceUtilities.hpp"
#include "FELevelSetMoments.hpp"
#include "OuterIntervalQuadrature.hpp"
#include <array>
#include <map>
#include "CutFemPreparedContext.hpp"

namespace FELevelSetMoments {

struct LineIntegrationStrip {
  long double left;
  long double right;
  long double lowerConstant;
  long double lowerSlope;
  long double upperConstant;
  long double upperSlope;
  bool fullBoundaryDelta = false;
  bool certifiedLowerZero = false;
  bool certifiedUpperZero = false;
};

struct LineIntegrationStatistics {
  unsigned long long fastSamples = 0;
  unsigned long long rootSolves = 0;
  unsigned long long singleRootSamples = 0;
  unsigned long long inactiveSkipped = 0;
};

struct LineIntegrationDiscriminant {
  std::vector<long double> coefficient;
  long double constant = 0;
};

struct GraphFiberRoots {
  long double positiveLeft, positiveRight, root, signedDerivative;
  bool boundary;
};

struct PhysicalGraphFiber {
  GraphFiberRoots geometry;
  PhysicalSurfaceMetric metric;
};

struct GraphStripSamples {
  std::map<long double,PhysicalGraphFiber> fibers;
  SharedSurfacePanel panel;
};

struct GraphStripGeometry {
  std::vector<long double> events, tangencies;
  std::map<std::array<long double,4>,GraphStripSamples> panels;
};

struct SharedGraphStrips {
  std::array<GraphStripGeometry,2> directions;
};

struct LineGraphContext {
  GraphStripSamples* samples = nullptr;
  const PreparedPlane* prepared = nullptr;
  unsigned preparedOrder = 1;
  const std::vector<std::vector<long double>>* coefficient = nullptr;
  unsigned degree = 0;
  LineIntegrationStrip strip = {0,1,0,0,1,0};
  bool computeDelta = true;
  bool inverseSurface = false;
  int normalComponent = -1;
  bool lowerZero = false;
  bool upperZero = false;
  const std::vector<std::vector<std::vector<long double>>>* gradient = nullptr;
  const std::array<long double,5>* coordinates = nullptr;
  int derivativeSign = 0;
  OuterInterval interval;
  LineIntegrationDiscriminant discriminant;
  unsigned* evaluations = nullptr;
};

class LineIntegration {
public:
  // Shared affine-fiber graph quadrature retains accepted moments during selective refinement.
  static LineIntegrationStatistics& Statistics() {
    static thread_local LineIntegrationStatistics data;
    return data;
  }

  static std::array<long double,5> MomentCoordinates(const std::array<long double,5>* coordinates, unsigned direction) {
    std::array<long double,5> result = coordinates == nullptr ? std::array<long double,5>{0,1,0,1,1} : *coordinates;
    if(direction == 2) { std::swap(result[0],result[2]); std::swap(result[1],result[3]); }
    return result;
  }

  static std::vector<unsigned char> MomentRequest(unsigned degree, const std::vector<unsigned char>* requested, unsigned direction, bool computeDelta) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    std::vector<unsigned char> result = requested == nullptr ? std::vector<unsigned char>(2U*count,1) : *requested;
    if(result.size() != 2U*count) throw std::invalid_argument("Invalid requested moment mask.");
    if(!computeDelta) std::fill(result.begin()+count,result.end(),0);
    if(direction == 2) for(unsigned i = 0; i <= degree; ++i) for(unsigned j = i+1U; j <= degree-i; ++j) for(unsigned pair = 0; pair < 2; ++pair) std::swap(result[pair*count+LineSliceUtilities::MomentIndex(degree,i,j)],result[pair*count+LineSliceUtilities::MomentIndex(degree,j,i)]);
    return result;
  }

  static bool BoundaryIsZero(const std::vector<std::vector<long double>>& coefficient, long double constant, long double slope) {
    const std::vector<long double> boundary = BoundaryPolynomial(coefficient,constant,slope);
    return LineSliceUtilities::PolynomialNorm(boundary) == 0.0L;
  }

  static long double QuadraticRootDerivative(const std::array<long double,3>& polynomial, long double mappedDerivative) {
    if(mappedDerivative >= 0) return mappedDerivative;
    const long double discriminant = std::fma(-4.0L*polynomial[2],polynomial[0],polynomial[1]*polynomial[1]);
    return discriminant >= 0 ? std::sqrt(discriminant) : -1.0L;
  }

  static GraphFiberRoots FindGraphFiber(const std::vector<std::vector<long double>>& coefficient, const LineIntegrationStrip& strip, long double x, bool lowerZero, bool upperZero, long double mappedDerivative, int derivativeSign, const PreparedPlane* prepared, unsigned preparedOrder) {
    const long double lower = StripLower(strip,x);
    const long double upper = StripUpper(strip,x);
    std::array<long double,3> polynomial = {0,0,0};
    if(prepared == nullptr || mappedDerivative >= 0) for(unsigned j = 0; j < 3; ++j) polynomial[j] = coefficient[2][j] == 0.0L ? coefficient[0][j]+coefficient[1][j]*x : std::fma(std::fma(coefficient[2][j],x,coefficient[1][j]),x,coefficient[0][j]);
    ++Statistics().rootSolves;
    const long double tolerance = 256.0L*std::numeric_limits<long double>::epsilon()*std::max(1.0L,std::max(std::abs(lower),std::abs(upper)));
    long double root = 0.0L;
    long double derivative = 0.0L;
    long double signedDerivative = 0.0L;
    unsigned number = 0;
    bool preparedPositive = false;
    if(prepared != nullptr && mappedDerivative < 0) {
      if(derivativeSign != 0) ++Statistics().singleRootSamples;
      const CutFemLevelSetData::LineAnalysis& analysis = prepared->Analyze(preparedOrder,x);
      if(analysis.identicallyZero) throw std::runtime_error("A graph projection reached an entire zero fiber.");
      preparedPositive = analysis.positiveIntervalCount != 0;
      for(unsigned n = 0; n < analysis.rootCount; ++n) { const long double candidate = prepared->RootCoordinate(preparedOrder,analysis.roots[n]); if(candidate >= lower-tolerance && candidate <= upper+tolerance) { root = std::max(lower,std::min(upper,candidate)); signedDerivative = prepared->WorkingDerivative(preparedOrder,analysis.roots[n]); derivative = std::abs(signedDerivative); ++number; } }
      if(number > 1) throw std::runtime_error("A graph projection contains multiple roots.");
      if(number == 1 && (lowerZero || upperZero)) root = lowerZero ? lower : upper;
    }
    else if(derivativeSign != 0) {
      ++Statistics().singleRootSamples;
      bool available = false;
      if(lowerZero || upperZero) {
        root = lowerZero ? lower : upper;
        derivative = std::abs(polynomial[2] == 0.0L ? polynomial[1] : std::fma(2.0L*polynomial[2],root,polynomial[1]));
        available = true;
      }
      else if(polynomial[2] == 0.0L) {
        if(polynomial[1] != 0.0L) { root = -polynomial[0]/polynomial[1]; derivative = std::abs(polynomial[1]); available = true; }
      }
      else {
        derivative = QuadraticRootDerivative(polynomial,mappedDerivative);
        if(derivative >= 0.0L) {
          const long double denominator = polynomial[1]+derivativeSign*derivative;
          root = derivativeSign*polynomial[1] >= 0.0L && denominator != 0.0L ? -2.0L*polynomial[0]/denominator : 0.5L*(-polynomial[1]+derivativeSign*derivative)/polynomial[2];
          available = true;
        }
      }
      if(available && root >= lower-tolerance && root <= upper+tolerance) { root = std::max(lower,std::min(upper,root)); number = 1; }
      signedDerivative = derivativeSign;
    }
    else {
      const long double quadraticDerivative = polynomial[2] == 0.0L ? -1.0L : QuadraticRootDerivative(polynomial,mappedDerivative);
      std::array<long double,2> roots;
      unsigned rootsCount = 0;
      if(lowerZero || upperZero) {
        roots[rootsCount++] = lowerZero ? lower : upper;
        if(polynomial[2] != 0.0L) {
          const long double other = -polynomial[1]/polynomial[2]-roots[0];
          if(other != roots[0]) roots[rootsCount++] = other;
        }
      }
      else if(polynomial[2] == 0.0L) {
        if(polynomial[1] != 0.0L) roots[rootsCount++] = -polynomial[0]/polynomial[1];
      }
      else {
        if(quadraticDerivative == 0.0L) roots[rootsCount++] = -polynomial[1]/(2.0L*polynomial[2]);
        else if(quadraticDerivative > 0.0L) {
          const long double numerator = -0.5L*(polynomial[1]+std::copysign(quadraticDerivative,polynomial[1]));
          roots[rootsCount++] = numerator/polynomial[2];
          roots[rootsCount++] = polynomial[0]/numerator;
        }
      }
      for(unsigned r = 0; r < rootsCount; ++r) if(roots[r] >= lower-tolerance && roots[r] <= upper+tolerance) { root = std::max(lower,std::min(upper,roots[r])); ++number; }
      if(number > 1) throw std::runtime_error("The cached graph branch is not a single root.");
      if(number == 1) {
        signedDerivative = polynomial[2] == 0.0L ? polynomial[1] : std::fma(2.0L*polynomial[2],root,polynomial[1]);
        derivative = polynomial[2] == 0.0L ? std::abs(polynomial[1]) : (quadraticDerivative >= 0.0L ? quadraticDerivative : std::abs(signedDerivative));
      }
    }
    long double positiveLeft = lower;
    long double positiveRight = lower;
    if(number == 1) {
      if(mappedDerivative >= 0.0L) derivative = mappedDerivative;
      if(!(derivative > 0.0L)) throw std::runtime_error("A singular fiber has no regular graph delta weight.");
      positiveLeft = signedDerivative > 0.0L ? root : lower;
      positiveRight = signedDerivative > 0.0L ? upper : root;
    }
    else if(prepared != nullptr && mappedDerivative < 0 ? preparedPositive : std::fma(std::fma(polynomial[2],0.5L*(lower+upper),polynomial[1]),0.5L*(lower+upper),polynomial[0]) > 0.0L) positiveRight = upper;
    return {positiveLeft,positiveRight,root,number == 1 ? std::copysign(derivative,signedDerivative) : 0,root == lower || root == upper};
  }

  static void AddPhysicalSurfaceMoments(unsigned degree, long double x, long double y, long double factor, const PhysicalSurfaceMetric& metric, const std::vector<unsigned char>& active, std::vector<long double>& values) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    const long double inverse = 1/metric.norm;
    const bool regularInverse = std::isfinite(inverse);
    std::array<long double*,4> output;
    std::array<long double,4> weights;
    unsigned size = 0;
    for(unsigned block = 1; block < 5; ++block) if(active[block*count]) {
      output[size] = values.data()+block*count;
      weights[size++] = block == 1 ? (regularInverse ? factor*inverse : factor/metric.norm) : factor*(regularInverse ? metric.gradient[block-2]*inverse : metric.gradient[block-2]/metric.norm);
    }
    long double xPower = 1;
    unsigned n = 0;
    for(unsigned i = 0; i <= degree; ++i) {
      long double monomial = xPower;
      for(unsigned j = 0; j <= degree-i; ++j, ++n) {
        for(unsigned b = 0; b < size; ++b) output[b][n] += weights[b]*monomial;
        monomial *= y;
      }
      xPower *= x;
    }
  }

  static void FastGraphSample(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, long double x, bool computeDelta, bool lowerZero, bool upperZero, std::vector<long double>& values, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>& active, long double mappedDerivative = -1.0L, const std::array<long double,5>* coordinates = nullptr, int derivativeSign = 0, bool inverseSurface = false, int normalComponent = -1, const PreparedPlane* prepared = nullptr, unsigned preparedOrder = 1, const PhysicalGraphFiber* shared = nullptr) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    const bool batch = active.size() == 5U*count;
    values.assign(active.size(),0.0L);
    const long double parentX = coordinates == nullptr ? x : (*coordinates)[0]+(*coordinates)[1]*x;
    const long double measure = coordinates == nullptr ? 1.0L : (*coordinates)[4];
    GraphFiberRoots local;
    if(shared == nullptr) local = FindGraphFiber(coefficient,strip,x,lowerZero,upperZero,mappedDerivative,derivativeSign,prepared,preparedOrder);
    const GraphFiberRoots& fiber = shared == nullptr ? local : shared->geometry;
    if(std::find(active.begin(),active.begin()+count,1) != active.begin()+count) {
      const long double parentLeft = coordinates == nullptr ? fiber.positiveLeft : (*coordinates)[2]+(*coordinates)[3]*fiber.positiveLeft;
      const long double parentRight = coordinates == nullptr ? fiber.positiveRight : (*coordinates)[2]+(*coordinates)[3]*fiber.positiveRight;
      long double rightPower = 1.0L;
      long double geometric = 0.0L;
      for(unsigned j = 0; j <= degree; ++j) {
        geometric = rightPower+parentLeft*geometric;
        values[j] = measure*(fiber.positiveRight-fiber.positiveLeft)*geometric/(j+1U);
        rightPower *= parentRight;
      }
    }
    if(batch && fiber.signedDerivative != 0) {
      const long double factor = measure*(fiber.boundary && !strip.fullBoundaryDelta ? 0.5L : 1.0L)/std::abs(fiber.signedDerivative);
      const long double parentRoot = (*coordinates)[2]+(*coordinates)[3]*fiber.root;
      AddPhysicalSurfaceMoments(degree,parentX,parentRoot,factor,prepared->EvaluateMetric(preparedOrder,x,fiber.root),active,values);
    }
    else if(!batch && computeDelta && fiber.signedDerivative != 0 && std::find(active.begin()+count,active.end(),1) != active.end()) {
      long double weight = measure*(fiber.boundary && !strip.fullBoundaryDelta ? 0.5L : 1.0L);
      if(shared != nullptr) weight *= shared->metric.Factor(inverseSurface,normalComponent)/std::abs(fiber.signedDerivative);
      else if(prepared != nullptr) weight *= prepared->Density(preparedOrder,x,fiber.root,fiber.signedDerivative);
      else { weight /= std::abs(fiber.signedDerivative); if(surfaceGradient != nullptr) weight *= SurfaceGradientFactor(*surfaceGradient,x,fiber.root,inverseSurface,normalComponent); }
      const long double parentRoot = coordinates == nullptr ? fiber.root : (*coordinates)[2]+(*coordinates)[3]*fiber.root;
      long double rootPower = 1.0L;
      for(unsigned j = 0; j <= degree; ++j) { values[count+j] = rootPower*weight; rootPower *= parentRoot; }
    }
    if(batch && !active[0]) return;
    long double xPower = 1.0L;
    for(unsigned i = 1; i <= degree; ++i) {
      xPower *= parentX;
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
        if(active[k]) values[k] = xPower*values[j];
        else ++Statistics().inactiveSkipped;
        if(!batch) { if(active[count+k]) values[count+k] = xPower*values[count+j]; else ++Statistics().inactiveSkipped; }
      }
    }
  }

  static LineIntegrationDiscriminant EndpointDiscriminant(const std::vector<std::vector<long double>>& coefficient, const OuterInterval& interval, const PreparedPlane* prepared = nullptr, unsigned preparedOrder = 1) {
    LineIntegrationDiscriminant result;
    if(interval.mapping == 0) return result;
    if(prepared != nullptr) { std::array<long double,5> shifted; prepared->ShiftedDiscriminant(preparedOrder,interval.endpoint,shifted); result.constant = interval.exactEndpoint ? 0 : shifted[0]; result.coefficient.assign(shifted.begin()+1,shifted.end()); return result; }
    result.constant = interval.exactEndpoint ? 0 : FiberDiscriminantValue(coefficient,interval.endpoint);
    std::vector<std::vector<long double>> field(3,std::vector<long double>(3));
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) field[2U-j][i] = coefficient[i][j];
    const std::vector<long double> polynomial = LineSliceUtilities::DiscriminantPolynomial(field);
    result.coefficient.assign(polynomial.size()-1U,0);
    for(unsigned i = 1; i < polynomial.size(); ++i) for(unsigned k = i; k < polynomial.size(); ++k) result.coefficient[i-1U] += ReferenceBinomial(k,i)*polynomial[k]*ReferencePower(interval.endpoint,k-i);
    return result;
  }

  static bool NearEndpoint(const std::vector<std::vector<long double>>& coefficient, long double point, int inward, long double width, const PreparedPlane* prepared = nullptr, unsigned preparedOrder = 1) {
    if(FiberPolynomial(coefficient,point)[2] == 0) return false;
    const long double value = prepared == nullptr ? FiberDiscriminantValue(coefficient,point) : prepared->Discriminant(preparedOrder,point);
    const long double derivative = inward*(prepared == nullptr ? FiberDiscriminantDerivative(coefficient,point) : prepared->DiscriminantDerivative(preparedOrder,point));
    return value >= 0 && derivative > 0 && value <= 0.25L*width*derivative;
  }

  static long double EndpointDerivative(const OuterInterval& interval, const LineIntegrationDiscriminant& discriminant, const OuterIntervalPoint& point) {
    const long double reduced = discriminant.constant/(point.parameter*point.parameter)+interval.mapping*interval.delta*LineSliceUtilities::PolynomialValue(discriminant.coefficient,point.offset);
    if(!(reduced > 0)) throw std::runtime_error("A transformed simple fold requires a positive discriminant.");
    return point.parameter*std::sqrt(reduced);
  }

  static const PhysicalGraphFiber& SharedGraphFiber(LineGraphContext& context, const OuterIntervalPoint& point) {
    std::map<long double,PhysicalGraphFiber>& fibers = context.samples->fibers;
    const std::map<long double,PhysicalGraphFiber>::iterator entry = fibers.lower_bound(point.parameter);
    if(entry != fibers.end() && entry->first == point.parameter) return entry->second;
    const long double derivative = context.interval.mapping == 0 ? -1 : EndpointDerivative(context.interval,context.discriminant,point);
    PhysicalGraphFiber fiber;
    fiber.geometry = FindGraphFiber(*context.coefficient,context.strip,point.coordinate,context.lowerZero,context.upperZero,derivative,context.derivativeSign,context.prepared,context.preparedOrder);
    if(fiber.geometry.signedDerivative != 0) fiber.metric = context.prepared->EvaluateMetric(context.preparedOrder,point.coordinate,fiber.geometry.root);
    return fibers.emplace_hint(entry,point.parameter,std::move(fiber))->second;
  }

  static void GraphSample(LineGraphContext& context, const OuterIntervalPoint& point, const std::vector<unsigned char>& active, std::vector<long double>& values, std::vector<long double>& errors, long double absoluteH, long double absoluteD, long double relative) {
    ++*context.evaluations;
    ++Statistics().fastSamples;
    const unsigned count = (context.degree+1U)*(context.degree+2U)/2U;
    const PhysicalGraphFiber* shared = context.samples == nullptr || point.precisionProbe || active.size() != 2U*count ? nullptr : &SharedGraphFiber(context,point);
    const long double derivative = shared != nullptr || context.interval.mapping == 0 ? -1 : EndpointDerivative(context.interval,context.discriminant,point);
    FastGraphSample(*context.coefficient,context.degree,context.strip,point.coordinate,context.computeDelta,context.lowerZero,context.upperZero,values,context.gradient,active,derivative,context.coordinates,context.derivativeSign,context.inverseSurface,context.normalComponent,context.prepared,context.preparedOrder,shared);
    errors.clear();
    if(point.jacobian != 1) for(unsigned n = 0; n < values.size(); ++n) if(active[n]) values[n] *= point.jacobian;
  }

  static bool Accept(const std::vector<long double>& values, const std::vector<long double>& errors, long double absoluteH, long double absoluteD, long double relative, std::vector<unsigned char>& active, std::vector<long double>& total, std::vector<long double>& totalErrors) {
    return OuterIntervalQuadrature::Accept(values,errors,absoluteH,absoluteD,relative,active,total,totalErrors,2);
  }

  static void GaussRule(unsigned number, std::vector<long double>& points, std::vector<long double>& weights) {
    points.resize(number);
    weights.resize(number);
    const long double pi = std::acos(-1.0L);
    for(unsigned k = 0; k < (number+1U)/2U; ++k) {
      long double root = std::cos(pi*(k+0.75L)/(number+0.5L));
      long double derivative = 0.0L;
      for(unsigned iteration = 0; iteration < 100; ++iteration) {
        long double previous = 1.0L;
        long double value = root;
        for(unsigned n = 2; n <= number; ++n) {
          const long double next = ((2.0L*n-1.0L)*root*value-(n-1.0L)*previous)/n;
          previous = value;
          value = next;
        }
        derivative = number*(root*value-previous)/(root*root-1.0L);
        const long double next = root-value/derivative;
        if(next == root || std::abs(next-root) <= 4.0L*std::numeric_limits<long double>::epsilon()) {
          root = next;
          break;
        }
        root = next;
      }
      points[k] = 0.5L*(1.0L-root);
      points[number-1U-k] = 0.5L*(1.0L+root);
      weights[k] = 1.0L/((1.0L-root*root)*derivative*derivative);
      weights[number-1U-k] = weights[k];
    }
  }

  struct FullStripRule {
    std::vector<long double> points;
    std::vector<long double> weights;
  };

  static std::array<FullStripRule,34> BuildFullStripRules() { std::array<FullStripRule,34> rules; for(unsigned count = 1; count < rules.size(); ++count) GaussRule(count,rules[count].points,rules[count].weights); return rules; }

  static const std::array<FullStripRule,34>& FullStripRules() { static const std::array<FullStripRule,34> rules = BuildFullStripRules(); return rules; }

  static void FullStripMoments(unsigned degree, const LineIntegrationStrip& strip, long double left, long double right, std::vector<long double>& values, std::vector<long double>& errors, const std::array<long double,5>* coordinates = nullptr, const std::vector<unsigned char>* requested = nullptr) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    if(requested != nullptr && std::find(requested->begin(),requested->begin()+count,1) == requested->begin()+count) return;
    const FullStripRule& rule = FullStripRules()[(degree+3U)/2U];
    const std::vector<long double>& points = rule.points;
    const std::vector<long double>& weights = rule.weights;
    std::vector<long double> fiber(degree+1U);
    for(unsigned q = 0; q < points.size(); ++q) {
      const long double x = left+(right-left)*points[q];
      const long double lower = StripLower(strip,x);
      const long double upper = StripUpper(strip,x);
      const long double parentX = coordinates == nullptr ? x : LineSliceUtilities::ParentCoordinate((*coordinates)[1],x,(*coordinates)[0]);
      const long double parentLower = coordinates == nullptr ? lower : LineSliceUtilities::ParentCoordinate((*coordinates)[3],lower,(*coordinates)[2]);
      const long double parentUpper = coordinates == nullptr ? upper : LineSliceUtilities::ParentCoordinate((*coordinates)[3],upper,(*coordinates)[2]);
      const long double weight = (right-left)*weights[q]*(coordinates == nullptr ? 1.0L : (*coordinates)[4]);
      long double xPower = 1.0L;
      long double geometric = 0.0L;
      long double power = 1.0L;
      for(unsigned j = 0; j <= degree; ++j) {
        geometric = power+parentLower*geometric;
        fiber[j] = (upper-lower)*geometric/(j+1U);
        power *= parentUpper;
      }
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) {
          const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
          if(requested != nullptr && (*requested)[k] == 0) continue;
          const long double term = weight*xPower*fiber[j];
          values[k] += term;
          errors[k] += (128.0L+32.0L*degree)*std::numeric_limits<long double>::epsilon()*std::abs(term);
        }
        xPower *= parentX;
      }
    }
  }

  static void IntegrateStrip(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, bool graph, bool mapped, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr, const std::array<long double,5>* coordinates = nullptr, const std::vector<unsigned char>* requested = nullptr, const OuterQuadratureRules& rules = OuterQuadratureRules(), bool inverseSurface = false, int normalComponent = -1, const PreparedPlane* prepared = nullptr, unsigned preparedOrder = 1) {
    rules.Validate();
    if(!graph) throw std::invalid_argument("The LINE outer graph path requires a certified graph.");
    GraphStripGeometry* shared = prepared != nullptr && prepared->graphs != nullptr ? &prepared->graphs->directions[preparedOrder-1] : nullptr;
    std::vector<long double> localEvents, localTangencies;
    if(shared == nullptr || shared->events.empty()) {
      CollectEvents(coefficient,strip,localEvents,&localTangencies,!mapped,prepared,preparedOrder);
      if(shared != nullptr) { shared->events.swap(localEvents); shared->tangencies.swap(localTangencies); }
    }
    const std::vector<long double>& events = shared == nullptr ? localEvents : shared->events;
    const std::vector<long double>& tangencies = shared == nullptr ? localTangencies : shared->tangencies;
    const long double intervalBudget = events.size() > 1 ? 1.0L/(events.size()-1U) : 1;
    const bool lowerZero = strip.certifiedLowerZero || (prepared == nullptr ? BoundaryIsZero(coefficient,strip.lowerConstant,strip.lowerSlope) : prepared->ZeroBoundary(preparedOrder,strip.lowerConstant,strip.lowerSlope));
    const bool upperZero = strip.certifiedUpperZero || (prepared == nullptr ? BoundaryIsZero(coefficient,strip.upperConstant,strip.upperSlope) : prepared->ZeroBoundary(preparedOrder,strip.upperConstant,strip.upperSlope));
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    std::vector<unsigned char> active = requested == nullptr ? std::vector<unsigned char>(2U*count,1) : *requested;
    if(!computeDelta) std::fill(active.begin()+count,active.end(),0);
    for(unsigned interval = 1; interval < events.size(); ++interval) {
      const long double left = events[interval-1U];
      const long double right = events[interval];
      if(!(right > left)) continue;
      const bool leftTangency = std::find(tangencies.begin(),tangencies.end(),left) != tangencies.end();
      const bool rightTangency = std::find(tangencies.begin(),tangencies.end(),right) != tangencies.end();
      const long double middle = left+0.5L*(right-left);
      const long double lower = StripLower(strip,middle);
      const long double upper = StripUpper(strip,middle);
      const std::vector<long double> polynomial = FiberPolynomial(coefficient,middle);
      unsigned number = 0;
      long double middleRoot = 0, preparedDerivative = 0;
      if(prepared != nullptr) {
        const CutFemLevelSetData::LineAnalysis& analysis = prepared->Analyze(preparedOrder,middle);
        for(unsigned n = 0; n < analysis.rootCount; ++n) { const long double root = prepared->RootCoordinate(preparedOrder,analysis.roots[n]); if(root >= lower && root <= upper) { middleRoot = root; preparedDerivative = prepared->WorkingDerivative(preparedOrder,analysis.roots[n]); ++number; } }
      }
      else { const std::vector<long double> roots = LineSliceUtilities::QuadraticRoots(polynomial[2],polynomial[1],polynomial[0]); for(unsigned n = 0; n < roots.size(); ++n) if(roots[n] >= lower && roots[n] <= upper) { middleRoot = roots[n]; ++number; } }
      if(lowerZero || upperZero) { number = 1; middleRoot = lowerZero ? lower : upper; }
      if(number == 0) {
        if(Evaluate(coefficient,middle,0.5L*(lower+upper)) > 0) FullStripMoments(degree,strip,left,right,values,errors,coordinates,&active);
        ++accepted;
        continue;
      }
      if(number != 1) throw std::runtime_error("The selected graph interval contains multiple roots.");
      const long double middleDerivative = prepared == nullptr ? std::fma(2.0L*polynomial[2],middleRoot,polynomial[1]) : preparedDerivative;
      const unsigned leftExponent = leftTangency || NearEndpoint(coefficient,left,1,right-left,prepared,preparedOrder) ? 2U : 1U;
      const unsigned rightExponent = rightTangency || NearEndpoint(coefficient,right,-1,right-left,prepared,preparedOrder) ? 2U : 1U;
      const std::vector<OuterInterval> plan = OuterIntervalQuadrature::Plan(left,right,leftExponent,rightExponent,2,leftTangency,rightTangency);
      for(unsigned part = 0; part < plan.size(); ++part) {
        LineGraphContext context;
        context.prepared = prepared;
        context.preparedOrder = preparedOrder;
        context.coefficient = &coefficient;
        context.degree = degree;
        context.strip = strip;
        context.computeDelta = computeDelta;
        context.lowerZero = lowerZero;
        context.upperZero = upperZero;
        context.gradient = surfaceGradient;
        context.inverseSurface = inverseSurface;
        context.normalComponent = normalComponent;
        context.coordinates = coordinates;
        context.derivativeSign = middleDerivative > 0 ? 1 : -1;
        context.interval = plan[part];
        if(shared != nullptr) {
          const OuterInterval& panel = plan[part];
          const std::array<long double,4> key = {panel.left,panel.right,static_cast<long double>(panel.mapping)*panel.exponent,static_cast<long double>(panel.exactEndpoint)};
          context.samples = &shared->panels[key];
        }
        context.discriminant = EndpointDiscriminant(coefficient,plan[part],prepared,preparedOrder);
        context.evaluations = &evaluations;
        const unsigned blocks = prepared == nullptr ? 0 : prepared->surfaceBlocks;
        OuterIntervalQuadrature::Adaptive(context,GraphSample,plan[part],plan[part].mapping == 0 ? plan[part].left : 0,plan[part].mapping == 0 ? plan[part].right : 1,intervalBudget*absoluteH/plan.size(),intervalBudget*absoluteD/plan.size(),relative,0,maximumDepth,maximumEvaluations,evaluations,accepted,2,active,values,errors,rules,blocks ? &context.samples->panel : nullptr,degree,normalComponent+2,blocks);
      }
    }
  }

  static bool HasVerticalFiber(const std::vector<std::vector<long double>>& coefficient, unsigned order) {
    std::vector<std::vector<long double>> fibers(3,std::vector<long double>(3,0.0L));
    for(unsigned j = 0; j < 3; ++j) {
      for(unsigned i = 0; i < 3; ++i) fibers[j][i] = order == 1 ? coefficient[i][j] : coefficient[j][i];
      LineSliceUtilities::Trim(fibers[j]);
      if(fibers[j].size() == 1 && fibers[j][0] != 0.0L) return false;
    }
    for(unsigned j = 0; j < 3; ++j) {
      if(fibers[j].size() < 2) continue;
      std::vector<long double> roots;
      LineSliceUtilities::PolynomialRoots(fibers[j],0.0L,1.0L,roots);
      for(unsigned r = 0; r < roots.size(); ++r) {
        bool common = true;
        for(unsigned k = 0; k < 3; ++k) {
          const long double residual = std::abs(LineSliceUtilities::PolynomialValue(fibers[k],roots[r]));
          const long double rounding = 128.0L*std::numeric_limits<long double>::epsilon()*LineSliceUtilities::PolynomialBound(fibers[k],roots[r]);
          common = common && residual <= rounding;
        }
        if(common) return true;
      }
      return false;
    }
    return false;
  }


  static long double Power(long double value, unsigned degree) {
    long double result = 1.0L;
    for(unsigned k = 0; k < degree; ++k) result *= value;
    return result;
  }

  static std::vector<long double> BoundaryPolynomial(const std::vector<std::vector<long double>>& coefficient, long double constant, long double slope) {
    std::vector<long double> result(5,0.0L);
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) {
        for(unsigned k = 0; k <= j; ++k) result[i+k] += coefficient[i][j]*LineSliceUtilities::Binomial(j,k)*Power(constant,j-k)*Power(slope,k);
      }
    }
    LineSliceUtilities::Trim(result);
    return result;
  }

  static long double Evaluate(const std::vector<std::vector<long double>>& coefficient, long double x, long double y) {
    long double value = 0.0L;
    for(unsigned i = 3; i > 0; --i) {
      value *= x;
      long double row = 0.0L;
      for(unsigned j = 3; j > 0; --j) row = std::fma(row,y,coefficient[i-1U][j-1U]);
      value += row;
    }
    return value;
  }

  static long double StripLower(const LineIntegrationStrip& strip, long double x) {
    return std::fma(strip.lowerSlope,x,strip.lowerConstant);
  }

  static long double StripUpper(const LineIntegrationStrip& strip, long double x) {
    return std::fma(strip.upperSlope,x,strip.upperConstant);
  }

  static std::vector<long double> FiberPolynomial(const std::vector<std::vector<long double>>& coefficient, long double x) {
    std::vector<long double> polynomial(3,0.0L);
    long double power = 1.0L;
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) polynomial[j] += coefficient[i][j]*power;
      power *= x;
    }
    return polynomial;
  }

  static std::array<long double,2> CompensatedSum(const std::array<long double,2>& first, const std::array<long double,2>& second) {
    const long double sum = first[0]+second[0];
    const long double part = sum-first[0];
    const long double error = (first[0]-(sum-part))+(second[0]-part)+first[1]+second[1];
    const long double value = sum+error;
    const long double tail = value-sum;
    return {value,(sum-(value-tail))+(error-tail)};
  }

  static std::array<long double,2> CompensatedProduct(const std::array<long double,2>& first, const std::array<long double,2>& second) {
    const long double product = first[0]*second[0];
    const long double error = std::fma(first[1],second[1],std::fma(first[0],second[1],std::fma(first[1],second[0],std::fma(first[0],second[0],-product))));
    return CompensatedSum({product,0},{error,0});
  }

  // Refines projected events against the original polynomial without cancellation in its discriminant.
  static std::array<std::array<long double,2>,3> CompensatedFiberPolynomial(const std::vector<std::vector<long double>>& coefficient, long double x) {
    std::array<std::array<long double,2>,3> polynomial = {{{0,0},{0,0},{0,0}}};
    for(unsigned j = 0; j < 3; ++j) for(unsigned i = 3; i > 0; --i) polynomial[j] = CompensatedSum(CompensatedProduct(polynomial[j],{x,0}),{coefficient[i-1U][j],0});
    return polynomial;
  }

  static long double FiberDiscriminantValue(const std::vector<std::vector<long double>>& coefficient, long double x) {
    const std::array<std::array<long double,2>,3> polynomial = CompensatedFiberPolynomial(coefficient,x);
    const std::array<long double,2> value = CompensatedSum(CompensatedProduct(polynomial[1],polynomial[1]),CompensatedProduct({-4,0},CompensatedProduct(polynomial[2],polynomial[0])));
    return value[0]+value[1];
  }

  static long double BoundaryValue(const std::vector<std::vector<long double>>& coefficient, long double x, long double constant, long double slope) {
    const std::array<std::array<long double,2>,3> polynomial = CompensatedFiberPolynomial(coefficient,x);
    const std::array<long double,2> point = CompensatedSum({constant,0},CompensatedProduct({slope,0},{x,0}));
    const std::array<long double,2> value = CompensatedSum(CompensatedProduct(CompensatedSum(CompensatedProduct(polynomial[2],point),polynomial[1]),point),polynomial[0]);
    return value[0]+value[1];
  }

  static long double RefineBoundaryEvent(const std::vector<std::vector<long double>>& coefficient, long double root, long double left, long double right, long double constant, long double slope) {
    const std::vector<long double> derivative = LineSliceUtilities::PolynomialDerivative(BoundaryPolynomial(coefficient,constant,slope));
    long double current = root;
    for(unsigned iteration = 0; iteration < 12; ++iteration) {
      const long double value = BoundaryValue(coefficient,current,constant,slope);
      const long double gradient = LineSliceUtilities::PolynomialValue(derivative,current);
      if(gradient == 0 || !std::isfinite(value) || !std::isfinite(gradient)) break;
      const long double next = current-value/gradient;
      if(!(next >= left && next <= right) || !std::isfinite(next) || next == current) break;
      current = next;
    }
    return current;
  }

  static long double FiberDiscriminantBound(const std::vector<std::vector<long double>>& coefficient, long double x) {
    long double bound[3] = {0.0L,0.0L,0.0L};
    for(unsigned j = 0; j < 3; ++j) for(unsigned i = 3; i > 0; --i) bound[j] = std::fma(bound[j],std::abs(x),std::abs(coefficient[i-1U][j]));
    return bound[1]*bound[1]+4.0L*bound[2]*bound[0];
  }

  static long double FiberDiscriminantDerivative(const std::vector<std::vector<long double>>& coefficient, long double x) {
    long double a = 0.0L;
    long double b = 0.0L;
    long double c = 0.0L;
    long double da = 0.0L;
    long double db = 0.0L;
    long double dc = 0.0L;
    long double power = 1.0L;
    for(unsigned i = 0; i < 3; ++i) {
      a += coefficient[i][2]*power;
      b += coefficient[i][1]*power;
      c += coefficient[i][0]*power;
      if(i > 0) {
        const long double derivativePower = i == 1 ? 1.0L : x;
        da += static_cast<long double>(i)*coefficient[i][2]*derivativePower;
        db += static_cast<long double>(i)*coefficient[i][1]*derivativePower;
        dc += static_cast<long double>(i)*coefficient[i][0]*derivativePower;
      }
      power *= x;
    }
    return 2.0L*b*db-4.0L*(da*c+a*dc);
  }

  static long double RefineDiscriminantEvent(const std::vector<std::vector<long double>>& coefficient, long double root, long double left, long double right, const PreparedPlane* prepared = nullptr, unsigned preparedOrder = 1) {
    long double current = root;
    for(unsigned iteration = 0; iteration < 12; ++iteration) {
      const long double value = prepared == nullptr ? FiberDiscriminantValue(coefficient,current) : prepared->Discriminant(preparedOrder,current);
      const long double derivative = prepared == nullptr ? FiberDiscriminantDerivative(coefficient,current) : prepared->DiscriminantDerivative(preparedOrder,current);
      if(derivative == 0.0L || !std::isfinite(value) || !std::isfinite(derivative)) break;
      const long double next = current-value/derivative;
      if(!(next >= left && next <= right) || !std::isfinite(next)) break;
      if(next == current) break;
      current = next;
    }
    return current;
  }

  static void CollectEvents(const std::vector<std::vector<long double>>& coefficient, const LineIntegrationStrip& strip, std::vector<long double>& events, std::vector<long double>* tangencies = nullptr, bool certifiedGraph = false, const PreparedPlane* prepared = nullptr, unsigned preparedOrder = 1) {
    events = {strip.left,strip.right};
    std::vector<long double> lower = prepared == nullptr ? BoundaryPolynomial(coefficient,strip.lowerConstant,strip.lowerSlope) : prepared->BoundaryPolynomial(preparedOrder,strip.lowerConstant,strip.lowerSlope);
    std::vector<long double> upper = prepared == nullptr ? BoundaryPolynomial(coefficient,strip.upperConstant,strip.upperSlope) : prepared->BoundaryPolynomial(preparedOrder,strip.upperConstant,strip.upperSlope);
    std::vector<long double> boundaryEvents[2];
    LineSliceUtilities::PolynomialRoots(lower,strip.left,strip.right,boundaryEvents[0]);
    LineSliceUtilities::PolynomialRoots(upper,strip.left,strip.right,boundaryEvents[1]);
    std::vector<long double> a(3,0.0L);
    std::vector<long double> b(3,0.0L);
    std::vector<long double> c(3,0.0L);
    for(unsigned i = 0; i < 3; ++i) {
      a[i] = coefficient[i][2];
      b[i] = coefficient[i][1];
      c[i] = coefficient[i][0];
    }
    const std::vector<long double> discriminant = prepared == nullptr ? LineSliceUtilities::PolynomialSum(LineSliceUtilities::PolynomialProduct(b,b),LineSliceUtilities::PolynomialProduct(a,c),-4.0L) : prepared->DiscriminantPolynomial(preparedOrder);
    std::vector<long double> discriminantEvents;
    if(!certifiedGraph && (prepared == nullptr || !prepared->data->HasNodalPolynomialSource() || !prepared->data->GetCompensatedDiscriminantRoots(PreparedPlane::InnerAxis(preparedOrder),strip.left,strip.right,discriminantEvents))) LineSliceUtilities::PolynomialRoots(discriminant,strip.left,strip.right,discriminantEvents);
    for(unsigned k = 0; k < discriminantEvents.size(); ++k) {
      const long double event = RefineDiscriminantEvent(coefficient,discriminantEvents[k],strip.left,strip.right,prepared,preparedOrder);
      discriminantEvents[k] = event;
      events.push_back(event);
      if(tangencies != nullptr) tangencies->push_back(event);
    }
    for(unsigned side = 0; side < 2; ++side) for(unsigned k = 0; k < boundaryEvents[side].size(); ++k) {
      long double event = RefineBoundaryEvent(coefficient,boundaryEvents[side][k],strip.left,strip.right,side == 0 ? strip.lowerConstant : strip.upperConstant,side == 0 ? strip.lowerSlope : strip.upperSlope);
      const long double boundary = side == 0 ? StripLower(strip,event) : StripUpper(strip,event);
      const std::vector<long double> polynomial = FiberPolynomial(coefficient,event);
      const long double derivative = std::fma(2.0L*polynomial[2],boundary,polynomial[1]);
      const long double derivativeBound = 2.0L*std::abs(polynomial[2]*boundary)+std::abs(polynomial[1]);
      const long double padding = 512.0L*std::numeric_limits<long double>::epsilon();
      if(std::abs(derivative) <= padding*std::max(1.0L,derivativeBound)) {
        for(unsigned r = 0; r < discriminantEvents.size(); ++r) {
          if(std::abs(event-discriminantEvents[r]) <= padding*std::max(1.0L,std::abs(event))) {
            event = discriminantEvents[r];
            break;
          }
        }
      }
      events.push_back(event);
    }
    LineSliceUtilities::SortDistinct(events);
  }


};

}



#endif
