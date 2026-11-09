#ifndef FEMUS_LINE_DOMAIN_INTEGRATION_HPP
#define FEMUS_LINE_DOMAIN_INTEGRATION_HPP

#include "AdaptiveLineSliceMoments.hpp"
#include "FELevelSetMoments.hpp"

namespace FELevelSetMoments {

struct LineIntegrationStrip {
  long double left;
  long double right;
  long double lowerConstant;
  long double lowerSlope;
  long double upperConstant;
  long double upperSlope;
};

struct LineIntegrationStatistics {
  unsigned long long fastSamples = 0;
  unsigned long long fullSamples = 0;
  unsigned long long taylorPanels = 0;
  unsigned long long rootSolves = 0;
  unsigned long long inactiveSkipped = 0;
  unsigned long long fallbackComponents = 0;
};

class LineDomainIntegration {
public:
  // Shared affine-fiber integration retains accepted moments through every fallback stage.
  static LineIntegrationStatistics& Statistics() {
    static thread_local LineIntegrationStatistics data;
    return data;
  }

  static bool BoundaryIsZero(const std::vector<std::vector<long double>>& coefficient, long double constant, long double slope) {
    const std::vector<long double> boundary = BoundaryPolynomial(coefficient,constant,slope);
    return LineSliceUtilities::PolynomialNorm(boundary) == 0.0L;
  }

  static void Sample(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, long double x, bool computeDelta, bool graph, bool lowerZero, bool upperZero, std::vector<long double>& values, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>& active) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    values.assign(2U*count,0.0L);
    const long double lower = StripLower(strip,x);
    const long double upper = StripUpper(strip,x);
    const std::vector<long double> polynomial = FiberPolynomial(coefficient,x);
    std::vector<long double> roots;
    if(lowerZero) {
      roots.push_back(lower);
      if(polynomial[2] != 0.0L) {
        const long double other = -polynomial[1]/polynomial[2]-lower;
        if(other != lower) roots.push_back(other);
      }
    }
    else if(upperZero) {
      roots.push_back(upper);
      if(polynomial[2] != 0.0L) {
        const long double other = -polynomial[1]/polynomial[2]-upper;
        if(other != upper) roots.push_back(other);
      }
    }
    else roots = LineSliceUtilities::QuadraticRoots(polynomial[2],polynomial[1],polynomial[0]);
    ++Statistics().rootSolves;
    std::sort(roots.begin(),roots.end());
    std::vector<long double> inside;
    std::vector<long double> cuts(1,lower);
    for(unsigned r = 0; r < roots.size(); ++r) {
      if(roots[r] < lower || roots[r] > upper) continue;
      inside.push_back(roots[r]);
      if(roots[r] > lower && roots[r] < upper) cuts.push_back(roots[r]);
    }
    cuts.push_back(upper);
    if(graph && inside.size() != 1) throw std::runtime_error("The cached graph branch is not a single root.");
    std::vector<unsigned char> needH(degree+1U,0);
    std::vector<unsigned char> needD(degree+1U,0);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
      const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
      needH[j] = needH[j] || active[k] != 0;
      needD[j] = needD[j] || (computeDelta && active[count+k] != 0);
    }
    std::vector<long double> h(degree+1U,0.0L);
    std::vector<long double> d(degree+1U,0.0L);
    for(unsigned segment = 1; segment < cuts.size(); ++segment) {
      const long double a = cuts[segment-1U];
      const long double b = cuts[segment];
      const long double middle = a+0.5L*(b-a);
      const long double sign = std::fma(std::fma(polynomial[2],middle,polynomial[1]),middle,polynomial[0]);
      if(sign <= 0.0L) continue;
      long double power = 1.0L;
      long double geometric = 0.0L;
      for(unsigned j = 0; j <= degree; ++j) {
        geometric = power+a*geometric;
        if(needH[j]) h[j] += (b-a)*geometric/(j+1U);
        power *= b;
      }
    }
    if(computeDelta) for(unsigned r = 0; r < inside.size(); ++r) {
      const long double root = inside[r];
      const long double derivative = std::abs(std::fma(2.0L*polynomial[2],root,polynomial[1]));
      if(!(derivative > 0.0L)) throw std::runtime_error("A singular fiber has no regular delta weight.");
      long double weight = (root == lower || root == upper ? 0.5L : 1.0L)/derivative;
      if(surfaceGradient != nullptr) weight *= SurfaceGradientNorm(*surfaceGradient,x,root);
      long double power = 1.0L;
      for(unsigned j = 0; j <= degree; ++j) {
        if(needD[j]) d[j] += power*weight;
        power *= root;
      }
    }
    long double xPower = 1.0L;
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
        if(active[k]) values[k] = xPower*h[j];
        else ++Statistics().inactiveSkipped;
        if(active[count+k]) values[count+k] = xPower*d[j];
        else ++Statistics().inactiveSkipped;
      }
      xPower *= x;
    }
  }

  static void FastGraphSample(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, long double x, bool computeDelta, bool lowerZero, bool upperZero, std::vector<long double>& values, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>& active) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    values.assign(2U*count,0.0L);
    const long double lower = StripLower(strip,x);
    const long double upper = StripUpper(strip,x);
    const std::vector<long double> polynomial = FiberPolynomial(coefficient,x);
    std::vector<long double> roots;
    if(lowerZero) {
      roots.push_back(lower);
      if(polynomial[2] != 0.0L) {
        const long double other = -polynomial[1]/polynomial[2]-lower;
        if(other != lower) roots.push_back(other);
      }
    }
    else if(upperZero) {
      roots.push_back(upper);
      if(polynomial[2] != 0.0L) {
        const long double other = -polynomial[1]/polynomial[2]-upper;
        if(other != upper) roots.push_back(other);
      }
    }
    else roots = LineSliceUtilities::QuadraticRoots(polynomial[2],polynomial[1],polynomial[0]);
    ++Statistics().rootSolves;
    const long double tolerance = 256.0L*std::numeric_limits<long double>::epsilon()*std::max(1.0L,std::max(std::abs(lower),std::abs(upper)));
    long double root = 0.0L;
    unsigned number = 0;
    for(unsigned r = 0; r < roots.size(); ++r) {
      if(roots[r] >= lower-tolerance && roots[r] <= upper+tolerance) {
        root = std::max(lower,std::min(upper,roots[r]));
        ++number;
      }
    }
    if(number > 1) throw std::runtime_error("The cached graph branch is not a single root.");
    long double positiveLeft = lower;
    long double positiveRight = lower;
    long double derivative = 0.0L;
    if(number == 1) {
      const long double signedDerivative = std::fma(2.0L*polynomial[2],root,polynomial[1]);
      derivative = std::abs(signedDerivative);
      if(!(derivative > 0.0L)) throw std::runtime_error("A singular fiber has no regular graph delta weight.");
      positiveLeft = signedDerivative > 0.0L ? root : lower;
      positiveRight = signedDerivative > 0.0L ? upper : root;
    }
    else if(Evaluate(coefficient,x,0.5L*(lower+upper)) > 0.0L) {
      positiveLeft = lower;
      positiveRight = upper;
    }
    std::vector<long double> h(degree+1U,0.0L);
    std::vector<long double> d(degree+1U,0.0L);
    long double leftPower = 1.0L;
    long double rightPower = 1.0L;
    for(unsigned j = 0; j <= degree; ++j) {
      leftPower *= positiveLeft;
      rightPower *= positiveRight;
      h[j] = (rightPower-leftPower)/(j+1U);
    }
    if(computeDelta && number == 1) {
      long double weight = (root == lower || root == upper ? 0.5L : 1.0L)/derivative;
      if(surfaceGradient != nullptr) weight *= SurfaceGradientNorm(*surfaceGradient,x,root);
      long double rootPower = 1.0L;
      for(unsigned j = 0; j <= degree; ++j) {
        d[j] = rootPower*weight;
        rootPower *= root;
      }
    }
    long double xPower = 1.0L;
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
        if(active[k]) values[k] = xPower*h[j];
        else ++Statistics().inactiveSkipped;
        if(active[count+k]) values[count+k] = xPower*d[j];
        else ++Statistics().inactiveSkipped;
      }
      xPower *= x;
    }
  }

  static long double MappedPoint(long double left, long double right, long double parameter, bool mapped) {
    if(!mapped) return left+(right-left)*parameter;
    const long double angle = 0.5L*std::acos(-1.0L)*parameter;
    const long double s = parameter <= 0.5L ? std::sin(angle) : std::cos(angle);
    return parameter <= 0.5L ? left+(right-left)*s*s : right-(right-left)*s*s;
  }

  static bool GaussKronrodPanel(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, long double left, long double right, long double parameterLeft, long double parameterRight, bool mapped, bool graph, bool computeDelta, bool lowerZero, bool upperZero, const std::vector<unsigned char>& active, unsigned maximumEvaluations, unsigned& evaluations, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient) {
    static const long double nodes[8] = {0.991455371120812639206854697526329L,0.949107912342758524526189684047851L,0.864864423359769072789712788640926L,0.741531185599394439863864773280788L,0.586087235467691130294144838258730L,0.405845151377397166906606412076961L,0.207784955007898467600689403773245L,0.0L};
    static const long double kronrod[8] = {0.022935322010529224963732008058970L,0.063092092629978553290700663189204L,0.104790010322250183839876322541518L,0.140653259715525918745189590510238L,0.169004726639267902826583426598550L,0.190350578064785409913256402421014L,0.204432940075298892414161999234649L,0.209482141084727828012999174891714L};
    static const long double gauss[4] = {0.129484966168869693270611432679082L,0.279705391489276667901467771423780L,0.381830050505118944950369775488975L,0.417959183673469387755102040816327L};
    const unsigned size = active.size();
    values.assign(size,0.0L);
    errors.assign(size,std::numeric_limits<long double>::infinity());
    std::vector<long double> low(size,0.0L);
    std::vector<long double> absolute(size,0.0L);
    std::vector<std::vector<long double>> samples(15);
    const long double middle = 0.5L*(parameterLeft+parameterRight);
    const long double halfWidth = 0.5L*(parameterRight-parameterLeft);
    const long double pi = std::acos(-1.0L);
    try {
      for(unsigned node = 0; node < 8; ++node) for(unsigned side = 0; side < (node == 7 ? 1U : 2U); ++side) {
        if(evaluations >= maximumEvaluations) throw std::runtime_error("LINE evaluation budget exhausted.");
        ++evaluations;
        if(graph) ++Statistics().fastSamples;
        else ++Statistics().fullSamples;
        const unsigned sample = 2U*node+side;
        const long double parameter = middle+(side == 0 ? -1.0L : 1.0L)*halfWidth*nodes[node];
        const long double point = MappedPoint(left,right,parameter,mapped);
        const long double jacobian = mapped ? 0.5L*pi*(right-left)*std::sin(pi*parameter) : right-left;
        if(graph) FastGraphSample(coefficient,degree,strip,point,computeDelta,lowerZero,upperZero,samples[sample],surfaceGradient,active);
        else Sample(coefficient,degree,strip,point,computeDelta,false,lowerZero,upperZero,samples[sample],surfaceGradient,active);
        for(unsigned k = 0; k < size; ++k) {
          if(!active[k]) continue;
          samples[sample][k] *= jacobian;
          if(!std::isfinite(samples[sample][k])) throw std::runtime_error("Nonfinite LINE sample.");
          values[k] += kronrod[node]*samples[sample][k];
          absolute[k] += kronrod[node]*std::abs(samples[sample][k]);
          if(node&1U) low[k] += gauss[node/2U]*samples[sample][k];
        }
      }
    }
    catch(const std::runtime_error&) {
      std::fill(errors.begin(),errors.end(),std::numeric_limits<long double>::infinity());
      return false;
    }
    for(unsigned k = 0; k < size; ++k) {
      if(!active[k]) continue;
      long double variation = 0.0L;
      for(unsigned sample = 0; sample < 15; ++sample) variation += kronrod[sample/2U]*std::abs(samples[sample][k]-0.5L*values[k]);
      long double estimate = std::abs(values[k]-low[k]);
      if(variation > 0.0L && estimate > 0.0L) {
        const long double ratio = 200.0L*estimate/variation;
        estimate = ratio >= 1.0L ? variation : variation*ratio*std::sqrt(ratio);
      }
      errors[k] = halfWidth*std::max(estimate,128.0L*std::numeric_limits<long double>::epsilon()*absolute[k]);
      values[k] *= halfWidth;
    }
    return true;
  }

  static bool Accept(const std::vector<long double>& values, const std::vector<long double>& errors, long double absoluteH, long double absoluteD, long double relative, std::vector<unsigned char>& active, std::vector<long double>& total, std::vector<long double>& totalError) {
    bool remaining = false;
    const unsigned count = active.size()/2U;
    for(unsigned k = 0; k < active.size(); ++k) {
      if(!active[k]) continue;
      const long double tolerance = (k < count ? absoluteH : absoluteD)+relative*std::abs(values[k]);
      if(std::isfinite(values[k]) && std::isfinite(errors[k]) && errors[k] >= 0.0L && errors[k] <= tolerance) {
        total[k] += values[k];
        totalError[k] += errors[k];
        active[k] = 0;
      }
      else remaining = true;
    }
    return remaining;
  }

  static void Adaptive(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, long double left, long double right, long double parameterLeft, long double parameterRight, bool mapped, bool graph, bool useTaylor, bool computeDelta, bool lowerZero, bool upperZero, long double absoluteH, long double absoluteD, long double relative, unsigned depth, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, std::vector<long double>& total, std::vector<long double>& totalError, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>& requested) {
    std::vector<unsigned char> active(requested);
    std::vector<long double> values;
    std::vector<long double> errors;
    const long double fraction = (right-left)*(parameterRight-parameterLeft);
    const bool sampled = GaussKronrodPanel(coefficient,degree,strip,left,right,parameterLeft,parameterRight,mapped,graph,computeDelta,lowerZero,upperZero,active,maximumEvaluations,evaluations,values,errors,surfaceGradient);
    if(!Accept(values,errors,absoluteH*fraction,absoluteD*fraction,relative,active,total,totalError)) {
      ++accepted;
      return;
    }
    const long double middle = 0.5L*(parameterLeft+parameterRight);
    const unsigned fastDepthLimit = mapped ? maximumDepth : std::min(maximumDepth,8U);
    if(graph && sampled && depth < fastDepthLimit && middle != parameterLeft && middle != parameterRight) {
      Adaptive(coefficient,degree,strip,left,right,parameterLeft,middle,mapped,true,useTaylor,computeDelta,lowerZero,upperZero,absoluteH,absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError,surfaceGradient,active);
      Adaptive(coefficient,degree,strip,left,right,middle,parameterRight,mapped,true,useTaylor,computeDelta,lowerZero,upperZero,absoluteH,absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError,surfaceGradient,active);
      return;
    }
    if(graph && useTaylor) {
      ++Statistics().taylorPanels;
      for(unsigned k = 0; k < active.size(); ++k) Statistics().fallbackComponents += active[k] != 0;
      try {
        TaylorPanel(coefficient,degree,strip,MappedPoint(left,right,parameterLeft,mapped),MappedPoint(left,right,parameterRight,mapped),computeDelta,values,errors,surfaceGradient,&active);
        if(!Accept(values,errors,absoluteH*fraction,absoluteD*fraction,relative,active,total,totalError)) {
          ++accepted;
          return;
        }
      }
      catch(const std::runtime_error&) {}
    }
    if(graph) {
      GaussKronrodPanel(coefficient,degree,strip,left,right,parameterLeft,parameterRight,mapped,false,computeDelta,lowerZero,upperZero,active,maximumEvaluations,evaluations,values,errors,surfaceGradient);
      if(!Accept(values,errors,absoluteH*fraction,absoluteD*fraction,relative,active,total,totalError)) {
        ++accepted;
        return;
      }
    }
    if(depth >= maximumDepth || middle == parameterLeft || middle == parameterRight || evaluations >= maximumEvaluations) throw std::runtime_error("LINE domain integration did not reach the requested tolerance.");
    Adaptive(coefficient,degree,strip,left,right,parameterLeft,middle,mapped,false,false,computeDelta,lowerZero,upperZero,absoluteH,absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError,surfaceGradient,active);
    Adaptive(coefficient,degree,strip,left,right,middle,parameterRight,mapped,false,false,computeDelta,lowerZero,upperZero,absoluteH,absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError,surfaceGradient,active);
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

  static void FullStripMoments(unsigned degree, const LineIntegrationStrip& strip, long double left, long double right, std::vector<long double>& values, std::vector<long double>& errors) {
    std::vector<long double> points;
    std::vector<long double> weights;
    GaussRule((degree+3U)/2U,points,weights);
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    for(unsigned q = 0; q < points.size(); ++q) {
      const long double x = left+(right-left)*points[q];
      const long double lower = StripLower(strip,x);
      const long double upper = StripUpper(strip,x);
      const long double weight = (right-left)*weights[q];
      long double xPower = 1.0L;
      std::vector<long double> fiber(degree+1U);
      long double geometric = 0.0L;
      long double power = 1.0L;
      for(unsigned j = 0; j <= degree; ++j) {
        geometric = power+lower*geometric;
        fiber[j] = (upper-lower)*geometric/(j+1U);
        power *= upper;
      }
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) {
          const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
          const long double term = weight*xPower*fiber[j];
          values[k] += term;
          errors[k] += (128.0L+32.0L*degree)*std::numeric_limits<long double>::epsilon()*std::abs(term);
        }
        xPower *= x;
      }
    }
  }

  static void IntegrateStrip(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, bool graph, bool mapped, bool useTaylor, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
    std::vector<long double> events;
    CollectEvents(coefficient,strip,events);
    const bool lowerZero = BoundaryIsZero(coefficient,strip.lowerConstant,strip.lowerSlope);
    const bool upperZero = BoundaryIsZero(coefficient,strip.upperConstant,strip.upperSlope);
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    const std::vector<unsigned char> active(2U*count,1);
    for(unsigned interval = 1; interval < events.size(); ++interval) {
      const long double left = events[interval-1U];
      const long double right = events[interval];
      if(!(right > left)) continue;
      const long double middle = left+0.5L*(right-left);
      const long double lower = StripLower(strip,middle);
      const long double upper = StripUpper(strip,middle);
      const std::vector<long double> polynomial = FiberPolynomial(coefficient,middle);
      const std::vector<long double> roots = LineSliceUtilities::QuadraticRoots(polynomial[2],polynomial[1],polynomial[0]);
      unsigned number = 0;
      for(unsigned r = 0; r < roots.size(); ++r) if(roots[r] >= lower && roots[r] <= upper) ++number;
      if(number == 0 && !lowerZero && !upperZero) {
        if(Evaluate(coefficient,middle,0.5L*(lower+upper)) > 0.0L) FullStripMoments(degree,strip,left,right,values,errors);
        ++accepted;
        continue;
      }
      Adaptive(coefficient,degree,strip,left,right,0.0L,1.0L,mapped,graph && number == 1,useTaylor,computeDelta,lowerZero,upperZero,absoluteH,absoluteD,relative,0,maximumDepth,maximumEvaluations,evaluations,accepted,values,errors,surfaceGradient,active);
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

  static void IntegrateCell(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned order, bool triangle, bool computeDelta, long double derivativeBound, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
    if(computeDelta && HasVerticalFiber(coefficient,order)) {
      LineSliceDirection(coefficient,degree,triangle,order,computeDelta,false,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,0,nullptr,surfaceGradient);
      if(order == 2) {
        const unsigned count = (degree+1U)*(degree+2U)/2U;
        for(unsigned i = 0; i <= degree; ++i) for(unsigned j = i+1U; j <= degree-i; ++j) {
          const unsigned a = LineSliceUtilities::MomentIndex(degree,i,j);
          const unsigned b = LineSliceUtilities::MomentIndex(degree,j,i);
          for(unsigned pair = 0; pair < 2; ++pair) {
            std::swap(values[pair*count+a],values[pair*count+b]);
            std::swap(errors[pair*count+a],errors[pair*count+b]);
          }
        }
      }
      return;
    }
    std::vector<std::vector<long double>> working(coefficient);
    std::vector<std::vector<std::vector<long double>>> gradient;
    if(order == 2) {
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = i+1U; j < 3; ++j) std::swap(working[i][j],working[j][i]);
      if(surfaceGradient != nullptr) {
        SurfaceGradientTranspose(*surfaceGradient,gradient);
        surfaceGradient = &gradient;
      }
    }
    const LineIntegrationStrip strip = {0.0L,1.0L,0.0L,triangle && order == 2 ? 1.0L : 0.0L,triangle && order == 1 ? 0.0L : 1.0L,triangle && order == 1 ? 1.0L : 0.0L};
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    values.assign(2U*count,0.0L);
    errors.assign(2U*count,0.0L);
    IntegrateStrip(working,degree,strip,true,derivativeBound == 0.0L,true,computeDelta,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,evaluations,accepted,values,errors,surfaceGradient);
    if(order == 2) for(unsigned i = 0; i <= degree; ++i) for(unsigned j = i+1U; j <= degree-i; ++j) {
      const unsigned a = LineSliceUtilities::MomentIndex(degree,i,j);
      const unsigned b = LineSliceUtilities::MomentIndex(degree,j,i);
      for(unsigned pair = 0; pair < 2; ++pair) {
        std::swap(values[pair*count+a],values[pair*count+b]);
        std::swap(errors[pair*count+a],errors[pair*count+b]);
      }
    }
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

  static void CollectEvents(const std::vector<std::vector<long double>>& coefficient, const LineIntegrationStrip& strip, std::vector<long double>& events) {
    events = {strip.left,strip.right};
    std::vector<long double> lower = BoundaryPolynomial(coefficient,strip.lowerConstant,strip.lowerSlope);
    std::vector<long double> upper = BoundaryPolynomial(coefficient,strip.upperConstant,strip.upperSlope);
    LineSliceUtilities::PolynomialRoots(lower,strip.left,strip.right,events);
    LineSliceUtilities::PolynomialRoots(upper,strip.left,strip.right,events);
    std::vector<long double> a(3,0.0L);
    std::vector<long double> b(3,0.0L);
    std::vector<long double> c(3,0.0L);
    for(unsigned i = 0; i < 3; ++i) {
      a[i] = coefficient[i][2];
      b[i] = coefficient[i][1];
      c[i] = coefficient[i][0];
    }
    std::vector<long double> discriminant = LineSliceUtilities::PolynomialSum(LineSliceUtilities::PolynomialProduct(b,b),LineSliceUtilities::PolynomialProduct(a,c),-4.0L);
    LineSliceUtilities::PolynomialRoots(discriminant,strip.left,strip.right,events);
    LineSliceUtilities::SortDistinct(events);
  }

  static void TaylorPanel(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, long double left, long double right, bool computeDelta, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>* active) {
    const unsigned order = 32;
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    const unsigned size = 2U*count;
    const long double thetaLeft = 0.0L;
    const long double thetaRight = 0.5L*std::acos(-1.0L);
    const long double thetaMiddle = 0.5L*(thetaLeft+thetaRight);
    const long double thetaHalfWidth = 0.5L*(thetaRight-thetaLeft);
    std::vector<long double> sine;
    std::vector<long double> cosine;
    SineCosineSeries(thetaMiddle,thetaHalfWidth,order,sine,cosine);
    std::vector<long double> x = SeriesProduct(sine,sine);
    const long double width = right-left;
    for(unsigned k = 0; k <= order; ++k) x[k] *= width;
    x[0] += left;
    std::vector<long double> jacobian = SeriesProduct(sine,cosine);
    for(unsigned k = 0; k <= order; ++k) jacobian[k] *= 2.0L*width*thetaHalfWidth;
    std::vector<long double> aPolynomial(3,0.0L);
    std::vector<long double> bPolynomial(3,0.0L);
    std::vector<long double> cPolynomial(3,0.0L);
    for(unsigned i = 0; i < 3; ++i) {
      aPolynomial[i] = coefficient[i][2];
      bPolynomial[i] = coefficient[i][1];
      cPolynomial[i] = coefficient[i][0];
    }
    std::vector<long double> a = SeriesPolynomial(aPolynomial,x);
    std::vector<long double> b = SeriesPolynomial(bPolynomial,x);
    std::vector<long double> c = SeriesPolynomial(cPolynomial,x);
    std::vector<long double> root(order+1U,0.0L);
    const bool quadratic = LineSliceUtilities::PolynomialNorm(aPolynomial) != 0.0L;
    if(!quadratic) {
      root = SeriesQuotient(c,b);
      for(unsigned k = 0; k <= order; ++k) root[k] = -root[k];
    }
    else {
      if(a[0] == 0.0L) throw std::runtime_error("Triangle Taylor quadratic coefficient vanishes at the expansion point.");
      std::vector<long double> discriminant = SeriesProduct(b,b);
      std::vector<long double> term = SeriesProduct(a,c);
      for(unsigned k = 0; k <= order; ++k) discriminant[k] -= 4.0L*term[k];
      std::vector<long double> squareRoot = SeriesSquareRoot(discriminant);
      const long double lowerMiddle = StripLower(strip,0.5L*(left+right));
      const long double upperMiddle = StripUpper(strip,0.5L*(left+right));
      const long double candidateMinus = (-b[0]-squareRoot[0])/(2.0L*a[0]);
      const long double candidatePlus = (-b[0]+squareRoot[0])/(2.0L*a[0]);
      const long double sign = candidateMinus >= lowerMiddle && candidateMinus <= upperMiddle ? -1.0L : (candidatePlus >= lowerMiddle && candidatePlus <= upperMiddle ? 1.0L : 0.0L);
      if(sign == 0.0L) throw std::runtime_error("Triangle Taylor graph root is outside the strip.");
      std::vector<long double> numerator(b);
      std::vector<long double> denominator(a);
      for(unsigned k = 0; k <= order; ++k) {
        numerator[k] = -numerator[k]+sign*squareRoot[k];
        denominator[k] *= 2.0L;
      }
      root = SeriesQuotient(numerator,denominator);
    }
    std::vector<long double> lowerSeries = x;
    std::vector<long double> upperSeries = x;
    for(unsigned k = 0; k <= order; ++k) {
      lowerSeries[k] *= strip.lowerSlope;
      upperSeries[k] *= strip.upperSlope;
    }
    lowerSeries[0] += strip.lowerConstant;
    upperSeries[0] += strip.upperConstant;
    std::vector<long double> gy = SeriesProduct(a,root);
    for(unsigned k = 0; k <= order; ++k) gy[k] = 2.0L*gy[k]+b[k];
    const long double signedDerivative = gy[0];
    if(signedDerivative == 0.0L) throw std::runtime_error("Triangle Taylor graph has a singular derivative.");
    if(signedDerivative < 0.0L) for(unsigned k = 0; k <= order; ++k) gy[k] = -gy[k];
    std::vector<long double> deltaWeight;
    if(computeDelta) {
      deltaWeight = SeriesQuotient(jacobian,gy);
      const bool boundaryRoot = BoundaryIsZero(coefficient,strip.lowerConstant,strip.lowerSlope) || BoundaryIsZero(coefficient,strip.upperConstant,strip.upperSlope);
      if(boundaryRoot) for(unsigned k = 0; k <= order; ++k) deltaWeight[k] *= 0.5L;
      if(surfaceGradient != nullptr) deltaWeight = SeriesProduct(deltaWeight,SurfaceGradientSeriesNorm(*surfaceGradient,x,root));
    }
    std::vector<std::vector<long double>> lowerPowers(degree+2U,SeriesConstant(1.0L,order));
    std::vector<std::vector<long double>> upperPowers(degree+2U,SeriesConstant(1.0L,order));
    std::vector<std::vector<long double>> rootPowers(degree+2U,SeriesConstant(1.0L,order));
    for(unsigned j = 1; j <= degree+1U; ++j) {
      lowerPowers[j] = SeriesProduct(lowerPowers[j-1U],lowerSeries);
      upperPowers[j] = SeriesProduct(upperPowers[j-1U],upperSeries);
      rootPowers[j] = SeriesProduct(rootPowers[j-1U],root);
    }
    std::vector<std::vector<long double>> xPowers(degree+1U,SeriesConstant(1.0L,order));
    for(unsigned i = 1; i <= degree; ++i) xPowers[i] = SeriesProduct(xPowers[i-1U],x);
    values.assign(size,0.0L);
    errors.assign(size,0.0L);
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned index = LineSliceUtilities::MomentIndex(degree,i,j);
        if(active == nullptr || (*active)[index] != 0) {
          std::vector<long double> fiber(order+1U,0.0L);
          const std::vector<long double>& first = signedDerivative > 0.0L ? upperPowers[j+1U] : rootPowers[j+1U];
          const std::vector<long double>& second = signedDerivative > 0.0L ? rootPowers[j+1U] : lowerPowers[j+1U];
          for(unsigned k = 0; k <= order; ++k) fiber[k] = (first[k]-second[k])/(j+1U);
          fiber = SeriesProduct(fiber,jacobian);
          bool converged = true;
          SeriesProductIntegral(xPowers[i],fiber,values[index],errors[index],converged);
          if(!converged) errors[index] = std::numeric_limits<long double>::infinity();
        }
        if(computeDelta && (active == nullptr || (*active)[count+index] != 0)) {
          std::vector<long double> fiber = SeriesProduct(rootPowers[j],deltaWeight);
          bool converged = true;
          SeriesProductIntegral(xPowers[i],fiber,values[count+index],errors[count+index],converged);
          if(!converged) errors[count+index] = std::numeric_limits<long double>::infinity();
        }
      }
    }
  }
};

}

#endif
