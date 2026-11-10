#ifndef FEMUS_LINE_DOMAIN_INTEGRATION_HPP
#define FEMUS_LINE_DOMAIN_INTEGRATION_HPP

#include "LineSliceUtilities.hpp"
#include "FELevelSetMoments.hpp"
#include <array>

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
  unsigned long long fullSamples = 0;
  unsigned long long rootSolves = 0;
  unsigned long long singleRootSamples = 0;
  unsigned long long inactiveSkipped = 0;
  unsigned long long mappedPanels = 0;
  unsigned long long regularGraphPanels = 0;
  unsigned long long pairPanels = 0;
  long double worstPairRatio = 0.0L;
  unsigned worstPairIndex = std::numeric_limits<unsigned>::max();
  long double worstPairValue = 0.0L;
  long double worstPairError = 0.0L;
  long double worstPairTolerance = 0.0L;
  long double worstPairLeft = 0.0L;
  long double worstPairRight = 0.0L;
};

struct LineIntegrationDiscriminant {
  std::vector<long double> coefficient;
  bool leftZero = false;
  bool rightZero = false;
};

class LineDomainIntegration {
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

  static void FastGraphSample(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, long double x, bool computeDelta, bool lowerZero, bool upperZero, std::vector<long double>& values, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>& active, long double mappedDerivative = -1.0L, const std::array<long double,5>* coordinates = nullptr, int derivativeSign = 0) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    values.assign(2U*count,0.0L);
    const long double parentX = coordinates == nullptr ? x : (*coordinates)[0]+(*coordinates)[1]*x;
    const long double measure = coordinates == nullptr ? 1.0L : (*coordinates)[4];
    const long double lower = StripLower(strip,x);
    const long double upper = StripUpper(strip,x);
    std::array<long double,3> polynomial;
    for(unsigned j = 0; j < 3; ++j) polynomial[j] = coefficient[2][j] == 0.0L ? coefficient[0][j]+coefficient[1][j]*x : std::fma(std::fma(coefficient[2][j],x,coefficient[1][j]),x,coefficient[0][j]);
    ++Statistics().rootSolves;
    const long double tolerance = 256.0L*std::numeric_limits<long double>::epsilon()*std::max(1.0L,std::max(std::abs(lower),std::abs(upper)));
    long double root = 0.0L;
    long double derivative = 0.0L;
    long double signedDerivative = 0.0L;
    unsigned number = 0;
    if(derivativeSign != 0) {
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
        const long double discriminant = std::fma(-4.0L*polynomial[2],polynomial[0],polynomial[1]*polynomial[1]);
        if(discriminant >= 0.0L) {
          derivative = std::sqrt(discriminant);
          const long double denominator = polynomial[1]+derivativeSign*derivative;
          root = derivativeSign*polynomial[1] >= 0.0L && denominator != 0.0L ? -2.0L*polynomial[0]/denominator : 0.5L*(-polynomial[1]+derivativeSign*derivative)/polynomial[2];
          available = true;
        }
      }
      if(available && root >= lower-tolerance && root <= upper+tolerance) { root = std::max(lower,std::min(upper,root)); number = 1; }
      signedDerivative = derivativeSign;
    }
    else {
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
        const long double discriminant = std::fma(-4.0L*polynomial[2],polynomial[0],polynomial[1]*polynomial[1]);
        if(discriminant == 0.0L) roots[rootsCount++] = -polynomial[1]/(2.0L*polynomial[2]);
        else if(discriminant > 0.0L) {
          const long double numerator = -0.5L*(polynomial[1]+std::copysign(std::sqrt(discriminant),polynomial[1]));
          roots[rootsCount++] = numerator/polynomial[2];
          roots[rootsCount++] = polynomial[0]/numerator;
        }
      }
      for(unsigned r = 0; r < rootsCount; ++r) if(roots[r] >= lower-tolerance && roots[r] <= upper+tolerance) { root = std::max(lower,std::min(upper,roots[r])); ++number; }
      if(number > 1) throw std::runtime_error("The cached graph branch is not a single root.");
      if(number == 1) {
        signedDerivative = polynomial[2] == 0.0L ? polynomial[1] : std::fma(2.0L*polynomial[2],root,polynomial[1]);
        const long double discriminant = polynomial[2] == 0.0L ? polynomial[1]*polynomial[1] : std::fma(-4.0L*polynomial[2],polynomial[0],polynomial[1]*polynomial[1]);
        derivative = polynomial[2] == 0.0L ? std::abs(polynomial[1]) : (discriminant >= 0.0L ? std::sqrt(discriminant) : std::abs(signedDerivative));
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
    else if(std::fma(std::fma(polynomial[2],0.5L*(lower+upper),polynomial[1]),0.5L*(lower+upper),polynomial[0]) > 0.0L) positiveRight = upper;
    const long double parentLeft = coordinates == nullptr ? positiveLeft : (*coordinates)[2]+(*coordinates)[3]*positiveLeft;
    const long double parentRight = coordinates == nullptr ? positiveRight : (*coordinates)[2]+(*coordinates)[3]*positiveRight;
    long double rightPower = 1.0L;
    long double geometric = 0.0L;
    for(unsigned j = 0; j <= degree; ++j) {
      geometric = rightPower+parentLeft*geometric;
      values[j] = measure*(positiveRight-positiveLeft)*geometric/(j+1U);
      rightPower *= parentRight;
    }
    if(computeDelta && number == 1) {
      long double weight = measure*((root == lower || root == upper) && !strip.fullBoundaryDelta ? 0.5L : 1.0L)/derivative;
      if(surfaceGradient != nullptr) weight *= SurfaceGradientNorm(*surfaceGradient,x,root);
      const long double parentRoot = coordinates == nullptr ? root : (*coordinates)[2]+(*coordinates)[3]*root;
      long double rootPower = 1.0L;
      for(unsigned j = 0; j <= degree; ++j) { values[count+j] = rootPower*weight; rootPower *= parentRoot; }
    }
    long double xPower = 1.0L;
    for(unsigned i = 1; i <= degree; ++i) {
      xPower *= parentX;
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
        if(active[k]) values[k] = xPower*values[j];
        else ++Statistics().inactiveSkipped;
        if(active[count+k]) values[count+k] = xPower*values[count+j];
        else ++Statistics().inactiveSkipped;
      }
    }
  }

  static long double MappedPoint(long double left, long double right, long double parameter, bool mapped) {
    if(!mapped) return left+(right-left)*parameter;
    const long double angle = 0.5L*std::acos(-1.0L)*parameter;
    const long double s = parameter <= 0.5L ? std::sin(angle) : std::cos(angle);
    return parameter <= 0.5L ? left+(right-left)*s*s : right-(right-left)*s*s;
  }

  static LineIntegrationDiscriminant MappedDiscriminant(const std::vector<std::vector<long double>>& coefficient, long double left, long double right, bool leftCandidate = true, bool rightCandidate = true) {
    LineIntegrationDiscriminant result;
    std::vector<std::vector<long double>> shifted(3,std::vector<long double>(3,0.0L));
    const long double width = right-left;
    for(unsigned j = 0; j < 3; ++j) {
      shifted[0][j] = std::fma(std::fma(coefficient[2][j],left,coefficient[1][j]),left,coefficient[0][j]);
      shifted[1][j] = width*std::fma(2.0L*coefficient[2][j],left,coefficient[1][j]);
      shifted[2][j] = width*width*coefficient[2][j];
    }
    result.coefficient.assign(5,0.0L);
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) {
      result.coefficient[i+j] = std::fma(shifted[i][1],shifted[j][1],result.coefficient[i+j]);
      result.coefficient[i+j] = std::fma(-4.0L*shifted[i][2],shifted[j][0],result.coefficient[i+j]);
    }
    const long double firstBound = FiberDiscriminantBound(coefficient,left);
    const long double lastBound = FiberDiscriminantBound(coefficient,right);
    const long double padding = 64.0L*std::numeric_limits<long double>::epsilon();
    result.leftZero = leftCandidate && std::abs(result.coefficient[0]) <= padding*firstBound;
    result.rightZero = rightCandidate && std::abs(LineSliceUtilities::PolynomialValue(result.coefficient,1.0L)) <= padding*lastBound;
    LineSliceUtilities::Trim(result.coefficient);
    if(result.leftZero && result.coefficient.size() > 1) result.coefficient.erase(result.coefficient.begin());
    else result.leftZero = false;
    if(result.rightZero && result.coefficient.size() > 1) {
      std::vector<long double> quotient(result.coefficient.size()-1U,0.0L);
      quotient.back() = -result.coefficient.back();
      for(unsigned k = quotient.size()-1U; k > 0; --k) quotient[k-1U] = quotient[k]-result.coefficient[k];
      result.coefficient.swap(quotient);
    }
    else result.rightZero = false;
    return result;
  }

  static long double MappedDerivativeMagnitude(const LineIntegrationDiscriminant& discriminant, long double parameter) {
    if(!discriminant.leftZero && !discriminant.rightZero) return -1.0L;
    const long double angle = 0.5L*std::acos(-1.0L)*parameter;
    const long double sine = std::sin(angle);
    const long double cosine = std::cos(angle);
    const long double point = parameter <= 0.5L ? sine*sine : 1.0L-cosine*cosine;
    long double value = LineSliceUtilities::PolynomialValue(discriminant.coefficient,point);
    if(discriminant.leftZero) value *= sine*sine;
    if(discriminant.rightZero) value *= cosine*cosine;
    return value >= 0.0L ? std::sqrt(value) : -1.0L;
  }

  static bool GaussKronrodPanel(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, long double left, long double right, long double parameterLeft, long double parameterRight, bool mapped, bool graph, bool computeDelta, bool lowerZero, bool upperZero, const std::vector<unsigned char>& active, unsigned maximumEvaluations, unsigned& evaluations, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const LineIntegrationDiscriminant* cachedDiscriminant = nullptr, const std::array<long double,5>* coordinates = nullptr, int derivativeSign = 0) {
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
    const LineIntegrationDiscriminant localDiscriminant = mapped && cachedDiscriminant == nullptr ? MappedDiscriminant(coefficient,left,right) : LineIntegrationDiscriminant();
    const LineIntegrationDiscriminant& discriminant = cachedDiscriminant == nullptr ? localDiscriminant : *cachedDiscriminant;
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
        const long double mappedDerivative = mapped ? MappedDerivativeMagnitude(discriminant,parameter) : -1.0L;
        FastGraphSample(coefficient,degree,strip,point,computeDelta,lowerZero,upperZero,samples[sample],surfaceGradient,active,mappedDerivative,coordinates,derivativeSign);
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

  static bool GraphGaussPairPanel(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, long double left, long double right, long double parameterLeft, long double parameterRight, bool mapped, bool computeDelta, bool lowerZero, bool upperZero, const std::vector<unsigned char>& active, unsigned maximumEvaluations, unsigned& evaluations, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const LineIntegrationDiscriminant* cachedDiscriminant = nullptr, const std::array<long double,5>* coordinates = nullptr, int derivativeSign = 0) {
    const unsigned size = active.size();
    values.assign(size,0.0L);
    errors.assign(size,std::numeric_limits<long double>::infinity());
    std::vector<long double> low(size,0.0L);
    std::vector<long double> absolute(size,0.0L);
    std::vector<long double> lowPoints;
    std::vector<long double> lowWeights;
    std::vector<long double> highPoints;
    std::vector<long double> highWeights;
    GaussRule(16,lowPoints,lowWeights);
    GaussRule(32,highPoints,highWeights);
    const long double parameterWidth = parameterRight-parameterLeft;
    const long double pi = std::acos(-1.0L);
    const LineIntegrationDiscriminant localDiscriminant = mapped && cachedDiscriminant == nullptr ? MappedDiscriminant(coefficient,left,right) : LineIntegrationDiscriminant();
    const LineIntegrationDiscriminant& discriminant = cachedDiscriminant == nullptr ? localDiscriminant : *cachedDiscriminant;
    try {
      for(unsigned pass = 0; pass < 2; ++pass) {
        const std::vector<long double>& points = pass == 0 ? lowPoints : highPoints;
        const std::vector<long double>& weights = pass == 0 ? lowWeights : highWeights;
        std::vector<long double>& result = pass == 0 ? low : values;
        for(unsigned q = 0; q < points.size(); ++q) {
          if(evaluations >= maximumEvaluations) throw std::runtime_error("LINE mapped graph evaluation budget exhausted.");
          ++evaluations;
          ++Statistics().fastSamples;
          const long double parameter = parameterLeft+parameterWidth*points[q];
          const long double point = MappedPoint(left,right,parameter,mapped);
          const long double jacobian = mapped ? 0.5L*pi*(right-left)*std::sin(pi*parameter) : right-left;
          const long double weight = parameterWidth*weights[q]*jacobian;
          std::vector<long double> sample;
          const long double mappedDerivative = mapped ? MappedDerivativeMagnitude(discriminant,parameter) : -1.0L;
          FastGraphSample(coefficient,degree,strip,point,computeDelta,lowerZero,upperZero,sample,surfaceGradient,active,mappedDerivative,coordinates,derivativeSign);
          for(unsigned k = 0; k < size; ++k) {
            if(!active[k]) continue;
            const long double term = weight*sample[k];
            result[k] += term;
            if(pass == 1) absolute[k] += std::abs(term);
          }
        }
      }
    }
    catch(const std::runtime_error&) {
      std::fill(errors.begin(),errors.end(),std::numeric_limits<long double>::infinity());
      return false;
    }
    for(unsigned k = 0; k < size; ++k) {
      if(!active[k]) continue;
      errors[k] = std::max(std::abs(values[k]-low[k]),256.0L*std::numeric_limits<long double>::epsilon()*absolute[k]);
    }
    return true;
  }

  static void Adaptive(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, long double left, long double right, long double parameterLeft, long double parameterRight, bool mapped, bool graph, bool computeDelta, bool lowerZero, bool upperZero, long double absoluteH, long double absoluteD, long double relative, unsigned depth, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, std::vector<long double>& total, std::vector<long double>& totalError, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>& requested, const LineIntegrationDiscriminant* cachedDiscriminant = nullptr, const std::array<long double,5>* coordinates = nullptr, int derivativeSign = 0) {
    const LineIntegrationDiscriminant localDiscriminant = mapped && cachedDiscriminant == nullptr ? MappedDiscriminant(coefficient,left,right) : LineIntegrationDiscriminant();
    if(cachedDiscriminant == nullptr) cachedDiscriminant = &localDiscriminant;
    std::vector<unsigned char> active(requested);
    std::vector<long double> values;
    std::vector<long double> errors;
    const long double parameterFraction = parameterRight-parameterLeft;
    const long double fractionH = parameterFraction;
    const long double fractionD = parameterFraction;
    if(graph) {
      if(mapped) ++Statistics().mappedPanels;
      else ++Statistics().regularGraphPanels;
    }
    const bool sampled = GaussKronrodPanel(coefficient,degree,strip,left,right,parameterLeft,parameterRight,mapped,graph,computeDelta,lowerZero,upperZero,active,maximumEvaluations,evaluations,values,errors,surfaceGradient,cachedDiscriminant,coordinates,derivativeSign);
    if(!Accept(values,errors,absoluteH*fractionH,absoluteD*fractionD,relative,active,total,totalError)) {
      ++accepted;
      return;
    }
    if(graph && sampled && (coordinates == nullptr || mapped || degree > 5)) {
      ++Statistics().pairPanels;
      const bool paired = GraphGaussPairPanel(coefficient,degree,strip,left,right,parameterLeft,parameterRight,mapped,computeDelta,lowerZero,upperZero,active,maximumEvaluations,evaluations,values,errors,surfaceGradient,cachedDiscriminant,coordinates,derivativeSign);
      if(paired) {
        const unsigned pairCount = active.size()/2U;
        for(unsigned k = 0; k < active.size(); ++k) {
          if(!active[k]) continue;
          const long double tolerance = (k < pairCount ? absoluteH*fractionH : absoluteD*fractionD)+relative*std::abs(values[k]);
          if(tolerance > 0.0L && std::isfinite(errors[k])) {
            const long double ratio = errors[k]/tolerance;
            if(ratio > Statistics().worstPairRatio) {
              Statistics().worstPairRatio = ratio;
              Statistics().worstPairIndex = k;
              Statistics().worstPairValue = values[k];
              Statistics().worstPairError = errors[k];
              Statistics().worstPairTolerance = tolerance;
              Statistics().worstPairLeft = MappedPoint(left,right,parameterLeft,mapped);
              Statistics().worstPairRight = MappedPoint(left,right,parameterRight,mapped);
            }
          }
        }
        if(!Accept(values,errors,absoluteH*fractionH,absoluteD*fractionD,relative,active,total,totalError)) {
          ++accepted;
          return;
        }
      }
    }
    const long double middle = 0.5L*(parameterLeft+parameterRight);
    const unsigned fastDepthLimit = maximumDepth;
    if(graph && sampled && depth < fastDepthLimit && middle != parameterLeft && middle != parameterRight) {
      Adaptive(coefficient,degree,strip,left,right,parameterLeft,middle,mapped,true,computeDelta,lowerZero,upperZero,absoluteH,absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError,surfaceGradient,active,cachedDiscriminant,coordinates,derivativeSign);
      Adaptive(coefficient,degree,strip,left,right,middle,parameterRight,mapped,true,computeDelta,lowerZero,upperZero,absoluteH,absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError,surfaceGradient,active,cachedDiscriminant,coordinates,derivativeSign);
      return;
    }
    if(depth >= maximumDepth || middle == parameterLeft || middle == parameterRight || evaluations >= maximumEvaluations) throw std::runtime_error("LINE domain integration did not reach the requested tolerance.");
    Adaptive(coefficient,degree,strip,left,right,parameterLeft,middle,mapped,graph,computeDelta,lowerZero,upperZero,absoluteH,absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError,surfaceGradient,active,cachedDiscriminant,coordinates,derivativeSign);
    Adaptive(coefficient,degree,strip,left,right,middle,parameterRight,mapped,graph,computeDelta,lowerZero,upperZero,absoluteH,absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError,surfaceGradient,active,cachedDiscriminant,coordinates,derivativeSign);
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

  static void FullStripMoments(unsigned degree, const LineIntegrationStrip& strip, long double left, long double right, std::vector<long double>& values, std::vector<long double>& errors, const std::array<long double,5>* coordinates = nullptr, const std::vector<unsigned char>* requested = nullptr) {
    std::vector<long double> points;
    std::vector<long double> weights;
    GaussRule((degree+3U)/2U,points,weights);
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    for(unsigned q = 0; q < points.size(); ++q) {
      const long double x = left+(right-left)*points[q];
      const long double lower = StripLower(strip,x);
      const long double upper = StripUpper(strip,x);
      const long double parentX = coordinates == nullptr ? x : std::fma((*coordinates)[1],x,(*coordinates)[0]);
      const long double parentLower = coordinates == nullptr ? lower : std::fma((*coordinates)[3],lower,(*coordinates)[2]);
      const long double parentUpper = coordinates == nullptr ? upper : std::fma((*coordinates)[3],upper,(*coordinates)[2]);
      const long double weight = (right-left)*weights[q]*(coordinates == nullptr ? 1.0L : (*coordinates)[4]);
      long double xPower = 1.0L;
      std::vector<long double> fiber(degree+1U);
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

  static void IntegrateStrip(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, bool graph, bool mapped, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr, const std::array<long double,5>* coordinates = nullptr, const std::vector<unsigned char>* requested = nullptr) {
    if(!(strip.right > strip.left)) throw std::invalid_argument("A LINE integration interval must have positive width.");
    std::vector<long double> events;
    std::vector<long double> tangencies;
    CollectEvents(coefficient,strip,events,&tangencies,graph && !mapped);
    const long double intervalBudget = events.size() > 1 ? 1.0L/static_cast<long double>(events.size()-1U) : 1.0L;
    const bool lowerZero = strip.certifiedLowerZero || BoundaryIsZero(coefficient,strip.lowerConstant,strip.lowerSlope);
    const bool upperZero = strip.certifiedUpperZero || BoundaryIsZero(coefficient,strip.upperConstant,strip.upperSlope);
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    std::vector<unsigned char> active = requested == nullptr ? std::vector<unsigned char>(2U*count,1) : *requested;
    if(!computeDelta) std::fill(active.begin()+count,active.end(),0);
    for(unsigned interval = 1; interval < events.size(); ++interval) {
      const long double left = events[interval-1U];
      const long double right = events[interval];
      if(!(right > left)) continue;
      const bool leftTangency = std::find(tangencies.begin(),tangencies.end(),left) != tangencies.end();
      const bool rightTangency = std::find(tangencies.begin(),tangencies.end(),right) != tangencies.end();
      const bool intervalMapped = mapped || leftTangency || rightTangency;
      const long double middle = left+0.5L*(right-left);
      const long double lower = StripLower(strip,middle);
      const long double upper = StripUpper(strip,middle);
      const std::vector<long double> polynomial = FiberPolynomial(coefficient,middle);
      const std::vector<long double> roots = LineSliceUtilities::QuadraticRoots(polynomial[2],polynomial[1],polynomial[0]);
      unsigned number = 0;
      long double middleRoot = 0.0L;
      for(unsigned r = 0; r < roots.size(); ++r) if(roots[r] >= lower && roots[r] <= upper) { middleRoot = roots[r]; ++number; }
      if(graph && (lowerZero || upperZero)) { number = 1; middleRoot = lowerZero ? lower : upper; }
      const long double middleDerivative = graph && !intervalMapped && number == 1 ? std::fma(2.0L*polynomial[2],middleRoot,polynomial[1]) : 0.0L;
      const int derivativeSign = middleDerivative > 0.0L ? 1 : (middleDerivative < 0.0L ? -1 : 0);
      if(number == 0 && !lowerZero && !upperZero) {
        if(Evaluate(coefficient,middle,0.5L*(lower+upper)) > 0.0L) FullStripMoments(degree,strip,left,right,values,errors,coordinates,&active);
        ++accepted;
        continue;
      }
      const LineIntegrationDiscriminant discriminant = intervalMapped ? MappedDiscriminant(coefficient,left,right,leftTangency,rightTangency) : LineIntegrationDiscriminant();
      Adaptive(coefficient,degree,strip,left,right,0.0L,1.0L,intervalMapped,graph && number == 1,computeDelta,lowerZero,upperZero,intervalBudget*absoluteH,intervalBudget*absoluteD,relative,0,maximumDepth,maximumEvaluations,evaluations,accepted,values,errors,surfaceGradient,active,&discriminant,coordinates,derivativeSign);
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

  static long double FiberDiscriminantValue(const std::vector<std::vector<long double>>& coefficient, long double x) {
    const std::vector<long double> polynomial = FiberPolynomial(coefficient,x);
    return std::fma(-4.0L*polynomial[2],polynomial[0],polynomial[1]*polynomial[1]);
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

  static long double RefineDiscriminantEvent(const std::vector<std::vector<long double>>& coefficient, long double root, long double left, long double right) {
    long double current = root;
    for(unsigned iteration = 0; iteration < 12; ++iteration) {
      const long double value = FiberDiscriminantValue(coefficient,current);
      const long double derivative = FiberDiscriminantDerivative(coefficient,current);
      if(derivative == 0.0L || !std::isfinite(value) || !std::isfinite(derivative)) break;
      const long double next = current-value/derivative;
      if(!(next >= left && next <= right) || !std::isfinite(next)) break;
      if(next == current) break;
      current = next;
    }
    return current;
  }

  static void CollectEvents(const std::vector<std::vector<long double>>& coefficient, const LineIntegrationStrip& strip, std::vector<long double>& events, std::vector<long double>* tangencies = nullptr, bool certifiedGraph = false) {
    events = {strip.left,strip.right};
    std::vector<long double> lower = BoundaryPolynomial(coefficient,strip.lowerConstant,strip.lowerSlope);
    std::vector<long double> upper = BoundaryPolynomial(coefficient,strip.upperConstant,strip.upperSlope);
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
    std::vector<long double> discriminant = LineSliceUtilities::PolynomialSum(LineSliceUtilities::PolynomialProduct(b,b),LineSliceUtilities::PolynomialProduct(a,c),-4.0L);
    std::vector<long double> discriminantEvents;
    if(!certifiedGraph) LineSliceUtilities::PolynomialRoots(discriminant,strip.left,strip.right,discriminantEvents);
    for(unsigned k = 0; k < discriminantEvents.size(); ++k) {
      const long double event = RefineDiscriminantEvent(coefficient,discriminantEvents[k],strip.left,strip.right);
      discriminantEvents[k] = event;
      events.push_back(event);
      if(tangencies != nullptr) tangencies->push_back(event);
    }
    for(unsigned side = 0; side < 2; ++side) for(unsigned k = 0; k < boundaryEvents[side].size(); ++k) {
      long double event = boundaryEvents[side][k];
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
