from pathlib import Path
import re

root = Path(__file__).resolve().parents[2]
folder = root / 'src/quadrature/cutfem_ls'

def function_span(text, name):
    match = re.search(r'^  static [^\n]*\b' + re.escape(name) + r'\(', text, re.M)
    if match is None:
        raise RuntimeError('Missing function ' + name)
    start = match.start()
    opening = text.index('{', match.start())
    depth = 1
    end = opening + 1
    while depth:
        if text[end] == '{':
            depth += 1
        elif text[end] == '}':
            depth -= 1
        end += 1
    return start, end

def extract(text, name):
    start, end = function_span(text, name)
    return text[start:end]

def replace_function(text, name, new):
    start, end = function_span(text, name)
    return text[:start] + new + text[end:]

def once(text, old, new):
    if old not in text:
        raise RuntimeError('Missing patch anchor: ' + old[:100])
    return text.replace(old, new, 1)

tri = (folder / 'TriangleIntegration.hpp').read_text()
fe = (folder / 'FELevelSetMoments.hpp').read_text()
fe = once(fe, 'SeriesProductIntegral(xPower, positiveSlice[j], values[index], errors[index], converged);', 'SeriesProductIntegral(xPower, positiveSlice[j], values[index], errors[index], converged);\n          if(!converged) errors[index] = std::numeric_limits<long double>::infinity();')
fe = once(fe, 'SeriesProductIntegral(xPower, deltaSlice[j], values[count + index], errors[count + index], converged);', 'SeriesProductIntegral(xPower, deltaSlice[j], values[count + index], errors[count + index], converged);\n          if(!converged) errors[count + index] = std::numeric_limits<long double>::infinity();')
(folder / 'FELevelSetMoments.hpp').write_text(fe)

triangle_moment = r'''  static long double TriangleMoment(const std::vector<TriangleDomainVertex>& vertices, unsigned firstDegree, unsigned secondDegree, long double* estimatedError = nullptr) {
    const unsigned degree = firstDegree+secondDegree;
    std::vector<long double> coefficient(1,1.0L);
    std::vector<long double> magnitude(1,1.0L);
    for(unsigned n = 1; n <= degree; ++n) {
      std::vector<long double> next((n+1U)*(n+2U)/2U,0.0L);
      std::vector<long double> bound(next.size(),0.0L);
      for(unsigned a = 0; a <= n; ++a) {
        for(unsigned b = 0; b <= n-a; ++b) {
          const unsigned c = n-a-b;
          const unsigned target = a*(2U*n+3U-a)/2U+b;
          const unsigned power[3] = {a,b,c};
          for(unsigned v = 0; v < 3; ++v) {
            if(power[v] == 0) continue;
            const unsigned previousA = a-(v == 0 ? 1U : 0U);
            const unsigned previousB = b-(v == 1 ? 1U : 0U);
            const unsigned source = previousA*(2U*n+1U-previousA)/2U+previousB;
            const long double coordinate = n <= firstDegree ? vertices[v].x : vertices[v].y;
            const long double multiplier = coordinate*static_cast<long double>(power[v])/n;
            next[target] += multiplier*coefficient[source];
            bound[target] += std::abs(multiplier)*magnitude[source];
          }
        }
      }
      coefficient.swap(next);
      magnitude.swap(bound);
    }
    const long double determinant = std::abs(std::fma(vertices[1].x-vertices[0].x,vertices[2].y-vertices[0].y,-(vertices[2].x-vertices[0].x)*(vertices[1].y-vertices[0].y)));
    const long double factor = determinant/((degree+1.0L)*(degree+2.0L));
    long double value = 0.0L;
    long double correction = 0.0L;
    long double absolute = 0.0L;
    for(unsigned k = 0; k < coefficient.size(); ++k) {
      const long double term = coefficient[k]-correction;
      const long double next = value+term;
      correction = (next-value)-term;
      value = next;
      absolute += magnitude[k];
    }
    if(estimatedError != nullptr) *estimatedError = (32.0L+16.0L*degree+4.0L*coefficient.size())*std::numeric_limits<long double>::epsilon()*factor*absolute;
    return factor*value;
  }'''
tri = replace_function(tri, 'TriangleMoment', triangle_moment)
tri = once(tri, 'values[index] = TriangleMoment(vertices,i,j);\n        errors[index] = 32.0L*std::numeric_limits<long double>::epsilon()*std::abs(values[index]);', 'values[index] = TriangleMoment(vertices,i,j,&errors[index]);')
tri = once(tri, '      selectedDirection = 0;\n      return;', '      selectedDirection = 0;\n      CheckErrors(values,errors,degree,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance);\n      return;')
tri = once(tri, 'if(std::abs(determinant) > 256.0L*std::numeric_limits<long double>::epsilon()) triangles.push_back(triangle);', 'if(determinant != 0.0L) triangles.push_back(triangle);')
tri = once(tri, '    point.y = first.y+parameter*(second.y-first.y);\n    return point;', '    point.y = first.y+parameter*(second.y-first.y);\n    if(axis == 0) point.x = coordinate;\n    else point.y = coordinate;\n    return point;')

# Move domain-independent sample/Taylor mathematics into LINE; retain only geometry in TRI.
helpers = []
for name in ['Power', 'BoundaryPolynomial', 'Evaluate', 'StripLower', 'StripUpper', 'FiberPolynomial', 'CollectEvents']:
    helpers.append(extract(tri, name).replace('TriangleDomainStrip', 'LineIntegrationStrip'))
taylor = extract(tri, 'GraphTaylorInterval').replace('TriangleDomainStrip', 'LineIntegrationStrip')
body_start = taylor.index('    const unsigned order = 32;')
body_end = taylor.index('    std::vector<unsigned char> failed(size,0);')
taylor = '  static void TaylorPanel(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const LineIntegrationStrip& strip, long double left, long double right, bool computeDelta, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>* active) {\n' + taylor[body_start:body_end] + '  }'
taylor = taylor.replace('std::vector<long double> values(size,0.0L);', 'values.assign(size,0.0L);').replace('std::vector<long double> errors(size,0.0L);', 'errors.assign(size,0.0L);')
# A Taylor fiber touching a whole strip boundary has a half-delta ownership weight.
taylor = taylor.replace('deltaWeight = SeriesQuotient(jacobian,gy);', 'deltaWeight = SeriesQuotient(jacobian,gy);\n      const bool boundaryRoot = BoundaryIsZero(coefficient,strip.lowerConstant,strip.lowerSlope) || BoundaryIsZero(coefficient,strip.upperConstant,strip.upperSlope);\n      if(boundaryRoot) for(unsigned k = 0; k <= order; ++k) deltaWeight[k] *= 0.5L;')

shared = r'''#ifndef FEMUS_LINE_DOMAIN_INTEGRATION_HPP
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
        Sample(coefficient,degree,strip,point,computeDelta,graph,lowerZero,upperZero,samples[sample],surfaceGradient,active);
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
    if(graph && sampled && depth < std::min(maximumDepth,8U) && middle != parameterLeft && middle != parameterRight) {
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

  static void IntegrateCell(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned order, bool triangle, bool computeDelta, long double derivativeBound, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
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
'''
shared += '\n\n'.join(helpers) + '\n\n' + taylor + '\n};\n\n}\n\n#endif\n'
(folder / 'LineDomainIntegration.hpp').write_text(shared)
line = (folder / 'LineIntegration.hpp').read_text()
line = once(line, '#include "FELevelSetMoments.hpp"', '#include "FELevelSetMoments.hpp"\n#include "LineDomainIntegration.hpp"')
line = once(line, 'class LineIntegration {', 'class LineIntegration : public LineDomainIntegration {')
(folder / 'LineIntegration.hpp').write_text(line)

# TRI now contributes domain geometry only, never a second adaptive integrator.
start = tri.index('struct TriangleDomainStrip {')
end = tri.index('};', start)+2
tri = tri[:start] + tri[end:]
tri = tri.replace('TriangleDomainStrip', 'LineIntegrationStrip')
start, unused = function_span(tri, 'BoundaryPolynomial')
end, unused = function_span(tri, 'IntegrateDirection')
tri = tri[:start]+tri[end:]
old_direction = extract(tri, 'IntegrateDirection')
opening = old_direction.index('    for(unsigned stripIndex = 0;')
ending = old_direction.index('    if(direction == 1)')
new_direction = old_direction[:opening] + '    for(unsigned stripIndex = 0; stripIndex < strips.size(); ++stripIndex) LineIntegration::IntegrateStrip(workingCoefficient,degree,strips[stripIndex],requireGraph,mappedGraph,useTaylor,computeDelta,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance,maximumDepth,maximumEvaluations,evaluations,accepted,raw,rawError,gradient);\n' + old_direction[ending:]
tri = replace_function(tri, 'IntegrateDirection', new_direction)
(folder / 'TriangleIntegration.hpp').write_text(tri)

quad = (folder / 'Quad9DirectionSelection.hpp').read_text()
start = quad.index('  if(!useLine) {', quad.index('inline void Quad9SelectedOneDirection('))
end = quad.index('  if(useLine) {', start)
quad = quad[:start] + '  if(!useLine) LineIntegration::IntegrateCell(coefficient,degree,order,false,computeDelta,derivativeBound,absoluteTolerance,absoluteTolerance*scale,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,surfaceGradient);\n' + quad[end:]
# A second direction uses the independent, non-split path rather than recomputing the same split tree.
needle = 'inline void Quad9SelectedReferenceOneDirection('
start = quad.index(needle)
end = quad.index(' {', start)
quad = quad[:end] + quad[end:]
signature = quad[start:end]
quad = quad[:start] + signature[:-1] + ', bool allowSplit = true)' + quad[end:]
quad = quad.replace('if(derivativeBound == 0.0L && !fullBoundaryDelta) Quad9SelectedSplitCanonical(', 'if(allowSplit && derivativeBound == 0.0L && !fullBoundaryDelta) Quad9SelectedSplitCanonical(')
quad = quad.replace('const unsigned quadrantOrder = Quad9SelectedOrder(normalized,lower,upper);', 'const unsigned automaticQuadrantOrder = Quad9SelectedOrder(normalized,lower,upper);\n        const unsigned quadrantOrder = allowSplit ? automaticQuadrantOrder : order;')
needle = 'trial[candidate],trialError[candidate],evaluations,trialAccepted[candidate],useTaylor,surfaceGradient);'
quad = once(quad, needle, 'trial[candidate],trialError[candidate],evaluations,trialAccepted[candidate],candidate == 1 && checkBoth ? false : useTaylor,surfaceGradient,!(checkBoth && candidate == 1));')
(folder / 'Quad9DirectionSelection.hpp').write_text(quad)

selected = (folder / 'Tri7DirectionSelection.hpp').read_text()
start = selected.index('  if(!line) {', selected.index('inline void Tri7SelectedOneDirection('))
end = selected.index('  if(line) {', start)
selected = selected[:start] + '  if(!line) LineIntegration::IntegrateCell(coefficient,degree,order,true,computeDelta,derivativeBound,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,raw,rawError,evaluations,accepted,surfaceGradient);\n' + selected[end:]
selected = selected.replace('integrationDirection == 0 && validationMode != 2 && !fullBoundaryDelta', 'integrationDirection == 0 && !fullBoundaryDelta')
needle = '      TriangleIntegration::Integrate(normalized,degree,vertices,computeDelta,chartAbsoluteTolerance,chartAbsoluteTolerance*scale,chartRelativeTolerance,maximumDepth,maximumEvaluations,0,10,values,errors,evaluations,accepted,selected,referenceSurface ? &referenceGradient : nullptr);'
verification = r'''
      if(validationMode != 0) {
        bool verified = false;
        std::string failure = "TRI split and independent direction did not agree.";
        for(unsigned candidate = 0; candidate < 3 && !verified; ++candidate) {
          std::vector<long double> independent;
          std::vector<long double> independentError;
          unsigned independentAccepted = 0;
          unsigned independentEvaluations = 0;
          ++attempts;
          try {
            Tri7SelectedOneDirection(normalized,scale,degree,ranked[candidate],lower[ranked[candidate]-1U],computeDelta,fullBoundaryDelta,chartAbsoluteTolerance,chartRelativeTolerance,maximumDepth,maximumEvaluations,zeroBoundary,false,independent,independentError,independentEvaluations,independentAccepted,referenceSurface ? &referenceGradient : nullptr);
            verified = true;
            for(unsigned k = 0; k < values.size(); ++k) {
              const long double estimate = std::max(std::abs(values[k]-independent[k]),std::max(errors[k],independentError[k]));
              const long double tolerance = chartAbsoluteTolerance*(k < count ? 1.0L : scale)+chartRelativeTolerance*std::abs(values[k]);
              if(!std::isfinite(estimate) || estimate > tolerance) verified = false;
            }
            if(verified) {
              for(unsigned k = 0; k < values.size(); ++k) errors[k] = std::max(std::abs(values[k]-independent[k]),std::max(errors[k],independentError[k]));
              compared = ranked[candidate];
              evaluations += independentEvaluations;
            }
          }
          catch(const std::runtime_error& exception) { failure = exception.what(); }
        }
        if(!verified) throw std::runtime_error(failure);
      }'''
selected = once(selected, needle, needle + verification)
(folder / 'Tri7DirectionSelection.hpp').write_text(selected)

for name in ['cutfem-rebuild.yml', 'quad-tri-benchmark.yml']:
    path = root / '.github/workflows' / name
    text = path.read_text().replace('ref: rifat-cutfem-rebuild', 'ref: ${{ github.sha }}')
    if name == 'quad-tri-benchmark.yml':
        text = text.replace("      - '.github/workflows/quad-tri-benchmark.yml'", "      - '.github/workflows/quad-tri-benchmark.yml'\n      - 'src/quadrature/cutfem_ls/**'")
        text = text.replace('run: python3 candidate/unittests/testCutFemWeightLSNative/compare_quad_tri_benchmark.py baseline_a.csv candidate_a.csv candidate_b.csv baseline_b.csv | tee benchmark_summary.txt', 'run: |\n          set -o pipefail\n          python3 candidate/unittests/testCutFemWeightLSNative/compare_quad_tri_benchmark.py baseline_a.csv candidate_a.csv candidate_b.csv baseline_b.csv | tee benchmark_summary.txt')
    path.write_text(text)
print('Review source repairs applied to the pinned source tree.')
