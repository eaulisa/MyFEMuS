#ifndef FEMUS_OUTER_INTERVAL_QUADRATURE_HPP
#define FEMUS_OUTER_INTERVAL_QUADRATURE_HPP

#include "FELevelSetMoments.hpp"
#include <array>
#include "GaussKronrodRules.hpp"

namespace FELevelSetMoments {

class OuterQuadratureResolutionError : public std::runtime_error {
public:
  explicit OuterQuadratureResolutionError(const std::string& message) : std::runtime_error(message) {}
};

struct OuterIntervalStatistics {
  unsigned long long regularPanels = 0;
  unsigned long long transformedPanels = 0;
  unsigned long long emergency15 = 0;
  unsigned long long emergency7 = 0;
  unsigned long long emergency3 = 0;
  unsigned long long subdivisions = 0;
  unsigned long long coordinateCollapses = 0;
  unsigned long long resolutionLimits = 0;
  unsigned long long childResolutionFailures = 0;
  unsigned long long regularIntervals = 0;
  unsigned long long oneSingularIntervals = 0;
  unsigned long long twoSingularIntervals = 0;
  unsigned long long higherOrderMaps = 0;
  unsigned long long samples = 0;
  unsigned long long precisionSamples = 0;
  unsigned long long acceptedMoments = 0;
  unsigned long long skippedMoments = 0;
  std::array<unsigned long long,15> regularRulePanels = {};
  std::array<unsigned long long,15> transformedRulePanels = {};
  std::array<unsigned long long,15> emergencyRulePanels = {};
};

struct OuterInterval {
  long double left = 0;
  long double right = 1;
  long double endpoint = 0;
  long double delta = 0;
  int mapping = 0;
  unsigned exponent = 1;
  bool exactEndpoint = true;
};

struct OuterIntervalPoint {
  long double parameter = 0;
  long double coordinate = 0;
  long double offset = 0;
  long double jacobian = 1;
};

class OuterIntervalQuadrature {
public:
  static OuterIntervalStatistics& Statistics(unsigned dimension) {
    static thread_local std::array<OuterIntervalStatistics,2> data;
    return data[dimension-2U];
  }

  static std::vector<OuterInterval> Plan(long double left, long double right, unsigned leftExponent, unsigned rightExponent, unsigned dimension, bool leftExact = true, bool rightExact = true) {
    if(!(right > left) || leftExponent == 0 || rightExponent == 0) throw std::invalid_argument("Invalid outer integration interval.");
    std::vector<OuterInterval> result;
    const bool leftSingular = leftExponent > 1;
    const bool rightSingular = rightExponent > 1;
    if(!leftSingular && !rightSingular) {
      ++Statistics(dimension).regularIntervals;
      result.push_back({left,right,0,0,0,1,true});
    }
    else if(leftSingular && rightSingular) {
      ++Statistics(dimension).twoSingularIntervals;
      const long double middle = left+0.5L*(right-left);
      if(!(middle > left && middle < right)) throw OuterQuadratureResolutionError("An outer midpoint is below coordinate resolution.");
      result.push_back({left,middle,left,middle-left,1,leftExponent,leftExact});
      result.push_back({middle,right,right,right-middle,-1,rightExponent,rightExact});
    }
    else {
      ++Statistics(dimension).oneSingularIntervals;
      result.push_back({left,right,leftSingular ? left : right,right-left,leftSingular ? 1 : -1,leftSingular ? leftExponent : rightExponent,leftSingular ? leftExact : rightExact});
    }
    for(unsigned i = 0; i < result.size(); ++i) if(result[i].exponent > 2) ++Statistics(dimension).higherOrderMaps;
    return result;
  }

  static OuterIntervalPoint Point(const OuterInterval& interval, long double parameter) {
    OuterIntervalPoint result;
    result.parameter = parameter;
    result.coordinate = parameter;
    if(interval.mapping != 0) {
      const long double power = ReferencePower(parameter,interval.exponent);
      result.offset = interval.mapping*interval.delta*power;
      result.coordinate = std::fma(interval.mapping*interval.delta,power,interval.endpoint);
      result.jacobian = interval.exponent*interval.delta*ReferencePower(parameter,interval.exponent-1U);
    }
    return result;
  }

  static bool Resolved(const GaussKronrodRule& rule, const OuterInterval& interval, long double left, long double right) {
    if(interval.mapping == 0) return true;
    const long double parameter = left+0.5L*(right-left)-0.5L*(right-left)*rule.nodes[rule.size-1U];
    return Point(interval,parameter).coordinate != interval.endpoint;
  }

  template<unsigned Rule> static bool Resolved(const OuterInterval& interval, long double left, long double right) {
    return Resolved(GaussKronrodRules::Get(Rule),interval,left,right);
  }

  static bool Accept(const std::vector<long double>& values, const std::vector<long double>& errors, long double absoluteH, long double absoluteD, long double relative, std::vector<unsigned char>& active, std::vector<long double>& total, std::vector<long double>& totalErrors, unsigned dimension) {
    bool remaining = false;
    const unsigned count = active.size()/2U;
    for(unsigned n = 0; n < active.size(); ++n) {
      if(!active[n]) continue;
      const long double tolerance = (n < count ? absoluteH : absoluteD)+relative*std::abs(values[n]);
      if(std::isfinite(values[n]) && std::isfinite(errors[n]) && errors[n] >= 0 && errors[n] <= tolerance) {
        total[n] += values[n];
        totalErrors[n] += errors[n];
        active[n] = 0;
        ++Statistics(dimension).acceptedMoments;
      }
      else remaining = true;
    }
    return remaining;
  }

  template<class Context> static bool Panel(const GaussKronrodRule& rule, Context& context, void (*sample)(Context&, const OuterIntervalPoint&, const std::vector<unsigned char>&, std::vector<long double>&, std::vector<long double>&, long double, long double, long double), const OuterInterval& interval, long double left, long double right, const std::vector<unsigned char>& active, long double absoluteH, long double absoluteD, long double relative, unsigned maximumEvaluations, unsigned& evaluations, unsigned dimension, bool precisionCheck, std::vector<long double>& values, std::vector<long double>& errors) {
    const unsigned long long previousWarnings = QuadratureWarningCount();
    const long double middle = left+0.5L*(right-left);
    const long double half = 0.5L*(right-left);
    const long double physicalWidth = interval.mapping == 0 ? right-left : interval.delta*(ReferencePower(right,interval.exponent)-ReferencePower(left,interval.exponent));
    values.assign(active.size(),0);
    errors.assign(active.size(),0);
    std::vector<long double> low(active.size(),0);
    std::vector<long double> norm(active.size(),0);
    for(unsigned q = 0; q < rule.size; ++q) for(unsigned side = 0; side < (q == 0 ? 1U : 2U); ++side) {
      const long double parameter = middle+(side == 0 ? -half : half)*rule.nodes[q];
      const OuterIntervalPoint point = Point(interval,parameter);
      if(interval.mapping != 0 && point.coordinate == interval.endpoint) throw OuterQuadratureResolutionError("A transformed outer point collapsed onto its endpoint.");
      std::vector<long double> value;
      std::vector<long double> error;
      ++Statistics(dimension).samples;
      sample(context,point,active,value,error,absoluteH/physicalWidth,absoluteD/physicalWidth,relative);
      if(precisionCheck) {
        OuterIntervalPoint adjacent = point;
        if(interval.mapping == 0) {
          adjacent.coordinate = std::nextafter(point.coordinate,left);
          if(adjacent.coordinate == left) adjacent.coordinate = std::nextafter(point.coordinate,right);
          adjacent.parameter = adjacent.coordinate;
        }
        else {
          adjacent.coordinate = std::nextafter(point.coordinate,interval.endpoint);
          if(adjacent.coordinate == interval.endpoint) adjacent.coordinate = std::nextafter(point.coordinate,interval.mapping > 0 ? std::numeric_limits<long double>::infinity() : -std::numeric_limits<long double>::infinity());
          adjacent.offset = adjacent.coordinate-interval.endpoint;
          if(adjacent.coordinate == interval.endpoint) throw OuterQuadratureResolutionError("An adjacent outer point collapsed onto its endpoint.");
          adjacent.parameter = std::pow(std::abs(adjacent.offset)/interval.delta,1.0L/interval.exponent);
        }
          std::vector<long double> adjacentValue;
        std::vector<long double> adjacentError;
        ++Statistics(dimension).samples;
        ++Statistics(dimension).precisionSamples;
        sample(context,adjacent,active,adjacentValue,adjacentError,absoluteH/physicalWidth,absoluteD/physicalWidth,relative);
        for(unsigned n = 0; n < active.size(); ++n) if(active[n]) error[n] += std::abs(value[n]-adjacentValue[n])+adjacentError[n];
      }
      const long double weight = half*rule.weights[q];
      const long double gaussWeight = q%2U == (rule.gauss%2U == 0 ? 1U : 0U) ? half*rule.gaussWeights[q/2U] : 0;
      for(unsigned n = 0; n < active.size(); ++n) {
        if(!active[n]) { ++Statistics(dimension).skippedMoments; continue; }
        values[n] += weight*value[n];
        low[n] += gaussWeight*value[n];
        errors[n] += weight*error[n];
        norm[n] += weight*std::abs(value[n]);
      }
    }
    for(unsigned n = 0; n < active.size(); ++n) if(active[n]) errors[n] += std::max(std::abs(values[n]-low[n]),256*std::numeric_limits<long double>::epsilon()*norm[n]);
    return QuadratureWarningCount() == previousWarnings;
  }

  template<class Context> static void Adaptive(Context& context, void (*sample)(Context&, const OuterIntervalPoint&, const std::vector<unsigned char>&, std::vector<long double>&, std::vector<long double>&, long double, long double, long double), const OuterInterval& interval, long double left, long double right, long double absoluteH, long double absoluteD, long double relative, unsigned depth, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, unsigned dimension, const std::vector<unsigned char>& requested, std::vector<long double>& total, std::vector<long double>& totalErrors, const OuterQuadratureRules& rules = OuterQuadratureRules()) {
    QuadratureWarningScope warningScope(dimension == 2 ? "2D outer quadrature" : "3D outer quadrature");
    std::vector<unsigned char> active = requested;
    std::vector<long double> values, errors;
    std::vector<long double> best(active.size(),std::numeric_limits<long double>::quiet_NaN());
    std::vector<long double> bestErrors(active.size(),std::numeric_limits<long double>::infinity());
    const GaussKronrodRule& selected = GaussKronrodRules::Get(rules.Select(dimension,interval.mapping != 0));
    const bool budgetLimited = evaluations >= maximumEvaluations || maximumEvaluations-evaluations < selected.kronrod;
    const GaussKronrodRule& primary = budgetLimited ? GaussKronrodRules::Get(3) : selected;
    bool reduce = interval.mapping != 0 && !Resolved(primary,interval,left,right);
    if(!reduce) {
      if(interval.mapping == 0) { ++Statistics(dimension).regularPanels; ++Statistics(dimension).regularRulePanels[primary.index]; }
      else { ++Statistics(dimension).transformedPanels; ++Statistics(dimension).transformedRulePanels[primary.index]; }
      try { const bool childAccuracy = Panel(primary,context,sample,interval,left,right,active,absoluteH,absoluteD,relative,maximumEvaluations,evaluations,dimension,false,values,errors); KeepBestQuadrature(values,errors,active,best,bestErrors); if(!childAccuracy) { ++Statistics(dimension).childResolutionFailures; reduce = true; } }
      catch(const QuadratureAccuracyRetry&) { ++Statistics(dimension).childResolutionFailures; reduce = true; }
      catch(const OuterQuadratureResolutionError&) { ++Statistics(dimension).childResolutionFailures; reduce = true; }
    }
    if(reduce) {
      ++Statistics(dimension).coordinateCollapses;
      for(unsigned lower = primary.gauss/2U; lower > 0; lower /= 2U) {
        const GaussKronrodRule& coarse = GaussKronrodRules::Get(2U*lower+1U);
        if(!Resolved(coarse,interval,left,right)) continue;
        ++Statistics(dimension).emergencyRulePanels[coarse.index];
        if(coarse.kronrod == 15) ++Statistics(dimension).emergency15;
        else if(coarse.kronrod == 7) ++Statistics(dimension).emergency7;
        else if(coarse.kronrod == 3) ++Statistics(dimension).emergency3;
        try { const bool childAccuracy = Panel(coarse,context,sample,interval,left,right,active,absoluteH,absoluteD,relative,maximumEvaluations,evaluations,dimension,true,values,errors); KeepBestQuadrature(values,errors,active,best,bestErrors); if(!childAccuracy) { ++Statistics(dimension).childResolutionFailures; continue; } }
        catch(const QuadratureAccuracyRetry&) { ++Statistics(dimension).childResolutionFailures; continue; }
        catch(const OuterQuadratureResolutionError&) { ++Statistics(dimension).childResolutionFailures; continue; }
        if(!Accept(values,errors,absoluteH,absoluteD,relative,active,total,totalErrors,dimension)) { ++accepted; return; }
      }
      ++Statistics(dimension).resolutionLimits;
      for(unsigned n = 0; n < active.size(); ++n) if(active[n] && !std::isfinite(best[n])) throw OuterQuadratureResolutionError("No complete finite outer rule is available at coordinate resolution.");
      RetainUncertifiedQuadrature("Outer quadrature","child quadrature accuracy or coordinate precision prevents the requested tolerance",best,bestErrors,absoluteH,absoluteD,relative,active,total,totalErrors);
      ++accepted;
      return;
    }
    if(!Accept(values,errors,absoluteH,absoluteD,relative,active,total,totalErrors,dimension)) { ++accepted; return; }
    std::vector<unsigned char> precisionLimited(active.size(),0);
    bool pending = false;
    for(unsigned n = 0; n < active.size(); ++n) if(active[n]) {
      if(bestErrors[n] <= 512*std::numeric_limits<long double>::epsilon()*std::abs(best[n])) { precisionLimited[n] = 1; active[n] = 0; }
      else pending = true;
    }
    RetainUncertifiedQuadrature("Outer quadrature","floating-point error floor prevents requested tolerance",best,bestErrors,absoluteH,absoluteD,relative,precisionLimited,total,totalErrors);
    if(!pending) { ++accepted; return; }
    const long double middle = left+0.5L*(right-left);
    if(budgetLimited || evaluations >= maximumEvaluations || depth >= maximumDepth || middle == left || middle == right) {
      const std::string problem = budgetLimited || evaluations >= maximumEvaluations ? "evaluation budget exhausted before tolerance" : (depth >= maximumDepth ? "maximum depth reached before tolerance" : "coordinate precision exhausted before tolerance");
      RetainUncertifiedQuadrature("Outer quadrature",problem,best,bestErrors,absoluteH,absoluteD,relative,active,total,totalErrors);
      ++accepted;
      return;
    }
    ++Statistics(dimension).subdivisions;
    const std::vector<long double> savedTotal = total;
    const std::vector<long double> savedErrors = totalErrors;
    const unsigned savedAccepted = accepted;
    try {
      Adaptive(context,sample,interval,left,middle,0.5L*absoluteH,0.5L*absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,dimension,active,total,totalErrors,rules);
      Adaptive(context,sample,interval,middle,right,0.5L*absoluteH,0.5L*absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,dimension,active,total,totalErrors,rules);
    }
    catch(const OuterQuadratureResolutionError&) {
      for(unsigned n = 0; n < active.size(); ++n) if(active[n]) { total[n] = savedTotal[n]; totalErrors[n] = savedErrors[n]; }
      accepted = savedAccepted+1U;
      RetainUncertifiedQuadrature("Outer quadrature","child coordinate resolution exhausted; retained the complete parent rule",best,bestErrors,absoluteH,absoluteD,relative,active,total,totalErrors);
    }
  }

};

}

#endif
