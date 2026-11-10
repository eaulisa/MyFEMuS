#ifndef FEMUS_OUTER_INTERVAL_QUADRATURE_HPP
#define FEMUS_OUTER_INTERVAL_QUADRATURE_HPP

#include "FELevelSetMoments.hpp"
#include <array>
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace FELevelSetMoments {

struct OuterIntervalStatistics {
  unsigned long long regularPanels = 0;
  unsigned long long transformedPanels = 0;
  unsigned long long emergency7 = 0;
  unsigned long long emergency3 = 0;
  unsigned long long subdivisions = 0;
  unsigned long long coordinateCollapses = 0;
  unsigned long long resolutionLimits = 0;
  unsigned long long regularIntervals = 0;
  unsigned long long oneSingularIntervals = 0;
  unsigned long long twoSingularIntervals = 0;
  unsigned long long higherOrderMaps = 0;
  unsigned long long samples = 0;
  unsigned long long precisionSamples = 0;
  unsigned long long acceptedMoments = 0;
  unsigned long long skippedMoments = 0;
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
      if(!(middle > left && middle < right)) throw std::runtime_error("An outer midpoint is below coordinate resolution.");
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

  template<unsigned Rule> static bool Resolved(const OuterInterval& interval, long double left, long double right) {
    if(interval.mapping == 0) return true;
    const long double parameter = left+0.5L*(right-left)-0.5L*(right-left)*boost::math::quadrature::gauss_kronrod<long double,Rule>::abscissa().back();
    const long double coordinate = Point(interval,parameter).coordinate;
    return coordinate != interval.endpoint;
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

  template<unsigned Rule, class Context> static void Panel(Context& context, void (*sample)(Context&, const OuterIntervalPoint&, const std::vector<unsigned char>&, std::vector<long double>&, std::vector<long double>&, long double, long double, long double), const OuterInterval& interval, long double left, long double right, const std::vector<unsigned char>& active, long double absoluteH, long double absoluteD, long double relative, unsigned maximumEvaluations, unsigned& evaluations, unsigned dimension, bool precisionCheck, std::vector<long double>& values, std::vector<long double>& errors) {
    const long double middle = left+0.5L*(right-left);
    const long double half = 0.5L*(right-left);
    const long double physicalWidth = interval.mapping == 0 ? right-left : interval.delta*(ReferencePower(right,interval.exponent)-ReferencePower(left,interval.exponent));
    values.assign(active.size(),0);
    errors.assign(active.size(),0);
    std::vector<long double> low(active.size(),0);
    std::vector<long double> norm(active.size(),0);
    for(unsigned q = 0; q < boost::math::quadrature::gauss_kronrod<long double,Rule>::abscissa().size(); ++q) for(unsigned side = 0; side < (q == 0 ? 1U : 2U); ++side) {
      const long double parameter = middle+(side == 0 ? -half : half)*boost::math::quadrature::gauss_kronrod<long double,Rule>::abscissa()[q];
      const OuterIntervalPoint point = Point(interval,parameter);
      if(interval.mapping != 0 && point.coordinate == interval.endpoint) throw std::runtime_error("A transformed outer point collapsed onto its endpoint.");
      if(evaluations >= maximumEvaluations) throw std::runtime_error("Outer quadrature evaluation budget exhausted.");
      std::vector<long double> value;
      std::vector<long double> error;
      ++Statistics(dimension).samples;
      sample(context,point,active,value,error,absoluteH/physicalWidth,absoluteD/physicalWidth,relative);
      if(precisionCheck) {
        OuterIntervalPoint adjacent = point;
        adjacent.coordinate = std::nextafter(point.coordinate,interval.endpoint);
        if(adjacent.coordinate == interval.endpoint) adjacent.coordinate = std::nextafter(point.coordinate,interval.mapping > 0 ? std::numeric_limits<long double>::infinity() : -std::numeric_limits<long double>::infinity());
        adjacent.offset = adjacent.coordinate-interval.endpoint;
        if(adjacent.coordinate == interval.endpoint) throw std::runtime_error("An adjacent outer point collapsed onto its endpoint.");
        adjacent.parameter = std::pow(std::abs(adjacent.offset)/interval.delta,1.0L/interval.exponent);
        if(evaluations >= maximumEvaluations) throw std::runtime_error("Outer quadrature evaluation budget exhausted.");
        std::vector<long double> adjacentValue;
        std::vector<long double> adjacentError;
        ++Statistics(dimension).samples;
        ++Statistics(dimension).precisionSamples;
        sample(context,adjacent,active,adjacentValue,adjacentError,absoluteH/physicalWidth,absoluteD/physicalWidth,relative);
        for(unsigned n = 0; n < active.size(); ++n) if(active[n]) error[n] += std::abs(value[n]-adjacentValue[n])+adjacentError[n];
      }
      const long double weight = half*boost::math::quadrature::gauss_kronrod<long double,Rule>::weights()[q];
      const long double gaussWeight = q%2U == 0 ? half*boost::math::quadrature::gauss<long double,(Rule-1U)/2U>::weights()[q/2U] : 0;
      for(unsigned n = 0; n < active.size(); ++n) {
        if(!active[n]) { ++Statistics(dimension).skippedMoments; continue; }
        values[n] += weight*value[n];
        low[n] += gaussWeight*value[n];
        errors[n] += weight*error[n];
        norm[n] += weight*std::abs(value[n]);
      }
    }
    for(unsigned n = 0; n < active.size(); ++n) if(active[n]) errors[n] += std::max(std::abs(values[n]-low[n]),256*std::numeric_limits<long double>::epsilon()*norm[n]);
  }

  template<class Context> static void Adaptive(Context& context, void (*sample)(Context&, const OuterIntervalPoint&, const std::vector<unsigned char>&, std::vector<long double>&, std::vector<long double>&, long double, long double, long double), const OuterInterval& interval, long double left, long double right, long double absoluteH, long double absoluteD, long double relative, unsigned depth, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, unsigned dimension, const std::vector<unsigned char>& requested, std::vector<long double>& total, std::vector<long double>& totalErrors) {
    std::vector<unsigned char> active = requested;
    std::vector<long double> values;
    std::vector<long double> errors;
    if(interval.mapping != 0 && !Resolved<15>(interval,left,right)) {
      ++Statistics(dimension).coordinateCollapses;
      if(Resolved<7>(interval,left,right)) {
        ++Statistics(dimension).emergency7;
        Panel<7>(context,sample,interval,left,right,active,absoluteH,absoluteD,relative,maximumEvaluations,evaluations,dimension,true,values,errors);
        if(!Accept(values,errors,absoluteH,absoluteD,relative,active,total,totalErrors,dimension)) { ++accepted; return; }
      }
      if(Resolved<3>(interval,left,right)) {
        ++Statistics(dimension).emergency3;
        Panel<3>(context,sample,interval,left,right,active,absoluteH,absoluteD,relative,maximumEvaluations,evaluations,dimension,true,values,errors);
        if(!Accept(values,errors,absoluteH,absoluteD,relative,active,total,totalErrors,dimension)) { ++accepted; return; }
      }
      ++Statistics(dimension).resolutionLimits;
      throw std::runtime_error("Outer tolerance cannot be satisfied at the available coordinate precision.");
    }
    if(interval.mapping == 0) ++Statistics(dimension).regularPanels; else ++Statistics(dimension).transformedPanels;
    Panel<15>(context,sample,interval,left,right,active,absoluteH,absoluteD,relative,maximumEvaluations,evaluations,dimension,false,values,errors);
    if(!Accept(values,errors,absoluteH,absoluteD,relative,active,total,totalErrors,dimension)) { ++accepted; return; }
    const long double middle = left+0.5L*(right-left);
    if(depth >= maximumDepth || middle == left || middle == right) {
      throw std::runtime_error("Outer adaptive quadrature did not reach tolerance.");
    }
    ++Statistics(dimension).subdivisions;
    Adaptive(context,sample,interval,left,middle,0.5L*absoluteH,0.5L*absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,dimension,active,total,totalErrors);
    Adaptive(context,sample,interval,middle,right,0.5L*absoluteH,0.5L*absoluteD,relative,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,dimension,active,total,totalErrors);
  }
};

}

#endif
