#ifndef FEMUS_QUADRATURE_WARNINGS_HPP
#define FEMUS_QUADRATURE_WARNINGS_HPP

#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace FELevelSetMoments {

class QuadratureAccuracyRetry : public std::runtime_error {
public:
  explicit QuadratureAccuracyRetry(const std::string& problem) : std::runtime_error(problem) {}
};

inline unsigned& QuadratureAccuracyRetryDepth() { static thread_local unsigned depth = 0; return depth; }

class QuadratureAccuracyRetryScope {
public:
  explicit QuadratureAccuracyRetryScope(bool enabled) : _previous(QuadratureAccuracyRetryDepth()) { QuadratureAccuracyRetryDepth() = enabled ? _previous+1U : 0U; }
  ~QuadratureAccuracyRetryScope() { QuadratureAccuracyRetryDepth() = _previous; }
  QuadratureAccuracyRetryScope(const QuadratureAccuracyRetryScope&) = delete;
  QuadratureAccuracyRetryScope& operator=(const QuadratureAccuracyRetryScope&) = delete;
private:
  unsigned _previous;
};

inline unsigned long long& QuadratureWarningCount() { static thread_local unsigned long long count = 0; return count; }
inline unsigned& QuadratureWarningDepth() { static thread_local unsigned depth = 0; return depth; }
inline std::string& QuadratureWarningContext() { static thread_local std::string context; return context; }
inline std::map<std::string,std::vector<long double>>& QuadratureWarningMessages() { static thread_local std::map<std::string,std::vector<long double>> messages; return messages; }

inline void FlushQuadratureWarnings() {
  if(QuadratureWarningMessages().empty()) return;
  std::ostringstream message;
  message << "CutFemWeight_LS warning [" << QuadratureWarningContext() << "]";
  for(const std::pair<const std::string,std::vector<long double>>& entry : QuadratureWarningMessages()) message << "; " << entry.first << " (failed checks=" << entry.second[0] << ", worst moment=" << entry.second[2] << ", estimated error=" << entry.second[3] << ", requested tolerance=" << entry.second[4] << ')';
  message << "; returning the best available quadrature; the reported accuracy checks were not satisfied.\n";
  std::cerr << message.str();
  QuadratureWarningMessages().clear();
}

inline bool QuadratureExceptionActive() {
#if __cplusplus >= 201703L
  return std::uncaught_exceptions() != 0;
#else
  return std::uncaught_exception();
#endif
}

class QuadratureWarningScope {
public:
  explicit QuadratureWarningScope(const std::string& context) {
    if(QuadratureWarningDepth()++ == 0) { QuadratureWarningContext() = context; QuadratureWarningMessages().clear(); }
  }
  ~QuadratureWarningScope() {
    if(--QuadratureWarningDepth() == 0) { if(!QuadratureExceptionActive()) FlushQuadratureWarnings(); else QuadratureWarningMessages().clear(); }
  }
  QuadratureWarningScope(const QuadratureWarningScope&) = delete;
  QuadratureWarningScope& operator=(const QuadratureWarningScope&) = delete;
};

// Aggregates repeated accuracy failures into one warning for each top-level quadrature request.
inline void WarnQuadrature(const std::string& stage, const std::string& problem, unsigned moment, long double error, long double tolerance) {
  if(QuadratureAccuracyRetryDepth() != 0) throw QuadratureAccuracyRetry(stage+": "+problem);
  ++QuadratureWarningCount();
  std::vector<long double>& record = QuadratureWarningMessages()[stage+": "+problem];
  if(record.empty()) record.assign(5,0);
  const long double estimate = std::isfinite(error) && error >= 0 ? error : std::numeric_limits<long double>::infinity();
  const long double ratio = tolerance > 0 ? estimate/tolerance : std::numeric_limits<long double>::infinity();
  if(record[0] == 0 || ratio > record[1]) { record[1] = ratio; record[2] = moment; record[3] = estimate; record[4] = tolerance; }
  ++record[0];
  if(QuadratureWarningDepth() == 0) { QuadratureWarningContext() = stage; FlushQuadratureWarnings(); }
}

inline void CheckQuadratureTolerance(const std::string& stage, const std::string& problem, const std::vector<long double>& values, const std::vector<long double>& errors, long double absoluteH, long double absoluteD, long double relative, const std::vector<unsigned char>* active = nullptr) {
  if(errors.size() != values.size() || (active != nullptr && active->size() != values.size())) throw std::invalid_argument("Invalid quadrature error-check dimensions.");
  const unsigned count = values.size()/2U;
  for(unsigned n = 0; n < values.size(); ++n) if(active == nullptr || (*active)[n]) {
    if(!std::isfinite(values[n])) throw std::runtime_error(stage+": a quadrature moment is not finite.");
    const long double tolerance = (n < count ? absoluteH : absoluteD)+relative*std::abs(values[n]);
    if(!std::isfinite(errors[n]) || errors[n] < 0 || errors[n] > tolerance) WarnQuadrature(stage,problem,n,errors[n],tolerance);
  }
}

inline void KeepBestQuadrature(const std::vector<long double>& values, const std::vector<long double>& errors, const std::vector<unsigned char>& active, std::vector<long double>& best, std::vector<long double>& bestErrors) {
  for(unsigned n = 0; n < active.size(); ++n) if(active[n] && std::isfinite(values[n])) {
    const long double error = std::isfinite(errors[n]) && errors[n] >= 0 ? errors[n] : std::numeric_limits<long double>::infinity();
    if(!std::isfinite(best[n]) || error < bestErrors[n]) { best[n] = values[n]; bestErrors[n] = error; }
  }
}

inline void RetainUncertifiedQuadrature(const std::string& stage, const std::string& problem, const std::vector<long double>& best, const std::vector<long double>& bestErrors, long double absoluteH, long double absoluteD, long double relative, std::vector<unsigned char>& active, std::vector<long double>& total, std::vector<long double>& errors) {
  const unsigned count = active.size()/2U;
  for(unsigned n = 0; n < active.size(); ++n) if(active[n]) {
    if(!std::isfinite(best[n])) throw std::runtime_error(stage+": no complete finite quadrature estimate is available.");
    total[n] += best[n];
    errors[n] += bestErrors[n];
    active[n] = 0;
    WarnQuadrature(stage,problem,n,bestErrors[n],(n < count ? absoluteH : absoluteD)+relative*std::abs(best[n]));
  }
}

}

#endif
