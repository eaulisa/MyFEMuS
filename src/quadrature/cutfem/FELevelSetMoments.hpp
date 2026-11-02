#ifndef FE_LEVEL_SET_MOMENTS_HPP
#define FE_LEVEL_SET_MOMENTS_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

// Native moments use [0,1] for LINE3, [0,1]^2 for QUAD9, and 0<=y<=x<=1 for TRI7; legacy entry points retain FEMuS-coordinate monomials.

namespace FELevelSetMoments {

inline long double Power(long double x, unsigned n) {
  long double value = 1.0L;
  while (n != 0) {
    if ((n & 1U) != 0) {
      value *= x;
    }
    x *= x;
    n >>= 1U;
  }
  return value;
}

inline long double Binomial(unsigned n, unsigned k) {
  if (k > n) {
    return 0.0L;
  }
  k = std::min(k, n - k);
  long double value = 1.0L;
  for (unsigned j = 1; j <= k; ++j) {
    value *= static_cast<long double>(n - k + j) / j;
  }
  return value;
}

inline long double Heaviside(long double x) {
  return (x > 0.0L) ? 1.0L : ((x < 0.0L) ? 0.0L : 0.5L);
}

inline long double MonomialIntegral(long double left, long double right, unsigned degree) {
  long double sum = 0.0L;
  long double correction = 0.0L;
  for (unsigned k = 0; k <= degree; ++k) {
    long double term = Power(left, k) * Power(right, degree - k) - correction;
    long double next = sum + term;
    correction = (next - sum) - term;
    sum = next;
  }
  return (right - left) * sum / static_cast<long double>(degree + 1U);
}

inline void Trim(std::vector<long double>& polynomial) {
  while (polynomial.size() > 1 && polynomial.back() == 0.0L) {
    polynomial.pop_back();
  }
}

inline long double PolynomialValue(const std::vector<long double>& polynomial, long double x) {
  long double value = 0.0L;
  for (unsigned k = static_cast<unsigned>(polynomial.size()); k > 0; --k) {
    value = std::fma(value, x, polynomial[k - 1U]);
  }
  return value;
}

inline long double PolynomialBound(const std::vector<long double>& polynomial, long double x) {
  long double value = 0.0L;
  for (unsigned k = static_cast<unsigned>(polynomial.size()); k > 0; --k) {
    value = value * std::abs(x) + std::abs(polynomial[k - 1U]);
  }
  return value;
}

inline long double PolynomialNorm(const std::vector<long double>& polynomial) {
  long double value = 0.0L;
  for (unsigned k = 0; k < polynomial.size(); ++k) {
    value = std::max(value, std::abs(polynomial[k]));
  }
  return value;
}

inline std::vector<long double> PolynomialDerivative(const std::vector<long double>& polynomial) {
  if (polynomial.size() <= 1) {
    return std::vector<long double>(1, 0.0L);
  }
  std::vector<long double> result(polynomial.size() - 1U, 0.0L);
  for (unsigned k = 1; k < polynomial.size(); ++k) {
    result[k - 1U] = k * polynomial[k];
  }
  Trim(result);
  return result;
}

inline std::vector<long double> PolynomialProduct(const std::vector<long double>& first, const std::vector<long double>& second) {
  std::vector<long double> result(first.size() + second.size() - 1U, 0.0L);
  for (unsigned i = 0; i < first.size(); ++i) {
    for (unsigned j = 0; j < second.size(); ++j) {
      result[i + j] += first[i] * second[j];
    }
  }
  Trim(result);
  return result;
}

inline std::vector<long double> PolynomialSum(const std::vector<long double>& first, const std::vector<long double>& second, long double multiplier = 1.0L) {
  std::vector<long double> result(std::max(first.size(), second.size()), 0.0L);
  for (unsigned k = 0; k < first.size(); ++k) {
    result[k] += first[k];
  }
  for (unsigned k = 0; k < second.size(); ++k) {
    result[k] += multiplier * second[k];
  }
  Trim(result);
  return result;
}

inline std::vector<long double> DivideLinear(const std::vector<long double>& polynomial, long double root) {
  if (polynomial.size() <= 1) {
    throw std::runtime_error("Invalid polynomial deflation.");
  }
  std::vector<long double> result(polynomial.size() - 1U, 0.0L);
  result.back() = polynomial.back();
  for (unsigned k = static_cast<unsigned>(result.size()) - 1U; k > 0; --k) {
    result[k - 1U] = polynomial[k] + root * result[k];
  }
  Trim(result);
  return result;
}

inline void SortDistinct(std::vector<long double>& values) {
  std::sort(values.begin(), values.end());
  unsigned count = 0;
  for (unsigned k = 0; k < values.size(); ++k) {
    long double tolerance = 16.0L * std::numeric_limits<long double>::epsilon() * std::max(1.0L, std::abs(values[k]));
    if (count == 0 || values[k] - values[count - 1U] > tolerance) {
      values[count++] = values[k];
    }
  }
  values.resize(count);
}

inline std::vector<long double> QuadraticRoots(long double a, long double b, long double c) {
  std::vector<long double> result;
  if (a == 0.0L) {
    if (b != 0.0L) {
      result.push_back(-c / b);
    }
    return result;
  }
  long double discriminant = std::fma(-4.0L * a, c, b * b);
  if (discriminant < 0.0L) {
    return result;
  }
  if (discriminant == 0.0L) {
    result.push_back(-b / (2.0L * a));
    return result;
  }
  long double squareRoot = std::sqrt(discriminant);
  long double numerator = -0.5L * (b + std::copysign(squareRoot, b));
  result.push_back(numerator / a);
  result.push_back(c / numerator);
  std::sort(result.begin(), result.end());
  return result;
}

inline void PolynomialRoots(const std::vector<long double>& input, long double left, long double right, std::vector<long double>& roots) {
  std::vector<long double> polynomial(input);
  Trim(polynomial);
  if (polynomial.size() <= 1) {
    return;
  }
  if (polynomial.size() <= 3) {
    std::vector<long double> candidates = QuadraticRoots((polynomial.size() == 3) ? polynomial[2] : 0.0L, polynomial[1], polynomial[0]);
    for (unsigned k = 0; k < candidates.size(); ++k) {
      if (candidates[k] >= left && candidates[k] <= right) {
        roots.push_back(candidates[k]);
      }
    }
    return;
  }
  std::vector<long double> derivative = PolynomialDerivative(polynomial);
  std::vector<long double> critical;
  PolynomialRoots(derivative, left, right, critical);
  critical.push_back(left);
  critical.push_back(right);
  SortDistinct(critical);
  for (unsigned k = 0; k < critical.size(); ++k) {
    long double value = PolynomialValue(polynomial, critical[k]);
    long double tolerance = 64.0L * std::numeric_limits<long double>::epsilon() * PolynomialBound(polynomial, critical[k]);
    if (std::abs(value) <= tolerance) {
      roots.push_back(critical[k]);
    }
  }
  for (unsigned k = 1; k < critical.size(); ++k) {
    long double a = critical[k - 1U];
    long double b = critical[k];
    long double fa = PolynomialValue(polynomial, a);
    long double fb = PolynomialValue(polynomial, b);
    if (fa == 0.0L || fb == 0.0L || (fa > 0.0L) == (fb > 0.0L)) {
      continue;
    }
    for (unsigned iteration = 0; iteration < 180; ++iteration) {
      long double middle = a + 0.5L * (b - a);
      if (middle == a || middle == b) {
        break;
      }
      long double fm = PolynomialValue(polynomial, middle);
      if (fm == 0.0L) {
        a = middle;
        b = middle;
        break;
      }
      if ((fm > 0.0L) == (fa > 0.0L)) {
        a = middle;
        fa = fm;
      }
      else {
        b = middle;
        fb = fm;
      }
    }
    roots.push_back(a + 0.5L * (b - a));
  }
}

inline void Integrate1DInternal(long double a, long double b, long double c, unsigned degree, long double left, long double right, bool computeDelta, bool fullBoundaryDelta, std::vector<long double>& positive, std::vector<long double>& delta) {
  positive.assign(degree + 1U, 0.0L);
  delta.assign(degree + 1U, 0.0L);
  if (!(right > left)) {
    throw std::invalid_argument("The 1D interval must have positive length.");
  }
  if (a == 0.0L && b == 0.0L) {
    if (c == 0.0L && computeDelta) {
      throw std::domain_error("The delta of an identically zero field is undefined.");
    }
    for (unsigned j = 0; j <= degree; ++j) {
      positive[j] = Heaviside(c) * MonomialIntegral(left, right, j);
    }
    return;
  }
  std::vector<long double> roots = QuadraticRoots(a, b, c);
  long double discriminant = std::fma(-4.0L * a, c, b * b);
  if (a != 0.0L && discriminant == 0.0L && !roots.empty() && roots[0] >= left && roots[0] <= right && computeDelta) {
    throw std::domain_error("A repeated 1D zero has no regular delta moment.");
  }
  std::vector<long double> cuts(1, left);
  for (unsigned k = 0; k < roots.size(); ++k) {
    if (roots[k] > left && roots[k] < right) {
      cuts.push_back(roots[k]);
    }
  }
  cuts.push_back(right);
  SortDistinct(cuts);
  for (unsigned k = 1; k < cuts.size(); ++k) {
    long double middle = cuts[k - 1U] + 0.5L * (cuts[k] - cuts[k - 1U]);
    if (std::fma(a, middle * middle, std::fma(b, middle, c)) > 0.0L) {
      for (unsigned j = 0; j <= degree; ++j) {
        positive[j] += MonomialIntegral(cuts[k - 1U], cuts[k], j);
      }
    }
  }
  if (!computeDelta) {
    return;
  }
  for (unsigned k = 0; k < roots.size(); ++k) {
    long double root = roots[k];
    long double endpointTolerance = 32.0L * std::numeric_limits<long double>::epsilon() * std::max(1.0L, std::max(std::abs(left), std::abs(right)));
    if (std::fma(a, left * left, std::fma(b, left, c)) == 0.0L && std::abs(root - left) <= endpointTolerance) {
      root = left;
    }
    if (std::fma(a, right * right, std::fma(b, right, c)) == 0.0L && std::abs(root - right) <= endpointTolerance) {
      root = right;
    }
    if (root < left || root > right) {
      continue;
    }
    long double derivative = std::abs(2.0L * a * root + b);
    if (derivative == 0.0L) {
      throw std::domain_error("A repeated 1D zero has no regular delta moment.");
    }
    long double weight = (!fullBoundaryDelta && (root == left || root == right)) ? 0.5L : 1.0L;
    long double power = 1.0L;
    for (unsigned j = 0; j <= degree; ++j) {
      delta[j] += weight * power / derivative;
      power *= root;
    }
  }
}

inline void Integrate1D(double a, double b, double c, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, double left = 0.0, double right = 1.0, bool computeDelta = true, bool fullBoundaryDelta = false) {
  if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(left) || !std::isfinite(right) || degree == std::numeric_limits<unsigned>::max()) {
    throw std::invalid_argument("Invalid 1D input.");
  }
  long double scale = std::max(std::abs(static_cast<long double>(a)), std::max(std::abs(static_cast<long double>(b)), std::abs(static_cast<long double>(c))));
  if (scale == 0.0L) {
    scale = 1.0L;
  }
  std::vector<long double> positiveInternal;
  std::vector<long double> deltaInternal;
  Integrate1DInternal(a / scale, b / scale, c / scale, degree, left, right, computeDelta, fullBoundaryDelta, positiveInternal, deltaInternal);
  positive.resize(degree + 1U);
  delta.resize(degree + 1U);
  for (unsigned j = 0; j <= degree; ++j) {
    positive[j] = static_cast<double>(positiveInternal[j]);
    delta[j] = static_cast<double>(deltaInternal[j] / scale);
    if (!std::isfinite(positive[j]) || !std::isfinite(delta[j])) {
      throw std::overflow_error("A 1D moment exceeds double precision range.");
    }
  }
}

inline std::vector<long double> SeriesConstant(long double value, unsigned order) {
  std::vector<long double> result(order + 1U, 0.0L);
  result[0] = value;
  return result;
}

inline std::vector<long double> SeriesProduct(const std::vector<long double>& first, const std::vector<long double>& second) {
  std::vector<long double> result(first.size(), 0.0L);
  unsigned firstDegree = static_cast<unsigned>(first.size())-1U;
  unsigned secondDegree = static_cast<unsigned>(second.size())-1U;
  while(firstDegree != 0 && first[firstDegree] == 0.0L) --firstDegree;
  while(secondDegree != 0 && second[secondDegree] == 0.0L) --secondDegree;
  unsigned degree = std::min(static_cast<unsigned>(result.size())-1U,firstDegree+secondDegree);
  for (unsigned k = 0; k <= degree; ++k) {
    unsigned begin = k > secondDegree ? k-secondDegree : 0;
    unsigned end = std::min(k,firstDegree);
    for (unsigned j = begin; j <= end; ++j) {
      result[k] += first[j] * second[k - j];
    }
  }
  return result;
}

inline std::vector<long double> SeriesQuotient(const std::vector<long double>& numerator, const std::vector<long double>& denominator) {
  if (denominator[0] == 0.0L) {
    throw std::runtime_error("A Taylor chart has a zero denominator.");
  }
  std::vector<long double> result(numerator.size(), 0.0L);
  for (unsigned k = 0; k < result.size(); ++k) {
    long double value = numerator[k];
    for (unsigned j = 1; j <= k; ++j) {
      value -= denominator[j] * result[k - j];
    }
    result[k] = value / denominator[0];
  }
  return result;
}

inline std::vector<long double> SeriesSquareRoot(const std::vector<long double>& input) {
  if (!(input[0] > 0.0L)) {
    throw std::runtime_error("A Taylor square-root chart is not positive.");
  }
  std::vector<long double> result(input.size(), 0.0L);
  result[0] = std::sqrt(input[0]);
  for (unsigned k = 1; k < result.size(); ++k) {
    long double value = input[k];
    for (unsigned j = 1; j < k; ++j) {
      value -= result[j] * result[k - j];
    }
    result[k] = value / (2.0L * result[0]);
  }
  return result;
}

inline std::vector<long double> SeriesPolynomial(const std::vector<long double>& polynomial, const std::vector<long double>& argument) {
  std::vector<long double> result(argument.size(), 0.0L);
  if(polynomial.empty()) return result;
  result[0] = polynomial.back();
  for (unsigned k = static_cast<unsigned>(polynomial.size())-1U; k > 0; --k) {
    result = SeriesProduct(result, argument);
    result[0] += polynomial[k - 1U];
  }
  return result;
}

inline std::vector<long double> SeriesPower(const std::vector<long double>& argument, unsigned exponent) {
  std::vector<long double> result = SeriesConstant(1.0L, static_cast<unsigned>(argument.size()) - 1U);
  std::vector<long double> factor(argument);
  while (exponent != 0) {
    if ((exponent & 1U) != 0) {
      result = SeriesProduct(result, factor);
    }
    exponent >>= 1U;
    if (exponent != 0) {
      factor = SeriesProduct(factor, factor);
    }
  }
  return result;
}

inline void SineCosineSeries(long double middle, long double halfWidth, unsigned order, std::vector<long double>& sine, std::vector<long double>& cosine) {
  sine.assign(order + 1U, 0.0L);
  cosine.assign(order + 1U, 0.0L);
  sine[0] = std::sin(middle);
  cosine[0] = std::cos(middle);
  for (unsigned k = 1; k <= order; ++k) {
    sine[k] = halfWidth * cosine[k - 1U] / k;
    cosine[k] = -halfWidth * sine[k - 1U] / k;
  }
}

inline void SeriesIntegral(const std::vector<long double>& coefficients, long double& value, long double& error, bool& converged) {
  value = 0.0L;
  long double correction = 0.0L;
  long double norm = 0.0L;
  long double tail = 0.0L;
  long double preceding = 0.0L;
  unsigned order = static_cast<unsigned>(coefficients.size()) - 1U;
  unsigned block = std::max(8U, order / 6U);
  for (unsigned k = 0; k <= order; k += 2U) {
    long double term = 2.0L * coefficients[k] / (k + 1U);
    long double corrected = term - correction;
    long double next = value + corrected;
    correction = (next - value) - corrected;
    value = next;
    norm += std::abs(term);
    if (k + block > order) {
      tail += std::abs(term);
    }
    else if (k + 2U * block > order) {
      preceding += std::abs(term);
    }
  }
  long double rounding = 64.0L * std::numeric_limits<long double>::epsilon() * norm;
  error = 8.0L * tail + rounding;
  converged = std::isfinite(value) && std::isfinite(error) && (tail <= 0.5L * preceding || tail <= rounding);
}

inline void SeriesProductIntegral(const std::vector<long double>& first, const std::vector<long double>& second, long double& value, long double& error, bool& converged) {
  value = 0.0L;
  long double correction = 0.0L;
  long double norm = 0.0L;
  long double tail = 0.0L;
  long double preceding = 0.0L;
  unsigned order = static_cast<unsigned>(first.size())-1U;
  unsigned block = std::max(8U,order/6U);
  unsigned firstDegree = order;
  unsigned secondDegree = order;
  while(firstDegree != 0 && first[firstDegree] == 0.0L) --firstDegree;
  while(secondDegree != 0 && second[secondDegree] == 0.0L) --secondDegree;
  for(unsigned k = 0; k <= order; k += 2U) {
    long double coefficient = 0.0L;
    unsigned begin = k > secondDegree ? k-secondDegree : 0;
    unsigned end = std::min(k,firstDegree);
    for(unsigned j = begin; j <= end; ++j) coefficient += first[j]*second[k-j];
    long double term = 2.0L*coefficient/(k+1U);
    long double corrected = term-correction;
    long double next = value+corrected;
    correction = (next-value)-corrected;
    value = next;
    norm += std::abs(term);
    if(k+block > order) tail += std::abs(term);
    else if(k+2U*block > order) preceding += std::abs(term);
  }
  long double rounding = 64.0L*std::numeric_limits<long double>::epsilon()*norm;
  error = 8.0L*tail+rounding;
  converged = std::isfinite(value) && std::isfinite(error) && (tail <= 0.5L*preceding || tail <= rounding);
}

inline unsigned MomentIndex(unsigned degree, unsigned i, unsigned j) {
  return i * (2U * degree + 3U - i) / 2U + j;
}

inline long double LowerBoundary(long double x, unsigned triangle, const std::vector<long double>& domain) {
  return (triangle == 3U) ? x : domain[2];
}

inline long double UpperBoundary(long double x, unsigned triangle, const std::vector<long double>& domain) {
  return (triangle == 1U) ? 1.0L - x : ((triangle == 2U) ? x : domain[3]);
}

inline long double FullStripMoment(long double left, long double right, unsigned i, unsigned j, unsigned triangle, const std::vector<long double>& domain) {
  if (!triangle) {
    return MonomialIntegral(left, right, i) * MonomialIntegral(domain[2], domain[3], j);
  }
  if (triangle == 2U) {
    return MonomialIntegral(left, right, i + j + 1U) / static_cast<long double>(j + 1U);
  }
  if (triangle == 3U) {
    return (MonomialIntegral(left, right, i) - MonomialIntegral(left, right, i + j + 1U)) / static_cast<long double>(j + 1U);
  }
  long double width = right - left;
  long double result = 0.0L;
  for (unsigned p = 0; p <= i; ++p) {
    for (unsigned q = 0; q <= j + 1U; ++q) {
      long double term = Binomial(i, p) * Binomial(j + 1U, q) * Power(left, i - p) * Power(1.0L - right, j + 1U - q) * Power(width, p + q);
      term /= static_cast<long double>(p + q + 1U) * Binomial(p + q, p);
      result += term;
    }
  }
  return width * result / static_cast<long double>(j + 1U);
}

inline std::vector<long double> BoundaryPolynomial(const std::vector<std::vector<long double>>& field, bool top, unsigned triangle, const std::vector<long double>& domain) {
  std::vector<long double> boundary(1, top ? domain[3] : domain[2]);
  if (top && triangle == 1U) {
    boundary.assign(2, 0.0L);
    boundary[0] = 1.0L;
    boundary[1] = -1.0L;
  }
  if ((top && triangle == 2U) || (!top && triangle == 3U)) {
    boundary.assign(2, 0.0L);
    boundary[1] = 1.0L;
  }
  std::vector<long double> result = PolynomialSum(field[2], PolynomialProduct(field[1], boundary));
  result = PolynomialSum(result, PolynomialProduct(field[0], PolynomialProduct(boundary, boundary)));
  Trim(result);
  return result;
}

inline std::vector<long double> DiscriminantPolynomial(const std::vector<std::vector<long double>>& field) {
  return PolynomialSum(PolynomialProduct(field[1], field[1]), PolynomialProduct(field[0], field[2]), -4.0L);
}

inline bool NearZero(const std::vector<long double>& polynomial, long double x) {
  long double scale = PolynomialBound(polynomial, x);
  return std::abs(PolynomialValue(polynomial, x)) <= 128.0L * std::numeric_limits<long double>::epsilon() * scale;
}

inline void CheckSingularZeros(const std::vector<std::vector<long double>>& field, const std::vector<long double>& events, unsigned triangle, const std::vector<long double>& domain) {
  std::vector<long double> da = PolynomialDerivative(field[0]);
  std::vector<long double> db = PolynomialDerivative(field[1]);
  std::vector<long double> dc = PolynomialDerivative(field[2]);
  for (unsigned k = 0; k < events.size(); ++k) {
    long double x = events[k];
    long double low = LowerBoundary(x, triangle, domain);
    long double high = UpperBoundary(x, triangle, domain);
    long double a = PolynomialValue(field[0], x);
    long double b = PolynomialValue(field[1], x);
    long double c = PolynomialValue(field[2], x);
    if (NearZero(field[0], x) && NearZero(field[1], x) && NearZero(field[2], x)) {
      long double dx2 = PolynomialValue(da, x);
      long double dx1 = PolynomialValue(db, x);
      long double dx0 = PolynomialValue(dc, x);
      std::vector<long double> roots = QuadraticRoots(dx2, dx1, dx0);
      if (dx2 == 0.0L && dx1 == 0.0L && dx0 == 0.0L && high >= low) {
        throw std::domain_error("Singular or numerically unresolved zero of the level set; request positive moments only to omit delta moments.");
      }
      for (unsigned r = 0; r < roots.size(); ++r) {
        if (roots[r] >= low && roots[r] <= high) {
          throw std::domain_error("The level set has a singular zero in the closed reference cell.");
        }
      }
      continue;
    }
    if (a == 0.0L) {
      continue;
    }
    long double y = -b / (2.0L * a);
    if (y < low || y > high) {
      continue;
    }
    long double value = std::fma(a, y * y, std::fma(b, y, c));
    long double bound = std::abs(a * y * y) + std::abs(b * y) + std::abs(c);
    long double gx = std::fma(PolynomialValue(da, x), y * y, std::fma(PolynomialValue(db, x), y, PolynomialValue(dc, x)));
    long double gxBound = PolynomialBound(da, x) * y * y + PolynomialBound(db, x) * std::abs(y) + PolynomialBound(dc, x);
    long double tolerance = 256.0L * std::numeric_limits<long double>::epsilon();
    if (std::abs(value) <= tolerance * bound && std::abs(gx) <= tolerance * gxBound) {
      throw std::domain_error("Singular or numerically unresolved zero of the level set; request positive moments only to omit delta moments.");
    }
  }
}

inline void IntervalBranches(const std::vector<std::vector<long double>>& field, long double x, unsigned triangle, const std::vector<long double>& domain, bool bottomZero, bool topZero, bool computeDelta, std::vector<int>& signs, std::vector<int>& boundary, std::vector<long double>& positiveSegments, long double& uniformPhase) {
  signs.clear();
  boundary.clear();
  positiveSegments.clear();
  long double a = PolynomialValue(field[0], x);
  long double b = PolynomialValue(field[1], x);
  long double c = PolynomialValue(field[2], x);
  long double low = LowerBoundary(x, triangle, domain);
  long double high = UpperBoundary(x, triangle, domain);
  long double discriminant = std::fma(-4.0L * a, c, b * b);
  std::vector<long double> roots = QuadraticRoots(a, b, c);
  if (a != 0.0L && discriminant == 0.0L) {
    if (!roots.empty() && roots[0] >= low && roots[0] <= high && computeDelta) {
      throw std::domain_error("The level set contains a repeated zero curve.");
    }
    uniformPhase = Heaviside(a);
    return;
  }
  long double tolerance = 256.0L * std::numeric_limits<long double>::epsilon() * std::max(1.0L, std::max(std::abs(low), std::abs(high)));
  for (unsigned k = 0; k < roots.size(); ++k) {
    long double root = roots[k];
    int tag = 0;
    if (bottomZero && std::abs(root - low) <= tolerance) {
      root = low;
      tag = -1;
    }
    if (topZero && std::abs(root - high) <= tolerance) {
      root = high;
      tag = 1;
    }
    if (root < low || root > high) {
      continue;
    }
    if (root == low && tag == 0) {
      tag = -1;
    }
    if (root == high && tag == 0) {
      tag = 1;
    }
    roots[k] = root;
    signs.push_back((a == 0.0L) ? 0 : ((2.0L * a * root + b > 0.0L) ? 1 : -1));
    boundary.push_back(tag);
  }
  std::vector<long double> cuts(1, low);
  for (unsigned k = 0; k < roots.size(); ++k) {
    if (roots[k] > low && roots[k] < high) {
      cuts.push_back(roots[k]);
    }
  }
  cuts.push_back(high);
  SortDistinct(cuts);
  for (unsigned k = 1; k < cuts.size(); ++k) {
    long double y = cuts[k - 1U] + 0.5L * (cuts[k] - cuts[k - 1U]);
    bool isPositive = false;
    if (a != 0.0L && discriminant > 0.0L) {
      bool between = y > roots.front() && y < roots.back();
      isPositive = between ? (a < 0.0L) : (a > 0.0L);
    }
    else if (b != 0.0L) {
      isPositive = b * (y + c / b) > 0.0L;
    }
    else {
      isPositive = std::fma(a, y * y, c) > 0.0L;
    }
    positiveSegments.push_back(isPositive ? 1.0L : 0.0L);
  }
  long double y = low + 0.5L * (high - low);
  uniformPhase = Heaviside(std::fma(a, y * y, std::fma(b, y, c)));
}

inline void FactorEndpoint(std::vector<long double>& polynomial, long double root, unsigned& multiplicity) {
  multiplicity = 0;
  while (polynomial.size() > 1 && NearZero(polynomial, root)) {
    polynomial = DivideLinear(polynomial, root);
    ++multiplicity;
  }
}

inline bool SliceTaylor(const std::vector<std::vector<long double>>& field, const std::vector<long double>& reducedDiscriminant, unsigned leftMultiplicity, unsigned rightMultiplicity, long double xLeft, long double xRight, long double thetaLeft, long double thetaRight, const std::vector<int>& signs, const std::vector<int>& boundary, const std::vector<long double>& positiveSegments, unsigned degree, unsigned order, unsigned triangle, const std::vector<long double>& domain, bool computePositive, bool computeDelta, bool fullBoundaryDelta, std::vector<long double>& values, std::vector<long double>& errors, bool fullDeltaChart = false) {
  long double middle = 0.5L * (thetaLeft + thetaRight);
  long double halfWidth = 0.5L * (thetaRight - thetaLeft);
  long double width = xRight - xLeft;
  std::vector<long double> sine;
  std::vector<long double> cosine;
  SineCosineSeries(middle, halfWidth, order, sine, cosine);
  std::vector<long double> x = SeriesProduct(sine, sine);
  for (unsigned k = 0; k <= order; ++k) {
    x[k] *= width;
  }
  x[0] += xLeft;
  std::vector<long double> jacobian = SeriesProduct(sine, cosine);
  for (unsigned k = 0; k <= order; ++k) {
    jacobian[k] *= 2.0L * width * halfWidth;
  }
  std::vector<long double> lower = SeriesConstant(domain[2], order);
  std::vector<long double> upper = SeriesConstant(domain[3], order);
  if (triangle == 1U) {
    for (unsigned k = 0; k <= order; ++k) {
      upper[k] = -x[k];
    }
    upper[0] += 1.0L;
  }
  if (triangle == 2U) {
    upper = x;
  }
  if (triangle == 3U) {
    lower = x;
  }
  std::vector<long double> a = SeriesPolynomial(field[0], x);
  std::vector<long double> b = SeriesPolynomial(field[1], x);
  std::vector<long double> c = SeriesPolynomial(field[2], x);
  bool quadratic = PolynomialNorm(field[0]) != 0.0L;
  std::vector<long double> squareRoot(order + 1U, 0.0L);
  std::vector<long double> reducedSquareRoot;
  if (quadratic) {
    reducedSquareRoot = SeriesSquareRoot(SeriesPolynomial(reducedDiscriminant,x));
    squareRoot = reducedSquareRoot;
    squareRoot = SeriesProduct(squareRoot, SeriesPower(sine, leftMultiplicity));
    squareRoot = SeriesProduct(squareRoot, SeriesPower(cosine, rightMultiplicity));
    long double factor = std::pow(width, 0.5L * static_cast<long double>(leftMultiplicity + rightMultiplicity));
    for (unsigned k = 0; k <= order; ++k) {
      squareRoot[k] *= factor;
    }
  }
  std::vector<std::vector<long double>> roots;
  std::vector<std::vector<long double>> weights;
  std::vector<long double> da,db,dc,commonDeltaWeight;
  if(computeDelta && !fullDeltaChart) {
    da = SeriesPolynomial(PolynomialDerivative(field[0]),x);
    db = SeriesPolynomial(PolynomialDerivative(field[1]),x);
    dc = SeriesPolynomial(PolynomialDerivative(field[2]),x);
  }
  if(computeDelta && fullDeltaChart) {
    if(quadratic && leftMultiplicity <= 1 && rightMultiplicity <= 1) {
      commonDeltaWeight = SeriesConstant(2.0L*width*halfWidth/std::pow(width,0.5L*(leftMultiplicity+rightMultiplicity)),order);
      if(leftMultiplicity == 0) commonDeltaWeight = SeriesProduct(commonDeltaWeight,sine);
      if(rightMultiplicity == 0) commonDeltaWeight = SeriesProduct(commonDeltaWeight,cosine);
      commonDeltaWeight = SeriesQuotient(commonDeltaWeight,reducedSquareRoot);
    }
    else {
      std::vector<long double> derivative = quadratic ? squareRoot : b;
      if(!quadratic && derivative[0] < 0.0L) {
        for(unsigned k = 0; k <= order; ++k) derivative[k] = -derivative[k];
      }
      commonDeltaWeight = SeriesQuotient(jacobian,derivative);
    }
  }
  for (unsigned r = 0; r < signs.size(); ++r) {
    std::vector<long double> root(order + 1U, 0.0L);
    if (boundary[r] != 0) {
      root = (boundary[r] < 0) ? lower : upper;
    }
    else if (!quadratic) {
      root = SeriesQuotient(c, b);
      for (unsigned k = 0; k <= order; ++k) {
        root[k] = -root[k];
      }
    }
    else if (static_cast<long double>(signs[r]) * b[0] > 0.0L && std::abs(-b[0] + signs[r] * squareRoot[0]) < 0.25L * (std::abs(b[0]) + squareRoot[0])) {
      std::vector<long double> denominator(b);
      std::vector<long double> numerator(c);
      for (unsigned k = 0; k <= order; ++k) {
        denominator[k] += signs[r] * squareRoot[k];
        numerator[k] *= -2.0L;
      }
      root = SeriesQuotient(numerator, denominator);
    }
    else {
      std::vector<long double> numerator(b);
      std::vector<long double> denominator(a);
      for (unsigned k = 0; k <= order; ++k) {
        numerator[k] = -numerator[k] + signs[r] * squareRoot[k];
        denominator[k] *= 2.0L;
      }
      root = SeriesQuotient(numerator, denominator);
    }
    roots.push_back(root);
    if (computeDelta) {
      if(fullDeltaChart) {
        std::vector<long double> weight(commonDeltaWeight);
        if(!fullBoundaryDelta && boundary[r] != 0) {
          for(unsigned k = 0; k <= order; ++k) weight[k] *= 0.5L;
        }
        weights.push_back(weight);
        continue;
      }
      std::vector<long double> gy = quadratic ? squareRoot : b;
      if (!quadratic && gy[0] < 0.0L) {
        for (unsigned k = 0; k <= order; ++k) {
          gy[k] = -gy[k];
        }
      }
      std::vector<long double> gx = SeriesProduct(da, SeriesProduct(root, root));
      std::vector<long double> term = SeriesProduct(db, root);
      for (unsigned k = 0; k <= order; ++k) {
        gx[k] += term[k] + dc[k];
      }
      std::vector<long double> denominator = SeriesProduct(gx, gx);
      term = SeriesProduct(gy, gy);
      for (unsigned k = 0; k <= order; ++k) {
        denominator[k] += term[k];
      }
      std::vector<long double> weight = SeriesQuotient(SeriesProduct(jacobian, gy), denominator);
      if (!fullBoundaryDelta && boundary[r] != 0) {
        for (unsigned k = 0; k <= order; ++k) {
          weight[k] *= 0.5L;
        }
      }
      weights.push_back(weight);
    }
  }
  std::vector<std::vector<long double>> cuts(1, lower);
  for (unsigned r = 0; r < roots.size(); ++r) {
    if (boundary[r] == 0) {
      cuts.push_back(roots[r]);
    }
  }
  cuts.push_back(upper);
  if (cuts.size() != positiveSegments.size() + 1U) {
    throw std::runtime_error("Unresolved slice topology.");
  }
  std::vector<std::vector<long double>> positivePowers(cuts.size(), SeriesConstant(1.0L, order));
  std::vector<std::vector<long double>> rootPowers(roots.size(), SeriesConstant(1.0L, order));
  std::vector<std::vector<long double>> positiveSlice(degree + 1U, SeriesConstant(0.0L, order));
  std::vector<std::vector<long double>> deltaSlice(degree + 1U, SeriesConstant(0.0L, order));
  for (unsigned j = 0; j <= degree; ++j) {
    if (computePositive) {
      for (unsigned r = 0; r < cuts.size(); ++r) {
        positivePowers[r] = SeriesProduct(positivePowers[r], cuts[r]);
      }
      for (unsigned r = 0; r < positiveSegments.size(); ++r) {
        if (positiveSegments[r] != 0.0L) {
          for (unsigned k = 0; k <= order; ++k) {
            positiveSlice[j][k] += (positivePowers[r + 1U][k] - positivePowers[r][k]) / static_cast<long double>(j + 1U);
          }
        }
      }
      positiveSlice[j] = SeriesProduct(positiveSlice[j], jacobian);
    }
    if (computeDelta) {
      for (unsigned r = 0; r < roots.size(); ++r) {
        std::vector<long double> term = SeriesProduct(rootPowers[r], weights[r]);
        for (unsigned k = 0; k <= order; ++k) {
          deltaSlice[j][k] += term[k];
        }
        rootPowers[r] = SeriesProduct(rootPowers[r], roots[r]);
      }
    }
  }
  unsigned count = (degree + 1U) * (degree + 2U) / 2U;
  values.assign(2U * count, 0.0L);
  errors.assign(2U * count, 0.0L);
  bool allConverged = true;
  std::vector<long double> xPower = SeriesConstant(1.0L, order);
  for (unsigned i = 0; i <= degree; ++i) {
    for (unsigned j = 0; j <= degree - i; ++j) {
      unsigned index = MomentIndex(degree, i, j);
      bool converged = true;
      if (computePositive) {
        SeriesProductIntegral(xPower,positiveSlice[j],values[index],errors[index],converged);
        allConverged = allConverged && converged;
      }
      if (computeDelta) {
        SeriesProductIntegral(xPower,deltaSlice[j],values[count+index],errors[count+index],converged);
        allConverged = allConverged && converged;
      }
    }
    xPower = SeriesProduct(xPower, x);
  }
  return allConverged;
}

inline void IntegrateTaylorInterval(const std::vector<std::vector<long double>>& field, const std::vector<long double>& reducedDiscriminant, unsigned leftMultiplicity, unsigned rightMultiplicity, long double xLeft, long double xRight, long double thetaLeft, long double thetaRight, const std::vector<int>& signs, const std::vector<int>& boundary, const std::vector<long double>& positiveSegments, unsigned degree, unsigned order, unsigned triangle, const std::vector<long double>& domain, bool computePositive, bool computeDelta, bool fullBoundaryDelta, long double absolutePositive, long double absoluteDelta, long double relativeTolerance, unsigned depth, unsigned maximumDepth, std::vector<long double>& sum, std::vector<long double>& sumError, unsigned& acceptedIntervals, bool fullDeltaChart = false) {
  std::vector<long double> values;
  std::vector<long double> errors;
  bool accepted = SliceTaylor(field, reducedDiscriminant, leftMultiplicity, rightMultiplicity, xLeft, xRight, thetaLeft, thetaRight, signs, boundary, positiveSegments, degree, order,triangle,domain,computePositive,computeDelta,fullBoundaryDelta,values,errors,fullDeltaChart);
  unsigned count = (degree + 1U) * (degree + 2U) / 2U;
  for (unsigned k = 0; k < values.size(); ++k) {
    long double tolerance = ((k < count) ? absolutePositive : absoluteDelta) + relativeTolerance * std::abs(values[k]);
    accepted = accepted && errors[k] <= tolerance;
  }
  if (accepted) {
    for (unsigned k = 0; k < values.size(); ++k) {
      sum[k] += values[k];
      sumError[k] += errors[k];
    }
    ++acceptedIntervals;
    return;
  }
  long double middle = 0.5L * (thetaLeft + thetaRight);
  if (depth >= maximumDepth || middle == thetaLeft || middle == thetaRight) {
    throw std::runtime_error("Taylor recurrence did not reach the requested tolerance; increase maximumDepth or relax the tolerance.");
  }
  IntegrateTaylorInterval(field,reducedDiscriminant,leftMultiplicity,rightMultiplicity,xLeft,xRight,thetaLeft,middle,signs,boundary,positiveSegments,degree,order,triangle,domain,computePositive,computeDelta,fullBoundaryDelta,0.5L*absolutePositive,0.5L*absoluteDelta,relativeTolerance,depth+1U,maximumDepth,sum,sumError,acceptedIntervals,fullDeltaChart);
  IntegrateTaylorInterval(field,reducedDiscriminant,leftMultiplicity,rightMultiplicity,xLeft,xRight,middle,thetaRight,signs,boundary,positiveSegments,degree,order,triangle,domain,computePositive,computeDelta,fullBoundaryDelta,0.5L*absolutePositive,0.5L*absoluteDelta,relativeTolerance,depth+1U,maximumDepth,sum,sumError,acceptedIntervals,fullDeltaChart);
}

inline void IntegrateDirection(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned triangle, const std::vector<long double>& domain, bool computePositive, bool computeDelta, bool fullBoundaryDelta, long double absoluteTolerance, long double deltaTolerance, long double relativeTolerance, unsigned order, unsigned maximumDepth, std::vector<long double>& result, std::vector<long double>& resultError, unsigned& acceptedIntervals, bool fullDeltaChart = false) {
  std::vector<std::vector<long double>> field(3, std::vector<long double>(3, 0.0L));
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      field[2U - j][i] = coefficient[i][j];
    }
  }
  for (unsigned k = 0; k < field.size(); ++k) {
    Trim(field[k]);
  }
  std::vector<long double> discriminant = DiscriminantPolynomial(field);
  std::vector<long double> bottom = BoundaryPolynomial(field, false, triangle, domain);
  std::vector<long double> top = BoundaryPolynomial(field, true, triangle, domain);
  std::vector<long double> events(1, domain[0]);
  events.push_back(domain[1]);
  PolynomialRoots(discriminant, domain[0], domain[1], events);
  PolynomialRoots(field[0], domain[0], domain[1], events);
  PolynomialRoots(field[1], domain[0], domain[1], events);
  PolynomialRoots(field[2], domain[0], domain[1], events);
  PolynomialRoots(bottom, domain[0], domain[1], events);
  PolynomialRoots(top, domain[0], domain[1], events);
  SortDistinct(events);
  if (computeDelta) {
    CheckSingularZeros(field, events, triangle, domain);
  }
  unsigned count = (degree + 1U) * (degree + 2U) / 2U;
  result.assign(2U * count, 0.0L);
  resultError.assign(2U * count, 0.0L);
  for (unsigned k = 1; k < events.size(); ++k) {
    long double left = events[k - 1U];
    long double right = events[k];
    if (!(right > left)) {
      continue;
    }
    long double x = left + 0.5L * (right - left);
    std::vector<int> signs;
    std::vector<int> boundary;
    std::vector<long double> positiveSegments;
    long double uniformPhase = 0.0L;
    IntervalBranches(field, x, triangle, domain, PolynomialNorm(bottom) == 0.0L, PolynomialNorm(top) == 0.0L, computeDelta, signs, boundary, positiveSegments, uniformPhase);
    if (signs.empty()) {
      if (computePositive && uniformPhase != 0.0L) {
        for (unsigned i = 0; i <= degree; ++i) {
          for (unsigned j = 0; j <= degree - i; ++j) {
            result[MomentIndex(degree, i, j)] += uniformPhase * FullStripMoment(left, right, i, j, triangle, domain);
          }
        }
      }
      continue;
    }
    std::vector<long double> reduced(discriminant);
    unsigned leftMultiplicity = 0;
    unsigned rightMultiplicity = 0;
    if (PolynomialNorm(field[0]) != 0.0L) {
      FactorEndpoint(reduced, left, leftMultiplicity);
      FactorEndpoint(reduced, right, rightMultiplicity);
      if ((rightMultiplicity & 1U) != 0) {
        for (unsigned p = 0; p < reduced.size(); ++p) {
          reduced[p] = -reduced[p];
        }
      }
      if (!(PolynomialValue(reduced, x) > 0.0L)) {
        throw std::runtime_error("Unable to resolve the discriminant chart.");
      }
    }
    long double fraction = (right - left) / (domain[1] - domain[0]);
    IntegrateTaylorInterval(field,reduced,leftMultiplicity,rightMultiplicity,left,right,0.0L,0.5L*std::acos(-1.0L),signs,boundary,positiveSegments,degree,order,triangle,domain,computePositive,computeDelta,fullBoundaryDelta,absoluteTolerance*fraction,deltaTolerance*fraction,relativeTolerance,0,maximumDepth,result,resultError,acceptedIntervals,fullDeltaChart);
  }
}

inline bool PreferReverseOrder(const std::vector<std::vector<long double>>& coefficient, const std::vector<long double>& domain) {
  long double radiusX = std::max(std::abs(domain[0]), std::abs(domain[1]));
  long double radiusY = std::max(std::abs(domain[2]), std::abs(domain[3]));
  long double derivativeX = 0.0L;
  long double derivativeY = 0.0L;
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      if (i != 0) {
        derivativeX += i * std::abs(coefficient[i][j]) * Power(radiusX, i - 1U) * Power(radiusY, j);
      }
      if (j != 0) {
        derivativeY += j * std::abs(coefficient[i][j]) * Power(radiusX, i) * Power(radiusY, j - 1U);
      }
    }
  }
  return derivativeX > derivativeY;
}

inline bool RegularDeltaProjection(const std::vector<std::vector<long double>>& coefficient, const std::vector<long double>& domain, unsigned order) {
  std::vector<std::vector<long double>> fibers(3,std::vector<long double>(3,0.0L));
  long double left = domain[order == 1 ? 0 : 2];
  long double right = domain[order == 1 ? 1 : 3];
  std::vector<long double> candidates = {left,right};
  bool selected = false;
  long double norm = 0.0L;
  long double radius = std::max(1.0L,std::max(std::abs(left),std::abs(right)));
  for(unsigned j = 0; j < 3; ++j) {
    for(unsigned i = 0; i < 3; ++i) fibers[j][i] = order == 1 ? coefficient[i][j] : coefficient[j][i];
    Trim(fibers[j]);
    norm += PolynomialBound(fibers[j],radius);
    if(PolynomialNorm(fibers[j]) != 0.0L) {
      PolynomialRoots(fibers[j],left,right,candidates);
      selected = true;
    }
  }
  if(!selected) return false;
  long double tolerance = 32.0L*std::sqrt(std::numeric_limits<long double>::epsilon())*norm;
  for(unsigned k = 0; k < candidates.size(); ++k) {
    bool common = true;
    for(unsigned j = 0; j < 3; ++j) common = common && std::abs(PolynomialValue(fibers[j],candidates[k])) <= tolerance;
    if(common) return false;
  }
  return true;
}

inline void IntegrateNormalizedCell(const std::vector<std::vector<long double>>& normalized, unsigned degree, unsigned triangle, const std::vector<long double>& domain, bool computePositive, bool computeDelta, bool fullBoundaryDelta, long double absolutePositive, long double absoluteDelta, long double relativeTolerance, unsigned seriesOrder, unsigned maximumDepth, unsigned integrationOrder, std::vector<long double>& result, std::vector<long double>& errors, unsigned& acceptedIntervals) {
  if (integrationOrder > 2U) {
    throw std::invalid_argument("integrationOrder must be 0 (automatic), 1 (y first), or 2 (x first).");
  }
  unsigned count = (degree + 1U) * (degree + 2U) / 2U;
  bool reversePositive = computePositive && ((integrationOrder == 2U) || (integrationOrder == 0U && PreferReverseOrder(normalized, domain)));
  result.assign(2U * count, 0.0L);
  errors.assign(2U * count, 0.0L);
  if(computeDelta) {
    unsigned intervalsBefore = acceptedIntervals;
    unsigned selected = integrationOrder == 0 ? (PreferReverseOrder(normalized,domain) ? 2U : 1U) : integrationOrder;
    if(!RegularDeltaProjection(normalized,domain,selected)) {
      if(integrationOrder == 0 && RegularDeltaProjection(normalized,domain,3U-selected)) selected = 3U-selected;
      else selected = 0;
    }
    try {
      if(selected == 1) {
        IntegrateDirection(normalized,degree,triangle,domain,computePositive,true,fullBoundaryDelta,absolutePositive,absoluteDelta,relativeTolerance,seriesOrder,maximumDepth,result,errors,acceptedIntervals,true);
        return;
      }
      if(selected == 2) {
        std::vector<std::vector<long double>> transposed(3,std::vector<long double>(3));
        for(unsigned i = 0; i < 3; ++i) {
          for(unsigned j = 0; j < 3; ++j) transposed[i][j] = normalized[j][i];
        }
        std::vector<long double> swapped(domain);
        std::swap(swapped[0],swapped[2]);
        std::swap(swapped[1],swapped[3]);
        unsigned shape = triangle == 2 ? 3U : (triangle == 3 ? 2U : triangle);
        std::vector<long double> value,error;
        IntegrateDirection(transposed,degree,shape,swapped,computePositive,true,fullBoundaryDelta,absolutePositive,absoluteDelta,relativeTolerance,seriesOrder,maximumDepth,value,error,acceptedIntervals,true);
        for(unsigned i = 0; i <= degree; ++i) {
          for(unsigned j = 0; j <= degree-i; ++j) {
            unsigned target = MomentIndex(degree,i,j);
            unsigned source = MomentIndex(degree,j,i);
            result[target] = value[source];
            result[count+target] = value[count+source];
            errors[target] = error[source];
            errors[count+target] = error[count+source];
          }
        }
        return;
      }
    }
    catch(const std::runtime_error&) {
      std::fill(result.begin(),result.end(),0.0L);
      std::fill(errors.begin(),errors.end(),0.0L);
      acceptedIntervals = intervalsBefore;
    }
  }
  if ((computePositive && !reversePositive) || computeDelta) {
    IntegrateDirection(normalized, degree, triangle, domain, computePositive && !reversePositive, computeDelta, fullBoundaryDelta, absolutePositive, 0.5L * absoluteDelta, relativeTolerance, seriesOrder, maximumDepth, result, errors, acceptedIntervals);
  }
  if (computeDelta || reversePositive) {
    std::vector<std::vector<long double>> transposed(3, std::vector<long double>(3, 0.0L));
    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j) {
        transposed[i][j] = normalized[j][i];
      }
    }
    std::vector<long double> swappedDomain(domain);
    std::swap(swappedDomain[0], swappedDomain[2]);
    std::swap(swappedDomain[1], swappedDomain[3]);
    unsigned swappedTriangle = (triangle == 2U) ? 3U : ((triangle == 3U) ? 2U : triangle);
    std::vector<long double> reverse;
    std::vector<long double> reverseError;
    IntegrateDirection(transposed, degree, swappedTriangle, swappedDomain, reversePositive, computeDelta, fullBoundaryDelta, absolutePositive, 0.5L * absoluteDelta, relativeTolerance, seriesOrder, maximumDepth, reverse, reverseError, acceptedIntervals);
    for (unsigned i = 0; i <= degree; ++i) {
      for (unsigned j = 0; j <= degree - i; ++j) {
        unsigned original = MomentIndex(degree, i, j);
        unsigned exchanged = MomentIndex(degree, j, i);
        if (reversePositive) {
          result[original] = reverse[exchanged];
          errors[original] = reverseError[exchanged];
        }
        if (computeDelta) {
          result[count + original] += reverse[count + exchanged];
          errors[count + original] += reverseError[count + exchanged];
        }
      }
    }
  }
}

inline void IntegrateCell(const std::vector<std::vector<double>>& coefficient, unsigned degree, unsigned triangle, const std::vector<long double>& domain, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance, double relativeTolerance, bool computeDelta, bool fullBoundaryDelta, unsigned seriesOrder, unsigned maximumDepth, std::vector<double>* estimatedError, unsigned integrationOrder = 0) {
  if (coefficient.size() != 3 || triangle > 3U || domain.size() != 4 || degree > std::numeric_limits<unsigned>::max() - 2U || !(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance) || seriesOrder < 24 || seriesOrder > 256 || maximumDepth > 60 || integrationOrder > 2U) {
    throw std::invalid_argument("Invalid cell integration input.");
  }
  unsigned long long countCheck = (static_cast<unsigned long long>(degree) + 1ULL) * (static_cast<unsigned long long>(degree) + 2ULL) / 2ULL;
  if (countCheck > std::numeric_limits<unsigned>::max() / 2U) {
    throw std::length_error("The requested moment table is too large.");
  }
  for (unsigned k = 0; k < domain.size(); ++k) {
    if (!std::isfinite(domain[k])) {
      throw std::invalid_argument("Invalid reference domain.");
    }
  }
  if (!(domain[1] > domain[0]) || !(domain[3] > domain[2])) {
    throw std::invalid_argument("The reference domain must have positive size.");
  }
  std::vector<std::vector<long double>> normalized(3, std::vector<long double>(3, 0.0L));
  long double scale = 0.0L;
  for (unsigned i = 0; i < 3; ++i) {
    if (coefficient[i].size() != 3) {
      throw std::invalid_argument("QUAD9 coefficients must be a 3 by 3 matrix.");
    }
    for (unsigned j = 0; j < 3; ++j) {
      if (!std::isfinite(coefficient[i][j])) {
        throw std::invalid_argument("A level-set coefficient is not finite.");
      }
      scale = std::max(scale, std::abs(static_cast<long double>(coefficient[i][j])));
    }
  }
  if (scale == 0.0L) {
    if (computeDelta) {
      throw std::domain_error("The delta of an identically zero field is undefined.");
    }
    scale = 1.0L;
  }
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      normalized[i][j] = coefficient[i][j] / scale;
    }
  }
  std::vector<long double> result;
  std::vector<long double> errors;
  unsigned acceptedIntervals = 0;
  IntegrateNormalizedCell(normalized, degree, triangle, domain, true, computeDelta, fullBoundaryDelta, absoluteTolerance, absoluteTolerance * scale, relativeTolerance, seriesOrder, maximumDepth, integrationOrder, result, errors, acceptedIntervals);
  unsigned count = static_cast<unsigned>(countCheck);
  positive.resize(degree + 1U);
  delta.resize(degree + 1U);
  long double maximumPositiveError = 0.0L;
  long double maximumDeltaError = 0.0L;
  for (unsigned i = 0; i <= degree; ++i) {
    positive[i].resize(degree - i + 1U);
    delta[i].resize(degree - i + 1U);
    for (unsigned j = 0; j <= degree - i; ++j) {
      unsigned index = MomentIndex(degree, i, j);
      positive[i][j] = static_cast<double>(result[index]);
      delta[i][j] = computeDelta ? static_cast<double>(result[count + index] / scale) : 0.0;
      if (!std::isfinite(positive[i][j]) || !std::isfinite(delta[i][j])) {
        throw std::overflow_error("A cell moment exceeds double precision range.");
      }
      maximumPositiveError = std::max(maximumPositiveError, errors[index]);
      maximumDeltaError = std::max(maximumDeltaError, errors[count + index] / scale);
    }
  }
  if (estimatedError != nullptr) {
    estimatedError->resize(3);
    (*estimatedError)[0] = static_cast<double>(maximumPositiveError);
    (*estimatedError)[1] = static_cast<double>(maximumDeltaError);
    (*estimatedError)[2] = acceptedIntervals;
  }
}

inline void IntegrateQuad9(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, double xLeft = -1.0, double xRight = 1.0, double yBottom = -1.0, double yTop = 1.0, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned seriesOrder = 32, unsigned maximumDepth = 24, std::vector<double>* estimatedError = nullptr) {
  std::vector<long double> domain(4, 0.0L);
  domain[0] = xLeft;
  domain[1] = xRight;
  domain[2] = yBottom;
  domain[3] = yTop;
  IntegrateCell(coefficient, degree, false, domain, positive, delta, absoluteTolerance, relativeTolerance, computeDelta, fullBoundaryDelta, seriesOrder, maximumDepth, estimatedError);
}

inline void IntegrateTri7(const std::vector<double>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned seriesOrder = 32, unsigned maximumDepth = 24, std::vector<double>* estimatedError = nullptr) {
  if (coefficient.size() != 7) {
    throw std::invalid_argument("TRI7 requires coefficients of 1,x,y,x*x,x*y,y*y,x*y*(x+y).");
  }
  std::vector<std::vector<double>> polynomial(3, std::vector<double>(3, 0.0));
  polynomial[0][0] = coefficient[0];
  polynomial[1][0] = coefficient[1];
  polynomial[0][1] = coefficient[2];
  polynomial[2][0] = coefficient[3];
  polynomial[1][1] = coefficient[4];
  polynomial[0][2] = coefficient[5];
  polynomial[2][1] = coefficient[6];
  polynomial[1][2] = coefficient[6];
  std::vector<long double> domain(4, 0.0L);
  domain[1] = 1.0L;
  domain[3] = 1.0L;
  IntegrateCell(polynomial, degree, true, domain, positive, delta, absoluteTolerance, relativeTolerance, computeDelta, fullBoundaryDelta, seriesOrder, maximumDepth, estimatedError);
}

inline void NormalizeCanonicalPolynomial(const std::vector<std::vector<double>>& coefficient, std::vector<std::vector<long double>>& normalized, long double& scale) {
  if (coefficient.size() != 3) {
    throw std::invalid_argument("The coefficient matrix must be 3 by 3.");
  }
  scale = 0.0L;
  for (unsigned i = 0; i < 3; ++i) {
    if (coefficient[i].size() != 3) {
      throw std::invalid_argument("The coefficient matrix must be 3 by 3.");
    }
    for (unsigned j = 0; j < 3; ++j) {
      if (!std::isfinite(coefficient[i][j])) {
        throw std::invalid_argument("A level-set coefficient is not finite.");
      }
      scale = std::max(scale, std::abs(static_cast<long double>(coefficient[i][j])));
    }
  }
  if (scale == 0.0L) {
    scale = 1.0L;
  }
  normalized.assign(3, std::vector<long double>(3, 0.0L));
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      normalized[i][j] = coefficient[i][j] / scale;
    }
  }
}

inline long double CanonicalFullMoment(unsigned i, unsigned j, bool triangle) {
  return triangle ? 1.0L / (static_cast<long double>(j + 1U) * (i + j + 2U)) : 1.0L / (static_cast<long double>(i + 1U) * (j + 1U));
}

inline void CanonicalBoundaryMoments(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool triangle, std::vector<long double>& right, std::vector<long double>& top) {
  std::vector<long double> unused;
  long double a = coefficient[0][2] + coefficient[1][2] + coefficient[2][2];
  long double b = coefficient[0][1] + coefficient[1][1] + coefficient[2][1];
  long double c = coefficient[0][0] + coefficient[1][0] + coefficient[2][0];
  Integrate1DInternal(a, b, c, degree, 0.0L, 1.0L, false, false, right, unused);
  top.assign(degree + 1U, 0.0L);
  if (!triangle) {
    a = coefficient[2][0] + coefficient[2][1] + coefficient[2][2];
    b = coefficient[1][0] + coefficient[1][1] + coefficient[1][2];
    c = coefficient[0][0] + coefficient[0][1] + coefficient[0][2];
    Integrate1DInternal(a, b, c, degree, 0.0L, 1.0L, false, false, top, unused);
  }
}

inline void CanonicalBernstein(const std::vector<std::vector<long double>>& coefficient, long double bernstein[3][3]) {
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      bernstein[i][j] = 0.0L;
      for(unsigned r = 0; r <= i; ++r) {
        for(unsigned s = 0; s <= j; ++s) bernstein[i][j] += coefficient[r][s]*Binomial(i,r)*Binomial(j,s)/(Binomial(2,r)*Binomial(2,s));
      }
    }
  }
}

inline unsigned CanonicalGraphOrder(const std::vector<std::vector<long double>>& coefficient, unsigned requestedOrder, int& phase, long double* derivativeBound = nullptr) {
  long double b[3][3];
  CanonicalBernstein(coefficient,b);
  long double norm = 0.0L;
  long double minimum = b[0][0];
  long double maximum = minimum;
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      norm += std::abs(coefficient[i][j]);
      minimum = std::min(minimum,b[i][j]);
      maximum = std::max(maximum,b[i][j]);
    }
  }
  long double padding = 512.0L*std::numeric_limits<long double>::epsilon()*norm;
  long double conditioning = 32.0L*std::sqrt(std::numeric_limits<long double>::epsilon())*norm;
  phase = minimum > padding ? 1 : (maximum < -padding ? -1 : 0);
  unsigned selected = 0;
  long double selectedBound = 0.0L;
  for(unsigned order = 1; order <= 2; ++order) {
    if(requestedOrder != 0 && requestedOrder != order) continue;
    long double derivativeMinimum = std::numeric_limits<long double>::infinity();
    long double derivativeMaximum = -derivativeMinimum;
    for(unsigned i = 0; i < 3; ++i) {
      long double low = order == 1 ? b[i][0] : b[0][i];
      long double middle = order == 1 ? b[i][1] : b[1][i];
      long double high = order == 1 ? b[i][2] : b[2][i];
      const long double derivative[3] = {2.0L*(middle-low),high-low,2.0L*(high-middle)};
      for(unsigned j = 0; j < 3; ++j) {
        derivativeMinimum = std::min(derivativeMinimum,derivative[j]);
        derivativeMaximum = std::max(derivativeMaximum,derivative[j]);
      }
    }
    long double bound = 0.0L;
    if(derivativeMinimum > conditioning) bound = derivativeMinimum-padding;
    if(derivativeMaximum < -conditioning) bound = -derivativeMaximum-padding;
    if(bound > selectedBound) {
      selectedBound = bound;
      selected = order;
    }
  }
  if(derivativeBound != nullptr) *derivativeBound = selectedBound;
  return selected;
}

inline void CanonicalGraphSample(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned order, long double point, bool computeDelta, std::vector<long double>& values, bool triangle = false) {
  long double polynomial[3] = {};
  for(unsigned k = 0; k < 3; ++k) {
    if(order == 1) polynomial[k] = std::fma(std::fma(coefficient[2][k],point,coefficient[1][k]),point,coefficient[0][k]);
    else polynomial[k] = std::fma(std::fma(coefficient[k][2],point,coefficient[k][1]),point,coefficient[k][0]);
  }
  std::vector<long double> roots = QuadraticRoots(polynomial[2],polynomial[1],polynomial[0]);
  long double lower = triangle && order == 2 ? point : 0.0L;
  long double upper = triangle && order == 1 ? point : 1.0L;
  long double root = -1.0L;
  for(unsigned k = 0; k < roots.size(); ++k) {
    if(roots[k] > lower && roots[k] < upper) {
      if(root >= 0.0L) throw std::runtime_error("A monotone graph has more than one interior root.");
      root = roots[k];
    }
  }
  if(root < 0.0L) throw std::runtime_error("A monotone graph root was not resolved.");
  long double derivative = std::abs(std::fma(2.0L*polynomial[2],root,polynomial[1]));
  if(!(derivative > 0.0L)) throw std::runtime_error("A monotone graph derivative was not resolved.");
  std::vector<long double> positive(degree+1U),delta(degree+1U);
  bool above = std::fma(polynomial[2],lower*lower,std::fma(polynomial[1],lower,polynomial[0])) < 0.0L;
  long double left = above ? root : lower;
  long double right = above ? upper : root;
  long double power = 1.0L;
  long double rightPower = 1.0L;
  long double geometric = 0.0L;
  for(unsigned j = 0; j <= degree; ++j) {
    geometric = rightPower+left*geometric;
    delta[j] = computeDelta ? power/derivative : 0.0L;
    power *= root;
    rightPower *= right;
    positive[j] = (right-left)*geometric/(j+1U);
  }
  unsigned count = (degree+1U)*(degree+2U)/2U;
  values.assign(2U*count,0.0L);
  power = 1.0L;
  for(unsigned i = 0; i <= degree; ++i) {
    for(unsigned j = 0; j <= degree-i; ++j) {
      unsigned index = order == 1 ? MomentIndex(degree,i,j) : MomentIndex(degree,j,i);
      values[index] = power*positive[j];
      values[count+index] = power*delta[j];
    }
    power *= point;
  }
}

inline void IntegrateCanonicalGraphInterval(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned order, long double left, long double right, bool computeDelta, long double absoluteH, long double absoluteDelta, long double relativeTolerance, unsigned depth, unsigned maximumDepth, std::vector<long double>& total, std::vector<long double>& totalError, unsigned& accepted, bool triangle = false) {
  const long double nodes[8] = {0.991455371120812639206854697526329L,0.949107912342758524526189684047851L,0.864864423359769072789712788640926L,0.741531185599394439863864773280788L,0.586087235467691130294144838258730L,0.405845151377397166906606412076961L,0.207784955007898467600689403773245L,0.0L};
  const long double kronrod[8] = {0.022935322010529224963732008058970L,0.063092092629978553290700663189204L,0.104790010322250183839876322541518L,0.140653259715525918745189590510238L,0.169004726639267902826583426598550L,0.190350578064785409913256402421014L,0.204432940075298892414161999234649L,0.209482141084727828012999174891714L};
  const long double gauss[4] = {0.129484966168869693270611432679082L,0.279705391489276667901467771423780L,0.381830050505118944950369775488975L,0.417959183673469387755102040816327L};
  unsigned count = (degree+1U)*(degree+2U)/2U;
  std::vector<long double> high(2U*count,0.0L),low(2U*count,0.0L),absolute(2U*count,0.0L);
  std::vector<std::vector<long double>> samples;
  std::vector<long double> sampleWeights;
  long double middle = 0.5L*(left+right);
  long double halfWidth = 0.5L*(right-left);
  for(unsigned node = 0; node < 8; ++node) {
    unsigned sides = node == 7 ? 1 : 2;
    for(unsigned side = 0; side < sides; ++side) {
      long double point = middle+(side == 0 ? -1.0L : 1.0L)*halfWidth*nodes[node];
      std::vector<long double> value;
      CanonicalGraphSample(coefficient,degree,order,point,computeDelta,value,triangle);
      for(unsigned k = 0; k < value.size(); ++k) {
        high[k] += kronrod[node]*value[k];
        absolute[k] += kronrod[node]*std::abs(value[k]);
        if((node&1U) != 0) low[k] += gauss[node/2U]*value[k];
      }
      samples.push_back(value);
      sampleWeights.push_back(kronrod[node]);
    }
  }
  bool converged = true;
  std::vector<long double> error(high.size());
  for(unsigned k = 0; k < high.size(); ++k) {
    long double variation = 0.0L;
    for(unsigned n = 0; n < samples.size(); ++n) variation += sampleWeights[n]*std::abs(samples[n][k]-0.5L*high[k]);
    long double estimate = std::abs(high[k]-low[k]);
    if(variation > 0.0L && estimate > 0.0L) estimate = variation*std::min(1.0L,std::pow(200.0L*estimate/variation,1.5L));
    estimate = std::max(estimate,128.0L*std::numeric_limits<long double>::epsilon()*absolute[k]);
    error[k] = halfWidth*estimate;
    high[k] *= halfWidth;
    long double tolerance = (k < count ? absoluteH : absoluteDelta)*(right-left)+relativeTolerance*std::abs(high[k]);
    converged = converged && error[k] <= tolerance;
  }
  if(converged) {
    for(unsigned k = 0; k < high.size(); ++k) {
      total[k] += high[k];
      totalError[k] += error[k];
    }
    ++accepted;
    return;
  }
  if(depth >= maximumDepth || middle == left || middle == right) throw std::runtime_error("Monotone graph integration did not reach the requested tolerance.");
  IntegrateCanonicalGraphInterval(coefficient,degree,order,left,middle,computeDelta,absoluteH,absoluteDelta,relativeTolerance,depth+1U,maximumDepth,total,totalError,accepted,triangle);
  IntegrateCanonicalGraphInterval(coefficient,degree,order,middle,right,computeDelta,absoluteH,absoluteDelta,relativeTolerance,depth+1U,maximumDepth,total,totalError,accepted,triangle);
}

inline bool IntegrateCanonicalGraph(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned order, bool triangle, bool computeDelta, long double derivativeBound, long double absoluteH, long double absoluteDelta, long double relativeTolerance, unsigned maximumDepth, std::vector<long double>& total, std::vector<long double>& totalError, unsigned& accepted) {
  const std::vector<long double> domain = {0.0L,1.0L,0.0L,1.0L};
  std::vector<std::vector<long double>> field(3,std::vector<long double>(3));
  long double norm = 0.0L;
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      field[2U-j][i] = order == 1 ? coefficient[i][j] : coefficient[j][i];
      norm += std::abs(coefficient[i][j]);
    }
  }
  unsigned shape = triangle ? (order == 1 ? 2U : 3U) : 0U;
  std::vector<std::vector<long double>> edges = {BoundaryPolynomial(field,false,shape,domain),BoundaryPolynomial(field,true,shape,domain)};
  std::vector<long double> events = {0.0L,1.0L};
  long double rootError = 0.0L;
  long double rounding = 512.0L*std::numeric_limits<long double>::epsilon();
  for(unsigned side = 0; side < 2; ++side) {
    Trim(edges[side]);
    if(edges[side].size() > 3 || PolynomialNorm(edges[side]) <= rounding*norm) return false;
    if(edges[side].size() == 3) {
      long double a = edges[side][2], b = edges[side][1], c = edges[side][0];
      long double discriminant = std::fma(-4.0L*a,c,b*b);
      if(std::abs(discriminant) <= rounding*(b*b+4.0L*std::abs(a*c))) return false;
    }
    std::vector<long double> roots;
    PolynomialRoots(edges[side],0.0L,1.0L,roots);
    std::vector<long double> derivative = PolynomialDerivative(edges[side]);
    for(unsigned k = 0; k < roots.size(); ++k) {
      long double error = std::abs(PolynomialValue(edges[side],roots[k]))+rounding*PolynomialBound(edges[side],roots[k]);
      long double slope = std::abs(PolynomialValue(derivative,roots[k]))-rounding*PolynomialBound(derivative,roots[k]);
      long double quadratic = edges[side].size() == 3 ? std::abs(edges[side][2]) : 0.0L;
      if(!(slope > 0.0L) || slope*slope <= 8.0L*quadratic*error) return false;
      rootError += 2.0L*error/slope;
      events.push_back(roots[k]);
    }
  }
  if(rootError > 0.125L*absoluteH || (computeDelta && rootError/derivativeBound > 0.125L*absoluteDelta)) return false;
  SortDistinct(events);
  unsigned count = (degree+1U)*(degree+2U)/2U;
  for(unsigned interval = 1; interval < events.size(); ++interval) {
    long double left = events[interval-1U], right = events[interval];
    if(!(right > left)) continue;
    long double middle = 0.5L*(left+right);
    long double bottom = PolynomialValue(edges[0],middle);
    long double top = PolynomialValue(edges[1],middle);
    if(bottom == 0.0L || top == 0.0L) return false;
    if((bottom > 0.0L) == (top > 0.0L)) {
      if(bottom > 0.0L) {
        for(unsigned i = 0; i <= degree; ++i) {
          for(unsigned j = 0; j <= degree-i; ++j) {
            unsigned target = order == 1 ? MomentIndex(degree,i,j) : MomentIndex(degree,j,i);
            total[target] += FullStripMoment(left,right,i,j,shape,domain);
          }
        }
      }
      ++accepted;
    }
    else IntegrateCanonicalGraphInterval(coefficient,degree,order,left,right,computeDelta,0.875L*absoluteH,0.875L*absoluteDelta,relativeTolerance,0,maximumDepth,total,totalError,accepted,triangle);
  }
  for(unsigned k = 0; k < count; ++k) {
    totalError[k] += rootError;
    if(computeDelta) totalError[count+k] += rootError/derivativeBound;
  }
  return true;
}

inline void IntegrateCanonicalPolynomial(const std::vector<std::vector<double>>& coefficient, unsigned degree, bool triangle, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance, double relativeTolerance, bool computeDelta, bool fullBoundaryDelta, unsigned seriesOrder, unsigned maximumDepth, unsigned integrationOrder, bool useBoundaryRecurrence, std::vector<double>* estimatedError, std::vector<long double>* momentErrors = nullptr) {
  if (degree > std::numeric_limits<unsigned>::max() - 6U || !(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance) || seriesOrder < 24U || seriesOrder > 256U || maximumDepth > 60U || integrationOrder > 2U) {
    throw std::invalid_argument("Invalid canonical integration input.");
  }
  std::vector<std::vector<long double>> normalized;
  long double scale = 1.0L;
  NormalizeCanonicalPolynomial(coefficient, normalized, scale);
  unsigned extraDegree = 0;
  long double coefficientWeight = 0.0L;
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      if (i + j != 2U && normalized[i][j] != 0.0L) {
        extraDegree = std::max(extraDegree, i + j);
        coefficientWeight += std::abs((2.0L - i - j) * normalized[i][j]);
      }
    }
  }
  unsigned internalDegree = useBoundaryRecurrence ? degree + extraDegree : degree;
  unsigned long long check = (static_cast<unsigned long long>(internalDegree) + 1ULL) * (static_cast<unsigned long long>(internalDegree) + 2ULL) / 2ULL;
  if (check > std::numeric_limits<unsigned>::max() / 2U) {
    throw std::length_error("The requested moment table is too large.");
  }
  std::vector<long double> domain = {0.0L, 1.0L, 0.0L, 1.0L};
  unsigned shape = triangle ? 2U : 0U;
  unsigned count = (degree + 1U) * (degree + 2U) / 2U;
  std::vector<long double> values(2U * count, 0.0L);
  std::vector<long double> valueError(2U * count, 0.0L);
  std::vector<long double> high;
  std::vector<long double> highError;
  unsigned acceptedIntervals = 0;
  bool usedRecurrence = false;
  bool completed = false;
  int phase = 0;
  long double graphDerivative = 0.0L;
  unsigned graphOrder = CanonicalGraphOrder(normalized,integrationOrder,phase,&graphDerivative);
  if(phase != 0) {
    if(phase > 0) {
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) values[MomentIndex(degree,i,j)] = CanonicalFullMoment(i,j,triangle);
      }
    }
    completed = true;
  }
  else if(!useBoundaryRecurrence && graphOrder != 0) {
    try {
      completed = IntegrateCanonicalGraph(normalized,degree,graphOrder,triangle,computeDelta,graphDerivative,absoluteTolerance,static_cast<long double>(absoluteTolerance)*scale,relativeTolerance,maximumDepth,values,valueError,acceptedIntervals);
    }
    catch(const std::runtime_error&) {
    }
    if(!completed) {
      std::fill(values.begin(),values.end(),0.0L);
      std::fill(valueError.begin(),valueError.end(),0.0L);
      acceptedIntervals = 0;
    }
  }
  if (!completed && useBoundaryRecurrence) {
    try {
      long double deltaTolerance = std::min(static_cast<long double>(absoluteTolerance) * scale, absoluteTolerance / (1.0L + coefficientWeight));
      IntegrateNormalizedCell(normalized, internalDegree, shape, domain, false, true, false, absoluteTolerance, deltaTolerance, relativeTolerance / 8.0L, seriesOrder, maximumDepth, integrationOrder, high, highError, acceptedIntervals);
      usedRecurrence = true;
    }
    catch (const std::domain_error&) {
      if (computeDelta) {
        throw;
      }
    }
    catch (const std::runtime_error&) {
      usedRecurrence = false;
    }
  }
  if (!completed && !usedRecurrence) {
    IntegrateNormalizedCell(normalized, degree, shape, domain, true, computeDelta, fullBoundaryDelta, absoluteTolerance, absoluteTolerance * scale, relativeTolerance, seriesOrder, maximumDepth, integrationOrder, values, valueError, acceptedIntervals);
  }
  else if(!completed) {
    unsigned highCount = (internalDegree + 1U) * (internalDegree + 2U) / 2U;
    std::vector<long double> right;
    std::vector<long double> top;
    CanonicalBoundaryMoments(normalized, degree, triangle, right, top);
    bool recomputePositive = false;
    for (unsigned i = 0; i <= degree; ++i) {
      for (unsigned j = 0; j <= degree - i; ++j) {
        unsigned index = MomentIndex(degree, i, j);
        long double value = right[j] + top[i];
        long double norm = std::abs(right[j]) + std::abs(top[i]);
        long double error = 0.0L;
        long double correction = 0.0L;
        for (unsigned r = 0; r < 3; ++r) {
          for (unsigned s = 0; s < 3; ++s) {
            long double factor = (2.0L - r - s) * normalized[r][s];
            if (factor != 0.0L) {
              unsigned source = highCount + MomentIndex(internalDegree, i + r, j + s);
              long double term = factor * high[source];
              long double corrected = term - correction;
              long double next = value + corrected;
              correction = (next - value) - corrected;
              value = next;
              norm += std::abs(term);
              error += std::abs(factor) * highError[source];
            }
          }
        }
        values[index] = value / (i + j + 2U);
        valueError[index] = (error + 128.0L * std::numeric_limits<long double>::epsilon() * norm) / (i + j + 2U);
        long double tolerance = absoluteTolerance + relativeTolerance * std::abs(values[index]);
        long double full = CanonicalFullMoment(i, j, triangle);
        recomputePositive = recomputePositive || valueError[index] > tolerance || values[index] < -tolerance || values[index] > full + tolerance;
        if (computeDelta) {
          unsigned source = highCount + MomentIndex(internalDegree, i, j);
          values[count + index] = high[source];
          valueError[count + index] = highError[source];
        }
      }
    }
    if (recomputePositive) {
      std::vector<long double> direct;
      std::vector<long double> directError;
      IntegrateNormalizedCell(normalized, degree, shape, domain, true, false, false, absoluteTolerance, absoluteTolerance * scale, relativeTolerance, seriesOrder, maximumDepth, integrationOrder, direct, directError, acceptedIntervals);
      for (unsigned k = 0; k < count; ++k) {
        values[k] = direct[k];
        valueError[k] = directError[k];
      }
    }
    if (computeDelta && fullBoundaryDelta) {
      std::vector<long double> full;
      std::vector<long double> fullError;
      IntegrateNormalizedCell(normalized, degree, shape, domain, false, true, true, absoluteTolerance, absoluteTolerance * scale, relativeTolerance, seriesOrder, maximumDepth, integrationOrder, full, fullError, acceptedIntervals);
      for (unsigned k = count; k < 2U * count; ++k) {
        values[k] = full[k];
        valueError[k] = fullError[k];
      }
    }
  }
  positive.resize(degree + 1U);
  delta.resize(degree + 1U);
  long double positiveError = 0.0L;
  long double deltaError = 0.0L;
  for (unsigned i = 0; i <= degree; ++i) {
    positive[i].resize(degree - i + 1U);
    delta[i].resize(degree - i + 1U);
    for (unsigned j = 0; j <= degree - i; ++j) {
      unsigned index = MomentIndex(degree, i, j);
      positive[i][j] = static_cast<double>(values[index]);
      delta[i][j] = computeDelta ? static_cast<double>(values[count + index] / scale) : 0.0;
      if (!std::isfinite(positive[i][j]) || !std::isfinite(delta[i][j])) {
        throw std::overflow_error("A canonical moment exceeds double precision range.");
      }
      positiveError = std::max(positiveError, valueError[index]);
      if (computeDelta) {
        deltaError = std::max(deltaError, valueError[count + index] / scale);
      }
    }
  }
  if(momentErrors != nullptr) {
    *momentErrors = valueError;
    for(unsigned k = count; k < 2U*count; ++k) {
      (*momentErrors)[k] = computeDelta ? valueError[k]/scale : 0.0L;
    }
  }
  if (estimatedError != nullptr) {
    estimatedError->assign(3, 0.0);
    (*estimatedError)[0] = static_cast<double>(positiveError);
    (*estimatedError)[1] = static_cast<double>(deltaError);
    (*estimatedError)[2] = acceptedIntervals;
  }
}

inline void Quad9CoefficientsToUnitSquare(const std::vector<std::vector<double>>& coefficient, std::vector<std::vector<double>>& transformed, double xLeft = -1.0, double xRight = 1.0, double yBottom = -1.0, double yTop = 1.0) {
  if (!(xRight > xLeft) || !(yTop > yBottom) || !std::isfinite(xLeft) || !std::isfinite(xRight) || !std::isfinite(yBottom) || !std::isfinite(yTop)) {
    throw std::invalid_argument("Invalid QUAD9 coordinate map.");
  }
  std::vector<std::vector<long double>> checked;
  long double scale = 1.0L;
  NormalizeCanonicalPolynomial(coefficient, checked, scale);
  std::vector<std::vector<double>> result(3, std::vector<double>(3, 0.0));
  long double widthX = static_cast<long double>(xRight) - xLeft;
  long double widthY = static_cast<long double>(yTop) - yBottom;
  for (unsigned r = 0; r < 3; ++r) {
    for (unsigned s = 0; s < 3; ++s) {
      long double value = 0.0L;
      for (unsigned i = r; i < 3; ++i) {
        for (unsigned j = s; j < 3; ++j) {
          value += static_cast<long double>(coefficient[i][j]) * Binomial(i, r) * Binomial(j, s) * Power(xLeft, i - r) * Power(yBottom, j - s);
        }
      }
      value *= Power(widthX, r) * Power(widthY, s);
      result[r][s] = static_cast<double>(value);
      if (!std::isfinite(result[r][s])) {
        throw std::overflow_error("A transformed QUAD9 coefficient exceeds double precision range.");
      }
    }
  }
  transformed.swap(result);
}

inline void Tri7CoefficientsToReflectedTriangle(const std::vector<double>& coefficient, std::vector<double>& transformed) {
  if (coefficient.size() != 7) {
    throw std::invalid_argument("TRI7 requires seven coefficients.");
  }
  for (unsigned k = 0; k < 7; ++k) {
    if (!std::isfinite(coefficient[k])) {
      throw std::invalid_argument("A TRI7 coefficient is not finite.");
    }
  }
  std::vector<long double> value(7, 0.0L);
  value[0] = static_cast<long double>(coefficient[0]) + coefficient[1] + coefficient[3];
  value[1] = -static_cast<long double>(coefficient[1]) - 2.0L * coefficient[3];
  value[2] = static_cast<long double>(coefficient[2]) + coefficient[4] + coefficient[6];
  value[3] = coefficient[3];
  value[4] = -static_cast<long double>(coefficient[4]) - 2.0L * coefficient[6];
  value[5] = static_cast<long double>(coefficient[5]) + coefficient[6];
  value[6] = coefficient[6];
  transformed.resize(7);
  for (unsigned k = 0; k < 7; ++k) {
    transformed[k] = static_cast<double>(value[k]);
    if (!std::isfinite(transformed[k])) {
      throw std::overflow_error("A reflected TRI7 coefficient exceeds double precision range.");
    }
  }
}

inline void IntegrateQuad9Unit(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned seriesOrder = 32, unsigned maximumDepth = 24, unsigned integrationOrder = 0, bool useBoundaryRecurrence = false, std::vector<double>* estimatedError = nullptr) {
  IntegrateCanonicalPolynomial(coefficient, degree, false, positive, delta, absoluteTolerance, relativeTolerance, computeDelta, fullBoundaryDelta, seriesOrder, maximumDepth, integrationOrder, useBoundaryRecurrence, estimatedError);
}

inline void IntegrateTri7Reflected(const std::vector<double>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned seriesOrder = 32, unsigned maximumDepth = 24, unsigned integrationOrder = 0, bool useBoundaryRecurrence = false, std::vector<double>* estimatedError = nullptr) {
  if (coefficient.size() != 7) {
    throw std::invalid_argument("Reflected TRI7 requires coefficients of 1,x,y,x*x,x*y,y*y,x*y*(x-y).");
  }
  std::vector<std::vector<double>> polynomial(3, std::vector<double>(3, 0.0));
  polynomial[0][0] = coefficient[0];
  polynomial[1][0] = coefficient[1];
  polynomial[0][1] = coefficient[2];
  polynomial[2][0] = coefficient[3];
  polynomial[1][1] = coefficient[4];
  polynomial[0][2] = coefficient[5];
  polynomial[2][1] = coefficient[6];
  polynomial[1][2] = -coefficient[6];
  IntegrateCanonicalPolynomial(polynomial, degree, true, positive, delta, absoluteTolerance, relativeTolerance, computeDelta, fullBoundaryDelta, seriesOrder, maximumDepth, integrationOrder, useBoundaryRecurrence, estimatedError);
}

inline void IntegrateQuad9Canonical(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, double xLeft = -1.0, double xRight = 1.0, double yBottom = -1.0, double yTop = 1.0, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned seriesOrder = 32, unsigned maximumDepth = 24, unsigned integrationOrder = 0, bool useBoundaryRecurrence = false, std::vector<double>* estimatedError = nullptr) {
  std::vector<std::vector<double>> transformed;
  Quad9CoefficientsToUnitSquare(coefficient, transformed, xLeft, xRight, yBottom, yTop);
  IntegrateQuad9Unit(transformed, degree, positive, delta, absoluteTolerance, relativeTolerance, computeDelta, fullBoundaryDelta, seriesOrder, maximumDepth, integrationOrder, useBoundaryRecurrence, estimatedError);
}

inline void IntegrateTri7Canonical(const std::vector<double>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned seriesOrder = 32, unsigned maximumDepth = 24, unsigned integrationOrder = 0, bool useBoundaryRecurrence = false, std::vector<double>* estimatedError = nullptr) {
  std::vector<double> transformed;
  Tri7CoefficientsToReflectedTriangle(coefficient, transformed);
  IntegrateTri7Reflected(transformed, degree, positive, delta, absoluteTolerance, relativeTolerance, computeDelta, fullBoundaryDelta, seriesOrder, maximumDepth, integrationOrder, useBoundaryRecurrence, estimatedError);
}

inline void Integrate1DCanonical(double a, double b, double c, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, double left = -1.0, double right = 1.0, bool computeDelta = true, bool fullBoundaryDelta = false) {
  if (!(right > left) || !std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(left) || !std::isfinite(right) || degree == std::numeric_limits<unsigned>::max()) {
    throw std::invalid_argument("Invalid 1D canonical input.");
  }
  long double width = static_cast<long double>(right) - left;
  long double aa = static_cast<long double>(a) * width * width;
  long double bb = (2.0L * a * left + b) * width;
  long double cc = std::fma(static_cast<long double>(a), static_cast<long double>(left) * left, std::fma(static_cast<long double>(b), left, static_cast<long double>(c)));
  long double scale = std::max(std::abs(aa), std::max(std::abs(bb), std::abs(cc)));
  if (scale == 0.0L) {
    scale = 1.0L;
  }
  std::vector<long double> h;
  std::vector<long double> d;
  Integrate1DInternal(aa / scale, bb / scale, cc / scale, degree, 0.0L, 1.0L, computeDelta, fullBoundaryDelta, h, d);
  positive.resize(degree + 1U);
  delta.resize(degree + 1U);
  for (unsigned k = 0; k <= degree; ++k) {
    positive[k] = static_cast<double>(h[k]);
    delta[k] = static_cast<double>(d[k] / scale);
    if (!std::isfinite(positive[k]) || !std::isfinite(delta[k])) {
      throw std::overflow_error("A canonical 1D moment exceeds double precision range.");
    }
  }
}

inline std::vector<long double> SolveInterpolation(std::vector<std::vector<long double>> matrix, std::vector<long double> values) {
  unsigned size = static_cast<unsigned>(values.size());
  for (unsigned column = 0; column < size; ++column) {
    unsigned pivot = column;
    for (unsigned row = column + 1U; row < size; ++row) {
      if (std::abs(matrix[row][column]) > std::abs(matrix[pivot][column])) {
        pivot = row;
      }
    }
    if (matrix[pivot][column] == 0.0L) {
      throw std::invalid_argument("Singular interpolation nodes.");
    }
    std::swap(matrix[pivot], matrix[column]);
    std::swap(values[pivot], values[column]);
    for (unsigned row = column + 1U; row < size; ++row) {
      long double multiplier = matrix[row][column] / matrix[column][column];
      matrix[row][column] = 0.0L;
      for (unsigned k = column + 1U; k < size; ++k) {
        matrix[row][k] -= multiplier * matrix[column][k];
      }
      values[row] -= multiplier * values[column];
    }
  }
  std::vector<long double> solution(size, 0.0L);
  for (unsigned row = size; row > 0; --row) {
    unsigned index = row - 1U;
    long double value = values[index];
    for (unsigned column = index + 1U; column < size; ++column) {
      value -= matrix[index][column] * solution[column];
    }
    solution[index] = value / matrix[index][index];
  }
  return solution;
}

inline void Line3CoefficientsFromNodalValues(const std::vector<double>& nodalValues, double& a, double& b, double& c, double left = 0.0, double right = 1.0) {
  if (nodalValues.size() != 3 || !(right > left) || !std::isfinite(left) || !std::isfinite(right)) {
    throw std::invalid_argument("LINE3 requires values at the left endpoint, right endpoint, and midpoint.");
  }
  for (unsigned k = 0; k < 3; ++k) {
    if (!std::isfinite(nodalValues[k])) {
      throw std::invalid_argument("A nodal value is not finite.");
    }
  }
  long double width = static_cast<long double>(right) - left;
  long double quadratic = 2.0L * (static_cast<long double>(nodalValues[0]) + nodalValues[1] - 2.0L * nodalValues[2]);
  long double linear = 4.0L * nodalValues[2] - 3.0L * nodalValues[0] - nodalValues[1];
  long double aa = quadratic / (width * width);
  long double bb = linear / width - 2.0L * aa * left;
  long double cc = nodalValues[0] - linear * left / width + aa * left * left;
  a = static_cast<double>(aa);
  b = static_cast<double>(bb);
  c = static_cast<double>(cc);
  if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c)) {
    throw std::overflow_error("A LINE3 coefficient exceeds double precision range.");
  }
}

inline void Quad9UnitCoefficientsFromNodalValues(const std::vector<double>& nodalValues, std::vector<std::vector<double>>& coefficient) {
  if (nodalValues.size() != 9) {
    throw std::invalid_argument("QUAD9 requires nine nodal values in FEMuS order.");
  }
  for (unsigned k = 0; k < 9; ++k) {
    if (!std::isfinite(nodalValues[k])) {
      throw std::invalid_argument("A nodal value is not finite.");
    }
  }
  const unsigned node[3][3] = {{0, 3, 7}, {1, 2, 5}, {4, 6, 8}};
  const long double interpolation[3][3] = {{1.0L, 0.0L, 0.0L}, {-3.0L, -1.0L, 4.0L}, {2.0L, 2.0L, -4.0L}};
  long double intermediate[3][3] = {};
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned sampleY = 0; sampleY < 3; ++sampleY) {
      for (unsigned sampleX = 0; sampleX < 3; ++sampleX) {
        intermediate[i][sampleY] += interpolation[i][sampleX] * nodalValues[node[sampleX][sampleY]];
      }
    }
  }
  std::vector<std::vector<double>> result(3, std::vector<double>(3, 0.0));
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      long double value = 0.0L;
      for (unsigned sampleY = 0; sampleY < 3; ++sampleY) {
        value += interpolation[j][sampleY] * intermediate[i][sampleY];
      }
      result[i][j] = static_cast<double>(value);
      if (!std::isfinite(result[i][j])) {
        throw std::overflow_error("A unit-square QUAD9 coefficient exceeds double precision range.");
      }
    }
  }
  coefficient.swap(result);
}

inline void Tri7ReflectedCoefficientsFromNodalValues(const std::vector<double>& nodalValues, std::vector<double>& coefficient) {
  if (nodalValues.size() != 7) {
    throw std::invalid_argument("TRI7 requires seven nodal values in FEMuS order.");
  }
  for (unsigned k = 0; k < 7; ++k) {
    if (!std::isfinite(nodalValues[k])) {
      throw std::invalid_argument("A nodal value is not finite.");
    }
  }
  const long double interpolation[7][7] = {{0.0L, 1.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L}, {-1.0L, -3.0L, 0.0L, 4.0L, 0.0L, 0.0L, 0.0L}, {1.0L, 0.0L, -1.0L, -4.0L, 4.0L, 0.0L, 0.0L}, {2.0L, 2.0L, 0.0L, -4.0L, 0.0L, 0.0L, 0.0L}, {-1.0L, 3.0L, 3.0L, -8.0L, -16.0L, -8.0L, 27.0L}, {-1.0L, -3.0L, -1.0L, 12.0L, 12.0L, 8.0L, -27.0L}, {-3.0L, -3.0L, -3.0L, 12.0L, 12.0L, 12.0L, -27.0L}};
  std::vector<double> result(7, 0.0);
  for (unsigned i = 0; i < 7; ++i) {
    long double value = 0.0L;
    for (unsigned k = 0; k < 7; ++k) {
      value += interpolation[i][k] * nodalValues[k];
    }
    result[i] = static_cast<double>(value);
    if (!std::isfinite(result[i])) {
      throw std::overflow_error("A reflected TRI7 coefficient exceeds double precision range.");
    }
  }
  coefficient.swap(result);
}

inline void Quad9CoefficientsFromNodalValues(const std::vector<double>& nodalValues, std::vector<std::vector<double>>& coefficient, double xLeft = -1.0, double xRight = 1.0, double yBottom = -1.0, double yTop = 1.0) {
  if (nodalValues.size() != 9 || !(xRight > xLeft) || !(yTop > yBottom) || !std::isfinite(xLeft) || !std::isfinite(xRight) || !std::isfinite(yBottom) || !std::isfinite(yTop)) {
    throw std::invalid_argument("Invalid QUAD9 nodal input.");
  }
  std::vector<std::vector<long double>> nodes(9, std::vector<long double>(2, 0.0L));
  nodes[0] = {xLeft, yBottom};
  nodes[1] = {xRight, yBottom};
  nodes[2] = {xRight, yTop};
  nodes[3] = {xLeft, yTop};
  nodes[4] = {0.5L * (xLeft + xRight), yBottom};
  nodes[5] = {xRight, 0.5L * (yBottom + yTop)};
  nodes[6] = {0.5L * (xLeft + xRight), yTop};
  nodes[7] = {xLeft, 0.5L * (yBottom + yTop)};
  nodes[8] = {0.5L * (xLeft + xRight), 0.5L * (yBottom + yTop)};
  std::vector<std::vector<long double>> matrix(9, std::vector<long double>(9, 0.0L));
  std::vector<long double> values(9, 0.0L);
  for (unsigned row = 0; row < 9; ++row) {
    if (!std::isfinite(nodalValues[row])) {
      throw std::invalid_argument("A nodal value is not finite.");
    }
    values[row] = nodalValues[row];
    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j) {
        matrix[row][3U * i + j] = Power(nodes[row][0], i) * Power(nodes[row][1], j);
      }
    }
  }
  std::vector<long double> solution = SolveInterpolation(matrix, values);
  coefficient.assign(3, std::vector<double>(3, 0.0));
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      coefficient[i][j] = static_cast<double>(solution[3U * i + j]);
    }
  }
}

inline void Tri7CoefficientsFromNodalValues(const std::vector<double>& nodalValues, std::vector<double>& coefficient) {
  if (nodalValues.size() != 7) {
    throw std::invalid_argument("TRI7 requires seven nodal values.");
  }
  std::vector<std::vector<long double>> nodes = {{0.0L, 0.0L}, {1.0L, 0.0L}, {0.0L, 1.0L}, {0.5L, 0.0L}, {0.5L, 0.5L}, {0.0L, 0.5L}, {1.0L / 3.0L, 1.0L / 3.0L}};
  std::vector<std::vector<long double>> matrix(7, std::vector<long double>(7, 0.0L));
  std::vector<long double> values(7, 0.0L);
  for (unsigned row = 0; row < 7; ++row) {
    if (!std::isfinite(nodalValues[row])) {
      throw std::invalid_argument("A nodal value is not finite.");
    }
    long double x = nodes[row][0];
    long double y = nodes[row][1];
    matrix[row] = {1.0L, x, y, x * x, x * y, y * y, x * y * (x + y)};
    values[row] = nodalValues[row];
  }
  std::vector<long double> solution = SolveInterpolation(matrix, values);
  coefficient.resize(7);
  for (unsigned k = 0; k < 7; ++k) {
    coefficient[k] = static_cast<double>(solution[k]);
  }
}

}

#ifdef FE_LEVEL_SET_MOMENTS_DEMO
int main() {
  try {
    std::vector<double> linePositive;
    std::vector<double> lineDelta;
    FELevelSetMoments::Integrate1D(1.0, -1.0, 0.16, 4, linePositive, lineDelta);
    std::vector<double> quadNodes = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0};
    std::vector<std::vector<double>> coefficient;
    FELevelSetMoments::Quad9UnitCoefficientsFromNodalValues(quadNodes, coefficient);
    std::vector<std::vector<double>> positive;
    std::vector<std::vector<double>> delta;
    FELevelSetMoments::IntegrateQuad9Unit(coefficient, 4, positive, delta);
    std::cout << std::setprecision(16) << "1D H[0] = " << linePositive[0] << ", delta[0] = " << lineDelta[0] << '\n';
    std::cout << "QUAD9 canonical H[0][0] = " << positive[0][0] << ", delta[0][0] = " << delta[0][0] << '\n';
    std::vector<double> triangleNodes = {-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 0.8};
    std::vector<double> triangleCoefficient;
    FELevelSetMoments::Tri7ReflectedCoefficientsFromNodalValues(triangleNodes, triangleCoefficient);
    FELevelSetMoments::IntegrateTri7Reflected(triangleCoefficient, 4, positive, delta);
    std::cout << "TRI7 H[0][0] = " << positive[0][0] << ", delta[0][0] = " << delta[0][0] << '\n';
  }
  catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
#endif

#endif
