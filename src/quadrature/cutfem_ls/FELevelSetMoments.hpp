#ifndef FE_LEVEL_SET_MOMENTS_HPP
#define FE_LEVEL_SET_MOMENTS_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
#include "ReferenceSurfaceMoments.hpp"
#include "ReferenceElementLS.hpp"

// Normalized coordinates below are local integration charts; nodal entries use FEMuS reference elements.

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

  inline std::vector<long double> SeriesConstant(long double value, unsigned order) {
    std::vector<long double> result(order + 1U, 0.0L);
    result[0] = value;
    return result;
  }

  inline std::vector<long double> SeriesProduct(const std::vector<long double>& first, const std::vector<long double>& second) {
    std::vector<long double> result(first.size(), 0.0L);
    unsigned firstDegree = static_cast<unsigned>(first.size()) - 1U;
    unsigned secondDegree = static_cast<unsigned>(second.size()) - 1U;
    while(firstDegree != 0 && first[firstDegree] == 0.0L) --firstDegree;
    while(secondDegree != 0 && second[secondDegree] == 0.0L) --secondDegree;
    unsigned degree = std::min(static_cast<unsigned>(result.size()) - 1U, firstDegree + secondDegree);
    for (unsigned k = 0; k <= degree; ++k) {
      unsigned begin = k > secondDegree ? k - secondDegree : 0;
      unsigned end = std::min(k, firstDegree);
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
    for (unsigned k = static_cast<unsigned>(polynomial.size()) - 1U; k > 0; --k) {
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

  inline std::vector<long double> SurfaceGradientSeriesNorm(const std::vector<std::vector<std::vector<long double>>>& gradient, const std::vector<long double>& first, const std::vector<long double>& second) {
    const unsigned order = static_cast<unsigned>(first.size())-1U;
    std::vector<long double> squared(order+1U,0.0L);
    for(unsigned component = 0; component < gradient.size(); ++component) {
      std::vector<long double> value(order+1U,0.0L);
      for(unsigned i = static_cast<unsigned>(gradient[component].size()); i > 0; --i) {
        value = SeriesProduct(value,first);
        const std::vector<long double> row = SeriesPolynomial(gradient[component][i-1U],second);
        for(unsigned k = 0; k <= order; ++k) value[k] += row[k];
      }
      const std::vector<long double> term = SeriesProduct(value,value);
      for(unsigned k = 0; k <= order; ++k) squared[k] += term[k];
    }
    return SeriesSquareRoot(squared);
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

  inline void SeriesProductIntegral(const std::vector<long double>& first, const std::vector<long double>& second, long double& value, long double& error, bool& converged) {
    value = 0.0L;
    long double correction = 0.0L;
    long double norm = 0.0L;
    long double tail = 0.0L;
    long double preceding = 0.0L;
    unsigned order = static_cast<unsigned>(first.size()) - 1U;
    unsigned block = std::max(8U, order / 6U);
    unsigned firstDegree = order;
    unsigned secondDegree = order;
    while(firstDegree != 0 && first[firstDegree] == 0.0L) --firstDegree;
    while(secondDegree != 0 && second[secondDegree] == 0.0L) --secondDegree;
    for(unsigned k = 0; k <= order; k += 2U) {
      long double coefficient = 0.0L;
      unsigned begin = k > secondDegree ? k - secondDegree : 0;
      unsigned end = std::min(k, firstDegree);
      for(unsigned j = begin; j <= end; ++j) coefficient += first[j] * second[k - j];
      long double term = 2.0L * coefficient / (k + 1U);
      long double corrected = term - correction;
      long double next = value + corrected;
      correction = (next - value) - corrected;
      value = next;
      norm += std::abs(term);
      if(k + block > order) tail += std::abs(term);
      else if(k + 2U * block > order) preceding += std::abs(term);
    }
    long double rounding = 64.0L * std::numeric_limits<long double>::epsilon() * norm;
    error = 8.0L * tail + rounding;
    converged = std::isfinite(value) && std::isfinite(error) && (tail <= 0.5L * preceding || tail <= rounding);
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

  inline bool SliceTaylor(const std::vector<std::vector<long double>>& field, const std::vector<long double>& reducedDiscriminant, unsigned leftMultiplicity, unsigned rightMultiplicity, long double xLeft, long double xRight, long double thetaLeft, long double thetaRight, const std::vector<int>& signs, const std::vector<int>& boundary, const std::vector<long double>& positiveSegments, unsigned degree, unsigned order, unsigned triangle, const std::vector<long double>& domain, bool computePositive, bool computeDelta, bool fullBoundaryDelta, std::vector<long double>& values, std::vector<long double>& errors, bool fullDeltaChart = false, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr, const std::vector<unsigned char>* active = nullptr) {
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
      reducedSquareRoot = SeriesSquareRoot(SeriesPolynomial(reducedDiscriminant, x));
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
    std::vector<long double> da, db, dc, commonDeltaWeight;
    if(computeDelta && !fullDeltaChart) {
      da = SeriesPolynomial(PolynomialDerivative(field[0]), x);
      db = SeriesPolynomial(PolynomialDerivative(field[1]), x);
      dc = SeriesPolynomial(PolynomialDerivative(field[2]), x);
    }
    if(computeDelta && fullDeltaChart) {
      if(quadratic && leftMultiplicity <= 1 && rightMultiplicity <= 1) {
        commonDeltaWeight = SeriesConstant(2.0L * width * halfWidth / std::pow(width, 0.5L * (leftMultiplicity + rightMultiplicity)), order);
        if(leftMultiplicity == 0) commonDeltaWeight = SeriesProduct(commonDeltaWeight, sine);
        if(rightMultiplicity == 0) commonDeltaWeight = SeriesProduct(commonDeltaWeight, cosine);
        commonDeltaWeight = SeriesQuotient(commonDeltaWeight, reducedSquareRoot);
      }
      else {
        std::vector<long double> derivative = quadratic ? squareRoot : b;
        if(!quadratic && derivative[0] < 0.0L) {
          for(unsigned k = 0; k <= order; ++k) derivative[k] = -derivative[k];
        }
        commonDeltaWeight = SeriesQuotient(jacobian, derivative);
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
          if(surfaceGradient != nullptr) weight = SeriesProduct(weight,SurfaceGradientSeriesNorm(*surfaceGradient,x,root));
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
        if(surfaceGradient != nullptr) weight = SeriesProduct(weight,SurfaceGradientSeriesNorm(*surfaceGradient,x,root));
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
    const unsigned count = (degree + 1U) * (degree + 2U) / 2U;
    std::vector<unsigned char> positiveNeeded(degree + 1U, 0);
    std::vector<unsigned char> deltaNeeded(degree + 1U, 0);
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned index = MomentIndex(degree,i,j);
        positiveNeeded[j] = positiveNeeded[j] || (computePositive && (active == nullptr || (*active)[index] != 0));
        deltaNeeded[j] = deltaNeeded[j] || (computeDelta && (active == nullptr || (*active)[count+index] != 0));
      }
    }
    std::vector<std::vector<long double>> positivePowers(cuts.size(), SeriesConstant(1.0L, order));
    std::vector<std::vector<long double>> rootPowers(roots.size(), SeriesConstant(1.0L, order));
    std::vector<std::vector<long double>> positiveSlice(degree + 1U, SeriesConstant(0.0L, order));
    std::vector<std::vector<long double>> deltaSlice(degree + 1U, SeriesConstant(0.0L, order));
    for (unsigned j = 0; j <= degree; ++j) {
      for (unsigned r = 0; r < cuts.size(); ++r) positivePowers[r] = SeriesProduct(positivePowers[r], cuts[r]);
      if (positiveNeeded[j] != 0) {
        for (unsigned r = 0; r < positiveSegments.size(); ++r) {
          if (positiveSegments[r] != 0.0L) {
            for (unsigned k = 0; k <= order; ++k) {
              positiveSlice[j][k] += (positivePowers[r + 1U][k] - positivePowers[r][k]) / static_cast<long double>(j + 1U);
            }
          }
        }
        positiveSlice[j] = SeriesProduct(positiveSlice[j], jacobian);
      }
      for (unsigned r = 0; r < roots.size(); ++r) {
        if (deltaNeeded[j] != 0) {
          std::vector<long double> term = SeriesProduct(rootPowers[r], weights[r]);
          for (unsigned k = 0; k <= order; ++k) deltaSlice[j][k] += term[k];
        }
        rootPowers[r] = SeriesProduct(rootPowers[r], roots[r]);
      }
    }
    values.assign(2U * count, 0.0L);
    errors.assign(2U * count, 0.0L);
    bool allConverged = true;
    std::vector<long double> xPower = SeriesConstant(1.0L, order);
    for (unsigned i = 0; i <= degree; ++i) {
      for (unsigned j = 0; j <= degree - i; ++j) {
        unsigned index = MomentIndex(degree, i, j);
        bool converged = true;
        if (computePositive && (active == nullptr || (*active)[index] != 0)) {
          SeriesProductIntegral(xPower, positiveSlice[j], values[index], errors[index], converged);
          if(!converged) errors[index] = std::numeric_limits<long double>::infinity();
          allConverged = allConverged && converged;
        }
        converged = true;
        if (computeDelta && (active == nullptr || (*active)[count + index] != 0)) {
          SeriesProductIntegral(xPower, deltaSlice[j], values[count + index], errors[count + index], converged);
          if(!converged) errors[count + index] = std::numeric_limits<long double>::infinity();
          allConverged = allConverged && converged;
        }
      }
      xPower = SeriesProduct(xPower, x);
    }
    return allConverged;
  }

  inline void IntegrateTaylorInterval(const std::vector<std::vector<long double>>& field, const std::vector<long double>& reducedDiscriminant, unsigned leftMultiplicity, unsigned rightMultiplicity, long double xLeft, long double xRight, long double thetaLeft, long double thetaRight, const std::vector<int>& signs, const std::vector<int>& boundary, const std::vector<long double>& positiveSegments, unsigned degree, unsigned order, unsigned triangle, const std::vector<long double>& domain, bool computePositive, bool computeDelta, bool fullBoundaryDelta, long double absolutePositive, long double absoluteDelta, long double relativeTolerance, unsigned depth, unsigned maximumDepth, std::vector<long double>& sum, std::vector<long double>& sumError, unsigned& acceptedIntervals, bool fullDeltaChart = false, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr, const std::vector<unsigned char>* active = nullptr) {
    std::vector<long double> values;
    std::vector<long double> errors;
    SliceTaylor(field, reducedDiscriminant, leftMultiplicity, rightMultiplicity, xLeft, xRight, thetaLeft, thetaRight, signs, boundary, positiveSegments, degree, order, triangle, domain, computePositive, computeDelta, fullBoundaryDelta, values, errors, fullDeltaChart, surfaceGradient, active);
    const unsigned count = (degree + 1U) * (degree + 2U) / 2U;
    std::vector<unsigned char> failed(values.size(), 0);
    bool anyFailed = false;
    for (unsigned k = 0; k < values.size(); ++k) {
      if(active != nullptr && (*active)[k] == 0) continue;
      const long double tolerance = ((k < count) ? absolutePositive : absoluteDelta) + relativeTolerance * std::abs(values[k]);
      if(errors[k] <= tolerance) {
        sum[k] += values[k];
        sumError[k] += errors[k];
      }
      else {
        failed[k] = 1;
        anyFailed = true;
      }
    }
    if (!anyFailed) {
      ++acceptedIntervals;
      return;
    }
    long double middle = 0.5L * (thetaLeft + thetaRight);
    if (depth >= maximumDepth || middle == thetaLeft || middle == thetaRight) {
      throw std::runtime_error("Taylor recurrence did not reach the requested tolerance; increase maximumDepth or relax the tolerance.");
    }
    IntegrateTaylorInterval(field, reducedDiscriminant, leftMultiplicity, rightMultiplicity, xLeft, xRight, thetaLeft, middle, signs, boundary, positiveSegments, degree, order, triangle, domain, computePositive, computeDelta, fullBoundaryDelta, 0.5L * absolutePositive, 0.5L * absoluteDelta, relativeTolerance, depth + 1U, maximumDepth, sum, sumError, acceptedIntervals, fullDeltaChart, surfaceGradient, &failed);
    IntegrateTaylorInterval(field, reducedDiscriminant, leftMultiplicity, rightMultiplicity, xLeft, xRight, middle, thetaRight, signs, boundary, positiveSegments, degree, order, triangle, domain, computePositive, computeDelta, fullBoundaryDelta, 0.5L * absolutePositive, 0.5L * absoluteDelta, relativeTolerance, depth + 1U, maximumDepth, sum, sumError, acceptedIntervals, fullDeltaChart, surfaceGradient, &failed);
  }

  inline void IntegrateDirection(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned triangle, const std::vector<long double>& domain, bool computePositive, bool computeDelta, bool fullBoundaryDelta, long double absoluteTolerance, long double deltaTolerance, long double relativeTolerance, unsigned order, unsigned maximumDepth, std::vector<long double>& result, std::vector<long double>& resultError, unsigned& acceptedIntervals, bool fullDeltaChart = false, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
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
      IntegrateTaylorInterval(field, reduced, leftMultiplicity, rightMultiplicity, left, right, 0.0L, 0.5L * std::acos(-1.0L), signs, boundary, positiveSegments, degree, order, triangle, domain, computePositive, computeDelta, fullBoundaryDelta, absoluteTolerance * fraction, deltaTolerance * fraction, relativeTolerance, 0, maximumDepth, result, resultError, acceptedIntervals, fullDeltaChart, surfaceGradient);
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
    std::vector<std::vector<long double>> fibers(3, std::vector<long double>(3, 0.0L));
    long double left = domain[order == 1 ? 0 : 2];
    long double right = domain[order == 1 ? 1 : 3];
    std::vector<long double> candidates = {left, right};
    bool selected = false;
    long double norm = 0.0L;
    long double radius = std::max(1.0L, std::max(std::abs(left), std::abs(right)));
    for(unsigned j = 0; j < 3; ++j) {
      for(unsigned i = 0; i < 3; ++i) fibers[j][i] = order == 1 ? coefficient[i][j] : coefficient[j][i];
      Trim(fibers[j]);
      norm += PolynomialBound(fibers[j], radius);
      if(PolynomialNorm(fibers[j]) != 0.0L) {
        PolynomialRoots(fibers[j], left, right, candidates);
        selected = true;
      }
    }
    if(!selected) return false;
    long double tolerance = 32.0L * std::sqrt(std::numeric_limits<long double>::epsilon()) * norm;
    for(unsigned k = 0; k < candidates.size(); ++k) {
      bool common = true;
      for(unsigned j = 0; j < 3; ++j) common = common && std::abs(PolynomialValue(fibers[j], candidates[k])) <= tolerance;
      if(common) return false;
    }
    return true;
  }

  inline void IntegrateNormalizedCell(const std::vector<std::vector<long double>>& normalized, unsigned degree, unsigned triangle, const std::vector<long double>& domain, bool computePositive, bool computeDelta, bool fullBoundaryDelta, long double absolutePositive, long double absoluteDelta, long double relativeTolerance, unsigned seriesOrder, unsigned maximumDepth, unsigned integrationOrder, std::vector<long double>& result, std::vector<long double>& errors, unsigned& acceptedIntervals, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
    if (integrationOrder > 2U) {
      throw std::invalid_argument("integrationOrder must be 0 (automatic), 1 (y first), or 2 (x first).");
    }
    unsigned count = (degree + 1U) * (degree + 2U) / 2U;
    bool reversePositive = computePositive && ((integrationOrder == 2U) || (integrationOrder == 0U && PreferReverseOrder(normalized, domain)));
    result.assign(2U * count, 0.0L);
    errors.assign(2U * count, 0.0L);
    if(computeDelta) {
      unsigned intervalsBefore = acceptedIntervals;
      unsigned selected = integrationOrder == 0 ? (PreferReverseOrder(normalized, domain) ? 2U : 1U) : integrationOrder;
      if(!RegularDeltaProjection(normalized, domain, selected)) {
        if(integrationOrder == 0 && RegularDeltaProjection(normalized, domain, 3U - selected)) selected = 3U - selected;
        else selected = 0;
      }
      try {
        if(selected == 1) {
          IntegrateDirection(normalized, degree, triangle, domain, computePositive, true, fullBoundaryDelta, absolutePositive, absoluteDelta, relativeTolerance, seriesOrder, maximumDepth, result, errors, acceptedIntervals, true, surfaceGradient);
          return;
        }
        if(selected == 2) {
          std::vector<std::vector<long double>> transposed(3, std::vector<long double>(3));
          for(unsigned i = 0; i < 3; ++i) {
            for(unsigned j = 0; j < 3; ++j) transposed[i][j] = normalized[j][i];
          }
          std::vector<long double> swapped(domain);
          std::swap(swapped[0], swapped[2]);
          std::swap(swapped[1], swapped[3]);
          unsigned shape = triangle == 2 ? 3U : (triangle == 3 ? 2U : triangle);
          std::vector<long double> value, error;
          std::vector<std::vector<std::vector<long double>>> transposedGradient;
          if(surfaceGradient != nullptr) SurfaceGradientTranspose(*surfaceGradient,transposedGradient);
          IntegrateDirection(transposed, degree, shape, swapped, computePositive, true, fullBoundaryDelta, absolutePositive, absoluteDelta, relativeTolerance, seriesOrder, maximumDepth, value, error, acceptedIntervals, true, surfaceGradient == nullptr ? nullptr : &transposedGradient);
          for(unsigned i = 0; i <= degree; ++i) {
            for(unsigned j = 0; j <= degree - i; ++j) {
              unsigned target = MomentIndex(degree, i, j);
              unsigned source = MomentIndex(degree, j, i);
              result[target] = value[source];
              result[count + target] = value[count + source];
              errors[target] = error[source];
              errors[count + target] = error[count + source];
            }
          }
          return;
        }
      }
      catch(const std::runtime_error&) {
        std::fill(result.begin(), result.end(), 0.0L);
        std::fill(errors.begin(), errors.end(), 0.0L);
        acceptedIntervals = intervalsBefore;
      }
    }
    if ((computePositive && !reversePositive) || computeDelta) {
      IntegrateDirection(normalized, degree, triangle, domain, computePositive && !reversePositive, computeDelta, fullBoundaryDelta, absolutePositive, 0.5L * absoluteDelta, relativeTolerance, seriesOrder, maximumDepth, result, errors, acceptedIntervals, false, surfaceGradient);
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
      std::vector<std::vector<std::vector<long double>>> transposedGradient;
      if(surfaceGradient != nullptr) SurfaceGradientTranspose(*surfaceGradient,transposedGradient);
      IntegrateDirection(transposed, degree, swappedTriangle, swappedDomain, reversePositive, computeDelta, fullBoundaryDelta, absolutePositive, 0.5L * absoluteDelta, relativeTolerance, seriesOrder, maximumDepth, reverse, reverseError, acceptedIntervals, false, surfaceGradient == nullptr ? nullptr : &transposedGradient);
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

  inline void CanonicalBernstein(const std::vector<std::vector<long double>>& coefficient, long double bernstein[3][3]) {
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) {
        bernstein[i][j] = 0.0L;
        for(unsigned r = 0; r <= i; ++r) {
          for(unsigned s = 0; s <= j; ++s) bernstein[i][j] += coefficient[r][s] * Binomial(i, r) * Binomial(j, s) / (Binomial(2, r) * Binomial(2, s));
        }
      }
    }
  }

  inline unsigned CanonicalGraphOrder(const std::vector<std::vector<long double>>& coefficient, unsigned requestedOrder, int& phase, long double* derivativeBound = nullptr) {
    long double b[3][3];
    CanonicalBernstein(coefficient, b);
    long double norm = 0.0L;
    long double minimum = b[0][0];
    long double maximum = minimum;
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) {
        norm += std::abs(coefficient[i][j]);
        minimum = std::min(minimum, b[i][j]);
        maximum = std::max(maximum, b[i][j]);
      }
    }
    long double padding = 512.0L * std::numeric_limits<long double>::epsilon() * norm;
    long double conditioning = 32.0L * std::sqrt(std::numeric_limits<long double>::epsilon()) * norm;
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
        const long double derivative[3] = {2.0L * (middle - low), high - low, 2.0L * (high - middle)};
        for(unsigned j = 0; j < 3; ++j) {
          derivativeMinimum = std::min(derivativeMinimum, derivative[j]);
          derivativeMaximum = std::max(derivativeMaximum, derivative[j]);
        }
      }
      long double bound = 0.0L;
      if(derivativeMinimum > conditioning) bound = derivativeMinimum - padding;
      if(derivativeMaximum < -conditioning) bound = -derivativeMaximum - padding;
      if(bound > selectedBound) {
        selectedBound = bound;
        selected = order;
      }
    }
    if(derivativeBound != nullptr) *derivativeBound = selectedBound;
    return selected;
  }


  inline void CanonicalGraphSample(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned order, long double point, bool computeDelta, std::vector<long double>& values, bool triangle = false, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr, const std::vector<unsigned char>* active = nullptr) {
  long double polynomial[3];
  for(unsigned k = 0; k < 3; ++k) {
    const long double quadratic = order == 1 ? coefficient[2][k] : coefficient[k][2];
    const long double linear = order == 1 ? coefficient[1][k] : coefficient[k][1];
    const long double constant = order == 1 ? coefficient[0][k] : coefficient[k][0];
    if(quadratic == 0.0L) polynomial[k] = linear == 0.0L ? constant : std::fma(linear, point, constant);
    else polynomial[k] = std::fma(std::fma(quadratic, point, linear), point, constant);
  }

  long double roots[2];
    unsigned numberOfRoots = 0;
    if(polynomial[2] == 0.0L) {
      if(polynomial[1] != 0.0L) roots[numberOfRoots++] = -polynomial[0] / polynomial[1];
    }
    else {
      const long double discriminant = std::fma(-4.0L * polynomial[2], polynomial[0], polynomial[1] * polynomial[1]);
      if(discriminant == 0.0L) roots[numberOfRoots++] = -polynomial[1] / (2.0L * polynomial[2]);
      else if(discriminant > 0.0L) {
        const long double numerator = -0.5L * (polynomial[1] + std::copysign(std::sqrt(discriminant), polynomial[1]));
        roots[numberOfRoots++] = numerator / polynomial[2];
        roots[numberOfRoots++] = polynomial[0] / numerator;
      }
    }

    const long double lower = triangle && order == 2 ? point : 0.0L;
    const long double upper = triangle && order == 1 ? point : 1.0L;
    long double root = -1.0L;
    for(unsigned k = 0; k < numberOfRoots; ++k) {
      if(roots[k] > lower && roots[k] < upper) {
        if(root >= 0.0L) throw std::runtime_error("A monotone graph has more than one interior root.");
        root = roots[k];
      }
    }
    if(root < 0.0L) throw std::runtime_error("A monotone graph root was not resolved.");

    const long double signedDerivative = polynomial[2] == 0.0L ? polynomial[1] : std::fma(2.0L * polynomial[2], root, polynomial[1]);
    const long double derivative = std::abs(signedDerivative);
    if(!(derivative > 0.0L)) throw std::runtime_error("A monotone graph derivative was not resolved.");

    const long double surfaceNorm = computeDelta && surfaceGradient != nullptr ? SurfaceGradientNorm(*surfaceGradient,order == 1 ? point : root,order == 1 ? root : point) : 1.0L;
    const bool above = signedDerivative > 0.0L;
    const long double left = above ? root : lower;
    const long double right = above ? upper : root;
    const long double width = right - left;
    const unsigned count = (degree + 1U) * (degree + 2U) / 2U;
    values.assign(2U * count, 0.0L);

    long double rootPower = 1.0L;
    long double rightPower = 1.0L;
    long double geometric = 0.0L;
    for(unsigned j = 0; j <= degree; ++j) {
      geometric = rightPower + left * geometric;
      const unsigned index = order == 1 ? MomentIndex(degree, 0, j) : MomentIndex(degree, j, 0);
      values[index] = width * geometric / (j + 1U);
      values[count + index] = computeDelta ? surfaceNorm * rootPower / derivative : 0.0L;
      if(computeDelta) rootPower *= root;
      rightPower *= right;
    }

    long double pointPower = 1.0L;
    for(unsigned i = 1; i <= degree; ++i) {
      pointPower *= point;
      for(unsigned j = 0; j <= degree - i; ++j) {
        const unsigned source = order == 1 ? MomentIndex(degree, 0, j) : MomentIndex(degree, j, 0);
        const unsigned target = order == 1 ? MomentIndex(degree, i, j) : MomentIndex(degree, j, i);
        if(active != nullptr && (*active)[target] == 0 && (*active)[count + target] == 0) continue;
        values[target] = pointPower * values[source];
        values[count + target] = computeDelta ? pointPower * values[count + source] : 0.0L;
      }
    }
  }



  inline void IntegrateCanonicalGraphInterval(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned order, long double left, long double right, bool computeDelta, long double absoluteH, long double absoluteDelta, long double relativeTolerance, unsigned depth, unsigned maximumDepth, std::vector<long double>& total, std::vector<long double>& totalError, unsigned& accepted, bool triangle = false, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr, const std::vector<unsigned char>* active = nullptr) {
    static const long double nodes[8] = {0.991455371120812639206854697526329L, 0.949107912342758524526189684047851L, 0.864864423359769072789712788640926L, 0.741531185599394439863864773280788L, 0.586087235467691130294144838258730L, 0.405845151377397166906606412076961L, 0.207784955007898467600689403773245L, 0.0L};
    static const long double kronrod[8] = {0.022935322010529224963732008058970L, 0.063092092629978553290700663189204L, 0.104790010322250183839876322541518L, 0.140653259715525918745189590510238L, 0.169004726639267902826583426598550L, 0.190350578064785409913256402421014L, 0.204432940075298892414161999234649L, 0.209482141084727828012999174891714L};
    static const long double gauss[4] = {0.129484966168869693270611432679082L, 0.279705391489276667901467771423780L, 0.381830050505118944950369775488975L, 0.417959183673469387755102040816327L};

    const unsigned count = (degree + 1U) * (degree + 2U) / 2U;
    const unsigned size = 2U * count;
    const long double width = right - left;
    const long double middle = 0.5L * (left + right);
    const long double halfWidth = 0.5L * width;

    std::vector<long double> high(size, 0.0L);
    std::vector<long double> low(size, 0.0L);
    std::vector<long double> absolute(size, 0.0L);
    std::vector<long double> samples[15];

    for(unsigned node = 0; node < 8; ++node) {
      const unsigned sides = node == 7 ? 1U : 2U;
      for(unsigned side = 0; side < sides; ++side) {
        const unsigned sample = 2U * node + side;
        const long double point = middle + (side == 0 ? -1.0L : 1.0L) * halfWidth * nodes[node];

        CanonicalGraphSample(coefficient, degree, order, point, computeDelta, samples[sample], triangle, surfaceGradient, active);

        for(unsigned k = 0; k < size; ++k) {
          if(active != nullptr && (*active)[k] == 0) continue;
          high[k] += kronrod[node] * samples[sample][k];
          absolute[k] += kronrod[node] * std::abs(samples[sample][k]);
          if((node & 1U) != 0) low[k] += gauss[node / 2U] * samples[sample][k];
        }
      }
    }

    std::vector<long double> error(size, 0.0L);
    std::vector<unsigned char> failed(size, 0);
    bool anyFailed = false;

    for(unsigned k = 0; k < size; ++k) {
      if(active != nullptr && (*active)[k] == 0) continue;
      const long double mean = 0.5L * high[k];
      long double variation = 0.0L;

      for(unsigned sample = 0; sample < 15; ++sample) {
        variation += kronrod[sample / 2U] * std::abs(samples[sample][k] - mean);
      }

      long double estimate = std::abs(high[k] - low[k]);

      if(variation > 0.0L && estimate > 0.0L) {
        const long double ratio = 200.0L * estimate / variation;
        estimate = ratio >= 1.0L ? variation : variation * ratio * std::sqrt(ratio);
      }

      estimate = std::max(estimate, 128.0L * std::numeric_limits<long double>::epsilon() * absolute[k]);
      error[k] = halfWidth * estimate;
      high[k] *= halfWidth;

      const long double tolerance = (k < count ? absoluteH : absoluteDelta) * width + relativeTolerance * std::abs(high[k]);
      if(error[k] <= tolerance) {
        total[k] += high[k];
        totalError[k] += error[k];
      }
      else {
        failed[k] = 1;
        anyFailed = true;
      }
    }

    if(!anyFailed) {
      ++accepted;
      return;
    }

    if(depth >= maximumDepth || middle == left || middle == right) throw std::runtime_error("Monotone graph integration did not reach the requested tolerance.");

    IntegrateCanonicalGraphInterval(coefficient, degree, order, left, middle, computeDelta, absoluteH, absoluteDelta, relativeTolerance, depth + 1U, maximumDepth, total, totalError, accepted, triangle, surfaceGradient, &failed);
    IntegrateCanonicalGraphInterval(coefficient, degree, order, middle, right, computeDelta, absoluteH, absoluteDelta, relativeTolerance, depth + 1U, maximumDepth, total, totalError, accepted, triangle, surfaceGradient, &failed);
  }



  inline bool IntegrateCanonicalGraph(const std::vector<std::vector<long double>>& coefficient, unsigned degree, unsigned order, bool triangle, bool computeDelta, long double derivativeBound, long double absoluteH, long double absoluteDelta, long double relativeTolerance, unsigned maximumDepth, std::vector<long double>& total, std::vector<long double>& totalError, unsigned& accepted, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
    const std::vector<long double> domain = {0.0L, 1.0L, 0.0L, 1.0L};
    std::vector<std::vector<long double>> field(3, std::vector<long double>(3));
    long double norm = 0.0L;
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) {
        field[2U - j][i] = order == 1 ? coefficient[i][j] : coefficient[j][i];
        norm += std::abs(coefficient[i][j]);
      }
    }
    unsigned shape = triangle ? (order == 1 ? 2U : 3U) : 0U;
    std::vector<std::vector<long double>> edges = {BoundaryPolynomial(field, false, shape, domain), BoundaryPolynomial(field, true, shape, domain)};
    std::vector<long double> events = {0.0L, 1.0L};
    long double rootError = 0.0L;
    long double rounding = 512.0L * std::numeric_limits<long double>::epsilon();
    for(unsigned side = 0; side < 2; ++side) {
      Trim(edges[side]);
      if(edges[side].size() > 3 || PolynomialNorm(edges[side]) <= rounding * norm) return false;
      if(edges[side].size() == 3) {
        long double a = edges[side][2], b = edges[side][1], c = edges[side][0];
        long double discriminant = std::fma(-4.0L * a, c, b * b);
        if(std::abs(discriminant) <= rounding * (b * b + 4.0L * std::abs(a * c))) return false;
      }
      std::vector<long double> roots;
      PolynomialRoots(edges[side], 0.0L, 1.0L, roots);
      std::vector<long double> derivative = PolynomialDerivative(edges[side]);
      for(unsigned k = 0; k < roots.size(); ++k) {
        long double error = std::abs(PolynomialValue(edges[side], roots[k])) + rounding * PolynomialBound(edges[side], roots[k]);
        long double slope = std::abs(PolynomialValue(derivative, roots[k])) - rounding * PolynomialBound(derivative, roots[k]);
        long double quadratic = edges[side].size() == 3 ? std::abs(edges[side][2]) : 0.0L;
        if(!(slope > 0.0L) || slope * slope <= 8.0L * quadratic * error) return false;
        rootError += 2.0L * error / slope;
        events.push_back(roots[k]);
      }
    }
    const long double surfaceBound = surfaceGradient == nullptr ? 1.0L : SurfaceGradientBound(*surfaceGradient);
    if(rootError > 0.125L * absoluteH || (computeDelta && rootError * surfaceBound / derivativeBound > 0.125L * absoluteDelta)) return false;
    SortDistinct(events);
    unsigned count = (degree + 1U) * (degree + 2U) / 2U;
    for(unsigned interval = 1; interval < events.size(); ++interval) {
      long double left = events[interval - 1U], right = events[interval];
      if(!(right > left)) continue;
      long double middle = 0.5L * (left + right);
      long double bottom = PolynomialValue(edges[0], middle);
      long double top = PolynomialValue(edges[1], middle);
      if(bottom == 0.0L || top == 0.0L) return false;
      if((bottom > 0.0L) == (top > 0.0L)) {
        if(bottom > 0.0L) {
          for(unsigned i = 0; i <= degree; ++i) {
            for(unsigned j = 0; j <= degree - i; ++j) {
              unsigned target = order == 1 ? MomentIndex(degree, i, j) : MomentIndex(degree, j, i);
              total[target] += FullStripMoment(left, right, i, j, shape, domain);
            }
          }
        }
        ++accepted;
      }
      else IntegrateCanonicalGraphInterval(coefficient, degree, order, left, right, computeDelta, 0.875L * absoluteH, 0.875L * absoluteDelta, relativeTolerance, 0, maximumDepth, total, totalError, accepted, triangle, surfaceGradient);
    }
    for(unsigned k = 0; k < count; ++k) {
      totalError[k] += rootError;
      if(computeDelta) totalError[count + k] += rootError * surfaceBound / derivativeBound;
    }
    return true;
  }

  inline void Line3ReferenceCoefficientsFromNodalValues(const std::vector<double>& nodalValues, double& a, double& b, double& c) {
    if(nodalValues.size() != 3) throw std::invalid_argument("LINE3 requires endpoint and midpoint values in FEMuS order.");
    for(unsigned i = 0; i < 3; ++i) if(!std::isfinite(nodalValues[i])) throw std::invalid_argument("A nodal value is not finite.");
    a = static_cast<double>(0.5L*(static_cast<long double>(nodalValues[0])+nodalValues[1])-nodalValues[2]);
    b = static_cast<double>(0.5L*(static_cast<long double>(nodalValues[1])-nodalValues[0]));
    c = nodalValues[2];
    if(!std::isfinite(a) || !std::isfinite(b)) throw std::overflow_error("A LINE3 coefficient exceeds double range.");
  }

}

#endif
