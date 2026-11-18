#ifndef FE_LEVEL_SET_MOMENTS_HPP
#define FE_LEVEL_SET_MOMENTS_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
#include "ReferenceSurfaceMoments.hpp"
#include "QuadratureWarnings.hpp"
#include "ReferenceElementLS.hpp"

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

}

#endif
