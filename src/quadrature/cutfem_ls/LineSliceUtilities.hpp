#ifndef LINE_SLICE_UTILITIES_HPP
#define LINE_SLICE_UTILITIES_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "ReferenceSurfaceMoments.hpp"

// Shared polynomial and root utilities support the retained LINE integration kernels.

namespace LineSliceUtilities {

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

inline long double ParentCoordinate(long double scale, long double coordinate, long double origin) { return scale == 1 && origin == 0 ? coordinate : std::fma(scale,coordinate,origin); }

inline long double PolynomialValue(const std::vector<long double>& polynomial, long double x) {
  if (polynomial.empty()) return 0.0L;
  long double value = polynomial.back();
  for (unsigned k = static_cast<unsigned>(polynomial.size())-1U; k > 0; --k) {
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
    long double tolerance = 16.0L*std::numeric_limits<long double>::epsilon()*std::max(std::abs(values[k]),count == 0 ? 0.0L : std::abs(values[count-1U]));
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

template <class Type> inline void NormalizeCanonicalPolynomial(const std::vector<std::vector<Type>>& coefficient, std::vector<std::vector<long double>>& normalized, long double& scale) {
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
  int exponent = 0;
  std::frexp(scale,&exponent);
  scale = std::ldexp(1.0L,exponent-1);
  normalized.assign(3, std::vector<long double>(3, 0.0L));
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      normalized[i][j] = coefficient[i][j] / scale;
    }
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

inline unsigned MomentIndex(unsigned degree, unsigned i, unsigned j) {
  return i * (2U * degree + 3U - i) / 2U + j;
}

inline std::vector<long double> DiscriminantPolynomial(const std::vector<std::vector<long double>>& field) {
  return PolynomialSum(PolynomialProduct(field[1], field[1]), PolynomialProduct(field[0], field[2]), -4.0L);
}

}

#endif
