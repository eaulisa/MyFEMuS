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




}

#endif
