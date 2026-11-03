#ifndef ADAPTIVE_LINE_SLICE_MOMENTS_HPP
#define ADAPTIVE_LINE_SLICE_MOMENTS_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Paired LINE3/QUAD9/TRI7 H/delta moments use analytic segment formulas and adaptive Gauss-Kronrod integration without Taylor integration or divergence recurrence.

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

inline void FactorEndpoint(std::vector<long double>& polynomial, long double root, unsigned& multiplicity) {
  multiplicity = 0;
  while (polynomial.size() > 1 && NearZero(polynomial, root)) {
    polynomial = DivideLinear(polynomial, root);
    ++multiplicity;
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

inline long double CanonicalFullMoment(unsigned i, unsigned j, bool triangle) {
  return triangle ? 1.0L / (static_cast<long double>(j + 1U) * (i + j + 2U)) : 1.0L / (static_cast<long double>(i + 1U) * (j + 1U));
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

inline void RationalSegmentMoments(long double a, long double b, long double c, long double left, long double right, unsigned degree, long double tolerance, unsigned depth, std::vector<long double>& values, std::vector<long double>& errors) {
  values.assign(degree+1U,0.0L);
  errors.assign(degree+1U,0.0L);
  if(!(right > left)) return;
  long double middle = 0.5L*(left+right);
  long double halfWidth = 0.5L*(right-left);
  long double constant = std::fma(a,middle*middle,std::fma(b,middle,c));
  if(!(constant > 0.0L) || !std::isfinite(constant)) throw std::domain_error("A rational segment has a nonregular denominator.");
  long double quadratic = a*halfWidth*halfWidth/constant;
  long double linear = (2.0L*a*middle+b)*halfWidth/constant;
  long double ratio = std::abs(quadratic)+std::abs(linear);
  if(ratio > 0.5L) {
    if(depth >= 48) throw std::runtime_error("A rational segment did not reach its analytic convergence radius.");
    std::vector<long double> other;
    std::vector<long double> otherError;
    RationalSegmentMoments(a,b,c,left,middle,degree,0.5L*tolerance,depth+1U,values,errors);
    RationalSegmentMoments(a,b,c,middle,right,degree,0.5L*tolerance,depth+1U,other,otherError);
    for(unsigned j = 0; j <= degree; ++j) {
      values[j] += other[j];
      errors[j] += otherError[j];
    }
    return;
  }
  std::vector<std::vector<long double>> monomial(degree+1U);
  for(unsigned j = 0; j <= degree; ++j) {
    monomial[j].resize(j+1U);
    for(unsigned k = 0; k <= j; ++k) monomial[j][k] = Binomial(j,k)*Power(middle,j-k)*Power(halfWidth,k);
  }
  std::vector<long double> term(1,1.0L);
  long double nextPower = ratio;
  long double magnitude = 0.0L;
  for(unsigned iteration = 0; iteration < 96; ++iteration) {
    for(unsigned j = 0; j <= degree; ++j) {
      long double contribution = 0.0L;
      for(unsigned k = 0; k < term.size(); ++k) {
        for(unsigned p = 0; p <= j; ++p) {
          if((k+p)%2U == 0) contribution += term[k]*monomial[j][p]*2.0L/(k+p+1U);
        }
      }
      contribution *= halfWidth/constant;
      values[j] += contribution;
      magnitude = std::max(magnitude,std::abs(values[j]));
    }
    long double remainder = 2.0L*halfWidth*nextPower/(constant*(1.0L-ratio));
    long double roundoff = 128.0L*std::numeric_limits<long double>::epsilon()*magnitude*(iteration+1U);
    if(remainder <= std::max(tolerance,roundoff)) {
      for(unsigned j = 0; j <= degree; ++j) errors[j] = remainder+roundoff;
      return;
    }
    term = PolynomialProduct(term,{0.0L,-linear,-quadratic});
    nextPower *= ratio;
  }
  throw std::runtime_error("A rational segment did not resolve its analytic remainder.");
}


}

namespace FELevelSetMoments {

inline long double LineSlicePolynomial(const std::vector<long double>& polynomial, long double x, long double offset = 0.0L) {
  long double high = 0.0L, low = 0.0L, derivative = 0.0L;
  for(unsigned k = static_cast<unsigned>(polynomial.size()); k > 0; --k) {
    derivative = std::fma(derivative,x,high);
    long double product = high*x;
    long double productError = std::fma(high,x,-product);
    long double sum = product+polynomial[k-1U];
    long double second = sum-product;
    long double sumError = (product-(sum-second))+(polynomial[k-1U]-second);
    low = std::fma(low,x,productError+sumError);
    high = sum;
  }
  return high+std::fma(derivative,offset,low);
}

inline long double LineSliceEventOffset(long double point, const std::vector<std::vector<long double>>& polynomials) {
  if(point == 0.0L || point == 1.0L) return 0.0L;
  long double offset = 0.0L;
  long double best = std::numeric_limits<long double>::infinity();
  for(unsigned p = 0; p < polynomials.size(); ++p) {
    if(polynomials[p].size() <= 1 || !LineSliceUtilities::NearZero(polynomials[p],point)) continue;
    std::vector<long double> derivative = LineSliceUtilities::PolynomialDerivative(polynomials[p]);
    long double slope = LineSliceUtilities::PolynomialValue(derivative,point);
    if(slope == 0.0L) continue;
    long double correction = -LineSlicePolynomial(polynomials[p],point)/slope;
    if(std::abs(correction) > 1024.0L*std::numeric_limits<long double>::epsilon()) continue;
    long double residual = std::abs(LineSlicePolynomial(polynomials[p],point,correction))/std::abs(slope);
    if(residual < best) { best = residual; offset = correction; }
  }
  return offset;
}

inline void LineSliceSegmentMoments(long double a, long double b, long double c, unsigned degree, long double lower, long double upper, bool computeDelta, bool fullBoundaryDelta, long double discriminant, long double jacobian, long double deltaJacobian, bool lowerZero, bool upperZero, std::vector<long double>& positive, std::vector<long double>& delta) {
  positive.assign(degree+1U,0.0L);
  delta.assign(degree+1U,0.0L);
  if(!(upper > lower)) return;
  long double roots[2] = {0.0L,0.0L};
  unsigned numberOfRoots = 0;
  long double derivative = 0.0L;
  if(a != 0.0L) {
    if(discriminant > 0.0L) {
      derivative = std::sqrt(discriminant);
      long double numerator = -0.5L*(b+std::copysign(derivative,b));
      roots[0] = numerator/a;
      roots[1] = numerator != 0.0L ? c/numerator : -b/a-roots[0];
      if(lowerZero || upperZero) {
        long double endpoint = lowerZero ? lower : upper;
        if(std::abs(roots[0]-endpoint) <= std::abs(roots[1]-endpoint)) roots[0] = endpoint;
        else roots[1] = endpoint;
        if(lowerZero && upperZero) { roots[0] = lower; roots[1] = upper; }
      }
      if(roots[1] < roots[0]) std::swap(roots[0],roots[1]);
      numberOfRoots = 2;
    }
    else if(discriminant == 0.0L) {
      roots[0] = -b/(2.0L*a);
      numberOfRoots = 1;
      if(computeDelta && roots[0] >= lower && roots[0] <= upper) throw std::domain_error("A repeated segment zero requires a different slice or endpoint treatment.");
    }
  }
  else if(b != 0.0L) {
    roots[0] = lowerZero ? lower : (upperZero ? upper : -c/b);
    numberOfRoots = 1;
    derivative = std::abs(b);
  }
  else if(c == 0.0L && computeDelta) throw std::domain_error("An entire zero segment requires a different slicing direction.");
  long double cuts[4] = {lower,upper,0.0L,0.0L};
  unsigned numberOfCuts = 2;
  for(unsigned r = 0; r < numberOfRoots; ++r) {
    if(roots[r] > lower && roots[r] < upper) cuts[numberOfCuts++] = roots[r];
  }
  for(unsigned k = 1; k < numberOfCuts; ++k) {
    long double value = cuts[k];
    unsigned position = k;
    while(position > 0 && cuts[position-1U] > value) {
      cuts[position] = cuts[position-1U];
      --position;
    }
    cuts[position] = value;
  }
  for(unsigned segment = 1; segment < numberOfCuts; ++segment) {
    long double left = cuts[segment-1U], right = cuts[segment];
    if(!(right > left)) continue;
    long double middle = 0.5L*(left+right);
    long double phase;
    if(a != 0.0L && numberOfRoots == 2) phase = ((middle > roots[0] && middle < roots[1]) == (a < 0.0L)) ? 1.0L : 0.0L;
    else phase = LineSliceUtilities::Heaviside(std::fma(a,middle*middle,std::fma(b,middle,c)));
    if(phase == 0.0L) continue;
    long double rightPower = 1.0L, geometric = 0.0L;
    for(unsigned j = 0; j <= degree; ++j) {
      geometric = rightPower+left*geometric;
      positive[j] += phase*jacobian*(right-left)*geometric/(j+1U);
      rightPower *= right;
    }
  }
  if(!computeDelta) return;
  for(unsigned r = 0; r < numberOfRoots; ++r) {
    long double root = roots[r];
    if(root < lower || root > upper) continue;
    if(!(derivative > 0.0L)) throw std::domain_error("The segment delta moment is undefined.");
    long double weight = std::isfinite(deltaJacobian) ? deltaJacobian : jacobian/derivative;
    if(!fullBoundaryDelta && (root == lower || root == upper)) weight *= 0.5L;
    long double power = 1.0L;
    for(unsigned j = 0; j <= degree; ++j) {
      delta[j] += weight*power;
      power *= root;
    }
  }
}

inline bool LineSliceZeroFiber(const std::vector<std::vector<long double>>& field) {
  std::vector<long double> candidates = {0.0L,1.0L};
  for(unsigned k = 0; k < 3; ++k) LineSliceUtilities::PolynomialRoots(field[k],0.0L,1.0L,candidates);
  for(unsigned r = 0; r < candidates.size(); ++r) {
    bool zero = true;
    for(unsigned k = 0; k < 3; ++k) zero = zero && LineSliceUtilities::NearZero(field[k],candidates[r]);
    if(zero) return true;
  }
  return false;
}

inline std::vector<long double> LineSliceFlatFibers(const std::vector<std::vector<long double>>& field) {
  long double norm = 0.0L;
  for(unsigned term = 0; term < 3; ++term) norm = std::max(norm,LineSliceUtilities::PolynomialNorm(field[term]));
  unsigned selected = 3;
  for(unsigned term = 0; term < 3; ++term) {
    if(LineSliceUtilities::PolynomialNorm(field[term]) != 0.0L && (selected == 3 || field[term].size() < field[selected].size())) selected = term;
  }
  std::vector<long double> candidates;
  std::vector<long double> roots;
  if(selected == 3) return roots;
  LineSliceUtilities::PolynomialRoots(field[selected],0.0L,1.0L,candidates);
  for(unsigned n = 0; n < candidates.size(); ++n) {
    bool common = true;
    for(unsigned term = 0; term < 3; ++term) common = common && std::abs(LineSliceUtilities::PolynomialValue(field[term],candidates[n])) <= 2048.0L*std::numeric_limits<long double>::epsilon()*std::max(norm,LineSliceUtilities::PolynomialBound(field[term],candidates[n]));
    if(common) roots.push_back(candidates[n]);
  }
  LineSliceUtilities::SortDistinct(roots);
  return roots;
}

inline void LineSliceSample(const std::vector<std::vector<long double>>& field, const std::vector<long double>& flat, const std::vector<long double>& reduced, unsigned leftMultiplicity, unsigned rightMultiplicity, long double left, long double right, long double leftOffset, long double rightOffset, long double parameter, bool linear, unsigned degree, unsigned shape, bool computeDelta, bool fullBoundaryDelta, bool lowerZero, bool upperZero, std::vector<long double>& values) {
  long double width = (right-left)+(rightOffset-leftOffset);
  long double sine = linear ? 0.0L : std::sin(parameter);
  long double cosine = linear ? 0.0L : std::cos(parameter);
  long double fromLeft = width*(linear ? parameter : sine*sine);
  long double fromRight = width*(linear ? 1.0L-parameter : cosine*cosine);
  long double x = fromLeft <= fromRight ? left+fromLeft : right-fromRight;
  long double first = fromLeft <= fromRight ? left : right;
  long double second = fromLeft <= fromRight ? fromLeft : -fromRight;
  long double rounded = x-first;
  long double coordinateError = (first-(x-rounded))+(second-rounded)+(fromLeft <= fromRight ? leftOffset : rightOffset);
  long double jacobian = linear ? width : 2.0L*width*sine*cosine;
  long double a = LineSliceUtilities::PolynomialValue(field[0],x);
  long double b = LineSliceUtilities::PolynomialValue(field[1],x);
  long double c = LineSliceUtilities::PolynomialValue(field[2],x);
  if(std::abs(a)+std::abs(b)+std::abs(c) < 1.0e-5L) {
    a = LineSlicePolynomial(field[0],x,coordinateError);
    b = LineSlicePolynomial(field[1],x,coordinateError);
    c = LineSlicePolynomial(field[2],x,coordinateError);
  }
  long double factor = 1.0L;
  for(unsigned n = 0; n < flat.size(); ++n) {
    long double distance = x-flat[n]+coordinateError;
    if(flat[n] == left) distance = fromLeft+leftOffset;
    if(flat[n] == right) distance = -fromRight+rightOffset;
    factor *= distance;
  }
  if(factor == 0.0L) throw std::runtime_error("A factored line sample reached its zero fiber.");
  if(factor < 0.0L) {
    a = -a;
    b = -b;
    c = -c;
  }
  long double discriminant = std::fma(-4.0L*a,c,b*b);
  long double deltaJacobian = std::numeric_limits<long double>::quiet_NaN();
  if(a != 0.0L && !reduced.empty()) {
    long double reducedValue = LineSliceUtilities::PolynomialValue(reduced,x);
    if(std::abs(reducedValue) < 1.0e-5L) reducedValue = LineSlicePolynomial(reduced,x,coordinateError);
    discriminant = reducedValue*LineSliceUtilities::Power(fromLeft,leftMultiplicity)*LineSliceUtilities::Power(fromRight,rightMultiplicity);
    if(computeDelta && reducedValue > 0.0L && !linear && leftMultiplicity <= 1 && rightMultiplicity <= 1) {
      deltaJacobian = 2.0L*width/std::sqrt(reducedValue*LineSliceUtilities::Power(width,leftMultiplicity+rightMultiplicity));
      if(leftMultiplicity == 0) deltaJacobian *= sine;
      if(rightMultiplicity == 0) deltaJacobian *= cosine;
    }
  }
  long double lower = shape == 3 ? x : 0.0L;
  long double upper = shape == 2 ? x : 1.0L;
  std::vector<long double> h,d;
  LineSliceSegmentMoments(a,b,c,degree,lower,upper,computeDelta,fullBoundaryDelta,discriminant,jacobian,deltaJacobian,lowerZero,upperZero,h,d);
  if(computeDelta) {
    for(unsigned j = 0; j <= degree; ++j) d[j] /= std::abs(factor);
  }
  unsigned count = (degree+1U)*(degree+2U)/2U;
  values.assign(2U*count,0.0L);
  long double power = 1.0L;
  for(unsigned i = 0; i <= degree; ++i) {
    for(unsigned j = 0; j <= degree-i; ++j) {
      unsigned index = LineSliceUtilities::MomentIndex(degree,i,j);
      values[index] = power*h[j];
      values[count+index] = power*d[j];
    }
    power *= x;
  }
}

inline void LineSliceAdaptiveInterval(const std::vector<std::vector<long double>>& field, const std::vector<long double>& flat, const std::vector<long double>& reduced, unsigned leftMultiplicity, unsigned rightMultiplicity, long double left, long double right, long double leftOffset, long double rightOffset, long double parameterLeft, long double parameterRight, bool linear, unsigned degree, unsigned shape, bool computeDelta, bool fullBoundaryDelta, bool lowerZero, bool upperZero, long double absoluteH, long double absoluteDelta, long double relativeTolerance, unsigned depth, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, std::vector<long double>& total, std::vector<long double>& totalError) {
  const long double nodes[8] = {0.991455371120812639206854697526329L,0.949107912342758524526189684047851L,0.864864423359769072789712788640926L,0.741531185599394439863864773280788L,0.586087235467691130294144838258730L,0.405845151377397166906606412076961L,0.207784955007898467600689403773245L,0.0L};
  const long double kronrod[8] = {0.022935322010529224963732008058970L,0.063092092629978553290700663189204L,0.104790010322250183839876322541518L,0.140653259715525918745189590510238L,0.169004726639267902826583426598550L,0.190350578064785409913256402421014L,0.204432940075298892414161999234649L,0.209482141084727828012999174891714L};
  const long double gauss[4] = {0.129484966168869693270611432679082L,0.279705391489276667901467771423780L,0.381830050505118944950369775488975L,0.417959183673469387755102040816327L};
  long double middle = 0.5L*(parameterLeft+parameterRight);
  long double halfWidth = 0.5L*(parameterRight-parameterLeft);
  unsigned count = (degree+1U)*(degree+2U)/2U;
  std::vector<long double> high(2U*count,0.0L), low(2U*count,0.0L), absolute(2U*count,0.0L);
  std::vector<std::vector<long double>> samples(15);
  long double sampleWeights[15] = {};
  unsigned numberOfSamples = 0;
  for(unsigned node = 0; node < 8; ++node) {
    for(unsigned side = 0; side < (node == 7 ? 1U : 2U); ++side) {
      if(evaluations >= maximumEvaluations) throw std::runtime_error("The line-slice evaluation budget was exhausted.");
      ++evaluations;
      long double point = middle+(side == 0 ? -1.0L : 1.0L)*halfWidth*nodes[node];
      LineSliceSample(field,flat,reduced,leftMultiplicity,rightMultiplicity,left,right,leftOffset,rightOffset,point,linear,degree,shape,computeDelta,fullBoundaryDelta,lowerZero,upperZero,samples[numberOfSamples]);
      const std::vector<long double>& value = samples[numberOfSamples];
      sampleWeights[numberOfSamples++] = kronrod[node];
      for(unsigned k = 0; k < value.size(); ++k) {
        if(!std::isfinite(value[k])) throw std::runtime_error("A line-slice sample exceeds floating point range.");
        high[k] += kronrod[node]*value[k];
        absolute[k] += kronrod[node]*std::abs(value[k]);
        if((node&1U) != 0) low[k] += gauss[node/2U]*value[k];
      }
    }
  }
  long double fraction = std::max(1.0e-2L,(right-left)*(parameterRight-parameterLeft)/(linear ? 1.0L : 0.5L*std::acos(-1.0L)));
  std::vector<long double> error(high.size());
  bool converged = true;
  for(unsigned k = 0; k < high.size(); ++k) {
    long double variation = 0.0L;
    for(unsigned n = 0; n < numberOfSamples; ++n) variation += sampleWeights[n]*std::abs(samples[n][k]-0.5L*high[k]);
    long double estimate = std::abs(high[k]-low[k]);
    if(variation > 0.0L && estimate > 0.0L) estimate = variation*std::min(1.0L,std::pow(200.0L*estimate/variation,1.5L));
    error[k] = halfWidth*std::max(estimate,128.0L*std::numeric_limits<long double>::epsilon()*absolute[k]);
    high[k] *= halfWidth;
    long double tolerance = (k < count ? absoluteH : absoluteDelta)*fraction+relativeTolerance*std::abs(high[k]);
    converged = converged && error[k] <= tolerance;
  }
  if(converged) {
    for(unsigned k = 0; k < high.size(); ++k) { total[k] += high[k]; totalError[k] += error[k]; }
    ++accepted;
    return;
  }
  if(depth >= maximumDepth || middle == parameterLeft || middle == parameterRight) throw std::runtime_error("Line-slice quadrature did not reach the requested tolerance.");
  LineSliceAdaptiveInterval(field,flat,reduced,leftMultiplicity,rightMultiplicity,left,right,leftOffset,rightOffset,parameterLeft,middle,linear,degree,shape,computeDelta,fullBoundaryDelta,lowerZero,upperZero,absoluteH,absoluteDelta,relativeTolerance,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError);
  LineSliceAdaptiveInterval(field,flat,reduced,leftMultiplicity,rightMultiplicity,left,right,leftOffset,rightOffset,middle,parameterRight,linear,degree,shape,computeDelta,fullBoundaryDelta,lowerZero,upperZero,absoluteH,absoluteDelta,relativeTolerance,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError);
}

inline void LineSliceDirection(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool triangle, unsigned order, bool computeDelta, bool fullBoundaryDelta, long double absoluteH, long double absoluteDelta, long double relativeTolerance, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned guaranteedMonotoneOrder = 0, const std::vector<unsigned>* zeroBoundary = nullptr) {
  std::vector<std::vector<long double>> field(3,std::vector<long double>(3,0.0L));
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) field[2U-j][i] = order == 1 ? coefficient[i][j] : coefficient[j][i];
  }
  for(unsigned k = 0; k < 3; ++k) LineSliceUtilities::Trim(field[k]);
  int phase = 0;
  bool linear = order == guaranteedMonotoneOrder || LineSliceUtilities::CanonicalGraphOrder(coefficient,order,phase) != 0;
  std::vector<long double> flat = LineSliceFlatFibers(field);
  if(computeDelta && zeroBoundary != nullptr) {
    bool entire = triangle ? (*zeroBoundary)[order == 1 ? 1U : 0U] != 0 : ((*zeroBoundary)[order == 1 ? 0U : 2U] != 0 || (*zeroBoundary)[order == 1 ? 1U : 3U] != 0);
    if(entire) {
      if(triangle) flat.push_back(order == 1 ? 1.0L : 0.0L);
      else {
        if((*zeroBoundary)[order == 1 ? 0U : 2U] != 0) flat.push_back(0.0L);
        if((*zeroBoundary)[order == 1 ? 1U : 3U] != 0) flat.push_back(1.0L);
      }
    }
  }
  LineSliceUtilities::SortDistinct(flat);
  std::vector<std::vector<long double>> originalField = field;
  for(unsigned n = 0; n < flat.size(); ++n) {
    for(unsigned term = 0; term < 3; ++term) {
      if(field[term].size() > 1) field[term] = LineSliceUtilities::DivideLinear(field[term],flat[n]);
      else field[term].assign(1,0.0L);
    }
  }
  unsigned shape = triangle ? (order == 1 ? 2U : 3U) : 0U;
  const std::vector<long double> domain = {0.0L,1.0L,0.0L,1.0L};
  std::vector<long double> discriminant = LineSliceUtilities::DiscriminantPolynomial(field);
  std::vector<long double> bottom = LineSliceUtilities::BoundaryPolynomial(field,false,shape,domain);
  std::vector<long double> top = LineSliceUtilities::BoundaryPolynomial(field,true,shape,domain);
  if(zeroBoundary != nullptr) {
    unsigned lowerIndex = triangle ? (order == 1 ? 0U : 2U) : (order == 1 ? 2U : 0U);
    unsigned upperIndex = triangle ? (order == 1 ? 2U : 1U) : (order == 1 ? 3U : 1U);
    if((*zeroBoundary)[lowerIndex] != 0) bottom.assign(1,0.0L);
    if((*zeroBoundary)[upperIndex] != 0) top.assign(1,0.0L);
  }
  std::vector<long double> events = {0.0L,1.0L};
  if(!linear) {
    LineSliceUtilities::PolynomialRoots(discriminant,0.0L,1.0L,events);
    for(unsigned k = 0; k < 3; ++k) LineSliceUtilities::PolynomialRoots(field[k],0.0L,1.0L,events);
  }
  LineSliceUtilities::PolynomialRoots(bottom,0.0L,1.0L,events);
  LineSliceUtilities::PolynomialRoots(top,0.0L,1.0L,events);
  for(unsigned n = 0; n < flat.size(); ++n) events.push_back(flat[n]);
  LineSliceUtilities::SortDistinct(events);
  const std::vector<std::vector<long double>> eventPolynomials = {discriminant,field[0],field[1],field[2],bottom,top};
  std::vector<long double> eventOffsets(events.size(),0.0L);
  for(unsigned k = 0; k < events.size(); ++k) eventOffsets[k] = LineSliceEventOffset(events[k],eventPolynomials);
  for(unsigned n = 0; n < flat.size(); ++n) {
    for(unsigned k = 0; k < events.size(); ++k) {
      if(events[k] == flat[n]) eventOffsets[k] = 0.0L;
    }
  }
  if(computeDelta && !linear) LineSliceUtilities::CheckSingularZeros(originalField,events,shape,domain);
  bool lowerZero = LineSliceUtilities::PolynomialNorm(bottom) == 0.0L;
  bool upperZero = LineSliceUtilities::PolynomialNorm(top) == 0.0L;
  unsigned count = (degree+1U)*(degree+2U)/2U;
  values.assign(2U*count,0.0L);
  errors.assign(2U*count,0.0L);
  for(unsigned n = 0; computeDelta && n < flat.size(); ++n) {
    long double root = flat[n];
    long double lower = LineSliceUtilities::LowerBoundary(root,shape,domain);
    long double upper = LineSliceUtilities::UpperBoundary(root,shape,domain);
    if(!(upper > lower)) continue;
    long double a = LineSliceUtilities::PolynomialValue(LineSliceUtilities::PolynomialDerivative(originalField[0]),root);
    long double b = LineSliceUtilities::PolynomialValue(LineSliceUtilities::PolynomialDerivative(originalField[1]),root);
    long double c = LineSliceUtilities::PolynomialValue(LineSliceUtilities::PolynomialDerivative(originalField[2]),root);
    std::vector<long double> singular = LineSliceUtilities::QuadraticRoots(a,b,c);
    for(unsigned r = 0; r < singular.size(); ++r) {
      if(singular[r] >= lower && singular[r] <= upper) throw std::domain_error("An entire zero segment has a critical transverse derivative.");
    }
    long double sign = std::fma(a,0.25L*(lower+upper)*(lower+upper),std::fma(b,0.5L*(lower+upper),c)) > 0.0L ? 1.0L : -1.0L;
    std::vector<long double> moments;
    std::vector<long double> momentError;
    LineSliceUtilities::RationalSegmentMoments(sign*a,sign*b,sign*c,lower,upper,degree,0.002L*absoluteDelta,0,moments,momentError);
    long double boundary = !fullBoundaryDelta && (root == 0.0L || root == 1.0L) ? 0.5L : 1.0L;
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        unsigned index = LineSliceUtilities::MomentIndex(degree,i,j);
        values[count+index] += boundary*LineSliceUtilities::Power(root,i)*moments[j];
        errors[count+index] += boundary*LineSliceUtilities::Power(root,i)*momentError[j];
      }
    }
  }
  for(unsigned interval = 1; interval < events.size(); ++interval) {
    long double left = events[interval-1U], right = events[interval];
    if(!(right > left)) continue;
    long double middle = 0.5L*(left+right);
    std::vector<int> signs,boundary;
    std::vector<long double> segments;
    long double uniform = 0.0L;
    long double factor = 1.0L;
    for(unsigned n = 0; n < flat.size(); ++n) factor *= middle-flat[n];
    std::vector<std::vector<long double>> signedField = field;
    if(factor < 0.0L) {
      for(unsigned term = 0; term < 3; ++term) {
        for(unsigned k = 0; k < signedField[term].size(); ++k) signedField[term][k] = -signedField[term][k];
      }
    }
    LineSliceUtilities::IntervalBranches(signedField,middle,shape,domain,lowerZero,upperZero,computeDelta,signs,boundary,segments,uniform);
    if(signs.empty() && uniform != 0.5L && !(computeDelta && (lowerZero || upperZero))) {
      if(uniform > 0.0L) {
        for(unsigned i = 0; i <= degree; ++i) {
          for(unsigned j = 0; j <= degree-i; ++j) values[LineSliceUtilities::MomentIndex(degree,i,j)] += LineSliceUtilities::FullStripMoment(left,right,i,j,shape,domain);
        }
      }
      ++accepted;
      continue;
    }
    std::vector<long double> reduced;
    unsigned leftMultiplicity = 0, rightMultiplicity = 0;
    if(!linear && LineSliceUtilities::PolynomialNorm(field[0]) != 0.0L) {
      reduced = discriminant;
      if(!linear) {
        LineSliceUtilities::FactorEndpoint(reduced,left,leftMultiplicity);
        LineSliceUtilities::FactorEndpoint(reduced,right,rightMultiplicity);
      }
      if((rightMultiplicity&1U) != 0) {
        for(unsigned k = 0; k < reduced.size(); ++k) reduced[k] = -reduced[k];
      }
      if(!(LineSliceUtilities::PolynomialValue(reduced,middle) > 0.0L)) {
        if(LineSliceUtilities::PolynomialNorm(discriminant) == 0.0L && !computeDelta) reduced.clear();
        else throw std::runtime_error("A line-slice discriminant chart could not be resolved.");
      }
    }
    bool intervalLinear = linear && leftMultiplicity == 0 && rightMultiplicity == 0;
    LineSliceAdaptiveInterval(field,flat,reduced,leftMultiplicity,rightMultiplicity,left,right,eventOffsets[interval-1U],eventOffsets[interval],0.0L,intervalLinear ? 1.0L : 0.5L*std::acos(-1.0L),intervalLinear,degree,shape,computeDelta,fullBoundaryDelta,lowerZero,upperZero,absoluteH,absoluteDelta,relativeTolerance,0,maximumDepth,maximumEvaluations,evaluations,accepted,values,errors);
  }
}


inline void LineSliceTriangleChart(const std::vector<std::vector<long double>>& input, unsigned chart, std::vector<std::vector<long double>>& output) {
  output = input;
  if(chart == 0) return;
  const long double c0 = input[0][0];
  const long double c1 = input[1][0];
  const long double c2 = input[0][1];
  const long double c3 = input[2][0];
  const long double c4 = input[1][1];
  const long double c5 = input[0][2];
  const long double c6 = input[2][1];
  if(chart == 1) {
    output[0][0] = c0+c1+c2+c3+c4+c5;
    output[1][0] = -c1-c2-2.0L*(c3+c4+c5);
    output[0][1] = c1+2.0L*c3+c4+c6;
    output[2][0] = c3+c4+c5;
    output[1][1] = -2.0L*c3-c4-2.0L*c6;
    output[0][2] = c3+c6;
  }
  else {
    output[0][0] = c0+c1+c3;
    output[1][0] = c2+c4+c6;
    output[0][1] = -c1-c2-2.0L*c3-c4-c6;
    output[2][0] = c5-c6;
    output[1][1] = -c4-2.0L*c5;
    output[0][2] = c3+c4+c5+c6;
  }
}

inline void LineSliceTriangleDerivativeBounds(const std::vector<std::vector<long double>>& coefficient, long double& lower, long double& upper) {
  const long double c2 = coefficient[0][1];
  const long double c4 = coefficient[1][1];
  const long double c5 = coefficient[0][2];
  const long double c6 = coefficient[2][1];
  const long double bernstein[6] = {c2,c2+0.5L*c4,c2+0.5L*c4+c5,c2+c4+c6,c2+c4+c5,c2+c4+2.0L*c5-c6};
  long double minimum = bernstein[0];
  long double maximum = bernstein[0];
  long double norm = 0.0L;
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) norm += std::abs(coefficient[i][j]);
  }
  for(unsigned k = 1; k < 6; ++k) {
    minimum = std::min(minimum,bernstein[k]);
    maximum = std::max(maximum,bernstein[k]);
  }
  const long double padding = 2048.0L*std::numeric_limits<long double>::epsilon()*norm;
  lower = minimum > padding ? minimum-padding : (maximum < -padding ? -maximum-padding : 0.0L);
  upper = std::max(std::abs(minimum),std::abs(maximum))+padding;
}

inline void LineSliceTriangleRestoreMoments(unsigned degree, unsigned chart, const std::vector<long double>& source, const std::vector<long double>& sourceError, std::vector<long double>& values, std::vector<long double>& errors) {
  if(chart == 0) {
    values = source;
    errors = sourceError;
    return;
  }
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  values.assign(2U*count,0.0L);
  errors.assign(2U*count,0.0L);
  for(unsigned i = 0; i <= degree; ++i) {
    for(unsigned j = 0; j <= degree-i; ++j) {
      const unsigned target = LineSliceUtilities::MomentIndex(degree,i,j);
      long double absoluteSum[2] = {0.0L,0.0L};
      long double correction[2] = {0.0L,0.0L};
      unsigned terms = 0;
      for(unsigned a = 0; a <= (chart == 1 ? i : j); ++a) {
        const unsigned bound = chart == 1 ? i+j-a : i;
        for(unsigned b = 0; b <= bound; ++b) {
          const unsigned first = chart == 1 ? b : a;
          const unsigned second = chart == 1 ? a : j-a+b;
          const unsigned index = LineSliceUtilities::MomentIndex(degree,first,second);
          long double multiplier = chart == 1 ? LineSliceUtilities::Binomial(i,a)*LineSliceUtilities::Binomial(i+j-a,b) : LineSliceUtilities::Binomial(j,a)*LineSliceUtilities::Binomial(i,b);
          if((chart == 1 ? b : j-a+b)%2U != 0) multiplier = -multiplier;
          for(unsigned pair = 0; pair < 2; ++pair) {
            const unsigned destination = pair*count+target;
            const long double term = multiplier*source[pair*count+index];
            const long double adjusted = term-correction[pair];
            const long double next = values[destination]+adjusted;
            correction[pair] = (next-values[destination])-adjusted;
            values[destination] = next;
            errors[destination] += std::abs(multiplier)*sourceError[pair*count+index];
            absoluteSum[pair] += std::abs(term);
          }
          ++terms;
        }
      }
      for(unsigned pair = 0; pair < 2; ++pair) errors[pair*count+target] += (8.0L+2.0L*terms)*std::numeric_limits<long double>::epsilon()*absoluteSum[pair];
    }
  }
}

inline void LineSliceTriangleDirections(const std::vector<std::vector<long double>>& normalized, long double scale, unsigned degree, bool computeDelta, bool fullBoundaryDelta, long double absoluteTolerance, long double relativeTolerance, unsigned maximumDepth, unsigned integrationOrder, unsigned maximumEvaluations, const std::vector<unsigned>* zeroBoundary, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned& selected) {
  std::vector<std::vector<std::vector<long double>>> charts(3);
  long double lower[3];
  long double upper[3];
  unsigned order[3] = {0,1,2};
  for(unsigned chart = 0; chart < 3; ++chart) {
    LineSliceTriangleChart(normalized,chart,charts[chart]);
    LineSliceTriangleDerivativeBounds(charts[chart],lower[chart],upper[chart]);
  }
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = i+1U; j < 3; ++j) {
      const bool better = lower[order[j]] > lower[order[i]] || (lower[order[j]] == lower[order[i]] && upper[order[j]] > upper[order[i]]);
      if(better) std::swap(order[i],order[j]);
    }
  }
  if(integrationOrder != 0) {
    for(unsigned i = 0; i < 3; ++i) {
      if(order[i] == integrationOrder-1U) std::swap(order[0],order[i]);
    }
  }
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  std::string failure = "No triangle direction converged.";
  for(unsigned attempt = 0; attempt < 3; ++attempt) {
    const unsigned chart = order[attempt];
    std::vector<unsigned> boundary;
    if(zeroBoundary != nullptr) {
      boundary.resize(3);
      for(unsigned side = 0; side < 3; ++side) boundary[side] = (*zeroBoundary)[(side+3U-chart)%3U];
    }
    const long double amplification = chart == 0 ? 1.0L : LineSliceUtilities::Power(chart == 1 ? 3.0L : 2.0L,degree);
    std::vector<long double> trial;
    std::vector<long double> trialError;
    unsigned trialAccepted = 0;
    try {
      LineSliceDirection(charts[chart],degree,true,1,computeDelta,fullBoundaryDelta,absoluteTolerance/amplification,absoluteTolerance*scale/amplification,relativeTolerance/amplification,maximumDepth,maximumEvaluations,trial,trialError,evaluations,trialAccepted,lower[chart] > 0.0L ? 1U : 0U,zeroBoundary == nullptr ? nullptr : &boundary);
      LineSliceTriangleRestoreMoments(degree,chart,trial,trialError,values,errors);
      for(unsigned k = 0; k < values.size(); ++k) {
        const long double tolerance = (k < count ? absoluteTolerance : absoluteTolerance*scale)+relativeTolerance*std::abs(values[k]);
        if(!std::isfinite(values[k]) || !std::isfinite(errors[k]) || errors[k] > tolerance) throw std::runtime_error("Triangle chart conversion exceeds the requested error budget.");
      }
      accepted = trialAccepted;
      selected = chart+1U;
      return;
    }
    catch(const std::runtime_error& exception) {
      failure = exception.what();
    }
  }
  throw std::runtime_error(failure);
}

inline void IntegrateNormalizedByLineSlices(const std::vector<std::vector<long double>>& normalized, long double scale, unsigned degree, bool triangle, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned maximumDepth = 28, unsigned integrationOrder = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr, const std::vector<unsigned>* zeroBoundary = nullptr) {
  if(degree > 10000U || !(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance) || maximumDepth > 60 || integrationOrder > (triangle ? 3U : 2U) || maximumEvaluations < 15) throw std::invalid_argument("Invalid adaptive line-slice input.");
  unsigned long long size = (static_cast<unsigned long long>(degree)+1ULL)*(degree+2ULL)/2ULL;
  if(size > std::numeric_limits<unsigned>::max()/2U) throw std::length_error("The moment table is too large.");
  if(zeroBoundary != nullptr && zeroBoundary->size() != (triangle ? 3U : 4U)) throw std::invalid_argument("Invalid nodal boundary flags.");
  unsigned count = static_cast<unsigned>(size);
  std::vector<long double> values(2U*count,0.0L),errors(2U*count,0.0L);
  int phase = 0;
  unsigned monotoneOrder = LineSliceUtilities::CanonicalGraphOrder(normalized,0,phase);
  if(computeDelta && zeroBoundary != nullptr) {
    for(unsigned side = 0; side < zeroBoundary->size(); ++side) {
      if((*zeroBoundary)[side] != 0) phase = 0;
    }
  }
  unsigned evaluations = 0, accepted = 0, selected = 0;
  if(phase != 0) {
    if(phase > 0) {
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) values[LineSliceUtilities::MomentIndex(degree,i,j)] = LineSliceUtilities::CanonicalFullMoment(i,j,triangle);
      }
    }
  }
  else if(triangle && normalized[2][2] == 0.0L && normalized[2][1] == -normalized[1][2]) {
    LineSliceTriangleDirections(normalized,scale,degree,computeDelta,fullBoundaryDelta,absoluteTolerance,relativeTolerance,maximumDepth,integrationOrder,maximumEvaluations,zeroBoundary,values,errors,evaluations,accepted,selected);
  }
  else {
    if(integrationOrder > 2) throw std::invalid_argument("The third triangle direction requires the enriched TRI7 polynomial space.");
    bool completed = false;
    std::string failure;
    const std::vector<long double> domain = {0.0L,1.0L,0.0L,1.0L};
    unsigned preferred = integrationOrder == 0 ? (monotoneOrder != 0 ? monotoneOrder : (LineSliceUtilities::PreferReverseOrder(normalized,domain) ? 2U : 1U)) : integrationOrder;
    for(unsigned attempt = 0; attempt < 2 && !completed; ++attempt) {
      unsigned order = attempt == 0 ? preferred : 3U-preferred;
      accepted = 0;
      try {
        LineSliceDirection(normalized,degree,triangle,order,computeDelta,fullBoundaryDelta,absoluteTolerance,static_cast<long double>(absoluteTolerance)*scale,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,monotoneOrder,zeroBoundary);
        if(order == 2) {
          std::vector<long double> copy(values),copyError(errors);
          for(unsigned i = 0; i <= degree; ++i) {
            for(unsigned j = 0; j <= degree-i; ++j) {
              unsigned target = LineSliceUtilities::MomentIndex(degree,i,j), source = LineSliceUtilities::MomentIndex(degree,j,i);
              values[target] = copy[source]; values[count+target] = copy[count+source];
              errors[target] = copyError[source]; errors[count+target] = copyError[count+source];
            }
          }
        }
        for(unsigned k = 0; k < values.size(); ++k) {
          long double tolerance = (k < count ? absoluteTolerance : static_cast<long double>(absoluteTolerance)*scale)+relativeTolerance*std::abs(values[k]);
          if(errors[k] > tolerance) throw std::runtime_error("The accumulated line-slice error estimate exceeds tolerance.");
        }
        completed = true;
        selected = order;
      }
      catch(const std::runtime_error& exception) { failure = exception.what(); }
    }
    if(!completed) throw std::runtime_error(failure);
  }
  positive.resize(degree+1U);
  delta.resize(degree+1U);
  long double maxH = 0.0L, maxD = 0.0L;
  for(unsigned i = 0; i <= degree; ++i) {
    positive[i].resize(degree-i+1U);
    delta[i].resize(degree-i+1U);
    for(unsigned j = 0; j <= degree-i; ++j) {
      unsigned index = LineSliceUtilities::MomentIndex(degree,i,j);
      long double h = values[index], d = values[count+index]/scale;
      if(!std::isfinite(h) || !std::isfinite(d) || std::abs(d) > std::numeric_limits<double>::max()) throw std::overflow_error("A returned line-slice moment exceeds double range.");
      positive[i][j] = static_cast<double>(h);
      delta[i][j] = static_cast<double>(d);
      maxH = std::max(maxH,errors[index]);
      maxD = std::max(maxD,errors[count+index]/scale);
    }
  }
  if(estimatedError != nullptr) *estimatedError = {static_cast<double>(maxH),static_cast<double>(maxD),static_cast<double>(accepted),static_cast<double>(evaluations),static_cast<double>(selected)};
}

inline void IntegrateCanonicalByLineSlices(const std::vector<std::vector<double>>& coefficient, unsigned degree, bool triangle, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned maximumDepth = 28, unsigned integrationOrder = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr, const std::vector<unsigned>* zeroBoundary = nullptr) {
  std::vector<std::vector<long double>> normalized;
  long double scale;
  LineSliceUtilities::NormalizeCanonicalPolynomial(coefficient,normalized,scale);
  IntegrateNormalizedByLineSlices(normalized,scale,degree,triangle,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta,maximumDepth,integrationOrder,maximumEvaluations,estimatedError,zeroBoundary);
}

inline void IntegrateQuad9UnitByLineSlices(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned maximumDepth = 28, unsigned integrationOrder = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr) {
  IntegrateCanonicalByLineSlices(coefficient,degree,false,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta,maximumDepth,integrationOrder,maximumEvaluations,estimatedError);
}

inline void IntegrateQuad9NodalByLineSlices(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false) {
  std::vector<std::vector<double>> coefficient;
  LineSliceUtilities::Quad9UnitCoefficientsFromNodalValues(nodalValues,coefficient);
  std::vector<unsigned> zeroBoundary = {static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[3] == 0.0 && nodalValues[7] == 0.0),static_cast<unsigned>(nodalValues[1] == 0.0 && nodalValues[2] == 0.0 && nodalValues[5] == 0.0),static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[1] == 0.0 && nodalValues[4] == 0.0),static_cast<unsigned>(nodalValues[3] == 0.0 && nodalValues[2] == 0.0 && nodalValues[6] == 0.0)};
  IntegrateCanonicalByLineSlices(coefficient,degree,false,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta,28,0,200000,nullptr,&zeroBoundary);
}

inline void IntegrateTri7ReflectedByLineSlices(const std::vector<double>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned maximumDepth = 28, unsigned integrationOrder = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr, const std::vector<unsigned>* zeroBoundary = nullptr) {
  if(coefficient.size() != 7) throw std::invalid_argument("TRI7 requires seven coefficients.");
  std::vector<std::vector<double>> polynomial(3,std::vector<double>(3,0.0));
  polynomial[0][0] = coefficient[0]; polynomial[1][0] = coefficient[1]; polynomial[0][1] = coefficient[2];
  polynomial[2][0] = coefficient[3]; polynomial[1][1] = coefficient[4]; polynomial[0][2] = coefficient[5];
  polynomial[2][1] = coefficient[6]; polynomial[1][2] = -coefficient[6];
  IntegrateCanonicalByLineSlices(polynomial,degree,true,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta,maximumDepth,integrationOrder,maximumEvaluations,estimatedError,zeroBoundary);
}

inline void IntegrateTri7NodalByLineSlices(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false) {
  std::vector<double> coefficient;
  LineSliceUtilities::Tri7ReflectedCoefficientsFromNodalValues(nodalValues,coefficient);
  std::vector<unsigned> zeroBoundary = {static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[1] == 0.0 && nodalValues[3] == 0.0),static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[2] == 0.0 && nodalValues[5] == 0.0),static_cast<unsigned>(nodalValues[1] == 0.0 && nodalValues[2] == 0.0 && nodalValues[4] == 0.0)};
  IntegrateTri7ReflectedByLineSlices(coefficient,degree,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta,28,0,200000,nullptr,&zeroBoundary);
}

inline void IntegrateTri7ReflectedByLineSlicesLongDouble(const std::vector<long double>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned maximumDepth = 28, unsigned integrationOrder = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr, const std::vector<unsigned>* zeroBoundary = nullptr) {
  if(coefficient.size() != 7) throw std::invalid_argument("TRI7 requires seven coefficients.");
  long double scale = 0.0L;
  for(unsigned term = 0; term < 7; ++term) {
    if(!std::isfinite(coefficient[term])) throw std::invalid_argument("A TRI7 coefficient is not finite.");
    scale = std::max(scale,std::abs(coefficient[term]));
  }
  if(scale == 0.0L) scale = 1.0L;
  int exponent = 0;
  std::frexp(scale,&exponent);
  scale = std::ldexp(1.0L,exponent-1);
  std::vector<std::vector<long double>> normalized(3,std::vector<long double>(3,0.0L));
  normalized[0][0] = coefficient[0]/scale;
  normalized[1][0] = coefficient[1]/scale;
  normalized[0][1] = coefficient[2]/scale;
  normalized[2][0] = coefficient[3]/scale;
  normalized[1][1] = coefficient[4]/scale;
  normalized[0][2] = coefficient[5]/scale;
  normalized[2][1] = coefficient[6]/scale;
  normalized[1][2] = -coefficient[6]/scale;
  IntegrateNormalizedByLineSlices(normalized,scale,degree,true,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta,maximumDepth,integrationOrder,maximumEvaluations,estimatedError,zeroBoundary);
}

inline void IntegrateQuad9CanonicalByLineSlices(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, double xLeft = -1.0, double xRight = 1.0, double yBottom = -1.0, double yTop = 1.0, bool computeDelta = true, bool fullBoundaryDelta = false) {
  std::vector<std::vector<double>> transformed;
  LineSliceUtilities::Quad9CoefficientsToUnitSquare(coefficient,transformed,xLeft,xRight,yBottom,yTop);
  IntegrateQuad9UnitByLineSlices(transformed,degree,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta);
}

inline void IntegrateTri7CanonicalByLineSlices(const std::vector<double>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false) {
  std::vector<double> transformed;
  LineSliceUtilities::Tri7CoefficientsToReflectedTriangle(coefficient,transformed);
  IntegrateTri7ReflectedByLineSlices(transformed,degree,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta);
}

inline void IntegrateSegmentByLineSlices(double a, double b, double c, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, double left = 0.0, double right = 1.0, bool computeDelta = true, bool fullBoundaryDelta = false, bool forceLeftZero = false, bool forceRightZero = false) {
  if(!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(left) || !std::isfinite(right) || !(right > left) || degree > 10000U) throw std::invalid_argument("Invalid analytic segment input.");
  long double scale = std::max(std::abs(static_cast<long double>(a)),std::max(std::abs(static_cast<long double>(b)),std::abs(static_cast<long double>(c))));
  if(scale == 0.0L) scale = 1.0L;
  int exponent = 0;
  std::frexp(scale,&exponent);
  scale = std::ldexp(1.0L,exponent-1);
  long double qa = a/scale, qb = b/scale, qc = c/scale;
  std::vector<long double> h,d;
  bool lowerZero = forceLeftZero || std::fma(static_cast<long double>(a),static_cast<long double>(left)*left,std::fma(static_cast<long double>(b),static_cast<long double>(left),static_cast<long double>(c))) == 0.0L;
  bool upperZero = forceRightZero || std::fma(static_cast<long double>(a),static_cast<long double>(right)*right,std::fma(static_cast<long double>(b),static_cast<long double>(right),static_cast<long double>(c))) == 0.0L;
  LineSliceSegmentMoments(qa,qb,qc,degree,left,right,computeDelta,fullBoundaryDelta,std::fma(-4.0L*qa,qc,qb*qb),1.0L,std::numeric_limits<long double>::quiet_NaN(),lowerZero,upperZero,h,d);
  positive.resize(degree+1U);
  delta.resize(degree+1U);
  for(unsigned j = 0; j <= degree; ++j) {
    positive[j] = static_cast<double>(h[j]);
    delta[j] = static_cast<double>(d[j]/scale);
    if(!std::isfinite(positive[j]) || !std::isfinite(delta[j])) throw std::overflow_error("An analytic segment moment exceeds double range.");
  }
}

inline void IntegrateLine3NodalByLineSlices(const std::vector<double>& nodalValues, unsigned degree, std::vector<double>& positive, std::vector<double>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false) {
  if(!(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance)) throw std::invalid_argument("Invalid line tolerance.");
  double a,b,c;
  LineSliceUtilities::Line3CoefficientsFromNodalValues(nodalValues,a,b,c);
  IntegrateSegmentByLineSlices(a,b,c,degree,positive,delta,0.0,1.0,computeDelta,fullBoundaryDelta,nodalValues[0] == 0.0,nodalValues[1] == 0.0);
}

}

#endif
