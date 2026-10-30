#ifndef FEMUS_SQUARE_PARABOLA_MOMENT_HPP
#define FEMUS_SQUARE_PARABOLA_MOMENT_HPP

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <boost/multiprecision/cpp_bin_float.hpp>

template <class Type> Type SquareParabolaAbs(const Type &value) {
  if(value < Type(0)) return -value;
  return value;
}

template <class Type> Type SquareParabolaPower(const Type &value, const unsigned &exponent) {
  Type result = 1;
  for(unsigned i = 0; i < exponent; i++) result *= value;
  return result;
}

template <class Type> Type SquareParabolaBinomial(const unsigned &n, const unsigned &r) {
  unsigned count = std::min(r, n - r);
  Type result = 1;
  for(unsigned i = 1; i <= count; i++) result = result * Type(n - count + i) / Type(i);
  return result;
}

template <class Type> std::vector<Type> SquareParabolaMultiply(const std::vector<Type> &first, const std::vector<Type> &second) {
  unsigned firstDegree = static_cast<unsigned>(first.size()) - 1;
  unsigned secondDegree = static_cast<unsigned>(second.size()) - 1;
  std::vector<Type> result(firstDegree + secondDegree + 1, Type(0));
  for(unsigned i = 0; i <= firstDegree; i++) {
    for(unsigned j = 0; j <= secondDegree; j++) {
      Type factor = SquareParabolaBinomial<Type>(firstDegree, i) * SquareParabolaBinomial<Type>(secondDegree, j) / SquareParabolaBinomial<Type>(firstDegree + secondDegree, i + j);
      result[i + j] += factor * first[i] * second[j];
    }
  }
  return result;
}

template <class Type> std::vector<Type> SquareParabolaPolynomialPower(const std::vector<Type> &polynomial, const unsigned &exponent) {
  std::vector<Type> result(1, Type(1));
  for(unsigned i = 0; i < exponent; i++) result = SquareParabolaMultiply<Type>(result, polynomial);
  return result;
}

template <class Type> Type SquareParabolaIntegrate(const std::vector<Type> &polynomial, const unsigned &m, const Type &left, const Type &right) {
  std::vector<Type> coordinate = {left, right};
  std::vector<Type> power = SquareParabolaPolynomialPower<Type>(coordinate, m);
  std::vector<Type> product = SquareParabolaMultiply<Type>(power, polynomial);
  Type result = 0;
  for(unsigned i = 0; i < product.size(); i++) result += product[i];
  return (right - left) * result / Type(product.size());
}

template <class Type> void SquareParabolaAddRoot(const Type &root, std::vector<Type> &points) {
  if(root >= Type(0) && root <= Type(1)) points.push_back(root);
}

template <class Type> void SquareParabolaRoots(const Type &k, const Type &b, const Type &d, std::vector<Type> &points) {
  if(k == Type(0)) {
    if(b != Type(0)) SquareParabolaAddRoot<Type>(Type(-d / b), points);
    return;
  }
  if(d == Type(0)) {
    SquareParabolaAddRoot<Type>(Type(0), points);
    SquareParabolaAddRoot<Type>(Type(-b / k), points);
    return;
  }
  if(k + b + d == Type(0)) {
    SquareParabolaAddRoot<Type>(Type(1), points);
    SquareParabolaAddRoot<Type>(Type(d / k), points);
    return;
  }
  if(b * b - Type(4) * k * d == Type(0)) {
    SquareParabolaAddRoot<Type>(Type(-b / (Type(2) * k)), points);
    return;
  }
  Type scale = std::max(SquareParabolaAbs<Type>(k), std::max(SquareParabolaAbs<Type>(b), SquareParabolaAbs<Type>(d)));
  Type aa = k / scale;
  Type bb = b / scale;
  Type dd = d / scale;
  Type discriminant = bb * bb - Type(4) * aa * dd;
  if(discriminant < Type(0)) return;
  if(discriminant == Type(0)) {
    SquareParabolaAddRoot<Type>(Type(-bb / (Type(2) * aa)), points);
    return;
  }
  Type squareRoot = sqrt(discriminant);
  if(bb < Type(0)) squareRoot = -squareRoot;
  Type q = -(bb + squareRoot) / Type(2);
  SquareParabolaAddRoot<Type>(Type(q / aa), points);
  SquareParabolaAddRoot<Type>(Type(dd / q), points);
}

template <class Type> std::vector<Type> SquareParabolaPolynomial(const Type &k, const Type &b, const Type &d, const Type &left, const Type &right) {
  std::vector<Type> result(3, Type(0));
  result[0] = (k * left + b) * left + d;
  result[1] = result[0] + (right - left) * (Type(2) * k * left + b) / Type(2);
  result[2] = (k * right + b) * right + d;
  return result;
}

template <class Type> Type SquareParabolaCanonicalMoment(const int &s, const unsigned &m, const unsigned &n, const Type &k, const Type &b, const Type &c, const Type &d) {
  if(s < -1) throw std::invalid_argument("SquareParabolaMoment requires s >= -1");
  if(k == Type(0) && b == Type(0) && c == Type(0)) {
    if(s == -1) {
      if(d == Type(0)) throw std::domain_error("Delta of an identically zero level set");
      return Type(0);
    }
    Type result = 0;
    if(d > Type(0)) {
      result = 1;
      for(int r = 1; r <= s; r++) result *= d / Type(r);
    }
    else if(d == Type(0) && s == 0) result = Type(1) / Type(2);
    return result / Type(m + 1) / Type(n + 1);
  }
  if(c == Type(0) && s == -1) {
    if(k != Type(0) && b * b - Type(4) * k * d == Type(0)) {
      Type root = -b / (Type(2) * k);
      if(root >= Type(0) && root <= Type(1)) throw std::domain_error("Delta at a repeated level-set root");
      return Type(0);
    }
    std::vector<Type> roots;
    SquareParabolaRoots<Type>(k, b, d, roots);
    std::sort(roots.begin(), roots.end());
    roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
    Type result = 0;
    for(unsigned i = 0; i < roots.size(); i++) {
      Type derivative = Type(2) * k * roots[i] + b;
      if(derivative == Type(0)) throw std::domain_error("Delta at a repeated level-set root");
      Type contribution = SquareParabolaPower<Type>(roots[i], m) / SquareParabolaAbs<Type>(derivative) / Type(n + 1);
      if(roots[i] == Type(0) || roots[i] == Type(1)) contribution /= Type(2);
      result += contribution;
    }
    return result;
  }
  std::vector<Type> points = {Type(0), Type(1)};
  SquareParabolaRoots<Type>(k, b, d, points);
  SquareParabolaRoots<Type>(k, b, Type(d + c), points);
  if(k != Type(0)) SquareParabolaAddRoot<Type>(Type(-b / (Type(2) * k)), points);
  std::sort(points.begin(), points.end());
  points.erase(std::unique(points.begin(), points.end()), points.end());
  Type result = 0;
  unsigned order = (s >= 0) ? static_cast<unsigned>(s) : 0;
  for(unsigned interval = 0; interval + 1 < points.size(); interval++) {
    Type left = points[interval];
    Type right = points[interval + 1];
    Type midpoint = (left + right) / Type(2);
    Type bottom = (k * midpoint + b) * midpoint + d;
    Type top = bottom + c;
    if(s == -1) {
      if((bottom > Type(0) && top > Type(0)) || (bottom < Type(0) && top < Type(0))) continue;
      std::vector<Type> height = SquareParabolaPolynomial<Type>(k, b, d, left, right);
      for(unsigned i = 0; i < height.size(); i++) height[i] = std::max(Type(0), std::min(Type(1), Type(-height[i] / c)));
      std::vector<Type> polynomial = SquareParabolaPolynomialPower<Type>(height, n);
      Type factor = Type(1) / SquareParabolaAbs<Type>(c);
      if(k == Type(0) && b == Type(0) && (d == Type(0) || d + c == Type(0))) factor /= Type(2);
      result += factor * SquareParabolaIntegrate<Type>(polynomial, m, left, right);
      continue;
    }
    if(bottom >= Type(0) && top >= Type(0)) {
      std::vector<Type> lowerPolynomial = SquareParabolaPolynomial<Type>(k, b, d, left, right);
      std::vector<Type> upperPolynomial = SquareParabolaPolynomial<Type>(k, b, Type(d + c), left, right);
      for(unsigned i = 0; i < lowerPolynomial.size(); i++) lowerPolynomial[i] = std::max(Type(0), lowerPolynomial[i]);
      for(unsigned i = 0; i < upperPolynomial.size(); i++) upperPolynomial[i] = std::max(Type(0), upperPolynomial[i]);
      Type factor = 1;
      for(unsigned r = 1; r <= order + 1; r++) factor /= Type(n + r);
      for(unsigned r = 0; r <= order; r++) {
        std::vector<Type> lowerPower = SquareParabolaPolynomialPower<Type>(lowerPolynomial, order - r);
        std::vector<Type> upperPower = SquareParabolaPolynomialPower<Type>(upperPolynomial, r);
        std::vector<Type> polynomial = SquareParabolaMultiply<Type>(lowerPower, upperPower);
        result += factor * SquareParabolaIntegrate<Type>(polynomial, m, left, right);
        factor *= Type(n + r + 1) / Type(r + 1);
      }
      continue;
    }
    if(bottom <= Type(0) && top <= Type(0)) continue;
    std::vector<Type> height = SquareParabolaPolynomial<Type>(k, b, d, left, right);
    std::vector<Type> remainder(3, Type(0));
    for(unsigned i = 0; i < height.size(); i++) {
      height[i] = std::max(Type(0), std::min(Type(1), Type(-height[i] / c)));
      remainder[i] = Type(1) - height[i];
    }
    if(c < Type(0)) {
      Type factor = Type(1) / Type(n + 1);
      for(unsigned r = 1; r <= order; r++) factor *= -c / Type(n + r + 1);
      std::vector<Type> polynomial = SquareParabolaPolynomialPower<Type>(height, n + order + 1);
      result += factor * SquareParabolaIntegrate<Type>(polynomial, m, left, right);
    }
    else {
      Type factor = 1;
      for(unsigned r = 1; r <= order; r++) factor *= c / Type(r);
      for(unsigned r = 0; r <= n; r++) {
        std::vector<Type> lowerPower = SquareParabolaPolynomialPower<Type>(height, n - r);
        std::vector<Type> upperPower = SquareParabolaPolynomialPower<Type>(remainder, order + r + 1);
        std::vector<Type> polynomial = SquareParabolaMultiply<Type>(lowerPower, upperPower);
        Type coefficient = factor * SquareParabolaBinomial<Type>(n, r) / Type(order + r + 1);
        result += coefficient * SquareParabolaIntegrate<Type>(polynomial, m, left, right);
      }
    }
  }
  return result;
}

template <class Type> Type SquareParabolaMoment(const int &s, const std::vector<unsigned> &m, const std::vector<Type> &a, const Type &d, const bool &xSquared) {
  if(m.size() != 2 || a.size() != 3) throw std::invalid_argument("SquareParabolaMoment requires m={i,j} and a={k,b,c}");
  Type k = Type(4) * a[0];
  Type b = Type(2) * a[1] - Type(4) * a[0];
  Type c = Type(2) * a[2];
  Type constant = a[0] - a[1] - a[2] + d;
  unsigned i = xSquared ? m[0] : m[1];
  unsigned j = xSquared ? m[1] : m[0];
  return SquareParabolaCanonicalMoment<Type>(s, i, j, k, b, c, constant);
}

template <class Type> Type SquareParabolaMoment(const int &s, const std::vector<unsigned> &m, const std::vector<Type> &aX, const Type &dX, const std::vector<Type> &aY, const Type &dY, const bool &xSquared) {
  if(xSquared) return SquareParabolaMoment<Type>(s, m, aX, dX, true);
  return SquareParabolaMoment<Type>(s, m, aY, dY, false);
}

#endif
