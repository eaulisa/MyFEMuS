#ifndef FEMUS_ENDPOINT_PROJECTION_HPP
#define FEMUS_ENDPOINT_PROJECTION_HPP

#include "Hex27LevelSetMoments.hpp"

namespace FELevelSetMoments {

struct EndpointProjectionEvent {
  long double point = 0;
  unsigned exponent = 2;

};

class EndpointProjection {
public:
  // Projection polynomials keep every coefficient, including the quartic discriminant.
  static std::vector<std::vector<std::vector<long double>>> Field(const std::vector<long double>& coefficient, unsigned outer, unsigned inner) {
    const unsigned middle = 3U-outer-inner;
    std::vector<std::vector<std::vector<long double>>> field(3,std::vector<std::vector<long double>>(3,std::vector<long double>(3,0)));
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) {
      unsigned powers[3] = {0,0,0};
      powers[inner] = i;
      powers[middle] = j;
      powers[outer] = k;
      field[i][j][k] = coefficient[Hex27Index(powers[0],powers[1],powers[2])];
    }
    return field;
  }

  static std::vector<std::vector<long double>> Derivative(const std::vector<std::vector<long double>>& polynomial) {
    std::vector<std::vector<long double>> result(std::max(1U,static_cast<unsigned>(polynomial.size())-1U),std::vector<long double>(1,0));
    for(unsigned i = 1; i < polynomial.size(); ++i) { result[i-1U] = polynomial[i]; for(unsigned j = 0; j < result[i-1U].size(); ++j) result[i-1U][j] *= i; }
    return result;
  }

  static std::vector<long double> BoundaryValue(const std::vector<std::vector<long double>>& polynomial, bool top, bool triangle) {
    std::vector<long double> result(1,0);
    std::vector<long double> power(1,1);
    for(unsigned i = 0; i < polynomial.size(); ++i) {
      if(top || i == 0) result = PolynomialSum(result,PolynomialProduct(polynomial[i],power));
      if(top && triangle) power = PolynomialProduct(power,{1,-1});
    }
    return result;
  }

  static std::vector<std::vector<long double>> Face(const std::vector<std::vector<std::vector<long double>>>& field, bool top, bool triangle) {
    std::vector<std::vector<long double>> result(5,std::vector<long double>(1,0));
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) {
      if(!top && i != 0) continue;
      for(unsigned a = 0; a <= (top && triangle ? i : 0U); ++a) {
        std::vector<long double> power(1,1);
        for(unsigned k = 0; k < (top && triangle ? i-a : 0U); ++k) power = PolynomialProduct(power,{1,-1});
        long double factor = ReferenceBinomial(i,a)*(a%2U == 0 ? 1.0L : -1.0L);
        if(!triangle) factor = 1;
        result[j+a] = PolynomialSum(result[j+a],PolynomialProduct(field[i][j],power),factor);
      }
    }
    Hex27TrimBivariate(result);
    return result;
  }

  static std::vector<std::vector<long double>> Discriminant(const std::vector<std::vector<std::vector<long double>>>& field) {
    std::vector<std::vector<long double>> result(5,std::vector<long double>(1,0));
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) result[i+j] = PolynomialSum(result[i+j],PolynomialSum(PolynomialProduct(field[1][i],field[1][j]),PolynomialProduct(field[0][i],field[2][j]),-4));
    Hex27TrimBivariate(result);
    return result;
  }

  static std::vector<std::vector<long double>> FaceNormalDerivative(const std::vector<std::vector<std::vector<long double>>>& field, bool top, bool triangle) {
    std::vector<std::vector<std::vector<long double>>> derivative(3,std::vector<std::vector<long double>>(3,std::vector<long double>(3,0)));
    for(unsigned i = 0; i < 2; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) derivative[i][j][k] = (i+1U)*field[i+1U][j][k];
    return Face(derivative,top,triangle);
  }

  // A definite face critical point is isolated when the parent normal derivative is nonzero.
  static bool RegularFaceTangency(const std::vector<std::vector<long double>>& constraint, const std::vector<std::vector<long double>>& normalDerivative, long double point, long double height) {
    long double normal = 0;
    long double normalBound = 0;
    for(unsigned i = 0; i < normalDerivative.size(); ++i) {
      normal += PolynomialValue(normalDerivative[i],height)*ReferencePower(point,i);
      normalBound += PolynomialBound(normalDerivative[i],height)*ReferencePower(std::abs(point),i);
    }
    if(std::abs(normal) <= 65536*std::numeric_limits<long double>::epsilon()*normalBound) return false;
    long double first = 0;
    long double mixed = 0;
    long double second = 0;
    for(unsigned i = 0; i < constraint.size(); ++i) {
      if(i > 1) first += i*(i-1U)*PolynomialValue(constraint[i],height)*ReferencePower(point,i-2U);
      if(i > 0) mixed += i*PolynomialValue(PolynomialDerivative(constraint[i]),height)*ReferencePower(point,i-1U);
      second += PolynomialValue(PolynomialDerivative(PolynomialDerivative(constraint[i])),height)*ReferencePower(point,i);
    }
    const long double determinant = std::fma(first,second,-mixed*mixed);
    return determinant > 65536*std::numeric_limits<long double>::epsilon()*(std::abs(first*second)+mixed*mixed);
  }

  static bool Relevant(const std::vector<std::vector<long double>>& polynomial, long double height, bool triangle) {
    std::vector<long double> fiber(polynomial.size());
    long double norm = 0;
    for(unsigned i = 0; i < polynomial.size(); ++i) { fiber[i] = PolynomialValue(polynomial[i],height); norm += PolynomialBound(polynomial[i],height); }
    const long double upper = triangle ? 1-height : 1;
    const long double padding = 8192*std::numeric_limits<long double>::epsilon()*norm;
    if(std::abs(PolynomialValue(fiber,0)) <= padding || std::abs(PolynomialValue(fiber,upper)) <= padding) return true;
    std::vector<long double> critical;
    PolynomialRoots(PolynomialDerivative(fiber),0,upper,critical);
    for(unsigned i = 0; i < critical.size(); ++i) if(std::abs(PolynomialValue(fiber,critical[i])) <= padding) return true;
    return false;
  }

  static void AddConstraint(const std::vector<std::vector<long double>>& constraint, bool triangle, std::vector<long double>& events, bool regularGraphBoundary = false) {
    std::vector<std::vector<long double>> polynomial = constraint;
    Hex27TrimBivariate(polynomial);
    std::vector<long double> candidates;
    for(unsigned i = 0; i < polynomial.size(); ++i) Hex27AddEventRoots(polynomial[i],candidates);
    Hex27AddEventRoots(BoundaryValue(polynomial,false,triangle),candidates);
    Hex27AddEventRoots(BoundaryValue(polynomial,true,triangle),candidates);
    if(polynomial.size() > 1) {
      const std::vector<long double> resultant = Hex27Resultant(polynomial,Derivative(polynomial));
      if(PolynomialNorm(resultant) == 0) {
        const bool repeatedQuadraticTrace = regularGraphBoundary && polynomial.size() == 3 && PolynomialNorm(PolynomialSum(PolynomialProduct(polynomial[1],polynomial[1]),PolynomialProduct(polynomial[0],polynomial[2]),-4)) == 0;
        if(!repeatedQuadraticTrace) throw std::runtime_error("Unsupported repeated projection factor; endpoint topology is not certified.");
      }
      Hex27AddEventRoots(resultant,candidates);
    }
    for(unsigned i = 0; i < candidates.size(); ++i) if(Relevant(polynomial,candidates[i],triangle)) events.push_back(candidates[i]);
  }

  static unsigned CommonExponent(unsigned first, unsigned second) {
    unsigned a = first;
    unsigned b = second;
    while(b != 0) { const unsigned remainder = a%b; a = b; b = remainder; }
    return first/a*second;
  }

  static unsigned ContactExponent(const std::vector<std::vector<long double>>& constraint, long double height, bool triangle, bool regularGraphBoundary = false, const std::vector<std::vector<long double>>* normalDerivative = nullptr) {
    std::vector<long double> polynomial(constraint.size());
    for(unsigned i = 0; i < constraint.size(); ++i) polynomial[i] = PolynomialValue(constraint[i],height);
    Trim(polynomial);
    if(PolynomialNorm(polynomial) == 0) return 1;
    std::vector<long double> roots;
    PolynomialRoots(polynomial,0,triangle ? 1-height : 1,roots);
    std::vector<long double> stationary;
    PolynomialRoots(PolynomialDerivative(polynomial),0,triangle ? 1-height : 1,stationary);
    for(unsigned n = 0; n < stationary.size(); ++n) if(std::abs(PolynomialValue(polynomial,stationary[n])) <= 65536*std::numeric_limits<long double>::epsilon()*PolynomialBound(polynomial,stationary[n])) roots.push_back(stationary[n]);
    SortDistinct(roots);
    unsigned exponent = 1;
    for(unsigned n = 0; n < roots.size(); ++n) {
      std::vector<long double> derivative = PolynomialDerivative(polynomial);
      unsigned multiplicity = 1;
      while(derivative.size() > 1 && std::abs(PolynomialValue(derivative,roots[n])) <= 65536*std::numeric_limits<long double>::epsilon()*PolynomialBound(derivative,roots[n])) { ++multiplicity; derivative = PolynomialDerivative(derivative); }
      if(multiplicity < 2) continue;
      long double transverse = 0;
      for(unsigned i = 0; i < constraint.size(); ++i) transverse += PolynomialValue(PolynomialDerivative(constraint[i]),height)*ReferencePower(roots[n],i);
      long double norm = 0;
      for(unsigned i = 0; i < constraint.size(); ++i) norm += PolynomialNorm(constraint[i]);
      if(!regularGraphBoundary && std::abs(transverse) <= 65536*std::numeric_limits<long double>::epsilon()*norm && (normalDerivative == nullptr || !RegularFaceTangency(constraint,*normalDerivative,roots[n],height))) throw std::runtime_error("Unsupported non-transverse projected endpoint; no regular-fold convergence claim is available.");
      if(multiplicity > 4) throw std::runtime_error("Unsupported projected contact multiplicity.");
      exponent = CommonExponent(exponent,CommonExponent(2,multiplicity));
    }
    return exponent;
  }

  static std::vector<EndpointProjectionEvent> BuildDirection(const std::vector<long double>& coefficient, bool triangle, unsigned outer, unsigned certifiedInner, bool graphOnly) {
    std::vector<long double> events = {0,1};
    std::vector<std::vector<std::vector<long double>>> constraints;
    std::vector<std::vector<std::vector<long double>>> normalDerivatives;
    for(unsigned inner = 0; inner < 3; ++inner) if(inner != outer && (certifiedInner == 3 || inner == certifiedInner)) {
      const std::vector<std::vector<std::vector<long double>>> field = Field(coefficient,outer,inner);
      for(unsigned side = 0; side < 2; ++side) {
        const std::vector<std::vector<long double>> face = Face(field,side != 0,triangle);
        AddConstraint(face,triangle,events,graphOnly);
        constraints.push_back(face);
        normalDerivatives.push_back(FaceNormalDerivative(field,side != 0,triangle));
      }
      if(graphOnly) continue;
      bool quadratic = false;
      bool linear = false;
      for(unsigned i = 0; i < 3; ++i) { quadratic = quadratic || PolynomialNorm(field[2][i]) != 0; linear = linear || PolynomialNorm(field[1][i]) != 0; }
      const std::vector<std::vector<long double>> discriminant = quadratic ? Discriminant(field) : (linear ? field[1] : field[0]);
      if(!quadratic && linear) Hex27AddEventRoots(Hex27Resultant(field[0],field[1]),events);
      if(quadratic && discriminant.size() == 1 && PolynomialNorm(discriminant[0]) == 0) throw std::runtime_error("A repeated fiber root has no regular delta moment.");
      AddConstraint(discriminant,triangle,events);
      constraints.push_back(discriminant);
      normalDerivatives.push_back({});
    }
    SortDistinct(events);
    std::vector<EndpointProjectionEvent> result;
    for(unsigned n = 0; n < events.size(); ++n) {
      unsigned exponent = graphOnly ? 1U : 2U;
      for(unsigned i = 0; i < constraints.size(); ++i) exponent = CommonExponent(exponent,ContactExponent(constraints[i],events[n],triangle,graphOnly,normalDerivatives[i].empty() ? nullptr : &normalDerivatives[i]));
      result.push_back({events[n],exponent});
    }
    return result;
  }
  static std::vector<EndpointProjectionEvent> Build(const std::vector<long double>& coefficient, bool triangle, unsigned outer, unsigned inner = 3, bool graphOnly = false) {
    std::runtime_error failure("No supported endpoint projection direction.");
    for(unsigned direction = 0; direction < 3; ++direction) if(direction != outer && (inner == 3 || inner == direction)) {
      try { return BuildDirection(coefficient,triangle,outer,direction,graphOnly); }
      catch(const std::runtime_error& error) { failure = std::runtime_error(error.what()); }
    }
    throw failure;
  }


};

}

#endif
