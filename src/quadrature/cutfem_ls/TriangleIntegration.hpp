#ifndef FEMUS_TRIANGLE_INTEGRATION_HPP
#define FEMUS_TRIANGLE_INTEGRATION_HPP

#include "LineIntegration.hpp"

namespace FELevelSetMoments {

struct TriangleDomainVertex {
  long double x;
  long double y;
};

struct TriangleDomainStrip {
  long double left;
  long double right;
  long double lowerConstant;
  long double lowerSlope;
  long double upperConstant;
  long double upperSlope;
};

class TriangleIntegration {
public:
  // Integrates the unchanged parent level-set polynomial over an arbitrary triangular subregion.
  static void Integrate(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const std::vector<TriangleDomainVertex>& vertices, bool computeDelta, long double absoluteTolerance, long double absoluteDeltaTolerance, long double relativeTolerance, unsigned maximumDepth, unsigned maximumEvaluations, unsigned splitDepth, unsigned maximumSplitDepth, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned& selectedDirection, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
    if(vertices.size() != 3 || coefficient.size() != 3) throw std::invalid_argument("Triangle integration requires three vertices and a 3 by 3 coefficient tensor.");
    for(unsigned i = 0; i < 3; ++i) if(coefficient[i].size() != 3) throw std::invalid_argument("Triangle integration requires a 3 by 3 coefficient tensor.");
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    values.assign(2U*count,0.0L);
    errors.assign(2U*count,0.0L);
    long double minimum = 0.0L;
    long double maximum = 0.0L;
    PolynomialBounds(coefficient,vertices,minimum,maximum);
    const long double norm = PolynomialNorm(coefficient);
    const long double padding = 1024.0L*std::numeric_limits<long double>::epsilon()*std::max(1.0L,norm);
    if(minimum > padding || maximum < -padding) {
      if(minimum > padding) FillFullTriangle(vertices,degree,values,errors);
      ++accepted;
      selectedDirection = 0;
      return;
    }
    long double lower[2] = {0.0L,0.0L};
    long double upper[2] = {0.0L,0.0L};
    DerivativeBounds(coefficient,vertices,0,lower[0],upper[0]);
    DerivativeBounds(coefficient,vertices,1,lower[1],upper[1]);
    unsigned direction = lower[1] > lower[0] || (lower[1] == lower[0] && upper[1] > upper[0]) ? 1U : 2U;
    const long double selectedLower = direction == 1 ? lower[1] : lower[0];
    if(selectedLower == 0.0L && splitDepth < maximumSplitDepth) {
      std::vector<std::vector<TriangleDomainVertex>> children;
      if(!SplitAtDerivativeChange(coefficient,vertices,children)) {
        std::vector<TriangleDomainVertex> first;
        std::vector<TriangleDomainVertex> second;
        SplitLongestEdge(vertices,first,second);
        children.push_back(first);
        children.push_back(second);
      }
      std::vector<long double> correction(values.size(),0.0L);
      unsigned commonDirection = 0;
      for(unsigned child = 0; child < children.size(); ++child) {
        std::vector<long double> childValues;
        std::vector<long double> childErrors;
        unsigned childDirection = 0;
        const long double fraction = 1.0L/static_cast<long double>(children.size());
        Integrate(coefficient,degree,children[child],computeDelta,fraction*absoluteTolerance,fraction*absoluteDeltaTolerance,relativeTolerance,maximumDepth,maximumEvaluations,splitDepth+1U,maximumSplitDepth,childValues,childErrors,evaluations,accepted,childDirection,surfaceGradient);
        if(child == 0) commonDirection = childDirection;
        else if(commonDirection != childDirection) commonDirection = 3U;
        for(unsigned k = 0; k < values.size(); ++k) {
          const long double adjusted = childValues[k]-correction[k];
          const long double next = values[k]+adjusted;
          correction[k] = (next-values[k])-adjusted;
          values[k] = next;
          errors[k] += childErrors[k]+8.0L*std::numeric_limits<long double>::epsilon()*std::abs(childValues[k]);
        }
      }
      selectedDirection = commonDirection;
      CheckErrors(values,errors,degree,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance);
      return;
    }
    selectedDirection = direction;
    if(selectedLower > 0.0L) {
      try {
        IntegrateDirection(coefficient,degree,vertices,direction,computeDelta,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,surfaceGradient,true,true);
      }
      catch(const std::runtime_error&) {
        IntegrateDirection(coefficient,degree,vertices,direction,computeDelta,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,surfaceGradient,false,false);
      }
    }
    else IntegrateDirection(coefficient,degree,vertices,direction,computeDelta,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,surfaceGradient,false,false);
    CheckErrors(values,errors,degree,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance);
  }

private:
  static long double Power(long double value, unsigned degree) {
    long double result = 1.0L;
    for(unsigned k = 0; k < degree; ++k) result *= value;
    return result;
  }

  static long double PolynomialNorm(const std::vector<std::vector<long double>>& coefficient) {
    long double result = 0.0L;
    for(unsigned i = 0; i < coefficient.size(); ++i) for(unsigned j = 0; j < coefficient[i].size(); ++j) result += std::abs(coefficient[i][j]);
    return result;
  }

  static void PowerBounds(long double lower, long double upper, unsigned degree, long double& minimum, long double& maximum) {
    if(degree == 0) {
      minimum = 1.0L;
      maximum = 1.0L;
      return;
    }
    long double a = Power(lower,degree);
    long double b = Power(upper,degree);
    minimum = std::min(a,b);
    maximum = std::max(a,b);
    if((degree&1U) == 0 && lower <= 0.0L && upper >= 0.0L) minimum = 0.0L;
  }

  static void BoundingBox(const std::vector<TriangleDomainVertex>& vertices, long double& xmin, long double& xmax, long double& ymin, long double& ymax) {
    xmin = xmax = vertices[0].x;
    ymin = ymax = vertices[0].y;
    for(unsigned k = 1; k < 3; ++k) {
      xmin = std::min(xmin,vertices[k].x);
      xmax = std::max(xmax,vertices[k].x);
      ymin = std::min(ymin,vertices[k].y);
      ymax = std::max(ymax,vertices[k].y);
    }
  }

  static void PolynomialBounds(const std::vector<std::vector<long double>>& coefficient, const std::vector<TriangleDomainVertex>& vertices, long double& minimum, long double& maximum) {
    long double xmin;
    long double xmax;
    long double ymin;
    long double ymax;
    BoundingBox(vertices,xmin,xmax,ymin,ymax);
    minimum = 0.0L;
    maximum = 0.0L;
    for(unsigned i = 0; i < 3; ++i) {
      long double xminPower;
      long double xmaxPower;
      PowerBounds(xmin,xmax,i,xminPower,xmaxPower);
      for(unsigned j = 0; j < 3; ++j) {
        if(coefficient[i][j] == 0.0L) continue;
        long double yminPower;
        long double ymaxPower;
        PowerBounds(ymin,ymax,j,yminPower,ymaxPower);
        const long double products[4] = {xminPower*yminPower,xminPower*ymaxPower,xmaxPower*yminPower,xmaxPower*ymaxPower};
        long double termMinimum = products[0];
        long double termMaximum = products[0];
        for(unsigned k = 1; k < 4; ++k) {
          termMinimum = std::min(termMinimum,products[k]);
          termMaximum = std::max(termMaximum,products[k]);
        }
        if(coefficient[i][j] >= 0.0L) {
          minimum += coefficient[i][j]*termMinimum;
          maximum += coefficient[i][j]*termMaximum;
        }
        else {
          minimum += coefficient[i][j]*termMaximum;
          maximum += coefficient[i][j]*termMinimum;
        }
      }
    }
  }

  static void DerivativeBounds(const std::vector<std::vector<long double>>& coefficient, const std::vector<TriangleDomainVertex>& vertices, unsigned axis, long double& lower, long double& upper) {
    std::vector<std::vector<long double>> derivative(3,std::vector<long double>(3,0.0L));
    if(axis == 0) {
      for(unsigned i = 1; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) derivative[i-1U][j] = i*coefficient[i][j];
    }
    else {
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = 1; j < 3; ++j) derivative[i][j-1U] = j*coefficient[i][j];
    }
    long double minimum = 0.0L;
    long double maximum = 0.0L;
    PolynomialBounds(derivative,vertices,minimum,maximum);
    const long double magnitude = std::max(std::abs(minimum),std::abs(maximum));
    const long double padding = 64.0L*std::sqrt(std::numeric_limits<long double>::epsilon())*std::max(1.0L,magnitude);
    lower = minimum > padding ? minimum-padding : (maximum < -padding ? -maximum-padding : 0.0L);
    upper = magnitude+padding;
  }

  static long double Coordinate(const TriangleDomainVertex& vertex, unsigned axis) {
    return axis == 0 ? vertex.x : vertex.y;
  }

  static long double DerivativeValue(const std::vector<std::vector<long double>>& coefficient, unsigned axis, const TriangleDomainVertex& vertex) {
    long double value = 0.0L;
    if(axis == 0) {
      for(unsigned i = 1; i < 3; ++i) {
        for(unsigned j = 0; j < 3; ++j) value += static_cast<long double>(i)*coefficient[i][j]*Power(vertex.x,i-1U)*Power(vertex.y,j);
      }
    }
    else {
      for(unsigned i = 0; i < 3; ++i) {
        for(unsigned j = 1; j < 3; ++j) value += static_cast<long double>(j)*coefficient[i][j]*Power(vertex.x,i)*Power(vertex.y,j-1U);
      }
    }
    return value;
  }

  static TriangleDomainVertex InterpolateToCoordinate(const TriangleDomainVertex& first, const TriangleDomainVertex& second, unsigned axis, long double coordinate) {
    const long double firstCoordinate = Coordinate(first,axis);
    const long double secondCoordinate = Coordinate(second,axis);
    const long double parameter = (coordinate-firstCoordinate)/(secondCoordinate-firstCoordinate);
    TriangleDomainVertex point;
    point.x = first.x+parameter*(second.x-first.x);
    point.y = first.y+parameter*(second.y-first.y);
    return point;
  }

  static void ClipByCoordinate(const std::vector<TriangleDomainVertex>& input, unsigned axis, long double coordinate, bool lowerSide, std::vector<TriangleDomainVertex>& output) {
    output.clear();
    if(input.empty()) return;
    TriangleDomainVertex previous = input.back();
    bool previousInside = lowerSide ? Coordinate(previous,axis) <= coordinate : Coordinate(previous,axis) >= coordinate;
    for(unsigned index = 0; index < input.size(); ++index) {
      const TriangleDomainVertex current = input[index];
      const bool currentInside = lowerSide ? Coordinate(current,axis) <= coordinate : Coordinate(current,axis) >= coordinate;
      if(currentInside != previousInside && Coordinate(current,axis) != Coordinate(previous,axis)) output.push_back(InterpolateToCoordinate(previous,current,axis,coordinate));
      if(currentInside) output.push_back(current);
      previous = current;
      previousInside = currentInside;
    }
  }

  static void TriangulatePolygon(const std::vector<TriangleDomainVertex>& polygon, std::vector<std::vector<TriangleDomainVertex>>& triangles) {
    if(polygon.size() < 3) return;
    for(unsigned index = 1; index+1U < polygon.size(); ++index) {
      std::vector<TriangleDomainVertex> triangle = {polygon[0],polygon[index],polygon[index+1U]};
      const long double determinant = (triangle[1].x-triangle[0].x)*(triangle[2].y-triangle[0].y)-(triangle[2].x-triangle[0].x)*(triangle[1].y-triangle[0].y);
      if(std::abs(determinant) > 256.0L*std::numeric_limits<long double>::epsilon()) triangles.push_back(triangle);
    }
  }

  static bool SplitAtDerivativeChange(const std::vector<std::vector<long double>>& coefficient, const std::vector<TriangleDomainVertex>& vertices, std::vector<std::vector<TriangleDomainVertex>>& children) {
    unsigned bestAxis = 2;
    long double bestCoordinate = 0.0L;
    long double bestScore = std::numeric_limits<long double>::infinity();
    for(unsigned axis = 0; axis < 2; ++axis) {
      long double derivative[3];
      for(unsigned vertex = 0; vertex < 3; ++vertex) derivative[vertex] = DerivativeValue(coefficient,axis,vertices[vertex]);
      std::vector<long double> crossings;
      for(unsigned edge = 0; edge < 3; ++edge) {
        const unsigned next = (edge+1U)%3U;
        if(derivative[edge] == 0.0L) crossings.push_back(Coordinate(vertices[edge],axis));
        if(derivative[edge]*derivative[next] < 0.0L) {
          const long double parameter = derivative[edge]/(derivative[edge]-derivative[next]);
          crossings.push_back(Coordinate(vertices[edge],axis)+parameter*(Coordinate(vertices[next],axis)-Coordinate(vertices[edge],axis)));
        }
      }
      if(crossings.size() < 2) continue;
      long double minimumCoordinate = Coordinate(vertices[0],axis);
      long double maximumCoordinate = minimumCoordinate;
      for(unsigned vertex = 1; vertex < 3; ++vertex) {
        minimumCoordinate = std::min(minimumCoordinate,Coordinate(vertices[vertex],axis));
        maximumCoordinate = std::max(maximumCoordinate,Coordinate(vertices[vertex],axis));
      }
      const long double width = maximumCoordinate-minimumCoordinate;
      if(!(width > 0.0L)) continue;
      long double coordinate = 0.0L;
      for(unsigned crossing = 0; crossing < crossings.size(); ++crossing) coordinate += crossings[crossing];
      coordinate /= static_cast<long double>(crossings.size());
      const long double margin = 0.05L*width;
      coordinate = std::max(minimumCoordinate+margin,std::min(maximumCoordinate-margin,coordinate));
      long double spread = 0.0L;
      for(unsigned crossing = 0; crossing < crossings.size(); ++crossing) spread = std::max(spread,std::abs(crossings[crossing]-coordinate));
      const long double score = spread/width;
      if(score < bestScore) {
        bestScore = score;
        bestAxis = axis;
        bestCoordinate = coordinate;
      }
    }
    if(bestAxis > 1) return false;
    std::vector<TriangleDomainVertex> lowerPolygon;
    std::vector<TriangleDomainVertex> upperPolygon;
    ClipByCoordinate(vertices,bestAxis,bestCoordinate,true,lowerPolygon);
    ClipByCoordinate(vertices,bestAxis,bestCoordinate,false,upperPolygon);
    children.clear();
    TriangulatePolygon(lowerPolygon,children);
    TriangulatePolygon(upperPolygon,children);
    return children.size() >= 2;
  }

  static void SplitLongestEdge(const std::vector<TriangleDomainVertex>& vertices, std::vector<TriangleDomainVertex>& first, std::vector<TriangleDomainVertex>& second) {
    unsigned edgeA = 0;
    unsigned edgeB = 1;
    unsigned other = 2;
    long double longest = -1.0L;
    for(unsigned a = 0; a < 3; ++a) {
      const unsigned b = (a+1U)%3U;
      const unsigned c = 3U-a-b;
      const long double dx = vertices[a].x-vertices[b].x;
      const long double dy = vertices[a].y-vertices[b].y;
      const long double length = dx*dx+dy*dy;
      if(length > longest) {
        longest = length;
        edgeA = a;
        edgeB = b;
        other = c;
      }
    }
    TriangleDomainVertex middle = {0.5L*(vertices[edgeA].x+vertices[edgeB].x),0.5L*(vertices[edgeA].y+vertices[edgeB].y)};
    first = {vertices[edgeA],middle,vertices[other]};
    second = {middle,vertices[edgeB],vertices[other]};
  }

  static long double TriangleMoment(const std::vector<TriangleDomainVertex>& vertices, unsigned firstDegree, unsigned secondDegree) {
    const long double ax = vertices[1].x-vertices[0].x;
    const long double bx = vertices[2].x-vertices[0].x;
    const long double ay = vertices[1].y-vertices[0].y;
    const long double by = vertices[2].y-vertices[0].y;
    const long double determinant = std::abs(ax*by-bx*ay);
    long double result = 0.0L;
    for(unsigned a = 0; a <= firstDegree; ++a) {
      for(unsigned b = 0; b <= firstDegree-a; ++b) {
        const unsigned constantX = firstDegree-a-b;
        const long double firstMultiplier = LineSliceUtilities::Binomial(firstDegree,a)*LineSliceUtilities::Binomial(firstDegree-a,b)*Power(vertices[0].x,constantX)*Power(ax,a)*Power(bx,b);
        for(unsigned c = 0; c <= secondDegree; ++c) {
          for(unsigned d = 0; d <= secondDegree-c; ++d) {
            const unsigned constantY = secondDegree-c-d;
            const long double secondMultiplier = LineSliceUtilities::Binomial(secondDegree,c)*LineSliceUtilities::Binomial(secondDegree-c,d)*Power(vertices[0].y,constantY)*Power(ay,c)*Power(by,d);
            result += firstMultiplier*secondMultiplier*ReferenceTriangleMoment(a+c,b+d);
          }
        }
      }
    }
    return determinant*result;
  }

  static void FillFullTriangle(const std::vector<TriangleDomainVertex>& vertices, unsigned degree, std::vector<long double>& values, std::vector<long double>& errors) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned index = LineSliceUtilities::MomentIndex(degree,i,j);
        values[index] = TriangleMoment(vertices,i,j);
        errors[index] = 32.0L*std::numeric_limits<long double>::epsilon()*std::abs(values[index]);
      }
    }
    for(unsigned k = count; k < 2U*count; ++k) {
      values[k] = 0.0L;
      errors[k] = 0.0L;
    }
  }

  static void LineEquation(const TriangleDomainVertex& first, const TriangleDomainVertex& second, long double& constant, long double& slope) {
    if(first.x == second.x) throw std::runtime_error("A vertical edge cannot define an inner strip boundary.");
    slope = (second.y-first.y)/(second.x-first.x);
    constant = first.y-slope*first.x;
  }

  static void AddStrip(const TriangleDomainVertex& first, const TriangleDomainVertex& second, const TriangleDomainVertex& longFirst, const TriangleDomainVertex& longSecond, std::vector<TriangleDomainStrip>& strips) {
    const long double left = std::min(first.x,second.x);
    const long double right = std::max(first.x,second.x);
    if(!(right > left)) return;
    long double firstConstant;
    long double firstSlope;
    long double longConstant;
    long double longSlope;
    LineEquation(first,second,firstConstant,firstSlope);
    LineEquation(longFirst,longSecond,longConstant,longSlope);
    const long double middle = 0.5L*(left+right);
    const long double firstValue = firstConstant+firstSlope*middle;
    const long double longValue = longConstant+longSlope*middle;
    TriangleDomainStrip strip;
    strip.left = left;
    strip.right = right;
    if(firstValue <= longValue) {
      strip.lowerConstant = firstConstant;
      strip.lowerSlope = firstSlope;
      strip.upperConstant = longConstant;
      strip.upperSlope = longSlope;
    }
    else {
      strip.lowerConstant = longConstant;
      strip.lowerSlope = longSlope;
      strip.upperConstant = firstConstant;
      strip.upperSlope = firstSlope;
    }
    strips.push_back(strip);
  }

  static bool VertexLess(const TriangleDomainVertex& first, const TriangleDomainVertex& second) {
    if(first.x != second.x) return first.x < second.x;
    return first.y < second.y;
  }

  static void BuildStrips(const std::vector<TriangleDomainVertex>& input, std::vector<TriangleDomainStrip>& strips) {
    std::vector<TriangleDomainVertex> vertices(input);
    std::sort(vertices.begin(),vertices.end(),VertexLess);
    strips.clear();
    if(vertices[0].x == vertices[2].x) return;
    if(vertices[0].x != vertices[1].x) AddStrip(vertices[0],vertices[1],vertices[0],vertices[2],strips);
    if(vertices[1].x != vertices[2].x) AddStrip(vertices[1],vertices[2],vertices[0],vertices[2],strips);
  }

  static std::vector<long double> BoundaryPolynomial(const std::vector<std::vector<long double>>& coefficient, long double constant, long double slope) {
    std::vector<long double> result(5,0.0L);
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) {
        for(unsigned k = 0; k <= j; ++k) result[i+k] += coefficient[i][j]*LineSliceUtilities::Binomial(j,k)*Power(constant,j-k)*Power(slope,k);
      }
    }
    LineSliceUtilities::Trim(result);
    return result;
  }

  static long double Evaluate(const std::vector<std::vector<long double>>& coefficient, long double x, long double y) {
    long double value = 0.0L;
    for(unsigned i = 3; i > 0; --i) {
      value *= x;
      long double row = 0.0L;
      for(unsigned j = 3; j > 0; --j) row = std::fma(row,y,coefficient[i-1U][j-1U]);
      value += row;
    }
    return value;
  }

  static long double StripLower(const TriangleDomainStrip& strip, long double x) {
    return std::fma(strip.lowerSlope,x,strip.lowerConstant);
  }

  static long double StripUpper(const TriangleDomainStrip& strip, long double x) {
    return std::fma(strip.upperSlope,x,strip.upperConstant);
  }

  static std::vector<long double> FiberPolynomial(const std::vector<std::vector<long double>>& coefficient, long double x) {
    std::vector<long double> polynomial(3,0.0L);
    long double power = 1.0L;
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) polynomial[j] += coefficient[i][j]*power;
      power *= x;
    }
    return polynomial;
  }

  static long double StripMoment(const TriangleDomainStrip& strip, long double left, long double right, unsigned i, unsigned j) {
    std::vector<long double> upper(j+2U,0.0L);
    std::vector<long double> lower(j+2U,0.0L);
    for(unsigned k = 0; k <= j+1U; ++k) {
      upper[k] = LineSliceUtilities::Binomial(j+1U,k)*Power(strip.upperConstant,j+1U-k)*Power(strip.upperSlope,k);
      lower[k] = LineSliceUtilities::Binomial(j+1U,k)*Power(strip.lowerConstant,j+1U-k)*Power(strip.lowerSlope,k);
    }
    long double result = 0.0L;
    for(unsigned k = 0; k < upper.size(); ++k) result += (upper[k]-lower[k])*LineSliceUtilities::MonomialIntegral(left,right,i+k)/(j+1U);
    return result;
  }

  static void CollectEvents(const std::vector<std::vector<long double>>& coefficient, const TriangleDomainStrip& strip, std::vector<long double>& events) {
    events = {strip.left,strip.right};
    std::vector<long double> lower = BoundaryPolynomial(coefficient,strip.lowerConstant,strip.lowerSlope);
    std::vector<long double> upper = BoundaryPolynomial(coefficient,strip.upperConstant,strip.upperSlope);
    LineSliceUtilities::PolynomialRoots(lower,strip.left,strip.right,events);
    LineSliceUtilities::PolynomialRoots(upper,strip.left,strip.right,events);
    std::vector<long double> a(3,0.0L);
    std::vector<long double> b(3,0.0L);
    std::vector<long double> c(3,0.0L);
    for(unsigned i = 0; i < 3; ++i) {
      a[i] = coefficient[i][2];
      b[i] = coefficient[i][1];
      c[i] = coefficient[i][0];
    }
    std::vector<long double> discriminant = LineSliceUtilities::PolynomialSum(LineSliceUtilities::PolynomialProduct(b,b),LineSliceUtilities::PolynomialProduct(a,c),-4.0L);
    LineSliceUtilities::PolynomialRoots(discriminant,strip.left,strip.right,events);
    LineSliceUtilities::SortDistinct(events);
  }

  static void GraphSample(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const TriangleDomainStrip& strip, long double x, bool computeDelta, std::vector<long double>& values, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>* active) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    values.assign(2U*count,0.0L);
    const long double lower = StripLower(strip,x);
    const long double upper = StripUpper(strip,x);
    std::vector<long double> polynomial = FiberPolynomial(coefficient,x);
    std::vector<long double> roots = LineSliceUtilities::QuadraticRoots(polynomial[2],polynomial[1],polynomial[0]);
    long double root = 0.0L;
    unsigned number = 0;
    const long double tolerance = 256.0L*std::numeric_limits<long double>::epsilon()*std::max(1.0L,std::max(std::abs(lower),std::abs(upper)));
    for(unsigned k = 0; k < roots.size(); ++k) {
      if(roots[k] >= lower-tolerance && roots[k] <= upper+tolerance) {
        root = std::max(lower,std::min(upper,roots[k]));
        ++number;
      }
    }
    if(number > 1) throw std::runtime_error("A local triangle graph fiber has multiple roots.");
    long double positiveLeft = lower;
    long double positiveRight = lower;
    long double derivative = 0.0L;
    if(number == 1) {
      derivative = std::abs(std::fma(2.0L*polynomial[2],root,polynomial[1]));
      if(!(derivative > 0.0L)) throw std::runtime_error("A local triangle graph has a singular fiber derivative.");
      const long double signedDerivative = std::fma(2.0L*polynomial[2],root,polynomial[1]);
      positiveLeft = signedDerivative > 0.0L ? root : lower;
      positiveRight = signedDerivative > 0.0L ? upper : root;
    }
    else {
      const long double middle = 0.5L*(lower+upper);
      if(Evaluate(coefficient,x,middle) > 0.0L) {
        positiveLeft = lower;
        positiveRight = upper;
      }
    }
    long double xPower = 1.0L;
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned index = LineSliceUtilities::MomentIndex(degree,i,j);
        if(active == nullptr || (*active)[index] != 0) values[index] = xPower*LineSliceUtilities::MonomialIntegral(positiveLeft,positiveRight,j);
        if(computeDelta && number == 1 && (active == nullptr || (*active)[count+index] != 0)) {
          long double weight = 1.0L/derivative;
          if(surfaceGradient != nullptr) weight *= SurfaceGradientNorm(*surfaceGradient,x,root);
          values[count+index] = xPower*Power(root,j)*weight;
        }
      }
      xPower *= x;
    }
  }

  static void GraphGaussKronrodInterval(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const TriangleDomainStrip& strip, long double left, long double right, bool computeDelta, long double absoluteTolerance, long double absoluteDeltaTolerance, long double relativeTolerance, unsigned depth, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, std::vector<long double>& total, std::vector<long double>& totalError, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>* active = nullptr) {
    static const long double nodes[8] = {0.991455371120812639206854697526329L,0.949107912342758524526189684047851L,0.864864423359769072789712788640926L,0.741531185599394439863864773280788L,0.586087235467691130294144838258730L,0.405845151377397166906606412076961L,0.207784955007898467600689403773245L,0.0L};
    static const long double kronrod[8] = {0.022935322010529224963732008058970L,0.063092092629978553290700663189204L,0.104790010322250183839876322541518L,0.140653259715525918745189590510238L,0.169004726639267902826583426598550L,0.190350578064785409913256402421014L,0.204432940075298892414161999234649L,0.209482141084727828012999174891714L};
    static const long double gauss[4] = {0.129484966168869693270611432679082L,0.279705391489276667901467771423780L,0.381830050505118944950369775488975L,0.417959183673469387755102040816327L};
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    const unsigned size = 2U*count;
    const long double middle = 0.5L*(left+right);
    const long double halfWidth = 0.5L*(right-left);
    std::vector<long double> high(size,0.0L);
    std::vector<long double> low(size,0.0L);
    std::vector<long double> absolute(size,0.0L);
    std::vector<std::vector<long double>> samples(15);
    unsigned sampleNumber = 0;
    for(unsigned node = 0; node < 8; ++node) {
      const unsigned sides = node == 7 ? 1U : 2U;
      for(unsigned side = 0; side < sides; ++side) {
        if(evaluations >= maximumEvaluations) throw std::runtime_error("Triangle graph evaluation budget exhausted.");
        ++evaluations;
        const long double point = middle+(side == 0 ? -1.0L : 1.0L)*halfWidth*nodes[node];
        GraphSample(coefficient,degree,strip,point,computeDelta,samples[sampleNumber],surfaceGradient,active);
        for(unsigned k = 0; k < size; ++k) {
          if(active != nullptr && (*active)[k] == 0) continue;
          high[k] += kronrod[node]*samples[sampleNumber][k];
          absolute[k] += kronrod[node]*std::abs(samples[sampleNumber][k]);
          if((node&1U) != 0) low[k] += gauss[node/2U]*samples[sampleNumber][k];
        }
        ++sampleNumber;
      }
    }
    std::vector<unsigned char> failed(size,0);
    bool anyFailed = false;
    for(unsigned k = 0; k < size; ++k) {
      if(active != nullptr && (*active)[k] == 0) continue;
      long double variation = 0.0L;
      for(unsigned n = 0; n < sampleNumber; ++n) variation += kronrod[n/2U]*std::abs(samples[n][k]-0.5L*high[k]);
      long double estimate = std::abs(high[k]-low[k]);
      if(variation > 0.0L && estimate > 0.0L) estimate = variation*std::min(1.0L,std::pow(200.0L*estimate/variation,1.5L));
      const long double error = halfWidth*std::max(estimate,128.0L*std::numeric_limits<long double>::epsilon()*absolute[k]);
      const long double value = halfWidth*high[k];
      const long double tolerance = (k < count ? absoluteTolerance : absoluteDeltaTolerance)*(right-left)+relativeTolerance*std::abs(value);
      if(error <= tolerance) {
        total[k] += value;
        totalError[k] += error;
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
    if(depth >= maximumDepth || middle == left || middle == right) throw std::runtime_error("Triangle fast graph GK did not reach tolerance.");
    GraphGaussKronrodInterval(coefficient,degree,strip,left,middle,computeDelta,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError,surfaceGradient,&failed);
    GraphGaussKronrodInterval(coefficient,degree,strip,middle,right,computeDelta,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError,surfaceGradient,&failed);
  }

  static void GraphTaylorInterval(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const TriangleDomainStrip& strip, long double left, long double right, bool computeDelta, long double absoluteTolerance, long double absoluteDeltaTolerance, long double relativeTolerance, unsigned depth, unsigned maximumDepth, std::vector<long double>& total, std::vector<long double>& totalError, unsigned& accepted, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>* active = nullptr) {
    const unsigned order = 32;
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    const unsigned size = 2U*count;
    const long double thetaLeft = 0.0L;
    const long double thetaRight = 0.5L*std::acos(-1.0L);
    const long double thetaMiddle = 0.5L*(thetaLeft+thetaRight);
    const long double thetaHalfWidth = 0.5L*(thetaRight-thetaLeft);
    std::vector<long double> sine;
    std::vector<long double> cosine;
    SineCosineSeries(thetaMiddle,thetaHalfWidth,order,sine,cosine);
    std::vector<long double> x = SeriesProduct(sine,sine);
    const long double width = right-left;
    for(unsigned k = 0; k <= order; ++k) x[k] *= width;
    x[0] += left;
    std::vector<long double> jacobian = SeriesProduct(sine,cosine);
    for(unsigned k = 0; k <= order; ++k) jacobian[k] *= 2.0L*width*thetaHalfWidth;
    std::vector<long double> aPolynomial(3,0.0L);
    std::vector<long double> bPolynomial(3,0.0L);
    std::vector<long double> cPolynomial(3,0.0L);
    for(unsigned i = 0; i < 3; ++i) {
      aPolynomial[i] = coefficient[i][2];
      bPolynomial[i] = coefficient[i][1];
      cPolynomial[i] = coefficient[i][0];
    }
    std::vector<long double> a = SeriesPolynomial(aPolynomial,x);
    std::vector<long double> b = SeriesPolynomial(bPolynomial,x);
    std::vector<long double> c = SeriesPolynomial(cPolynomial,x);
    std::vector<long double> root(order+1U,0.0L);
    const bool quadratic = LineSliceUtilities::PolynomialNorm(aPolynomial) != 0.0L;
    if(!quadratic) {
      root = SeriesQuotient(c,b);
      for(unsigned k = 0; k <= order; ++k) root[k] = -root[k];
    }
    else {
      if(a[0] == 0.0L) throw std::runtime_error("Triangle Taylor quadratic coefficient vanishes at the expansion point.");
      std::vector<long double> discriminant = SeriesProduct(b,b);
      std::vector<long double> term = SeriesProduct(a,c);
      for(unsigned k = 0; k <= order; ++k) discriminant[k] -= 4.0L*term[k];
      std::vector<long double> squareRoot = SeriesSquareRoot(discriminant);
      const long double lowerMiddle = StripLower(strip,0.5L*(left+right));
      const long double upperMiddle = StripUpper(strip,0.5L*(left+right));
      const long double candidateMinus = (-b[0]-squareRoot[0])/(2.0L*a[0]);
      const long double candidatePlus = (-b[0]+squareRoot[0])/(2.0L*a[0]);
      const long double sign = candidateMinus >= lowerMiddle && candidateMinus <= upperMiddle ? -1.0L : (candidatePlus >= lowerMiddle && candidatePlus <= upperMiddle ? 1.0L : 0.0L);
      if(sign == 0.0L) throw std::runtime_error("Triangle Taylor graph root is outside the strip.");
      std::vector<long double> numerator(b);
      std::vector<long double> denominator(a);
      for(unsigned k = 0; k <= order; ++k) {
        numerator[k] = -numerator[k]+sign*squareRoot[k];
        denominator[k] *= 2.0L;
      }
      root = SeriesQuotient(numerator,denominator);
    }
    std::vector<long double> lowerSeries = x;
    std::vector<long double> upperSeries = x;
    for(unsigned k = 0; k <= order; ++k) {
      lowerSeries[k] *= strip.lowerSlope;
      upperSeries[k] *= strip.upperSlope;
    }
    lowerSeries[0] += strip.lowerConstant;
    upperSeries[0] += strip.upperConstant;
    std::vector<long double> gy = SeriesProduct(a,root);
    for(unsigned k = 0; k <= order; ++k) gy[k] = 2.0L*gy[k]+b[k];
    const long double signedDerivative = gy[0];
    if(signedDerivative == 0.0L) throw std::runtime_error("Triangle Taylor graph has a singular derivative.");
    if(signedDerivative < 0.0L) for(unsigned k = 0; k <= order; ++k) gy[k] = -gy[k];
    std::vector<long double> deltaWeight;
    if(computeDelta) {
      deltaWeight = SeriesQuotient(jacobian,gy);
      if(surfaceGradient != nullptr) deltaWeight = SeriesProduct(deltaWeight,SurfaceGradientSeriesNorm(*surfaceGradient,x,root));
    }
    std::vector<std::vector<long double>> lowerPowers(degree+2U,SeriesConstant(1.0L,order));
    std::vector<std::vector<long double>> upperPowers(degree+2U,SeriesConstant(1.0L,order));
    std::vector<std::vector<long double>> rootPowers(degree+2U,SeriesConstant(1.0L,order));
    for(unsigned j = 1; j <= degree+1U; ++j) {
      lowerPowers[j] = SeriesProduct(lowerPowers[j-1U],lowerSeries);
      upperPowers[j] = SeriesProduct(upperPowers[j-1U],upperSeries);
      rootPowers[j] = SeriesProduct(rootPowers[j-1U],root);
    }
    std::vector<std::vector<long double>> xPowers(degree+1U,SeriesConstant(1.0L,order));
    for(unsigned i = 1; i <= degree; ++i) xPowers[i] = SeriesProduct(xPowers[i-1U],x);
    std::vector<long double> values(size,0.0L);
    std::vector<long double> errors(size,0.0L);
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned index = LineSliceUtilities::MomentIndex(degree,i,j);
        if(active == nullptr || (*active)[index] != 0) {
          std::vector<long double> fiber(order+1U,0.0L);
          const std::vector<long double>& first = signedDerivative > 0.0L ? upperPowers[j+1U] : rootPowers[j+1U];
          const std::vector<long double>& second = signedDerivative > 0.0L ? rootPowers[j+1U] : lowerPowers[j+1U];
          for(unsigned k = 0; k <= order; ++k) fiber[k] = (first[k]-second[k])/(j+1U);
          fiber = SeriesProduct(fiber,jacobian);
          bool converged = true;
          SeriesProductIntegral(xPowers[i],fiber,values[index],errors[index],converged);
          if(!converged) errors[index] = std::numeric_limits<long double>::infinity();
        }
        if(computeDelta && (active == nullptr || (*active)[count+index] != 0)) {
          std::vector<long double> fiber = SeriesProduct(rootPowers[j],deltaWeight);
          bool converged = true;
          SeriesProductIntegral(xPowers[i],fiber,values[count+index],errors[count+index],converged);
          if(!converged) errors[count+index] = std::numeric_limits<long double>::infinity();
        }
      }
    }
    std::vector<unsigned char> failed(size,0);
    bool anyFailed = false;
    for(unsigned k = 0; k < size; ++k) {
      if(active != nullptr && (*active)[k] == 0) continue;
      const long double tolerance = (k < count ? absoluteTolerance : absoluteDeltaTolerance)+relativeTolerance*std::abs(values[k]);
      if(errors[k] <= tolerance) {
        total[k] += values[k];
        totalError[k] += errors[k];
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
    const long double middle = 0.5L*(left+right);
    if(depth >= maximumDepth || middle == left || middle == right) throw std::runtime_error("Triangle Taylor graph did not reach tolerance.");
    GraphTaylorInterval(coefficient,degree,strip,left,middle,computeDelta,0.5L*absoluteTolerance,0.5L*absoluteDeltaTolerance,relativeTolerance,depth+1U,maximumDepth,total,totalError,accepted,surfaceGradient,&failed);
    GraphTaylorInterval(coefficient,degree,strip,middle,right,computeDelta,0.5L*absoluteTolerance,0.5L*absoluteDeltaTolerance,relativeTolerance,depth+1U,maximumDepth,total,totalError,accepted,surfaceGradient,&failed);
  }

  static void FullGaussKronrodSample(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const TriangleDomainStrip& strip, long double x, bool computeDelta, std::vector<long double>& values, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>* active) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    values.assign(2U*count,0.0L);
    const long double lower = StripLower(strip,x);
    const long double upper = StripUpper(strip,x);
    std::vector<long double> polynomial = FiberPolynomial(coefficient,x);
    std::vector<long double> roots = LineSliceUtilities::QuadraticRoots(polynomial[2],polynomial[1],polynomial[0]);
    std::vector<long double> cuts = {lower,upper};
    for(unsigned k = 0; k < roots.size(); ++k) if(roots[k] > lower && roots[k] < upper) cuts.push_back(roots[k]);
    LineSliceUtilities::SortDistinct(cuts);
    long double xPower = 1.0L;
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned index = LineSliceUtilities::MomentIndex(degree,i,j);
        if(active == nullptr || (*active)[index] != 0) {
          long double positive = 0.0L;
          for(unsigned segment = 1; segment < cuts.size(); ++segment) {
            const long double middle = 0.5L*(cuts[segment-1U]+cuts[segment]);
            if(Evaluate(coefficient,x,middle) > 0.0L) positive += LineSliceUtilities::MonomialIntegral(cuts[segment-1U],cuts[segment],j);
          }
          values[index] = xPower*positive;
        }
        if(computeDelta && (active == nullptr || (*active)[count+index] != 0)) {
          long double delta = 0.0L;
          for(unsigned r = 0; r < roots.size(); ++r) {
            if(roots[r] < lower || roots[r] > upper) continue;
            const long double derivative = std::abs(std::fma(2.0L*polynomial[2],roots[r],polynomial[1]));
            if(!(derivative > 0.0L)) throw std::runtime_error("Triangle full GK found a singular interface point.");
            long double weight = 1.0L/derivative;
            if(surfaceGradient != nullptr) weight *= SurfaceGradientNorm(*surfaceGradient,x,roots[r]);
            delta += Power(roots[r],j)*weight;
          }
          values[count+index] = xPower*delta;
        }
      }
      xPower *= x;
    }
  }

  static void FullGaussKronrodInterval(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const TriangleDomainStrip& strip, long double left, long double right, bool computeDelta, long double absoluteTolerance, long double absoluteDeltaTolerance, long double relativeTolerance, unsigned depth, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, std::vector<long double>& total, std::vector<long double>& totalError, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, const std::vector<unsigned char>* active = nullptr) {
    static const long double nodes[8] = {0.991455371120812639206854697526329L,0.949107912342758524526189684047851L,0.864864423359769072789712788640926L,0.741531185599394439863864773280788L,0.586087235467691130294144838258730L,0.405845151377397166906606412076961L,0.207784955007898467600689403773245L,0.0L};
    static const long double kronrod[8] = {0.022935322010529224963732008058970L,0.063092092629978553290700663189204L,0.104790010322250183839876322541518L,0.140653259715525918745189590510238L,0.169004726639267902826583426598550L,0.190350578064785409913256402421014L,0.204432940075298892414161999234649L,0.209482141084727828012999174891714L};
    static const long double gauss[4] = {0.129484966168869693270611432679082L,0.279705391489276667901467771423780L,0.381830050505118944950369775488975L,0.417959183673469387755102040816327L};
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    const unsigned size = 2U*count;
    const long double middle = 0.5L*(left+right);
    const long double halfWidth = 0.5L*(right-left);
    std::vector<long double> high(size,0.0L);
    std::vector<long double> low(size,0.0L);
    std::vector<long double> absolute(size,0.0L);
    std::vector<std::vector<long double>> samples(15);
    unsigned sampleNumber = 0;
    for(unsigned node = 0; node < 8; ++node) {
      const unsigned sides = node == 7 ? 1U : 2U;
      for(unsigned side = 0; side < sides; ++side) {
        if(evaluations >= maximumEvaluations) throw std::runtime_error("Triangle full GK evaluation budget exhausted.");
        ++evaluations;
        const long double point = middle+(side == 0 ? -1.0L : 1.0L)*halfWidth*nodes[node];
        FullGaussKronrodSample(coefficient,degree,strip,point,computeDelta,samples[sampleNumber],surfaceGradient,active);
        for(unsigned k = 0; k < size; ++k) {
          if(active != nullptr && (*active)[k] == 0) continue;
          high[k] += kronrod[node]*samples[sampleNumber][k];
          absolute[k] += kronrod[node]*std::abs(samples[sampleNumber][k]);
          if((node&1U) != 0) low[k] += gauss[node/2U]*samples[sampleNumber][k];
        }
        ++sampleNumber;
      }
    }
    std::vector<unsigned char> failed(size,0);
    bool anyFailed = false;
    for(unsigned k = 0; k < size; ++k) {
      if(active != nullptr && (*active)[k] == 0) continue;
      long double variation = 0.0L;
      for(unsigned n = 0; n < sampleNumber; ++n) variation += kronrod[n/2U]*std::abs(samples[n][k]-0.5L*high[k]);
      long double estimate = std::abs(high[k]-low[k]);
      if(variation > 0.0L && estimate > 0.0L) estimate = variation*std::min(1.0L,std::pow(200.0L*estimate/variation,1.5L));
      const long double error = halfWidth*std::max(estimate,128.0L*std::numeric_limits<long double>::epsilon()*absolute[k]);
      const long double value = halfWidth*high[k];
      const long double tolerance = (k < count ? absoluteTolerance : absoluteDeltaTolerance)*(right-left)+relativeTolerance*std::abs(value);
      if(error <= tolerance) {
        total[k] += value;
        totalError[k] += error;
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
    if(depth >= maximumDepth || middle == left || middle == right) throw std::runtime_error("Triangle full GK did not reach tolerance.");
    FullGaussKronrodInterval(coefficient,degree,strip,left,middle,computeDelta,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError,surfaceGradient,&failed);
    FullGaussKronrodInterval(coefficient,degree,strip,middle,right,computeDelta,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance,depth+1U,maximumDepth,maximumEvaluations,evaluations,accepted,total,totalError,surfaceGradient,&failed);
  }

  static void IntegrateDirection(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const std::vector<TriangleDomainVertex>& inputVertices, unsigned direction, bool computeDelta, long double absoluteTolerance, long double absoluteDeltaTolerance, long double relativeTolerance, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, bool requireGraph, bool useTaylor) {
    std::vector<std::vector<long double>> workingCoefficient = coefficient;
    std::vector<TriangleDomainVertex> vertices = inputVertices;
    std::vector<std::vector<std::vector<long double>>> workingGradient;
    const std::vector<std::vector<std::vector<long double>>>* gradient = surfaceGradient;
    if(direction == 2) {
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = i+1U; j < 3; ++j) std::swap(workingCoefficient[i][j],workingCoefficient[j][i]);
      for(unsigned k = 0; k < vertices.size(); ++k) std::swap(vertices[k].x,vertices[k].y);
      if(surfaceGradient != nullptr) {
        SurfaceGradientTranspose(*surfaceGradient,workingGradient);
        gradient = &workingGradient;
      }
    }
    std::vector<TriangleDomainStrip> strips;
    BuildStrips(vertices,strips);
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    std::vector<long double> raw(2U*count,0.0L);
    std::vector<long double> rawError(2U*count,0.0L);
    std::vector<long double> correction(2U*count,0.0L);
    for(unsigned stripIndex = 0; stripIndex < strips.size(); ++stripIndex) {
      std::vector<long double> events;
      CollectEvents(workingCoefficient,strips[stripIndex],events);
      for(unsigned event = 1; event < events.size(); ++event) {
        const long double left = events[event-1U];
        const long double right = events[event];
        if(!(right > left)) continue;
        const long double middle = 0.5L*(left+right);
        const long double lower = StripLower(strips[stripIndex],middle);
        const long double upper = StripUpper(strips[stripIndex],middle);
        std::vector<long double> polynomial = FiberPolynomial(workingCoefficient,middle);
        std::vector<long double> roots = LineSliceUtilities::QuadraticRoots(polynomial[2],polynomial[1],polynomial[0]);
        unsigned rootsInside = 0;
        for(unsigned r = 0; r < roots.size(); ++r) rootsInside += roots[r] >= lower && roots[r] <= upper ? 1U : 0U;
        if(rootsInside == 0) {
          if(Evaluate(workingCoefficient,middle,0.5L*(lower+upper)) > 0.0L) {
            for(unsigned i = 0; i <= degree; ++i) {
              for(unsigned j = 0; j <= degree-i; ++j) {
                const unsigned index = LineSliceUtilities::MomentIndex(degree,i,j);
                raw[index] += StripMoment(strips[stripIndex],left,right,i,j);
              }
            }
          }
          ++accepted;
          continue;
        }
        if(requireGraph && rootsInside != 1) throw std::runtime_error("Triangle child is not a graph in the selected direction.");
        std::vector<long double> local(2U*count,0.0L);
        std::vector<long double> localError(2U*count,0.0L);
        bool completed = false;
        if(requireGraph) {
          try {
            GraphGaussKronrodInterval(workingCoefficient,degree,strips[stripIndex],left,right,computeDelta,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance,0,maximumDepth,maximumEvaluations,evaluations,accepted,local,localError,gradient);
            completed = true;
          }
          catch(const std::runtime_error&) {
            if(useTaylor) {
              std::fill(local.begin(),local.end(),0.0L);
              std::fill(localError.begin(),localError.end(),0.0L);
              try {
                GraphTaylorInterval(workingCoefficient,degree,strips[stripIndex],left,right,computeDelta,absoluteTolerance*(right-left),absoluteDeltaTolerance*(right-left),relativeTolerance,0,maximumDepth,local,localError,accepted,gradient);
                completed = true;
              }
              catch(const std::runtime_error&) {}
            }
          }
        }
        if(!completed) {
          std::fill(local.begin(),local.end(),0.0L);
          std::fill(localError.begin(),localError.end(),0.0L);
          FullGaussKronrodInterval(workingCoefficient,degree,strips[stripIndex],left,right,computeDelta,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance,0,maximumDepth,maximumEvaluations,evaluations,accepted,local,localError,gradient);
        }
        for(unsigned k = 0; k < raw.size(); ++k) {
          const long double adjusted = local[k]-correction[k];
          const long double next = raw[k]+adjusted;
          correction[k] = (next-raw[k])-adjusted;
          raw[k] = next;
          rawError[k] += localError[k];
        }
      }
    }
    if(direction == 1) {
      values.swap(raw);
      errors.swap(rawError);
      return;
    }
    values.assign(2U*count,0.0L);
    errors.assign(2U*count,0.0L);
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned target = LineSliceUtilities::MomentIndex(degree,i,j);
        const unsigned source = LineSliceUtilities::MomentIndex(degree,j,i);
        values[target] = raw[source];
        values[count+target] = raw[count+source];
        errors[target] = rawError[source];
        errors[count+target] = rawError[count+source];
      }
    }
  }

  static void CheckErrors(const std::vector<long double>& values, const std::vector<long double>& errors, unsigned degree, long double absoluteTolerance, long double absoluteDeltaTolerance, long double relativeTolerance) {
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    for(unsigned k = 0; k < values.size(); ++k) {
      const long double tolerance = (k < count ? absoluteTolerance : absoluteDeltaTolerance)+relativeTolerance*std::abs(values[k]);
      if(!std::isfinite(values[k]) || !std::isfinite(errors[k]) || errors[k] < 0.0L || errors[k] > tolerance) throw std::runtime_error(k < count ? "Triangle volume moment error exceeds tolerance." : "Triangle interface moment error exceeds tolerance.");
    }
  }
};

}

#endif
