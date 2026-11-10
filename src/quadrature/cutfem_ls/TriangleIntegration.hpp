#ifndef FEMUS_TRIANGLE_INTEGRATION_HPP
#define FEMUS_TRIANGLE_INTEGRATION_HPP

#include "LineIntegration.hpp"

namespace FELevelSetMoments {

struct TriangleDomainVertex {
  long double x;
  long double y;
};



class TriangleIntegration {
public:
  // Integrates the unchanged parent level-set polynomial over an arbitrary triangular subregion.
  static void Integrate(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const std::vector<TriangleDomainVertex>& vertices, bool computeDelta, long double absoluteTolerance, long double absoluteDeltaTolerance, long double relativeTolerance, unsigned maximumDepth, unsigned maximumEvaluations, unsigned splitDepth, unsigned maximumSplitDepth, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned& selectedDirection, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr, const std::array<long double,5>* coordinates = nullptr, const std::vector<unsigned char>* requested = nullptr, unsigned certifiedDirection = 0) {
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
      if(minimum > padding) FillFullTriangle(vertices,degree,values,errors,coordinates,requested);
      ++accepted;
      selectedDirection = 0;
      CheckErrors(values,errors,degree,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance);
      return;
    }
    long double lower[2] = {0.0L,0.0L};
    long double upper[2] = {0.0L,0.0L};
    if(certifiedDirection == 0) { DerivativeBounds(coefficient,vertices,0,lower[0],upper[0]); DerivativeBounds(coefficient,vertices,1,lower[1],upper[1]); }
    else lower[certifiedDirection == 1 ? 1 : 0] = 1.0L;
    unsigned direction = lower[1] > lower[0] || (lower[1] == lower[0] && upper[1] > upper[0]) ? 1U : 2U;
    const long double selectedLower = direction == 1 ? lower[1] : lower[0];
    const unsigned weakDirection = selectedLower == 0.0L ? WeakDerivativeDirection(coefficient,vertices) : 0U;
    if(weakDirection != 0) direction = weakDirection;
    if(selectedLower == 0.0L && weakDirection == 0 && splitDepth < maximumSplitDepth) {
      ++LineIntegration::Statistics().triangleSplitCells;
      std::vector<std::vector<TriangleDomainVertex>> children;
      if(!SplitAtDerivativeChange(coefficient,vertices,children)) {
        std::vector<TriangleDomainVertex> first;
        std::vector<TriangleDomainVertex> second;
        SplitLongestEdge(vertices,first,second);
        children.push_back(first);
        children.push_back(second);
      }
      const unsigned acceptedBeforeChildren = accepted;
      std::runtime_error lastFailure("Triangle child integration did not reach tolerance.");
      for(unsigned budgetPass = 0; budgetPass < 2; ++budgetPass) {
        values.assign(2U*count,0.0L);
        errors.assign(2U*count,0.0L);
        std::vector<long double> correction(values.size(),0.0L);
        unsigned commonDirection = 0;
        accepted = acceptedBeforeChildren;
        bool completedChildren = true;
        try {
          for(unsigned child = 0; child < children.size(); ++child) {
            std::vector<long double> childValues;
            std::vector<long double> childErrors;
            unsigned childDirection = 0;
            const long double fraction = budgetPass == 0 ? 1.0L : 1.0L/static_cast<long double>(children.size());
            Integrate(coefficient,degree,children[child],computeDelta,fraction*absoluteTolerance,fraction*absoluteDeltaTolerance,relativeTolerance,maximumDepth,maximumEvaluations,splitDepth+1U,maximumSplitDepth,childValues,childErrors,evaluations,accepted,childDirection,surfaceGradient,coordinates,requested,certifiedDirection);
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
          CheckErrors(values,errors,degree,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance);
        }
        catch(const std::runtime_error& exception) {
          completedChildren = false;
          lastFailure = std::runtime_error(exception.what());
        }
        if(completedChildren) {
          selectedDirection = commonDirection;
          return;
        }
      }
      throw lastFailure;
    }
    selectedDirection = direction;
    if(selectedLower > 0.0L || weakDirection != 0) {
      try {
        IntegrateDirection(coefficient,degree,vertices,direction,computeDelta,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,surfaceGradient,true,true,weakDirection != 0,coordinates,requested);
      }
      catch(const std::runtime_error&) {
        IntegrateDirection(coefficient,degree,vertices,direction,computeDelta,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,surfaceGradient,false,false,false,coordinates,requested);
      }
    }
    else IntegrateDirection(coefficient,degree,vertices,direction,computeDelta,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,surfaceGradient,false,false,false,coordinates,requested);
    CheckErrors(values,errors,degree,absoluteTolerance,absoluteDeltaTolerance,relativeTolerance);
  }

private:
  static unsigned WeakDerivativeDirection(const std::vector<std::vector<long double>>& coefficient, const std::vector<TriangleDomainVertex>& vertices) {
    unsigned selected = 0;
    long double selectedMagnitude = 0.0L;
    for(unsigned axis = 0; axis < 2; ++axis) {
      std::vector<std::vector<long double>> derivative(3,std::vector<long double>(3,0.0L));
      long double norm = 0.0L;
      if(axis == 0) {
        for(unsigned i = 1; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) {
          derivative[i-1U][j] = static_cast<long double>(i)*coefficient[i][j];
          norm += std::abs(derivative[i-1U][j]);
        }
      }
      else {
        for(unsigned i = 0; i < 3; ++i) for(unsigned j = 1; j < 3; ++j) {
          derivative[i][j-1U] = static_cast<long double>(j)*coefficient[i][j];
          norm += std::abs(derivative[i][j-1U]);
        }
      }
      long double minimum = 0.0L;
      long double maximum = 0.0L;
      PolynomialBounds(derivative,vertices,minimum,maximum);
      const long double padding = 1024.0L*std::numeric_limits<long double>::epsilon()*std::max(1.0L,norm);
      const bool nonnegative = minimum >= -padding && maximum > padding;
      const bool nonpositive = maximum <= padding && minimum < -padding;
      if(nonnegative || nonpositive) {
        const long double magnitude = std::max(std::abs(minimum),std::abs(maximum));
        if(magnitude > selectedMagnitude) {
          selectedMagnitude = magnitude;
          selected = axis == 0 ? 2U : 1U;
        }
      }
    }
    return selected;
  }

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
    std::vector<std::vector<std::vector<long double>>> local;
    SurfaceGradientAffine({coefficient},vertices[0].x,vertices[1].x-vertices[0].x,vertices[2].x-vertices[0].x,vertices[0].y,vertices[1].y-vertices[0].y,vertices[2].y-vertices[0].y,local);
    unsigned degree = 0;
    for(unsigned i = 0; i < local[0].size(); ++i) for(unsigned j = 0; j < local[0][i].size(); ++j) if(local[0][i][j] != 0) degree = std::max(degree,i+j);
    minimum = std::numeric_limits<long double>::infinity();
    maximum = -minimum;
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
      long double value = 0;
      for(unsigned a = 0; a <= i && a < local[0].size(); ++a) for(unsigned b = 0; b <= j && b < local[0][a].size(); ++b) {
        long double multiplier = 1;
        for(unsigned n = 0; n < a; ++n) multiplier *= static_cast<long double>(i-n)/(degree-n);
        for(unsigned n = 0; n < b; ++n) multiplier *= static_cast<long double>(j-n)/(degree-a-n);
        value = std::fma(multiplier,local[0][a][b],value);
      }
      minimum = std::min(minimum,value);
      maximum = std::max(maximum,value);
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
      for(unsigned i = 1; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) value += static_cast<long double>(i)*coefficient[i][j]*Power(vertex.x,i-1U)*Power(vertex.y,j);
    }
    else {
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = 1; j < 3; ++j) value += static_cast<long double>(j)*coefficient[i][j]*Power(vertex.x,i)*Power(vertex.y,j-1U);
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
    if(axis == 0) point.x = coordinate;
    else point.y = coordinate;
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
      if(determinant != 0.0L) triangles.push_back(triangle);
    }
  }

  static bool DerivativeAffineFactor(const std::vector<std::vector<long double>>& coefficient, const std::vector<TriangleDomainVertex>& vertices, unsigned axis, std::array<long double,3>& factor) {
    std::vector<std::vector<long double>> derivative(3,std::vector<long double>(3,0));
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) if((axis == 0 ? i : j) != 0) derivative[i-(axis == 0 ? 1U : 0U)][j-(axis == 1 ? 1U : 0U)] = (axis == 0 ? i : j)*coefficient[i][j];
    unsigned first = 3;
    unsigned second = 3;
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) if(derivative[i][j] != 0) { first = std::min(first,i); second = std::min(second,j); }
    if(first == 3) return false;
    for(unsigned d = 0; d < 2; ++d) if((d == 0 ? first : second)%2U != 0) {
      long double minimum = Coordinate(vertices[0],d);
      long double maximum = minimum;
      for(unsigned n = 1; n < 3; ++n) { minimum = std::min(minimum,Coordinate(vertices[n],d)); maximum = std::max(maximum,Coordinate(vertices[n],d)); }
      if(minimum < 0 && maximum > 0) { if(d == 0) first = 0; else second = 0; }
    }
    factor = {0,0,0};
    for(unsigned i = first; i < 3; ++i) for(unsigned j = second; j < 3; ++j) if(derivative[i][j] != 0) {
      if(i+j-first-second > 1) return false;
      factor[i > first ? 1U : (j > second ? 2U : 0U)] = derivative[i][j];
    }
    return factor[1] != 0 || factor[2] != 0;
  }

  static void ClipByAffineLine(const std::vector<TriangleDomainVertex>& input, const std::array<long double,3>& line, bool lowerSide, std::vector<TriangleDomainVertex>& output) {
    output.clear();
    TriangleDomainVertex previous = input.back();
    long double previousValue = std::fma(line[1],previous.x,std::fma(line[2],previous.y,line[0]));
    bool previousInside = lowerSide ? previousValue <= 0 : previousValue >= 0;
    for(unsigned n = 0; n < input.size(); ++n) {
      const TriangleDomainVertex current = input[n];
      const long double value = std::fma(line[1],current.x,std::fma(line[2],current.y,line[0]));
      const bool inside = lowerSide ? value <= 0 : value >= 0;
      if(inside != previousInside) {
        const long double parameter = previousValue/(previousValue-value);
        output.push_back({std::fma(parameter,current.x-previous.x,previous.x),std::fma(parameter,current.y-previous.y,previous.y)});
      }
      if(inside) output.push_back(current);
      previous = current; previousValue = value; previousInside = inside;
    }
  }

  static bool SplitAtAffineDerivative(const std::vector<std::vector<long double>>& coefficient, const std::vector<TriangleDomainVertex>& vertices, std::vector<std::vector<TriangleDomainVertex>>& children) {
    for(unsigned axis = 0; axis < 2; ++axis) {
      std::array<long double,3> line;
      if(!DerivativeAffineFactor(coefficient,vertices,axis,line)) continue;
      long double minimum = std::numeric_limits<long double>::infinity();
      long double maximum = -minimum;
      for(unsigned n = 0; n < vertices.size(); ++n) { const long double value = std::fma(line[1],vertices[n].x,std::fma(line[2],vertices[n].y,line[0])); minimum = std::min(minimum,value); maximum = std::max(maximum,value); }
      const long double padding = 1024*std::numeric_limits<long double>::epsilon()*std::max(1.0L,std::max(std::abs(minimum),std::abs(maximum)));
      if(!(minimum < -padding && maximum > padding)) continue;
      std::vector<TriangleDomainVertex> lower;
      std::vector<TriangleDomainVertex> upper;
      ClipByAffineLine(vertices,line,true,lower);
      ClipByAffineLine(vertices,line,false,upper);
      children.clear();
      TriangulatePolygon(lower,children);
      TriangulatePolygon(upper,children);
      if(children.size() >= 2) return true;
    }
    return false;
  }

  static bool SplitAtDerivativeChange(const std::vector<std::vector<long double>>& coefficient, const std::vector<TriangleDomainVertex>& vertices, std::vector<std::vector<TriangleDomainVertex>>& children) {
    if(SplitAtAffineDerivative(coefficient,vertices,children)) return true;
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

  static long double TriangleMoment(const std::vector<TriangleDomainVertex>& vertices, unsigned firstDegree, unsigned secondDegree, long double* estimatedError = nullptr) {
    const unsigned degree = firstDegree+secondDegree;
    std::vector<long double> coefficient(1,1.0L);
    std::vector<long double> magnitude(1,1.0L);
    for(unsigned n = 1; n <= degree; ++n) {
      std::vector<long double> next((n+1U)*(n+2U)/2U,0.0L);
      std::vector<long double> bound(next.size(),0.0L);
      for(unsigned a = 0; a <= n; ++a) {
        for(unsigned b = 0; b <= n-a; ++b) {
          const unsigned c = n-a-b;
          const unsigned target = a*(2U*n+3U-a)/2U+b;
          const unsigned power[3] = {a,b,c};
          for(unsigned v = 0; v < 3; ++v) {
            if(power[v] == 0) continue;
            const unsigned previousA = a-(v == 0 ? 1U : 0U);
            const unsigned previousB = b-(v == 1 ? 1U : 0U);
            const unsigned source = previousA*(2U*n+1U-previousA)/2U+previousB;
            const long double coordinate = n <= firstDegree ? vertices[v].x : vertices[v].y;
            const long double multiplier = coordinate*static_cast<long double>(power[v])/n;
            next[target] += multiplier*coefficient[source];
            bound[target] += std::abs(multiplier)*magnitude[source];
          }
        }
      }
      coefficient.swap(next);
      magnitude.swap(bound);
    }
    const long double determinant = std::abs(std::fma(vertices[1].x-vertices[0].x,vertices[2].y-vertices[0].y,-(vertices[2].x-vertices[0].x)*(vertices[1].y-vertices[0].y)));
    const long double factor = determinant/((degree+1.0L)*(degree+2.0L));
    long double value = 0.0L;
    long double correction = 0.0L;
    long double absolute = 0.0L;
    for(unsigned k = 0; k < coefficient.size(); ++k) {
      const long double term = coefficient[k]-correction;
      const long double next = value+term;
      correction = (next-value)-term;
      value = next;
      absolute += magnitude[k];
    }
    if(estimatedError != nullptr) *estimatedError = (32.0L+16.0L*degree+4.0L*coefficient.size())*std::numeric_limits<long double>::epsilon()*factor*absolute;
    return factor*value;
  }

  static void FillFullTriangle(const std::vector<TriangleDomainVertex>& vertices, unsigned degree, std::vector<long double>& values, std::vector<long double>& errors, const std::array<long double,5>* coordinates = nullptr, const std::vector<unsigned char>* requested = nullptr) {
    std::vector<TriangleDomainVertex> parent = vertices;
    if(coordinates != nullptr) for(unsigned n = 0; n < parent.size(); ++n) { parent[n].x = std::fma((*coordinates)[1],parent[n].x,(*coordinates)[0]); parent[n].y = std::fma((*coordinates)[3],parent[n].y,(*coordinates)[2]); }
    const long double measure = coordinates == nullptr ? 1.0L : (*coordinates)[4]/std::abs((*coordinates)[1]*(*coordinates)[3]);
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        const unsigned index = LineSliceUtilities::MomentIndex(degree,i,j);
        if(requested != nullptr && (*requested)[index] == 0) continue;
        values[index] = measure*TriangleMoment(parent,i,j,&errors[index]);
        errors[index] *= measure;
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

  static void AddStrip(const TriangleDomainVertex& first, const TriangleDomainVertex& second, const TriangleDomainVertex& longFirst, const TriangleDomainVertex& longSecond, std::vector<LineIntegrationStrip>& strips) {
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
    LineIntegrationStrip strip;
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

  static void BuildStrips(const std::vector<TriangleDomainVertex>& input, std::vector<LineIntegrationStrip>& strips) {
    std::vector<TriangleDomainVertex> vertices(input);
    std::sort(vertices.begin(),vertices.end(),VertexLess);
    strips.clear();
    if(vertices[0].x == vertices[2].x) return;
    if(vertices[0].x != vertices[1].x) AddStrip(vertices[0],vertices[1],vertices[0],vertices[2],strips);
    if(vertices[1].x != vertices[2].x) AddStrip(vertices[1],vertices[2],vertices[0],vertices[2],strips);
  }

  static void IntegrateDirection(const std::vector<std::vector<long double>>& coefficient, unsigned degree, const std::vector<TriangleDomainVertex>& inputVertices, unsigned direction, bool computeDelta, long double absoluteTolerance, long double absoluteDeltaTolerance, long double relativeTolerance, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, bool requireGraph, bool useTaylor, bool mappedGraph, const std::array<long double,5>* coordinates, const std::vector<unsigned char>* requested) {
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
    const std::array<long double,5> parent = LineIntegration::MomentCoordinates(coordinates,direction);
    const std::vector<unsigned char> active = LineIntegration::MomentRequest(degree,requested,direction,computeDelta);
    std::vector<LineIntegrationStrip> strips;
    BuildStrips(vertices,strips);
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    std::vector<long double> raw(2U*count,0.0L);
    std::vector<long double> rawError(2U*count,0.0L);
    std::vector<long double> correction(2U*count,0.0L);
    for(unsigned stripIndex = 0; stripIndex < strips.size(); ++stripIndex) LineIntegration::IntegrateStrip(workingCoefficient,degree,strips[stripIndex],requireGraph,mappedGraph,useTaylor,computeDelta,absoluteTolerance/strips.size(),absoluteDeltaTolerance/strips.size(),relativeTolerance,maximumDepth,maximumEvaluations,evaluations,accepted,raw,rawError,gradient,&parent,&active);
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
