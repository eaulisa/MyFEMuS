#ifndef FEMUS_PREPARED_ENDPOINT_PROJECTION_HPP
#define FEMUS_PREPARED_ENDPOINT_PROJECTION_HPP

#include "Hex27LevelSetMoments.hpp"
#include "CutFemPreparedContext.hpp"
#include "GaussKronrodRules.hpp"

namespace FELevelSetMoments {

struct EndpointProjectionEvent {
  long double point = 0;
  unsigned exponent = 1;
  bool logarithmic = false;
};

struct EndpointCriticalPoint {
  long double height;
  std::vector<long double> point;
  bool saddle;
};

struct EndpointEventPlan {
  std::vector<EndpointProjectionEvent> events;
  std::string failure;
};

struct EndpointAxisProjection {
  bool analyzed = false;
  bool hasSaddle = false;
  const std::vector<std::vector<long double>>* normalMetric = nullptr;
  unsigned long long metricRevision = 0;
  std::array<bool,3> dominant = {};
  std::vector<EndpointCriticalPoint> critical;
  std::array<std::vector<long double>,3> projectedRoots;
  std::array<bool,3> repeated = {};
  std::array<std::array<EndpointEventPlan,3>,2> plans;
};

class EndpointProjection {
  GeomElType geometry = LINE;
  std::vector<long double> coefficients;
  std::array<EndpointAxisProjection,3> axes;
  unsigned faceRule = 0;
  unsigned faceOuter = 3;

  static void ScoreFaceFiber(const CutFemLevelSetData& field, const CutFemLevelSetData& face, unsigned root, unsigned other, long double coordinate, std::vector<long double>& point, std::array<long double,3>& alignment) {
    point[other] = coordinate;
    CutFemLevelSetData::LineAnalysis line;
    face.AnalyzeLocalFiber(root < other ? 0U : 1U,coordinate,line);
    for(unsigned i = 0; i < line.rootCount; ++i) {
      point[root] = line.roots[i].coordinate;
      std::array<long double,3> gradient;
      field.GetGradient(point,gradient);
      const long double scale = std::max(std::abs(gradient[0]),std::max(std::abs(gradient[1]),std::abs(gradient[2])));
      if(scale == 0) continue;
      for(long double& g : gradient) g /= scale;
      const long double norm2 = gradient[0]*gradient[0]+gradient[1]*gradient[1]+gradient[2]*gradient[2];
      for(unsigned d = 0; d < 3; ++d) alignment[d] = std::max(alignment[d],gradient[d]*gradient[d]/norm2);
    }
  }

  // Minimize the largest squared reference-normal projection on the HEX face traces.
  static unsigned FindHexOuter(const CutFemLevelSetData& field, unsigned order) {
    for(unsigned d = 0; d < 3; ++d) if(!field.DependsOn(d)) return d;
    std::array<long double,3> lower, upper;
    for(unsigned d = 0; d < 3; ++d) {
      const std::array<long double,2>& bounds = field.GetGraphCertificate(d).bounds;
      lower[d] = std::max(0.0L,std::max(bounds[0],-bounds[1]));
      upper[d] = std::max(std::abs(bounds[0]),std::abs(bounds[1]));
    }
    for(unsigned d = 0; d < 3; ++d) if(upper[d] <= std::min(lower[(d+1)%3],lower[(d+2)%3])) return d;
    std::array<long double,3> alignment = {-1,-1,-1};
    const GaussKronrodRule& rule = GaussKronrodRules::Get(order);
    CutFemLevelSetData face;
    std::vector<long double> point(3);
    for(unsigned fixed = 0; fixed < 3; ++fixed) for(unsigned side = 0; side < 2; ++side) {
      point[fixed] = side;
      field.Restrict(fixed,side,face,false);
      const std::array<long double,2> bounds = face.GetBounds();
      if(bounds[0] > 0 || bounds[1] < 0) continue;
      for(unsigned root = 0; root < 3; ++root) if(root != fixed) {
        const unsigned other = 3-fixed-root;
        ScoreFaceFiber(field,face,root,other,0.5L,point,alignment);
        for(unsigned k = 1; k < rule.size; ++k) {
          ScoreFaceFiber(field,face,root,other,(1-rule.nodes[k])/2,point,alignment);
          ScoreFaceFiber(field,face,root,other,(1+rule.nodes[k])/2,point,alignment);
        }
        if(fixed < other) for(unsigned edge = 0; edge < 2; ++edge) ScoreFaceFiber(field,face,root,other,edge,point,alignment);
      }
    }
    return alignment[0] < 0 ? 3U : static_cast<unsigned>(std::min_element(alignment.begin(),alignment.end())-alignment.begin());
  }

public:
  unsigned SharedHexOuter(const CutFemLevelSetData& field, unsigned order) {
    if(faceRule != order) { faceOuter = FindHexOuter(field,order); faceRule = order; }
    return faceOuter;
  }

  void Bind(const CutFemLevelSetData& data) {
    const std::vector<long double>& normalized = data.GetNormalizedPowerCoefficients();
    if(geometry == data.GetGeometry() && coefficients == normalized) return;
    faceRule = 0;
    geometry = data.GetGeometry();
    coefficients = normalized;
    for(EndpointAxisProjection& axis : axes) {
      axis.analyzed = false;
      for(std::array<EndpointEventPlan,3>& plans : axis.plans) for(EndpointEventPlan& plan : plans) { plan.events.clear(); plan.failure.clear(); }
    }
  }

  static void AddRoots(const std::vector<long double>& polynomial, long double lower, std::vector<long double>& events) { if(PolynomialNorm(polynomial) != 0) PolynomialRoots(polynomial,lower,1,events); }

  static std::vector<std::vector<std::vector<long double>>> PreparedField(const std::vector<long double>& coefficient, unsigned outer, unsigned inner, const CutFemLevelSetData* data) {
    if(data == nullptr) return Field(coefficient,outer,inner);
    const CutFemLevelSetData::FiberProjection& projection = data->GetFiberProjection(inner);
    const unsigned first = inner == 0 ? 1U : 0U;
    std::vector<std::vector<std::vector<long double>>> result(3,std::vector<std::vector<long double>>(3,std::vector<long double>(3)));
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) result[i][j][k] = projection.coefficients[i][first == outer ? 3*k+j : 3*j+k];
    return result;
  }

  static std::vector<std::vector<long double>> PreparedDiscriminant(const std::vector<std::vector<std::vector<long double>>>& field, unsigned outer, unsigned inner, const CutFemLevelSetData* data) {
    if(data == nullptr) return Discriminant(field);
    const CutFemLevelSetData::FiberProjection& projection = data->GetFiberProjection(inner);
    const unsigned first = inner == 0 ? 1U : 0U;
    std::vector<std::vector<long double>> result(5,std::vector<long double>(5));
    for(unsigned j = 0; j < 5; ++j) for(unsigned k = 0; k < 5; ++k) result[j][k] = projection.discriminant[first == outer ? 5*k+j : 5*j+k];
    Hex27TrimBivariate(result);
    return result;
  }

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

  static std::vector<long double> BoundaryValue(const std::vector<std::vector<long double>>& polynomial, bool top, bool triangle, long double lower) {
    std::vector<long double> result(1,0);
    std::vector<long double> power(1,1);
    for(unsigned i = 0; i < polynomial.size(); ++i) {
      result = PolynomialSum(result,PolynomialProduct(polynomial[i],power));
      power = PolynomialProduct(power,top ? (triangle ? std::vector<long double>{1,-1} : std::vector<long double>{1}) : std::vector<long double>{lower});
    }
    return result;
  }

  static std::vector<std::vector<long double>> Face(const std::vector<std::vector<std::vector<long double>>>& field, bool top, bool triangle, bool middleBoundary = true, bool outerBoundary = true, long double lower = 0) {
    std::vector<std::vector<long double>> result(5,std::vector<long double>(1,0));
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) {
      for(unsigned a = 0; a <= (top && triangle && middleBoundary ? i : 0U); ++a) {
        std::vector<long double> power(1,1);
        for(unsigned k = 0; k < (top && triangle && outerBoundary ? i-a : 0U); ++k) power = PolynomialProduct(power,{1,-1});
        long double factor = ReferenceBinomial(i,a)*(a%2U == 0 ? 1.0L : -1.0L);
        if(!triangle) factor = 1;
        if(!top) factor = ReferencePower(lower,i);
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

  static std::vector<std::vector<long double>> FaceNormalDerivative(const std::vector<std::vector<std::vector<long double>>>& field, bool top, bool triangle, bool middleBoundary = true, bool outerBoundary = true, long double lower = 0) {
    std::vector<std::vector<std::vector<long double>>> derivative(3,std::vector<std::vector<long double>>(3,std::vector<long double>(3,0)));
    for(unsigned i = 0; i < 2; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) derivative[i][j][k] = (i+1U)*field[i+1U][j][k];
    return Face(derivative,top,triangle,middleBoundary,outerBoundary,lower);
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

  static bool Relevant(const std::vector<std::vector<long double>>& polynomial, long double height, bool triangle, long double lower) {
    std::vector<long double> fiber(polynomial.size());
    long double norm = 0;
    for(unsigned i = 0; i < polynomial.size(); ++i) { fiber[i] = PolynomialValue(polynomial[i],height); norm += PolynomialBound(polynomial[i],height); }
    const long double upper = triangle ? 1-height : 1;
    const long double padding = 8192*std::numeric_limits<long double>::epsilon()*norm;
    if(std::abs(PolynomialValue(fiber,lower)) <= padding || std::abs(PolynomialValue(fiber,upper)) <= padding) return true;
    std::vector<long double> critical;
    PolynomialRoots(PolynomialDerivative(fiber),lower,upper,critical);
    for(unsigned i = 0; i < critical.size(); ++i) if(std::abs(PolynomialValue(fiber,critical[i])) <= padding) return true;
    return false;
  }

  static void AddConstraint(const std::vector<std::vector<long double>>& constraint, bool triangle, std::vector<long double>& events, bool regularGraphBoundary = false, long double middleLower = 0, long double outerLower = 0, const std::vector<long double>* projectedRoots = nullptr, bool repeated = false) {
    std::vector<std::vector<long double>> polynomial = constraint;
    Hex27TrimBivariate(polynomial);
    std::vector<long double> candidates;
    for(unsigned i = 0; i < polynomial.size(); ++i) AddRoots(polynomial[i],outerLower,candidates);
    AddRoots(BoundaryValue(polynomial,false,triangle,middleLower),outerLower,candidates);
    AddRoots(BoundaryValue(polynomial,true,triangle,middleLower),outerLower,candidates);
    if(polynomial.size() > 1) {
      if(projectedRoots != nullptr) candidates.insert(candidates.end(),projectedRoots->begin(),projectedRoots->end());
      else {
        const std::vector<long double> resultant = Hex27Resultant(polynomial,Derivative(polynomial));
        repeated = PolynomialNorm(resultant) == 0;
        AddRoots(resultant,outerLower,candidates);
      }
      if(repeated) {
        const bool repeatedQuadraticTrace = regularGraphBoundary && polynomial.size() == 3 && PolynomialNorm(PolynomialSum(PolynomialProduct(polynomial[1],polynomial[1]),PolynomialProduct(polynomial[0],polynomial[2]),-4)) == 0;
        if(!repeatedQuadraticTrace) throw std::runtime_error("Unsupported repeated projection factor; endpoint topology is not certified.");
      }
    }
    for(unsigned i = 0; i < candidates.size(); ++i) if(Relevant(polynomial,candidates[i],triangle,middleLower)) events.push_back(candidates[i]);
  }

  static unsigned CommonExponent(unsigned first, unsigned second) {
    unsigned a = first;
    unsigned b = second;
    while(b != 0) { const unsigned remainder = a%b; a = b; b = remainder; }
    return first/a*second;
  }

  static long double FoldVertex(const std::vector<std::vector<std::vector<long double>>>& field, long double middle, long double height, bool triangle, bool middleBoundary = true, bool outerBoundary = true, long double lower = 0, bool verifyZero = false) {
    long double a = 0;
    long double b = 0;
    for(unsigned i = 0; i < field[0].size(); ++i) { a += PolynomialValue(field[2][i],height)*ReferencePower(middle,i); b += PolynomialValue(field[1][i],height)*ReferencePower(middle,i); }
    if(a == 0) return std::numeric_limits<long double>::infinity();
    const long double root = -b/(2*a);
    const long double padding = 65536*std::numeric_limits<long double>::epsilon();
    if(!(root >= lower-padding && root <= (triangle ? 1-(outerBoundary ? height : 0)-(middleBoundary ? middle : 0) : 1)+padding)) return std::numeric_limits<long double>::infinity();
    if(verifyZero) {
      long double c = 0;
      long double bound = 0;
      for(unsigned j = 0; j < field[0].size(); ++j) {
        c += PolynomialValue(field[0][j],height)*ReferencePower(middle,j);
        for(unsigned i = 0; i < 3; ++i) bound += PolynomialBound(field[i][j],height)*ReferencePower(middle,j)*ReferencePower(std::abs(root),i);
      }
      if(std::abs(std::fma(root,std::fma(a,root,b),c)) > padding*bound) return std::numeric_limits<long double>::infinity();
    }
    return root;
  }

  static bool Interior(GeomElType geometry, const std::vector<long double>& point, long double padding) {
    for(long double x : point) if(!(x > padding && x < 1-padding)) return false;
    if(geometry == TET) return point[0]+point[1]+point[2] < 1-padding;
    return geometry != WEDGE || point[0]+point[1] < 1-padding;
  }

  // Refine G = G_u = G_v = 0 and classify regular interior critical points of the slicing height.
  static bool RefineCriticalPoint(const CutFemLevelSetData& data, unsigned outer, unsigned u, unsigned v, std::vector<long double>& point, bool& saddle) {
    const long double tolerance = 512*std::numeric_limits<long double>::epsilon();
    const long double scale = data.GetNormalizationScale();
    for(unsigned iteration = 0; iteration < 16; ++iteration) {
      if(!Interior(data.GetGeometry(),point,0)) return false;
      std::array<long double,3> gradient;
      data.GetGradient(point,gradient);
      const long double value = data.Evaluate(point)/scale, gu = gradient[u]/scale, gv = gradient[v]/scale, gt = gradient[outer]/scale;
      const long double uu = data.EvaluateSecondDerivative(u,u,point)/scale, uv = data.EvaluateSecondDerivative(u,v,point)/scale, vv = data.EvaluateSecondDerivative(v,v,point)/scale;
      const long double determinant = std::fma(uu,vv,-uv*uv);
      if(std::abs(determinant) <= tolerance*(std::abs(uu*vv)+uv*uv)) return false;
      const long double ut = data.EvaluateSecondDerivative(u,outer,point)/scale, vt = data.EvaluateSecondDerivative(v,outer,point)/scale;
      const long double du = (vv*gu-uv*gv)/determinant, dv = (uu*gv-uv*gu)/determinant;
      const long double su = (vv*ut-uv*vt)/determinant, sv = (uu*vt-uv*ut)/determinant;
      const long double slope = gt-gu*su-gv*sv;
      if(std::abs(slope) <= tolerance*(std::abs(gt)+std::abs(gu*su)+std::abs(gv*sv))) return false;
      const long double dt = (value-gu*du-gv*dv)/slope;
      const long double step = std::max({std::abs(du-su*dt),std::abs(dv-sv*dt),std::abs(dt)});
      const bool converged = step <= tolerance && std::max({std::abs(value),std::abs(gu),std::abs(gv)}) <= tolerance;
      point[u] -= du-su*dt;
      point[v] -= dv-sv*dt;
      point[outer] -= dt;
      if(converged) {
        if(!Interior(data.GetGeometry(),point,tolerance)) return false;
        saddle = determinant < 0;
        return std::abs(gt) > tolerance;
      }
    }
    return false;
  }

  static void RecordCriticalPoint(EndpointAxisProjection& axis, unsigned outer, const std::vector<long double>& point, bool saddle, long double& height) {
    const long double padding = 512*std::numeric_limits<long double>::epsilon();
    height = point[outer];
    for(const EndpointCriticalPoint& previous : axis.critical) if(std::abs(previous.height-height) <= padding) {
      height = previous.height;
      if(std::abs(previous.point[0]-point[0]) <= padding && std::abs(previous.point[1]-point[1]) <= padding && std::abs(previous.point[2]-point[2]) <= padding) return;
    }
    axis.critical.push_back({height,point,saddle});
    axis.hasSaddle = axis.hasSaddle || saddle;
  }

  // Near a squared discriminant, seed Newton from the transverse bilinear part and refine on the full field.
  static void AddBilinearSeeds(const CutFemLevelSetData& data, unsigned outer, unsigned u, unsigned v, const std::vector<std::vector<std::vector<long double>>>& field, EndpointAxisProjection& axis) {
    long double curvature = 0;
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) if(i == 2 || j == 2) curvature += PolynomialNorm(field[i][j]);
    if(curvature == 0 || curvature > std::sqrt(std::numeric_limits<long double>::epsilon())*PolynomialNorm(field[1][1])) return;
    std::vector<long double> heights;
    AddRoots(PolynomialSum(PolynomialProduct(field[1][1],field[0][0]),PolynomialProduct(field[1][0],field[0][1]),-1),0,heights);
    for(long double height : heights) {
      const long double mixed = PolynomialValue(field[1][1],height);
      if(mixed == 0) continue;
      std::vector<long double> point(3);
      point[outer] = height;
      point[u] = -PolynomialValue(field[0][1],height)/mixed;
      point[v] = -PolynomialValue(field[1][0],height)/mixed;
      bool saddle = false;
      if(RefineCriticalPoint(data,outer,u,v,point,saddle)) RecordCriticalPoint(axis,outer,point,saddle,height);
    }
  }

  EndpointAxisProjection& CriticalPoints(const CutFemLevelSetData& data, unsigned outer, unsigned long long& searches) {
    EndpointAxisProjection& axis = axes[outer];
    if(axis.analyzed) return axis;
    ++searches;
    axis.critical.clear();
    axis.hasSaddle = false;
    axis.normalMetric = nullptr;
    const bool tet = geometry == TET;
    const bool wedge = geometry == WEDGE;
    for(unsigned u = 0; u < 3; ++u) if(u != outer) {
      const unsigned v = 3-outer-u;
      const std::vector<std::vector<std::vector<long double>>> field = PreparedField(data.GetNormalizedPowerCoefficients(),outer,u,&data);
      if(u < v) AddBilinearSeeds(data,outer,u,v,field,axis);
      bool quadratic = false;
      for(const std::vector<long double>& row : field[2]) if(PolynomialNorm(row) != 0) quadratic = true;
      const std::vector<std::vector<long double>> trace = quadratic ? PreparedDiscriminant(field,outer,u,&data) : field[1];
      const std::vector<std::vector<long double>> stationary = quadratic ? Derivative(trace) : trace;
      const std::vector<long double> resultant = quadratic ? Hex27Resultant(trace,stationary) : Hex27Resultant(field[0],field[1]);
      axis.repeated[u] = PolynomialNorm(resultant) == 0;
      std::vector<long double>& candidates = axis.projectedRoots[u];
      candidates.clear();
      AddRoots(resultant,0,candidates);
      const bool limitedMiddle = tet || (wedge && outer < 2 && v < 2);
      for(long double& candidate : candidates) {
        const long double height = candidate;
        bool matched = false;
        std::vector<long double> fiber(stationary.size()), roots;
        for(unsigned j = 0; j < fiber.size(); ++j) fiber[j] = PolynomialValue(stationary[j],height);
        PolynomialRoots(fiber,0,limitedMiddle ? 1-height : 1,roots);
        const std::vector<std::vector<long double>>& constraint = quadratic ? trace : field[0];
        std::vector<long double> section(constraint.size());
        long double bound = 0;
        for(unsigned j = 0; j < section.size(); ++j) { section[j] = PolynomialValue(constraint[j],height); bound += PolynomialBound(constraint[j],height); }
        for(long double middle : roots) {
          std::vector<long double> point(3);
          point[v] = middle;
          point[outer] = height;
          if(quadratic) point[u] = FoldVertex(field,middle,height,tet || (wedge && u < 2),tet || v < 2,tet || outer < 2);
          else {
            const long double duv = PolynomialValue(field[1][1],height)+2*middle*PolynomialValue(field[1][2],height);
            if(duv == 0) continue;
            point[u] = -(PolynomialValue(field[0][1],height)+2*middle*PolynomialValue(field[0][2],height))/duv;
          }
          bool saddle = false;
          if(!RefineCriticalPoint(data,outer,u,v,point,saddle)) continue;
          long double refined;
          RecordCriticalPoint(axis,outer,point,saddle,refined);
          const bool projectedRoot = std::abs(PolynomialValue(section,middle)) <= 8192*std::numeric_limits<long double>::epsilon()*bound;
          if(projectedRoot && !matched) { candidate = refined; matched = true; }
        }
      }
    }
    axis.analyzed = true;
    return axis;
  }

  static unsigned ContactExponent(const std::vector<std::vector<long double>>& constraint, long double height, bool triangle, bool regularGraphBoundary = false, const std::vector<std::vector<long double>>* normalDerivative = nullptr, const std::vector<std::vector<std::vector<long double>>>* foldField = nullptr, bool wedge = false, bool middleBoundary = true, bool outerBoundary = true, long double innerLower = 0, long double middleLower = 0) {
    std::vector<long double> polynomial(constraint.size());
    for(unsigned i = 0; i < constraint.size(); ++i) polynomial[i] = PolynomialValue(constraint[i],height);
    Trim(polynomial);
    if(PolynomialNorm(polynomial) == 0) {
      if(foldField == nullptr) return 1;
      const long double upper = triangle ? 1-height : 1;
      for(const long double middle : {middleLower,0.5L*(middleLower+upper),upper}) if(std::isfinite(FoldVertex(*foldField,middle,height,triangle || wedge,middleBoundary,outerBoundary,innerLower,true))) {
        long double transverse = 0;
        long double norm = 0;
        for(unsigned i = 0; i < constraint.size(); ++i) { transverse += PolynomialValue(PolynomialDerivative(constraint[i]),height)*ReferencePower(middle,i); norm += PolynomialNorm(constraint[i]); }
        if(std::abs(transverse) <= 65536*std::numeric_limits<long double>::epsilon()*norm) throw std::runtime_error("Unsupported non-transverse zero projection fiber.");
        return 2;
      }
      return 1;
    }
    std::vector<long double> roots;
    PolynomialRoots(polynomial,middleLower,triangle ? 1-height : 1,roots);
    std::vector<long double> stationary;
    PolynomialRoots(PolynomialDerivative(polynomial),middleLower,triangle ? 1-height : 1,stationary);
    const long double contactPadding = 65536*std::numeric_limits<long double>::epsilon()*PolynomialNorm(polynomial);
    for(unsigned n = 0; n < stationary.size(); ++n) if(std::abs(PolynomialValue(polynomial,stationary[n])) <= contactPadding) roots.push_back(stationary[n]);
    for(const long double boundary : {middleLower,triangle ? 1-height : 1.0L}) if(std::abs(PolynomialValue(polynomial,boundary)) <= contactPadding) roots.push_back(boundary);
    SortDistinct(roots);
    unsigned exponent = 1;
    for(unsigned n = 0; n < roots.size(); ++n) {
      if(foldField != nullptr && !std::isfinite(FoldVertex(*foldField,roots[n],height,triangle || wedge,middleBoundary,outerBoundary,innerLower,true))) continue;
      std::vector<long double> derivative = PolynomialDerivative(polynomial);
      unsigned multiplicity = 1;
      while(derivative.size() > 1 && std::abs(PolynomialValue(derivative,roots[n])) <= 65536*std::numeric_limits<long double>::epsilon()*PolynomialNorm(derivative)) { ++multiplicity; derivative = PolynomialDerivative(derivative); }
      if(multiplicity < 2) {
        const long double upper = triangle ? 1-height : 1;
        const long double padding = 65536*std::numeric_limits<long double>::epsilon();
        if(foldField != nullptr && (roots[n] <= middleLower+padding || roots[n] >= upper-padding)) exponent = CommonExponent(exponent,2);
        continue;
      }
      long double transverse = 0;
      for(unsigned i = 0; i < constraint.size(); ++i) transverse += PolynomialValue(PolynomialDerivative(constraint[i]),height)*ReferencePower(roots[n],i);
      long double norm = 0;
      for(unsigned i = 0; i < constraint.size(); ++i) norm += PolynomialNorm(constraint[i]);
      if(std::abs(transverse) <= 65536*std::numeric_limits<long double>::epsilon()*norm && normalDerivative != nullptr && RegularFaceTangency(constraint,*normalDerivative,roots[n],height)) continue;
      if(!regularGraphBoundary && std::abs(transverse) <= 65536*std::numeric_limits<long double>::epsilon()*norm) throw std::runtime_error("Unsupported non-transverse projected endpoint; no regular-fold convergence claim is available.");
      if(multiplicity > 4) throw std::runtime_error("Unsupported projected contact multiplicity.");
      exponent = CommonExponent(exponent,CommonExponent(2,multiplicity));
    }
    return exponent;
  }

  static std::vector<EndpointProjectionEvent> BuildDirection(const std::vector<long double>& coefficient, bool triangle, unsigned outer, unsigned certifiedInner, bool graphOnly, bool wedge, const CutFemLevelSetData* prepared, const EndpointAxisProjection* projection) {
    const long double outerLower = 0;
    const long double innerLower = 0;
    std::vector<long double> events = {outerLower,1};
    const unsigned middle = 3U-outer-certifiedInner;
    const long double middleLower = 0;
    const bool projectedTriangle = triangle || (wedge && outer < 2 && middle < 2);
    const bool middleBoundary = triangle || (wedge && certifiedInner < 2 && middle < 2);
    const bool outerBoundary = triangle || (wedge && certifiedInner < 2 && outer < 2);
    std::vector<std::vector<std::vector<long double>>> constraints;
    std::vector<std::vector<std::vector<long double>>> normalDerivatives;
    std::vector<std::vector<std::vector<std::vector<long double>>>> foldFields;
    for(unsigned inner = 0; inner < 3; ++inner) if(inner != outer && (certifiedInner == 3 || inner == certifiedInner)) {
      const std::vector<std::vector<std::vector<long double>>> field = PreparedField(coefficient,outer,inner,prepared);
      for(unsigned side = 0; side < 2; ++side) {
        const std::vector<std::vector<long double>> face = Face(field,side != 0,triangle || wedge,middleBoundary,outerBoundary,innerLower);
        AddConstraint(face,projectedTriangle,events,graphOnly,middleLower,outerLower);
        constraints.push_back(face);
        normalDerivatives.push_back(FaceNormalDerivative(field,side != 0,triangle || wedge,middleBoundary,outerBoundary,innerLower));
        foldFields.push_back({});
      }
      if(graphOnly) continue;
      bool quadratic = false;
      bool linear = false;
      for(unsigned i = 0; i < 3; ++i) { quadratic = quadratic || PolynomialNorm(field[2][i]) != 0; linear = linear || PolynomialNorm(field[1][i]) != 0; }
      const std::vector<std::vector<long double>> discriminant = quadratic ? PreparedDiscriminant(field,outer,inner,prepared) : (linear ? field[1] : field[0]);
      if(!quadratic && linear) events.insert(events.end(),projection->projectedRoots[inner].begin(),projection->projectedRoots[inner].end());
      if(quadratic && discriminant.size() == 1 && PolynomialNorm(discriminant[0]) == 0) throw std::runtime_error("A repeated fiber root has no regular delta moment.");
      AddConstraint(discriminant,projectedTriangle,events,false,middleLower,outerLower,quadratic ? &projection->projectedRoots[inner] : nullptr,quadratic && projection->repeated[inner]);
      constraints.push_back(discriminant);
      normalDerivatives.push_back({});
      foldFields.push_back(quadratic ? field : std::vector<std::vector<std::vector<long double>>>());
    }
    if(projection != nullptr) for(const EndpointCriticalPoint& critical : projection->critical) events.push_back(critical.height);
    SortDistinct(events);
    std::vector<EndpointProjectionEvent> result;
    for(unsigned n = 0; n < events.size(); ++n) {
      bool logarithmic = false;
      if(projection != nullptr && projection->hasSaddle) for(const EndpointCriticalPoint& critical : projection->critical) if(critical.saddle && std::abs(critical.height-events[n]) <= 16*std::numeric_limits<long double>::epsilon()) { logarithmic = true; break; }
      unsigned exponent = 1;
      if(!logarithmic) for(unsigned i = 0; i < constraints.size(); ++i) exponent = CommonExponent(exponent,ContactExponent(constraints[i],events[n],projectedTriangle,graphOnly,normalDerivatives[i].empty() ? nullptr : &normalDerivatives[i],foldFields[i].empty() ? nullptr : &foldFields[i],wedge,middleBoundary,outerBoundary,innerLower,middleLower));
      result.push_back({events[n],exponent,logarithmic});
    }
    return result;
  }
  const std::vector<EndpointProjectionEvent>& Events(const CutFemLevelSetData& data, unsigned outer, unsigned inner, bool graphOnly, unsigned long long& searches, unsigned long long& builds) {
    EndpointAxisProjection& axis = axes[outer];
    std::string failure;
    for(unsigned direction = 0; direction < 3; ++direction) if(direction != outer && (inner == 3 || inner == direction)) {
      EndpointEventPlan& plan = axis.plans[graphOnly][direction];
      if(plan.events.empty() && plan.failure.empty()) {
        ++builds;
        const EndpointAxisProjection* critical = graphOnly ? nullptr : &CriticalPoints(data,outer,searches);
        try { plan.events = BuildDirection(coefficients,geometry == TET,outer,direction,graphOnly,geometry == WEDGE,&data,critical); }
        catch(const std::runtime_error& error) { plan.failure = error.what(); }
      }
      if(!plan.events.empty()) return plan.events;
      failure = plan.failure;
    }
    throw std::runtime_error(failure.empty() ? "No supported endpoint projection direction." : failure);
  }


};

}



#endif
