#include "include/CutFemLevelSetData.hpp"
#include "Basis.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>

unsigned long long checks = 0;
long double maximumScaledError = 0;
std::string currentCheck;

void Check(bool condition, const std::string& description) { ++checks; if(!condition) throw std::runtime_error(currentCheck+": "+description); }

void CheckNear(long double actual, long double expected, long double tolerance, const std::string& description) {
  ++checks;
  const long double error = std::abs(actual-expected)/(1+std::abs(expected));
  maximumScaledError = std::max(maximumScaledError,error);
  if(!std::isfinite(actual) || error > tolerance) { std::ostringstream message; message << std::setprecision(21) << currentCheck << ": " << description << ", actual=" << actual << ", expected=" << expected << ", scaled error=" << error; throw std::runtime_error(message.str()); }
}

std::unique_ptr<femus::basis> MakeBasis(GeomElType geometry) {
  if(geometry == LINE) return std::unique_ptr<femus::basis>(new femus::LineBiquadratic);
  if(geometry == QUAD) return std::unique_ptr<femus::basis>(new femus::QuadBiquadratic);
  if(geometry == TRI) return std::unique_ptr<femus::basis>(new femus::TriBiquadratic);
  if(geometry == HEX) return std::unique_ptr<femus::basis>(new femus::HexBiquadratic);
  if(geometry == TET) return std::unique_ptr<femus::basis>(new femus::TetBiquadratic);
  return std::unique_ptr<femus::basis>(new femus::WedgeBiquadratic);
}

std::string GeometryName(GeomElType geometry) { const std::string names[6] = {"HEX27","TET15","WEDGE21","QUAD9","TRI7","LINE3"}; return names[static_cast<unsigned>(geometry)]; }
unsigned Dimension(GeomElType geometry) { return geometry == LINE ? 1U : (geometry == QUAD || geometry == TRI ? 2U : 3U); }

std::vector<long double> RandomPoint(GeomElType geometry, std::mt19937& random) {
  std::uniform_real_distribution<double> uniform(0.01,0.99);
  std::vector<long double> point(Dimension(geometry));
  if(geometry == TRI || geometry == TET || geometry == WEDGE) { point[0] = uniform(random); point[1] = (1-point[0])*uniform(random); if(geometry == TET) point[2] = (1-point[0]-point[1])*uniform(random); if(geometry == WEDGE) point[2] = 2*uniform(random)-1; }
  else for(long double& coordinate : point) coordinate = 2*uniform(random)-1;
  return point;
}

long double ReferenceValue(const femus::basis& basis, const std::vector<double>& nodal, const std::vector<long double>& point, int derivative = -1) {
  double native[3] = {0,0,0};
  for(unsigned axis = 0; axis < point.size(); ++axis) native[axis] = static_cast<double>(point[axis]);
  long double result = 0;
  for(unsigned n = 0; n < nodal.size(); ++n) result += static_cast<long double>(nodal[n])*(derivative < 0 ? basis.eval_phi(basis.GetIND(n),native) : basis.eval_dphidxyz(static_cast<unsigned>(derivative),basis.GetIND(n),native));
  return result;
}

void TestReconstruction(unsigned cases, unsigned points, std::mt19937& random) {
  const GeomElType elements[6] = {LINE,QUAD,TRI,HEX,TET,WEDGE};
  std::uniform_real_distribution<double> value(-2,2);
  for(GeomElType geometry : elements) {
    currentCheck = "reconstruction "+GeometryName(geometry);
    std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    for(unsigned test = 0; test < cases; ++test) {
      std::vector<double> nodal(basis->_nc);
      for(double& coefficient : nodal) coefficient = value(random);
      CutFemLevelSetData data(geometry,nodal);
      Check(data.GetMaximumDepth() == 1,"lower default depth");
      for(unsigned n = 0; n < nodal.size(); ++n) { const double* coordinate = basis->GetXcoarse(n); std::vector<long double> point(Dimension(geometry)); for(unsigned axis = 0; axis < point.size(); ++axis) point[axis] = coordinate[axis]; CheckNear(data.Evaluate(point),nodal[n],2e-13L,"nodal interpolation"); }
      for(unsigned sample = 0; sample < points; ++sample) {
        const std::vector<long double> point = RandomPoint(geometry,random);
        const long double actual = data.Evaluate(point);
        CheckNear(actual,ReferenceValue(*basis,nodal,point),5e-13L,"actual FEMuS basis");
        CheckNear(data.EvaluateBernstein(point),actual,3e-15L,"Bernstein/power agreement");
        CheckNear(CutFemLevelSetData::EvaluatePolynomial(data.GetNormalizedPowerCoefficients(),data.GetDimension(),2,point)*data.GetNormalizationScale(),actual,3e-18L,"normalization");
        const std::array<long double,2> bounds = data.GetBounds();
        Check(actual >= bounds[0] && actual <= bounds[1],"Bernstein range bound");
        const std::vector<long double> gradient = data.GetGradient(point);
        const std::vector<long double> rootGradient = data.GetRootGradient(point);
        for(unsigned axis = 0; axis < point.size(); ++axis) { CheckNear(gradient[axis],ReferenceValue(*basis,nodal,point,static_cast<int>(axis)),2e-12L,"FEM derivative"); CheckNear(rootGradient[axis],gradient[axis],0,"root gradient"); const std::array<long double,2> derivativeBounds = data.GetDerivativeBounds(axis); Check(gradient[axis] >= derivativeBounds[0] && gradient[axis] <= derivativeBounds[1],"derivative range bound"); }
      }
      const std::array<unsigned long long,4> before = data.GetWorkCounts();
      data.GetPowerCoefficients();
      data.GetBernsteinCoefficients();
      data.GetNormalizedPowerCoefficients();
      data.GetGradient(RandomPoint(geometry,random));
      Check(data.GetWorkCounts() == before,"prepared polynomial reuse");
      data.SetContext(4,123456,789);
      Check(data.GetWorkCounts() == before,"context changes preserve coefficients");
      data.Reset(geometry,nodal,0);
      Check(data.GetWorkCounts()[0] == before[0] && data.GetWorkCounts()[1] == before[1],"same nodal data preserve coefficients");
      Check(data.GetMeshLevel() == 4 && data.GetElementId() == 123456 && data.GetNodalVersion() == 789,"context metadata");
    }
  }
}

void TestRestrictions(unsigned cases, unsigned points, std::mt19937& random) {
  const GeomElType elements[5] = {QUAD,TRI,HEX,TET,WEDGE};
  std::uniform_real_distribution<double> value(-1,1), position(0.02,0.98);
  for(GeomElType geometry : elements) {
    currentCheck = "restriction "+GeometryName(geometry);
    std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    for(unsigned test = 0; test < cases; ++test) {
      std::vector<double> nodal(basis->_nc);
      for(double& coefficient : nodal) coefficient = value(random);
      CutFemLevelSetData parent(geometry,nodal);
      parent.SetContext(2,34,56);
      CutFemLevelSetData child, line;
      for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) {
        const bool simplexAxis = geometry == TRI || geometry == TET || (geometry == WEDGE && axis < 2);
        const long double coordinate = simplexAxis ? position(random) : 2*position(random)-1;
        parent.Restrict(axis,coordinate,child);
        Check(child.GetDimension()+1 == parent.GetDimension(),"lower integration dimension");
        Check(child.GetRootDimension() == parent.GetRootDimension(),"retained root dimension");
        Check(child.GetNodalValues().empty() && child.GetWorkCounts()[0] == 0 && child.GetWorkCounts()[1] == 0,"coefficient-based child construction");
        Check(child.GetMeshLevel() == 2 && child.GetElementId() == 34 && child.GetNodalVersion() == 56,"slice context");
        for(unsigned sample = 0; sample < points; ++sample) {
          const std::vector<long double> point = RandomPoint(child.GetGeometry(),random);
          const std::vector<long double> parentPoint = child.MapToParent(point);
          CheckNear(parentPoint[axis],coordinate,0,"fixed coordinate");
          CheckNear(child.Evaluate(point),parent.Evaluate(parentPoint),3e-15L,"restricted power polynomial");
          CheckNear(child.EvaluateBernstein(point),child.Evaluate(point),8e-15L,"direct Bernstein restriction");
          const std::vector<long double> parentGradient = parent.GetGradient(parentPoint);
          const std::vector<long double> childGradient = child.GetGradient(point);
          const std::vector<long double> rootGradient = child.GetRootGradient(point);
          for(unsigned original = 0; original < parent.GetDimension(); ++original) CheckNear(rootGradient[original],parentGradient[original],4e-15L,"all original derivative components");
          for(unsigned local = 0; local < child.GetDimension(); ++local) { long double expected = 0; for(unsigned original = 0; original < parent.GetDimension(); ++original) expected += parentGradient[original]*child.GetParentMap()[original][local]; CheckNear(childGradient[local],expected,5e-15L,"reference-map chain rule"); }
          if(child.GetDimension() == 2) for(unsigned secondAxis = 0; secondAxis < 2; ++secondAxis) {
            const bool triangle = child.GetGeometry() == TRI;
            const long double secondCoordinate = triangle ? 0.37L : -0.23L;
            child.Restrict(secondAxis,secondCoordinate,line);
            const std::vector<long double> linePoint = {0.29L};
            const std::vector<long double> originalPoint = line.MapToRoot(linePoint);
            Check(line.GetGeometry() == LINE && line.GetRootDimension() == 3,"recursive same-type restriction");
            CheckNear(line.Evaluate(linePoint),parent.Evaluate(originalPoint),7e-15L,"recursive root-coordinate evaluation");
            CheckNear(line.EvaluateBernstein(linePoint),line.Evaluate(linePoint),1e-14L,"recursive Bernstein evaluation");
            const std::vector<long double> originalGradient = parent.GetGradient(originalPoint);
            const std::vector<long double> lineRootGradient = line.GetRootGradient(linePoint);
            for(unsigned original = 0; original < 3; ++original) CheckNear(lineRootGradient[original],originalGradient[original],8e-15L,"three components survive to LINE");
          }
        }
      }
      Check(parent.GetWorkCounts()[0] == 1 && parent.GetWorkCounts()[1] == 1,"one parent reconstruction for all slices");
    }
  }
}

long double ReferenceSecondDerivative(const femus::basis& basis, const std::vector<double>& nodal, const std::vector<long double>& point, unsigned first, unsigned second) {
  double native[3] = {0,0,0};
  for(unsigned axis = 0; axis < point.size(); ++axis) native[axis] = static_cast<double>(point[axis]);
  long double result = 0;
  for(unsigned n = 0; n < nodal.size(); ++n) {
    const int* index = basis.GetIND(n);
    double value;
    if(first == second) value = first == 0 ? basis.eval_d2phidx2(index,native) : (first == 1 ? basis.eval_d2phidy2(index,native) : basis.eval_d2phidz2(index,native));
    else if(first+second == 1) value = basis.eval_d2phidxdy(index,native);
    else if(first+second == 2) value = basis.eval_d2phidzdx(index,native);
    else value = basis.eval_d2phidydz(index,native);
    result += static_cast<long double>(nodal[n])*value;
  }
  return result;
}

void TestSecondDerivatives(unsigned cases, std::mt19937& random) {
  const GeomElType elements[6] = {LINE,QUAD,TRI,HEX,TET,WEDGE};
  std::uniform_real_distribution<double> value(-1,1);
  for(GeomElType geometry : elements) {
    currentCheck = "second derivatives "+GeometryName(geometry);
    std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    for(unsigned test = 0; test < cases; ++test) {
      std::vector<double> nodal(basis->_nc);
      for(double& coefficient : nodal) coefficient = value(random);
      CutFemLevelSetData data(geometry,nodal);
      const std::vector<long double> point = RandomPoint(geometry,random);
      for(unsigned first = 0; first < data.GetDimension(); ++first) for(unsigned second = 0; second < data.GetDimension(); ++second) { const long double actual = data.EvaluateSecondDerivative(first,second,point); CheckNear(actual,ReferenceSecondDerivative(*basis,nodal,point,first,second),3e-12L,"actual FEM Hessian"); CheckNear(actual,data.EvaluateSecondDerivative(second,first,point),0,"symmetric derivative cache"); }
      const std::array<unsigned long long,4> before = data.GetWorkCounts();
      for(unsigned first = 0; first < data.GetDimension(); ++first) for(unsigned second = 0; second < data.GetDimension(); ++second) data.EvaluateSecondDerivative(first,second,point);
      Check(data.GetWorkCounts() == before,"Hessian cache reuse");
    }
  }
}

std::vector<double> MakeField(GeomElType geometry, unsigned field, long double sign = 1, long double radius = 1e-6L) {
  std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
  std::vector<double> nodal(basis->_nc);
  const unsigned dimension = Dimension(geometry);
  const bool simplex = geometry == TRI || geometry == TET || geometry == WEDGE;
  const long double center[3] = {simplex ? 0.13L : 0.17L,simplex ? 0.17L : 0.23L,0.19L};
  for(unsigned n = 0; n < nodal.size(); ++n) {
    const double* point = basis->GetXcoarse(n);
    long double value = 1;
    if(field == 1) value = point[0]-(simplex ? 0.25L : 0.2L);
    if(field == 2) { value = -radius*radius; for(unsigned axis = 0; axis < dimension; ++axis) value += (point[axis]-center[axis])*(point[axis]-center[axis]); }
    if(field == 3) value = (point[0]-(simplex ? 0.5L : 0))*(point[0]-(simplex ? 0.5L : 0))+0.03L;
    if(field == 4) value = point[0]+(simplex ? 0 : 1);
    if(field == 5) value = (point[0]-(simplex ? 0.5L : 0))*(point[0]-(simplex ? 0.5L : 0));
    if(field == 6) value = 0;
    if(field == 7) value = dimension == 1 ? point[0] : point[0]-point[1];
    nodal[n] = static_cast<double>(sign*value);
  }
  return nodal;
}

void TestClassification() {
  const GeomElType elements[6] = {LINE,QUAD,TRI,HEX,TET,WEDGE};
  for(GeomElType geometry : elements) {
    currentCheck = "classification "+GeometryName(geometry);
    for(int sign = -1; sign <= 1; sign += 2) {
      CutFemLevelSetData constant(geometry,MakeField(geometry,0,sign));
      Check(constant.GetPhase() == sign,"constant phase");
      CutFemLevelSetData crossing(geometry,MakeField(geometry,1,sign));
      Check(crossing.GetPhase() == 0 && crossing.GetClassification() == "crossing","nodal crossing");
      Check(!crossing.HasPowerCoefficients() && !crossing.HasBernsteinCoefficients(),"nodal signs checked first");
      CutFemLevelSetData contact(geometry,MakeField(geometry,4,sign));
      Check(contact.GetPhase() == 0 && contact.GetClassification() == "contact","nodal contact");
      Check(!contact.HasPowerCoefficients() && !contact.HasBernsteinCoefficients(),"contact fast path");
      const std::vector<double> pocketValues = MakeField(geometry,2,sign);
      for(double value : pocketValues) Check(sign*value > 0,"small cut has uniform nodal signs");
      CutFemLevelSetData pocket(geometry,pocketValues,0);
      Check(pocket.IsCutCandidate(),"small interior cut survives a single Bernstein test");
      pocket.SetMaximumDepth(12);
      Check(pocket.IsCutCandidate(),"small interior cut survives subdivision");
      CutFemLevelSetData positive(geometry,MakeField(geometry,3,sign),0);
      Check(positive.IsCutCandidate(),"mixed Bernstein coefficients remain uncertified");
      const std::array<unsigned long long,4> before = positive.GetWorkCounts();
      positive.SetMaximumDepth(12);
      Check(positive.GetPhase() == sign,"subdivision removes a false alarm");
      Check(positive.GetWorkCounts() == before,"depth change reuses parent coefficients");
      CutFemLevelSetData tangent(geometry,MakeField(geometry,5,sign),2);
      Check(tangent.IsCutCandidate(),"tangent contact remains a candidate");
    }
    CutFemLevelSetData zero(geometry,MakeField(geometry,6));
    Check(zero.IsCutCandidate() && zero.GetClassification() == "contact","identically zero classification");
    CheckNear(zero.GetNormalizationScale(),1,0,"zero normalization");
    const std::vector<long double> center(Dimension(geometry),0);
    CheckNear(zero.Evaluate(center),0,0,"zero polynomial");
    CutFemLevelSetData plane(geometry,MakeField(geometry,1));
    Check(plane.GetGraphDirection() == 0,"certified graph direction");
  }
  currentCheck = "classification zero Bernstein coefficient";
  CutFemLevelSetData positive(LINE,{2,2,1},0);
  Check(positive.IsCutCandidate(),"a zero control is not proof of a crossing");
  positive.SetMaximumDepth(1);
  Check(positive.GetPhase() == 1,"positive polynomial certified after one split");
}

void TestBoundarySlices() {
  const GeomElType elements[5] = {QUAD,TRI,HEX,TET,WEDGE};
  for(GeomElType geometry : elements) {
    currentCheck = "boundary slices "+GeometryName(geometry);
    CutFemLevelSetData parent(geometry,MakeField(geometry,1));
    CutFemLevelSetData child, line;
    for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) {
      const bool simplexAxis = geometry == TRI || geometry == TET || (geometry == WEDGE && axis < 2);
      const long double endpoints[4] = {simplexAxis ? 0.0L : -1.0L,1.0L,std::nextafter(1.0L,0.0L),simplexAxis ? 1e-15L : -1+1e-15L};
      for(long double coordinate : endpoints) {
        parent.Restrict(axis,coordinate,child);
        const std::vector<long double> point(child.GetDimension(),0);
        CheckNear(child.Evaluate(point),parent.Evaluate(child.MapToParent(point)),5e-17L,"endpoint power restriction");
        CheckNear(child.EvaluateBernstein(point),child.Evaluate(point),5e-16L,"endpoint Bernstein restriction");
        const bool collapsed = simplexAxis && coordinate == 1;
        Check(child.IsDegenerate() == collapsed,"degenerate slice flag");
        long double expected = 1;
        if(geometry == TET) expected = (1-coordinate)*(1-coordinate);
        if(geometry == TRI || (geometry == WEDGE && axis < 2)) expected = (1-coordinate)/2;
        CheckNear(child.GetParentMeasureScale(),expected,0,"slice measure");
        if(child.GetDimension() == 2) { child.Restrict(0,child.GetGeometry() == TRI ? 1.0L : -1.0L,line); const std::vector<long double> linePoint = {0.3L}; CheckNear(line.Evaluate(linePoint),parent.Evaluate(line.MapToRoot(linePoint)),2e-16L,"recursive endpoint restriction"); CheckNear(line.EvaluateBernstein(linePoint),line.Evaluate(linePoint),8e-16L,"recursive endpoint Bernstein restriction"); }
      }
    }
  }
}

void TestComponents(unsigned cases, std::mt19937& random) {
  const GeomElType elements[5] = {QUAD,TRI,HEX,TET,WEDGE};
  std::uniform_real_distribution<double> value(-0.1,0.1);
  for(GeomElType geometry : elements) {
    currentCheck = "external surface components "+GeometryName(geometry);
    CutFemLevelSetData parent(geometry,MakeField(geometry,1));
    for(unsigned degree = 2; degree <= 6; ++degree) {
      unsigned count = 1;
      for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) count *= degree+1;
      std::vector<std::vector<long double>> components(3,std::vector<long double>(count)), output;
      for(std::vector<long double>& component : components) for(long double& coefficient : component) coefficient = value(random);
      CutFemLevelSetData child;
      for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) {
        const bool simplexAxis = geometry == TRI || geometry == TET || (geometry == WEDGE && axis < 2);
        const long double coordinate = simplexAxis ? 0.31L : -0.27L;
        parent.Restrict(axis,coordinate,child);
        parent.RestrictComponents(axis,coordinate,degree,components,output);
        Check(output.size() == 3,"all physical components retained");
        for(unsigned sample = 0; sample < cases; ++sample) {
          const std::vector<long double> point = RandomPoint(child.GetGeometry(),random);
          const std::vector<long double> parentPoint = child.MapToParent(point);
          for(unsigned component = 0; component < 3; ++component) CheckNear(CutFemLevelSetData::EvaluatePolynomial(output[component],child.GetDimension(),degree,point),CutFemLevelSetData::EvaluatePolynomial(components[component],parent.GetDimension(),degree,parentPoint),1e-14L,"full-degree component contraction");
        }
      }
    }
  }
}

void TestOwnershipAndInvalidation() {
  currentCheck = "ownership and invalidation";
  std::vector<double> nodal = MakeField(HEX,1);
  std::unique_ptr<femus::basis> basis = MakeBasis(HEX);
  for(unsigned n = 0; n < nodal.size(); ++n) { const double* coordinate = basis->GetXcoarse(n); nodal[n] = coordinate[0]*coordinate[0]*coordinate[2]+coordinate[1]*coordinate[2]*coordinate[2]+0.25*coordinate[0]-0.125; }
  CutFemLevelSetData parent(HEX,nodal), child;
  const std::vector<long double> point = {0.11L,0.22L,0.33L};
  const long double original = parent.Evaluate(point);
  nodal.assign(27,42);
  CheckNear(parent.Evaluate(point),original,0,"owned nodal snapshot");
  parent.Restrict(2,0.17L,child);
  const std::vector<long double> childPoint = {0.13L,-0.29L};
  const long double saved = child.Evaluate(childPoint);
  const std::vector<long double> savedGradient = child.GetRootGradient(childPoint);
  CutFemLevelSetData copied = parent;
  parent.Reset(HEX,nodal);
  CheckNear(parent.Evaluate(point),42,0,"new nodal state");
  CheckNear(copied.Evaluate(point),original,0,"independent copy");
  CheckNear(child.Evaluate(childPoint),saved,0,"slice snapshot survives parent reset");
  Check(child.GetRootGradient(childPoint) == savedGradient,"slice gradient survives parent reset");
  CutFemLevelSetData reusable;
  copied.Restrict(2,-0.4L,reusable);
  const long double* powerStorage = reusable.GetPowerCoefficients().data();
  const long double* bernsteinStorage = reusable.GetBernsteinCoefficients().data();
  const long double* derivativeStorage = reusable.GetDerivativeCoefficients(0).data();
  reusable.EvaluateSecondDerivative(0,0,childPoint);
  const std::uint64_t previous = reusable.GetRevision();
  copied.Restrict(2,0.6L,reusable);
  Check(reusable.GetRevision() == previous+1,"slice revision changes");
  Check(reusable.GetPowerCoefficients().data() == powerStorage && reusable.GetBernsteinCoefficients().data() == bernsteinStorage,"coefficient output storage reused");
  Check(reusable.GetDerivativeCoefficients(0).data() == derivativeStorage,"derivative storage reused after invalidation");
  CheckNear(reusable.EvaluateSecondDerivative(0,0,childPoint),copied.EvaluateSecondDerivative(0,0,reusable.MapToRoot(childPoint)),0,"Hessian invalidated and rebuilt");
  Check(reusable.GetWorkCounts()[0] == 0 && reusable.GetWorkCounts()[1] == 0,"no nodal or Bernstein reconversion for slices");
}

void TestPrecision() {
  currentCheck = "precision and line roots";
  Check(std::numeric_limits<long double>::digits >= std::numeric_limits<double>::digits,"working precision");
  CutFemLevelSetData cancellation(LINE,{1e16,1,-1e16});
  if(std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits) CheckNear(cancellation.Evaluate({1}),1,0,"extended precision retains cancellation residual");
  const double radii[5] = {0.2,1e-3,1e-6,1e-12,1e-100};
  for(double radius : radii) {
    const double squared = radius*radius;
    const std::vector<double> nodal = {1-squared,1-squared,-squared};
    CutFemLevelSetData line(LINE,nodal);
    const std::vector<long double> roots = line.GetLineRoots();
    Check(roots.size() == 2,"two regular roots");
    const boost::multiprecision::cpp_bin_float_quad a = (boost::multiprecision::cpp_bin_float_quad(nodal[0])+boost::multiprecision::cpp_bin_float_quad(nodal[1]))/2-boost::multiprecision::cpp_bin_float_quad(nodal[2]);
    const boost::multiprecision::cpp_bin_float_quad reference = sqrt(-boost::multiprecision::cpp_bin_float_quad(nodal[2])/a);
    const long double expected = static_cast<long double>(reference);
    CheckNear(roots[0]/expected,-1,3e-18L,"negative root versus multiprecision");
    CheckNear(roots[1]/expected,1,3e-18L,"positive root versus multiprecision");
  }
  CutFemLevelSetData tangent(LINE,{1,1,0});
  Check(tangent.GetLineRoots() == std::vector<long double>{0},"repeated root");
  CutFemLevelSetData endpoints(LINE,{0,0,-1});
  Check(endpoints.GetLineRoots() == std::vector<long double>({-1,1}),"boundary roots");
  CutFemLevelSetData empty(LINE,{1,1,1});
  Check(empty.GetLineRoots().empty(),"constant has no roots");
  const double scales[5] = {1e-200,1e-40,1,1e40,1e200};
  for(GeomElType geometry : {LINE,QUAD,TRI,HEX,TET,WEDGE}) {
    const std::vector<double> base = MakeField(geometry,1);
    for(double scale : scales) {
      std::vector<double> nodal = base;
      for(double& value : nodal) value *= scale;
      CutFemLevelSetData scaled(geometry,nodal);
      Check(scaled.GetPhase() == 0,"scale-invariant crossing classification");
      const std::vector<long double> point(Dimension(geometry),0);
      const long double expected = ReferenceValue(*MakeBasis(geometry),nodal,point);
      CheckNear(scaled.Evaluate(point)/scale,expected/scale,2e-13L,"scaled input reconstruction");
      CheckNear(scaled.EvaluateBernstein(point)/scale,scaled.Evaluate(point)/scale,3e-15L,"scaled Bernstein reconstruction");
    }
  }
}

void TestSliceRoundoffMargins() {
  currentCheck = "inherited slice roundoff margins";
  const GeomElType elements[5] = {QUAD,TRI,HEX,TET,WEDGE};
  for(GeomElType geometry : elements) {
    const long double offset = geometry == TRI || geometry == TET || geometry == WEDGE ? 0.25L : 0.2L;
    CutFemLevelSetData parent(geometry,MakeField(geometry,1)), child, line;
    parent.Restrict(0,offset,child);
    Check(child.IsCutCandidate(),"cancellation slice keeps the parent uncertainty");
    const std::array<long double,2> bounds = child.GetBounds();
    Check(bounds[0] <= 0 && bounds[1] >= 0,"slice bounds retain a near-zero interval");
    if(child.GetDimension() == 2) { child.Restrict(0,0.17L,line); const std::array<long double,2> lineBounds = line.GetBounds(); Check(line.IsCutCandidate() && lineBounds[0] <= 0 && lineBounds[1] >= 0,"recursive slices retain the uncertainty"); }
  }
  CutFemLevelSetData asymmetric(LINE,{2,1e-30,1e-30});
  const std::vector<long double> roots = asymmetric.GetLineRoots();
  Check(roots.size() == 2,"asymmetric roots");
  const boost::multiprecision::cpp_bin_float_quad c(1e-30);
  const boost::multiprecision::cpp_bin_float_quad a = 1-c/2;
  const boost::multiprecision::cpp_bin_float_quad b = -1+c/2;
  const boost::multiprecision::cpp_bin_float_quad discriminant = sqrt(b*b-4*a*c);
  const boost::multiprecision::cpp_bin_float_quad q = (-b+discriminant)/2;
  CheckNear(roots[0]/static_cast<long double>(c/q),1,3e-18L,"small root avoids subtraction cancellation");
  CheckNear(roots[1],static_cast<long double>(q/a),3e-18L,"large asymmetric root");
}

void TestGKRestrictions(std::mt19937& random) {
  currentCheck = "nonuniform Gauss-Kronrod slices";
  const GeomElType elements[5] = {QUAD,TRI,HEX,TET,WEDGE};
  const std::vector<long double> abscissa(boost::math::quadrature::gauss_kronrod<long double,15>::abscissa().begin(),boost::math::quadrature::gauss_kronrod<long double,15>::abscissa().end());
  for(GeomElType geometry : elements) {
    CutFemLevelSetData parent(geometry,MakeField(geometry,2));
    CutFemLevelSetData child;
    for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) for(long double node : abscissa) for(int sign = -1; sign <= 1; sign += 2) {
      const bool simplexAxis = geometry == TRI || geometry == TET || (geometry == WEDGE && axis < 2);
      const long double coordinate = simplexAxis ? (1+sign*node)/2 : sign*node;
      parent.Restrict(axis,coordinate,child);
      const std::vector<long double> point = RandomPoint(child.GetGeometry(),random);
      CheckNear(child.Evaluate(point),parent.Evaluate(child.MapToParent(point)),4e-17L,"GK point coefficient restriction");
      CheckNear(child.EvaluateBernstein(point),child.Evaluate(point),3e-16L,"GK point Bernstein restriction");
    }
    Check(parent.GetWorkCounts()[0] == 1 && parent.GetWorkCounts()[1] == 1,"GK sweeps reuse parent coefficients");
  }
}

void TestInvalidInputs() {
  currentCheck = "input validation";
  CutFemLevelSetData empty;
  bool rejected = false;
  try { empty.GetPhase(); } catch(const std::logic_error&) { rejected = true; }
  Check(rejected && !empty.IsInitialized(),"uninitialized access rejected");
  rejected = false;
  try { CutFemLevelSetData invalid(POINT,{}); } catch(const std::invalid_argument&) { rejected = true; }
  Check(rejected,"unsupported geometry rejected");
  const GeomElType elements[6] = {LINE,QUAD,TRI,HEX,TET,WEDGE};
  for(GeomElType geometry : elements) {
    const std::vector<double> nodal = MakeField(geometry,1);
    CutFemLevelSetData data(geometry,nodal);
    const std::vector<long double> point(Dimension(geometry),0);
    const long double saved = data.Evaluate(point);
    std::vector<double> invalid = nodal;
    invalid.pop_back();
    rejected = false;
    try { data.Reset(geometry,invalid); } catch(const std::invalid_argument&) { rejected = true; }
    Check(rejected,"incorrect node count rejected");
    CheckNear(data.Evaluate(point),saved,0,"failed reset preserves previous state");
    for(double value : {std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) { invalid = nodal; invalid[0] = value; rejected = false; try { data.Reset(geometry,invalid); } catch(const std::invalid_argument&) { rejected = true; } Check(rejected,"nonfinite nodal input rejected"); }
    rejected = false;
    try { data.SetMaximumDepth(33); } catch(const std::invalid_argument&) { rejected = true; }
    Check(rejected,"depth limit validated");
    rejected = false;
    try { data.Evaluate({}); } catch(const std::invalid_argument&) { rejected = true; }
    Check(rejected,"evaluation dimension validated");
    std::vector<long double> outside = point;
    outside[0] = 2;
    rejected = false;
    try { data.Evaluate(outside); } catch(const std::out_of_range&) { rejected = true; }
    Check(rejected,"reference domain validated");
    rejected = false;
    try { data.GetDerivativeCoefficients(data.GetDimension()); } catch(const std::out_of_range&) { rejected = true; }
    Check(rejected,"derivative axis validated");
    rejected = false;
    try { data.Restrict(0,0,data); } catch(const std::invalid_argument&) { rejected = true; }
    Check(rejected,"parent-child alias rejected");
    CutFemLevelSetData child;
    rejected = false;
    try { data.Restrict(data.GetDimension(),0,child); } catch(const std::exception&) { rejected = true; }
    Check(rejected,"restriction dimension validated");
    if(data.GetDimension() > 1) {
      rejected = false;
      try { data.Restrict(0,2,child); } catch(const std::out_of_range&) { rejected = true; }
      Check(rejected,"restriction coordinate validated");
      rejected = false;
      try { data.Restrict(0,std::numeric_limits<long double>::quiet_NaN(),child); } catch(const std::invalid_argument&) { rejected = true; }
      Check(rejected,"nonfinite slice rejected");
      std::vector<std::vector<long double>> components(3,std::vector<long double>(1));
      rejected = false;
      try { data.RestrictComponents(0,0,2,components,components); } catch(const std::invalid_argument&) { rejected = true; }
      Check(rejected,"component alias rejected");
      std::vector<std::vector<long double>> output;
      rejected = false;
      try { data.RestrictComponents(0,0,2,components,output); } catch(const std::invalid_argument&) { rejected = true; }
      Check(rejected,"component coefficient count validated");
    }
    else { rejected = false; try { data.Restrict(0,0,child); } catch(const std::logic_error&) { rejected = true; } Check(rejected,"LINE termination validated"); }
  }
  CutFemLevelSetData zero(LINE,{0,0,0});
  rejected = false;
  try { zero.GetLineRoots(); } catch(const std::domain_error&) { rejected = true; }
  Check(rejected,"undefined isolated roots rejected");
  rejected = false;
  try { CutFemLevelSetData::EvaluatePolynomial({1,2},1,2,{0}); } catch(const std::invalid_argument&) { rejected = true; }
  Check(rejected,"external polynomial layout validated");
}

std::array<long double,2> IntegratePrepared(CutFemLevelSetData& data, const std::vector<long double>& points, const std::vector<long double>& weights) {
  std::array<long double,2> result = {0,0};
  if(data.GetDimension() == 1) {
    const std::vector<long double> roots = data.GetLineRoots();
    std::vector<long double> cuts = {-1};
    for(long double root : roots) if(root > -1 && root < 1) cuts.push_back(root);
    cuts.push_back(1);
    for(unsigned i = 1; i < cuts.size(); ++i) if(data.Evaluate({(cuts[i-1]+cuts[i])/2}) > 0) result[0] += cuts[i]-cuts[i-1];
    for(long double root : roots) { const std::vector<long double> gradient = data.GetRootGradient({root}); long double norm = 0; for(long double value : gradient) norm = std::hypot(norm,value); const long double derivative = std::abs(data.EvaluateDerivative(0,{root})); if(derivative > 0) result[1] += (root == -1 || root == 1 ? 0.5L : 1.0L)*norm/derivative; }
    return result;
  }
  const unsigned axis = data.GetDimension()-1;
  const bool triangle = data.GetGeometry() == TRI;
  const bool tetrahedron = data.GetGeometry() == TET;
  const std::vector<long double> panels = triangle ? std::vector<long double>{0,0.5L,1} : (tetrahedron ? std::vector<long double>{0,1} : std::vector<long double>{-1,1});
  CutFemLevelSetData child;
  for(unsigned panel = 1; panel < panels.size(); ++panel) {
    const long double origin = (panels[panel]+panels[panel-1])/2;
    const long double scale = (panels[panel]-panels[panel-1])/2;
    for(unsigned n = 0; n < points.size(); ++n) {
      data.Restrict(axis,std::fma(scale,points[n],origin),child);
      const std::array<long double,2> inner = IntegratePrepared(child,points,weights);
      const long double weight = weights[n]*scale*child.GetParentMeasureScale();
      for(unsigned component = 0; component < 2; ++component) result[component] += weight*inner[component];
    }
  }
  return result;
}

void TestNestedIntegration() {
  currentCheck = "nested volume and surface integration";
  const std::vector<long double> abscissa(boost::math::quadrature::gauss_kronrod<long double,15>::abscissa().begin(),boost::math::quadrature::gauss_kronrod<long double,15>::abscissa().end());
  const std::vector<long double> kronrodWeights(boost::math::quadrature::gauss_kronrod<long double,15>::weights().begin(),boost::math::quadrature::gauss_kronrod<long double,15>::weights().end());
  std::vector<long double> points, weights;
  for(unsigned i = 0; i < abscissa.size(); ++i) { points.push_back(abscissa[i]); weights.push_back(kronrodWeights[i]); if(abscissa[i] != 0) { points.push_back(-abscissa[i]); weights.push_back(kronrodWeights[i]); } }
  const GeomElType elements[6] = {LINE,QUAD,TRI,HEX,TET,WEDGE};
  const long double volumes[6] = {2,4,0.5L,8,1.0L/6,1};
  const long double surfaces[6] = {1,2*std::sqrt(2.0L),1/std::sqrt(2.0L),4*std::sqrt(2.0L),1/(2*std::sqrt(2.0L)),std::sqrt(2.0L)};
  for(unsigned element = 0; element < 6; ++element) {
    CutFemLevelSetData data(elements[element],MakeField(elements[element],7));
    const std::array<long double,2> integrals = IntegratePrepared(data,points,weights);
    CheckNear(integrals[0],volumes[element]/2,3e-15L,"half-domain volume "+GeometryName(elements[element]));
    CheckNear(integrals[1],surfaces[element],3e-14L,"interface measure "+GeometryName(elements[element]));
    Check(data.GetWorkCounts()[0] == 1,"one nodal reconstruction through nested integration");
  }
}

template<class Real> Real BenchmarkValue(const std::vector<Real>& coefficient, unsigned dimension, const std::vector<Real>& point) {
  Real result = 0;
  if(dimension == 1) return (coefficient[2]*point[0]+coefficient[1])*point[0]+coefficient[0];
  for(unsigned i = 3; i > 0; --i) { Real row = 0; for(unsigned j = 3; j > 0; --j) { Real line; if(dimension == 2) line = coefficient[3*(i-1)+j-1]; else line = (coefficient[(3*(i-1)+j-1)*3+2]*point[2]+coefficient[(3*(i-1)+j-1)*3+1])*point[2]+coefficient[(3*(i-1)+j-1)*3]; row = j == 3 ? line : row*point[1]+line; } result = i == 3 ? row : result*point[0]+row; }
  return result;
}

boost::multiprecision::cpp_bin_float_quad QuadReference(const std::vector<long double>& coefficient, unsigned dimension, const std::vector<double>& point) {
  boost::multiprecision::cpp_bin_float_quad result = 0;
  for(unsigned n = 0; n < coefficient.size(); ++n) { unsigned index = n; std::array<unsigned,3> exponent = {0,0,0}; for(unsigned axis = dimension; axis > 0; --axis) { exponent[axis-1] = index%3; index /= 3; } boost::multiprecision::cpp_bin_float_quad term(coefficient[n]); for(unsigned axis = 0; axis < dimension; ++axis) for(unsigned power = 0; power < exponent[axis]; ++power) term *= boost::multiprecision::cpp_bin_float_quad(point[axis]); result += term; }
  return result;
}

void TestPrecisionPaths(bool benchmark, std::mt19937& random) {
  currentCheck = "double/long-double polynomial arithmetic";
  const GeomElType elements[6] = {LINE,QUAD,TRI,HEX,TET,WEDGE};
  std::uniform_real_distribution<double> value(-1,1);
  for(GeomElType geometry : elements) {
    std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    std::vector<double> nodal(basis->_nc);
    for(double& coefficient : nodal) coefficient = value(random);
    CutFemLevelSetData data(geometry,nodal);
    const std::vector<long double>& coefficient = data.GetPowerCoefficients();
    const std::vector<double> doubleCoefficient(coefficient.begin(),coefficient.end());
    std::vector<std::vector<double>> doublePoints;
    std::vector<std::vector<long double>> extendedPoints;
    long double doubleError = 0, extendedError = 0;
    for(unsigned sample = 0; sample < 64; ++sample) {
      const std::vector<long double> point = RandomPoint(geometry,random);
      doublePoints.push_back(std::vector<double>(point.begin(),point.end()));
      extendedPoints.push_back(std::vector<long double>(doublePoints.back().begin(),doublePoints.back().end()));
      const boost::multiprecision::cpp_bin_float_quad reference = QuadReference(coefficient,data.GetDimension(),doublePoints.back());
      const double doubleValue = BenchmarkValue(doubleCoefficient,data.GetDimension(),doublePoints.back());
      const long double extendedValue = BenchmarkValue(coefficient,data.GetDimension(),extendedPoints.back());
      doubleError = std::max(doubleError,static_cast<long double>(abs(boost::multiprecision::cpp_bin_float_quad(doubleValue)-reference)/(1+abs(reference))));
      extendedError = std::max(extendedError,static_cast<long double>(abs(boost::multiprecision::cpp_bin_float_quad(extendedValue)-reference)/(1+abs(reference))));
      CheckNear(extendedValue,static_cast<long double>(reference),2e-16L,"extended evaluation versus multiprecision");
      CheckNear(doubleValue,static_cast<long double>(reference),3e-13L,"double evaluation versus multiprecision");
    }
    if(!benchmark) continue;
    volatile double doubleChecksum = 0;
    volatile long double extendedChecksum = 0;
    const std::chrono::steady_clock::time_point doubleStart = std::chrono::steady_clock::now();
    for(unsigned iteration = 0; iteration < 2000; ++iteration) for(const std::vector<double>& point : doublePoints) doubleChecksum += BenchmarkValue(doubleCoefficient,data.GetDimension(),point);
    const double doubleSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-doubleStart).count();
    const std::chrono::steady_clock::time_point extendedStart = std::chrono::steady_clock::now();
    for(unsigned iteration = 0; iteration < 2000; ++iteration) for(const std::vector<long double>& point : extendedPoints) extendedChecksum += BenchmarkValue(coefficient,data.GetDimension(),point);
    const double extendedSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-extendedStart).count();
    Check(std::isfinite(doubleChecksum) && std::isfinite(extendedChecksum),"benchmark checksums");
    std::cout << std::setprecision(8) << "precision " << GeometryName(geometry) << ": double_error=" << doubleError << ", extended_error=" << extendedError << ", double_seconds=" << doubleSeconds << ", extended_seconds=" << extendedSeconds << '\n';
  }
}

void PhysicalMap(const std::vector<long double>& point, bool nonAffine, std::vector<long double>& physical, std::vector<std::vector<long double>>& jacobian) {
  const unsigned dimension = point.size();
  const long double linear[3][3] = {{1,0.13L,0.04L},{0.03L,1,0.08L},{-0.06L,0.02L,1}};
  const long double mixed[3][4] = {{0.07L,0.03L,0.05L,0.04L},{-0.04L,0.06L,-0.03L,0.02L},{0.03L,-0.02L,0.04L,-0.03L}};
  const unsigned exponents[4][3] = {{1,1,0},{1,0,1},{0,1,1},{1,1,1}};
  physical.assign(dimension,0);
  jacobian.assign(dimension,std::vector<long double>(dimension,0));
  for(unsigned component = 0; component < dimension; ++component) {
    for(unsigned axis = 0; axis < dimension; ++axis) { physical[component] += linear[component][axis]*point[axis]; jacobian[component][axis] = linear[component][axis]; }
    if(!nonAffine) continue;
    for(unsigned term = 0; term < 4; ++term) {
      if(dimension == 2 && exponents[term][2] != 0) continue;
      long double value = mixed[component][term];
      for(unsigned axis = 0; axis < dimension; ++axis) if(exponents[term][axis]) value *= point[axis];
      physical[component] += value;
      for(unsigned derivative = 0; derivative < dimension; ++derivative) if(exponents[term][derivative]) { value = mixed[component][term]; for(unsigned axis = 0; axis < dimension; ++axis) if(axis != derivative && exponents[term][axis]) value *= point[axis]; jacobian[component][derivative] += value; }
    }
  }
}

void TestPhysicalSphere(unsigned cases, std::mt19937& random) {
  currentCheck = "physical circle/sphere on affine and multilinear geometry";
  const long double center[3] = {0.13L,-0.17L,0.11L};
  for(GeomElType geometry : {QUAD,HEX}) for(unsigned mode = 0; mode < 2; ++mode) {
    std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    std::vector<double> nodal(basis->_nc);
    for(unsigned n = 0; n < nodal.size(); ++n) { const double* native = basis->GetXcoarse(n); std::vector<long double> point(Dimension(geometry)), physical; std::vector<std::vector<long double>> jacobian; for(unsigned axis = 0; axis < point.size(); ++axis) point[axis] = native[axis]; PhysicalMap(point,mode != 0,physical,jacobian); long double value = 0.35L*0.35L; for(unsigned axis = 0; axis < point.size(); ++axis) value -= (physical[axis]-center[axis])*(physical[axis]-center[axis]); nodal[n] = static_cast<double>(value); }
    CutFemLevelSetData data(geometry,nodal), child;
    for(unsigned sample = 0; sample < cases*4; ++sample) {
      const std::vector<long double> point = RandomPoint(geometry,random);
      std::vector<long double> physical;
      std::vector<std::vector<long double>> jacobian;
      PhysicalMap(point,mode != 0,physical,jacobian);
      long double expected = 0.35L*0.35L;
      for(unsigned component = 0; component < point.size(); ++component) expected -= (physical[component]-center[component])*(physical[component]-center[component]);
      CheckNear(data.Evaluate(point),expected,8e-15L,"physical sphere pullback");
      CheckNear(data.EvaluateBernstein(point),expected,1e-14L,"physical sphere Bernstein coefficients");
      const std::vector<long double> gradient = data.GetGradient(point);
      for(unsigned axis = 0; axis < point.size(); ++axis) { long double derivative = 0; for(unsigned component = 0; component < point.size(); ++component) derivative -= 2*(physical[component]-center[component])*jacobian[component][axis]; CheckNear(gradient[axis],derivative,2e-14L,"physical sphere pullback gradient"); }
      data.Restrict(data.GetDimension()-1,point.back(),child);
      const std::vector<long double> childPoint(point.begin(),point.end()-1);
      CheckNear(child.Evaluate(childPoint),expected,8e-15L,"physical sphere slice");
    }
  }
}

int main(int argc, char** argv) {
  try {
    unsigned cases = 64, points = 32;
    bool benchmark = false;
    for(int i = 1; i < argc; ++i) { const std::string option = argv[i]; if((option == "--cases" || option == "--points") && i+1 < argc) { const std::string text = argv[++i]; if(text.empty() || text.find_first_not_of("0123456789") != std::string::npos) throw std::invalid_argument("Test counts must be positive integers."); const unsigned long value = std::stoul(text); if(value == 0 || value > 10000) throw std::invalid_argument("Test counts must lie between 1 and 10000."); if(option == "--cases") cases = static_cast<unsigned>(value); else points = static_cast<unsigned>(value); } else if(option == "--benchmark") benchmark = true; else if(option == "--help") { std::cout << "MultiPhase_ex47 [--cases N] [--points N] [--benchmark]\n"; return 0; } else throw std::invalid_argument("Unknown or incomplete argument: "+option); }
    std::mt19937 random(470013);
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    TestReconstruction(cases,points,random);
    TestRestrictions(cases,points,random);
    TestSecondDerivatives(cases,random);
    TestClassification();
    TestBoundarySlices();
    TestComponents(cases,random);
    TestOwnershipAndInvalidation();
    TestPrecision();
    TestSliceRoundoffMargins();
    TestGKRestrictions(random);
    TestInvalidInputs();
    TestNestedIntegration();
    TestPrecisionPaths(benchmark,random);
    TestPhysicalSphere(cases,random);
    std::cout << std::setprecision(8) << "ex47: " << checks << " checks passed; maximum scaled error=" << maximumScaledError << "; seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count() << '\n';
    return 0;
  }
  catch(const std::exception& error) { std::cerr << "ex47 failed after " << checks << " checks: " << error.what() << '\n'; return 1; }
}
