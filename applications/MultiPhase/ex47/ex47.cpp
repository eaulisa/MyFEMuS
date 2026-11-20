#include "CutFemLevelSetData.hpp"
#include "CutFemWeight_LS.hpp"
#include "Basis.hpp"
#include "ReferenceSurfaceMoments.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <cstdlib>
#include <new>

#ifdef EX47_TRACK_ALLOCATIONS
unsigned long long heapAllocations = 0;
void* operator new(std::size_t bytes) { void* result = std::malloc(bytes == 0 ? 1 : bytes); if(result == nullptr) throw std::bad_alloc(); ++heapAllocations; return result; }
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept { std::free(pointer); }
#endif

unsigned long long checks = 0;
long double maximumScaledError = 0;
std::string currentCheck;

void Check(bool condition, const std::string& description) {
  ++checks;
  if(!condition) throw std::runtime_error(currentCheck + ": " + description);
}

void CheckNear(long double actual, long double expected, long double tolerance, const std::string& description) {
  ++checks;
  const long double error = std::abs(actual - expected) / (1 + std::abs(expected));
  maximumScaledError = std::max(maximumScaledError, error);
  if(!std::isfinite(actual) || error > tolerance) {
    std::ostringstream message;
    message << std::setprecision(21) << currentCheck << ": " << description << ", actual=" << actual << ", expected=" << expected << ", scaled error=" << error;
    throw std::runtime_error(message.str());
  }
}

std::unique_ptr<femus::basis> MakeBasis(GeomElType geometry) {
  if(geometry == LINE) return std::unique_ptr<femus::basis>(new femus::LineBiquadratic);
  if(geometry == QUAD) return std::unique_ptr<femus::basis>(new femus::QuadBiquadratic);
  if(geometry == TRI) return std::unique_ptr<femus::basis>(new femus::TriBiquadratic);
  if(geometry == HEX) return std::unique_ptr<femus::basis>(new femus::HexBiquadratic);
  if(geometry == TET) return std::unique_ptr<femus::basis>(new femus::TetBiquadratic);
  return std::unique_ptr<femus::basis>(new femus::WedgeBiquadratic);
}

std::string GeometryName(GeomElType geometry) {
  const std::string names[6] = {"HEX27", "TET15", "WEDGE21", "QUAD9", "TRI7", "LINE3"};
  return names[static_cast<unsigned>(geometry)];
}
unsigned Dimension(GeomElType geometry) {
  return geometry == LINE ? 1U : (geometry == QUAD || geometry == TRI ? 2U : 3U);
}

bool TensorAxis(const femus::basis& basis, unsigned axis) {
  return dynamic_cast<const femus::LineBiquadratic*>(&basis) != nullptr || dynamic_cast<const femus::QuadBiquadratic*>(&basis) != nullptr || dynamic_cast<const femus::HexBiquadratic*>(&basis) != nullptr || (dynamic_cast<const femus::WedgeBiquadratic*>(&basis) != nullptr && axis == 2);
}

long double UnitNode(const femus::basis& basis, unsigned node, unsigned axis) { return TensorAxis(basis,axis) ? (basis.GetXcoarse(node)[axis]+1)/2 : basis.GetXcoarse(node)[axis]; }

std::vector<long double> RandomPoint(GeomElType geometry, std::mt19937& random) {
  std::uniform_real_distribution<double> uniform(0.01, 0.99);
  std::vector<long double> point(Dimension(geometry));
  if(geometry == TRI || geometry == TET || geometry == WEDGE) {
    point[0] = uniform(random);
    point[1] = (1 - point[0]) * uniform(random);
    if(geometry == TET) point[2] = (1 - point[0] - point[1]) * uniform(random);
    if(geometry == WEDGE) point[2] = uniform(random);
  }
  else for(long double& coordinate : point) coordinate = uniform(random);
  return point;
}

long double ReferenceValue(const femus::basis& basis, const std::vector<double>& nodal, const std::vector<long double>& point, int derivative = -1) {
  double native[3] = {0, 0, 0};
  for(unsigned axis = 0; axis < point.size(); ++axis) native[axis] = static_cast<double>(TensorAxis(basis,axis) ? 2*point[axis]-1 : point[axis]);
  long double result = 0;
  for(unsigned n = 0; n < nodal.size(); ++n) result += static_cast<long double>(nodal[n]) * (derivative < 0 ? basis.eval_phi(basis.GetIND(n), native) : basis.eval_dphidxyz(static_cast<unsigned>(derivative), basis.GetIND(n), native));
  return derivative >= 0 && TensorAxis(basis,static_cast<unsigned>(derivative)) ? 2*result : result;
}

void TestReconstruction(unsigned cases, unsigned points, std::mt19937& random) {
  const GeomElType elements[6] = {LINE, QUAD, TRI, HEX, TET, WEDGE};
  std::uniform_real_distribution<double> value(-2, 2);
  for(GeomElType geometry : elements) {
    currentCheck = "reconstruction " + GeometryName(geometry);
    std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    for(unsigned test = 0; test < cases; ++test) {
      std::vector<double> nodal(basis->_nc);
      for(double& coefficient : nodal) coefficient = value(random);
      CutFemLevelSetData data(geometry, nodal);
      Check(data.GetMaximumDepth() == 1, "lower default depth");
      for(unsigned n = 0; n < nodal.size(); ++n) {
        const double* coordinate = basis->GetXcoarse(n);
        std::vector<long double> point(Dimension(geometry));
        for(unsigned axis = 0; axis < point.size(); ++axis) point[axis] = UnitNode(*basis,n,axis);
        CheckNear(data.Evaluate(point), nodal[n], 2e-13L, "nodal interpolation");
      }
      for(unsigned sample = 0; sample < points; ++sample) {
        const std::vector<long double> point = RandomPoint(geometry, random);
        const long double actual = data.Evaluate(point);
        CheckNear(actual, ReferenceValue(*basis, nodal, point), 5e-13L, "actual FEMuS basis");
        CheckNear(data.EvaluateBernstein(point), actual, 3e-15L, "Bernstein/power agreement");
        CheckNear(CutFemLevelSetData::EvaluatePolynomial(data.GetNormalizedPowerCoefficients(), data.GetDimension(), 2, point)*data.GetNormalizationScale(), actual, 3e-18L, "normalization");
        const std::array<long double, 2> bounds = data.GetBounds();
        Check(actual >= bounds[0] && actual <= bounds[1], "Bernstein range bound");
        const std::vector<long double> gradient = data.GetGradient(point);
        const std::vector<long double> rootGradient = data.GetRootGradient(point);
        for(unsigned axis = 0; axis < point.size(); ++axis) {
          CheckNear(gradient[axis], ReferenceValue(*basis, nodal, point, static_cast<int>(axis)), 2e-12L, "FEM derivative");
          CheckNear(rootGradient[axis], gradient[axis], 0, "root gradient");
          const std::array<long double, 2> derivativeBounds = data.GetDerivativeBounds(axis);
          Check(gradient[axis] >= derivativeBounds[0] && gradient[axis] <= derivativeBounds[1], "derivative range bound");
        }
      }
      const std::array<unsigned long long, 4> before = data.GetWorkCounts();
      data.GetPowerCoefficients();
      data.GetBernsteinCoefficients();
      data.GetNormalizedPowerCoefficients();
      data.GetGradient(RandomPoint(geometry, random));
      Check(data.GetWorkCounts() == before, "prepared polynomial reuse");
      data.SetContext(4, 123456, 789);
      Check(data.GetWorkCounts() == before, "context changes preserve coefficients");
      data.Reset(geometry, nodal, 0);
      Check(data.GetWorkCounts()[0] == before[0] && data.GetWorkCounts()[1] == before[1], "same nodal data preserve coefficients");
      Check(data.GetMeshLevel() == 4 && data.GetElementId() == 123456 && data.GetNodalVersion() == 789, "context metadata");
    }
  }
}

void TestRestrictions(unsigned cases, unsigned points, std::mt19937& random) {
  const GeomElType elements[5] = {QUAD, TRI, HEX, TET, WEDGE};
  std::uniform_real_distribution<double> value(-1, 1), position(0.02, 0.98);
  for(GeomElType geometry : elements) {
    currentCheck = "restriction " + GeometryName(geometry);
    std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    for(unsigned test = 0; test < cases; ++test) {
      std::vector<double> nodal(basis->_nc);
      for(double& coefficient : nodal) coefficient = value(random);
      CutFemLevelSetData parent(geometry, nodal);
      parent.SetContext(2, 34, 56);
      CutFemLevelSetData child, line;
      for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) {
        const bool simplexAxis = geometry == TRI || geometry == TET || (geometry == WEDGE && axis < 2);
        const long double coordinate = position(random);
        parent.Restrict(axis, coordinate, child);
        Check(child.GetDimension() + 1 == parent.GetDimension(), "lower integration dimension");
        Check(child.GetRootDimension() == parent.GetRootDimension(), "retained root dimension");
        Check(child.GetNodalValues().empty() && child.GetWorkCounts()[0] == 0 && child.GetWorkCounts()[1] == 0, "coefficient-based child construction");
        Check(child.GetMeshLevel() == 2 && child.GetElementId() == 34 && child.GetNodalVersion() == 56, "slice context");
        for(unsigned sample = 0; sample < points; ++sample) {
          const std::vector<long double> point = RandomPoint(child.GetGeometry(), random);
          const std::vector<long double> parentPoint = child.MapToParent(point);
          CheckNear(parentPoint[axis], coordinate, 0, "fixed coordinate");
          CheckNear(child.Evaluate(point), parent.Evaluate(parentPoint), 3e-15L, "restricted power polynomial");
          CheckNear(child.EvaluateBernstein(point), child.Evaluate(point), 8e-15L, "direct Bernstein restriction");
          const std::vector<long double> parentGradient = parent.GetGradient(parentPoint);
          const std::vector<long double> childGradient = child.GetGradient(point);
          const std::vector<long double> rootGradient = child.GetRootGradient(point);
          for(unsigned original = 0; original < parent.GetDimension(); ++original) CheckNear(rootGradient[original], parentGradient[original], 4e-15L, "all original derivative components");
          for(unsigned local = 0; local < child.GetDimension(); ++local) {
            long double expected = 0;
            for(unsigned original = 0; original < parent.GetDimension(); ++original) expected += parentGradient[original] * child.GetParentMap()[original][local];
            CheckNear(childGradient[local], expected, 5e-15L, "reference-map chain rule");
          }
          if(child.GetDimension() == 2) for(unsigned secondAxis = 0; secondAxis < 2; ++secondAxis) {
              const bool triangle = child.GetGeometry() == TRI;
              const long double secondCoordinate = 0.37L;
              child.Restrict(secondAxis, secondCoordinate, line);
              const std::vector<long double> linePoint = {0.29L};
              const std::vector<long double> originalPoint = line.MapToRoot(linePoint);
              Check(line.GetGeometry() == LINE && line.GetRootDimension() == 3, "recursive same-type restriction");
              CheckNear(line.Evaluate(linePoint), parent.Evaluate(originalPoint), 7e-15L, "recursive root-coordinate evaluation");
              CheckNear(line.EvaluateBernstein(linePoint), line.Evaluate(linePoint), 1e-14L, "recursive Bernstein evaluation");
              const std::vector<long double> originalGradient = parent.GetGradient(originalPoint);
              const std::vector<long double> lineRootGradient = line.GetRootGradient(linePoint);
              for(unsigned original = 0; original < 3; ++original) CheckNear(lineRootGradient[original], originalGradient[original], 8e-15L, "three components survive to LINE");
            }
        }
      }
      Check(parent.GetWorkCounts()[0] == 1 && parent.GetWorkCounts()[1] == 1, "one parent reconstruction for all slices");
    }
  }
}

long double ReferenceSecondDerivative(const femus::basis& basis, const std::vector<double>& nodal, const std::vector<long double>& point, unsigned first, unsigned second) {
  double native[3] = {0, 0, 0};
  for(unsigned axis = 0; axis < point.size(); ++axis) native[axis] = static_cast<double>(TensorAxis(basis,axis) ? 2*point[axis]-1 : point[axis]);
  long double result = 0;
  for(unsigned n = 0; n < nodal.size(); ++n) {
    const int* index = basis.GetIND(n);
    double value;
    if(first == second) value = first == 0 ? basis.eval_d2phidx2(index, native) : (first == 1 ? basis.eval_d2phidy2(index, native) : basis.eval_d2phidz2(index, native));
    else if(first + second == 1) value = basis.eval_d2phidxdy(index, native);
    else if(first + second == 2) value = basis.eval_d2phidzdx(index, native);
    else value = basis.eval_d2phidydz(index, native);
    result += static_cast<long double>(nodal[n]) * value;
  }
  return result*(TensorAxis(basis,first) ? 2 : 1)*(TensorAxis(basis,second) ? 2 : 1);
}

void TestSecondDerivatives(unsigned cases, std::mt19937& random) {
  const GeomElType elements[6] = {LINE, QUAD, TRI, HEX, TET, WEDGE};
  std::uniform_real_distribution<double> value(-1, 1);
  for(GeomElType geometry : elements) {
    currentCheck = "second derivatives " + GeometryName(geometry);
    std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    for(unsigned test = 0; test < cases; ++test) {
      std::vector<double> nodal(basis->_nc);
      for(double& coefficient : nodal) coefficient = value(random);
      CutFemLevelSetData data(geometry, nodal);
      const std::vector<long double> point = RandomPoint(geometry, random);
      for(unsigned first = 0; first < data.GetDimension(); ++first) for(unsigned second = 0; second < data.GetDimension(); ++second) {
          const long double actual = data.EvaluateSecondDerivative(first, second, point);
          CheckNear(actual, ReferenceSecondDerivative(*basis, nodal, point, first, second), 3e-12L, "actual FEM Hessian");
          CheckNear(actual, data.EvaluateSecondDerivative(second, first, point), 0, "symmetric derivative cache");
        }
      const std::array<unsigned long long, 4> before = data.GetWorkCounts();
      for(unsigned first = 0; first < data.GetDimension(); ++first) for(unsigned second = 0; second < data.GetDimension(); ++second) data.EvaluateSecondDerivative(first, second, point);
      Check(data.GetWorkCounts() == before, "Hessian cache reuse");
    }
  }
}

std::vector<double> MakeField(GeomElType geometry, unsigned field, long double sign = 1, long double radius = 1e-6L) {
  std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
  std::vector<double> nodal(basis->_nc);
  const unsigned dimension = Dimension(geometry);
  const bool simplex = geometry == TRI || geometry == TET || geometry == WEDGE;
  const long double center[3] = {simplex ? 0.13L : 0.17L, simplex ? 0.17L : 0.23L, 0.19L};
  for(unsigned n = 0; n < nodal.size(); ++n) {
    const double* point = basis->GetXcoarse(n);
    long double value = 1;
    if(field == 1) value = point[0] - (simplex ? 0.25L : 0.2L);
    if(field == 2) {
      value = -radius * radius;
      for(unsigned axis = 0; axis < dimension; ++axis) value += (point[axis] - center[axis]) * (point[axis] - center[axis]);
    }
    if(field == 3) value = (point[0] - (simplex ? 0.5L : 0)) * (point[0] - (simplex ? 0.5L : 0)) + 0.03L;
    if(field == 4) value = point[0] + (simplex ? 0 : 1);
    if(field == 5) value = (point[0] - (simplex ? 0.5L : 0)) * (point[0] - (simplex ? 0.5L : 0));
    if(field == 6) value = 0;
    if(field == 7) value = dimension == 1 ? point[0] : point[0] - point[1];
    nodal[n] = static_cast<double>(sign * value);
  }
  return nodal;
}

void TestClassification() {
  const GeomElType elements[6] = {LINE, QUAD, TRI, HEX, TET, WEDGE};
  for(GeomElType geometry : elements) {
    currentCheck = "classification " + GeometryName(geometry);
    for(int sign = -1; sign <= 1; sign += 2) {
      CutFemLevelSetData constant(geometry, MakeField(geometry, 0, sign));
      Check(constant.GetPhase() == sign, "constant phase");
      CutFemLevelSetData crossing(geometry, MakeField(geometry, 1, sign));
      Check(crossing.GetPhase() == 0 && crossing.GetClassification() == "crossing", "nodal crossing");
      Check(!crossing.HasPowerCoefficients() && !crossing.HasBernsteinCoefficients(), "nodal signs checked first");
      CutFemLevelSetData contact(geometry, MakeField(geometry, 4, sign));
      Check(contact.GetPhase() == 0 && contact.GetClassification() == "contact", "nodal contact");
      Check(!contact.HasPowerCoefficients() && !contact.HasBernsteinCoefficients(), "contact fast path");
      const std::vector<double> pocketValues = MakeField(geometry, 2, sign);
      for(double value : pocketValues) Check(sign * value > 0, "small cut has uniform nodal signs");
      CutFemLevelSetData pocket(geometry, pocketValues, 0);
      Check(pocket.IsCutCandidate(), "small interior cut survives a single Bernstein test");
      pocket.SetMaximumDepth(12);
      Check(pocket.IsCutCandidate(), "small interior cut survives subdivision");
      CutFemLevelSetData positive(geometry, MakeField(geometry, 3, sign), 0);
      Check(positive.IsCutCandidate(), "mixed Bernstein coefficients remain uncertified");
      const std::array<unsigned long long, 4> before = positive.GetWorkCounts();
      positive.SetMaximumDepth(12);
      Check(positive.GetPhase() == sign, "subdivision removes a false alarm");
      Check(positive.GetWorkCounts() == before, "depth change reuses parent coefficients");
      CutFemLevelSetData tangent(geometry, MakeField(geometry, 5, sign), 2);
      Check(tangent.IsCutCandidate(), "tangent contact remains a candidate");
    }
    CutFemLevelSetData zero(geometry, MakeField(geometry, 6));
    Check(zero.IsCutCandidate() && zero.GetClassification() == "contact", "identically zero classification");
    CheckNear(zero.GetNormalizationScale(), 1, 0, "zero normalization");
    const std::vector<long double> center(Dimension(geometry), 0);
    CheckNear(zero.Evaluate(center), 0, 0, "zero polynomial");
    CutFemLevelSetData plane(geometry, MakeField(geometry, 1));
    Check(plane.GetGraphDirection() == 0, "certified graph direction");
  }
  currentCheck = "classification zero Bernstein coefficient";
  CutFemLevelSetData positive(LINE, {2, 2, 1}, 0);
  Check(positive.IsCutCandidate(), "a zero control is not proof of a crossing");
  positive.SetMaximumDepth(1);
  Check(positive.GetPhase() == 1, "positive polynomial certified after one split");
}

void TestBoundarySlices() {
  const GeomElType elements[5] = {QUAD, TRI, HEX, TET, WEDGE};
  for(GeomElType geometry : elements) {
    currentCheck = "boundary slices " + GeometryName(geometry);
    CutFemLevelSetData parent(geometry, MakeField(geometry, 1));
    CutFemLevelSetData child, line;
    for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) {
      const bool simplexAxis = geometry == TRI || geometry == TET || (geometry == WEDGE && axis < 2);
      const long double endpoints[4] = {0.0L, 1.0L, std::nextafter(1.0L, 0.0L), 1e-15L};
      for(long double coordinate : endpoints) {
        parent.Restrict(axis, coordinate, child);
        const std::vector<long double> point(child.GetDimension(), 0);
        CheckNear(child.Evaluate(point), parent.Evaluate(child.MapToParent(point)), 5e-17L, "endpoint power restriction");
        CheckNear(child.EvaluateBernstein(point), child.Evaluate(point), 5e-16L, "endpoint Bernstein restriction");
        const bool collapsed = simplexAxis && coordinate == 1;
        Check(child.IsDegenerate() == collapsed, "degenerate slice flag");
        long double expected = 1;
        if(geometry == TET) expected = (1 - coordinate) * (1 - coordinate);
        if(geometry == TRI || (geometry == WEDGE && axis < 2)) expected = (1 - coordinate);
        CheckNear(child.GetParentMeasureScale(), expected, 0, "slice measure");
        if(child.GetDimension() == 2) {
          child.Restrict(0, child.GetGeometry() == TRI ? 1.0L : 0.0L, line);
          const std::vector<long double> linePoint = {0.3L};
          CheckNear(line.Evaluate(linePoint), parent.Evaluate(line.MapToRoot(linePoint)), 2e-16L, "recursive endpoint restriction");
          CheckNear(line.EvaluateBernstein(linePoint), line.Evaluate(linePoint), 8e-16L, "recursive endpoint Bernstein restriction");
        }
      }
    }
  }
}

void TestComponents(unsigned cases, std::mt19937& random) {
  const GeomElType elements[5] = {QUAD, TRI, HEX, TET, WEDGE};
  std::uniform_real_distribution<double> value(-0.1, 0.1);
  for(GeomElType geometry : elements) {
    currentCheck = "external surface components " + GeometryName(geometry);
    CutFemLevelSetData parent(geometry, MakeField(geometry, 1));
    for(unsigned degree = 2; degree <= 6; ++degree) {
      unsigned count = 1;
      for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) count *= degree + 1;
      std::vector<std::vector<long double>> components(3, std::vector<long double>(count)), output;
      for(std::vector<long double>& component : components) for(long double& coefficient : component) coefficient = value(random);
      CutFemLevelSetData child;
      for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) {
        const bool simplexAxis = geometry == TRI || geometry == TET || (geometry == WEDGE && axis < 2);
        const long double coordinate = 0.31L;
        parent.Restrict(axis, coordinate, child);
        parent.RestrictComponents(axis, coordinate, degree, components, output);
        Check(output.size() == 3, "all physical components retained");
        for(unsigned sample = 0; sample < cases; ++sample) {
          const std::vector<long double> point = RandomPoint(child.GetGeometry(), random);
          const std::vector<long double> parentPoint = child.MapToParent(point);
          for(unsigned component = 0; component < 3; ++component) CheckNear(CutFemLevelSetData::EvaluatePolynomial(output[component], child.GetDimension(), degree, point), CutFemLevelSetData::EvaluatePolynomial(components[component], parent.GetDimension(), degree, parentPoint), 1e-14L, "full-degree component contraction");
        }
      }
    }
  }
}

void TestOwnershipAndInvalidation() {
  currentCheck = "ownership and invalidation";
  std::vector<double> nodal = MakeField(HEX, 1);
  std::unique_ptr<femus::basis> basis = MakeBasis(HEX);
  for(unsigned n = 0; n < nodal.size(); ++n) {
    const double* coordinate = basis->GetXcoarse(n);
    nodal[n] = coordinate[0] * coordinate[0] * coordinate[2] + coordinate[1] * coordinate[2] * coordinate[2] + 0.25 * coordinate[0] - 0.125;
  }
  CutFemLevelSetData parent(HEX, nodal), child;
  const std::vector<long double> point = {0.11L, 0.22L, 0.33L};
  const long double original = parent.Evaluate(point);
  nodal.assign(27, 42);
  CheckNear(parent.Evaluate(point), original, 0, "owned nodal snapshot");
  parent.Restrict(2, 0.17L, child);
  const std::vector<long double> childPoint = {0.13L, 0.29L};
  const long double saved = child.Evaluate(childPoint);
  const std::vector<long double> savedGradient = child.GetRootGradient(childPoint);
  CutFemLevelSetData copied = parent;
  parent.Reset(HEX, nodal);
  CheckNear(parent.Evaluate(point), 42, 0, "new nodal state");
  CheckNear(copied.Evaluate(point), original, 0, "independent copy");
  CheckNear(child.Evaluate(childPoint), saved, 0, "slice snapshot survives parent reset");
  Check(child.GetRootGradient(childPoint) == savedGradient, "slice gradient survives parent reset");
  CutFemLevelSetData reusable;
  copied.Restrict(2, 0.3L, reusable);
  const long double* powerStorage = reusable.GetPowerCoefficients().data();
  const long double* bernsteinStorage = reusable.GetBernsteinCoefficients().data();
  const long double* derivativeStorage = reusable.GetDerivativeCoefficients(0).data();
  reusable.EvaluateSecondDerivative(0, 0, childPoint);
  const std::uint64_t previous = reusable.GetRevision();
  copied.Restrict(2, 0.6L, reusable);
  Check(reusable.GetRevision() == previous + 1, "slice revision changes");
  Check(reusable.GetPowerCoefficients().data() == powerStorage && reusable.GetBernsteinCoefficients().data() == bernsteinStorage, "coefficient output storage reused");
  Check(reusable.GetDerivativeCoefficients(0).data() == derivativeStorage, "derivative storage reused after invalidation");
  CheckNear(reusable.EvaluateSecondDerivative(0, 0, childPoint), copied.EvaluateSecondDerivative(0, 0, reusable.MapToRoot(childPoint)), 0, "Hessian invalidated and rebuilt");
  Check(reusable.GetWorkCounts()[0] == 0 && reusable.GetWorkCounts()[1] == 0, "no nodal or Bernstein reconversion for slices");
}

void TestPrecision() {
  currentCheck = "precision and line roots";
  Check(std::numeric_limits<long double>::digits >= std::numeric_limits<double>::digits, "working precision");
  CutFemLevelSetData cancellation(LINE, {1e16, 1, -1e16});
  if(std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits) CheckNear(cancellation.Evaluate({1}), 1, 0, "extended precision retains cancellation residual");
  const double radii[5] = {0.2, 1e-3, 1e-6, 1e-12, 1e-100};
  for(double radius : radii) {
    const double squared = radius * radius;
    const std::vector<double> nodal = {1 - squared, 1 - squared, -squared};
    CutFemLevelSetData line(LINE, nodal);
    const std::vector<long double> roots = line.GetLineRoots();
    Check(roots.size() == 2, "two regular roots");
    const boost::multiprecision::cpp_bin_float_quad a = (boost::multiprecision::cpp_bin_float_quad(nodal[0]) + boost::multiprecision::cpp_bin_float_quad(nodal[1])) / 2 - boost::multiprecision::cpp_bin_float_quad(nodal[2]);
    const boost::multiprecision::cpp_bin_float_quad reference = sqrt(-boost::multiprecision::cpp_bin_float_quad(nodal[2]) / a);
    const long double expected = static_cast<long double>(reference);
    CheckNear(roots[0],0.5L*(1-expected),3e-18L,"lower unit root versus multiprecision");
    CheckNear(roots[1],0.5L*(1+expected),3e-18L,"upper unit root versus multiprecision");
  }
  CutFemLevelSetData tangent(LINE, {1, 1, 0});
  Check(tangent.GetLineRoots() == std::vector<long double> {0.5L}, "repeated root");
  CutFemLevelSetData endpoints(LINE, {0, 0, -1});
  Check(endpoints.GetLineRoots() == std::vector<long double>({0, 1}), "boundary roots");
  CutFemLevelSetData empty(LINE, {1, 1, 1});
  Check(empty.GetLineRoots().empty(), "constant has no roots");
  const double scales[5] = {1e-200, 1e-40, 1, 1e40, 1e200};
  for(GeomElType geometry : {
        LINE, QUAD, TRI, HEX, TET, WEDGE
      }) {
    const std::vector<double> base = MakeField(geometry, 1);
    for(double scale : scales) {
      std::vector<double> nodal = base;
      for(double& value : nodal) value *= scale;
      CutFemLevelSetData scaled(geometry, nodal);
      Check(scaled.GetPhase() == 0, "scale-invariant crossing classification");
      const std::vector<long double> point(Dimension(geometry), 0);
      const long double expected = ReferenceValue(*MakeBasis(geometry), nodal, point);
      CheckNear(scaled.Evaluate(point) / scale, expected / scale, 2e-13L, "scaled input reconstruction");
      CheckNear(scaled.EvaluateBernstein(point) / scale, scaled.Evaluate(point) / scale, 3e-15L, "scaled Bernstein reconstruction");
    }
  }
}

void TestSliceRoundoffMargins() {
  currentCheck = "inherited slice roundoff margins";
  const GeomElType elements[5] = {QUAD, TRI, HEX, TET, WEDGE};
  for(GeomElType geometry : elements) {
    const long double offset = geometry == TRI || geometry == TET || geometry == WEDGE ? 0.25L : 0.6L;
    CutFemLevelSetData parent(geometry, MakeField(geometry, 1)), child, line;
    parent.Restrict(0, offset, child);
    Check(child.IsCutCandidate(), "cancellation slice keeps the parent uncertainty");
    const std::array<long double, 2> bounds = child.GetBounds();
    Check(bounds[0] <= 0 && bounds[1] >= 0, "slice bounds retain a near-zero interval");
    if(child.GetDimension() == 2) {
      child.Restrict(0, 0.17L, line);
      const std::array<long double, 2> lineBounds = line.GetBounds();
      Check(line.IsCutCandidate() && lineBounds[0] <= 0 && lineBounds[1] >= 0, "recursive slices retain the uncertainty");
    }
  }
  CutFemLevelSetData asymmetric(LINE, {2, 1e-30, 1e-30});
  const std::vector<long double> roots = asymmetric.GetLineRoots();
  Check(roots.size() == 2, "asymmetric roots");
  const boost::multiprecision::cpp_bin_float_quad c(1e-30);
  const boost::multiprecision::cpp_bin_float_quad a = 1 - c / 2;
  const boost::multiprecision::cpp_bin_float_quad b = -1 + c / 2;
  const boost::multiprecision::cpp_bin_float_quad discriminant = sqrt(b * b - 4 * a * c);
  const boost::multiprecision::cpp_bin_float_quad q = (-b + discriminant) / 2;
  CheckNear(roots[0],0.5L*(1+static_cast<long double>(c/q)),3e-18L,"small native root mapped to unit coordinates");
  CheckNear(roots[1], 0.5L*(1+static_cast<long double>(q/a)), 3e-18L, "large asymmetric root");
}

void TestGKRestrictions(std::mt19937& random) {
  currentCheck = "nonuniform Gauss-Kronrod slices";
  const GeomElType elements[5] = {QUAD, TRI, HEX, TET, WEDGE};
  const std::vector<long double> abscissa(boost::math::quadrature::gauss_kronrod<long double, 15>::abscissa().begin(), boost::math::quadrature::gauss_kronrod<long double, 15>::abscissa().end());
  for(GeomElType geometry : elements) {
    CutFemLevelSetData parent(geometry, MakeField(geometry, 2));
    CutFemLevelSetData child;
    for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) for(long double node : abscissa) for(int sign = -1; sign <= 1; sign += 2) {
          const bool simplexAxis = geometry == TRI || geometry == TET || (geometry == WEDGE && axis < 2);
          const long double coordinate = (1 + sign * node) / 2;
          parent.Restrict(axis, coordinate, child);
          const std::vector<long double> point = RandomPoint(child.GetGeometry(), random);
          CheckNear(child.Evaluate(point), parent.Evaluate(child.MapToParent(point)), 4e-17L, "GK point coefficient restriction");
          CheckNear(child.EvaluateBernstein(point), child.Evaluate(point), 3e-16L, "GK point Bernstein restriction");
        }
    Check(parent.GetWorkCounts()[0] == 1 && parent.GetWorkCounts()[1] == 1, "GK sweeps reuse parent coefficients");
  }
}

void TestInvalidInputs() {
  currentCheck = "input validation";
  CutFemLevelSetData empty;
  bool rejected = false;
  try {
    empty.GetPhase();
  }
  catch(const std::logic_error&) {
    rejected = true;
  }
  Check(rejected && !empty.IsInitialized(), "uninitialized access rejected");
  rejected = false;
  try {
    CutFemLevelSetData invalid(POINT, {});
  }
  catch(const std::invalid_argument&) {
    rejected = true;
  }
  Check(rejected, "unsupported geometry rejected");
  const GeomElType elements[6] = {LINE, QUAD, TRI, HEX, TET, WEDGE};
  for(GeomElType geometry : elements) {
    const std::vector<double> nodal = MakeField(geometry, 1);
    CutFemLevelSetData data(geometry, nodal);
    const std::vector<long double> point(Dimension(geometry), 0);
    const long double saved = data.Evaluate(point);
    std::vector<double> invalid = nodal;
    invalid.pop_back();
    rejected = false;
    try {
      data.Reset(geometry, invalid);
    }
    catch(const std::invalid_argument&) {
      rejected = true;
    }
    Check(rejected, "incorrect node count rejected");
    CheckNear(data.Evaluate(point), saved, 0, "failed reset preserves previous state");
    for(double value : {
          std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()
        }) {
      invalid = nodal;
      invalid[0] = value;
      rejected = false;
      try {
        data.Reset(geometry, invalid);
      }
      catch(const std::invalid_argument&) {
        rejected = true;
      }
      Check(rejected, "nonfinite nodal input rejected");
    }
    rejected = false;
    try {
      data.SetMaximumDepth(33);
    }
    catch(const std::invalid_argument&) {
      rejected = true;
    }
    Check(rejected, "depth limit validated");
    rejected = false;
    try {
      data.Evaluate({});
    }
    catch(const std::invalid_argument&) {
      rejected = true;
    }
    Check(rejected, "evaluation dimension validated");
    std::vector<long double> outside = point;
    outside[0] = 2;
    rejected = false;
    try {
      data.Evaluate(outside);
    }
    catch(const std::out_of_range&) {
      rejected = true;
    }
    Check(rejected, "reference domain validated");
    rejected = false;
    try {
      data.GetDerivativeCoefficients(data.GetDimension());
    }
    catch(const std::out_of_range&) {
      rejected = true;
    }
    Check(rejected, "derivative axis validated");
    rejected = false;
    try {
      data.Restrict(0, 0, data);
    }
    catch(const std::invalid_argument&) {
      rejected = true;
    }
    Check(rejected, "parent-child alias rejected");
    CutFemLevelSetData child;
    rejected = false;
    try {
      data.Restrict(data.GetDimension(), 0, child);
    }
    catch(const std::exception&) {
      rejected = true;
    }
    Check(rejected, "restriction dimension validated");
    if(data.GetDimension() > 1) {
      rejected = false;
      try {
        data.Restrict(0, 2, child);
      }
      catch(const std::out_of_range&) {
        rejected = true;
      }
      Check(rejected, "restriction coordinate validated");
      rejected = false;
      try {
        data.Restrict(0, std::numeric_limits<long double>::quiet_NaN(), child);
      }
      catch(const std::invalid_argument&) {
        rejected = true;
      }
      Check(rejected, "nonfinite slice rejected");
      std::vector<std::vector<long double>> components(3, std::vector<long double>(1));
      rejected = false;
      try {
        data.RestrictComponents(0, 0, 2, components, components);
      }
      catch(const std::invalid_argument&) {
        rejected = true;
      }
      Check(rejected, "component alias rejected");
      std::vector<std::vector<long double>> output;
      rejected = false;
      try {
        data.RestrictComponents(0, 0, 2, components, output);
      }
      catch(const std::invalid_argument&) {
        rejected = true;
      }
      Check(rejected, "component coefficient count validated");
    }
    else {
      rejected = false;
      try {
        data.Restrict(0, 0, child);
      }
      catch(const std::logic_error&) {
        rejected = true;
      }
      Check(rejected, "LINE termination validated");
    }
  }
  CutFemLevelSetData zero(LINE, {0, 0, 0});
  rejected = false;
  try {
    zero.GetLineRoots();
  }
  catch(const std::domain_error&) {
    rejected = true;
  }
  Check(rejected, "undefined isolated roots rejected");
  rejected = false;
  try {
    CutFemLevelSetData::EvaluatePolynomial({1, 2}, 1, 2, {0});
  }
  catch(const std::invalid_argument&) {
    rejected = true;
  }
  Check(rejected, "external polynomial layout validated");
}

std::array<long double, 2> IntegratePrepared(CutFemLevelSetData& data, const std::vector<long double>& points, const std::vector<long double>& weights) {
  std::array<long double, 2> result = {0, 0};
  if(data.GetDimension() == 1) {
    const std::vector<long double> roots = data.GetLineRoots();
    std::vector<long double> cuts = {0};
    for(long double root : roots) if(root > 0 && root < 1) cuts.push_back(root);
    cuts.push_back(1);
    for(unsigned i = 1; i < cuts.size(); ++i) if(data.Evaluate({(cuts[i - 1] + cuts[i]) / 2}) > 0) result[0] += cuts[i] - cuts[i - 1];
    for(long double root : roots) {
      const std::vector<long double> gradient = data.GetRootGradient({root});
      long double norm = 0;
      for(long double value : gradient) norm = std::hypot(norm, value);
      const long double derivative = std::abs(data.EvaluateDerivative(0, {root}));
      if(derivative > 0) result[1] += (root == 0 || root == 1 ? 0.5L : 1.0L) * norm / derivative;
    }
    return result;
  }
  const unsigned axis = data.GetDimension() - 1;
  const bool triangle = data.GetGeometry() == TRI;
  const bool tetrahedron = data.GetGeometry() == TET;
  const std::vector<long double> panels = triangle ? std::vector<long double> {0, 0.5L, 1} :
                                          (tetrahedron ? std::vector<long double> {0, 1} : std::vector<long double> {0, 1});
  CutFemLevelSetData child;
  for(unsigned panel = 1; panel < panels.size(); ++panel) {
    const long double origin = (panels[panel] + panels[panel - 1]) / 2;
    const long double scale = (panels[panel] - panels[panel - 1]) / 2;
    for(unsigned n = 0; n < points.size(); ++n) {
      data.Restrict(axis, std::fma(scale, points[n], origin), child);
      const std::array<long double, 2> inner = IntegratePrepared(child, points, weights);
      const long double weight = weights[n] * scale * child.GetParentMeasureScale();
      for(unsigned component = 0; component < 2; ++component) result[component] += weight * inner[component];
    }
  }
  return result;
}

void TestNestedIntegration() {
  currentCheck = "nested volume and surface integration";
  const std::vector<long double> abscissa(boost::math::quadrature::gauss_kronrod<long double, 15>::abscissa().begin(), boost::math::quadrature::gauss_kronrod<long double, 15>::abscissa().end());
  const std::vector<long double> kronrodWeights(boost::math::quadrature::gauss_kronrod<long double, 15>::weights().begin(), boost::math::quadrature::gauss_kronrod<long double, 15>::weights().end());
  std::vector<long double> points, weights;
  for(unsigned i = 0; i < abscissa.size(); ++i) {
    points.push_back(abscissa[i]);
    weights.push_back(kronrodWeights[i]);
    if(abscissa[i] != 0) {
      points.push_back(-abscissa[i]);
      weights.push_back(kronrodWeights[i]);
    }
  }
  const GeomElType elements[6] = {LINE, QUAD, TRI, HEX, TET, WEDGE};
  const long double volumes[6] = {1, 1, 0.5L, 1, 1.0L / 6, 0.5L};
  const long double surfaces[6] = {1, std::sqrt(2.0L), 1 / std::sqrt(2.0L), std::sqrt(2.0L), 1 / (2 * std::sqrt(2.0L)), std::sqrt(2.0L)/2};
  for(unsigned element = 0; element < 6; ++element) {
    CutFemLevelSetData data(elements[element], MakeField(elements[element], 7));
    const std::array<long double, 2> integrals = IntegratePrepared(data, points, weights);
    CheckNear(integrals[0], volumes[element] / 2, 3e-15L, "half-domain volume " + GeometryName(elements[element]));
    CheckNear(integrals[1], surfaces[element], 3e-14L, "interface measure " + GeometryName(elements[element]));
    Check(data.GetWorkCounts()[0] == 1, "one nodal reconstruction through nested integration");
  }
}

template<class Real> Real BenchmarkValue(const std::vector<Real>& coefficient, unsigned dimension, const std::vector<Real>& point) {
  Real result = 0;
  if(dimension == 1) return (coefficient[2] * point[0] + coefficient[1]) * point[0] + coefficient[0];
  for(unsigned i = 3; i > 0; --i) {
    Real row = 0;
    for(unsigned j = 3; j > 0; --j) {
      Real line;
      if(dimension == 2) line = coefficient[3 * (i - 1) + j - 1];
      else line = (coefficient[(3 * (i - 1) + j - 1) * 3 + 2] * point[2] + coefficient[(3 * (i - 1) + j - 1) * 3 + 1]) * point[2] + coefficient[(3 * (i - 1) + j - 1) * 3];
      row = j == 3 ? line : row * point[1] + line;
    }
    result = i == 3 ? row : result * point[0] + row;
  }
  return result;
}

boost::multiprecision::cpp_bin_float_quad QuadReference(const std::vector<long double>& coefficient, unsigned dimension, const std::vector<double>& point) {
  boost::multiprecision::cpp_bin_float_quad result = 0;
  for(unsigned n = 0; n < coefficient.size(); ++n) {
    unsigned index = n;
    std::array<unsigned, 3> exponent = {0, 0, 0};
    for(unsigned axis = dimension; axis > 0; --axis) {
      exponent[axis - 1] = index % 3;
      index /= 3;
    }
    boost::multiprecision::cpp_bin_float_quad term(coefficient[n]);
    for(unsigned axis = 0; axis < dimension; ++axis) for(unsigned power = 0; power < exponent[axis]; ++power) term *= boost::multiprecision::cpp_bin_float_quad(point[axis]);
    result += term;
  }
  return result;
}

void TestPrecisionPaths(bool benchmark, std::mt19937& random) {
  currentCheck = "double/long-double polynomial arithmetic";
  const GeomElType elements[6] = {LINE, QUAD, TRI, HEX, TET, WEDGE};
  std::uniform_real_distribution<double> value(-1, 1);
  for(GeomElType geometry : elements) {
    std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    std::vector<double> nodal(basis->_nc);
    for(double& coefficient : nodal) coefficient = value(random);
    CutFemLevelSetData data(geometry, nodal);
    const std::vector<long double>& coefficient = data.GetPowerCoefficients();
    const std::vector<double> doubleCoefficient(coefficient.begin(), coefficient.end());
    std::vector<std::vector<double>> doublePoints;
    std::vector<std::vector<long double>> extendedPoints;
    long double doubleError = 0, extendedError = 0;
    for(unsigned sample = 0; sample < 64; ++sample) {
      const std::vector<long double> point = RandomPoint(geometry, random);
      doublePoints.push_back(std::vector<double>(point.begin(), point.end()));
      extendedPoints.push_back(std::vector<long double>(doublePoints.back().begin(), doublePoints.back().end()));
      const boost::multiprecision::cpp_bin_float_quad reference = QuadReference(coefficient, data.GetDimension(), doublePoints.back());
      const double doubleValue = BenchmarkValue(doubleCoefficient, data.GetDimension(), doublePoints.back());
      const long double extendedValue = BenchmarkValue(coefficient, data.GetDimension(), extendedPoints.back());
      doubleError = std::max(doubleError, static_cast<long double>(abs(boost::multiprecision::cpp_bin_float_quad(doubleValue) - reference) / (1 + abs(reference))));
      extendedError = std::max(extendedError, static_cast<long double>(abs(boost::multiprecision::cpp_bin_float_quad(extendedValue) - reference) / (1 + abs(reference))));
      CheckNear(extendedValue, static_cast<long double>(reference), 2e-16L, "extended evaluation versus multiprecision");
      CheckNear(doubleValue, static_cast<long double>(reference), 3e-13L, "double evaluation versus multiprecision");
    }
    if(!benchmark) continue;
    volatile double doubleChecksum = 0;
    volatile long double extendedChecksum = 0;
    const std::chrono::steady_clock::time_point doubleStart = std::chrono::steady_clock::now();
    for(unsigned iteration = 0; iteration < 2000; ++iteration) for(const std::vector<double>& point : doublePoints) doubleChecksum += BenchmarkValue(doubleCoefficient, data.GetDimension(), point);
    const double doubleSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - doubleStart).count();
    const std::chrono::steady_clock::time_point extendedStart = std::chrono::steady_clock::now();
    for(unsigned iteration = 0; iteration < 2000; ++iteration) for(const std::vector<long double>& point : extendedPoints) extendedChecksum += BenchmarkValue(coefficient, data.GetDimension(), point);
    const double extendedSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - extendedStart).count();
    Check(std::isfinite(doubleChecksum) && std::isfinite(extendedChecksum), "benchmark checksums");
    std::cout << std::setprecision(8) << "precision " << GeometryName(geometry) << ": double_error=" << doubleError << ", extended_error=" << extendedError << ", double_seconds=" << doubleSeconds << ", extended_seconds=" << extendedSeconds << '\n';
  }
}

void PhysicalMap(const std::vector<long double>& point, bool nonAffine, std::vector<long double>& physical, std::vector<std::vector<long double>>& jacobian) {
  const unsigned dimension = point.size();
  const long double linear[3][3] = {{1, 0.13L, 0.04L}, {0.03L, 1, 0.08L}, {-0.06L, 0.02L, 1}};
  const long double mixed[3][4] = {{0.07L, 0.03L, 0.05L, 0.04L}, {-0.04L, 0.06L, -0.03L, 0.02L}, {0.03L, -0.02L, 0.04L, -0.03L}};
  const unsigned exponents[4][3] = {{1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
  physical.assign(dimension, 0);
  jacobian.assign(dimension, std::vector<long double>(dimension, 0));
  for(unsigned component = 0; component < dimension; ++component) {
    for(unsigned axis = 0; axis < dimension; ++axis) {
      physical[component] += linear[component][axis] * point[axis];
      jacobian[component][axis] = linear[component][axis];
    }
    if(!nonAffine) continue;
    for(unsigned term = 0; term < 4; ++term) {
      if(dimension == 2 && exponents[term][2] != 0) continue;
      long double value = mixed[component][term];
      for(unsigned axis = 0; axis < dimension; ++axis) if(exponents[term][axis]) value *= point[axis];
      physical[component] += value;
      for(unsigned derivative = 0; derivative < dimension; ++derivative) if(exponents[term][derivative]) {
          value = mixed[component][term];
          for(unsigned axis = 0; axis < dimension; ++axis) if(axis != derivative && exponents[term][axis]) value *= point[axis];
          jacobian[component][derivative] += value;
        }
    }
  }
}

void TestPhysicalSphere(unsigned cases, std::mt19937& random) {
  currentCheck = "physical circle/sphere on affine and multilinear geometry";
  const long double center[3] = {0.13L, -0.17L, 0.11L};
  for(GeomElType geometry : {
        QUAD, HEX
      }) for(unsigned mode = 0; mode < 2; ++mode) {
      std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
      std::vector<double> nodal(basis->_nc);
      for(unsigned n = 0; n < nodal.size(); ++n) {
        const double* native = basis->GetXcoarse(n);
        std::vector<long double> point(Dimension(geometry)), physical;
        std::vector<std::vector<long double>> jacobian;
        for(unsigned axis = 0; axis < point.size(); ++axis) point[axis] = UnitNode(*basis,n,axis);
        PhysicalMap(point, mode != 0, physical, jacobian);
        long double value = 0.35L * 0.35L;
        for(unsigned axis = 0; axis < point.size(); ++axis) value -= (physical[axis] - center[axis]) * (physical[axis] - center[axis]);
        nodal[n] = static_cast<double>(value);
      }
      CutFemLevelSetData data(geometry, nodal), child;
      for(unsigned sample = 0; sample < cases * 4; ++sample) {
        const std::vector<long double> point = RandomPoint(geometry, random);
        std::vector<long double> physical;
        std::vector<std::vector<long double>> jacobian;
        PhysicalMap(point, mode != 0, physical, jacobian);
        long double expected = 0.35L * 0.35L;
        for(unsigned component = 0; component < point.size(); ++component) expected -= (physical[component] - center[component]) * (physical[component] - center[component]);
        CheckNear(data.Evaluate(point), expected, 8e-15L, "physical sphere pullback");
        CheckNear(data.EvaluateBernstein(point), expected, 1e-14L, "physical sphere Bernstein coefficients");
        const std::vector<long double> gradient = data.GetGradient(point);
        for(unsigned axis = 0; axis < point.size(); ++axis) {
          long double derivative = 0;
          for(unsigned component = 0; component < point.size(); ++component) derivative -= 2 * (physical[component] - center[component]) * jacobian[component][axis];
          CheckNear(gradient[axis], derivative, 2e-14L, "physical sphere pullback gradient");
        }
        data.Restrict(data.GetDimension() - 1, point.back(), child);
        const std::vector<long double> childPoint(point.begin(), point.end() - 1);
        CheckNear(child.Evaluate(childPoint), expected, 8e-15L, "physical sphere slice");
      }
    }
}

void SurfaceTestMap(const std::vector<long double>& point, unsigned mode, std::vector<long double>& physical, std::vector<std::vector<long double>>& jacobian) {
  PhysicalMap(point,mode != 0,physical,jacobian);
  if(mode != 2) return;
  long double product = 1;
  for(long double coordinate : point) product *= coordinate*coordinate;
  for(unsigned component = 0; component < point.size(); ++component) {
    const long double coefficient = 0.006L*(component+1);
    physical[component] += coefficient*product;
    for(unsigned axis = 0; axis < point.size(); ++axis) { long double derivative = 2*coefficient*point[axis]; for(unsigned other = 0; other < point.size(); ++other) if(other != axis) derivative *= point[other]*point[other]; jacobian[component][axis] += derivative; }
  }
}

void PrepareSurfaceTestPolynomials(const CutFemLevelSetData& data, unsigned mode, unsigned& degree, std::vector<std::vector<long double>>& components) {
  const unsigned dimension = data.GetDimension();
  std::unique_ptr<femus::basis> basis = MakeBasis(data.GetGeometry());
  std::vector<std::vector<double>> coordinates(dimension,std::vector<double>(basis->_nc));
  for(unsigned n = 0; n < coordinates[0].size(); ++n) { const double* native = basis->GetXcoarse(n); std::vector<long double> point(dimension), physical; std::vector<std::vector<long double>> jacobian; for(unsigned axis = 0; axis < dimension; ++axis) point[axis] = UnitNode(*basis,n,axis); SurfaceTestMap(point,mode,physical,jacobian); for(unsigned component = 0; component < dimension; ++component) coordinates[component][n] = static_cast<double>(physical[component]); }
  std::vector<std::vector<long double>> coordinatePolynomials(dimension);
  for(unsigned component = 0; component < dimension; ++component) { CutFemLevelSetData coordinateData(data.GetGeometry(),coordinates[component]); coordinatePolynomials[component] = coordinateData.GetPowerCoefficients(); }
  components.resize(dimension);
  for(unsigned axis = 0; axis < dimension; ++axis) components[axis] = data.GetDerivativeCoefficients(axis);
  if(dimension == 3) { FELevelSetMoments::SurfaceGradientCofactor3D(components,coordinatePolynomials,1); degree = FELevelSetMoments::SurfacePolynomialCubeSize(static_cast<unsigned>(components[0].size()))-1; }
  else if(dimension == 2) {
    std::vector<std::vector<std::vector<long double>>> gradient(2,std::vector<std::vector<long double>>(3,std::vector<long double>(3))), mapping = gradient;
    for(unsigned component = 0; component < 2; ++component) for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) { gradient[component][i][j] = components[component][3*i+j]; mapping[component][i][j] = coordinatePolynomials[component][3*i+j]; }
    FELevelSetMoments::SurfaceGradientCofactor2D(gradient,mapping,1,1);
    const unsigned side = static_cast<unsigned>(gradient[0].size());
    degree = side-1;
    for(unsigned component = 0; component < 2; ++component) { components[component].resize(side*side); for(unsigned i = 0; i < side; ++i) for(unsigned j = 0; j < side; ++j) components[component][side*i+j] = gradient[component][i][j]; }
  }
  else degree = 2;
}

void CheckSelectiveSurfaceFactors(const CutFemLevelSetData& line, long double root, const std::array<long double,3>& cofactorGradient, const std::array<std::array<long double,3>,3>& cofactorMatrix, const std::array<long double,4>& scalar, const std::array<long double,3>& normal, const std::array<long double,4>& matrixScalar, const std::array<long double,3>& matrixNormal) {
  CheckNear(line.GetLineDeltaFactor(root),scalar[0],0,"selective delta");
  CheckNear(line.GetLineReferenceSurfaceFactor(root),scalar[1],0,"selective reference surface");
  CheckNear(line.GetLinePhysicalSurfaceFactor(root,cofactorGradient),scalar[2],0,"selective physical surface");
  CheckNear(line.GetLineTangentialSurfaceFactor(root,cofactorGradient),scalar[3],0,"selective tangential surface");
  CheckNear(line.GetLinePhysicalSurfaceFactorFromCofactor(root,cofactorMatrix),matrixScalar[2],0,"selective matrix physical surface");
  CheckNear(line.GetLineTangentialSurfaceFactorFromCofactor(root,cofactorMatrix),matrixScalar[3],0,"selective matrix tangential surface");
  std::array<long double,3> selectedNormals, selectedMatrixNormals, alias = cofactorGradient;
  line.GetLineNormalSurfaceFactors(root,cofactorGradient,selectedNormals);
  line.GetLineNormalSurfaceFactorsFromCofactor(root,cofactorMatrix,selectedMatrixNormals);
  line.GetLineNormalSurfaceFactors(root,alias,alias);
  for(unsigned component = 0; component < 3; ++component) { CheckNear(selectedNormals[component],normal[component],0,"selected normal family"); CheckNear(selectedMatrixNormals[component],matrixNormal[component],0,"selected matrix normal family"); CheckNear(alias[component],normal[component],0,"selected normal input/output alias"); }
  for(unsigned component = 0; component < line.GetRootDimension(); ++component) { CheckNear(line.GetLineNormalSurfaceFactor(root,cofactorGradient,component),normal[component],0,"single normal component"); CheckNear(line.GetLineNormalSurfaceFactorFromCofactor(root,cofactorMatrix,component),matrixNormal[component],0,"single matrix normal component"); }
}

void TestFastQuadratureInterfaces(unsigned cases, std::mt19937& random) {
  currentCheck = "fast quadrature interfaces and shared surface factors";
  for(GeomElType geometry : {LINE,QUAD,TRI,HEX,TET,WEDGE}) {
    CutFemLevelSetData data(geometry,MakeField(geometry,1));
    for(unsigned sample = 0; sample < cases; ++sample) {
      const std::vector<long double> point = RandomPoint(geometry,random);
      std::array<long double,3> gradient, rootGradient;
      data.GetGradient(point,gradient);
      data.GetRootGradient(point,rootGradient);
      const std::vector<long double> expected = data.GetGradient(point);
      for(unsigned axis = 0; axis < 3; ++axis) { CheckNear(gradient[axis],axis < data.GetDimension() ? expected[axis] : 0,0,"allocation-free gradient"); CheckNear(rootGradient[axis],gradient[axis],0,"allocation-free root gradient"); }
    }
    const unsigned modes = geometry == QUAD || geometry == HEX ? 3U : 2U;
    for(unsigned mode = 0; mode < modes; ++mode) {
      unsigned degree = 0;
      std::vector<std::vector<long double>> components, childComponents, lineComponents;
      PrepareSurfaceTestPolynomials(data,mode,degree,components);
      CutFemLevelSetData child, line;
      std::uniform_real_distribution<double> position(0.1,0.3);
      for(unsigned sample = 0; sample < cases; ++sample) {
        if(data.GetDimension() == 3) { const long double coordinate = position(random); data.Restrict(2,coordinate,child); data.RestrictComponents(2,coordinate,degree,components,childComponents); const long double second = position(random); child.Restrict(1,second,line); child.RestrictComponents(1,second,degree,childComponents,lineComponents); }
        else if(data.GetDimension() == 2) { const long double coordinate = position(random); data.Restrict(1,coordinate,line); data.RestrictComponents(1,coordinate,degree,components,lineComponents); }
        else { line = data; lineComponents = components; }
        std::array<long double,2> roots;
        const unsigned count = line.GetLineRoots(roots);
        Check(count == 1,"regular plane slice has one root");
        Check(line.GetLineRoots() == std::vector<long double>(roots.begin(),roots.begin()+count),"fixed root output agrees with compatibility interface");
        const long double root = roots[0];
        CheckNear(line.EvaluateLine(root),line.Evaluate({root}),0,"scalar line evaluation");
        const std::vector<long double> parentPoint = line.MapToRoot({root});
        const std::vector<long double> gradient = data.GetGradient(parentPoint);
        std::vector<long double> physical;
        std::vector<std::vector<long double>> jacobian;
        SurfaceTestMap(parentPoint,mode,physical,jacobian);
        std::array<long double,3> cofactor = {0,0,0}, expected = {0,0,0}, rootGradient;
        std::array<std::array<long double,3>,3> cofactorMatrix = {};
        for(unsigned component = 0; component < data.GetDimension(); ++component) cofactor[component] = CutFemLevelSetData::EvaluatePolynomial(lineComponents[component],1,degree,{root});
        if(data.GetDimension() == 3) for(unsigned a = 0; a < 3; ++a) for(unsigned b = 0; b < 3; ++b) cofactorMatrix[a][b] = jacobian[(a+1)%3][(b+1)%3]*jacobian[(a+2)%3][(b+2)%3]-jacobian[(a+1)%3][(b+2)%3]*jacobian[(a+2)%3][(b+1)%3];
        else if(data.GetDimension() == 2) { cofactorMatrix[0][0] = jacobian[1][1]; cofactorMatrix[0][1] = -jacobian[1][0]; cofactorMatrix[1][0] = -jacobian[0][1]; cofactorMatrix[1][1] = jacobian[0][0]; }
        else cofactorMatrix[0][0] = 1;
        for(unsigned a = 0; a < data.GetDimension(); ++a) for(unsigned b = 0; b < data.GetDimension(); ++b) expected[a] += cofactorMatrix[a][b]*gradient[b];
        long double referenceNorm = 0, physicalNorm = 0;
        for(unsigned component = 0; component < data.GetDimension(); ++component) { CheckNear(cofactor[component],expected[component],2e-14L,"prepared cofactor polynomial through recursive slices"); referenceNorm = std::hypot(referenceNorm,gradient[component]); physicalNorm = std::hypot(physicalNorm,expected[component]); }
        const long double delta = line.GetLineSurfaceData(root,rootGradient);
        CheckNear(delta,1/std::abs(line.EvaluateDerivative(0,{root})),2e-18L,"shared line delta density");
        for(unsigned component = 0; component < data.GetDimension(); ++component) CheckNear(rootGradient[component],gradient[component],4e-15L,"all original surface-gradient components");
        std::array<long double,4> scalar;
        std::array<long double,3> normal;
        line.GetLineAllSurfaceFactors(root,cofactor,scalar,normal);
        CheckNear(scalar[0],delta,0,"delta option");
        CheckNear(scalar[1],referenceNorm*delta,5e-15L,"reference surface option");
        CheckNear(scalar[2],physicalNorm*delta,3e-14L,"physical surface option");
        CheckNear(scalar[3],delta/physicalNorm,3e-14L,"inverse-metric tangential option");
        for(unsigned component = 0; component < 3; ++component) CheckNear(normal[component],delta*expected[component]/physicalNorm,3e-14L,"three signed normal options");
        std::array<long double,4> matrixScalar;
        std::array<long double,3> matrixNormal;
        line.GetLineAllSurfaceFactorsFromCofactor(root,cofactorMatrix,matrixScalar,matrixNormal);
        for(unsigned option = 0; option < 4; ++option) CheckNear(matrixScalar[option],scalar[option],3e-14L,"matrix and prepared polynomial surface options agree");
        for(unsigned component = 0; component < 3; ++component) CheckNear(matrixNormal[component],normal[component],3e-14L,"matrix and prepared polynomial normal options agree");
        CheckSelectiveSurfaceFactors(line,root,cofactor,cofactorMatrix,scalar,normal,matrixScalar,matrixNormal);
        line.GetLineAllSurfaceFactors(root,cofactor,scalar,cofactor);
        Check(cofactor == normal,"surface input/output storage alias is supported");
      }
    }
  }
  std::unique_ptr<femus::basis> basis = MakeBasis(HEX);
  std::vector<double> nodal(basis->_nc);
  for(unsigned n = 0; n < nodal.size(); ++n) { const double* point = basis->GetXcoarse(n); nodal[n] = (point[2]+0.25)*point[0]+(0.25-point[2])*point[1]; }
  CutFemLevelSetData parent(HEX,nodal), reusable, snapshot;
  parent.Restrict(2,0.375L,reusable);
  Check(reusable.GetGraphDirection() == 1 && reusable.GetGraphDirection() == 1,"cached graph certificate");
  Check(reusable.GetWorkCounts()[2] == 0,"graph bounds need no derivative reconstruction");
  parent.Restrict(2,0.625L,reusable);
  Check(reusable.GetGraphDirection() == 0,"graph certificate invalidated with slice state");
  { CutFemLevelSetData temporary(HEX,nodal); temporary.Restrict(2,0.625L,snapshot); }
  CheckNear(snapshot.Evaluate({0.6L,0.3L}),0.1L,2e-18L,"slice survives destruction of its parent");
  CutFemLevelSetData tangent(LINE,{1,1,0}), regular(LINE,{-1,1,0});
  std::array<long double,3> gradient, normal;
  std::array<long double,4> scalar;
  bool rejected = false;
  try { tangent.GetLineSurfaceData(0.5L,gradient); } catch(const std::domain_error&) { rejected = true; }
  Check(rejected,"singular graph root requires endpoint handling");
  rejected = false;
  try { regular.GetLineAllSurfaceFactors(0,{0,0,0},scalar,normal); } catch(const std::domain_error&) { rejected = true; }
  Check(rejected,"zero surface metric rejected");
  rejected = false;
  try { regular.GetLineAllSurfaceFactors(0,{std::numeric_limits<long double>::infinity(),0,0},scalar,normal); } catch(const std::invalid_argument&) { rejected = true; }
  Check(rejected,"nonfinite surface components rejected");
}

void TestSelectiveSurfacePreparation(unsigned cases, std::mt19937& random) {
  currentCheck = "selective gradient and geometry preparation";
  for(GeomElType geometry : {LINE,QUAD,TRI,HEX,TET,WEDGE}) {
    CutFemLevelSetData parent(geometry,MakeField(geometry,1)), full, bare, fullLine, bareLine;
    for(unsigned sample = 0; sample < cases; ++sample) {
      if(parent.GetDimension() > 1) {
        const std::array<unsigned long long,4> before = parent.GetWorkCounts();
        parent.Restrict(parent.GetDimension()-1,0.17L,bare,false);
        Check(parent.GetWorkCounts()[2] == before[2],"omitted root derivatives are not prepared");
        parent.Restrict(parent.GetDimension()-1,0.17L,full);
        Check(!bare.HasRootGradientContext() && full.HasRootGradientContext(),"explicit slice dependencies");
        const std::vector<long double> point = RandomPoint(bare.GetGeometry(),random);
        CheckNear(bare.Evaluate(point),full.Evaluate(point),0,"identical selective slice polynomial");
        Check(bare.GetBernsteinCoefficients() == full.GetBernsteinCoefficients() && bare.GetBounds() == full.GetBounds() && bare.GetPhase() == full.GetPhase(),"identical selective slice certification");
        Check(bare.MapToRoot(point) == full.MapToRoot(point),"selective slice preserves maps");
        bool rejected = false;
        try { bare.GetRootGradient(point); } catch(const std::logic_error&) { rejected = true; }
        Check(rejected,"omitted original-gradient context is explicit");
        if(bare.GetDimension() == 2) {
          full.Restrict(1,0.19L,fullLine);
          bare.Restrict(1,0.19L,bareLine);
          Check(!bareLine.HasRootGradientContext(),"slice dependencies propagate recursively");
          const std::uint64_t revision = bareLine.GetRevision();
          rejected = false;
          try { bare.Restrict(1,0.19L,bareLine,true); } catch(const std::logic_error&) { rejected = true; }
          Check(rejected && bareLine.GetRevision() == revision,"unavailable context cannot be silently restored");
        }
        else { fullLine = full; bareLine = bare; }
      }
      else { fullLine = parent; bareLine = parent; }
      std::array<long double,2> roots;
      Check(fullLine.GetLineRoots(roots) == 1,"plane root for selective preparation");
      std::array<long double,2> selectedRoots;
      Check(bareLine.GetLineRoots(selectedRoots) == 1 && selectedRoots == roots,"identical roots without original-gradient preparation");
      const long double root = roots[0];
      const std::array<long double,3> cofactorGradient = {1,0.2L,-0.1L};
      std::array<long double,4> scalar;
      std::array<long double,3> normal, selected;
      fullLine.GetLineAllSurfaceFactors(root,cofactorGradient,scalar,normal);
      const std::array<unsigned long long,4> before = bareLine.GetWorkCounts();
      CheckNear(bareLine.GetLineDeltaFactor(root),scalar[0],0,"delta without original gradients or geometry");
      CheckNear(bareLine.GetLinePhysicalSurfaceFactor(root,cofactorGradient),scalar[2],0,"physical density without original gradients");
      CheckNear(bareLine.GetLineTangentialSurfaceFactor(root,cofactorGradient),scalar[3],0,"tangential density without original gradients");
      bareLine.GetLineNormalSurfaceFactors(root,cofactorGradient,selected);
      Check(selected == normal,"normal family without original gradients");
      for(unsigned component = 0; component < bareLine.GetRootDimension(); ++component) CheckNear(bareLine.GetLineNormalSurfaceFactor(root,cofactorGradient,component),normal[component],0,"normal component without original gradients");
      Check(bareLine.GetWorkCounts() == before,"unrequested derivative work remains absent");
      if(parent.GetDimension() > 1) {
        bool rejected = false;
        try { bareLine.GetLineReferenceSurfaceFactor(root); } catch(const std::logic_error&) { rejected = true; }
        Check(rejected,"reference mode requires retained original gradients");
        parent.Restrict(parent.GetDimension()-1,0.17L,bare,true);
        Check(bare.HasRootGradientContext(),"context can be restored from its prepared ancestor");
      }
    }
    if(parent.GetDimension() > 1) for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) { parent.Restrict(axis,0.17L,bare,false); parent.Restrict(axis,0.17L,full); const std::vector<long double> point = RandomPoint(bare.GetGeometry(),random); CheckNear(bare.Evaluate(point),full.Evaluate(point),0,"selective restriction on every axis"); Check(bare.GetGradient(point) == full.GetGradient(point) && bare.GetGraphDirection() == full.GetGraphDirection(),"local derivatives and graph certification remain available"); }
  }
  CutFemLevelSetData line(LINE,{-1,1,0}), tangent(LINE,{1,1,0}), saved;
  const std::array<unsigned long long,4> before = line.GetWorkCounts();
  CheckNear(line.GetLineDeltaFactor(0.5L),0.5L,0,"delta needs only the line derivative");
  Check(line.GetWorkCounts()[2] == before[2],"delta does not prepare even local derivative coefficients");
  long double savedDelta = 0;
  { CutFemLevelSetData temporary(HEX,MakeField(HEX,1)), child; temporary.Restrict(2,0.17L,child,false); child.Restrict(1,0.19L,saved); savedDelta = saved.GetLineDeltaFactor(0.2L); }
  CheckNear(saved.GetLineDeltaFactor(0.2L),savedDelta,0,"selective snapshot survives parent destruction");
  bool rejected = false;
  try { tangent.GetLineDeltaFactor(0.5L); } catch(const std::domain_error&) { rejected = true; }
  Check(rejected,"selective delta still rejects a singular root");
  rejected = false;
  try { line.GetLinePhysicalSurfaceFactor(0,{0,0,0}); } catch(const std::domain_error&) { rejected = true; }
  Check(rejected,"selective physical metric validated");
  rejected = false;
  try { line.GetLineNormalSurfaceFactor(0,{1,0,0},1); } catch(const std::out_of_range&) { rejected = true; }
  Check(rejected,"selected normal component validated");
  rejected = false;
  try { line.GetLineTangentialSurfaceFactor(0,{std::numeric_limits<long double>::quiet_NaN(),0,0}); } catch(const std::invalid_argument&) { rejected = true; }
  Check(rejected,"selected metric rejects nonfinite components");
}

void TestSelectiveSurfaceScaling() {
  currentCheck = "selective surface factors across scales and small cuts";
  const std::array<std::array<long double,3>,3> identity = {{{1,0,0},{0,1,0},{0,0,1}}};
  for(double scale : {1e-200,1.0,1e200}) for(double radius : {0.2,1e-6,1e-12}) {
    const std::vector<double> nodal = {scale*(1-radius*radius),scale*(1-radius*radius),-scale*radius*radius};
    CutFemLevelSetData line(LINE,nodal);
    std::array<long double,2> roots;
    Check(line.GetLineRoots(roots) == 2,"scaled small cut retains both roots");
    const boost::multiprecision::cpp_bin_float_quad a = (boost::multiprecision::cpp_bin_float_quad(nodal[0])+boost::multiprecision::cpp_bin_float_quad(nodal[1]))/2-boost::multiprecision::cpp_bin_float_quad(nodal[2]);
    const boost::multiprecision::cpp_bin_float_quad exactRoot = sqrt(-boost::multiprecision::cpp_bin_float_quad(nodal[2])/a);
    const long double exactDelta = static_cast<long double>(1/(4*a*exactRoot));
    for(long double root : roots) {
      std::array<long double,3> gradient, normal, matrixNormal;
      std::array<long double,4> scalar, matrixScalar;
      line.GetLineSurfaceData(root,gradient);
      line.GetLineAllSurfaceFactors(root,gradient,scalar,normal);
      line.GetLineAllSurfaceFactorsFromCofactor(root,identity,matrixScalar,matrixNormal);
      CheckSelectiveSurfaceFactors(line,root,gradient,identity,scalar,normal,matrixScalar,matrixNormal);
      CheckNear(line.GetLineDeltaFactor(root)/exactDelta,1,5e-18L,"scaled delta versus multiprecision");
      CheckNear(line.GetLineReferenceSurfaceFactor(root),1,5e-18L,"scaled reference surface is invariant");
      CheckNear(line.GetLinePhysicalSurfaceFactor(root,gradient),1,5e-18L,"scaled physical surface is invariant");
      CheckNear(line.GetLineTangentialSurfaceFactor(root,gradient)/(exactDelta*exactDelta),1,5e-18L,"scaled inverse metric versus multiprecision");
    }
  }
}

long double QuadratureKernelSample(CutFemLevelSetData& parent, unsigned operation, unsigned sample, const std::vector<std::vector<long double>>& points, const std::vector<long double>& positions, CutFemLevelSetData& child, CutFemLevelSetData& line, const std::vector<std::vector<long double>>& components, std::vector<std::vector<long double>>& output) {
  if(operation == 0) return parent.Evaluate(points[sample]);
  if(operation == 1) return parent.GetRootGradient(points[sample])[0];
  if(operation == 2) return parent.EvaluateBernstein(points[sample]);
  if(operation == 7) { const std::vector<long double> roots = parent.GetLineRoots(); return roots.empty() ? 0 : roots[0]; }
  if(operation == 5) { parent.RestrictComponents(parent.GetDimension()-1,positions[sample],4,components,output); return output[0][0]; }
  parent.Restrict(parent.GetDimension()-1,positions[sample],child);
  if(operation == 6) return child.GetGraphDirection();
  if(operation == 3 || child.GetDimension() == 1) return child.GetPowerCoefficients()[0];
  child.Restrict(child.GetDimension()-1,child.GetGeometry() == TRI ? 0.37L : 0.23L,line);
  return line.GetPowerCoefficients()[0];
}

void BenchmarkQuadratureKernels() {
  std::mt19937 random(470014);
  std::uniform_real_distribution<double> value(-0.25,0.25);
  const GeomElType elements[6] = {LINE,QUAD,TRI,HEX,TET,WEDGE};
  const std::string names[8] = {"value","root_gradient","bernstein","restrict","nested_restrict","components_degree4","restrict_graph","line_roots"};
  const unsigned iterations = 20000;
  for(GeomElType geometry : elements) {
    CutFemLevelSetData parent(geometry,MakeField(geometry,2,1,0.29L)), child, line;
    std::vector<std::vector<long double>> points, output;
    std::vector<long double> positions;
    for(unsigned sample = 0; sample < 64; ++sample) { points.push_back(RandomPoint(geometry,random)); positions.push_back(points.back().back()); }
    unsigned count = 1;
    for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) count *= 5;
    std::vector<std::vector<long double>> components(3,std::vector<long double>(count));
    for(std::vector<long double>& component : components) for(long double& coefficient : component) coefficient = value(random);
    for(unsigned operation = 0; operation < 8; ++operation) {
      if(parent.GetDimension() == 1 && operation >= 3 && operation != 7) continue;
      if(parent.GetDimension() != 1 && operation == 7) continue;
      if(parent.GetDimension() < 3 && operation == 4) continue;
      volatile long double checksum = 0;
      for(unsigned sample = 0; sample < 64; ++sample) checksum += QuadratureKernelSample(parent,operation,sample,points,positions,child,line,components,output);
#ifdef EX47_TRACK_ALLOCATIONS
      const unsigned long long initialAllocations = heapAllocations;
#endif
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      for(unsigned iteration = 0; iteration < iterations; ++iteration) checksum += QuadratureKernelSample(parent,operation,iteration%64,points,positions,child,line,components,output);
      const double nanoseconds = 1e9*std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/iterations;
#ifdef EX47_TRACK_ALLOCATIONS
      const unsigned long long allocated = heapAllocations-initialAllocations;
      if(operation != 1 && operation != 7) Check(allocated == 0,"prepared quadrature kernel performs no heap allocations");
#endif
      Check(std::isfinite(checksum),"quadrature benchmark checksum");
      std::cout << std::setprecision(10) << "kernel " << GeometryName(geometry) << ' ' << names[operation] << " ns=" << nanoseconds << " checksum=" << checksum;
#ifdef EX47_TRACK_ALLOCATIONS
      std::cout << " allocations=" << allocated;
#endif
      std::cout << '\n';
    }
    if(parent.GetDimension() == 3) { parent.Restrict(2,positions[0],child); child.Restrict(1,child.GetGeometry() == TRI ? 0.17L : 0.13L,line); }
    else if(parent.GetDimension() == 2) parent.Restrict(1,positions[0],line);
    else line = parent;
    const std::string fastNames[4] = {"root_gradient_array","line_value_scalar","line_roots_array","line_surface_all"};
    std::array<long double,3> gradient, normal;
    std::array<long double,4> scalar;
    std::array<long double,2> roots;
    const unsigned rootCount = line.GetLineRoots(roots);
    const long double coordinate = rootCount > 0 ? roots[0] : 0.91L;
    line.GetLineAllSurfaceFactors(coordinate,{1,0.2L,0.3L},scalar,normal);
    for(unsigned operation = 0; operation < 4; ++operation) {
      volatile long double checksum = 0;
#ifdef EX47_TRACK_ALLOCATIONS
      const unsigned long long initialAllocations = heapAllocations;
#endif
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      for(unsigned iteration = 0; iteration < iterations; ++iteration) {
        if(operation == 0) { parent.GetRootGradient(points[iteration%64],gradient); checksum += gradient[0]; }
        else if(operation == 1) checksum += line.EvaluateLine(coordinate);
        else if(operation == 2) checksum += line.GetLineRoots(roots);
        else { line.GetLineAllSurfaceFactors(coordinate,{1,0.2L,0.3L},scalar,normal); checksum += scalar[2]+normal[0]; }
      }
      const double nanoseconds = 1e9*std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/iterations;
#ifdef EX47_TRACK_ALLOCATIONS
      const unsigned long long allocated = heapAllocations-initialAllocations;
      Check(allocated == 0,"fixed-output quadrature kernel performs no heap allocations");
#endif
      Check(std::isfinite(checksum),"fast quadrature benchmark checksum");
      std::cout << std::setprecision(10) << "kernel " << GeometryName(geometry) << ' ' << fastNames[operation] << " ns=" << nanoseconds << " checksum=" << checksum;
#ifdef EX47_TRACK_ALLOCATIONS
      std::cout << " allocations=" << allocated;
#endif
      std::cout << '\n';
    }
  }
}

long double RunSelectiveSurfaceKernel(unsigned operation, unsigned iterations, const CutFemLevelSetData& full, const CutFemLevelSetData& bare, const std::array<long double,64>& points, const std::array<std::array<long double,3>,64>& gradients, const std::array<std::array<long double,3>,3>& matrix) {
  volatile long double checksum = 0;
  std::array<long double,4> scalar;
  std::array<long double,3> normal;
  switch(operation) {
    case 0: for(unsigned i = 0; i < iterations; ++i) { full.GetLineAllSurfaceFactors(points[i%64],gradients[i%64],scalar,normal); checksum += scalar[0]; } break;
    case 1: for(unsigned i = 0; i < iterations; ++i) checksum += bare.GetLineDeltaFactor(points[i%64]); break;
    case 2: for(unsigned i = 0; i < iterations; ++i) { full.GetLineAllSurfaceFactors(points[i%64],gradients[i%64],scalar,normal); checksum += scalar[1]; } break;
    case 3: for(unsigned i = 0; i < iterations; ++i) checksum += full.GetLineReferenceSurfaceFactor(points[i%64]); break;
    case 4: for(unsigned i = 0; i < iterations; ++i) { full.GetLineAllSurfaceFactors(points[i%64],gradients[i%64],scalar,normal); checksum += scalar[2]; } break;
    case 5: for(unsigned i = 0; i < iterations; ++i) checksum += bare.GetLinePhysicalSurfaceFactor(points[i%64],gradients[i%64]); break;
    case 6: for(unsigned i = 0; i < iterations; ++i) { full.GetLineAllSurfaceFactors(points[i%64],gradients[i%64],scalar,normal); checksum += scalar[3]; } break;
    case 7: for(unsigned i = 0; i < iterations; ++i) checksum += bare.GetLineTangentialSurfaceFactor(points[i%64],gradients[i%64]); break;
    case 8: for(unsigned i = 0; i < iterations; ++i) { full.GetLineAllSurfaceFactors(points[i%64],gradients[i%64],scalar,normal); checksum += normal[0]; } break;
    case 9: for(unsigned i = 0; i < iterations; ++i) checksum += bare.GetLineNormalSurfaceFactor(points[i%64],gradients[i%64],0); break;
    case 10: for(unsigned i = 0; i < iterations; ++i) { full.GetLineAllSurfaceFactors(points[i%64],gradients[i%64],scalar,normal); checksum += normal[0]+normal[1]+normal[2]; } break;
    case 11: for(unsigned i = 0; i < iterations; ++i) { bare.GetLineNormalSurfaceFactors(points[i%64],gradients[i%64],normal); checksum += normal[0]+normal[1]+normal[2]; } break;
    case 12: for(unsigned i = 0; i < iterations; ++i) { full.GetLineAllSurfaceFactorsFromCofactor(points[i%64],matrix,scalar,normal); checksum += scalar[2]; } break;
    case 13: for(unsigned i = 0; i < iterations; ++i) checksum += full.GetLinePhysicalSurfaceFactorFromCofactor(points[i%64],matrix); break;
    case 14: for(unsigned i = 0; i < iterations; ++i) { full.GetLineAllSurfaceFactorsFromCofactor(points[i%64],matrix,scalar,normal); checksum += scalar[3]; } break;
    case 15: for(unsigned i = 0; i < iterations; ++i) checksum += full.GetLineTangentialSurfaceFactorFromCofactor(points[i%64],matrix); break;
    case 16: for(unsigned i = 0; i < iterations; ++i) { full.GetLineAllSurfaceFactorsFromCofactor(points[i%64],matrix,scalar,normal); checksum += normal[0]; } break;
    case 17: for(unsigned i = 0; i < iterations; ++i) checksum += full.GetLineNormalSurfaceFactorFromCofactor(points[i%64],matrix,0); break;
    case 18: for(unsigned i = 0; i < iterations; ++i) { full.GetLineAllSurfaceFactorsFromCofactor(points[i%64],matrix,scalar,normal); checksum += normal[0]+normal[1]+normal[2]; } break;
    case 19: for(unsigned i = 0; i < iterations; ++i) { full.GetLineNormalSurfaceFactorsFromCofactor(points[i%64],matrix,normal); checksum += normal[0]+normal[1]+normal[2]; } break;
    default: throw std::logic_error("Invalid surface benchmark operation.");
  }
  return checksum;
}

long double RunSelectiveRestrictionKernel(unsigned iterations, const CutFemLevelSetData& parent, bool retainRootGradient, CutFemLevelSetData& child, CutFemLevelSetData& line) {
  volatile long double checksum = 0;
  for(unsigned i = 0; i < iterations; ++i) { parent.Restrict(parent.GetDimension()-1,0.1L+0.002L*(i%64),child,retainRootGradient); if(child.GetDimension() == 2) { child.Restrict(1,0.17L,line); checksum += line.GetPowerCoefficients()[0]; } else checksum += child.GetPowerCoefficients()[0]; }
  return checksum;
}

void BenchmarkSelectiveSurfaceKernels() {
  const std::string names[22] = {"all_delta","delta","all_reference","reference","all_physical","physical","all_tangential","tangential","all_normal0","normal0","all_normals","normals","all_matrix_physical","matrix_physical","all_matrix_tangential","matrix_tangential","all_matrix_normal0","matrix_normal0","all_matrix_normals","matrix_normals","restrict_full","restrict_selected"};
  const unsigned iterations = 20000;
  for(GeomElType geometry : {LINE,QUAD,TRI,HEX,TET,WEDGE}) {
    CutFemLevelSetData parent(geometry,MakeField(geometry,2,1,0.29L)), fullChild, bareChild, full, bare;
    if(parent.GetDimension() == 3) { parent.Restrict(2,0.19L,fullChild); parent.Restrict(2,0.19L,bareChild,false); const long double y = geometry == TET ? 0.17L/0.81L : (geometry == WEDGE ? 0.17L : 0.23L); fullChild.Restrict(1,y,full); bareChild.Restrict(1,y,bare); }
    else if(parent.GetDimension() == 2) { parent.Restrict(1,geometry == TRI ? 0.17L : 0.23L,full); parent.Restrict(1,geometry == TRI ? 0.17L : 0.23L,bare,false); }
    else { full = parent; bare = parent; }
    std::array<long double,2> roots;
    Check(full.GetLineRoots(roots) > 0,"surface benchmark has an actual intersection");
    const std::array<std::array<long double,3>,3> matrix = {{{1,0.1L,0.03L},{0.05L,1,-0.02L},{0.1L,0.04L,1}}};
    std::array<long double,64> points;
    std::array<std::array<long double,3>,64> gradients = {};
    for(unsigned sample = 0; sample < points.size(); ++sample) { points[sample] = roots[0]+0.02L*(static_cast<long double>(sample)-31.5L)/64; std::array<long double,3> gradient; full.GetRootGradient({points[sample]},gradient); for(unsigned a = 0; a < parent.GetDimension(); ++a) for(unsigned b = 0; b < parent.GetDimension(); ++b) gradients[sample][a] += matrix[a][b]*gradient[b]; }
    CutFemLevelSetData restrictionChild, restrictionLine;
    long double preceding = 0;
    for(unsigned operation = 0; operation < (parent.GetDimension() == 1 ? 20U : 22U); ++operation) {
      if(operation < 20) RunSelectiveSurfaceKernel(operation,128,full,bare,points,gradients,matrix);
      else RunSelectiveRestrictionKernel(128,parent,operation == 20,restrictionChild,restrictionLine);
#ifdef EX47_TRACK_ALLOCATIONS
      const unsigned long long before = heapAllocations;
#endif
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      const long double checksum = operation < 20 ? RunSelectiveSurfaceKernel(operation,iterations,full,bare,points,gradients,matrix) : RunSelectiveRestrictionKernel(iterations,parent,operation == 20,restrictionChild,restrictionLine);
      const double nanoseconds = 1e9*std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/iterations;
#ifdef EX47_TRACK_ALLOCATIONS
      const unsigned long long allocated = heapAllocations-before;
      Check(allocated == 0,"selective surface kernels allocate no heap memory");
#endif
      if(operation%2 == 0) preceding = checksum;
      else CheckNear(checksum,preceding,0,"selective and all-measure benchmark checksums agree");
      std::cout << std::setprecision(10) << "surface " << GeometryName(geometry) << ' ' << names[operation] << " ns=" << nanoseconds << " checksum=" << checksum;
#ifdef EX47_TRACK_ALLOCATIONS
      std::cout << " allocations=" << allocated;
#endif
      std::cout << '\n';
    }
  }
}

long double maximumRootDensityError = 0;
long double maximumFoldDensityError = 0;
long double maximumRoundedRootError = 0;

void TestPreparedLineAnalysis() {
  currentCheck = "prepared line topology and stable root densities";
  CutFemLevelSetData positive(LINE,{1,1,1}), negative(LINE,{-1,-1,-1}), zero(LINE,{0,0,0});
  Check(positive.GetLineAnalysis().positiveIntervalCount == 1 && positive.GetLineAnalysis().positiveIntervals[0] == std::array<long double,2>{{0,1}},"positive constant interval");
  Check(negative.GetLineAnalysis().positiveIntervalCount == 0 && negative.GetLineAnalysis().rootCount == 0,"negative constant interval");
  Check(zero.GetLineAnalysis().identicallyZero && zero.IsIdenticallyZero(),"zero fiber is explicit metadata");
  CutFemLevelSetData quadratic(LINE,{0.75,0.75,-0.25});
  const CutFemLevelSetData::LineAnalysis& analysis = quadratic.GetLineAnalysis();
  Check(analysis.degree == 2 && analysis.rootCount == 2 && analysis.positiveIntervalCount == 2,"quadratic topology");
  Check(analysis.positiveIntervals[0] == std::array<long double,2>{{0,0.25L}} && analysis.positiveIntervals[1] == std::array<long double,2>{{0.75L,1}},"positive quadratic branches");
  Check(analysis.roots[0].derivative == -2 && analysis.roots[1].derivative == 2,"stable signed root slopes");
  Check(&analysis == &quadratic.GetLineAnalysis(),"line analysis is retained");
  CutFemLevelSetData repeated(LINE,{1,1,0}), endpoints(LINE,{0,0,1}), linear(LINE,{-1.25,0.75,-0.25});
  Check(repeated.GetLineAnalysis().rootCount == 1 && repeated.GetLineAnalysis().roots[0].multiplicity == 2,"double-root multiplicity");
  Check(endpoints.GetLineAnalysis().roots[0].boundary == -1 && endpoints.GetLineAnalysis().roots[1].boundary == 1,"endpoint root tags");
  Check(linear.GetLineAnalysis().degree == 1 && linear.GetLineAnalysis().positiveIntervals[0] == std::array<long double,2>{{0.625L,1}},"linear topology");
  for(unsigned centerIndex = 1; centerIndex <= 64; ++centerIndex) for(unsigned exponent = 6; exponent <= 10; ++exponent) {
    const long double center = 0.123456789L+0.01L*centerIndex, radius = std::pow(10.0L,-static_cast<int>(exponent));
    CutFemLevelSetData line(LINE,{static_cast<double>((-1-center)*(-1-center)-radius*radius),static_cast<double>((1-center)*(1-center)-radius*radius),static_cast<double>(center*center-radius*radius)});
    const std::vector<long double>& power = line.GetPowerCoefficients();
    const boost::multiprecision::cpp_bin_float_100 a(power[2]), b(power[1]), c(power[0]);
    const boost::multiprecision::cpp_bin_float_100 discriminant = b*b-4*a*c;
    const CutFemLevelSetData::LineAnalysis& roots = line.GetLineAnalysis();
    if(discriminant <= 0) { Check(roots.rootCount == (discriminant == 0 ? 1U : 0U),"near-double root count versus multiprecision"); continue; }
    Check(roots.rootCount == 2,"translated near-double roots survive");
    const long double expected = static_cast<long double>(1/sqrt(discriminant));
    for(unsigned n = 0; n < roots.rootCount; ++n) {
      const long double actual = line.GetLineDeltaFactor(roots.roots[n]);
      const long double error = std::abs(actual/expected-1);
      maximumRootDensityError = std::max(maximumRootDensityError,error);
      maximumRoundedRootError = std::max(maximumRoundedRootError,std::abs(1/(std::abs(2*power[2]*roots.roots[n].coordinate+power[1])*expected)-1));
      Check(error < 3e-18L,"root density versus 100-digit retained-polynomial reference");
      CheckNear(line.GetLineDeltaFactor(roots.roots[n].coordinate),actual,0,"coordinate compatibility reuses the stable root slope");
      CheckNear(line.GetLineReferenceSurfaceFactor(roots.roots[n]),1,3e-18L,"translated LINE reference density preserves cancellation between gradient and slope");
      const std::array<std::array<long double,3>,3> identity = {{{1,0,0},{0,1,0},{0,0,1}}};
      CheckNear(line.GetLinePhysicalSurfaceFactorFromCofactor(roots.roots[n],identity),1,3e-18L,"translated LINE physical density uses consistent stable root gradients");
    }
  }
  quadratic.Reset(LINE,{-1,-1,-1});
  Check(quadratic.GetLineAnalysis().rootCount == 0 && quadratic.GetLineAnalysis().positiveIntervalCount == 0,"reset invalidates line topology");
}

std::vector<boost::multiprecision::cpp_bin_float_100> ReferenceNodalPower(GeomElType geometry, const std::vector<double>& nodal) {
  const unsigned size = Dimension(geometry) == 1 ? 3U : (Dimension(geometry) == 2 ? 9U : 27U);
  std::vector<boost::multiprecision::cpp_bin_float_100> power(size,0);
  std::vector<double> basis(nodal.size(),0);
  for(unsigned node = 0; node < nodal.size(); ++node) {
    basis[node] = 1;
    CutFemLevelSetData unitBasis(geometry,basis);
    const std::vector<long double>& coefficient = unitBasis.GetPowerCoefficients();
    for(unsigned term = 0; term < size; ++term) power[term] += boost::multiprecision::cpp_bin_float_100(nodal[node])*boost::multiprecision::cpp_bin_float_100(coefficient[term]);
    basis[node] = 0;
  }
  return power;
}

std::array<boost::multiprecision::cpp_bin_float_100,3> ReferenceFoldPolynomial(const std::vector<boost::multiprecision::cpp_bin_float_100>& power, GeomElType geometry, long double height, const boost::multiprecision::cpp_bin_float_100& coordinate) {
  std::array<boost::multiprecision::cpp_bin_float_100,3> values;
  const std::vector<boost::multiprecision::cpp_bin_float_100> coefficient(power.begin(),power.end());
  for(unsigned sample = 0; sample < 3; ++sample) {
    const boost::multiprecision::cpp_bin_float_100 t = boost::multiprecision::cpp_bin_float_100(sample)/2;
    std::vector<boost::multiprecision::cpp_bin_float_100> point(Dimension(geometry));
    point[0] = coordinate;
    point[1] = t;
    if(geometry == TET) { point[0] *= 1-boost::multiprecision::cpp_bin_float_100(height); point[1] = (1-boost::multiprecision::cpp_bin_float_100(height))*(1-coordinate)*t; }
    else if(geometry == WEDGE) point[1] = (1-coordinate)*t;
    if(point.size() == 3) point[2] = height;
    values[sample] = BenchmarkValue(coefficient,point.size(),point);
  }
  return {values[0],-3*values[0]+4*values[1]-values[2],2*values[0]-4*values[1]+2*values[2]};
}

void TestCompensatedFoldRestrictions() {
  currentCheck = "last-ULP fold restrictions from retained parent polynomials";
  for(GeomElType geometry : {QUAD,HEX,TET,WEDGE}) {
    const std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    std::vector<double> nodal(basis->_nc);
    for(unsigned n = 0; n < nodal.size(); ++n) {
      const double* point = basis->GetXcoarse(n);
      const long double x = point[0]-(geometry == QUAD ? 0.17L : 0.35L), y = point[1]-(geometry == QUAD ? 0.23L : 0.22L), z = geometry == QUAD ? 0 : point[2]-0.16L;
      nodal[n] = static_cast<double>(x*x+y*y+z*z-(geometry == QUAD ? 0.29L*0.29L : 0.12L*0.12L));
    }
    CutFemLevelSetData parent(geometry,nodal), planar, line;
    const std::vector<boost::multiprecision::cpp_bin_float_100> power = ReferenceNodalPower(geometry,nodal);
    const long double height = geometry == QUAD ? 0 : (geometry == TET ? 0.16L : 0.58L);
    if(geometry != QUAD) parent.Restrict(2,height,planar,false);
    boost::multiprecision::cpp_bin_float_100 lower = 0, upper = geometry == QUAD ? 0.585L : (geometry == HEX ? 0.675L : (geometry == TET ? 0.35L/(1-height) : 0.35L));
    for(unsigned iteration = 0; iteration < 220; ++iteration) {
      const boost::multiprecision::cpp_bin_float_100 middle = (lower+upper)/2;
      const std::array<boost::multiprecision::cpp_bin_float_100,3> coefficient = ReferenceFoldPolynomial(power,geometry,height,middle);
      if(coefficient[1]*coefficient[1]-4*coefficient[2]*coefficient[0] > 0) upper = middle;
      else lower = middle;
    }
    long double coordinate = static_cast<long double>((lower+upper)/2);
    for(unsigned n = 0; n < 5; ++n) coordinate = std::nextafter(coordinate,-std::numeric_limits<long double>::infinity());
    unsigned positive = 0;
    for(unsigned sample = 0; sample < 24; ++sample) {
      if(geometry == QUAD) parent.Restrict(0,coordinate,line,false);
      else planar.Restrict(0,coordinate,line,false);
      const std::array<boost::multiprecision::cpp_bin_float_100,3> coefficient = ReferenceFoldPolynomial(power,geometry,height,boost::multiprecision::cpp_bin_float_100(coordinate));
      const boost::multiprecision::cpp_bin_float_100 discriminant = coefficient[1]*coefficient[1]-4*coefficient[2]*coefficient[0];
      const CutFemLevelSetData::LineAnalysis& analysis = line.GetLineAnalysis();
      Check(analysis.compensated,"near-fold slice takes the compensated path");
      Check(analysis.rootCount == (discriminant > 0 ? 2U : 0U),"near-fold topology versus the original parent polynomial");
      if(geometry != QUAD) {
        const boost::multiprecision::cpp_bin_float_100 scale = planar.GetGeometry() == TRI ? (1-boost::multiprecision::cpp_bin_float_100(coordinate)) : boost::multiprecision::cpp_bin_float_100(1), normalization(planar.GetNormalizationScale());
        const long double expected = static_cast<long double>(discriminant/(scale*scale*normalization*normalization));
        Check(std::abs(planar.EvaluateFiberDiscriminant(1,coordinate)/expected-1) < 4e-16L,"3D slice discriminant refinement retains the original parent information");
      }
      if(discriminant > 0) {
        ++positive;
        const long double expected = static_cast<long double>(1/sqrt(discriminant));
        for(unsigned n = 0; n < analysis.rootCount; ++n) {
          const long double error = std::abs(line.GetLineDeltaFactor(analysis.roots[n])/expected-1);
          maximumFoldDensityError = std::max(maximumFoldDensityError,error);
          Check(error < 2e-16L,"fold slope retains compensated parent contraction accuracy");
        }
      }
      coordinate = std::nextafter(coordinate,std::numeric_limits<long double>::infinity());
    }
    Check(positive > 0 && positive < 24,"fold sweep crosses both phases");
    const CutFemLevelSetData saved = line;
    const CutFemLevelSetData::LineAnalysis savedAnalysis = saved.GetLineAnalysis();
    parent.Reset(geometry,std::vector<double>(nodal.size(),1));
    Check(saved.GetLineAnalysis().rootCount == savedAnalysis.rootCount,"compensated slice owns its source after parent reset");
  }
  const std::unique_ptr<femus::basis> basis = MakeBasis(QUAD);
  std::vector<double> nodal(9);
  for(unsigned n = 0; n < nodal.size(); ++n) { const double* point = basis->GetXcoarse(n); nodal[n] = (point[0]-0.125)*(point[0]-0.125)+point[1]*point[1]-0.375*0.375; }
  CutFemLevelSetData centered(QUAD,nodal), line;
  for(long double coordinate : {0.75L-1e-4L,0.75L-1e-8L,0.75L-1e-12L,0.75L-1e-16L,std::nextafter(0.75L,0.0L)}) {
    centered.Restrict(0,coordinate,line,false);
    const std::array<boost::multiprecision::cpp_bin_float_100,3> coefficient = ReferenceFoldPolynomial(ReferenceNodalPower(QUAD,nodal),QUAD,0,boost::multiprecision::cpp_bin_float_100(coordinate));
    const long double expected = static_cast<long double>(1/sqrt(coefficient[1]*coefficient[1]-4*coefficient[2]*coefficient[0]));
    const CutFemLevelSetData::LineAnalysis& analysis = line.GetLineAnalysis();
    Check(analysis.compensated && analysis.rootCount == 2,"centered fold recovers cancellation in the contracted constant coefficient");
    for(unsigned n = 0; n < 2; ++n) Check(std::abs(line.GetLineDeltaFactor(analysis.roots[n])/expected-1) < 3e-18L,"centered-fold density versus retained-parent reference");
  }
  for(unsigned n = 0; n < nodal.size(); ++n) { const double* point = basis->GetXcoarse(n); nodal[n] = (point[0]-0.125)*(point[0]-0.375)*point[1]*point[1]-1e-200; }
  CutFemLevelSetData small(QUAD,nodal);
  const std::vector<boost::multiprecision::cpp_bin_float_100> power = ReferenceNodalPower(QUAD,nodal);
  long double coordinate = 0.5625L;
  for(unsigned sample = 0; sample < 8; ++sample) {
    coordinate = std::nextafter(coordinate,0.0L);
    small.Restrict(0,coordinate,line,false);
    const std::array<boost::multiprecision::cpp_bin_float_100,3> coefficient = ReferenceFoldPolynomial(power,QUAD,0,boost::multiprecision::cpp_bin_float_100(coordinate));
    const CutFemLevelSetData::LineAnalysis& analysis = line.GetLineAnalysis();
    Check(analysis.compensated && analysis.rootCount == 2,"a rounded-away quadratic coefficient cannot hide a tiny cut");
    const boost::multiprecision::cpp_bin_float_100 x(coordinate), native = 2*x-1, exactA = 4*(native-boost::multiprecision::cpp_bin_float_100(0.125L))*(native-boost::multiprecision::cpp_bin_float_100(0.375L));
    const long double expected = static_cast<long double>(1/sqrt(-4*exactA*boost::multiprecision::cpp_bin_float_100(nodal[8])));
    for(unsigned n = 0; n < 2; ++n) Check(std::abs(line.GetLineDeltaFactor(analysis.roots[n])/expected-1) < 3e-18L,"tiny-cut slope versus retained-parent reference");
  }
}

long double BoundaryDistance(GeomElType geometry, unsigned boundary, const std::vector<long double>& point) {
  if(geometry == TRI || geometry == TET || (geometry == WEDGE && boundary < 3)) return boundary == 0 ? 1-point[0]-point[1]-(geometry == TET ? point[2] : 0) : point[boundary-1];
  const unsigned axis = geometry == WEDGE ? 2U : boundary/2, side = geometry == WEDGE ? boundary-3 : boundary%2;
  return point[axis]-(side == 0 ? 0 : 1);
}

long double EmbeddingMeasure(const std::vector<std::vector<long double>>& map) {
  if(map[0].size() == 1) { long double square = 0; for(const std::vector<long double>& row : map) square += row[0]*row[0]; return std::sqrt(square); }
  std::array<long double,3> normal;
  for(unsigned axis = 0; axis < 3; ++axis) normal[axis] = map[(axis+1)%3][0]*map[(axis+2)%3][1]-map[(axis+2)%3][0]*map[(axis+1)%3][1];
  return std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
}

void TestBoundaryCertificates(std::mt19937& random) {
  currentCheck = "exact zero traces on every boundary";
  for(GeomElType geometry : {LINE,QUAD,TRI,HEX,TET,WEDGE}) {
    const std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    CutFemLevelSetData count(geometry,std::vector<double>(basis->_nc,1));
    for(unsigned boundary = 0; boundary < count.GetBoundaryCount(); ++boundary) for(double scale : {1e-200,1.0,1e200}) {
      std::vector<double> nodal(basis->_nc);
      for(unsigned n = 0; n < nodal.size(); ++n) {
        const double* coordinate = basis->GetXcoarse(n);
        std::vector<long double> point(Dimension(geometry));
        for(unsigned axis = 0; axis < point.size(); ++axis) point[axis] = UnitNode(*basis,n,axis);
        nodal[n] = std::abs(BoundaryDistance(geometry,boundary,point)) < 1e-14L ? 0 : scale*(0.5+0.125*n);
      }
      CutFemLevelSetData data(geometry,nodal), face;
      Check(!data.IsIdenticallyZero() && data.IsZeroOnBoundary(boundary),"nodal trace zeros are retained exactly");
      Check(data.GetZeroBoundaryMask() == (1U<<boundary),"only the prescribed boundary is zero");
      if(geometry == LINE) { const CutFemLevelSetData::LineAnalysis& analysis = data.GetLineAnalysis(); bool found = false; for(unsigned n = 0; n < analysis.rootCount; ++n) found = found || analysis.roots[n].boundary == (boundary == 0 ? -1 : 1); Check(found,"zero endpoint is tagged in line analysis"); continue; }
      data.RestrictBoundary(boundary,face);
      Check(face.IsIdenticallyZero() && face.GetZeroBoundaryMask() == (1U<<face.GetBoundaryCount())-1,"boundary restriction preserves the exact zero trace");
      if(face.GetDimension() == 2) CheckNear(face.EvaluateFiberDiscriminant(0,0.0L),0,0,"zero boundary traces keep zero discriminants");
      for(unsigned sample = 0; sample < 8; ++sample) {
        const std::vector<long double> point = RandomPoint(face.GetGeometry(),random), parentPoint = face.MapToParent(point);
        CheckNear(BoundaryDistance(geometry,boundary,parentPoint),0,3e-18L,"boundary chart lies on the requested face");
        CheckNear(face.Evaluate(point),0,0,"certified zero trace has zero scalar coefficients");
        const std::vector<long double> gradient = face.GetRootGradient(point), expected = data.GetGradient(face.MapToRoot(point));
        for(unsigned component = 0; component < gradient.size(); ++component) CheckNear(gradient[component]/scale,expected[component]/scale,3e-14L,"boundary retains every original gradient component");
      }
      CheckNear(face.GetParentMeasureScale(),EmbeddingMeasure(face.GetParentMap()),2e-18L,"boundary embedding measure includes sloping faces");
    }
    CutFemLevelSetData tiny(geometry,std::vector<double>(basis->_nc,1e-300));
    Check(!tiny.IsIdenticallyZero() && tiny.GetZeroBoundaryMask() == 0,"small nonzero values are never treated as exact trace zeros");
  }
}

void CheckPreparedChart(const CutFemLevelSetData& parent, const CutFemLevelSetData& chart, unsigned points, std::mt19937& random) {
  Check(chart.GetRootDimension() == parent.GetRootDimension() && chart.GetParentMeasureScale() == 1,"chart metadata");
  for(unsigned sample = 0; sample < points; ++sample) {
    const std::vector<long double> point = RandomPoint(chart.GetGeometry(),random), parentPoint = chart.MapToParent(point);
    CheckNear(chart.Evaluate(point),parent.Evaluate(parentPoint),2e-15L,"chart preserves the polynomial");
    CheckNear(chart.EvaluateBernstein(point),chart.Evaluate(point),3e-15L,"chart preserves the Bernstein controls");
    const std::vector<long double> rootPoint = chart.MapToRoot(point), expectedRootPoint = parent.MapToRoot(parentPoint);
    for(unsigned axis = 0; axis < rootPoint.size(); ++axis) CheckNear(rootPoint[axis],expectedRootPoint[axis],2e-18L,"chart composes the original coordinate map");
    const std::vector<long double> gradient = chart.GetGradient(point), parentGradient = parent.GetGradient(parentPoint), rootGradient = chart.GetRootGradient(point), expectedRootGradient = parent.GetRootGradient(parentPoint);
    for(unsigned axis = 0; axis < gradient.size(); ++axis) { long double expected = 0; for(unsigned row = 0; row < parentGradient.size(); ++row) expected += parentGradient[row]*chart.GetParentMap()[row][axis]; CheckNear(gradient[axis],expected,3e-15L,"chart gradient follows the chain rule"); }
    for(unsigned axis = 0; axis < rootGradient.size(); ++axis) CheckNear(rootGradient[axis],expectedRootGradient[axis],3e-15L,"chart preserves original component order and signs");
  }
}

void TestPreparedCharts(unsigned cases, unsigned points, std::mt19937& random) {
  currentCheck = "supported axis and barycentric charts";
  std::uniform_real_distribution<double> value(-1,1);
  for(GeomElType geometry : {LINE,QUAD,TRI,HEX,TET,WEDGE}) {
    const std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    for(unsigned test = 0; test < std::min(cases,32U); ++test) {
      std::vector<double> nodal(basis->_nc);
      for(double& coefficient : nodal) coefficient = value(random);
      CutFemLevelSetData parent(geometry,nodal), chart, child;
      parent.SetContext(2,47,13);
      parent.GetZeroBoundaryMask();
      std::array<unsigned,3> axes = {0,1,2};
      do {
        if(geometry == WEDGE && axes[2] != 2) continue;
        parent.PermuteAxes(axes,chart);
        CheckPreparedChart(parent,chart,std::min(points,8U),random);
        Check(chart.GetMeshLevel() == 2 && chart.GetElementId() == 47 && chart.GetNodalVersion() == 13,"chart context is retained");
        if(chart.GetDimension() > 1) {
          const unsigned axis = chart.GetDimension()-1;
          const bool simplex = geometry == TRI || geometry == TET || (geometry == WEDGE && axis < 2);
          chart.Restrict(axis,0.2L,child);
          const std::vector<long double> point = RandomPoint(child.GetGeometry(),random);
          CheckNear(child.Evaluate(point),parent.Evaluate(child.MapToRoot(point)),3e-15L,"axis chart supports recursive restriction");
        }
      } while(std::next_permutation(axes.begin(),axes.begin()+parent.GetDimension()));
      if(geometry == TRI || geometry == TET || geometry == WEDGE) {
        std::array<unsigned,4> vertices = {0,1,2,3};
        const unsigned count = geometry == TET ? 4U : 3U;
        do { parent.PermuteSimplexVertices(vertices,chart); CheckPreparedChart(parent,chart,std::min(points,8U),random); } while(std::next_permutation(vertices.begin(),vertices.begin()+count));
        if(geometry != TRI) {
          parent.Restrict(geometry == TET ? 0U : 2U,0.2L,child);
          child.PermuteSimplexVertices({1,0,2,3},chart);
          CheckPreparedChart(child,chart,std::min(points,8U),random);
        }
      }
    }
  }
  CutFemLevelSetData wedge(WEDGE,MakeField(WEDGE,2)), quad, swapped, line;
  wedge.Restrict(0,0.2L,quad);
  quad.PermuteAxes({1,0,2},swapped);
  Check(swapped.GetBernsteinOrders() == std::array<unsigned,3>{{2,3,0}},"swapped wedge slice retains the elevated Bernstein degree");
  for(unsigned axis = 0; axis < 2; ++axis) { swapped.Restrict(axis,0.3L,line); for(unsigned sample = 0; sample < 8; ++sample) { const std::vector<long double> point = RandomPoint(LINE,random); CheckNear(line.EvaluateBernstein(point),wedge.Evaluate(line.MapToRoot(point)),3e-15L,"both directions of an elevated swapped quad restrict correctly"); } }
  wedge.Restrict(2,0.1L,quad,false);
  quad.PermuteSimplexVertices({1,0,2,3},swapped);
  Check(!swapped.HasRootGradientContext(),"charts inherit omitted original-gradient context");
  const std::vector<long double> point = RandomPoint(TRI,random);
  const long double expected = swapped.Evaluate(point);
  wedge.Reset(WEDGE,std::vector<double>(21,1));
  CheckNear(swapped.Evaluate(point),expected,0,"charted slices own their state after parent reset");
}

void TestProjectionAndGraphData(unsigned cases, std::mt19937& random) {
  currentCheck = "projection algebra, graph certificates, and unit export";
  std::uniform_real_distribution<double> value(-1,1);
  for(GeomElType geometry : {LINE,QUAD,TRI,HEX,TET,WEDGE}) {
    const std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    for(unsigned test = 0; test < cases; ++test) {
      std::vector<double> nodal(basis->_nc);
      for(double& coefficient : nodal) coefficient = value(random);
      CutFemLevelSetData data(geometry,nodal);
      const std::vector<long double> point = RandomPoint(geometry,random);
      std::vector<long double> unit = point;
      CheckNear(CutFemLevelSetData::EvaluatePolynomial(data.GetUnitPowerCoefficients(),data.GetDimension(),2,unit),data.Evaluate(point),3e-15L,"native-to-unit coefficient export");
      CheckNear(data.GetUnitToReferenceMeasureScale(),geometry == LINE || geometry == QUAD || geometry == HEX ? static_cast<long double>(1U<<data.GetDimension()) : (geometry == WEDGE ? 2 : 1),0,"unit-to-native measure scale");
      for(unsigned axis = 0; axis < data.GetDimension(); ++axis) {
        const CutFemLevelSetData::GraphCertificate& graph = data.GetGraphCertificate(axis);
        const long double derivative = data.EvaluateDerivative(axis,point);
        Check(derivative >= graph.bounds[0] && derivative <= graph.bounds[1],"graph bounds enclose sampled derivatives");
        if(graph.strict) Check(graph.sign*derivative > 0 && graph.zeroFiber == -1,"strict graph certificate");
        bool depends = false;
        unsigned stride = 1;
        for(unsigned next = axis+1; next < data.GetDimension(); ++next) stride *= 3;
        for(unsigned n = 0; n < data.GetPowerCoefficients().size(); ++n) depends = depends || ((n/stride)%3 != 0 && data.GetPowerCoefficients()[n] != 0);
        Check(depends == data.DependsOn(axis),"dependence follows retained nonzero monomials");
        if(data.GetDimension() == 1) continue;
        const CutFemLevelSetData::FiberProjection& projection = data.GetFiberProjection(axis);
        Check(&projection == &data.GetFiberProjection(axis) && projection.degree == 4,"projection is retained with its full degree");
        std::vector<long double> outer;
        for(unsigned direction = 0; direction < data.GetDimension(); ++direction) if(direction != axis) outer.push_back(point[direction]);
        std::array<boost::multiprecision::cpp_bin_float_100,3> coefficient;
        for(unsigned power = 0; power < 3; ++power) coefficient[power] = BenchmarkValue(std::vector<boost::multiprecision::cpp_bin_float_100>(projection.coefficients[power].begin(),projection.coefficients[power].end()),projection.dimension,std::vector<boost::multiprecision::cpp_bin_float_100>(outer.begin(),outer.end()));
        const long double expected = static_cast<long double>(coefficient[1]*coefficient[1]-4*coefficient[2]*coefficient[0]);
        CheckNear(data.EvaluateFiberDiscriminant(axis,outer),expected,2e-18L,"compensated discriminant versus 100-digit reference");
        CheckNear(CutFemLevelSetData::EvaluatePolynomial(projection.discriminant,projection.dimension,4,outer),expected,2e-16L,"full projection polynomial");
        long double reconstructed = 0;
        for(unsigned power = 3; power > 0; --power) reconstructed = reconstructed*point[axis]+static_cast<long double>(coefficient[power-1]);
        CheckNear(reconstructed*data.GetNormalizationScale(),data.Evaluate(point),2e-15L,"fiber coefficients reconstruct the parent polynomial");
        if(data.GetDimension() == 2) {
          std::vector<long double> boundary;
          data.GetFiberBoundaryPolynomial(axis,0.3L,-0.2L,boundary);
          std::vector<long double> boundaryPoint = point;
          boundaryPoint[axis] = 0.3L-0.2L*outer[0];
          CheckNear(CutFemLevelSetData::EvaluatePolynomial(boundary,1,4,outer)*data.GetNormalizationScale(),CutFemLevelSetData::EvaluatePolynomial(data.GetPowerCoefficients(),2,2,boundaryPoint),2e-15L,"affine boundary polynomial");
          std::array<long double,5> shifted;
          data.GetShiftedFiberDiscriminant(axis,outer[0],shifted);
          CheckNear(shifted[0],expected,2e-18L,"shifted discriminant retains compensated endpoint value");
          const long double offset = 0.01L;
          CheckNear(CutFemLevelSetData::EvaluatePolynomial(std::vector<long double>(shifted.begin(),shifted.end()),1,4,{offset}),data.EvaluateFiberDiscriminant(axis,{outer[0]+offset}),2e-16L,"shifted discriminant polynomial");
        }
      }
    }
  }
  const std::unique_ptr<femus::basis> basis = MakeBasis(QUAD);
  std::vector<double> nodal(9);
  for(unsigned n = 0; n < 9; ++n) { const double* point = basis->GetXcoarse(n); nodal[n] = point[0]*((point[1]-0.2)*(point[1]-0.2)+0.03); }
  CutFemLevelSetData tight(QUAD,nodal);
  Check(tight.GetDerivativeBounds(0)[0] < 0 && tight.GetGraphCertificate(0).bounds[0] > 0.029L && tight.GetGraphDirection() == 0,"tight quadratic extrema certify a graph missed by Bernstein");
  const std::unique_ptr<femus::basis> hexBasis = MakeBasis(HEX);
  std::vector<double> hexNodes(hexBasis->_nc);
  for(unsigned axis = 0; axis < 3; ++axis) for(double shift : {0.5,0.0,-0.03125}) for(int sign : {-1,1}) {
    for(unsigned n = 0; n < hexNodes.size(); ++n) {
      const double* point = hexBasis->GetXcoarse(n);
      const double u = point[(axis+1)%3]-0.25, v = point[(axis+2)%3]+0.25;
      hexNodes[n] = sign*(point[axis]-0.125)*(u*u+v*v+shift);
    }
    const CutFemLevelSetData hex(HEX,hexNodes);
    const CutFemLevelSetData::GraphCertificate& graph = hex.GetGraphCertificate(axis);
    Check(graph.strict == (shift > 0),"HEX strict graph requires a sign-definite derivative");
    if(shift > 0) Check(graph.sign == sign && hex.GetGraphDirection() == static_cast<int>(axis),"HEX finds the monotone axis for either orientation");
  }
  for(unsigned n = 0; n < 9; ++n) { const double* point = basis->GetXcoarse(n); nodal[n] = point[0]*(point[1]+1); }
  CutFemLevelSetData zeroFiber(QUAD,nodal);
  Check(zeroFiber.GetGraphCertificate(0).sign == 1 && !zeroFiber.GetGraphCertificate(0).strict && zeroFiber.GetGraphCertificate(0).zeroFiber == 1,"weak graph detects an entire zero fiber");
  Check(zeroFiber.GetGraphDirection(true) == -1,"zero fibers cannot become regular graphs");
  for(unsigned n = 0; n < 9; ++n) { const double* point = basis->GetXcoarse(n); nodal[n] = (point[0]+1)*(point[0]+1)+0.25*point[1]*point[1]; }
  CutFemLevelSetData weak(QUAD,nodal);
  Check(weak.GetGraphCertificate(0).sign == 1 && !weak.GetGraphCertificate(0).strict && weak.GetGraphCertificate(0).zeroFiber == -1,"weak graph excludes zero fibers");
  Check(weak.GetGraphDirection() == -1 && weak.GetGraphDirection(true) == 0,"weak directions require an explicit request");
  Check(weak.GetReflectionParity(1) == 1 && weak.GetReflectionParity(0) == 0,"native reflection parity");
  for(unsigned n = 0; n < 9; ++n) { const double* point = basis->GetXcoarse(n); nodal[n] = point[0]*point[0]*point[1]*point[1]; }
  CutFemLevelSetData quartic(QUAD,nodal);
  std::vector<long double> boundary;
  quartic.GetFiberBoundaryPolynomial(1,0.5L,1,boundary);
  CheckNear(boundary[4]*quartic.GetNormalizationScale(),16,0,"projected boundary quartic coefficient is retained");
  CutFemLevelSetData child;
  quartic.Restrict(0,0.2L,child);
  child.GetLineAnalysis();
  const CutFemLevelSetData savedLine = child;
  child.GetUnitPowerCoefficients();
  quartic.Reset(QUAD,std::vector<double>(9,1));
  quartic.Restrict(0,0.2L,child);
  Check(child.GetLineAnalysis().rootCount == 0 && !child.DependsOn(0),"reused child invalidates line and dependence caches");
  Check(savedLine.GetLineAnalysis().rootCount == 1,"copied line analysis survives reuse of the original workspace");
  CheckNear(child.GetUnitPowerCoefficients()[0],1,0,"reused child invalidates unit coefficients");
  CutFemLevelSetData spatial(HEX,MakeField(HEX,2,1,0.29L)), slice;
  spatial.Restrict(2,0.1L,slice);
  const std::vector<long double> savedProjection = slice.GetFiberProjection(0).discriminant;
  const CutFemLevelSetData::GraphCertificate savedGraph = slice.GetGraphCertificate(0);
  const CutFemLevelSetData copied = slice;
  spatial.Restrict(2,0.2L,slice);
  slice.GetFiberProjection(0);
  slice.GetGraphCertificate(0);
  Check(copied.GetFiberProjection(0).discriminant == savedProjection,"copied projection caches survive reuse of the original workspace");
  Check(slice.GetFiberProjection(0).discriminant != savedProjection,"reused projection follows the new slice");
  Check(copied.GetGraphCertificate(0).bounds == savedGraph.bounds,"copied graph certificates survive reuse of the original workspace");
  CutFemLevelSetData partial = tight;
  const CutFemLevelSetData::GraphCertificate firstGraph = partial.GetGraphCertificate(0);
  const CutFemLevelSetData partialCopy = partial;
  partial.GetGraphCertificate(1);
  Check(partial.GetGraphCertificate(0).bounds == firstGraph.bounds && partial.GetGraphCertificate(0).strict == firstGraph.strict,"preparing another graph axis preserves existing detached cache entries");
  Check(partialCopy.GetGraphCertificate(0).bounds == firstGraph.bounds,"partial cache copies remain independent");
}

void TestExternalChartComponents(std::mt19937& random) {
  currentCheck = "full-degree external component chart transforms";
  std::uniform_real_distribution<double> value(-0.1,0.1);
  for(GeomElType geometry : {QUAD,TRI,HEX,TET,WEDGE}) {
    CutFemLevelSetData parent(geometry,MakeField(geometry,1)), chart;
    if(geometry == TRI || geometry == TET || geometry == WEDGE) parent.PermuteSimplexVertices({1,0,2,3},chart);
    else parent.PermuteAxes(geometry == QUAD ? std::array<unsigned,3>{{1,0,2}} : std::array<unsigned,3>{{2,0,1}},chart);
    for(unsigned degree : {2U,3U}) {
      unsigned count = 1, outputDegree = 0;
      for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) count *= degree+1;
      std::vector<std::vector<long double>> components(5,std::vector<long double>(count)), output;
      for(std::vector<long double>& component : components) for(long double& coefficient : component) coefficient = value(random);
      chart.TransformParentComponents(degree,components,outputDegree,output);
      Check(output.size() == components.size() && outputDegree >= degree,"chart preserves every supplied component and its degree");
      if(geometry == TRI || geometry == TET || geometry == WEDGE) Check(outputDegree == 2*degree,"barycentric mixing retains the increased coordinate degree");
      for(unsigned sample = 0; sample < 16; ++sample) {
        const std::vector<long double> point = RandomPoint(geometry,random), parentPoint = chart.MapToParent(point);
        for(unsigned component = 0; component < components.size(); ++component) CheckNear(CutFemLevelSetData::EvaluatePolynomial(output[component],chart.GetDimension(),outputDegree,point),CutFemLevelSetData::EvaluatePolynomial(components[component],parent.GetDimension(),degree,parentPoint),3e-15L,"external polynomial chart composition");
      }
    }
  }
}

void TestRootSurfaceInterfaces() {
  currentCheck = "selected densities from prepared root records";
  for(GeomElType geometry : {LINE,QUAD,TRI,HEX,TET,WEDGE}) {
    CutFemLevelSetData parent(geometry,MakeField(geometry,1)), planar, line;
    if(parent.GetDimension() == 3) { parent.Restrict(2,0.1L,planar); planar.Restrict(1,0.2L,line); }
    else if(parent.GetDimension() == 2) parent.Restrict(1,0.2L,line);
    else line = parent;
    const CutFemLevelSetData::LineAnalysis& analysis = line.GetLineAnalysis();
    Check(analysis.rootCount == 1,"regular prepared root for every element");
    const CutFemLevelSetData::LineRoot& root = analysis.roots[0];
    std::array<std::array<long double,3>,3> matrix = {{{1.2L,0.1L,0.2L},{0.2L,1.1L,0.1L},{0.1L,0.2L,1.3L}}};
    const std::array<long double,3> h = {1,0.2L,0.3L};
    std::array<long double,4> scalar, coordinateScalar, matrixScalar;
    std::array<long double,3> normal, coordinateNormal, matrixNormal, gradient;
    line.GetLineAllSurfaceFactors(root,h,scalar,normal);
    line.GetLineAllSurfaceFactors(root.coordinate,h,coordinateScalar,coordinateNormal);
    Check(scalar == coordinateScalar && normal == coordinateNormal,"prepared-root batch compatibility");
    line.GetLineAllSurfaceFactorsFromCofactor(root,matrix,matrixScalar,matrixNormal);
    line.GetLineAllSurfaceFactorsFromCofactor(root.coordinate,matrix,coordinateScalar,coordinateNormal);
    Check(matrixScalar == coordinateScalar && matrixNormal == coordinateNormal,"prepared-root matrix compatibility");
    CheckNear(line.GetLineDeltaFactor(root),scalar[0],0,"prepared-root delta");
    CheckNear(line.GetLineSurfaceData(root,gradient),scalar[0],0,"prepared-root gradient data");
    CheckNear(line.GetLineReferenceSurfaceFactor(root),scalar[1],0,"prepared-root reference measure");
    CheckNear(line.GetLinePhysicalSurfaceFactor(root,h),scalar[2],0,"prepared-root physical measure");
    CheckNear(line.GetLineTangentialSurfaceFactor(root,h),scalar[3],0,"prepared-root tangential measure");
    CheckNear(line.GetLinePhysicalSurfaceFactorFromCofactor(root,matrix),matrixScalar[2],0,"prepared-root matrix physical measure");
    CheckNear(line.GetLineTangentialSurfaceFactorFromCofactor(root,matrix),matrixScalar[3],0,"prepared-root matrix tangential measure");
    for(unsigned component = 0; component < line.GetRootDimension(); ++component) { CheckNear(line.GetLineNormalSurfaceFactor(root,h,component),normal[component],0,"prepared-root normal component"); CheckNear(line.GetLineNormalSurfaceFactorFromCofactor(root,matrix,component),matrixNormal[component],0,"prepared-root matrix normal component"); }
    line.GetLineNormalSurfaceFactors(root,h,coordinateNormal);
    Check(coordinateNormal == normal,"prepared-root normal family");
    line.GetLineNormalSurfaceFactorsFromCofactor(root,matrix,coordinateNormal);
    Check(coordinateNormal == matrixNormal,"prepared-root matrix normal family");
    if(geometry == TRI || geometry == TET || geometry == WEDGE) {
      CutFemLevelSetData chart;
      parent.PermuteSimplexVertices({1,0,2,3},chart);
      if(chart.GetDimension() == 3) { chart.Restrict(2,0.1L,planar); planar.Restrict(0,0.1L,line); }
      else chart.Restrict(0,0.1L,line);
      const CutFemLevelSetData::LineAnalysis& chartAnalysis = line.GetLineAnalysis();
      Check(chartAnalysis.rootCount == 1,"regular root in an oblique barycentric chart");
      const CutFemLevelSetData::LineRoot& chartRoot = chartAnalysis.roots[0];
      line.GetLineSurfaceData(chartRoot,gradient);
      const std::vector<long double> native = line.MapToRoot({chartRoot.coordinate}), expected = parent.GetGradient(native);
      long double directional = 0;
      for(unsigned component = 0; component < parent.GetDimension(); ++component) { CheckNear(gradient[component],expected[component],3e-15L,"oblique-root gradient retains original components"); }
      const std::vector<long double> first = line.MapToRoot({0}), second = line.MapToRoot({1});
      for(unsigned component = 0; component < parent.GetDimension(); ++component) directional += (second[component]-first[component])*gradient[component];
      CheckNear(directional,chartRoot.derivative,3e-18L,"oblique-root gradient agrees with the retained directional slope");
    }
  }
}

void TestLocalFiberAnalysis(unsigned cases, std::mt19937& random) {
  currentCheck = "local unit fibers without LINE construction";
  std::uniform_real_distribution<double> value(-1,1), position(0.01,0.9);
  for(GeomElType geometry : {QUAD,TRI,HEX,TET,WEDGE}) {
    const std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    for(unsigned test = 0; test < cases; ++test) {
      std::vector<double> nodal(basis->_nc);
      for(double& node : nodal) node = value(random);
      CutFemLevelSetData root(geometry,nodal), slice, local, line;
      if(root.GetDimension() == 3) root.Restrict(2,position(random),slice,false);
      else slice = root;
      local.ResetPowerCoefficients(slice.GetGeometry(),slice.GetPowerCoefficients());
      for(unsigned axis = 0; axis < 2; ++axis) for(long double fixed : {0.0L,0.25L,0.5L,0.75L,0.99L}) {
        const unsigned long long restrictions = local.GetWorkCounts()[3];
        CutFemLevelSetData::LineAnalysis actual;
        local.AnalyzeLocalFiber(axis,fixed,actual);
        Check(local.GetWorkCounts()[3] == restrictions,"direct fiber analysis creates no restricted object");
        local.Restrict(1-axis,fixed,line,false);
        const CutFemLevelSetData::LineAnalysis& expected = line.GetLineAnalysis();
        const long double width = local.GetGeometry() == TRI ? 1-fixed : 1;
        Check(actual.rootCount == expected.rootCount && actual.positiveIntervalCount == expected.positiveIntervalCount && actual.identicallyZero == expected.identicallyZero,"direct local topology agrees with independently restricted LINE analysis");
        for(unsigned n = 0; n < actual.rootCount; ++n) {
          CheckNear(actual.roots[n].coordinate,width*expected.roots[n].coordinate,3e-16L,"fiber roots use parent unit coordinates");
          CheckNear(actual.roots[n].derivative,expected.roots[n].derivative/width,3e-15L,"fiber slopes use the parent unit coordinate");
          Check(actual.roots[n].boundary == expected.roots[n].boundary,"fiber boundary tags");
        }
        for(unsigned n = 0; n < actual.positiveIntervalCount; ++n) for(unsigned side = 0; side < 2; ++side) CheckNear(actual.positiveIntervals[n][side],width*expected.positiveIntervals[n][side],3e-16L,"fiber positive intervals");
      }
    }
  }
  CutFemLevelSetData zero(QUAD,std::vector<double>(9,0));
  CutFemLevelSetData::LineAnalysis analysis;
  zero.AnalyzeLocalFiber(0,0.3L,analysis);
  Check(analysis.identicallyZero && analysis.rootCount == 0,"zero local fibers are explicit");
}

void TestPreparedDataInvalidInputs() {
  currentCheck = "prepared-data interface validation";
  CutFemLevelSetData quad(QUAD,MakeField(QUAD,1)), tri(TRI,MakeField(TRI,1)), wedge(WEDGE,MakeField(WEDGE,1)), child;
  unsigned rejected = 0;
  try { quad.PermuteAxes({0,0,2},child); } catch(const std::invalid_argument&) { ++rejected; }
  try { wedge.PermuteAxes({0,2,1},child); } catch(const std::invalid_argument&) { ++rejected; }
  try { quad.PermuteAxes({1,0,2},quad); } catch(const std::invalid_argument&) { ++rejected; }
  try { tri.PermuteSimplexVertices({1,1,2,3},child); } catch(const std::invalid_argument&) { ++rejected; }
  try { quad.PermuteSimplexVertices({1,0,2,3},child); } catch(const std::logic_error&) { ++rejected; }
  try { quad.RestrictBoundary(4,child); } catch(const std::out_of_range&) { ++rejected; }
  try { quad.GetLineAnalysis(); } catch(const std::logic_error&) { ++rejected; }
  try { quad.GetFiberProjection(2); } catch(const std::out_of_range&) { ++rejected; }
  try { quad.GetZeroFiberStatus(2); } catch(const std::out_of_range&) { ++rejected; }
  try { quad.GetReflectionParity(2); } catch(const std::out_of_range&) { ++rejected; }
  Check(rejected == 10,"invalid charts, boundaries, and analysis requests are rejected");
  tri.PermuteSimplexVertices({1,0,2,3},child);
  std::vector<std::vector<long double>> components(1,std::vector<long double>(9,1)), output;
  unsigned degree = 0;
  bool aliasRejected = false;
  try { child.TransformParentComponents(2,components,degree,components); } catch(const std::invalid_argument&) { aliasRejected = true; }
  Check(aliasRejected,"external chart transform rejects aliased input");
  components[0].assign(18*18,0);
  components[0].back() = 1;
  bool degreeRejected = false;
  try { child.TransformParentComponents(17,components,degree,output); } catch(const std::invalid_argument&) { degreeRejected = true; }
  Check(degreeRejected,"unsupported expanded degrees are rejected without truncation");
}

long double DataReuseSample(const CutFemLevelSetData& planar, unsigned operation, long double coordinate, CutFemLevelSetData& line) {
  if(operation < 3) planar.Restrict(1,coordinate,line,false);
  if(operation == 0) return line.GetPowerCoefficients()[0];
  std::array<long double,2> roots;
  const unsigned count = line.GetLineRoots(roots);
  if(operation == 1 || operation == 3) return count == 0 ? 0 : count+roots[0];
  long double value = 0;
  for(unsigned n = 0; n < count; ++n) value += line.GetLinePhysicalSurfaceFactor(roots[n],{1,0.2L,0.3L});
  return value;
}

void BenchmarkDataReuse() {
  std::cout << "data_record bytes=" << sizeof(CutFemLevelSetData) << '\n';
  const std::string names[5] = {"slice","slice_roots","slice_physical","cached_roots","cached_physical"};
  const unsigned iterations = 20000;
  for(GeomElType geometry : {LINE,QUAD,TRI,HEX,TET,WEDGE}) {
    CutFemLevelSetData parent(geometry,MakeField(geometry,2,1,0.29L)), planar, line;
    if(parent.GetDimension() == 3) parent.Restrict(2,0.19L,planar,false);
    else if(parent.GetDimension() == 2) planar = parent;
    else line = parent;
    for(unsigned operation = parent.GetDimension() == 1 ? 3U : 0U; operation < 5; ++operation) {
      volatile long double checksum = 0;
      for(unsigned sample = 0; sample < 128; ++sample) checksum += DataReuseSample(planar,operation,0.15L+0.001L*(sample%64),line);
#ifdef EX47_TRACK_ALLOCATIONS
      const unsigned long long before = heapAllocations;
#endif
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      for(unsigned sample = 0; sample < iterations; ++sample) checksum += DataReuseSample(planar,operation,0.15L+0.001L*(sample%64),line);
      const double nanoseconds = 1e9*std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/iterations;
#ifdef EX47_TRACK_ALLOCATIONS
      const unsigned long long allocated = heapAllocations-before;
      Check(allocated == 0,"prepared slice/root reuse performs no warmed heap allocation");
#endif
      Check(std::isfinite(checksum),"prepared-data benchmark checksum");
      std::cout << std::setprecision(10) << "data_reuse " << GeometryName(geometry) << ' ' << names[operation] << " ns=" << nanoseconds << " checksum=" << checksum;
#ifdef EX47_TRACK_ALLOCATIONS
      std::cout << " allocations=" << allocated;
#endif
      std::cout << '\n';
    }
  }
}

long double PreparedAlgebraSample(const CutFemLevelSetData& parent, const CutFemLevelSetData& line, unsigned operation, long double coordinate, CutFemLevelSetData& chart, const std::vector<std::vector<long double>>& components, std::vector<std::vector<long double>>& output) {
  if(operation == 0) return parent.GetZeroBoundaryMask();
  if(operation == 1) return parent.DependsOn(0)+parent.GetReflectionParity(0);
  if(operation == 2) return parent.GetGraphCertificate(0).bounds[0];
  if(operation == 3) return parent.GetUnitPowerCoefficients()[0];
  if(operation == 4) return parent.GetFiberProjection(0).discriminant[0];
  if(operation == 5) return parent.EvaluateFiberDiscriminant(0,coordinate);
  if(operation == 6) { parent.PermuteAxes({1,0,2},chart); return chart.GetPowerCoefficients()[0]; }
  if(operation == 7) { parent.PermuteSimplexVertices({1,0,2,3},chart); return chart.GetPowerCoefficients()[0]; }
  if(operation == 8) { unsigned degree = 0; chart.TransformParentComponents(3,components,degree,output); return output[0][0]; }
  const CutFemLevelSetData::LineAnalysis& analysis = line.GetLineAnalysis();
  return analysis.rootCount == 0 ? 0 : line.GetLineDeltaFactor(analysis.roots[0]);
}

void BenchmarkPreparedAlgebra() {
  const std::string names[10] = {"boundary_cached","dependence_parity_cached","graph_cached","unit_cached","projection_cached","discriminant_compensated","axis_chart","barycentric_chart","external_chart_degree3","root_record_delta"};
  for(GeomElType geometry : {QUAD,TRI,HEX,TET,WEDGE}) {
    CutFemLevelSetData parent(geometry,MakeField(geometry,2,1,0.29L)), chart, planar, line;
    if(parent.GetDimension() == 3) { parent.Restrict(2,0.19L,planar,false); planar.Restrict(1,0.17L,line,false); }
    else parent.Restrict(1,0.17L,line,false);
    unsigned size = 1;
    for(unsigned axis = 0; axis < parent.GetDimension(); ++axis) size *= 4;
    std::vector<std::vector<long double>> components(3,std::vector<long double>(size,0)), output;
    for(unsigned component = 0; component < 3; ++component) for(unsigned n = 0; n < size; ++n) components[component][n] = 0.001L*(1+component+n);
    for(unsigned operation = 0; operation < 10; ++operation) {
      if(operation == 5 && parent.GetDimension() != 2) continue;
      if(operation == 7 && geometry != TRI && geometry != TET && geometry != WEDGE) continue;
      if(operation == 8) { if(geometry == TRI || geometry == TET || geometry == WEDGE) parent.PermuteSimplexVertices({1,0,2,3},chart); else parent.PermuteAxes({1,0,2},chart); }
      const unsigned iterations = operation == 8 ? 128U : 20000U;
      volatile long double checksum = 0;
      for(unsigned sample = 0; sample < 64; ++sample) checksum += PreparedAlgebraSample(parent,line,operation,0.15L+0.001L*(sample%64),chart,components,output);
#ifdef EX47_TRACK_ALLOCATIONS
      const unsigned long long before = heapAllocations;
#endif
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      for(unsigned sample = 0; sample < iterations; ++sample) checksum += PreparedAlgebraSample(parent,line,operation,0.15L+0.001L*(sample%64),chart,components,output);
      const double nanoseconds = 1e9*std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/iterations;
#ifdef EX47_TRACK_ALLOCATIONS
      const unsigned long long allocated = heapAllocations-before;
      Check(allocated == 0,"prepared algebra performs no warmed heap allocation");
#endif
      Check(std::isfinite(checksum),"prepared algebra benchmark checksum");
      std::cout << std::setprecision(10) << "prepared_algebra " << GeometryName(geometry) << ' ' << names[operation] << " ns=" << nanoseconds << " checksum=" << checksum;
#ifdef EX47_TRACK_ALLOCATIONS
      std::cout << " allocations=" << allocated;
#endif
      std::cout << '\n';
    }
  }
}

template<class Weight> void GetWeightPair(Weight& quadrature, unsigned mode, const std::vector<double>& nodal, const std::vector<std::vector<double>>& coordinates, std::vector<long double>& positive, std::vector<long double>& surface, double tolerance = 1e-10) {
  if(mode == 0 || mode == 1) quadrature.GetWeightsFromNodalValues(nodal,positive,surface,tolerance,tolerance,mode == 1);
  else if(mode == 2) quadrature.GetPhysicalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,surface,tolerance,tolerance);
  else if(mode == 3) quadrature.GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,surface,tolerance,tolerance);
  else { std::vector<std::vector<long double>> normal; quadrature.GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,normal,tolerance,tolerance); surface = normal[0]; }
}

void ReferenceLineWeightMoments(const std::vector<double>& nodal, unsigned degree, unsigned mode, std::vector<long double>& positive, std::vector<long double>& surface) {
  const boost::multiprecision::cpp_bin_float_100 a = (boost::multiprecision::cpp_bin_float_100(nodal[0])+boost::multiprecision::cpp_bin_float_100(nodal[1]))/2-boost::multiprecision::cpp_bin_float_100(nodal[2]), b = (boost::multiprecision::cpp_bin_float_100(nodal[1])-boost::multiprecision::cpp_bin_float_100(nodal[0]))/2, c(nodal[2]);
  std::vector<boost::multiprecision::cpp_bin_float_100> roots, cuts(1,-1);
  const boost::multiprecision::cpp_bin_float_100 discriminant = b*b-4*a*c;
  if(a == 0 && b != 0) roots.push_back(-c/b);
  else if(a != 0 && discriminant > 0) { roots.push_back((-b-sqrt(discriminant))/(2*a)); roots.push_back((-b+sqrt(discriminant))/(2*a)); }
  std::sort(roots.begin(),roots.end());
  for(const boost::multiprecision::cpp_bin_float_100& root : roots) if(root > -1 && root < 1) cuts.push_back(root);
  cuts.push_back(1);
  positive.assign(degree+1,0);
  surface.assign(degree+1,0);
  for(unsigned power = 0; power <= degree; ++power) {
    boost::multiprecision::cpp_bin_float_100 h = 0, d = 0;
    for(unsigned interval = 1; interval < cuts.size(); ++interval) { const boost::multiprecision::cpp_bin_float_100 middle = (cuts[interval-1]+cuts[interval])/2; if((a*middle+b)*middle+c > 0) h += (pow((cuts[interval]+1)/2,power+1)-pow((cuts[interval-1]+1)/2,power+1))/(power+1); }
    for(const boost::multiprecision::cpp_bin_float_100& root : roots) if(root >= -1 && root <= 1) {
      const boost::multiprecision::cpp_bin_float_100 derivative = abs(2*a*root+b), factor = root == -1 || root == 1 ? 0.25 : 0.5;
      d += factor*pow((root+1)/2,power)*(mode == 0 ? 1/derivative : (mode == 3 ? 1/(derivative*derivative) : boost::multiprecision::cpp_bin_float_100(1)));
    }
    positive[power] = static_cast<long double>(2*h);
    surface[power] = static_cast<long double>(2*d);
  }
}

template<class Weight> std::vector<long double> FittedLineMoments(const Weight& quadrature, const std::vector<long double>& weights, unsigned degree) {
  std::vector<long double> result(degree+1,0);
  for(unsigned point = 0; point < quadrature.GetGaussQuadraturePointNumber(); ++point) {
    const long double x = (static_cast<long double>(quadrature.GetGaussCoordinatePointer(0)[point])+1)/2, weight = quadrature.GetGaussWeightPointer()[point]*weights[point];
    long double power = 1;
    for(unsigned i = 0; i <= degree; ++i) { result[i] += weight*power; power *= x; }
  }
  return result;
}

void TestLineWeightIntegration() {
  currentCheck = "1D native fitted-weight integration";
  const std::vector<std::vector<double>> fields = {{1,1,1},{-1,-1,-1},{-1,1,0},{0,1,0.5},{-1,0,-0.5},{0,0,-1},{0.75,0.75,-0.25},{-0.75,-0.75,0.25},{1,1,-1e-24},{1e-100,1e-100,-2.5e-101},{1e100,1e100,-2.5e99}};
  const std::vector<std::vector<double>> coordinates = {{-0.85,1.35,0.2}};
  for(unsigned degree : {2U,5U,10U}) {
    CutFemWeight_LS<long double,boost::multiprecision::cpp_bin_float_quad> referenceFitter(LINE,degree,"legendre");
    CutFemWeight_LS<long double,boost::multiprecision::cpp_bin_float_quad> prepared(LINE,degree,"legendre");
    long double newError = 0, fittingFloor = 0, rawError = 0;
    for(const std::vector<double>& nodal : fields) for(unsigned mode = 0; mode < 4; ++mode) {
      std::vector<long double> newH, newD, expectedH, expectedD;
      GetWeightPair(prepared,mode,nodal,coordinates,newH,newD);
      ReferenceLineWeightMoments(nodal,degree,mode,expectedH,expectedD);
      std::vector<long double> rawH, rawD, exactH = expectedH, exactD = expectedD, floorH, floorD;
      prepared.GetLastIntegratedMoments(rawH,rawD);
      for(unsigned n = 0; n <= degree; ++n) { exactH[n] /= 2; exactD[n] /= 2; rawError = std::max(rawError,std::max(std::abs(2*rawH[n]-expectedH[n])/(1+std::abs(expectedH[n])),std::abs(2*rawD[n]-expectedD[n])/(1+std::abs(expectedD[n])))); CheckNear(2*rawH[n],expectedH[n],3e-17L,"LINE raw volume moments versus multiprecision"); CheckNear(2*rawD[n],expectedD[n],3e-17L,"LINE raw surface moments versus multiprecision"); }
      referenceFitter.GetWeightsFromExternalMoments(exactH,exactD,floorH,floorD);
      const std::vector<long double> floorHMoments = FittedLineMoments(referenceFitter,floorH,degree), floorDMoments = FittedLineMoments(referenceFitter,floorD,degree);
      const std::vector<long double> newHMoments = FittedLineMoments(prepared,newH,degree), newDMoments = FittedLineMoments(prepared,newD,degree);
      for(unsigned n = 0; n <= degree; ++n) {
        newError = std::max(newError,std::max(std::abs(newHMoments[n]-expectedH[n])/(1+std::abs(expectedH[n])),std::abs(newDMoments[n]-expectedD[n])/(1+std::abs(expectedD[n]))));
        const long double hFloor = std::abs(floorHMoments[n]-expectedH[n])/(1+std::abs(expectedH[n])), dFloor = std::abs(floorDMoments[n]-expectedD[n])/(1+std::abs(expectedD[n]));
        fittingFloor = std::max(fittingFloor,std::max(hFloor,dFloor));
        CheckNear(newHMoments[n],expectedH[n],2*hFloor+5e-14L,"prepared LINE volume residual respects the fixed fitting floor");
        CheckNear(newDMoments[n],expectedD[n],2*dFloor+5e-14L,"prepared LINE surface residual respects the fixed fitting floor");
      }
      Check(prepared.GetPreparedLevelSetData().GetWorkCounts()[0] == prepared.GetPreparedLevelSetData().GetNodalVersion(),"one nodal expansion per distinct LINE state");
    }
    std::cout << std::setprecision(8) << "weight_check LINE3 q=" << degree << " prepared_error=" << newError << " raw_error=" << rawError << " fitting_floor=" << fittingFloor << '\n';
  }
}

template<class Algebra> void GetNodalFixtureWeightPair(CutFemWeight_LS<long double,Algebra>& quadrature, unsigned mode, const CutFemLevelSetData& data, const std::vector<std::vector<double>>& coordinates, std::vector<long double>& positive, std::vector<long double>& surface) { GetWeightPair(quadrature,mode,data.GetNodalValues(),coordinates,positive,surface); }

template<class Algebra> void GetFixtureWeightPair(CutFemWeight_LS<long double,Algebra>& quadrature, unsigned mode, const CutFemLevelSetData& data, const std::vector<std::vector<double>>& coordinates, std::vector<long double>& positive, std::vector<long double>& surface) {
  if(mode == 0 || mode == 1) quadrature.GetWeightsFromLevelSetData(data,positive,surface,1e-10,1e-10,mode == 1);
  else if(mode == 2) quadrature.GetPhysicalSurfaceWeightsFromLevelSetData(data,coordinates,positive,surface,1e-10,1e-10);
  else if(mode == 3) quadrature.GetTangentialSurfaceWeightsFromLevelSetData(data,coordinates,positive,surface,1e-10,1e-10);
  else { std::vector<std::vector<long double>> normal; quadrature.GetNormalSurfaceWeightsFromLevelSetData(data,coordinates,positive,normal,1e-10,1e-10); surface = normal[0]; }
}

template<class Weight> double TimeLineWeights(Weight& quadrature, unsigned mode, const std::array<CutFemLevelSetData,2>& data, const std::vector<std::vector<double>>& coordinates, unsigned iterations, long double& checksum, bool nodalInput = false) {
  std::vector<long double> positive, surface;
  for(unsigned sample = 0; sample < (data[0].GetDimension() == 3 ? 2U : 16U); ++sample) { quadrature.ClearCache(); if(nodalInput) GetNodalFixtureWeightPair(quadrature,mode,data[sample%2],coordinates,positive,surface); else GetFixtureWeightPair(quadrature,mode,data[sample%2],coordinates,positive,surface); }
  volatile long double total = 0;
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  for(unsigned sample = 0; sample < iterations; ++sample) { quadrature.ClearCache(); if(nodalInput) GetNodalFixtureWeightPair(quadrature,mode,data[sample%2],coordinates,positive,surface); else GetFixtureWeightPair(quadrature,mode,data[sample%2],coordinates,positive,surface); total += positive[0]+surface[0]; }
  checksum = total;
  return 1e6*std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/iterations;
}

void BenchmarkLineWeights() {
  const std::array<CutFemLevelSetData,2> data = {{CutFemLevelSetData(LINE,{0.75,0.75,-0.25}),CutFemLevelSetData(LINE,{1.5,0.5,0})}};
  const std::vector<std::vector<double>> coordinates = {{-0.85,1.35,0.2}};
  for(unsigned degree : {2U,5U,10U}) {
    CutFemWeight_LS<long double,long double> nodalRule(LINE,degree,"legendre");
    CutFemWeight_LS<long double,long double> prepared(LINE,degree,"legendre");
    for(unsigned mode = 0; mode < 4; ++mode) {
      std::array<double,5> oldTime, newTime;
      long double oldChecksum = 0, newChecksum = 0;
      for(unsigned repeat = 0; repeat < 5; ++repeat) {
        if(repeat%2 == 0) { oldTime[repeat] = TimeLineWeights(nodalRule,mode,data,coordinates,2000,oldChecksum,true); newTime[repeat] = TimeLineWeights(prepared,mode,data,coordinates,2000,newChecksum); }
        else { newTime[repeat] = TimeLineWeights(prepared,mode,data,coordinates,2000,newChecksum); oldTime[repeat] = TimeLineWeights(nodalRule,mode,data,coordinates,2000,oldChecksum,true); }
      }
      std::sort(oldTime.begin(),oldTime.end());
      std::sort(newTime.begin(),newTime.end());
      CheckNear(newChecksum,oldChecksum,degree == 10 ? 2e-10L : 1e-12L,"LINE benchmark checksums");
      std::cout << std::setprecision(8) << "weight_speed LINE3 q=" << degree << " mode=" << mode << " nodal_us=" << oldTime[2] << " prepared_us=" << newTime[2] << " nodal_over_prepared=" << oldTime[2]/newTime[2] << '\n';
    }
  }
}

void WeightBasisPowers(GeomElType geometry, const std::array<long double,3>& point, unsigned degree, std::array<std::vector<long double>,3>& powers) {
  const bool tensor = geometry == LINE || geometry == QUAD || geometry == HEX;
  for(unsigned axis = 0; axis < 3; ++axis) { const long double coordinate = tensor || (geometry == WEDGE && axis == 2) ? (point[axis]+1)/2 : point[axis]; powers[axis].resize(degree+1); powers[axis][0] = 1; for(unsigned n = 1; n <= degree; ++n) powers[axis][n] = powers[axis][n-1]*coordinate; }
}

template<class Weight> std::vector<long double> FittedWeightMoments(const Weight& quadrature, GeomElType geometry, unsigned dimension, const std::vector<long double>& weights, unsigned degree) {
  const unsigned count = dimension == 2 ? (degree+1)*(degree+2)/2 : (degree+1)*(degree+2)*(degree+3)/6;
  std::vector<long double> result(count,0);
  std::array<std::vector<long double>,3> powers;
  for(unsigned n = 0; n < quadrature.GetGaussQuadraturePointNumber(); ++n) {
    std::array<long double,3> point = {0,0,0};
    for(unsigned axis = 0; axis < dimension; ++axis) point[axis] = quadrature.GetGaussCoordinatePointer(axis)[n];
    WeightBasisPowers(geometry,point,degree,powers);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= (dimension == 3 ? degree-i-j : 0U); ++k) result[FELevelSetMoments::UnitMomentIndex(dimension,i,j,k)] += quadrature.GetGaussWeightPointer()[n]*weights[n]*powers[0][i]*powers[1][j]*powers[2][k];
  }
  return result;
}

void MakeWeightBall(GeomElType geometry, std::vector<double>& nodal, std::vector<std::vector<double>>& coordinates, std::array<long double,3>& center, long double& radius) {
  const std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
  const unsigned dimension = geometry == QUAD || geometry == TRI ? 2U : 3U;
  center = geometry == QUAD || geometry == HEX ? std::array<long double,3>{-0.25L,-0.25L,-0.25L} : std::array<long double,3>{0.25L,0.25L,geometry == WEDGE ? 0.0L : 0.125L};
  radius = geometry == TET ? 0.0625L : 0.125L;
  nodal.assign(basis->_nc,0);
  coordinates.assign(dimension,nodal);
  for(unsigned n = 0; n < nodal.size(); ++n) {
    const double* point = basis->GetXcoarse(n);
    long double value = -radius*radius;
    for(unsigned axis = 0; axis < dimension; ++axis) value += (point[axis]-center[axis])*(point[axis]-center[axis]);
    nodal[n] = value;
    coordinates[0][n] = point[0]+0.05*point[0]*point[1];
    coordinates[1][n] = point[1]+0.04*point[0]*point[1];
    if(dimension == 3) coordinates[2][n] = point[2]+0.03*point[0]*point[2]+0.02*point[1]*point[2];
  }
}

long double CenteredBallMoment(unsigned dimension, long double radius, unsigned i, unsigned j, unsigned k) {
  if(i%2 != 0 || j%2 != 0 || k%2 != 0) return 0;
  long double value = std::pow(radius,dimension+i+j+k)*std::tgamma((i+1)/2.0L)*std::tgamma((j+1)/2.0L)/std::tgamma(1+(dimension+i+j+k)/2.0L);
  if(dimension == 3) value *= std::tgamma((k+1)/2.0L);
  return value;
}

void ReferenceBallWeightMoments(GeomElType geometry, unsigned degree, unsigned mode, const std::array<long double,3>& center, long double radius, std::vector<long double>& positive, std::vector<long double>& surface) {
  const unsigned dimension = geometry == QUAD || geometry == TRI ? 2U : 3U, count = dimension == 2 ? (degree+1)*(degree+2)/2 : (degree+1)*(degree+2)*(degree+3)/6;
  const long double pi = std::acos(-1.0L);
  positive.assign(count,0);
  surface.assign(count,0);
  const std::array<long double,3> unitCenter = geometry == QUAD || geometry == HEX ? std::array<long double,3>{(center[0]+1)/2,(center[1]+1)/2,(center[2]+1)/2} : (geometry == WEDGE ? std::array<long double,3>{center[0],center[1],(center[2]+1)/2} : center);
  for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= (dimension == 3 ? degree-i-j : 0U); ++k) {
    long double full = geometry == QUAD ? 4.0L/((i+1)*(j+1)) : (geometry == HEX ? 8.0L/((i+1)*(j+1)*(k+1)) : (geometry == TRI ? FELevelSetMoments::ReferenceTriangleMoment(i,j) : (geometry == WEDGE ? 2*FELevelSetMoments::ReferenceTriangleMoment(i,j)/(k+1) : std::tgamma(i+1.0L)*std::tgamma(j+1.0L)*std::tgamma(k+1.0L)/std::tgamma(i+j+k+4.0L))));
    long double ball = 0;
    for(unsigned a = 0; a <= i; ++a) for(unsigned b = 0; b <= j; ++b) for(unsigned c = 0; c <= k; ++c) { const long double scale = std::pow(0.5L,(geometry == QUAD || geometry == HEX ? a+b+c : (geometry == WEDGE ? c : 0U))); ball += FELevelSetMoments::ReferenceBinomial(i,a)*FELevelSetMoments::ReferenceBinomial(j,b)*FELevelSetMoments::ReferenceBinomial(k,c)*std::pow(unitCenter[0],i-a)*std::pow(unitCenter[1],j-b)*std::pow(unitCenter[2],k-c)*scale*CenteredBallMoment(dimension,radius,a,b,c); }
    positive[FELevelSetMoments::UnitMomentIndex(dimension,i,j,k)] = full-ball;
  }
  const unsigned angles = 128;
  const std::vector<long double>& nodes = boost::math::quadrature::gauss<long double,60>::abscissa();
  const std::vector<long double>& weights = boost::math::quadrature::gauss<long double,60>::weights();
  std::array<std::vector<long double>,3> powers;
  const unsigned latitudes = dimension == 2 ? 1U : 60U;
  for(unsigned latitude = 0; latitude < latitudes; ++latitude) {
    const long double z = dimension == 2 ? 0 : (latitude%2 == 0 ? nodes[latitude/2] : -nodes[latitude/2]), wz = dimension == 2 ? 1 : weights[latitude/2], radial = std::sqrt(1-z*z);
    for(unsigned angle = 0; angle < angles; ++angle) {
      const long double theta = 2*pi*angle/angles;
      const std::array<long double,3> normal = {radial*std::cos(theta),radial*std::sin(theta),z}, point = {center[0]+radius*normal[0],center[1]+radius*normal[1],center[2]+radius*normal[2]};
      const long double j00 = 1+0.05L*point[1], j01 = 0.05L*point[0], j10 = 0.04L*point[1], j11 = 1+0.04L*point[0], j20 = dimension == 3 ? 0.03L*point[2] : 0, j21 = dimension == 3 ? 0.02L*point[2] : 0, j22 = dimension == 3 ? 1+0.03L*point[0]+0.02L*point[1] : 1;
      const std::array<long double,3> h = {j22*(j11*normal[0]-j10*normal[1])+(j10*j21-j11*j20)*normal[2],j22*(-j01*normal[0]+j00*normal[1])+(j01*j20-j00*j21)*normal[2],(j00*j11-j01*j10)*normal[2]};
      const long double norm = std::sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]);
      const long double density = mode == 0 ? (dimension == 2 ? 0.5L : radius/2) : (mode == 1 ? std::pow(radius,dimension-1) : (mode == 2 ? std::pow(radius,dimension-1)*norm : (mode == 3 ? std::pow(radius,static_cast<int>(dimension)-3)/(4*norm) : std::pow(radius,dimension-2)*h[mode-4]/(2*norm))));
      const long double weight = 2*pi*wz*density/angles;
      WeightBasisPowers(geometry,point,degree,powers);
      for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) for(unsigned k = 0; k <= (dimension == 3 ? degree-i-j : 0U); ++k) surface[FELevelSetMoments::UnitMomentIndex(dimension,i,j,k)] += weight*powers[0][i]*powers[1][j]*powers[2][k];
    }
  }
}

void TestBallWeightIntegration(unsigned dimension) {
  const std::vector<GeomElType> elements = dimension == 2 ? std::vector<GeomElType>{QUAD,TRI} : std::vector<GeomElType>{HEX,TET,WEDGE};
  for(GeomElType geometry : elements) for(unsigned degree : {2U,5U,10U}) {
    if(dimension == 3 && geometry != HEX && degree > 5) continue;
    currentCheck = std::string(GeometryName(geometry))+" analytic ball weights degree "+std::to_string(degree);
    CutFemWeight_LS<long double,boost::multiprecision::cpp_bin_float_quad> referenceFitter(geometry,degree,"legendre");
    CutFemWeight_LS<long double,boost::multiprecision::cpp_bin_float_quad> prepared(geometry,degree,"legendre");
    std::vector<double> nodal;
    std::vector<std::vector<double>> coordinates;
    std::array<long double,3> center;
    long double radius = 0, rawError = 0, newError = 0, fittingFloor = 0;
    MakeWeightBall(geometry,nodal,coordinates,center,radius);
    for(unsigned mode = 0; mode < 4; ++mode) {
      std::vector<long double> newH, newD, exactH, exactD, rawH, rawD, floorH, floorD;
      ReferenceBallWeightMoments(geometry,degree,mode,center,radius,exactH,exactD);
      GetWeightPair(prepared,mode,nodal,coordinates,newH,newD);
      prepared.GetLastIntegratedMoments(rawH,rawD);
      const long double measure = geometry == QUAD ? 4 : (geometry == HEX ? 8 : (geometry == WEDGE ? 2 : 1));
      std::vector<long double> unitH = exactH, unitD = exactD;
      for(unsigned n = 0; n < rawH.size(); ++n) { unitH[n] /= measure; unitD[n] /= measure; rawError = std::max(rawError,std::max(std::abs(measure*rawH[n]-exactH[n])/(1+std::abs(exactH[n])),std::abs(measure*rawD[n]-exactD[n])/(1+std::abs(exactD[n])))); CheckNear(measure*rawH[n],exactH[n],2e-9L,"raw positive ball moment"); CheckNear(measure*rawD[n],exactD[n],2e-9L,"raw surface ball moment"); }
      referenceFitter.GetWeightsFromExternalMoments(unitH,unitD,floorH,floorD);
      const std::vector<long double> newHM = FittedWeightMoments(prepared,geometry,dimension,newH,degree), newDM = FittedWeightMoments(prepared,geometry,dimension,newD,degree), floorHM = FittedWeightMoments(referenceFitter,geometry,dimension,floorH,degree), floorDM = FittedWeightMoments(referenceFitter,geometry,dimension,floorD,degree);
      for(unsigned n = 0; n < rawH.size(); ++n) {
        const long double hf = std::abs(floorHM[n]-exactH[n])/(1+std::abs(exactH[n])), df = std::abs(floorDM[n]-exactD[n])/(1+std::abs(exactD[n]));
        fittingFloor = std::max(fittingFloor,std::max(hf,df));
        newError = std::max(newError,std::max(std::abs(newHM[n]-exactH[n])/(1+std::abs(exactH[n])),std::abs(newDM[n]-exactD[n])/(1+std::abs(exactD[n]))));
        CheckNear(newHM[n],exactH[n],2*hf+1e-8L,"fitted positive ball moment");
        CheckNear(newDM[n],exactD[n],2*df+1e-8L,"fitted surface ball moment");
      }
    }
    std::cout << std::setprecision(8) << "weight_check " << GeometryName(geometry) << " q=" << degree << " prepared_error=" << newError << " raw_error=" << rawError << " fitting_floor=" << fittingFloor << std::endl;
    Check(prepared.GetPreparedLevelSetData().GetWorkCounts()[0] == 1,"all surface modes reuse the parent polynomial reconstruction");
  }
}

void BenchmarkElementWeights(unsigned dimension) {
  const std::vector<GeomElType> elements = dimension == 2 ? std::vector<GeomElType>{QUAD,TRI} : std::vector<GeomElType>{HEX,TET,WEDGE};
  for(GeomElType geometry : elements) for(unsigned degree : {2U,5U}) {
    std::vector<double> nodal;
    std::vector<std::vector<double>> coordinates;
    std::array<long double,3> center;
    long double radius = 0;
    MakeWeightBall(geometry,nodal,coordinates,center,radius);
    std::vector<double> second = nodal;
    const std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    for(unsigned n = 0; n < second.size(); ++n) second[n] += 0.03125L*(basis->GetXcoarse(n)[0]-center[0])+0.000244140625L;
    std::array<CutFemLevelSetData,2> data = {{CutFemLevelSetData(geometry,nodal),CutFemLevelSetData(geometry,second)}};
    CutFemWeight_LS<long double,long double> nodalRule(geometry,degree,"legendre");
    CutFemWeight_LS<long double,long double> prepared(geometry,degree,"legendre");
    for(unsigned mode = 0; mode < (dimension == 3 ? 5U : 4U); ++mode) {
      std::array<double,5> oldTime, newTime;
      long double oldChecksum = 0, newChecksum = 0;
      const unsigned iterations = dimension == 2 ? 16U : 2U;
      for(unsigned repeat = 0; repeat < 5; ++repeat) { oldTime[repeat] = TimeLineWeights(nodalRule,mode,data,coordinates,iterations,oldChecksum,true); newTime[repeat] = TimeLineWeights(prepared,mode,data,coordinates,iterations,newChecksum); }
      std::sort(oldTime.begin(),oldTime.end());
      std::sort(newTime.begin(),newTime.end());
      CheckNear(newChecksum,oldChecksum,2e-7L,"element benchmark checksums");
      std::cout << std::setprecision(8) << "weight_speed " << GeometryName(geometry) << " q=" << degree << " mode=" << mode << " nodal_us=" << oldTime[2] << " prepared_us=" << newTime[2] << " nodal_over_prepared=" << oldTime[2]/newTime[2] << std::endl;
    }
  }
}

void TestHighDegreePreparedMoments() {
  for(GeomElType geometry : {HEX,TET,WEDGE}) {
    currentCheck = std::string(GeometryName(geometry))+" degree 13 prepared moments";
    std::vector<double> nodal;
    std::vector<std::vector<double>> coordinates;
    std::array<long double,3> center;
    long double radius = 0, maximumError = 0;
    MakeWeightBall(geometry,nodal,coordinates,center,radius);
    const CutFemLevelSetData data(geometry,nodal);
    FELevelSetMoments::PreparedMomentWorkspace workspace;
    for(unsigned mode = 0; mode < 4; ++mode) {
      std::vector<long double> values, errors, exactH, exactD;
      FELevelSetMoments::IntegrateSpatial(data,13,1e-10,1e-10,mode != 0,mode >= 2 ? &coordinates : nullptr,0,1,FELevelSetMoments::OuterQuadratureRules(),mode == 3,-1,true,workspace,values,errors);
      if(mode == 1 && geometry == WEDGE) Check(workspace.surface.gradient == workspace.referenceMetric,"WEDGE reference bookkeeping uses the selected metric exactly once");
      if(mode >= 2) Check(workspace.surface.gradient == workspace.metric,"physical bookkeeping retains the selected cofactor metric");
      ReferenceBallWeightMoments(geometry,13,mode,center,radius,exactH,exactD);
      const long double measure = geometry == HEX ? 8 : (geometry == WEDGE ? 2 : 1);
      const long double surfaceFactor = mode == 0 ? 1 : (mode >= 2 ? (mode == 3 ? measure : 1/measure) : (geometry == HEX ? 0.5L : 1.0L));
      for(unsigned n = 0; n < exactH.size(); ++n) { CheckNear(measure*values[n],exactH[n],2e-9L,"high-degree raw volume moment"); CheckNear(measure*surfaceFactor*values[exactH.size()+n],exactD[n],2e-9L,"high-degree raw surface moment"); maximumError = std::max(maximumError,std::max(std::abs(measure*values[n]-exactH[n])/(1+std::abs(exactH[n])),std::abs(measure*surfaceFactor*values[exactH.size()+n]-exactD[n])/(1+std::abs(exactD[n])))); }
    }
    std::cout << "weight_raw " << GeometryName(geometry) << " q=13 error=" << maximumError << std::endl;
  }
}

void TestElementWeightOptions(unsigned dimension) {
  const std::vector<GeomElType> elements = dimension == 2 ? std::vector<GeomElType>{QUAD,TRI} : std::vector<GeomElType>{HEX,TET,WEDGE};
  for(GeomElType geometry : elements) {
    currentCheck = std::string(GeometryName(geometry))+" prepared options and enriched basis";
    CutFemWeight_LS<long double,boost::multiprecision::cpp_bin_float_quad> nodalRule(geometry,2,"legendre");
    CutFemWeight_LS<long double,boost::multiprecision::cpp_bin_float_quad> prepared(geometry,2,"legendre");
    const std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
    std::vector<double> nodal;
    std::vector<std::vector<double>> coordinates;
    std::array<long double,3> center;
    long double radius = 0;
    MakeWeightBall(geometry,nodal,coordinates,center,radius);
    for(unsigned field = 0; field < 4; ++field) {
      for(unsigned n = 0; n < nodal.size(); ++n) { const double* x = basis->GetXcoarse(n); nodal[n] = field == 0 ? 1 : (field == 1 ? -1 : x[0]-0.25+0.02*x[0]*x[0]+0.03*x[1]+(dimension == 3 ? 0.02*x[2]*x[1] : 0)+(field == 3 ? 0.0001*std::sin(3.0*n+0.25) : 0)); }
      CutFemLevelSetData data(geometry,nodal);
      const unsigned directions = field < 2 ? 0 : (geometry == QUAD ? 2U : (geometry == TRI ? 3U : (geometry == TET ? 4U : 6U)));
      for(unsigned direction = 0; direction <= directions; ++direction) {
        nodalRule.SetIntegrationDirection(direction);
        prepared.SetIntegrationDirection(direction);
        for(unsigned mode = 0; mode < 4; ++mode) {
          std::vector<long double> h, d, oldH, oldD;
          GetWeightPair(nodalRule,mode,nodal,coordinates,oldH,oldD);
          GetFixtureWeightPair(prepared,mode,data,coordinates,h,d);
          const std::vector<long double> oldHM = FittedWeightMoments(nodalRule,geometry,dimension,oldH,2), oldDM = FittedWeightMoments(nodalRule,geometry,dimension,oldD,2), newHM = FittedWeightMoments(prepared,geometry,dimension,h,2), newDM = FittedWeightMoments(prepared,geometry,dimension,d,2);
          for(unsigned n = 0; n < newHM.size(); ++n) { CheckNear(newHM[n],oldHM[n],2e-8L,"directional volume moment agreement"); CheckNear(newDM[n],oldDM[n],2e-8L,"directional surface moment agreement"); }
        }
      }
    }
    prepared.SetIntegrationDirection(0);
    nodalRule.SetIntegrationDirection(0);
    MakeWeightBall(geometry,nodal,coordinates,center,radius);
    CutFemLevelSetData ball(geometry,nodal);
    for(unsigned rule : {7U,31U}) {
      FELevelSetMoments::OuterQuadratureRules rules;
      rules.regular2D = rules.tangent2D = rules.regular3D = rules.tangent3D = rule;
      prepared.SetOuterQuadratureRules(rules);
      prepared.SetDirectionValidation(rule == 31 ? 2U : 0U);
      std::vector<long double> h, d, exactH, exactD;
      GetFixtureWeightPair(prepared,2,ball,coordinates,h,d);
      ReferenceBallWeightMoments(geometry,2,2,center,radius,exactH,exactD);
      const std::vector<long double> fitted = FittedWeightMoments(prepared,geometry,dimension,d,2);
      for(unsigned n = 0; n < fitted.size(); ++n) CheckNear(fitted[n],exactD[n],3e-9L,"configured GK surface moments");
    }
    prepared.SetDirectionValidation(1);
    const unsigned long long preparations = prepared.GetPreparedStatistics().geometryPreparations;
    std::vector<long double> h, d;
    prepared.ClearCache();
    GetFixtureWeightPair(prepared,3,ball,coordinates,h,d);
    Check(prepared.GetPreparedStatistics().geometryPreparations == preparations,"geometry preparation survives rules, measure, and final-weight cache changes");
    coordinates[0][0] += 1e-5;
    GetFixtureWeightPair(prepared,2,ball,coordinates,h,d);
    Check(prepared.GetPreparedStatistics().geometryPreparations == preparations+1,"geometry change invalidates the physical metric");
    if(dimension == 3) {
      MakeWeightBall(geometry,nodal,coordinates,center,radius);
      CutFemLevelSetData sphere(geometry,nodal);
      std::vector<std::vector<long double>> normal, oldNormal;
      nodalRule.GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,h,oldNormal,1e-10,1e-10);
      prepared.GetNormalSurfaceWeightsFromLevelSetData(sphere,coordinates,h,normal,1e-10,1e-10);
      if(geometry == TET) { prepared.SetIntegrationDirection(1); prepared.GetNormalSurfaceWeightsFromLevelSetData(sphere,coordinates,h,normal,1e-10,1e-10); }
      for(unsigned component = 0; component < 3; ++component) {
        std::vector<long double> exactH, exactD;
        ReferenceBallWeightMoments(geometry,2,4+component,center,radius,exactH,exactD);
        const std::vector<long double> fitted = FittedWeightMoments(prepared,geometry,dimension,normal[component],2), baseline = FittedWeightMoments(nodalRule,geometry,dimension,oldNormal[component],2);
        for(unsigned n = 0; n < fitted.size(); ++n) { CheckNear(fitted[n],exactD[n],3e-9L,"analytic signed-normal moment"); CheckNear(fitted[n],baseline[n],3e-9L,"nodalRule/prepared signed-normal moment"); }
      }
    }
    std::cout << "weight_options " << GeometryName(geometry) << " passed" << std::endl;
  }
}

void TestSmallPlanarWeights() {
  currentCheck = "native tiny planar cuts";
  CutFemWeight_LS<long double,boost::multiprecision::cpp_bin_float_quad> quadrature(QUAD,2,"legendre");
  for(long double radius : {1e-4L,1e-8L,1e-12L}) {
    std::vector<double> nodal(9);
    const std::unique_ptr<femus::basis> basis = MakeBasis(QUAD);
    for(unsigned n = 0; n < nodal.size(); ++n) { const double* x = basis->GetXcoarse(n); nodal[n] = x[0]*x[0]+x[1]*x[1]-radius*radius; }
    const CutFemLevelSetData data(QUAD,nodal);
    std::vector<long double> h, d;
    quadrature.GetWeightsFromLevelSetData(data,h,d,1e-10,1e-10,false);
    quadrature.GetLastIntegratedMoments(h,d);
    CheckNear(4*d[0],std::acos(-1.0L),2e-8L,"tiny-circle delta measure");
    CheckNear(4*h[0],4-std::acos(-1.0L)*radius*radius,2e-10L,"tiny-circle volume");
  }
}

void TestHexSaddleProjection() {
  currentCheck = "HEX27 tangential weights with saddle points in all three projections";
  const std::vector<double> nodal = {-0.00013635390802039783,-0.00064460436066085297,-0.0011683214673039054,0.00010500267777512871,-0.00030361125340942513,0.00069214235526596853,0.0015287449843700938,1.1934081100784282e-05,-8.0977860379127177e-05,-0.00096900392774501663,-0.00063573883900853751,0.00014999234775047162,-0.00015277686883369306,0.0013024057266662083,0.0011367187800165118,-0.00038773668000858705,-0.00065389955078982878,0.00024418757758534298,-8.51184712337801e-05,0.00017688451792260501,6.7841713189283255e-05,8.7424330826183446e-05,0.00015990330642173079,-0.00016207556802104204,-0.00034120090070308828,0.00048832599574893606,0.00021377421661503265};
  const std::unique_ptr<femus::basis> basis = MakeBasis(HEX);
  std::vector<std::vector<double>> coordinates(3,std::vector<double>(nodal.size()));
  for(unsigned n = 0; n < nodal.size(); ++n) for(unsigned axis = 0; axis < 3; ++axis) coordinates[axis][n] = basis->GetXcoarse(n)[axis]/1280;
  CutFemWeight_LS<long double,long double> quadrature(HEX,5,"legendre");
  std::vector<long double> positive, surface;
  const unsigned long long warnings = FELevelSetMoments::QuadratureWarningCount();
  const unsigned long long logPanels = FELevelSetMoments::OuterIntervalQuadrature::Statistics(3).logarithmicPanels;
  quadrature.GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,surface,1e-8,1e-8);
  Check(FELevelSetMoments::OuterIntervalQuadrature::Statistics(3).logarithmicPanels > logPanels,"interior saddles use logarithmic quadrature");
  Check(FELevelSetMoments::QuadratureWarningCount() == warnings,"saddle projection meets the requested tolerance");
  const std::vector<long double> h = FittedWeightMoments(quadrature,HEX,3,positive,5), d = FittedWeightMoments(quadrature,HEX,3,surface,5);
  const unsigned xy = FELevelSetMoments::UnitMomentIndex(3,1,1,0), z5 = FELevelSetMoments::UnitMomentIndex(3,0,0,5);
  // Independent fixed-order graph quadrature on element 69612, in unit coordinates.
  CheckNear(h[0]/8,0.58600081172474183L,1e-8L,"positive volume");
  CheckNear(d[0]/8,6333983265621.4111L,1e-8L,"tangential mass");
  CheckNear(h[xy]/8,0.15367998129365718L,1e-8L,"mixed positive moment");
  CheckNear(d[xy]/8,677258555576.12524L,1e-8L,"mixed tangential moment");
  CheckNear(h[z5]/8,0.12996399549813537L,1e-8L,"degree-five positive moment");
  CheckNear(d[z5]/8,199241442653.80463L,1e-8L,"degree-five tangential moment");

  std::vector<std::vector<long double>> normals;
  quadrature.GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,normals,1e-8,1e-8);
  const std::array<long double,3> expectedNormal = {168.65300890836528L,271.57008413395181L,534.28540021517574L};
  for(unsigned axis = 0; axis < 3; ++axis) CheckNear(FittedWeightMoments(quadrature,HEX,3,normals[axis],5)[0]/8,expectedNormal[axis],1e-8L,"signed normal mass");
  Check(FELevelSetMoments::QuadratureWarningCount() == warnings,"signed saddle projections meet the requested tolerance");

  currentCheck = "HEX27 whole-element outer direction";
  const std::vector<double> twoSaddles = {0.0011137319967982392,-0.00045001059519084603,-0.00097515682293665851,-0.00044996036212267804,0.0010944453954665851,-0.00046172092017875941,-0.00094752956408125589,-0.00046126626498387268,0.00043708252339941863,-0.00047416194807737297,-0.00047334447336950652,0.00043753463971559842,0.0004744843356212895,-0.00043373395702033028,-0.00043289383593618681,0.00047494675803111722,0.00094471009657230037,-0.00072743005977729839,-0.0009410772621390492,-0.00072732978923546558,0.00050921245406149401,-0.00045414287132016754,-0.00045327738702507608,0.00050912104630967091,0.00060012500935635962,0.00020820008447767745,0.00014636136428263074};
  const CutFemLevelSetData data(HEX,twoSaddles);
  FELevelSetMoments::PreparedMomentWorkspace workspace;
  std::vector<long double> values, errors;
  FELevelSetMoments::IntegrateSpatial(data,5,1e-8/64,1e-8,true,&coordinates,0,0,FELevelSetMoments::OuterQuadratureRules(),true,-1,true,workspace,values,errors);
  Check(workspace.surface.outerCounts == std::array<unsigned,3>{0,0,1},"element 23611 chooses z on its first attempt");
  Check(FELevelSetMoments::QuadratureWarningCount() == warnings,"selected outer axis meets tolerance");

  currentCheck = "HEX27 volume shared between normal and tangential weights";
  const std::vector<double> cachedVolumeField = {-0.0005452486121344316,-0.00016548584445942978,-0.0030631252117336158,-0.0037788307366267565,0.000394281025773052,-0.0010904698364475848,-0.00056801995789561559,0.000168276568626631,-0.00015360608204576716,-0.0012665497605746428,-0.00364177636619117,-0.001930954284192552,-0.00047884748321712855,-0.00073595872719249967,5.7853037213847363e-05,0.00017413329937618788,6.0737625104915975e-06,-0.00060871392795185022,-0.00098422816453096781,-0.0028786082169287614,9.0673555847626557e-05,-0.00021841635783366822,-0.0023006657090184823,-0.00092636562858408817,-0.001657415319868125,-7.8944039770497695e-05,-0.00059925034101617012};
  quadrature.GetNormalSurfaceWeightsFromNodalGeometry(cachedVolumeField,coordinates,positive,normals,1e-8,1e-8);
  const unsigned xxz = FELevelSetMoments::UnitMomentIndex(3,2,0,1);
  CheckNear(FittedWeightMoments(quadrature,HEX,3,positive,5)[xxz]/8,0.0032484132034037356L,1e-8L,"signed-normal volume moment");
  quadrature.GetTangentialSurfaceWeightsFromNodalGeometry(cachedVolumeField,coordinates,positive,surface,1e-8,1e-8);
  CheckNear(FittedWeightMoments(quadrature,HEX,3,positive,5)[xxz]/8,0.0032484132034037356L,1e-8L,"volume moment after a signed-normal request");
  Check(FELevelSetMoments::QuadratureWarningCount() == warnings,"reused positive weights meet tolerance");
}

void TestCenteredNormalWeights() {
  currentCheck = "centered sphere signed-normal parity";
  CutFemWeight_LS<long double,boost::multiprecision::cpp_bin_float_quad> quadrature(HEX,2,"legendre");
  const std::unique_ptr<femus::basis> basis = MakeBasis(HEX);
  const long double radius = 0.25L, expected = std::acos(-1.0L)*radius*radius/3;
  std::vector<double> nodal(basis->_nc);
  std::vector<std::vector<double>> coordinates(3,nodal);
  for(unsigned n = 0; n < nodal.size(); ++n) { const double* x = basis->GetXcoarse(n); nodal[n] = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]-radius*radius; for(unsigned axis = 0; axis < 3; ++axis) coordinates[axis][n] = x[axis]; }
  for(unsigned sign = 0; sign < 2; ++sign) {
    if(sign != 0) for(double& value : nodal) value = -value;
    std::vector<long double> h;
    std::vector<std::vector<long double>> normal;
    quadrature.GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,h,normal,1e-10,1e-10);
    for(unsigned component = 0; component < 3; ++component) {
      const std::vector<long double> moments = FittedWeightMoments(quadrature,HEX,3,normal[component],2);
      CheckNear(moments[0],0,2e-12L,"centered normal has zero mean");
      for(unsigned axis = 0; axis < 3; ++axis) { unsigned powers[3] = {0,0,0}; powers[axis] = 1; CheckNear(moments[FELevelSetMoments::UnitMomentIndex(3,powers[0],powers[1],powers[2])],axis == component ? (sign == 0 ? expected : -expected) : 0,2e-12L,"normal parity and orientation"); }
    }
  }
}

void TestReviewedChartContracts(std::mt19937& random) {
  currentCheck = "reviewed chart contracts";
  for(GeomElType geometry : {TRI,TET,WEDGE}) {
    const unsigned dimension = Dimension(geometry), count = dimension == 2 ? 9U : 27U;
    std::vector<long double> power(count,0);
    power[dimension == 2 ? 8U : (geometry == TET ? 24U : 21U)] = 1;
    CutFemLevelSetData imported, chart;
    imported.ResetPowerCoefficients(geometry,power);
    imported.PermuteSimplexVertices({0,2,1,3},chart);
    for(unsigned sample = 0; sample < 30; ++sample) {
      const std::vector<long double> point = RandomPoint(geometry,random), parent = chart.MapToParent(point);
      CheckNear(chart.Evaluate(point),imported.Evaluate(parent),3e-17L,"imported chart preserves every power coefficient");
      CheckNear(chart.EvaluateBernstein(point),chart.Evaluate(point),3e-17L,"imported chart Bernstein and powers agree");
      const std::vector<long double> expected = imported.GetRootGradient(parent), gradient = chart.GetRootGradient(point);
      for(unsigned axis = 0; axis < dimension; ++axis) CheckNear(gradient[axis],expected[axis],3e-17L,"imported chart retains original gradient components");
    }
    bool rejected = false;
    try { imported.PermuteSimplexVertices({1,0,2,3},chart); } catch(const std::domain_error&) { rejected = true; }
    Check(rejected && !chart.IsInitialized(),"unsupported chart degree is rejected without a partial polynomial");
    for(unsigned n = 0; n < count; ++n) {
      const unsigned i = dimension == 2 ? n/3 : n/9, j = dimension == 2 ? n%3 : (n/3)%3, k = dimension == 2 ? 0U : n%3;
      power[n] = i+j+k <= 2 ? static_cast<long double>(static_cast<int>(n%5)-2)/7 : 0;
    }
    imported.ResetPowerCoefficients(geometry,power);
    imported.PermuteSimplexVertices({1,0,2,3},chart);
    for(unsigned sample = 0; sample < 30; ++sample) { const std::vector<long double> point = RandomPoint(geometry,random); CheckNear(chart.Evaluate(point),imported.Evaluate(chart.MapToParent(point)),3e-17L,"imported affine composition preserves complete quadratics"); }
    for(unsigned test = 0; test < 20; ++test) {
      std::unique_ptr<femus::basis> basis = MakeBasis(geometry);
      std::vector<double> nodal(basis->_nc);
      for(unsigned node = 0; node < nodal.size(); ++node) nodal[node] = std::sin(0.37*(1+node+test));
      CutFemLevelSetData native(geometry,nodal), first, second;
      native.PermuteAxes({1,0,2},first);
      first.PermuteSimplexVertices({1,0,2,3},second);
      for(unsigned sample = 0; sample < 10; ++sample) { const std::vector<long double> point = RandomPoint(geometry,random); CheckNear(second.Evaluate(point),native.Evaluate(second.MapToRoot(point)),3e-13L,"repeated native charts retain FEM-space closure"); CheckNear(second.EvaluateBernstein(point),second.Evaluate(point),3e-13L,"repeated native chart representations agree"); }
      if(dimension == 3) {
        native.Restrict(2,0.17L,first);
        first.PermuteSimplexVertices({1,0,2,3},second);
        for(unsigned sample = 0; sample < 10; ++sample) { const std::vector<long double> point = RandomPoint(second.GetGeometry(),random); CheckNear(second.Evaluate(point),native.Evaluate(second.MapToRoot(point)),3e-13L,"native slices retain chart closure"); }
      }
    }
  }
  CutFemLevelSetData data(QUAD,MakeField(QUAD,2)), slice;
  std::array<long double,5> shifted = {19,19,19,19,19};
  data.PermuteAxes({1,0,2},slice);
  bool rejected = false;
  try { slice.GetLocalShiftedFiberDiscriminant(0,std::numeric_limits<long double>::quiet_NaN(),shifted); } catch(const std::invalid_argument&) { rejected = true; }
  Check(rejected && shifted == std::array<long double,5>{{19,19,19,19,19}},"invalid local shifted discriminants preserve outputs");
  std::vector<long double> roots;
  rejected = false;
  try { data.GetCompensatedDiscriminantRoots(0,0.8L,0.2L,roots); } catch(const std::invalid_argument&) { rejected = true; }
  Check(rejected,"compensated root bounds are ordered");
}

void TestReviewedChartPrecision() {
  currentCheck = "chart source precision";
  for(GeomElType geometry : {QUAD,HEX,TET,WEDGE}) for(double minimum : {-1e-24,-1e-80,-1e-200}) {
    if((geometry == TET || geometry == WEDGE) && minimum != -1e-24) continue;
    std::vector<double> nodal(MakeBasis(geometry)->_nc,1);
    nodal.back() = minimum;
    CutFemLevelSetData parent(geometry,nodal), chart, direct, transformed, planar;
    if(geometry == TET || geometry == WEDGE) parent.PermuteSimplexVertices({1,0,2,3},chart);
    else parent.PermuteAxes({1,0,2},chart);
    Check(chart.HasNodalPolynomialSource(),"a native chart retains its own nodal source");
    const long double height = geometry == TET ? 0.25L : 0.5L, fixed = geometry == TET || geometry == WEDGE ? 1.0L/3 : 0.5L;
    if(Dimension(geometry) == 3) { parent.Restrict(2,height,planar,false); planar.Restrict(0,fixed,direct,false); chart.Restrict(2,height,planar,false); planar.Restrict(0,fixed,transformed,false); }
    else { parent.Restrict(0,fixed,direct,false); chart.Restrict(0,fixed,transformed,false); }
    const CutFemLevelSetData::LineAnalysis& expected = direct.GetLineAnalysis();
    const CutFemLevelSetData::LineAnalysis& actual = transformed.GetLineAnalysis();
    Check(expected.rootCount == 2 && actual.rootCount == 2,"chart permutations preserve both small-cut roots");
    for(unsigned n = 0; n < 2; ++n) { Check(actual.roots[n].multiplicity == 1,"chart permutations retain simple roots"); CheckNear(transformed.GetLineDeltaFactor(actual.roots[n])/direct.GetLineDeltaFactor(expected.roots[n]),1,3e-15L,"chart root densities preserve retained nodal precision"); }
    if(geometry == HEX) {
      parent.Restrict(2,0.5L,planar,false);
      planar.PermuteAxes({1,0,2},chart);
      chart.Restrict(0,0.5L,transformed,false);
      Check(transformed.GetLineAnalysis().rootCount == 2,"a chart of a slice preserves sub-ULP nodal information");
      for(unsigned n = 0; n < 2; ++n) CheckNear(transformed.GetLineDeltaFactor(transformed.GetLineAnalysis().roots[n])/direct.GetLineDeltaFactor(expected.roots[n]),1,3e-15L,"slice-chart densities match direct restrictions");
    }
    const CutFemLevelSetData saved = chart;
    parent.Reset(geometry,std::vector<double>(nodal.size(),2));
    if(geometry == TET || geometry == WEDGE) parent.PermuteSimplexVertices({1,0,2,3},chart);
    else parent.PermuteAxes({1,0,2},chart);
    if(saved.GetDimension() == 3) { saved.Restrict(2,height,planar,false); planar.Restrict(0,fixed,transformed,false); }
    else saved.Restrict(0,fixed,transformed,false);
    Check(transformed.GetLineAnalysis().rootCount == 2,"chart source snapshots survive output reuse");
    chart.Reset(geometry,std::vector<double>(nodal.size(),3));
    Check(saved.HasNodalPolynomialSource(),"chart sources survive parent and output reset");
  }
}

void TestReviewedWeightState() {
  currentCheck = "reviewed weight state";
  CutFemWeight_LS<long double,boost::multiprecision::cpp_bin_float_quad> rule(LINE,2,"legendre");
  std::vector<long double> h,d,originalH,originalD,afterH,afterD;
  rule.GetWeightsFromNodalValues({-0.3,0.7,0.2},h,d);
  rule.GetLastIntegratedMoments(originalH,originalD);
  std::vector<long double> invalid = originalD, changed = originalH;
  changed[0] += 1;
  invalid.back() = std::numeric_limits<long double>::quiet_NaN();
  bool rejected = false;
  try { rule.GetWeightsFromExternalMoments(changed,invalid,h,d); } catch(const std::invalid_argument&) { rejected = true; }
  Check(rejected,"invalid external moments are rejected before state mutation");
  rule.GetLastIntegratedMoments(afterH,afterD);
  Check(afterH == originalH && afterD == originalD,"invalid external fitting preserves integrated moments");
  rule.GetWeightsFromExternalMoments(originalH,originalD,h,d);
  rejected = false;
  try { rule.GetLastIntegratedMoments(afterH,afterD); } catch(const std::logic_error&) { rejected = true; }
  Check(rejected,"external fitting cannot masquerade as completed integration");
  rule.GetWeightsFromNodalValues({-0.3,0.7,0.2},h,d);
  rule.GetLastIntegratedMoments(afterH,afterD);
  Check(afterH == originalH && afterD == originalD,"integration state is restored after external fitting");
  rule.GetWeightsFromNodalValues({-0.3,0.7,0.2},h,d,1e-12,1e-12,true);
  rule.GetLastIntegratedMoments(originalH,originalD);
  rule.GetWeightsFromNodalValues({-0.3,0.7,0.2},h,d);
  rule.GetLastIntegratedMoments(afterH,afterD);
  Check(afterH == originalH && afterD == originalD,"cache hits retain the last actually integrated pair");
}

int main(int argc, char** argv) {
  try {
    unsigned cases = 64, points = 32;
    bool benchmark = false;
    bool quadratureBenchmark = false;
    bool surfaceBenchmark = false;
    bool dataBenchmark = false;
    bool algebraBenchmark = false;
    bool weightTests = false;
    bool weightBenchmark = false;
    unsigned weightDimension = 0;
    for(int i = 1; i < argc; ++i) {
      const std::string option = argv[i];
      if((option == "--cases" || option == "--points") && i + 1 < argc) {
        const std::string text = argv[++i];
        if(text.empty() || text.find_first_not_of("0123456789") != std::string::npos) throw std::invalid_argument("Test counts must be positive integers.");
        const unsigned long value = std::stoul(text);
        if(value == 0 || value > 10000) throw std::invalid_argument("Test counts must lie between 1 and 10000.");
        if(option == "--cases") cases = static_cast<unsigned>(value);
        else points = static_cast<unsigned>(value);
      }
      else if(option == "--benchmark") benchmark = true;
      else if(option == "--benchmark-quadrature") quadratureBenchmark = true;
      else if(option == "--benchmark-surface") surfaceBenchmark = true;
      else if(option == "--benchmark-data") dataBenchmark = true;
      else if(option == "--benchmark-algebra") algebraBenchmark = true;
      else if(option == "--test-weights") weightTests = true;
      else if(option == "--benchmark-weights") weightBenchmark = true;
      else if(option == "--weight-dimension" && i+1 < argc) { const std::string value = argv[++i]; if(value.empty() || value.find_first_not_of("0123456789") != std::string::npos) throw std::invalid_argument("Invalid weight-test dimension."); weightDimension = static_cast<unsigned>(std::stoul(value)); if(weightDimension < 1 || weightDimension > 3) throw std::invalid_argument("Weight-test dimension must be one, two, or three."); }
      else if(option == "--help") {
        std::cout << "MultiPhase_ex47 [--cases N] [--points N] [--benchmark] [--benchmark-quadrature] [--benchmark-surface] [--benchmark-data] [--benchmark-algebra] [--test-weights] [--benchmark-weights] [--weight-dimension 1|2|3]\n";
        return 0;
      }
      else throw std::invalid_argument("Unknown or incomplete argument: " + option);
    }
    if(quadratureBenchmark) { BenchmarkQuadratureKernels(); return 0; }
    if(surfaceBenchmark) { BenchmarkSelectiveSurfaceKernels(); return 0; }
    if(dataBenchmark) { BenchmarkDataReuse(); return 0; }
    if(algebraBenchmark) { BenchmarkPreparedAlgebra(); return 0; }
    if(weightTests) { TestReviewedWeightState(); if(weightDimension == 0 || weightDimension == 1) TestLineWeightIntegration(); if(weightDimension == 0 || weightDimension == 2) { TestBallWeightIntegration(2); TestElementWeightOptions(2); TestSmallPlanarWeights(); } if(weightDimension == 0 || weightDimension == 3) { TestBallWeightIntegration(3); TestElementWeightOptions(3); TestHexSaddleProjection(); TestCenteredNormalWeights(); TestHighDegreePreparedMoments(); } std::cout << "weight checks passed: " << checks << '\n'; return 0; }
    if(weightBenchmark) { if(weightDimension == 0 || weightDimension == 1) BenchmarkLineWeights(); if(weightDimension == 0 || weightDimension == 2) BenchmarkElementWeights(2); if(weightDimension == 0 || weightDimension == 3) BenchmarkElementWeights(3); return 0; }
    std::mt19937 random(470013);
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    TestReconstruction(cases, points, random);
    TestRestrictions(cases, points, random);
    TestSecondDerivatives(cases, random);
    TestClassification();
    TestBoundarySlices();
    TestComponents(cases, random);
    TestOwnershipAndInvalidation();
    TestPrecision();
    TestSliceRoundoffMargins();
    TestGKRestrictions(random);
    TestInvalidInputs();
    TestNestedIntegration();
    TestPrecisionPaths(benchmark, random);
    TestPhysicalSphere(cases, random);
    TestFastQuadratureInterfaces(cases,random);
    TestSelectiveSurfacePreparation(cases,random);
    TestSelectiveSurfaceScaling();
    TestPreparedLineAnalysis();
    TestCompensatedFoldRestrictions();
    TestBoundaryCertificates(random);
    TestPreparedCharts(cases,points,random);
    TestProjectionAndGraphData(cases,random);
    TestExternalChartComponents(random);
    TestRootSurfaceInterfaces();
    TestLocalFiberAnalysis(cases,random);
    TestPreparedDataInvalidInputs();
    TestReviewedChartContracts(random);
    TestReviewedChartPrecision();
    std::cout << std::setprecision(8) << "root accuracy: rounded-root density=" << maximumRoundedRootError << "; retained-slope density=" << maximumRootDensityError << "; last-ULP fold density=" << maximumFoldDensityError << '\n';
    std::cout << std::setprecision(8) << "ex47: " << checks << " checks passed; maximum scaled error=" << maximumScaledError << "; seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << '\n';
    return 0;
  }
  catch(const std::exception& error) {
    std::cerr << "ex47 failed after " << checks << " checks: " << error.what() << '\n';
    return 1;
  }
}
