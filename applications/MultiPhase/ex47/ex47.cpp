#include "include/CutFemLevelSetData.hpp"
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

std::vector<long double> RandomPoint(GeomElType geometry, std::mt19937& random) {
  std::uniform_real_distribution<double> uniform(0.01, 0.99);
  std::vector<long double> point(Dimension(geometry));
  if(geometry == TRI || geometry == TET || geometry == WEDGE) {
    point[0] = uniform(random);
    point[1] = (1 - point[0]) * uniform(random);
    if(geometry == TET) point[2] = (1 - point[0] - point[1]) * uniform(random);
    if(geometry == WEDGE) point[2] = 2 * uniform(random) - 1;
  }
  else for(long double& coordinate : point) coordinate = 2 * uniform(random) - 1;
  return point;
}

long double ReferenceValue(const femus::basis& basis, const std::vector<double>& nodal, const std::vector<long double>& point, int derivative = -1) {
  double native[3] = {0, 0, 0};
  for(unsigned axis = 0; axis < point.size(); ++axis) native[axis] = static_cast<double>(point[axis]);
  long double result = 0;
  for(unsigned n = 0; n < nodal.size(); ++n) result += static_cast<long double>(nodal[n]) * (derivative < 0 ? basis.eval_phi(basis.GetIND(n), native) : basis.eval_dphidxyz(static_cast<unsigned>(derivative), basis.GetIND(n), native));
  return result;
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
        for(unsigned axis = 0; axis < point.size(); ++axis) point[axis] = coordinate[axis];
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
        const long double coordinate = simplexAxis ? position(random) : 2 * position(random) - 1;
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
              const long double secondCoordinate = triangle ? 0.37L : -0.23L;
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
  for(unsigned axis = 0; axis < point.size(); ++axis) native[axis] = static_cast<double>(point[axis]);
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
  return result;
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
      const long double endpoints[4] = {simplexAxis ? 0.0L : -1.0L, 1.0L, std::nextafter(1.0L, 0.0L), simplexAxis ? 1e-15L : -1 + 1e-15L};
      for(long double coordinate : endpoints) {
        parent.Restrict(axis, coordinate, child);
        const std::vector<long double> point(child.GetDimension(), 0);
        CheckNear(child.Evaluate(point), parent.Evaluate(child.MapToParent(point)), 5e-17L, "endpoint power restriction");
        CheckNear(child.EvaluateBernstein(point), child.Evaluate(point), 5e-16L, "endpoint Bernstein restriction");
        const bool collapsed = simplexAxis && coordinate == 1;
        Check(child.IsDegenerate() == collapsed, "degenerate slice flag");
        long double expected = 1;
        if(geometry == TET) expected = (1 - coordinate) * (1 - coordinate);
        if(geometry == TRI || (geometry == WEDGE && axis < 2)) expected = (1 - coordinate) / 2;
        CheckNear(child.GetParentMeasureScale(), expected, 0, "slice measure");
        if(child.GetDimension() == 2) {
          child.Restrict(0, child.GetGeometry() == TRI ? 1.0L : -1.0L, line);
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
        const long double coordinate = simplexAxis ? 0.31L : -0.27L;
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
  const std::vector<long double> childPoint = {0.13L, -0.29L};
  const long double saved = child.Evaluate(childPoint);
  const std::vector<long double> savedGradient = child.GetRootGradient(childPoint);
  CutFemLevelSetData copied = parent;
  parent.Reset(HEX, nodal);
  CheckNear(parent.Evaluate(point), 42, 0, "new nodal state");
  CheckNear(copied.Evaluate(point), original, 0, "independent copy");
  CheckNear(child.Evaluate(childPoint), saved, 0, "slice snapshot survives parent reset");
  Check(child.GetRootGradient(childPoint) == savedGradient, "slice gradient survives parent reset");
  CutFemLevelSetData reusable;
  copied.Restrict(2, -0.4L, reusable);
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
    CheckNear(roots[0] / expected, -1, 3e-18L, "negative root versus multiprecision");
    CheckNear(roots[1] / expected, 1, 3e-18L, "positive root versus multiprecision");
  }
  CutFemLevelSetData tangent(LINE, {1, 1, 0});
  Check(tangent.GetLineRoots() == std::vector<long double> {0}, "repeated root");
  CutFemLevelSetData endpoints(LINE, {0, 0, -1});
  Check(endpoints.GetLineRoots() == std::vector<long double>({-1, 1}), "boundary roots");
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
    const long double offset = geometry == TRI || geometry == TET || geometry == WEDGE ? 0.25L : 0.2L;
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
  CheckNear(roots[0] / static_cast<long double>(c / q), 1, 3e-18L, "small root avoids subtraction cancellation");
  CheckNear(roots[1], static_cast<long double>(q / a), 3e-18L, "large asymmetric root");
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
          const long double coordinate = simplexAxis ? (1 + sign * node) / 2 : sign * node;
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
    std::vector<long double> cuts = {-1};
    for(long double root : roots) if(root > -1 && root < 1) cuts.push_back(root);
    cuts.push_back(1);
    for(unsigned i = 1; i < cuts.size(); ++i) if(data.Evaluate({(cuts[i - 1] + cuts[i]) / 2}) > 0) result[0] += cuts[i] - cuts[i - 1];
    for(long double root : roots) {
      const std::vector<long double> gradient = data.GetRootGradient({root});
      long double norm = 0;
      for(long double value : gradient) norm = std::hypot(norm, value);
      const long double derivative = std::abs(data.EvaluateDerivative(0, {root}));
      if(derivative > 0) result[1] += (root == -1 || root == 1 ? 0.5L : 1.0L) * norm / derivative;
    }
    return result;
  }
  const unsigned axis = data.GetDimension() - 1;
  const bool triangle = data.GetGeometry() == TRI;
  const bool tetrahedron = data.GetGeometry() == TET;
  const std::vector<long double> panels = triangle ? std::vector<long double> {0, 0.5L, 1} :
                                          (tetrahedron ? std::vector<long double> {0, 1} : std::vector<long double> {-1, 1});
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
  const long double volumes[6] = {2, 4, 0.5L, 8, 1.0L / 6, 1};
  const long double surfaces[6] = {1, 2 * std::sqrt(2.0L), 1 / std::sqrt(2.0L), 4 * std::sqrt(2.0L), 1 / (2 * std::sqrt(2.0L)), std::sqrt(2.0L)};
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
        for(unsigned axis = 0; axis < point.size(); ++axis) point[axis] = native[axis];
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
  for(unsigned n = 0; n < coordinates[0].size(); ++n) { const double* native = basis->GetXcoarse(n); std::vector<long double> point(dimension), physical; std::vector<std::vector<long double>> jacobian; for(unsigned axis = 0; axis < dimension; ++axis) point[axis] = native[axis]; SurfaceTestMap(point,mode,physical,jacobian); for(unsigned component = 0; component < dimension; ++component) coordinates[component][n] = static_cast<double>(physical[component]); }
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
  parent.Restrict(2,-0.25L,reusable);
  Check(reusable.GetGraphDirection() == 1 && reusable.GetGraphDirection() == 1,"cached graph certificate");
  Check(reusable.GetWorkCounts()[2] == 0,"graph bounds need no derivative reconstruction");
  parent.Restrict(2,0.25L,reusable);
  Check(reusable.GetGraphDirection() == 0,"graph certificate invalidated with slice state");
  { CutFemLevelSetData temporary(HEX,nodal); temporary.Restrict(2,0.25L,snapshot); }
  CheckNear(snapshot.Evaluate({0.2L,-0.4L}),0.1L,2e-18L,"slice survives destruction of its parent");
  CutFemLevelSetData tangent(LINE,{1,1,0}), regular(LINE,{-1,1,0});
  std::array<long double,3> gradient, normal;
  std::array<long double,4> scalar;
  bool rejected = false;
  try { tangent.GetLineSurfaceData(0,gradient); } catch(const std::domain_error&) { rejected = true; }
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
  CheckNear(line.GetLineDeltaFactor(0),1,0,"delta needs only the line derivative");
  Check(line.GetWorkCounts()[2] == before[2],"delta does not prepare even local derivative coefficients");
  long double savedDelta = 0;
  { CutFemLevelSetData temporary(HEX,MakeField(HEX,1)), child; temporary.Restrict(2,0.17L,child,false); child.Restrict(1,0.19L,saved); savedDelta = saved.GetLineDeltaFactor(0.2L); }
  CheckNear(saved.GetLineDeltaFactor(0.2L),savedDelta,0,"selective snapshot survives parent destruction");
  bool rejected = false;
  try { tangent.GetLineDeltaFactor(0); } catch(const std::domain_error&) { rejected = true; }
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
    const long double exactDelta = static_cast<long double>(1/(2*a*exactRoot));
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
  child.Restrict(child.GetDimension()-1,child.GetGeometry() == TRI ? 0.37L : -0.23L,line);
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
    if(parent.GetDimension() == 3) { parent.Restrict(2,positions[0],child); child.Restrict(1,child.GetGeometry() == TRI ? 0.17L : -0.13L,line); }
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

int main(int argc, char** argv) {
  try {
    unsigned cases = 64, points = 32;
    bool benchmark = false;
    bool quadratureBenchmark = false;
    bool surfaceBenchmark = false;
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
      else if(option == "--help") {
        std::cout << "MultiPhase_ex47 [--cases N] [--points N] [--benchmark] [--benchmark-quadrature] [--benchmark-surface]\n";
        return 0;
      }
      else throw std::invalid_argument("Unknown or incomplete argument: " + option);
    }
    if(quadratureBenchmark) { BenchmarkQuadratureKernels(); return 0; }
    if(surfaceBenchmark) { BenchmarkSelectiveSurfaceKernels(); return 0; }
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
    std::cout << std::setprecision(8) << "ex47: " << checks << " checks passed; maximum scaled error=" << maximumScaledError << "; seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << '\n';
    return 0;
  }
  catch(const std::exception& error) {
    std::cerr << "ex47 failed after " << checks << " checks: " << error.what() << '\n';
    return 1;
  }
}
