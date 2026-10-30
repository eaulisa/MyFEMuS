#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

struct BurstingBubblePoint {
  double r = 0.;
  double z = 0.;
};

struct BurstingBubbleCircle {
  double r = 0.;
  double z = 0.;
  double radius = 0.;
};

struct BurstingBubbleConfig {
  double bondNumber = 0.1172;
  unsigned pointLess = 200u;
  double outerRadius = 10.;
  double freeSurfaceHeight = 0.;
  unsigned filletArcPoints = 64u;
  double simplificationTolerance = 1.e-5;
};

struct BurstingBubbleGeometry {
  std::vector<BurstingBubblePoint> interface;
  BurstingBubbleCircle fillet;
  double curvatureBond = 0.;
  double effectiveBond = 0.;
  double effectiveRadius = 0.;
};

namespace BurstingBubbleDetail {

  constexpr double pi = 3.141592653589793238462643383279502884;
  constexpr double fStep = 1.e-3;
  constexpr double gStep = 1.e-4;
  constexpr double tailStep = 2.e-4;
  constexpr double tailStopTolerance = 1.e-14;
  constexpr double bondRelativeTolerance = 2.e-4;
  constexpr unsigned maxIterations = 100u;

  struct ShapeData {
    std::vector<BurstingBubblePoint> bottom;
    std::vector<BurstingBubblePoint> upper;
    std::vector<double> upperSlope;
    std::vector<double> phi;
    std::vector<BurstingBubblePoint> tail;
    unsigned triplePosition = 0u;
    unsigned bPosition = 0u;
    double phiTriple = 0.;
    double L = 0.;
    double effectiveRadius = 0.;
    double effectiveBond = 0.;
  };

  inline double Sqr(const double x) {
    return x * x;
  }

  inline double Distance(const BurstingBubblePoint& a,
                         const BurstingBubblePoint& b) {
    return std::hypot(a.r - b.r, a.z - b.z);
  }

  inline double Fz(const double r,
                   const double z,
                   const double dzdr,
                   const double Bo) {
    return -dzdr * (1. + dzdr * dzdr) / r
           + (Bo * z + 2.) * std::pow(1. + dzdr * dzdr, 1.5);
  }

  inline double Gx(const double z,
                   const double r,
                   const double drdz,
                   const double Bo) {
    return (1. + drdz * drdz) / r
           - (Bo * z + 2.) * std::pow(1. + drdz * drdz, 1.5);
  }

  inline double TailCurvature(const double r,
                              const double z,
                              const double dzdr,
                              const double Bo,
                              const double L) {
    return Bo * (z - L) * std::pow(1. + dzdr * dzdr, 1.5)
           - dzdr * (1. + dzdr * dzdr) / r;
  }

  template <class Function>
  inline void RK4SecondOrderStep(double& t,
                                 double& y,
                                 double& dy,
                                 const double h,
                                 const Function& rhs) {
    const double k1 = rhs(t, y, dy);

    const double k2 = rhs(t + 0.5 * h,
                          y + 0.5 * h * dy,
                          dy + 0.5 * h * k1);

    const double k3 = rhs(t + 0.5 * h,
                          y + 0.5 * h * dy + 0.5 * h * h * k2,
                          dy + 0.5 * h * k2);

    const double k4 = rhs(t + h,
                          y + h * dy + 0.5 * h * h * k2,
                          dy + h * k3);

    const double yNew = y + h * dy + h * h / 6. * (k1 + k2 + k3);
    const double dyNew = dy + h / 6. * (k1 + 2. * k2 + 2. * k3 + k4);

    t += h;
    y = yNew;
    dy = dyNew;
  }

  inline unsigned FindPhiLowerIndex(const std::vector<double>& phi,
                                    const double target) {
    int bestIndex = -1;
    double bestValue = -std::numeric_limits<double>::max();

    for(unsigned i = 0u; i < phi.size(); ++i) {
      if(phi[i] < target && phi[i] > bestValue) {
        bestValue = phi[i];
        bestIndex = static_cast<int>(i);
      }
    }

    if(bestIndex < 0 || static_cast<unsigned>(bestIndex + 1) >= phi.size()) {
      throw std::runtime_error(
        "BurstingBubble: unable to bracket triple-point angle");
    }

    return static_cast<unsigned>(bestIndex);
  }

  inline unsigned FindNearestBubblePoint(const BurstingBubblePoint& A,
                                         const std::vector<BurstingBubblePoint>& upper,
                                         const unsigned lastIndex) {
    if(upper.empty()) {
      throw std::runtime_error("BurstingBubble: empty upper branch");
    }

    const unsigned end = std::min(lastIndex,
                                  static_cast<unsigned>(upper.size() - 1u));

    unsigned best = 0u;
    double bestDistance = std::numeric_limits<double>::max();

    for(unsigned i = 0u; i <= end; ++i) {
      const double d = Distance(A, upper[i]);
      if(d < bestDistance) {
        bestDistance = d;
        best = i;
      }
    }

    return best;
  }

  inline double RevolutionVolumeSegment(const BurstingBubblePoint& a,
                                        const BurstingBubblePoint& b) {
    const double meanRadius = 0.5 * (a.r + b.r);
    return pi * meanRadius * meanRadius * (b.z - a.z);
  }

  inline double SphericalCapVolume(const BurstingBubblePoint& p,
                                   const BurstingBubblePoint& previous) {
    const double alpha = p.r - previous.r;
    const double beta = p.z - previous.z;

    if(std::abs(beta) <= 1.e-14) {
      throw std::runtime_error(
        "BurstingBubble: degenerate spherical-cap construction");
    }

    const double c = -alpha * p.r - beta * p.z;
    const double zCenter = -c / beta;
    const double radius = std::hypot(p.r, p.z - zCenter);
    const double theta0 = std::atan(-alpha / beta);
    const unsigned n = 1000u;
    const double dtheta = (0.5 * pi - theta0) / static_cast<double>(n - 1u);

    double theta = theta0;
    BurstingBubblePoint oldPoint = {
      radius * std::cos(theta),
      radius * std::sin(theta) + zCenter
    };

    double volume = 0.;

    for(unsigned i = 1u; i < n; ++i) {
      theta += dtheta;

      const BurstingBubblePoint newPoint = {
        radius * std::cos(theta),
        radius * std::sin(theta) + zCenter
      };

      volume += RevolutionVolumeSegment(oldPoint, newPoint);
      oldPoint = newPoint;
    }

    return volume;
  }

  inline std::vector<BurstingBubblePoint>
  IntegrateTail(const double r0,
                const double z0,
                const double slope0,
                const double Bo,
                const double L,
                const double outerRadius,
                double& minimumZ) {
    std::vector<BurstingBubblePoint> tail;

    double r = r0;
    double z = z0;
    double slope = slope0;

    tail.push_back({r, z});
    minimumZ = z;

    while(r < outerRadius) {
      RK4SecondOrderStep(
        r,
        z,
        slope,
        tailStep,
      [Bo, L](const double rr, const double zz, const double ss) {
        return TailCurvature(rr, zz, ss, Bo, L);
      });

      if(!std::isfinite(r) || !std::isfinite(z) || !std::isfinite(slope)) {
        break;
      }

      tail.push_back({r, z});
      minimumZ = std::min(minimumZ, z);

      if(std::abs(slope) > 1000.) {
        break;
      }

      if(tail.size() > 1000000u) {
        throw std::runtime_error("BurstingBubble: tail integration did not stop");
      }
    }

    return tail;
  }

  inline ShapeData SolveShape(const double Bo,
                              const unsigned pointLess,
                              const double outerRadius) {
    if(!(Bo > 0.)) {
      throw std::runtime_error("BurstingBubble: curvature Bond must be positive");
    }

    ShapeData shape;

    double r = 1.e-8;
    double z = 0.;
    double dzdr = 0.;

    shape.bottom.push_back({0., 0.});

    while(dzdr < 1.) {
      RK4SecondOrderStep(
        r,
        z,
        dzdr,
        fStep,
      [Bo](const double rr, const double zz, const double ss) {
        return Fz(rr, zz, ss, Bo);
      });

      if(!std::isfinite(r) || !std::isfinite(z) || !std::isfinite(dzdr)) {
        throw std::runtime_error("BurstingBubble: lower branch diverged");
      }

      shape.bottom.push_back({r, z});

      if(shape.bottom.size() > 100000u) {
        throw std::runtime_error("BurstingBubble: lower branch did not stop");
      }
    }

    double drdz = 1.;

    shape.upper.push_back({r, z});
    shape.upperSlope.push_back(1.);
    shape.phi.push_back(45.);

    while(drdz >= -1000.) {
      RK4SecondOrderStep(
        z,
        r,
        drdz,
        gStep,
      [Bo](const double zz, const double rr, const double ss) {
        return Gx(zz, rr, ss, Bo);
      });

      if(!std::isfinite(r) || !std::isfinite(z) || !std::isfinite(drdz)) {
        break;
      }

      double dzdrLocal = 0.;
      if(std::abs(drdz) <= 1.e-15) {
        dzdrLocal = drdz >= 0. ? 1.e15 : -1.e15;
      }
      else {
        dzdrLocal = 1. / drdz;
      }

      double phi = std::atan(dzdrLocal) * 180. / pi;
      if(phi < 0.) {
        phi += 180.;
      }

      shape.upper.push_back({r, z});
      shape.upperSlope.push_back(dzdrLocal);
      shape.phi.push_back(phi);

      if(shape.upper.size() > 100000u) {
        throw std::runtime_error("BurstingBubble: upper branch did not stop");
      }
    }

    if(shape.upper.size() < 3u) {
      throw std::runtime_error("BurstingBubble: upper branch is too short");
    }

    double phiMin = 90.;
    double phiMax = 180.;
    double phiC = 0.5 * (phiMin + phiMax);

    unsigned position = FindPhiLowerIndex(shape.phi, phiC);

    double rInit = 0.5 * (shape.upper[position].r + shape.upper[position + 1u].r);
    double zInit = 0.5 * (shape.upper[position].z + shape.upper[position + 1u].z);
    double slopeInit = 0.5 * (shape.upperSlope[position]
                              + shape.upperSlope[position + 1u]);

    double rInitOld = 0.;
    double zInitOld = 0.;
    double L = 0.;
    std::vector<BurstingBubblePoint> tail;

    unsigned tailIteration = 0u;

    while(std::abs(rInit - rInitOld) > tailStopTolerance
          && std::abs(zInit - zInitOld) > tailStopTolerance) {
      if(++tailIteration > maxIterations) {
        throw std::runtime_error("BurstingBubble: triple-point search did not converge");
      }

      const double R = rInit / std::sin(phiC * pi / 180.);
      L = 2. / Bo * (2. / R - 1.);

      double minimumZ = zInit;
      tail = IntegrateTail(rInit,
                           zInit,
                           slopeInit,
                           Bo,
                           L,
                           outerRadius,
                           minimumZ);

      if(tail.size() < 2u) {
        throw std::runtime_error("BurstingBubble: tail integration failed");
      }

      if(minimumZ < L) {
        phiMin = phiC;
        phiC = 0.5 * (phiC + phiMax);
      }
      else {
        phiMax = phiC;
        phiC = 0.5 * (phiC + phiMin);
      }

      rInitOld = rInit;
      zInitOld = zInit;

      position = FindPhiLowerIndex(shape.phi, phiC);

      const double denominator = shape.phi[position + 1u] - shape.phi[position];
      if(std::abs(denominator) <= 1.e-15) {
        throw std::runtime_error("BurstingBubble: degenerate phi interpolation");
      }

      const double ratio = (phiC - shape.phi[position]) / denominator;

      rInit = shape.upper[position].r
              + ratio * (shape.upper[position + 1u].r - shape.upper[position].r);

      zInit = shape.upper[position].z
              + ratio * (shape.upper[position + 1u].z - shape.upper[position].z);

      slopeInit = shape.upperSlope[position]
                  + ratio * (shape.upperSlope[position + 1u]
                             - shape.upperSlope[position]);
    }

    if(tail.empty()) {
      throw std::runtime_error("BurstingBubble: no tail was generated");
    }

    if(pointLess >= tail.size()) {
      throw std::runtime_error("BurstingBubble: pointLess exceeds tail size");
    }

    shape.tail = std::move(tail);
    shape.triplePosition = position;
    shape.phiTriple = phiC;
    shape.L = L;

    const BurstingBubblePoint A = shape.tail[pointLess];

    if(position == 0u) {
      throw std::runtime_error("BurstingBubble: invalid triple-point index");
    }

    shape.bPosition = FindNearestBubblePoint(A, shape.upper, position - 1u);

    double bubbleVolume = SphericalCapVolume(shape.upper[position],
                          shape.upper[position - 1u]);

    for(unsigned i = 1u; i < shape.bottom.size(); ++i) {
      BurstingBubblePoint a = shape.bottom[i - 1u];
      BurstingBubblePoint b = shape.bottom[i];
      a.z -= L;
      b.z -= L;
      bubbleVolume += RevolutionVolumeSegment(a, b);
    }

    BurstingBubblePoint previous = shape.bottom.back();
    previous.z -= L;

    for(unsigned i = 0u; i <= shape.bPosition; ++i) {
      BurstingBubblePoint current = shape.upper[i];
      current.z -= L;
      bubbleVolume += RevolutionVolumeSegment(previous, current);
      previous = current;
    }

    if(!(bubbleVolume > 0.) || !std::isfinite(bubbleVolume)) {
      throw std::runtime_error("BurstingBubble: invalid bubble volume");
    }

    shape.effectiveRadius = std::cbrt(3. * bubbleVolume / (4. * pi));
    shape.effectiveBond = Bo * shape.effectiveRadius * shape.effectiveRadius;

    return shape;
  }

  inline double FindCurvatureBond(const double targetBond,
                                  const double outerRadius) {
    if(!(targetBond > 0.)) {
      throw std::runtime_error("BurstingBubble: target Bond must be positive");
    }

    if(targetBond > 30.) {
      throw std::runtime_error("BurstingBubble: target Bond is above supported range");
    }

    double BoMin = 0.;
    double BoMax = 1000.;
    double BoInput = 0.5 * (BoMin + BoMax);

    for(unsigned iteration = 0u; iteration < maxIterations; ++iteration) {
      const ShapeData shape = SolveShape(BoInput, 0u, outerRadius);

      const double relativeError = std::abs(shape.effectiveBond - targetBond)
                                   / targetBond;

      if(relativeError <= bondRelativeTolerance) {
        return BoInput;
      }

      if(shape.effectiveBond > targetBond) {
        BoMax = BoInput;
      }
      else {
        BoMin = BoInput;
      }

      BoInput = 0.5 * (BoMin + BoMax);
    }

    throw std::runtime_error("BurstingBubble: Bond conversion did not converge");
  }

  inline BurstingBubbleCircle FindFilletCircle(const BurstingBubblePoint& B,
      const BurstingBubblePoint& A,
      const double slopeAtB) {
    const double norm = std::sqrt(1. + slopeAtB * slopeAtB);

    const double denominator = 2. * ((B.r - A.r) * slopeAtB
                                     + A.z - B.z);

    if(std::abs(denominator) <= 1.e-14) {
      throw std::runtime_error("BurstingBubble: singular fillet construction");
    }

    const double signedRadius = norm
                                * (Sqr(A.r - B.r) + Sqr(A.z - B.z))
                                / denominator;

    BurstingBubbleCircle circle;
    circle.r = B.r - signedRadius * slopeAtB / norm;
    circle.z = B.z + signedRadius / norm;
    circle.radius = std::abs(signedRadius);
    return circle;
  }

  inline double NormalizePositiveAngle(double angle) {
    const double twoPi = 2. * pi;
    while(angle < 0.) angle += twoPi;
    while(angle >= twoPi) angle -= twoPi;
    return angle;
  }

  inline std::vector<BurstingBubblePoint>
  BuildFilletArc(const BurstingBubblePoint& previousB,
                 const BurstingBubblePoint& B,
                 const BurstingBubblePoint& A,
                 const BurstingBubbleCircle& circle,
                 const unsigned nPoints) {
    if(nPoints < 2u) {
      return {B, A};
    }

    const double thetaB = std::atan2(B.z - circle.z, B.r - circle.r);
    const double thetaA = std::atan2(A.z - circle.z, A.r - circle.r);

    const double radiusB = std::hypot(B.r - circle.r, B.z - circle.z);
    if(radiusB <= 1.e-14) {
      throw std::runtime_error("BurstingBubble: invalid fillet radius at B");
    }

    const BurstingBubblePoint incoming = {
      B.r - previousB.r,
      B.z - previousB.z
    };

    const double incomingNorm = std::hypot(incoming.r, incoming.z);
    if(incomingNorm <= 1.e-14) {
      throw std::runtime_error("BurstingBubble: degenerate incoming tangent at B");
    }

    const BurstingBubblePoint tangentCCW = {
      -(B.z - circle.z) / radiusB,
      (B.r - circle.r) / radiusB
    };

    const BurstingBubblePoint tangentCW = {
      -tangentCCW.r,
      -tangentCCW.z
    };

    const double dotCCW = tangentCCW.r * incoming.r + tangentCCW.z * incoming.z;
    const double dotCW = tangentCW.r * incoming.r + tangentCW.z * incoming.z;

    double delta = 0.;

    if(dotCCW >= dotCW) {
      delta = NormalizePositiveAngle(thetaA - thetaB);
    }
    else {
      delta = -NormalizePositiveAngle(thetaB - thetaA);
    }

    std::vector<BurstingBubblePoint> arc;
    arc.reserve(nPoints);

    for(unsigned i = 0u; i < nPoints; ++i) {
      const double s = static_cast<double>(i) / static_cast<double>(nPoints - 1u);
      const double theta = thetaB + s * delta;
      arc.push_back({circle.r + circle.radius * std::cos(theta),
                     circle.z + circle.radius * std::sin(theta)});
    }

    arc.front() = B;
    arc.back() = A;
    return arc;
  }

  inline double PointSegmentDistanceSquared(const BurstingBubblePoint& p,
      const BurstingBubblePoint& a,
      const BurstingBubblePoint& b) {
    const double dr = b.r - a.r;
    const double dz = b.z - a.z;
    const double lengthSquared = dr * dr + dz * dz;

    if(lengthSquared <= 1.e-30) {
      return Sqr(p.r - a.r) + Sqr(p.z - a.z);
    }

    double t = ((p.r - a.r) * dr + (p.z - a.z) * dz) / lengthSquared;
    t = std::max(0., std::min(1., t));

    const double qr = a.r + t * dr;
    const double qz = a.z + t * dz;
    return Sqr(p.r - qr) + Sqr(p.z - qz);
  }

  inline std::vector<BurstingBubblePoint>
  SimplifyPolyline(const std::vector<BurstingBubblePoint>& points,
                   const double tolerance) {
    if(points.size() <= 2u || tolerance <= 0.) {
      return points;
    }

    std::vector<unsigned char> keep(points.size(), 0u);
    keep.front() = 1u;
    keep.back() = 1u;

    std::vector<std::pair<std::size_t, std::size_t>> stack;
    stack.push_back({0u, points.size() - 1u});

    const double toleranceSquared = tolerance * tolerance;

    while(!stack.empty()) {
      const auto range = stack.back();
      stack.pop_back();

      const std::size_t first = range.first;
      const std::size_t last = range.second;

      if(last <= first + 1u) {
        continue;
      }

      double maxDistanceSquared = -1.;
      std::size_t maxIndex = first;

      for(std::size_t i = first + 1u; i < last; ++i) {
        const double d2 = PointSegmentDistanceSquared(points[i],
                          points[first],
                          points[last]);
        if(d2 > maxDistanceSquared) {
          maxDistanceSquared = d2;
          maxIndex = i;
        }
      }

      if(maxDistanceSquared > toleranceSquared) {
        keep[maxIndex] = 1u;
        stack.push_back({first, maxIndex});
        stack.push_back({maxIndex, last});
      }
    }

    std::vector<BurstingBubblePoint> simplified;
    simplified.reserve(points.size());

    for(std::size_t i = 0u; i < points.size(); ++i) {
      if(keep[i]) {
        simplified.push_back(points[i]);
      }
    }

    return simplified;
  }

  inline double DistanceToInterface(const BurstingBubblePoint& p,
                                    const std::vector<BurstingBubblePoint>& interface) {
    if(interface.size() < 2u) {
      throw std::runtime_error("BurstingBubble: interface contains fewer than two points");
    }

    double distanceSquared = std::numeric_limits<double>::max();

    for(std::size_t i = 0u; i + 1u < interface.size(); ++i) {
      distanceSquared = std::min(distanceSquared,
                                 PointSegmentDistanceSquared(p,
                                     interface[i],
                                     interface[i + 1u]));
    }

    return std::sqrt(distanceSquared);
  }

  inline bool IsGasSide(const BurstingBubblePoint& p,
                        const std::vector<BurstingBubblePoint>& interface) {
    if(interface.size() < 2u) {
      return false;
    }

    double zTop = p.z + 1.;
    for(const BurstingBubblePoint& q : interface) {
      zTop = std::max(zTop, q.z + 2.);
    }

    std::vector<BurstingBubblePoint> polygon = interface;
    polygon.push_back({interface.back().r, zTop});
    polygon.push_back({0., zTop});

    bool inside = false;
    const std::size_t n = polygon.size();

    for(std::size_t i = 0u, j = n - 1u; i < n; j = i++) {
      const BurstingBubblePoint& a = polygon[i];
      const BurstingBubblePoint& b = polygon[j];

      if((a.z > p.z) == (b.z > p.z)) {
        continue;
      }

      const double rCross = a.r
                            + (p.z - a.z) * (b.r - a.r) / (b.z - a.z);

      if(p.r < rCross) {
        inside = !inside;
      }
    }

    return inside;
  }

} // namespace BurstingBubbleDetail

inline BurstingBubbleGeometry
ComputeBurstingBubbleGeometry(const BurstingBubbleConfig& cfg) {
  using namespace BurstingBubbleDetail;

  if(!(cfg.bondNumber > 0.)) {
    throw std::runtime_error(
      "ComputeBurstingBubbleGeometry: bondNumber must be positive");
  }

  if(!(cfg.outerRadius > 0.)) {
    throw std::runtime_error(
      "ComputeBurstingBubbleGeometry: outerRadius must be positive");
  }

  const double curvatureBond = FindCurvatureBond(cfg.bondNumber,
                               cfg.outerRadius);

  const ShapeData shape = SolveShape(curvatureBond,
                                     cfg.pointLess,
                                     cfg.outerRadius);

  if(shape.bPosition == 0u) {
    throw std::runtime_error(
      "ComputeBurstingBubbleGeometry: B point has no previous bubble point");
  }

  const BurstingBubblePoint Araw = shape.tail[cfg.pointLess];
  const BurstingBubblePoint Braw = shape.upper[shape.bPosition];

  BurstingBubblePoint A = {Araw.r, Araw.z - shape.L};
  BurstingBubblePoint B = {Braw.r, Braw.z - shape.L};

  BurstingBubbleCircle fillet = FindFilletCircle(B,
                                A,
                                shape.upperSlope[shape.bPosition]);

  const double scale = 1. / shape.effectiveRadius;

  auto ScalePoint = [&](const BurstingBubblePoint & p) {
    return BurstingBubblePoint{
      p.r * scale,
      p.z * scale + cfg.freeSurfaceHeight
    };
  };

  BurstingBubbleGeometry geometry;
  geometry.curvatureBond = curvatureBond;
  geometry.effectiveBond = shape.effectiveBond;
  geometry.effectiveRadius = shape.effectiveRadius;

  geometry.fillet = {
    fillet.r * scale,
    fillet.z * scale + cfg.freeSurfaceHeight,
    fillet.radius * scale
  };

  for(const BurstingBubblePoint& p : shape.bottom) {
    const BurstingBubblePoint shifted = {p.r, p.z - shape.L};
    const BurstingBubblePoint q = ScalePoint(shifted);

    if(geometry.interface.empty()
        || Distance(geometry.interface.back(), q) > 1.e-14) {
      geometry.interface.push_back(q);
    }
  }

  for(unsigned i = 0u; i <= shape.bPosition; ++i) {
    const BurstingBubblePoint shifted = {
      shape.upper[i].r,
      shape.upper[i].z - shape.L
    };

    const BurstingBubblePoint q = ScalePoint(shifted);

    if(geometry.interface.empty()
        || Distance(geometry.interface.back(), q) > 1.e-14) {
      geometry.interface.push_back(q);
    }
  }

  const BurstingBubblePoint scaledA = ScalePoint(A);
  const BurstingBubblePoint scaledB = ScalePoint(B);

  const BurstingBubblePoint previousB = geometry.interface.size() >= 2u
                                        ? geometry.interface[geometry.interface.size() - 2u]
                                        : geometry.interface.back();

  const std::vector<BurstingBubblePoint> arc = BuildFilletArc(previousB,
      scaledB,
      scaledA,
      geometry.fillet,
      cfg.filletArcPoints);

  for(std::size_t i = 1u; i < arc.size(); ++i) {
    if(Distance(geometry.interface.back(), arc[i]) > 1.e-14) {
      geometry.interface.push_back(arc[i]);
    }
  }

  for(unsigned i = cfg.pointLess; i < shape.tail.size(); ++i) {
    const BurstingBubblePoint shifted = {
      shape.tail[i].r,
      shape.tail[i].z - shape.L
    };

    const BurstingBubblePoint q = ScalePoint(shifted);

    if(Distance(geometry.interface.back(), q) > 1.e-14) {
      geometry.interface.push_back(q);
    }
  }

  if(geometry.interface.back().r < cfg.outerRadius) {
    geometry.interface.push_back({cfg.outerRadius,
                                  geometry.interface.back().z});
  }

  geometry.interface = SimplifyPolyline(geometry.interface,
                                        cfg.simplificationTolerance);

  return geometry;
}

inline double BurstingBubbleLevelSet(const std::vector<double>& x,
                                     const BurstingBubbleGeometry& geometry) {
  using namespace BurstingBubbleDetail;

  if(x.size() != 2u && x.size() != 3u) {
    throw std::runtime_error(
      "BurstingBubbleLevelSet: only 2D axisymmetric or 3D coordinates are supported");
  }

  const double r = x.size() == 2u
                   ? std::abs(x[0])
                   : std::hypot(x[0], x[1]);

  const double z = x.size() == 2u ? x[1] : x[2];

  const BurstingBubblePoint p = {r, z};
  const double distance = DistanceToInterface(p, geometry.interface);

  // Same sign convention used by Psi::BUBBLES_FREE_SURFACE:
  // psi > 0 in the gas, psi < 0 in the liquid.
  return IsGasSide(p, geometry.interface) ? distance : -distance;
}
