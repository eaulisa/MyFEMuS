#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if defined(__GNUC__) || defined(__clang__)
  #define FEMUS_RESTRICT __restrict__
#else
  #define FEMUS_RESTRICT
#endif

class RungeKutta {
public:
  enum class VelKind : unsigned {
    RisingBubble = 0,
    Rotation     = 1,
    Vortex       = 2
  };

  // Moves from (time - dt) -> time. Input Xp is at old time.
  static inline void rkForward(std::vector<std::vector<double>>& Xp,
                               const double time,
                               const double dt,
                               const VelKind kind)
  {
    if (dt == 0.0) return;
    rk4_inplace(Xp, time - dt, dt, kind);
  }

  // Moves from time -> (time - dt). Input Xp is at new time.
  static inline void rkBackward(std::vector<std::vector<double>>& Xp,
                                const double time,
                                const double dt,
                                const VelKind kind)
  {
    if (dt == 0.0) return;
    rk4_inplace(Xp, time, -dt, kind);
  }

  // ---------------- Specific velocities ----------------

  static inline std::array<double, 2>
  RisingBubble2D(const double x, const double y, const double time) noexcept
  {
    const double Lx = 1.0, Ly = 1.0, A = 0.1;
    const double pi = M_PI, T = 8.0;
    const double cost = std::cos(pi * time / T);

    const double s1 = std::sin(pi * (x + 0.5 * Lx) / Lx);
    const double X  = (s1 * s1) * std::sin(2.0 * pi * x / Lx);

    const double sy = std::sin(pi * y / Ly);
    const double cy = std::cos(pi * y / Ly);
    const double Y  = sy * sy;
    const double dYdy = 2.0 * sy * cy * (pi / Ly);

    const double a = pi / Lx;
    const double s = std::sin(a * (x + 0.5 * Lx));
    const double c = std::cos(a * (x + 0.5 * Lx));
    const double Xprime =
      2.0 * a * s * c * std::sin(2.0 * a * x) +
      (s * s) * 2.0 * a * std::cos(2.0 * a * x);

    const double u = -A * X * dYdy * cost;
    const double v =  A * Xprime * Y * cost;
    return {u, v};
  }

  static inline std::array<double, 3>
  RisingBubble3D(const double x, const double y, const double z, const double time) noexcept
  {
    const double Lx = 1.0, Ly = 1.0, Lz = 1.0, A = 0.1;
    const double pi = M_PI, T = 4.0, cost = std::cos(pi * time / T);
    const double ax = pi / Lx, ay = pi / Ly, az = pi / Lz;

    const double sxh = std::sin(ax * (x + 0.5 * Lx));
    const double cxh = std::cos(ax * (x + 0.5 * Lx));
    const double X   = (sxh * sxh) * std::sin(2.0 * ax * x);
    const double Xp  = 2.0 * ax * sxh * cxh * std::sin(2.0 * ax * x)
                     + (sxh * sxh) * 2.0 * ax * std::cos(2.0 * ax * x);

    const double syh = std::sin(ay * (y + 0.5 * Ly));
    const double cyh = std::cos(ay * (y + 0.5 * Ly));
    const double Y   = (syh * syh) * std::sin(2.0 * ay * y);
    const double Yp  = 2.0 * ay * syh * cyh * std::sin(2.0 * ay * y)
                     + (syh * syh) * 2.0 * ay * std::cos(2.0 * ay * y);

    const double sz = std::sin(az * z), cz = std::cos(az * z);
    const double Z  = sz * sz;
    const double Zp = 2.0 * sz * cz * az;

    const double u = -A * X * Zp * cost;
    const double v = -A * Y * Zp * cost;
    const double w =  A * (Xp + Yp) * Z * cost;
    return {u, v, w};
  }

  static inline std::array<double, 2>
  RotationVel2D(const double x, const double y, const double /*time*/) noexcept
  { return {y, -x}; }

  static inline std::array<double, 3>
  RotationVel3D(const double x, const double y, const double z, const double /*time*/) noexcept
  { return {y - z, z - x, x - y}; }

  static inline std::array<double, 2>
  VortexVel2D(double x, double y, const double time) noexcept
  {
    const double T = 8.0;
    x += 0.5; y += 0.5;
    const double sx = std::sin(M_PI * x);
    const double cx = std::cos(M_PI * x);
    const double sy = std::sin(M_PI * y);
    const double cy = std::cos(M_PI * y);
    const double cosT = std::cos(M_PI * time / T);

    const double u = -2.0 * (sx*sx) * (sy*cy) * cosT;
    const double v =  2.0 * (sx*cx) * (sy*sy) * cosT;
    return {u, v};
  }

  static inline std::array<double, 3>
  VortexVel3D(double x, double y, double z, const double time) noexcept
  {
    x += 0.5; y += 0.5; z += 0.5;

    double sx, cx, sy, cy, sz, cz;
  #if defined(__GLIBC__) || defined(__gnu_linux__)
    ::sincos(M_PI * x, &sx, &cx);
    ::sincos(M_PI * y, &sy, &cy);
    ::sincos(M_PI * z, &sz, &cz);
  #else
    sx = std::sin(M_PI * x); cx = std::cos(M_PI * x);
    sy = std::sin(M_PI * y); cy = std::cos(M_PI * y);
    sz = std::sin(M_PI * z); cz = std::cos(M_PI * z);
  #endif

    const double T = 4.0;
    const double cosT = std::cos(M_PI * (time / T));
    const double s2x = sx * sx, s2y = sy * sy, s2z = sz * sz;
    const double sxcx = sx * cx, sycy = sy * cy, szcz = sz * cz;
    const double k = 2.0 * cosT;
    return { k*s2x*(-sycy + szcz), k*s2y*(-szcz + sxcx), k*s2z*(-sxcx + sycy) };
  }

private:
  using VelEvalFn2D = void(*)(double, double, double, double&, double&) noexcept;
  using VelEvalFn3D = void(*)(double, double, double, double, double&, double&, double&) noexcept;

  // wrappers returning by reference (avoid array temporaries in inner loop)
  static inline void evalRB2D(double x,double y,double t,double& u,double& v) noexcept {
    const auto a = RisingBubble2D(x,y,t); u=a[0]; v=a[1];
  }
  static inline void evalRot2D(double x,double y,double t,double& u,double& v) noexcept {
    const auto a = RotationVel2D(x,y,t); u=a[0]; v=a[1];
  }
  static inline void evalVor2D(double x,double y,double t,double& u,double& v) noexcept {
    const auto a = VortexVel2D(x,y,t); u=a[0]; v=a[1];
  }

  static inline void evalRB3D(double x,double y,double z,double t,double& u,double& v,double& w) noexcept {
    const auto a = RisingBubble3D(x,y,z,t); u=a[0]; v=a[1]; w=a[2];
  }
  static inline void evalRot3D(double x,double y,double z,double t,double& u,double& v,double& w) noexcept {
    const auto a = RotationVel3D(x,y,z,t); u=a[0]; v=a[1]; w=a[2];
  }
  static inline void evalVor3D(double x,double y,double z,double t,double& u,double& v,double& w) noexcept {
    const auto a = VortexVel3D(x,y,z,t); u=a[0]; v=a[1]; w=a[2];
  }

  static inline VelEvalFn2D select2D(const VelKind kind)
  {
    switch (kind) {
      case VelKind::RisingBubble: return &evalRB2D;
      case VelKind::Rotation:     return &evalRot2D;
      case VelKind::Vortex:       return &evalVor2D;
      default: throw std::runtime_error("RungeKutta: unknown VelKind (2D)");
    }
  }

  static inline VelEvalFn3D select3D(const VelKind kind)
  {
    switch (kind) {
      case VelKind::RisingBubble: return &evalRB3D;
      case VelKind::Rotation:     return &evalRot3D;
      case VelKind::Vortex:       return &evalVor3D;
      default: throw std::runtime_error("RungeKutta: unknown VelKind (3D)");
    }
  }

  static inline void rk4_inplace(std::vector<std::vector<double>>& Xp,
                                 const double t0,
                                 const double dt,
                                 const VelKind kind)
  {
    const std::size_t dim = Xp.size();
    if (dim != 2 && dim != 3) throw std::runtime_error("RungeKutta::rk4_inplace: Xp.size() must be 2 or 3");
    if (Xp[0].empty()) return;

    const std::size_t nPts = Xp[0].size();
    for (std::size_t a = 1; a < dim; ++a) {
      if (Xp[a].size() != nPts) throw std::runtime_error("RungeKutta::rk4_inplace: inconsistent Xp[a].size()");
    }

    // raw pointers (speed)
    double* FEMUS_RESTRICT x = Xp[0].data();
    double* FEMUS_RESTRICT y = Xp[1].data();
    double* FEMUS_RESTRICT z = (dim == 3) ? Xp[2].data() : nullptr;

    const double halfdt = 0.5 * dt;
    const double sixth  = dt / 6.0;

    if (dim == 2) {
      const VelEvalFn2D vel = select2D(kind);

      for (std::size_t i = 0; i < nPts; ++i) {
        const double x0 = x[i];
        const double y0 = y[i];

        double k1x, k1y;
        double k2x, k2y;
        double k3x, k3y;
        double k4x, k4y;

        vel(x0, y0, t0, k1x, k1y);

        vel(x0 + halfdt*k1x, y0 + halfdt*k1y, t0 + halfdt, k2x, k2y);
        vel(x0 + halfdt*k2x, y0 + halfdt*k2y, t0 + halfdt, k3x, k3y);
        vel(x0 + dt*k3x,     y0 + dt*k3y,     t0 + dt,     k4x, k4y);

        x[i] = x0 + sixth * (k1x + 2.0*k2x + 2.0*k3x + k4x);
        y[i] = y0 + sixth * (k1y + 2.0*k2y + 2.0*k3y + k4y);
      }
      return;
    }

    // dim == 3
    const VelEvalFn3D vel = select3D(kind);

    for (std::size_t i = 0; i < nPts; ++i) {
      const double x0 = x[i];
      const double y0 = y[i];
      const double z0 = z[i];

      double k1x, k1y, k1z;
      double k2x, k2y, k2z;
      double k3x, k3y, k3z;
      double k4x, k4y, k4z;

      vel(x0, y0, z0, t0, k1x, k1y, k1z);

      vel(x0 + halfdt*k1x, y0 + halfdt*k1y, z0 + halfdt*k1z, t0 + halfdt, k2x, k2y, k2z);
      vel(x0 + halfdt*k2x, y0 + halfdt*k2y, z0 + halfdt*k2z, t0 + halfdt, k3x, k3y, k3z);
      vel(x0 + dt*k3x,     y0 + dt*k3y,     z0 + dt*k3z,     t0 + dt,     k4x, k4y, k4z);

      x[i] = x0 + sixth * (k1x + 2.0*k2x + 2.0*k3x + k4x);
      y[i] = y0 + sixth * (k1y + 2.0*k2y + 2.0*k3y + k4y);
      z[i] = z0 + sixth * (k1z + 2.0*k2z + 2.0*k3z + k4z);
    }
  }
};

#undef FEMUS_RESTRICT
