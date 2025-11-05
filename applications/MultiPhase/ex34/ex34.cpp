
// main.cpp
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <variant>
#include <vector>
#include <gperftools/profiler.h>


#include "Mollifier.hpp"
#include "Reinitializer.hpp"
#include "OctTree.hpp"
#include "FieldAdvection.hpp"
#include "FieldAdvectionWithAnalyticVelocity.hpp"

using namespace fem;

template<std::size_t DIM> using Pt = std::array<double, DIM>;

// ---------------- Level sets ----------------
template<std::size_t DIM> struct Psi;

template<> struct Psi<2> {
  double xc{0.0}, yc{0.0};
  double sigma{0.18016836131796748}, r{0.15}, delta{0.0}, eps{1. / 512};
  double operator()(double x, double y) const {
    // const double rr = (x - xc) * (x - xc) + (y - yc) * (y - yc);
    // return std::exp((r * r - rr) / (sigma * sigma)) - 1.0 + delta;
    Mollifier m(eps);
    const double d = r - sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
    return m.SigmoidC1(d);
  }
};

template<> struct Psi<3> {
  double xc{0.0}, yc{0.0}, zc{0.0};
  double sigma{0.18016836131796748}, r{0.15}, delta{0.0}, eps{1. / 128};
  double operator()(double x, double y, double z) const {
    const double dx = x - xc, dy = y - yc, dz = z - zc;
    //const double rr = dx * dx + dy * dy + dz * dz;
    //return std::exp((r * r - rr) / (sigma * sigma)) - 1.0 + delta;

    Mollifier m(eps);
    const double d = r - sqrt(dx * dx + dy * dy + dz * dz);
    return m.SigmoidC1(d);

  }
};

// ---------------- Analytic velocities ----------------

// Rising-bubble-like incompressible velocity field in [-Lx/2,Lx/2] x [0,Ly]
static inline std::array<double, 2>
RisingBubble2D(double x, double y, double time) {
  double Lx = 1, Ly = 1, A = .1;
  const double pi = M_PI, T = 8.;
  const double cost = std::cos(M_PI * time / T);
  const double s1 = std::sin(pi * (x + 0.5 * Lx) / Lx);
  const double X  = (s1 * s1) * std::sin(2.0 * pi * x / Lx);
  const double Y  = std::pow(std::sin(pi * y / Ly), 2);
  const double dYdy = 2.0 * std::sin(pi * y / Ly) * std::cos(pi * y / Ly) * (pi / Ly);
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

// Rising-bubble-like incompressible velocity field (3D)
// Domain: x∈[-Lx/2,Lx/2], y∈[-Ly/2,Ly/2], z∈[0,Lz]
static inline std::array<double, 3>
RisingBubble3D(double x, double y, double z, double time) {
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
  const double Z  = sz * sz, Zp = 2.0 * sz * cz * az;

  const double u = -A * X * Zp * cost;
  const double v = -A * Y * Zp * cost;
  const double w =  A * (Xp + Yp) * Z * cost;
  return {u, v, w};
}

// Rotational velocities (your original)
static inline std::array<double, 2>
RotationVel2D(double x, double y, double time) noexcept {
  return {y, -x};
}

static inline std::array<double, 3>
RotationVel3D(double x, double y, double z, double time) noexcept {
  return {y - z, z - x, x - y};
}

static inline std::array<double, 2>
VortexVel2D(double x, double y, double time) noexcept {
  const double T = 8.0; x += 0.5; y += 0.5;
  double u = -2.*std::sin(M_PI * x) * std::sin(M_PI * x)
             * std::sin(M_PI * y) * std::cos(M_PI * y)
             * std::cos(M_PI * time / T);
  double v =  2.*std::sin(M_PI * x) * std::cos(M_PI * x)
              * std::sin(M_PI * y) * std::sin(M_PI * y)
              * std::cos(M_PI * time / T);
  return {u, v};
}

static inline std::array<double, 3>
VortexVel3D(double x, double y, double z, double time) noexcept {
  x += 0.5; y += 0.5; z += 0.5;
  double sx, cx, sy, cy, sz, cz;
  ::sincos(M_PI * x, &sx, &cx);
  ::sincos(M_PI * y, &sy, &cy);
  ::sincos(M_PI * z, &sz, &cz);
  const double cosT = std::cos(M_PI * (time * 0.25)); // your original factor
  const double s2x = sx * sx, s2y = sy * sy, s2z = sz * sz;
  const double sxcx = sx * cx, sycy = sy * cy, szcz = sz * cz;
  const double k = 2.0 * cosT;
  return { k*s2x*(-sycy + szcz), k*s2y*(-szcz + sxcx), k*s2z*(-sxcx + sycy) };
}


// ---------------- Geometry builders ----------------
template<std::size_t DIM> struct Box;

// 2D
template<> struct Box<2> {
  static std::array<std::array<double, 9>, 2>
  coords(const std::array<double, 2>& minCorner = {-0.5, -0.5},
         const std::array<double, 2>& maxCorner = {+0.5, +0.5}) {
    std::array<std::array<double, 9>, 2> X{};
    const double ax = minCorner[0], bx = maxCorner[0];
    const double ay = minCorner[1], by = maxCorner[1];
    const double xm = 0.5 * (ax + bx), ym = 0.5 * (ay + by);
    X[0][0] = ax; X[1][0] = ay;  X[0][1] = bx; X[1][1] = ay;
    X[0][2] = bx; X[1][2] = by;  X[0][3] = ax; X[1][3] = by;
    X[0][4] = xm; X[1][4] = ay;  X[0][5] = bx; X[1][5] = ym;
    X[0][6] = xm; X[1][6] = by;  X[0][7] = ax; X[1][7] = ym;
    X[0][8] = xm; X[1][8] = ym;
    return X;
  }
};
// 3D
template<> struct Box<3> {
  static std::array<std::array<double, 27>, 3>
  coords(const std::array<double, 3>& minCorner = {-0.5, -0.5, -0.5},
         const std::array<double, 3>& maxCorner = {+0.5, +0.5, +0.5}) {
    std::array<std::array<double, 27>, 3> X{};
    const double ax = minCorner[0], bx = maxCorner[0];
    const double ay = minCorner[1], by = maxCorner[1];
    const double az = minCorner[2], bz = maxCorner[2];
    const double xm = 0.5 * (ax + bx), ym = 0.5 * (ay + by), zm = 0.5 * (az + bz);
    auto set = [&](int i, double x, double y, double z) {
      X[0][i] = x;
      X[1][i] = y;
      X[2][i] = z;
    };
    set(0, ax, ay, az); set(1, bx, ay, az); set(2, bx, by, az); set(3, ax, by, az);
    set(4, ax, ay, bz); set(5, bx, ay, bz); set(6, bx, by, bz); set(7, ax, by, bz);
    set(8, xm, ay, az); set(9, bx, ym, az); set(10, xm, by, az); set(11, ax, ym, az);
    set(12, xm, ay, bz); set(13, bx, ym, bz); set(14, xm, by, bz); set(15, ax, ym, bz);
    set(16, ax, ay, zm); set(17, bx, ay, zm); set(18, bx, by, zm); set(19, ax, by, zm);
    set(20, xm, ay, zm); set(21, bx, ym, zm); set(22, xm, by, zm); set(23, ax, ym, zm);
    set(24, xm, ym, az); set(25, xm, ym, bz); set(26, xm, ym, zm);
    return X;
  }
};

template<std::size_t DIM> struct Ball;
// 2D circle

template<> struct Ball<2> {
  static std::array<std::array<double, 9>, 2>
  data(const double R = 0.5, const double cx = 0.0, const double cy = 0.0) {
    const double inv_sqrt2 = 0.70710678118654752440084436210485; // 1/sqrt(2)
    const double c = R * inv_sqrt2;

    std::array<std::array<double, 9>, 2> XY{{
        // x row [0][j]
        { -c,  +c,  +c,  -c,  0.0, +R,  0.0, -R,  0.0 },
        // y row [1][j]
        { -c,  -c,  +c,  +c, -R,   0.0, +R,   0.0, 0.0 }
      }};

    // shift by center (cx, cy)
    for (int j = 0; j < 9; ++j) {
      XY[0][j] += cx;
      XY[1][j] += cy;
    }
    return XY;
  }
};

//3D sphere
template<> struct Ball<3> {
  static std::array<std::array<double, 27>, 3>
  data(const double R = 0.5, const double cx = 0.0, const double cy = 0.0, const double cz = 0.0) {
    const double inv_sqrt3 = 0.57735026918962576451; // 1/sqrt(3)
    const double inv_sqrt2 = 0.70710678118654752440; // 1/sqrt(2)
    const double a = R * inv_sqrt3; // corners on diagonals
    const double b = R * inv_sqrt2; // edge midpoints

    std::array<std::array<double, 27>, 3> XYZ{};

    // ---- corners (0..7) ----
    // Morton-like: z-, then z+
    const int sgn[8][3] = {
      {-1, -1, -1}, {+1, -1, -1}, {+1, +1, -1}, {-1, +1, -1},
      {-1, -1, +1}, {+1, -1, +1}, {+1, +1, +1}, {-1, +1, +1}
    };
    for (int k = 0; k < 8; ++k) {
      XYZ[0][k] = a * sgn[k][0];
      XYZ[1][k] = a * sgn[k][1];
      XYZ[2][k] = a * sgn[k][2];
    }

    // ---- z = -1 edges (8..11): bottom ring ----
    // order: 0-1, 1-2, 2-3, 3-0
    XYZ[0][ 8] =  0.0; XYZ[1][ 8] = -b;  XYZ[2][ 8] = -b;  // between (-a,-a,-a) and (+a,-a,-a)
    XYZ[0][ 9] = +b;  XYZ[1][ 9] =  0.0; XYZ[2][ 9] = -b;  // between (+a,-a,-a) and (+a,+a,-a)
    XYZ[0][10] = 0.0; XYZ[1][10] = +b;  XYZ[2][10] = -b;  // between (+a,+a,-a) and (-a,+a,-a)
    XYZ[0][11] = -b;  XYZ[1][11] = 0.0; XYZ[2][11] = -b;  // between (-a,+a,-a) and (-a,-a,-a)

    // ---- z = +1 edges (12..15): top ring ----
    // order: 4-5, 5-6, 6-7, 7-4
    XYZ[0][12] =  0.0; XYZ[1][12] = -b;  XYZ[2][12] = +b;  // between (-a,-a,+a) and (+a,-a,+a)
    XYZ[0][13] = +b;  XYZ[1][13] =  0.0; XYZ[2][13] = +b;  // between (+a,-a,+a) and (+a,+a,+a)
    XYZ[0][14] =  0.0; XYZ[1][14] = +b;  XYZ[2][14] = +b;  // between (+a,+a,+a) and (-a,+a,+a)
    XYZ[0][15] = -b;  XYZ[1][15] =  0.0; XYZ[2][15] = +b;  // between (-a,+a,+a) and (-a,-a,+a)

    // ---- vertical edges (16..19): 0-4, 1-5, 2-6, 3-7 ----
    XYZ[0][16] = -b;  XYZ[1][16] = -b;  XYZ[2][16] = 0.0;  // 0-4
    XYZ[0][17] = +b;  XYZ[1][17] = -b;  XYZ[2][17] = 0.0;  // 1-5
    XYZ[0][18] = +b;  XYZ[1][18] = +b;  XYZ[2][18] = 0.0;  // 2-6
    XYZ[0][19] = -b;  XYZ[1][19] = +b;  XYZ[2][19] = 0.0;  // 3-7

    // ---- face centers (20..25): (y=-1, x=+1, y=+1, x=-1, z=-1, z=+1) ----
    XYZ[0][20] =  0.0; XYZ[1][20] = -R;  XYZ[2][20] =  0.0; // y = -1 (−Y)
    XYZ[0][21] = +R;   XYZ[1][21] =  0.0; XYZ[2][21] =  0.0; // x = +1 (+X)
    XYZ[0][22] =  0.0; XYZ[1][22] = +R;  XYZ[2][22] =  0.0; // y = +1 (+Y)
    XYZ[0][23] = -R;   XYZ[1][23] =  0.0; XYZ[2][23] =  0.0; // x = -1 (−X)
    XYZ[0][24] =  0.0; XYZ[1][24] =  0.0; XYZ[2][24] = -R;  // z = -1 (bottom)
    XYZ[0][25] =  0.0; XYZ[1][25] =  0.0; XYZ[2][25] = +R;  // z = +1 (top)

    // ---- cell center (26) ----
    XYZ[0][26] = 0.0; XYZ[1][26] = 0.0; XYZ[2][26] = 0.0;

    // translate by (cx, cy, cz)
    for (int j = 0; j < 27; ++j) {
      XYZ[0][j] += cx;
      XYZ[1][j] += cy;
      XYZ[2][j] += cz;
    }

    return XYZ;
  }
};


// ---------------- Funnel geometry builder ----------------
template<std::size_t DIM> struct Funnel;

// 2D specialization
template<>
struct Funnel<2> {
  static std::array<std::array<double, 9>, 2>
  data(double a, double b, double h, double k, double y0, double y1) {
    std::array<std::array<double, 9>, 2> X{};

    // Compute y mid
    const double ym = 0.5 * (y0 + y1);

    // Compute x positions from hyperbola for left/right boundaries
    auto x_from_y = [&](double y, double sign) {
      const double term = 1.0 + std::pow((y - k) / b, 2);
      return h + sign * a * std::sqrt(term);
    };

    // Left (−) and right (+) x at bottom, mid, top
    const double xl0 = x_from_y(y0, -1.0);
    const double xr0 = x_from_y(y0, +1.0);
    const double xlm = x_from_y(ym, -1.0);
    const double xrm = x_from_y(ym, +1.0);
    const double xl1 = x_from_y(y1, -1.0);
    const double xr1 = x_from_y(y1, +1.0);

    // Fill Q9 nodes (same order as Box<2>)
    X[0][0] = xl0;  X[1][0] = y0;  // bottom-left
    X[0][1] = xr0;  X[1][1] = y0;  // bottom-right
    X[0][2] = xr1;  X[1][2] = y1;  // top-right
    X[0][3] = xl1;  X[1][3] = y1;  // top-left

    X[0][4] = h;     X[1][4] = y0;  // mid-bottom
    X[0][5] = xrm;   X[1][5] = ym;  // mid-right
    X[0][6] = h;     X[1][6] = y1;  // mid-top
    X[0][7] = xlm;   X[1][7] = ym;  // mid-left

    X[0][8] = 0.5 * (xlm + xrm);  X[1][8] = ym;  // center (approx average midsection)

    return X;
  }
};

// 3D specialization -----------------------------------------------------------
template<>
struct Funnel<3> {
  static std::array<std::array<double, 27>, 3>
  data(double a, double b, double c,
       double h, double k, double l,
       double z0, double z1) {
    std::array<std::array<double, 27>, 3> X{};

    // --- Midplane ---
    const double zm = 0.5 * (z0 + z1);

    // --- Hyperbolic scale factor ---
    auto scale = [&](double z) {
      return std::sqrt(1.0 + std::pow((z - l) / c, 2));
    };

    const double s0 = scale(z0);
    const double sm = scale(zm);
    const double s1 = scale(z1);

    // --- Local semi-axes for each z plane ---
    const double ax0 = a * s0, ay0 = b * s0;
    const double axm = a * sm, aym = b * sm;
    const double ax1 = a * s1, ay1 = b * s1;

    auto set = [&](int i, double x, double y, double z) {
      X[0][i] = x;  X[1][i] = y;  X[2][i] = z;
    };

    // --- Corner nodes (bottom z=z0, top z=z1) -------------------------------
    set(0, h - ax0, k - ay0, z0);
    set(1, h + ax0, k - ay0, z0);
    set(2, h + ax0, k + ay0, z0);
    set(3, h - ax0, k + ay0, z0);

    set(4, h - ax1, k - ay1, z1);
    set(5, h + ax1, k - ay1, z1);
    set(6, h + ax1, k + ay1, z1);
    set(7, h - ax1, k + ay1, z1);

    // --- Midpoints on bottom face (z=z0) -----------------------------------
    set(8,  h,       k - ay0, z0); // front
    set(9,  h + ax0, k,       z0); // right
    set(10, h,       k + ay0, z0); // back
    set(11, h - ax0, k,       z0); // left

    // --- Midpoints on top face (z=z1) --------------------------------------
    set(12, h,       k - ay1, z1); // front
    set(13, h + ax1, k,       z1); // right
    set(14, h,       k + ay1, z1); // back
    set(15, h - ax1, k,       z1); // left

    // --- Corners on middle z plane (z=zm) ----------------------------------
    set(16, h - axm, k - aym, zm);
    set(17, h + axm, k - aym, zm);
    set(18, h + axm, k + aym, zm);
    set(19, h - axm, k + aym, zm);

    // --- Mid-edge nodes on middle layer (front→right→back→left) ------------
    set(20, h,       k - aym, zm); // front (y=-1)
    set(21, h + axm, k,       zm); // right (x=+1)
    set(22, h,       k + aym, zm); // back  (y=+1)
    set(23, h - axm, k,       zm); // left  (x=-1)

    // --- Face centers + element center -------------------------------------
    set(24, h, k, z0);  // bottom center
    set(25, h, k, z1);  // top center
    set(26, h, k, zm);  // mid center

    return X;
  }
};


// ---------------- Dimension traits (as you had) ----------------
template<std::size_t DIM> struct DimOps;

template<> struct DimOps<2> {
  typedef BasisT<2> BasisType;
  static BasisType highest_basis() {
    return BasisT<2>::Q9;
  }
  static void shift(Psi<2>& psi) {
    psi.xc = 0.0;
    psi.yc = 0.25;
  }
  static double evalPsi(const Psi<2>& psi, const Pt<2>& p) {
    return psi(p[0], p[1]);
  }
  static std::array<std::array<double, 9>, 2> box_geom(
    const std::array<double, 2>& minCorner = {-0.5, -0.5},
    const std::array<double, 2>& maxCorner = {+0.5, +0.5}) {
    return Box<2>::coords(minCorner, maxCorner);
  }
  static std::array<std::array<double, 9>, 2> ball_geom(
    const double R = 0.5, const std::array<double, 2> xc = {0., 0.}) {
    return Ball<2>::data(R, xc[0], xc[1]);
  }

  static std::array<std::array<double, 9>, 2> funnel_geom(
    const double a = 0.375, const double b = 0.5,
    const std::array<double, 2> xc = {0., 0.5},
    const double y0 = 0, const double y1 = 1) {
    return Funnel<2>::data(a, b, xc[0], xc[1], y0,  y1);
  }

};

template<> struct DimOps<3> {
  typedef BasisT<3> BasisType;
  static BasisType highest_basis() {
    return BasisT<3>::H27;
  }
  static void shift(Psi<3>& psi) {
    psi.xc = 0.0;
    psi.yc = 0.0;
    psi.zc = 0.25;
  }
  static double evalPsi(const Psi<3>& psi, const Pt<3>& p) {
    return psi(p[0], p[1], p[2]);
  }
  static std::array<std::array<double, 27>, 3> box_geom(
    const std::array<double, 3>& minCorner = {-0.5, -0.5, -0.5},
    const std::array<double, 3>& maxCorner = {+0.5, +0.5, +0.5}) {
    return Box<3>::coords(minCorner, maxCorner);
  }
  static std::array<std::array<double, 27>, 3> ball_geom(
    const double R = 0.5, const std::array<double, 3> xc = {0., 0., 0.}) {
    return Ball<3>::data(R, xc[0], xc[1], xc[2]);
  }

  static std::array<std::array<double, 27>, 3> funnel_geom(
    const double a = 0.375, const double b = 0.4, const double c = 0.5,
    const std::array<double, 3> xc = {0., 0., 0.5},
    const double z0 = 0, const double z1 = 1) {
    return Funnel<3>::data(a, b, c, xc[0], xc[1], xc[2], z0,  z1);
  }
};

// ---------------- Scenarios ----------------
enum class Scenario { RB2D, RB3D, VX2D, VX3D, ROT2D, ROT3D };
static inline Scenario parse_scenario(const std::string& s) {
  if (s == "RB2D"  || s == "rb2d")  return Scenario::RB2D;
  if (s == "RB3D"  || s == "rb3d")  return Scenario::RB3D;
  if (s == "VX2D" || s == "vx2d") return Scenario::VX2D;
  if (s == "VX3D" || s == "vx3d") return Scenario::VX3D;
  if (s == "ROT2D" || s == "rot2d") return Scenario::ROT2D;
  if (s == "ROT3D" || s == "rot3d") return Scenario::ROT3D;
  return Scenario::VX2D;
}


// ======================= Helpers for DIM-dependent types =======================

// Velocity pointer types (match your advectors: scalar args)
template<std::size_t DIM> struct VelT;
template<> struct VelT<2> {
  using type = std::array<double, 2> (*)(double, double, double);
};
template<> struct VelT<3> {
  using type = std::array<double, 3> (*)(double, double, double, double);
};
template<std::size_t DIM> using Vel = typename VelT<DIM>::type;

// Shape-cache types used by OctTree refine/coarsen
template<std::size_t DIM> struct ShapeCacheT;
template<> struct ShapeCacheT<2> {
  using type = std::vector<std::array<double, 9>>;
};
template<> struct ShapeCacheT<3> {
  using type = std::vector<std::array<double, 27>>;
};
template<std::size_t DIM> using ShapeCache = typename ShapeCacheT<DIM>::type;

// Domain per (DIM, Scenario)
template<std::size_t DIM>
static inline std::pair<std::array<double, DIM>, std::array<double, DIM>>
make_domain(Scenario s) {
  if constexpr(DIM == 2) {
    switch (s) {
    case Scenario::RB2D: return { { -0.5,  0.0 }, { +0.5, 1.0 } };
    case Scenario::VX2D: return { { -0.5, -0.5 }, { +0.5, +0.5 } };
    case Scenario::ROT2D: return { { -0.5, -0.5 }, { +0.5, +0.5 } };
    default:             return { { -0.5, -0.5 }, { +0.5, +0.5 } };
    }
  }
  else {
    switch (s) {
    case Scenario::RB3D: return { { -0.5, -0.5, 0.0 }, { +0.5, +0.5, 1.0 } };
    case Scenario::VX3D: return { { -0.5, -0.5, -0.5 }, { +0.5, +0.5, +0.5 } };
    case Scenario::ROT3D: return { { -0.5, -0.5, -0.5 }, { +0.5, +0.5, +0.5 } };
    default:             return { { -0.5, -0.5, -0.5 }, { +0.5, +0.5, +0.5 } };
    }
  }
}

// Velocity selector (function pointer with scalar args)
template<std::size_t DIM>
static inline Vel<DIM> make_velocity(Scenario s) {
  if constexpr(DIM == 2) {
    if (s == Scenario::RB2D) return &RisingBubble2D;
    else if (s == Scenario::VX2D) return &VortexVel2D;
    else return &RotationVel2D; // default ROT2D
  }
  else {
    if (s == Scenario::RB3D) return &RisingBubble3D;
    else if (s == Scenario::VX3D) return &VortexVel3D;
    else return &RotationVel3D; // default ROT3D
  }
}

// Reconcile scenario to DIM (keeps CLI flexible)
template<std::size_t DIM>
static inline Scenario reconcile_dim(Scenario s) {
  if constexpr(DIM == 2) {
    if (s == Scenario::RB3D)  return Scenario::RB2D;
    if (s == Scenario::VX3D) return Scenario::VX2D;
    if (s == Scenario::ROT3D) return Scenario::ROT2D;
    return s;
  }
  else {
    if (s == Scenario::RB2D)  return Scenario::RB3D;
    if (s == Scenario::VX2D) return Scenario::VX3D;
    if (s == Scenario::ROT2D) return Scenario::ROT3D;
    return s;
  }
}

template<std::size_t DIM>
static inline double period_for(Scenario s) {
  if constexpr(DIM == 2) {
    if (s == Scenario::ROT2D) return 2. * M_PI;
    else return 8.0;
  }
  else {
    if (s == Scenario::ROT3D) return 2. * M_PI / sqrt(3.);
    else return 4.0;
  }
}

// DIM-generic VTU vector write (phi must already exist as field 0)
template<std::size_t DIM>
static inline void write_vtu_frame_generic(OctTree<DIM>& ot,
                                           unsigned max_depth, unsigned k, double time,
                                           Vel<DIM> evalVelocity,
                                           BasisT<DIM> bHi,
                                           const char* stem) {
  static const std::array<const char*, 3> names = {"u", "v", "w"};
  std::array<u32, DIM> fidv{};
  for (size_t d = 0; d < DIM; ++d) fidv[d] = ot.add_field(bHi, names[d]);
  ot.resize_fields_to_nodes();

  std::array<Field*, DIM> f{};
  for (size_t d = 0; d < DIM; ++d) f[d] = &ot.field(fidv[d]);

  const BasisT<DIM> bs = to_basis<DIM>(f[0]->basis_id);
  const auto& nodes = ot.basis_nodes(bs);
  for (size_t d = 0; d < DIM; ++d) f[d]->nodal.resize(nodes.size());

  for (size_t gid = 0; gid < nodes.size(); ++gid) {
    const auto& p = nodes[gid].physical;
    std::array<double, DIM> v{};
    if constexpr(DIM == 2) v = evalVelocity(p[0], p[1], time);
    else                    v = evalVelocity(p[0], p[1], p[2], time);
    for (size_t d = 0; d < DIM; ++d) f[d]->nodal[gid] = v[d];
  }

  std::vector<std::variant<u32, std::vector<u32>>> groups;
  groups.emplace_back(0u); // phi
  groups.emplace_back(std::vector<u32>(fidv.begin(), fidv.begin() + DIM)); // velocity

  const std::string filename =
    std::string("./output/") + stem + "." + std::to_string(max_depth) + "." + std::to_string(k) + ".vtu";
  ot.write_binary_vtu(filename, bHi, groups, false);
}

// ======================== Unified run<DIM> (exact signatures) ========================
template<std::size_t DIM>
std::pair<double, double> run(int /*argc*/, char** /*argv*/, unsigned nSteps, Scenario scenario, bool vtu, bool pprof, const u32 max_depth, const bool box_domain, const unsigned reinit) {
  if (pprof) ProfilerStart((DIM == 3) ? "profiling_3d.prof" : "profiling_2d.prof");

  const u32 maxDepth = max_depth;
  const u32 minDepth = (DIM == 3 ? 2u : 3u);
  const bool allowDrop = true;

  std::cout << (DIM == 3 ? "Octree 3D" : "Octree 2D")
            << " config: maxDepth=" << maxDepth
            << " minDepth=" << minDepth
            << " allowDrop=" << (allowDrop ? 1 : 0) << "\n";

  // Select scenario/domain/velocity per DIM
  scenario = reconcile_dim<DIM>(scenario);

  Vel<DIM> evalVelocity = make_velocity<DIM>(scenario);

  OctTree<DIM> ot(maxDepth, minDepth);
  ot.set_allow_coarsen_below_min(allowDrop);

  const auto [minCorner, maxCorner] = make_domain<DIM>(scenario);
  auto X = (box_domain) ?
           DimOps<DIM>::box_geom(minCorner, maxCorner) :
           ((scenario == Scenario::RB2D || scenario == Scenario::RB3D) ? DimOps<DIM>::funnel_geom() : DimOps<DIM>::ball_geom(0.65, fem::Point<DIM> {0}));

  ot.set_physical_coordinates(X);


  double eps = (DIM == 2) ? 1. / pow(2, std::max(maxDepth - 7u, 1u)) : 1. / pow(2, std::max(maxDepth - 4u, 1u)); //this parameter is fundamental, if I make it too small the convergence breaks.
  //double eps = 1. / pow(2, 5);

  // Level set
  Psi<DIM> psi;
  psi.eps = eps;

  DimOps<DIM>::shift(psi);

  const double tau_ref    = 2.0;
  const double tau_coarse = 1e-5;

  const typename DimOps<DIM>::BasisType bHi = DimOps<DIM>::highest_basis();

  // ---- EXACT signatures as in your original ----
  // refine:  (u32 idx, const vector<Pt<DIM>>& pts_s, const vector<Pt<DIM>>& pts_xyz, const ShapeCache<DIM>& cache)
  // coarsen: (u32 morton_or_id, u32 level, const vector<Pt<DIM>>& pts_s, const vector<Pt<DIM>>& pts_xyz, const ShapeCache<DIM>& cache)

  const auto refine1 = [&](u32 /*idx*/,
                           const std::vector<Pt<DIM>>& /*pts_s*/,
                           const std::vector<Pt<DIM>>& pts_xyz,
  const ShapeCache<DIM>& /*cache*/) -> bool {
    if (pts_xyz.empty()) return false;
    double v0 = DimOps<DIM>::evalPsi(psi, pts_xyz[0]), mn = v0, mx = v0;
    for (size_t i = 1; i < pts_xyz.size(); ++i) {
      double v = DimOps<DIM>::evalPsi(psi, pts_xyz[i]);
      mn = std::min(mn, v);
      mx = std::max(mx, v);
    }
    return (mn <= 0.0 && mx >= 0.0) || (mx - mn > tau_ref);
  };

  const auto coarsen1 = [&](u32 /*id_or_morton*/, u32 level,
                            const std::vector<Pt<DIM>>& /*pts_s*/,
                            const std::vector<Pt<DIM>>& pts_xyz,
  const ShapeCache<DIM>& /*cache*/) -> bool {
    if (level <= minDepth + 1u || level == maxDepth || pts_xyz.empty()) return false;
    double v0 = DimOps<DIM>::evalPsi(psi, pts_xyz[0]), mn = v0, mx = v0;
    for (size_t i = 1; i < pts_xyz.size(); ++i) {
      double v = DimOps<DIM>::evalPsi(psi, pts_xyz[i]);
      mn = std::min(mn, v);
      mx = std::max(mx, v);
    }
    return (mn > +tau_coarse) || (mx < -tau_coarse);
  };

  std::size_t changed = ot.adapt_cycle(coarsen1, refine1, /*max_passes=*/10);
  std::cout << "Adapt-cycle changed " << changed << " cells\n";
  std::cout << "Leaves after cycle: " << ot.leaf_count() << "\n";

  // Scalar field "phi"
  const u32 fid = ot.add_field(bHi, "phi");
  ot.resize_fields_to_nodes();
  {
    Field& fld = ot.field(fid);
    const auto bs = to_basis<DIM>(fld.basis_id);
    const auto& nodes = ot.basis_nodes(bs);
    fld.nodal.resize(nodes.size());
    for (size_t gid = 0; gid < nodes.size(); ++gid) {
      const auto& p = nodes[gid].physical;
      fld.nodal[gid] = DimOps<DIM>::evalPsi(psi, p);
    }
  }


  Mollifier m(eps);

  Reinitializer<DIM> reinitializer(&ot, fid, [&m](double x) noexcept {
    return m.SigmoidC1(x);
  }, true /*projection flag*/);

  // if (reinit) {
  //   reinitializer.compute_signed_distance();
  // }
  // First VTU frame
  if (vtu) {
    const char* stem = (DIM == 2 ? "element_adaptive2d" : "element_adaptive3d");
    write_vtu_frame_generic(ot, max_depth, 0, 0.0, evalVelocity, bHi, stem);
  }

  // Time-stepping (dt from fixed nIter=320)
  const unsigned nIter = 320;
  const double period = period_for<DIM>(scenario);
  const double dt = period / static_cast<double>(nIter);

  OctTree<DIM> ot0 = ot;
  double N_cells_sum = 0.0;
  double N_cells_max = 0.0;
  unsigned n_samples = 0;
  const u32 fid_t0 = fid;
  OctTree<DIM> ot1(maxDepth, minDepth);

  std::vector<Pt<DIM>> markers;
  if (reinit){
    markers.clear();
    markers = reinitializer.collect_markers(0. /*marker density*/, 3 /*min segments*/);
  }
  
  for (u32 k = 1; k <= nSteps; ++k) {
    const double time = k * dt;

    std::vector<Pt<DIM>> leftOld, stayedNew;
    if (reinit) fem::advect_physical_markers_forward_analytic(ot, time, dt, evalVelocity, markers, leftOld, stayedNew);
    else //TODO (check Samuele New Advection Functions)
    fem::advect_interface_markers_forward_analytic(ot, fid, time, dt, evalVelocity, leftOld, stayedNew);

    ot1.reset(false, false);
    ot1.set_allow_coarsen_below_min(allowDrop);
    ot1.set_physical_coordinates(X);
    ot1.refine_to_contain_points(stayedNew, ot1.max_depth());
    ot1.enlarge_top_layer_and_enforce_balance();

    const u32 fid0 = fid;
    const u32 fid1 = ot1.add_field(bHi, "phi");

    fem::advect_nodes_backward_and_transport_field_analytic(ot, fid0, time, dt, evalVelocity, ot1, fid1);

    u32 num_coarsened = 0;
    if (reinit) {
      //num_coarsened = ot1.coarsen_only_cycle_safe(fid0, {tau_coarse, 1.0e-6}, ot);
    }
    else {
      //num_coarsened = ot1.coarsen_only_cycle_safe(fid0, {tau_coarse}, ot);
    }
    std::cout << "\x1b[1A" << "\x1b[2K"   // cursor up 1, and erase entire line
              << "\x1b[1A" << "\x1b[2K"   // cursor up 1, and erase entire line
              << "\r"        // return to column 1
              << std::flush;
    std::cout << "Coarsened " << num_coarsened << " leaves.\n";

    using std::swap; swap(ot, ot1);

    if (reinit) {
      int marker_density = (k % reinit == 0) ? 20. : 0.;
      markers.clear();
      markers = reinitializer.collect_markers(marker_density, 3 /*min segments*/);
    }

    // --- Statistiche sulle celle ---
    const double n_cells = static_cast<double>(ot.leaf_count());
    N_cells_sum += n_cells;
    N_cells_max = std::max(N_cells_max, n_cells);
    n_samples++;


    if (reinit && (k % reinit == 0)) {
      reinitializer.compute_signed_distance();
    }

    if (vtu) {
      const char* stem = (DIM == 2 ? "element_adaptive2d" : "element_adaptive3d");
      write_vtu_frame_generic(ot, max_depth, k, time, evalVelocity, bHi, stem);
      // reinitializer.write_markers_csv(
      //   (DIM == 2 ? "./output_markers/markers2d." : "./output_markers/markers3d.") +
      //   std::to_string(max_depth) + "." + std::to_string(k) + ".csv");
    }

  }

  // --- Calcolo di h_min ---
  const double domain_length = maxCorner[0] - minCorner[0];  // vale 1.0 per te
  const double h_min = domain_length / std::pow(2.0, static_cast<double>(maxDepth));
  const double N_cells_mean = (n_samples > 0) ? (N_cells_sum / n_samples) : 0.0;

  std::cout << "Statistiche celle: max=" << N_cells_max
            << "  mean=" << N_cells_mean
            << "  h_min=" << h_min << "\n";

  // Union mesh VTU
  OctTree<DIM> ot3(ot0, 0, ot, 0);
  if (vtu) {
    const std::string filename =
      std::string("./output/") + (DIM == 2 ? "element_adaptive_union2d.vtu"
                                  : "element_adaptive_union3d.vtu");
    ot3.write_binary_vtu(filename, bHi, false);
    std::cout << "Printing " << filename << "\n";
  }

// --- Compute integral of sign(phi) ---
  // double I = compute_sign_integral(ot3, 0);
  // double total = (DIM == 2 ? 4.0 : 8.0); // dominio [-1,1]^DIM
  // double V_neg = 0.5 * (total - I);
  // double V_pos = total - V_neg;
  // std::cout << "Integral of sign(phi): " << I << "\n";
  // std::cout << "Volume(phi<0): " << V_neg
  //           << "  Volume(phi>0): " << V_pos << "\n";

  // --- Compute relative mass and geometrical errors ---
  const u32 fid_t1 = fid;  // id del campo φ al tempo t1

  auto [Em, Eg] = compute_mass_and_geom_errors_overlay(ot3, ot0, fid_t0, ot, fid_t1);

  std::cout << "---------------------------------------------\n";
  std::cout << "Relative mass error  E_m = " << Em << "\n";
  std::cout << "Geometrical error     E_g = " << Eg << "\n";
  std::cout << "---------------------------------------------\n";

  std::string suffix = (DIM == 3 ? "_3D" : "_2D");
  std::ofstream fout("errors_vs_depth" + suffix + ".dat", std::ios::app);

  if (fout) {
    fout << std::setw(3) << max_depth << " "
         << std::scientific << std::setprecision(10)
         << std::setw(15) << Em << " "
         << std::setw(15) << Eg << " "
         << std::scientific << std::setprecision(3)
         << std::setw(14) << N_cells_max << " "
         << std::setw(14) << N_cells_mean << " "
         << std::setw(14) << h_min << "\n";
  }
  else {
    std::cerr << "⚠️ Error: can't open errors_vs_depth.dat for writing.\n";
  }

  if (pprof) ProfilerStop();
  return {Em, Eg};
}


// ---------------- CLI dispatcher ----------------
static void print_usage(const char* prog) {
  std::cerr
      << "Usage: " << prog << "\n"
      << "  [-d|--dim 2|3]\n"
      << "  [-n|--niter N]\n"
      << "  [--scenario RB2D|RB3D|VX2D|VX3D|ROT2D|ROT3D]\n"
      << "  [--max_depth M]\n"
      << "  [--profiling|--no-profiling]\n"
      << "  [--vtu-output|--no-vtu-output]\n"
      << "  [--box_domain|--curved_domain]\n"
      << "  [--reinit R]\n"
      << "\n"
      << "Options:\n"
      << "  -d, --dim           Dimension (2 or 3). Default: 2\n"
      << "  -n, --niter         Number of time steps to execute (nSteps). Default: 320\n"
      << "      --scenario      Flow/velocity setup. Default reconciled from dim; default is VX2D.\n"
      << "                      Valid: RB2D, RB3D, VX2D, VX3D, ROT2D, ROT3D\n"
      << "      --max_depth     Maximum tree depth (u32). Default: 8\n"
      << "      --profiling     Enable gperftools CPU profiling\n"
      << "      --no_profiling  Disable profiling (default)\n"
      << "      --vtu           Write VTU frames (default)\n"
      << "      --no_vtu        Do not write VTU\n"
      << "      --box_domain    Use box geometry (default)\n"
      << "      --ball_domain   Use curved domain (circle or hyperbola in 2D, sphere or funnel in 3D depending on the scenario)\n"
      << "      --reinit        Number of time steps beween reinitializations; deafult is 0, i.e., no reinitialization\n"
      << "\n"
      << "Notes:\n"
      << "  --niter controls the number of STEPS executed (nSteps),\n"
      << "  while dt is computed with a fixed nIter=320.\n";
}


int main(int argc, char** argv) {
  int dim = 2;            // default dimension
  unsigned nSteps = 320;  // default steps
  bool pprof = false;
  bool vtu = true;
  bool box_domain = true;
  unsigned reinit = 0; //(default false)
  Scenario scenario = Scenario::VX2D; // default; reconciled with dim inside run

  unsigned max_depth = 8;
  unsigned delta_depth = 1;

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "-d" || a == "--dim") {
      if (i + 1 >= argc) {
        print_usage(argv[0]);
        return 1;
      }
      dim = std::atoi(argv[++i]);
    }
    else if (a == "-n" || a == "--niter") {
      if (i + 1 >= argc) {
        print_usage(argv[0]);
        return 1;
      }
      long v = std::strtol(argv[++i], nullptr, 10);
      if (v < 0) {
        std::cerr << "Error: --niter must be non-negative\n";
        return 1;
      }
      nSteps = static_cast<unsigned>(v);
    }
    else if (a == "--scenario") {
      if (i + 1 >= argc) {
        print_usage(argv[0]);
        return 1;
      }
      scenario = parse_scenario(argv[++i]);
    }
    else if (a == "-h" || a == "--help") {
      print_usage(argv[0]); return 0;
    }
    else if (a == "--profiling") {
      pprof = true;
    }
    else if (a == "--no_profiling") {
      pprof = false;
    }
    else if (a == "--vtu") {
      vtu = true;
    }
    else if (a == "--no_vtu") {
      vtu = false;
    }
    else if (a == "--max_depth") {
      if (i + 1 >= argc) {
        print_usage(argv[0]);
        return 1;
      }
      max_depth = std::atoi(argv[++i]);
    }
    else if (a == "--delta_depth") {
      if (i + 1 >= argc) {
        print_usage(argv[0]);
        return 1;
      }
      delta_depth = std::atoi(argv[++i]);
    }
    else if (a == "--curved_domain") {
      box_domain = false;
    }
    else if (a == "--box_domain") {
      box_domain = true;
    }
    if (a == "--reinit") {
      if (i + 1 >= argc) {
        print_usage(argv[0]);
        return 1;
      }
      reinit = std::atoi(argv[++i]);
    }
  }

  std::vector<std::pair<double, double>> Er(delta_depth);
  std::vector<clock_t> Time(delta_depth);
  clock_t start_time;

  for (unsigned i = 0; i < delta_depth; ++i) {
    switch (dim) {
    case 2:
      start_time = clock();
      Er[i] = run<2>(argc, argv, nSteps, scenario, vtu, pprof, max_depth + i, box_domain, reinit);
      Time[i] = clock() - start_time;
      break;
    case 3:
      start_time = clock();
      Er[i] = run<3>(argc, argv, nSteps, scenario, vtu, pprof, max_depth + i, box_domain, reinit);
      Time[i] = clock() - start_time;
      break;
    default:
      std::cerr << "Error: DIM must be 2 or 3 (got " << dim << ")\n";
      print_usage(argv[0]);
      return 2;
    }
  }

  std::cout << "Max_Depth\tMass_Error\tGeometric_Error\tCompt_Time(s)" << std::endl;
  for (unsigned i = 0; i < delta_depth; i++) {
    std::cout << max_depth + i << "\t\t" << Er[i].first << "\t" << Er[i].second << "\t"
              /*      */ << static_cast<double>(Time[i]) / CLOCKS_PER_SEC << std::endl;
    if (i + 1 < delta_depth) {
      std::cout << "conv.\t\t" << log(Er[i].first / Er[i + 1].first) / log(2.) << "\t\t"
                /*      */ << log(Er[i].second / Er[i + 1].second) / log(2.) << "\t\t"
                /*      */ << log((double)Time[i + 1] / (double)Time[i]) / log(2.) << std::endl;
    }
  }
  if (delta_depth > 1) {
    std::cout << "\naver. conv. \t" << log(Er[0].first / Er[delta_depth - 1].first) / ((delta_depth - 1) * log(2.))
              << "\t\t" << log(Er[0].second / Er[delta_depth - 1].second) / ((delta_depth - 1) * log(2.))
              << "\t\t" << log((double)Time[delta_depth - 1] / (double)Time[0]) / ((delta_depth - 1) * log(2.)) << std::endl;
  }
  return 0;

}

