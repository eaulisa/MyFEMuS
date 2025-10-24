
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

#include "OctTree.hpp"
#include "FieldAdvection3D.hpp"
#include "FieldAdvectionWithAnalyticVelocity3D.hpp"

using namespace fem;

template<std::size_t DIM> using Pt = std::array<double, DIM>;

// ---------------- Level sets ----------------
template<std::size_t DIM> struct Psi;

template<> struct Psi<2> {
  double xc{0.0}, yc{0.0};
  double sigma{0.18016836131796748}, r{0.15}, delta{0.0};
  double operator()(double x, double y) const {
    const double rr = (x - xc) * (x - xc) + (y - yc) * (y - yc);
    return std::exp((r * r - rr) / (sigma * sigma)) - 1.0 + delta;
  }
};

template<> struct Psi<3> {
  double xc{0.0}, yc{0.0}, zc{0.0};
  double sigma{0.18016836131796748}, r{0.15}, delta{0.0};
  double operator()(double x, double y, double z) const {
    const double dx = x - xc, dy = y - yc, dz = z - zc;
    const double rr = dx * dx + dy * dy + dz * dz;
    return std::exp((r * r - rr) / (sigma * sigma)) - 1.0 + delta;
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
rotVel2D(double x, double y, double time) noexcept {
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
rotVel3D(double x, double y, double z, double time) noexcept {
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
template<std::size_t DIM> struct Geom;

// 2D
template<> struct Geom<2> {
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
template<> struct Geom<3> {
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
  static std::array<std::array<double, 9>, 2> geom(
    const std::array<double, 2>& minCorner = {-0.5, -0.5},
    const std::array<double, 2>& maxCorner = {+0.5, +0.5}) {
    return Geom<2>::coords(minCorner, maxCorner);
  }
  static double period() {
    return 8.0;
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
  static std::array<std::array<double, 27>, 3> geom(
    const std::array<double, 3>& minCorner = {-0.5, -0.5, -0.5},
    const std::array<double, 3>& maxCorner = {+0.5, +0.5, +0.5}) {
    return Geom<3>::coords(minCorner, maxCorner);
  }
  static double period() {
    return 4.0;  // keep your original choice
  }
};

// ---------------- Scenarios ----------------
enum class Scenario { RB2D, RB3D, ROT2D, ROT3D };
static inline Scenario parse_scenario(const std::string& s) {
  if (s == "RB2D"  || s == "rb2d")  return Scenario::RB2D;
  if (s == "RB3D"  || s == "rb3d")  return Scenario::RB3D;
  if (s == "ROT2D" || s == "rot2d") return Scenario::ROT2D;
  if (s == "ROT3D" || s == "rot3d") return Scenario::ROT3D;
  return Scenario::RB2D;
}

// ---------------- The dimensioned runner ----------------
template<std::size_t DIM>
int run(int /*argc*/, char** /*argv*/, unsigned nSteps, Scenario scenario, bool vtu, bool pprof, const u32 max_depth) {
  if (pprof) ProfilerStart((DIM == 3) ? "profiling_3d.prof" : "profiling_2d.prof");

  const u32 maxDepth = max_depth;
  const u32 minDepth = (DIM == 3 ? 2u : 3u);
  const bool allowDrop = true;

  std::cout << (DIM == 3 ? "Octree 3D" : "Octree 2D")
            << " config: maxDepth=" << maxDepth
            << " minDepth=" << minDepth
            << " allowDrop=" << (allowDrop ? 1 : 0) << "\n";

  OctTree<DIM> ot(maxDepth, minDepth);
  ot.set_allow_coarsen_below_min(allowDrop);

  // ---- Scenario selection (geometry + velocity) ----
  std::array<double, DIM> minCorner{};
  std::array<double, DIM> maxCorner{};

  // velocity function pointers
  if constexpr(DIM == 2) {
    using Vel2 = std::array<double, 2> (*)(double, double, double);
    Vel2 evalVelocity = nullptr;

    // reconcile scenario with dimension
    if (scenario == Scenario::RB3D || scenario == Scenario::ROT3D) scenario = Scenario::ROT2D;

    switch (scenario) {
    case Scenario::RB2D:
      minCorner = { -0.5, 0.0 };
      maxCorner = { +0.5, 1.0 };
      evalVelocity = &RisingBubble2D;
      break;
    case Scenario::ROT2D:
      minCorner = { -0.5, -0.5 };
      maxCorner = { +0.5, +0.5 };
      evalVelocity = &rotVel2D;
      break;
    default: // fallback
      minCorner = { -0.5, -0.5 };
      maxCorner = { +0.5, +0.5 };
      evalVelocity = &rotVel2D;
      break;
    }

    auto X = DimOps<DIM>::geom(minCorner, maxCorner);
    ot.set_physical_coordinates(X);

    Psi<DIM> psi0, psi1 = psi0;
    DimOps<DIM>::shift(psi1);

    const double tau_ref    = 2.0;
    const double tau_coarse = 1e-5;
    const u32    min_level  = 0;

    const auto refine1  = [&](u32, const std::vector<Pt<2>>&, const std::vector<Pt<2>>& pts_xyz,
    const std::vector<std::array<double, 9>>&) -> bool {
      if (pts_xyz.empty()) return false;
      double v0 = psi1(pts_xyz[0][0], pts_xyz[0][1]), mn = v0, mx = v0;
      for (size_t i = 1; i < pts_xyz.size(); ++i) {
        double v = psi1(pts_xyz[i][0], pts_xyz[i][1]);
        mn = std::min(mn, v);
        mx = std::max(mx, v);
      }
      return (mn <= 0.0 && mx >= 0.0) || (mx - mn > tau_ref);
    };
    const auto coarsen1 = [&](u32, u32 level,
                              const std::vector<Pt<2>>&,
                              const std::vector<Pt<2>>& pts_xyz,
    const std::vector<std::array<double, 9>>&) -> bool {
      if (level <= min_level || pts_xyz.empty()) return false;
      double v0 = psi1(pts_xyz[0][0], pts_xyz[0][1]), mn = v0, mx = v0;
      for (size_t i = 1; i < pts_xyz.size(); ++i) {
        double v = psi1(pts_xyz[i][0], pts_xyz[i][1]);
        mn = std::min(mn, v);
        mx = std::max(mx, v);
      }
      return (mn > +tau_coarse) || (mx < -tau_coarse);
    };

    const typename DimOps<DIM>::BasisType bHi = DimOps<DIM>::highest_basis();
    std::size_t changed = ot.adapt_cycle(coarsen1, refine1, /*max_passes=*/10);
    std::cout << "Adapt-cycle changed " << changed << " cells\n";
    std::cout << "Leaves after cycle: " << ot.leaf_count() << "\n";

    const u32 fid = ot.add_field(bHi, "phi");
    ot.resize_fields_to_nodes();

    Field& fld = ot.field(fid);
    const auto bs = to_basis<DIM>(fld.basis_id);
    const auto& nodes = ot.basis_nodes(bs);
    fld.nodal.resize(nodes.size());
    for (size_t gid = 0; gid < nodes.size(); ++gid) {
      const auto& p = nodes[gid].physical;
      fld.nodal[gid] = DimOps<DIM>::evalPsi(psi1, p);
    }

    // --- VTU frame ---
    auto write_vtu_frame = [&](int k, double time) {
      std::array<const char*, 3> names = {"u", "v", "w"};
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
        const auto v = evalVelocity(p[0], p[1], time);
        for (size_t d = 0; d < DIM; ++d) f[d]->nodal[gid] = v[d];
      }

      std::vector<std::variant<u32, std::vector<u32>>> groups;
      groups.emplace_back(0u); // phi
      groups.emplace_back(std::vector<u32>(fidv.begin(), fidv.begin() + DIM)); // vector

      const std::string filename =
        std::string("./output/element_adaptive2d.") + std::to_string(k) + ".vtu";
      ot.write_binary_vtu(filename, bHi, groups, false);
    };

    if (vtu) write_vtu_frame(0, 0.);

    // time-stepping
    const unsigned nIter = 320;
    const double period = (scenario == Scenario::ROT2D ? 8.0 : 8.0);
    const double dt = period / static_cast<double>(nIter);

    OctTree<DIM> ot0 = ot;
    OctTree<DIM> ot1(maxDepth, minDepth);

    for (u32 k = 1; k <= nSteps; ++k) {
      const double time = k * dt;

      std::vector<Pt<DIM>> leftOld, stayedNew;
      fem::advect_markers_forward_analytic(ot, time, dt, evalVelocity, leftOld, stayedNew);

      ot1.reset(false, false);
      ot1.set_allow_coarsen_below_min(allowDrop);
      ot1.set_physical_coordinates(X);
      ot1.refine_to_contain_points(stayedNew, ot1.max_depth());

      const u32 fid0 = fid;
      const u32 fid1 = ot1.add_field(bHi, "phi");

      fem::advect_nodes_backward_and_transport_field_analytic(ot, fid0, time, dt, evalVelocity, ot1, fid1);

      const u32 num_coarsened = ot1.coarsen_only_cycle_safe(fid0, tau_coarse, ot);
      std::cout << "Coarsened " << num_coarsened << " leaves.\n";

      using std::swap; swap(ot, ot1);
      if (vtu) write_vtu_frame(k, time);
    }

    OctTree<DIM> ot3(ot0, 0, ot, 0);
    if (vtu) {
      std::string filename = std::string("./output/element_adaptive_union2d.vtu");
      ot3.write_binary_vtu(filename, bHi, false);
      std::cout << "Printing " << filename << "\n";
    }
  }
  else {   // DIM == 3
    using Vel3 = std::array<double, 3> (*)(double, double, double, double);
    Vel3 evalVelocity = nullptr;

    if (scenario == Scenario::RB2D || scenario == Scenario::ROT2D) scenario = Scenario::ROT3D;

    switch (scenario) {
    case Scenario::RB3D:
      minCorner = { -0.5, -0.5, 0.0 };
      maxCorner = { +0.5, +0.5, 1.0 };
      evalVelocity = &RisingBubble3D;
      break;
    case Scenario::ROT3D:
      minCorner = { -0.5, -0.5, -0.5 };
      maxCorner = { +0.5, +0.5, +0.5 };
      evalVelocity = &rotVel3D;
      break;
    default:
      minCorner = { -0.5, -0.5, -0.5 };
      maxCorner = { +0.5, +0.5, +0.5 };
      evalVelocity = &RisingBubble3D;
      break;
    }

    auto X = DimOps<DIM>::geom(minCorner, maxCorner);
    ot.set_physical_coordinates(X);

    Psi<DIM> psi0, psi1 = psi0;
    DimOps<DIM>::shift(psi1);

    const double tau_ref    = 2.0;
    const double tau_coarse = 1e-5;
    const u32    min_level  = 0;

    const auto refine1  = [&](u32, const std::vector<Pt<3>>&,
                              const std::vector<Pt<3>>& pts_xyz,
    const std::vector<std::array<double, 27>>&) -> bool {
      if (pts_xyz.empty()) return false;
      double v0 = psi1(pts_xyz[0][0], pts_xyz[0][1], pts_xyz[0][2]), mn = v0, mx = v0;
      for (size_t i = 1; i < pts_xyz.size(); ++i) {
        double v = psi1(pts_xyz[i][0], pts_xyz[i][1], pts_xyz[i][2]);
        mn = std::min(mn, v);
        mx = std::max(mx, v);
      }
      return (mn <= 0.0 && mx >= 0.0) || (mx - mn > tau_ref);
    };
    const auto coarsen1 = [&](u32, u32 level,
                              const std::vector<Pt<3>>&,
                              const std::vector<Pt<3>>& pts_xyz,
    const std::vector<std::array<double, 27>>&) -> bool {
      if (level <= min_level || pts_xyz.empty()) return false;
      double v0 = psi1(pts_xyz[0][0], pts_xyz[0][1], pts_xyz[0][2]), mn = v0, mx = v0;
      for (size_t i = 1; i < pts_xyz.size(); ++i) {
        double v = psi1(pts_xyz[i][0], pts_xyz[i][1], pts_xyz[i][2]);
        mn = std::min(mn, v);
        mx = std::max(mx, v);
      }
      return (mn > +tau_coarse) || (mx < -tau_coarse);
    };

    const typename DimOps<DIM>::BasisType bHi = DimOps<DIM>::highest_basis();
    std::size_t changed = ot.adapt_cycle(coarsen1, refine1, /*max_passes=*/10);
    std::cout << "Adapt-cycle changed " << changed << " cells\n";
    std::cout << "Leaves after cycle: " << ot.leaf_count() << "\n";

    const u32 fid = ot.add_field(bHi, "phi");
    ot.resize_fields_to_nodes();

    Field& fld = ot.field(fid);
    const auto bs = to_basis<DIM>(fld.basis_id);
    const auto& nodes = ot.basis_nodes(bs);
    fld.nodal.resize(nodes.size());
    for (size_t gid = 0; gid < nodes.size(); ++gid) {
      const auto& p = nodes[gid].physical;
      fld.nodal[gid] = DimOps<DIM>::evalPsi(psi1, p);
    }

    auto write_vtu_frame = [&](int k, double time) {
      std::array<const char*, 3> names = {"u", "v", "w"};
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
        const auto v = evalVelocity(p[0], p[1], p[2], time);
        for (size_t d = 0; d < DIM; ++d) f[d]->nodal[gid] = v[d];
      }

      std::vector<std::variant<u32, std::vector<u32>>> groups;
      groups.emplace_back(0u); // phi
      groups.emplace_back(std::vector<u32>(fidv.begin(), fidv.begin() + DIM)); // vector

      const std::string filename =
        std::string("./output/element_adaptive3d.") + std::to_string(k) + ".vtu";
      ot.write_binary_vtu(filename, bHi, groups, false);
    };

    if (vtu) write_vtu_frame(0, 0.);

    const unsigned nIter = 320;
    const double period = (scenario == Scenario::ROT3D ? 4.0 : 4.0);
    const double dt = period / static_cast<double>(nIter);

    OctTree<DIM> ot0 = ot;
    OctTree<DIM> ot1(maxDepth, minDepth);

    for (u32 k = 1; k <= nSteps; ++k) {
      const double time = k * dt;

      std::vector<Pt<DIM>> leftOld, stayedNew;
      fem::advect_markers_forward_analytic(ot, time, dt, evalVelocity, leftOld, stayedNew);

      ot1.reset(false, false);
      ot1.set_allow_coarsen_below_min(allowDrop);
      ot1.set_physical_coordinates(X);
      ot1.refine_to_contain_points(stayedNew, ot1.max_depth());

      const u32 fid0 = fid;
      const u32 fid1 = ot1.add_field(bHi, "phi");

      fem::advect_nodes_backward_and_transport_field_analytic(ot, fid0, time, dt, evalVelocity, ot1, fid1);

      const u32 num_coarsened = ot1.coarsen_only_cycle_safe(fid0, tau_coarse, ot);
      std::cout << "Coarsened " << num_coarsened << " leaves.\n";

      using std::swap; swap(ot, ot1);
      if (vtu) write_vtu_frame(k, time);
    }

    OctTree<DIM> ot3(ot0, 0, ot, 0);
    if (vtu) {
      std::string filename = std::string("./output/element_adaptive_union3d.vtu");
      ot3.write_binary_vtu(filename, bHi, false);
      std::cout << "Printing " << filename << "\n";
    }
  }

  if (pprof) ProfilerStop();
  return 0;
}

// ---------------- CLI dispatcher ----------------
static void print_usage(const char* prog) {
  std::cerr << "Usage: " << prog
            << " [-d|--dim 2|3] [-n|--niter N] [--scenario RB2D|RB3D|ROT2D|ROT3D]\n"
            << "  --niter controls the number of STEPS executed (nSteps),\n"
            << "  while dt is computed with a fixed nIter=320.\n";
}

int main(int argc, char** argv) {
  int dim = 2;            // default dimension
  unsigned nSteps = 320;  // default steps
  bool pprof = false;
  bool vtu = true;
  Scenario scenario = Scenario::ROT2D; // default; reconciled with dim inside run

  unsigned max_depth = 8;

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
      if (v <= 0) {
        std::cerr << "Error: --niter must be positive\n";
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
    else if (a == "--no-profiling") {
      pprof = false;
    }
    else if (a == "--vtu-output") {
      vtu = true;
    }
    else if (a == "--no-vtu-output") {
      vtu = false;
    }
    if (a == "--max_depth") {
      if (i + 1 >= argc) {
        print_usage(argv[0]);
        return 1;
      }
      max_depth = std::atoi(argv[++i]);
    }
  }

  switch (dim) {
  case 2: return run<2>(argc, argv, nSteps, scenario, vtu, pprof, max_depth);
  case 3: return run<3>(argc, argv, nSteps, scenario, vtu, pprof, max_depth);
  default:
    std::cerr << "Error: DIM must be 2 or 3 (got " << dim << ")\n";
    print_usage(argv[0]);
    return 2;
  }
}
