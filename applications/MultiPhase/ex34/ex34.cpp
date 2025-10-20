// main.cpp
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#ifdef USE_GPERFTOOLS
#include <gperftools/profiler.h>
#endif

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
static inline std::array<double, 2>
rotVel2D(double x, double y, double time) noexcept {
  const double T = 8.0; x += 0.5; y += 0.5;
  double u = -2.*std::sin(M_PI * x) * std::sin(M_PI * x) * std::sin(M_PI * y) * std::cos(M_PI * y) * std::cos(M_PI * time / T);
  double v =  2.*std::sin(M_PI * x) * std::cos(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * y) * std::cos(M_PI * time / T);
  return {u, v};
}

static inline std::array<double, 3>
rotVel3D(double x, double y, double z, double time) noexcept {
  x += 0.5; y += 0.5; z += 0.5;
  double sx, cx, sy, cy, sz, cz;
  ::sincos(M_PI * x, &sx, &cx);
  ::sincos(M_PI * y, &sy, &cy);
  ::sincos(M_PI * z, &sz, &cz);
  const double cosT = std::cos(M_PI * (time * 0.25));
  const double s2x = sx * sx, s2y = sy * sy, s2z = sz * sz;
  const double sxcx = sx * cx, sycy = sy * cy, szcz = sz * cz;
  const double k = 2.0 * cosT;
  return { k*s2x*(-sycy + szcz), k*s2y*(-szcz + sxcx), k*s2z*(-sxcx + sycy) };
}

// ---------------- Geometry builders ----------------
template<std::size_t DIM> struct Geom;

template<> struct Geom<2> {
  static std::array<std::array<double, 9>, 2> coords(double a = -0.5, double b = +0.5) {
    std::array<std::array<double, 9>, 2> X{};
    const double xm = 0, ym = 0;
    X[0][0] = a; X[1][0] = a; X[0][1] = b; X[1][1] = a; X[0][2] = b; X[1][2] = b; X[0][3] = a; X[1][3] = b;
    X[0][4] = xm; X[1][4] = a; X[0][5] = b;  X[1][5] = ym; X[0][6] = xm; X[1][6] = b; X[0][7] = a; X[1][7] = ym;
    X[0][8] = xm; X[1][8] = ym;
    return X;
  }
};

template<> struct Geom<3> {
  static std::array<std::array<double, 27>, 3> coords(double a = -0.5, double b = +0.5) {
    std::array<std::array<double, 27>, 3> X{};
    auto set = [&](int i, double xi, double eta, double z) {
      X[0][i] = xi;
      X[1][i] = eta;
      X[2][i] = z;
    };
    const double xm = 0, ym = 0, zm = 0;
    set(0, a, a, a); set(1, b, a, a); set(2, b, b, a); set(3, a, b, a);
    set(4, a, a, b); set(5, b, a, b); set(6, b, b, b); set(7, a, b, b);
    set(8, xm, a, a); set(9, b, ym, a); set(10, xm, b, a); set(11, a, ym, a);
    set(12, xm, a, b); set(13, b, ym, b); set(14, xm, b, b); set(15, a, ym, b);
    set(16, a, a, zm); set(17, b, a, zm); set(18, b, b, zm); set(19, a, b, zm);
    set(20, xm, a, zm); set(21, b, ym, zm); set(22, xm, b, zm); set(23, a, ym, zm);
    set(24, xm, ym, a); set(25, xm, ym, b); set(26, xm, ym, zm);
    return X;
  }
};

// ---------------- Dimension traits ----------------
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

  template<class PsiFunc>
  static auto make_refine_pred(const PsiFunc& psi, double tau_refine, u32 /*cap*/) {
    return [&, tau_refine](u32, const std::vector<Pt<2>>&,
                           const std::vector<Pt<2>>& pts_xyz,
    const std::vector<std::array<double, 9>>&) -> bool {
      if (pts_xyz.empty()) return false;
      double v0 = psi(pts_xyz[0][0], pts_xyz[0][1]), mn = v0, mx = v0;
      for (size_t i = 1; i < pts_xyz.size(); ++i) {
        double v = psi(pts_xyz[i][0], pts_xyz[i][1]);
        mn = std::min(mn, v);
        mx = std::max(mx, v);
      }
      return (mn <= 0.0 && mx >= 0.0) || (mx - mn > tau_refine);
    };
  }
  template<class PsiFunc>
  static auto make_coarsen_pred(const PsiFunc& psi, double tau_coarse, u32 min_level) {
    return [&, tau_coarse, min_level](u32, u32 level,
                                      const std::vector<Pt<2>>&,
                                      const std::vector<Pt<2>>& pts_xyz,
    const std::vector<std::array<double, 9>>&) -> bool {
      if (level <= min_level || pts_xyz.empty()) return false;
      double v0 = psi(pts_xyz[0][0], pts_xyz[0][1]), mn = v0, mx = v0;
      for (size_t i = 1; i < pts_xyz.size(); ++i) {
        double v = psi(pts_xyz[i][0], pts_xyz[i][1]);
        mn = std::min(mn, v);
        mx = std::max(mx, v);
      }
      return (mn > +tau_coarse) || (mx < -tau_coarse);
    };
  }

  template<class Tree>
  static void advect_markers(const Tree& ot, double time, double dt,
                             std::vector<Pt<2>>& leftOld, std::vector<Pt<2>>& stayedNew) {
    fem::advect_markers_forward_analytic(ot, time, dt, rotVel2D, leftOld, stayedNew);
  }
  template<class Tree>
  static void advect_nodes(const Tree& ot, u32 fid0, double time, double dt,
                           Tree& ot1, u32 fid1) {
    fem::advect_nodes_backward_and_transport_field_analytic(ot, fid0, time, dt, rotVel2D, ot1, fid1);
  }

  static std::array<std::array<double, 9>, 2> geom() {
    return Geom<2>::coords();
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
    psi.yc = 0.25;
    psi.zc = 0.0;
  }
  static double evalPsi(const Psi<3>& psi, const Pt<3>& p) {
    return psi(p[0], p[1], p[2]);
  }

  template<class PsiFunc>
  static auto make_refine_pred(const PsiFunc& psi, double tau_refine, u32 /*cap*/) {
    return [&, tau_refine](u32, const std::vector<Pt<3>>&,
                           const std::vector<Pt<3>>& pts_xyz,
    const std::vector<std::array<double, 27>>&) -> bool {
      if (pts_xyz.empty()) return false;
      double v0 = psi(pts_xyz[0][0], pts_xyz[0][1], pts_xyz[0][2]), mn = v0, mx = v0;
      for (size_t i = 1; i < pts_xyz.size(); ++i) {
        double v = psi(pts_xyz[i][0], pts_xyz[i][1], pts_xyz[i][2]);
        mn = std::min(mn, v);
        mx = std::max(mx, v);
      }
      return (mn <= 0.0 && mx >= 0.0) || (mx - mn > tau_refine);
    };
  }
  template<class PsiFunc>
  static auto make_coarsen_pred(const PsiFunc& psi, double tau_coarse, u32 min_level) {
    return [&, tau_coarse, min_level](u32, u32 level,
                                      const std::vector<Pt<3>>&,
                                      const std::vector<Pt<3>>& pts_xyz,
    const std::vector<std::array<double, 27>>&) -> bool {
      if (level <= min_level || pts_xyz.empty()) return false;
      double v0 = psi(pts_xyz[0][0], pts_xyz[0][1], pts_xyz[0][2]), mn = v0, mx = v0;
      for (size_t i = 1; i < pts_xyz.size(); ++i) {
        double v = psi(pts_xyz[i][0], pts_xyz[i][1], pts_xyz[i][2]);
        mn = std::min(mn, v);
        mx = std::max(mx, v);
      }
      return (mn > +tau_coarse) || (mx < -tau_coarse);
    };
  }

  template<class Tree>
  static void advect_markers(const Tree& ot, double time, double dt,
                             std::vector<Pt<3>>& leftOld, std::vector<Pt<3>>& stayedNew) {
    fem::advect_markers_forward_analytic(ot, time, dt, rotVel3D, leftOld, stayedNew);
  }
  template<class Tree>
  static void advect_nodes(const Tree& ot, u32 fid0, double time, double dt,
                           Tree& ot1, u32 fid1) {
    fem::advect_nodes_backward_and_transport_field_analytic(ot, fid0, time, dt, rotVel3D, ot1, fid1);
  }

  static std::array<std::array<double, 27>, 3> geom() {
    return Geom<3>::coords();
  }
  static double period() {
    return 4.0;
  }
};

// ---------------- The dimensioned runner ----------------
template<std::size_t DIM>
int run(int /*argc*/, char** /*argv*/, unsigned nSteps /* from CLI or default */) {
#ifdef USE_GPERFTOOLS
  ProfilerStart(DIM == 3 ? "profiling_3d.prof" : "profiling_2d.prof");
#endif

  const u32 maxDepth = (DIM == 3 ? 8u : 12u);
  const u32 minDepth = (DIM == 3 ? 2u : 3u);
  const bool allowDrop = true;

  std::cout << (DIM == 3 ? "Octree 3D" : "Octree 2D")
            << " config: maxDepth=" << maxDepth
            << " minDepth=" << minDepth
            << " allowDrop=" << (allowDrop ? 1 : 0) << "\n";

  OctTree<DIM> ot(maxDepth, minDepth);
  ot.set_allow_coarsen_below_min(allowDrop);

  auto X = DimOps<DIM>::geom();
  ot.set_physical_coordinates(X);

  Psi<DIM> psi0, psi1 = psi0;
  DimOps<DIM>::shift(psi1);

  const double tau_ref    = 2.0;
  const double tau_coarse = 1e-5;
  const u32    min_level  = 0;

  const auto refine1  = DimOps<DIM>::make_refine_pred(psi1, tau_ref, maxDepth);
  const auto coarsen1 = DimOps<DIM>::make_coarsen_pred(psi1, tau_coarse, min_level);

  const typename DimOps<DIM>::BasisType bHi = DimOps<DIM>::highest_basis();
  std::size_t changed = ot.adapt_cycle(coarsen1, refine1, bHi, /*max_passes=*/10);
  std::cout << "Adapt-cycle changed " << changed << " cells\n";
  std::cout << "Leaves after cycle: " << ot.leaf_count() << "\n";

  const u32 fid = ot.add_field(bHi);
  ot.resize_fields_to_nodes();

  Field& fld = ot.field(fid);
  const auto bs = to_basis<DIM>(fld.basis_id);
  const auto& nodes = ot.basis_nodes(bs);
  fld.nodal.resize(nodes.size());
  for (size_t gid = 0; gid < nodes.size(); ++gid) {
    const auto& p = nodes[gid].physical;
    fld.nodal[gid] = DimOps<DIM>::evalPsi(psi1, p);
  }

  std::string filename = std::string("./output/element_adaptive")
                         + (DIM == 3 ? "3d." : "2d.")
                         + std::to_string(0) + ".vtu";
  ot.write_binary_vtu(filename, fid, "u", false);
  std::cout << "Printing " << filename << "\n";

  // nIter is FIXED to 320 to define dt.
  const unsigned nIter = 320;
  const double period = DimOps<DIM>::period();
  const double dt = period / static_cast<double>(nIter);

  OctTree<DIM> ot0 = ot;

  OctTree<DIM> ot1(maxDepth, minDepth);



  for (u32 k = 1; k <= nSteps; ++k) {
    ot1.reset(false, false);
    ot1.set_allow_coarsen_below_min(allowDrop);
    ot1.set_physical_coordinates(X);

    const double time = k * dt;

    std::vector<Pt<DIM>> leftOld, stayedNew;
    DimOps<DIM>::advect_markers(ot, time, dt, leftOld, stayedNew);

    ot1.refine_to_contain_points(stayedNew, ot1.max_depth());

    const u32 fid0 = fid;
    const u32 fid1 = ot1.add_field(bHi);

    DimOps<DIM>::advect_nodes(ot, fid0, time, dt, ot1, fid1);

    const u32 num_coarsened = ot1.coarsen_only_cycle_safe(fid0, tau_coarse, ot);
    std::cout << "Coarsened " << num_coarsened << " leaves.\n";

    using std::swap;
    swap(ot, ot1);

    filename = std::string("./output/element_adaptive")
               + (DIM == 3 ? "3d." : "2d.")
               + std::to_string(k) + ".vtu";
    ot.write_binary_vtu(filename, fid, "u", false);
    std::cout << "Printing " << filename << "\n";
  }


  OctTree<DIM> ot3(ot0, 0, ot, 0);
  filename = std::string("./output/element_adaptive_union")
             + (DIM == 3 ? "3d." : "2d.")
             + std::to_string(0) + ".vtu";
  ot3.write_binary_vtu(filename, 0, "u", false);
  std::cout << "Printing " << filename << "\n";

  filename = std::string("./output/element_adaptive_union")
             + (DIM == 3 ? "3d." : "2d.")
             + std::to_string(1) + ".vtu";
  ot3.write_binary_vtu(filename, 1, "u", false);
  std::cout << "Printing " << filename << "\n";

#ifdef USE_GPERFTOOLS
  ProfilerStop();
#endif
  return 0;
}

// ---------------- CLI dispatcher ----------------
static void print_usage(const char* prog) {
  std::cerr << "Usage: " << prog << " [-d|--dim 2|3] [-n|--niter N]\n"
            << "  --niter controls the number of STEPS executed (nSteps),\n"
            << "  while dt is computed with a fixed nIter=320.\n";
}

int main(int argc, char** argv) {
  int dim = 3;            // default dimension
  unsigned nSteps = 320;  // default steps if --niter not provided

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
    else if (a == "-h" || a == "--help") {
      print_usage(argv[0]); return 0;
    }
  }

  switch (dim) {
  case 2: return run<2>(argc, argv, nSteps);
  case 3: return run<3>(argc, argv, nSteps);
  default:
    std::cerr << "Error: DIM must be 2 or 3 (got " << dim << ")\n";
    print_usage(argv[0]);
    return 2;
  }
}
