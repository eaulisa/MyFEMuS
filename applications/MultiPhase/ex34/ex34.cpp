#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <gperftools/profiler.h>

#include "HexTree3D.hpp"
//#include "FieldAdvection3D.hpp"
//#include "FieldAdvectionWithAnalyticVelocity3D.hpp"

using namespace fem;

// ---------------- Level set (example, 3D) ----------------
struct Psi3D {
  double xc{0.0}, yc{0.0}, zc{0.0};
  double sigma{0.18016836131796748}, r{0.15}, delta{0.0};

  double operator()(double x, double y, double z) const {
    const double dx = x - xc, dy = y - yc, dz = z - zc;
    const double rr = dx*dx + dy*dy + dz*dz;
    return std::exp((r*r - rr) / (sigma*sigma)) - 1.0 + delta;
  }
};

// ---------------- Refine predicate (3D) ----------------
template <class PsiFunc>
auto make_refine_pred_3d(const PsiFunc& psi, double tau_refine, fem::u32 max_level_cap) {
  return [&, tau_refine, max_level_cap](fem::u32 /*leaf*/,
                                        const std::vector<std::array<double,3>>& /*pts_s*/,
                                        const std::vector<std::array<double,3>>& pts_xyz,
                                        const std::vector<std::array<double,27>>& /*Nvals*/) -> bool {
    if (pts_xyz.empty()) return false;

    double v0 = psi(pts_xyz[0][0], pts_xyz[0][1], pts_xyz[0][2]);
    double mn = v0, mx = v0;
    for (size_t i = 1; i < pts_xyz.size(); ++i) {
      const auto& p = pts_xyz[i];
      double v = psi(p[0], p[1], p[2]);
      mn = std::min(mn, v);
      mx = std::max(mx, v);
    }
    const bool sign_change = (mn <= 0.0 && mx >= 0.0);
    const double var = mx - mn;
    (void)max_level_cap; // enforced by the tree anyway
    return sign_change || (var > tau_refine);
  };
}

// ---------------- Coarsen predicate (robust, 3D) ----------------
template <class PsiFunc>
auto make_coarsen_pred_3d(const PsiFunc& psi, double tau_coarse, u32 min_level) {
  return [&, tau_coarse, min_level](u32 /*parent*/, u32 level,
                                    const std::vector<std::array<double,3>>& /*pts_s*/,
                                    const std::vector<std::array<double,3>>& pts_xyz,
                                    const std::vector<std::array<double,27>>& /*Nvals*/) -> bool {
    if (level <= min_level) return false;
    if (pts_xyz.empty()) return false;

    double v0 = psi(pts_xyz[0][0], pts_xyz[0][1], pts_xyz[0][2]);
    double mn = v0, mx = v0;
    for (size_t i = 1; i < pts_xyz.size(); ++i) {
      const auto& p = pts_xyz[i];
      double v = psi(p[0], p[1], p[2]);
      mn = std::min(mn, v);
      mx = std::max(mx, v);
    }
    return (mn > +tau_coarse) || (mx < -tau_coarse);
  };
}

// ---------------- Analytic velocity (3D: TG-like, w=0 by default) ----------------
auto rotVel3D = [](double x, double y, double z, double time) -> std::array<double,3> {
  (void)z; // easy extension: add a z-structure if you want
  double T = 8.0;
  x += 0.5; y += 0.5;
  double u = -2. * std::sin(M_PI*x)*std::sin(M_PI*x) * std::sin(M_PI*y)*std::cos(M_PI*y) * std::cos(M_PI*time/T);
  double v =  2. * std::sin(M_PI*x)*std::cos(M_PI*x) * std::sin(M_PI*y)*std::sin(M_PI*y) * std::cos(M_PI*time/T);
  double w = 0.0;
  return {u, v, w};
};

int main() {
  //ProfilerStart("profiling_3d.prof");

  // -------- Explicit control --------
  const u32 maxDepth   = 6;   // absolute cap on tree depth (<=20 is supported by your impl)
  const u32 minDepth   = 2;    // baseline depth
  const bool allowDrop = true; // allow coarsening below minDepth where safe?

  std::cout << "Octree config: maxDepth=" << maxDepth
            << " minDepth=" << minDepth
            << " allowDrop=" << (allowDrop ? 1 : 0) << "\n";

  // Construct the octree
  OctTree3D ot(maxDepth, minDepth);

  // --- Coarse HEX27 geometry: identity map on [-0.5,0.5]^3 ---
  double X[27], Y[27], Z[27];

  // Parent nodes in your H27 order:
  // corners 0..7
  // bottom edges 8..11, top edges 12..15, vertical edges 16..19
  // face centers 20..25, cell center 26
  auto add = [&](int i, double xi, double eta, double z) {
    X[i]=xi; Y[i]=eta; Z[i]=z;
  };

  const double a=-0.5,b=+0.5, xm=0.0, ym=0.0, zm=0.0;

  // corners
  add(0,a,a,a); add(1,b,a,a); add(2,b,b,a); add(3,a,b,a);
  add(4,a,a,b); add(5,b,a,b); add(6,b,b,b); add(7,a,b,b);
  // z=-1 face edges
  add( 8,xm,a,a); add( 9,b,ym,a); add(10,xm,b,a); add(11,a,ym,a);
  // z=+1 face edges
  add(12,xm,a,b); add(13,b,ym,b); add(14,xm,b,b); add(15,a,ym,b);
  // vertical edges
  add(16,a,a,zm); add(17,b,a,zm); add(18,b,b,zm); add(19,a,b,zm);
  // face centers
  add(20,xm,a,zm); add(21,b,ym,zm); add(22,xm,b,zm); add(23,a,ym,zm);
  add(24,xm,ym,a); add(25,xm,ym,b);
  // cell center
  add(26,xm,ym,zm);

  ot.set_physical_hex27(X,Y,Z);
  ot.set_allow_coarsen_below_min(allowDrop);

  // Level set and thresholds
  Psi3D psi0;

  const double tau_ref = 2.0;  // refine tolerance
  const u32 fid = ot.add_field(Basis::H27);

  // Slightly shifted level set
  Psi3D psi1 = psi0;
  psi1.xc = 0.0;
  psi1.yc = 0.25;
  psi1.zc = 0.0;

  const double tau_coarse = 1e-5;
  const u32    min_level  = 0;

  auto refine1  = make_refine_pred_3d(psi1, tau_ref, maxDepth);
  auto coarsen1 = make_coarsen_pred_3d(psi1, tau_coarse, min_level);

  std::size_t changed = ot.adapt_cycle(coarsen1, refine1, Basis::H27, /*max_passes=*/10);
  std::cout << "Adapt-cycle changed " << changed << " cells\n";
  std::cout << "Leaves after cycle: " << ot.leaf_count() << "\n";

  // --- Resize field storage to match current global node set ---
  ot.resize_fields_to_nodes();

  Field& fld = ot.field(fid);
  const Basis bs = fld.basis;

  // Unique global FEM nodes with physical coords
  const auto& nodes = ot.basis_nodes(bs);

  // Initialize nodal values from psi1 at node physical coordinates
  fld.nodal.resize(nodes.size());
  for (size_t gid = 0; gid < nodes.size(); ++gid) {
    const auto& xyz = nodes[gid].physical;
    fld.nodal[gid] = psi1(xyz[0], xyz[1], xyz[2]);
  }

  std::string filename = "./output/element_adaptive3d." + std::to_string(0) + ".vtu";
  ot.write_binary_vtu(filename, fid, "u", /*cell_centered=*/false);
  std::cout << "Printing " << filename << "\n";

  // --- Time loop (analytic advection) ---
  double period = 8.0;
  unsigned nIterations = 320;
  double dt = period / nIterations;

  OctTree3D ot1(ot.max_depth());

  // for (u32 k = 1; k <= nIterations; ++k) {
  //   // coarse-node reference marker set for Hex27 in parent space
  //   std::vector<std::array<double,3>> vOld = {
  //     {+1,-1,-1},{+1,+1,-1},{-1,+1,-1},{-1,-1,-1},
  //     {+1,-1,+1},{+1,+1,+1},{-1,+1,+1},{-1,-1,+1},
  //     {+1, 0,-1},{ 0,+1,-1},{-1, 0,-1},{ 0,-1,-1},
  //     {+1, 0,+1},{ 0,+1,+1},{-1, 0,+1},{ 0,-1,+1},
  //     {-1,-1, 0},{+1,-1, 0},{+1,+1, 0},{-1,+1, 0},
  //     { 0,-1, 0},{+1, 0, 0},{ 0,+1, 0},{-1, 0, 0},
  //     { 0, 0,-1},{ 0, 0,+1},{ 0, 0, 0}
  //   };
  //   std::vector<std::array<double,3>> vNew = vOld; // placeholder
  //
  //   ot1.reset(false, false);
  //   ot1.set_allow_coarsen_below_min(allowDrop);
  //   ot1.set_physical_hex27(X, Y, Z);
  //
  //   double time = k * dt;
  //
  //   std::vector<std::array<double,3>> leftOld, stayedNew;
  //   // Forward trace analytic markers
  //   fem::advect_markers_forward_analytic_3d(ot, time, dt, rotVel3D, leftOld, stayedNew);
  //
  //   // Refine ot1 to contain all traced markers
  //   ot1.refine_to_contain_points(stayedNew, ot1.max_depth());
  //
  //   u32 fid0 = fid; // advected scalar from ot
  //   u32 fid1 = ot1.add_field(Basis::H27);
  //
  //   // Backtrace element nodes and transport field (analytic 3D velocity)
  //   fem::advect_nodes_backward_and_transport_field_analytic_3d(ot, fid0, time, dt, rotVel3D, ot1, fid1);
  //
  //   // Conservative coarsen (safe), using snapshot
  //   u32 num_coarsened = ot1.coarsen_only_cycle_safe(fid0, tau_coarse, ot);
  //   std::cout << "Coarsened " << num_coarsened << " leaves.\n";
  //
  //   using std::swap;
  //   swap(ot, ot1);
  //
  //   filename = "./output/element_adaptive3d." + std::to_string(k) + ".vtu";
  //   ot.write_binary_vtu(filename, fid, "u", false);
  //   std::cout << "Printing " << filename << "\n";
  // }

  //ProfilerStop();
  return 0;
}



























// #include <algorithm>
// #include <array>
// #include <cmath>
// #include <iostream>
// #include <limits>
// #include <string>
// #include <vector>
// #include <gperftools/profiler.h>
//
// //#include "FieldAdvection.hpp"
// //#include "FieldAdvectionWithAnalyticVelocity.hpp"
//
// #include "HexTree3D.hpp"
// #include <iostream>
//
// using namespace  fem;
//
// int main(){
//   fem::OctTree3D T(6);
//   // unit cube mapping: parent == physical (for quick test)
//   double X[27],Y[27],Z[27];
//   // fill H27 grid points in the given order
//   const fem::Point3 P[27] = {
//     {-1,-1,-1},{ 1,-1,-1},{ 1, 1,-1},{-1, 1,-1},{-1,-1, 1},{ 1,-1, 1},{ 1, 1, 1},{-1, 1, 1},
//     { 0,-1,-1},{ 1, 0,-1},{ 0, 1,-1},{-1, 0,-1},{ 0,-1, 1},{ 1, 0, 1},{ 0, 1, 1},{-1, 0, 1},
//     {-1,-1, 0},{ 1,-1, 0},{ 1, 1, 0},{-1, 1, 0},{ 0,-1, 0},{ 1, 0, 0},{ 0, 1, 0},{-1, 0, 0},
//     { 0, 0,-1},{ 0, 0, 1},{ 0, 0, 0}
//   };
//   for(int i=0;i<27;++i){ X[i]=P[i][0]; Y[i]=P[i][1]; Z[i]=P[i][2]; }
//   T.set_physical_hex27(X,Y,Z);
//
//   // refine root once
//   T.refine_leaf_once(T.leaves()[0]);
//
//   // pick a point in parent space and find its leaf
//   fem::Point3 s{0.2, -0.1, 0.3};
//   auto leaf = T.locate_leaf_on_parent(s);
//   std::cout << "leaf idx: " << leaf << "\n";
//
//   // neighbor across +x
//   auto nb = T.neighbor_leaf_by_face_any(leaf, 1);
//   std::cout << "neighbor +x: " << nb << "\n";
//
//   // parent->physical and inverse map
//   auto x = T.parent_to_physical(s);
//   fem::Point3 back;
//   bool ok = T.inverse_map_hex27(x, back);
//   std::cout << "inverse ok? " << ok << "  back=("<<back[0]<<","<<back[1]<<","<<back[2]<<")\n";
// }


