#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <gperftools/profiler.h>

#include "QuadTree2D.hpp"
#include "FieldAdvection.hpp"
#include "FieldAdvectionWithAnalyticVelocity.hpp"

using namespace  fem;

// ---------------- Level set (example) ----------------
struct Psi {
  double xc{0.0}, yc{0.0}, sigma{0.18016836131796748}, r{0.15}, delta{0.0};
  double operator()(double x, double y) const {
    return (std::exp((r * r - (x - xc) * (x - xc) - (y - yc) * (y - yc)) / (sigma * sigma)) + 0 * std::exp((r * r - (x - xc) * (x - xc) - (y + yc) * (y + yc)) / (sigma * sigma)) - 1.0 + delta);
  }
};


template <class PsiFunc>
auto make_refine_pred(const PsiFunc &psi, double tau_refine, fem::u32 max_level_cap) {
  return [&, tau_refine, max_level_cap](fem::u32 /*leaf*/,
                                        const std::vector<std::array<double, 2>> & /*pts_xi*/,
                                        const std::vector<std::array<double, 2>> &pts_xy,
  const std::vector<std::array<double, 9>> & /*Nvals*/) -> bool {
    // guard
    if (pts_xy.empty()) return false;

    // Variation of psi over the probe points
    double v0 = psi(pts_xy[0][0], pts_xy[0][1]);
    double mn = v0, mx = v0;
    for (size_t i = 1; i < pts_xy.size(); ++i) {
      double v = psi(pts_xy[i][0], pts_xy[i][1]);
      mn = std::min(mn, v);
      mx = std::max(mx, v);
    }
    const double var = mx - mn;

    // If sign changes OR variation exceeds tolerance, request refine.
    // The max_level_cap is enforced by QuadTree2D::refine() anyway,
    // but we keep it here to avoid asking for refinement needlessly.
    (void)max_level_cap; // cap is enforced in refine()
    const bool sign_change = (mn <= 0.0 && mx >= 0.0);
    return sign_change || (var > tau_refine);
  };
}

// ---------------- Coarsen predicate (robust) ----------------
template<class PsiFunc>
auto make_coarsen_pred(const PsiFunc& psi, double tau_coarse, u32 min_level) {
  return [&, tau_coarse, min_level](u32 /*parent*/, u32 level,
                                    const std::vector<std::array<double, 2>>& /*pts_xi*/,
                                    const std::vector<std::array<double, 2>>& pts_xy,
  const std::vector<std::array<double, 9>>& /*Nvals*/) -> bool {
    if (level <= min_level) return false;
    if (pts_xy.empty()) return false;
    double v0 = psi(pts_xy[0][0], pts_xy[0][1]);
    double mn = v0, mx = v0;
    for (size_t i = 1; i < pts_xy.size(); ++i) {
      double v = psi(pts_xy[i][0], pts_xy[i][1]);
      mn = std::min(mn, v);
      mx = std::max(mx, v);
    }
    // Far from interface: strictly on one side with margin tau_coarse
    return (mn > +tau_coarse) || (mx < -tau_coarse);
  };
}


auto rotVel = [](double x, double y, double time) -> std::array<double, 2> {

  double T = 8.;
  x += 0.5;
  y += 0.5;

  double u = -2. * sin(M_PI * x) * sin(M_PI * x) * sin(M_PI * y) * cos(M_PI * y) * cos(M_PI * time / T);
  double v =  2. * sin(M_PI * x) * cos(M_PI * x) * sin(M_PI * y) * sin(M_PI * y) * cos(M_PI * time / T);



  // y += 0.25;
  //
  //
  // double u = sin(M_PI * 2 * x) * sin(M_PI * 2 * y) * cos(M_PI * time / T);
  // double v = cos(M_PI * 2 * x) * cos(M_PI * 2 * y) * cos(M_PI * time / T);



  return { u, v };
};

int main() {

  ProfilerStart("profiling.prof");
  // -------- Explicit control --------
  const u32 maxDepth   = 12;   // absolute cap on tree depth
  const u32 minDepth   = 3;    // baseline depth enforced by the class
  const bool allowDrop = true; // allow coarsening below minDepth where safe?

  std::cout << "Quadtree config: maxDepth=" << maxDepth
            << " minDepth=" << minDepth
            << " allowDrop=" << (allowDrop ? 1 : 0) << "\n";

  // Construct the quadtree with explicit minDepth
  QuadTree2D qt(maxDepth, minDepth);

  double X[9], Y[9];
  auto phys = [](double xi, double eta) {
    // identity map (distort later if you want)
    return std::array<double, 2> {xi, eta};
  };
  const double nodePE[9][2] = {
    {-0.5, -0.5}, {+0.5, -0.5}, {+0.5, +0.5}, {-0.5, +0.5}, {0, -0.5}, {+0.5, 0}, {0, +0.5}, {-0.5, 0}, {0, 0}
  };
  for (int a = 0; a < 9; ++a) {
    auto p = phys(nodePE[a][0], nodePE[a][1]);
    X[a] = p[0];
    Y[a] = p[1];
  }
  qt.set_physical_quad9(X, Y);

  // Soft floor policy (true = may coarsen below floor where coarsen_pred says it's safe)
  qt.set_allow_coarsen_below_min(allowDrop);

  // Level set and thresholds
  Psi psi0;

  const double tau_ref    = 2;  // refine if |psi| small or sign-crossing
  const u32 fid = qt.add_field(Basis::Q2_Quad9);

  // --- Now change Psi slightly and run coarsen+refine with hysteresis ---

  Psi psi1 = psi0;
  psi1.xc = 0.;      // small shift of interface
  psi1.yc = 0.25;      // small shift of interface

  const double tau_coarse = 1e-5;  // coarsen margin ( > tau_ref for hysteresis)
  const u32    min_level  = 0;     // optionally keep a minimum level

  auto refine1  = make_refine_pred(psi1, tau_ref, maxDepth);
  auto coarsen1 = make_coarsen_pred(psi1, tau_coarse, min_level);

  std::size_t changed = qt.adapt_cycle(coarsen1, refine1, Basis::Q2_Quad9, /*max_passes=*/10);
  std::cout << "Adapt-cycle changed " << changed << " cells\n";
  std::cout << "Leaves after cycle: " << qt.leaf_count() << "\n";

// --- Resize field storage to match current global node set ---
  qt.resize_fields_to_nodes();

  Field& fld = qt.field(fid);
  const Basis b = fld.basis;

// Get the unique global FEM nodes (with physical coords)
  const auto& nodes = qt.basis_nodes(b);

// Initialize nodal values from psi1 at node physical coordinates
  fld.nodal.resize(nodes.size());
  for (size_t gid = 0; gid < nodes.size(); ++gid) {
    const auto& xy = nodes[gid].physical;  // (x,y)
    fld.nodal[gid] = psi1(xy[0], xy[1]);
  }

  std::string filename = "./output/element_adaptive." + std::to_string(0) + ".vtu";
  qt.write_binary_vtu(filename, fid, "u", /*cell_centered=*/false);
  std::cout << "Printing " << filename << "\n";





  double period = 8;
  unsigned nIterations = 320;
  double dt = period / nIterations;

  QuadTree2D qt1(qt.max_depth());


  for (u32 k = 1; k <= 0 + 1 * nIterations; k++) {

    std::vector<std::array<double, 2>> vOld = {{+1, -1}, { +1, +1}, {-1, +1}, {-1, -1}, {+1, 0}, {0, +1}, {-1, 0,}, {0, -1}, {0, 0}};
    std::vector<std::array<double, 2>> vNew = {{+1, -1}, { +1, +1}, {-1, +1}, {-1, -1}, {+1, 0}, {0, +1}, {-1, 0,}, {0, -1}, {0, 0}};
    // fill vOld, vNew with {u,v} at coarse nodes

    qt1.reset(false, false);
    qt1.set_allow_coarsen_below_min(allowDrop);
    qt1.set_physical_quad9(X, Y);

    double time = k * dt;

    std::vector<std::array<double, 2>> leftOld, stayedNew;
    //fem::advect_markers_forward(qt, Basis::Q2_Quad9, vOld, vNew, dt, leftOld, stayedNew);
    fem::advect_markers_forward_analytic(qt, time, dt, rotVel, leftOld, stayedNew);


    qt1.refine_to_contain_points(stayedNew, qt1.max_depth());

    u32 fid0 = 0; // some scalar field id in qt0
    u32 fid1 = qt1.add_field(Basis::Q2_Quad9);
    //fem::advect_nodes_backward_and_transport_field(qt, fid0, Basis::Q2_Quad9, vOld, vNew, dt, qt1, fid1);
    fem::advect_nodes_backward_and_transport_field_analytic(qt, fid0, time, dt, rotVel, qt1, fid1);


    u32 num_coarsened = qt1.coarsen_only_cycle_safe(fid0, tau_coarse, qt);
    std::cout << "Coarsened " << num_coarsened << " leaves.\n";

    using std::swap;
    swap(qt, qt1);

    filename = "./output/element_adaptive." + std::to_string(k) + ".vtu";
    qt.write_binary_vtu(filename, fid, "u", false);
    std::cout << "Printing " << filename << "\n";
  }

  // filename = "./output/element_adaptive." + std::to_string(100) + ".vtu";
  // qt.write_binary_vtu(filename, fid, "u", false);
  // std::cout << "Printing " << filename << "\n";

  ProfilerStop();
  return 0;
}




