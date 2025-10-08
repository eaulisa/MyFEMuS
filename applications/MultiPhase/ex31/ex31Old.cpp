#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <gperftools/profiler.h>

#include "QuadTree2DOld.hpp"
#include "FieldAdvectionOld.hpp"
#include "FieldAdvectionWithAnalyticVelocityOld.hpp"

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




    // --- Re-size and refresh field coefficients for new mesh (from psi1) ---
  qt.resize_fields_to_leaves();
  const auto& L1 = qt.leaf_indices();
  for (u32 k = 0; k < (u32)L1.size(); ++k) {
    std::vector<std::array<double, 2>> xy9;
    qt.leaf_physical_nodes(Basis::Q2_Quad9, L1[k], xy9);
    double* c = qt.leaf_coeff_ptr(fid, k);
    for (int a = 0; a < 9; ++a) c[a] = psi1(xy9[a][0], xy9[a][1]);
  }

// --- Resize field storage to match current global node set ---
//   qt.resize_fields_to_nodes();
//
//   Field& fld = qt.field(fid);
//   const Basis b = fld.basis;
//
// // Get the unique global FEM nodes (with physical coords)
//   const auto& nodes = qt.basis_nodes(b);
//
// // Initialize nodal values from psi1 at node physical coordinates
//   fld.nodal.resize(nodes.size());
//   for (size_t gid = 0; gid < nodes.size(); ++gid) {
//     const auto& xy = nodes[gid].physical;  // (x,y)
//     fld.nodal[gid] = psi1(xy[0], xy[1]);
//   }

  std::string filename = "./output/element_adaptive." + std::to_string(0) + ".vtu";
  qt.write_binary_vtu(filename, fid, "u", /*cell_centered=*/false);
  std::cout << "Printing " << filename << "\n";





  double period = 8;
  unsigned nIterations = 320;
  double dt = period / nIterations;

  for (u32 k = 1; k <= nIterations; k++) {

    std::vector<std::array<double, 2>> vOld = {{+1, -1}, { +1, +1}, {-1, +1}, {-1, -1}, {+1, 0}, {0, +1}, {-1, 0,}, {0, -1}, {0, 0}};
    std::vector<std::array<double, 2>> vNew = {{+1, -1}, { +1, +1}, {-1, +1}, {-1, -1}, {+1, 0}, {0, +1}, {-1, 0,}, {0, -1}, {0, 0}};
    // fill vOld, vNew with {u,v} at coarse nodes




    double time = k * dt;

    std::vector<std::array<double, 2>> leftOld, stayedNew;
    //fem::advect_markers_forward(qt, Basis::Q2_Quad9, vOld, vNew, dt, leftOld, stayedNew);
    fem::advect_markers_forward_analytic(qt, time, dt, rotVel, leftOld, stayedNew);

    QuadTree2D qt1(qt.max_depth());
    qt1.set_allow_coarsen_below_min(allowDrop);

    qt1.set_physical_quad9(X, Y);
    qt1.refine_to_contain_points(stayedNew, qt1.max_depth());

    u32 fid0 = 0; // some scalar field id in qt0
    u32 fid1 = qt1.add_field(Basis::Q2_Quad9);
    //fem::advect_nodes_backward_and_transport_field(qt, fid0, Basis::Q2_Quad9, vOld, vNew, dt, qt1, fid1);
    fem::advect_nodes_backward_and_transport_field_analytic(qt, fid0, time, dt, rotVel, qt1, fid1);

    u32 num_coarsened = qt1.coarsen_only_cycle_safe(fid0, tau_coarse);
    std::cout << "Coarsened " << num_coarsened << " leaves.\n";

    std::swap(qt, qt1);

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












// // ---------------- Refine predicate based on a field ----------------
// inline auto make_refine_pred_field(const QuadTree2D& qt, u32 fid,
//                                    double tau_refine, u32 max_level_cap) {
//   return [&, fid, tau_refine, max_level_cap](u32 leaf,
//          const std::vector<std::array<double, 2>>& /*pts_xi*/,
//          const std::vector<std::array<double, 2>>& pts_xy,
//   const std::vector<std::array<double, 9>>& /*Nvals*/) -> bool {
//     if (pts_xy.empty()) return false;
//
//     double v0;
//     bool ok0 = qt.evaluate_physical(fid, pts_xy[0][0], pts_xy[0][1], v0);
//     if (!ok0) return false;
//     double mn = v0, mx = v0;
//
//     for (size_t i = 1; i < pts_xy.size(); ++i) {
//       double val;
//       bool ok = qt.evaluate_physical(fid, pts_xy[i][0], pts_xy[i][1], val);
//       if (!ok) continue;
//       mn = std::min(mn, val);
//       mx = std::max(mx, val);
//     }
//
//     double var = mx - mn;
//     bool sign_change = (mn <= 0.0 && mx >= 0.0);
//
//     (void)max_level_cap; // cap is enforced by QuadTree2D::refine
//     return sign_change || (var > tau_refine);
//   };
// }

// // ---------------- Coarsen predicate based on a field ----------------
// inline auto make_coarsen_pred_field(const QuadTree2D& qt, u32 fid,
//                                     double tau_coarse, u32 min_level)
// {
//   return [&, fid, tau_coarse, min_level](u32 /*parent*/, u32 level,
//                                          const std::vector<std::array<double,2>>& /*pts_xi*/,
//                                          const std::vector<std::array<double,2>>& pts_xy,
//                                          const std::vector<std::array<double,9>>& /*Nvals*/) -> bool
//   {
//     if (level <= min_level) return false;
//     if (pts_xy.empty()) return false;
//
//     double v0;
//     bool ok0 = qt.evaluate_physical(fid, pts_xy[0][0], pts_xy[0][1], v0);
//     if (!ok0) return false;
//     double mn = v0, mx = v0;
//
//     for (size_t i = 1; i < pts_xy.size(); ++i) {
//       double val;
//       bool ok = qt.evaluate_physical(fid, pts_xy[i][0], pts_xy[i][1], val);
//       if (!ok) continue;
//       mn = std::min(mn, val);
//       mx = std::max(mx, val);
//     }
//
//     // Far from interface: strictly on one side with margin tau_coarse
//     return (mn > +tau_coarse) || (mx < -tau_coarse);
//   };
// }













// int main() {
//   // -------- Explicit control --------
//   const u32 maxDepth   = 12;   // absolute cap on tree depth
//   const u32 minDepth   = 3;    // baseline depth enforced by the class
//   const bool allowDrop = true; // allow coarsening below minDepth where safe?
//
//   std::cout << "Quadtree config: maxDepth=" << maxDepth
//             << " minDepth=" << minDepth
//             << " allowDrop=" << (allowDrop ? 1 : 0) << "\n";
//
//   // Construct the quadtree with explicit minDepth
//   QuadTree2D qt(maxDepth, minDepth);
//
//   double X[9], Y[9];
//   auto phys = [](double xi, double eta) {
//     // identity map (distort later if you want)
//     return std::array<double, 2> {xi, eta};
//   };
//   const double nodePE[9][2] = {
//     {-1, -1}, {+1, -1}, {+1, +1}, {-1, +1}, {0, -1}, {+1, 0}, {0, +1}, {-1, 0}, {0, 0}
//   };
//   for (int a = 0; a < 9; ++a) {
//     auto p = phys(nodePE[a][0], nodePE[a][1]);
//     X[a] = p[0];
//     Y[a] = p[1];
//   }
//   qt.set_physical_quad9(X, Y);
//
//   // Soft floor policy (true = may coarsen below floor where coarsen_pred says it's safe)
//   qt.set_allow_coarsen_below_min(allowDrop);
//
//   // Level set and thresholds
//   Psi psi0;
//   // psi0.xc = 0.6;        // small shift of interface
//   // psi0.yc = 0.5;
//   const double tau_ref    = 2;  // refine if |psi| small or sign-crossing
//   // const double tau_coarse = 2e-3;  // coarsen if strictly away from interface
//   //
//   // auto refine0  = make_refine_pred(psi0, tau_ref, maxDepth);//make_refine_pred(psi0, tau_ref);
//   // auto coarsen0 = make_coarsen_pred(psi0, tau_coarse, /*min_level=*/0);
//
//   // // Example: a refine-only phase (still respects minDepth automatically)
//   // std::size_t r0 = qt.adapt_refine_until(refine0, Basis::Q2_Quad9, /*max_passes=*/10);
//   // std::cout << "Initial refine: " << r0 << " cells\n";
//   // std::cout << "Leaves after init: " << qt.leaf_count() << "\n";
//
//   // Example: a full coarsen+refine cycle phase
//   // std::size_t total = qt.adapt_cycle(coarsen0, refine0, Basis::Q2_Quad9, /*max_passes=*/8);
//   // std::cout << "Cycle changes: " << total << " cells\n";
//   // std::cout << "Leaves after cycle: " << qt.leaf_count() << "\n";
//
//   // --- Create a Q2 field and initialize coeffs from Psi0 at leaf nodes ---
//   const u32 fid = qt.add_field(Basis::Q2_Quad9);
//   //qt.resize_fields_to_leaves();
//   // const auto& L0 = qt.leaf_indices();
//   // for (u32 k = 0; k < (u32)L0.size(); ++k) {
//   //   std::vector<std::array<double, 2>> xy9;
//   //   qt.leaf_physical_nodes(Basis::Q2_Quad9, L0[k], xy9);
//   //   double* c = qt.leaf_coeff_ptr(fid, k);
//   //   for (int a = 0; a < 9; ++a) c[a] = psi0(xy9[a][0], xy9[a][1]);
//   // }
//   //
//   // qt.write_vtu("./output/element_adaptive.0.vtu", fid, "u", /*cell_centered=*/true);
//
//   // --- Now change Psi slightly and run coarsen+refine with hysteresis ---
//
//   for (u32 k = 0; k <= 100; k++) {
//
//     Psi psi1 = psi0;
//     psi1.xc = 0.5 + 0.1 * cos(2. * M_PI / 100. * k);      // small shift of interface
//     psi1.yc = 0.5 + 0.1 * sin(2. * M_PI / 100. * k);      // small shift of interface
//
//     const double tau_coarse = 2e-3;  // coarsen margin ( > tau_ref for hysteresis)
//     const u32    min_level  = 0;     // optionally keep a minimum level
//
//     auto refine1  = make_refine_pred(psi1, tau_ref, maxDepth);
//     auto coarsen1 = make_coarsen_pred(psi1, tau_coarse, min_level);
//
//     std::size_t changed = qt.adapt_cycle(coarsen1, refine1, Basis::Q2_Quad9, /*max_passes=*/10);
//     std::cout << "Adapt-cycle changed " << changed << " cells\n";
//     std::cout << "Leaves after cycle: " << qt.leaf_count() << "\n";
//
//     // --- Re-size and refresh field coefficients for new mesh (from psi1) ---
//     qt.resize_fields_to_leaves();
//     const auto& L1 = qt.leaf_indices();
//     for (u32 k = 0; k < (u32)L1.size(); ++k) {
//       std::vector<std::array<double, 2>> xy9;
//       qt.leaf_physical_nodes(Basis::Q2_Quad9, L1[k], xy9);
//       double* c = qt.leaf_coeff_ptr(fid, k);
//       for (int a = 0; a < 9; ++a) c[a] = psi1(xy9[a][0], xy9[a][1]);
//     }
//
//     std::string filename = "./output/element_adaptive." + std::to_string(k) + ".vtu";
//     qt.write_vtu(filename, fid, "u", /*cell_centered=*/false);
//     std::cout << "Printing " << filename << "\n";
//   }
//
//
//
//
//   // Fill coeffs somehow...
//   qt.gather_field(fid);
//   for (size_t g = 0; g < qt.global_coords.size(); ++g) {
//     auto [x, y] = qt.global_coords[g];
//     qt.global_field[g] = x + y;
//   }
//   qt.scatter_field(fid);
//   std::string filename = "./output/element_adaptive." + std::to_string(101) + ".vtu";
//   qt.write_vtu(filename, fid, "u", /*cell_centered=*/false);
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//   //ProfilerStop();
//
//   return 0;
//
// }










// for (u32 k = 1; k <= 5; k++) {
//
//   auto coords = qt.extract_node_coords_in_level_range(maxDepth, maxDepth, Basis::Q2_Quad9);
//
//   std::cout << "Collected " << coords.size() << " node coordinates\n";
//   for (auto& xy : coords) {
//     xy[0] -= 0.02;
//     xy[1] -= 0.02;
//   }
//
//   qt.reset();
//   qt.refine_to_contain_points(coords, maxDepth);
//
//   Psi psi1 = psi0;
//   psi1.xc = 0.5 - 0.02 * k;      // small shift of interface
//   psi1.yc = 0.5 - 0.02 * k;      // small shift of interface
//
//   const double tau_coarse = 2e-3;  // coarsen margin ( > tau_ref for hysteresis)
//   const u32    min_level  = 0;     // optionally keep a minimum level
//
//   auto refine1  = make_refine_pred(psi1, tau_ref, maxDepth);
//   auto coarsen1 = make_coarsen_pred(psi1, tau_coarse, min_level);
//
//   std::size_t changed = qt.adapt_cycle(coarsen1, refine1, Basis::Q2_Quad9, /*max_passes=*/10);
//   std::cout << "Adapt-cycle changed " << changed << " cells\n";
//   std::cout << "Leaves after cycle: " << qt.leaf_count() << "\n";
//
//   // --- Re-size and refresh field coefficients for new mesh (from psi1) ---
//   qt.resize_fields_to_leaves();
//   const auto& L1 = qt.leaf_indices();
//   for (u32 k = 0; k < (u32)L1.size(); ++k) {
//     std::vector<std::array<double, 2>> xy9;
//     qt.leaf_physical_nodes(Basis::Q2_Quad9, L1[k], xy9);
//     double* c = qt.leaf_coeff_ptr(fid, k);
//     for (int a = 0; a < 9; ++a) c[a] = psi1(xy9[a][0], xy9[a][1]);
//   }
//
//   std::string filename = "./output/element_adaptive." + std::to_string(k) + ".vtu";
//   qt.write_vtu(filename, fid, "u", /*cell_centered=*/false);
//   std::cout << "Printing " << filename << "\n";
// }
//
// //ProfilerStop();
//
// return 0;

//}



// #include "QuadTree2D.hpp"
// #include <iostream>
// #include <vector>
// #include <array>
// #include <cmath>
//
// int main() {
//     using namespace fem;
//
//     // --- 1. Build tree with max depth 6
//     QuadTree2D tree(10);
//
//     // --- 2. Geometry (parent [-1,1]^2 mapped to same physical square)
//     double X9[9] = {-1, 1, 1, -1, 0, 1, 0, -1, 0};
//     double Y9[9] = {-1,-1, 1,  1,-1, 0, 1,  0, 0};
//     tree.set_physical_quad9(X9, Y9);
//
//     // --- 3. Pick some points to force refinement around
//     std::vector<std::array<double,2>> pts = {
//         {0.2, 0.2},
//         {-0.5, 0.7},
//         {0.8, -0.4}
//     };
//     tree.refine_to_contain_points(pts, 10); // ensure depth >= 4
//
//     // --- 4. Add a fake field
//     u32 fid = tree.add_field(Basis::Q2_Quad9);
//     tree.resize_fields_to_leaves();
//
//     // --- 5. Fill with fake values (here: v = sin(pi*x)*cos(pi*y))
//     const auto& L = tree.leaf_indices();
//     for (u32 k = 0; k < (u32)L.size(); ++k) {
//         std::vector<std::array<double,2>> xy;
//         tree.leaf_physical_nodes(Basis::Q2_Quad9, L[k], xy);
//         double* coeffs = tree.leaf_coeff_ptr(fid, k);
//         for (int i = 0; i < 9; ++i) {
//             double x = xy[i][0], y = xy[i][1];
//             coeffs[i] = std::sin(M_PI * x) * std::cos(M_PI * y);
//         }
//     }
//
//     // --- 6. Dump to VTU
//     if (tree.write_vtu("fake_field.vtu", fid, "fake_field", false)) {
//         std::cout << "Wrote fake_field.vtu" << std::endl;
//     } else {
//         std::cerr << "Failed to write VTU" << std::endl;
//     }
//
//     return 0;
// }
