#ifndef __femus_classify_hpp__
#define __femus_classify_hpp__

#include <array>
#include <cmath>
#include <algorithm>

// ---------- Public API ----------
enum class ElemKind { Hex, Tet, Wedge, Quad, Tri, Line };

struct Options {
  double base_eps = 1.0e-12;     // unit-size base tolerance
  bool   normalize_by_normal = true;
  double tie_tol = 1e-12;        // relative tie tolerance for “equal” violations
};

struct Classify {
  bool inside = true;  // true if inside
  int  face  = -1;     // most-likely outside face index (0-based), or -1 if inside
  int  face2 = -1;     // second-most-likely face index (distinct from face), or -1 if none
};

// ---------- Internal ----------
namespace detail {

  struct Face {
    // s(x) = ax*x + ay*y + az*z + b ; "inside" means s(x) <= 0
    double ax = 0, ay = 0, az = 0, b = 0;
    int idx = -1;  // face index (user-specified order)
    double nrm() const {
      return std::sqrt(ax * ax + ay * ay + az * az);
    }
    double s(double x, double y = 0, double z = 0) const {
      return ax * x + ay * y + az * z + b;
    }
  };

  inline double scaled_eps(double base_eps, double x, double y = 0, double z = 0) {
    double m = std::max({1.0, std::fabs(x), std::fabs(y), std::fabs(z)});
    return base_eps * m;
  }

  inline bool greater_with_tol(double a, double b, double rel_tol) {
    // a > b by a relative margin
    double scale = std::max(1.0, std::fabs(b));
    return (a - b) > rel_tol * scale;
  }
  inline bool approx_equal(double a, double b, double rel_tol) {
    double scale = std::max({1.0, std::fabs(a), std::fabs(b)});
    return std::fabs(a - b) <= rel_tol * scale;
  }

  template <class FaceList>
  Classify classify_common(const FaceList& faces, double x, double y, double z, const Options& opt) {
    const double eps = scaled_eps(opt.base_eps, x, y, z);

    // Track top-2 violating faces (normalized)
    bool   any_out   = false;
    double best_v    = -1.0;
    int best_i    = -1;
    double second_v  = -1.0;
    int second_i  = -1;

    for (const auto& f : faces) {
      const double val  = f.s(x, y, z);
      const double denom = opt.normalize_by_normal ? std::max(f.nrm(), 1e-300) : 1.0;
      const double vn   = val / denom;     // distance-like
      const double viol = vn - eps;        // positive => outside

      if (viol <= 0.0) continue;           // inside w.r.t. this face

      if (!any_out) {
        any_out = true;
        best_v = viol;
        best_i = f.idx;
        continue;
      }

      // Compare against current best
      if (greater_with_tol(viol, best_v, opt.tie_tol)) {
        // New best; demote old best to second
        second_v = best_v;
        second_i = best_i;
        best_v   = viol;
        best_i   = f.idx;
      }
      else if (approx_equal(viol, best_v, opt.tie_tol)) {
        // Tie with best: choose smallest index as best; the other becomes candidate for second
        if (f.idx != best_i) {
          int lo = std::min(f.idx, best_i);
          int hi = std::max(f.idx, best_i);
          // Best is smaller index
          if (lo != best_i) {
            // swap: put current as best, old best becomes second
            second_v = best_v;
            second_i = best_i;
            best_i = lo;
            best_v = viol; // viol == best_v
          }
          else {
            // keep best as is; use current as second candidate
            if (second_i == -1 || greater_with_tol(viol, second_v, opt.tie_tol) ||
                (approx_equal(viol, second_v, opt.tie_tol) && f.idx < second_i)) {
              second_v = viol;
              second_i = f.idx;
            }
          }
          // If second equals best, fix it
          if (second_i == best_i) second_i = hi;
        }
      }
      else {
        // Between best and second
        if (second_i == -1 || greater_with_tol(viol, second_v, opt.tie_tol) ||
            (approx_equal(viol, second_v, opt.tie_tol) && f.idx < second_i)) {
          if (f.idx != best_i) {
            second_v = viol;
            second_i = f.idx;
          }
        }
      }
    }

    if (!any_out) return {true, -1, -1};

    // Ensure second is distinct and valid; if not, leave as -1
    if (second_i == best_i) second_i = -1;
    return {false, best_i, second_i};
  }

} // namespace detail

// ===================================================================
// Element 0: HEX [-1,1]^3
// Face order:
//  0: front  (y=-1)
//  1: right  (x=1)
//  2: back   (y=1)
//  3: left   (x=-1)
//  4: bottom (z=-1)
//  5: top    (z=1)
inline Classify classify_hex(const std::array<double, 3>& X, const Options& opt = {}) {
  using namespace detail;
  static const std::array<Face, 6> faces = {
    Face{  0, -1,  0, -1, 0},  // y >= -1
    Face{ +1,  0,  0, -1, 1},  // x <=  1
    Face{  0, +1,  0, -1, 2},  // y <=  1
    Face{ -1,  0,  0, -1, 3},  // x >= -1
    Face{  0,  0, -1, -1, 4},  // z >= -1
    Face{  0,  0, +1, -1, 5}   // z <=  1
  };
  return classify_common(faces, X[0], X[1], X[2], opt);
}

// ===================================================================
// Element 1: TET (0,0,0)-(1,0,0)-(0,1,0)-(0,0,1)
// Face order: 0: z=0 (bottom), 1: y=0 (front), 2: 1-x-y-z=0, 3: x=0 (left)
inline Classify classify_tet(const std::array<double, 3>& X, const Options& opt = {}) {
  using namespace detail;
  static const std::array<Face, 4> faces = {
    Face{  0,  0, -1,  0, 0},   // z >= 0
    Face{  0, -1,  0,  0, 1},   // y >= 0
    Face{ +1, +1, +1, -1, 2},   // x+y+z <= 1
    Face{ -1,  0,  0,  0, 3}    // x >= 0
  };
  return classify_common(faces, X[0], X[1], X[2], opt);
}

// ===================================================================
// Element 2: WEDGE (x>=0, y>=0, x+y<=1, -1<=z<=1)
// Face order: 0:y=0, 1:1-x-y=0, 2:x=0, 3:z=-1, 4:z=1
inline Classify classify_wedge(const std::array<double, 3>& X, const Options& opt = {}) {
  using namespace detail;
  static const std::array<Face, 5> faces = {
    Face{  0, -1,  0,  0, 0},   // y >= 0
    Face{ +1, +1,  0, -1, 1},   // x+y <= 1
    Face{ -1,  0,  0,  0, 2},   // x >= 0
    Face{  0,  0, -1, -1, 3},   // z >= -1
    Face{  0,  0, +1, -1, 4}    // z <= 1
  };
  return classify_common(faces, X[0], X[1], X[2], opt);
}

// ===================================================================
// Element 3: QUAD [-1,1]^2
// Edge order: 0:y=-1, 1:x=1, 2:y=1, 3:x=-1
inline Classify classify_quad(const std::array<double, 2>& X, const Options& opt = {}) {
  using namespace detail;
  static const std::array<Face, 4> faces = {
    Face{  0, -1, 0, -1, 0},   // y >= -1
    Face{ +1,  0, 0, -1, 1},   // x <=  1
    Face{  0, +1, 0, -1, 2},   // y <=  1
    Face{ -1,  0, 0, -1, 3}    // x >= -1
  };
  return classify_common(faces, X[0], X[1], 0.0, opt);
}

// ===================================================================
// Element 4: TRI (x>=0, y>=0, x+y<=1)
// Face order: 0:y=0, 1:1-x-y=0, 2:x=0
inline Classify classify_tri(const std::array<double, 2>& X, const Options& opt = {}) {
  using namespace detail;
  static const std::array<Face, 3> faces = {
    Face{  0, -1, 0,  0, 0},   // y >= 0
    Face{ +1, +1, 0, -1, 1},   // x+y <= 1
    Face{ -1,  0, 0,  0, 2}    // x >= 0
  };
  return classify_common(faces, X[0], X[1], 0.0, opt);
}

// ===================================================================
// Element 5: LINE [-1,1]
// Point order: 0:x=-1, 1:x=1
inline Classify classify_line(double x, const Options& opt = {}) {
  using namespace detail;
  static const std::array<Face, 2> faces = {
    Face{ -1, 0, 0, -1, 0},   // x >= -1
    Face{ +1, 0, 0, -1, 1}    // x <=  1
  };
  return classify_common(faces, x, 0.0, 0.0, opt);
}

// ---------- Dispatcher ----------
inline Classify classify_point(ElemKind kind, const double* x, const Options& opt = {}) {
  switch (kind) {
    case ElemKind::Hex:
      return classify_hex({x[0], x[1], x[2]}, opt);
    case ElemKind::Tet:
      return classify_tet({x[0], x[1], x[2]}, opt);
    case ElemKind::Wedge:
      return classify_wedge({x[0], x[1], x[2]}, opt);
    case ElemKind::Quad:
      return classify_quad({x[0], x[1]},     opt);
    case ElemKind::Tri:
      return classify_tri({x[0], x[1]},      opt);
    case ElemKind::Line:
      return classify_line(x[0],            opt);
    default:
      return {};
  }
}


inline Classify classify_point(unsigned elementType, const double* x, const Options& opt = {}) {
  switch (elementType) {
    case 0:
      return classify_hex({x[0], x[1], x[2]}, opt);
    case 1:
      return classify_tet({x[0], x[1], x[2]}, opt);
    case 2:
      return classify_wedge({x[0], x[1], x[2]}, opt);
    case 3:
      return classify_quad({x[0], x[1]},     opt);
    case 4:
      return classify_tri({x[0], x[1]},      opt);
    case 5:
      return classify_line(x[0],            opt);
    default:
      return {};
  }
}
#endif
