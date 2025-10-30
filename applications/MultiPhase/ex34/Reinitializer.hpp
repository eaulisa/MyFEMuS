#pragma once

#include <cstdint>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <cstddef>
#include <cassert>

#include "OctTree.hpp"
#include "nanoflann.hpp"

template <std::size_t DIM>
struct Geom_op {
  using Point  = std::array<double, DIM>;
  using Vector = std::array<double, DIM>;

  static inline Point add(const Point& a, const Vector& v) {
    Point r{};
    for (std::size_t d = 0; d < DIM; ++d) r[d] = a[d] + v[d];
    return r;
  }

  static inline Point sub(const Point& a, const Vector& v) {
    Point r{};
    for (std::size_t d = 0; d < DIM; ++d) r[d] = a[d] - v[d];
    return r;
  }

  static inline Vector mul(double c, const Vector& v) {
    Vector r{};
    for (std::size_t d = 0; d < DIM; ++d) r[d] = c * v[d];
    return r;
  }

  static inline Vector div(const Vector& v, double c) {
    Vector r{};
    for (std::size_t d = 0; d < DIM; ++d) r[d] = v[d] / c;
    return r;
  }

  static inline double distance2(const Point& a, const Point& b) {
    double s = 0.0;
    for (std::size_t d = 0; d < DIM; ++d) {
      const double t = a[d] - b[d];
      s += t * t;
    }
    return s;
  }

  static inline double distance(const Point& a, const Point& b) {
    return std::sqrt(distance2(a,b));
  }

  static inline double norm2(const Vector& v) {
    double s = 0.0;
    for (std::size_t d = 0; d < DIM; ++d) s += v[d] * v[d];
    return s;
  }

  static inline double norm(const Vector& v) {
    return std::sqrt(norm2(v));
  }
};

template <std::size_t DIM>
using Point  = typename Geom_op<DIM>::Point;

template <std::size_t DIM>
using Vector = typename Geom_op<DIM>::Vector;


namespace fem {
    using u32 = uint32_t;
    using u64 = uint64_t;

    template<std::size_t DIM>
    class OctTree;

    template<std::size_t DIM>
    class Reinitializer
    {
    public:
        explicit Reinitializer(OctTree<DIM>* tree_ptr, u32 fid_phi, bool flag = true, double density = 10.);

        // this routine reinitializes the level set field fid
        void compute_signed_distance();

        // this routine creates a csv with all markers
        void write_markers_csv(const std::string&) const;

    private:
        // this routine identifies cut cells and compute cell intersections with the interface
        bool leaf_is_cut(u32 leaf_pos, std::vector<Point<DIM>>* cell_intersections = nullptr);
        
        // this routine compute markers on the cell interface
        void compute_cell_markers(std::vector<Point<DIM>> cell_intersections);

        // this routine collect markers on the level set interface
        void collect_markers();

        // this routine project the cut cells nodes on to the interface
        void project_cut_cells_nodes();

        // this routine find the root on a intersected edge
        std::vector<double> edge_roots(double v0, double v1, double v2);

        // this routines compute leaf quantities
        Point<DIM>  evaluate_coord_on_leaf(const Point<DIM> P_local);
        double      evaluate_field_on_leaf(const Point<DIM> P_local);
        Vector<DIM> evaluate_gradient_on_leaf(const Point<DIM> P);


        OctTree<DIM> * tree = nullptr;
        std::vector<Point<DIM>> markers;
        std::vector<u32> cut_cells;
        std::vector<std::vector<double>> cut_cells_nodes_dist;
        u32 leaf_id;
        u32 fid;

        std::vector<Point<DIM>> local_coord;
        std::vector<double> local_values;

        bool proj_flag;
        double marker_density;
    };
}

// adaptor for nanoflann
template <std::size_t DIM>
struct PointCloud {
    using PointT = std::array<double, DIM>;
    std::vector<PointT> pts;

    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline double kdtree_distance(const double* query, const size_t idx, size_t /*dim*/) const {
        double acc = 0.0;
        const auto& p = pts[idx];
        for (std::size_t d = 0; d < DIM; ++d) {
            const double t = query[d] - p[d];
            acc += t * t;
        }
        return acc;
    }

    inline double kdtree_get_pt(const size_t idx, size_t dim) const {
        return pts[idx][dim];
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

// wrapper kdtree
template <std::size_t DIM>
class KDTree {
public:
    static_assert(DIM == 2 || DIM == 3, "KDTree supports only DIM=2 or 3");

    using Adaptor  = PointCloud<DIM>;
    using Metric   = nanoflann::L2_Simple_Adaptor<double, Adaptor>;
    using KDIndex  = nanoflann::KDTreeSingleIndexAdaptor<Metric, Adaptor, int(DIM)>; 

    KDTree(unsigned leaf_max = 10)
    : index_(int(DIM), cloud_, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max)) {}

    void build(const std::vector<std::array<double, DIM>>& points) {
        cloud_.pts = points;
        index_.buildIndex();
        built_ = true;
    }

    template <class PointLike>
    void build_from_any(const std::vector<PointLike>& points) {
        cloud_.pts.resize(points.size());
        for (size_t i = 0; i < points.size(); ++i) {
            for (std::size_t d = 0; d < DIM; ++d) cloud_.pts[i][d] = points[i][d];
        }
        index_.buildIndex();
        built_ = true;
    }

    void rebuild() {
        index_.buildIndex();
        built_ = true;
    }

    std::pair<std::vector<size_t>, std::vector<double>>
    knn(const std::array<double, DIM>& q, size_t k) const {
        assert(built_ && "KDTree: call build() before querying");
        std::vector<size_t> idx(k);
        std::vector<double> dist2(k);
        nanoflann::KNNResultSet<double> rs(k);
        rs.init(idx.data(), dist2.data());
        index_.findNeighbors(rs, q.data(), nanoflann::SearchParameters());
        return {idx, dist2}; // distanze al quadrato
    }

    std::vector<std::pair<size_t,double>>
    radius(const std::array<double, DIM>& q, double r,
           size_t maxMatches = std::numeric_limits<size_t>::max()) const
    {
        assert(built_ && "KDTree: call build() before querying");
        using IndexType    = typename KDIndex::IndexType;
        using DistanceType = typename KDIndex::DistanceType;

        std::vector< nanoflann::ResultItem<IndexType, DistanceType> > hits;
        nanoflann::SearchParameters params;
        const DistanceType r2 = static_cast<DistanceType>(r*r);

        index_.radiusSearch(q.data(), r2, hits, params);

        if (maxMatches != std::numeric_limits<size_t>::max() && hits.size() > maxMatches)
            hits.resize(maxMatches);

        std::vector<std::pair<size_t,double>> out;
        out.reserve(hits.size());
        for (const auto& h : hits)
            out.emplace_back(static_cast<size_t>(h.first), static_cast<double>(h.second)); // dist²
        return out;
    }

    const std::vector<std::array<double, DIM>>& points() const { return cloud_.pts; }

private:
    Adaptor cloud_;
    KDIndex index_;
    bool    built_{false};
};

// some routines for marker placement in 3D
namespace Tri_ord{
  constexpr double EPS_LEN2   = 1e-28; 
  constexpr double EPS_NORMAL = 1e-14;  

  template <std::size_t D>
  inline double dot(const Point<D>& a, const Vector<D>& b){
      double s=0.0; for (std::size_t i=0;i<D;++i) s+=a[i]*b[i]; return s;
  }
  template <std::size_t D>
  inline double dotv(const Vector<D>& a, const Vector<D>& b){
      double s=0.0; for (std::size_t i=0;i<D;++i) s+=a[i]*b[i]; return s;
  }
  inline Vector<3> cross(const Vector<3>& a, const Vector<3>& b){
      return { a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] };
  }
  template <std::size_t D>
  inline double norm2(const Vector<D>& a){ return dotv<D>(a,a); }
  template <std::size_t D>
  inline double norm(const Vector<D>& a){ return std::sqrt(norm2<D>(a)); }

  template <std::size_t D>
  inline Vector<D> to_vec(const Point<D>& a){ Vector<D> r; for (std::size_t i=0;i<D;++i) r[i]=a[i]; return r; }
  template <std::size_t D>
  inline Point<D> to_pt(const Vector<D>& a){ Point<D> r; for (std::size_t i=0;i<D;++i) r[i]=a[i]; return r; }

  struct Tri3 { Point<3> A,B,C; };

  auto tri_area2 = [](const Point<3>& A, const Point<3>& B, const Point<3>& C){
      Vector<3> u{B[0]-A[0], B[1]-A[1], B[2]-A[2]};
      Vector<3> v{C[0]-A[0], C[1]-A[1], C[2]-A[2]};
      Vector<3> w{ u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0] };
      return std::sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  };

  auto triangulate_poly = [&](const std::vector<Point<3>>& poly) {
      std::vector<Tri3> tris;
      const std::size_t m = poly.size();
      if (m < 3) return tris;

      if (m == 3) {
          tris.push_back({poly[0], poly[1], poly[2]});
          return tris;
      }
      if (m == 4) {
          auto d02 = (poly[0][0]-poly[2][0])*(poly[0][0]-poly[2][0]) +
                  (poly[0][1]-poly[2][1])*(poly[0][1]-poly[2][1]) +
                  (poly[0][2]-poly[2][2])*(poly[0][2]-poly[2][2]);
          auto d13 = (poly[1][0]-poly[3][0])*(poly[1][0]-poly[3][0]) +
                  (poly[1][1]-poly[3][1])*(poly[1][1]-poly[3][1]) +
                  (poly[1][2]-poly[3][2])*(poly[1][2]-poly[3][2]);
          if (d02 <= d13) {
              tris.push_back({poly[0], poly[1], poly[2]});
              tris.push_back({poly[0], poly[2], poly[3]});
          } else {
              tris.push_back({poly[1], poly[2], poly[3]});
              tris.push_back({poly[1], poly[3], poly[0]});
          }
          return tris;
      }

      Point<3> C{0,0,0};
      for (auto& p : poly) { C[0]+=p[0]; C[1]+=p[1]; C[2]+=p[2]; }
      C[0]/=m; C[1]/=m; C[2]/=m;

      tris.reserve(m);
      for (std::size_t i=0; i<m; ++i) {
          std::size_t j = (i+1)%m;
          if (tri_area2(C, poly[i], poly[j]) > 1e-16)
              tris.push_back({C, poly[i], poly[j]});
      }
      return tris;
  };


  auto sample_triangle = [&](const Tri3& T, int res, std::vector<Point<3>>& out) {
      for (int i=0; i<=res; ++i) {
          for (int j=0; j<=res-i; ++j) {
              double u = static_cast<double>(i)/res;
              double v = static_cast<double>(j)/res;
              double w = 1.0 - u - v;
              Point<3> P;
              for (int d=0; d<3; ++d)
                  P[d] = w*T.A[d] + u*T.B[d] + v*T.C[d];
              out.push_back(P);
          }
      }
  };



  template <std::size_t D>
  std::vector<Point<D>> dedup_nearby(const std::vector<Point<D>>& pts, double eps2 = EPS_LEN2) {
      if (pts.empty()) return {};
      std::vector<Point<D>> out; out.reserve(pts.size());
      for (const auto& p : pts) {
          bool dup=false;
          for (const auto& q : out) {
              double s=0.0; for (std::size_t i=0;i<D;++i){ double e=p[i]-q[i]; s+=e*e; }
              if (s < eps2) { dup=true; break; }
          }
          if (!dup) out.push_back(p);
      }
      return out;
  }

  inline std::vector<Point<3>> order_polygon_3d(const std::vector<Point<3>>& pts_in) {
      auto pts = dedup_nearby<3>(pts_in);
      if (pts.size() <= 2) return pts;

      Point<3> C{0.0,0.0,0.0};
      for (auto& p: pts){ C[0]+=p[0]; C[1]+=p[1]; C[2]+=p[2]; }
      C[0]/=pts.size(); C[1]/=pts.size(); C[2]/=pts.size();

      Vector<3> n{0.0,0.0,0.0};
      for (std::size_t i=0;i<pts.size();++i){
          const auto& p = pts[i];
          const auto& q = pts[(i+1)%pts.size()];
          n[0] += (p[1]-q[1])*(p[2]+q[2]);
          n[1] += (p[2]-q[2])*(p[0]+q[0]);
          n[2] += (p[0]-q[0])*(p[1]+q[1]);
      }
      double nn = norm<3>(n);
      if (nn < EPS_NORMAL) {
          Vector<3> ranges{0,0,0};
          Point<3>  minp = pts[0], maxp = pts[0];
          for (auto& p: pts){
              for (int d=0; d<3; ++d){ minp[d]=std::min(minp[d],p[d]); maxp[d]=std::max(maxp[d],p[d]); }
          }
          for (int d=0; d<3; ++d) ranges[d] = maxp[d]-minp[d];
          int drop = (ranges[0] > ranges[1]) ? (ranges[0] > ranges[2] ? 0 : 2) : (ranges[1] > ranges[2] ? 1 : 2);
          struct AngIdx { std::size_t i; double ang; };
          std::vector<AngIdx> ord; ord.reserve(pts.size());
          for (std::size_t i=0;i<pts.size();++i){
              double x = pts[i][(drop+1)%3]-C[(drop+1)%3];
              double y = pts[i][(drop+2)%3]-C[(drop+2)%3];
              ord.push_back({i, std::atan2(y,x)});
          }
          std::sort(ord.begin(), ord.end(), [](const AngIdx& a, const AngIdx& b){ return a.ang < b.ang; });
          std::vector<Point<3>> poly; poly.reserve(pts.size());
          for (auto& ai : ord) poly.push_back(pts[ai.i]);
          return poly;
      }
      n[0]/=nn; n[1]/=nn; n[2]/=nn;

      Vector<3> t{1.0,0.0,0.0};
      if (std::fabs(n[0]) > 0.9) t = {0.0,1.0,0.0};

      double tn = dotv<3>(t,n);
      Vector<3> u{ t[0]-tn*n[0], t[1]-tn*n[1], t[2]-tn*n[2] };
      double un = norm<3>(u); if (un < EPS_NORMAL) u = {0.0,1.0,0.0}, un = 1.0;
      u[0]/=un; u[1]/=un; u[2]/=un;

      Vector<3> v = cross(n,u);

      struct AngIdx { std::size_t i; double ang; };
      std::vector<AngIdx> ord; ord.reserve(pts.size());
      for (std::size_t i=0;i<pts.size();++i){
          Vector<3> d{ pts[i][0]-C[0], pts[i][1]-C[1], pts[i][2]-C[2] };
          double x = dotv<3>(d,u);
          double y = dotv<3>(d,v);
          ord.push_back({i, std::atan2(y,x)});
      }
      std::sort(ord.begin(), ord.end(), [](const AngIdx& a, const AngIdx& b){ return a.ang < b.ang; });

      std::vector<Point<3>> poly; poly.reserve(pts.size());
      for (auto& ai : ord) poly.push_back(pts[ai.i]);

      double area2 = 0.0;
      for (std::size_t i=0;i<poly.size();++i){
          Vector<3> di{ poly[i][0]-C[0], poly[i][1]-C[1], poly[i][2]-C[2] };
          Vector<3> dj{ poly[(i+1)%poly.size()][0]-C[0], poly[(i+1)%poly.size()][1]-C[1], poly[(i+1)%poly.size()][2]-C[2] };
          double xi = dotv<3>(di,u), yi = dotv<3>(di,v);
          double xj = dotv<3>(dj,u), yj = dotv<3>(dj,v);
          area2 += (xi*yj - xj*yi);
      }
      if (area2 < 0.0) std::reverse(poly.begin(), poly.end());

      return poly;
  }
}

#include "Reinitializer.tpp"