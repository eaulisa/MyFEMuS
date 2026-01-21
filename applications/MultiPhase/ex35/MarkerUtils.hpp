#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <locale>


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

  static inline double dot(const Vector& v, const Vector& w) {
    double d = 0.0;
    for (std::size_t idim = 0; idim < DIM; ++idim) d += v[idim]*w[idim];
    return d;
  }
};

template <std::size_t DIM>
using Point  = typename Geom_op<DIM>::Point;

template <std::size_t DIM>
using Vector = typename Geom_op<DIM>::Vector;

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


    static inline double dot3(const Point<3>& a, const Point<3>& b){
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    static inline void sub3(const Point<3>& a, const Point<3>& b, Point<3>& r){
        r[0]=a[0]-b[0]; r[1]=a[1]-b[1]; r[2]=a[2]-b[2];
    }
    static inline void madd3(Point<3>& x, const Point<3>& v, double s){
        x[0]+=s*v[0]; x[1]+=s*v[1]; x[2]+=s*v[2];
    }
    static inline void lerp3(const Point<3>& P, const Point<3>& Q, double a, Point<3>& R){
        R[0]=(1.0-a)*P[0]+a*Q[0]; R[1]=(1.0-a)*P[1]+a*Q[1]; R[2]=(1.0-a)*P[2]+a*Q[2];
    }
    static inline double norm3(const Point<3>& v){ return std::sqrt(dot3(v,v)); }
    static inline void cross3(const Point<3>& a, const Point<3>& b, Point<3>& c){
        c[0]=a[1]*b[2]-a[2]*b[1]; c[1]=a[2]*b[0]-a[0]*b[2]; c[2]=a[0]*b[1]-a[1]*b[0];
    }

    void sample_triangle_with_density(const Tri3& T,
                                    const double& den,
                                    const int& min_segments,
                                    std::vector<Point<3>>& out){
        Point<3> ABv, BCv, CAv;
        sub3(T.B, T.A, ABv);
        sub3(T.C, T.B, BCv);
        sub3(T.A, T.C, CAv);
        const double LAB = norm3(ABv);
        const double LBC = norm3(BCv);
        const double LCA = norm3(CAv);

        const double rho = std::max(0.0, den);
        const double s_line = (rho>0.0) ? 1.0/rho : 0.0;

        const int N_AB = std::max((int)std::lround(rho*LAB), min_segments+1);
        const int N_BC = std::max((int)std::lround(rho*LBC), min_segments+1);
        const int N_CA = std::max((int)std::lround(rho*LCA), min_segments+1);

        for (int i=0; i<=N_AB; ++i){
            double t = (N_AB==0)?0.5 : (double)i/N_AB;
            Point<3> P; lerp3(T.A, T.B, t, P);
            out.push_back(P);
        }
        for (int i=1; i<=N_BC; ++i){
            double t = (N_BC==0)?0.5 : (double)i/N_BC;
            Point<3> P; lerp3(T.B, T.C, t, P);
            out.push_back(P);
        }
        for (int i=1; i<=N_CA; ++i){
            double t = (N_CA==0)?0.5 : (double)i/N_CA;
            Point<3> P; lerp3(T.C, T.A, t, P);
            out.push_back(P);
        }

        const Point<3> *P, *Q, *R;
        Point<3> baseQP;
        double   Lbase;
        if (LBC >= LAB && LBC >= LCA) { P=&T.B; Q=&T.C; R=&T.A; Lbase=LBC; baseQP=BCv; }
        else if (LAB >= LBC && LAB >= LCA) { P=&T.A; Q=&T.B; R=&T.C; Lbase=LAB; baseQP=ABv; }
        else { P=&T.C; Q=&T.A; R=&T.B; Lbase=LCA; baseQP=CAv; }

        Point<3> RP, QP;
        sub3(*R, *P, RP);
        sub3(*Q, *P, QP); 
        Point<3> cr; cross3(QP, RP, cr);
        const double h = norm3(cr) / std::max(1e-30, Lbase);
        if (Lbase<=0.0 || h<=0.0) return;

        int nT = 0;
        if (s_line>0.0) nT = std::max(1, (int)std::ceil(h / s_line));
        else            nT = std::max(1, std::max(N_AB, std::max(N_BC, N_CA)));

        if (rho>0.0){
            const double sum_t = (nT>1) ? ( (nT-1)*nT*0.5 / nT ) : 0.0;
            const size_t est_internal = (size_t)std::max(0.0, rho * Lbase * sum_t);
            out.reserve(out.size() + est_internal + 8);
        } else {
            out.reserve(out.size() + (size_t)(nT * std::max({N_AB,N_BC,N_CA})) + 8);
        }

        Point<3> RtoP, RtoQ;
        sub3(*P, *R, RtoP);
        sub3(*Q, *R, RtoQ);

        for (int k=1; k<=nT-1; ++k){
            const double t = (double)k / nT;
            const double ell = t * Lbase;

            int nk = 1;
            if (rho>0.0) nk = std::max(1, (int)std::lround(ell * rho));
            else         nk = std::max(1, (int)std::lround((ell / Lbase) * std::max(N_AB, std::max(N_BC, N_CA))));

            Point<3> segDir; sub3(RtoQ, RtoP, segDir); 
            for (int j=0; j<nk; ++j){
                const double a = (nk==1) ? 0.5 : (double)j/(nk-1);
                Point<3> X = *R;
                madd3(X, RtoP, t);
                madd3(X, segDir, t*a);
                out.push_back(X);
            }
        }
    }



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

namespace MarkerUtils {

inline void writeMarkersCSV(const std::string& filepath,
                            const std::vector<std::vector<double>>& markers,
                            bool include_index = true)
{
    const std::size_t dim = markers.size();
    if (dim == 0) {
        throw std::runtime_error("writeMarkersCSV: markers has dim=0");
    }
    if (dim > 3) {
        throw std::runtime_error("writeMarkersCSV: dim > 3 not supported");
    }

    const std::size_t N = markers[0].size();
    for (std::size_t a = 1; a < dim; ++a) {
        if (markers[a].size() != N) {
            throw std::runtime_error("writeMarkersCSV: inconsistent sizes across dimensions");
        }
    }

    std::ofstream out(filepath);
    if (!out.is_open()) {
        throw std::runtime_error("writeMarkersCSV: cannot open file: " + filepath);
    }

    // Forza il separatore decimale '.' indipendentemente dalla locale di sistema
    out.imbue(std::locale::classic());
    out << std::setprecision(17);

    // Header
    if (include_index) out << "idx,";
    if (dim >= 1) out << "x";
    if (dim >= 2) out << ",y";
    if (dim >= 3) out << ",z";
    out << "\n";

    // Rows: una riga per marker
    for (std::size_t i = 0; i < N; ++i) {
        if (include_index) out << i << ",";
        for (std::size_t a = 0; a < dim; ++a) {
            out << markers[a][i];
            if (a + 1 < dim) out << ",";
        }
        out << "\n";
    }

    out.flush();
    if (!out) {
        throw std::runtime_error("writeMarkersCSV: error while writing file: " + filepath);
    }
}

} // namespace Marker3DUtils