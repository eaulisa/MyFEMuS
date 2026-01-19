#include "include/nanoflann.hpp"

#include <vector>
#include <array>
#include <cassert>
#include <stdexcept>
#include <limits>
#include <utility>

/// Adaptor for nanoflann with AoS layout: pts[idx][d], pts is array<3>
struct PointCloudAoS3
{
    int dim = 0;                                  // 2 or 3
    std::vector<std::array<double,3>> pts;        // pts[idx] = {x,y,z}

    std::size_t kdtree_get_point_count() const
    {
        return pts.size();
    }

    double kdtree_get_pt(const std::size_t idx, const std::size_t d) const
    {
        // nanoflann will only ask d in [0..dim-1]
        return pts[idx][d];
    }

    // Optional, but consistent
    double kdtree_distance(const double* query, const std::size_t idx, std::size_t /*unused*/) const
    {
        double acc = 0.0;
        const std::size_t D = static_cast<std::size_t>(dim);
        for (std::size_t d = 0; d < D; ++d)
        {
            const double t = query[d] - pts[idx][d];
            acc += t * t;
        }
        return acc;
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

class KDTreeRT
{
public:
    using Adaptor = PointCloudAoS3;
    using Metric  = nanoflann::L2_Simple_Adaptor<double, Adaptor>;
    using KDIndex = nanoflann::KDTreeSingleIndexAdaptor<Metric, Adaptor, -1 /* dynamic DIM */>;

    explicit KDTreeRT(int dim, std::size_t leaf_max = 10)
    : dim_(dim),
      cloud_{dim, {}},
      index_(dim, cloud_, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max))
    {
        if (dim_ != 2 && dim_ != 3)
            throw std::invalid_argument("KDTreeRT supports only dim=2 or dim=3");
    }

    int  dim()   const noexcept { return dim_; }
    bool built() const noexcept { return built_; }

    /// Build from AoS: points[idx] = array<double,3>
    void build(const std::vector<std::array<double,3>>& points)
    {
        cloud_.dim = dim_;
        cloud_.pts = points;

        index_.buildIndex();
        built_ = true;
    }

    /// If you update cloud_.pts externally, call rebuild()
    void rebuild()
    {
        if (cloud_.pts.empty())
            throw std::runtime_error("KDTreeRT::rebuild: empty point set");

        index_.buildIndex();
        built_ = true;
    }

    std::pair<std::vector<std::size_t>, std::vector<double>>
    knn(const double* q, std::size_t k) const
    {
        if (!built_) throw std::runtime_error("KDTreeRT::knn: call build() before querying");

        const std::size_t N = cloud_.pts.size();
        if (N == 0 || k == 0) return {{}, {}};
        if (k > N) k = N;

        std::vector<std::size_t> idx(k);
        std::vector<double> dist2(k);

        nanoflann::KNNResultSet<double> rs(k);
        rs.init(idx.data(), dist2.data());

        index_.findNeighbors(rs, q, nanoflann::SearchParameters());
        return {idx, dist2};
    }

    std::pair<std::vector<std::size_t>, std::vector<double>>
    knn(const std::array<double,3>& q, std::size_t k) const
    {
        // dim_=2: usa q[0], q[1]; dim_=3: usa q[0], q[1], q[2]
        return knn(q.data(), k);
    }

    std::vector<std::pair<std::size_t,double>>
    radius(const double* q, double r,
           std::size_t maxMatches = std::numeric_limits<std::size_t>::max()) const
    {
        if (!built_) throw std::runtime_error("KDTreeRT::radius: call build() before querying");

        using IndexType    = typename KDIndex::IndexType;
        using DistanceType = typename KDIndex::DistanceType;

        std::vector< nanoflann::ResultItem<IndexType, DistanceType> > hits;

        nanoflann::SearchParameters params;
        const DistanceType r2 = static_cast<DistanceType>(r * r);

        index_.radiusSearch(q, r2, hits, params);

        if (maxMatches != std::numeric_limits<std::size_t>::max() && hits.size() > maxMatches)
            hits.resize(maxMatches);

        std::vector<std::pair<std::size_t,double>> out;
        out.reserve(hits.size());
        for (const auto& h : hits)
            out.emplace_back(static_cast<std::size_t>(h.first), static_cast<double>(h.second));

        return out;
    }

    std::vector<std::pair<std::size_t,double>>
    radius(const std::array<double,3>& q, double r,
           std::size_t maxMatches = std::numeric_limits<std::size_t>::max()) const
    {
        return radius(q.data(), r, maxMatches);
    }

    const std::vector<std::array<double,3>>& points() const { return cloud_.pts; }

private:
    int     dim_;
    Adaptor cloud_;
    KDIndex index_;
    bool    built_{false};
};
