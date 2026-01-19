#pragma once

#include <cmath>
#include <cstddef>
#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include "Mesh.hpp"

#include "FemProjection.hpp"
#include "PointLocator.hpp"
#include "PointLocatorResult.hpp"
#include "Field.hpp"
#include "ElementTopology.hpp"
#include "KDTree.hpp"
#include "Mollifier.hpp"



class Reinit {
    public:
        using Vec = std::vector<double>;

        enum class Location : unsigned {
            Nodal     = 0,
            Elemental = 1
        };

        Reinit() = delete;         
        explicit Reinit(std::array<std::unique_ptr<FemProjection>, 6>& el_proj, std::vector<Field>& f, unsigned id, Mollifier m) noexcept : 
            _el_proj(el_proj),
            _field(f),
            _id(id),
            _mesh(_field.back().mesh()),
            _elTplgy(_mesh.elTplgy()),
            _elType(_mesh.elType()),
            _elLevel(_mesh.elLevel()),
            _X(_mesh.X()),
            _U(_field.back().getNodalById(_id)),
            _d(_mesh.dim()),
            _nEl(_mesh.numElements()),
            _nN(_mesh.numNodes()),
            _m(m)
            {
                _cut_region_nodes_cnt.resize(_nN);
            }

        ~Reinit() = default;

        //=================================================================================================
        // ===============================         MAIN FUNCTIONS          ===============================
        //=================================================================================================

        void computeInterfaceMarkers(std::vector<std::vector<double>>& global_markers, unsigned n_int) {

            _n_int = n_int;

            if (_d == 0) throw std::runtime_error("Reinit::computeInterfaceMarkers: mesh.dim()==0");

            global_markers.clear();
            global_markers.resize(_d);
            for (std::size_t a = 0; a < _d; ++a) {
                global_markers[a].clear();
            }

            for (std::size_t e = 0; e < _nEl; ++e) {

                std::vector<std::vector<double>> edge_roots;
                if (cellIsCut(e, edge_roots)) {

                    const auto& conn = _elTplgy[e];
                    for (unsigned i = 0; i < conn.size(); ++i) {
                        const unsigned gnode = conn[i];                      
                        _cut_region_nodes_cnt[gnode]++;
                    }

                    const auto local_markers = createLocalMarkers(edge_roots);

                    for (std::size_t a = 0; a < _d; a++)
                        global_markers[a].insert(global_markers[a].end(), local_markers[a].begin(), local_markers[a].end());

                }
            }
        }

        void reinitializeSignedDistance(const std::vector<std::vector<double>>& markers_soa) {

            collectCutRegionNodes();
            const auto converged_flag = projectPointsOnInterface(_cut_region_nodes);

            reinitFarField(markers_soa);
            reinitCutField(converged_flag);
        }

        //=================================================================================================
        // ============================         AUXILIARY FUNCTIONS          =============================
        //=================================================================================================

        bool cellIsCut(std::size_t e, std::vector<std::vector<double>>& edge_roots) {

            const auto& conn = _elTplgy[e];
            if (conn.empty()) return false;

            const auto& topo = ex35_topo::topologyFromTypeId(_elType[e]);

            const unsigned Nv = topo.nVertices;
            if (conn.size() < static_cast<std::size_t>(Nv + 1u)) {
                throw std::runtime_error("Reinit::cellIsCut: connectivity too small for element type");
            }

            // sign change on vertices only
            double umin = 0.0, umax = 0.0;
            bool first = true;
            for (unsigned i = 0; i < Nv; ++i) {
                const unsigned node = conn[i];
                if (node >= _nN) throw std::runtime_error("Reinit::cellIsCut: vertex node out of range");
                const double u = _U[node];
                if (first) {
                    umin = umax = u;
                    first = false;
                }
                else {
                    if (u < umin) umin = u;
                    if (u > umax) umax = u;
                }
            }

            if (!(umin <= 0.0 && umax >= 0.0)) 
                return false;
            else {
                if (conn.size() < static_cast<std::size_t>(topo.nNodes)) {
                    throw std::runtime_error("Reinit::cellIsCut: conn.size() < topo.nNodes (missing higher-order nodes?)");
                }

                edge_roots.clear();
                edge_roots.resize(_d);         
                for (std::size_t a = 0; a < _d; ++a) {
                    edge_roots[a].clear();   
                }

                _h_eff = 0.0;

                for (unsigned k = 0; k < topo.nEdges; ++k) {

                    const auto& E = topo.edges[k];

                    if (E.nodes.n != 3 || E.vertices.n != 2) 
                        throw std::runtime_error("Reinit::cellIsCut: expected quadratic edge (3 nodes, 2 vertices)");

                    unsigned lv0 = E.vertices[0];
                    unsigned lv1 = E.vertices[1];
                    unsigned lvm = E.nodes.p[2];

                    if (lv0 >= conn.size() || lv1 >= conn.size() || lvm >= conn.size())
                        throw std::runtime_error("Reinit::cellIsCut: local edge indices out of conn range");

                    unsigned g0 = conn[lv0];
                    unsigned g1 = conn[lv1];
                    unsigned gm = conn[lvm];

                    if (g0 >= _nN || g1 >= _nN || gm >= _nN)
                        throw std::runtime_error("Reinit::cellIsCut: global edge node out of range");

                    double h_eff_tmp = 0.0;
                    for (std::size_t a = 0; a < _d; ++a) {
                        const double dx = _X[a][g1] - _X[a][g0];
                        h_eff_tmp += dx * dx;
                    }
                    h_eff_tmp = std::sqrt(h_eff_tmp);
                    if (h_eff_tmp > _h_eff) _h_eff = h_eff_tmp;

                    double u0 = _U[g0];
                    double u1 = _U[g1];
                    double um = _U[gm];

                    const double a = 2.0 * (u0 + u1 - 2.0 * um);
                    const double b = 4.0 * um - 3.0 * u0 - u1;
                    const double c = u0;

                    double roots[2];
                    unsigned nRoots = 0;
                    const double epsT = 1e-12;

                    auto pushRootIfIn01 = [&](double t) {
                        if (t >= -epsT && t <= 1.0 + epsT) {
                            if (t < 0.0) t = 0.0;
                            if (t > 1.0) t = 1.0;
                            roots[nRoots++] = t;
                        }
                    };

                    if (std::abs(a) < 1e-16) {
                        if (std::abs(b) > 1e-16) {
                            const double t = -c / b;
                            pushRootIfIn01(t);
                        }
                    } else {
                        const double disc = b*b - 4.0*a*c;
                        if (disc >= -1e-16) {
                            const double sdisc = std::sqrt(std::max(0.0, disc));
                            const double inv2a = 1.0 / (2.0 * a);
                            pushRootIfIn01((-b - sdisc) * inv2a);
                            pushRootIfIn01((-b + sdisc) * inv2a);
                        }
                    }

                    if (nRoots == 0) continue;

                    for (unsigned iroot = 0; iroot < nRoots; iroot++) {

                        double t = roots[iroot];
                        const double N0 = 2.0 * (t - 0.5) * (t - 1.0);
                        const double Nm = 4.0 * t * (1.0 - t);
                        const double N1 = 2.0 * t * (t - 0.5);

                        const std::size_t ip = edge_roots[0].size();

                        for (std::size_t a = 0; a < _d; ++a) {
                            const double x0 = _X[a][g0];
                            const double xm = _X[a][gm];
                            const double x1 = _X[a][g1];
                            edge_roots[a].push_back(N0 * x0 + Nm * xm + N1 * x1);
                        }
                    }
                }
            }

            return true;
        }

        std::vector<std::vector<double>> createLocalMarkers(std::vector<std::vector<double>>& edge_roots) {
            std::vector<std::vector<double>> local_markers;

            if (_d == 1) {
                return edge_roots;
            }
            else if (_d == 2) {
                local_markers.clear();
                local_markers.resize(2);
                for (int a = 0; a < 2; a ++)
                    local_markers[a].clear();

                if (edge_roots[0].size() < 2) {
                    return local_markers;
                }
                if (edge_roots[0].size() == 2) {
                    double len = 0;
                    for (int a = 0; a < 2; a ++) 
                        len += (edge_roots[a][1] - edge_roots[a][0]) * (edge_roots[a][1] - edge_roots[a][0]);
                    len = std::sqrt(len);

                    int n_int_eff = static_cast<int>(std::ceil(len / _h_eff * _n_int )); 
                    n_int_eff = std::max(2, n_int_eff);

                    std::array<double,2> dist = {0,0};
                    for (int a = 0; a < 2; a ++)
                        dist[a] = (edge_roots[a][1] - edge_roots[a][0]) / static_cast<double>(n_int_eff);

                    for (int imarker = 0; imarker < n_int_eff + 1; imarker++) 
                        for (int a = 0; a < 2; a ++) 
                            local_markers[a].push_back(edge_roots[a][0] + dist[a]*imarker);
                }
                else {
                    std::runtime_error("Reinit::createLocalMarkers: multiple cuts not implemented yet");
                }
                
                return local_markers;
            }
            else if (_d == 3) {
                std::runtime_error("Reinit::createLocalMarkers: 3d not implemented yet");
            }
            else 
                std::runtime_error("Reinit::createLocalMarkers: unexpected dimension");
        }

        void collectCutRegionNodes() {
            _cut_region_nodes_idx.clear();
            _cut_region_nodes_idx.reserve(_cut_region_nodes_cnt.size());

            for (std::size_t gid = 0; gid < _cut_region_nodes_cnt.size(); ++gid) {
                if (_cut_region_nodes_cnt[gid] > 0) {
                    _cut_region_nodes_idx.push_back(gid);
                }
            }

            _cut_region_nodes.clear();
            _cut_region_nodes.resize(_d);

            for (int a = 0; a < _d; ++a) {
                _cut_region_nodes[a].clear();
                _cut_region_nodes[a].reserve(_cut_region_nodes_idx.size());

                for (std::size_t li = 0; li < _cut_region_nodes_idx.size(); ++li) {
                    const std::size_t gid = _cut_region_nodes_idx[li];
                    _cut_region_nodes[a].push_back(_X[a][gid]);
                }
            }
        }

        std::vector<uint8_t> projectPointsOnInterface(std::vector<std::vector<double>>& points) {
            PointLocator pl = PointLocator(_field[0].mesh(), .1);
            std::vector<PointLocatorResult> out, in;

            std::vector<double> p_values;       
            p_values.reserve(points[0].size());                
            std::vector<std::vector<double>> p_gradients(_d);
            for (int a = 0; a < _d; ++a)
                p_gradients[a].reserve(points[0].size());

            std::vector<uint8_t> converged_flag(points[0].size(), 0);

            const double eps_g2 = 1e-30;
            const double tol_phi = (_d <= 2) ? 1e-12 : 1e-8;
            const int max_iter = 30;

            std::vector<std::vector<double>> active_points(points);
            std::vector<std::size_t> active_idx(active_points[0].size());
            std::iota(active_idx.begin(), active_idx.end(), 0);

            std::size_t n_active = active_idx.size();
            int cnt = 0;
            unsigned it;

            for (it = 0; it < max_iter; ++it) {

                if (active_points.empty() || active_points[0].size() != active_idx.size())
                    throw std::runtime_error("Reinit: active_points/active_idx size mismatch");

                out.clear();
                pl.locateAll(out, active_points);
                for (unsigned l = 1; l <= _field.size()-1; l++) {
                    std::swap(in, out);
                    _field[l - 1].mesh().projectPointLocatorResultsToNextLevel(_field[l].mesh(), in, out);
                }
                if (out.size() != n_active) throw std::runtime_error("Reinit::reinitializeCutField: plr size != #marker points");

                _field.back().evalNodalAtLocatedPointsById(_id, out, _el_proj, p_values, p_gradients, /*outsideVal=*/0.0);

                if (p_values.size() != n_active) throw std::runtime_error("Reinit::reinitializeCutField: phi size mismatch");
                for (std::size_t a = 0; a < _d; ++a)
                    if (p_gradients[a].size() != n_active) throw std::runtime_error("Reinit::reinitializeCutField: grad size mismatch");

                std::vector<std::size_t> new_active_idx;
                std::vector<std::size_t> new_local_idx;
                new_active_idx.reserve(out.size());
                new_local_idx.reserve(out.size());

                for (std::size_t j = 0; j < out.size(); ++j) {
                    if (out[j].ok && std::abs(p_values[j]) > tol_phi) {
                        new_active_idx.push_back(active_idx[j]); 
                        new_local_idx.push_back(j);
                    }
                    else if (out[j].ok && std::abs(p_values[j]) <= tol_phi)
                        converged_flag[active_idx[j]] = 1;
                    else 
                        cnt++;
                }

                active_idx.swap(new_active_idx);
                n_active = active_idx.size();

                if (n_active == 0) break;

                std::vector<double> g2_values(n_active, 0.0);
                for (int a = 0; a < _d; ++a) {
                    for (std::size_t idx = 0; idx < n_active; ++idx) {
                        const std::size_t j = new_local_idx[idx];
                        g2_values[idx] += p_gradients[a][j] * p_gradients[a][j];
                    }
                }

                for (int a = 0; a < _d; ++a) {
                    for (std::size_t idx = 0; idx < n_active; ++idx) {
                        const std::size_t i = active_idx[idx];
                        const std::size_t j = new_local_idx[idx];
                        points[a][i] -= p_values[j] / (g2_values[idx] + eps_g2) * p_gradients[a][j];
                    }
                }

                for (int a = 0; a < _d; ++a) {
                    active_points[a].clear();
                    active_points[a].reserve(n_active);
                    for (std::size_t idx = 0; idx < n_active; ++idx) {
                        active_points[a].push_back(points[a][active_idx[idx]]);
                    }
                }
            }

            cnt += n_active;
            if (cnt > 0)
                std::cout<<"REINIT::projectPointsOnInterface : "<<cnt<<" points not converged."<<std::endl;

            return converged_flag;
        }

        void reinitCutField(std::vector<uint8_t> converged_flag) {
            std::vector<double> dist2(_cut_region_nodes_idx.size(), 0.0);
            for (std::size_t a = 0; a < _d; a++)
                for (std::size_t i = 0; i < _cut_region_nodes_idx.size(); i++) {
                    if (converged_flag[i] == 1) {
                        dist2[i] += (_cut_region_nodes[a][i] - _X[a][_cut_region_nodes_idx[i]]) * (_cut_region_nodes[a][i] - _X[a][_cut_region_nodes_idx[i]]);
                    }
                }
            for (std::size_t i = 0; i < _cut_region_nodes_idx.size(); i++) {
                if (converged_flag[i] == 1) {
                    double sign = (_U[_cut_region_nodes_idx[i]] > 0) ? 1.0 : -1.0;
                    double dist = sqrt(dist2[i]);
                    _U[_cut_region_nodes_idx[i]] = sign * _m.SigmoidC1(dist);
                }
            } 
        }

        void reinitFarField(const std::vector<std::vector<double>>& markers_soa) {
            const auto markers_aos = soaToAos(markers_soa);

            KDTreeRT kdtree(_d);
            kdtree.build(markers_aos);

            std::vector<std::size_t> nearest_marker(_nN, std::numeric_limits<std::size_t>::max());
            std::vector<double>      nearest_dist2(_nN, std::numeric_limits<double>::infinity());

            for (std::size_t i = 0; i < _nN; ++i)
            {
                std::array<double,3> p{0.0, 0.0, 0.0};
                for (std::size_t a = 0; a < static_cast<std::size_t>(_d); ++a)
                    p[a] = _X[a][i];

                auto knn_res = kdtree.knn(p, 1);
                const auto& idx   = knn_res.first;
                const auto& dist2 = knn_res.second;
                if (idx.empty()) continue; 

                const std::size_t j = idx[0];       
                nearest_marker[i] = j;
                nearest_dist2[i]  = dist2[0];

                const double dist = std::sqrt(dist2[0]);
                double sign = (_U[i] > 0) ? 1.0 : -1.0;

                _U[i] = sign * _m.SigmoidC1(dist);
            }
        }

        std::vector<std::array<double,3>> soaToAos(const std::vector<std::vector<double>>& points_soa) {
            std::size_t n_points = points_soa[0].size();
            std::vector<std::array<double,3>> points_aos(n_points);

            for (std::size_t a = 0; a < _d; a++) 
                for (std::size_t i = 0; i < n_points; i++)
                    points_aos[i][a] = points_soa[a][i];
            
            return points_aos;
        }

    private:
        std::array<std::unique_ptr<FemProjection>, 6>& _el_proj;
        std::vector<Field>& _field;
        const unsigned _id;

        const Mesh _mesh;
        const Mesh::ElemConnectivity& _elTplgy;
        const Mesh::ElemType& _elType;
        const Mesh::ElemLevel& _elLevel;
        const Mesh::Coordinates& _X;
        Vec& _U;
        const std::size_t _d;
        const std::size_t _nEl;
        const std::size_t _nN;
        const Mollifier _m;

        unsigned _n_int;
        double _h_eff;
        std::vector<int> _cut_region_nodes_cnt;
        std::vector<int> _cut_region_nodes_idx;
        std::vector<std::vector<double>> _cut_region_nodes;
        
};