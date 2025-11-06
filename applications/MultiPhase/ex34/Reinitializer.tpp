
namespace fem {

    template<std::size_t DIM>
    Reinitializer<DIM>::Reinitializer(OctTree<DIM>* tree_ptr, u32 fid_phi, std::function<double(double)> mollifier_, bool flag) 
        : tree(tree_ptr), fid(fid_phi), mollifier(std::move(mollifier_)), proj_flag(flag) {}

    template <std::size_t DIM>
    void Reinitializer<DIM>::write_markers_csv(const std::string& filename) const
    {
        std::ofstream out(filename);
        out.setf(std::ios::fixed, std::ios::floatfield);
        out << "x,y,z\n";
        for (const auto& p : markers) {
            out << std::setprecision(16) << p[0];
            for (unsigned int idim = 1; idim < DIM; idim++)
                out << "," << p[idim];
            if (DIM==2)
                out << "," << 0.;
            out << "\n";
        }
    }

    template <std::size_t DIM>
    bool Reinitializer<DIM>::leaf_is_cut(u32 leaf_pos, std::vector<Point<DIM>>* cell_intersections) 
    {
        // we check if the cell (leaf_id) is cut
        // if it is we compute its intersections with the interface

        // check for cut cells
        double mn = local_values[0], mx = local_values[0];
        for (int a = 1; a < local_values.size(); ++a) { 
            mn = std::min(mn, local_values[a]); 
            mx = std::max(mx, local_values[a]); 
        }

        const bool cut = (mn < 0.0) && (mx > 0.0);

        // if not cut clear and return false
        if (!cut) {
            if (cell_intersections) cell_intersections->clear();
            return false;
        }

        cell_intersections->clear();        
        
        if constexpr (DIM == 2)
        {
            cell_intersections->reserve(4); 

            static const int E[4][3]      = { {0,4,1}, {1,5,2}, {2,6,3}, {3,7,0} };
            static const int varDim[4]    = { 0,    1,   0,    1 };
            static const int varVersus[4] = { 1,    1,   -1,    -1 };
            static const double XiC[4]    = { 0.0,     +1.0,    0.0,     -1.0    };
            static const double EtaC[4]   = { -1.0,    0.0,     +1.0,    0.0     };

            std::map<int,std::vector<Point<DIM>>> edge_intersections_map;

            for (int e = 0; e < 4; ++e) 
            {
                const double vL = local_values[E[e][0]];
                const double vM = local_values[E[e][1]];
                const double vR = local_values[E[e][2]];
                auto roots = edge_roots(vL, vM, vR); 
                for (double s : roots) 
                {
                    double xi  = (varDim[e] == 0) ? s*varVersus[e] : XiC[e];
                    double eta = (varDim[e] == 1) ? s*varVersus[e] : EtaC[e];

                    Point<DIM> P = {xi, eta};
                    cell_intersections->push_back(P);
                }
            }

            if (cell_intersections->size() == 4) { // TODO CHECK
                const double epsb = 1e-12;

                auto edge_and_t = [&](const Point<2>& P){
                    const double xi = P[0], eta = P[1];
                    int e = -1; double t = 0.0;
                    if (std::abs(eta + 1.0) <= epsb)       { e = 0; t = 0.5*(xi  + 1.0); } 
                    else if (std::abs(xi  - 1.0) <= epsb)  { e = 1; t = 0.5*(eta + 1.0); } 
                    else if (std::abs(eta - 1.0) <= epsb)  { e = 2; t = 0.5*(-xi + 1.0);} 
                    else if (std::abs(xi  + 1.0) <= epsb)  { e = 3; t = 0.5*(-eta+ 1.0);} 
                    return std::pair<int,double>(e,t);
                };

                struct Node { Point<2> p; int e; double t; };
                std::array<Node,4> Q;
                for (int i=0; i<4; ++i) {
                    auto [e,t] = edge_and_t((*cell_intersections)[i]);
                    Q[i] = { (*cell_intersections)[i], e, t };
                }

                std::sort(Q.begin(), Q.end(), [](const Node& a, const Node& b){
                    if (a.e != b.e) return a.e < b.e;
                    return a.t < b.t;
                });

                const double f00 = local_values[0];
                const double f10 = local_values[1]; 
                const double f11 = local_values[2]; 
                const double f01 = local_values[3]; 
                const double S   = f00*f11 - f10*f01;

                std::array<Point<2>,4> ordered;

                if (S > 0.0) {
                    ordered[0] = Q[0].p; ordered[1] = Q[1].p;
                    ordered[2] = Q[2].p; ordered[3] = Q[3].p;
                } else if (S < 0.0) {
                    ordered[0] = Q[1].p; ordered[1] = Q[2].p;
                    ordered[2] = Q[3].p; ordered[3] = Q[0].p;
                } else {
                    const double fcc = local_values[8];
                    const bool like0011 = (fcc*f00 >= 0.0) && (fcc*f11 >= 0.0);
                    if (like0011) {
                        ordered[0] = Q[0].p; ordered[1] = Q[1].p;
                        ordered[2] = Q[2].p; ordered[3] = Q[3].p;
                    } else {
                        ordered[0] = Q[1].p; ordered[1] = Q[2].p;
                        ordered[2] = Q[3].p; ordered[3] = Q[0].p;
                    }
                }

                (*cell_intersections)[0] = ordered[0];
                (*cell_intersections)[1] = ordered[1];
                (*cell_intersections)[2] = ordered[2];
                (*cell_intersections)[3] = ordered[3];
            }
        }
        else 
        {
            static const int E[12][3] = {
                {0,  8, 1}, {1,  9, 2}, {2, 10, 3}, {3, 11, 0}, // bottom
                {0, 16, 4}, {1, 17, 5}, {2, 18, 6}, {3, 19, 7}, // verticals
                {4, 12, 5}, {5, 13, 6}, {6, 14, 7}, {7, 15, 4}  // top
            };

            static const int varDim[12]  = {
                0, 1, 0, 1,
                2, 2, 2, 2,
                0, 1, 0, 1
            };

            static const int varVersus[12]  = {
                +1, +1,  -1,  -1,
                +1, +1,  +1,  +1,
                +1, +1,  -1,  -1
            };

            static const double XiC[12]   = {  0, +1,  0, -1,  -1, +1, +1, -1,   0, +1,  0, -1 };
            static const double EtaC[12]  = { -1,  0, +1,  0,  -1, -1, +1, +1,  -1,  0, +1,  0 };
            static const double ZetaC[12] = { -1, -1, -1, -1,   0,  0,  0,  0,  +1, +1, +1, +1 };

            cell_intersections->reserve(5);
            for (int e = 0; e < 12; ++e) {
                const double vL = local_values[E[e][0]];
                const double vM = local_values[E[e][1]];
                const double vR = local_values[E[e][2]];
                auto roots = edge_roots(vL, vM, vR); 
                for (double s : roots) {
                    double xi   = (varDim[e] == 0) ? s*varVersus[e] : XiC[e];
                    double eta  = (varDim[e] == 1) ? s*varVersus[e] : EtaC[e];
                    double zeta = (varDim[e] == 2) ? s*varVersus[e] : ZetaC[e];

                    Point<DIM> P = {xi, eta, zeta};
                    cell_intersections->push_back(P);
                }
            }

            // TODO gestire doppie intersezioni
        }


        return true;
    }

    template <std::size_t DIM>
    void Reinitializer<DIM>::compute_cell_markers(std::vector<Point<DIM>> cell_intersections)
    {
        const int cell_density = marker_density * pow(2 , tree->max_depth() - tree->_tree_nodes[leaf_id].level);
        const int cell_min_segments = min_segments * pow(2 , tree->max_depth() - tree->_tree_nodes[leaf_id].level);
        std::vector<Point<DIM>> cell_markers;
        if constexpr (DIM==2)
        {
            if (cell_intersections.size() == 2) {
                const auto& A = cell_intersections[0];
                const auto& B = cell_intersections[1];

                int n_segments = static_cast<int>(std::ceil(Geom_op<DIM>::norm(Geom_op<DIM>::sub(B,A)) * cell_density));
                n_segments = std::max(cell_min_segments, n_segments);

                Vector<DIM> step = Geom_op<DIM>::div( Geom_op<DIM>::sub(B,A) , static_cast<double>(n_segments));

                for (int k = 0; k < n_segments; ++k) 
                    cell_markers.push_back(Geom_op<DIM>::add( A , Geom_op<DIM>::mul( static_cast<double>(k) , step) ) );
                cell_markers.push_back(B);
            }
            else if (cell_intersections.size() == 4) {
                for (int i = 0; i < 2; i++) {
                    const auto& A = cell_intersections[0+i*2];
                    const auto& B = cell_intersections[1+i*2];

                    int n_segments = static_cast<int>(std::ceil(Geom_op<DIM>::norm(Geom_op<DIM>::sub(B,A)) * cell_density));
                    n_segments = std::max(cell_min_segments, n_segments);

                    Vector<DIM> step = Geom_op<DIM>::div( Geom_op<DIM>::sub(B,A) , static_cast<double>(n_segments));

                    for (int k = 0; k < n_segments; ++k) 
                        cell_markers.push_back(Geom_op<DIM>::add( A , Geom_op<DIM>::mul( static_cast<double>(k) , step) ) );
                    cell_markers.push_back(B);
                }
            }
            

        }
        else if constexpr (DIM==3)
        {
            // we create a marker grid per cell along this "plane"
            
            assert(cell_intersections.size() <= 5 && "More than 5 intersections in cut cell");

            // we create a polygon with the cell intersections
            auto poly = Tri_ord::order_polygon_3d(cell_intersections);
            if (poly.size() < 3) return;

            // 1 / 2 / 3 triangles are created from the polygon 
            auto tris = Tri_ord::triangulate_poly(poly);
                
            // finally we place markers on each triangles
            for (const auto& tr : tris)
                Tri_ord::sample_triangle_with_density(tr, cell_density, cell_min_segments, cell_markers);

        }
        
        // and finally each marker is projected on the interface using leaf coordinates
        double tol = 1e-10;
        const int max_iters = 30;   

        std::vector<size_t> unmatched(cell_markers.size());
        std::iota(unmatched.begin(), unmatched.end(), 0);

        for (int it = 0; it < max_iters && !unmatched.empty(); ++it) 
        {
            std::vector<size_t> next;
            next.reserve(unmatched.size());

            for (size_t idx : unmatched) 
            {
                auto& p = cell_markers[idx]; 

                double value = evaluate_field_on_leaf(p);

                if (std::abs(value) <= tol) continue; 

                Vector<DIM> grad =  evaluate_gradient_on_leaf(p);

                p = Geom_op<DIM>::sub(p , Geom_op<DIM>::div(Geom_op<DIM>::mul(value , grad) , Geom_op<DIM>::norm2(grad)));

                next.push_back(idx);
            }

            unmatched.swap(next);
        }

        std::vector<Point<DIM>> accepted_markers;
        for (unsigned int i = 0; i < cell_markers.size(); i++){
            double value = evaluate_field_on_leaf(cell_markers[i]);
            if (std::abs(value) <= tol){
                accepted_markers.push_back(cell_markers[i]);
            }
        }
        cell_markers.swap(accepted_markers);
        
        for (unsigned int i = 0; i < cell_markers.size(); i++)
        {
            cell_markers[i] = evaluate_coord_on_leaf(cell_markers[i]);
        }

        markers.insert(markers.end(), cell_markers.begin(), cell_markers.end());
    }

    template <std::size_t DIM>
    Point<DIM> Reinitializer<DIM>::evaluate_coord_on_leaf(Point<DIM> P_local)
    {
        // we compute physical coordinates from leaf ones

        Point<DIM> P_global;
        for (unsigned int idim = 0; idim < DIM; idim++)
            P_global[idim] = 0;

        constexpr int Nn = (DIM == 2 ? 9 : 27);
        std::array<double, Nn> N{};   

        if constexpr (DIM == 2) {
            Shapes2D::Q9(P_local, N.data());
        } 
        else 
        {
            Shapes3D::H27(P_local, N.data());
        }
             
        for (unsigned int i = 0; i<local_coord.size(); i++)
            for (int idim = 0; idim < DIM; idim++)
            {
                P_global[idim] += local_coord[i][idim] * N[i];
            }
        
        return P_global;
    }

    template <std::size_t DIM>
    double Reinitializer<DIM>::evaluate_field_on_leaf(const Point<DIM> P)
    {
        constexpr int Nn = (DIM == 2 ? 9 : 27);
        std::array<double, Nn> N{};   

        if constexpr (DIM == 2) {
            Shapes2D::Q9(P, N.data());
        } 
        else 
        {
            Shapes3D::H27(P, N.data());
        }
        
        double value = 0;
        for (unsigned int i = 0; i<local_values.size(); i++)
        {
            value += local_values[i] * N[i];
        }

        return value;
    }

    template <std::size_t DIM>
    Vector<DIM> Reinitializer<DIM>::evaluate_gradient_on_leaf(const Point<DIM> P)
    {
        constexpr int Nn = (DIM == 2 ? 9 : 27);
        std::array<std::array<double, Nn> , DIM> dN{};   

        if constexpr (DIM == 2) {
            Shapes2D::Q9_dN(P, dN[0].data(), dN[1].data());
        } 
        else 
        {
            Shapes3D::H27_dN(P, dN[0].data(), dN[1].data(), dN[2].data());
        }
        
        Vector<DIM> grad;
        for (int idim = 0; idim < DIM; idim++)
            grad[idim] = 0;

        for (unsigned int i = 0; i<local_values.size(); i++)
            for (int idim = 0; idim < DIM; idim++)
                grad[idim] += local_values[i] * dN[idim][i];
        
        return grad;
    }


    template <std::size_t DIM>
    std::vector<double> Reinitializer<DIM>::edge_roots(double v0, double v1, double v2)
    {
        // we look for roots in the edge's quadratic polynomial
        // v0 v1 and v2 are consecutive in xi / eta

        const double A = 0.5*v0 - v1 + 0.5*v2;
        const double B = -0.5*v0 + 0.5*v2;
        const double C = v1;
        std::vector<double> r;
        const double eps = 1e-14;
        if (std::abs(A) < eps) 
        {
            if (std::abs(B) > eps) 
            {
                double s = -C/B;
                if (s >= -1.0-1e-12 && s <= 1.0+1e-12) r.push_back(std::max(-1.0, std::min(1.0, s)));
            }
        } 
        else 
        {
            double D = B*B - 4.0*A*C;
            if (D >= 0.0) 
            {
                double rt = std::sqrt(D);
                double s1 = (-B - rt) / (2.0*A);
                double s2 = (-B + rt) / (2.0*A);
                if (s1 >= -1.0-1e-12 && s1 <= 1.0+1e-12) r.push_back(std::max(-1.0, std::min(1.0, s1)));
                if (s2 >= -1.0-1e-12 && s2 <= 1.0+1e-12) r.push_back(std::max(-1.0, std::min(1.0, s2)));
                if (r.size()==2 && std::abs(r[0]-r[1])<1e-12) r.pop_back();
            }
        }
        return r;
    }

    

    template <std::size_t DIM>
    std::vector<Point<DIM>> Reinitializer<DIM>::collect_markers(int density, int min_segments_)
    {
        marker_density = density;
        min_segments = min_segments_;

        const auto& L = tree->leaves();
        const size_t numCells = L.size();

        auto& Field      = tree->field(fid);
        auto  Basis      = to_basis<DIM>(Field.basis_id);
        const auto& nodes = tree->basis_nodes(Basis);
        BasisRegistry<DIM> basis_registry = tree->_basisReg[(int)Basis];

        markers.clear();
        markers.reserve(L.size()*marker_density); 

        cut_cells.clear();
        cut_cells.reserve(L.size());

        for (int k = 0; k < static_cast<int>(numCells); k++) 
        {
            leaf_id = L[k];

            const u32 leaf_pos = (leaf_id < tree->_tree_nodes.size()) ? tree->_tree_nodes[leaf_id].node2leafIdx : npos32;

            const auto &conn = basis_registry.elem2glob[leaf_pos];
            
            local_coord.assign(conn.size(), Point<DIM>{}); 
            local_values.assign(conn.size(), 0.0);

            for (u32 a = 0; a < conn.size(); ++a) {
                u32 gid          = conn[a];                  
                local_coord[a]   = nodes[gid].physical;       
                local_values[a]  = Field.nodal[gid];
            }

            std::vector<Point<DIM>> cell_intersections;    

            if (leaf_is_cut(leaf_pos, &cell_intersections)) 
            {
                cut_cells.push_back(leaf_pos);

                compute_cell_markers(cell_intersections);
            }
            else if (tree->_tree_nodes[leaf_id].level == tree->max_depth()) {
                cut_cells.push_back(leaf_pos);
            }
        }

        return markers;
    }

    template <std::size_t DIM>
    void Reinitializer<DIM>::project_cut_cells_nodes()
    {
        // we project each cutcell node on the interface and calculate the distance 
        // the projection is done in parent coordinates

        const auto field = tree->field(fid);
        const auto basis = to_basis<DIM>(field.basis_id);
        
        const double tol      = 1e-10;
        const int    max_iter = 30;

        cut_cells_nodes_dist.clear();
        cut_cells_nodes_dist.resize(cut_cells.size());

        for (size_t i_cell = 0; i_cell < cut_cells.size(); ++i_cell)
        {
            std::vector<Point<DIM>> parent_coords;
            tree->extract_leaf_parent_coords(basis, cut_cells[i_cell], parent_coords);
            Point<DIM> P0 = tree->parent_to_physical(parent_coords[0]);
            Point<DIM> P1 = tree->parent_to_physical(parent_coords[1]);

            double h = Geom_op<DIM>::norm(Geom_op<DIM>::sub(P0,P1));

            cut_cells_nodes_dist[i_cell].resize(parent_coords.size());

            const double lam0    = 1e-12;  
            const double eps_g2  = 1e-12;  
            const int    max_bt  = 8;       
            const double smax    = h;  

            for (size_t i_node = 0; i_node < parent_coords.size(); ++i_node)
            {
                bool converged = false;

                Point<DIM> p = parent_coords[i_node];
                p = Geom_op<DIM>::add(p , Geom_op<DIM>::mul(1.e-4 , Geom_op<DIM>::sub(parent_coords[parent_coords.size()-1] , p)));
                Vector<DIM> grad_parent;
            
                for (int it = 0; it < max_iter; ++it)
                {
                    double value = 0.0;
                    if (!tree->evaluate_field_on_parent(fid, p, value)) {
                        break;
                    }
                    if (std::abs(value) <= tol) {
                        converged = true;
                        break;
                    }

                    for (int d = 0; d < DIM; ++d) grad_parent[d] = 0.0;
                    if (!tree->evaluate_gradient_on_parent(fid, p, grad_parent)) {
                        break;
                    }

                    const double g2   = Geom_op<DIM>::norm2(grad_parent);
                    const double denom = g2 + lam0;

                    Vector<DIM> step = Geom_op<DIM>::mul(value / denom, grad_parent);

                    double step_norm = Geom_op<DIM>::norm(step);
                    if (step_norm > smax && step_norm > 0.0) {
                        step = Geom_op<DIM>::mul(smax / step_norm, step);
                        step_norm = smax;
                    }

                    bool accepted = false;
                    Point<DIM> p_try;
                    double value_new = std::numeric_limits<double>::infinity();

                    double alpha = 1.0;
                    const double f0 = std::abs(value);

                    for (int bt = 0; bt < max_bt; ++bt) 
                    {
                        p_try = Geom_op<DIM>::sub(p, Geom_op<DIM>::mul(alpha, step));

                        if (tree->evaluate_field_on_parent(fid, p_try, value_new)) 
                        {
                            if ((std::abs(value_new) < f0 && value*value_new >= 0) || std::abs(value_new) < f0*0.5) 
                            { 
                                accepted = true;
                                break;
                            }
                        }

                        alpha *= 0.5; 
                    }

                    if (accepted) 
                    {
                        p = p_try;

                        if (std::abs(value_new) <= tol) {
                            converged = true;
                            break;
                        }
                    }
                    else 
                        break; // TODO add something
                }

                Point<DIM> p_physical = tree->parent_to_physical(p);
                Point<DIM> node_physical = tree->parent_to_physical(parent_coords[i_node]);

                cut_cells_nodes_dist[i_cell][i_node] = Geom_op<DIM>::norm( Geom_op<DIM>::sub( p_physical , node_physical )); 

            }
            
        }

    }

    template <std::size_t DIM>
    void Reinitializer<DIM>::compute_signed_distance()
    {

        if (proj_flag)
            project_cut_cells_nodes();

        auto& Field      = tree->field(fid);
        auto  Basis      = to_basis<DIM>(Field.basis_id);
        const auto& nodes = tree->basis_nodes(Basis);
        BasisRegistry<DIM> basis_registry = tree->_basisReg[(int)Basis];

        KDTree<DIM> kdtree;
        kdtree.build_from_any(markers);

        // for each node in the grid we look for the closest marker
        for (size_t i = 0; i < nodes.size(); ++i) 
        {
            Point<DIM> p_physical = nodes[i].physical;
            auto knn = kdtree.knn(p_physical, 1);
            double dist = 0.0;
            if (!knn.first.empty())
                dist = std::sqrt(knn.second[0]);
            else 
            {
                std::cout<<"knn search failed"<<std::endl;
                abort();
            }

            const double sign = (Field.nodal[i] > 0) ? 1. : -1.;
            Field.nodal[i] = sign * dist;
        }

        // for each cutcell node we substitute its projection
        // in this way we reduce tangential errors on the interface
        if (proj_flag)
            for (unsigned int i_cell = 0; i_cell < cut_cells.size(); i_cell++)
            {
                const u32 leaf_pos = cut_cells[i_cell];
                const auto &conn = basis_registry.elem2glob[leaf_pos];

                for (unsigned int i_node = 0; i_node < cut_cells_nodes_dist[i_cell].size(); i_node++)
                {   
                    const double sign = (Field.nodal[conn[i_node]] > 0) ? 1. : -1.;
                    const double new_dist = cut_cells_nodes_dist[i_cell][i_node];
                    const double old_dist = Field.nodal[conn[i_node]];
                    if (new_dist < fabs(old_dist))
                    {
                        Field.nodal[conn[i_node]] = sign * new_dist;
                    }
                }
            }

        if (mollifier)
            for (size_t i = 0; i < nodes.size(); ++i) {
                Field.nodal[i] = mollifier(Field.nodal[i]);
            }
    }

    
}
