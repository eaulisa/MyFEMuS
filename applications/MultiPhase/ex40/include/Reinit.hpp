#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "MultiLevelSolution.hpp"

#include "BBoxToIel.hpp"
#include "FemProjection.hpp"
#include "KDTree.hpp"
#include "LevelMarkers.hpp"
#include "Mollifier.hpp"
#include "MyVector.hpp"

using namespace femus;

class Reinit {
public:
  Reinit(std::string name, MultiLevelSolution &mlSol, double eps)
      : _name(name), _mlSol(mlSol), _m(eps) {
    _el_proj[0].reset(new Hex27Projection());
    _el_proj[1].reset(new Tet15Projection());
    _el_proj[2].reset(new Wedge21Projection());
    _el_proj[3].reset(new Quad9Projection());
    _el_proj[4].reset(new Tri7Projection());
    _el_proj[5].reset(new Line3Projection());

    MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();

    _mlSol1.Build(&mlMsh);
    _mlSol1.AddSolution("Psi", LAGRANGE, SECOND);
    _mlSol1.Initialize("All");
  }

  struct DofPair {
    unsigned meshDof;
    unsigned solDof;
  };

  void farFieldReinit(const std::vector<MyVector<double>> &markers) {

    int my_rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    MultiLevelMesh &mlMsh = *_mlSol.GetMultilevelMesh();
    const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
    Mesh &msh = *mlMsh.GetLevel(level);
    const unsigned dim = msh.GetDimension();
    Solution &solR = *_mlSol.GetLevel(level);
    Solution &solW = *_mlSol1.GetLevel(level);
    const unsigned mesh_proc = msh.processor_id();
    if (mesh_proc != my_rank)
      throw std::runtime_error("Mesh processor ID does not match MPI rank");

    const unsigned solIndex = _mlSol.GetIndex(_name.c_str());
    const unsigned solType = _mlSol.GetSolutionType(_name.c_str());
    const unsigned xType = 2u;
    int n_local = markers[0].end() - markers[0].begin();

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // compute query points
    auto &xv = msh._topology->_Sol;
    auto &solVecR = solR._Sol[solIndex];
    auto &solVecW = solW._Sol[solIndex];
    const unsigned offset = msh._elementOffset[my_rank];
    const unsigned offsetp1 = msh._elementOffset[my_rank + 1];
    const unsigned solOffset = msh._dofOffset[solType][my_rank];
    const unsigned solOffsetp1 = msh._dofOffset[solType][my_rank + 1];

    std::vector<DofPair> dof_pairs;

    for (unsigned iel = offset; iel < offsetp1; ++iel) {
      const unsigned nDof = msh.GetElementDofNumber(iel, solType);
      for (unsigned i = 0; i < nDof; ++i) {
        if (msh.GetSolutionDof(i, iel, solType) >= solOffset &&
            msh.GetSolutionDof(i, iel, solType) < solOffsetp1)
          dof_pairs.push_back({msh.GetSolutionDof(i, iel, xType),
                               msh.GetSolutionDof(i, iel, solType)});
      }
    }
    std::sort(dof_pairs.begin(), dof_pairs.end(),
              [](const DofPair &a, const DofPair &b) {
                return a.meshDof < b.meshDof;
              });
    dof_pairs.erase(std::unique(dof_pairs.begin(), dof_pairs.end(),
                                [](const DofPair &a, const DofPair &b) {
                                  return a.meshDof == b.meshDof;
                                }),
                    dof_pairs.end());

    const unsigned nq = dof_pairs.size();
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // compute bounding sphere of local points

    // baricenter of the local subdomain -> center of the sphere
    std::vector<double> c(dim, 0.);
    double r = 0.0;

    if (!dof_pairs.empty()) {
      // find farest point B from first dof A
      auto first_dof = dof_pairs.front().meshDof;

      double max_d2 = -1.0;
      size_t far_dof_B = first_dof;

      for (const auto idx : dof_pairs) {
        double d2 = 0.0;
        for (int idim = 0; idim < dim; idim++) {
          double diff = (*xv[idim])(idx.meshDof) - (*xv[idim])(first_dof);
          d2 += diff * diff;
        }
        if (d2 > max_d2) {
          max_d2 = d2;
          far_dof_B = idx.meshDof;
        }
      }

      // find farest point C from B
      max_d2 = -1.0;
      size_t far_dof_C = far_dof_B;
      for (const auto idx : dof_pairs) {
        double d2 = 0.0;
        for (int idim = 0; idim < dim; idim++) {
          double diff = (*xv[idim])(idx.meshDof) - (*xv[idim])(far_dof_B);
          d2 += diff * diff;
        }
        if (d2 > max_d2) {
          max_d2 = d2;
          far_dof_C = idx.meshDof;
        }
      }

      // initialize center as (B+C)/2 and radius as ||B-C||/2
      r = std::sqrt(max_d2) / 2.0;
      for (int idim = 0; idim < dim; idim++) {
        c[idim] = ((*xv[idim])(far_dof_B) + (*xv[idim])(far_dof_C)) * 0.5;
      }

      // correct sphere to contain all points
      for (const auto idx : dof_pairs) {
        double d2 = 0.0;
        for (int idim = 0; idim < dim; idim++) {
          double dist = (*xv[idim])(idx.meshDof) - c[idim];
          d2 += dist * dist;
        }

        if (d2 > r * r) {
          double d = std::sqrt(d2);
          double r_new = (r + d) * 0.5;
          double alpha = (r_new - r) / d;

          for (int idim = 0; idim < dim; idim++) {
            double dist = (*xv[idim])(idx.meshDof) - c[idim];
            c[idim] += dist * alpha;
          }
          r = r_new;
        }
      }
    }
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // gather relevant markers for this rank

    // Allgather bounding spheres on all ranks
    const int sphere_stride = dim + 1;
    std::vector<double> sphere_sendbuf(sphere_stride);
    for (int idim = 0; idim < dim; idim++) {
      sphere_sendbuf[idim] = c[idim];
    }
    sphere_sendbuf[dim] = r;

    std::vector<double> sphere_recvbuf(nprocs * sphere_stride);

    MPI_Allgather(sphere_sendbuf.data(), sphere_stride, MPI_DOUBLE,
                  sphere_recvbuf.data(), sphere_stride, MPI_DOUBLE,
                  MPI_COMM_WORLD);

    // Accessor lambda
    auto center = [&](int p, int d) -> double {
      return sphere_recvbuf[p * sphere_stride + d];
    };
    auto radius = [&](int p) -> double {
      return sphere_recvbuf[p * sphere_stride + dim];
    };

    // compute min (max distance) S of local markers from each sphere
    const unsigned b = markers[0].begin();
    std::vector<double> local_min_dist2(
        nprocs, std::numeric_limits<double>::infinity());

    for (int iproc = 0; iproc < nprocs; ++iproc) {
      double r_p = radius(iproc);
      for (int imarker = 0; imarker < n_local; ++imarker) {
        double d2 = 0.0;
        for (int idim = 0; idim < dim; idim++) {
          double dist = markers[idim][b + imarker] - center(iproc, idim);
          d2 += dist * dist;
        }
        local_min_dist2[iproc] = std::min(local_min_dist2[iproc], d2);
      }
    }

    std::vector<double> global_min_dist2(nprocs);
    MPI_Allreduce(local_min_dist2.data(), global_min_dist2.data(), nprocs,
                  MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    std::vector<double> S(nprocs);
    for (int p = 0; p < nprocs; ++p)
      S[p] = std::sqrt(global_min_dist2[p]) + radius(p);

    // select markers for each ranks
    //
    // we will send markers to the interested rank using MPI_Alltoallv
    //
    // int MPI_Alltoallv(const void *sendbuf, const int sendcounts[], const int
    // sdispls[], MPI_Datatype sendtype,
    //                   void *recvbuf, const int recvcounts[], const int
    //                   rdispls[], MPI_Datatype recvtype, MPI_Comm
    //                   MPI_COMM_WORLD)
    //
    // - sendbuf contains the data to send
    // - sendcounts[iproc] contain the number of elements to send from my_rank
    // to iproc
    // - sdispls[iproc] contain the displacement in sendbuf from which to take
    // the data to send to iproc
    //    : sendbuf[sdispls[iproc]:sdispls[iproc]+sendcounts[iproc]]
    // - recvbuf is the buffer to receive data from other ranks
    // - recvcounts[iproc] contain the number of elements that my_rank will
    // receive from iproc
    // - rdispls[iproc] contain the displacement in recvbuf at which to place
    // the incoming data from iproc
    //    : recvbuf[rdispls[iproc]:rdispls[iproc]+recvcounts[iproc]]
    //
    // we first compute sendcounts and recvcounts
    std::vector<int> sendcounts(nprocs, 0);

    for (int imarker = 0; imarker < n_local; ++imarker) {
      for (int iproc = 0; iproc < nprocs; iproc++) {
        double d2 = 0.0;
        for (int idim = 0; idim < dim; idim++) {
          double dist = markers[idim][b + imarker] - center(iproc, idim);
          d2 += dist * dist;
        }
        if (sqrt(d2) - radius(iproc) <= S[iproc]) {
          sendcounts[iproc] += dim;
        }
      }
    }

    // sendcounts[iproc] contain the number of double that my_rank will send to
    // iproc recvcounts[iproc] contain the number of double that my_rank will
    // receive from iproc MPI_Alltoall works as a transpose of the sendcounts
    // matrix
    std::vector<int> recvcounts(nprocs);
    MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT,
                 MPI_COMM_WORLD);

    // now we compute displacements for writing and riding the send and receive
    // buffers
    std::vector<int> sdispls(nprocs, 0), rdispls(nprocs, 0);
    for (int p = 1; p < nprocs; ++p) {
      sdispls[p] = sdispls[p - 1] + sendcounts[p - 1];
      rdispls[p] = rdispls[p - 1] + recvcounts[p - 1];
    }
    const int total_send = sdispls[nprocs - 1] + sendcounts[nprocs - 1];
    const int total_recv = rdispls[nprocs - 1] + recvcounts[nprocs - 1];

    // now we fill the send buffer
    std::vector<double> sendbuf(total_send);
    std::vector<int> fill_offset = sdispls;

    for (int imarker = 0; imarker < n_local; ++imarker) {
      for (int iproc = 0; iproc < nprocs; iproc++) {
        double d2 = 0.0;
        for (int idim = 0; idim < dim; idim++) {
          double dist = markers[idim][b + imarker] - center(iproc, idim);
          d2 += dist * dist;
        }
        if (sqrt(d2) - radius(iproc) <= S[iproc]) {
          for (int idim = 0; idim < dim; idim++) {
            sendbuf[fill_offset[iproc]++] = markers[idim][b + imarker];
          }
        }
      }
    }

    // now we can perform the alltoallv
    std::vector<double> recvbuf(total_recv);
    MPI_Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_DOUBLE,
                  recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_DOUBLE,
                  MPI_COMM_WORLD);

    // finally we reconstruct the received markers
    const int total_markers_recv = total_recv / dim;

    std::vector<std::vector<double>> local_markers(
        dim, std::vector<double>(total_markers_recv));
    for (int im = 0; im < total_markers_recv; ++im)
      for (int d = 0; d < dim; ++d)
        local_markers[d][im] = recvbuf[im * dim + d];
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Perform KDTree search
    const auto markers_aos = soaToAos(local_markers);

    KDTreeRT kdtree(dim);
    kdtree.build(markers_aos);

    std::vector<std::size_t> nearest_marker(
        nq, std::numeric_limits<std::size_t>::max());
    std::vector<double> nearest_dist2(nq,
                                      std::numeric_limits<double>::infinity());

    for (std::size_t i = 0; i < nq; ++i) {
      std::array<double, 3> p{0.0, 0.0, 0.0};
      for (std::size_t a = 0; a < static_cast<std::size_t>(dim); ++a)
        p[a] = (*xv[a])(dof_pairs[i].meshDof);

      auto knn_res = kdtree.knn(p, 1);
      const auto &idx = knn_res.first;
      const auto &dist2 = knn_res.second;

      if (idx.empty())
        throw std::runtime_error("empty search");

      const std::size_t j = idx[0];
      nearest_marker[i] = j;
      nearest_dist2[i] = dist2[0];

      const double dist = std::sqrt(dist2[0]);
      double sign = ((*solVecR)(dof_pairs[i].solDof) > 0) ? 1.0 : -1.0;

      solVecW->set(dof_pairs[i].solDof, dist * sign);
    }
    solVecW->close();
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  }

  void interfaceFieldReinit(BBoxToIel &bbox) {
    int my_rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    MultiLevelMesh &mlMsh = *_mlSol.GetMultilevelMesh();
    const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
    Mesh &msh = *mlMsh.GetLevel(level);
    auto &xv = msh._topology->_Sol;
    const unsigned dim = msh.GetDimension();
    Solution &solR = *_mlSol.GetLevel(level);
    Solution &solW = *_mlSol1.GetLevel(level);
    const unsigned mesh_proc = msh.processor_id();
    if (mesh_proc != my_rank)
      throw std::runtime_error("Mesh processor ID does not match MPI rank");

    const unsigned solIndex = _mlSol.GetIndex(_name.c_str());
    const unsigned solType = _mlSol.GetSolutionType(_name.c_str());
    const unsigned xType = 2u;

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Extract nodes at max level
    std::vector<MyVector<double>> X1;
    std::vector<unsigned> maxLevelDofOffset;
    maxLevelDofOffset.clear();
    maxLevelDofOffset.resize(nprocs);
    const auto dof_pairs =
        GetMaxLevelSolutionPoints(_mlSol, _name, X1, maxLevelDofOffset);
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // this set associates the ith active index to the original local position
    std::vector<unsigned> activeIds(maxLevelDofOffset[my_rank + 1] -
                                    maxLevelDofOffset[my_rank]);
    std::iota(activeIds.begin(), activeIds.end(), 0);

    // here we save projection info
    std::vector<std::vector<double>> converged_points;
    std::vector<bool> converged_flag(
        maxLevelDofOffset[my_rank + 1] - maxLevelDofOffset[my_rank], false);
    converged_points.resize(dim);
    for (int idim = 0; idim < dim; idim++)
      converged_points[idim].resize(maxLevelDofOffset[my_rank + 1] -
                                    maxLevelDofOffset[my_rank]);

    const unsigned nLevels = mlMsh.GetNumberOfLevels();
    const unsigned bboxLevels = nLevels - bbox.GetLevel();
    const unsigned level0 = nLevels - 1u;

    Mesh &msh0 = *mlMsh.GetLevel(level0);
    Solution &sol0 = *_mlSol.GetLevel(level0);

    const unsigned solIndex0 = _mlSol.GetIndex(_name.c_str());
    const unsigned solType0 = _mlSol.GetSolutionType(_name.c_str());

    auto &solVecR = solR._Sol[solIndex0];
    auto &solVecW = solW._Sol[solIndex0];

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Projection loop
    int it = 0;
    int max_it = 10;
    double tol = (dim == 2) ? 1.e-8 : 1.e-6;
    while (it++ < max_it) {

      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      // Inverse mapping of the X1 ste of active points
      LevelMarkers l0;
      std::vector<LevelMarkers> lX(bboxLevels);

      bbox.GetInverseMappingOnCoarseLevel(X1, l0, lX[0]);
      for (unsigned k = 1; k < bboxLevels; ++k) {
        bbox.Project(mlMsh, lX[k - 1], lX[k]);
      }
      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      LevelMarkers &lTop = lX.back();

      std::vector<MyVector<double>> &Xi = lTop.GetLocalCoordinates();
      MyVector<unsigned> &Iel = lTop.GetElements();

      std::vector<double> psiLocal;
      std::vector<std::vector<double>> pointLocal;
      psiLocal.resize(Iel.end() - Iel.begin(), 0.0);
      pointLocal.resize(dim);
      for (int idim = 0; idim < dim; idim++)
        pointLocal[idim].resize(Iel.end() - Iel.begin(), 0.0);

      std::vector<double> xi(dim);
      std::vector<double> phi;
      std::vector<double> dphi;

      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      // Compute the iterative projection for each active point
      for (unsigned ip = Iel.begin(); ip < Iel.end(); ++ip) {

        const unsigned iel = Iel[ip];
        short unsigned ielType = msh0.GetElementType(iel);

        for (unsigned k = 0; k < dim; ++k) {
          xi[k] = Xi[k][ip];
        }

        const unsigned nDof = msh0.GetElementDofNumber(iel, solType0);

        std::vector<std::vector<double>> iel_coord;
        iel_coord.resize(dim);
        for (int idim = 0; idim < dim; idim++) {
          iel_coord[idim].resize(nDof);
          for (unsigned j = 0; j < nDof; ++j) {
            const unsigned xdof = msh.GetSolutionDof(j, iel, xType);
            iel_coord[idim][j] = (*xv[idim])(xdof);
          }
        }

        dphi.resize(nDof * dim, 0.0);
        phi.resize(nDof, 0.0);

        double weight;
        _el_proj[ielType]->fem().Jacobian(iel_coord, xi, weight, phi, dphi);

        double value = 0.0;
        std::vector<double> grad(dim, 0.0);
        std::vector<double> point(dim, 0.0);
        for (unsigned j = 0; j < nDof; ++j) {
          const unsigned solDof = msh.GetSolutionDof(j, iel, solType);
          const unsigned xdof = msh.GetSolutionDof(j, iel, xType);
          value += phi[j] * (*solVecR)(solDof);
          for (int idim = 0; idim < dim; idim++) {
            grad[idim] += dphi[j * dim + idim] * (*solVecR)(solDof);
            point[idim] += phi[j] * (*xv[idim])(xdof);
          }
        }

        double grad_norm = 0;
        for (int idim = 0; idim < dim; idim++)
          grad_norm += grad[idim] * grad[idim];

        psiLocal[ip - Iel.begin()] = value;
        for (int idim = 0; idim < dim; idim++)
          pointLocal[idim][ip - Iel.begin()] =
              point[idim] - value * grad[idim] / grad_norm;
      }
      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      lTop.GetFields().resize(dim + 1);
      lTop.GetFields()[0].buildFromLocal(psiLocal);
      for (int idim = 0; idim < dim; idim++)
        lTop.GetFields()[idim + 1].buildFromLocal(pointLocal[idim]);

      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      // Project backward through marker hierarchy
      std::vector<std::vector<double>> Wfield_r;
      std::vector<std::vector<double>> Wfield_s;

      const unsigned nFields = dim + 1;
      const bool backward = true;

      for (int l = static_cast<int>(bboxLevels) - 1; l >= 1; --l) {
        lX[l].RebuildLocalFromField(Wfield_r, nFields, backward);
        lX[l - 1].SendLocalField(Wfield_r, Wfield_s);
        lX[l - 1].RebuildFieldFromLocal(Wfield_s, nFields, backward);
      }

      lX[0].RebuildLocalFromField(Wfield_r, nFields, backward);
      l0.SendLocalField(Wfield_r, Wfield_s);
      l0.RebuildFieldFromLocal(Wfield_s, nFields, backward);
      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      const MyVector<double> &psiProjected = l0.GetFields()[0];
      const std::vector<MyVector<double>> &pointsProjected =
          (dim == 2)
              ? std::vector<MyVector<double>>(
                    {l0.GetFields()[1], l0.GetFields()[2]})
              : std::vector<MyVector<double>>(
                    {l0.GetFields()[1], l0.GetFields()[2], l0.GetFields()[3]});

      const std::vector<bool> &isInsideDomain = l0.GetPointInsideDomain();

      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      // Check for convergence and update active points
      std::vector<unsigned> oldIds;
      std::vector<unsigned> activeIds_tmp;
      unsigned offset = psiProjected.begin();
      for (unsigned i = psiProjected.begin(); i < psiProjected.end(); ++i) {
        const double value = psiProjected[i];
        if (isInsideDomain[i - offset]) {
          if (fabs(value) > tol) {
            oldIds.push_back(i);
            activeIds_tmp.push_back(activeIds[i - psiProjected.begin()]);
          } else {
            converged_flag[activeIds[i - psiProjected.begin()]] = true;
            for (int idim = 0; idim < dim; idim++)
              converged_points[idim][activeIds[i - psiProjected.begin()]] =
                  pointsProjected[idim][i];
          }
        }
      }
      std::swap(activeIds, activeIds_tmp);

      int nLocal = oldIds.size();
      int nGlobal = 0.0;

      MPI_Allreduce(&nLocal, &nGlobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      std::cout << "[Reinit] projection step: " << it
                << " - Active dofs: " << nGlobal << std::endl;

      if (nGlobal == 0) {
        std::cout << "[Reinit] projection converged in " << it << " iterations"
                  << std::endl;
        break;
      }

      std::vector<std::vector<double>> Xloc_tmp(dim);
      for (unsigned k = 0; k < dim; k++)
        Xloc_tmp[k].reserve(nLocal);

      for (const unsigned i : oldIds) {
        for (unsigned k = 0; k < dim; k++)
          Xloc_tmp[k].push_back(pointsProjected[k][i]);
      }

      std::vector<MyVector<double>> X_tmp;
      X_tmp.resize(dim);
      for (unsigned k = 0; k < dim; ++k)
        X_tmp[k].buildFromLocal(Xloc_tmp[k]);

      std::swap(X_tmp, X1);
      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    }
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Set converged values
    for (int i = 0; i < converged_flag.size(); i++) {
      if (converged_flag[i]) {
        unsigned xdof = dof_pairs[i].first;
        unsigned soldof = dof_pairs[i].second;

        double sign = ((*solVecR)(soldof) > 0.) ? 1. : -1.;

        double dist = 0.;
        for (int idim = 0; idim < dim; idim++)
          dist += ((*xv[idim])(xdof)-converged_points[idim][i]) *
                  ((*xv[idim])(xdof)-converged_points[idim][i]);

        dist = sqrt(dist);

        if (dist < fabs((*solVecW)(soldof)))
          solVecW->set(soldof, dist * sign);
      }
    }
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  }

  std::vector<std::pair<unsigned, unsigned>>
  GetMaxLevelSolutionPoints(MultiLevelSolution &mlSol, const std::string &name,
                            std::vector<MyVector<double>> &X,
                            std::vector<unsigned> &maxLevelDofOffset) {
    MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
    const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
    Mesh &msh = *mlMsh.GetLevel(level);
    const unsigned dim = msh.GetDimension();
    const unsigned solType = mlSol.GetSolutionType(name.c_str());
    const unsigned xType = 2u;
    if (solType > xType)
      throw std::runtime_error(
          "coordinate FE space is too low-order for solType");

    const unsigned iproc = msh.processor_id();
    const unsigned nprocs = msh.n_processors();
    auto &xv = msh._topology->_Sol;

    const unsigned solOffset = msh._dofOffset[solType][iproc];
    const unsigned solOffsetp1 = msh._dofOffset[solType][iproc + 1];
    const unsigned elOffset = msh._elementOffset[iproc];
    const unsigned elOffsetp1 = msh._elementOffset[iproc + 1];

    std::map<unsigned, unsigned> xdof_to_sdof;
    for (unsigned iel = elOffset; iel < elOffsetp1; ++iel) {
      if (msh.el->GetElementLevel(iel) != level)
        continue;
      const unsigned nDof = msh.GetElementDofNumber(iel, solType);
      for (unsigned i = 0; i < nDof; ++i) {
        const unsigned sdof = msh.GetSolutionDof(i, iel, solType);
        const unsigned xdof = msh.GetSolutionDof(i, iel, xType);
        if (solOffset <= sdof && sdof < solOffsetp1)
          xdof_to_sdof[xdof] = sdof;
      }
    }

    std::vector<std::pair<unsigned, unsigned>> pairs(xdof_to_sdof.begin(),
                                                     xdof_to_sdof.end());

    const unsigned nLocal = pairs.size();

    maxLevelDofOffset.resize(nprocs + 1, 0);
    MPI_Allgather(&nLocal, 1, MPI_UNSIGNED, maxLevelDofOffset.data() + 1, 1,
                  MPI_UNSIGNED, MPI_COMM_WORLD);
    for (unsigned p = 1; p <= nprocs; ++p)
      maxLevelDofOffset[p] += maxLevelDofOffset[p - 1];

    std::vector<std::vector<double>> Xloc(dim);
    for (unsigned k = 0; k < dim; k++)
      Xloc[k].resize(nLocal);

    for (unsigned iq = 0; iq < nLocal; ++iq) {
      const unsigned xdof = pairs[iq].first;
      for (unsigned k = 0; k < dim; k++)
        Xloc[k][iq] = (*xv[k])(xdof);
    }

    X.resize(dim);
    for (unsigned k = 0; k < dim; ++k)
      X[k].buildFromLocal(Xloc[k]);

    return pairs;
  }

  std::vector<std::array<double, 3>>
  soaToAos(const std::vector<std::vector<double>> &points_soa) {
    std::size_t n_points = points_soa[0].size();
    std::vector<std::array<double, 3>> points_aos(n_points);

    for (std::size_t a = 0; a < points_soa.size(); a++)
      for (std::size_t i = 0; i < n_points; i++)
        points_aos[i][a] = points_soa[a][i];

    return points_aos;
  }

  void updateSolution() {

    int my_rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    MultiLevelMesh &mlMsh = *_mlSol.GetMultilevelMesh();
    const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
    Mesh &msh = *mlMsh.GetLevel(level);
    auto &xv = msh._topology->_Sol;
    const unsigned dim = msh.GetDimension();
    Solution &solR = *_mlSol.GetLevel(level);
    Solution &solW = *_mlSol1.GetLevel(level);
    const unsigned mesh_proc = msh.processor_id();
    if (mesh_proc != my_rank)
      throw std::runtime_error("Mesh processor ID does not match MPI rank");

    const unsigned solIndex = _mlSol.GetIndex(_name.c_str());
    const unsigned solType = _mlSol.GetSolutionType(_name.c_str());

    auto &solVecR = solR._Sol[solIndex];
    auto &solVecW = solW._Sol[solIndex];

    const unsigned solOffset = msh._dofOffset[solType][my_rank];
    const unsigned solOffsetp1 = msh._dofOffset[solType][my_rank + 1];

    for (unsigned i = solOffset; i < solOffsetp1; i++) {
      solVecR->set(i, _m.SigmoidC1((*solVecW)(i)));
    }

    solVecR->close();
  }

private:
  const std::string _name;
  MultiLevelSolution &_mlSol;
  Mollifier _m;

  std::array<std::unique_ptr<FemProjection>, 6> _el_proj;

  MultiLevelSolution _mlSol1;
};