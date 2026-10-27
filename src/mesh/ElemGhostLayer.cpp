#include <algorithm>
#include <array>
#include <climits>
#include <iostream>
#include <map>
#include <mpi.h>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Elem.hpp"
#include "refinement/MeshRefinement.hpp"

namespace femus {

  static void ExchangeNodeBuffers(const std::vector<std::vector<unsigned>>& buffers_send, std::vector<std::vector<unsigned>>& buffers_recv, const unsigned iproc, const unsigned nprocs, const int tag) {
    std::vector<int> counts_send(nprocs, 0);
    std::vector<int> counts_recv(nprocs, 0);

    for (unsigned jproc = 0; jproc < nprocs; jproc++) {
      if (jproc == iproc) continue;
      if (buffers_send[jproc].size() > INT_MAX) MPI_Abort(MPI_COMM_WORLD, 1);
      counts_send[jproc] = static_cast<int>(buffers_send[jproc].size());
    }

    MPI_Alltoall(counts_send.data(), 1, MPI_INT, counts_recv.data(), 1, MPI_INT, MPI_COMM_WORLD);

    buffers_recv.clear();
    buffers_recv.resize(nprocs);
    buffers_recv[iproc] = buffers_send[iproc];

    std::vector<MPI_Request> requests(2 * nprocs, MPI_REQUEST_NULL);

    for (unsigned jproc = 0; jproc < nprocs; jproc++) {
      if (jproc == iproc) continue;

      if (counts_recv[jproc] > 0) {
        buffers_recv[jproc].resize(counts_recv[jproc]);
        MPI_Irecv(buffers_recv[jproc].data(), counts_recv[jproc], MPI_UNSIGNED, static_cast<int>(jproc), tag, MPI_COMM_WORLD, &requests[jproc]);
      }

      if (counts_send[jproc] > 0) {
        MPI_Isend(buffers_send[jproc].data(), counts_send[jproc], MPI_UNSIGNED, static_cast<int>(jproc), tag, MPI_COMM_WORLD, &requests[nprocs + jproc]);
      }
    }

    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
  }

  std::unordered_map<unsigned, std::vector<unsigned>>& elem::GetGhostElementLayer(const unsigned jproc) {
    return _ghostElementLayer[jproc];
  }

  const std::unordered_map<unsigned, std::vector<unsigned>>& elem::GetGhostElementLayer(const unsigned jproc) const {
    return _ghostElementLayer[jproc];
  }

  unsigned elem::GetGhostElementType(const unsigned jproc, const unsigned iel) const {
    return _ghostElementType[jproc].at(iel);
  }

  unsigned elem::GetGhostElementLevel(const unsigned jproc, const unsigned iel) const {
    return _ghostElementLevel[jproc].at(iel);
  }

  bool elem::GetIfGhostElementFatherHasBeenRefined(const unsigned jproc, const unsigned iel) const {
    return _ghostElementLevel[jproc].at(iel) == _level;
  }

  void elem::BuildElementNearElement(Mesh *msh, const bool updateElementNearFace) {
    const unsigned dim = msh->GetDimension();
    const unsigned nprocs = msh->n_processors();
    const unsigned iproc = msh->processor_id();

    _ghostElementLayer.clear();
    _ghostElementLayer.resize(nprocs);
    _ghostElementType.clear();
    _ghostElementType.resize(nprocs);
    _ghostElementLevel.clear();
    _ghostElementLevel.resize(nprocs);

    std::vector<std::vector<unsigned>> ghostElements_send(nprocs);

    {
      MyVector<unsigned> rowSize(_elementOffset, 1);

      std::unordered_map<unsigned, std::vector<unsigned>> nodeToElement;
      nodeToElement.reserve(rowSize.end() - rowSize.begin());

      std::vector<std::vector<unsigned>> ghostNodes_send(nprocs);

      for (unsigned iel = rowSize.begin(); iel < rowSize.end(); iel++) {
        for (unsigned i = 0; i < GetElementDofNumber(iel, 0); i++) {
          unsigned inode = GetElementDofIndex(iel, i);
          std::vector<unsigned>& elements = nodeToElement[inode];

          if (elements.empty()) {
            elements.reserve(2u << dim);
            unsigned jproc = msh->IsdomBisectionSearch(inode, 2);
            if (jproc != iproc) ghostNodes_send[jproc].push_back(inode);
          }

          elements.push_back(iel);
        }
      }

      std::vector<std::vector<unsigned>> ghost_send(nprocs);
      std::vector<std::vector<unsigned>> ghost_recv;

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        unsigned count_send = 0;

        for (unsigned i = 0; i < ghostNodes_send[jproc].size(); i++) {
          unsigned inode_send = ghostNodes_send[jproc][i];
          const std::vector<unsigned>& elements = nodeToElement.at(inode_send);
          if (count_send > static_cast<unsigned>(INT_MAX - 2) || elements.size() > static_cast<unsigned>(INT_MAX - 2) - count_send) MPI_Abort(MPI_COMM_WORLD, 1);
          count_send += 2 + static_cast<unsigned>(elements.size());
        }

        ghost_send[jproc].reserve(count_send);

        for (unsigned i = 0; i < ghostNodes_send[jproc].size(); i++) {
          unsigned inode_send = ghostNodes_send[jproc][i];
          const std::vector<unsigned>& elements = nodeToElement.at(inode_send);

          ghost_send[jproc].push_back(inode_send);
          ghost_send[jproc].push_back(static_cast<unsigned>(elements.size()));
          ghost_send[jproc].insert(ghost_send[jproc].end(), elements.begin(), elements.end());
        }
      }

      ExchangeNodeBuffers(ghost_send, ghost_recv, iproc, nprocs, 100);

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        unsigned pos_recv = 0;

        while (pos_recv < ghost_recv[jproc].size()) {
          unsigned inode_recv = ghost_recv[jproc][pos_recv++];
          unsigned count_recv = ghost_recv[jproc][pos_recv++];
          std::vector<unsigned>& elements = nodeToElement[inode_recv];

          elements.insert(elements.end(), ghost_recv[jproc].begin() + pos_recv, ghost_recv[jproc].begin() + pos_recv + count_recv);
          pos_recv += count_recv;
        }
      }

      for (std::unordered_map<unsigned, std::vector<unsigned>>::iterator it = nodeToElement.begin(); it != nodeToElement.end(); ++it) {
        std::vector<unsigned>& elements = it->second;
        std::sort(elements.begin(), elements.end());
        elements.erase(std::unique(elements.begin(), elements.end()), elements.end());
      }

      std::vector<std::vector<unsigned>> updated_send(nprocs);
      std::vector<std::vector<unsigned>> updated_recv;

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        unsigned pos_recv = 0;
        unsigned count_send = 0;

        while (pos_recv < ghost_recv[jproc].size()) {
          unsigned inode_recv = ghost_recv[jproc][pos_recv++];
          unsigned count_recv = ghost_recv[jproc][pos_recv++];
          pos_recv += count_recv;

          const std::vector<unsigned>& elements = nodeToElement.at(inode_recv);
          if (count_send > static_cast<unsigned>(INT_MAX - 2) || elements.size() > static_cast<unsigned>(INT_MAX - 2) - count_send) MPI_Abort(MPI_COMM_WORLD, 1);
          count_send += 2 + static_cast<unsigned>(elements.size());
        }

        updated_send[jproc].reserve(count_send);
        pos_recv = 0;

        while (pos_recv < ghost_recv[jproc].size()) {
          unsigned inode_recv = ghost_recv[jproc][pos_recv++];
          unsigned count_recv = ghost_recv[jproc][pos_recv++];
          pos_recv += count_recv;

          const std::vector<unsigned>& elements = nodeToElement.at(inode_recv);

          updated_send[jproc].push_back(inode_recv);
          updated_send[jproc].push_back(static_cast<unsigned>(elements.size()));
          updated_send[jproc].insert(updated_send[jproc].end(), elements.begin(), elements.end());
        }
      }

      ExchangeNodeBuffers(updated_send, updated_recv, iproc, nprocs, 101);

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        unsigned pos_recv = 0;

        while (pos_recv < updated_recv[jproc].size()) {
          unsigned inode_recv = updated_recv[jproc][pos_recv++];
          unsigned count_recv = updated_recv[jproc][pos_recv++];
          std::vector<unsigned>& elements = nodeToElement.at(inode_recv);

          elements.assign(updated_recv[jproc].begin() + pos_recv, updated_recv[jproc].begin() + pos_recv + count_recv);
          pos_recv += count_recv;
        }
      }

      std::vector<std::vector<unsigned>> elementNeighbors(rowSize.end() - rowSize.begin());

      for (unsigned iel = rowSize.begin(); iel < rowSize.end(); iel++) {
        std::vector<unsigned>& elements = elementNeighbors[iel - rowSize.begin()];
        elements.push_back(iel);

        for (unsigned i = 0; i < GetElementDofNumber(iel, 0); i++) {
          unsigned inode = GetElementDofIndex(iel, i);
          const std::vector<unsigned>& adjacentElements = nodeToElement.at(inode);

          for (unsigned j = 0; j < adjacentElements.size(); j++) {
            if (adjacentElements[j] != iel) elements.push_back(adjacentElements[j]);
          }
        }

        std::sort(elements.begin() + 1, elements.end());
        elements.erase(std::unique(elements.begin() + 1, elements.end()), elements.end());
        rowSize[iel] = static_cast<unsigned>(elements.size());

        for (unsigned j = 1; j < elements.size(); j++) {
          unsigned jproc = static_cast<unsigned>(std::upper_bound(_elementOffset.begin(), _elementOffset.end(), elements[j]) - _elementOffset.begin() - 1);
          if (jproc != iproc && (ghostElements_send[jproc].empty() || ghostElements_send[jproc].back() != iel)) ghostElements_send[jproc].push_back(iel);
        }
      }

      _elementNearElement.clear();
      _elementNearElement = MyMatrix<unsigned>(rowSize, UINT_MAX);

      for (unsigned iel = _elementNearElement.begin(); iel < _elementNearElement.end(); iel++) {
        const std::vector<unsigned>& elements = elementNeighbors[iel - rowSize.begin()];

        for (unsigned j = 0; j < elements.size(); j++) {
          _elementNearElement[iel][j] = elements[j];
        }
      }
      if(updateElementNearFace) {
        for (unsigned iel = rowSize.begin(); iel < rowSize.end(); iel++) {
          for (unsigned iface = 0; iface < GetElementFaceNumber(iel); iface++) {
            if (GetFaceElementIndex(iel, iface) != -1) continue;

            const std::vector<unsigned>* faceElements[3] = {nullptr, nullptr, nullptr};
            unsigned first = 0;

            for (unsigned k = 0; k < dim; k++) {
              unsigned inode = GetFaceVertexIndex(iel, iface, k);
              faceElements[k] = &nodeToElement.at(inode);
              if (faceElements[k]->size() < faceElements[first]->size()) first = k;
            }

            const std::vector<unsigned>& candidates = *faceElements[first];

            for (unsigned j = 0; j < candidates.size(); j++) {
              unsigned jel = candidates[j];
              if (jel == iel) continue;

              bool shared = true;

              for (unsigned k = 0; k < dim && shared; k++) {
                if (k != first) shared = std::binary_search(faceElements[k]->begin(), faceElements[k]->end(), jel);
              }

              if (shared) {
                SetFaceElementIndex(iel, iface, jel + 1u);
                break;
              }
            }
            if (_level == 0 && GetFaceElementIndex(iel, iface) == -1) SetFaceElementIndex(iel, iface, -101);
          }
        }
      }
    }

    std::vector<std::vector<unsigned>> connectivity_send(nprocs);
    std::vector<std::vector<unsigned>> connectivity_recv;

    for (unsigned jproc = 0; jproc < nprocs; jproc++) {
      unsigned count_send = 0;

      for (unsigned i = 0; i < ghostElements_send[jproc].size(); i++) {
        unsigned nDofs_send = GetElementDofNumber(ghostElements_send[jproc][i], 2);
        if (count_send > static_cast<unsigned>(INT_MAX - 4) || nDofs_send > static_cast<unsigned>(INT_MAX - 4) - count_send) MPI_Abort(MPI_COMM_WORLD, 1);
        count_send += 4 + nDofs_send;
      }

      connectivity_send[jproc].reserve(count_send);

      for (unsigned i = 0; i < ghostElements_send[jproc].size(); i++) {
        unsigned iel_send = ghostElements_send[jproc][i];
        unsigned nDofs_send = GetElementDofNumber(iel_send, 2);

        connectivity_send[jproc].push_back(iel_send);
        connectivity_send[jproc].push_back(GetElementType(iel_send));
        connectivity_send[jproc].push_back(GetElementLevel(iel_send));
        connectivity_send[jproc].push_back(nDofs_send);

        for (unsigned j = 0; j < nDofs_send; j++) connectivity_send[jproc].push_back(GetElementDofIndex(iel_send, j));
      }
    }

    ExchangeNodeBuffers(connectivity_send, connectivity_recv, iproc, nprocs, 102);

    for (unsigned jproc = 0; jproc < nprocs; jproc++) {
      unsigned pos_recv = 0;

      while (pos_recv < connectivity_recv[jproc].size()) {
        unsigned iel_recv = connectivity_recv[jproc][pos_recv++];
        unsigned ielType_recv = connectivity_recv[jproc][pos_recv++];
        unsigned ielLevel_recv = connectivity_recv[jproc][pos_recv++];
        unsigned nDofs_recv = connectivity_recv[jproc][pos_recv++];

        _ghostElementType[jproc][iel_recv] = ielType_recv;
        _ghostElementLevel[jproc][iel_recv] = ielLevel_recv;
        std::vector<unsigned>& dofs_recv = _ghostElementLayer[jproc][iel_recv];
        dofs_recv.resize(GetNVE(ielType_recv, 2), UINT_MAX);
        std::copy(connectivity_recv[jproc].begin() + pos_recv, connectivity_recv[jproc].begin() + pos_recv + nDofs_recv, dofs_recv.begin());
        pos_recv += nDofs_recv;
      }
    }

    _ghostElements_send.swap(ghostElements_send);
  }

  // ********************************************

  void elem::BuildEdgePoints(Mesh *msh) {

    unsigned nnodes = GetNodeNumber();
    const unsigned iproc = _iproc;
    const unsigned nprocs = _nprocs;

    if (_ghostElements_send.size() != nprocs || _ghostElementLayer.size() != nprocs || _ghostElementType.size() != nprocs || _ghostElementLevel.size() != nprocs) {
      std::cerr << "BuildEdgePoints: rank=" << iproc << " ghost layers have not been built" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    {
      unsigned counts_send[4] = {0, 0, 0, nnodes};
      unsigned pointOffset[4] = {0, 0, 0, 0};
      std::vector<unsigned> pointNodes_send;
      std::vector<std::vector<unsigned>> pointDofs_send(nprocs);

      for (unsigned kind = 0; kind < 2; kind++) {
        std::unordered_map<unsigned, std::unordered_map<unsigned, unsigned>> edgeIndex;
        std::map<std::array<unsigned, 5>, unsigned> faceIndex;
        std::vector<unsigned> pointElement;
        std::vector<unsigned> pointNode;
        std::vector<std::vector<unsigned>> pointDofs;

        for (unsigned iel = _elementOffset[iproc]; iel < _elementOffset[iproc + 1]; iel++) {
          unsigned ielt = GetElementType(iel);
          if (ielt >= 6) MPI_Abort(MPI_COMM_WORLD, 1);

          if (kind == 0) {
            for (unsigned inode = 0; inode < GetNVE(ielt, 2); inode++) {
              unsigned node = GetElementDofIndex(iel, inode);
              if ((inode < GetNVE(ielt, 0) && node == UINT_MAX) || (node != UINT_MAX && node >= nnodes)) {
                std::cerr << "BuildEdgePoints: rank=" << iproc << " invalid input element=" << iel << " slot=" << inode << " node=" << node << " nnodes=" << nnodes << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
              }
            }
          }

          if (!GetIfFatherHasBeenRefined(iel) || (kind == 1 && ielt >= 3)) continue;
          unsigned istart = GetNVE(ielt, kind);
          unsigned iend = kind == 0 ? GetNVE(ielt, 1) : istart + GetNFC(ielt, 1);

          for (unsigned inode = istart; inode < iend; inode++) {
            unsigned localPoint = inode - istart;
            unsigned nvertices = kind == 0 ? 2 : GetNFACENODES(ielt, localPoint, 0);
            std::array<unsigned, 5> vertices = {{nvertices, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX}};

            for (unsigned i = 0; i < nvertices; i++) {
              unsigned localVertex = kind == 0 ? edge2VerticesMapping[ielt][localPoint][i] : GetIG(ielt, localPoint, i);
              vertices[i + 1] = GetElementDofIndex(iel, localVertex);
            }

            std::sort(vertices.begin() + 1, vertices.begin() + 1 + nvertices);
            unsigned index = UINT_MAX;

            if (kind == 0) {
              std::unordered_map<unsigned, unsigned>& edges = edgeIndex[vertices[1]];
              std::unordered_map<unsigned, unsigned>::const_iterator it = edges.find(vertices[2]);
              if (it != edges.end()) index = it->second;
            }
            else {
              std::map<std::array<unsigned, 5>, unsigned>::const_iterator it = faceIndex.find(vertices);
              if (it != faceIndex.end()) index = it->second;
            }

            if (index == UINT_MAX) {
              if (pointElement.size() >= UINT_MAX) MPI_Abort(MPI_COMM_WORLD, 1);
              index = static_cast<unsigned>(pointElement.size());
              if (kind == 0) edgeIndex[vertices[1]][vertices[2]] = index;
              else faceIndex[vertices] = index;
              pointElement.push_back(iel);
              pointNode.push_back(UINT_MAX);
              pointDofs.emplace_back();
            }

            unsigned node = GetElementDofIndex(iel, inode);

            if (node != UINT_MAX) {
              if (pointNode[index] != UINT_MAX && pointNode[index] != node) {
                std::cerr << "BuildEdgePoints: rank=" << iproc << " conflicting shared point element=" << iel << " slot=" << inode << " node=" << node << " previous=" << pointNode[index] << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
              }
              pointNode[index] = node;
            }
            else {
              pointDofs[index].push_back(iproc);
              pointDofs[index].push_back(iel);
              pointDofs[index].push_back(inode);
            }
          }
        }

        for (unsigned jproc = 0; jproc < nprocs; jproc++) {
          for (std::unordered_map<unsigned, std::vector<unsigned>>::const_iterator it = _ghostElementLayer[jproc].begin(); it != _ghostElementLayer[jproc].end(); ++it) {
            unsigned iel = it->first;
            unsigned ielt = GetGhostElementType(jproc, iel);
            if (ielt >= 6) MPI_Abort(MPI_COMM_WORLD, 1);
            const std::vector<unsigned>& dofs = it->second;

            if (kind == 0) {
              if (dofs.size() != GetNVE(ielt, 2)) {
                std::cerr << "BuildEdgePoints: rank=" << iproc << " invalid ghost row owner=" << jproc << " element=" << iel << " size=" << dofs.size() << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
              }

              for (unsigned inode = 0; inode < dofs.size(); inode++) {
                unsigned node = dofs[inode];
                if ((inode < GetNVE(ielt, 0) && node == UINT_MAX) || (node != UINT_MAX && node >= nnodes)) {
                  std::cerr << "BuildEdgePoints: rank=" << iproc << " invalid ghost input owner=" << jproc << " element=" << iel << " slot=" << inode << " node=" << node << " nnodes=" << nnodes << std::endl;
                  MPI_Abort(MPI_COMM_WORLD, 1);
                }
              }
            }

            if (!GetIfGhostElementFatherHasBeenRefined(jproc, iel) || (kind == 1 && ielt >= 3)) continue;
            unsigned istart = GetNVE(ielt, kind);
            unsigned iend = kind == 0 ? GetNVE(ielt, 1) : istart + GetNFC(ielt, 1);

            for (unsigned inode = istart; inode < iend; inode++) {
              unsigned localPoint = inode - istart;
              unsigned nvertices = kind == 0 ? 2 : GetNFACENODES(ielt, localPoint, 0);
              std::array<unsigned, 5> vertices = {{nvertices, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX}};

              for (unsigned i = 0; i < nvertices; i++) {
                unsigned localVertex = kind == 0 ? edge2VerticesMapping[ielt][localPoint][i] : GetIG(ielt, localPoint, i);
                vertices[i + 1] = dofs[localVertex];
              }

              std::sort(vertices.begin() + 1, vertices.begin() + 1 + nvertices);
              unsigned index = UINT_MAX;

              if (kind == 0) {
                std::unordered_map<unsigned, std::unordered_map<unsigned, unsigned>>::const_iterator first = edgeIndex.find(vertices[1]);
                if (first == edgeIndex.end()) continue;
                std::unordered_map<unsigned, unsigned>::const_iterator second = first->second.find(vertices[2]);
                if (second != first->second.end()) index = second->second;
              }
              else {
                std::map<std::array<unsigned, 5>, unsigned>::const_iterator found = faceIndex.find(vertices);
                if (found != faceIndex.end()) index = found->second;
              }

              if (index == UINT_MAX) continue;
              pointElement[index] = std::min(pointElement[index], iel);
              unsigned node = dofs[inode];

              if (node != UINT_MAX) {
                if (pointNode[index] != UINT_MAX && pointNode[index] != node) {
                  std::cerr << "BuildEdgePoints: rank=" << iproc << " conflicting ghost point owner=" << jproc << " element=" << iel << " slot=" << inode << " node=" << node << " previous=" << pointNode[index] << std::endl;
                  MPI_Abort(MPI_COMM_WORLD, 1);
                }
                pointNode[index] = node;
              }
              else {
                pointDofs[index].push_back(jproc);
                pointDofs[index].push_back(iel);
                pointDofs[index].push_back(inode);
              }
            }
          }
        }

        for (unsigned i = 0; i < pointElement.size(); i++) {
          if (pointElement[i] < _elementOffset[iproc] || pointElement[i] >= _elementOffset[iproc + 1] || pointDofs[i].empty()) continue;
          if (pointNodes_send.size() >= UINT_MAX) MPI_Abort(MPI_COMM_WORLD, 1);
          unsigned index = static_cast<unsigned>(pointNodes_send.size());
          pointNodes_send.push_back(pointNode[i]);
          if (pointNode[i] == UINT_MAX) counts_send[kind]++;

          for (unsigned j = 0; j < pointDofs[i].size(); j += 3) {
            unsigned jproc = pointDofs[i][j];
            if (pointDofs_send[jproc].size() > static_cast<unsigned>(INT_MAX - 3)) MPI_Abort(MPI_COMM_WORLD, 1);
            pointDofs_send[jproc].push_back(pointDofs[i][j + 1]);
            pointDofs_send[jproc].push_back(pointDofs[i][j + 2]);
            pointDofs_send[jproc].push_back(index);
          }
        }

        pointOffset[kind + 1] = static_cast<unsigned>(pointNodes_send.size());
      }

      for (unsigned iel = _elementOffset[iproc]; iel < _elementOffset[iproc + 1]; iel++) {
        unsigned ielt = GetElementType(iel);
        if (!GetIfFatherHasBeenRefined(iel) || ielt >= 5) continue;
        unsigned inode = GetNVE(ielt, 2) - 1;
        if (GetElementDofIndex(iel, inode) != UINT_MAX) continue;
        if (pointNodes_send.size() >= UINT_MAX || pointDofs_send[iproc].size() > static_cast<unsigned>(INT_MAX - 3)) MPI_Abort(MPI_COMM_WORLD, 1);
        unsigned index = static_cast<unsigned>(pointNodes_send.size());
        pointNodes_send.push_back(UINT_MAX);
        counts_send[2]++;
        pointDofs_send[iproc].push_back(iel);
        pointDofs_send[iproc].push_back(inode);
        pointDofs_send[iproc].push_back(index);
      }

      pointOffset[3] = static_cast<unsigned>(pointNodes_send.size());
      std::vector<unsigned> counts_recv(4 * static_cast<std::size_t>(nprocs), 0);
      MPI_Allgather(counts_send, 4, MPI_UNSIGNED, counts_recv.data(), 4, MPI_UNSIGNED, MPI_COMM_WORLD);

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        if (counts_recv[4 * static_cast<std::size_t>(jproc) + 3] != counts_send[3]) {
          std::cerr << "BuildEdgePoints: rank=" << iproc << " inconsistent initial nnodes=" << counts_send[3] << " rank=" << jproc << " nnodes=" << counts_recv[4 * static_cast<std::size_t>(jproc) + 3] << std::endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      }

      for (unsigned kind = 0; kind < 3; kind++) {
        unsigned nextNode = nnodes;

        for (unsigned jproc = 0; jproc < nprocs; jproc++) {
          unsigned count_recv = counts_recv[4 * static_cast<std::size_t>(jproc) + kind];
          if (count_recv > UINT_MAX - nnodes) MPI_Abort(MPI_COMM_WORLD, 1);
          if (jproc < iproc) nextNode += count_recv;
          nnodes += count_recv;
        }

        for (unsigned i = pointOffset[kind]; i < pointOffset[kind + 1]; i++) {
          if (pointNodes_send[i] == UINT_MAX) pointNodes_send[i] = nextNode++;
        }
      }

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        for (unsigned pos_send = 0; pos_send < pointDofs_send[jproc].size(); pos_send += 3) {
          unsigned iel_send = pointDofs_send[jproc][pos_send];
          unsigned inode_send = pointDofs_send[jproc][pos_send + 1];
          unsigned node_send = pointNodes_send[pointDofs_send[jproc][pos_send + 2]];
          pointDofs_send[jproc][pos_send + 2] = node_send;

          if (jproc == iproc) SetElementDofIndex(iel_send, inode_send, node_send);
          else _ghostElementLayer[jproc].at(iel_send)[inode_send] = node_send;
        }
      }

      pointDofs_send[iproc].clear();
      std::vector<std::vector<unsigned>> pointDofs_recv;
      ExchangeNodeBuffers(pointDofs_send, pointDofs_recv, iproc, nprocs, 103);

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        for (unsigned pos_recv = 0; pos_recv < pointDofs_recv[jproc].size(); pos_recv += 3) {
          unsigned iel_recv = pointDofs_recv[jproc][pos_recv];
          unsigned inode_recv = pointDofs_recv[jproc][pos_recv + 1];
          unsigned node_recv = pointDofs_recv[jproc][pos_recv + 2];
          unsigned node = GetElementDofIndex(iel_recv, inode_recv);
          if (node != UINT_MAX && node != node_recv) MPI_Abort(MPI_COMM_WORLD, 1);
          SetElementDofIndex(iel_recv, inode_recv, node_recv);
        }
      }
    }

    for (unsigned iel = _elementOffset[iproc]; iel < _elementOffset[iproc + 1]; iel++) {
      for (unsigned inode = 0; inode < GetElementDofNumber(iel, 2); inode++) {
        unsigned node = GetElementDofIndex(iel, inode);
        if (node == UINT_MAX || node >= nnodes) {
          std::cerr << "BuildEdgePoints: rank=" << iproc << " incomplete owned connectivity element=" << iel << " level=" << GetElementLevel(iel) << " meshLevel=" << _level << " slot=" << inode << " node=" << node << std::endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      }
    }

    {
      std::vector<std::vector<unsigned>> connectivity_send(nprocs);
      std::vector<std::vector<unsigned>> connectivity_recv;

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        for (unsigned i = 0; i < _ghostElements_send[jproc].size(); i++) {
          unsigned iel_send = _ghostElements_send[jproc][i];
          if (!GetIfFatherHasBeenRefined(iel_send)) continue;
          unsigned nDofs_send = GetElementDofNumber(iel_send, 2);
          if (connectivity_send[jproc].size() > static_cast<unsigned>(INT_MAX - 2) - nDofs_send) MPI_Abort(MPI_COMM_WORLD, 1);
          connectivity_send[jproc].push_back(iel_send);
          connectivity_send[jproc].push_back(nDofs_send);

          for (unsigned inode = 0; inode < nDofs_send; inode++) connectivity_send[jproc].push_back(GetElementDofIndex(iel_send, inode));
        }
      }

      ExchangeNodeBuffers(connectivity_send, connectivity_recv, iproc, nprocs, 104);

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        unsigned pos_recv = 0;

        while (pos_recv < connectivity_recv[jproc].size()) {
          unsigned iel_recv = connectivity_recv[jproc][pos_recv++];
          unsigned nDofs_recv = connectivity_recv[jproc][pos_recv++];
          std::vector<unsigned>& dofs_recv = _ghostElementLayer[jproc].at(iel_recv);
          if (nDofs_recv != dofs_recv.size()) MPI_Abort(MPI_COMM_WORLD, 1);
          std::copy(connectivity_recv[jproc].begin() + pos_recv, connectivity_recv[jproc].begin() + pos_recv + nDofs_recv, dofs_recv.begin());
          pos_recv += nDofs_recv;
        }
      }
    }

    for (unsigned jproc = 0; jproc < nprocs; jproc++) {
      for (std::unordered_map<unsigned, std::vector<unsigned>>::const_iterator it = _ghostElementLayer[jproc].begin(); it != _ghostElementLayer[jproc].end(); ++it) {
        for (unsigned inode = 0; inode < it->second.size(); inode++) {
          unsigned node = it->second[inode];
          if (node == UINT_MAX || node >= nnodes) {
            std::cerr << "BuildEdgePoints: rank=" << iproc << " incomplete ghost connectivity owner=" << jproc << " element=" << it->first << " level=" << GetGhostElementLevel(jproc, it->first) << " meshLevel=" << _level << " slot=" << inode << " node=" << node << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
          }
        }
      }
    }

    SetNodeNumber(nnodes);
    msh->SetNumberOfNodes(nnodes);
  }

  // ********************************************

  static void DebugConnectivityError(const unsigned iproc, const std::string& message, unsigned long long& errors_send) {
    if (errors_send < 12) std::cerr << "[DebugCheckConnectivity rank " << iproc << "] " << message << std::endl;
    errors_send++;
  }

  static void DebugConnectivityStage(const unsigned iproc, const char* stage, const unsigned long long errors_send) {
    unsigned long long errors_recv = 0;
    MPI_Allreduce(&errors_send, &errors_recv, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (errors_recv != 0) {
      if (iproc == 0) std::cerr << "[DebugCheckConnectivity] failed during " << stage << ": " << errors_recv << " errors; at most 12 details per rank." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  static bool DebugConnectivityAppendAllowed(const std::vector<unsigned>& buffer, const unsigned count, const unsigned iproc, const char* stage, unsigned long long& errors_send) {
    if (count > static_cast<unsigned>(INT_MAX) || buffer.size() > static_cast<unsigned>(INT_MAX) - count) {
      DebugConnectivityError(iproc, std::string("MPI buffer exceeds INT_MAX during ") + stage, errors_send);
      return false;
    }
    return true;
  }

  static void DebugConnectivityCheckRow(const std::vector<unsigned>& row, const unsigned iel, const unsigned owner, const unsigned level, const unsigned required, const unsigned nnodes, const unsigned iproc, unsigned long long& errors_send) {
    for (unsigned inode = 0; inode < row.size(); inode++) {
      unsigned node = row[inode];
      if (node == UINT_MAX) {
        if (inode < required) {
          std::ostringstream message;
          message << (owner == iproc ? "owned" : "ghost") << " element=" << iel << " owner=" << owner << " level=" << level << " slot=" << inode << " remains UINT_MAX; required slots=[0," << required << ").";
          DebugConnectivityError(iproc, message.str(), errors_send);
        }
        continue;
      }
      if (node >= nnodes) {
        std::ostringstream message;
        message << (owner == iproc ? "owned" : "ghost") << " element=" << iel << " owner=" << owner << " level=" << level << " slot=" << inode << " node=" << node << " is outside [0," << nnodes << ").";
        DebugConnectivityError(iproc, message.str(), errors_send);
      }
      for (unsigned jnode = 0; jnode < inode; jnode++) {
        if (row[jnode] != node) continue;
        std::ostringstream message;
        message << (owner == iproc ? "owned" : "ghost") << " element=" << iel << " owner=" << owner << " level=" << level << " slots=" << jnode << "," << inode << " duplicate node=" << node << ".";
        DebugConnectivityError(iproc, message.str(), errors_send);
      }
    }
  }

  void elem::DebugCheckConnectivity(const unsigned nnodes, const unsigned connectivityType) {
    const unsigned iproc = _iproc;
    const unsigned nprocs = _nprocs;
    unsigned long long errors_send = 0;

    {
      bool offsetsValid = _elementOffset.size() == static_cast<std::size_t>(nprocs) + 1;
      unsigned values_send[5] = {nnodes, connectivityType, _level, offsetsValid ? _elementOffset[iproc] : UINT_MAX, offsetsValid ? _elementOffset[iproc + 1] : UINT_MAX};
      std::vector<unsigned> values_recv(5 * static_cast<std::size_t>(nprocs), 0);
      MPI_Allgather(values_send, 5, MPI_UNSIGNED, values_recv.data(), 5, MPI_UNSIGNED, MPI_COMM_WORLD);

      if (nnodes != values_recv[0] || connectivityType != values_recv[1] || _level != values_recv[2]) {
        std::ostringstream message;
        message << "rank settings disagree with rank 0: nnodes=" << nnodes << "/" << values_recv[0] << " connectivityType=" << connectivityType << "/" << values_recv[1] << " meshLevel=" << _level << "/" << values_recv[2] << ".";
        DebugConnectivityError(iproc, message.str(), errors_send);
      }
      if (connectivityType > 2) DebugConnectivityError(iproc, "connectivityType must be 0, 1, or 2.", errors_send);
      if (!offsetsValid) DebugConnectivityError(iproc, "_elementOffset does not contain nprocs + 1 entries.", errors_send);
      if (_ghostElementLayer.size() != nprocs || _ghostElementType.size() != nprocs || _ghostElementLevel.size() != nprocs) DebugConnectivityError(iproc, "ghost layer/type/level arrays do not contain nprocs entries; build the current ghost layer first.", errors_send);

      if (offsetsValid) {
        if (_elementOffset[0] != 0) DebugConnectivityError(iproc, "_elementOffset[0] must be zero.", errors_send);
        for (unsigned jproc = 0; jproc < nprocs; jproc++) {
          if (_elementOffset[jproc] > _elementOffset[jproc + 1] || _elementOffset[jproc] != values_recv[5 * static_cast<std::size_t>(jproc) + 3] || _elementOffset[jproc + 1] != values_recv[5 * static_cast<std::size_t>(jproc) + 4]) {
            std::ostringstream message;
            message << "element owner interval for rank=" << jproc << " is inconsistent: local=[" << _elementOffset[jproc] << "," << _elementOffset[jproc + 1] << ") reported=[" << values_recv[5 * static_cast<std::size_t>(jproc) + 3] << "," << values_recv[5 * static_cast<std::size_t>(jproc) + 4] << ").";
            DebugConnectivityError(iproc, message.str(), errors_send);
          }
        }
      }
      DebugConnectivityStage(iproc, "rank and storage validation", errors_send);
    }

    for (unsigned iel = _elementOffset[iproc]; iel < _elementOffset[iproc + 1]; iel++) {
      unsigned type = GetElementType(iel);
      if (type >= 6) {
        std::ostringstream message;
        message << "owned element=" << iel << " has invalid type=" << type << ".";
        DebugConnectivityError(iproc, message.str(), errors_send);
        continue;
      }
      std::vector<unsigned> row(GetNVE(type, 2));
      for (unsigned inode = 0; inode < row.size(); inode++) row[inode] = GetElementDofIndex(iel, inode);
      DebugConnectivityCheckRow(row, iel, iproc, GetElementLevel(iel), GetNVE(type, connectivityType), nnodes, iproc, errors_send);
    }

    for (unsigned jproc = 0; jproc < nprocs; jproc++) {
      if (_ghostElementType[jproc].size() != _ghostElementLayer[jproc].size() || _ghostElementLevel[jproc].size() != _ghostElementLayer[jproc].size()) {
        std::ostringstream message;
        message << "ghost metadata sizes disagree with connectivity layer for owner=" << jproc << ".";
        DebugConnectivityError(iproc, message.str(), errors_send);
      }
      for (std::unordered_map<unsigned, std::vector<unsigned>>::const_iterator it = _ghostElementLayer[jproc].begin(); it != _ghostElementLayer[jproc].end(); ++it) {
        unsigned iel = it->first;
        std::unordered_map<unsigned, unsigned>::const_iterator typeIt = _ghostElementType[jproc].find(iel);
        std::unordered_map<unsigned, unsigned>::const_iterator levelIt = _ghostElementLevel[jproc].find(iel);
        if (jproc == iproc || iel < _elementOffset[jproc] || iel >= _elementOffset[jproc + 1]) {
          std::ostringstream message;
          message << "ghost element=" << iel << " is in wrong owner layer=" << jproc << " with interval=[" << _elementOffset[jproc] << "," << _elementOffset[jproc + 1] << ").";
          DebugConnectivityError(iproc, message.str(), errors_send);
        }
        if (typeIt == _ghostElementType[jproc].end() || levelIt == _ghostElementLevel[jproc].end()) {
          std::ostringstream message;
          message << "ghost element=" << iel << " owner=" << jproc << " is missing type or level metadata.";
          DebugConnectivityError(iproc, message.str(), errors_send);
          continue;
        }
        unsigned type = typeIt->second;
        if (type >= 6) {
          std::ostringstream message;
          message << "ghost element=" << iel << " owner=" << jproc << " has invalid type=" << type << ".";
          DebugConnectivityError(iproc, message.str(), errors_send);
          continue;
        }
        if (it->second.size() != GetNVE(type, 2)) {
          std::ostringstream message;
          message << "ghost element=" << iel << " owner=" << jproc << " type=" << type << " has rowSize=" << it->second.size() << " expected=" << GetNVE(type, 2) << ".";
          DebugConnectivityError(iproc, message.str(), errors_send);
          continue;
        }
        DebugConnectivityCheckRow(it->second, iel, jproc, levelIt->second, GetNVE(type, connectivityType), nnodes, iproc, errors_send);
      }
    }
    DebugConnectivityStage(iproc, "owned and ghost connectivity validation", errors_send);

    std::vector<std::vector<unsigned>> requirements_send(nprocs);
    std::vector<std::vector<unsigned>> midpoint_send(nprocs);

    {
      std::map<std::pair<unsigned, unsigned>, std::vector<std::array<unsigned, 5>>> edges;

      {
        std::vector<std::vector<unsigned>> edges_send(nprocs);
        std::vector<std::vector<unsigned>> edges_recv;

        for (unsigned iel = _elementOffset[iproc]; iel < _elementOffset[iproc + 1]; iel++) {
          unsigned type = GetElementType(iel);
          unsigned first = GetNVE(type, 0);
          unsigned last = GetNVE(type, 1);

          for (unsigned inode = first; inode < last; inode++) {
            unsigned im = GetElementDofIndex(iel, edge2VerticesMapping[type][inode - first][0]);
            unsigned ip = GetElementDofIndex(iel, edge2VerticesMapping[type][inode - first][1]);
            if (ip < im) std::swap(im, ip);
            unsigned jproc = static_cast<unsigned>((static_cast<unsigned long long>(im) * 2654435761ULL + static_cast<unsigned long long>(ip) * 2246822519ULL) % nprocs);
            if (!DebugConnectivityAppendAllowed(edges_send[jproc], 6, iproc, "edge routing", errors_send)) continue;
            edges_send[jproc].insert(edges_send[jproc].end(), {im, ip, GetElementLevel(iel), GetElementDofIndex(iel, inode), iel, inode});
          }
        }
        DebugConnectivityStage(iproc, "edge routing buffer preparation", errors_send);
        ExchangeNodeBuffers(edges_send, edges_recv, iproc, nprocs, 105);

        for (unsigned jproc = 0; jproc < nprocs; jproc++) {
          if (edges_recv[jproc].size() % 6 != 0) {
            DebugConnectivityError(iproc, "malformed authoritative edge-record buffer.", errors_send);
            continue;
          }
          for (unsigned pos_recv = 0; pos_recv < edges_recv[jproc].size(); pos_recv += 6) {
            std::array<unsigned, 5> record = {{edges_recv[jproc][pos_recv + 2], edges_recv[jproc][pos_recv + 3], edges_recv[jproc][pos_recv + 4], edges_recv[jproc][pos_recv + 5], jproc}};
            edges[std::make_pair(edges_recv[jproc][pos_recv], edges_recv[jproc][pos_recv + 1])].push_back(record);
          }
        }
        DebugConnectivityStage(iproc, "edge record parsing", errors_send);
      }

      for (std::map<std::pair<unsigned, unsigned>, std::vector<std::array<unsigned, 5>>>::const_iterator it = edges.begin(); it != edges.end(); ++it) {
        std::map<unsigned, std::array<unsigned, 4>> levelMidpoint;
        std::vector<unsigned> owners;

        for (unsigned i = 0; i < it->second.size(); i++) {
          const std::array<unsigned, 5>& record = it->second[i];
          if (std::find(owners.begin(), owners.end(), record[4]) == owners.end()) owners.push_back(record[4]);
          if (record[1] == UINT_MAX) continue;
          std::map<unsigned, std::array<unsigned, 4>>::const_iterator known = levelMidpoint.find(record[0]);
          if (known == levelMidpoint.end()) {
            std::array<unsigned, 4> midpoint = {{record[1], record[2], record[3], record[4]}};
            levelMidpoint[record[0]] = midpoint;
          }
          else if (known->second[0] != record[1]) {
            std::ostringstream message;
            message << "edge=(" << it->first.first << "," << it->first.second << ") level=" << record[0] << " has different midpoint IDs: owner=" << known->second[3] << " element=" << known->second[1] << " slot=" << known->second[2] << " node=" << known->second[0] << "; owner=" << record[4] << " element=" << record[2] << " slot=" << record[3] << " node=" << record[1] << ".";
            DebugConnectivityError(iproc, message.str(), errors_send);
          }
          unsigned jproc = record[1] % nprocs;
          if (!DebugConnectivityAppendAllowed(midpoint_send[jproc], 7, iproc, "midpoint routing", errors_send)) continue;
          midpoint_send[jproc].insert(midpoint_send[jproc].end(), {record[1], it->first.first, it->first.second, record[2], record[3], record[0], record[4]});
        }

        for (unsigned i = 0; i < owners.size(); i++) {
          unsigned jproc = owners[i];
          for (unsigned j = 0; j < it->second.size(); j++) {
            const std::array<unsigned, 5>& record = it->second[j];
            if (record[4] == jproc) continue;
            if (!DebugConnectivityAppendAllowed(requirements_send[jproc], 5, iproc, "edge-neighbor routing", errors_send)) continue;
            requirements_send[jproc].insert(requirements_send[jproc].end(), {record[4], record[2], record[0], it->first.first, it->first.second});
          }
        }
      }
      DebugConnectivityStage(iproc, "global edge agreement", errors_send);
    }

    {
      std::vector<std::vector<unsigned>> requirements_recv;
      ExchangeNodeBuffers(requirements_send, requirements_recv, iproc, nprocs, 106);
      std::vector<std::vector<unsigned>>().swap(requirements_send);

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        if (requirements_recv[jproc].size() % 5 != 0) {
          DebugConnectivityError(iproc, "malformed edge-neighbor requirement buffer.", errors_send);
          continue;
        }
        for (unsigned pos_recv = 0; pos_recv < requirements_recv[jproc].size(); pos_recv += 5) {
          unsigned owner = requirements_recv[jproc][pos_recv];
          unsigned iel = requirements_recv[jproc][pos_recv + 1];
          if (owner >= nprocs || owner == iproc) {
            DebugConnectivityError(iproc, "invalid owner in edge-neighbor requirement buffer.", errors_send);
            continue;
          }
          if (_ghostElementLayer[owner].find(iel) == _ghostElementLayer[owner].end()) {
            std::ostringstream message;
            message << "missing ghost element=" << iel << " owner=" << owner << " level=" << requirements_recv[jproc][pos_recv + 2] << " sharing edge=(" << requirements_recv[jproc][pos_recv + 3] << "," << requirements_recv[jproc][pos_recv + 4] << ") with an owned element.";
            DebugConnectivityError(iproc, message.str(), errors_send);
          }
        }
      }
      DebugConnectivityStage(iproc, "edge-neighbor ghost completeness", errors_send);
    }

    {
      std::vector<std::vector<unsigned>> midpoint_recv;
      ExchangeNodeBuffers(midpoint_send, midpoint_recv, iproc, nprocs, 107);
      std::vector<std::vector<unsigned>>().swap(midpoint_send);
      std::map<unsigned, std::array<unsigned, 6>> nodes;

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        if (midpoint_recv[jproc].size() % 7 != 0) {
          DebugConnectivityError(iproc, "malformed midpoint-record buffer.", errors_send);
          continue;
        }
        for (unsigned pos_recv = 0; pos_recv < midpoint_recv[jproc].size(); pos_recv += 7) {
          unsigned node = midpoint_recv[jproc][pos_recv];
          std::array<unsigned, 6> record = {{midpoint_recv[jproc][pos_recv + 1], midpoint_recv[jproc][pos_recv + 2], midpoint_recv[jproc][pos_recv + 3], midpoint_recv[jproc][pos_recv + 4], midpoint_recv[jproc][pos_recv + 5], midpoint_recv[jproc][pos_recv + 6]}};
          std::map<unsigned, std::array<unsigned, 6>>::const_iterator known = nodes.find(node);
          if (known == nodes.end()) nodes[node] = record;
          else if (known->second[0] != record[0] || known->second[1] != record[1]) {
            std::ostringstream message;
            message << "midpoint node=" << node << " is reused for different edges: owner=" << known->second[5] << " element=" << known->second[2] << " slot=" << known->second[3] << " level=" << known->second[4] << " edge=(" << known->second[0] << "," << known->second[1] << "); owner=" << record[5] << " element=" << record[2] << " slot=" << record[3] << " level=" << record[4] << " edge=(" << record[0] << "," << record[1] << ").";
            DebugConnectivityError(iproc, message.str(), errors_send);
          }
        }
      }
      DebugConnectivityStage(iproc, "global midpoint uniqueness", errors_send);
    }

    {
      std::vector<std::vector<unsigned>> ghosts_send(nprocs);
      std::vector<std::vector<unsigned>> ghosts_recv;

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        for (std::unordered_map<unsigned, std::vector<unsigned>>::const_iterator it = _ghostElementLayer[jproc].begin(); it != _ghostElementLayer[jproc].end(); ++it) {
          unsigned count = static_cast<unsigned>(it->second.size());
          if (!DebugConnectivityAppendAllowed(ghosts_send[jproc], 4 + count, iproc, "ghost-owner comparison", errors_send)) continue;
          ghosts_send[jproc].insert(ghosts_send[jproc].end(), {it->first, _ghostElementType[jproc].at(it->first), _ghostElementLevel[jproc].at(it->first), count});
          ghosts_send[jproc].insert(ghosts_send[jproc].end(), it->second.begin(), it->second.end());
        }
      }
      DebugConnectivityStage(iproc, "ghost-owner comparison buffer preparation", errors_send);
      ExchangeNodeBuffers(ghosts_send, ghosts_recv, iproc, nprocs, 108);

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        unsigned pos_recv = 0;
        while (pos_recv < ghosts_recv[jproc].size()) {
          if (ghosts_recv[jproc].size() - pos_recv < 4) {
            DebugConnectivityError(iproc, "truncated ghost-owner comparison header.", errors_send);
            break;
          }
          unsigned iel = ghosts_recv[jproc][pos_recv++];
          unsigned type = ghosts_recv[jproc][pos_recv++];
          unsigned level = ghosts_recv[jproc][pos_recv++];
          unsigned count = ghosts_recv[jproc][pos_recv++];
          if (count > ghosts_recv[jproc].size() - pos_recv) {
            DebugConnectivityError(iproc, "truncated ghost-owner comparison row.", errors_send);
            break;
          }
          if (iel < _elementOffset[iproc] || iel >= _elementOffset[iproc + 1]) {
            std::ostringstream message;
            message << "rank=" << jproc << " claims this rank owns ghost element=" << iel << " outside [" << _elementOffset[iproc] << "," << _elementOffset[iproc + 1] << ").";
            DebugConnectivityError(iproc, message.str(), errors_send);
            pos_recv += count;
            continue;
          }
          unsigned actualType = GetElementType(iel);
          unsigned actualLevel = GetElementLevel(iel);
          unsigned actualCount = GetNVE(actualType, 2);
          if (type != actualType || level != actualLevel || count != actualCount) {
            std::ostringstream message;
            message << "rank=" << jproc << " ghost metadata differs from owner=" << iproc << " element=" << iel << ": type=" << type << "/" << actualType << " level=" << level << "/" << actualLevel << " rowSize=" << count << "/" << actualCount << ".";
            DebugConnectivityError(iproc, message.str(), errors_send);
          }
          else {
            for (unsigned inode = 0; inode < count; inode++) {
              unsigned ghostNode = ghosts_recv[jproc][pos_recv + inode];
              unsigned ownedNode = GetElementDofIndex(iel, inode);
              if (ghostNode == ownedNode) continue;
              std::ostringstream message;
              message << "stale ghost on rank=" << jproc << " element=" << iel << " owner=" << iproc << " level=" << level << " slot=" << inode << " ghostNode=" << ghostNode << " ownedNode=" << ownedNode << ".";
              DebugConnectivityError(iproc, message.str(), errors_send);
            }
          }
          pos_recv += count;
        }
      }
      DebugConnectivityStage(iproc, "ghost-owner connectivity agreement", errors_send);
    }

    if (iproc == 0) std::cout << "[DebugCheckConnectivity] passed: connectivityType=" << connectivityType << " nnodes=" << nnodes << " elements=" << _elementOffset.back() << " ranks=" << nprocs << "." << std::endl;
  }

  void elem::DebugCheckNodeLayout(Mesh *msh) {
    const unsigned iproc = _iproc;
    const unsigned nprocs = _nprocs;
    unsigned long long errors_send = 0;

    if (msh == nullptr) {
      DebugConnectivityError(iproc, "type-2 node layout: Mesh pointer is null.", errors_send);
    }
    else {
      if (msh->el != this) DebugConnectivityError(iproc, "type-2 node layout: Mesh::el is not this elem object.", errors_send);
      if (msh->processor_id() != iproc || msh->n_processors() != nprocs) DebugConnectivityError(iproc, "type-2 node layout: Mesh and elem process information differs.", errors_send);
      if (msh->_elementOffset.size() != static_cast<std::size_t>(nprocs) + 1) DebugConnectivityError(iproc, "type-2 node layout: Mesh::_elementOffset must contain nprocs + 1 entries.", errors_send);
      if (msh->_dofOffset[2].size() != static_cast<std::size_t>(nprocs) + 1) DebugConnectivityError(iproc, "type-2 node layout: Mesh::_dofOffset[2] must contain nprocs + 1 entries.", errors_send);
      if (msh->_ownSize[2].size() != nprocs) DebugConnectivityError(iproc, "type-2 node layout: Mesh::_ownSize[2] must contain nprocs entries.", errors_send);
    }
    if (_elementOffset.size() != static_cast<std::size_t>(nprocs) + 1) DebugConnectivityError(iproc, "type-2 node layout: elem::_elementOffset must contain nprocs + 1 entries.", errors_send);
    DebugConnectivityStage(iproc, "type-2 node layout storage", errors_send);

    if (_elementOffset != msh->_elementOffset) DebugConnectivityError(iproc, "type-2 node layout: Mesh and elem element offsets differ.", errors_send);
    if (_elementOffset[0] != 0 || msh->_dofOffset[2][0] != 0) DebugConnectivityError(iproc, "type-2 node layout: element and node offsets must start at zero.", errors_send);
    if (_elementOffset.back() != GetElementNumber() || _elementOffset.back() != msh->GetNumberOfElements()) DebugConnectivityError(iproc, "type-2 node layout: element counts disagree with the final element offset.", errors_send);

    unsigned nnodes = msh->GetNumberOfNodes();
    if (nnodes != GetNodeNumber() || nnodes != msh->_dofOffset[2].back()) {
      std::ostringstream message;
      message << "type-2 node layout: node totals disagree: Mesh=" << nnodes << " elem=" << GetNodeNumber() << " final type-2 offset=" << msh->_dofOffset[2].back() << ". Complete node renumbering and distribute its metadata before assembly.";
      DebugConnectivityError(iproc, message.str(), errors_send);
    }
    if (_elementType.begin() != _elementOffset[iproc] || _elementType.end() != _elementOffset[iproc + 1] || _elementDof.begin() != _elementOffset[iproc] || _elementDof.end() != _elementOffset[iproc + 1]) DebugConnectivityError(iproc, "type-2 node layout: owned element storage does not match this process's element interval.", errors_send);

    for (unsigned jproc = 0; jproc < nprocs; jproc++) {
      if (_elementOffset[jproc] > _elementOffset[jproc + 1] || msh->_dofOffset[2][jproc] > msh->_dofOffset[2][jproc + 1]) {
        std::ostringstream message;
        message << "type-2 node layout: decreasing element or node offsets at rank=" << jproc << ".";
        DebugConnectivityError(iproc, message.str(), errors_send);
      }
      else if (msh->_ownSize[2][jproc] != msh->_dofOffset[2][jproc + 1] - msh->_dofOffset[2][jproc]) {
        std::ostringstream message;
        message << "type-2 node layout: _ownSize[2] disagrees with the node interval at rank=" << jproc << ".";
        DebugConnectivityError(iproc, message.str(), errors_send);
      }
    }
    DebugConnectivityStage(iproc, "type-2 node layout metadata", errors_send);

    {
      unsigned values_send[6] = {nnodes, GetElementNumber(), _elementOffset[iproc], _elementOffset[iproc + 1], msh->_dofOffset[2][iproc], msh->_dofOffset[2][iproc + 1]};
      std::vector<unsigned> values_recv(6 * static_cast<std::size_t>(nprocs), 0);
      MPI_Allgather(values_send, 6, MPI_UNSIGNED, values_recv.data(), 6, MPI_UNSIGNED, MPI_COMM_WORLD);

      for (unsigned jproc = 0; jproc < nprocs; jproc++) {
        std::size_t pos_recv = 6 * static_cast<std::size_t>(jproc);
        if (nnodes != values_recv[pos_recv] || GetElementNumber() != values_recv[pos_recv + 1] || _elementOffset[jproc] != values_recv[pos_recv + 2] || _elementOffset[jproc + 1] != values_recv[pos_recv + 3] || msh->_dofOffset[2][jproc] != values_recv[pos_recv + 4] || msh->_dofOffset[2][jproc + 1] != values_recv[pos_recv + 5]) {
          std::ostringstream message;
          message << "type-2 node layout: counts or ownership intervals disagree with authoritative metadata from rank=" << jproc << ".";
          DebugConnectivityError(iproc, message.str(), errors_send);
        }
      }
    }
    DebugConnectivityStage(iproc, "type-2 node layout agreement across ranks", errors_send);

    for (unsigned iel = _elementOffset[iproc]; iel < _elementOffset[iproc + 1]; iel++) {
      unsigned type = GetElementType(iel);
      if (type >= 6) {
        std::ostringstream message;
        message << "type-2 node layout: element=" << iel << " has invalid type=" << type << ".";
        DebugConnectivityError(iproc, message.str(), errors_send);
        continue;
      }
      unsigned nDofs = GetNVE(type, 2);
      if (_elementDof.end(iel) < nDofs) {
        std::ostringstream message;
        message << "type-2 node layout: element=" << iel << " has " << _elementDof.end(iel) << " stored entries; its full connectivity requires " << nDofs << ".";
        DebugConnectivityError(iproc, message.str(), errors_send);
        continue;
      }

      for (unsigned inode = 0; inode < nDofs; inode++) {
        unsigned node = GetElementDofIndex(iel, inode);
        if (node == UINT_MAX || node >= nnodes) {
          std::ostringstream message;
          message << "type-2 node layout: element=" << iel << " slot=" << inode << " node=" << node << " is " << (node == UINT_MAX ? "UINT_MAX" : "outside the final node range") << "; required range=[0," << nnodes << ").";
          DebugConnectivityError(iproc, message.str(), errors_send);
        }
        else if (node >= msh->_dofOffset[2][iproc + 1]) {
          unsigned owner = static_cast<unsigned>(std::upper_bound(msh->_dofOffset[2].begin(), msh->_dofOffset[2].end(), node) - msh->_dofOffset[2].begin() - 1);
          std::ostringstream message;
          message << "type-2 node layout: element=" << iel << " slot=" << inode << " node=" << node << " belongs to higher rank=" << owner << "; this rank's node interval=[" << msh->_dofOffset[2][iproc] << "," << msh->_dofOffset[2][iproc + 1] << "). Existing FEMuS sparsity construction requires local or lower-rank node ownership.";
          DebugConnectivityError(iproc, message.str(), errors_send);
        }
      }
    }
    DebugConnectivityStage(iproc, "type-2 node layout connectivity", errors_send);
    if (iproc == 0) std::cout << "[DebugCheckNodeLayout] passed: type-2 node layout, nnodes=" << nnodes << " elements=" << GetElementNumber() << " ranks=" << nprocs << "; lower-order and system DOF maps were not checked." << std::endl;
  }

}
