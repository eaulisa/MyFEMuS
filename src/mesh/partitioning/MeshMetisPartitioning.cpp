/*=========================================================================

 Program: FEMUS
 Module: MeshMetisPartitioning
 Authors: Simone Bnà, Eugenio Aulisa

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MeshMetisPartitioning.hpp"
#include "Mesh.hpp"
#include "FemusConfig.hpp"

#ifdef HAVE_METIS
#include "metis.h"
#endif

#ifdef HAVE_PARMETIS
#include <parmetis.h>
#endif

//C++ include
#include <iostream>


namespace femus {

  using std::cout;
  using std::endl;


  MeshMetisPartitioning::MeshMetisPartitioning(Mesh& mesh) : MeshPartitioning(mesh) {

  }

//------------------------------------------------------------------------------------------------------
void MeshMetisPartitioning::DoPartition(std::vector<unsigned>& partition, const bool& AMR) {

  idx_t nnodes = static_cast<idx_t>(_mesh.GetNumberOfNodes());
  idx_t nelem  = static_cast<idx_t>(_mesh.GetNumberOfElements());

  // If partition already has one value for each element,
  // use it as the inherited/father partition.
  const bool useInitialPartition =
      (partition.size() == static_cast<unsigned>(nelem));

  std::vector<idx_t> epart(nelem);


  //========================================================================
  // Trivial cases
  //========================================================================

  if(_nprocs == 1) {

    for(idx_t i = 0; i < nelem; ++i) {
      epart[i] = 0;
    }
  }

  else if(_nprocs > static_cast<unsigned>(nelem)) {

    std::cerr
        << "Error in MeshMetisPartitioning::DoPartition: "
        << "number of processes "
        << _nprocs
        << " is greater than number of elements "
        << nelem
        << std::endl;

    abort();
  }

  else if(_nprocs == static_cast<unsigned>(nelem)) {

    for(idx_t i = 0; i < nelem; ++i) {
      epart[i] = i;
    }
  }

  else {

#ifndef HAVE_METIS

    std::cerr
        << "Fatal error: METIS library was not found."
        << std::endl;

    exit(1);

#else

    //========================================================================
    // Build METIS mesh connectivity
    //========================================================================

    unsigned eind_size =
        _mesh.el->GetElementNumber("Hex")      * NVE[0][2]
      + _mesh.el->GetElementNumber("Tet")      * NVE[1][2]
      + _mesh.el->GetElementNumber("Wedge")    * NVE[2][2]
      + _mesh.el->GetElementNumber("Quad")     * NVE[3][2]
      + _mesh.el->GetElementNumber("Triangle") * NVE[4][2]
      + _mesh.el->GetElementNumber("Line")     * NVE[5][2];

    std::vector<idx_t> eptr(nelem + 1);
    std::vector<idx_t> eind(eind_size);

    eptr[0] = 0;

    unsigned counter = 0;

    for(idx_t iel = 0; iel < nelem; ++iel) {

      unsigned ndofs =
          _mesh.el->GetElementDofNumber(
              static_cast<unsigned>(iel), 2);

      eptr[iel + 1] =
          eptr[iel] + static_cast<idx_t>(ndofs);

      for(unsigned inode = 0; inode < ndofs; ++inode) {

        eind[counter] =
            static_cast<idx_t>(
                _mesh.el->GetElementDofIndex(
                    static_cast<unsigned>(iel),
                    inode
                )
            );

        ++counter;
      }
    }


    idx_t ncommon =
        (AMR || _mesh.GetDimension() == 1)
        ? static_cast<idx_t>(1)
        : static_cast<idx_t>(_mesh.GetDimension() + 1);


    //========================================================================
    // EXISTING PARTITION:
    //
    // Adaptive repartitioning starting from the father partition.
    //========================================================================

    if(useInitialPartition) {

#ifndef HAVE_PARMETIS

      std::cerr
          << "Fatal error: an initial partition was supplied, "
          << "but ParMETIS was not found."
          << std::endl;

      exit(1);

#else

      int iproc = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

      idx_t nparts =
          static_cast<idx_t>(_nprocs);


      //======================================================================
      // Build the original element dual graph
      //======================================================================

      idx_t numflag = 0;

      idx_t* xadjGlobal   = NULL;
      idx_t* adjncyGlobal = NULL;

      int err = METIS_MeshToDual(
          &nelem,
          &nnodes,
          eptr.data(),
          eind.data(),
          &ncommon,
          &numflag,
          &xadjGlobal,
          &adjncyGlobal
      );

      if(err != METIS_OK) {

        std::cerr
            << "METIS_MeshToDual failed with error "
            << err
            << std::endl;

        exit(1);
      }


      //======================================================================
      // Validate inherited partition and count elements in each partition
      //======================================================================

      std::vector<idx_t> initialCount(nparts, 0);

      for(idx_t iel = 0; iel < nelem; ++iel) {

        if(partition[iel] >= _nprocs) {

          std::cerr
              << "Invalid inherited partition for element "
              << iel
              << ": "
              << partition[iel]
              << std::endl;

          abort();
        }

        initialCount[partition[iel]]++;
      }


      if(iproc == 0) {

        std::cout << std::endl;
        std::cout << "Inherited partition:" << std::endl;

        for(idx_t p = 0; p < nparts; ++p) {

          std::cout
              << "  partition "
              << p
              << " : "
              << initialCount[p]
              << " elements"
              << std::endl;
        }
      }


      //======================================================================
      // Construct a NEW ParMETIS numbering grouped by inherited partition.
      //
      // New numbering:
      //
      //   [ all partition 0 elements ]
      //   [ all partition 1 elements ]
      //   ...
      //
      // Therefore processor p stores exactly the vertices belonging to
      // inherited partition p.
      //======================================================================

      std::vector<idx_t> vtxdist(_nprocs + 1);

      vtxdist[0] = 0;

      for(unsigned p = 0; p < _nprocs; ++p) {

        vtxdist[p + 1] =
            vtxdist[p] + initialCount[p];
      }


      std::vector<idx_t> oldToNew(nelem);
      std::vector<idx_t> newToOld(nelem);

      std::vector<idx_t> nextPosition(_nprocs);

      for(unsigned p = 0; p < _nprocs; ++p) {

        nextPosition[p] = vtxdist[p];
      }


      for(idx_t oldElement = 0;
          oldElement < nelem;
          ++oldElement) {

        unsigned p =
            partition[oldElement];

        idx_t newElement =
            nextPosition[p]++;

        oldToNew[oldElement] =
            newElement;

        newToOld[newElement] =
            oldElement;
      }


      //======================================================================
      // Rebuild the global dual graph using the new numbering
      //======================================================================

      std::vector<idx_t> xadjReordered(nelem + 1);

      xadjReordered[0] = 0;


      for(idx_t newElement = 0;
          newElement < nelem;
          ++newElement) {

        idx_t oldElement =
            newToOld[newElement];

        idx_t degree =
            xadjGlobal[oldElement + 1]
            -
            xadjGlobal[oldElement];

        xadjReordered[newElement + 1] =
            xadjReordered[newElement]
            +
            degree;
      }


      idx_t totalAdjacency =
          xadjReordered[nelem];

      std::vector<idx_t> adjncyReordered(totalAdjacency);


      for(idx_t newElement = 0;
          newElement < nelem;
          ++newElement) {

        idx_t oldElement =
            newToOld[newElement];

        idx_t destination =
            xadjReordered[newElement];


        for(idx_t jj = xadjGlobal[oldElement];
            jj < xadjGlobal[oldElement + 1];
            ++jj) {

          idx_t oldNeighbor =
              adjncyGlobal[jj];

          adjncyReordered[destination++] =
              oldToNew[oldNeighbor];
        }
      }


      //======================================================================
      // Extract this processor's local CSR graph.
      //
      // Because of the permutation above:
      //
      // processor p stores exactly inherited partition p.
      //======================================================================

      idx_t firstElement =
          vtxdist[iproc];

      idx_t lastElement =
          vtxdist[iproc + 1];

      idx_t nLocalElements =
          lastElement - firstElement;


      idx_t firstAdjacency =
          xadjReordered[firstElement];

      idx_t lastAdjacency =
          xadjReordered[lastElement];

      idx_t nLocalAdjacency =
          lastAdjacency - firstAdjacency;


      std::vector<idx_t> xadjLocal(
          nLocalElements + 1);

      std::vector<idx_t> adjncyLocal(
          nLocalAdjacency);


      for(idx_t i = 0;
          i <= nLocalElements;
          ++i) {

        xadjLocal[i] =
            xadjReordered[firstElement + i]
            -
            firstAdjacency;
      }


      for(idx_t i = 0;
          i < nLocalAdjacency;
          ++i) {

        adjncyLocal[i] =
            adjncyReordered[firstAdjacency + i];
      }


      //======================================================================
      // Initial local partition.
      //
      // In COUPLED mode ParMETIS obtains the initial partition from
      // vtxdist / processor ownership.
      //
      // Set localPart = iproc anyway because it is also the output array.
      //======================================================================

      std::vector<idx_t> localPart(
          nLocalElements,
          static_cast<idx_t>(iproc)
      );


      //======================================================================
      // AdaptiveRepart parameters
      //======================================================================

      idx_t wgtflag = 0;
      numflag       = 0;
      idx_t ncon    = 1;


      std::vector<real_t> tpwgts(nparts);

      for(idx_t p = 0; p < nparts; ++p) {

        tpwgts[p] =
            static_cast<real_t>(1.0)
            /
            static_cast<real_t>(nparts);
      }


      //----------------------------------------------------------------------
      // Allow 10% load imbalance.
      //
      // We deliberately do NOT force an extremely tight balance because
      // preserving the inherited subdomains is important here.
      //----------------------------------------------------------------------

      real_t ubvec[1];

      ubvec[0] =
          static_cast<real_t>(1.10);


      //----------------------------------------------------------------------
      // Equal migration cost for all elements
      //----------------------------------------------------------------------

      std::vector<idx_t> vsize(
          nLocalElements,
          static_cast<idx_t>(1)
      );


      //----------------------------------------------------------------------
      // IMPORTANT:
      //
      // Small itr strongly favors keeping elements in their original
      // partition.
      //
      // itr = 1000 was much too aggressive for this purpose.
      //----------------------------------------------------------------------

      real_t itr =
          static_cast<real_t>(1.0);


      //----------------------------------------------------------------------
      // COUPLED:
      //
      // Initial partition is the processor distribution itself.
      //
      // Because we reordered the graph, processor p currently stores
      // exactly inherited partition p.
      //----------------------------------------------------------------------

      idx_t parmetisOptions[4];

      parmetisOptions[0] = 1;
      parmetisOptions[1] = 0;
      parmetisOptions[2] = 0;
      parmetisOptions[3] = PARMETIS_PSR_COUPLED;


      idx_t edgecut = 0;

      MPI_Comm comm =
          MPI_COMM_WORLD;


      //======================================================================
      // Adaptive repartitioning
      //======================================================================

      err = ParMETIS_V3_AdaptiveRepart(
          vtxdist.data(),
          xadjLocal.data(),
          adjncyLocal.data(),
          NULL,                       // vwgt
          vsize.data(),               // migration cost
          NULL,                       // adjwgt
          &wgtflag,
          &numflag,
          &ncon,
          &nparts,
          tpwgts.data(),
          ubvec,
          &itr,
          parmetisOptions,
          &edgecut,
          localPart.data(),
          &comm
      );


      if(err != METIS_OK) {

        std::cerr
            << "ParMETIS_V3_AdaptiveRepart failed with error "
            << err
            << std::endl;

        METIS_Free(xadjGlobal);
        METIS_Free(adjncyGlobal);

        exit(1);
      }


      //======================================================================
      // Count how many elements ParMETIS moved
      //======================================================================

      idx_t localMoved = 0;

      for(idx_t i = 0;
          i < nLocalElements;
          ++i) {

        if(localPart[i] !=
           static_cast<idx_t>(iproc)) {

          ++localMoved;
        }
      }


      idx_t globalMoved = 0;

      MPI_Allreduce(
          &localMoved,
          &globalMoved,
          1,
          IDX_T,
          MPI_SUM,
          MPI_COMM_WORLD
      );


      //======================================================================
      // Gather the repartitioned result in REORDERED numbering
      //======================================================================

      std::vector<int> recvCounts(_nprocs);
      std::vector<int> displacements(_nprocs);


      for(unsigned p = 0;
          p < _nprocs;
          ++p) {

        recvCounts[p] =
            static_cast<int>(
                vtxdist[p + 1]
                -
                vtxdist[p]
            );

        displacements[p] =
            static_cast<int>(
                vtxdist[p]
            );
      }


      std::vector<idx_t> reorderedPart(nelem);


      MPI_Allgatherv(
          localPart.data(),
          static_cast<int>(nLocalElements),
          IDX_T,
          reorderedPart.data(),
          recvCounts.data(),
          displacements.data(),
          IDX_T,
          MPI_COMM_WORLD
      );


      //======================================================================
      // Convert partition back to original FEMuS element numbering
      //======================================================================

      for(idx_t oldElement = 0;
          oldElement < nelem;
          ++oldElement) {

        idx_t newElement =
            oldToNew[oldElement];

        epart[oldElement] =
            reorderedPart[newElement];
      }


      //======================================================================
      // CONNECTIVITY CLEANUP
      //
      // ParMETIS does not provide a METIS_OPTION_CONTIG equivalent here.
      //
      // For every partition:
      //
      //   1. Find its connected components in the ORIGINAL dual graph.
      //   2. Keep its largest component.
      //   3. Reassign each smaller island to an adjacent partition.
      //
      // The target partition is chosen primarily by the number of shared
      // dual edges. Balance is used as a secondary criterion.
      //======================================================================

      std::vector<idx_t> partitionSize(
          nparts,
          static_cast<idx_t>(0)
      );


      for(idx_t iel = 0;
          iel < nelem;
          ++iel) {

        partitionSize[epart[iel]]++;
      }


      const double idealSize =
          static_cast<double>(nelem)
          /
          static_cast<double>(nparts);

      const idx_t preferredMaximum =
          static_cast<idx_t>(
              1.10 * idealSize + 0.999999
          );


      idx_t cleanupMoved = 0;


      for(idx_t p = 0;
          p < nparts;
          ++p) {

        //--------------------------------------------------------------------
        // Find all connected components belonging to partition p
        //--------------------------------------------------------------------

        std::vector<char> visited(
            nelem,
            static_cast<char>(0)
        );

        std::vector< std::vector<idx_t> > components;


        for(idx_t start = 0;
            start < nelem;
            ++start) {

          if(epart[start] != p ||
             visited[start]) {

            continue;
          }


          std::vector<idx_t> component;
          std::vector<idx_t> stack;

          stack.push_back(start);
          visited[start] = 1;


          while(!stack.empty()) {

            idx_t current =
                stack.back();

            stack.pop_back();

            component.push_back(current);


            for(idx_t jj = xadjGlobal[current];
                jj < xadjGlobal[current + 1];
                ++jj) {

              idx_t neighbor =
                  adjncyGlobal[jj];

              if(epart[neighbor] == p &&
                 !visited[neighbor]) {

                visited[neighbor] = 1;
                stack.push_back(neighbor);
              }
            }
          }


          components.push_back(component);
        }


        //--------------------------------------------------------------------
        // Already connected
        //--------------------------------------------------------------------

        if(components.size() <= 1) {
          continue;
        }


        //--------------------------------------------------------------------
        // Keep largest component as the core of partition p
        //--------------------------------------------------------------------

        unsigned largestComponent = 0;

        for(unsigned c = 1;
            c < components.size();
            ++c) {

          if(components[c].size() >
             components[largestComponent].size()) {

            largestComponent = c;
          }
        }


        //--------------------------------------------------------------------
        // Every other component is an island.
        //--------------------------------------------------------------------

        for(unsigned c = 0;
            c < components.size();
            ++c) {

          if(c == largestComponent) {
            continue;
          }


          const std::vector<idx_t>& island =
              components[c];


          //------------------------------------------------------------------
          // Count how many dual edges connect this island to every other
          // partition.
          //------------------------------------------------------------------

          std::vector<idx_t> boundaryEdges(
              nparts,
              static_cast<idx_t>(0)
          );


          for(unsigned ii = 0;
              ii < island.size();
              ++ii) {

            idx_t element =
                island[ii];


            for(idx_t jj = xadjGlobal[element];
                jj < xadjGlobal[element + 1];
                ++jj) {

              idx_t neighbor =
                  adjncyGlobal[jj];

              idx_t q =
                  epart[neighbor];


              if(q != p) {

                boundaryEdges[q]++;
              }
            }
          }


          //------------------------------------------------------------------
          // First try to find a target that:
          //
          //   - shares a face/edge with the island
          //   - does not exceed preferredMaximum
          //
          // Among them choose the one with largest shared boundary.
          //------------------------------------------------------------------

          idx_t target = -1;
          idx_t bestBoundary = -1;


          for(idx_t q = 0;
              q < nparts;
              ++q) {

            if(q == p ||
               boundaryEdges[q] == 0) {

              continue;
            }


            idx_t resultingSize =
                partitionSize[q]
                +
                static_cast<idx_t>(island.size());


            if(resultingSize <= preferredMaximum) {

              if(target < 0 ||
                 boundaryEdges[q] > bestBoundary ||
                 (boundaryEdges[q] == bestBoundary &&
                  partitionSize[q] < partitionSize[target])) {

                target = q;
                bestBoundary = boundaryEdges[q];
              }
            }
          }


          //------------------------------------------------------------------
          // If no balanced target exists, connectivity has priority.
          //
          // Choose the neighboring partition with the largest shared
          // boundary.
          //------------------------------------------------------------------

          if(target < 0) {

            bestBoundary = -1;


            for(idx_t q = 0;
                q < nparts;
                ++q) {

              if(q == p ||
                 boundaryEdges[q] == 0) {

                continue;
              }


              if(target < 0 ||
                 boundaryEdges[q] > bestBoundary ||
                 (boundaryEdges[q] == bestBoundary &&
                  partitionSize[q] < partitionSize[target])) {

                target = q;
                bestBoundary = boundaryEdges[q];
              }
            }
          }


          //------------------------------------------------------------------
          // A disconnected global dual graph could theoretically give an
          // island with no adjacent partition.
          //------------------------------------------------------------------

          if(target < 0) {

            if(iproc == 0) {

              std::cerr
                  << "Warning: unable to reconnect component of partition "
                  << p
                  << " containing "
                  << island.size()
                  << " elements."
                  << std::endl;
            }

            continue;
          }


          //------------------------------------------------------------------
          // Move the complete connected island.
          //------------------------------------------------------------------

          for(unsigned ii = 0;
              ii < island.size();
              ++ii) {

            idx_t element =
                island[ii];

            epart[element] =
                target;

            ++cleanupMoved;
          }


          partitionSize[p] -=
              static_cast<idx_t>(island.size());

          partitionSize[target] +=
              static_cast<idx_t>(island.size());
        }
      }


      //======================================================================
      // Verify connectivity after cleanup
      //======================================================================

      bool connectivityOK = true;


      for(idx_t p = 0;
          p < nparts;
          ++p) {

        idx_t first = -1;

        for(idx_t iel = 0;
            iel < nelem;
            ++iel) {

          if(epart[iel] == p) {

            first = iel;
            break;
          }
        }


        if(first < 0) {

          if(iproc == 0) {

            std::cerr
                << "Warning: partition "
                << p
                << " is empty after connectivity cleanup."
                << std::endl;
          }

          continue;
        }


        std::vector<char> visited(
            nelem,
            static_cast<char>(0)
        );

        std::vector<idx_t> stack;

        stack.push_back(first);
        visited[first] = 1;

        idx_t count = 0;


        while(!stack.empty()) {

          idx_t current =
              stack.back();

          stack.pop_back();

          ++count;


          for(idx_t jj = xadjGlobal[current];
              jj < xadjGlobal[current + 1];
              ++jj) {

            idx_t neighbor =
                adjncyGlobal[jj];

            if(epart[neighbor] == p &&
               !visited[neighbor]) {

              visited[neighbor] = 1;
              stack.push_back(neighbor);
            }
          }
        }


        if(count != partitionSize[p]) {

          connectivityOK = false;

          if(iproc == 0) {

            std::cerr
                << "WARNING: partition "
                << p
                << " remains disconnected. "
                << count
                << " / "
                << partitionSize[p]
                << " elements are connected to its main component."
                << std::endl;
          }
        }
      }


      //======================================================================
      // Diagnostics
      //======================================================================

      if(iproc == 0) {

        std::cout << std::endl;

        std::cout
            << "ParMETIS adaptive repartitioning:"
            << std::endl;

        std::cout
            << "  edgecut            = "
            << edgecut
            << std::endl;

        std::cout
            << "  ParMETIS moved     = "
            << globalMoved
            << " / "
            << nelem
            << " elements"
            << std::endl;

        std::cout
            << "  cleanup moved      = "
            << cleanupMoved
            << " elements"
            << std::endl;

        std::cout
            << "  connected          = "
            << (connectivityOK ? "YES" : "NO")
            << std::endl;

        std::cout
            << "Final partition:"
            << std::endl;


        for(idx_t p = 0;
            p < nparts;
            ++p) {

          std::cout
              << "  partition "
              << p
              << " : "
              << partitionSize[p]
              << " elements"
              << std::endl;
        }

        std::cout << std::endl;
      }


      //======================================================================
      // METIS allocated these arrays
      //======================================================================

      METIS_Free(xadjGlobal);
      METIS_Free(adjncyGlobal);

#endif // HAVE_PARMETIS
    }


    //========================================================================
    // NO INITIAL PARTITION:
    //
    // Original METIS behavior.
    //========================================================================

    else {

      std::vector<idx_t> npart(nnodes);

      idx_t objval;

      idx_t options[METIS_NOPTIONS];

      METIS_SetDefaultOptions(options);

      options[METIS_OPTION_NUMBERING] = 0;
      options[METIS_OPTION_DBGLVL]    = 0;
      options[METIS_OPTION_CTYPE]     = METIS_CTYPE_SHEM;
      options[METIS_OPTION_PTYPE]     = METIS_PTYPE_KWAY;
      options[METIS_OPTION_IPTYPE]    = METIS_IPTYPE_RANDOM;
      options[METIS_OPTION_CONTIG]    = 0;
      options[METIS_OPTION_MINCONN]   = 1;
      options[METIS_OPTION_NITER]     = 10;
      options[METIS_OPTION_UFACTOR]   = 100;


      idx_t nparts =
          static_cast<idx_t>(_nprocs);


      int err = METIS_PartMeshDual(
          &nelem,
          &nnodes,
          eptr.data(),
          eind.data(),
          NULL,
          NULL,
          &ncommon,
          &nparts,
          NULL,
          options,
          &objval,
          epart.data(),
          npart.data()
      );


      if(err == METIS_OK) {

        std::cout
            << " METIS PARTITIONING IS OK "
            << std::endl;
      }

      else if(err == METIS_ERROR_INPUT) {

        std::cerr
            << " METIS_ERROR_INPUT "
            << std::endl;

        exit(1);
      }

      else if(err == METIS_ERROR_MEMORY) {

        std::cerr
            << " METIS_ERROR_MEMORY "
            << std::endl;

        exit(2);
      }

      else {

        std::cerr
            << " METIS_GENERIC_ERROR "
            << std::endl;

        exit(3);
      }
    }


    std::vector<unsigned> MaterialElementCounter =
        _mesh.el->GetMaterialElementCounter();

#endif // HAVE_METIS
  }


  //========================================================================
  // Return partition to FEMuS
  //========================================================================

  partition.resize(
      static_cast<unsigned>(epart.size())
  );


  for(unsigned i = 0;
      i < static_cast<unsigned>(epart.size());
      ++i) {

    partition[i] =
        static_cast<unsigned>(epart[i]);
  }
}


//------------------------------------------------------------------------------------------------------
//   void MeshMetisPartitioning::DoPartition(std::vector <unsigned>& partition, const bool& AMR) {
//
//     int nnodes = _mesh.GetNumberOfNodes();
//     int nelem = _mesh.GetNumberOfElements();
//
//     std::vector <int> epart(nelem);
//
//     if(_nprocs == 1) {
//       //serial computation
//       for(unsigned i = 0; i < nelem; i++) {
//         epart[i] = 0;
//       }
//     }
//     else if(_nprocs > nelem) {
//       std::cout << "Error In MeshMetis::DoPartition, the number of processes " << _nprocs
//                 << " is greater than the number of elements " << nelem << std::endl;
//       abort();
//     }
//     else if(_nprocs == nelem) {
//       for(unsigned i = 0; i < nelem; i++) {
//         epart[i] = i;
//       }
//     }
//     else {
//
// #ifndef HAVE_METIS
//       std::cerr << "Fatal error: Metis library was not found. Metis partioning algorithm cannot be called!" << std::endl;
//       exit(1);
// #endif
//
//       unsigned eind_size = _mesh.el->GetElementNumber("Hex") * NVE[0][2]      + _mesh.el->GetElementNumber("Tet") * NVE[1][2]
//                            + _mesh.el->GetElementNumber("Wedge") * NVE[2][2]    + _mesh.el->GetElementNumber("Quad") * NVE[3][2]
//                            + _mesh.el->GetElementNumber("Triangle") * NVE[4][2] + _mesh.el->GetElementNumber("Line") * NVE[5][2];
//
//       vector < idx_t > eptr(nelem + 1);
//       vector < idx_t > eind(eind_size);
//
//       vector < int > npart(nnodes);
//
//       idx_t objval;
//       idx_t options[METIS_NOPTIONS];
//
//       METIS_SetDefaultOptions(options);
//
//       options[METIS_OPTION_NUMBERING] = 0;
//       options[METIS_OPTION_DBGLVL]   = 0;
//       options[METIS_OPTION_CTYPE]    = METIS_CTYPE_SHEM;
//       options[METIS_OPTION_PTYPE]    = METIS_PTYPE_KWAY;
//       options[METIS_OPTION_IPTYPE]   = METIS_IPTYPE_RANDOM;
//       options[METIS_OPTION_CONTIG]   = 0;
//       options[METIS_OPTION_MINCONN]  = 1;
//       options[METIS_OPTION_NITER]    = 10;
//       options[METIS_OPTION_UFACTOR]  = 100;
//
//       eptr[0] = 0;
//       unsigned counter = 0;
//       for(unsigned iel = 0; iel < nelem; iel++) {
//         unsigned ndofs = _mesh.el->GetElementDofNumber(iel, 2);
//         eptr[iel + 1] = eptr[iel] + ndofs;
//         for(unsigned inode = 0; inode < ndofs; inode++) {
//           eind[counter] = _mesh.el->GetElementDofIndex(iel, inode);
//           counter++;
//         }
//       }
//
//
//       int ncommon = (AMR || _mesh.GetDimension() == 1) ? 1 : _mesh.GetDimension() + 1;
//
//       //I call the Mesh partioning function of Metis library (output is epart(own elem) and npart (own nodes))
//       int err = METIS_PartMeshDual(&nelem, &nnodes, &eptr[0], &eind[0], NULL, NULL, &ncommon, &_nprocs, NULL, options, &objval, &epart[0], &npart[0]);
//
//       if(err == METIS_OK) {
//         std::cout << " METIS PARTITIONING IS OK " << std::endl;
//       }
//       else if(err == METIS_ERROR_INPUT) {
//         cout << " METIS_ERROR_INPUT " << endl;
//         exit(1);
//       }
//       else if(err == METIS_ERROR_MEMORY) {
//         cout << " METIS_ERROR_MEMORY " << endl;
//         exit(2);
//       }
//       else {
//         cout << " METIS_GENERIC_ERROR " << endl;
//         exit(3);
//       }
//
//       std::vector<unsigned> MaterialElementCounter = _mesh.el->GetMaterialElementCounter();
//
// //         for (unsigned i = 0 ; i < MaterialElementCounter.size(); i++)   std::cout << MaterialElementCounter[i] << " ";
// //         std::cout << std::endl;
//
//
//     }
//
//     partition.resize(epart.size());
//     for(unsigned i = 0; i < epart.size(); i++) {
//       partition[i] = epart[i];
//     }
//
//     return;
//   }






















  void MeshMetisPartitioning::DoPartition(std::vector <unsigned>& partition, const Mesh& meshc) {
    partition.resize(_mesh.GetNumberOfElements());
    unsigned refIndex = _mesh.GetRefIndex();
    for(int isdom = 0; isdom < _nprocs; isdom++) {
      for(unsigned iel = meshc._elementOffset[isdom]; iel < meshc._elementOffset[isdom + 1]; iel++) {
        for(unsigned j = 0; j < refIndex; j++) {
          partition[ iel * refIndex + j ] = isdom;
        }
      }
    }
  }

}

