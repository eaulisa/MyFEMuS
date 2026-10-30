/*=========================================================================

 Program: FEMUS
 Module: Elem
 Authors: Eugenio Aulisa, Sara Calandrini, Giacomo Capodaglio

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <assert.h>
#include <cfloat>

#include "Elem.hpp"
#include "GeomElTypeEnum.hpp"

#include <unordered_map>
#include <unordered_set>

#include "MeshRefinement.hpp"

namespace femus {

  using std::cout;
  using std::endl;

  /**
   * This constructor allocates the memory for the \textit{coarsest elem}
   **/
  elem::elem(const unsigned& other_nel) {
    _coarseElem = NULL;

    _level = 0;

    _nelt[0] = _nelt[1] = _nelt[2] = _nelt[3] = _nelt[4] = _nelt[5] = 0;
    _nel = other_nel;

    _elementType.resize(_nel);
    _elementGroup.resize(_nel);
    _elementMaterial.resize(_nel);
    _elementLevel.resize(_nel, _level);

    _elementDof.resize(_nel, NVE[0][2], UINT_MAX);
    _elementNearFace.resize(_nel, NFC[0][1], -1);

  }

  void elem::ShrinkToFit() {

    _elementDof.shrinkToFit(UINT_MAX);

    MyVector <unsigned> rowSize(_nel);
    for (unsigned iel = 0; iel < _nel; iel++) {
      unsigned ielType = GetElementType(iel);
      rowSize[iel] = NFC[ielType][1];
    }
    _elementNearFace.shrinkToFit(rowSize);

  }
  /**
   * This constructor allocates the memory for the \textit{finer elem}
   * starting from the parameters of the \textit{coarser elem}
   **/
  elem::elem(elem* elc, const unsigned refindex, const std::vector < double >& coarseAmrVector) {
    _coarseElem = elc;

    _level = elc->_level + 1;

    _nelt[0] = _nelt[1] = _nelt[2] = _nelt[3] = _nelt[4] = _nelt[5] = 0;
    _nel = elc->GetRefinedElementNumber() * refindex; //refined
    _nel += elc->GetElementNumber() - elc->GetRefinedElementNumber(); // + non-refined;

    _elementType.resize(_nel);
    _elementGroup.resize(_nel);
    _elementMaterial.resize(_nel);
    _elementLevel.resize(_nel, _level);

    //**************************
    MyVector <unsigned> rowSizeElDof(_nel);
    MyVector <unsigned> rowSizeElNearFace(_nel);
    unsigned jel = 0;
    for (unsigned isdom = 0; isdom < elc->_nprocs; isdom++) {
      elc->_elementType.broadcast(isdom);
      for (unsigned iel = elc->_elementType.begin(); iel < elc->_elementType.end(); iel++) {
        short unsigned elType = elc->_elementType[iel];
        int increment = 1;
        if (static_cast < short unsigned >(coarseAmrVector[iel] + 0.25) == 1) {
          increment = NRE[elType];
        }
        for (unsigned j = 0; j < increment; j++) {
          rowSizeElDof[jel + j] += NVE[elType][2];
          rowSizeElNearFace[jel + j] += NFC[elType][1];
        }
        jel += increment;
      }
      elc->_elementType.clearBroadcast();
    }
    _elementDof = MyMatrix <unsigned> (rowSizeElDof);
    _elementNearFace = MyMatrix <int> (rowSizeElNearFace, -1);

    rowSizeElDof.clear();
    rowSizeElNearFace.clear();

  }

  elem::elem(Mesh* mshc) {

    elem* elc = mshc->el;
    const unsigned refindex = mshc->GetRefIndex();
    auto amr = mshc->_topology->_Sol[mshc->GetAmrIndex()];

    _coarseElem = elc;
    _iproc = elc->_iproc;
    _nprocs = elc->_nprocs;
    _level = elc->_level + 1;
    _ngroup = elc->_ngroup;
    _nelr = 0;
    _nvt = 0;

    unsigned count_send = 0;
    unsigned nelt_send[6] = {0, 0, 0, 0, 0, 0};
    unsigned material_send[3] = {0, 0, 0};

    for (unsigned iel = elc->_elementType.begin(); iel < elc->_elementType.end(); iel++) {
      short unsigned elType = elc->_elementType[iel];
      short unsigned elMaterial = elc->_elementMaterial[iel];
      unsigned increment = (static_cast<short unsigned>((*amr)(iel) + 0.25) == 1) ? refindex : 1;
      count_send += increment;
      nelt_send[elType] += increment;
      if (elMaterial == 2) material_send[0] += increment;
      else if (elMaterial == 3) material_send[1] += increment;
      else material_send[2] += increment;
    }

    _elementOwned = count_send;
    _elementOffset.assign(_nprocs + 1, 0);
    _materialElementCounter.assign(3, 0);

    MPI_Allgather(&count_send, 1, MPI_UNSIGNED, _elementOffset.data() + 1, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    MPI_Allreduce(nelt_send, _nelt, 6, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(material_send, _materialElementCounter.data(), 3, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

    for (unsigned p = 0; p < _nprocs; p++) _elementOffset[p + 1] += _elementOffset[p];

    _nel = _elementOffset[_nprocs];

    _elementType.resize(_elementOffset);
    _elementGroup.resize(_elementOffset);
    _elementMaterial.resize(_elementOffset);
    _elementLevel.resize(_elementOffset);

    MyVector<unsigned> rowSizeElDof(_elementOffset);
    MyVector<unsigned> rowSizeElNearFace(_elementOffset);

    unsigned jel = rowSizeElDof.begin();
    for (unsigned iel = elc->_elementType.begin(); iel < elc->_elementType.end(); iel++) {
      bool refined = (static_cast<short unsigned>((*amr)(iel) + 0.25) == 1);
      short unsigned elType = elc->_elementType[iel];
      short unsigned elGroup = elc->_elementGroup[iel];
      short unsigned elMaterial = elc->_elementMaterial[iel];
      short unsigned elLevel = refined ? static_cast<short unsigned>(elc->_elementLevel[iel] + 1) : elc->_elementLevel[iel];
      unsigned increment = refined ? refindex : 1;

      for (unsigned j = 0; j < increment; j++) {
        _elementType[jel + j] = elType;
        _elementGroup[jel + j] = elGroup;
        _elementMaterial[jel + j] = elMaterial;
        _elementLevel[jel + j] = elLevel;
        rowSizeElDof[jel + j] = NVE[elType][2];
        rowSizeElNearFace[jel + j] = NFC[elType][1];
        elc->SetChildElement(iel, j, jel + j);
      }

      jel += increment;
    }

    _elementDof = MyMatrix<unsigned>(rowSizeElDof, UINT_MAX);
    _elementNearFace = MyMatrix<int>(rowSizeElNearFace, -1);

    jel = rowSizeElDof.begin();
    for (unsigned iel = elc->_elementType.begin(); iel < elc->_elementType.end(); iel++) {
      bool refined = (static_cast<short unsigned>((*amr)(iel) + 0.25) == 1);
      short unsigned elType = elc->_elementType[iel];
      unsigned increment = refined ? refindex : 1;

      if (refined) {
        for (unsigned j = 0; j < increment; j++) {
          for (unsigned inode = 0; inode < NVE[elType][0]; inode++) {
            unsigned jDof = mshc->_finiteElement[elType][0]->GetBasis()->GetFine2CoarseVertexMapping(j, inode);
            _elementDof[jel + j][inode] = elc->_elementDof[iel][jDof];
          }
        }
      }
      else {
        for (unsigned inode = 0; inode < NVE[elType][2]; inode++) _elementDof[jel][inode] = elc->_elementDof[iel][inode];
      }

      for (unsigned iface = 0; iface < NFC[elType][1]; iface++) {
        int value = elc->_elementNearFace[iel][iface];
        if (value < -1) {
          if (refined) {
            for (unsigned jface = 0; jface < refindex / 2; jface++) {
              unsigned child = coarse2FineFaceMapping[elType][iface][jface][0];
              unsigned face = coarse2FineFaceMapping[elType][iface][jface][1];
              _elementNearFace[jel + child][face] = value;
            }
          }
          else {
            _elementNearFace[jel][iface] = value;
          }
        }
      }

      jel += increment;
    }

    rowSizeElDof.clear();
    rowSizeElNearFace.clear();

    SetNodeNumber(elc->GetNodeNumber());
  }

  // elem::elem(elem* elc, const unsigned, const MyVector<double>& coarseAmrVector) {
  //   _coarseElem = elc;
  //   _iproc = elc->_iproc;
  //   _nprocs = elc->_nprocs;
  //   _level = elc->_level + 1;
  //   _ngroup = elc->_ngroup;
  //   _nelr = 0;
  //   _nvt = 0;
  //
  //   unsigned count_send = 0;
  //   unsigned nelt_send[6] = {0, 0, 0, 0, 0, 0};
  //
  //   for (unsigned iel = elc->_elementType.begin(); iel < elc->_elementType.end(); iel++) {
  //     short unsigned elType = elc->_elementType[iel];
  //     unsigned increment = (static_cast<short unsigned>(coarseAmrVector[iel] + 0.25) == 1) ? NRE[elType] : 1;
  //     count_send += increment;
  //     nelt_send[elType] += increment;
  //   }
  //
  //   _elementOwned = count_send;
  //   _elementOffset.assign(_nprocs + 1, 0);
  //
  //   MPI_Allgather(&count_send, 1, MPI_UNSIGNED, _elementOffset.data() + 1, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
  //   MPI_Allreduce(nelt_send, _nelt, 6, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  //
  //   for (unsigned p = 0; p < _nprocs; p++) _elementOffset[p + 1] += _elementOffset[p];
  //
  //   _nel = _elementOffset[_nprocs];
  //
  //   _elementType.resize(_elementOffset);
  //   _elementGroup.resize(_elementOffset);
  //   _elementMaterial.resize(_elementOffset);
  //   _elementLevel.resize(_elementOffset);
  //
  //   MyVector<unsigned> rowSizeElDof(_elementOffset);
  //   MyVector<unsigned> rowSizeElNearFace(_elementOffset);
  //
  //   unsigned jel = rowSizeElDof.begin();
  //   for (unsigned iel = elc->_elementType.begin(); iel < elc->_elementType.end(); iel++) {
  //     bool refined = (static_cast<short unsigned>(coarseAmrVector[iel] + 0.25) == 1);
  //     short unsigned elType = elc->_elementType[iel];
  //     short unsigned elGroup = elc->_elementGroup[iel];
  //     short unsigned elMaterial = elc->_elementMaterial[iel];
  //     short unsigned elLevel = refined ? static_cast<short unsigned>(_level) : elc->_elementLevel[iel];
  //     unsigned increment = refined ? NRE[elType] : 1;
  //
  //     for (unsigned j = 0; j < increment; j++) {
  //       _elementType[jel + j] = elType;
  //       _elementGroup[jel + j] = elGroup;
  //       _elementMaterial[jel + j] = elMaterial;
  //       _elementLevel[jel + j] = elLevel;
  //       rowSizeElDof[jel + j] = NVE[elType][2];
  //       rowSizeElNearFace[jel + j] = NFC[elType][1];
  //     }
  //
  //     jel += increment;
  //   }
  //
  //   _elementDof = MyMatrix<unsigned>(rowSizeElDof);
  //   _elementNearFace = MyMatrix<int>(rowSizeElNearFace, -1);
  //
  //   jel = rowSizeElDof.begin();
  //   for (unsigned iel = elc->_elementType.begin(); iel < elc->_elementType.end(); iel++) {
  //     bool refined = (static_cast<short unsigned>(coarseAmrVector[iel] + 0.25) == 1);
  //     short unsigned elType = elc->_elementType[iel];
  //
  //     if (!refined) {
  //       for (unsigned iface = 0; iface < NFC[elType][1]; iface++) {
  //         int face = elc->_elementNearFace[iel][iface];
  //         if (face < -1) _elementNearFace[jel][iface] = face;
  //       }
  //     }
  //
  //     jel += refined ? NRE[elType] : 1;
  //   }
  //
  //   rowSizeElDof.clear();
  //   rowSizeElNearFace.clear();
  // }

  /*
   elem::elem(elem* elc, const unsigned, const MyVector<double>& coarseAmrVector) {
     _coarseElem = elc;
     _level = elc->_level + 1;

     _nelt[0] = _nelt[1] = _nelt[2] = _nelt[3] = _nelt[4] = _nelt[5] = 0;

     unsigned count_send = 0;
     for (unsigned iel = elc->_elementType.begin(); iel < elc->_elementType.end(); iel++) {
       short unsigned elType = elc->_elementType[iel];
       count_send += (static_cast<short unsigned>(coarseAmrVector[iel] + 0.25) == 1) ? NRE[elType] : 1;
     }

     _elementOwned = count_send;

     _elementOffset.assign(elc->_nprocs + 1, 0);
     MPI_Allgather(&count_send, 1, MPI_UNSIGNED, _elementOffset.data() + 1, 1, MPI_UNSIGNED, MPI_COMM_WORLD);

     for (unsigned p = 0; p < elc->_nprocs; p++) _elementOffset[p + 1] += _elementOffset[p];

     _nel = _elementOffset[elc->_nprocs];

     _elementType.resize(_elementOffset);
     _elementGroup.resize(_elementOffset);
     _elementMaterial.resize(_elementOffset);
     _elementLevel.resize(_elementOffset, _level);

     MyVector<unsigned> rowSizeElDof(_elementOffset);
     MyVector<unsigned> rowSizeElNearFace(_elementOffset);

     unsigned jel = rowSizeElDof.begin();
     for (unsigned iel = elc->_elementType.begin(); iel < elc->_elementType.end(); iel++) {
       short unsigned elType = elc->_elementType[iel];
       unsigned increment = 1;
       if (static_cast<short unsigned>(coarseAmrVector[iel] + 0.25) == 1) increment = NRE[elType];

       for (unsigned j = 0; j < increment; j++) {
         rowSizeElDof[jel + j] = NVE[elType][2];
         rowSizeElNearFace[jel + j] = NFC[elType][1];
       }

       jel += increment;
     }

     _elementDof = MyMatrix<unsigned>(rowSizeElDof);
     _elementNearFace = MyMatrix<int>(rowSizeElNearFace, -1);

     rowSizeElDof.clear();
     rowSizeElNearFace.clear();
   }*/

  void elem::ReorderMeshElements(const std::vector < unsigned >& elementMapping) {

    //BEGIN reordering _elementType
    MyVector <short unsigned> tempElementType;
    tempElementType = _elementType;
    for (unsigned iel = 0; iel < _nel; iel++) {
      _elementType[ elementMapping [iel] ]  = tempElementType[iel] ;
    }
    tempElementType.clear();
    //END reordering _elementType

    //BEGIN reordering _elementLevel
    if (_level != 0) {
      MyVector <short unsigned> tempElementLevel;
      tempElementLevel = _elementLevel;
      for (unsigned iel = 0; iel < _nel; iel++) {
        _elementLevel[elementMapping [iel]] = tempElementLevel[iel];
      }
      tempElementLevel.clear();
    }
    //END reordering _elementCanBeRefined

    //BEGIN reordering _elementGroup and _elementMaterial

    MyVector <short unsigned> tempElementGroup;
    MyVector <short unsigned> tempElementMaterial;
    tempElementGroup = _elementGroup;
    tempElementMaterial = _elementMaterial;

    for (unsigned iel = 0; iel < _nel; iel++) {
      _elementGroup[elementMapping[iel]]   = tempElementGroup[iel];
      _elementMaterial[elementMapping[iel]] = tempElementMaterial[iel] ;
    }
    tempElementGroup.clear();
    tempElementMaterial.clear();

    //END reordering _elementGroup and _elementMaterial

    //BEGIN reordering _elementDof (rows)
    MyVector <unsigned> rowSize(_nel, 0);
    for (unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      rowSize[elementMapping[i]] = _elementDof.size(i);
    }

    MyMatrix <unsigned> tmpElDof = _elementDof;
    _elementDof = MyMatrix <unsigned>(rowSize, 0);
    for (unsigned i = tmpElDof.begin(); i < tmpElDof.end(); i++) {
      for (unsigned j = tmpElDof.begin(i); j < tmpElDof.end(i); j++) {
        _elementDof[elementMapping [i]][j] = tmpElDof[i][j];
      }
    }
    tmpElDof.clear();
    //END reordering OF _elementDof

    //BEGIN reordering _elementNearFace (rows)
    for (unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      rowSize[elementMapping[i]] = _elementNearFace.size(i);
    }
    MyMatrix <int> tmpElNearFace = _elementNearFace;
    _elementNearFace = MyMatrix <int>(rowSize, 0);
    for (unsigned i = tmpElNearFace.begin(); i < tmpElNearFace.end(); i++) {
      for (unsigned j = tmpElNearFace.begin(i); j < tmpElNearFace.end(i); j++) {
        _elementNearFace[elementMapping [i]][j] = tmpElNearFace[i][j];
      }
    }
    tmpElNearFace.clear();
    //END reordering _elementNearFace

    //BEGIN reordering _childElementDof (columns) on coarse level
    if (_level != 0) {
      for (unsigned i = _coarseElem->_childElem.begin(); i < _coarseElem->_childElem.end(); i++) {
        for (unsigned j = _coarseElem->_childElem.begin(i); j < _coarseElem->_childElem.end(i); j++) {
          _coarseElem->_childElem[i][j] =  elementMapping[ _coarseElem->_childElem[i][j]];
        }
      }
    }
    //END reordering _childElementDof
  }

  void elem::ReorderMeshNodes(const std::vector < unsigned >& nodeMapping) {
    for (unsigned i = _elementDof.begin(); i < _elementDof.end(); i++) {
      for (unsigned j = _elementDof.begin(i); j < _elementDof.end(i); j++) {
        _elementDof[i][j] =  nodeMapping[_elementDof[i][j]];
      }
    }
  }

  elem::~elem() {
  }

  void elem::DeleteElementNearVertex() {
    _elementNearVertex.clear();
  }

  /**
   * Return the number of vertices(type=0) + midpoints(type=1) + facepoints(type=2) + interiorpoits(type=2)
   **/
  unsigned elem::GetElementDofNumber(const unsigned& iel, const unsigned& type) {
    return NVE[_elementType[iel]][type];
  }

  /**
   * Set the local->global node number
   **/
  void elem::SetElementDofIndex(const unsigned& iel, const unsigned& inode, const unsigned& value) {
    _elementDof[iel][inode] = value;
  }

  /**
   * Return the local->global face node index
   **/
  unsigned elem::GetFaceVertexIndex(const unsigned& iel, const unsigned& iface, const unsigned& inode) {
    return _elementDof[iel][ig[_elementType[iel]][iface][inode]];
  }

  /**
   * Return the total node number
   **/
  unsigned elem::GetNodeNumber()const {
    return _nvt;
  }

  /**
   * Set the total node number
   **/
  void elem::SetNodeNumber(const unsigned& value) {
    _nvt = value;
  }

  /**
   * Return the total number of elements
   **/
  unsigned elem::GetElementNumber(const char* name) const {
    if (!strcmp(name, "All")) {
      return _nel;
    }
    unsigned i;
    i = this->GetIndex(name);
    return _nelt[i];
  }

  /**
   * Add value to the total number of the element
   **/
  void elem::AddToElementNumber(const unsigned& value, const char name[]) {
    unsigned i;
    i = GetIndex(name);
    _nelt[i] += value;
  }
  void elem::AddToElementNumber(const unsigned& value, short unsigned ielt) {
    _nelt[ielt] += value;
  }
  void elem::BroadcastElementNumber() {
    MPI_Bcast(_nelt, 6, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  }
  unsigned elem::GetElementFaceNumber(const unsigned& iel, const unsigned& type) {
    return NFC[ _elementType[iel] ][type];
  }

  unsigned elem::GetFaceRangeStart(const unsigned &ielt) const {
    return FACERANGE[ielt][0];
  }

  unsigned elem::GetFaceRangeEnd(const unsigned &ielt) const {
    return FACERANGE[ielt][1];
  }

  /**
   * Return the global adiacent-to-face element number + (-1 is and AMR face) and boundary indexes -2,-3,-4,... -> 1, 2, 3,...
   **/
  int elem::GetFaceElementIndex(const unsigned& iel, const unsigned& iface) {
    return _elementNearFace[iel][iface];
  }

  int elem::GetBoundaryIndex(const unsigned& iel, const unsigned& iface) {
    return  -(GetFaceElementIndex(iel, iface) + 1);
  }

  /**
   * Set the global adiacent-to-face element number
   **/
  void elem::SetFaceElementIndex(const unsigned& iel, const unsigned& iface, const int& value) {
    _elementNearFace[iel][iface] = value;
  }

  /**
   * Return element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
   **/
  short unsigned elem::GetElementType(const unsigned& iel) {
    return _elementType[iel];
  }

  /**
   * Set element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
   **/
  void elem::SetElementType(const unsigned& iel, const short unsigned& value) {
    _elementType[iel] = value;
  }

  /**
   * Return element group
   **/
  short unsigned elem::GetElementGroup(const unsigned& iel) {
    return _elementGroup[iel];
  }

  /**
   * Set element group
   **/
  void elem::SetElementGroup(const unsigned& iel, const short unsigned& value) {
    _elementGroup[iel] = value;
  }

  /**
   * Set element Material
  **/
  void elem::SetElementMaterial(const unsigned& iel, const short unsigned& value) {
    _elementMaterial[iel] = value;
  }

  /**
   * Return element material
   **/
  short unsigned elem::GetElementMaterial(const unsigned& iel) {
    return _elementMaterial[iel];
  }

  /**
   * Return element group number
   **/
  unsigned elem::GetElementGroupNumber() const {
    return _ngroup;
  }

  /**
   * Set element group
   **/
  void elem::SetElementGroupNumber(const unsigned& value) {
    _ngroup = value;
  }

  /**
   * Set the memory storage and initialize nve and kvtel (node->element vectors)
   **/

  void elem::BuildElementNearElement() {
    MyVector < unsigned > rowSize(_elementOffset, 1);
    for (unsigned iel = rowSize.begin(); iel < rowSize.end(); iel++) {
      std::map< unsigned, bool> elements;
      for (unsigned i = 0; i < GetElementDofNumber(iel, 0); i++) {
        unsigned inode = GetElementDofIndex(iel, i);
        for (unsigned j = _elementNearVertex.begin(inode); j < _elementNearVertex.end(inode); j++) {
          elements[_elementNearVertex[inode][j]] = true;
        }
      }
      rowSize[iel] = elements.size();
    }
    _elementNearElement = MyMatrix <unsigned> (rowSize, UINT_MAX);
    for (unsigned iel = _elementNearElement.begin(); iel < _elementNearElement.end(); iel++) {
      std::map< unsigned, bool> elements;
      _elementNearElement[iel][0] = iel;
      for (unsigned i = 0; i < GetElementDofNumber(iel, 0); i++) {
        unsigned inode = GetElementDofIndex(iel, i);
        for (unsigned j = _elementNearVertex.begin(inode); j < _elementNearVertex.end(inode); j++) {
          if (_elementNearVertex[inode][j] != iel) {
            elements[_elementNearVertex[inode][j]] = true;
          }
        }
      }
      unsigned j = 1;
      for (std::map<unsigned, bool>::iterator it = elements.begin(); it != elements.end(); it++, j++) {
        _elementNearElement[iel][j] = it->first;
      }
    }

  }

// #include <mpi.h>
// #include <algorithm>
// #include <climits>
// #include <unordered_map>
// #include <vector>

  // void ExchangeNodeBuffers(const std::vector<std::vector<unsigned>>& buffers_send, std::vector<std::vector<unsigned>>& buffers_recv, const unsigned iproc, const unsigned nprocs, const int tag) {
  //   std::vector<int> counts_send(nprocs, 0);
  //   std::vector<int> counts_recv(nprocs, 0);
  //
  //   for (unsigned jproc = 0; jproc < nprocs; jproc++) {
  //     if (jproc == iproc) continue;
  //     if (buffers_send[jproc].size() > INT_MAX) MPI_Abort(MPI_COMM_WORLD, 1);
  //     counts_send[jproc] = static_cast<int>(buffers_send[jproc].size());
  //   }
  //
  //   MPI_Alltoall(counts_send.data(), 1, MPI_INT, counts_recv.data(), 1, MPI_INT, MPI_COMM_WORLD);
  //
  //   buffers_recv.clear();
  //   buffers_recv.resize(nprocs);
  //   buffers_recv[iproc] = buffers_send[iproc];
  //
  //   std::vector<MPI_Request> requests(2 * nprocs, MPI_REQUEST_NULL);
  //
  //   for (unsigned jproc = 0; jproc < nprocs; jproc++) {
  //     if (jproc == iproc) continue;
  //
  //     if (counts_recv[jproc] > 0) {
  //       buffers_recv[jproc].resize(counts_recv[jproc]);
  //       MPI_Irecv(buffers_recv[jproc].data(), counts_recv[jproc], MPI_UNSIGNED, static_cast<int>(jproc), tag, MPI_COMM_WORLD, &requests[jproc]);
  //     }
  //
  //     if (counts_send[jproc] > 0) {
  //       MPI_Isend(buffers_send[jproc].data(), counts_send[jproc], MPI_UNSIGNED, static_cast<int>(jproc), tag, MPI_COMM_WORLD, &requests[nprocs + jproc]);
  //     }
  //   }
  //
  //   MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
  // }
  //
  // void elem::BuildElementNearElement(Mesh *msh) {
  //   unsigned dim = msh->GetDimension();
  //   unsigned nprocs = msh->n_processors();
  //   const unsigned iproc = msh->processor_id();
  //
  //   MyVector<unsigned> rowSize(_elementOffset, 1);
  //
  //   std::unordered_map<unsigned, std::vector<unsigned>> nodeToElement;
  //   nodeToElement.reserve(rowSize.end() - rowSize.begin()); //TODO get this infor from msh
  //
  //   std::vector<std::vector<unsigned>> ghostNodes_send(nprocs);
  //
  //   for (unsigned iel = rowSize.begin(); iel < rowSize.end(); iel++) {
  //     for (unsigned i = 0; i < GetElementDofNumber(iel, 0); i++) {
  //       unsigned inode = GetElementDofIndex(iel, i);
  //       std::vector<unsigned>& elements = nodeToElement[inode];
  //       if (elements.empty()) {
  //         elements.reserve(2u << dim);
  //         unsigned jproc = msh->IsdomBisectionSearch(inode, 2);
  //
  //         if (jproc != iproc) ghostNodes_send[jproc].push_back(inode);
  //       }
  //       elements.push_back(iel);
  //     }
  //   }
  //
  //   std::vector<std::vector<unsigned>> ghost_send(nprocs);
  //   std::vector<std::vector<unsigned>> ghost_recv;
  //
  //   for (unsigned jproc = 0; jproc < nprocs; jproc++) {
  //     unsigned count_send = 0;
  //
  //     for (unsigned i = 0; i < ghostNodes_send[jproc].size(); i++) {
  //       unsigned inode = ghostNodes_send[jproc][i];
  //       const std::vector<unsigned>& elements = nodeToElement.at(inode);
  //       if (count_send > static_cast<unsigned>(INT_MAX - 2) || elements.size() > static_cast<unsigned>(INT_MAX - 2) - count_send) MPI_Abort(MPI_COMM_WORLD, 1);
  //       count_send += 2 + static_cast<unsigned>(elements.size());
  //     }
  //
  //     ghost_send[jproc].reserve(count_send);
  //
  //     for (unsigned i = 0; i < ghostNodes_send[jproc].size(); i++) {
  //       unsigned inode = ghostNodes_send[jproc][i];
  //       const std::vector<unsigned>& elements = nodeToElement.at(inode);
  //
  //       ghost_send[jproc].push_back(inode);
  //       ghost_send[jproc].push_back(static_cast<unsigned>(elements.size()));
  //       ghost_send[jproc].insert(ghost_send[jproc].end(), elements.begin(), elements.end());
  //     }
  //   }
  //
  //   ExchangeNodeBuffers(ghost_send, ghost_recv, iproc, nprocs, 100);
  //
  //   for (unsigned jproc = 0; jproc < nprocs; jproc++) {
  //     unsigned pos_recv = 0;
  //
  //     while (pos_recv < ghost_recv[jproc].size()) {
  //       unsigned inode = ghost_recv[jproc][pos_recv++];
  //       unsigned count_recv = ghost_recv[jproc][pos_recv++];
  //       std::vector<unsigned>& elements = nodeToElement[inode];
  //
  //       elements.insert(elements.end(), ghost_recv[jproc].begin() + pos_recv, ghost_recv[jproc].begin() + pos_recv + count_recv);
  //       pos_recv += count_recv;
  //     }
  //   }
  //
  //   for (std::unordered_map<unsigned, std::vector<unsigned>>::iterator it = nodeToElement.begin(); it != nodeToElement.end(); ++it) {
  //     std::vector<unsigned>& elements = it->second;
  //     std::sort(elements.begin(), elements.end());
  //     elements.erase(std::unique(elements.begin(), elements.end()), elements.end());
  //   }
  //
  //   std::vector<std::vector<unsigned>> updated_send(nprocs);
  //   std::vector<std::vector<unsigned>> updated_recv;
  //
  //   for (unsigned jproc = 0; jproc < nprocs; jproc++) {
  //     unsigned pos_recv = 0;
  //     unsigned count_send = 0;
  //
  //     while (pos_recv < ghost_recv[jproc].size()) {
  //       unsigned inode = ghost_recv[jproc][pos_recv++];
  //       unsigned count_recv = ghost_recv[jproc][pos_recv++];
  //       pos_recv += count_recv;
  //
  //       const std::vector<unsigned>& elements = nodeToElement.at(inode);
  //       if (count_send > static_cast<unsigned>(INT_MAX - 2) || elements.size() > static_cast<unsigned>(INT_MAX - 2) - count_send) MPI_Abort(MPI_COMM_WORLD, 1);
  //       count_send += 2 + static_cast<unsigned>(elements.size());
  //     }
  //
  //     updated_send[jproc].reserve(count_send);
  //     pos_recv = 0;
  //
  //     while (pos_recv < ghost_recv[jproc].size()) {
  //       unsigned inode = ghost_recv[jproc][pos_recv++];
  //       unsigned count_recv = ghost_recv[jproc][pos_recv++];
  //       pos_recv += count_recv;
  //
  //       const std::vector<unsigned>& elements = nodeToElement.at(inode);
  //
  //       updated_send[jproc].push_back(inode);
  //       updated_send[jproc].push_back(static_cast<unsigned>(elements.size()));
  //       updated_send[jproc].insert(updated_send[jproc].end(), elements.begin(), elements.end());
  //     }
  //   }
  //
  //   ExchangeNodeBuffers(updated_send, updated_recv, iproc, nprocs, 101);
  //
  //   for (unsigned jproc = 0; jproc < nprocs; jproc++) {
  //     unsigned pos_recv = 0;
  //
  //     while (pos_recv < updated_recv[jproc].size()) {
  //       unsigned inode = updated_recv[jproc][pos_recv++];
  //       unsigned count_recv = updated_recv[jproc][pos_recv++];
  //       std::vector<unsigned>& elements = nodeToElement.at(inode);
  //
  //       elements.assign(updated_recv[jproc].begin() + pos_recv, updated_recv[jproc].begin() + pos_recv + count_recv);
  //       pos_recv += count_recv;
  //     }
  //   }
  //
  //   std::vector<std::vector<unsigned>> elementNeighbors(rowSize.end() - rowSize.begin());
  //
  //   for (unsigned iel = rowSize.begin(); iel < rowSize.end(); iel++) {
  //     std::vector<unsigned>& elements = elementNeighbors[iel - rowSize.begin()];
  //     elements.push_back(iel);
  //
  //     for (unsigned i = 0; i < GetElementDofNumber(iel, 0); i++) {
  //       unsigned inode = GetElementDofIndex(iel, i);
  //       const std::vector<unsigned>& adjacentElements = nodeToElement.at(inode);
  //
  //       for (unsigned j = 0; j < adjacentElements.size(); j++) {
  //         if (adjacentElements[j] != iel) elements.push_back(adjacentElements[j]);
  //       }
  //     }
  //
  //     std::sort(elements.begin() + 1, elements.end());
  //     elements.erase(std::unique(elements.begin() + 1, elements.end()), elements.end());
  //     rowSize[iel] = static_cast<unsigned>(elements.size());
  //   }
  //
  //   _elementNearElement = MyMatrix<unsigned>(rowSize, UINT_MAX);
  //
  //   for (unsigned iel = _elementNearElement.begin(); iel < _elementNearElement.end(); iel++) {
  //     const std::vector<unsigned>& elements = elementNeighbors[iel - rowSize.begin()];
  //
  //     for (unsigned j = 0; j < elements.size(); j++) {
  //       _elementNearElement[iel][j] = elements[j];
  //     }
  //   }
  //
  //   for (unsigned iel = rowSize.begin(); iel < rowSize.end(); iel++) {
  //     for (unsigned iface = 0; iface < GetElementFaceNumber(iel); iface++) {
  //       if ( GetFaceElementIndex(iel, iface) != -1 )continue;
  //
  //       const std::vector<unsigned>* faceElements[3] = {nullptr, nullptr, nullptr};
  //       unsigned first = 0;
  //
  //       for (unsigned k = 0; k < dim; k++) {
  //         unsigned inode = GetFaceVertexIndex(iel, iface, k);
  //         faceElements[k] = &nodeToElement.at(inode);
  //         if (faceElements[k]->size() < faceElements[first]->size()) first = k;
  //       }
  //
  //       const std::vector<unsigned>& candidates = *faceElements[first];
  //
  //       for (unsigned j = 0; j < candidates.size(); j++) {
  //         unsigned jel = candidates[j];
  //         if (jel == iel) continue;
  //
  //         bool shared = true;
  //
  //         for (unsigned k = 0; k < dim && shared; k++) {
  //           if (k != first) shared = std::binary_search(faceElements[k]->begin(), faceElements[k]->end(), jel);
  //         }
  //
  //         if (shared) {
  //           SetFaceElementIndex(iel, iface, jel + 1u);
  //           break;
  //         }
  //       }
  //
  //       if (_level == 0 && GetFaceElementIndex(iel, iface) == -1) {
  //         SetFaceElementIndex(iel, iface, -101);
  //       }
  //     }
  //   }
  //
  // }

  void elem::BuildElementNearVertex() {
    MyVector <unsigned> rowSize(_nvt, 0);
    for (unsigned iel = 0; iel < _nel; iel++) {
      for (unsigned inode = 0; inode < GetElementDofNumber(iel, 0); inode++) {
        rowSize[GetElementDofIndex(iel, inode)]++;
      }
    }
    _elementNearVertex = MyMatrix <unsigned> (rowSize, _nel);
    for (unsigned iel = 0; iel < _nel; iel++) {
      for (unsigned inode = 0; inode < GetElementDofNumber(iel, 0); inode++) {
        unsigned irow = GetElementDofIndex(iel, inode);
        unsigned j = 0;
        while (_nel != _elementNearVertex[irow][j]) j++;
        _elementNearVertex[irow][j] = iel;
      }
    }
  }

  /**
   * Return the number of elements which have a given vertex
   **/
  unsigned elem::GetElementNearVertexNumber(const unsigned & inode) {
    return _elementNearVertex.size(inode);
  }

  /**
   * Return the element index for the given i-node in the j-position with 0<=j<nve(i)
   **/
  unsigned elem::GetElementNearVertex(const unsigned & inode, const unsigned & j) {
    return _elementNearVertex[inode][j];
  }

  /**
   * return the index 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
   **/
  unsigned elem::GetIndex(const char name[]) const {
    unsigned index = 0;
    if (!strcmp(name, "Hex")) {
      index = 0;
    }
    else if (!strcmp(name, "Tet")) {
      index = 1;
    }
    else if (!strcmp(name, "Wedge")) {
      index = 2;
    }
    else if (!strcmp(name, "Quad")) {
      index = 3;
    }
    else if (!strcmp(name, "Triangle")) {
      index = 4;
    }
    else if (!strcmp(name, "Line")) {
      index = 5;
    }
    else {
      cout << "error! invalid Element Shape in elem::GetIndex(...)" << endl;
      exit(0);
    }
    return index;
  }

  void elem::AllocateChildrenElement(const unsigned & refindex, Mesh * msh) {
    MyVector <unsigned> rowSize(_elementOffset, 0);
    for (unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      rowSize[i] = (msh->GetRefinedElementIndex(i) == 1) ? refindex : 1;
    }
    _childElem = MyMatrix <unsigned> (rowSize, 0);

    for (unsigned i = rowSize.begin(); i < rowSize.end(); i++) {
      unsigned elementType = msh->GetElementType(i);
      rowSize[i] = (msh->GetRefinedElementIndex(i) == 1) ? refindex * NVE[elementType][2] : NVE[elementType][2];
    }
    _childElemDof = MyMatrix <unsigned> (rowSize, 0);
  }

  void elem::SetChildElementDof(elem * elf) {

    for (unsigned iel = _childElemDof.begin(); iel < _childElemDof.end(); iel++) {
      // unsigned nDofs = GetElementDofNumber(iel, 2);
      unsigned nDofs = _childElemDof.size(iel) / _childElem.size(iel);
      for (unsigned j = _childElemDof.begin(iel); j < _childElemDof.end(iel); j++) {
        unsigned ielf = GetChildElement(iel, j / nDofs);
        _childElemDof[iel][j] = elf->GetElementDofIndex(ielf, j % nDofs);
      }
    }
  }

  void elem::SetChildElementDof(Mesh* meshf) {
    const int tag_request = 200;
    const int tag_values = 201;

    std::vector<unsigned> requests_send;
    std::vector<unsigned> values_recv;
    unsigned count_recv = 0;

    if (_iproc != 0) {
      for (unsigned iel = _childElemDof.begin(); iel < _childElemDof.end(); iel++) {
        unsigned nDofs = GetElementDofNumber(iel, 2);
        unsigned rowDofs = _childElemDof.end(iel) - _childElemDof.begin(iel);

        if (rowDofs > static_cast<unsigned>(INT_MAX) - count_recv) MPI_Abort(MPI_COMM_WORLD, 1);
        count_recv += rowDofs;

        for (unsigned j = _childElemDof.begin(iel); j < _childElemDof.end(iel); j += nDofs) {
          unsigned ielf = GetChildElement(iel, j / nDofs);
          requests_send.push_back(ielf);
          requests_send.push_back(nDofs);
        }
      }

      values_recv.resize(count_recv);
    }

    if (requests_send.size() > INT_MAX) MPI_Abort(MPI_COMM_WORLD, 1);

    int count_send = static_cast<int>(requests_send.size());
    std::vector<int> counts_recv(_nprocs, 0);

    MPI_Gather(&count_send, 1, MPI_INT, counts_recv.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (_iproc == 0) {
      elem* elf = meshf->el;

      std::vector<std::vector<unsigned>> requests_recv(_nprocs);
      std::vector<std::vector<unsigned>> values_send(_nprocs);
      std::vector<MPI_Request> communication_recv(_nprocs, MPI_REQUEST_NULL);
      std::vector<MPI_Request> communication_send(_nprocs, MPI_REQUEST_NULL);

      for (unsigned jproc = 1; jproc < _nprocs; jproc++) {
        if (counts_recv[jproc] > 0) {
          requests_recv[jproc].resize(counts_recv[jproc]);
          MPI_Irecv(requests_recv[jproc].data(), counts_recv[jproc], MPI_UNSIGNED, static_cast<int>(jproc), tag_request, MPI_COMM_WORLD, &communication_recv[jproc]);
        }
      }

      for (unsigned iel = _childElemDof.begin(); iel < _childElemDof.end(); iel++) {
        unsigned nDofs = GetElementDofNumber(iel, 2);

        for (unsigned j = _childElemDof.begin(iel); j < _childElemDof.end(iel); j++) {
          unsigned ielf = GetChildElement(iel, j / nDofs);
          _childElemDof[iel][j] = elf->GetElementDofIndex(ielf, j % nDofs);
        }
      }

      MPI_Waitall(static_cast<int>(_nprocs), communication_recv.data(), MPI_STATUSES_IGNORE);

      for (unsigned jproc = 1; jproc < _nprocs; jproc++) {
        if (counts_recv[jproc] == 0) continue;

        unsigned count_send = 0;

        for (unsigned j = 0; j < requests_recv[jproc].size(); j += 2) {
          count_send += requests_recv[jproc][j + 1];
        }

        values_send[jproc].resize(count_send);
        unsigned pos_send = 0;

        for (unsigned j = 0; j < requests_recv[jproc].size(); j += 2) {
          unsigned ielf = requests_recv[jproc][j];
          unsigned nDofs = requests_recv[jproc][j + 1];

          for (unsigned k = 0; k < nDofs; k++) {
            values_send[jproc][pos_send++] = elf->GetElementDofIndex(ielf, k);
          }
        }

        MPI_Isend(values_send[jproc].data(), static_cast<int>(count_send), MPI_UNSIGNED, static_cast<int>(jproc), tag_values, MPI_COMM_WORLD, &communication_send[jproc]);
      }

      MPI_Waitall(static_cast<int>(_nprocs), communication_send.data(), MPI_STATUSES_IGNORE);
    }
    else {
      if (count_send > 0) {
        MPI_Request communication[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

        MPI_Irecv(values_recv.data(), static_cast<int>(count_recv), MPI_UNSIGNED, 0, tag_values, MPI_COMM_WORLD, &communication[0]);
        MPI_Isend(requests_send.data(), count_send, MPI_UNSIGNED, 0, tag_request, MPI_COMM_WORLD, &communication[1]);

        MPI_Waitall(2, communication, MPI_STATUSES_IGNORE);
      }

      unsigned pos_recv = 0;

      for (unsigned iel = _childElemDof.begin(); iel < _childElemDof.end(); iel++) {
        for (unsigned j = _childElemDof.begin(iel); j < _childElemDof.end(iel); j++) {
          _childElemDof[iel][j] = values_recv[pos_recv++];
        }
      }
    }
  }

  unsigned elem::GetChildElementDof(const unsigned & iel, const unsigned & i0, const unsigned i1) {
    return _childElemDof[iel][i0 * GetElementDofNumber(iel, 2) + i1];
  }

  void elem::SetChildElement(const unsigned & iel, const unsigned & json, const unsigned & value) {
    _childElem[iel][json] = value;
  }

  unsigned elem::GetChildElement(const unsigned & iel, const unsigned & json) {
    return _childElem[iel][json];
  }

  const unsigned elem::GetNVE(const unsigned & elementType, const unsigned & doftype) const {
    return NVE[elementType][doftype];
  }

  const unsigned elem::GetNFACENODES(const unsigned & elementType, const unsigned & jface, const unsigned & dof) const {
    return NFACENODES[elementType][jface][dof];
  }

  const unsigned elem::GetNFC(const unsigned & elementType, const unsigned & type) const {
    return NFC[elementType][type];
  }

  const unsigned elem::GetIG(const unsigned & elementType, const unsigned & iface, const unsigned & jnode) const {
    return ig[elementType][iface][jnode];
  }

  unsigned elem::GetElementDofIndex(const unsigned & iel, const unsigned & inode) {
    return _elementDof[iel][inode];
  };

  const unsigned FELT[6][2] = {{3, 3}, {4, 4}, {3, 4}, {5, 5}, {5, 5}, {6, 6}};
  unsigned elem::GetFaceType(const unsigned & ielt, const unsigned & jface) {
    return FELT[ielt][jface >= NFC[ielt][0]];
  }

  void elem::GetAMRRestriction(Mesh * msh) {

    std::vector < std::map < unsigned,  std::map < unsigned, double  > > > &restriction = msh->GetAmrRestrictionMap();
    restriction.resize(3);

    std::vector < std::map < unsigned, bool > > &interfaceSolidMark = msh->GetAmrSolidMark();
    interfaceSolidMark.resize(3);

    std::vector < MyVector<unsigned> > interfaceElement; // iel = interfaceElement[ilevel][i]
    std::vector < MyMatrix<unsigned> > interfaceLocalDof; // ldof = interfaceElement[ilevel][i][j]
    std::vector < std::vector < MyMatrix<unsigned> > > interfaceDof; // gdof = interfaceDof[solType][ilevel][i][j]
    std::vector < std::vector < MyMatrix<unsigned> > > levelInterfaceSolidMark; // FSIdof = levelInterfaceSolidMark[solType][ilevel][i][j]
    std::vector < std::vector < MyMatrix< double > > > interfaceNodeCoordinates; // x[dim] = interfaceNodeCoordinates[ilevel][dim][i][j]

    interfaceElement.resize(_level + 1);
    interfaceLocalDof.resize(_level + 1);
    interfaceDof.resize(3);
    levelInterfaceSolidMark.resize(3);
    for (unsigned i = 0; i < 3; i++) {
      interfaceDof[i].resize(_level + 1);
      levelInterfaceSolidMark[i].resize(_level + 1);
    }
    interfaceNodeCoordinates.resize(_level + 1);
    unsigned dim = msh->GetDimension();

    const double __amr_total_start = MPI_Wtime();
    double __amr_t0 = 0.0;

    double __amr_t_interface_build = 0.0;
    double __amr_t_bounding_boxes  = 0.0;

    double __amr_t_jmin[3]                = {0.0, 0.0, 0.0};
    double __amr_t_candidate_exchange[3]  = {0.0, 0.0, 0.0};
    double __amr_t_restriction_scan[3]    = {0.0, 0.0, 0.0};
    double __amr_t_bbx_lookup[3]          = {0.0, 0.0, 0.0};
    double __amr_t_candidate_loop[3]      = {0.0, 0.0, 0.0};
    double __amr_t_same_node[3]           = {0.0, 0.0, 0.0};
    double __amr_t_geometry_filter[3]     = {0.0, 0.0, 0.0};
    double __amr_t_projection_inverse[3]  = {0.0, 0.0, 0.0};
    double __amr_t_eval_store[3]          = {0.0, 0.0, 0.0};
    double __amr_t_parallel_closure[3]    = {0.0, 0.0, 0.0};
    double __amr_t_genealogy[3]           = {0.0, 0.0, 0.0};
    double __amr_t_solidmark_bcast[3]     = {0.0, 0.0, 0.0};

    unsigned long long __amr_cnt_received_points[3]      = {0ull, 0ull, 0ull};
    unsigned long long __amr_cnt_candidate_elements[3]    = {0ull, 0ull, 0ull};
    unsigned long long __amr_cnt_same_node_skip[3]        = {0ull, 0ull, 0ull};
    unsigned long long __amr_cnt_geometry_pass[3]         = {0ull, 0ull, 0ull};
    unsigned long long __amr_cnt_inverse_calls[3]         = {0ull, 0ull, 0ull};
    unsigned long long __amr_cnt_inside_domain[3]         = {0ull, 0ull, 0ull};
    unsigned long long __amr_cnt_stored_entries[3]        = {0ull, 0ull, 0ull};

    __amr_t0 = MPI_Wtime();

    for (unsigned ilevel = 0; ilevel <= _level; ilevel++) {
      //BEGIN interface element search
      interfaceElement[ilevel] = MyVector <unsigned> (_elementOwned);
      unsigned counter = 0;
      for (unsigned i = _elementLevel.begin(); i < _elementLevel.end(); i++) {
        if (ilevel == _elementLevel[i]) {
          for (unsigned j = _elementNearFace.begin(i); j < _elementNearFace.end(i); j++) {
            if (-1 == _elementNearFace[i][j]) {
              interfaceElement[ilevel][counter] = i;
              counter++;
              break;
            }
          }
        }
      }
      interfaceElement[ilevel].resize(counter);
      interfaceElement[ilevel].stack();
      //END interface element search

      //BEGIN interface node search
      std::vector< unsigned > offset = interfaceElement[ilevel].getOffset();
      interfaceLocalDof[ilevel] = MyMatrix <unsigned>(offset, NVE[0][2], UINT_MAX);
      for (unsigned i = interfaceElement[ilevel].begin(); i < interfaceElement[ilevel].end(); i++) {
        unsigned iel =  interfaceElement[ilevel][i];
        std::map <unsigned, bool> ldofs;
        for (unsigned jface = _elementNearFace.begin(iel); jface < _elementNearFace.end(iel); jface++) {
          if (-1 == _elementNearFace[iel][jface]) {
            for (unsigned k = 0; k < GetNFACENODES(GetElementType(iel), jface, 2); k++) {
              unsigned index = GetIG(GetElementType(iel), jface, k);
              ldofs[index] = true;
            }
          }
        }
        unsigned j = 0;
        for (std::map<unsigned, bool>::iterator it = ldofs.begin(); it != ldofs.end(); it++) {
          interfaceLocalDof[ilevel][i][j] = it->first;
          j++;
        }
      }
      interfaceLocalDof[ilevel].shrinkToFit(UINT_MAX);
      //END interface node search

      //BEGIN interface node dof global search, one for each soltype
      MyVector <unsigned> rowSize = interfaceLocalDof[ilevel].getRowSize();
      for (unsigned soltype = 0; soltype < 3; soltype++) {
        interfaceDof[soltype][ilevel] = MyMatrix< unsigned > (rowSize, UINT_MAX);
        levelInterfaceSolidMark[soltype][ilevel] = MyMatrix< unsigned > (rowSize, UINT_MAX);
        for (unsigned i = interfaceLocalDof[ilevel].begin(); i < interfaceLocalDof[ilevel].end(); i++) {
          unsigned iel = interfaceElement[ilevel][i];
          unsigned counter = 0;
          for (unsigned j = interfaceLocalDof[ilevel].begin(i); j < interfaceLocalDof[ilevel].end(i); j++) {
            unsigned jloc = interfaceLocalDof[ilevel][i][j];
            if (jloc < GetElementDofNumber(iel, soltype)) {
              unsigned jdof  = msh->GetSolutionDof(jloc, iel, soltype);
              interfaceDof[soltype][ilevel][i][counter] = jdof;
              unsigned jdof2  = msh->GetSolutionDof(jloc, iel, 2);
              levelInterfaceSolidMark[soltype][ilevel][i][counter] = msh->GetSolidMark(jdof2);
              counter++;
            }
            else {
              break;
            }
          }
        }
        interfaceDof[soltype][ilevel].shrinkToFit(UINT_MAX);
        levelInterfaceSolidMark[soltype][ilevel].shrinkToFit(UINT_MAX);
      }
      //END interface node dof global search, one for each soltype

      //BEGIN interface node coordinates
      interfaceNodeCoordinates[ilevel].resize(dim);
      for (unsigned k = 0; k < dim; k++) {
        interfaceNodeCoordinates[ilevel][k] = MyMatrix< double > (rowSize, 0.);
      }
      for (unsigned i = interfaceLocalDof[ilevel].begin(); i < interfaceLocalDof[ilevel].end(); i++) {
        unsigned iel = interfaceElement[ilevel][i];
        for (unsigned j = interfaceLocalDof[ilevel].begin(i); j < interfaceLocalDof[ilevel].end(i); j++) {
          unsigned jnode = interfaceLocalDof[ilevel][i][j];
          unsigned xDof  = msh->GetSolutionDof(jnode, iel, 2);
          for (unsigned k = 0; k < dim; k++) {
            interfaceNodeCoordinates[ilevel][k][i][j] = (*msh->_topology->_Sol[k])(xDof);
          }
        }
      }
      //END interface node coordinates search
    }

    __amr_t_interface_build += MPI_Wtime() - __amr_t0;

    const unsigned nLevels = _level + 1u;

    std::vector<std::vector<double*> > xMin(_nprocs);
    std::vector<std::vector<double*> > xMax(_nprocs);

    std::vector<double> xMinMemory(_nprocs * nLevels * dim);
    std::vector<double> xMaxMemory(_nprocs * nLevels * dim);

    std::vector<std::vector<std::vector<unsigned>>> bbx_to_elements;
    std::vector<std::vector<unsigned>> bbxN;

    __amr_t0 = MPI_Wtime();
    GetAMRInterfaceBoundingBoxes(msh, interfaceElement, xMin, xMax, xMinMemory, xMaxMemory, bbx_to_elements, bbxN);
    __amr_t_bounding_boxes += MPI_Wtime() - __amr_t0;

    for (unsigned soltype = 0; soltype < 3; soltype++) {

      std::vector<double> xl(dim);
      std::vector<double> phi;
      std::vector<std::vector<double>> gradPhi(dim);
      std::vector < std::vector <double > > xv(dim);
      std::vector <double> xc(dim);
      std::vector < std::vector< double > > xe(dim);
      std::vector <double> xi(dim);
      std::vector < std::vector < std::vector <double > > > aP(3);

      std::vector<unsigned> jMinLevel(nLevels, UINT_MAX);
      std::vector<unsigned> jMaxLevel(nLevels, 0);
      std::vector<unsigned> jSizeLevel(nLevels, 0);

      __amr_t0 = MPI_Wtime();

      for (unsigned level = 0; level < nLevels; level++) {
        for (unsigned k = interfaceDof[soltype][level].begin(); k < interfaceDof[soltype][level].end(); k++) {
          for (unsigned l = interfaceDof[soltype][level].begin(k); l < interfaceDof[soltype][level].end(k); l++) {
            const unsigned ldof = interfaceDof[soltype][level][k][l];

            if (ldof < jMinLevel[level]) jMinLevel[level] = ldof;
            if (ldof > jMaxLevel[level]) jMaxLevel[level] = ldof;
          }
        }

        if (jMaxLevel[level] >= jMinLevel[level]) {
          jSizeLevel[level] = jMaxLevel[level] - jMinLevel[level] + 1u;
        }
      }

      __amr_t_jmin[soltype] += MPI_Wtime() - __amr_t0;

      std::vector<std::vector<unsigned>> jDof_r(_nprocs);
      std::vector<std::vector<unsigned>> jSolidMark_r(_nprocs);
      std::vector<std::vector<std::vector<double>>> jCoordinate_r(_nprocs);

      for(unsigned lproc = 0u; lproc < _nprocs; lproc++) {
        jCoordinate_r[lproc].resize(dim);
      }

      for (int ilevel = 0; ilevel < _level; ilevel++) {
        for (int jlevel = ilevel + 1; jlevel <= _level; jlevel++) {

          const unsigned jMin  = jMinLevel[jlevel];
          const unsigned jSize = jSizeLevel[jlevel];

          __amr_t0 = MPI_Wtime();

          GetAMRRestrictionCandidateDataForLevelPair(
            static_cast<unsigned>(ilevel),                // ilevel
            dim,                                          // dim
            jMin,                                         // jMin
            jSize,                                        // jSize
            interfaceDof[soltype][jlevel],                // interfaceDofJ
            levelInterfaceSolidMark[soltype][jlevel],     // levelInterfaceSolidMarkJ
            interfaceNodeCoordinates[jlevel],             // interfaceNodeCoordinatesJ
            xMin,                                         // xMin
            xMax,                                         // xMax
            jDof_r,                                       // output
            jSolidMark_r,                                 // output
            jCoordinate_r                                 // output
          );

          __amr_t_candidate_exchange[soltype] += MPI_Wtime() - __amr_t0;

          __amr_t0 = MPI_Wtime();

          for (unsigned lproc = 0; lproc < _nprocs; lproc++) {

            for (unsigned kl = 0; kl < jDof_r[lproc].size(); kl++) {

              __amr_cnt_received_points[soltype]++;

              const unsigned ldof = jDof_r[lproc][kl];

              for (unsigned d = 0; d < dim; d++) {
                xl[d] = jCoordinate_r[lproc][d][kl];
              }

              unsigned bbxIdx = 0u;

              const double __amr_t_bbx0 = MPI_Wtime();

              const bool pointInsideGrid =
                GetAMRBBoxIndex(ilevel,
                                dim,
                                xl,
                                xMin,
                                xMax,
                                bbxN,
                                bbxIdx);

              __amr_t_bbx_lookup[soltype] += MPI_Wtime() - __amr_t_bbx0;

              if (!pointInsideGrid) continue;
              if (bbxIdx >= bbx_to_elements[ilevel].size()) continue;

              const std::vector<unsigned>& candidateElements =
                bbx_to_elements[ilevel][bbxIdx];

              const double __amr_t_cand0 = MPI_Wtime();

              for (unsigned ee = 0; ee < candidateElements.size(); ee++) {

                __amr_cnt_candidate_elements[soltype]++;

                const unsigned iel = interfaceElement[ilevel][candidateElements[ee]];

                const unsigned i = candidateElements[ee];

                const double __amr_t_same0 = MPI_Wtime();

                bool sameNode = false;

                const unsigned nElementDofs = GetElementDofNumber(iel, soltype);

                for (unsigned j = 0; j < nElementDofs; j++) {
                  const unsigned jdof = msh->GetSolutionDof(j, iel, soltype);

                  if (jdof == ldof) {
                    sameNode = true;
                    break;
                  }
                }

                __amr_t_same_node[soltype] += MPI_Wtime() - __amr_t_same0;

                if (sameNode) {
                  __amr_cnt_same_node_skip[soltype]++;
                  continue;
                }

                const double __amr_t_geom0 = MPI_Wtime();

                short unsigned ielType = _elementType[iel];
                basis* base = msh->GetBasis(ielType, soltype);

                msh->GetElementNodeCoordinates(xv, iel);

                double r2;
                GetConvexHullSphereRadiousSquare(xv, xc, r2, 0.01);

                double d2 = 0.0;

                for (unsigned d = 0; d < dim; d++) {
                  d2 += (xl[d] - xc[d]) * (xl[d] - xc[d]);
                }

                if (d2 > r2) {
                  __amr_t_geometry_filter[soltype] += MPI_Wtime() - __amr_t_geom0;
                  continue;
                }

                GetBoundingBox(xv, xe, 0.01);

                bool insideHull = true;

                for (unsigned d = 0; d < dim; d++) {
                  if (xl[d] < xe[d][0] || xl[d] > xe[d][1]) {
                    insideHull = false;
                    break;
                  }
                }

                __amr_t_geometry_filter[soltype] += MPI_Wtime() - __amr_t_geom0;

                if (!insideHull) continue;

                __amr_cnt_geometry_pass[soltype]++;

                const double __amr_t_inv0 = MPI_Wtime();

                for (unsigned jtype = 0; jtype < 3; jtype++) {
                  ProjectNodalToPolynomialCoefficients(aP[jtype], xv, ielType, jtype);
                }

                GetClosestPointInReferenceElement(xv, xl, ielType, xi);

                GetInverseMapping_fast(2u,
                                       ielType,
                                       aP,
                                       xl,
                                       xi,
                                       100u,
                                       phi,
                                       gradPhi);

                __amr_t_projection_inverse[soltype] += MPI_Wtime() - __amr_t_inv0;
                __amr_cnt_inverse_calls[soltype]++;

                const bool insideDomain =
                  CheckIfPointIsInsideReferenceDomain(xi, ielType, 0.0001);

                if (!insideDomain) continue;

                __amr_cnt_inside_domain[soltype]++;

                const double __amr_t_store0 = MPI_Wtime();

                for (unsigned j = interfaceDof[soltype][ilevel].begin(i);
                     j < interfaceDof[soltype][ilevel].end(i);
                     j++) {

                  const unsigned jloc = interfaceLocalDof[ilevel][i][j];

                  const double value = base->eval_phi(jloc, xi);

                  if (fabs(value) >= 1.0e-10) {
                    const unsigned jdof = interfaceDof[soltype][ilevel][i][j];

                    auto& rowJ = restriction[soltype][jdof];

                    auto it = rowJ.find(jdof);

                    if (it == rowJ.end()) {
                      rowJ.emplace(jdof, 1.0);
                      interfaceSolidMark[soltype][jdof] =
                        levelInterfaceSolidMark[soltype][ilevel][i][j];
                    }

                    rowJ[ldof] = value;

                    auto& rowL = restriction[soltype][ldof];
                    rowL[ldof] = 10.0;

                    interfaceSolidMark[soltype][ldof] = jSolidMark_r[lproc][kl];
                    __amr_cnt_stored_entries[soltype]++;
                  }
                }

                __amr_t_eval_store[soltype] += MPI_Wtime() - __amr_t_store0;
              }

              __amr_t_candidate_loop[soltype] += MPI_Wtime() - __amr_t_cand0;
            }
          }
        }
      }

      __amr_t_restriction_scan[soltype] += MPI_Wtime() - __amr_t0;

      // __amr_t0 = MPI_Wtime();
      //
      // unsigned counter = 1;
      // while (counter != 0) {
      //   counter = 0;
      //
      //   //BEGIN  saving the restriction object in parallel vectors and matrices
      //
      //   MyVector <unsigned> rowSize(restriction[soltype].size(), 0);
      //   unsigned cnt1 = 0;
      //   for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
      //     rowSize[cnt1] = restriction[soltype][it1->first].size();
      //     cnt1++;
      //   }
      //   rowSize.stack();
      //
      //   std::vector< unsigned > offset = rowSize.getOffset();
      //
      //   MyVector <unsigned> masterNode(offset);
      //   MyMatrix <unsigned> slaveNodes(rowSize);
      //   MyMatrix <double> slaveNodesValues(rowSize);
      //
      //   cnt1 = 0;
      //   for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
      //     masterNode[offset[_iproc] + cnt1] = it1->first;
      //     unsigned cnt2 = 0;
      //     for (std::map<unsigned, double> ::iterator it2 = restriction[soltype][it1->first].begin(); it2 != restriction[soltype][it1->first].end(); it2++) {
      //       slaveNodes[offset[_iproc] + cnt1][cnt2] = it2->first;
      //       slaveNodesValues[offset[_iproc] + cnt1][cnt2] = it2->second;
      //       cnt2++;
      //     }
      //     cnt1++;
      //   }
      //   //END  saving the restriction object in parallel vectors and matrices
      //
      //
      //   //BEGIN filling the restriction object with infos coming form the parallel vectors and matrices
      //   unsigned solutionOffset = msh->_dofOffset[soltype][_iproc];
      //   unsigned solutionOffsetp1 = msh->_dofOffset[soltype][_iproc + 1];
      //   for (unsigned lproc = 0; lproc < _nprocs; lproc++) {
      //     masterNode.broadcast(lproc);
      //     slaveNodes.broadcast(lproc);
      //     slaveNodesValues.broadcast(lproc);
      //     for (unsigned i = slaveNodes.begin(); i < slaveNodes.end(); i++) {
      //       unsigned inode = masterNode[i];
      //       if (inode >= solutionOffset && inode < solutionOffsetp1 && // inode belongs to _iproc
      //           restriction[soltype].find(inode) == restriction[soltype].end()) { // but inode is not set as master node of _iproc
      //         counter++;
      //         for (unsigned j = slaveNodes.begin(i); j < slaveNodes.end(i); j++) { //copy information for lproc to _iproc
      //           unsigned jnode = slaveNodes[i][j];
      //           restriction[soltype][inode][jnode] = slaveNodesValues[i][j];
      //         }
      //       }
      //       else { // either inode does not belong to _iproc or it was already defined as a master for _iproc
      //         if (restriction[soltype].find(inode) != restriction[soltype].end()) { // inode is already defined as master node in _iproc (either it does or does not belong to _iproc)
      //           for (unsigned j = slaveNodes.begin(i); j < slaveNodes.end(i); j++) { // loop on all the columns of restriction[lproc][inode]
      //             unsigned jnode = slaveNodes[i][j];
      //             double value = slaveNodesValues[i][j];
      //             if (inode != jnode || value > 5.) { // if off-diagonal or hanging node for lproc
      //               restriction[soltype][inode][jnode] =  value;
      //             }
      //             if (restriction[soltype].find(jnode) == restriction[soltype].end()) { // if jnode is not yet a master node for _iproc
      //               counter++;
      //               for (unsigned k = masterNode.begin(); k < masterNode.end(); k++) {
      //                 if (masterNode[k] == jnode) { // and if jnode is also a master node for lproc
      //                   for (unsigned l = slaveNodes.begin(k); l < slaveNodes.end(k); l++) {
      //                     unsigned lnode = slaveNodes[k][l];
      //                     restriction[soltype][jnode][lnode] = slaveNodesValues[k][l]; //copy the rule of lptoc into _iproc
      //                   }
      //                   break;
      //                 }
      //               }
      //             }
      //           }
      //         }
      //       }
      //     }
      //     masterNode.clearBroadcast();
      //     slaveNodes.clearBroadcast();
      //     slaveNodesValues.clearBroadcast();
      //   }
      //   //END filling the restriction object with infos coming form the parallel vectors and matrices
      //
      //
      //   unsigned globalCounter = 0u;
      //
      //   MPI_Allreduce(&counter,
      //                 &globalCounter,
      //                 1,
      //                 MPI_UNSIGNED,
      //                 MPI_SUM,
      //                 MPI_COMM_WORLD);
      //
      //   counter = globalCounter;
      // }
      unsigned counter = 1;

      while (counter != 0) {
        counter = 0;

        // ------------------------------------------------------------
        // Save restriction object into parallel vectors/matrices
        // ------------------------------------------------------------

        MyVector<unsigned> rowSize(restriction[soltype].size(), 0);

        unsigned cnt1 = 0;
        for (std::map<unsigned, std::map<unsigned, double>>::iterator it1 =
               restriction[soltype].begin();
             it1 != restriction[soltype].end();
             ++it1) {

          rowSize[cnt1] = it1->second.size();
          cnt1++;
        }

        rowSize.stack();

        std::vector<unsigned> offset = rowSize.getOffset();

        MyVector<unsigned> masterNode(offset);
        MyMatrix<unsigned> slaveNodes(rowSize);
        MyMatrix<double> slaveNodesValues(rowSize);

        cnt1 = 0;
        for (std::map<unsigned, std::map<unsigned, double>>::iterator it1 =
               restriction[soltype].begin();
             it1 != restriction[soltype].end();
             ++it1) {

          const unsigned inode = it1->first;
          const std::map<unsigned, double>& row = it1->second;

          const unsigned rowIndex = offset[_iproc] + cnt1;

          masterNode[rowIndex] = inode;

          unsigned cnt2 = 0;
          for (std::map<unsigned, double>::const_iterator it2 = row.begin();
               it2 != row.end();
               ++it2) {

            slaveNodes[rowIndex][cnt2]       = it2->first;
            slaveNodesValues[rowIndex][cnt2] = it2->second;
            cnt2++;
          }

          cnt1++;
        }

        // ------------------------------------------------------------
        // Fill restriction object with information from broadcast data
        // ------------------------------------------------------------

        const unsigned solutionOffset   = msh->_dofOffset[soltype][_iproc];
        const unsigned solutionOffsetp1 = msh->_dofOffset[soltype][_iproc + 1];

        for (unsigned lproc = 0; lproc < _nprocs; lproc++) {

          masterNode.broadcast(lproc);
          slaveNodes.broadcast(lproc);
          slaveNodesValues.broadcast(lproc);

          // Build lookup: master node id -> row index in received data
          std::unordered_map<unsigned, unsigned> masterIndex;

          for (unsigned k = masterNode.begin(); k < masterNode.end(); k++) {
            masterIndex[masterNode[k]] = k;
          }

          for (unsigned i = slaveNodes.begin(); i < slaveNodes.end(); i++) {

            const unsigned inode = masterNode[i];

            std::map<unsigned, std::map<unsigned, double>>::iterator inodeIt =
                  restriction[soltype].find(inode);

            if (inode >= solutionOffset &&
                inode < solutionOffsetp1 &&
                inodeIt == restriction[soltype].end()) {

              counter++;

              std::map<unsigned, double>& rowInode = restriction[soltype][inode];

              for (unsigned j = slaveNodes.begin(i); j < slaveNodes.end(i); j++) {
                const unsigned jnode = slaveNodes[i][j];
                rowInode[jnode] = slaveNodesValues[i][j];
              }
            }
            else {

              if (inodeIt != restriction[soltype].end()) {

                std::map<unsigned, double>& rowInode = inodeIt->second;

                for (unsigned j = slaveNodes.begin(i); j < slaveNodes.end(i); j++) {

                  const unsigned jnode = slaveNodes[i][j];
                  const double value   = slaveNodesValues[i][j];

                  if (inode != jnode || value > 5.0) {
                    rowInode[jnode] = value;
                  }

                  if (restriction[soltype].find(jnode) == restriction[soltype].end()) {

                    auto mit = masterIndex.find(jnode);

                    if (mit != masterIndex.end()) {

                      counter++;

                      const unsigned k = mit->second;

                      std::map<unsigned, double>& rowJnode =
                        restriction[soltype][jnode];

                      for (unsigned l = slaveNodes.begin(k);
                           l < slaveNodes.end(k);
                           l++) {

                        const unsigned lnode = slaveNodes[k][l];
                        rowJnode[lnode] = slaveNodesValues[k][l];
                      }
                    }
                  }
                }
              }
            }
          }

          masterNode.clearBroadcast();
          slaveNodes.clearBroadcast();
          slaveNodesValues.clearBroadcast();
        }

        unsigned globalCounter = 0u;

        MPI_Allreduce(&counter,
                      &globalCounter,
                      1,
                      MPI_UNSIGNED,
                      MPI_SUM,
                      MPI_COMM_WORLD);

        counter = globalCounter;
      }
      //delete pvector;

      __amr_t_parallel_closure[soltype] += MPI_Wtime() - __amr_t0;

      __amr_t0 = MPI_Wtime();

      std::vector<std::vector < unsigned > > genealogy;
      std::vector<std::vector < double > > heredity;
      std::vector< unsigned > index;
      std::map < unsigned,  std::map < unsigned, double  > >  restrictionCopy = restriction[soltype];

      for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restrictionCopy.begin(); it1 != restrictionCopy.end(); it1++) { // loop all over master, hanging and master+hanging nodes

        genealogy.resize(1);
        heredity.resize(1);
        index.resize(1);

        genealogy[0].resize(1);
        heredity[0].resize(1);

        unsigned inode = it1->first;

        if (restrictionCopy[inode][inode] < 5.) { // only if a real master node

          // initialize master node genealogy and heredity at level 0
          genealogy[0][0] = inode;
          heredity[0][0] = 1.;
          index[0] = 0;

          restriction[soltype][inode].clear();
          restriction[soltype][inode][inode] = 1.;

          unsigned level = 1;
          while (level > 0) {

            // initialize master node genealogy and heredity at genemeric level

            unsigned father = genealogy[level - 1][index[level - 1]];

            genealogy.resize(level + 1);
            genealogy[level].reserve(restrictionCopy[ father ].size() - 1);
            genealogy[level].resize(0);

            heredity.resize(level + 1);
            heredity[level].reserve(restrictionCopy[ father ].size() - 1);
            heredity[level].resize(0);

            index.resize(level + 1);
            index[level] = 0;

            unsigned cnt  = 0;
            for (std::map <unsigned, double>::iterator it2 = restrictionCopy[ father ].begin(); it2 != restrictionCopy[ father ].end(); it2++) { // loop on all the father sons
              unsigned son = it2->first;
              bool alreadyFound = false;
              for (unsigned klevel = 0; klevel < level; klevel++) { // check if the son is in the previous genealogy
                for (unsigned k = 0; k < genealogy[klevel].size(); k++) {
                  if (genealogy[klevel][k] == son) alreadyFound = true;
                }
              }
              if (!alreadyFound) {  // if never found add the the restionction value in the master node line and zero the hanging node line
                genealogy[level].resize(genealogy[level].size() + 1);
                heredity[level].resize(heredity[level].size() + 1);

                genealogy[level][cnt] = son;
                heredity[level][cnt] = it2->second * heredity[level - 1][index[level - 1]];

                restriction[soltype][inode][son] += heredity[level][cnt];
                cnt++;

                restriction[soltype][son].clear();
                restriction[soltype][son][son] = 0.;
              }
            }

            if (cnt > 0) {
              level++;
            }
            else {
              bool test = true;
              while (test && level > 0) {
                index[level - 1]++;
                test = false;
                if (index[level - 1] == genealogy[level - 1].size()) {
                  level--;
                  test = true;
                }
              }
            }
          }
        }
        else {
          restriction[soltype][inode].clear();
          restriction[soltype][inode][inode] = 0.;
        }
      }

      __amr_t_genealogy[soltype] += MPI_Wtime() - __amr_t0;

      __amr_t0 = MPI_Wtime();

      MyVector <unsigned> InterfaceSolidMarkNode(interfaceSolidMark[soltype].size());
      MyVector <short unsigned> InterfaceSolidMarkValue(interfaceSolidMark[soltype].size());

      unsigned cnt = 0;
      for (std::map<unsigned, bool >::iterator it = interfaceSolidMark[soltype].begin(); it != interfaceSolidMark[soltype].end(); it++) {
        InterfaceSolidMarkNode[cnt] = it->first;
        InterfaceSolidMarkValue[cnt] = it->second;
        cnt++;
      }
      InterfaceSolidMarkNode.stack();
      InterfaceSolidMarkValue.stack();

      for (unsigned lproc = 0; lproc < _nprocs; lproc++) {
        InterfaceSolidMarkNode.broadcast(lproc);
        InterfaceSolidMarkValue.broadcast(lproc);
        for (unsigned i = InterfaceSolidMarkNode.begin(); i < InterfaceSolidMarkNode.end(); i++) {
          unsigned jnode = InterfaceSolidMarkNode[i];
          if (restriction[soltype].find(jnode) != restriction[soltype].end()) {
            interfaceSolidMark[soltype][jnode] = InterfaceSolidMarkValue[i];
          }
        }
        InterfaceSolidMarkNode.clearBroadcast();
        InterfaceSolidMarkValue.clearBroadcast();
      }

      __amr_t_solidmark_bcast[soltype] += MPI_Wtime() - __amr_t0;
    }

    const double __amr_total = MPI_Wtime() - __amr_total_start;

    if (false && _iproc == 0u) {
      std::cout << "\n========== GetAMRRestriction timing ==========" << std::endl;
      std::cout << "total                              : " << __amr_total << std::endl;
      std::cout << "interface build                    : " << __amr_t_interface_build << std::endl;
      std::cout << "bounding boxes + grid              : " << __amr_t_bounding_boxes << std::endl;

      double __amr_accounted =
        __amr_t_interface_build + __amr_t_bounding_boxes;

      for (unsigned s = 0u; s < 3u; ++s) {
        const double subtotal =
          __amr_t_jmin[s] +
          __amr_t_candidate_exchange[s] +
          __amr_t_restriction_scan[s] +
          __amr_t_parallel_closure[s] +
          __amr_t_genealogy[s] +
          __amr_t_solidmark_bcast[s];

        __amr_accounted += subtotal;

        std::cout << "---- soltype " << s << " ----" << std::endl;
        std::cout << "jMin/jSize scan                    : " << __amr_t_jmin[s] << std::endl;
        std::cout << "candidate exchange method          : " << __amr_t_candidate_exchange[s] << std::endl;
        std::cout << "restriction scan total             : " << __amr_t_restriction_scan[s] << std::endl;
        std::cout << "  bbox lookup subset               : " << __amr_t_bbx_lookup[s] << std::endl;
        std::cout << "  candidate-element loop subset    : " << __amr_t_candidate_loop[s] << std::endl;
        std::cout << "  same-node check subset           : " << __amr_t_same_node[s] << std::endl;
        std::cout << "  geometry filter subset           : " << __amr_t_geometry_filter[s] << std::endl;
        std::cout << "  projection/inverse subset        : " << __amr_t_projection_inverse[s] << std::endl;
        std::cout << "  eval/store subset                : " << __amr_t_eval_store[s] << std::endl;
        std::cout << "parallel closure                   : " << __amr_t_parallel_closure[s] << std::endl;
        std::cout << "genealogy flattening               : " << __amr_t_genealogy[s] << std::endl;
        std::cout << "solid-mark broadcast               : " << __amr_t_solidmark_bcast[s] << std::endl;
        std::cout << "received candidate points          : " << __amr_cnt_received_points[s] << std::endl;
        std::cout << "candidate elements tested          : " << __amr_cnt_candidate_elements[s] << std::endl;
        std::cout << "same-node skips                    : " << __amr_cnt_same_node_skip[s] << std::endl;
        std::cout << "passed geometry filter             : " << __amr_cnt_geometry_pass[s] << std::endl;
        std::cout << "inverse mapping calls              : " << __amr_cnt_inverse_calls[s] << std::endl;
        std::cout << "inside reference domain            : " << __amr_cnt_inside_domain[s] << std::endl;
        std::cout << "stored entries                     : " << __amr_cnt_stored_entries[s] << std::endl;
      }

      std::cout << "accounted subtotal                 : " << __amr_accounted << std::endl;
      std::cout << "unaccounted / local overhead       : " << (__amr_total - __amr_accounted) << std::endl;
      std::cout << "==============================================\n" << std::endl;
    }
  }
  /*
    void elem::GetAMRRestriction(Mesh *msh) {

      std::vector < std::map < unsigned,  std::map < unsigned, double  > > > &restriction = msh->GetAmrRestrictionMap();
      restriction.resize(3);

      std::vector < std::map < unsigned, bool > > &interfaceSolidMark = msh->GetAmrSolidMark();
      interfaceSolidMark.resize(3);

      std::vector < MyVector<unsigned> > interfaceElement; // iel = interfaceElement[ilevel][i]
      std::vector < MyMatrix<unsigned> > interfaceLocalDof; // ldof = interfaceElement[ilevel][i][j]
      std::vector < std::vector < MyMatrix<unsigned> > > interfaceDof; // gdof = interfaceDof[solType][ilevel][i][j]
      std::vector < std::vector < MyMatrix<unsigned> > > levelInterfaceSolidMark; // FSIdof = levelInterfaceSolidMark[solType][ilevel][i][j]
      std::vector < std::vector < MyMatrix< double > > > interfaceNodeCoordinates; // x[dim] = interfaceNodeCoordinates[ilevel][dim][i][j]

      interfaceElement.resize(_level + 1);
      interfaceLocalDof.resize(_level + 1);
      interfaceDof.resize(3);
      levelInterfaceSolidMark.resize(3);
      for (unsigned i = 0; i < 3; i++) {
        interfaceDof[i].resize(_level + 1);
        levelInterfaceSolidMark[i].resize(_level + 1);
      }
      interfaceNodeCoordinates.resize(_level + 1);
      unsigned dim = msh->GetDimension();

      for (unsigned ilevel = 0; ilevel <= _level; ilevel++) {
        //BEGIN interface element search
        interfaceElement[ilevel] = MyVector <unsigned> (_elementOwned);
        unsigned counter = 0;
        for (unsigned i = _elementLevel.begin(); i < _elementLevel.end(); i++) {
          if (ilevel == _elementLevel[i]) {
            for (unsigned j = _elementNearFace.begin(i); j < _elementNearFace.end(i); j++) {
              if (-1 == _elementNearFace[i][j]) {
                interfaceElement[ilevel][counter] = i;
                counter++;
                break;
              }
            }
          }
        }
        interfaceElement[ilevel].resize(counter);
        interfaceElement[ilevel].stack();
        //END interface element search

        //BEGIN interface node search
        std::vector< unsigned > offset = interfaceElement[ilevel].getOffset();
        interfaceLocalDof[ilevel] = MyMatrix <unsigned>(offset, NVE[0][2], UINT_MAX);
        for (unsigned i = interfaceElement[ilevel].begin(); i < interfaceElement[ilevel].end(); i++) {
          unsigned iel =  interfaceElement[ilevel][i];
          std::map <unsigned, bool> ldofs;
          for (unsigned jface = _elementNearFace.begin(iel); jface < _elementNearFace.end(iel); jface++) {
            if (-1 == _elementNearFace[iel][jface]) {
              for (unsigned k = 0; k < GetNFACENODES(GetElementType(iel), jface, 2); k++) {
                unsigned index = GetIG(GetElementType(iel), jface, k);
                ldofs[index] = true;
              }
            }
          }
          unsigned j = 0;
          for (std::map<unsigned, bool>::iterator it = ldofs.begin(); it != ldofs.end(); it++) {
            interfaceLocalDof[ilevel][i][j] = it->first;
            j++;
          }
        }
        interfaceLocalDof[ilevel].shrinkToFit(UINT_MAX);
        //END interface node search

        //BEGIN interface node dof global search, one for each soltype
        MyVector <unsigned> rowSize = interfaceLocalDof[ilevel].getRowSize();
        for (unsigned soltype = 0; soltype < 3; soltype++) {
          interfaceDof[soltype][ilevel] = MyMatrix< unsigned > (rowSize, UINT_MAX);
          levelInterfaceSolidMark[soltype][ilevel] = MyMatrix< unsigned > (rowSize, UINT_MAX);
          for (unsigned i = interfaceLocalDof[ilevel].begin(); i < interfaceLocalDof[ilevel].end(); i++) {
            unsigned iel = interfaceElement[ilevel][i];
            unsigned counter = 0;
            for (unsigned j = interfaceLocalDof[ilevel].begin(i); j < interfaceLocalDof[ilevel].end(i); j++) {
              unsigned jloc = interfaceLocalDof[ilevel][i][j];
              if (jloc < GetElementDofNumber(iel, soltype)) {
                unsigned jdof  = msh->GetSolutionDof(jloc, iel, soltype);
                interfaceDof[soltype][ilevel][i][counter] = jdof;
                unsigned jdof2  = msh->GetSolutionDof(jloc, iel, 2);
                levelInterfaceSolidMark[soltype][ilevel][i][counter] = msh->GetSolidMark(jdof2);
                counter++;
              }
              else {
                break;
              }
            }
          }
          interfaceDof[soltype][ilevel].shrinkToFit(UINT_MAX);
          levelInterfaceSolidMark[soltype][ilevel].shrinkToFit(UINT_MAX);
        }
        //END interface node dof global search, one for each soltype

        //BEGIN interface node coordinates
        interfaceNodeCoordinates[ilevel].resize(dim);
        for (unsigned k = 0; k < dim; k++) {
          interfaceNodeCoordinates[ilevel][k] = MyMatrix< double > (rowSize, 0.);
        }
        for (unsigned i = interfaceLocalDof[ilevel].begin(); i < interfaceLocalDof[ilevel].end(); i++) {
          unsigned iel = interfaceElement[ilevel][i];
          for (unsigned j = interfaceLocalDof[ilevel].begin(i); j < interfaceLocalDof[ilevel].end(i); j++) {
            unsigned jnode = interfaceLocalDof[ilevel][i][j];
            unsigned xDof  = msh->GetSolutionDof(jnode, iel, 2);
            for (unsigned k = 0; k < dim; k++) {
              interfaceNodeCoordinates[ilevel][k][i][j] = (*msh->_topology->_Sol[k])(xDof);
            }
          }
        }
        //END interface node coordinates search
      }


      const unsigned nLevels = _level + 1u;

      std::vector<std::vector<double*> > xMin(_nprocs);
      std::vector<std::vector<double*> > xMax(_nprocs);

      std::vector<double> xMinMemory(_nprocs * nLevels * dim);
      std::vector<double> xMaxMemory(_nprocs * nLevels * dim);


      std::vector<std::vector<std::vector<unsigned>>> bbx_to_elements;
      std::vector<std::vector<unsigned>> bbxN;

      GetAMRInterfaceBoundingBoxes(msh, interfaceElement, xMin, xMax, xMinMemory, xMaxMemory, bbx_to_elements, bbxN);

      for (unsigned soltype = 0; soltype < 3; soltype++) {

        std::vector<double> xl(dim);
        std::vector<double> phi;
        std::vector<std::vector<double>> gradPhi(dim);
        std::vector < std::vector <double > > xv(dim);
        std::vector <double> xc(dim);
        std::vector < std::vector< double > > xe(dim);
        std::vector <double> xi(dim);
        std::vector < std::vector < std::vector <double > > > aP(3);

        std::vector<unsigned> jMinLevel(nLevels, UINT_MAX);
        std::vector<unsigned> jMaxLevel(nLevels, 0);
        std::vector<unsigned> jSizeLevel(nLevels, 0);

        for (unsigned level = 0; level < nLevels; level++) {
          for (unsigned k = interfaceDof[soltype][level].begin(); k < interfaceDof[soltype][level].end(); k++) {
            for (unsigned l = interfaceDof[soltype][level].begin(k); l < interfaceDof[soltype][level].end(k); l++) {
              const unsigned ldof = interfaceDof[soltype][level][k][l];

              if (ldof < jMinLevel[level]) jMinLevel[level] = ldof;
              if (ldof > jMaxLevel[level]) jMaxLevel[level] = ldof;
            }
          }

          if (jMaxLevel[level] >= jMinLevel[level]) {
            jSizeLevel[level] = jMaxLevel[level] - jMinLevel[level] + 1u;
          }
        }

        std::vector<std::vector<unsigned>> jDof_r(_nprocs);
        std::vector<std::vector<unsigned>> jSolidMark_r(_nprocs);
        std::vector<std::vector<std::vector<double>>> jCoordinate_r(_nprocs);

        for(unsigned lproc = 0u; lproc < _nprocs; lproc++) {
          jCoordinate_r[lproc].resize(dim);
        }

        for (int ilevel = 0; ilevel < _level; ilevel++) {
          for (int jlevel = ilevel + 1; jlevel <= _level; jlevel++) {

            const unsigned jMin  = jMinLevel[jlevel];
            const unsigned jSize = jSizeLevel[jlevel];

            GetAMRRestrictionCandidateDataForLevelPair(
              static_cast<unsigned>(ilevel),                // ilevel
              dim,                                          // dim
              jMin,                                         // jMin
              jSize,                                        // jSize
              interfaceDof[soltype][jlevel],                // interfaceDofJ
              levelInterfaceSolidMark[soltype][jlevel],     // levelInterfaceSolidMarkJ
              interfaceNodeCoordinates[jlevel],             // interfaceNodeCoordinatesJ
              xMin,                                         // xMin
              xMax,                                         // xMax
              jDof_r,                                       // output
              jSolidMark_r,                                 // output
              jCoordinate_r                                 // output
            );

            for (unsigned lproc = 0; lproc < _nprocs; lproc++) {

              for (unsigned kl = 0; kl < jDof_r[lproc].size(); kl++) {

                const unsigned ldof = jDof_r[lproc][kl];

                for (unsigned d = 0; d < dim; d++) {
                  xl[d] = jCoordinate_r[lproc][d][kl];
                }

                unsigned bbxIdx = 0u;

                const bool pointInsideGrid =
                  GetAMRBBoxIndex(ilevel,
                                  dim,
                                  xl,
                                  xMin,
                                  xMax,
                                  bbxN,
                                  bbxIdx);

                if (!pointInsideGrid) continue;
                if (bbxIdx >= bbx_to_elements[ilevel].size()) continue;

                const std::vector<unsigned>& candidateElements =
                  bbx_to_elements[ilevel][bbxIdx];

                for (unsigned ee = 0; ee < candidateElements.size(); ee++) {

                  const unsigned iel = interfaceElement[ilevel][candidateElements[ee]];

                  const unsigned i = candidateElements[ee];

                  bool sameNode = false;

                  const unsigned nElementDofs = GetElementDofNumber(iel, soltype);

                  for (unsigned j = 0; j < nElementDofs; j++) {
                    const unsigned jdof = msh->GetSolutionDof(j, iel, soltype);

                    if (jdof == ldof) {
                      sameNode = true;
                      break;
                    }
                  }

                  if (sameNode) continue;

                  short unsigned ielType = _elementType[iel];
                  basis* base = msh->GetBasis(ielType, soltype);

                  msh->GetElementNodeCoordinates(xv, iel);

                  double r2;
                  GetConvexHullSphereRadiousSquare(xv, xc, r2, 0.01);

                  double d2 = 0.0;

                  for (unsigned d = 0; d < dim; d++) {
                    d2 += (xl[d] - xc[d]) * (xl[d] - xc[d]);
                  }

                  if (d2 > r2) continue;

                  GetBoundingBox(xv, xe, 0.01);

                  bool insideHull = true;

                  for (unsigned d = 0; d < dim; d++) {
                    if (xl[d] < xe[d][0] || xl[d] > xe[d][1]) {
                      insideHull = false;
                      break;
                    }
                  }

                  if (!insideHull) continue;

                  for (unsigned jtype = 0; jtype < 3; jtype++) {
                    ProjectNodalToPolynomialCoefficients(aP[jtype], xv, ielType, jtype);
                  }

                  GetClosestPointInReferenceElement(xv, xl, ielType, xi);

                  GetInverseMapping_fast(2u,
                                         ielType,
                                         aP,
                                         xl,
                                         xi,
                                         100u,
                                         phi,
                                         gradPhi);

                  const bool insideDomain =
                    CheckIfPointIsInsideReferenceDomain(xi, ielType, 0.0001);

                  if (!insideDomain) continue;

                  for (unsigned j = interfaceDof[soltype][ilevel].begin(i);
                       j < interfaceDof[soltype][ilevel].end(i);
                       j++) {

                    const unsigned jloc = interfaceLocalDof[ilevel][i][j];

                    const double value = base->eval_phi(jloc, xi);

                    if (fabs(value) >= 1.0e-10) {
                      const unsigned jdof = interfaceDof[soltype][ilevel][i][j];

                      auto& rowJ = restriction[soltype][jdof];

                      auto it = rowJ.find(jdof);

                      if (it == rowJ.end()) {
                        rowJ.emplace(jdof, 1.0);
                        interfaceSolidMark[soltype][jdof] =
                          levelInterfaceSolidMark[soltype][ilevel][i][j];
                      }

                      rowJ[ldof] = value;

                      auto& rowL = restriction[soltype][ldof];
                      rowL[ldof] = 10.0;

                      interfaceSolidMark[soltype][ldof] = jSolidMark_r[lproc][kl];
                    }
                  }
                }
              }
            }
          }
        }










        NumericVector* pvector;
        pvector = NumericVector::build().release();
        pvector->init(_nprocs, 1, false, AUTOMATIC);

        unsigned counter = 1;
        while (counter != 0) {
          counter = 0;

          //BEGIN  saving the restriction object in parallel vectors and matrices

          MyVector <unsigned> rowSize(restriction[soltype].size(), 0);
          unsigned cnt1 = 0;
          for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
            rowSize[cnt1] = restriction[soltype][it1->first].size();
            cnt1++;
          }
          rowSize.stack();

          std::vector< unsigned > offset = rowSize.getOffset();

          MyVector <unsigned> masterNode(offset);
          MyMatrix <unsigned> slaveNodes(rowSize);
          MyMatrix <double> slaveNodesValues(rowSize);

          cnt1 = 0;
          for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
            masterNode[offset[_iproc] + cnt1] = it1->first;
            unsigned cnt2 = 0;
            for (std::map<unsigned, double> ::iterator it2 = restriction[soltype][it1->first].begin(); it2 != restriction[soltype][it1->first].end(); it2++) {
              slaveNodes[offset[_iproc] + cnt1][cnt2] = it2->first;
              slaveNodesValues[offset[_iproc] + cnt1][cnt2] = it2->second;
              cnt2++;
            }
            cnt1++;
          }
          //END  saving the restriction object in parallel vectors and matrices


          //BEGIN filling the restriction object with infos coming form the parallel vectors and matrices
          unsigned solutionOffset = msh->_dofOffset[soltype][_iproc];
          unsigned solutionOffsetp1 = msh->_dofOffset[soltype][_iproc + 1];
          for (unsigned lproc = 0; lproc < _nprocs; lproc++) {
            masterNode.broadcast(lproc);
            slaveNodes.broadcast(lproc);
            slaveNodesValues.broadcast(lproc);
            for (unsigned i = slaveNodes.begin(); i < slaveNodes.end(); i++) {
              unsigned inode = masterNode[i];
              if (inode >= solutionOffset && inode < solutionOffsetp1 && // inode belongs to _iproc
                  restriction[soltype].find(inode) == restriction[soltype].end()) { // but inode is not set as master node of _iproc
                counter++;
                for (unsigned j = slaveNodes.begin(i); j < slaveNodes.end(i); j++) { //copy information for lproc to _iproc
                  unsigned jnode = slaveNodes[i][j];
                  restriction[soltype][inode][jnode] = slaveNodesValues[i][j];
                }
              }
              else { // either inode does not belong to _iproc or it was already defined as a master for _iproc
                if (restriction[soltype].find(inode) != restriction[soltype].end()) { // inode is already defined as master node in _iproc (either it does or does not belong to _iproc)
                  for (unsigned j = slaveNodes.begin(i); j < slaveNodes.end(i); j++) { // loop on all the columns of restriction[lproc][inode]
                    unsigned jnode = slaveNodes[i][j];
                    double value = slaveNodesValues[i][j];
                    if (inode != jnode || value > 5.) { // if off-diagonal or hanging node for lproc
                      restriction[soltype][inode][jnode] =  value;
                    }
                    if (restriction[soltype].find(jnode) == restriction[soltype].end()) { // if jnode is not yet a master node for _iproc
                      counter++;
                      for (unsigned k = masterNode.begin(); k < masterNode.end(); k++) {
                        if (masterNode[k] == jnode) { // and if jnode is also a master node for lproc
                          for (unsigned l = slaveNodes.begin(k); l < slaveNodes.end(k); l++) {
                            unsigned lnode = slaveNodes[k][l];
                            restriction[soltype][jnode][lnode] = slaveNodesValues[k][l]; //copy the rule of lptoc into _iproc
                          }
                          break;
                        }
                      }
                    }
                  }
                }
              }
            }
            masterNode.clearBroadcast();
            slaveNodes.clearBroadcast();
            slaveNodesValues.clearBroadcast();
          }
          //END filling the restriction object with infos coming form the parallel vectors and matrices


          pvector->set(_iproc, counter);
          pvector->close();
          counter = static_cast <unsigned>(floor(pvector->l1_norm() + 0.5));
        }
        delete pvector;

        std::vector<std::vector < unsigned > > genealogy;
        std::vector<std::vector < double > > heredity;
        std::vector< unsigned > index;
        std::map < unsigned,  std::map < unsigned, double  > >  restrictionCopy = restriction[soltype];

        for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restrictionCopy.begin(); it1 != restrictionCopy.end(); it1++) { // loop all over master, hanging and master+hanging nodes

          genealogy.resize(1);
          heredity.resize(1);
          index.resize(1);

          genealogy[0].resize(1);
          heredity[0].resize(1);

          unsigned inode = it1->first;

          if (restrictionCopy[inode][inode] < 5.) { // only if a real master node

            // initialize master node genealogy and heredity at level 0
            genealogy[0][0] = inode;
            heredity[0][0] = 1.;
            index[0] = 0;

            restriction[soltype][inode].clear();
            restriction[soltype][inode][inode] = 1.;

            unsigned level = 1;
            while (level > 0) {

              // initialize master node genealogy and heredity at genemeric level

              unsigned father = genealogy[level - 1][index[level - 1]];

              genealogy.resize(level + 1);
              genealogy[level].reserve(restrictionCopy[ father ].size() - 1);
              genealogy[level].resize(0);

              heredity.resize(level + 1);
              heredity[level].reserve(restrictionCopy[ father ].size() - 1);
              heredity[level].resize(0);

              index.resize(level + 1);
              index[level] = 0;

              unsigned cnt  = 0;
              for (std::map <unsigned, double>::iterator it2 = restrictionCopy[ father ].begin(); it2 != restrictionCopy[ father ].end(); it2++) { // loop on all the father sons
                unsigned son = it2->first;
                bool alreadyFound = false;
                for (unsigned klevel = 0; klevel < level; klevel++) { // check if the son is in the previous genealogy
                  for (unsigned k = 0; k < genealogy[klevel].size(); k++) {
                    if (genealogy[klevel][k] == son) alreadyFound = true;
                  }
                }
                if (!alreadyFound) {  // if never found add the the restionction value in the master node line and zero the hanging node line
                  genealogy[level].resize(genealogy[level].size() + 1);
                  heredity[level].resize(heredity[level].size() + 1);

                  genealogy[level][cnt] = son;
                  heredity[level][cnt] = it2->second * heredity[level - 1][index[level - 1]];

                  restriction[soltype][inode][son] += heredity[level][cnt];
                  cnt++;

                  restriction[soltype][son].clear();
                  restriction[soltype][son][son] = 0.;
                }
              }

              if (cnt > 0) {
                level++;
              }
              else {
                bool test = true;
                while (test && level > 0) {
                  index[level - 1]++;
                  test = false;
                  if (index[level - 1] == genealogy[level - 1].size()) {
                    level--;
                    test = true;
                  }
                }
              }
            }
          }
          else {
            restriction[soltype][inode].clear();
            restriction[soltype][inode][inode] = 0.;
          }
        }

        MyVector <unsigned> InterfaceSolidMarkNode(interfaceSolidMark[soltype].size());
        MyVector <short unsigned> InterfaceSolidMarkValue(interfaceSolidMark[soltype].size());

        unsigned cnt = 0;
        for (std::map<unsigned, bool >::iterator it = interfaceSolidMark[soltype].begin(); it != interfaceSolidMark[soltype].end(); it++) {
          InterfaceSolidMarkNode[cnt] = it->first;
          InterfaceSolidMarkValue[cnt] = it->second;
          cnt++;
        }
        InterfaceSolidMarkNode.stack();
        InterfaceSolidMarkValue.stack();

        for (unsigned lproc = 0; lproc < _nprocs; lproc++) {
          InterfaceSolidMarkNode.broadcast(lproc);
          InterfaceSolidMarkValue.broadcast(lproc);
          for (unsigned i = InterfaceSolidMarkNode.begin(); i < InterfaceSolidMarkNode.end(); i++) {
            unsigned jnode = InterfaceSolidMarkNode[i];
            if (restriction[soltype].find(jnode) != restriction[soltype].end()) {
              interfaceSolidMark[soltype][jnode] = InterfaceSolidMarkValue[i];
            }
          }
          InterfaceSolidMarkNode.clearBroadcast();
          InterfaceSolidMarkValue.clearBroadcast();
        }
      }
    }*/

  void elem::GetAMRInterfaceBoundingBoxes(
    Mesh * msh,
    const std::vector<MyVector<unsigned>>& interfaceElement,
    std::vector<std::vector<double*>>& xMin,
    std::vector<std::vector<double*>>& xMax,
    std::vector<double>& xMinMemory,
    std::vector<double>& xMaxMemory,
    std::vector<std::vector<std::vector<unsigned>>>& bbx_to_elements,
    std::vector<std::vector<unsigned>>& bbxN
  ) {
    const unsigned dim = msh->GetDimension();
    const unsigned nLevels = _level + 1u;

    std::vector<double> xMinLocal(nLevels * dim, DBL_MAX);
    std::vector<double> xMaxLocal(nLevels * dim, -DBL_MAX);

    std::vector<double*> xMinL(nLevels);
    std::vector<double*> xMaxL(nLevels);

    for (unsigned level = 0; level < nLevels; level++) {
      xMinL[level] = xMinLocal.data() + level * dim;
      xMaxL[level] = xMaxLocal.data() + level * dim;
    }

    // ------------------------------------------------------------
    // Compute local level bounding boxes
    // ------------------------------------------------------------
    for (unsigned ilevel = 0; ilevel < nLevels; ilevel++) {
      for (unsigned i = interfaceElement[ilevel].begin();
           i < interfaceElement[ilevel].end();
           i++) {

        const unsigned iel = interfaceElement[ilevel][i];
        const unsigned ndofs = GetElementDofNumber(iel, 2u);

        for (unsigned j = 0; j < ndofs; j++) {
          const unsigned xdof = msh->GetSolutionDof(j, iel, 2u);

          for (unsigned d = 0; d < dim; d++) {
            const double x = (*msh->_topology->_Sol[d])(xdof);

            if (x < xMinL[ilevel][d]) xMinL[ilevel][d] = x;
            if (x > xMaxL[ilevel][d]) xMaxL[ilevel][d] = x;
          }
        }
      }
    }

    // ------------------------------------------------------------
    // Enlarge each valid local level bounding box
    // ------------------------------------------------------------
    for (unsigned level = 0; level < nLevels; level++) {
      bool hasPoints = true;

      for (unsigned d = 0; d < dim; d++) {
        if (xMinL[level][d] == DBL_MAX || xMaxL[level][d] == -DBL_MAX) {
          hasPoints = false;
          break;
        }
      }

      if (!hasPoints) continue;

      for (unsigned d = 0; d < dim; d++) {
        const double length = xMaxL[level][d] - xMinL[level][d];
        const double delta = 0.01 * length + 1.0e-14;

        xMinL[level][d] -= delta;
        xMaxL[level][d] += delta;
      }
    }

    // ------------------------------------------------------------
    // Gather bounding boxes from all processors
    // ------------------------------------------------------------
    xMinMemory.resize(_nprocs * nLevels * dim);
    xMaxMemory.resize(_nprocs * nLevels * dim);

    MPI_Allgather(xMinLocal.data(), nLevels * dim, MPI_DOUBLE,
                  xMinMemory.data(), nLevels * dim, MPI_DOUBLE,
                  MPI_COMM_WORLD);

    MPI_Allgather(xMaxLocal.data(), nLevels * dim, MPI_DOUBLE,
                  xMaxMemory.data(), nLevels * dim, MPI_DOUBLE,
                  MPI_COMM_WORLD);

    xMin.resize(_nprocs);
    xMax.resize(_nprocs);

    for (unsigned proc = 0; proc < _nprocs; proc++) {
      xMin[proc].resize(nLevels);
      xMax[proc].resize(nLevels);

      for (unsigned level = 0; level < nLevels; level++) {
        xMin[proc][level] =
          xMinMemory.data() + (proc * nLevels + level) * dim;

        xMax[proc][level] =
          xMaxMemory.data() + (proc * nLevels + level) * dim;
      }
    }

    // ------------------------------------------------------------
    // Build local bounding-box grid:
    // bbx_to_elements[level][flat_box_idx] = list of elements
    // bbxN[level][d] = number of boxes in direction d
    // ------------------------------------------------------------
    bbx_to_elements.clear();
    bbx_to_elements.resize(nLevels);

    bbxN.clear();
    bbxN.resize(nLevels, std::vector<unsigned>(dim, 1u));

    for (unsigned level = 0; level < nLevels; level++) {
      bool hasPoints = true;

      for (unsigned d = 0; d < dim; d++) {
        if (xMinL[level][d] == DBL_MAX || xMaxL[level][d] == -DBL_MAX) {
          hasPoints = false;
          break;
        }
      }

      if (!hasPoints) {
        bbx_to_elements[level].clear();
        continue;
      }

      std::vector<double> minElemSize(dim, DBL_MAX);

      // ----------------------------------------------------------
      // Compute minimum element size in each coordinate direction
      // ----------------------------------------------------------
      for (unsigned i = interfaceElement[level].begin();
           i < interfaceElement[level].end();
           i++) {

        const unsigned iel = interfaceElement[level][i];
        const unsigned ndofs = GetElementDofNumber(iel, 2u);

        std::vector<double> eMin(dim, DBL_MAX);
        std::vector<double> eMax(dim, -DBL_MAX);

        for (unsigned j = 0; j < ndofs; j++) {
          const unsigned xdof = msh->GetSolutionDof(j, iel, 2u);

          for (unsigned d = 0; d < dim; d++) {
            const double x = (*msh->_topology->_Sol[d])(xdof);

            if (x < eMin[d]) eMin[d] = x;
            if (x > eMax[d]) eMax[d] = x;
          }
        }

        for (unsigned d = 0; d < dim; d++) {
          const double h = eMax[d] - eMin[d];

          if (h > 0.0 && h < minElemSize[d]) {
            minElemSize[d] = h;
          }
        }
      }

      // ----------------------------------------------------------
      // Choose number of grid boxes.
      // Box size is approximately half of the minimum element size.
      // ----------------------------------------------------------
      unsigned totalBoxes = 1u;

      for (unsigned d = 0; d < dim; d++) {
        const double domainLength = xMaxL[level][d] - xMinL[level][d];

        if (domainLength <= 0.0 || minElemSize[d] == DBL_MAX || minElemSize[d] <= 0.0) {
          bbxN[level][d] = 1u;
        }
        else {
          const double targetBoxSize = 0.5 * minElemSize[d];

          unsigned nBoxes =
            static_cast<unsigned>(std::ceil(domainLength / targetBoxSize));

          if (nBoxes < 1u) nBoxes = 1u;

          bbxN[level][d] = nBoxes;
        }

        totalBoxes *= bbxN[level][d];
      }

      bbx_to_elements[level].clear();
      bbx_to_elements[level].resize(totalBoxes);

      std::vector<double> boxSize(dim, 1.0);

      for (unsigned d = 0; d < dim; d++) {
        boxSize[d] = (xMaxL[level][d] - xMinL[level][d]) /
                     static_cast<double>(bbxN[level][d]);

        if (boxSize[d] <= 0.0) boxSize[d] = 1.0;
      }

      // ----------------------------------------------------------
      // Assign each element to all grid boxes overlapped
      // ----------------------------------------------------------
      for (unsigned i = interfaceElement[level].begin();
           i < interfaceElement[level].end();
           i++) {

        const unsigned iel = interfaceElement[level][i];
        const unsigned ndofs = GetElementDofNumber(iel, 2u);

        std::vector<double> eMin(dim, DBL_MAX);
        std::vector<double> eMax(dim, -DBL_MAX);

        for (unsigned j = 0; j < ndofs; j++) {
          const unsigned xdof = msh->GetSolutionDof(j, iel, 2u);

          for (unsigned d = 0; d < dim; d++) {
            const double x = (*msh->_topology->_Sol[d])(xdof);

            if (x < eMin[d]) eMin[d] = x;
            if (x > eMax[d]) eMax[d] = x;
          }
        }

        std::vector<unsigned> iMin(dim, 0u);
        std::vector<unsigned> iMax(dim, 0u);

        for (unsigned d = 0; d < dim; d++) {
          int a = static_cast<int>(std::floor((eMin[d] - xMinL[level][d]) / boxSize[d]));
          int b = static_cast<int>(std::floor((eMax[d] - xMinL[level][d]) / boxSize[d]));

          if (a < 0) a = 0;
          if (b < 0) b = 0;

          if (a >= static_cast<int>(bbxN[level][d])) a = static_cast<int>(bbxN[level][d]) - 1;
          if (b >= static_cast<int>(bbxN[level][d])) b = static_cast<int>(bbxN[level][d]) - 1;

          iMin[d] = static_cast<unsigned>(a);
          iMax[d] = static_cast<unsigned>(b);
        }

        if (dim == 1u) {
          for (unsigned ix = iMin[0]; ix <= iMax[0]; ix++) {
            const unsigned idx = ix;
            bbx_to_elements[level][idx].push_back(i);
          }
        }
        else if (dim == 2u) {
          const unsigned nx = bbxN[level][0];
          const unsigned ny = bbxN[level][1];

          for (unsigned iy = iMin[1]; iy <= iMax[1]; iy++) {
            for (unsigned ix = iMin[0]; ix <= iMax[0]; ix++) {
              const unsigned idx = ix + nx * iy;
              bbx_to_elements[level][idx].push_back(i);
            }
          }
        }
        else if (dim == 3u) {
          const unsigned nx = bbxN[level][0];
          const unsigned ny = bbxN[level][1];
          const unsigned nz = bbxN[level][2];

          for (unsigned iz = iMin[2]; iz <= iMax[2]; iz++) {
            for (unsigned iy = iMin[1]; iy <= iMax[1]; iy++) {
              for (unsigned ix = iMin[0]; ix <= iMax[0]; ix++) {
                const unsigned idx = ix + nx * (iy + ny * iz);
                bbx_to_elements[level][idx].push_back(i);
              }
            }
          }
        }
        else {
          std::cout << "GetAMRInterfaceBoundingBoxes: unsupported dimension "
                    << dim << std::endl;
          abort();
        }
      }
    }
  }

  inline bool elem::GetAMRBBoxIndex(
    const unsigned ilevel,
    const unsigned dim,
    const std::vector<double>& xl,
    const std::vector<std::vector<double*>>& xMin,
    const std::vector<std::vector<double*>>& xMax,
    const std::vector<std::vector<unsigned>>& bbxN,
    unsigned & bbxIdx
  ) const {
    bbxIdx = 0u;

    if (dim == 1u) {
      const unsigned nx = bbxN[ilevel][0];

      const double xmin = xMin[_iproc][ilevel][0];
      const double xmax = xMax[_iproc][ilevel][0];

      const double hx = (xmax - xmin) / static_cast<double>(nx);

      if (xl[0] < xmin || xl[0] > xmax || hx <= 0.0) return false;

      int ix = static_cast<int>(std::floor((xl[0] - xmin) / hx));

      if (ix < 0) ix = 0;
      if (ix >= static_cast<int>(nx)) ix = static_cast<int>(nx) - 1;

      bbxIdx = static_cast<unsigned>(ix);
      return true;
    }

    if (dim == 2u) {
      const unsigned nx = bbxN[ilevel][0];
      const unsigned ny = bbxN[ilevel][1];

      const double xmin = xMin[_iproc][ilevel][0];
      const double xmax = xMax[_iproc][ilevel][0];

      const double ymin = xMin[_iproc][ilevel][1];
      const double ymax = xMax[_iproc][ilevel][1];

      const double hx = (xmax - xmin) / static_cast<double>(nx);
      const double hy = (ymax - ymin) / static_cast<double>(ny);

      if (xl[0] < xmin || xl[0] > xmax ||
          xl[1] < ymin || xl[1] > ymax ||
          hx <= 0.0 || hy <= 0.0) {
        return false;
      }

      int ix = static_cast<int>(std::floor((xl[0] - xmin) / hx));
      int iy = static_cast<int>(std::floor((xl[1] - ymin) / hy));

      if (ix < 0) ix = 0;
      if (iy < 0) iy = 0;

      if (ix >= static_cast<int>(nx)) ix = static_cast<int>(nx) - 1;
      if (iy >= static_cast<int>(ny)) iy = static_cast<int>(ny) - 1;

      bbxIdx = static_cast<unsigned>(ix) +
               nx * static_cast<unsigned>(iy);

      return true;
    }

    if (dim == 3u) {
      const unsigned nx = bbxN[ilevel][0];
      const unsigned ny = bbxN[ilevel][1];
      const unsigned nz = bbxN[ilevel][2];

      const double xmin = xMin[_iproc][ilevel][0];
      const double xmax = xMax[_iproc][ilevel][0];

      const double ymin = xMin[_iproc][ilevel][1];
      const double ymax = xMax[_iproc][ilevel][1];

      const double zmin = xMin[_iproc][ilevel][2];
      const double zmax = xMax[_iproc][ilevel][2];

      const double hx = (xmax - xmin) / static_cast<double>(nx);
      const double hy = (ymax - ymin) / static_cast<double>(ny);
      const double hz = (zmax - zmin) / static_cast<double>(nz);

      if (xl[0] < xmin || xl[0] > xmax ||
          xl[1] < ymin || xl[1] > ymax ||
          xl[2] < zmin || xl[2] > zmax ||
          hx <= 0.0 || hy <= 0.0 || hz <= 0.0) {
        return false;
      }

      int ix = static_cast<int>(std::floor((xl[0] - xmin) / hx));
      int iy = static_cast<int>(std::floor((xl[1] - ymin) / hy));
      int iz = static_cast<int>(std::floor((xl[2] - zmin) / hz));

      if (ix < 0) ix = 0;
      if (iy < 0) iy = 0;
      if (iz < 0) iz = 0;

      if (ix >= static_cast<int>(nx)) ix = static_cast<int>(nx) - 1;
      if (iy >= static_cast<int>(ny)) iy = static_cast<int>(ny) - 1;
      if (iz >= static_cast<int>(nz)) iz = static_cast<int>(nz) - 1;

      bbxIdx = static_cast<unsigned>(ix) +
               nx * (static_cast<unsigned>(iy) +
                     ny * static_cast<unsigned>(iz));

      return true;
    }

    abort();
  }

  void elem::GetAMRRestrictionCandidateDataForLevelPair(
    const unsigned ilevel,
    const unsigned dim,
    const unsigned jMin,
    const unsigned jSize,
    MyMatrix<unsigned>& interfaceDofJ,
    MyMatrix<unsigned>& levelInterfaceSolidMarkJ,
    std::vector<MyMatrix<double>>& interfaceNodeCoordinatesJ,
    const std::vector<std::vector<double*>>& xMin,
    const std::vector<std::vector<double*>>& xMax,
    std::vector<std::vector<unsigned>>& jDof_r,
    std::vector<std::vector<unsigned>>& jSolidMark_r,
    std::vector<std::vector<std::vector<double>>>& jCoordinate_r
  ) {
    std::vector<std::vector<unsigned>> jDof_s(_nprocs);
    std::vector<std::vector<unsigned>> jSolidMark_s(_nprocs);
    std::vector<std::vector<std::vector<double>>> jCoordinate_s(_nprocs);

    std::vector<unsigned> size_s(_nprocs, 0u);
    std::vector<unsigned> size_r(_nprocs, 0u);

    for (unsigned lproc = 0u; lproc < _nprocs; lproc++) {
      jCoordinate_s[lproc].resize(dim);
    }

    jDof_r.assign(_nprocs, std::vector<unsigned>());
    jSolidMark_r.assign(_nprocs, std::vector<unsigned>());
    jCoordinate_r.assign(_nprocs, std::vector<std::vector<double>>(dim));

    for (unsigned lproc = 0u; lproc < _nprocs; lproc++) {

      jDof_s[lproc].reserve(jSize);
      jSolidMark_s[lproc].reserve(jSize);

      for (unsigned d = 0u; d < dim; d++) {
        jCoordinate_s[lproc][d].reserve(jSize);
      }

      if (jSize == 0u) continue;

      bool validBox = true;

      for (unsigned d = 0u; d < dim; d++) {
        if (xMin[lproc][ilevel][d] == DBL_MAX ||
            xMax[lproc][ilevel][d] == -DBL_MAX) {
          validBox = false;
          break;
        }
      }

      if (!validBox) continue;

      std::vector<char> visited(jSize, 0);

      for (unsigned k = interfaceDofJ.begin(); k < interfaceDofJ.end(); k++) {
        for (unsigned l = interfaceDofJ.begin(k); l < interfaceDofJ.end(k); l++) {

          const unsigned ldof = interfaceDofJ[k][l];
          const unsigned idx  = ldof - jMin;

          if (visited[idx]) continue;

          bool insideLproc = true;

          for (unsigned d = 0u; d < dim; d++) {
            const double x = interfaceNodeCoordinatesJ[d][k][l];

            if (x < xMin[lproc][ilevel][d] ||
                x > xMax[lproc][ilevel][d]) {
              insideLproc = false;
              break;
            }
          }

          if (!insideLproc) continue;

          visited[idx] = 1;

          jDof_s[lproc].push_back(ldof);
          jSolidMark_s[lproc].push_back(levelInterfaceSolidMarkJ[k][l]);

          for (unsigned d = 0u; d < dim; d++) {
            jCoordinate_s[lproc][d].push_back(interfaceNodeCoordinatesJ[d][k][l]);
          }
        }
      }

      size_s[lproc] = static_cast<unsigned>(jDof_s[lproc].size());
    }

    MPI_Alltoall(size_s.data(), 1, MPI_UNSIGNED,
                 size_r.data(), 1, MPI_UNSIGNED,
                 MPI_COMM_WORLD);

    for (unsigned lproc = 0u; lproc < _nprocs; ++lproc) {
      jDof_r[lproc].resize(size_r[lproc]);
      jSolidMark_r[lproc].resize(size_r[lproc]);

      for (unsigned d = 0u; d < dim; d++) {
        jCoordinate_r[lproc][d].resize(size_r[lproc]);
      }
    }

    std::vector<MPI_Request> reqs;
    reqs.reserve((2u + dim) * 2u * _nprocs);

    for (int p = 0; p < static_cast<int>(_nprocs); ++p) {

      if (p == static_cast<int>(_iproc)) {
        jDof_r[_iproc]       = jDof_s[_iproc];
        jSolidMark_r[_iproc] = jSolidMark_s[_iproc];

        for (unsigned d = 0u; d < dim; d++) {
          jCoordinate_r[_iproc][d] = jCoordinate_s[_iproc][d];
        }

        continue;
      }

      if (size_r[p] == 0u) continue;

      MPI_Request req;

      MPI_Irecv(jDof_r[p].data(),
                static_cast<int>(size_r[p]),
                MPI_UNSIGNED,
                p,
                200,
                MPI_COMM_WORLD,
                &req);
      reqs.push_back(req);

      MPI_Irecv(jSolidMark_r[p].data(),
                static_cast<int>(size_r[p]),
                MPI_UNSIGNED,
                p,
                201,
                MPI_COMM_WORLD,
                &req);
      reqs.push_back(req);

      for (unsigned d = 0u; d < dim; d++) {
        MPI_Irecv(jCoordinate_r[p][d].data(),
                  static_cast<int>(size_r[p]),
                  MPI_DOUBLE,
                  p,
                  300 + static_cast<int>(d),
                  MPI_COMM_WORLD,
                  &req);
        reqs.push_back(req);
      }
    }

    for (int p = 0; p < static_cast<int>(_nprocs); ++p) {

      if (p == static_cast<int>(_iproc)) continue;
      if (size_s[p] == 0u) continue;

      MPI_Request req;

      MPI_Isend(jDof_s[p].data(),
                static_cast<int>(size_s[p]),
                MPI_UNSIGNED,
                p,
                200,
                MPI_COMM_WORLD,
                &req);
      reqs.push_back(req);

      MPI_Isend(jSolidMark_s[p].data(),
                static_cast<int>(size_s[p]),
                MPI_UNSIGNED,
                p,
                201,
                MPI_COMM_WORLD,
                &req);
      reqs.push_back(req);

      for (unsigned d = 0u; d < dim; d++) {
        MPI_Isend(jCoordinate_s[p][d].data(),
                  static_cast<int>(size_s[p]),
                  MPI_DOUBLE,
                  p,
                  300 + static_cast<int>(d),
                  MPI_COMM_WORLD,
                  &req);
        reqs.push_back(req);
      }
    }

    if (!reqs.empty()) {
      MPI_Waitall(static_cast<int>(reqs.size()),
                  reqs.data(),
                  MPI_STATUSES_IGNORE);
    }
  }

  void elem::GetAMRRestrictionOld(Mesh * msh) {

    std::vector < std::map < unsigned,  std::map < unsigned, double  > > > &restriction = msh->GetAmrRestrictionMap();
    restriction.resize(3);

    std::vector < std::map < unsigned, bool > > &interfaceSolidMark = msh->GetAmrSolidMark();
    interfaceSolidMark.resize(3);

    std::vector < MyVector<unsigned> > interfaceElement; // iel = interfaceElement[ilevel][i]
    std::vector < MyMatrix<unsigned> > interfaceLocalDof; // ldof = interfaceElement[ilevel][i][j]
    std::vector < std::vector < MyMatrix<unsigned> > > interfaceDof; // gdof = interfaceDof[solType][ilevel][i][j]
    std::vector < std::vector < MyMatrix<unsigned> > > levelInterfaceSolidMark; // FSIdof = levelInterfaceSolidMark[solType][ilevel][i][j]
    std::vector < std::vector < MyMatrix< double > > > interfaceNodeCoordinates; // x[dim] = interfaceNodeCoordinates[ilevel][dim][i][j]

    interfaceElement.resize(_level + 1);
    interfaceLocalDof.resize(_level + 1);
    interfaceDof.resize(3);
    levelInterfaceSolidMark.resize(3);
    for (unsigned i = 0; i < 3; i++) {
      interfaceDof[i].resize(_level + 1);
      levelInterfaceSolidMark[i].resize(_level + 1);
    }
    interfaceNodeCoordinates.resize(_level + 1);
    unsigned dim = msh->GetDimension();

    for (unsigned ilevel = 0; ilevel <= _level; ilevel++) {
      //BEGIN interface element search
      interfaceElement[ilevel] = MyVector <unsigned> (_elementOwned);
      unsigned counter = 0;
      for (unsigned i = _elementLevel.begin(); i < _elementLevel.end(); i++) {
        if (ilevel == _elementLevel[i]) {
          for (unsigned j = _elementNearFace.begin(i); j < _elementNearFace.end(i); j++) {
            if (-1 == _elementNearFace[i][j]) {
              interfaceElement[ilevel][counter] = i;
              counter++;
              break;
            }
          }
        }
      }
      interfaceElement[ilevel].resize(counter);
      interfaceElement[ilevel].stack();
      //END interface element search

      //BEGIN interface node search
      std::vector< unsigned > offset = interfaceElement[ilevel].getOffset();
      interfaceLocalDof[ilevel] = MyMatrix <unsigned>(offset, NVE[0][2], UINT_MAX);
      for (unsigned i = interfaceElement[ilevel].begin(); i < interfaceElement[ilevel].end(); i++) {
        unsigned iel =  interfaceElement[ilevel][i];
        std::map <unsigned, bool> ldofs;
        for (unsigned jface = _elementNearFace.begin(iel); jface < _elementNearFace.end(iel); jface++) {
          if (-1 == _elementNearFace[iel][jface]) {
            for (unsigned k = 0; k < GetNFACENODES(GetElementType(iel), jface, 2); k++) {
              unsigned index = GetIG(GetElementType(iel), jface, k);
              ldofs[index] = true;
            }
          }
        }
        unsigned j = 0;
        for (std::map<unsigned, bool>::iterator it = ldofs.begin(); it != ldofs.end(); it++) {
          interfaceLocalDof[ilevel][i][j] = it->first;
          j++;
        }
      }
      interfaceLocalDof[ilevel].shrinkToFit(UINT_MAX);
      //END interface node search

      //BEGIN interface node dof global search, one for each soltype
      MyVector <unsigned> rowSize = interfaceLocalDof[ilevel].getRowSize();

      for (unsigned soltype = 0; soltype < 3; soltype++) {
        interfaceDof[soltype][ilevel] = MyMatrix< unsigned > (rowSize, UINT_MAX);
        levelInterfaceSolidMark[soltype][ilevel] = MyMatrix< unsigned > (rowSize, UINT_MAX);
        for (unsigned i = interfaceLocalDof[ilevel].begin(); i < interfaceLocalDof[ilevel].end(); i++) {
          unsigned iel = interfaceElement[ilevel][i];
          unsigned counter = 0;
          for (unsigned j = interfaceLocalDof[ilevel].begin(i); j < interfaceLocalDof[ilevel].end(i); j++) {
            unsigned jloc = interfaceLocalDof[ilevel][i][j];
            if (jloc < GetElementDofNumber(iel, soltype)) {
              unsigned jdof  = msh->GetSolutionDof(jloc, iel, soltype);
              interfaceDof[soltype][ilevel][i][counter] = jdof;
              unsigned jdof2  = msh->GetSolutionDof(jloc, iel, 2);
              levelInterfaceSolidMark[soltype][ilevel][i][counter] = msh->GetSolidMark(jdof2);
              counter++;
            }
            else {
              break;
            }
          }
        }
        interfaceDof[soltype][ilevel].shrinkToFit(UINT_MAX);
        levelInterfaceSolidMark[soltype][ilevel].shrinkToFit(UINT_MAX);
      }
      //END interface node dof global search, one for each soltype

      //BEGIN interface node coordinates
      interfaceNodeCoordinates[ilevel].resize(dim);
      for (unsigned k = 0; k < dim; k++) {
        interfaceNodeCoordinates[ilevel][k] = MyMatrix< double > (rowSize, 0.);
      }
      for (unsigned i = interfaceLocalDof[ilevel].begin(); i < interfaceLocalDof[ilevel].end(); i++) {
        unsigned iel = interfaceElement[ilevel][i];
        for (unsigned j = interfaceLocalDof[ilevel].begin(i); j < interfaceLocalDof[ilevel].end(i); j++) {
          unsigned jnode = interfaceLocalDof[ilevel][i][j];
          unsigned xDof  = msh->GetSolutionDof(jnode, iel, 2);
          for (unsigned k = 0; k < dim; k++) {
            interfaceNodeCoordinates[ilevel][k][i][j] = (*msh->_topology->_Sol[k])(xDof);
          }
        }
      }
      //END interface node coordinates search
    }

    const unsigned nLevels = _level + 1u;

    std::vector<double> xMinLocal(nLevels * dim, DBL_MAX);
    std::vector<double> xMaxLocal(nLevels * dim, -DBL_MAX);

    // pointer view: local [level][d]
    std::vector<double*> xMinL(nLevels);
    std::vector<double*> xMaxL(nLevels);

    for (unsigned level = 0; level < nLevels; level++) {
      xMinL[level] = xMinLocal.data() + level * dim;
      xMaxL[level] = xMaxLocal.data() + level * dim;
    }

    for (int ilevel = 0; ilevel < nLevels; ilevel++) {
      for (unsigned i = interfaceElement[ilevel].begin(); i < interfaceElement[ilevel].end(); i++) { //i-level element loop
        unsigned iel = interfaceElement[ilevel][i];
        for (unsigned j = 0; j < GetElementDofNumber(iel, 2u); j++) { // i-level i-elem node loop
          unsigned xdof  = msh->GetSolutionDof(j, iel, 2u);
          for (unsigned d = 0; d < dim; d++) {
            double x = (*msh->_topology->_Sol[d])(xdof);
            if(x < xMinL[ilevel][d]) xMinL[ilevel][d] = x;
            if(x > xMaxL[ilevel][d]) xMaxL[ilevel][d] = x;
          }
        }
      }
    }

    // Enlarge each valid local level bounding box by 1% on both sides.
    // Skip levels with no points, which still have DBL_MAX / -DBL_MAX.
    for (unsigned level = 0; level < nLevels; level++) {
      bool hasPoints = true;

      for (unsigned d = 0; d < dim; d++) {
        if (xMinL[level][d] == DBL_MAX || xMaxL[level][d] == -DBL_MAX) {
          hasPoints = false;
          break;
        }
      }

      if (!hasPoints) continue;

      for (unsigned d = 0; d < dim; d++) {
        const double length = xMaxL[level][d] - xMinL[level][d];
        const double delta = 0.01 * length;

        xMinL[level][d] -= delta;
        xMaxL[level][d] += delta;
      }
    }

    std::vector<double> xMinMemory(_nprocs * nLevels * dim);
    std::vector<double> xMaxMemory(_nprocs * nLevels * dim);

    MPI_Allgather(xMinLocal.data(), nLevels * dim, MPI_DOUBLE,
                  xMinMemory.data(), nLevels * dim, MPI_DOUBLE,
                  MPI_COMM_WORLD);

    MPI_Allgather(xMaxLocal.data(), nLevels * dim, MPI_DOUBLE,
                  xMaxMemory.data(), nLevels * dim, MPI_DOUBLE,
                  MPI_COMM_WORLD);

    // pointer view: global [proc][level][d]
    std::vector<std::vector<double*> > xMin(_nprocs);
    std::vector<std::vector<double*> > xMax(_nprocs);

    for (unsigned proc = 0; proc < _nprocs; proc++) {
      xMin[proc].resize(nLevels);
      xMax[proc].resize(nLevels);

      for (unsigned level = 0; level < nLevels; level++) {
        xMin[proc][level] =
          xMinMemory.data() + (proc * nLevels + level) * dim;

        xMax[proc][level] =
          xMaxMemory.data() + (proc * nLevels + level) * dim;
      }
    }

    for (unsigned soltype = 0; soltype < 3; soltype++) {

      //std::vector<bool> jlevelNodes;
      std::vector<bool> ilevelNodes;

      std::vector<double> xl(dim);
      std::vector<double> phi;
      std::vector<std::vector<double>> gradPhi(dim);
      std::vector < std::vector <double > > xv(dim);
      std::vector <double> xc(dim);
      std::vector < std::vector< double > > xe(dim);
      std::vector <double> xi(dim);
      std::vector < std::vector < std::vector <double > > > aP(3);

      /////////////////////////////

      std::vector<unsigned> jMinLevel(nLevels, UINT_MAX);
      std::vector<unsigned> jMaxLevel(nLevels, 0);
      std::vector<unsigned> jSizeLevel(nLevels, 0);
      std::vector<std::vector<unsigned>> jDof_s(_nprocs);
      std::vector<std::vector<unsigned>> jSolidMark_s(_nprocs);
      std::vector<std::vector<std::vector<double>>> jCoordinate_s(_nprocs);

      std::vector<std::vector<unsigned>> jDof_r(_nprocs);
      std::vector<std::vector<unsigned>> jSolidMark_r(_nprocs);
      std::vector<std::vector<std::vector<double>>> jCoordinate_r(_nprocs);

      std::vector<unsigned> size_s(_nprocs, 0);
      std::vector<unsigned> size_r(_nprocs, 0);

      for(unsigned lproc = 0u; lproc < _nprocs; lproc++) {
        jCoordinate_s[lproc].resize(dim);
        jCoordinate_r[lproc].resize(dim);
      }

      for (unsigned level = 0; level < nLevels; level++) {
        for (unsigned k = interfaceDof[soltype][level].begin(); k < interfaceDof[soltype][level].end(); k++) {
          for (unsigned l = interfaceDof[soltype][level].begin(k); l < interfaceDof[soltype][level].end(k); l++) {
            const unsigned ldof = interfaceDof[soltype][level][k][l];

            if (ldof < jMinLevel[level]) jMinLevel[level] = ldof;
            if (ldof > jMaxLevel[level]) jMaxLevel[level] = ldof;
          }
        }

        if (jMaxLevel[level] >= jMinLevel[level]) {
          jSizeLevel[level] = jMaxLevel[level] - jMinLevel[level] + 1u;
        }
      }

      for (int ilevel = 0; ilevel < _level; ilevel++) {
        for (int jlevel = ilevel + 1; jlevel <= _level; jlevel++) {

          const unsigned jMin  = jMinLevel[jlevel];
          const unsigned jSize = jSizeLevel[jlevel];

          for (unsigned lproc = 0; lproc < _nprocs; lproc++) {

            jDof_s[lproc].reserve(jSize);
            jSolidMark_s[lproc].reserve(jSize);
            jDof_s[lproc].clear();
            jSolidMark_s[lproc].clear();
            for(unsigned d = 0; d < dim; d++) {
              jCoordinate_s[lproc][d].reserve(jSize);
              jCoordinate_s[lproc][d].clear();
            }

            size_s[lproc] = 0u;

            if(jSize == 0) continue;

            // Skip processors with no valid ilevel bounding box.
            bool validBox = true;
            for (unsigned d = 0; d < dim; d++) {
              if (xMin[lproc][ilevel][d] == DBL_MAX || xMax[lproc][ilevel][d] == -DBL_MAX) {
                validBox = false;
                break;
              }
            }

            if (!validBox) continue;

            // One visited array per lproc, so each ldof is added at most once for this lproc.
            std::vector<char> visited(jSize, 0);

            for (unsigned k = interfaceDof[soltype][jlevel].begin(); k < interfaceDof[soltype][jlevel].end(); k++) {
              for (unsigned l = interfaceDof[soltype][jlevel].begin(k); l < interfaceDof[soltype][jlevel].end(k); l++) {
                const unsigned ldof = interfaceDof[soltype][jlevel][k][l];
                const unsigned idx = ldof - jMin;
                if (visited[idx]) continue;
                bool insideLproc = true;
                for (unsigned d = 0; d < dim; d++) {
                  const double x = interfaceNodeCoordinates[jlevel][d][k][l];
                  if (x < xMin[lproc][ilevel][d] || x > xMax[lproc][ilevel][d]) {
                    insideLproc = false;
                    break;
                  }
                }

                if (insideLproc) {

                  visited[idx] = 1;

                  size_s[lproc]++;

                  jDof_s[lproc].push_back(ldof);
                  jSolidMark_s[lproc].push_back(levelInterfaceSolidMark[soltype][jlevel][k][l]);

                  for(unsigned d = 0u; d < dim; d++) {
                    jCoordinate_s[lproc][d].push_back(interfaceNodeCoordinates[jlevel][d][k][l]);
                  }

                }
              }
            }
          }

          MPI_Alltoall(size_s.data(), 1, MPI_UNSIGNED,
                       size_r.data(), 1, MPI_UNSIGNED,
                       MPI_COMM_WORLD);

          // Resize receive buffers
          for(unsigned lproc = 0; lproc < _nprocs; ++lproc) {
            jDof_r[lproc].resize(size_r[lproc]);
            jSolidMark_r[lproc].resize(size_r[lproc]);
            for(unsigned d = 0u; d < dim; d++) {
              jCoordinate_r[lproc][d].resize(size_r[lproc]);
            }
          }

          // Nonblocking communication for jDof, jSolidMark, and jCoordinate.
          // Tags:
          //   200 = jDof
          //   201 = jSolidMark
          //   300 + d = jCoordinate[d]

          std::vector<MPI_Request> reqs;
          reqs.reserve((2u + dim) * 2u * _nprocs);

          // -----------------------------
          // Nonblocking receives first
          // -----------------------------
          for (int p = 0; p < static_cast<int>(_nprocs); ++p) {

            if (p == static_cast<int>(_iproc)) {
              jDof_r[_iproc]       = jDof_s[_iproc];
              jSolidMark_r[_iproc] = jSolidMark_s[_iproc];

              for (unsigned d = 0u; d < dim; d++) {
                jCoordinate_r[_iproc][d] = jCoordinate_s[_iproc][d];
              }
            }
            else {
              if (size_r[p] > 0u) {
                MPI_Request req;

                MPI_Irecv(jDof_r[p].data(),
                          static_cast<int>(size_r[p]),
                          MPI_UNSIGNED,
                          p,
                          200,
                          MPI_COMM_WORLD,
                          &req);
                reqs.push_back(req);

                MPI_Irecv(jSolidMark_r[p].data(),
                          static_cast<int>(size_r[p]),
                          MPI_UNSIGNED,
                          p,
                          201,
                          MPI_COMM_WORLD,
                          &req);
                reqs.push_back(req);

                for (unsigned d = 0u; d < dim; d++) {
                  MPI_Irecv(jCoordinate_r[p][d].data(),
                            static_cast<int>(size_r[p]),
                            MPI_DOUBLE,
                            p,
                            300 + static_cast<int>(d),
                            MPI_COMM_WORLD,
                            &req);
                  reqs.push_back(req);
                }
              }
            }
          }

          // -----------------------------
          // Nonblocking sends
          // -----------------------------
          for (int p = 0; p < static_cast<int>(_nprocs); ++p) {

            if (p != static_cast<int>(_iproc)) {
              if (size_s[p] > 0u) {
                MPI_Request req;

                MPI_Isend(jDof_s[p].data(),
                          static_cast<int>(size_s[p]),
                          MPI_UNSIGNED,
                          p,
                          200,
                          MPI_COMM_WORLD,
                          &req);
                reqs.push_back(req);

                MPI_Isend(jSolidMark_s[p].data(),
                          static_cast<int>(size_s[p]),
                          MPI_UNSIGNED,
                          p,
                          201,
                          MPI_COMM_WORLD,
                          &req);
                reqs.push_back(req);

                for (unsigned d = 0u; d < dim; d++) {
                  MPI_Isend(jCoordinate_s[p][d].data(),
                            static_cast<int>(size_s[p]),
                            MPI_DOUBLE,
                            p,
                            300 + static_cast<int>(d),
                            MPI_COMM_WORLD,
                            &req);
                  reqs.push_back(req);
                }
              }
            }
          }

          // -----------------------------
          // Complete communication
          // -----------------------------
          if (!reqs.empty()) {
            MPI_Waitall(static_cast<int>(reqs.size()),
                        reqs.data(),
                        MPI_STATUSES_IGNORE);
          }

          //////////////////////////

          for (unsigned lproc = 0; lproc < _nprocs; lproc++) {
            // interfaceDof[soltype][jlevel].broadcast(lproc);
            // levelInterfaceSolidMark[soltype][jlevel].broadcast(lproc);
            // for (unsigned d = 0; d < dim; d++) {
            //   interfaceNodeCoordinates[jlevel][d].broadcast(lproc);
            // }

            unsigned ilevelMinDof = UINT_MAX;
            unsigned ilevelMaxDof = 0;
            for (unsigned i = interfaceDof[soltype][ilevel].begin(); i < interfaceDof[soltype][ilevel].end(); i++) { //i-level element loop
              unsigned iel = interfaceElement[ilevel][i];
              for (unsigned j = 0; j < GetElementDofNumber(iel, soltype); j++) { // i-level i-elem node loop
                unsigned jdof  = msh->GetSolutionDof(j, iel, soltype);
                if (jdof < ilevelMinDof) ilevelMinDof = jdof;
                if (jdof > ilevelMaxDof) ilevelMaxDof = jdof;
              }
            }
            if (ilevelMaxDof > ilevelMinDof) ilevelNodes.resize(ilevelMaxDof - ilevelMinDof + 1u);
            else ilevelNodes.clear();

            unsigned jlevelMinDof = UINT_MAX;
            unsigned jlevelMaxDof = 0;
            for (unsigned k = interfaceDof[soltype][jlevel].begin(); k < interfaceDof[soltype][jlevel].end(); k++) { //j-level element loop
              for (unsigned l = interfaceDof[soltype][jlevel].begin(k); l < interfaceDof[soltype][jlevel].end(k); l++) { // j-level k-elem node loop
                unsigned ldof = interfaceDof[soltype][jlevel][k][l];
                if (ldof < jlevelMinDof) jlevelMinDof = ldof;
                if (ldof > jlevelMaxDof) jlevelMaxDof = ldof;
              }
            }
            // if (jlevelMaxDof > jlevelMinDof) jlevelNodes.resize(jlevelMaxDof - jlevelMinDof + 1u);
            // else jlevelNodes.clear();

            for (unsigned i = interfaceDof[soltype][ilevel].begin(); i < interfaceDof[soltype][ilevel].end(); i++) { //i-level element loop

              //candidateNodes.clear();
              std::fill(ilevelNodes.begin(), ilevelNodes.end(), false);
              // std::fill(jlevelNodes.begin(), jlevelNodes.end(), true);

              //elementNodes.clear();

              bool aPIsInitialized = false;

              unsigned iel = interfaceElement[ilevel][i];
              short unsigned ielType = _elementType[iel];
              basis* base = msh->GetBasis(ielType, soltype);

              for (unsigned j = 0; j < GetElementDofNumber(iel, soltype); j++) { // i-level i-elem node loop
                unsigned jdof  = msh->GetSolutionDof(j, iel, soltype);
                ilevelNodes[jdof - ilevelMinDof] = true;
              }

              msh->GetElementNodeCoordinates(xv, iel);
              unsigned ndofs = xv[0].size();
              double r2;
              GetConvexHullSphereRadiousSquare(xv, xc, r2, 0.01);

              GetBoundingBox(xv, xe, 0.01);

              // for (unsigned k = interfaceDof[soltype][jlevel].begin(); k < interfaceDof[soltype][jlevel].end(); k++) { //j-level element loop
              //   for (unsigned l = interfaceDof[soltype][jlevel].begin(k); l < interfaceDof[soltype][jlevel].end(k); l++) { // j-level k-elem node loop

              for(unsigned kl = 0; kl < jDof_r[lproc].size(); kl++ ) {
                //unsigned ldof = interfaceDof[soltype][jlevel][k][l];
                unsigned ldof = jDof_r[lproc][kl];

                if (ldof < ilevelMinDof || ldof > ilevelMaxDof || ilevelNodes[ldof - ilevelMinDof] == false) {

                  //if (jlevelNodes[ldof - jlevelMinDof] == true) {
                  double d2 = 0.;
                  for (int d = 0; d < dim; d++) {
                    //xl[d] = interfaceNodeCoordinates[jlevel][d][k][l];
                    xl[d] = jCoordinate_r[lproc][d][kl];
                    d2 += (xl[d] - xc[d]) * (xl[d] - xc[d]);
                  }
                  bool insideHull = true;
                  if (d2 > r2) {
                    insideHull = false;
                    //jlevelNodes[ldof - jlevelMinDof] = false;
                    continue;
                  }
                  for (unsigned d = 0; d < dim; d++) {
                    if (xl[d] < xe[d][0] || xl[d] > xe[d][1]) {
                      insideHull = false;
                      break;
                    }
                  }
                  if (insideHull) {

                    if (!aPIsInitialized) {
                      aPIsInitialized = true;
                      for (unsigned jtype = 0; jtype < 3; jtype++) {
                        ProjectNodalToPolynomialCoefficients(aP[jtype], xv, ielType, jtype) ;
                      }
                    }

                    {
                      GetClosestPointInReferenceElement(xv, xl, ielType, xi);
                      GetInverseMapping_fast(2u, ielType, aP, xl, xi, 100u, phi, gradPhi);
                    }

                    bool insideDomain = CheckIfPointIsInsideReferenceDomain(xi, ielType, 0.0001);
                    // if (insideDomain) {
                    //   for (unsigned j = interfaceDof[soltype][ilevel].begin(i); j < interfaceDof[soltype][ilevel].end(i); j++) {
                    //     unsigned jloc = interfaceLocalDof[ilevel][i][j];
                    //
                    //     double value = base->eval_phi(jloc, xi);
                    //
                    //     if (fabs(value) >= 1.0e-10) {
                    //       unsigned jdof = interfaceDof[soltype][ilevel][i][j];
                    //
                    //       auto& rowJ = restriction[soltype][jdof];
                    //
                    //       auto it = rowJ.find(jdof);
                    //       if (it == rowJ.end()) {
                    //         rowJ.emplace(jdof, 1.);
                    //         interfaceSolidMark[soltype][jdof] = levelInterfaceSolidMark[soltype][ilevel][i][j];
                    //       }
                    //
                    //       rowJ[ldof] = value;
                    //
                    //       auto& rowL = restriction[soltype][ldof];
                    //       rowL[ldof] = 10.;
                    //
                    //       interfaceSolidMark[soltype][ldof] = levelInterfaceSolidMark[soltype][jlevel][k][l];
                    //       jlevelNodes[ldof - jlevelMinDof] = true;
                    //     }
                    //   }
                    // }

                    if (insideDomain) {
                      for (unsigned j = interfaceDof[soltype][ilevel].begin(i); j < interfaceDof[soltype][ilevel].end(i); j++) {
                        unsigned jloc = interfaceLocalDof[ilevel][i][j];

                        double value = base->eval_phi(jloc, xi);

                        if (fabs(value) >= 1.0e-10) {
                          unsigned jdof = interfaceDof[soltype][ilevel][i][j];

                          auto& rowJ = restriction[soltype][jdof];

                          auto it = rowJ.find(jdof);
                          if (it == rowJ.end()) {
                            rowJ.emplace(jdof, 1.);
                            interfaceSolidMark[soltype][jdof] = levelInterfaceSolidMark[soltype][ilevel][i][j];
                          }

                          rowJ[ldof] = value;

                          auto& rowL = restriction[soltype][ldof];
                          rowL[ldof] = 10.;

                          // interfaceSolidMark[soltype][ldof] = levelInterfaceSolidMark[soltype][jlevel][k][l];
                          interfaceSolidMark[soltype][ldof] = jSolidMark_r[lproc][kl];
                        }
                      }
                    }

                    //jlevelNodes[ldof - jlevelMinDof] = false;
                  }
                  // }
                  // else {
                  //   jlevelNodes[ldof - jlevelMinDof] = false;
                  // }
                }
              }
              // }//end j-level element loop
            }
            // interfaceDof[soltype][jlevel].clearBroadcast();
            // levelInterfaceSolidMark[soltype][jlevel].clearBroadcast();
            // for (unsigned d = 0; d < dim; d++) {
            //   interfaceNodeCoordinates[jlevel][d].clearBroadcast();
            // }
          }
        }
      }

      NumericVector* pvector;
      pvector = NumericVector::build().release();
      pvector->init(_nprocs, 1, false, AUTOMATIC);

      unsigned counter = 1;
      while (counter != 0) {
        counter = 0;

        //BEGIN  saving the restriction object in parallel vectors and matrices

        MyVector <unsigned> rowSize(restriction[soltype].size(), 0);
        unsigned cnt1 = 0;
        for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
          rowSize[cnt1] = restriction[soltype][it1->first].size();
          cnt1++;
        }
        rowSize.stack();

        std::vector< unsigned > offset = rowSize.getOffset();

        MyVector <unsigned> masterNode(offset);
        MyMatrix <unsigned> slaveNodes(rowSize);
        MyMatrix <double> slaveNodesValues(rowSize);

        cnt1 = 0;
        for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
          masterNode[offset[_iproc] + cnt1] = it1->first;
          unsigned cnt2 = 0;
          for (std::map<unsigned, double> ::iterator it2 = restriction[soltype][it1->first].begin(); it2 != restriction[soltype][it1->first].end(); it2++) {
            slaveNodes[offset[_iproc] + cnt1][cnt2] = it2->first;
            slaveNodesValues[offset[_iproc] + cnt1][cnt2] = it2->second;
            cnt2++;
          }
          cnt1++;
        }
        //END  saving the restriction object in parallel vectors and matrices

        //BEGIN filling the restriction object with infos coming form the parallel vectors and matrices
        unsigned solutionOffset = msh->_dofOffset[soltype][_iproc];
        unsigned solutionOffsetp1 = msh->_dofOffset[soltype][_iproc + 1];
        for (unsigned lproc = 0; lproc < _nprocs; lproc++) {
          masterNode.broadcast(lproc);
          slaveNodes.broadcast(lproc);
          slaveNodesValues.broadcast(lproc);
          for (unsigned i = slaveNodes.begin(); i < slaveNodes.end(); i++) {
            unsigned inode = masterNode[i];
            if (inode >= solutionOffset && inode < solutionOffsetp1 && // inode belongs to _iproc
                restriction[soltype].find(inode) == restriction[soltype].end()) { // but inode is not set as master node of _iproc
              counter++;
              for (unsigned j = slaveNodes.begin(i); j < slaveNodes.end(i); j++) { //copy information for lproc to _iproc
                unsigned jnode = slaveNodes[i][j];
                restriction[soltype][inode][jnode] = slaveNodesValues[i][j];
              }
            }
            else { // either inode does not belong to _iproc or it was already defined as a master for _iproc
              if (restriction[soltype].find(inode) != restriction[soltype].end()) { // inode is already defined as master node in _iproc (either it does or does not belong to _iproc)
                for (unsigned j = slaveNodes.begin(i); j < slaveNodes.end(i); j++) { // loop on all the columns of restriction[lproc][inode]
                  unsigned jnode = slaveNodes[i][j];
                  double value = slaveNodesValues[i][j];
                  if (inode != jnode || value > 5.) { // if off-diagonal or hanging node for lproc
                    restriction[soltype][inode][jnode] =  value;
                  }
                  if (restriction[soltype].find(jnode) == restriction[soltype].end()) { // if jnode is not yet a master node for _iproc
                    counter++;
                    for (unsigned k = masterNode.begin(); k < masterNode.end(); k++) {
                      if (masterNode[k] == jnode) { // and if jnode is also a master node for lproc
                        for (unsigned l = slaveNodes.begin(k); l < slaveNodes.end(k); l++) {
                          unsigned lnode = slaveNodes[k][l];
                          restriction[soltype][jnode][lnode] = slaveNodesValues[k][l]; //copy the rule of lptoc into _iproc
                        }
                        break;
                      }
                    }
                  }
                }
              }
            }
          }
          masterNode.clearBroadcast();
          slaveNodes.clearBroadcast();
          slaveNodesValues.clearBroadcast();
        }
        //END filling the restriction object with infos coming form the parallel vectors and matrices

        pvector->set(_iproc, counter);
        pvector->close();
        counter = static_cast <unsigned>(floor(pvector->l1_norm() + 0.5));
      }
      delete pvector;

      //       for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
      //         unsigned inode = it1->first;
      //         if (restriction[soltype][inode][inode] > 5.) {
      //           if (restriction[soltype][inode].size() > 1) {
      //             for (std::map<unsigned, std::map<unsigned, double> >::iterator it2 = restriction[soltype].begin(); it2 != restriction[soltype].end(); it2++) {
      //               unsigned jnode = it2->first;
      //               if (jnode != inode && restriction[soltype][jnode].find(inode) != restriction[soltype][jnode].end()) {
      //                 double value =  restriction[soltype][jnode][inode];
      //                 for (std::map<unsigned, double> ::iterator it3 = restriction[soltype][inode].begin(); it3 != restriction[soltype][inode].end(); it3++) {
      //                   unsigned knode = it3->first;
      //                   if (knode != inode) {
      //                     restriction[soltype][jnode][knode] = it3->second * value;
      //                   }
      //                 }
      //               }
      //             }
      //           }
      //           restriction[soltype][inode].clear();
      //           restriction[soltype][inode][inode] = 0.;
      //         }
      //       }

      std::vector<std::vector < unsigned > > genealogy;
      std::vector<std::vector < double > > heredity;
      std::vector< unsigned > index;
      std::map < unsigned,  std::map < unsigned, double  > >  restrictionCopy = restriction[soltype];

      for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restrictionCopy.begin(); it1 != restrictionCopy.end(); it1++) { // loop all over master, hanging and master+hanging nodes

        genealogy.resize(1);
        heredity.resize(1);
        index.resize(1);

        genealogy[0].resize(1);
        heredity[0].resize(1);

        unsigned inode = it1->first;

        if (restrictionCopy[inode][inode] < 5.) { // only if a real master node

          // initialize master node genealogy and heredity at level 0
          genealogy[0][0] = inode;
          heredity[0][0] = 1.;
          index[0] = 0;

          restriction[soltype][inode].clear();
          restriction[soltype][inode][inode] = 1.;

          unsigned level = 1;
          while (level > 0) {

            // initialize master node genealogy and heredity at genemeric level

            unsigned father = genealogy[level - 1][index[level - 1]];

            genealogy.resize(level + 1);
            genealogy[level].reserve(restrictionCopy[ father ].size() - 1);
            genealogy[level].resize(0);

            heredity.resize(level + 1);
            heredity[level].reserve(restrictionCopy[ father ].size() - 1);
            heredity[level].resize(0);

            index.resize(level + 1);
            index[level] = 0;

            unsigned cnt  = 0;
            for (std::map <unsigned, double>::iterator it2 = restrictionCopy[ father ].begin(); it2 != restrictionCopy[ father ].end(); it2++) { // loop on all the father sons
              unsigned son = it2->first;
              bool alreadyFound = false;
              for (unsigned klevel = 0; klevel < level; klevel++) { // check if the son is in the previous genealogy
                for (unsigned k = 0; k < genealogy[klevel].size(); k++) {
                  if (genealogy[klevel][k] == son) alreadyFound = true;
                }
              }
              if (!alreadyFound) {  // if never found add the the restionction value in the master node line and zero the hanging node line
                genealogy[level].resize(genealogy[level].size() + 1);
                heredity[level].resize(heredity[level].size() + 1);

                genealogy[level][cnt] = son;
                heredity[level][cnt] = it2->second * heredity[level - 1][index[level - 1]];

                restriction[soltype][inode][son] += heredity[level][cnt];
                cnt++;

                restriction[soltype][son].clear();
                restriction[soltype][son][son] = 0.;
              }
            }

            if (cnt > 0) {
              level++;
            }
            else {
              bool test = true;
              while (test && level > 0) {
                index[level - 1]++;
                test = false;
                if (index[level - 1] == genealogy[level - 1].size()) {
                  level--;
                  test = true;
                }
              }
            }
          }
        }
        else {
          restriction[soltype][inode].clear();
          restriction[soltype][inode][inode] = 0.;
        }
      }

      MyVector <unsigned> InterfaceSolidMarkNode(interfaceSolidMark[soltype].size());
      MyVector <short unsigned> InterfaceSolidMarkValue(interfaceSolidMark[soltype].size());

      unsigned cnt = 0;
      for (std::map<unsigned, bool >::iterator it = interfaceSolidMark[soltype].begin(); it != interfaceSolidMark[soltype].end(); it++) {
        InterfaceSolidMarkNode[cnt] = it->first;
        InterfaceSolidMarkValue[cnt] = it->second;
        cnt++;
      }
      InterfaceSolidMarkNode.stack();
      InterfaceSolidMarkValue.stack();

      for (unsigned lproc = 0; lproc < _nprocs; lproc++) {
        InterfaceSolidMarkNode.broadcast(lproc);
        InterfaceSolidMarkValue.broadcast(lproc);
        for (unsigned i = InterfaceSolidMarkNode.begin(); i < InterfaceSolidMarkNode.end(); i++) {
          unsigned jnode = InterfaceSolidMarkNode[i];
          if (restriction[soltype].find(jnode) != restriction[soltype].end()) {
            interfaceSolidMark[soltype][jnode] = InterfaceSolidMarkValue[i];
          }
        }
        InterfaceSolidMarkNode.clearBroadcast();
        InterfaceSolidMarkValue.clearBroadcast();
      }

    }

  }

} //end namespace femus
