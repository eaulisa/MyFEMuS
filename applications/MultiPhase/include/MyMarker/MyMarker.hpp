/*=========================================================================

 Program: FEMuS
 Module: Marker
 Authors: Eugenio Aulisa and Giacomo Capodaglio

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_ism_BarbiMarker_hpp__
#define __femus_ism_BarbiMarker_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MarkerTypeEnum.hpp"
#include "ParallelObject.hpp"
#include "Mesh.hpp"
#include "vector"
#include "map"
#include "MyVector.hpp"

namespace femus {

  class MyMarker : public ParallelObject {
    public:

      MyMarker() {};

      MyMarker(const std::vector < double > &x, Solution *sol, const unsigned &solType,
               const unsigned &elem = UINT_MAX, const double &s1 = 0.) {

        _x = x;
        _solType = solType;
        _dim = sol->GetMesh()->GetDimension();

        if(elem == UINT_MAX) { //parallel search
          ParallelElementSearch(true, UINT_MAX, sol, s1);
        }
        else { //try first a serial search starting from the given elem

          _elem = elem;
          unsigned previousElem = elem;
          _mproc = GetMarkerProc(sol); //based on _elem we identify the process the search should start from

          unsigned preMproc = _mproc;
          if(_iproc == preMproc) {
            // careful this function can change _mproc, only in _iproc = preMproc, if _iproc believes that the marker does not belong to it
            SerialElementSearch(previousElem, sol, s1);
          }
          // previousElem and _elem, changed in _iproc =  preMproc, and are broadcast again
          MPI_Bcast(& _elem, 1, MPI_UNSIGNED, preMproc, PETSC_COMM_WORLD);
          MPI_Bcast(& previousElem, 1, MPI_UNSIGNED, preMproc, PETSC_COMM_WORLD);

          if(_elem != UINT_MAX) { // if the search in preMproc did not bring us outside the domain
            _mproc = GetMarkerProc(sol); //
            if(_mproc != preMproc) {  //if the search moved outside _preProc domain we call the global search
              // this is a parallel wrapper to serial search
              ParallelSerialElementSearch(previousElem, preMproc, sol, s1);
            }
          }
          else {// the particle is outside the domain
            //then check with a real parallel search
            ParallelElementSearch(true, UINT_MAX, sol, s1);
          }
        }

//         if(_iproc == _mproc) {
//           if(_elem != UINT_MAX) {
//             std::vector < std::vector < std::vector < std::vector < double > > > >aX;
//             FindLocalCoordinates(_solType, aX, true, sol, s1);
//           }
//         }
//         else {
//           std::vector < double > ().swap(_x);
//         }

        if(_iproc != _mproc) {
          std::vector < double > ().swap(_x);
        }

      };

      unsigned GetElement() {
        return _elem;
      }

      bool ParallelElementSearch(const std::vector < double > &x, Solution *sol, const unsigned & solType,
                                 const unsigned &elem = UINT_MAX, const double &s1 = 0.) {
        _x = x;
        _solType = solType;
        _dim = sol->GetMesh()->GetDimension();

        if(elem == UINT_MAX) { //parallel search
          ParallelElementSearch(true, UINT_MAX, sol, s1);
        }
        else { //try first a serial search starting from the given elem
          _elem = elem;
          unsigned previousElem = elem;
          _mproc = GetMarkerProc(sol); //based on _elem we identify the process the search should start from

          //std::cout <<previousElem<<" "<< _elem <<" "<< _mproc<<"\t";

          unsigned preMproc = _mproc;
          if(_iproc == preMproc) {
            // careful this function can change _mproc, only in _iproc = preMproc, if _iproc believes that the marker does not belong to it
            SerialElementSearch(previousElem, sol, s1);
          }
          // previousElem and _elem, changed in _iproc =  preMproc, and are broadcast again
          MPI_Bcast(& _elem, 1, MPI_UNSIGNED, preMproc, PETSC_COMM_WORLD);
          MPI_Bcast(& previousElem, 1, MPI_UNSIGNED, preMproc, PETSC_COMM_WORLD);

          //std::cout << previousElem<<" "<< _elem <<" "<< _mproc<<"\t";

          if(_elem != UINT_MAX) { // if the search in preMproc did not bring us outside the domain
            _mproc = GetMarkerProc(sol); //
            if(_mproc != preMproc) {  //if the search moved outside _preProc domain we call the global search
              // this is a parallel wrapper to serial search
              ParallelSerialElementSearch(previousElem, preMproc, sol, s1);
              //std::cout << _elem <<" "<< _mproc<<"\n";
            }
          }
          else {// the particle is outside the domain
            //then check with a real parallel search
            ParallelElementSearch(true, UINT_MAX, sol, s1);
          }
        }
        if(_iproc != _mproc) {
          std::vector < double > ().swap(_x);
        }
        return (_elem == UINT_MAX) ? false : true;
      }

      bool ParallelElementSearchWithInverseMapping(const std::vector < double > &x, Solution *sol, const unsigned & solType,
          const unsigned &elem = UINT_MAX, const double &s1 = 0.) {
        bool elemFound = ParallelElementSearch(x, sol, solType, elem, s1);
        if(_iproc == _mproc && elemFound) {
          FindLocalCoordinates(_solType, _aX, true, sol, s1);
        }
        return elemFound;
      }


      bool SerialElementSearch(const std::vector < double > &x, Solution *sol, const unsigned & solType,
                               const unsigned &elem, const double &s1 = 0.) {

        _x = x;
        _solType = solType;
        _dim = sol->GetMesh()->GetDimension();
        _elem = elem;
        _mproc = GetMarkerProc(sol);

        if(_iproc != _mproc) {
          return false;
        }
        else {
          unsigned previousElem = elem;
          return SerialElementSearch(previousElem, sol, s1);
        }
      }

      bool SerialElementSearchWithInverseMapping(const std::vector < double > &x, Solution *sol, const unsigned & solType,
          const unsigned &elem, const double &s1 = 0.) {
        bool elementFound = SerialElementSearch(x, sol, solType, elem, s1);
        if(elementFound == false) return false;
        else {
          FindLocalCoordinates(_solType, _aX, true, sol, s1);
          return true;
        }

      }

      void ClearElement() {
        _elem = UINT_MAX;
      }


      std::vector<double> GetIprocLocalCoordinates(const std::vector < double > &x, Solution *sol, const unsigned & solType, const unsigned &elem, const double &s1 = 0.) {
        _x = x;
        _solType = solType;
        _dim = sol->GetMesh()->GetDimension();
        bool newElement = (_elem != elem) ? true : false;
        _elem = elem;
        _mproc = GetMarkerProc(sol);
        FindLocalCoordinates(_solType, _aX, newElement, sol, s1);
        return _xi;
      }


      std::vector<double> GetIprocLocalCoordinates() {
        return _xi;
      }

      unsigned GetProc() {
        return _mproc;
      }

    private:
      double GetMeshCoordinates(Solution *sol, const unsigned &k, const unsigned &i, const double &s) {
        if(!sol->GetIfFSI()) {
          return (*sol->GetMesh()->_topology->_Sol[k])(i);
        }
        else {
          const char varname[3][3] = {"DX", "DY", "DZ"};
          unsigned solIndex = sol->GetIndex(&varname[k][0]);
          return (*sol->GetMesh()->_topology->_Sol[k])(i)
                 + (1. - s) * (*sol->_SolOld[solIndex])(i)
                 + s * (*sol->_Sol[solIndex])(i);
        }
      }

      unsigned GetMarkerProc(Solution *sol) {
        _mproc = (_elem == UINT_MAX) ? 0 : sol->GetMesh()->IsdomBisectionSearch(_elem, 3);
        return _mproc;
      }

      bool SerialElementSearch(unsigned & initialElem, Solution * sol, const double & s);
      void ParallelElementSearch(const bool & useInitialSearch, const unsigned & initialElem, Solution * sol, const double & s);
      void ParallelSerialElementSearch(unsigned & previousElem, const unsigned & previousMproc, Solution * sol, const double & s);

      void FindLocalCoordinates(const unsigned & solVType, std::vector < std::vector < std::vector < std::vector < double > > > > &aX,
                                const bool & pcElemUpdate, Solution * sol, const double & s);

      unsigned GetNextElement2D(const unsigned & iel, const unsigned & previousElem, Solution * sol, const double & s);
      unsigned GetNextElement3D(const unsigned & iel, const unsigned & previousElem, Solution * sol, const double & s);
      unsigned GetNextElement3D(const unsigned & iel, const std::vector< unsigned > &searchHistory, Solution * sol, const double & s);
      int FastForward(const unsigned & currentElem, const unsigned & previousElem, Solution * sol, const double & s);

      std::vector < double > _x; // global coordinates
      std::vector < double > _xi; // local coordinates
      unsigned _solType; // FEM interpolation: linear, serendipity, or quadratic

      unsigned _elem; // element that owns the marker, UINT_MAX if outside the domain
      unsigned _dim;
      unsigned _mproc; //process that owns the marker


      // support vectors
      std::vector < std::vector < std::vector < std::vector < double > > > >_aX;
      static const double _localCentralNode[6][3];
      static const double _a[4][4][4];
      static const double _b[4][4];
      static const double _c[4][4];
  };
  
} //end namespace femus

#include "./MyMarker.cpp"

#endif


