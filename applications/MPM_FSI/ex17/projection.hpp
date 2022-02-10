#include "Marker.hpp"

class Projection {
  public:
    Projection(MultiLevelSolution *mlSolM, MultiLevelSolution *mlSolB) {
      _mlSolM = mlSolM;
      _mlSolB = mlSolB;
      unsigned levelB = _mlSolB->_mlMesh->GetNumberOfLevels() - 1;
      Mesh *mshB   = _mlSolB->_mlMesh->GetLevel(levelB);
      _iproc  = mshB->processor_id();
      _nprocs = mshB->n_processors();
      _ielb.resize(_nprocs);
      for(unsigned kp = 0; kp < _nprocs; kp++) _ielb[kp].resize(0);
      _ielMInitialized = false;
    }
    ~Projection() {};

    void Init();
    void Allocate(const unsigned & dim, const unsigned & iproc, std::vector < unsigned > &pPntCnt);

    void SetNewmarkParameters(const double &beta, const double &gamma, const double &DT) {
      _beta = beta;
      _gamma = gamma;
      _DT = DT;
    }

    void FromMarkerToBackground();
    void FromBackgroundToMarker();

    const std::vector<std::vector<unsigned> > & GetIel() {
      return _ielb;
    }
    const std::vector < std::vector < std::vector <double > > > & GetXi() {
      return _xib;
    }
    const std::vector<std::vector<unsigned> > & GetMtype() {
      return _mtypeb;
    }
    const std::vector<std::vector<double> > & GetWeight() {
      return _weightb;
    }
    const std::vector < std::vector < std::vector <double > > > & GetX() {
      return _Xb;
    }
    const std::vector < std::vector < std::vector <double > > > & GetV() {
      return _Vb;
    }
    const std::vector < std::vector < std::vector <double > > > & GetA() {
      return _Ab;
    }
    const std::vector < std::vector<std::vector<double> > > & GetD() {
      return _Db;
    }
    const std::vector < std::vector<std::vector<double> > > & GetN() {
      return _Nb;
    }
    const std::vector < std::vector < std::vector < std::vector <double > > > > & GetGradD() {
      return _gradDb;
    }
    
    const std::vector < std::vector < std::vector < std::vector <double > > > > & GetF() {
      return _Fb;
    }

  private:

    double _gamma;
    double _beta;
    double _DT;

    unsigned _dim;
    unsigned _nprocs;
    unsigned _iproc;

    MultiLevelSolution *_mlSolM, *_mlSolB;

    std::vector < unsigned > _map; // sending coordinates

    std::vector < std::vector < std::vector <double > > > _Xm; // marker coordinates
    std::vector < std::vector < std::vector <double > > > _Vm; // marker velocity
    std::vector < std::vector < std::vector <double > > > _Am; // marker acceleration
    std::vector < std::vector<std::vector<double> > > _Dm;     // marker displacement
    std::vector < std::vector<std::vector<double> > > _Nm;     // marker normal
    std::vector < std::vector < std::vector < std::vector <double > > > > _gradDm; // marker grad D
    std::vector < std::vector < std::vector < std::vector <double > > > > _Fm; // marker grad D
    std::vector<std::vector<double> > _weightm; // marker weight
    std::vector<std::vector<unsigned> > _ielm; // marker elements
    std::vector<std::vector<unsigned> > _mtypem; // marker mtype

    std::vector < std::vector < std::vector <double > > > _xib; // background coordinates
    std::vector < std::vector < std::vector <double > > > _Xb; // background coordinates
    std::vector < std::vector < std::vector <double > > > _Vb; // background velocity
    std::vector < std::vector < std::vector <double > > > _Ab; // background acceleration
    std::vector < std::vector<std::vector<double> > > _Db;     // background displacement
    std::vector < std::vector<std::vector<double> > > _Nb;     // background normal
    std::vector < std::vector < std::vector < std::vector <double > > > > _gradDb; // background grad D
    std::vector < std::vector < std::vector < std::vector <double > > > > _Fb;
    std::vector<std::vector<double> > _weightb; // background weight
    std::vector<std::vector<unsigned> > _ielb; // background elements
    std::vector<std::vector<unsigned> > _mtypeb; // mbackground mtype

    bool _ielMInitialized;

    const char _Dname[3][3] = {"DX", "DY", "DZ"};
    const char _Vname[3][3] = {"VX", "VY", "VZ"};
    const char _Aname[3][3] = {"AX", "AY", "AZ"};
    const char _Nname[3][3] = {"NX", "NY", "NZ"};
    const char _gradDname[3][3][4] = {{"DXx", "DXy", "DXz"}, {"DYx", "DYy", "DYz"}, {"DZx", "DZy", "DZz"}};
    const char _Fname[3][3][4] = {{"F11", "F12", "F13"}, {"F21", "F22", "F23"}, {"F31", "F32", "F33"}};

};

////////////////////////////////////

void Projection::Allocate(const unsigned & dim, const unsigned & iproc, std::vector < unsigned > &pPntCnt) {
  _dim = dim;
  _iproc = iproc;
  _nprocs = pPntCnt.size() - 1;
  _Xm.resize(_nprocs);
  _Xb.resize(_nprocs);
  _xib.resize(_nprocs);
  _Vm.resize(_nprocs);
  _Vb.resize(_nprocs);
  _Am.resize(_nprocs);
  _Ab.resize(_nprocs);
  _Dm.resize(_nprocs);
  _Db.resize(_nprocs);

  _Nm.resize(_nprocs);
  _Nb.resize(_nprocs);

  _gradDm.resize(_nprocs);
  _gradDb.resize(_nprocs);

  _Fm.resize(_nprocs);
  _Fb.resize(_nprocs);

  _weightm.resize(_nprocs);
  _weightb.resize(_nprocs);

  _ielm.resize(_nprocs);
  _ielb.resize(_nprocs);
  _mtypem.resize(_nprocs);
  _mtypeb.resize(_nprocs);
  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    _Xm[kproc].resize(_dim);
    _Xb[kproc].resize(_dim);
    _xib[kproc].resize(_dim);
    _Vm[kproc].resize(_dim);
    _Vb[kproc].resize(_dim);
    _Am[kproc].resize(_dim);
    _Ab[kproc].resize(_dim);
    _Dm[kproc].resize(_dim);
    _Db[kproc].resize(_dim);
    _Nm[kproc].resize(_dim);
    _Nb[kproc].resize(_dim);
    _gradDm[kproc].resize(_dim);
    _gradDb[kproc].resize(_dim);
    _Fm[kproc].resize(_dim);
    _Fb[kproc].resize(_dim);

    unsigned size = pPntCnt[kproc + 1] - pPntCnt[kproc]; // size of the iproc marker grid belonging to the jproc background grid
    for(unsigned k = 0; k < _dim; k++) {
      _Xm[kproc][k].resize(size);
      _Vm[kproc][k].resize(size);
      _Am[kproc][k].resize(size);
      _Dm[kproc][k].resize(size);
      _Nm[kproc][k].resize(size);
      _gradDm[kproc][k].resize(_dim);
      _gradDb[kproc][k].resize(_dim);
      _Fm[kproc][k].resize(_dim);
      _Fb[kproc][k].resize(_dim);
      for(unsigned l = 0; l < _dim; l++) {
        _gradDm[kproc][k][l].resize(size);
        _Fm[kproc][k][l].resize(size);
      }
    }
    _weightm[kproc].resize(size);
    _ielm[kproc].resize(size);
    _mtypem[kproc].resize(size);
  }

  std::vector < MPI_Request > reqsSend(_nprocs) ;
  std::vector < MPI_Request > reqsRecv(_nprocs) ;

  std::vector<unsigned > vsize(_nprocs);
  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    MPI_Irecv(&vsize[kproc], 1, MPI_UNSIGNED, kproc, 0, MPI_COMM_WORLD, &reqsRecv[kproc]);
  }

  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    unsigned size = _ielm[kproc].size();
    MPI_Isend(&size, 1, MPI_UNSIGNED, kproc, 0, MPI_COMM_WORLD, &reqsSend[kproc]);
  }

  MPI_Status status;
  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    MPI_Wait(&reqsRecv[kproc], &status);
    MPI_Wait(&reqsSend[kproc], &status);
  }

  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    _weightb[kproc].resize(vsize[kproc]);
    _ielb[kproc].resize(vsize[kproc]);
    _mtypeb[kproc].resize(vsize[kproc]);
    for(unsigned k = 0; k < _dim; k++) {
      _Xb[kproc][k].resize(vsize[kproc]);
      _xib[kproc][k].resize(vsize[kproc]);
      _Vb[kproc][k].resize(vsize[kproc]);
      _Ab[kproc][k].resize(vsize[kproc]);
      _Db[kproc][k].resize(vsize[kproc]);
      _Nb[kproc][k].resize(vsize[kproc]);
      for(unsigned l = 0; l < _dim; l++) {
        _gradDb[kproc][k][l].resize(vsize[kproc]);
        _Fb[kproc][k][l].resize(vsize[kproc]);
      }
    }
  }
}

////////////////////////////////////

void Projection::Init() {

  unsigned levelM = _mlSolM->_mlMesh->GetNumberOfLevels() - 1;
  Solution *solM  = _mlSolM->GetSolutionLevel(levelM);
  Mesh     *mshM   = _mlSolM->_mlMesh->GetLevel(levelM);

  unsigned levelB = _mlSolB->_mlMesh->GetNumberOfLevels() - 1;
  Solution *solB  = _mlSolB->GetSolutionLevel(levelB);
  Mesh     *mshB   = _mlSolB->_mlMesh->GetLevel(levelB);


  unsigned iproc  = mshM->processor_id();
  unsigned nprocs  = mshM->n_processors();

  const unsigned dim = mshM->GetDimension();

  std::vector < unsigned > DIdx(dim);
  for(unsigned k = 0; k < dim; k++) {
    DIdx[k] = _mlSolM->GetIndex(&_Dname[k][0]);
  }

  unsigned ielIdx = _mlSolM->GetIndex("iel");

  unsigned solTypeM = 2;

  if(!_ielMInitialized) {
    std::vector <double> xm(dim);
    for(unsigned kproc = 0; kproc < nprocs; kproc++) {
      for(unsigned k = mshM->_dofOffset[solTypeM][kproc]; k < mshM->_dofOffset[solTypeM][kproc + 1]; k++) {

        xm.assign(dim, 0.);
        if(iproc == kproc) {
          for(unsigned j = 0; j < dim; j++) {
            xm[j] = (*mshM->_topology->_Sol[j])(k) + (*solM->_Sol[DIdx[j]])(k);
          }
        }
        MPI_Bcast(xm.data(), xm.size(), MPI_DOUBLE, kproc, PETSC_COMM_WORLD);

        Marker *gp = new Marker(xm, VOLUME, solB, solTypeM);
        unsigned kel = gp->GetMarkerElement();

        unsigned kelold = (*solM->_Sol[ielIdx])(k);

        if(kel == UINT_MAX) std::cout << "error";

        solM->_Sol[ielIdx]->set(k, kel);

        delete gp;

      }
    }
    solM->_Sol[ielIdx]->close();
  }

  clock_t time = clock();
  unsigned dof0 = mshM->_dofOffset[solTypeM][iproc];
  unsigned dof1 = mshM->_dofOffset[solTypeM][iproc + 1];
  unsigned size = dof1 - dof0; // number of markers iproc owns

  std::vector < unsigned > pSize(nprocs, 0);
  std::vector < unsigned > pPntCnt(nprocs + 1, 0); //process Marker counter
  _map.resize(size);

  std::vector<unsigned> ielp(size);
  std::vector<unsigned*> vec(ielp.size());

  for(unsigned i = 0; i < ielp.size(); i++) {
    ielp[i] = static_cast<unsigned>((*solM->_Sol[ielIdx])(dof0 + i));
    vec[i] = &ielp[i];
    unsigned kproc = mshB->IsdomBisectionSearch(ielp[i], 3);
    pSize[kproc]++;
  }

  std::sort(vec.begin(), vec.end(), [](const unsigned * a, const unsigned * b) {
    return *a < *b;
  });

  for(unsigned i = 0; i < _map.size(); i++) {
    _map[i] =  static_cast<unsigned>(vec[i] - &ielp[0]);
  }

  for(unsigned kproc = 1; kproc <= nprocs; kproc++) {
    pPntCnt[kproc] = pPntCnt[kproc - 1] + pSize[kproc - 1];
  }

//   for(unsigned i = 0; i < _map.size(); i++) _map[i] = i;
//   for(unsigned i = 0; i < _map.size(); i++) {
//     unsigned iel = (*solM->_Sol[ielIdx])(dof0 + _map[i]);
//     for(unsigned j = i + 1; j < _map.size(); j++) {
//       unsigned jel = (*solM->_Sol[ielIdx])(dof0 + _map[j]);
//       if(jel < iel) {
//         iel = jel;
//         std::swap(_map[j], _map[i]);
//         //_map[j] = _map[i];
//         //_map[i] = mapi;
//       }
//     }
//     unsigned kproc = mshB->IsdomBisectionSearch(iel, 3);
//     pPntCnt[kproc + 1] = i + 1;
//
//     //std::cout << (*solM->_Sol[ielIdx])(dof0 + _map[i]) << " ";
//
//     if((*solM->_Sol[ielIdx])(dof0 + _map[i]) != (*solM->_Sol[ielIdx])(dof0 + map1[i])) std::cout<< _map[i] << " " << map1[i]<<"   ";
//   }
//   for(unsigned kproc = 1; kproc <= nprocs; kproc++) {
//     if(pPntCnt[kproc] == 0)  pPntCnt[kproc] = pPntCnt[kproc - 1];
//   }
  std::cout << "sort time " << " = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl << std::flush;

  time = clock();
  this->Allocate(dim, iproc, pPntCnt);
  std::cout << "allocate time " << " = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl << std::flush;
}

////////////////////////////////////

void Projection::FromMarkerToBackground() {

  clock_t time = clock();
  this->Init(); // creates the mapping and allocate memory
  std::cout << "init time " << " = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl << std::flush;

  unsigned levelM = _mlSolM->_mlMesh->GetNumberOfLevels() - 1;
  Solution *solM  = _mlSolM->GetSolutionLevel(levelM);
  Mesh     *mshM   = _mlSolM->_mlMesh->GetLevel(levelM);

  std::vector < unsigned > DIdx(_dim);
  std::vector < unsigned > VIdx(_dim);
  std::vector < unsigned > AIdx(_dim);
  std::vector < unsigned > NIdx(_dim);
  std::vector < std::vector < unsigned > > gradDIdx(_dim);
  std::vector < std::vector < unsigned > > FIdx(_dim);
  for(unsigned k = 0; k < _dim; k++) {
    gradDIdx[k].resize(_dim);
    FIdx[k].resize(_dim);
    DIdx[k] = _mlSolM->GetIndex(&_Dname[k][0]);
    VIdx[k] = _mlSolM->GetIndex(&_Vname[k][0]);
    AIdx[k] = _mlSolM->GetIndex(&_Aname[k][0]);
    NIdx[k] = _mlSolM->GetIndex(&_Nname[k][0]);
    for(unsigned l = 0; l < _dim; l++) {
      gradDIdx[k][l] = _mlSolM->GetIndex(&_gradDname[k][l][0]);
      FIdx[k][l] = _mlSolM->GetIndex(&_Fname[k][l][0]);
    }
  }

  unsigned weightIdx = _mlSolM->GetIndex("weight");
  unsigned ielIdx = _mlSolM->GetIndex("iel");
  unsigned mtypeIdx = _mlSolM->GetIndex("mtype");

  unsigned solTypeM = 2;

  unsigned offset0 = mshM->_dofOffset[solTypeM][_iproc];
  unsigned offset1 = mshM->_dofOffset[solTypeM][_iproc + 1];
  unsigned size = offset1 - offset0;

  unsigned jproc = 0; // process on the background grid to which the i marker belongs
  unsigned nj = _ielm[0].size(); //offset on the send vectors

  for(unsigned i = 0; i < size; i++) { // local marker loop
    unsigned mapi = offset0 + _map[i]; //global vector mapped index

    unsigned j;
    for(unsigned kproc = jproc; kproc < _nprocs; kproc++) {
      if(i < nj)  {
        jproc = kproc;
        j = (i + _ielm[kproc].size()) - nj ;
        break;
      }
      else {
        nj += _ielm[kproc + 1].size();
      }
    }
    _weightm[jproc][j] = (*solM->_Sol[weightIdx])(mapi);
    _ielm[jproc][j] = static_cast <unsigned>(floor((*solM->_Sol[ielIdx])(mapi) + 0.1));
    _mtypem[jproc][j] = static_cast <unsigned>(floor((*solM->_Sol[mtypeIdx])(mapi) + 0.1));;
    for(unsigned k = 0; k < _dim; k++) {
      _Xm[jproc][k][j] = (*mshM->_topology->_Sol[k])(mapi) + (*solM->_Sol[DIdx[k]])(mapi);
      _Vm[jproc][k][j] = (*solM->_Sol[VIdx[k]])(mapi);
      _Am[jproc][k][j] = (*solM->_Sol[AIdx[k]])(mapi);
      _Nm[jproc][k][j] = (*solM->_Sol[NIdx[k]])(mapi);
      for(unsigned l = 0; l < _dim; l++) {
        _gradDm[jproc][k][l][j] = (*solM->_Sol[gradDIdx[k][l]])(mapi);
        _Fm[jproc][k][l][j] = (*solM->_Sol[FIdx[k][l]])(mapi);
      }
    }
  }

  std::vector<std::vector < MPI_Request >> reqsSend(_nprocs) ;
  std::vector<std::vector < MPI_Request >> reqsRecv(_nprocs) ;

  unsigned nOfMessages = 3 + 4 * _dim + 2 * _dim * _dim;
  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    reqsRecv[kproc].resize(nOfMessages);
    int cnt = 0;
    for(unsigned k = 0; k < _dim; k++) {
      MPI_Irecv(_Xb[kproc][k].data(), _Xb[kproc][k].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsRecv[kproc][cnt]);
      cnt++;
      MPI_Irecv(_Vb[kproc][k].data(), _Vb[kproc][k].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsRecv[kproc][cnt]);
      cnt++;
      MPI_Irecv(_Ab[kproc][k].data(), _Ab[kproc][k].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsRecv[kproc][cnt]);
      cnt++;
      MPI_Irecv(_Nb[kproc][k].data(), _Nb[kproc][k].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsRecv[kproc][cnt]);
      cnt++;
      for(unsigned l = 0; l < _dim; l++) {
        MPI_Irecv(_gradDb[kproc][k][l].data(), _gradDb[kproc][k][l].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsRecv[kproc][cnt]);
        cnt++;
        MPI_Irecv(_Fb[kproc][k][l].data(), _Fb[kproc][k][l].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsRecv[kproc][cnt]);
        cnt++;
      }
    }
    MPI_Irecv(_weightb[kproc].data(), _weightb[kproc].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsRecv[kproc][cnt]);
    cnt++;
    MPI_Irecv(_ielb[kproc].data(), _ielb[kproc].size(), MPI_UNSIGNED, kproc, cnt, MPI_COMM_WORLD, &reqsRecv[kproc][cnt]);
    cnt++;
    MPI_Irecv(_mtypeb[kproc].data(), _mtypeb[kproc].size(), MPI_UNSIGNED, kproc, cnt, MPI_COMM_WORLD, &reqsRecv[kproc][cnt]);
    cnt++;
  }

  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    reqsSend[kproc].resize(nOfMessages);
    int cnt = 0;
    for(unsigned k = 0; k < _dim; k++) {
      MPI_Isend(_Xm[kproc][k].data(), _Xm[kproc][k].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsSend[kproc][cnt]);
      cnt++;
      MPI_Isend(_Vm[kproc][k].data(), _Vm[kproc][k].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsSend[kproc][cnt]);
      cnt++;
      MPI_Isend(_Am[kproc][k].data(), _Am[kproc][k].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsSend[kproc][cnt]);
      cnt++;
      MPI_Isend(_Nm[kproc][k].data(), _Nm[kproc][k].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsSend[kproc][cnt]);
      cnt++;
      for(unsigned l = 0; l < _dim; l++) {
        MPI_Isend(_gradDm[kproc][k][l].data(), _gradDm[kproc][k][l].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsSend[kproc][cnt]);
        cnt++;
        MPI_Isend(_Fm[kproc][k][l].data(), _Fm[kproc][k][l].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsSend[kproc][cnt]);
        cnt++;
      }
    }
    MPI_Isend(_weightm[kproc].data(), _weightm[kproc].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsSend[kproc][cnt]);
    cnt++;
    MPI_Isend(_ielm[kproc].data(), _ielm[kproc].size(), MPI_UNSIGNED, kproc, cnt, MPI_COMM_WORLD, &reqsSend[kproc][cnt]);
    cnt++;
    MPI_Isend(_mtypem[kproc].data(), _mtypem[kproc].size(), MPI_UNSIGNED, kproc, cnt, MPI_COMM_WORLD, &reqsSend[kproc][cnt]);
    cnt++;
  }

  MPI_Status status;
  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    for(unsigned j = 0; j < nOfMessages; j++) {
      MPI_Wait(&reqsRecv[kproc][j], &status);
      MPI_Wait(&reqsSend[kproc][j], &status);
    }
  }

  //Flag the elements on the background grid and the get inverse mapping of the particles immersed in the background mesh
  unsigned levelB = _mlSolB->_mlMesh->GetNumberOfLevels() - 1;
  Solution *solB  = _mlSolB->GetSolutionLevel(levelB);
  Mesh     *mshB  = _mlSolB->_mlMesh->GetLevel(levelB);
  unsigned eflagIdx = _mlSolB->GetIndex("eflag");
  unsigned nflagIdx = _mlSolB->GetIndex("nflag");
  unsigned nodeType = _mlSolB->GetSolutionType(nflagIdx);
  unsigned solType = _mlSolB->GetSolutionType("DX");

  solB->_Sol[eflagIdx]->zero();
  solB->_Sol[nflagIdx]->zero();

  std::vector < std::vector < std::vector <double > > > aP(3);
  std::vector<std::vector<double>> vx(_dim);
  std::vector <double> x(_dim);
  std::vector <double> xi(_dim);
  std::vector < unsigned > im(_nprocs, 0);
  for(int iel = mshB->_elementOffset[_iproc]; iel < mshB->_elementOffset[_iproc + 1]; iel++) {
    bool ielIsInitialized = false;
    short unsigned ielType = mshB->GetElementType(iel);

    for(unsigned kp = 0; kp < _nprocs; kp++) { //loop on all the markers in ielB element, projected from all the processes of the marker mesh

      im[kp] = 0;
      while(im[kp] < _ielb[kp].size() && _ielb[kp][im[kp]] < iel) {
        im[kp]++;
      }

      while(im[kp] < _ielb[kp].size() && iel == _ielb[kp][im[kp]]) {
        if(_mtypeb[kp][im[kp]] == 1 && (*solB->_Sol[eflagIdx])(iel) != 1) {
          solB->_Sol[eflagIdx]->set(iel, 1);
          unsigned nDofu  = mshB->GetElementDofNumber(iel, nodeType);  // number of solution element dofs
          for(unsigned i = 0; i < nDofu; i++) {
            unsigned idof = mshB->GetSolutionDof(i, iel, nodeType);
            solB->_Sol[nflagIdx]->set(idof, 1);
          }
        }
        else if(_mtypeb[kp][im[kp]] == 2 && (*solB->_Sol[eflagIdx])(iel) != 1) {
          solB->_Sol[eflagIdx]->set(iel, 2);
        }
        else if(_mtypeb[kp][im[kp]] == 0 && (*solB->_Sol[eflagIdx])(iel) == 0) {
          solB->_Sol[eflagIdx]->set(iel, 0.5);
        }

        if(!ielIsInitialized) {
          ielIsInitialized = true;
          unsigned nDofs = mshB->GetElementDofNumber(iel, solType);
          for(unsigned k = 0; k < _dim; k++) {
            vx[k].resize(nDofs);
            for(unsigned i = 0; i < nDofs; i++) {
              unsigned idofX = mshB->GetSolutionDof(i, iel, 2);
              vx[k][i] = (*mshB->_topology->_Sol[k])(idofX);
            }
          }
          for(unsigned jtype = 0; jtype <= solType; jtype++) {
            ProjectNodalToPolynomialCoefficients(aP[jtype], vx, ielType, jtype) ;
          }
        }


        for(unsigned k = 0; k < _dim; k++) {
          x[k] = _Xb[kp][k][im[kp]];
        }

        GetClosestPointInReferenceElement(vx, x, ielType, xi);
        bool inverseMapping = GetInverseMapping(solType, ielType, aP, x, xi, 100);
        if(!inverseMapping) {
          std::cout << "InverseMapping failed at " << iel << " " << im[kp] << std::endl;
        }
        for(unsigned k = 0; k < _dim; k++) {
          _xib[kp][k][im[kp]] = xi[k];
        }
        im[kp]++;
      }
    }
  }
  solB->_Sol[eflagIdx]->close();
  solB->_Sol[nflagIdx]->close();
}

void Projection::FromBackgroundToMarker() {

  unsigned levelB = _mlSolB->_mlMesh->GetNumberOfLevels() - 1;
  Solution *solB  = _mlSolB->GetSolutionLevel(levelB);
  Mesh     *mshB  = _mlSolB->_mlMesh->GetLevel(levelB);

  const char Dname[3][3] = {"DX", "DY", "DZ"};
  const char Vname[3][3] = {"VX", "VY", "VZ"};
  const char Aname[3][3] = {"AX", "AY", "AZ"};

  std::vector < unsigned > DIdxB(_dim);
  for(unsigned k = 0; k < _dim; k++) {
    DIdxB[k] = _mlSolM->GetIndex(&Dname[k][0]);
  }

  unsigned solTypeB = _mlSolB->GetSolutionType(DIdxB[0]);

  std::vector<std::vector<double>> vx(_dim);
  std::vector<std::vector<double>> D(_dim);
  std::vector<double> phi;
  std::vector<double> gradPhi;
  double weight;

  std::vector < unsigned > im(_nprocs, 0);
  std::vector<double> xi(_dim);
  for(int iel = mshB->_elementOffset[_iproc]; iel < mshB->_elementOffset[_iproc + 1]; iel++) {
    bool ielIsInitialized = false;
    short unsigned ielType = mshB->GetElementType(iel);
    unsigned nDofs;
    for(unsigned kp = 0; kp < _nprocs; kp++) {
      im[kp] = 0;
      while(im[kp] < _ielb[kp].size() && _ielb[kp][im[kp]] < iel) {
        im[kp]++;
      }

      while(im[kp] < _ielb[kp].size() && iel == _ielb[kp][im[kp]]) {
        if(!ielIsInitialized) {
          ielIsInitialized = true;
          nDofs = mshB->GetElementDofNumber(iel, solTypeB);
          for(unsigned k = 0; k < _dim; k++) {
            vx[k].resize(nDofs);
            D[k].resize(nDofs);
            for(unsigned i = 0; i < nDofs; i++) {
              unsigned idofX = mshB->GetSolutionDof(i, iel, 2);
              vx[k][i] = (*mshB->_topology->_Sol[k])(idofX);
              unsigned idof = mshB->GetSolutionDof(i, iel, solTypeB);
              D[k][i] = (*solB->_Sol[DIdxB[k]])(idof);
            }
          }
        }
        for(unsigned k = 0; k < _dim; k++) {
          xi[k] = _xib[kp][k][im[kp]];
        }
        mshB->_finiteElement[ielType][solTypeB]->Jacobian(vx, xi, weight, phi, gradPhi);

        for(unsigned k = 0; k < _dim; k++) {
          _Db[kp][k][im[kp]] = 0.;
          for(unsigned i = 0; i < nDofs; i++) {
            _Db[kp][k][im[kp]] += phi[i] * D[k][i];
          }
        }


        std::vector<double> Fold(pow(_dim, 2));
        for(unsigned k = 0; k < _dim; k++) {
          for(unsigned j = 0; j < _dim; j++) {
            Fold[k * _dim + j] = _Fb[kp][k][j][im[kp]];
          }
        }

        std::vector<double> Fnew(pow(_dim, 2));
        for(unsigned k = 0; k < _dim; k++) {
          for(unsigned j = 0; j < _dim; j++) {
            Fnew[k * _dim + j] = (k == j) ? 1. : 0.;
            for(unsigned i = 0; i < nDofs; i++) {
              Fnew[k * _dim + j] +=  gradPhi[i * _dim + j] * D[k][i];
            }
          }
        }

        for(unsigned k = 0; k < _dim; k++) {
          for(unsigned j = 0; j < _dim; j++) {
            _Fb[kp][k][j][im[kp]] = 0.;
            for(unsigned l = 0; l < _dim; l++) {
              _Fb[kp][k][j][im[kp]] += Fnew[k * _dim + l] * Fold[l * _dim + j];
            }
          }
        }

        im[kp]++;
      }
    }
  }
  for(unsigned k = 0; k < _dim; k++) { //reset the background displacement
    solB->_Sol[DIdxB[k]]->zero();
  }


  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    for(unsigned jp = 0; jp < _nprocs; jp++) {
      unsigned np;
      if(_iproc == kproc) np = _ielb[jp].size();
      MPI_Bcast(&np, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);
      for(unsigned im = 0; im < np; im++) {
        unsigned iel;
        std::vector < double > xp(_dim);
        if(_iproc == kproc) {
          iel = _ielb[jp][im];
          for(unsigned k = 0; k < _dim; k++) {
            xp[k] = _Xb[jp][k][im] + _Db[jp][k][im];
          }
        }
        MPI_Bcast(&iel, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);
        MPI_Bcast(xp.data(), xp.size(), MPI_DOUBLE, kproc, PETSC_COMM_WORLD);

        Marker gp = Marker(xp, VOLUME, solB, 2, iel);
        if(_iproc == kproc) {
          _ielb[jp][im] = gp.GetMarkerElement();
          if(_ielb[jp][im] == UINT_MAX) std::cout << "error";
        }
      }
    }
  }



  unsigned solTypeM = 2;

  std::vector<std::vector < MPI_Request >> reqsSend(_nprocs) ;
  std::vector<std::vector < MPI_Request >> reqsRecv(_nprocs) ;

  unsigned nOfMessages = _dim + _dim * _dim + 1;
  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    reqsRecv[kproc].resize(nOfMessages);
    int cnt = 0;
    for(unsigned k = 0; k < _dim; k++) {
      MPI_Irecv(_Dm[kproc][k].data(), _Dm[kproc][k].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsRecv[kproc][cnt]);
      cnt++;
      for(unsigned j = 0; j < _dim; j++) {
        MPI_Irecv(_Fm[kproc][k][j].data(), _Fm[kproc][k][j].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsRecv[kproc][cnt]);
        cnt++;
      }
    }
    MPI_Irecv(_ielm[kproc].data(), _ielm[kproc].size(), MPI_UNSIGNED, kproc, cnt, MPI_COMM_WORLD, &reqsRecv[kproc][cnt]);
    cnt++;
  }

  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    reqsSend[kproc].resize(nOfMessages);
    int cnt = 0;
    for(unsigned k = 0; k < _dim; k++) {
      MPI_Isend(_Db[kproc][k].data(), _Db[kproc][k].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsSend[kproc][cnt]);
      cnt++;
      for(unsigned j = 0; j < _dim; j++) {
        MPI_Isend(_Fb[kproc][k][j].data(), _Fb[kproc][k][j].size(), MPI_DOUBLE, kproc, cnt, MPI_COMM_WORLD, &reqsSend[kproc][cnt]);
        cnt++;
      }
    }
    MPI_Isend(_ielb[kproc].data(), _ielb[kproc].size(), MPI_UNSIGNED, kproc, cnt, MPI_COMM_WORLD, &reqsSend[kproc][cnt]);
    cnt++;
  }

  MPI_Status status;
  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    for(unsigned j = 0; j < nOfMessages; j++) {
      MPI_Wait(&reqsRecv[kproc][j], &status);
      MPI_Wait(&reqsSend[kproc][j], &status);
    }
  }

 
  
  unsigned levelM = _mlSolM->_mlMesh->GetNumberOfLevels() - 1;
  Solution *solM  = _mlSolM->GetSolutionLevel(levelM);
  Mesh     *mshM   = _mlSolM->_mlMesh->GetLevel(levelM);

  std::vector < unsigned > DIdxM(_dim);
  std::vector < unsigned > VIdxM(_dim);
  std::vector < unsigned > AIdxM(_dim);
  std::vector < std::vector < unsigned > > FIdxM(_dim);
  for(unsigned k = 0; k < _dim; k++) {
    FIdxM[k].resize(_dim);  
    DIdxM[k] = _mlSolM->GetIndex(&Dname[k][0]);
    VIdxM[k] = _mlSolM->GetIndex(&Vname[k][0]);
    AIdxM[k] = _mlSolM->GetIndex(&Aname[k][0]);
    for(unsigned j = 0; j < _dim; j++) {
      FIdxM[k][j] = _mlSolM->GetIndex(&_Fname[k][j][0]);
    }
  }
  unsigned ielIdx = _mlSolM->GetIndex("iel");

  unsigned offset0 = mshM->_dofOffset[solTypeM][_iproc];
  unsigned i = 0;
  for(unsigned jproc = 0; jproc < _nprocs; jproc++) {
    for(unsigned j = 0; j < _Dm[jproc][0].size(); j++) {
      unsigned mapi = offset0 + _map[i];

      solM->_Sol[ielIdx]->set(mapi, _ielm[jproc][j]);

      for(unsigned k = 0; k < _dim; k++) {

        double Xold = (*solM->_Sol[DIdxM[k]])(mapi);
        double Vold = (*solM->_Sol[VIdxM[k]])(mapi);
        double Aold = (*solM->_Sol[AIdxM[k]])(mapi);

        double Dnew = _Dm[jproc][k][j];

        double Xnew = Xold + Dnew;
        double Anew = Dnew / (_beta * _DT * _DT) - Vold / (_beta * _DT) - Aold * (0.5 - _beta) / _beta;
        double Vnew = Vold + (Aold * (1. - _gamma) + Anew * _gamma) * _DT;

//         std::cout.precision(12);
//         if(mapi == 0) std::cout << Dnew << " " << Vnew << " " << Anew << " " << Vold << " "<< Aold <<std::endl;

        solM->_Sol[DIdxM[k]]->set(mapi, Xnew);
        solM->_Sol[VIdxM[k]]->set(mapi, Vnew);
        solM->_Sol[AIdxM[k]]->set(mapi, Anew);

        for(unsigned l = 0; l < _dim; l++) {
          solM->_Sol[FIdxM[k][l]]->set(mapi, _Fm[jproc][k][l][j]);
        }

      }
      i++;
    }
  }

  solM->_Sol[ielIdx]->close();
  for(unsigned k = 0; k < _dim; k++) {
    solM->_Sol[DIdxM[k]]->close();
    solM->_Sol[VIdxM[k]]->close();
    solM->_Sol[AIdxM[k]]->close();
    for(unsigned j = 0; j < _dim; j++) {
      solM->_Sol[FIdxM[k][j]]->close();
    }
  }

  _ielMInitialized = true;
  UpdateMeshQuantities(_mlSolM);

}

