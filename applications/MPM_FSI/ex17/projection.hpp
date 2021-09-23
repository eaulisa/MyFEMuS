

class Projection {
  public:
    Projection(MultiLevelSolution *mlSolM, MultiLevelSolution *mlSolB) {
      _mlSolM = mlSolM;
      _mlSolB = mlSolB;
    }

    void Init();
    void Allocate(const unsigned & dim, const unsigned & iproc, std::vector < unsigned > &pPntCnt);

    void SetNewmarkParameters(const double &beta,const double &gamma,const double &DT){
      _beta = beta;
      _gamma = gamma;
      _DT = DT;
    }
    
    void FromMarkerToBackground();
    void FromBackgroundToMarker();
    void FakeMovement();
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
    std::vector<std::vector<double> > _weightm; // marker weight
    std::vector<std::vector<unsigned> > _ielm; // marker elements
    std::vector<std::vector<unsigned> > _mtypem; // marker mtype

    std::vector < std::vector < std::vector <double > > > _Xb; // background coordinates
    std::vector < std::vector < std::vector <double > > > _Vb; // background velocity
    std::vector < std::vector < std::vector <double > > > _Ab; // background acceleration
    std::vector < std::vector<std::vector<double> > > _Db;     // background displacement
    std::vector < std::vector<std::vector<double> > > _Nb;     // background normal
    std::vector < std::vector < std::vector < std::vector <double > > > > _gradDb; // background grad D
    std::vector<std::vector<double> > _weightb; // background weight
    std::vector<std::vector<unsigned> > _ielb; // background elements
    std::vector<std::vector<unsigned> > _mtypeb; // mbackground mtype

};

////////////////////////////////////

void Projection::Allocate(const unsigned & dim, const unsigned & iproc, std::vector < unsigned > &pPntCnt) {
  _dim = dim;
  _iproc = iproc;
  _nprocs = pPntCnt.size() - 1;
  _Xm.resize(_nprocs);
  _Xb.resize(_nprocs);
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

  _weightm.resize(_nprocs);
  _weightb.resize(_nprocs);

  _ielm.resize(_nprocs);
  _ielb.resize(_nprocs);
  _mtypem.resize(_nprocs);
  _mtypeb.resize(_nprocs);
  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    _Xm[kproc].resize(_dim);
    _Xb[kproc].resize(_dim);
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
    unsigned size = pPntCnt[kproc + 1] - pPntCnt[kproc]; // size of the iproc marker grid belonging to the jproc background grid
    for(unsigned k = 0; k < _dim; k++) {
      _Xm[kproc][k].resize(size);
      _Vm[kproc][k].resize(size);
      _Am[kproc][k].resize(size);
      _Dm[kproc][k].resize(size);
      _Nm[kproc][k].resize(size);
      _gradDm[kproc][k].resize(_dim);
      _gradDb[kproc][k].resize(_dim);
      for(unsigned l = 0; l < _dim; l++) {
        _gradDm[kproc][k][l].resize(size);
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
      _Vb[kproc][k].resize(vsize[kproc]);
      _Ab[kproc][k].resize(vsize[kproc]);
      _Db[kproc][k].resize(vsize[kproc]);
      _Nb[kproc][k].resize(vsize[kproc]);
      for(unsigned l = 0; l < _dim; l++) {
        _gradDb[kproc][k][l].resize(vsize[kproc]);
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
  const char Dname[3][3] = {"DX", "DY", "DZ"};

  std::vector < unsigned > DIdx(dim);
  for(unsigned k = 0; k < dim; k++) {
    DIdx[k] = _mlSolM->GetIndex(&Dname[k][0]);
  }

  unsigned ielIdx = _mlSolM->GetIndex("iel");

  unsigned solType = 2;

  std::vector <double> xm(dim);
  for(unsigned kproc = 0; kproc < nprocs; kproc++) {
    for(unsigned k = mshM->_dofOffset[solType][kproc]; k < mshM->_dofOffset[solType][kproc + 1]; k++) {

      xm.assign(dim, 0.);
      if(iproc == kproc) {
        for(unsigned j = 0; j < dim; j++) {
          xm[j] = (*mshM->_topology->_Sol[j])(k) + (*solM->_Sol[DIdx[j]])(k);
        }
      }
      MPI_Bcast(xm.data(), xm.size(), MPI_DOUBLE, kproc, PETSC_COMM_WORLD);

      Marker *gp = new Marker(xm, VOLUME, solB, solType);
      unsigned kel = gp->GetMarkerElement();

      solM->_Sol[ielIdx]->set(k, kel);

      delete gp;

    }
  }

  solM->_Sol[ielIdx]->close();

  unsigned dof0 = mshM->_dofOffset[solType][iproc];
  unsigned dof1 = mshM->_dofOffset[solType][iproc + 1];
  unsigned size = dof1 - dof0; // number of markers iproc owns

  std::vector < unsigned > pPntCnt(nprocs + 1, 0); //process Marker counter
  //std::vector<unsigned> &map = projection.GetMap();
  _map.resize(size);

  unsigned iel0 = UINT_MAX - 1;

  for(unsigned i = 0; i < _map.size(); i++) _map[i] = i;
  for(unsigned i = 0; i < _map.size(); i++) {
    unsigned iel = (*solM->_Sol[ielIdx])(dof0 + _map[i]);
    for(unsigned j = i + 1; j < _map.size(); j++) {
      unsigned jel = (*solM->_Sol[ielIdx])(dof0 + _map[j]);
      if(jel < iel) {
        iel = jel;
        unsigned mapi = _map[j];
        _map[j] = _map[i];
        _map[i] = mapi;
      }
    }
    unsigned kproc = mshB->IsdomBisectionSearch(iel, 3);
    pPntCnt[kproc + 1] = i + 1;
  }
  for(unsigned kproc = 1; kproc <= nprocs; kproc++) {
    if(pPntCnt[kproc] == 0)  pPntCnt[kproc] = pPntCnt[kproc - 1];
  }
  this->Allocate(dim, iproc, pPntCnt);
}

////////////////////////////////////

void Projection::FromMarkerToBackground() {

  unsigned levelM = _mlSolM->_mlMesh->GetNumberOfLevels() - 1;
  Solution *solM  = _mlSolM->GetSolutionLevel(levelM);
  Mesh     *mshM   = _mlSolM->_mlMesh->GetLevel(levelM);

  const char Dname[3][3] = {"DX", "DY", "DZ"};
  const char Vname[3][3] = {"VX", "VY", "VZ"};
  const char Aname[3][3] = {"AX", "AY", "AZ"};
  const char Nname[3][3] = {"NX", "NY", "NZ"};
  const char gradDname[3][3][4] = {{"DXx", "DXy", "DXz"}, {"DYx", "DYy", "DYz"}, {"DZx", "DZy", "DZz"}};

  std::vector < unsigned > DIdx(_dim);
  std::vector < unsigned > VIdx(_dim);
  std::vector < unsigned > AIdx(_dim);
  std::vector < unsigned > NIdx(_dim);
  std::vector < std::vector < unsigned > > gradDIdx(_dim);
  for(unsigned k = 0; k < _dim; k++) {
    gradDIdx[k].resize(_dim);
    DIdx[k] = _mlSolM->GetIndex(&Dname[k][0]);
    VIdx[k] = _mlSolM->GetIndex(&Vname[k][0]);
    AIdx[k] = _mlSolM->GetIndex(&Aname[k][0]);
    NIdx[k] = _mlSolM->GetIndex(&Nname[k][0]);
    for(unsigned l = 0; l < _dim; l++) {
      gradDIdx[k][l] = _mlSolM->GetIndex(&gradDname[k][l][0]);
    }
  }

  unsigned weightIdx = _mlSolM->GetIndex("weight");
  unsigned ielIdx = _mlSolM->GetIndex("iel");
  unsigned mtypeIdx = _mlSolM->GetIndex("mtype");

  unsigned solType = 2;

  unsigned offset0 = mshM->_dofOffset[solType][_iproc];
  unsigned offset1 = mshM->_dofOffset[solType][_iproc + 1];
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
      }
    }
  }

  std::vector<std::vector < MPI_Request >> reqsSend(_nprocs) ;
  std::vector<std::vector < MPI_Request >> reqsRecv(_nprocs) ;

  unsigned nOfMessages = 3 + 4 * _dim + _dim * _dim;
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

}

void Projection::FromBackgroundToMarker() {

  unsigned levelM = _mlSolM->_mlMesh->GetNumberOfLevels() - 1;
  Solution *solM  = _mlSolM->GetSolutionLevel(levelM);
  Mesh     *mshM   = _mlSolM->_mlMesh->GetLevel(levelM);

  const char Dname[3][3] = {"DX", "DY", "DZ"};
  const char Vname[3][3] = {"VX", "VY", "VZ"};
  const char Aname[3][3] = {"AX", "AY", "AZ"};

  std::vector < unsigned > DIdx(_dim);
  std::vector < unsigned > VIdx(_dim);
  std::vector < unsigned > AIdx(_dim);
  for(unsigned k = 0; k < _dim; k++) {
    DIdx[k] = _mlSolM->GetIndex(&Dname[k][0]);
    VIdx[k] = _mlSolM->GetIndex(&Vname[k][0]);
    AIdx[k] = _mlSolM->GetIndex(&Aname[k][0]);
  }
//
//   unsigned ielIdx = _mlSolM->GetIndex("iel");
//
  unsigned solType = 2;
//
//   unsigned offset0 = mshM->_dofOffset[solType][_iproc];
//   unsigned offset1 = mshM->_dofOffset[solType][_iproc + 1];
//   unsigned size = offset1 - offset0;
//
//   unsigned jproc = 0; // process on the background grid to which the i marker belongs
//   unsigned nj = _ielm[0].size(); //offset on the send vectors
//
//   for(unsigned i = 0; i < size; i++) { // local marker loop
//     unsigned mapi = offset0 + _map[i]; //global vector mapped index
//
//     unsigned j;
//     for(unsigned kproc = jproc; kproc < _nprocs; kproc++) {
//       if(i < nj)  {
//         jproc = kproc;
//         j = (i + _ielm[kproc].size()) - nj ;
//         break;
//       }
//       else {
//         nj += _ielm[kproc + 1].size();
//       }
//     }
//     _ielm[jproc][j] = (*solM->_Sol[ielIdx])(mapi);
//     for(unsigned k = 0; k < _dim; k++) {
//       _Xm[jproc][k][j] = (*mshM->_topology->_Sol[k])(mapi) + (*solM->_Sol[DIdx[k]])(mapi);
//     }
//   }

  std::vector<std::vector < MPI_Request >> reqsSend(_nprocs) ;
  std::vector<std::vector < MPI_Request >> reqsRecv(_nprocs) ;

  unsigned nOfMessages = _dim;
  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    reqsRecv[kproc].resize(nOfMessages);
    for(unsigned k = 0; k < _dim; k++) {
      MPI_Irecv(_Dm[kproc][k].data(), _Dm[kproc][k].size(), MPI_DOUBLE, kproc, k, MPI_COMM_WORLD, &reqsRecv[kproc][k]);
    }
  }

  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    reqsSend[kproc].resize(nOfMessages);
    for(unsigned k = 0; k < _dim; k++) {
      MPI_Isend(_Db[kproc][k].data(), _Db[kproc][k].size(), MPI_DOUBLE, kproc, k, MPI_COMM_WORLD, &reqsSend[kproc][k]);
    }
  }

  MPI_Status status;
  for(unsigned kproc = 0; kproc < _nprocs; kproc++) {
    for(unsigned j = 0; j < nOfMessages; j++) {
      MPI_Wait(&reqsRecv[kproc][j], &status);
      MPI_Wait(&reqsSend[kproc][j], &status);
    }
  }


  unsigned offset0 = mshM->_dofOffset[solType][_iproc];
  unsigned i = 0;
  for(unsigned jproc = 0; jproc < _nprocs; jproc++) {
    for(unsigned j = 0; j < _Dm[jproc][0].size(); j++) {
      unsigned mapi = offset0 + _map[i];
      for(unsigned k = 0; k < _dim; k++) {

        double D = _Dm[jproc][k][j];
        double Xold = (*solM->_Sol[DIdx[k]])(mapi);
        double Vold = (*solM->_Sol[VIdx[k]])(mapi);
        double Aold = (*solM->_Sol[AIdx[k]])(mapi);

        double Xnew = Xold + D;
        double Anew = D / (_beta * _DT * _DT) - Vold / (_beta * _DT) - Aold * (0.5 - _beta) / _beta;
        double Vnew = Vold + (Aold * (1. - _gamma) + Anew * _gamma) * _DT;
        
        //std::cout<<Vnew <<" ";

        solM->_Sol[DIdx[k]]->set(mapi, Xnew);
        solM->_Sol[VIdx[k]]->set(mapi, Vnew);
        solM->_Sol[AIdx[k]]->set(mapi, Anew);

      }
      i++;
    }
  }
  for(unsigned k = 0; k < _dim; k++) {
    solM->_Sol[DIdx[k]]->close();
    solM->_Sol[VIdx[k]]->close();
    solM->_Sol[AIdx[k]]->close();
  }
}


void Projection::FakeMovement() {

  for(unsigned jproc = 0; jproc < _nprocs; jproc++) {
    for(unsigned j = 0; j < _Db[jproc][0].size(); j++) {

      double x = _Xb[jproc][0][j];
      _Db[jproc][0][j] = 0.1 * ((1. - x / 5.) * x / 5 * 0.5 + x / 5 * x / 5. * (-0.5));


      _Db[jproc][1][j] = -0.1;
    }
  }








}
