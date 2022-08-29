

const double WEIGHT[6][27] = {
  {
    1. / 64., 1. / 64., 1. / 64., 1. / 64.,
    1. / 64., 1. / 64., 1. / 64., 1. / 64.,
    1. / 32., 1. / 32., 1. / 32., 1. / 32.,
    1. / 32., 1. / 32., 1. / 32., 1. / 32.,
    1. / 32., 1. / 32., 1. / 32., 1. / 32.,
    1. / 16., 1. / 16., 1. / 16., 1. / 16., 1. / 16., 1. / 16., 1. / 8
  }, //hex
  {
    1. / 32., 1. / 32., 1. / 32., 1. / 32.,
    7. / 48., 7. / 48., 7. / 48., 7. / 48., 7. / 48., 7. / 48.
  }, //tet
  {
    1. / 48., 1. / 48., 1. / 48., 1. / 48., 1. / 48., 1. / 48.,
    1. / 16., 1. / 16., 1. / 16., 1. / 16., 1. / 16., 1. / 16.,
    1. / 24., 1. / 24., 1. / 24., 1. / 8., 1. / 8., 1. / 8.
  }, //wedge
  {
    0.0625, 0.0625, 0.0625, 0.0625,
    0.125, 0.125, 0.125, 0.125, 0.25
  },
  {
    1. / 12., 1. / 12., 1. / 12.,
    0.25, 0.25, 0.25, 0.
  },
  {0.25, 0.25, .5}
} ;

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int facename, const double t) {
  bool test = (1 == facename)? 1 : 0; //dirichlet
  value = 0.;
  return test;
}

double InitVariableDX(const std::vector < double >& x) {
  double um = 0.2;
  double  value = (1. - x[0] / 5.) * x[0] / 5 * 0.5 + x[0] / 5 * x[0] / 5. * (-0.5);
  return value;
}

double InitVariableDY(const std::vector < double >& x) {
  //  double um = 0.2;
  //  double  value = (1. - x[0] / 5.) * x[0] / 5. * 0.5 + x[0] / 5. * x[0] / 5. * (-0.5);
  double  value = -0.1;
  return value;
}

double InitVariableVY(const std::vector < double >& x) {
  double  value = -0.1;
  return value;
}

void FlagElements(MultiLevelMesh& mlMesh, const unsigned &layers);

void InitializeMarkerVariables(MultiLevelSolution &mlSol);
void UpdateMeshQuantities(MultiLevelSolution *mlSol);
void BuildInvariants(MultiLevelSolution& mlSol);



void InitializeMarkerVariables(MultiLevelSolution &mlSol) {

  unsigned dim = mlSol._mlMesh->GetDimension();

  FEOrder femOrder = SECOND;

  //add variables to mlSol
  mlSol.AddSolution("DX", LAGRANGE, femOrder, 0, true);
  mlSol.AddSolution("DY", LAGRANGE, femOrder, 0, true);
  if(dim == 3) mlSol.AddSolution("DZ", LAGRANGE, femOrder, 0, true);

  mlSol.AddSolution("VX", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("VY", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("VZ", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("AX", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("AY", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("AZ", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("DXx", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("DXy", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("DXz", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("DYx", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("DYy", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("DYz", LAGRANGE, femOrder, 0, false);

  if(dim == 3) mlSol.AddSolution("DZx", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("DZy", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("DZz", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("F11", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("F12", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("F13", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("F21", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("F22", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("F23", LAGRANGE, femOrder, 0, false);

  if(dim == 3) mlSol.AddSolution("F31", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("F32", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("F33", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("mtype", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("NX", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("NY", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("NZ", LAGRANGE, femOrder, 0, false);



  mlSol.AddSolution("weight", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("area", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("iel", LAGRANGE, femOrder, 0, false); // for each node it stores the background element
  mlSol.AddSolution("dist", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("kernel", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("d2max", LAGRANGE, femOrder, 0, false); // for each node it stores the maximum distance2 in the support


  mlSol.Initialize("All");
  //mlSol.Initialize("DX", InitVariableDX);
  //mlSol.Initialize("VY", InitVariableVY);

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("DX", "Steady");
  if(dim > 1) mlSol.GenerateBdc("DY", "Steady");
  if(dim > 2) mlSol.GenerateBdc("DZ", "Steady");

  BuildInvariants(mlSol);
}

void BuildInvariants(MultiLevelSolution& mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel(level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level);
  unsigned iproc  = msh->processor_id();

  const unsigned dim = msh->GetDimension();

  unsigned solType = 2;

  std::vector < double > phi;
  std::vector < double > gradPhi;
  std::vector <std::vector < double> > vx(dim);

  std::vector < unsigned> idof;

  unsigned d2maxIdx = mlSol.GetIndex("d2max");
  unsigned mtypeIdx = mlSol.GetIndex("mtype");

  const char Dname[3][3] = {"DX", "DY", "DZ"};

  std::vector < unsigned > DIdx(dim);
  for(unsigned k = 0; k < dim; k++) {
    DIdx[k] = sol->GetIndex(&Dname[k][0]);
    //(*sol->_Bdc[DIdx[k]]) = 2.;
  }

  sol->_Sol[d2maxIdx]->zero();
  sol->_Sol[mtypeIdx]->zero();

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned ielMat = msh->GetElementMaterial(iel);

    short unsigned ielt = msh->GetElementType(iel);
    unsigned nDofs = msh->GetElementDofNumber(iel, solType);
    for(unsigned  k = 0; k < dim; k++) {
      vx[k].resize(nDofs);
    }
    idof.resize(nDofs);

    for(unsigned i = 0; i < nDofs; i++) {
      idof[i] = msh->GetSolutionDof(i, iel, solType);
      for(unsigned  k = 0; k < dim; k++) {
        vx[k][i] = (*msh->_topology->_Sol[k])(idof[i]);
      }
    }

    if(ielMat == 4) {
      for(unsigned i = 0; i < nDofs; i++) {
        double d2max = (*sol->_Sol[d2maxIdx])(idof[i]);
        for(unsigned j = 0; j < nDofs; j++) {
          double d2maxj = 0.;
          for(unsigned k = 0; k < dim; k++) {
            d2maxj += (vx[k][i] - vx[k][j]) * (vx[k][i] - vx[k][j]);
          }
          if(d2maxj > d2max) {
            d2max = d2maxj;
            sol->_Sol[d2maxIdx]->set(idof[i], d2max);
          }
        }
        sol->_Sol[mtypeIdx]->set(idof[i], 2);

        for(unsigned k = 0; k < dim; k++) {
          sol->_Bdc[DIdx[k]]->set(idof[i], 0.);
        }

      }
    }
  }
  sol->_Sol[d2maxIdx]->closeWithMaxValue();
  sol->_Sol[mtypeIdx]->close();

  for(unsigned k = 0; k < dim; k++) {
    sol->_Bdc[DIdx[k]]->close();
  }

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned ielMat = msh->GetElementMaterial(iel);
    if(ielMat == 2) {

      unsigned nDofs = msh->GetElementDofNumber(iel, solType);

      for(unsigned i = 0; i < nDofs; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solType);
        if((*sol->_Sol[mtypeIdx])(idof) == 2) {
          sol->_Sol[mtypeIdx]->set(idof, 1);
        }
      }
    }
  }
  sol->_Sol[mtypeIdx]->close();


  const char Fname[3][3][4] = {{"F11", "F12", "F13"}, {"F21", "F22", "F23"}, {"F31", "F32", "F33"}};
  std::vector < std::vector < unsigned > > FIdx(dim);
  for(unsigned k = 0; k < dim; k++) {
    FIdx[k].resize(dim);
    for(unsigned j = 0; j < dim; j++) {
      FIdx[k][j] = mlSol.GetIndex(&Fname[k][j][0]);
    }
  }
  for(unsigned k = 0; k < dim; k++) {
    for(unsigned j = 0; j < dim; j++) {
      (*sol->_Sol[FIdx[k][j]]) = (j == k) ?  1. : 0.;
    }
  }
  for(unsigned k = 0; k < dim; k++) {
    for(unsigned j = 0; j < dim; j++) {
      sol->_Sol[FIdx[k][j]]->close();
    }
  }
}

void UpdateMeshQuantities(MultiLevelSolution *mlSol) {

  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol->GetSolutionLevel(level);
  Mesh     *msh   = mlSol->_mlMesh->GetLevel(level);
  unsigned iproc  = msh->processor_id();
  unsigned nprocs  = msh->n_processors();

  const unsigned dim = msh->GetDimension();

  unsigned solType = 2;

  std::vector < double > phi;
  std::vector < double > gradPhi;
  std::vector <std::vector < double> > vx(dim);
  std::vector <std::vector < double> > vxD(dim);
  std::vector <std::vector < double> > solD(dim);

  std::vector < unsigned> idof;
  std::vector < double > d2max;
  std::vector < unsigned > mtype;

  // solution and coordinate variables
  const char Dname[3][3] = {"DX", "DY", "DZ"};
  const char gradDname[3][3][4] = {{"DXx", "DXy", "DXz"}, {"DYx", "DYy", "DYz"}, {"DZx", "DZy", "DZz"}};

  std::vector < unsigned > DIdx(dim);
  std::vector < std::vector < unsigned > > gradDIdx(dim);
  for(unsigned k = 0; k < dim; k++) {
    DIdx[k] = mlSol->GetIndex(&Dname[k][0]);
    gradDIdx[k].resize(dim);
    for(unsigned j = 0; j < dim; j++) {
      gradDIdx[k][j] = mlSol->GetIndex(&gradDname[k][j][0]);
    }
  }
  unsigned weightIdx = mlSol->GetIndex("weight");
  unsigned kernelIdx = mlSol->GetIndex("kernel");
  unsigned d2maxIdx = mlSol->GetIndex("d2max");
  unsigned mtypeIdx = mlSol->GetIndex("mtype");

  //return;

  sol->_Sol[weightIdx]->zero();
  sol->_Sol[kernelIdx]->zero();
  for(unsigned k = 0; k < dim; k++) {
    for(unsigned j = 0; j < dim; j++) {
      sol->_Sol[gradDIdx[k][j]]->zero();
    }
  }

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = msh->GetElementType(iel);
    unsigned ielMat = msh->GetElementMaterial(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);
    for(unsigned  k = 0; k < dim; k++) {
      solD[k].resize(nDofs);
      vx[k].resize(nDofs);
      vxD[k].resize(nDofs);
    }
    idof.resize(nDofs);
    d2max.resize(nDofs);
    mtype.resize(nDofs);

    for(unsigned i = 0; i < nDofs; i++) {
      idof[i] = msh->GetSolutionDof(i, iel, solType);
      mtype[i] = (*sol->_Sol[mtypeIdx])(idof[i]);
      d2max[i] = (*sol->_Sol[d2maxIdx])(idof[i]);
      for(unsigned  k = 0; k < dim; k++) {
        solD[k][i] = (*sol->_Sol[DIdx[k]])(idof[i]);
        vx[k][i] = (*msh->_topology->_Sol[k])(idof[i]);
        vxD[k][i] = vx[k][i] + solD[k][i];
      }
    }

    for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {
      double jac;
      msh->_finiteElement[ielt][solType]->Jacobian(vx, ig, jac, phi, gradPhi);

      std::vector< std::vector< double > > gradSolDg(dim);
      std::vector< double > xg(dim, 0.);
      for(unsigned  k = 0; k < dim; k++) { //solution
        gradSolDg[k].assign(dim, 0.);
        for(unsigned i = 0; i < nDofs; i++) { //node
          xg[k] += vx[k][i] * phi[i];
          for(unsigned j = 0; j < dim; j++) { // derivative
            gradSolDg[k][j] += solD[k][i] * gradPhi[dim * i + j];
          }
        }
      }
      for(unsigned i = 0; i < nDofs; i++) { //node

        if(ielMat == 4) { // || mtype[i] == 0) {
          double d2 = 0.;
          for(unsigned  k = 0; k < dim; k++) { //solution
            d2 += (vx[k][i] - xg[k]) * (vx[k][i] - xg[k]);
          }
          if(d2max[i] < d2) std::cout << "e";

          sol->_Sol[kernelIdx]->add(idof[i], jac * pow(d2max[i] - d2, 2));
          for(unsigned  k = 0; k < dim; k++) {
            for(unsigned j = 0; j < dim; j++) {
              sol->_Sol[gradDIdx[k][j]]->add(idof[i], gradSolDg[k][j] * jac * pow(d2max[i] - d2, 2));
            }
          }
        }
      }
    }
    double ielArea = 0.;
    for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {
      double jac;
      msh->_finiteElement[ielt][solType]->Jacobian(vxD, ig, jac, phi, gradPhi);
      ielArea += jac;
    }
    for(unsigned i = 0; i < nDofs; i++) {
      sol->_Sol[weightIdx]->add(idof[i], ielArea * WEIGHT[ielt][i]);
    }
  }
  sol->_Sol[weightIdx]->close();
  sol->_Sol[kernelIdx]->close();
  for(unsigned k = 0; k < dim; k++) {
    for(unsigned j = 0; j < dim; j++) {
      sol->_Sol[gradDIdx[k][j]]->close();
    }
  }

  for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    double kernel = (*sol->_Sol[kernelIdx])(i);
    if(kernel != 0) {
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned j = 0; j < dim; j++) {
          double value = (*sol->_Sol[gradDIdx[k][j]])(i);
          sol->_Sol[gradDIdx[k][j]]->set(i, value / kernel);
        }
      }
    }
  }

  for(unsigned k = 0; k < dim; k++) {
    for(unsigned j = 0; j < dim; j++) {
      sol->_Sol[gradDIdx[k][j]]->close();
    }
  }

  //END bulk markers
  const char Nname[3][4] = {"NX", "NY", "NZ"};

  std::vector < unsigned > NIdx(dim);
  for(unsigned k = 0; k < dim; k++) {
    NIdx[k] = mlSol->GetIndex(&Nname[k][0]);
  }

  unsigned areaIdx = mlSol->GetIndex("area");
  sol->_Sol[areaIdx]->zero();
  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[NIdx[k]]->zero();
  }

  const MyVector< short unsigned> elMat = msh->el->GetElementMaterial();
  std::vector < short unsigned > elMatLoc;
  elMat.localize(elMatLoc);

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned ielMat = msh->GetElementMaterial(iel);
    if(ielMat == 4) {
      unsigned nDofs = msh->GetElementDofNumber(iel, solType);
      for(unsigned  k = 0; k < dim; k++) {
        vx[k].resize(nDofs);
        solD[k].resize(nDofs);
      }
      d2max.resize(nDofs);
      idof.resize(nDofs);
      for(unsigned i = 0; i < nDofs; i++) {
        idof[i] = msh->GetSolutionDof(i, iel, solType);
        d2max[i] = (*sol->_Sol[d2maxIdx])(idof[i]);
        for(unsigned  k = 0; k < dim; k++) {
          solD[k][i] = (*sol->_Sol[DIdx[k]])(idof[i]);
          vx[k][i] = (*msh->_topology->_Sol[k])(idof[i]) + solD[k][i];
        }
      }
      for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
        int jel = msh->el->GetFaceElementIndex(iel, iface) - 1;
        if(jel >= 0) { // iface is not a boundary of the domain

          unsigned jelMat = elMatLoc[jel];
          if(ielMat != jelMat) { //iel and jel are on the FSI interface

            const unsigned faceGeom = msh->GetElementFaceType(iel, iface);
            unsigned faceDofs = msh->GetElementFaceDofNumber(iel, iface, solType);
            std::vector  < std::vector  <  double > > fvx(dim);    // A matrix holding the face coordinates rowwise.
            std::vector  < std::vector  <  double > > fsolD(dim);
            std::vector  <  double > fd2max;
            for(int k = 0; k < dim; k++) {
              fvx[k].resize(faceDofs);
              fsolD[k].resize(faceDofs);
            }
            fd2max.resize(faceDofs);

            std::vector<unsigned> fdof(faceDofs);
            for(unsigned i = 0; i < faceDofs; i++) {
              unsigned inode = msh->GetLocalFaceVertexIndex(iel, iface, i);    // face-to-element local node mapping.
              fdof[i] = idof[inode];
              fd2max[i] = d2max[inode];
              for(unsigned k = 0; k < dim; k++) {
                fvx[k][i] =  vx[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
                fsolD[k][i] = solD[k][inode];
              }
            }

            double faceArea = 0.;
            for(unsigned ig = 0; ig < msh->_finiteElement[faceGeom][solType]->GetGaussPointNumber(); ig++) {
              double jac;
              std::vector < double> normal;
              msh->_finiteElement[faceGeom][solType]->JacobianSur(fvx, ig, jac, phi, gradPhi, normal);
              faceArea += jac;

              std::vector< double > xg(dim - 1, 0.);
              for(unsigned  k = 0; k < dim - 1; k++) { //solution
                for(unsigned i = 0; i < faceDofs; i++) { //node
                  xg[k] += (fvx[k][i] - fsolD[k][i]) * phi[i];
                }
              }

              for(unsigned i = 0; i < faceDofs; i++) { //node
                double d2 = 0.;
                for(unsigned  k = 0; k < dim - 1; k++) { //solution
                  d2 += ((fvx[k][i] - fsolD[k][i]) - xg[k]) * ((fvx[k][i] - fsolD[k][i]) - xg[k]);
                }

                for(unsigned  k = 0; k < dim; k++) {
                  sol->_Sol[NIdx[k]]->add(fdof[i], normal[k] * jac * (fd2max[i] - d2) * (fd2max[i] - d2));
                }
              }
            }
            for(unsigned i = 0; i < faceDofs; i++) {
              sol->_Sol[areaIdx]->add(fdof[i], faceArea * WEIGHT[faceGeom][i]);
            }
          }
        }
      }
    }
  }
  sol->_Sol[areaIdx]->close();
  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[NIdx[k]]->close();
  }

  for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    double area = (*sol->_Sol[areaIdx])(i);
    if(area != 0) {
      double norm = 0;
      for(unsigned k = 0; k < dim; k++) {
        double value = (*sol->_Sol[NIdx[k]])(i);
        norm += value * value;
      }
      norm = sqrt(norm);
      for(unsigned k = 0; k < dim; k++) {
        double value = (*sol->_Sol[NIdx[k]])(i) * area / norm;
        sol->_Sol[NIdx[k]]->set(i, value);
      }
    }
  }

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[NIdx[k]]->close();
  }
}

void FlagElements(MultiLevelMesh & mlMesh, const unsigned & layers) {

  unsigned level = mlMesh.GetNumberOfLevels() - 1;
  Mesh *msh   = mlMesh.GetLevel(level);
  unsigned iproc  = msh->processor_id();
  unsigned nprocs  = msh->n_processors();
  const unsigned dim = msh->GetDimension();

  const MyVector< short unsigned> elMat = msh->el->GetElementMaterial();
  std::vector < short unsigned > elMatLoc;
  elMat.localize(elMatLoc);

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    if(msh->el->GetIfElementCanBeRefined(iel)) {
      bool refine = 0;
      unsigned ielMat = msh->GetElementMaterial(iel);

      for(unsigned j = 0; j < msh->el->GetElementNearElementSize(iel, 1); j++) { //loop all over the neighboring elements
        int jel = msh->el->GetElementNearElement(iel, j);

        if(jel >= 0) { // jel is not outside the domain
          unsigned jelMat = elMatLoc[jel];
          if(ielMat != jelMat) { //iel and jel are on the FSI interface
            refine = true;
            break;
          }
        }
      }
      if(refine)  {
        msh->_topology->_Sol[msh->GetAmrIndex()]->set(iel, 1.);
      }
    }
  }
  msh->_topology->_Sol[msh->GetAmrIndex()]->close();

  for(unsigned k = 1; k < layers; k++) {

    std::vector<double> amrIdxLocal;
    msh->_topology->_Sol[msh->GetAmrIndex()]->localize_to_all(amrIdxLocal);

    for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
      if(msh->el->GetIfElementCanBeRefined(iel) && (*msh->_topology->_Sol[msh->GetAmrIndex()])(iel) == 0.) {
        bool refine = 0;

        for(unsigned j = 0; j < msh->el->GetElementNearElementSize(iel, 1); j++) {
          int jel = msh->el->GetElementNearElement(iel, j);
          if(jel >= 0) {
            //unsigned jproc = msh->IsdomBisectionSearch(jel, 3);
            if(amrIdxLocal[jel] == 1) {
              //if(iproc == jproc && (*msh->_topology->_Sol[msh->GetAmrIndex()])(jel) == 1.) { // iface is not a boundary of the domain
              refine = true;
              break;
            }
          }
        }
        if(refine)  {
          msh->_topology->_Sol[msh->GetAmrIndex()]->set(iel, 0.5);
        }
      }
    }
    for(unsigned iel = 0; iel < msh->el->GetElementNumber(); iel++) {
      if((*msh->_topology->_Sol[msh->GetAmrIndex()])(iel) == 0.5) {
        msh->_topology->_Sol[msh->GetAmrIndex()]->set(iel, 1.);
      }
    }
    msh->_topology->_Sol[msh->GetAmrIndex()]->close();
  }

}


void AssembleMarkerStructure(MultiLevelProblem& ml_prob) {

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  //pointers and references

  NonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<NonLinearImplicitSystem> ("Marker");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel(level);     // pointer to the solution (level) object

  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];  // pointer to the equation (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)
  SparseMatrix* myKK = myLinEqSolver->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* myRES =  myLinEqSolver->_RES;  // pointer to the global residual vector object in pdeSys (level)

  myKK->zero();
  myRES->zero();


  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = msh->GetDimension();

  // data
  unsigned iproc  = msh->processor_id();
  unsigned nprocs  = msh->n_processors();

  vector< vector< adept::adouble > > solD(dim);   // background displacement at n+1


  vector< double > rhs;    // local redidual vector
  vector< vector< adept::adouble > > aResD(dim);     // local redidual vector
  vector < double > Jac;

  std::vector <unsigned> sysDofsAll;

  vector < double > phi;

  vector < double > gradPhiHat;
  vector < adept::adouble> gradPhi;     // phi_x

  vector <vector < double> > vxHat(dim);
  vector <vector < adept::adouble> > vx(dim);

  double weightHat;
  adept::adouble weight;

  std::cout.precision(10);

  //variable-name handling
  const char varname[10][5] = {"DX", "DY", "DZ"};

  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexPdeD(dim);
  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol->GetIndex(&varname[ivar][0]);
    indexPdeD[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);       // DX, DY, DZ
  }
  unsigned solType = mlSol->GetSolutionType(&varname[0][0]);
  unsigned meshType = 2;

  start_time = clock();

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned ielMat = msh->GetElementMaterial(iel);

    short unsigned ielt = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofsAll = dim * nDofs;

    // resize local arrays
    sysDofsAll.resize(nDofsAll);

    for(unsigned  k = 0; k < dim; k++) {
      solD[k].resize(nDofs);
      aResD[k].assign(nDofs, 0.);
      vxHat[k].resize(nDofs);
      vx[k].resize(nDofs);
    }

    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solType);
      for(unsigned  k = 0; k < dim; k++) {
        solD[k][i] = (*mysolution->_Sol[indexSolD[k]])(idof); //t_{n+1} -t_n
        sysDofsAll[i + k * nDofs] = myLinEqSolver->GetSystemDof(indexSolD[k], indexPdeD[k], i, iel);
      }
    }

    if(ielMat == 2) {
      // start a new recording of all the operations involving adept::adouble variables
      s.new_recording();

      for(unsigned i = 0; i < nDofs; i++) {
        unsigned idofX = msh->GetSolutionDof(i, iel, meshType);
        for(unsigned  k = 0; k < dim; k++) {
          vxHat[k][i] = (*msh->_topology->_Sol[k])(idofX); // undeformed background configuration at t_{n}
          vx[k][i]  = vxHat[k][i] + solD[k][i]; // deformed background configuration at alpha_f/theta
        }
      }

      // *** Gauss point loop ***
      for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

        msh->_finiteElement[ielt][solType]->Jacobian(vxHat, ig, weightHat, phi, gradPhiHat);
        msh->_finiteElement[ielt][solType]->Jacobian(vx,    ig, weight,    phi, gradPhi);

        vector < vector < adept::adouble > > gradSolDgHat(dim, vector < adept::adouble >(dim, 0.));
        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned j = 0; j < dim; j++) {
            for(unsigned  k = 0; k < dim; k++) {
              gradSolDgHat[k][j] += gradPhiHat[i * dim + j] * solD[k][i]; //gradient of new solution with respect to deformed reference configuration
            }
          }
        }

        adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned k = 0; k < dim; k++) {
            F[j][k] += gradSolDgHat[j][k];
          }
        }

        adept::adouble J_Hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                                - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

        adept::adouble B[3][3];
        for(unsigned i = 0; i < 3; i++) {
          for(int j = 0; j < 3; j++) {
            B[i][j] = 0.;
            for(unsigned k = 0; k < 3; k++) {
              //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
              B[i][j] += F[i][k] * F[j][k];
            }
          }
        }

        adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];
        adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
        adept::adouble sigma[3][3];

        double E = 1.;
        double nu = 0.4;

        double mu = E / (2. * (1. + nu));
        double lambda = (E * nu) / ((1. + nu) * (1. - 2.*nu));

        for(unsigned j = 0; j < 3; j++) {
          for(unsigned k = 0; k < 3; k++) {
            sigma[j][k] = lambda * log(J_Hat) / J_Hat * Id2th[j][k] + mu / J_Hat * (B[j][k] - Id2th[j][k]);    // alternative formulation
          }
        }
        //END computation of the Cauchy Stress
        for(unsigned i = 0; i < nDofs; i++) {//Auxiliary Equations
          for(unsigned k = 0.; k < dim; k++) {
            adept::adouble cauchy = 0.;
            for(unsigned j = 0.; j < dim; j++) {
              cauchy += sigma[k][j] * gradPhi[i * dim + j] ;
            }
            aResD[k][i] += cauchy * weight;
          }
        }
      } // end gauss point loop
    }
    //copy the value of the adept::adoube aRes in double Res and store them in RES
    rhs.resize(nDofsAll);   //resize

    for(int i = 0; i < nDofs; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        rhs[ i +  k * nDofs ] = -aResD[k][i].value();
      }
    }
    myRES->add_vector_blocked(rhs, sysDofsAll);
    Jac.resize(nDofsAll * nDofsAll);
    // define the dependent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&aResD[k][0], nDofs);
    }
    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solD[k][0], nDofs);
    }

    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0], true);
    myKK->add_matrix_blocked(Jac, sysDofsAll, sysDofsAll);

    s.clear_independents();
    s.clear_dependents();

  }

  myRES->close();
  myKK->close();

//   PetscViewer    viewer1;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1800, 1800, &viewer1);
//   PetscObjectSetName((PetscObject) viewer1, "FSI matrix");
//   PetscViewerPushFormat(viewer1, PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast< PetscMatrix* >(myKK))->mat(), viewer1);
// 
//   double a;
//   std::cin >> a;

// *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
// ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}




