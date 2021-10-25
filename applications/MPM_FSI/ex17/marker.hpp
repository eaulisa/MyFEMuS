
#include "Marker.hpp"
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
  bool test = 1; //dirichlet
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
  mlSol.AddSolution("DX", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("DY", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("DZ", LAGRANGE, femOrder, 0, false);

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

  mlSol.AddSolution("mtype", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("NX", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("NY", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("NZ", LAGRANGE, femOrder, 0, false);



  mlSol.AddSolution("weight", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("area", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("iel", LAGRANGE, femOrder, 0, false); // for each node it stores the maximum distance2 in the support
  mlSol.AddSolution("dist", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("kernel", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("d2max", LAGRANGE, femOrder, 0, false); // for each node it stores the maximum distance2 in the support


  mlSol.Initialize("All");
  //mlSol.Initialize("DX", InitVariableDX);
  //mlSol.Initialize("VY", InitVariableVY);

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

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
      if(ielMat == 4) sol->_Sol[mtypeIdx]->set(idof[i], 2);
    }
  }
  sol->_Sol[d2maxIdx]->closeWithMaxValue();
  sol->_Sol[mtypeIdx]->close();

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned ielMat = msh->GetElementMaterial(iel);
    if(ielMat == 2) {

      short unsigned ielt = msh->GetElementType(iel);
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
      vx[k].resize(nDofs);
      solD[k].resize(nDofs);
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
      }
    }

    double ielArea = 0.;
    for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {
      double jac;
      msh->_finiteElement[ielt][solType]->Jacobian(vx, ig, jac, phi, gradPhi);
      ielArea += jac;


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

        if(ielMat == 4){// || mtype[i] == 0) {
          double d2 = 0.;
          for(unsigned  k = 0; k < dim; k++) { //solution
            d2 += (vx[k][i] - xg[k]) * (vx[k][i] - xg[k]);
          }

          sol->_Sol[kernelIdx]->add(idof[i], jac * (d2max[i] - d2) * (d2max[i] - d2));
          for(unsigned  k = 0; k < dim; k++) {
            for(unsigned j = 0; j < dim; j++) {
              sol->_Sol[gradDIdx[k][j]]->add(idof[i], gradSolDg[k][j] * jac * (d2max[i] - d2) * (d2max[i] - d2));
            }
          }
        }
      }


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




