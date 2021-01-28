#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "MultiLevelSolution.hpp"

#include "MeshRefinement.hpp"

const double WEIGHT[6][27] = {
  {}, //hex
  {}, //tet
  {}, //wedge
  {0.0625, 0.0625, 0.0625, 0.0625, 0.125, 0.125, 0.125, 0.125, 0.25 }, 
  {1./12., 1./12., 1./12., 0.25, 0.25, 0.25, 0.}, 
  {0.25, 0.25, .5}  
} ;


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int facename, const double t) {
  bool test = 1; //dirichlet
  value = 0.;
  return test;
}

void BuildMarkers(MultiLevelMesh& mlMesh);

void FlagElements(MultiLevelMesh& mlMesh, const unsigned &layer);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  mlMsh.ReadCoarseMesh("../input/beam.neu", "fifth", scalingFactor);

  unsigned numberOfRefinement = 1;

  for(unsigned i = 0; i < numberOfRefinement; i++) {
    FlagElements(mlMsh, 2);
    mlMsh.AddAMRMeshLevel();
  }

  mlMsh.EraseCoarseLevels(numberOfRefinement);

  unsigned dim = mlMsh.GetDimension();

  FEOrder femOrder = FIRST;

  MultiLevelSolution mlSol(&mlMsh);

  //add variables to mlSol
  mlSol.AddSolution("u", LAGRANGE, femOrder, 2);

  mlSol.Initialize("All");
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("u", "Steady");

  BuildMarkers(mlMsh);

  //******* Print solution *******
  mlSol.SetWriter(VTK);
  std::vector<std::string> mov_vars;
  std::vector<std::string> print_vars;
  print_vars.push_back("All");
  mlSol.GetWriter()->SetDebugOutput(false);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);

  return 0;

} //end main




void BuildMarkers(MultiLevelMesh& mlMesh) {

  unsigned level = mlMesh.GetNumberOfLevels() - 1;
  Mesh *msh   = mlMesh.GetLevel(level);
  unsigned iproc  = msh->processor_id();
  const unsigned dim = msh->GetDimension();

  unsigned solType = 2;
  unsigned xpSize = msh->el->GetNodeNumber();
  std::vector < std::vector < double > > xp(xpSize);

  for(unsigned i = 0; i < xpSize; i++) {
    xp[i].resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      xp[i][k] = (*msh->_topology->_Sol[k])(i);
    }
  }


  std::vector < double > weight(xpSize, 0.);
  std::vector < double > phi;
  std::vector < double > gradPhi;


  std::vector <std::vector < double> > vx(dim);
  std::vector<unsigned> dof;

  double area = 0.;

  for(unsigned iel = 0; iel < msh->el->GetElementNumber(); iel++) {

    short unsigned ielt = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);
    dof.resize(nDofs);
    for(unsigned  k = 0; k < dim; k++) {
      vx[k].resize(nDofs);
    }

    for(unsigned i = 0; i < nDofs; i++) {
      dof[i] = msh->GetSolutionDof(i, iel, solType);
    }

    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idofX = msh->GetSolutionDof(i, iel, solType);
      for(unsigned  k = 0; k < dim; k++) {
        vx[k][i] = (*msh->_topology->_Sol[k])(idofX);
      }
    }

    double ielArea = 0.;
    for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {
      double jac;
      msh->_finiteElement[ielt][solType]->Jacobian(vx, ig, jac, phi, gradPhi);
      ielArea += jac;
    }


    for(unsigned i = 0; i < nDofs; i++) {
      weight[msh->GetSolutionDof(i, iel, solType)] += ielArea * WEIGHT[ielt][i];
      area += ielArea * WEIGHT[ielt][i];
    }
    

  }

  std::cout << area << " " << 2 * 4.99 + 0.5 * M_PI * 1. * 1.;

  std::vector < std::vector < std::vector < double > > >  bulkPoints(1);
  bulkPoints[0] = xp;
  PrintLine(DEFAULT_OUTPUTDIR, "bulk", bulkPoints, 0);

}

void FlagElements(MultiLevelMesh& mlMesh, const unsigned &layers) {

  unsigned level = mlMesh.GetNumberOfLevels() - 1;
  Mesh *msh   = mlMesh.GetLevel(level);
  unsigned iproc  = msh->processor_id();
  const unsigned dim = msh->GetDimension();

  for(unsigned iel = 0; iel < msh->el->GetElementNumber(); iel++) {
    if(msh->el->GetIfElementCanBeRefined(iel)) {
      bool refine = 0;
      unsigned ielGrup = msh->GetElementGroup(iel);

      for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
        int jel = msh->el->GetFaceElementIndex(iel, iface) - 1;
        if(jel >= 0) { // iface is not a boundary of the domain
          unsigned jelGrup = msh->GetElementGroup(jel);
          if(ielGrup != jelGrup) { //iel and jel are on the FSI interface
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
    for(unsigned iel = 0; iel < msh->el->GetElementNumber(); iel++) {
      if(msh->el->GetIfElementCanBeRefined(iel) && (*msh->_topology->_Sol[msh->GetAmrIndex()])(iel) == 0.) {
        bool refine = 0;
        for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          int jel = msh->el->GetFaceElementIndex(iel, iface) - 1;
          if(jel >= 0 && (*msh->_topology->_Sol[msh->GetAmrIndex()])(jel) == 1.) { // iface is not a boundary of the domain
            refine = true;
            break;
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
