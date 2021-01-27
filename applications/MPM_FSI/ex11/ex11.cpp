#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "MultiLevelSolution.hpp"


const double WEIGHT[6][27] = {
    {}, //hex
    {}, //tet
    {}, //wedge
    {0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 1. }, //quad sum is 4
    {1./24., 1./24., 1./24., 1./8., 1./8., 1./8., 0.}, //tri sun 1/2
    {0.5, 0.5, 1.}  //line sum 2
} ;


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int facename, const double t) {
  bool test = 1; //dirichlet
  value = 0.;
  return test;
}

void BuildMarkers(MultiLevelMesh& mlMesh);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfUniformLevels = 2; //for refinement in 3D
  //unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;


  mlMsh.ReadCoarseMesh("../input/beam.neu", "fifth", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);

  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);
  numberOfUniformLevels = 1;

  unsigned dim = mlMsh.GetDimension();

  FEOrder femOrder = FIRST;

  MultiLevelSolution mlSol(&mlMsh);
  // add variables to mlSol
  mlSol.AddSolution("u", LAGRANGE, femOrder, 2);

  mlSol.Initialize("All");
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("u", "Steady");

  BuildMarkers(mlMsh);

  // ******* Print solution *******
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


    for(unsigned i = 0; i < nDofs; i++) {
      std::vector < double > xi(dim);
      
      for(unsigned k = 0; k < dim; k++) {
        xi[k] = *(msh->_finiteElement[ielt][solType]->GetBasis()->GetXcoarse(i) + k);
      }
      
      
      double jac;
      msh->_finiteElement[ielt][solType]->Jacobian(vx, xi, jac, phi, gradPhi);
      weight[msh->GetSolutionDof(i, iel, solType)] += jac * WEIGHT[ielt][i];
      area += jac * WEIGHT[ielt][i];
    }

   

  }

  std::cout << area << " " << 2 * 4.99 + 0.5 * M_PI * 1. * 1.; 
  
  std::vector < std::vector < std::vector < double > > >  bulkPoints(1);
  bulkPoints[0] = xp;
  PrintLine(DEFAULT_OUTPUTDIR, "bulk", bulkPoints, 0);

}
