

void InitializeBackgroundVariables(MultiLevelSolution &mlSol) {

  unsigned dim = mlSol._mlMesh->GetDimension();

  FEOrder femOrder = FIRST;

  //add variables to mlSol
  mlSol.AddSolution("DX", LAGRANGE, femOrder, 2);
  mlSol.AddSolution("DY", LAGRANGE, femOrder, 2);
  if(dim == 3) mlSol.AddSolution("DZ", LAGRANGE, femOrder, 2);
  mlSol.AddSolution("VX", LAGRANGE, femOrder, 2);
  mlSol.AddSolution("VY", LAGRANGE, femOrder, 2);
  if(dim == 3) mlSol.AddSolution("VZ", LAGRANGE, femOrder, 2);
  
   mlSol.AddSolution("P", LAGRANGE, FIRST, 2);

  mlSol.AddSolution("eflag", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("nflag", LAGRANGE, femOrder, 0, false);
  
  mlSol.Initialize("All");
  mlSol.Initialize("DX", InitVariableDX);
  mlSol.Initialize("DY", InitVariableDY);

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

}


