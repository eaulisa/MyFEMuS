

void InitializeBackgroundVariables(MultiLevelSolution &mlSol) {

  unsigned dim = mlSol._mlMesh->GetDimension();

  FEOrder femOrder = FIRST;

  //add variables to mlSol
  mlSol.AddSolution("DX", LAGRANGE, femOrder, 0, false);
  mlSol.AddSolution("DY", LAGRANGE, femOrder, 0, false);
  if(dim == 3) mlSol.AddSolution("DZ", LAGRANGE, femOrder, 0, false);

  mlSol.Initialize("All");
  //mlSol.Initialize("DX", InitVariableDX);
  //mlSol.Initialize("DY", InitVariableDY);

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

}


