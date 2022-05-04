void AssembleGhostPenaltyP(MultiLevelProblem& ml_prob, const bool &omega1) {

  clock_t start_time;

  //pointers and references

  LinearImplicitSystem& mlPdeSys  = ml_prob.get_system<LinearImplicitSystem> ("Nitsche");   // pointer to the linear implicit system named "Poisson"
  const unsigned  level = mlPdeSys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel(level);     // pointer to the solution (level) object

  //Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* myLinEqSolver = mlPdeSys._LinSolver[level];  // pointer to the equation (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)
  SparseMatrix* myKK = myLinEqSolver->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* myRES =  myLinEqSolver->_RES;  // pointer to the global residual vector object in pdeSys (level)

  MatSetOption((static_cast< PetscMatrix* >(myKK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);


  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = msh->GetDimension();
  const unsigned dim2 = 3 * (dim - 1);//Dimension of the second order derivative
  // this is independent. In 2D-> xx,yy,xy , In 3D-> xx,yy,zz,xy,xz,yz
  // data
  unsigned iproc  = msh->processor_id();
  unsigned nprocs  = msh->n_processors();//Number of processors that are running 

  //quantities for iel will have index1
  //quantities for jel will have index2
  //Only the bool omega1 refers to subdomain 1 and 2 in the Nietsche coupling

  //solution of one element and neighboring element
  vector< adept::adouble > solu1; 
  vector< adept::adouble > solu2; 

  vector< adept::adouble > aRes1;     // local redidual vector
  vector< adept::adouble > aRes2;     // local redidual vector

  vector< double > rhs; // local redidual vector
  vector < double > Jac;

  std::vector <unsigned> sysDofs;


  double weight;
  std::vector < double > phi;
  std::vector < double> gradPhi;

  
  double weight1;
  std::vector < double > phi1;
  std::vector < double> gradPhi1;
  std::vector < double> nablaPhi1;

  double weight2;
  std::vector < double > phi2;
  std::vector < double> gradPhi2;
  std::vector < double> nablaPhi2;

  //****coordinates of two neighboring elements?
  vector <vector < double> > vx1(dim);
  vector <vector < double> > vx2(dim);

  //reading parameters for fluid FEM domain
  double alpha = (omega1) ? alpha1 : alpha2;


  std::cout.precision(10);

  //variable-name handling
  const char varname[2][3] = {"u1", "u2"};//vector of strings. two entries and atmost 3 characters

  //if omega1 is true, it's going to give you u1. if it is false, it is u2.
  unsigned indexSol = mlSol->GetIndex(&varname[!omega1][0]);
  unsigned indexPde = mlPdeSys.GetSolPdeIndex(&varname[!omega1][0]);
  unsigned solType = mlSol->GetSolutionType(&varname[0][0]);

  unsigned eflagIndex = mlSol->GetIndex("eflag");//whether you have a cut or not. if eflag=1 we have a cut cell.

  start_time = clock();

  std::vector < std::vector < std::vector <double > > > aP1(3);
  std::vector < std::vector < std::vector <double > > > aP2(3);


  //flagmark
  //loop on the element on the mesh
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    //if we have a cut element go in
    unsigned eFlag1 = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(iel) + 0.25));
    if(eFlag1 == 1) {
      //in cut element we need to know what type of cut element we have and number of nodes we have.
      short unsigned ielt1 = msh->GetElementType(iel);
      unsigned nDofs1 = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
      // resize local arrays
      sysDofs.resize(nDofs1);
      solu1.resize(nDofs1);
      for(unsigned  k = 0; k < dim; k++) {
        vx1[k].resize(nDofs1);
      }

      for(unsigned i = 0; i < nDofs1; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solType); //local to global mapping for solution solType
        solu1[i] = (*mysolution->_Sol[indexSol])(idof);
        sysDofs[i] = myLinEqSolver->GetSystemDof(indexSol, indexPde, i, iel); //local to global mapping for the PDE
      }
      //coordinates of iel element
      for(unsigned i = 0; i < nDofs1; i++) {
        unsigned idofX = msh->GetSolutionDof(i, iel, 2);
        for(unsigned  k = 0; k < dim; k++) {
          vx1[k][i] = (*msh->_topology->_Sol[k])(idofX);
        }
      }

      bool aP1IsInitialized = false;
      //looping on faces of element iel
      for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
        int jel = el->GetFaceElementIndex(iel, iface) - 1;//jel is the neighboring element
        if(jel >= 0) { // iface is not a boundary of the domain.
          //if we are on the boundary jel is negative number.we ignore that case
            
          unsigned jproc = msh->IsdomBisectionSearch(jel, 3);
          //jel is piece wise discontinuous FEM it is 3. 0=Lagrange Linear, 1=Lagrange quadratic, 2=Lagrage biquadratic, 3= constat piecewise discontinuous linear. 
          if(jproc == iproc) {//iel and jel belongs to the same processor
            //since we are ine the same processor, we can exract flag jel
            unsigned eFlag2 = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(jel) + 0.25));
            
             //if you are in omega1, do ghost cell in 0,1 , if you are in omega2 do 1,2 with jel>iel (as we want to do integral one type)
            if(eFlag2 == !omega1 * 2 || (eFlag2 == 1 && jel > iel)) {

              short unsigned ielt2 = msh->GetElementType(jel);
              unsigned nDofs2 = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs
              // resize local arrays
              sysDofs.resize(nDofs1 + nDofs2);
              solu2.resize(nDofs2);

              for(unsigned  k = 0; k < dim; k++) {
                vx2[k].resize(nDofs2);
              }
              for(unsigned i = 0; i < nDofs2; i++) {
                unsigned idof = msh->GetSolutionDof(i, jel, solType);
                solu2[i] = (*mysolution->_Sol[indexSol])(idof);
                sysDofs[nDofs1 + i] = myLinEqSolver->GetSystemDof(indexSol, indexPde, i, jel);
              }
              for(unsigned i = 0; i < nDofs2; i++) {
                unsigned idofX = msh->GetSolutionDof(i, jel, 2);
                for(unsigned  k = 0; k < dim; k++) {
                  vx2[k][i] = (*msh->_topology->_Sol[k])(idofX);
                }
              }
              //now we are on iel and on a particular face and corresponding jel. we have all the information.
              aRes1.assign(nDofs1, 0.);
              aRes2.assign(nDofs2, 0.);

              s.new_recording();
               
              //gathering information on the boundary
              const unsigned faceGeom = msh->GetElementFaceType(iel, iface);
              unsigned faceDofs = msh->GetElementFaceDofNumber(iel, iface, solType);
              std::vector  < std::vector  <  double > > faceVx(dim);    // A matrix holding the face coordinates rowwise.
              for(int k = 0; k < dim; k++) {
                faceVx[k].resize(faceDofs);
              }
              for(unsigned i = 0; i < faceDofs; i++) {
                unsigned inode = msh->GetLocalFaceVertexIndex(iel, iface, i);    // face-to-element local node mapping.
                for(unsigned k = 0; k < dim; k++) {
                  faceVx[k][i] =  vx1[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
                }
              }

              double h = sqrt(dim - 1) * sqrt((faceVx[0][0] - faceVx[0][1]) * (faceVx[0][0] - faceVx[0][1]) +
                                              (faceVx[1][0] - faceVx[1][1]) * (faceVx[1][0] - faceVx[1][1])) ;
              double h2 = h * h;
              double h3 = h * h * h;


              if(!aP1IsInitialized) { //build the basis 1,x,y,z... corresponding to the solution type
                for(unsigned jtype = 0; jtype < solType + 1; jtype++) {
                  ProjectNodalToPolynomialCoefficients(aP1[jtype], vx1, ielt1, jtype);
                }
                aP1IsInitialized = true;
              }

              for(unsigned jtype = 0; jtype < solType + 1; jtype++) {
                ProjectNodalToPolynomialCoefficients(aP2[jtype], vx2, ielt2, jtype);
              }
              //**gauss loop on the face
              for(unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solType]->GetGaussPointNumber(); ig++) {

                std::vector < double> normal;//**normal is going from iel to jel
                msh->_finiteElement[faceGeom][solType]->JacobianSur(faceVx, ig, weight, phi, gradPhi, normal);//**here gradient for the test function is completely fake. important phi on the edge. we use phi to get physical gauss coordinates

                std::vector< double > xg(dim, 0.); // physical coordinates of the face Gauss point
                for(unsigned i = 0; i < faceDofs; i++) {
                  for(unsigned k = 0; k < dim; k++) {
                    xg[k] += phi[i] * faceVx[k][i];
                  }
                }
                //**now we can do the inverse mapping
                std::vector <double> xi1;//local coordinates of the face gauss point with respect to iel
                GetClosestPointInReferenceElement(vx1, xg, ielt1, xi1);

                bool inverseMapping = GetInverseMapping(solType, ielt1, aP1, xg, xi1, 100);
                if(!inverseMapping) {
                  std::cout << "InverseMapping1 failed at " << iel << " " << jel << " " << iface << std::endl;
                }

                std::vector <double> xi2;//local coordinates of the face gauss point with respect to jel
                GetClosestPointInReferenceElement(vx2, xg, ielt2, xi2);

                inverseMapping = GetInverseMapping(solType, ielt2, aP2, xg, xi2, 100);
                if(!inverseMapping) {
                  std::cout << "InverseMapping2 failed at " << iel << " " << jel << " " << iface << std::endl;
                }

                //** having xi1 and xi2, now we can get gradient of laplace and the test function
                msh->_finiteElement[ielt1][solType]->Jacobian(vx1, xi1, weight1, phi1, gradPhi1, nablaPhi1);
                msh->_finiteElement[ielt2][solType]->Jacobian(vx2, xi2, weight2, phi2, gradPhi2, nablaPhi2);

                //**now you have surface are and normal. now we can assemble penalty term.

                adept::adouble gradSolu1DotN = 0.;//**gradient of the solution one in normal direction
                adept::adouble hessSolu1DotN = 0.;//** Hessian of the solution in the normal direction
                for(unsigned i = 0; i < nDofs1; i++) {
                  for(unsigned J = 0; J < dim; J++) {

                    //**dot product between the gradient and the normal
                    //**gradient of phi in the direction j multiply by normal j
                    gradSolu1DotN += solu1[i] * gradPhi1[i * dim + J] * normal[J];
                    for(unsigned K = 0; K < dim; K++) {
                      //2D xx, yy, xy
                      //3D xx, yy, zz, xy, yz ,zx
                      unsigned L;
                      if(J == K) L = J;
                      else if(1 == J + K) L = dim;     // xy
                      else if(2 == J + K) L = dim + 2; // xz
                      else if(3 == J + K) L = dim + 1; // yz
                      hessSolu1DotN += solu1[i] * normal[J] * nablaPhi1[i * dim2 + L] * normal[K];
                    }
                  }
                }
//** same thing should be done for element 2
                adept::adouble gradSolu2DotN = 0.;
                adept::adouble hessSolu2DotN = 0.;
                for(unsigned i = 0; i < nDofs2; i++) {
                  for(unsigned J = 0; J < dim; J++) {
                    gradSolu2DotN += solu2[i] * gradPhi2[i * dim + J] * normal[J];
                    for(unsigned K = 0; K < dim; K++) {
                      //2D xx, yy, xy
                      //3D xx, yy, zz, xy, yz ,zx
                      unsigned L;
                      if(J == K) L = J;
                      else if(1 == J + K) L = dim;     // xy
                      else if(2 == J + K) L = dim + 2; // xz
                      else if(3 == J + K) L = dim + 1; // yz
                      hessSolu2DotN += solu2[i] * normal[J] * nablaPhi2[i * dim2 + L] * normal[K];
                    }
                  }
                }


                double C1 = 0.05 * alpha;  // [alpha];

                for(unsigned i = 0; i < nDofs1; i++) {

                  for(unsigned J = 0; J < dim; J++) {
                    aRes1[i] +=  C1 * h * gradPhi1[i * dim + J] * normal[J] * (gradSolu1DotN - gradSolu2DotN) * weight;



                    for(unsigned K = 0; K < dim; K++) {

                      unsigned L;
                      if(J == K) L = J;
                      else if(1 == J + K) L = dim;     // xy
                      else if(2 == J + K) L = dim + 2; // xz
                      else if(3 == J + K) L = dim + 1; // yz

                      aRes1[i] += C1 * h3 * normal[J] * nablaPhi1[i * dim2 + L] * normal[K] * (hessSolu1DotN - hessSolu2DotN) * weight;
                    }
                  }
                }

                for(unsigned i = 0; i < nDofs2; i++) {

                  for(unsigned J = 0; J < dim; J++) {
                    aRes2[i] +=  -C1 * h * gradPhi2[i * dim + J] * normal[J] * (gradSolu1DotN - gradSolu2DotN) * weight;

                    for(unsigned K = 0; K < dim; K++) {

                      unsigned L;
                      if(J == K) L = J;
                      else if(1 == J + K) L = dim;     // xy
                      else if(2 == J + K) L = dim + 2; // xz
                      else if(3 == J + K) L = dim + 1; // yz

                      aRes2[i] +=  -C1 * h3 * normal[J] * nablaPhi2[i * dim2 + L] * normal[K] * (hessSolu1DotN - hessSolu2DotN) * weight;
                    }
                  }
                }


              }

              //copy the value of the adept::adoube aRes in double Res and store them in RES
              rhs.resize(nDofs1 + nDofs2);   //resize
              for(int i = 0; i < nDofs1; i++) {
                rhs[ i ] = -aRes1[i].value();
              }

              for(int i = 0; i < nDofs2; i++) {
                rhs[nDofs1 + i] = -aRes2[i].value();
              }
              myRES->add_vector_blocked(rhs, sysDofs);

              s.dependent(&aRes1[0], nDofs1);
              s.dependent(&aRes2[0], nDofs2);
              s.independent(&solu1[0], nDofs1);
              s.independent(&solu2[0], nDofs2);

              Jac.resize((nDofs1 + nDofs2) * (nDofs1 + nDofs2));
              s.jacobian(&Jac[0], true);
              myKK->add_matrix_blocked(Jac, sysDofs, sysDofs);

              s.clear_independents();
              s.clear_dependents();

            }
          }
        }
      }
    }
  }

  //flagmark
  if(nprocs > 1) {
    for(unsigned kproc = 0; kproc < nprocs; kproc++) { // not scalable
      for(int iel = msh->_elementOffset[kproc]; iel < msh->_elementOffset[kproc + 1]; iel++) {
        //**who owns eFlag1? iproc. one of kproc is iproc. 
        unsigned eFlag1;

        if(iproc == kproc) {
          eFlag1 = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(iel) + 0.25));
        }
        //**broadcasting to all other processors eFlag1. they are expecting to recieve from kproc.
        MPI_Bcast(&eFlag1, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);
        
        //**now every processor knows eFlag1
        if(eFlag1 == 1) {
          //**only kproc knows the iformation about the faces.
          unsigned nFaces;
          if(iproc == kproc) {
            nFaces = msh->GetElementFaceNumber(iel);
          }
          MPI_Bcast(&nFaces, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);
          
          //**looping on the faces 
          for(unsigned iface = 0; iface < nFaces; iface++) {
            
            int jel;
            if(iproc == kproc) {
              jel = el->GetFaceElementIndex(iel, iface) - 1;
            }
            MPI_Bcast(&jel, 1, MPI_INT, kproc, PETSC_COMM_WORLD);

            if(jel >= 0) { // iface is not a boundary of the domain
            //**then somebody should owns jel. then we check to which process jel belogs to 
              unsigned jproc = msh->IsdomBisectionSearch(jel, 3);  // return  jproc for piece-wise constant discontinuous type (3)
              
              //**if jproc is different from kproc and then if iproc == kproc || iproc == jproc 
              if(jproc != kproc && (iproc == kproc || iproc == jproc)) {
              //**after this only jproc and kproc is working. So hereafter if we have exchanges we have to exchange only with in jproc and kproc.
                unsigned eFlag2;
                //**information about eFlag2 has jproc
                if(iproc == jproc) {
                  eFlag2 = static_cast <unsigned>(floor((*mysolution->_Sol[eflagIndex])(jel) + 0.25));
                  MPI_Send(&eFlag2, 1, MPI_UNSIGNED, kproc, 0, PETSC_COMM_WORLD);//**MPI_Send done by jproc to kproc
                }
                else if(iproc == kproc) {
                  MPI_Recv(&eFlag2, 1, MPI_UNSIGNED, jproc, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                }

                if(eFlag2 == !omega1 * 2 || (eFlag2 == 1 && jel > iel)) {
                //**both jproc and kproc going forward
                  
                  //std::cout << "I am " << iel << " on " << kproc << " talking with " << jel << " on " << jproc << std::endl;

                  short unsigned ielt1;
                  short unsigned ielt2;
                  //now jproc send all the information to iproc. and iproc does the essembly.
                  if(iproc == kproc) {
                    ielt1 = msh->GetElementType(iel);
                    MPI_Recv(&ielt2, 1, MPI_UNSIGNED_SHORT, jproc, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                  }
                  else if(iproc == jproc) {
                    ielt2 = msh->GetElementType(jel);
                    MPI_Send(&ielt2, 1, MPI_UNSIGNED_SHORT, kproc, 0, PETSC_COMM_WORLD);
                  }
                  

                  unsigned nDofs1;
                  unsigned nDofs2 = el->GetNVE(ielt2, solType);
                  solu2.resize(nDofs2);

                  for(unsigned  k = 0; k < dim; k++) {
                    vx2[k].resize(nDofs2);
                  }

                  std::vector < MPI_Request > reqs(2 + dim);
                  if(iproc == kproc) {

                    nDofs1 = el->GetNVE(ielt1, solType);
                    sysDofs.resize(nDofs1 + nDofs2);
                    solu1.resize(nDofs1);

                    for(unsigned  k = 0; k < dim; k++) {

                      vx1[k].resize(nDofs1);
                    }

                    for(unsigned i = 0; i < nDofs1; i++) {
                      unsigned idof = msh->GetSolutionDof(i, iel, solType);
                      solu1[i] = (*mysolution->_Sol[indexSol])(idof);
                      sysDofs[i] = myLinEqSolver->GetSystemDof(indexSol, indexPde, i, iel);
                    }

                    for(unsigned i = 0; i < nDofs1; i++) {
                      unsigned idofX = msh->GetSolutionDof(i, iel, 2);
                      for(unsigned  k = 0; k < dim; k++) {
                        vx1[k][i] = (*msh->_topology->_Sol[k])(idofX);
                      }
                    }
                    MPI_Irecv(solu2.data(), solu2.size(), MPI_DOUBLE, jproc, 0, PETSC_COMM_WORLD, &reqs[0]);
                    MPI_Irecv(&sysDofs[nDofs1], nDofs2, MPI_UNSIGNED, jproc, 1, PETSC_COMM_WORLD,  &reqs[1]);
                    for(unsigned k = 0; k < dim; k++) {
                      MPI_Irecv(vx2[k].data(), vx2[k].size(), MPI_DOUBLE, jproc, 2 + k, PETSC_COMM_WORLD, &reqs[2 + k]);
                    }
                  }
                  else if(iproc == jproc) {
                    sysDofs.resize(nDofs2);
                    for(unsigned i = 0; i < nDofs2; i++) {
                      unsigned idof = msh->GetSolutionDof(i, jel, solType);
                      solu2[i] = (*mysolution->_Sol[indexSol])(idof);
                      sysDofs[i] = myLinEqSolver->GetSystemDof(indexSol, indexPde, i, jel);
                    }

                    for(unsigned i = 0; i < nDofs2; i++) {
                      unsigned idofX = msh->GetSolutionDof(i, jel, 2);
                      for(unsigned  k = 0; k < dim; k++) {
                        vx2[k][i] = (*msh->_topology->_Sol[k])(idofX);
                      }
                    }

                    MPI_Isend(solu2.data(), solu2.size(), MPI_DOUBLE, kproc, 0, PETSC_COMM_WORLD, &reqs[0]);
                    MPI_Isend(sysDofs.data(), sysDofs.size(), MPI_UNSIGNED, kproc, 1, PETSC_COMM_WORLD, &reqs[1]);
                    for(unsigned k = 0; k < dim; k++) {
                      MPI_Isend(vx2[k].data(), vx2[k].size(), MPI_DOUBLE, kproc, 2 + k, PETSC_COMM_WORLD, &reqs[2 + k]);
                    }
                  }

                  MPI_Status status;
                  for(unsigned m = 0; m < 2 + dim; m++) {
                    MPI_Wait(&reqs[m], &status);
                  }

                  if(iproc == kproc) {
                      
                    //**now we are on iel and on a particular face and corresponding jel. we have all the information.
                    aRes1.assign(nDofs1, 0.);
                    aRes2.assign(nDofs2, 0.);

                    s.new_recording();
                    //**gathering information on the boundary
                    const unsigned faceGeom = msh->GetElementFaceType(iel, iface);
                    unsigned faceDofs = msh->GetElementFaceDofNumber(iel, iface, solType);
                    std::vector  < std::vector  <  double > > faceVx(dim);    // A matrix holding the face coordinates rowwise.
                    for(int k = 0; k < dim; k++) {
                      faceVx[k].resize(faceDofs);
                    }
                    for(unsigned i = 0; i < faceDofs; i++) {
                      unsigned inode = msh->GetLocalFaceVertexIndex(iel, iface, i);    // face-to-element local node mapping.
                      for(unsigned k = 0; k < dim; k++) {
                        faceVx[k][i] =  vx1[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
                      }
                    }

                    double h = sqrt(dim - 1) * sqrt((faceVx[0][0] - faceVx[0][1]) * (faceVx[0][0] - faceVx[0][1]) +
                                                    (faceVx[1][0] - faceVx[1][1]) * (faceVx[1][0] - faceVx[1][1])) ;
                    double h2 = h * h;
                    double h3 = h * h * h;


                    //if(!aP1IsInitialized) { //build the basis 1,x,y,z... corresponding to the solution type
                      for(unsigned jtype = 0; jtype < solType + 1; jtype++) {
                        ProjectNodalToPolynomialCoefficients(aP1[jtype], vx1, ielt1, jtype);
                      }
                      //aP1IsInitialized = true;
                    //}

                    for(unsigned jtype = 0; jtype < solType + 1; jtype++) {
                      ProjectNodalToPolynomialCoefficients(aP2[jtype], vx2, ielt2, jtype);
                    }
                    //**gauss loop on the face
                    for(unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solType]->GetGaussPointNumber(); ig++) {

                      std::vector < double> normal;//**normal is going from iel to jel
                      msh->_finiteElement[faceGeom][solType]->JacobianSur(faceVx, ig, weight, phi, gradPhi, normal);//**here gradient for the test function is completely fake. important phi on the edge. we use phi to get physical gauss coordinates

                      std::vector< double > xg(dim, 0.); // physical coordinates of the face Gauss point
                      for(unsigned i = 0; i < faceDofs; i++) {
                        for(unsigned k = 0; k < dim; k++) {
                          xg[k] += phi[i] * faceVx[k][i];
                        }
                      }
                      //**now we can do the inverse mapping
                      std::vector <double> xi1;//local coordinates of the face gauss point with respect to iel
                      GetClosestPointInReferenceElement(vx1, xg, ielt1, xi1);

                      bool inverseMapping = GetInverseMapping(solType, ielt1, aP1, xg, xi1, 100);
                      if(!inverseMapping) {
                        std::cout << "InverseMapping1 failed at " << iel << " " << jel << " " << iface << std::endl;
                      }

                      std::vector <double> xi2;//local coordinates of the face gauss point with respect to jel
                      GetClosestPointInReferenceElement(vx2, xg, ielt2, xi2);

                      inverseMapping = GetInverseMapping(solType, ielt2, aP2, xg, xi2, 100);
                      if(!inverseMapping) {
                        std::cout << "InverseMapping2 failed at " << iel << " " << jel << " " << iface << std::endl;
                      }

                      //** having xi1 and xi2, now we can get gradient of laplace and the test function
                      msh->_finiteElement[ielt1][solType]->Jacobian(vx1, xi1, weight1, phi1, gradPhi1, nablaPhi1);
                      msh->_finiteElement[ielt2][solType]->Jacobian(vx2, xi2, weight2, phi2, gradPhi2, nablaPhi2);

                      //**now you have surface are and normal. now we can assemble penalty term.

                      adept::adouble gradSolu1DotN = 0.;//**gradient of the solution one in normal direction
                      adept::adouble hessSolu1DotN = 0.;//** Hessian of the solution in the normal direction
                      for(unsigned i = 0; i < nDofs1; i++) {
                        for(unsigned J = 0; J < dim; J++) {

                          //**dot product between the gradient and the normal
                          //**gradient of phi in the direction j multiply by normal j
                          gradSolu1DotN += solu1[i] * gradPhi1[i * dim + J] * normal[J];
                          for(unsigned K = 0; K < dim; K++) {
                            //2D xx, yy, xy
                            //3D xx, yy, zz, xy, yz ,zx
                            unsigned L;
                            if(J == K) L = J;
                            else if(1 == J + K) L = dim;     // xy
                            else if(2 == J + K) L = dim + 2; // xz
                            else if(3 == J + K) L = dim + 1; // yz
                            hessSolu1DotN += solu1[i] * normal[J] * nablaPhi1[i * dim2 + L] * normal[K];
                          }
                        }
                      }
//** same thing should be done for element 2
                      adept::adouble gradSolu2DotN = 0.;
                      adept::adouble hessSolu2DotN = 0.;
                      for(unsigned i = 0; i < nDofs2; i++) {
                        for(unsigned J = 0; J < dim; J++) {
                          gradSolu2DotN += solu2[i] * gradPhi2[i * dim + J] * normal[J];
                          for(unsigned K = 0; K < dim; K++) {
                            //2D xx, yy, xy
                            //3D xx, yy, zz, xy, yz ,zx
                            unsigned L;
                            if(J == K) L = J;
                            else if(1 == J + K) L = dim;     // xy
                            else if(2 == J + K) L = dim + 2; // xz
                            else if(3 == J + K) L = dim + 1; // yz
                            hessSolu2DotN += solu2[i] * normal[J] * nablaPhi2[i * dim2 + L] * normal[K];
                          }
                        }
                      }


                      double C1 = 0.05 * alpha;  // [alpha];

                      for(unsigned i = 0; i < nDofs1; i++) {

                        for(unsigned J = 0; J < dim; J++) {
                          aRes1[i] +=  C1 * h * gradPhi1[i * dim + J] * normal[J] * (gradSolu1DotN - gradSolu2DotN) * weight;



                          for(unsigned K = 0; K < dim; K++) {

                            unsigned L;
                            if(J == K) L = J;
                            else if(1 == J + K) L = dim;     // xy
                            else if(2 == J + K) L = dim + 2; // xz
                            else if(3 == J + K) L = dim + 1; // yz

                            aRes1[i] += C1 * h3 * normal[J] * nablaPhi1[i * dim2 + L] * normal[K] * (hessSolu1DotN - hessSolu2DotN) * weight;
                          }
                        }
                      }

                      for(unsigned i = 0; i < nDofs2; i++) {

                        for(unsigned J = 0; J < dim; J++) {
                          aRes2[i] +=  -C1 * h * gradPhi2[i * dim + J] * normal[J] * (gradSolu1DotN - gradSolu2DotN) * weight;

                          for(unsigned K = 0; K < dim; K++) {

                            unsigned L;
                            if(J == K) L = J;
                            else if(1 == J + K) L = dim;     // xy
                            else if(2 == J + K) L = dim + 2; // xz
                            else if(3 == J + K) L = dim + 1; // yz

                            aRes2[i] +=  -C1 * h3 * normal[J] * nablaPhi2[i * dim2 + L] * normal[K] * (hessSolu1DotN - hessSolu2DotN) * weight;
                          }
                        }
                      }


                    }

                    //copy the value of the adept::adoube aRes in double Res and store them in RES
                    rhs.resize(nDofs1 + nDofs2);   //resize
                    for(int i = 0; i < nDofs1; i++) {
                      rhs[ i ] = -aRes1[i].value();
                    }

                    for(int i = 0; i < nDofs2; i++) {
                      rhs[nDofs1 + i] = -aRes2[i].value();
                    }
                    myRES->add_vector_blocked(rhs, sysDofs);

                    s.dependent(&aRes1[0], nDofs1);
                    s.dependent(&aRes2[0], nDofs2);
                    s.independent(&solu1[0], nDofs1);
                    s.independent(&solu2[0], nDofs2);

                    Jac.resize((nDofs1 + nDofs2) * (nDofs1 + nDofs2));
                    s.jacobian(&Jac[0], true);
                    myKK->add_matrix_blocked(Jac, sysDofs, sysDofs);

                    s.clear_independents();
                    s.clear_dependents();

                  }
                }
              }
            }
          }
        }
      }
    }
  }
  // *************************************
  std::cout << "Ghost Penalty Assembly time = " << static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC << std::endl;

}




