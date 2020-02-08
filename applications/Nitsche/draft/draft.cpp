// I will use this file as a draft for my silly experiments

#include <iostream>
#include <vector>


void PrintArray(const std::vector<double>& array);
void PrintMat(std::vector< std::vector<double> >& M);

int main (int argc, char** args) {
    unsigned dim = 2;
    unsigned nDofD = 4;
//     std::vector <double> array = {1,2,3,4,5}; // std::vector provides real time allocation
//     std::cout << "Before resize, lenght is: " << array.size() << std::endl;
//     PrintArray(array);
//     std::cout << " "<<std::endl;
//     array.resize(7); // resize the lenght with 0 fill
//     std::cout << "After resize, lenght is now: " << array.size() << std::endl;
//     PrintArray(array);
//     std::vector<double> x2;
//     x2.assign(3,11);
//     PrintArray(x2);
//     std::cout << " "<<std::endl;
    std::vector<std::vector<double> > M(dim);
    std::cout << "Before resize " << M.size();
    PrintMat(M);
    for(int k = 0; k < dim; k++) {
      M[k].resize(nDofD);
      //aResD1[k].assign(nDofD, 0.);    //resize
    }
    std::cout << "After resize " << M.size() << std::endl;
    PrintMat(M);
    
    
//     for (int k = 0; k < M.size(); k++) {
//        M[k].resize (7);
//     }
//     std::cout << "After resize" ;
//     PrintMat(M);
//     
    
    
    

  return 0;
}



void PrintArray(const std::vector<double>& array)
{
    for(unsigned i = 0; i < array.size() ; i++){
        std::cout << array[i] << ' ';
    }
}


void PrintMat(std::vector< std::vector<double> >& M){
    
    for(unsigned i=0;i<M.size();i++){
        for(unsigned j=0;j<M[i].size();j++){
            std::cout<<M[i][j]<<" ";
        }
        std::cout<< std::endl;
    }
    
}




// void AssembleNitscheProblem_AD (MultiLevelProblem& ml_prob) {
//   //  ml_prob is the global object from/to where get/set all the data
//   //  level is the level of the PDE system to be assembled
//   //  levelMax is the Maximum level of the MultiLevelProblem
//   //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled
// 
//   // call the adept stack object
// 
// 
//   adept::Stack& s = FemusInit::_adeptStack;
// 
//   //  extract pointers to the several objects that we are going to use
// 
//   LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Nitsche");   // pointer to the linear implicit system named "Poisson"
//   const unsigned level = mlPdeSys->GetLevelToAssemble(); // We have different level of meshes. we assemble the problem on the specified one.
// 
//   Mesh*                    msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
//   elem*                     el = msh->el;  // pointer to the elem object in msh (level)
// 
//   MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
//   Solution*                sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object
// 
//   LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
//   SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
//   NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)
// 
//   const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
// 
//   unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
// 
//   //solution variable
//   unsigned solu1Index = mlSol->GetIndex ("u1");   // get the position of "u" in the ml_sol object
//   unsigned solu2Index = mlSol->GetIndex ("u2");   // get the position of "u" in the ml_sol object
//   unsigned soluType = mlSol->GetSolutionType (solu1Index);   // get the finite element type for "u"
// 
//   unsigned solu1PdeIndex;
//   solu1PdeIndex = mlPdeSys->GetSolPdeIndex ("u1");   // get the position of "u" in the pdeSys object
//   unsigned solu2PdeIndex;
//   solu2PdeIndex = mlPdeSys->GetSolPdeIndex ("u2");   // get the position of "u" in the pdeSys object
// 
//   vector < adept::adouble >  solu1; // local solution
//   vector < adept::adouble >  solu2; // local solution
// 
//   unsigned eflagIndex = mlSol->GetIndex ("eflag");
//   unsigned nflagIndex = mlSol->GetIndex ("nflag");
// 
//   vector < unsigned >  nodeFlag; // local solution
// 
//   vector < vector < double > > x (dim);   // local coordinates. x is now dim x m matrix.
//   unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
// 
//   vector <double> phi;  // local test function
//   vector <double> phi_x; // local test function first order partial derivatives
//   double weight; // gauss point weight
// 
//   vector< adept::adouble > aResu1; // local redidual vector
//   vector< adept::adouble > aResu2; // local redidual vector
// 
//   vector< unsigned > l2GMap; // local to global mapping
//   vector< double > Res; // local redidual vector
//   vector < double > Jac;
// 
//   KK->zero(); // Set to zero all the entries of the Global Matrix
//   RES->zero(); // Set to zero all the entries of the Global Residual
// 
//   std::vector < std::vector < std::vector <double > > > aP (3);
// 
//   for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
// 
//     short unsigned ielGeom = msh->GetElementType (iel);
//     unsigned nDofu  = msh->GetElementDofNumber (iel, soluType); // number of solution element dofs
// 
//     unsigned eFlag = static_cast <unsigned> (floor ( (*sol->_Sol[eflagIndex]) (iel) + 0.5));
// 
//     // resize local arrays
//     l2GMap.resize (2 * nDofu);
//     solu1.resize (nDofu);
//     solu2.resize (nDofu);
//     nodeFlag.resize (nDofu);
// 
//     for (int k = 0; k < dim; k++) {
//       x[k].resize (nDofu);
//     }
// 
//     aResu1.assign (nDofu, 0.);   //resize
//     aResu2.assign (nDofu, 0.);   //resize
// 
//     // local storage of global mapping and solution
//     for (unsigned i = 0; i < nDofu; i++) {
//       unsigned solDof = msh->GetSolutionDof (i, iel, soluType);   // global to global mapping between solution node and solution dof
//       solu1[i] = (*sol->_Sol[solu1Index]) (solDof);                 // global extraction and local storage for the solution
//       solu2[i] = (*sol->_Sol[solu2Index]) (solDof);                 // global extraction and local storage for the solution
//       nodeFlag[i] = (*sol->_Sol[nflagIndex]) (solDof);
//       l2GMap[i] = pdeSys->GetSystemDof (solu1Index, solu1PdeIndex, i, iel);   // global to global mapping between solution node and pdeSys dof
//       l2GMap[nDofu + i] = pdeSys->GetSystemDof (solu2Index, solu2PdeIndex, i, iel);   // global to global mapping between solution node and pdeSys dof
//     }
// 
//     // local storage of coordinates
//     for (unsigned i = 0; i < nDofu; i++) {
//       unsigned xDof  = msh->GetSolutionDof (i, iel, 2);   // global to global mapping between coordinates node and coordinate dof
//       for (unsigned k = 0; k < dim; k++) {
//         x[k][i] = (*msh->_topology->_Sol[k]) (xDof);     // global extraction and local storage for the element coordinates
//       }
//     }
// 
// 
//     double alpha1 = .5;
//     double alpha2 = 3.;
// 
//     // start a new recording of all the operations involving adept::adouble variables
//     s.new_recording();
//     if (eFlag == 0 || eFlag == 2) {
//       // *** Element Gauss point loop ***
//       for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
//         // *** get gauss point weight, test function and test function partial derivatives ***
//         msh->_finiteElement[ielGeom][soluType]->Jacobian (x, ig, weight, phi, phi_x);
// 
//         // evaluate the solution, the solution derivatives and the coordinates in the gauss point
//         vector < adept::adouble > gradSolu1g (dim, 0.);
//         vector < adept::adouble > gradSolu2g (dim, 0.);
// 
//         for (unsigned i = 0; i < nDofu; i++) {
//           for (unsigned k = 0; k < dim; k++) {
//             gradSolu1g[k] += phi_x[i * dim + k] * solu1[i];
//             gradSolu2g[k] += phi_x[i * dim + k] * solu2[i];
//           }
//         }
// 
//         // *** phi_i loop ***
//         for (unsigned i = 0; i < nDofu; i++) {
// 
//           adept::adouble graduGradphi1 = 0.;
//           adept::adouble graduGradphi2 = 0.;
// 
//           for (unsigned k = 0; k < dim; k++) {
//             graduGradphi1 += gradSolu1g[k] * phi_x[i * dim + k];
//             graduGradphi2 += gradSolu2g[k] * phi_x[i * dim + k];
//           }
// 
// 
//           if (eFlag == 0) {
//             aResu1[i] += (- phi[i] + alpha1 * graduGradphi1) * weight;
//             if (nodeFlag[i] != 1) {
//               aResu2[i] += (graduGradphi2) * weight;
//             }
//           }
//           else if (eFlag == 2) {
//             aResu2[i] += (- phi[i] + alpha2 * graduGradphi2) * weight;
//             if (nodeFlag[i] != 1) {
//               aResu1[i] += (graduGradphi1) * weight;
//             }
//           }
// 
//         } // end phi_i loop
//       } // end gauss point loop
//     }
// 
//     else {
// 
//       // *** Element Gauss point loop ***
//       
//       std::string elemGeometry;
//       if (dim == 1) elemGeometry = "line";
//       else if (dim == 2) elemGeometry = "quad";
//       else if (dim == 3) elemGeometry = "hex";
// 
//       Gauss gpBulk (elemGeometry.c_str(), "seventh");
// 
//       for (unsigned jtype = 0; jtype < soluType + 1; jtype++) {
//         ProjectNodalToPolynomialCoefficients (aP[jtype], x, ielGeom, jtype) ;
//       }
// 
//       //bulk1
//       std::vector <double> N (dim, 0.);
//       N[0] = -1.;
// 
//       std::vector <double> H (dim);
//       std::vector <double> XC (dim);
//       if (dim == 1) {
//         H[0] = fabs (x[0][1] - x[0][0]);
//         XC[0] = 0.5 * (x[0][1] + x[0][0]);
//       }
//       else if (dim == 2) {
//         H[0] = fabs (x[0][2] - x[0][0]);
//         H[1] = fabs (x[1][2] - x[1][0]);
// 
//         XC[0] = 0.5 * (x[0][2] + x[0][0]);
//         XC[1] = 0.5 * (x[1][2] + x[1][0]);
//       }
//       else if (dim == 3) {
//         H[0] = fabs (x[0][6] - x[0][0]);
//         H[1] = fabs (x[1][6] - x[1][0]);
//         H[2] = fabs (x[2][6] - x[2][0]);
// 
//         XC[0] = 0.5 * (x[0][6] + x[0][0]);
//         XC[1] = 0.5 * (x[1][6] + x[1][0]);
//         XC[2] = 0.5 * (x[2][6] + x[2][0]);
//       }
// 
// 
//       unsigned NGPB = gpBulk.GetGaussPointsNumber();
//       std::vector < const double * > xig (dim);
//       for (unsigned k = 0; k < dim; k++) {
//         xig[k] = gpBulk.GetGaussCoordinatePointer(k);
//       }
// 
//       for (unsigned ig = 0; ig < NGPB; ig++) {
//         std::vector < double > xg (dim);
//         xg[0] = (XC[0] - 0.25 * H[0]) + (xig[0][ig]) * 0.25 * H[0] ;
//         for (unsigned k = 1; k < dim; k++) {
//           xg[k] = XC[k] + (xig[k][ig]) * 0.5 * H[k] ;
//         }
// 
//         std::vector <double> xi (dim, 0.);
//         GetClosestPointInReferenceElement (x, xg, ielGeom, xi);
//         bool inverseMapping = GetInverseMapping (soluType, ielGeom, aP, xg, xi, 100);
// 
//         msh->_finiteElement[ielGeom][soluType]->Jacobian (x, xi, weight, phi, phi_x);
//         double Area = 0.5 * H[0];
//         double Area0 = 2.;
//         for (unsigned k = 1; k < dim; k++) {
//           Area *= H[k];
//           Area0 *= 2.;
//         }
//         weight = gpBulk.GetGaussWeight (ig) * Area / Area0;
// 
//         // evaluate the solution, the solution derivatives and the coordinates in the gauss point
//         vector < adept::adouble > gradSolu1g (dim, 0.);
//         for (unsigned i = 0; i < nDofu; i++) {
//           for (unsigned k = 0; k < dim; k++) {
//             gradSolu1g[k] += phi_x[i * dim + k] * solu1[i];
//           }
//         }
// 
//         // *** phi_i loop ***
//         for (unsigned i = 0; i < nDofu; i++) {
//           adept::adouble graduGradphi1 = 0.;
//           for (unsigned k = 0; k < dim; k++) {
//             graduGradphi1 += gradSolu1g[k] * phi_x[i * dim + k];
//           }
//           aResu1[i] += (-phi[i] + alpha1 * graduGradphi1) * weight;
//         }
//       }
// 
//       // interface
//       std::string InterfaceGeometry;
//       if (dim == 1) InterfaceGeometry = "point";
//       else if (dim == 2) InterfaceGeometry = "line";
//       else if (dim == 3) InterfaceGeometry = "quad";
// 
//       Gauss gpInt (elemGeometry.c_str(), "seventh");
//       unsigned NGPI = gpInt.GetGaussPointsNumber();
//       for (unsigned k = 1; k < dim; k++) {
//         xig[k] = gpInt.GetGaussCoordinatePointer(k - 1);
//       }
// 
//       for (unsigned ig = 0; ig < gpInt.GetGaussPointsNumber(); ig++) {
//         std::vector < double > xg (dim);
// 
//         xg[0] = XC[0];
// 
//         for (unsigned k = 1; k < dim; k++) {
//           xg[k] = XC[k] + xig[k][ig] * 0.5 * H[k] ;
//         }
// 
//         std::vector <double> xi (dim, 0.);
//         GetClosestPointInReferenceElement (x, xg, ielGeom, xi);
//         bool inverseMapping = GetInverseMapping (soluType, ielGeom, aP, xg, xi, 100);
// 
//         msh->_finiteElement[ielGeom][soluType]->Jacobian (x, xi, weight, phi, phi_x);
//         double Area = 1.;
//         double Area0 = 1.;
//         for (unsigned k = 1; k < dim; k++) {
//           Area *= H[k];
//           Area0 *= 2.;
//         }
//         weight = gpInt.GetGaussWeight (ig) * Area / Area0;
// 
//         // evaluate the solution, the solution derivatives and the coordinates in the gauss point
// 
//         double theta = 1.;
//         double gamma1 = 0.5;
//         double gamma2 = 0.5;
// 
//         adept::adouble solu1g  = 0.;
//         adept::adouble solu2g  = 0.;
//         adept::adouble alphaGradSoluDotN = 0.;
// 
//         for (unsigned i = 0; i < nDofu; i++) {
//           solu1g += phi[i] * solu1[i];
//           solu2g += phi[i] * solu2[i];
//           for (unsigned k = 0; k < dim; k++) {
//             alphaGradSoluDotN += (alpha1 * gamma1  * solu1[i] + alpha2 * gamma2 * solu2[i]) * phi_x[i * dim + k] * N[k];
//           }
//         }
// 
//         // *** phi_i loop ***
//         for (unsigned i = 0; i < nDofu; i++) {
// 
//           adept::adouble gradPhiDotN = 0.;
//           for (unsigned k = 0; k < dim; k++) {
//             gradPhiDotN += phi_x[i * dim + k] * N[k];
//           }
// 
//           aResu1[i] += (solu2g - solu1g) * (-phi[i] * theta - alpha1 * gamma1 * gradPhiDotN) * weight;
//           aResu1[i] += alphaGradSoluDotN * (+phi[i]) * weight;
//           aResu2[i] += (solu2g - solu1g) * (+phi[i] * theta - alpha2 * gamma2 * gradPhiDotN) * weight;
//           aResu2[i] += alphaGradSoluDotN * (-phi[i]) * weight;
//         } // end phi_i loop
// 
//       }
// 
//       //bulk2
//       for (unsigned k = 0; k < dim; k++) {
//         xig[k] = gpBulk.GetGaussCoordinatePointer(k);
//       }
// 
//       for (unsigned ig = 0; ig < gpBulk.GetGaussPointsNumber(); ig++) {
//         std::vector < double > xg (dim);
//         xg[0] = (XC[0] + 0.25 * H[0]) + xig[0][ig] * 0.25 * H[0] ;
//         for (unsigned k = 1; k < dim; k++) {
//           xg[k] = XC[k] + xig[k][ig] * 0.5 * H[k] ;
//         }
// 
//         std::vector <double> xi (dim, 0.);
//         GetClosestPointInReferenceElement (x, xg, ielGeom, xi);
//         bool inverseMapping = GetInverseMapping (soluType, ielGeom, aP, xg, xi, 100);
// 
//         msh->_finiteElement[ielGeom][soluType]->Jacobian (x, xi, weight, phi, phi_x);
//         double Area = 0.5 * H[0];
//         double Area0 = 2.;
//         for (unsigned k = 1; k < dim; k++) {
//           Area *= H[k];
//           Area0 *= 2.;
//         }
//         weight = gpBulk.GetGaussWeight (ig) * Area / Area0;
// 
//         vector < adept::adouble > gradSolu2g (dim, 0.);
//         for (unsigned i = 0; i < nDofu; i++) {
//           for (unsigned k = 0; k < dim; k++) {
//             gradSolu2g[k] += phi_x[i * dim + k] * solu2[i];
//           }
//         }
// 
//         // *** phi_i loop ***
//         for (unsigned i = 0; i < nDofu; i++) {
//           adept::adouble graduGradphi2 = 0.;
//           for (unsigned k = 0; k < dim; k++) {
//             graduGradphi2 += gradSolu2g[k] * phi_x[i * dim + k];
//           }
//           aResu2[i] += (-phi[i] + alpha2 * graduGradphi2) * weight;
//         }
//       }
//     }
// 
//     //copy the value of the adept::adoube aRes in double Res and store
//     Res.resize (2 * nDofu);   //resize
// 
//     for (int i = 0; i < nDofu; i++) {
//       Res[i] = - aResu1[i].value();
//       Res[nDofu + i] = - aResu2[i].value();
//     }
// 
//     RES->add_vector_blocked (Res, l2GMap);
// 
//     // define the dependent variables
//     s.dependent (&aResu1[0], nDofu);
//     s.dependent (&aResu2[0], nDofu);
// 
//     // define the independent variables
//     s.independent (&solu1[0], nDofu);
//     s.independent (&solu2[0], nDofu);
// 
//     // get the jacobian matrix (ordered by row major )
//     Jac.resize (2 * nDofu * 2 * nDofu);   //resize
//     s.jacobian (&Jac[0], true);
// 
//     //store K in the global matrix KK
//     KK->add_matrix_blocked (Jac, l2GMap, l2GMap);
// 
//     s.clear_independents();
//     s.clear_dependents();
// 
//   } //end element loop for each process
// 
//   RES->close();
//   KK->close();
// 
// //   PetscViewer    viewer;
// //   //PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
// //   //PetscObjectSetName ( (PetscObject) viewer, "FSI matrix");
// //   //PetscViewerPushFormat (viewer, PETSC_VIEWER_DRAW_LG);
// //   MatView ( (static_cast<PetscMatrix*> (KK))->mat(), viewer);
// //   VecView ( (static_cast<PetscVector*> (RES))->vec(), viewer);
// 
//   //double a;
//   //std::cin >> a;
// 
// 
// 
//   //***************** END ASSEMBLY *******************
// }


