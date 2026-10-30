#pragma once

void GetGPSparsityPattern (
  MultiLevelProblem& ml_prob,
  const std::string& CName,
  const std::string& SolName,
  const unsigned level,
  const unsigned levelC,
  const unsigned phase, // 0 -> Velocity, 1 -> P1, 2 -> P2
  vector < std::map < int, bool > >& BlgToMe_d,
  vector < std::map < int, bool > >& BlgToMe_o,
  std::map < int, std::map <int, bool > >& DnBlgToMe_o,
  std::map < int, std::map <int, bool > >& DnBlgToMe_d
) {

  MultiLevelSolution* mlSol = ml_prob._ml_sol;
  Solution* sol = mlSol->GetSolutionLevel(level);

  Mesh* msh = mlSol->GetMultilevelMesh()->GetLevel(level);
  elem* el = msh->el;

  unsigned iproc = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  const unsigned cIndex = mlSol->GetIndex(CName.c_str());

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("NS");
  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];
  const auto KKoffset = myLinEqSolver->GetOffset();
  const auto KKIndex = myLinEqSolver->GetIndex();
  int IndexStart = KKoffset[0][iproc];
  int IndexEnd  = KKoffset[KKIndex.size() - 1][iproc];
  int owned_dofs    = IndexEnd - IndexStart;

  unsigned indexSol = mlSol->GetIndex(SolName.c_str());
  unsigned indexPde = my_nnlin_impl_sys.GetSolPdeIndex(SolName.c_str());
  unsigned solType = mlSol->GetSolutionType(SolName.c_str());

  std::vector <unsigned> sysDofs1;
  std::vector <unsigned> sysDofs2;

  auto flag = [&](unsigned iellevel, double Ciel) -> bool {
    if (iellevel != levelC) return false;

    if (phase == 0) {
      return true;
    }
    else if (phase == 1) {
      return (Ciel > 0.1 );
    }
    else if (phase == 2) {
      return (Ciel < 0.9 );
    }
    else {
      throw std::runtime_error("WRONG PHASE FOR GETGPFACES()");
    }
  };

  //flagmark
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    double Ciel = (*sol->_Sol[cIndex])(iel);

    unsigned iel_level = msh->el->GetElementLevel(iel);

    if (flag(iel_level, Ciel)) {

      unsigned nDofs1 = msh->GetElementDofNumber(iel, solType);
      sysDofs1.resize(nDofs1);

      for(unsigned i = 0; i < nDofs1; i++) {
        sysDofs1[i] = myLinEqSolver->GetSystemDof(indexSol, indexPde, i, iel);
      }

      for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
        int jel = el->GetFaceElementIndex(iel, iface) - 1;
        if(jel >= 0) { // iface is not a boundary of the domain

          unsigned jproc = msh->IsdomBisectionSearch(jel, 3);

          if(jproc == iproc) {

            unsigned jel_level = msh->el->GetElementLevel(jel);

            double Cjel = (*sol->_Sol[cIndex])(jel);

            if (flag(jel_level, Cjel)) {

              unsigned nDofs2 = msh->GetElementDofNumber(jel, solType);
              sysDofs2.resize(nDofs2);

              for(unsigned i = 0; i < nDofs2; i++) {
                sysDofs2[i] = myLinEqSolver->GetSystemDof(indexSol, indexPde, i, jel);
              }

              for(int inode = 0; inode < nDofs1; inode++) {
                // identify the process the i-row belogns to
                int iiproc = 0;
                while(sysDofs1[inode] >= KKoffset[KKIndex.size() - 1][iiproc]) iiproc++;

                for(int jnode = 0; jnode < nDofs2; jnode++) {
                  // identify the process the j-column belogns to
                  int jjproc = 0;
                  while(sysDofs2[jnode] >= KKoffset[KKIndex.size() - 1][jjproc]) jjproc++;

                  if(iiproc == iproc) {  // i-row belogns to this proc
                    if(jjproc == iproc) {  // j-column belongs to this proc (diagonal)
                      BlgToMe_d[ sysDofs1[inode] - IndexStart ][ sysDofs2[jnode] - IndexStart ] = 1;
                    }
                    else { // j-column does not belong to this proc (off-diagonal)
                      BlgToMe_o[ sysDofs1[inode] - IndexStart ][ sysDofs2[jnode] ] = 1;
                    }
                  }
                  else { // i-row does not belong to this proc
                    if(iiproc != jjproc) {  // if diagonal
                      DnBlgToMe_o[sysDofs1[inode]][sysDofs2[jnode]] = 1;
                    }
                    else {  // if off-diagonal
                      DnBlgToMe_d[sysDofs1[inode]][sysDofs2[jnode]] = 1;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if(nprocs > 1) {
    for(unsigned kproc = 0; kproc < nprocs; kproc++) {
      for(int iel = msh->_elementOffset[kproc]; iel < msh->_elementOffset[kproc + 1]; iel++) {

        unsigned eFlag1 = 0;
        if(iproc == kproc) {
          double Ciel = (*sol->_Sol[cIndex])(iel);
          unsigned iel_level = msh->el->GetElementLevel(iel);
          eFlag1 = static_cast<unsigned>(flag(iel_level, Ciel));
        }
        MPI_Bcast(&eFlag1, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);

        if(eFlag1 > 0) {
          unsigned nFaces;
          if(iproc == kproc) {
            nFaces = msh->GetElementFaceNumber(iel);
          }
          MPI_Bcast(&nFaces, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);

          for(unsigned iface = 0; iface < nFaces; iface++) {

            int jel;
            if(iproc == kproc) {
              jel = el->GetFaceElementIndex(iel, iface) - 1;
            }
            MPI_Bcast(&jel, 1, MPI_INT, kproc, PETSC_COMM_WORLD);

            if(jel >= 0) { // iface is not a boundary of the domain
              unsigned jproc = msh->IsdomBisectionSearch(jel, 3);  // return  jproc for piece-wise constant discontinuous type (3)
              if(jproc != kproc && (iproc == kproc || iproc == jproc)) {

                unsigned eFlag2;
                if(iproc == jproc) {
                  double Cjel = (*sol->_Sol[cIndex])(jel);
                  unsigned jel_level = msh->el->GetElementLevel(jel);
                  eFlag2 = static_cast<unsigned>(flag(jel_level, Cjel));
                  MPI_Send(&eFlag2, 1, MPI_UNSIGNED, kproc, 0, PETSC_COMM_WORLD);
                }
                else if(iproc == kproc) {
                  MPI_Recv(&eFlag2, 1, MPI_UNSIGNED, jproc, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
                }

                if(eFlag2 > 0) {
                  short unsigned ielt1;
                  short unsigned ielt2;

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

                  sysDofs2.resize(nDofs2);
                  std::vector < MPI_Request > reqs(1);

                  if(iproc == kproc) {

                    nDofs1 = el->GetNVE(ielt1, solType);

                    sysDofs1.resize(nDofs1);
                    for(unsigned i = 0; i < nDofs1; i++) {
                      sysDofs1[i] = myLinEqSolver->GetSystemDof(indexSol, indexPde, i, iel);
                    }

                    MPI_Irecv(sysDofs2.data(), sysDofs2.size(), MPI_UNSIGNED, jproc, 0, PETSC_COMM_WORLD,  &reqs[0]);
                  }
                  else if(iproc == jproc) {
                    for(unsigned i = 0; i < nDofs2; i++) {
                      sysDofs2[i] = myLinEqSolver->GetSystemDof(indexSol, indexPde, i, jel);
                    }

                    MPI_Isend(sysDofs2.data(), sysDofs2.size(), MPI_UNSIGNED, kproc, 0, PETSC_COMM_WORLD, &reqs[0]);
                  }

                  MPI_Status status;
                  for(unsigned m = 0; m < 1; m++) {
                    MPI_Wait(&reqs[m], &status);
                  }

                  if(iproc == kproc) {
                    for(int inode = 0; inode < nDofs1; inode++) {
                      // identify the process the i-row belogns to
                      int iiproc = 0;
                      while(sysDofs1[inode] >= KKoffset[KKIndex.size() - 1][iiproc]) iiproc++;

                      for(int jnode = 0; jnode < nDofs2; jnode++) {
                        // identify the process the j-column belogns to
                        int jjproc = 0;
                        while(sysDofs2[jnode] >= KKoffset[KKIndex.size() - 1][jjproc]) jjproc++;

                        if(iiproc == iproc) {  // i-row belogns to this proc
                          if(jjproc == iproc) {  // j-column belongs to this proc (diagonal)
                            BlgToMe_d[ sysDofs1[inode] - IndexStart ][ sysDofs2[jnode] - IndexStart ] = 1;
                          }
                          else { // j-column does not belong to this proc (off-diagonal)
                            BlgToMe_o[ sysDofs1[inode] - IndexStart ][ sysDofs2[jnode] ] = 1;
                          }
                        }
                        else { // i-row does not belong to this proc
                          if(iiproc != jjproc) {  // if diagonal
                            DnBlgToMe_o[sysDofs1[inode]][sysDofs2[jnode]] = 1;
                          }
                          else {  // if off-diagonal
                            DnBlgToMe_d[sysDofs1[inode]][sysDofs2[jnode]] = 1;
                          }
                        }
                      }
                    }
                  }

                }
              }
            }
          }
        }
      }
    }
  }

}
