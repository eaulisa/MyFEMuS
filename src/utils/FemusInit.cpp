/*=========================================================================

 Program: FEMUS
 Module: FemusInit
 Authors: Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include "FemusInit.hpp"
#include "UqQuadratureTypeEnum.hpp"

namespace femus {

  adept::Stack FemusInit::_adeptStack;

  uq FemusInit::_uqHermite(UQ_HERMITE);
  uq FemusInit::_uqLegendre(UQ_LEGENDRE);

// =======================================================
/// This function initializes the libraries if it is parallel
  FemusInit::FemusInit(
    int & argc,            // integer program input
    char** & argv,         // char program input
    MPI_Comm comm_world_in // communicator for MPI direct
  ) {// ======================================================

#ifdef HAVE_PETSC

    int ierr = PetscInitialize (&argc, &argv, NULL, NULL);    CHKERRABORT(PETSC_COMM_WORLD, ierr);

#endif

#ifdef HAVE_MPI
    // redirect libMesh::out to nothing on all
    // other processors unless explicitly told
    // not to via the --keep-cout command-line argument.

    int world_size;       //The number of processes that were spawned.
    int world_rank;       //The rank of the current process.
    char processor_name[MPI_MAX_PROCESSOR_NAME];  //The name of the current process
    int name_len;
    // Get and set the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get and set the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get and set the name and length of the processor
    MPI_Get_processor_name(processor_name, &name_len);

    std::ofstream fout;

    for(unsigned i = 0; i < world_size; i++) {
      if(i == world_rank) {
        fout.open("ProcessorTable.txt", std::ios_base::app);
        fout << "processor " << processor_name << " rank " << world_rank;
        fout << " of " << world_size << " processors" << std::endl;
        fout.close();
      }
    }

    if ( world_rank != 0) {
    if ( i != 0) {
      std::cout.rdbuf(NULL);
    }
#endif

    std::cout << " FemusInit(): PETSC_COMM_WORLD initialized" << std::endl << std::endl << std::flush;

    return;
  }


  FemusInit::~FemusInit() {

#ifdef HAVE_PETSC
    PetscFinalize();
    std::cout << std::endl << " ~FemusInit(): PETSC_COMM_WORLD ends" << std::endl;
#endif

    return;
  }


} //end namespace femus


