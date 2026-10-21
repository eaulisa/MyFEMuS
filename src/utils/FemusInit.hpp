/*=========================================================================

 Program: FEMUS
 Module: FemusInit
 Authors: Simone Bnà, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_utils_FemusInit_hpp__
#define __femus_utils_FemusInit_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FemusConfig.hpp"

//use this class only with petsc-mpi
#ifdef HAVE_MPI
#include <mpi.h>  //For MPI_COMM_WORLD
#endif

// PETSC ----------------------
#ifdef HAVE_PETSC
// EXTERN_C_FOR_PETSC_BEGIN
# include <petsc.h>
# include <petscerror.h>
// EXTERN_C_FOR_PETSC_END
#endif

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// ========================================
/// Class FemusInit
// ========================================

#include "adept.h"
#include "uq.hpp"

namespace femus {

  class FemusInit {

    public:

      /// Constructor
      FemusInit(int & argc, char** & argv, MPI_Comm comm_world_in = MPI_COMM_WORLD);

      /// Destructor
      ~FemusInit();

      static adept::Stack _adeptStack;
      static uq _uqHermite;
      static uq _uqLegendre;

  };

  inline void PrintProcessMemory(const char* checkpoint, const unsigned iproc) {
    std::ifstream status("/proc/self/status");
    std::string line;
    unsigned long rss = 0;
    unsigned long peak = 0;

    if (!status) return;

    while (std::getline(status, line)) {
      std::istringstream entry(line);
      std::string key;
      entry >> key;
      if (key == "VmRSS:") entry >> rss;
      else if (key == "VmHWM:") entry >> peak;
    }

    std::ostringstream message;
    message << "[RAM] proc=" << iproc << " checkpoint=" << checkpoint << " RSS=" << rss / 1024.0 << " MiB peak=" << peak / 1024.0 << " MiB\n";
    std::cerr << message.str() << std::flush;
  }

  inline void Print(std::string s, bool dbg = false) {
    int iproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    for (unsigned i = 0; i < 30; i++)
      (dbg) ? std::cerr << s : std::cout << s;
    (dbg) ? std::cerr << " [" << iproc << "]" << std::endl : std::cout << std::endl;
  }

} //end namespace femus

#endif // end _femus_init ----------------------
