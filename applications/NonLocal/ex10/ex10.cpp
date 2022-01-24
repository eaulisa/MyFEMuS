
#include "mpi.h"
#include <iostream>

int main(int argc, char* argv[]) {

  // init Petsc-MPI communicator
  MPI_Init(&argc, &argv);
  int p;
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);   
  if(myid == 0){
    std::cout <<"The total number of Processes is " << p <<std::endl;
  }
  MPI_Finalize();

  return 0;

} //end main
