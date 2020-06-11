#include <iostream>
#include <math.h>
#include <stdio.h>
#include "mpi.h"


// MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
// MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,MPI_Comm comm, MPI_Status * status)


int main(int argc, char** argv) {

  int my_PE_num, numbertoreceive, numbertosend = 42;
  MPI_Status status;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_PE_num);

  if(my_PE_num == 0) {

    MPI_Recv(&numbertoreceive, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    std::cout << "NumberReceived = " << numbertoreceive << std::endl;
  }
  else {
    std::cout << "I am PE = " << my_PE_num << std::endl;
    MPI_Send(&numbertosend, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
  }
  MPI_Finalize();

}



