#include <iostream>
#include <math.h>
#include <stdio.h>
#include "mpi.h"
#include <vector>

// MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
// MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,MPI_Comm comm, MPI_Status * status)


int main(int argc, char** argv) {

//   int my_PE_num, numbertoreceive, numbertosend = 42;
//   MPI_Status status;
//
//   MPI_Init(&argc, &argv);
//
//   MPI_Comm_rank(MPI_COMM_WORLD, &my_PE_num);
//
//   if(my_PE_num == 0) {
//
//     MPI_Recv(&numbertoreceive, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//
//     std::cout << "NumberReceived = " << numbertoreceive << std::endl;
//   }
//   else {
//     std::cout << "I am PE = " << my_PE_num << std::endl;
//     MPI_Send(&numbertosend, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
//   }
//   MPI_Finalize();

  std::vector < std::vector < std::vector < double > > > T;
  unsigned N = 5;
  unsigned dim = 3;


  T.reserve(pow(4 * N,dim));
  for(unsigned i = 0; i < N; i++) {
    T[i].resize(dim - 1);
    for(unsigned k = 0; k < dim - 1; k++) {
      T[i][k].resize(dim, 0.);
    }
  }

  for(unsigned cnt = 0; cnt < N ; cnt++) {
    T[cnt][0][0] = 1.;
    T[cnt][0][1] = 2.;
    T[cnt][0][2] = 3.;

    T[cnt][1][0] = 4.;
    T[cnt][1][1] = 5.;
    T[cnt][1][2] = 6.;
    
  }

  std::cout << T[N-2][1][1] << " " << T.size() <<std::endl;


}



