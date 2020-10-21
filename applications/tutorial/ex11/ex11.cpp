#include <iostream>
#include <math.h>
#include <stdio.h>
#include "mpi.h"

// MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
// buf: starting adress of the buffer
// count: number of elements in that buffer
// dest: ID of the receiver
// tag: tag the messages
// MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,MPI_Comm comm, MPI_Status * status)


void SimpleMPI1();
void SimpleSR();
void SRMultiply();
int main(int argc, char** argv) {


  MPI_Init(&argc, &argv);

  SRMultiply();

  MPI_Finalize();


}



void SimpleMPI1() {
  int iproc;
  int numprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  std::cout << " This process # " << iproc << " / " << numprocs << std::endl;

}


// PE0 sends pi to others and they print once they get it.
void SimpleSR() {
  int iproc;
  int numprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  if(iproc == 0) {
    double pi = 3.14;
    for(int index = 1; index < numprocs ; index++) {
      MPI_Send(&pi, 1, MPI_DOUBLE, index, 0,  MPI_COMM_WORLD);
      std::cout << " Sent to : " << index << std::endl;
    }
  }
  else {
    // carefully understand processes are waiting iproc=0 to send the info
    double value;
    MPI_Status status;
    int ierr = MPI_Recv(&value, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    std::cout << " Received by : " << iproc << std::endl;
  }
}

// PE0 sends 10 and they multiply it by their ID and send back to 0 where
// its sum is calcuated and printed. Faster with broadcast

void SRMultiply() {
  
  int iproc;
  int numprocs;
  int result = 0;
  
  MPI_Status status;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  if(iproc == 0) {
    int number = 10;
    for(int index = 1; index < numprocs ; index++) {
      MPI_Send(&number, 1, MPI_DOUBLE, index, 0,  MPI_COMM_WORLD);
      std::cout << " Sent to : " << index << std::endl;
    }
  }
  else {
    int value;
    MPI_Recv(&value, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    result = value * iproc;
    std::cout << "Result on ID : " << iproc << " =  " << result << std::endl;
  }
  
  

  if(iproc == 0) {
    for(int index = 1; index < numprocs ; index++) {
      int value;
      MPI_Recv(&value, 1, MPI_DOUBLE, index, 0, MPI_COMM_WORLD, &status);
      result += value;
    }
    std::cout << " Total =  " << result << std::endl;
  }
  else {
    MPI_Send(&result, 1, MPI_DOUBLE, 0, 0,  MPI_COMM_WORLD);
  }


}



void DotProd(){
    std::vector<int> u = {2,4,5}; 
    std::vector<int> v = {4,3,6};
}

