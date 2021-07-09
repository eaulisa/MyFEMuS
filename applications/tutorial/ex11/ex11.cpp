#include <iostream>
#include <math.h>
#include <stdio.h>
#include "mpi.h"
#include <vector>
/*
MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)

buf: starting adress of the buffer
count: number of elements in that buffer
tag: tag the messages

MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,MPI_Comm comm, MPI_Status * status)
MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)

*/

// void SimpleMPI1();
// void SimpleSR();
// void SRMultiply();
// double f(double x);
// double Trap(double a, double b, int n);
// void ParallelTrap();
void LaplaceSol();
int main(int argc, char** argv) {

  MPI_Init(&argc, &argv);
  MPI_Status status;
  LaplaceSol();

  MPI_Finalize();


}




void LaplaceSol() {
    
  unsigned rows = 100;
  unsigned cols = 100;
  std::vector<std::vector<double>> T(rows);
  std::vector<std::vector<double>> Told(rows);

  for(unsigned i = 0; i < rows; i++) {
    T[i].resize(cols, 0.);
    Told[i].resize(cols, 0.);
  }

  //boundary values
  for(unsigned i = 0; i < rows ; i++) {
    Told[i][0] = 0.;
    Told[i][cols - 1] = 100. * i / rows;
  }

  for(unsigned j = 0; j < cols ; j++) {
    Told[0][j] = 0.;
    Told[rows - 1][j] = 100. * j / cols;
  }


  //main solver starts here
  unsigned iter = 1;
  unsigned MaxIter = 3000;
  double MaxTempError = 0.01;
  double dt = 50.;

  std::cout << "----- Iteration number---- " << " --------Temprature----" << std::endl;

  while(dt > MaxTempError && iter < MaxIter) {

    for(unsigned i = 1; i < rows - 1; i++) {
      for(unsigned j = 1; j < cols - 1; j++) {
        T[i][j] = 0.25 * (Told[i - 1][j] + Told[i + 1][j] + Told[i][j - 1] + Told[i][j + 1]);
      }
    }
    dt = 0.;

    for(unsigned i = 1; i < rows - 1; i++) {
      for(unsigned j = 1; j < cols - 1; j++) {
        dt = fmax(fabs(T[i][j] - Told[i][j]), dt);
        Told[i][j] = T[i][j];
      }
    }

    if( (iter % 100) == 0) {
      for(unsigned i = rows - 5; i < rows-1; i++) {
        std::cout << "-----" << iter << "-------- " << T[i][i] << std::endl;
      }

    }
    iter++;

  }


}



// Each process makes its own traprule and we add them up to get the final sum
/*void ParallelTrap(){
    double a = 0.;
    double b = 1.;
    int nTotal = 100;
    double h = (b - a) / nTotal;

    int iproc;
    int numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    int n_local = nTotal / numprocs;

    MPI_Bcast(&n_local, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double a_loc = a + iproc * n_local * h;
    double b_loc = a + (iproc + 1) * n_local * h;
    double local_sum = Trap(a_loc, b_loc, n_local);
    std::cout << "iproc = " << iproc << " local_sum = " << local_sum << std::endl;

    double result;
    MPI_Reduce(&local_sum, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(iproc==0){
        std::cout << "Result is = " << result << std::endl;
    }

}

double Trap(double a, double b, int n){

    double h = (b - a)/n;
    double sum = h / 2. * (f(a) + f(b));
    for(int i = 1; i <= n - 1 ; i++){
        sum += h * f(a + h * i);
    }
    return sum;

}


double f(double x){
    return x*x;
}
*/

// void SimpleMPI1() {
//   int iproc;
//   int numprocs;
//   MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
//   MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//   std::cout << " This process # " << iproc << " / " << numprocs << std::endl;
//
// }

// // PE0 sends pi to others and they print once they get it.
// void SimpleSR() {
//   int iproc;
//   int numprocs;
//   MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
//   MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//
//   if(iproc == 0) {
//     double pi = 3.14;
//     for(int index = 1; index < numprocs ; index++) {
//       MPI_Send(&pi, 1, MPI_DOUBLE, index, 0,  MPI_COMM_WORLD);
//       std::cout << " Sent to : " << index << std::endl;
//     }
//   }
//   else {
//     // carefully understand processes are waiting iproc=0 to send the info
//     double value;
//     MPI_Status status;
//     int ierr = MPI_Recv(&value, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
//     std::cout << " Received by : " << iproc << std::endl;
//   }
// }

// PE0 sends 10 and they multiply it by their ID and send back to 0 where
// its sum is calcuated and printed. Faster with broadcast

// void SRMultiply() {
//
//   int iproc;
//   int numprocs;
//   int result = 0;
//
//   MPI_Status status;
//   MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
//   MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//
//   if(iproc == 0) {
//     int number = 10;
//     for(int index = 1; index < numprocs ; index++) {
//       MPI_Send(&number, 1, MPI_DOUBLE, index, 0,  MPI_COMM_WORLD);
//       std::cout << " Sent to : " << index << std::endl;
//     }
//   }
//   else {
//     int value;
//     MPI_Recv(&value, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
//     result = value * iproc;
//     std::cout << "Result on ID : " << iproc << " =  " << result << std::endl;
//   }
//
//
//
//   if(iproc == 0) {
//     for(int index = 1; index < numprocs ; index++) {
//       int value;
//       MPI_Recv(&value, 1, MPI_DOUBLE, index, 0, MPI_COMM_WORLD, &status);
//       result += value;
//     }
//     std::cout << " Total =  " << result << std::endl;
//   }
//   else {
//     MPI_Send(&result, 1, MPI_DOUBLE, 0, 0,  MPI_COMM_WORLD);
//   }
//
//
// }


