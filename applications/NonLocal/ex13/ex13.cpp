#include <iostream>
#include <vector>
#include <omp.h>
#include <cmath>
#include "FemusInit.hpp"

using namespace femus;

int main(int argc, char** argv)  {

    FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

    const int n = 10000000;
    //std::vector<double> array(n);

    double array[n];

    // Initialize the array
    for (int i = 0; i < n; i++) {
        array[i] = 0;
    }
    double sum;
    // Parallelized for loop using OpenMP
    #pragma omp parallel for reduction(+:sum)
    //{
    //#pragma omp loop
    for (int i = 0; i < n; i++) {
        array[i] = sqrt(i) + 1.;
        sum += array[i];
        //Print thread information (optional)
        //std::cout << "Thread " << omp_get_thread_num() << " processed index " << i << std::endl;
    }
   // }

    // Print the array
    // std::cout << "Array contents: ";
    // double sum = 0.;
    // for (int i = 0; i < n; i++) {
    //   sum += array[i];
    // }
    std::cout << sum << std::endl;


    return 0;
}
