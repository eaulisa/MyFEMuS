
#include "FemusInit.hpp"
#include "adept.h"


using namespace femus;


int main(int argc, char** args) {

  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
  std::cout<<"Hello World!"<<std::endl;

  return 0;
}


