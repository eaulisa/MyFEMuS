

#include "petsc.h"


int main(int argc, char** argv) {

 
  int ierr = PetscInitialize (&argc, &argv, NULL, NULL);    
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = PetscFinalize();
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  return 0;

} //end main
