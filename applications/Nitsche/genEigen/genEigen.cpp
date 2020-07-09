
#include "slepceps.h"
#include <iostream>

int main(int argc, char** args) {

  SlepcInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);

  double a[6][6] = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 1.9247225877689356, -1.9247225877689356, -0.96236129388446778, 0.0, 0.96236129388446778 },
    {0.0, -1.9247225877689356, 1.9247225877689356, 0.96236129388446778, 0.0, -0.96236129388446778},
    {0.0, -0.96236129388446778, 0.96236129388446778, 0.48118064694223389, 0.0, -0.48118064694223389},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.96236129388446778, -0.96236129388446778, -0.48118064694223389, 0.0, 0.48118064694223389}
  };


  double b[6][6] = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.046728852764308028, -0.046728852764308028, -0.023364426382154014, 0.0, 0.023364426382154014},
    {0.0, -0.046728852764308028, 0.046728852764308028, 0.023364426382154014, 0.0, -0.023364426382154014 },
    {0.0, -0.023364426382154014, 0.023364426382154014, 0.011682213191077007, 0.0, -0.011682213191077007 },
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.023364426382154014, -0.023364426382154014, -0.011682213191077007, 0.0, 0.011682213191077007 }
  };

  Mat A, B;
  EPS eps;
  
  int idx[6] = {0, 1, 2, 3, 4, 5};
 
  MatCreateAIJ(PETSC_COMM_SELF,6,6,6,6,6,PETSC_NULL,0,PETSC_NULL,&A);
  MatCreateAIJ(PETSC_COMM_SELF,6,6,6,6,6,PETSC_NULL,0,PETSC_NULL,&B);
 
  MatSetValues(A, 6, &idx[0], 6, &idx[0], &a[0][0], INSERT_VALUES);
  MatSetValues(B, 6, &idx[0], 6, &idx[0], &b[0][0], INSERT_VALUES);
    
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

  double nrm;
  MatNorm(B, NORM_INFINITY, &nrm);
  //MatShift(B, 1.0e-10 * nrm);

  EPSCreate(PETSC_COMM_SELF, &eps);
  EPSSetOperators(eps, A, B);

  EPSSetType(eps, EPSLAPACK);
  
  EPSSetWhichEigenpairs(eps, EPS_LARGEST_REAL);
  EPSSetPurify(eps, PETSC_TRUE);

  EPSSetFromOptions(eps);

  EPSSolve(eps);

  double real;
  EPSGetEigenpair(eps, 0, &real, PETSC_NULL, PETSC_NULL, PETSC_NULL);
  std::cout << real << " " << std::endl;

  EPSDestroy(&eps);
  MatDestroy(&A);
  MatDestroy(&B);

  return 1;
}
