
#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

#include "slepceps.h"

unsigned lmax1 = 1; // consistency form 3 -> 7
const bool correctConstant = false;

#include "../ex20/include/nonlocal_assembly_adaptive.hpp"


void BuildI2(MultiLevelSolution & mlSol, const std::vector<double> &xg, const double & r);

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {

  bool dirichlet = true;

  value = 0.;

  for(unsigned k = 0; k < x.size(); k++) {
    // value +=  x[k] * x[k]; //consistency
    value +=  x[k] * x[k] * x[k]; //cubic
//   value +=  x[k] * x[k] * x[k] * x[k];//quartic
  }

  return dirichlet;
}

using namespace femus;

unsigned numberOfUniformLevels = 3; //consistency
//unsigned numberOfUniformLevels = 1; //cubic-quartic 2->6 //cubic Marta4Quad Tri Mix
//unsigned numberOfUniformLevels = 2; //cubic-quartic 2->4 mappa a 4->6 //cubic Marta4Fine


int main(int argc, char** argv) {

  if(argc == 3) {
    numberOfUniformLevels = atoi(argv[1]);

    lmax1 = atoi(argv[2]);
  }
  else if(argc == 2) {
    numberOfUniformLevels = atoi(argv[1]);
  }

  std::cout << "USING VARIABLES " << numberOfUniformLevels << "  " << lmax1 << std::endl;

  clock_t total_time = clock();

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;

  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;


  char fileName[100] = "../input/martaTest4.neu"; // good form 2->6 in serial but in parallel use martaTest4Fine
//   char fileName[100] = "../input/martaTest4Fine.neu"; // works till 144 nprocs +2
//   char fileName[100] = "../input/martaTest4Finer.neu"; // works till 144 nprocs +4
//   char fileName[100] = "../input/martaTest4Tri.neu";
// char fileName[100] = "../input/martaTest4Unstr.neu"; // works till 144 nprocs
  // char fileName[100] = "../input/salome/martaTest4QuadUnstr.med";
  // char fileName[100] = "../input/salome/martaTest4QuadUnstr2.med";
  // char fileName[100] = "../input/martaTest4-3D-tet.neu"; // works till 288 nprocs 0.2
  //char fileName[100] = "../input/martaTest4-3D.neu"; // works till 288 nprocs 0.2
  //char fileName[100] = "../input/martaTest4-3Dfine.neu"; // works till 576 and more nprocs +1 0.1

  mlMsh.ReadCoarseMesh(fileName, "fifth", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  unsigned dim = mlMsh.GetDimension();

  MultiLevelSolution mlSol(&mlMsh);


  // add variables to mlSol
  FEOrder femType = SERENDIPITY;
//   FEOrder femType = FIRST;

  std::vector < std::string > femTypeName = {"zero", "linear", "quadratic", "biquadratic"};

  mlSol.AddSolution("u", LAGRANGE,  femType, 0);
  mlSol.Initialize("All");


  std::vector<double> xg = {-0.40563508326896, -0.39436491673104};
  double r = 0.2;
  BuildI2(mlSol, xg, r);

  // ******* Print solution *******
  mlSol.SetWriter(VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back("All");
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, femTypeName[femType].c_str(), print_vars, 0);

  std::cout << std::endl << " total CPU time : " << std::setw(11) << std::setprecision(6) << std::fixed
            << static_cast<double>((clock() - total_time)) / CLOCKS_PER_SEC << " s" << std::endl;

  return 0;

} //end main

void BuildI2(MultiLevelSolution & mlSol, const std::vector<double> &xg2, const double & delta) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  Solution* sol  = mlSol.GetSolutionLevel(level);

  const unsigned  dim = msh->GetDimension();

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

  RefineElement *refineElement[6][3];

  unsigned soluIndex = mlSol.GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol.GetSolutionType(soluIndex);


  NonLocal *nonlocal;

  double eps = 0.;


  BallApproximation *ballAprx = new BallApproximation();

  if(dim == 3) {
    refineElement[0][0] = new RefineElement(lmax1, "hex", "linear", "third", "third", "legendre");
    refineElement[0][1] = new RefineElement(lmax1, "hex", "quadratic", "third", "third", "legendre");
    refineElement[0][2] = new RefineElement(lmax1, "hex", "biquadratic", "third", "third", "legendre");

    refineElement[1][0] = new RefineElement(lmax1, "tet", "linear", "third", "third", "legendre");
    refineElement[1][1] = new RefineElement(lmax1, "tet", "quadratic", "third", "third", "legendre");
    refineElement[1][2] = new RefineElement(lmax1, "tet", "biquadratic", "third", "third", "legendre");

    refineElement[0][soluType]->SetConstants(eps);
    refineElement[1][soluType]->SetConstants(eps);

    nonlocal = new NonLocalBall3D();

  }
  else if(dim == 2) {
    refineElement[3][0] = new RefineElement(lmax1, "quad", "linear", "fourth", "fourth", "legendre");
    refineElement[3][1] = new RefineElement(lmax1, "quad", "quadratic", "fourth", "fourth", "legendre");
    refineElement[3][2] = new RefineElement(lmax1, "quad", "biquadratic", "fourth", "fourth", "legendre");

    refineElement[4][0] = new RefineElement(lmax1, "tri", "linear", "fourth", "fourth", "legendre");
    refineElement[4][1] = new RefineElement(lmax1, "tri", "quadratic", "fourth", "fourth", "legendre");
    refineElement[4][2] = new RefineElement(lmax1, "tri", "biquadratic", "fourth", "fourth", "legendre");

    refineElement[3][soluType]->SetConstants(eps);
    refineElement[4][soluType]->SetConstants(eps);

    nonlocal = new NonLocalBall();
    //nonlocal = new NonLocalBox();
  }

  double localAreap = 0.;
  double localI2p = 0.;

  double localAreal = 0.;
  double localI2l = 0.;

  unsigned xType = 2.;

  std::vector < std::vector <double> >  xv1(dim);
  std::vector < std::vector <double> >  xv1l;


    double subarea_total = 0.0;
    double subI2_total = 0.0;

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    short unsigned ielGeom = msh->GetElementType(iel);

    const OctTreeElement &octTreeElement1 = refineElement[ielGeom][soluType]->GetOctTreeElement1();
    const OctTreeElement &octTreeElement1CF = refineElement[ielGeom][soluType]->GetOctTreeElement1CF();
    RefineElement &element1 = *refineElement[ielGeom][soluType];

    const unsigned &nDof1 = element1.GetNumberOfNodes();

    for(int k = 0; k < dim; k++) {
      xv1[k].resize(nDof1);
    }
    for(unsigned i = 0; i < nDof1; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        xv1[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    const elem_type *fem1 = element1.GetFem1();
    const elem_type *fem1CF = element1.GetFem1CF();

    const unsigned &ng1 = fem1->GetGaussPointNumber();
    std::vector<std::vector<double>> xg1(ng1, std::vector<double>(dim, 0));
    std::vector<double> weight1(ng1);
    for(unsigned ig = 0; ig < ng1; ig++) {
      const double *phi;
      fem1->GetGaussQuantities(xv1, ig, weight1[ig], phi);
      for(unsigned i = 0; i < nDof1; i++) {
        for(unsigned k = 0; k < dim; k++) {
          xg1[ig][k] += phi[i] * xv1[k][i];
        }
      }
    }


    const unsigned &ng1CF = fem1CF->GetGaussPointNumber();
    std::vector<std::vector<double>> xg1CF(ng1CF, std::vector<double>(dim, 0));
    std::vector<double>weight1CF(ng1CF);
    for(unsigned ig = 0; ig < ng1CF; ig++) {
      const double *phi;
      fem1CF->GetGaussQuantities(xv1, ig, weight1CF[ig], phi);
      for(unsigned i = 0; i < nDof1; i++) {
        for(unsigned k = 0; k < dim; k++) {
          xg1CF[ig][k] += phi[i] * xv1[k][i];
        }
      }
    }


    //BEGIN NEW STUFF

    std::vector < std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> > x1MinMax(dim);
    for(unsigned k = 0; k < dim; k++) {
      x1MinMax[k] = std::minmax_element(xv1[k].begin(), xv1[k].end());
    }

    bool coarseIntersectionTest = true;
    for(unsigned k = 0; k < dim; k++) {
      if((xg2[k]  - * (x1MinMax[k].second)) > delta || (*(x1MinMax[k].first) - xg2[k]) > delta) {  // this can be improved with the l2 norm
        coarseIntersectionTest = false;
        break;
      }
    }

    if(coarseIntersectionTest) {

      xv1l = xv1;
      for(unsigned k = 0; k < xv1l.size(); k++) xv1l[k].resize(element1.GetNumberOfLinearNodes());

      std::vector<double> a;
      double d;
      unsigned cut;
      ballAprx->GetNormal(element1.GetElementType(), xv1l, xg2, delta, a, d, cut);




      if(cut == 0) { //interior element
        double d2W1 = 0.;
        for(unsigned ig = 0; ig < ng1; ig++) {
          double d2 = 0.;
          for(unsigned k = 0; k < dim; k++) {
            d2 += (xg2[k] - xg1[ig][k]) * (xg2[k] - xg1[ig][k]);
          }
          localI2l += d2 * weight1[ig];
          localAreal += weight1[ig];

          d2W1 += d2 * weight1[ig];

          localI2p += d2 * weight1[ig];
          localAreap += weight1[ig];
        }
        //std::cout << iel << " " << d2W1 << std::endl;
      }
      else if(cut == 1) { //cut element
        element1.GetCutFem()->clear();

        std::vector<double> weightsTMP;
        element1.GetCDweight()->GetWeight(a, d, weightsTMP);

        // BEGIN Parabola integration
        bool twoInt = false;



        std::vector<double> A(6, 0.);
        A[0] = -1;
        A[1] = 0;
        A[2] = -1;
        A[3] = + 2 * xg2[0];
        A[4] = + 2 * xg2[1];
        A[5] = - xg2[0] * xg2[0] - xg2[1] * xg2[1] + delta * delta;

        std::vector<double> eqPolyWeight;
        element1.GetCDWeightPar()->GetWeight(xv1l, A, eqPolyWeight, twoInt);
        if(!twoInt) {
          std::cout << "not twoIntersections!\n";
          element1.GetCDweight()->GetWeight(a, d, eqPolyWeight);
        }

        double d2W1 = 0.;
//         cout<< " ng1cf = " << ng1CF << endl;
        for(unsigned ig = 0; ig < ng1CF; ig++) {
          double d2 = 0.;
          for(unsigned k = 0; k < dim; k++) {
            d2 += (xg2[k] - xg1CF[ig][k]) * (xg2[k] - xg1CF[ig][k]);
          }

          localI2l += d2 * weight1CF[ig] * weightsTMP[ig];
          localAreal += weight1CF[ig] * weightsTMP[ig];

          d2W1 += d2 * weight1CF[ig] * eqPolyWeight[ig];

          localI2p += d2 * weight1CF[ig] * eqPolyWeight[ig];
          localAreap += weight1CF[ig] * eqPolyWeight[ig];
        }

        std::cout << iel << " " << d2W1 << " , "<< std::endl;
//         std::cout << "     center" << xg2[0]<<" " <<xg2[1] << " delta "<<delta << std::endl;
/*
        for(unsigned i = 0; i < eqPolyWeight.size(); i++) {
            std::cout << eqPolyWeight[i] << " , ";
        }
        std::cout<<std::endl;*/

//         std::cout<< " weight1 size = " << weight1CF.size()<<std::endl;
//         for(unsigned i = 0; i < eqPolyWeight.size(); i++) {
//             std::cout << weight1CF[i] << " , ";
//         }
//         std::cout<<std::endl;

        int n = 10 ;

        // Loop through each vertex to find the min and max of x and y
        double min_x = xv1l[0][0];
        double max_x = xv1l[0][0];
        double min_y = xv1l[1][0];
        double max_y = xv1l[1][0];

        // Loop through each vertex to find the min and max of x and y
        for (int i = 1; i < 4; ++i) {  // Start from the second vertex (index 1)
            min_x = std::min(min_x, xv1l[0][i]);
            max_x = std::max(max_x, xv1l[0][i]);
            min_y = std::min(min_y, xv1l[1][i]);
            max_y = std::max(max_y, xv1l[1][i]);
        }

        // Here hx and hy should be equal. If not print a warning
        double hx = (max_x - min_x) / n;
        double hy = (max_y - min_y) / n;
        cout << " h = " <<hx << " "<< hy<<endl;
        if (fabs(hx-hy) > 0.0000000001)std::cout<<":::::WARNING::::: hx and hy are not equal "<< endl;
        double subarea = 0.0;
        double subI2 = 0.0;

        // Loop over each sub-element
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                // Calculate the middle point of the sub-element
                double xm = min_x + hx * (i + 0.5);
                double ym = min_y + hx * (j + 0.5);

                // Check the condition (xm - xc)^2 + (ym - yc)^2 - delta^2
                double check_value = (xm - xg2[0])*(xm - xg2[0]) + (ym - xg2[1])*(ym - xg2[1]) - delta*delta;
                if (check_value < 0) {
                    // Add to the area if the sign is negative
                    subarea += hx * hx;
                    subI2 += hx * hx * (xm * xm + ym * ym);
                }
            }
        }

        // Output the results
        std::cout << "sub Area: " << subarea << std::endl;
        std::cout << "sub I2: " << subI2 << std::endl;

        subarea_total += subarea;
        subI2_total += subI2;
      }
    }
  }

  std::cout << "sub Area total : " << subarea_total << std::endl;
  std::cout << "sub I2 total : " << subI2_total << std::endl;

  double areal = 0.;
  double areap = 0.;
  double area = M_PI * delta * delta;
  std::cout << "Area Analytic = " << area << std::endl;

  MPI_Allreduce(&localAreal, &areal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  std::cout << "Areal = " << areal << " " << area/areal << std::endl;

  MPI_Allreduce(&localAreap, &areap, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  std::cout << "Areap = " << areap <<" " << area/areap << std::endl;

  double I2l = 0.;
  double I2p = 0.;
  double I2 = 0.5 * M_PI * pow(delta, 4);
  std::cout << "I2 Analytic = " << I2 << std::endl;

  MPI_Allreduce(&localI2l, &I2l, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  std::cout << "I2l = " << I2l <<" "<< I2/I2l << std::endl;

  MPI_Allreduce(&localI2p, &I2p, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  std::cout << "I2p = " << I2p <<" "<< I2/I2p << std::endl;




  delete ballAprx;

  if(dim == 3) {
    delete refineElement[0][0];
    delete refineElement[0][1];
    delete refineElement[0][2];
    delete refineElement[1][0];
    delete refineElement[1][1];
    delete refineElement[1][2];
  }
  else if(dim == 2) {
    delete refineElement[3][0];
    delete refineElement[3][1];
    delete refineElement[3][2];
    delete refineElement[4][0];
    delete refineElement[4][1];
    delete refineElement[4][2];
  }


}
