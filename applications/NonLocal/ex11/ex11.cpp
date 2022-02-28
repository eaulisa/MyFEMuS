
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "CurrentElem.hpp"
#include "LinearImplicitSystem.hpp"

#include "PolynomialBases.hpp"

#include "CutFemIntegration.hpp"

#include <vector>
#include <cmath>
#include <iostream>

using namespace std;
using namespace femus;

#define N_UNIFORM_LEVELS  8
#define N_ERASED_LEVELS   0

#define EX_1       -1.
#define EX_2        1.
#define EY_1       -1.
#define EY_2        1.
#define N_X         4
#define N_Y         4

double InitialValueU(const std::vector < double >& x) {
  return 0. * x[0] * x[0];
}
bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;

//   if(facename == 1) {
//     dirichlet = true; //dirichlet
//     value = 0.;
//   }
//   else if(facename == 2) {
//     dirichlet = true; //dirichlet
//     value = 0.;
//   }

  return dirichlet;
}

void GetNormalQuad(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &a, double &d,  std::vector<double> &xm, std::vector<double> &b, double &db, unsigned &cut);

void GetNormalTri(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &a, double &d,  std::vector<double> &xm, std::vector<double> &b, double &db, unsigned &cut);


void SimpleNonlocalAssembly(MultiLevelProblem& ml_prob);


const elem_type *finiteElementQuad;

int main(int argc, char** argv) {

  typedef double TypeIO;
  typedef cpp_bin_float_oct TypeA;



  const std::string fe_quad_rule_1 = "seventh";
  const std::string fe_quad_rule_2 = "eighth";

  // ======= Init ========================
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  finiteElementQuad = new const elem_type_2D("quad", "linear", "fifth", "legendre");

  unsigned numberOfUniformLevels = N_UNIFORM_LEVELS;

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.GenerateCoarseBoxMesh(N_X, N_Y, 0, EX_1, EX_2, EY_1, EY_2, 0., 0., TRI6, fe_quad_rule_1.c_str());
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);

  /*

  char fileName[100] = "../input/martaTest4Unstr.neu"; // works till 144 nprocs
  mlMsh.ReadCoarseMesh(fileName, "fifth", scalingFactor);
  MPI_Barrier(MPI_COMM_WORLD);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);*/






  // erase all the coarse mesh levels
  const unsigned erased_levels = N_ERASED_LEVELS;
  mlMsh.EraseCoarseLevels(erased_levels);

  const unsigned level = N_UNIFORM_LEVELS - N_ERASED_LEVELS - 1;

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("u", LAGRANGE, SECOND, 2);


  mlSol.Initialize("All");
  mlSol.Initialize("u", InitialValueU);

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("All");

  MultiLevelProblem ml_prob(&mlSol);

  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("FracProblem");


  Mesh*                    msh = mlMsh.GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)

  unsigned xType = 2;

  const unsigned  dim = msh->GetDimension();

  std::vector < std::vector < double > > x1;

  FILE * fp;

  fp = fopen("lines.dat", "w");

  unsigned qM = 3;
  CutFemIntegral <TypeIO, TypeA> quad  = CutFemIntegral<TypeIO, TypeA >(QUAD, qM, "legendre");
  CutFemIntegral <TypeIO, TypeA> tri  = CutFemIntegral<TypeIO, TypeA >(TRI, qM, "legendre");


  std::vector<double> xg(2, 0);
  xg[0] -= 0.054;
  xg[1] -= 0.012;
  double R = 0.5;
  double CircArea = 0.;

  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned ielType = msh->GetElementType(iel);
    unsigned nDof = msh->GetElementDofNumber(iel, 0);  // number of coordinate linear element dofs
    x1.resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      x1[k].resize(nDof);
    }

    for(unsigned k = 0; k < dim; k++) {
      for(unsigned i = 0; i < nDof; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
        x1[k][i] = (*msh->_topology->_Sol[k])(xDof);  // global extraction and local storage for the element coordinates
      }
    }

    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> xm;
    double d;
    double db;
    unsigned cut;

    double h2;

    if(ielType == 3) {
      GetNormalQuad(x1, xg, R, a, d, xm, b, db, cut);
      h2 = ((EX_2 - EX_1) / (N_X * pow(2, N_UNIFORM_LEVELS - 1))) * ((EY_2 - EY_1) / (N_Y * pow(2, N_UNIFORM_LEVELS - 1)));
    }
    else {
      GetNormalTri(x1, xg, R, a, d, xm, b, db, cut);
      h2 = 0.5 * ((EX_2 - EX_1) / (N_X * pow(2, N_UNIFORM_LEVELS - 1))) * ((EY_2 - EY_1) / (N_Y * pow(2, N_UNIFORM_LEVELS - 1)));
    }
    if(cut == 1) {
      if(ielType == 3) {
        std::vector <TypeIO> weightCF;
        quad(qM, 0, b, db, weightCF);

        const double* weightG = quad.GetGaussWeightPointer();

        double sum = 0.;
        for(unsigned ig = 0; ig < weightCF.size(); ig++) {
          sum += weightG[ig] * weightCF[ig] ; //TODO use the correct quad rule!!!!!
        }
        CircArea += sum / 4. * h2;
      }
      else if(ielType == 4) {
        std::vector <TypeIO> weightCF;
        tri(qM, 0, b, db, weightCF);

        const double* weightG = tri.GetGaussWeightPointer();

        double sum = 0.;
        for(unsigned ig = 0; ig < weightCF.size(); ig++) {
          sum += weightG[ig] * weightCF[ig] ; //TODO use the correct quad rule!!!!!
        }
        CircArea += 2. * sum * h2;
      }
      
      
    }
    else if(cut == 0) {
      CircArea += h2; //TODO
    }

    /* trivial print for xmgrace */
    if(cut == 1) {
      double xx = xm[0] - 0.5 * a[1] * 0.1 * sqrt(h2);
      double yy = xm[1] + 0.5 * a[0] * 0.1 * sqrt(h2);
      fprintf(fp, "%f %f \n", xx, yy);
      xx = xm[0] + 0.5 * a[1];
      yy = xm[1] - 0.5 * a[0];
      fprintf(fp, "%f %f \n", xx, yy);
      fprintf(fp, "\n \n");
    }
  }
  std::cout.precision(14);
  std::cout << "AREA CIRCLE = " << CircArea << "  analytic value = " << M_PI * R * R << "\n";

  mlSol.SetWriter(VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back("All");
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "quadratic", print_vars, 0);


  /*Testing the function GetNormalQuad inside a nonlocal assembly-like function*/
  //SimpleNonlocalAssembly(ml_prob);

// //   /* Basic numerical tests for the circle - no mesh involved*/
//   std::vector < std::vector<double> > xva = {{1., 2., 4., 0.}, {0., 0., 4., 1.}};
//   std::vector < std::vector<double> > xvb = {{0., 4., 3., 2.}, {0., 0., 1., 1.}};
//   std::vector < std::vector<double> > xvc0 = {{0., 1., 2., 2.}, {0., 0., 1., 2.}};
//   std::vector < std::vector<double> > xvc1 = {{2., 0., 1., 2.}, {2., 0., 0., 1.}};
//   std::vector < std::vector<double> > xvc2 = {{2., 2., 0., 1.}, {1., 2., 0., 0.}};
//   std::vector < std::vector<double> > xvc3 = {{1., 2., 2., 0.}, {0., 1., 2., 0.}};
//   std::vector < std::vector<double> > xv0 = {{1., 2., 2., 1.}, {1., 1., 2., 2.}};
//   std::vector < std::vector<double> > xv1 = {{1., 1., 2., 2.}, {2., 1., 1., 2.}};
//   std::vector < std::vector<double> > xv2 = {{2., 1., 1., 2.}, {2., 2., 1., 1.}};
//   std::vector < std::vector<double> > xv3 = {{2., 2., 1., 1.}, {1., 2., 2., 1.}};
//   std::vector < std::vector<double> > xvr0 = {{1., 2., 3., 2.}, {1., 2., 4., 3.}};
//   std::vector < std::vector<double> > xvr1 = {{2., 1., 2., 3.}, {3., 1., 2., 4.}};
//   std::vector < std::vector<double> > xvr2 = {{3., 2., 1., 2.}, {4., 3., 1., 2.}};
//   std::vector < std::vector<double> > xvr3 = {{2., 3., 2., 1.}, {2., 4., 3., 1.}};
//   std::vector<double> xg(2, 0);
//   double R = 2.;
//   std::vector<double> a;
//   double d;
//   std::vector<double> xm;
//   std::vector<double> b;
//   double db;
//   unsigned  cut;
//
//   GetNormalQuad(xva, xg, R, a, d, xm, b, db, cut);
//   std::cout << std::endl;
//   GetNormalQuad(xvb, xg, R, a, d, xm, b, db, cut);
//   std::cout << std::endl;
//   GetNormalQuad(xvc0, xg, R, a, d, xm, b, db, cut);
//   std::cout << std::endl;
//   GetNormalQuad(xvc1, xg, R, a, d, xm, b, db, cut);
//   std::cout << std::endl;
//   GetNormalQuad(xvc2, xg, R, a, d, xm, b, db, cut);
//   std::cout << std::endl;
//   GetNormalQuad(xvc3, xg, R, a, d, xm, b, db, cut);
//   std::cout << std::endl;
//   GetNormalQuad(xv0, xg, R, a, d, xm, b, db, cut);
//   std::cout << std::endl;
//   GetNormalQuad(xv1, xg, R, a, d, xm, b, db, cut);
//   std::cout << std::endl;
//   GetNormalQuad(xv2, xg, R, a, d, xm, b, db, cut);
//   std::cout << std::endl;
//   GetNormalQuad(xv3, xg, R, a, d, xm, b, db, cut);
//   std::cout << std::endl;
//   GetNormalQuad(xvr0, xg, R, a, d, xm, b, db, cut);
//   std::cout << std::endl;
//   GetNormalQuad(xvr1, xg, R, a, d, xm, b, db, cut);
//   std::cout << std::endl;
//   GetNormalQuad(xvr2, xg, R, a, d, xm, b, db, cut);
//   std::cout << std::endl;
//   GetNormalQuad(xvr3, xg, R, a, d, xm, b, db, cut);

  fclose(fp);

  delete finiteElementQuad;
  return 1;
}


void GetNormalQuad(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &a, double &d,  std::vector<double> &xm, std::vector<double> &b, double &db, unsigned &cut) {

  const unsigned &dim =  xv.size();
  const unsigned &nve =  xv[0].size();

  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& x4 = xv[0][3];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];
  const double& y4 = xv[1][3];

  double hx = 0.5 * (fabs(x3 - x1) + fabs(x4 - x2));
  double hy = 0.5 * (fabs(y3 - y1) + fabs(y4 - y2));
  double h = sqrt(hx * hx + hy * hy);
  double eps = 1.0e-10 * h;

  std::vector<double> dist(nve, 0);
  std::vector<double> dist0(nve);
  unsigned cnt0 = 0;
  for(unsigned i = 0; i < nve; i++) {
    for(unsigned k = 0;  k < dim; k++) {
      dist[i] += (xv[k][i] - xg[k]) * (xv[k][i] - xg[k]);
    }
    dist[i] = sqrt(dist[i]) - R;

    if(fabs(dist[i]) < eps) {
      dist0[i] = (dist[i] < 0) ? -eps : eps;
      dist[i] = 0.;
      cnt0++;
    }
    else {
      dist0[i] = dist[i];
    }
  }

  if(cnt0 > 0) {
    unsigned cntp = 0;
    for(unsigned i = 0; i < nve; i++) {
      if(dist[i] > 0) cntp++;
      dist[i] = dist0[i];
    }
    if(cntp == 0) { // the element is inside the ball
      cut = 0;
      return;
    }
    else if(cntp == nve - cnt0) {  // the element in outside the ball
      cut = 2;
      return;
    }
  }

  std::vector <double> theta(2);
  unsigned cnt = 0;
  for(unsigned e = 0; e < nve; e++) {
    unsigned ep1 = (e + 1) % nve;
    if(dist[e] * dist[ep1] < 0) {
      double s = 0.5  * (1 + (dist[e] + dist[ep1]) / (dist[e] - dist[ep1]));
      theta[cnt] = atan2((1 - s) * xv[1][e] + s * xv[1][ep1]  - xg[1], (1 - s) * xv[0][e] + s * xv[0][ep1] - xg[0]) ;
      cnt++;
    }
  }

  if(cnt == 0) {
    if(dist[0] < 0) cut = 0; // cell inside the ball
    else cut = 2; // cell outside the ball
    return;
  }
  else {
    cut = 1;
    if(theta[0] > theta[1]) {
      std::swap(theta[0], theta[1]);
    }
    double DT = theta[1] - theta[0];
    if(DT > M_PI) {
      std::swap(theta[0], theta[1]);
      theta[1] += 2. * M_PI;
      DT = theta[1] - theta[0];
    }
    xm.resize(dim);

    d = R * sqrt(0.5 * DT / tan(0.5 * DT)) ;
    a.resize(dim);
    a[0] = -cos(theta[0] + 0.5 * DT);
    a[1] = -sin(theta[0] + 0.5 * DT);

    for(unsigned k = 0; k < dim; k++) {
      xm[k] = -a[k] * d + xg[k];
    }
    d += - a[0] * xg[0] - a[1] * xg[1]; //TODO

    std::cout.precision(14);

//     std::cout << "xm = " << xm[0] << " " << xm[1] << std::endl;
//     std::cout << "a = " << a[0] << " b = " << a[1] << " d = " << d << std::endl;

    std::vector<double> xi(dim);
    double &u = xi[0];
    double &v = xi[1];

    std::vector < std::vector < double > > J(2, std::vector<double>(2));

    double dx12 = x1 - x2;
    double dx34 = x3 - x4;
    double dy12 = y1 - y2;
    double dy34 = y3 - y4;
    double hu = dx34 * dy12 - dx12 * dy34;

    double dx14 = (x1 - x4);
    double dy23 = (y2 - y3);
    double dx23 = (x2 - x3);
    double dy14 = (y1 - y4);
    double hv = dx14 * dy23 - dx23 * dy14;

    double eps2 = 1.0e-10 * h * h;

    if(fabs(hu) > eps2) {//edges 1 and 3 are not parallel
      double gu = -x4 * y1 + x3 * y2 - x2 * y3 + x1 * y4;
      double f = xm[0] * (dy12 + dy34) - xm[1] * (dx12 + dx34);
      double fpgu = f + gu;

      double det = sqrt(hu * (- 2. * xm[0] * (dy14 + dy23)
                              + (2. * xm[1] - y3 - y4) * (x1 + x2)
                              - (2. * xm[1] - y1 - y2) * (x3 + x4))
                        + fpgu * fpgu);
      u = (fpgu + det) / hu;

      if(fabs(hv) > eps2) { //edges 2 and 4 are not parallel
        double gv = -x4 * y3 + x3 * y4 - x2 * y1 + x1 * y2;
        v = (f + gv - det) / hv;

        J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
        J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
        J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);
        J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);

      }
      else { //edges 2 and 4 are parallel
        //   std::cout << "2 and 4 are parallel\n";
        J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
        J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);

        v = (J[0][1] > eps) ?
            (0.25 * ((-1. + u) * (x1 + x4) - (1. + u) * (x3 + x2)) + xm[0]) / J[0][1] :
            (0.25 * ((-1. + u) * (y1 + y4) - (1. + u) * (y3 + y2)) + xm[1]) / J[1][1];

        J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
        J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);

      }
    }
    else if(fabs(hv) > eps2) {  //edges 1 and 3 are parallel, but edges 2 and 4 are not
      // std::cout << "1 and 3 are parallel\n";
      double f = xm[0] * (dy12 + dy34) - xm[1] * (dx12 + dx34);
      double gv = -x4 * y3 + x3 * y4 - x2 * y1 + x1 * y2;
      double fpgv = f + gv;

      double det = sqrt(hv * (- 2. * xm[0] * (dy12 - dy34)
                              + (2. * xm[1] - y2 - y3) * (x1 + x4)
                              - (2. * xm[1] - y1 - y4) * (x2 + x3))
                        +  fpgv * fpgv);

      v = (fpgv - det) / hv;

      J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
      J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);

      u = (fabs(J[0][0]) > eps) ?
          (0.25 * ((-1. + v) * (x1 + x2) - (1. + v) * (x3 + x4)) + xm[0]) / J[0][0] :
          (0.25 * ((-1. + v) * (y1 + y2) - (1. + v) * (y3 + y4)) + xm[1]) / J[1][0];

      J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
      J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);
    }
    else { //edges 1 and 3, and  edges 2 and 4 are parallel
      //   std::cout << "Romboid\n";
      std::vector<std::vector<unsigned> > idx = {{3, 1}, {0, 2}};

      double A[2][2] = {{-dy14, dy23}, {dy12, dy34}};
      double B[2][2] = {{dx14, -dx23}, {-dx12, -dx34}};

      for(unsigned k = 0; k < 2; k++) {
        double d[2];
        for(unsigned j = 0 ; j < 2; j++) {
          double Ckj = - A[k][j] * xv[0][idx[k][j]] - B[k][j] * xv[1][idx[k][j]];
          d[j] = (A[k][j] * xm[0] + B[k][j] * xm[1] + Ckj) / sqrt(A[k][j] * A[k][j] + B[k][j] * B[k][j]);
        }
        xi[k] = -1. + 2. * d[0] / (d[0] + d[1]);
      }

      J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
      J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
      J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);
      J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);
    }

    double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    std::vector < std::vector < double > > Ji = {{J[1][1] / det, -J[0][1] / det}, {-J[1][0] / det, J[0][0] / det}};

    b.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        b[k] += Ji[k][j] * a[j];
      }
    }
    double bNorm = sqrt(b[0] * b[0] + b[1] * b[1]);
    b[0] /= bNorm;
    b[1] /= bNorm;
    db = - b[0] * xi[0] - b[1] * xi[1];
//   std::cout << b[0] << " " << b[1] << " " << db << " " << std::endl;

//     // Old inverse mapping for comparison
//
//     std::vector <  std::vector < std::vector <double > > > aP(1);
//     short unsigned quad = 3;
//     unsigned linear = 0;
//     bool ielIsInitialized = false;
//     if(!ielIsInitialized) {
//       ielIsInitialized = true;
//       ProjectNodalToPolynomialCoefficients(aP[0], xv, quad, linear) ;
//     }
//
//     std::vector<double> xib(dim);
//     GetClosestPointInReferenceElement(xv, xm, quad, xib);
//     bool inverseMapping = GetInverseMapping(linear, quad, aP, xm, xib, 100);
//     if(!inverseMapping) {
//       std::cout << "InverseMapping failed" << std::endl;
//     }
//
//     //std::cout << xib[0] << " " << xib[1] << std::endl;
//
//     vector < vector < double > > Jac2;
//     vector < vector < double > > JacI2;
//     finiteElementQuad->GetJacobianMatrix(xv, xib, Jac2, JacI2);
//
//
//     //std::swap(JacI[0][1],JacI[1][0]);
//     //std::cout << Jac[0][0] << " " << Jac[0][1] << " " << Jac[1][0] << " " <<Jac[1][1] << std::endl;
//
//     std::vector <double> b2(dim, 0.);
//
//     for(unsigned k = 0; k < dim; k++) {
//       for(unsigned j = 0; j < dim; j++) {
//         b2[k] += JacI2[k][j] * a[j];
//       }
//     }
//     double b2Norm = sqrt(b2[0] * b2[0] + b2[1] * b2[1]);
//     b2[0] /= b2Norm;
//     b2[1] /= b2Norm;
//     double db2 = - b2[0] * xi[0] - b2[1] * xi[1];
//     std::cout << b2[0] << " " << b2[1] << " " << db2 << " " << std::endl;

  }

}

void GetNormalTri(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &a, double &d,  std::vector<double> &xm, std::vector<double> &b, double &db, unsigned &cut) {

  const unsigned &dim =  xv.size();
  const unsigned &nve =  xv[0].size();

  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];

  double hx = (fabs(x2 - x1) + fabs(x3 - x2) + fabs(x3 - x1)) / 3.;
  double hy = (fabs(y2 - y1) + fabs(y3 - y2) + fabs(y3 - y1)) / 3.;

  double h = sqrt(hx * hx + hy * hy);
  double eps = 1.0e-10 * h;

  std::vector<double> dist(nve, 0);
  std::vector<double> dist0(nve);
  unsigned cnt0 = 0;
  for(unsigned i = 0; i < nve; i++) {
    for(unsigned k = 0;  k < dim; k++) {
      dist[i] += (xv[k][i] - xg[k]) * (xv[k][i] - xg[k]);
    }
    dist[i] = sqrt(dist[i]) - R;

    if(fabs(dist[i]) < eps) {
      dist0[i] = (dist[i] < 0) ? -eps : eps;
      dist[i] = 0.;
      cnt0++;
    }
    else {
      dist0[i] = dist[i];
    }
  }

  if(cnt0 > 0) {
    unsigned cntp = 0;
    for(unsigned i = 0; i < nve; i++) {
      if(dist[i] > 0) cntp++;
      dist[i] = dist0[i];
    }
    if(cntp == 0) { // the element is inside the ball
      cut = 0;
      return;
    }
    else if(cntp == nve - cnt0) {  // the element in outside the ball
      cut = 2;
      return;
    }
  }

  std::vector <double> theta(2);
  unsigned cnt = 0;
  for(unsigned e = 0; e < nve; e++) {
    unsigned ep1 = (e + 1) % nve;
    if(dist[e] * dist[ep1] < 0) {
      double s = 0.5  * (1 + (dist[e] + dist[ep1]) / (dist[e] - dist[ep1]));
      theta[cnt] = atan2((1 - s) * xv[1][e] + s * xv[1][ep1]  - xg[1], (1 - s) * xv[0][e] + s * xv[0][ep1] - xg[0]) ;
      cnt++;
    }
  }

  if(cnt == 0) {
    if(dist[0] < 0) cut = 0; // cell inside the ball
    else cut = 2; // cell outside the ball
    return;
  }
  else {
    cut = 1;
    if(theta[0] > theta[1]) {
      std::swap(theta[0], theta[1]);
    }
    double DT = theta[1] - theta[0];
    if(DT > M_PI) {
      std::swap(theta[0], theta[1]);
      theta[1] += 2. * M_PI;
      DT = theta[1] - theta[0];
    }
    xm.resize(dim);

    d = R * sqrt(0.5 * DT / tan(0.5 * DT)) ;
    a.resize(dim);
    a[0] = -cos(theta[0] + 0.5 * DT);
    a[1] = -sin(theta[0] + 0.5 * DT);

    for(unsigned k = 0; k < dim; k++) {
      xm[k] = -a[k] * d + xg[k];
    }
    d += - a[0] * xg[0] - a[1] * xg[1]; //TODO

    std::cout.precision(14);

//     std::cout << "xm = " << xm[0] << " " << xm[1] << std::endl;
//     std::cout << "a = " << a[0] << " b = " << a[1] << " d = " << d << std::endl;

    std::vector<double> xi(dim);

    std::vector < std::vector < double > > J(2, std::vector<double>(2));
    J[0][0] = (-x1 + x2);
    J[0][1] = (-x1 + x3);

    J[1][0] = (-y1 + y2);
    J[1][1] = (-y1 + y3);

    double den = (x3 * y1 - x1 * y3 + x2 * J[1][1] - y2 * J[0][1]);

    xi[0] = (x3 * y1 - x1 * y3 + xm[0] * J[1][1] - xm[1] * J[0][1]) / den;
    xi[1] = (x1 * y2 - x2 * y1 - xm[0] * J[1][0] + xm[1] * J[0][0]) / den;


    double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    std::vector < std::vector < double > > Ji = {{J[1][1] / det, -J[0][1] / det}, {-J[1][0] / det, J[0][0] / det}};

    b.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        b[k] += Ji[k][j] * a[j];
      }
    }
    double bNorm = sqrt(b[0] * b[0] + b[1] * b[1]);
    b[0] /= bNorm;
    b[1] /= bNorm;
    db = - b[0] * xi[0] - b[1] * xi[1];
  }
}



void SimpleNonlocalAssembly(MultiLevelProblem & ml_prob) {

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("FracProblem");

  const unsigned level = N_UNIFORM_LEVELS - N_ERASED_LEVELS - 1;


  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension();
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  constexpr unsigned int space_dim = 3;

  unsigned    iproc = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  vector < vector < double > > x1(dim);    // local coordinates
  vector < vector < double > > x2(dim);    // local coordinates
  for(unsigned k = 0; k < dim; k++) {
    x1[k].reserve(maxSize);
    x2[k].reserve(maxSize);
  }

  std::vector < std::vector < double > >  Jac_qp(dim);
  std::vector < std::vector < double > >  JacI_qp(space_dim);
  double detJac_qp;

  CurrentElem < double > geom_element1(dim, msh);            // must be adept if the domain is moving, otherwise double
  CurrentElem < double > geom_element2(dim, msh);

  unsigned soluIndex = mlSol->GetIndex("u");
  unsigned solType   = mlSol->GetSolutionType(soluIndex);

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < std::vector < /*const*/ elem_type_templ_base< double, double > *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);

  FILE * fp1;

  fp1 = fopen("lines_qp.dat", "w");

  for(int kproc = 0; kproc < nprocs; kproc++) {
    for(int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++) {
      short unsigned ielGeom2;
      unsigned nDof2;
      unsigned nDofx2;
      unsigned nDofLin;
      unsigned nDofu2;
      //unsigned n_face;

      if(iproc == kproc) {
        ielGeom2 = msh->GetElementType(jel);
        nDof2  = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs
        nDofx2 = msh->GetElementDofNumber(jel, xType);    // number of coordinate element dofs
        nDofLin = 4; //TODO

      }

      MPI_Bcast(&ielGeom2, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDof2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDofx2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      //MPI_Bcast(&n_face, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);

      for(int k = 0; k < dim; k++) {
        x2[k].resize(nDofx2);
      }

      vector < double > phi;
      vector < double > phi_x;

      phi.reserve(maxSize);
      phi_x.reserve(maxSize * dim);


      // local storage of coordinates  #######################################
      if(iproc == kproc) {
        for(unsigned j = 0; j < nDofx2; j++) {
          unsigned xDof  = msh->GetSolutionDof(j, jel, xType);  // global to global mapping between coordinates node and coordinate dof
          for(unsigned k = 0; k < dim; k++) {
            x2[k][j] = (*msh->_topology->_Sol[k])(xDof);  // global extraction and local storage for the element coordinates
          }
        }
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& x2[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }

      if(iproc == kproc) {
        geom_element2.set_coords_at_dofs_and_geom_type(jel, xType);
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& geom_element2.get_coords_at_dofs()[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
      for(unsigned k = 0; k < space_dim; k++) {
        MPI_Bcast(& geom_element2.get_coords_at_dofs_3d()[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }

      const unsigned jgNumber = msh->_finiteElement[ielGeom2][solType]->GetGaussPointNumber();
//       const unsigned jgNumber = ml_prob.GetQuadratureRule(ielGeom2).GetGaussPointsNumber();

      vector < vector < double > > xg2(jgNumber);
      vector <double> weight2(jgNumber);
      vector < vector <double> > phi2(jgNumber);  // local test function

      for(unsigned jg = 0; jg < jgNumber; jg++) {

        msh->_finiteElement[ielGeom2][solType]->Jacobian(x2, jg, weight2[jg], phi2[jg], phi_x);

//         elem_all[ielGeom2][xType]->JacJacInv(/*x2*/geom_element2.get_coords_at_dofs_3d(), jg, Jac_qp, JacI_qp, detJac_qp, space_dim);
//         weight2[jg] = detJac_qp * ml_prob.GetQuadratureRule(ielGeom2).GetGaussWeightsPointer()[jg];
//         elem_all[ielGeom2][solType]->shape_funcs_current_elem(jg, JacI_qp, phi2[jg], phi_x /*boost::none*/, boost::none /*phi_u_xx*/, space_dim);




        xg2[jg].assign(dim, 0.);

        for(unsigned j = 0; j < nDof2; j++) {
          for(unsigned k = 0; k < dim; k++) {
            xg2[jg][k] += x2[k][j] * phi2[jg][j];
          }
        }
      }

      std::vector <int> bd_face(0);
      unsigned nFaces;

      for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        short unsigned ielGeom1 = msh->GetElementType(iel);
        unsigned nDof1  = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
        unsigned nDofx1 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

        for(int k = 0; k < dim; k++) {
          x1[k].resize(nDofx1);
        }

        // local storage of coordinates
        for(unsigned i = 0; i < nDofx1; i++) {
          unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
          for(unsigned k = 0; k < dim; k++) {
            x1[k][i] = (*msh->_topology->_Sol[k])(xDof);
          }
        }

        const unsigned igNumber = msh->_finiteElement[ielGeom1][solType]->GetGaussPointNumber();
        double weight1;
        vector < double > phi1;  // local test function


        for(unsigned ig = 0; ig < igNumber; ig++) {
          msh->_finiteElement[ielGeom1][solType]->Jacobian(x1, ig, weight1, phi1, phi_x);
          vector < double > xg1(dim, 0.);
          for(unsigned i = 0; i < nDof1; i++) {
            for(unsigned k = 0; k < dim; k++) {
              xg1[k] += x1[k][i] * phi1[i];
            }
          }
          double R = 0.5;
          std::vector<double> a;
          std::vector<double> b;
          double d;
          double db;
          unsigned cut;
          std::vector<double> xm;
          GetNormalQuad(x2, xg1, R, a, d, xm, b, db, cut);

          if(cut && iel == 110 && ig == 0) { //TODO

            double xx = xm[0] - 0.5 * a[1];
            double yy = xm[1] + 0.5 * a[0];
            fprintf(fp1, "%f %f \n", xx, yy);
            xx = xm[0] + 0.5 * a[1];
            yy = xm[1] - 0.5 * a[0];
            fprintf(fp1, "%f %f \n", xx, yy);
            fprintf(fp1, "\n \n");
          }

        }


      }
    }
  }
  fclose(fp1);
}




