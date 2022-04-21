
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "CurrentElem.hpp"
#include "LinearImplicitSystem.hpp"

#include "PolynomialBases.hpp"

#include "CutFemWeight.hpp"

#include "CDWeights.hpp"

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

void GetNormalTet(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &a, double &d,  std::vector<double> &xm, std::vector<double> &b, double &db, double &vol, unsigned &cut);

double getHeightPolyhedronSphereInt(const std::vector < std::vector <double> > &aN, const std::vector <double> &a, const std::vector <double> &xg, const double &R);


const elem_type *finiteElementQuad;

int main(int argc, char** argv) {

  std::vector<std::vector<double>> xt = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  std::vector< double > xg1 = {0.0, 0., 1.5};
  double R1 = 1.;
  std::vector<double> a10;
  double d10;
  std::vector<double> xm1;
  std::vector<double> b1;
  double db1;
  double vol1;
  unsigned cut1;

  GetNormalTet(xt, xg1, R1, a10, d10, xm1, b1, db1, vol1, cut1);


//   return 1;


  typedef double TypeIO;
  typedef cpp_bin_float_oct TypeA;
  

  unsigned qM = 3;
  CutFemWeight <TypeIO, TypeA> quad  = CutFemWeight<TypeIO, TypeA >(QUAD, qM, "legendre");
  CutFemWeight <TypeIO, TypeA> tri  = CutFemWeight<TypeIO, TypeA >(TRI, qM, "legendre");
  CutFemWeight <TypeIO, TypeA> tet  = CutFemWeight<TypeIO, TypeA >(TET, qM, "legendre");

  std::vector<double> weight1;
  tet.GetWeightWithMap(0, {-0.034878236872063, 0.0012179748700879, 0.9993908270191}, 0.096573056501712, weight1 );
  for(unsigned j = 0; j < 1; j++) {
    std::cout << weight1[j] << "\n ";
  }
  
  //abort();
  

  double dx = .05;
  double dt = 2.;
  CDWeightQUAD <TypeA> quadCD(qM, dx, dt);
  CDWeightTRI <TypeA> triCD(qM, dx, dt);
  CDWeightTET <TypeA> tetCD(qM, dx, dt);

  std::cout<<std::endl;
  
  double theta1 = 45;
  double phi1 = 30;
  std::vector<double> a1 = {cos(theta1 * M_PI / 180), -sin(theta1 * M_PI / 180)};
  std::vector<double> a2 = {cos(theta1 * M_PI / 180)* sin(phi1 * M_PI / 180), sin(theta1 * M_PI / 180) * sin(phi1 * M_PI / 180), cos(phi1 * M_PI / 180)};
  double d1 = -0.1 * sqrt(2);

 
//   quad.GetWeightWithMap(0, a1, d1, weight1);
//
//   for(unsigned j = 0; j < weight1.size(); j++) {
//     std::cout << weight1[j] << " ";
//   }
//   std::cout << std::endl;
  std::vector<double> weight;
//   quadCD.GetWeight(a1, d1, weight);
//
//   for(unsigned j = 0; j < weight.size(); j++) {
//     std::cout << weight[j] << " ";
//   }
//   std::cout << std::endl;
//
//
//   tri.GetWeightWithMap(0, a1, d1, weight1);
//   for(unsigned j = 0; j < weight1.size(); j++) {
//     std::cout << weight1[j] << " ";
//   }
//   std::cout << std::endl;
//   triCD.GetWeight(a1, d1, weight);
//   for(unsigned j = 0; j < weight.size(); j++) {
//     std::cout << weight[j] << " ";
//   }
//   std::cout << std::endl;

  tet.GetWeightWithMap(0, a2, d1, weight1);
  for(unsigned j = 0; j < weight1.size(); j++) {
    std::cout << weight1[j] << " ";
  }
  std::cout << std::endl;
  tetCD.GetWeight(a2, d1, weight);
  for(unsigned j = 0; j < weight.size(); j++) {
    std::cout << weight[j] << " ";
  }
  std::cout << std::endl;
// 
// 
// 
//   return 1;

  const std::string fe_quad_rule_1 = "seventh";
  const std::string fe_quad_rule_2 = "eighth";

// ======= Init ========================
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

//   finiteElementQuad = new const elem_type_2D("quad", "linear", "fifth", "legendre");

  unsigned numberOfUniformLevels = N_UNIFORM_LEVELS;

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;
//   mlMsh.GenerateCoarseBoxMesh(N_X, N_Y, 0, EX_1, EX_2, EY_1, EY_2, 0., 0., QUAD9, fe_quad_rule_1.c_str());
  mlMsh.GenerateCoarseBoxMesh(N_X, N_Y, 0, EX_1, EX_2, EY_1, EY_2, 0., 0., TRI6, fe_quad_rule_1.c_str());

//   mlMsh.ReadCoarseMesh("./input/cube_tet.neu", fe_quad_rule_1.c_str(), 1.);

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


  std::vector<double> xg(dim);
  xg[0] = -0.04852;
  xg[1] = -0.0017;
  if(dim == 3) xg[2] = 0.0320;
  double R = 0.25;
  double volumeBall = 0.;

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
        x1[k][(i + 2) % nDof] = (*msh->_topology->_Sol[k])(xDof); // global extraction and local storage for the element coordinates
      }
    }

    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> xm;
    double d;
    double db;
    unsigned cut;

    double vol;

    if(ielType == 3) {
      GetNormalQuad(x1, xg, R, a, d, xm, b, db, cut);
      vol = ((EX_2 - EX_1) / (N_X * pow(2, N_UNIFORM_LEVELS - 1))) * ((EY_2 - EY_1) / (N_Y * pow(2, N_UNIFORM_LEVELS - 1)));
    }
    else if(ielType == 4) {
      GetNormalTri(x1, xg, R, a, d, xm, b, db, cut);
      vol = 0.5 * ((EX_2 - EX_1) / (N_X * pow(2, N_UNIFORM_LEVELS - 1))) * ((EY_2 - EY_1) / (N_Y * pow(2, N_UNIFORM_LEVELS - 1)));
    }
    else if(ielType == 1) {
      GetNormalTet(x1, xg, R, a, d, xm, b, db, vol, cut);
      //std::cout << cut <<" ";
    }

    if(cut == 1) {
      bool wMap = 1;
      if(ielType == 3) {
        std::vector <TypeIO> weightCF;
        //quad.GetWeightWithMap(0, b, db, weightCF);

        quadCD.GetWeight(b, db, weightCF);

        const double* weightG = quad.GetGaussWeightPointer();

        double sum = 0.;
        for(unsigned ig = 0; ig < weightCF.size(); ig++) {
          sum += weightG[ig] * weightCF[ig] ; //TODO use the correct quad rule!!!!!
        }
        volumeBall += sum / 4. * vol;
      }
      else if(ielType == 4) {
        std::vector <TypeIO> weightCF;
        //tri.GetWeightWithMap(0, b, db, weightCF);

        triCD.GetWeight(b, db, weightCF);
        const double* weightG = tri.GetGaussWeightPointer();

        double sum = 0.;
        for(unsigned ig = 0; ig < weightCF.size(); ig++) {
          sum += weightG[ig] * weightCF[ig] ;
        }
        //std::cout << sum << " ";
        volumeBall += 2. * sum * vol;
      }
      else if(ielType == 1) {
        std::vector <TypeIO> weightCF;
        //tet.clear();
        //tet.GetWeightWithMap(0, b, db, weightCF);
        tetCD.GetWeight(b, db, weightCF);

        //std::cout<<"a"<<std::flush;
        
        const double* weightG = tet.GetGaussWeightPointer();

        double sum = 0.;
        for(unsigned ig = 0; ig < weightCF.size(); ig++) {
          sum += weightG[ig] * weightCF[ig] ;
        }
        volumeBall += 6. * sum * vol;

      }


    }
    else if(cut == 0) {
      volumeBall += vol;
    }

    /* trivial print for xmgrace */
//     if(cut == 1) {
//       double xx = xm[0] - a[1] * d;
//       double yy = xm[1] + a[0] * d;
//       fprintf(fp, "%f %f \n", xx, yy);
//       xx = xm[0] + a[1] * d;
//       yy = xm[1] - a[0] * d;
//       fprintf(fp, "%f %f \n", xx, yy);
//       fprintf(fp, "\n \n");
//     }
  }

  std::cout << "numnber of calls in QUAD " << quad.GetCounter() << std::endl;
  std::cout << "numnber of calls in TRI  " << tri.GetCounter() << std::endl;
  std::cout << "numnber of calls in TET  " << tet.GetCounter() << std::endl;


  std::cout.precision(14);
  std::cout << "volume of the Ball  in " << dim << "D = " << volumeBall << "  analytic value = " << ((dim == 2) ? M_PI * R * R : 4. / 3. * M_PI * pow(R, 3)) << "\n";

  mlSol.SetWriter(VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back("All");
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "quadratic", print_vars, 0);

  delete finiteElementQuad;
  return 1;
}


void GetNormalQuad(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double & R, std::vector<double> &a, double & d,  std::vector<double> &xm, std::vector<double> &b, double & db, unsigned & cut) {

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
  
  
  std::vector<double> A(2, 0.);
  std::vector<std::vector<double>> xe(2, std::vector<double>(4));
  double D = 0.;
  unsigned intMax = 2;
  unsigned nEdge = 0;
  unsigned cnt = 0;
  
  for(unsigned i = 0; i < nve; i++) {
    unsigned ip1 = (i + 1) % nve;
    A[0] = xv[1][ip1] - xv[1][i];
    A[1] = - xv[0][ip1] + xv[0][i];
    D = - A[0] * xv[0][i] - A[1] * xv[1][i];
   
    
    std::vector<double> inters(intMax, 0.);
    unsigned dir = (fabs(A[0]) > fabs(A[1]))? 1 : 0 ;
    unsigned dirp1 = (dir + 1)%2;
    
    double iMax = std::max(xv[dir][ip1], xv[dir][i]);
    double iMin = std::min(xv[dir][ip1], xv[dir][i]);
    
    double delta = ((A[0] * A[0] + A[1] * A[1])* R * R) - (D + A[0] * xg[0] + A[1] * xg[1]) * (D + A[0] * xg[0] + A[1] * xg[1]);
    a.resize(dim);
    if(delta > 0.){
      inters[0] = (- A[dir] * (D + A[dirp1] * xg[dirp1]) + A[dirp1] * (A[dirp1] * xg[dir] - sqrt(delta))) / (A[0] * A[0] + A[1] * A[1]);  
      inters[1] = (- A[dir] * (D + A[dirp1] * xg[dirp1]) + A[dirp1] * (A[dirp1] * xg[dir] + sqrt(delta))) / (A[0] * A[0] + A[1] * A[1]);  
      unsigned nInt = 0;
      unsigned jInt = 2;
      for(unsigned j = 0; j < intMax; j++){
        if(inters[j] < iMax && inters[j] > iMin) {
            nInt++;
            jInt = j;
        }
      }
      if(nInt == 1){   
          xe[dir][cnt] = inters[jInt];
          xe[dirp1][cnt] = ( - D - A[dir] * xe[dir][cnt] ) / A[dirp1];
          cnt++;
      }
    }
  }
  if(cnt == 0) cut = ( R * R - (xv[0][0] - xg[0]) * (xv[0][0] - xg[0]) - (xv[1][0] - xg[1]) * (xv[1][0] - xg[1]) > 0 ) ? 0 : 2;
  else if (cnt == 4) cut = 0;
  else if(cnt == 2){
    cut = 1;
    std::vector<double> theta(2);
    
    a[0] = xe[1][1] - xe[1][0] ;
    a[1] = - xe[0][1] + xe[0][0] ;
    
    xm.resize(2);
    xm[0] = 0.5 * (xe[0][0] + xe[0][1]);
    xm[1] = 0.5 * (xe[1][0] + xe[1][1]);
    
    double det = 0;
    for(unsigned k = 0; k < dim; k++) {
      det += a[k] * ( xg[k] - xm[k] );
    }
    double sign = (det>=0) ? 1. : -1.;
    
    double norm = sign * sqrt(a[0] * a[0] + a[1] * a[1]);
    a[0] /= norm;
    a[1] /= norm;
    
    theta[0] = atan2(xe[1][0]  - xg[1], xe[0][0] - xg[0]);
    theta[1] = atan2(xe[1][1]  - xg[1], xe[0][1] - xg[0]);
    
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

    double d2 = sqrt(pow(xm[0] - xg[0], 2) + pow(xm[1] - xg[1], 2));
    d = d2 * tan(0.5 * DT);
    
    std::cout << "xm = " << xm[0] << " " << xm[1] << std::endl;
    std::cout << "a = " << a[0] << " b = " << a[1] << std::endl;
    
    
  
  
  
//   for(unsigned i = 0; i < nve; i++) {
//     for(unsigned k = 0;  k < dim; k++) {
//       dist[i] += (xv[k][i] - xg[k]) * (xv[k][i] - xg[k]);
//     }
//     dist[i] = sqrt(dist[i]) - R;
// 
//     if(fabs(dist[i]) < eps) {
//       dist0[i] = (dist[i] < 0) ? -eps : eps;
//       dist[i] = 0.;
//       cnt0++;
//     }
//     else {
//       dist0[i] = dist[i];
//     }
//   }
// 
//   if(cnt0 > 0) {
//     unsigned cntp = 0;
//     for(unsigned i = 0; i < nve; i++) {
//       if(dist[i] > 0) cntp++;
//       dist[i] = dist0[i];
//     }
//     if(cntp == 0) { // the element is inside the ball
//       cut = 0;
//       return;
//     }
//     else if(cntp == nve - cnt0) {  // the element in outside the ball
//       cut = 2;
//       return;
//     }
//   }
// 
//   std::vector <double> theta(2);
//   unsigned cnt = 0;
//   for(unsigned e = 0; e < nve; e++) {
//     unsigned ep1 = (e + 1) % nve;
//     if(dist[e] * dist[ep1] < 0) {
//       double s = 0.5  * (1 + (dist[e] + dist[ep1]) / (dist[e] - dist[ep1]));
//       theta[cnt] = atan2((1 - s) * xv[1][e] + s * xv[1][ep1]  - xg[1], (1 - s) * xv[0][e] + s * xv[0][ep1] - xg[0]) ;
//       cnt++;
//     }
//   }
// 
//   if(cnt == 0) {
//     if(dist[0] < 0) cut = 0; // cell inside the ball
//     else cut = 2; // cell outside the ball
//     return;
//   }
//   else {
//     cut = 1;
//     if(theta[0] > theta[1]) {
//       std::swap(theta[0], theta[1]);
//     }
//     double DT = theta[1] - theta[0];
//     if(DT > M_PI) {
//       std::swap(theta[0], theta[1]);
//       theta[1] += 2. * M_PI;
//       DT = theta[1] - theta[0];
//     }
//     xm.resize(dim);
// 
//     d = R * sqrt(0.5 * DT / tan(0.5 * DT)) ;
//     a.resize(dim);
//     a[0] = -cos(theta[0] + 0.5 * DT);
//     a[1] = -sin(theta[0] + 0.5 * DT);
// 
//     for(unsigned k = 0; k < dim; k++) {
//       xm[k] = -a[k] * d + xg[k];
//     }
//     d += - a[0] * xg[0] - a[1] * xg[1]; //TODO
// 
//     double d2 = sqrt(pow(xm[0] - xg[0], 2) + pow(xm[1] - xg[1], 2));
//     d = d2 * tan(0.5 * DT);
// 
//     std::cout.precision(14);
// 
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

//     double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
//     std::vector < std::vector < double > > Ji = {{J[1][1] / det, -J[0][1] / det}, {-J[1][0] / det, J[0][0] / det}};

    b.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        b[k] += J[j][k] * a[j];
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
//         b2[k] += JacI2[j][k] * a[j];
//       }
//     }
//     double b2Norm = sqrt(b2[0] * b2[0] + b2[1] * b2[1]);
//     b2[0] /= b2Norm;
//     b2[1] /= b2Norm;
//     double db2 = - b2[0] * xi[0] - b2[1] * xi[1];
//     std::cout << b2[0] << " " << b2[1] << " " << db2 << " " << std::endl;

  }

}

void GetNormalTri(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double & R, std::vector<double> &a, double & d,  std::vector<double> &xm, std::vector<double> &b, double & db, unsigned & cut) {

  const unsigned &dim =  xv.size();
  const unsigned &nve =  xv[0].size();
  
  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];
  
  
  
  
  
  
//   std::vector<double> A(2, 0.);
//   std::vector<std::vector<double>> xe(2, std::vector<double>(4));
//   double D = 0.;
//   unsigned intMax = 2;
//   unsigned nEdge = 0;
//   unsigned cnt = 0;
//   
//   for(unsigned i = 0; i < nve; i++) {
//     unsigned ip1 = (i + 1) % nve;
//     A[0] = xv[1][ip1] - xv[1][i];
//     A[1] = - xv[0][ip1] + xv[0][i];
//     D = - A[0] * xv[0][i] - A[1] * xv[1][i];
//    
//     
//     std::vector<double> inters(intMax, 0.);
//     unsigned dir = (fabs(A[0]) > fabs(A[1]))? 1 : 0 ;
//     unsigned dirp1 = (dir + 1)%2;
//     
//     double iMax = std::max(xv[dir][ip1], xv[dir][i]);
//     double iMin = std::min(xv[dir][ip1], xv[dir][i]);
//     
//     double delta = ((A[0] * A[0] + A[1] * A[1])* R * R) - (D + A[0] * xg[0] + A[1] * xg[1]) * (D + A[0] * xg[0] + A[1] * xg[1]);
//     a.resize(dim);
//     if(delta > 0.){
//       inters[0] = (- A[dir] * (D + A[dirp1] * xg[dirp1]) + A[dirp1] * (A[dirp1] * xg[dir] - sqrt(delta))) / (A[0] * A[0] + A[1] * A[1]);  
//       inters[1] = (- A[dir] * (D + A[dirp1] * xg[dirp1]) + A[dirp1] * (A[dirp1] * xg[dir] + sqrt(delta))) / (A[0] * A[0] + A[1] * A[1]);  
//       unsigned nInt = 0;
//       unsigned jInt = 2;
//       for(unsigned j = 0; j < intMax; j++){
//         if(inters[j] < iMax && inters[j] > iMin) {
//             nInt++;
//             jInt = j;
//         }
//       }
//       if(nInt == 1){   
//           xe[dir][cnt] = inters[jInt];
//           xe[dirp1][cnt] = ( - D - A[dir] * xe[dir][cnt] ) / A[dirp1];
//           cnt++;
//       }
//     }
//   }
//   if(cnt == 0) cut = ( R * R - (xv[0][0] - xg[0]) * (xv[0][0] - xg[0]) - (xv[1][0] - xg[1]) * (xv[1][0] - xg[1]) > 0 ) ? 0 : 2;
//   else if (cnt == 4) cut = 0;
//   else if(cnt == 2){
//     cut = 1;
//     std::vector<double> theta(2);
//     
//     a[0] = xe[1][1] - xe[1][0] ;
//     a[1] = - xe[0][1] + xe[0][0] ;
//     
//     xm.resize(2);
//     xm[0] = 0.5 * (xe[0][0] + xe[0][1]);
//     xm[1] = 0.5 * (xe[1][0] + xe[1][1]);
//     
//     double det = 0;
//     for(unsigned k = 0; k < dim; k++) {
//       det += a[k] * ( xg[k] - xm[k] );
//     }
//     double sign = (det>=0) ? 1. : -1.;
//     
//     double norm = sign * sqrt(a[0] * a[0] + a[1] * a[1]);
//     a[0] /= norm;
//     a[1] /= norm;
//     
//     theta[0] = atan2(xe[1][0]  - xg[1], xe[0][0] - xg[0]);
//     theta[1] = atan2(xe[1][1]  - xg[1], xe[0][1] - xg[0]);
//     
//     if(theta[0] > theta[1]) {
//       std::swap(theta[0], theta[1]);
//     }
//     double DT = theta[1] - theta[0];
//     if(DT > M_PI) {
//       std::swap(theta[0], theta[1]);
//       theta[1] += 2. * M_PI;
//       DT = theta[1] - theta[0];
//     }
//     xm.resize(dim);
// 
//     d = R * sqrt(0.5 * DT / tan(0.5 * DT)) ;
//     a.resize(dim);
//     a[0] = -cos(theta[0] + 0.5 * DT);
//     a[1] = -sin(theta[0] + 0.5 * DT);
// 
//     for(unsigned k = 0; k < dim; k++) {
//       xm[k] = -a[k] * d + xg[k];
//     }
//     d += - a[0] * xg[0] - a[1] * xg[1]; //TODO
// 
//     double d2 = sqrt(pow(xm[0] - xg[0], 2) + pow(xm[1] - xg[1], 2));
//     d = d2 * tan(0.5 * DT);
//     
//     std::cout << "xm = " << xm[0] << " " << xm[1] << std::endl;
//     std::cout << "a = " << a[0] << " b = " << a[1] << std::endl;
    
    
    
    
    
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

    //std::cout << dist[i] << std::endl;

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

    //std::cout << "xm = " << xm[0] << " " << xm[1] << std::endl;
    //std::cout << "a = " << a[0] << " b = " << a[1] << " d = " << d << std::endl;

    double d2 = sqrt(pow(xm[0] - xg[0], 2) + pow(xm[1] - xg[1], 2));
    d = d2 * tan(0.5 * DT);

    std::cout.precision(14);



    std::vector<double> xi(dim);

    std::vector < std::vector < double > > J(2, std::vector<double>(2));
    J[0][0] = (-x1 + x2);
    J[0][1] = (-x1 + x3);

    J[1][0] = (-y1 + y2);
    J[1][1] = (-y1 + y3);

    double den = (x3 * y1 - x1 * y3 + x2 * J[1][1] - y2 * J[0][1]);

    xi[0] = (x3 * y1 - x1 * y3 + xm[0] * J[1][1] - xm[1] * J[0][1]) / den;
    xi[1] = (x1 * y2 - x2 * y1 - xm[0] * J[1][0] + xm[1] * J[0][0]) / den;


    //std::cout << xi[0] << " " << xi[1] << std::endl;


    b.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        b[k] += J[j][k] * a[j];
      }
    }
    double bNorm = sqrt(b[0] * b[0] + b[1] * b[1]);
    b[0] /= bNorm;
    b[1] /= bNorm;
    db = - b[0] * xi[0] - b[1] * xi[1];


    //std::cout << b[0] << " " << b[1] << " " << db << " " << std::endl;
  }
}








void GetNormalTet(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double & R, std::vector<double> &a, double & d,  std::vector<double> &xm, std::vector<double> &a2, double & d2, double & volume,  unsigned & cut) {

  const unsigned dim =  3;
  const unsigned nve =  4;

  //std::cout<<nve<<std::endl;

  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& x4 = xv[0][3];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];
  const double& y4 = xv[1][3];
  const double& z1 = xv[2][0];
  const double& z2 = xv[2][1];
  const double& z3 = xv[2][2];
  const double& z4 = xv[2][3];

  double hx = (fabs(x2 - x1) + fabs(x3 - x2) + fabs(x3 - x1) + fabs(x4 - x1) + fabs(x4 - x2) + fabs(x4 - x3)) / 6.;
  double hy = (fabs(y2 - y1) + fabs(y3 - y2) + fabs(y3 - y1) + fabs(y4 - y1) + fabs(y4 - y2) + fabs(y4 - y3)) / 6.;
  double hz = (fabs(z2 - z1) + fabs(z3 - z2) + fabs(z3 - z1) + fabs(z4 - z1) + fabs(z4 - z2) + fabs(z4 - z3)) / 6.;
  double h = sqrt(hx * hx + hy * hy + hz * hz);
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

  std::vector < std::vector <double> > y(4, std::vector<double>(dim));
  std::vector < unsigned > i0(4);
  unsigned cnt = 0;
  for(unsigned i = 0; i < nve - 1; i++) {
    for(unsigned j = i + 1; j < nve; j++) {
      if(dist[i] * dist[j] < 0) {
        double s = dist[i] / (dist[i] - dist[j]);
        for(unsigned k = 0; k < dim; k++) {
          y[cnt][k] = (1. - s) * xv[k][i] + s * xv[k][j];
        }
        i0[cnt] = (i + j) - (i == 0);
        cnt++;
      }
    }
  }

  if(cnt == 0) {
    if(dist[0] < 0) cut = 0; // cell inside the ball
    else cut = 2; // cell outside the ball
    return;
  }
  else {
    cut = 1;

    if(cnt == 4) {

      if((i0[0] == 0 && i0[1] == 1) || (i0[0] == 1 && i0[1] == 2)) {
        std::swap(y[2], y[3]);
        std::swap(i0[2], i0[3]);
      }
      else {
        std::swap(y[1], y[3]);
        std::swap(y[1], y[2]);

        std::swap(i0[1], i0[3]);
        std::swap(i0[1], i0[2]);
      }

      //std::cout << i0[0] << " " << i0[1] << " " << i0[2] << " " << i0[3] << std::endl;
    }

    std::vector <double> yg(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned i = 0; i < cnt; i++) {
        yg[k] += y[i][k];
      }
      yg[k] /= cnt;
    }

    a.resize(dim);

    std::vector < std::vector <double> > b(cnt, std::vector<double>(dim));
    for(unsigned k = 0; k < dim; k++) {
      a[k] =  yg[k] - xg[k];
      for(unsigned i = 0; i < cnt; i++) {
        b[i][k] = y[i][k] - xg[k];
      }
    }
    double an = 0.;
    std::vector <double> bn(cnt, 0);
    for(unsigned k = 0; k < dim; k++) {
      an += a[k] * a[k];
      for(unsigned i = 0; i < cnt; i++) {
        bn[i] += b[i][k] * b[i][k];
      }
    }
    an = sqrt(an);
    for(unsigned i = 0; i < cnt; i++) {
      bn[i] = sqrt(bn[i]);
    }

    for(unsigned k = 0; k < dim; k++) {
      a[k] /= an;
      for(unsigned i = 0; i < cnt; i++) {
        b[i][k] /= bn[i];
      }
    }




//     double phig = 0;
//     for(unsigned i = 0; i < cnt; i++) {
//       double phii = 0;
//       for(unsigned k = 0; k < dim; k++) {
//         phii += a[k] * b[i][k];
//       }
//       phii = acos(phii / (an * bn[i]));
//       phig += phii;
//     }
//     phig /= cnt;
    //double H = R * pow(2. * (1. - cos(phig)) / (tan(phig) * tan(phig)), 1. / 3.);

    double H = getHeightPolyhedronSphereInt(b, a, xg, R);
    //else H = R * pow(2. * cos(phig) * cos(phig) / (1. + cos(phig)), 1. / 3.);;



//     if(cnt == 3) H = getHeightPolyhedronSphereInt(bn, a, R);
//     else H = R * pow(2. * cos(phig) * cos(phig) / (1. + cos(phig)), 1. / 3.);

    // H = R * pow(2. * cos(phig) * cos(phig) / (1. + cos(phig)), 1. / 3.);


    xm.resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      xm[k] = xg[k] + a[k] * H;
    }

    // std::cout << "\nBBB " << H << " " << a[0] << " " << a[1] << " " << a[2] << " \n" << xm[0] << " " << xm[1] << " " << xm[2] << std::endl;


    /*
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
        d += - a[0] * xg[0] - a[1] * xg[1]; //TODO*/

    //std::cout << "xm = " << xm[0] << " " << xm[1] << std::endl;
    //std::cout << "a = " << a[0] << " b = " << a[1] << " d = " << d << std::endl;
    /*
        double d2 = sqrt(pow(xm[0] - xg[0], 2) + pow(xm[1] - xg[1], 2));
        d = d2 * tan(0.5 * DT);*/

    std::cout.precision(14);

    std::vector<double> xi(dim);

    std::vector < std::vector < double > > J(3, std::vector<double>(3));
    //std::vector < std::vector < double > > JI(3, std::vector<double>(3));
    J[0][0] = (-x1 + x2);
    J[0][1] = (-x1 + x3);
    J[0][2] = (-x1 + x4);

    J[1][0] = (-y1 + y2);
    J[1][1] = (-y1 + y3);
    J[1][2] = (-y1 + y4);

    J[2][0] = (-z1 + z2);
    J[2][1] = (-z1 + z3);
    J[2][2] = (-z1 + z4);

    double den =   J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
                   - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
                   + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

    volume = den / 6.;

    xi[0] = -(x3 * y4 * z1 - x3 * xm[1] * z1 - x1 * y4 * z3 + x1 * xm[1] * z3 - x3 * y1 * z4 + x1 * y3 * z4 - x1 * xm[1] * z4 + x3 * xm[1] * z4 +
              xm[0] * (y3 * z1 - y4 * z1 - y1 * z3 + y4 * z3 + y1 * z4 - y3 * z4) +
              x3 * y1 * xm[2] - x1 * y3 * xm[2] + x1 * y4 * xm[2] - x3 * y4 * xm[2] +
              x4 * (xm[1] * z1 + y1 * z3 - xm[1] * z3 - y1 * xm[2] + y3 * (-z1 + xm[2]))) / den;

    xi[1] = -(-(x2 * y4 * z1) + x2 * xm[1] * z1 + x1 * y4 * z2 - x1 * xm[1] * z2 + x2 * y1 * z4 - x1 * y2 * z4 + x1 * xm[1] * z4 - x2 * xm[1] * z4 +
              xm[0] * (-(y2 * z1) + y4 * z1 + y1 * z2 - y4 * z2 - y1 * z4 + y2 * z4) +
              (-(x2 * y1) + x1 * y2 - x1 * y4 + x2 * y4) * xm[2] +
              x4 * (-(xm[1] * z1) - y1 * z2 + xm[1] * z2 + y2 * (z1 - xm[2]) + y1 * xm[2])) / den;


    xi[2] = -(x2 * y3 * z1 - x2 * xm[1] * z1 - x1 * y3 * z2 + x1 * xm[1] * z2 - x2 * y1 * z3 + x1 * y2 * z3 - x1 * xm[1] * z3 + x2 * xm[1] * z3 +
              xm[0] * (y2 * z1 - y3 * z1 - y1 * z2 + y3 * z2 + y1 * z3 - y2 * z3) +
              x2 * y1 * xm[2] - x1 * y2 * xm[2] + x1 * y3 * xm[2] - x2 * y3 * xm[2] +
              x3 * (xm[1] * z1 + y1 * z2 - xm[1] * z2 - y1 * xm[2] + y2 * (-z1 + xm[2]))) / den;


    // std::cout << "AAA " << xi[0] << " " << xi[1] << " " << xi[2] << " " << den << std::endl;


    a2.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        a2[k] -= J[j][k] * a[j]; // this normal has to point toward the center of the ball, thus -=
      }
    }
    double bNorm = sqrt(a2[0] * a2[0] + a2[1] * a2[1] + a2[2] * a2[2]);
    a2[0] /= bNorm;
    a2[1] /= bNorm;
    a2[2] /= bNorm;
    d2 = - a2[0] * xi[0] - a2[1] * xi[1] - a2[2] * xi[2];


    // std::cout << a2[0] << " " << a2[1] << " " << a2[2] << " " << d2 << " " << std::endl;
  }
}

double getHeightPolyhedronSphereInt(const std::vector < std::vector <double> > &b, const std::vector <double> &a, const std::vector <double> &xg, const double &R) {
  const unsigned& cnt = b.size();
  if(b.size() < 3) {
    abort();
  }
  const unsigned& dim = b[0].size();

  std::vector < std::vector <double> > v(cnt, std::vector<double>(dim, 0.));
  for(unsigned i = 0; i < cnt; i++) {
    unsigned ip1 = (i + 1) % cnt;
    for(unsigned k = 0; k < dim; k++) {
      v[i][k] += b[ip1][k] - b[i][k];
    }
  }


  double S = - M_PI * (cnt - 2u);
  for(unsigned i = 0; i < cnt; i++) {
    double dotf = 0.;
    double dotb = 0.;
    unsigned im1 = (cnt + i - 1u) % cnt;
    for(unsigned k = 0; k < dim; k++) {
      dotf += v[i][k] * b[i][k];
      dotb += v[im1][k] * b[i][k];
    }
    double PfdotPb = 0.;
    double normPf = 0.;
    double normPb = 0.;
    for(unsigned k = 0; k < dim; k++) {
      double pf = v[i][k] - dotf * b[i][k];
      double pb = - v[im1][k] + dotb * b[i][k];
      PfdotPb += pf * pb;
      normPf += pf * pf;
      normPb += pb * pb;
    }
    normPf = sqrt(normPf);
    normPb = sqrt(normPb);
    S += acos(PfdotPb / (normPf * normPb));

  }

  std::vector < std::vector <double> > x(cnt, xg);

  for(unsigned i = 0; i < cnt; i++) {
    double h = 0.;
    for(unsigned k = 0; k < dim; k++) {
      h += b[i][k] * a[k];
    }
    h = 1. / h;
    for(unsigned k = 0; k < dim; k++) {
      x[i][k] += h * b[i][k];
    }
  }
  for(unsigned i = 1; i < cnt; i++) {
    for(unsigned k = 0; k < dim; k++) {
      x[i][k] -= x[0][k];
    }
  }
  x[0] = {0., 0., 0.};

  double A = 0.;
  for(unsigned i = 1; i < cnt - 1; i++) {
    A += 0.5 * sqrt((x[i][1] * x[i + 1][2] - x[i][2] * x[i + 1][1]) * (x[i][1] * x[i + 1][2] - x[i][2] * x[i + 1][1]) +
                    (x[i][2] * x[i + 1][0] - x[i][0] * x[i + 1][2]) * (x[i][2] * x[i + 1][0] - x[i][0] * x[i + 1][2]) +
                    (x[i][0] * x[i + 1][1] - x[i][1] * x[i + 1][0]) * (x[i][0] * x[i + 1][1] - x[i][1] * x[i + 1][0]));
  }

  return R * pow(S / A, 1. / 3.);
}


















