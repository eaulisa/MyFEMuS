/** tutorial/Ex1
 * This example shows how to:
 * initialize a femus application;
 * define the multilevel-mesh object mlMsh;
 * read from the file ./input/square.neu the coarse-level mesh and associate it to mlMsh;
 * add in mlMsh uniform refined level-meshes;
 * define the multilevel-solution object mlSol associated to mlMsh;
 * add in mlSol different types of finite element solution variables;
 * initialize the solution varables;
 * define vtk and gmv writer objects associated to mlSol;
 * print vtk and gmv binary-format files in ./output directory.
 **/

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "PetscMatrix.hpp"
#include "slepceps.h"

#include "LinearImplicitSystem.hpp"
#include "Marker.hpp"
#include "Line.hpp"

#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"

using namespace femus;
#include <boost/math/special_functions/sign.hpp>
#include "RefineElement.hpp"


double GetIntegral(const double &dMax, const unsigned &level, const unsigned &levelMin, const unsigned &levelMax,
                   RefineElement &refineElement, const unsigned iFather = 0);

double GetIntegral1(const double &eps, const unsigned &level, const unsigned &levelMin, const unsigned &levelMax,
                    RefineElement &refineElement, const unsigned ii = 0);

double radius;
double xc = -1;
double yc = -1;

double a0;
double a1;
double a3;
double a5;
double a7;
double a9;

void SetConstants(const double &eps) {
  a0 = 0.5; // 128./256.;
  a1 = pow(eps, -1.) * 1.23046875; // 315/256.;
  a3 = -pow(eps, -3.) * 1.640625; //420./256.;
  a5 = pow(eps, -5.) * 1.4765625; // 378./256.;
  a7 = -pow(eps, -7.) * 0.703125; // 180./256.;
  a9 = pow(eps, -9.) * 0.13671875; // 35./256.;
}

double GetDistance(const std::vector < double>  &x) {
  return radius - sqrt((x[0] - xc) * (x[0] - xc) + (x[1] - yc) * (x[1] - yc));
  //return -.5 - x[0] - x[1];
}

double GetIntegrand(const std::vector < double>  &x) {
  return (radius * radius) - ((x[0] - xc) * (x[0] - xc) + (x[1] - yc) * (x[1] - yc));
}

bool printMesh = false;
std::ofstream fout;

void PrintElement(const std::vector < std::vector < double> > &xv, const RefineElement &refineElement) {
  fout.open("mesh.txt", std::ios::app);
  double f;
  for(unsigned j = 0; j < refineElement.GetNumberOfLinearNodes(); j++) {
    f = GetIntegrand({xv[0][j], xv[1][j]});
    fout << xv[0][j] << " " << xv[1][j] << " " << f << std::endl;
  }
  f = GetIntegrand({xv[0][0], xv[1][0]});
  fout << std::endl;
  fout.close();
}

int main(int argc, char** args) {
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  unsigned dim = 2;

  double dMax;

  std::vector < std::vector <double> > xv;
  char geometry[] = "quad";

  if(!strcmp(geometry, "quad")) {
    radius = 1.5;
    xv = {{ -1., 1., 1., -1., 0., 1., 0., -1., 0.}, { -1., -1., 1., 1., -1., 0., 1., 0., 0.}};
    dMax = sqrt(pow(xv[0][2] - xv[0][0], 2) + pow(xv[1][2] - xv[1][0], 2));
  }
  else if(!strcmp(geometry, "tri")) {
    radius = 1.25;
    xv = {{ -1., 1., -1., 0., 0., -1., -1. / 3.}, { -1., -1., 1., -1., 0., 0., -1. / 3.}};
    dMax = sqrt(pow(xv[0][2] - xv[0][1], 2) + pow(xv[1][2] - xv[1][1], 2));
  }

  // dMax is the characteristic length on the coarse element/mesh
  // eps0 is the characteristic half-thickness for the unit step function on the coarse grid

  double eps0 = dMax * 0.025;

  RefineElement refineElement = RefineElement(geometry, "biquadratic", "seventh");
  
  for(unsigned k = 0; k < dim; k++) xv[k].resize(refineElement.GetNumberOfNodes());
  if(printMesh) {
    fout.open("mesh.txt");
    fout.close();
  }

  unsigned lmin = 0;
  unsigned lmax = 10;
  refineElement.InitElement(xv, lmax);

  std::cout.precision(14);

  double integral;
  double analyticIntegral =  M_PI * pow(radius, 4.) / 8.;

  std::clock_t c_start = std::clock();

  //for a given level max of refinement eps is the characteristic length really used for the unit step function: eps = eps0 * 0.5^lmax
  double eps = eps0 * pow(0.5, lmin);
  SetConstants(eps);
  integral = GetIntegral(eps, 0, lmin, lmin, refineElement);

  std::cout << "Computed Integral level " << lmin << " = " << integral << " Analytic Integral = " << analyticIntegral << " ";
  double Errorl0 = fabs(integral - analyticIntegral) / analyticIntegral;
  std::cout << "Relative Error level " << lmin << " = " << Errorl0 << std::endl;

  double Errorlm1 = Errorl0;
  for(unsigned l = lmin + 1; l < lmax; l++) {
    eps = eps0 * pow(0.5, l);
    SetConstants(eps);
    integral = GetIntegral(eps, 0, lmin, l, refineElement);
    double Errorl = fabs(integral - analyticIntegral) / analyticIntegral;
    std::cout << "Order of Convergence1 = " << log(Errorlm1 / Errorl) / log(2) << " ";
    std::cout << "Order of Convergence2 = " << log(Errorl0 / Errorl) / log(pow(2, l - lmin)) << std::endl;
    std::cout << "Computed Integral level " << l << " = " << integral << " Analytic Integral = " << analyticIntegral << " ";
    std::cout << "Relative Error level " << l + 1 << " = " << Errorl << std::endl;
    Errorlm1 = Errorl;
  }

  std::clock_t c_end = std::clock();
  long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

  c_start = std::clock();

  eps = eps0 * pow(0.5, lmin);
  SetConstants(eps);
  integral = GetIntegral1(eps, 0, lmin, lmin, refineElement);

  std::cout << "Computed Integral level " << lmin << " = " << integral << " Analytic Integral = " << analyticIntegral << " ";
  Errorl0 = fabs(integral - analyticIntegral) / analyticIntegral;
  std::cout << "Relative Error level " << lmin << " = " << Errorl0 << std::endl;

  Errorlm1 = Errorl0;
  for(unsigned l = lmin + 1; l < lmax; l++) {
    eps = eps0 * pow(0.5, l);
    SetConstants(eps);
    integral = GetIntegral1(eps, 0, lmin, l, refineElement);
    double Errorl = fabs(integral - analyticIntegral) / analyticIntegral;
    std::cout << "Order of Convergence1 = " << log(Errorlm1 / Errorl) / log(2) << " ";
    std::cout << "Order of Convergence2 = " << log(Errorl0 / Errorl) / log(pow(2, l - lmin)) << std::endl;
    std::cout << "Computed Integral level " << l << " = " << integral << " Analytic Integral = " << analyticIntegral << " ";
    std::cout << "Relative Error level " << l + 1 << " = " << Errorl << std::endl;
    Errorlm1 = Errorl;
  }

  c_end = std::clock();
  time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

  return 1;
}

double GetIntegral(const double &eps, const unsigned &level,
                   const unsigned &levelMin, const unsigned &levelMax,
                   RefineElement &refineElement, const unsigned ii) {

  double integral = 0.;
  const unsigned &numberOfNodes = refineElement.GetNumberOfNodes();
  const unsigned &numberOfChildren = refineElement.GetNumberOfChildren();
  const unsigned &dim = refineElement.GetDimension();

  const std::vector < std::vector <double> >  &xv = refineElement.GetNodeCoordinates(level, ii);

  bool oneNodeIsInside = true;
  bool oneNodeIsOutside = true;
  if(level < levelMax) {
    if(level < levelMin) {
    refine:
      refineElement.BuildElementProlongation(level, ii);
      for(unsigned i = 0; i < numberOfChildren; i++) {
        integral += GetIntegral(eps, level + 1, levelMin, levelMax, refineElement, i);
      }
    }
    else {
      oneNodeIsInside = false;
      oneNodeIsOutside = false;
      double factor = 1.;
      double d;
      std::vector< double > x3(3, 0.);
      for(unsigned j = 0; j < numberOfNodes; j++) {
        for(unsigned k = 0; k < dim; k++) {
          x3[k] = xv[k][j];
        }
        d = GetDistance({x3[0], x3[1], x3[2]});
        if(d > factor * eps) { // check if one node is inside thick interface
          if(oneNodeIsOutside) goto refine;
          oneNodeIsInside = true;
        }
        else if(d < -factor * eps) { // check if one node is outside thick interface
          if(oneNodeIsInside) goto refine;
          oneNodeIsOutside = true;
        }
        else { // node is inside layer
          goto refine;
        }
      }
      if(!oneNodeIsOutside) { // the entire element is inside the thick interface
        goto integrate;
      }
    }
  }
  else { // integration rule for interface elements
  integrate:

    const elem_type &finiteElement = refineElement.GetFEM();
    std::vector < double> xg(dim);
    std::vector < double> xiFg(dim);
    double f;
    double dg1;
    double dg2;
    double weight;
    const double *phiC;
    std::vector < double > phiF(numberOfNodes);
    double U;
    const std::vector < std::vector <double> >  &xiF = refineElement.GetNodeLocalCoordinates(level, ii);
    for(unsigned ig = 0; ig < finiteElement.GetGaussPointNumber(); ig++) {
      finiteElement.GetGaussQuantities(xv, ig, weight, phiC);
      xg.assign(dim, 0.);
      xiFg.assign(dim, 0.);
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned j = 0; j < numberOfNodes; j++) {
          xg[k] += xv[k][j] * phiC[j];
          xiFg[k] += xiF[k][j] * phiC[j];
        }
      }
      finiteElement.GetPhi(phiF, xiFg);
      
      f = GetIntegrand(xg);

      /* Regularized Heaviside Function from
       * Efficient adaptive integration of functions with sharp gradients
       * and cusps in n-dimensional parallelepipeds, sec 4.1 (1)
       * https://arxiv.org/abs/1202.5341
       */
      if(level == levelMax) { // any element at level l = lmax
        dg1 = GetDistance(xg);
        dg2 = dg1 * dg1;
        if(dg1 < -eps)
          U = 0.;
        else if(dg1 > eps) {
          U = 1.;
        }
        else {
          U = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
        }
        integral += U * f * weight;
      }
      else { // interior element at level < lmax
        integral += f * weight;
      }
    }
    if(printMesh) PrintElement(xv, refineElement);
  }

  return integral;
}






























double GetIntegral1(const double & eps, const unsigned & level,
                    const unsigned & levelMin, const unsigned & levelMax,
                    RefineElement & refineElement, const unsigned ii) {

  double integral = 0.;
  const unsigned numberOfNodes = refineElement.GetNumberOfNodes();
  const unsigned numberOfChildren = refineElement.GetNumberOfChildren();
  const unsigned &dim = refineElement.GetDimension();
  const std::vector < std::vector <double> >  &xv = refineElement.GetNodeCoordinates(level, ii);

  if(level < levelMax) {

    if(level < levelMin) {
      refineElement.BuildElementProlongation(level, ii);
      for(unsigned i = 0; i < numberOfChildren; i++) {
        integral += GetIntegral1(eps, level + 1, levelMin, levelMax, refineElement, i);
      }
    }
    else {
      double factor = 1.;
      bool oneNodeIsInside = false;
      bool oneNodeIsOutside = false;
      for(unsigned j = 0; j < numberOfNodes; j++) {
        double d = GetDistance({xv[0][j], xv[1][j]});
        if(d > -factor * eps) {
          oneNodeIsInside = true;
          if(!oneNodeIsOutside) { //loop on the left nodes
            for(unsigned jj = j; jj < numberOfNodes; jj++) {
              if(GetDistance({xv[0][jj], xv[1][jj]}) < factor * eps) {
                oneNodeIsOutside = true;
                break;
              }
            }
          }
        }
        if(d < factor * eps) {
          oneNodeIsOutside = true;
        }
      }
      if((oneNodeIsInside * oneNodeIsOutside)) {
        refineElement.BuildElementProlongation(level, ii);
        for(unsigned i = 0; i < numberOfChildren; i++) {
          integral += GetIntegral1(eps, level + 1, levelMin, levelMax, refineElement, i);
        }
      }
      else if(!oneNodeIsOutside) {
        const elem_type &finiteElement = refineElement.GetFEM();
        double weight;
        const double *phiF;
        std::vector < double > phiC(numberOfNodes);
        const std::vector < std::vector <double> >  &xi = refineElement.GetNodeLocalCoordinates(level, ii);
        for(unsigned ig = 0; ig < finiteElement.GetGaussPointNumber(); ig++) {
          finiteElement.GetGaussQuantities(xv, ig, weight, phiF);

          std::vector < double> xg(dim, 0.);
          std::vector < double> xig(dim, 0.);
          for(unsigned j = 0; j < numberOfNodes; j++) {
            for(unsigned k = 0; k < dim; k++) {
              xg[k] += xv[k][j] * phiF[j];
              xig[k] += xv[k][j] * phiF[j];
            }
          }
          finiteElement.GetPhi(phiC, xig);
          double f = GetIntegrand(xg);

          integral += f * weight;
        }
        if(printMesh) PrintElement(xv, refineElement);
      }
    }
  }
  else { // integration rule for interface elements
    const elem_type &finiteElement = refineElement.GetFEM();
    std::vector < double> xg(3, 0.);
    std::vector < double> xig(3, 0.);
    double f;
    double dg1;
    double dg2;
    double weight;
    const double *phiF;
    std::vector < double > phiC(numberOfNodes);
    double psi;
    const std::vector < std::vector <double> >  &xi = refineElement.GetNodeLocalCoordinates(level, ii);
    for(unsigned ig = 0; ig < finiteElement.GetGaussPointNumber(); ig++) {
      finiteElement.GetGaussQuantities(xv, ig, weight, phiF);

      std::fill(xg.begin(), xg.begin() + dim, 0.);
      std::fill(xig.begin(), xig.begin() + dim, 0.);
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned j = 0; j < numberOfNodes; j++) {
          xg[k] += xv[k][j] * phiF[j];
          xig[k] += xi[k][j] * phiF[j];
        }
      }
      finiteElement.GetPhi(phiC, xig);

      f = GetIntegrand(xg);
      dg1 = GetDistance({xg[0], xg[1], xg[2]});
      dg2 = dg1 * dg1;

      /* Regularized Heaviside Function from
       * Efficient adaptive integration of functions with sharp gradients
       * and cusps in n-dimensional parallelepipeds, 1v1435.2021:viXra {sec 4.1 (1)}
       * https://arxiv.org/abs/1202.5341
       */

      if(dg1 < -eps)
        psi = 0.;
      else if(dg1 > eps) {
        psi = 1.;
      }
      else {
        psi = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
      }
      integral += psi * f * weight;
    }

    if(printMesh) PrintElement(xv, refineElement);
  }

  return integral;
}

