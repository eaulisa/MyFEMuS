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

#include "RefineElement.hpp"


double GetIntegral(const double &dMax, const unsigned &level, const unsigned &levelMin, const unsigned &levelMax,
                   RefineElement &refineElement, const unsigned iFather = 0);

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



std::vector <double > dist;
const double *phi;
std::vector <double > phi_x;
double weight;

int main(int argc, char** args) {
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  unsigned dim = 2;

  double dMax;

  std::vector < std::vector <double> > xv(dim);

  char geometry[] = "tri";

  if(!strcmp(geometry, "quad")) {

    radius = 1.5;

    for(unsigned k = 0; k < dim; k++) xv[k].resize(9);

    xv[0][0] = -1.;
    xv[1][0] = -1.;

    xv[0][1] =  1.;
    xv[1][1] = -1.;

    xv[0][2] =  1.;
    xv[1][2] =  1.;

    xv[0][3] = -1.;
    xv[1][3] =  1.;

    xv[0][4] =  0.;
    xv[1][4] = -1.;

    xv[0][5] =  1.;
    xv[1][5] =  0.;

    xv[0][6] =  0.;
    xv[1][6] =  1.;

    xv[0][7] = -1.;
    xv[1][7] =  0.;

    xv[0][8] =  0.25;
    xv[1][8] =  -0.25;

    dMax = sqrt(pow(xv[0][2] - xv[0][0], 2) + pow(xv[1][2] - xv[1][0], 2));
  }
  else if(!strcmp(geometry, "tri")) {
    radius = 1.25;
    for(unsigned k = 0; k < dim; k++) xv[k].resize(7);

    xv[0][0] = -1.;
    xv[1][0] = -1.;

    xv[0][1] =  1.;
    xv[1][1] = -1.;

    xv[0][2] =  -1.;
    xv[1][2] =  1.;

    xv[0][3] =  0.;
    xv[1][3] = -1.;

    xv[0][4] =  0.;
    xv[1][4] =  0.;

    xv[0][5] = -1.;
    xv[1][5] =  0.;

    xv[0][6] =  -1. / 3.;
    xv[1][6] =  -1. / 3.;

    dMax = sqrt(pow(xv[0][2] - xv[0][1], 2) + pow(xv[1][2] - xv[1][1], 2));
  }

  dMax *= 0.025;

  std::cout << dMax << std::endl;

  RefineElement refineElement = RefineElement(geometry, "biquadratic", "seventh");
  const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &PMatrix = refineElement.GetProlongationMatrix();


  for(unsigned k = 0; k < dim; k++) xv[k].resize(refineElement.GetNumberOfNodes());
  if(printMesh) {
    fout.open("mesh.txt");
    fout.close();
  }

  unsigned lmin = 0;
  unsigned lmax = 6;
  refineElement.InitElement(xv, lmax);

  std::clock_t c_start = std::clock();

  std::cout.precision(14);

  double integral;
  double analyticIntegral =  M_PI * pow(radius, 4.) / 8.;
  double eps = dMax * pow(0.5, lmin);
  SetConstants(eps);
  integral = GetIntegral(eps, 0, lmin, lmin, refineElement);

  std::cout << "Computed Integral level " << lmin << " = " << integral << " Analytic Integral = " << analyticIntegral << " ";
  double Errorlm1 = fabs(integral - analyticIntegral) / analyticIntegral;
  std::cout << "Relative Error level " << lmin << " = " << Errorlm1 << std::endl;

  for(unsigned l = lmin + 1; l < lmax; l++) {
    eps = dMax * pow(0.5, l);
    SetConstants(eps);
    integral = GetIntegral(eps, 0, lmin, l, refineElement);
    double Errorl = fabs(integral - analyticIntegral) / analyticIntegral;
    std::cout << "Order of Convergence = " << log(Errorlm1 / Errorl) / log(2) << std::endl;
    std::cout << "Computed Integral level " << l << " = " << integral << " Analytic Integral = " << analyticIntegral << " ";
    std::cout << "Relative Error level " << l + 1 << " = " << Errorl << std::endl;
    Errorlm1 = Errorl;
  }

  std::clock_t c_end = std::clock();
  long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

  return 1;
}

double GetIntegral(const double &eps, const unsigned &level,
                   const unsigned &levelMin, const unsigned &levelMax,
                   RefineElement &refineElement, const unsigned ii) {

  double integral = 0.;
  const unsigned numberOfNodes = refineElement.GetNumberOfNodes();
  const unsigned numberOfChildren = refineElement.GetNumberOfChildren();
  const std::vector < std::vector <double> >  &xv = refineElement.GetElement(level, ii);

  if(level < levelMax) {

    if(level < levelMin) {
      refineElement.BuildElementProlongation(level, ii);
      for(unsigned i = 0; i < numberOfChildren; i++) {
        integral += GetIntegral(eps, level + 1, levelMin, levelMax, refineElement, i);
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
          integral += GetIntegral(eps, level + 1, levelMin, levelMax, refineElement, i);
        }
      }
      else if(!oneNodeIsOutside) {
        const elem_type &finiteElement = refineElement.GetFEM();

        for(unsigned ig = 0; ig < finiteElement.GetGaussPointNumber(); ig++) {
          finiteElement.GetGaussQuantities(xv, ig, weight);

          phi = finiteElement.GetPhi(ig);
          std::vector < double> xg(2, 0.);
          for(unsigned j = 0; j < numberOfNodes; j++) {
            xg[0] += xv[0][j] * phi[j];
            xg[1] += xv[1][j] * phi[j];
          }
          double f = GetIntegrand(xg);

          integral += f * weight;
        }
        if(printMesh) PrintElement(xv, refineElement);
      }
    }
  }
  else { // integration rule for interface elements
    const elem_type &finiteElement = refineElement.GetFEM();
    for(unsigned ig = 0; ig < finiteElement.GetGaussPointNumber(); ig++) {
      finiteElement.GetGaussQuantities(xv, ig, weight);

      phi = finiteElement.GetPhi(ig);
      std::vector < double> xg(2, 0.);

      for(unsigned j = 0; j < numberOfNodes; j++) {
        xg[0] += xv[0][j] * phi[j];
        xg[1] += xv[1][j] * phi[j];
      }

      double f = GetIntegrand(xg);
      double dg1 = GetDistance({xg[0], xg[1]});
      double dg2 = dg1 * dg1;
      double dg3 = dg1 * dg2;
      double dg5 = dg3 * dg2;
      double dg7 = dg5 * dg2;
      double dg9 = dg7 * dg2;

      /* Regularized Heaviside Function from
       * Efficient adaptive integration of functions with sharp gradients
       * and cusps in n-dimensional parallelepipeds, 1v1435.2021:viXra {sec 4.1 (1)}
       * https://arxiv.org/abs/1202.5341
       */
      double phi;
      if(dg1 < -eps)
        phi = 0.;
      else if(dg1 > eps) {
        phi = 1.;
      }
      else {
        phi = (a0 + a1 * dg1 + a3 * dg3 + a5 * dg5 + a7 * dg7 + a9 * dg9);
      }
      integral += phi * f * weight;
    }

    if(printMesh) PrintElement(xv, refineElement);
  }

  return integral;
}



