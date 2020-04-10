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


double GetArea(const double &dMax, const unsigned &level, const unsigned &levelMin, const unsigned &levelMax,
               RefineElement &refineElement, const unsigned iFather = 0);

double radius;
double xc = -1;
double yc = -1;

double GetDistance(const double &x, const double &y) {
  return radius - sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
}

bool printMesh = false;
std::ofstream fout;

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


  RefineElement refineElement = RefineElement(geometry, "biquadratic", "seventh");
  const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &PMatrix = refineElement.GetProlongationMatrix();

  
  for(unsigned k = 0; k < dim; k++) xv[k].resize(refineElement.GetNumberOfNodes());


  if(printMesh) {
    fout.open("mesh.txt");
    fout.close();
  }

  unsigned lmin = 0;
  unsigned lmax = 18;
  refineElement.InitElement(xv,lmax);

  std::clock_t c_start = std::clock();

  double Area;
  double AnalyticArea =  M_PI * radius * radius / 4.;
  Area = GetArea(dMax, 0, lmin, lmin, refineElement);
  std::cout << "Computed Area level " << lmin << " = " << Area << " Analytic Area = " << AnalyticArea << std::endl;
  double Error0 = fabs(Area - AnalyticArea) / AnalyticArea;
  std::cout << "Relative Error level " << lmin << " = " << Error0 << std::endl;

  for(unsigned l = lmin + 1; l < lmax; l++) {
    Area = GetArea(dMax, 0, lmin, l, refineElement);
    std::cout << "Computed Area level " << l << " = " << Area << " Analytic Area = " << AnalyticArea << std::endl;

    double Errorl = fabs(Area - AnalyticArea) / (M_PI * AnalyticArea);
    std::cout << "Order of Convergence = " << log(Error0 / Errorl) / log(pow(2., l - lmin)) << std::endl;
    std::cout << "Relative Error level " << l + 1 << " = " << Errorl << std::endl;
  }

  std::clock_t c_end = std::clock();

  long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

  return 1;
}

double GetArea(const double &dMax, const unsigned &level,
               const unsigned &levelMin, const unsigned &levelMax,
               RefineElement &refineElement, const unsigned ii) {

  double area = 0.;

  const unsigned numberOfNodes = refineElement.GetNumberOfNodes();
  const unsigned numberOfChildren = refineElement.GetNumberOfChildren();
  const std::vector < std::vector <double> >  &xv = refineElement.GetElement(level, ii);
  
  if(level < levelMax) {

    if(level < levelMin) {
      refineElement.BuildElementProlongation(level,ii);
      for(unsigned i = 0; i < numberOfChildren; i++) {
        area += GetArea(dMax / 2., level + 1, levelMin, levelMax, refineElement, i);
      }
    }
    else {
      bool oneNodeIsInside = false;
      bool oneNodeIsOutside = false;
      for(unsigned j = 0; j < numberOfNodes; j++) {
        double d = GetDistance(xv[0][j], xv[1][j]);
        if(d > -0.01 * dMax) {
          oneNodeIsInside = true;
          if(!oneNodeIsOutside) { //loop on the left nodes
            for(unsigned jj = j + 1; jj < numberOfNodes; jj++) {
              if(GetDistance(xv[0][jj], xv[1][jj]) < 0.01 * dMax) {
                oneNodeIsOutside = true;
                break;
              }
            }
          }
          break;
        }
        if(d < 0.01 * dMax) {
          oneNodeIsOutside = true;
        }
      }
      if((oneNodeIsInside * oneNodeIsOutside)) {
        refineElement.BuildElementProlongation(level, ii);
        for(unsigned i = 0; i < numberOfChildren; i++) {
          area += GetArea(dMax / 2., level + 1, levelMin, levelMax, refineElement, i);
        }
      }
      else if(!oneNodeIsOutside) {
        
        const elem_type &finiteElement = refineElement.GetFEM();
        for(unsigned ig = 0; ig < finiteElement.GetGaussPointNumber(); ig++) {
          finiteElement.GetGaussQuantities(xv, ig, weight);
          area += weight;
        }
        if(printMesh) {
          fout.open("mesh.txt", std::ios::app);
          for(unsigned j = 0; j < refineElement.GetNumberOfLinearNodes(); j++) {
            fout << level << " " << j << " " << xv[0][j] << " " << xv[1][j] << std::endl;
          }
          fout << level << " " << 0 << " " << xv[0][0] << " " << xv[1][0] << std::endl;
          fout << std::endl;
          fout.close();
        }
      }
    }
  }
  else { // integration rule for interface elements
    dist.resize(numberOfNodes);
    for(unsigned j = 0; j < numberOfNodes; j++) {
      dist[j] = GetDistance(xv[0][j], xv[1][j]);
    }

    double C1 = dMax / 50.;
    double C2 = 1. / (0.5 * M_PI + atan(dMax / C1));

    const elem_type &finiteElement = refineElement.GetFEM();
    for(unsigned ig = 0; ig < finiteElement.GetGaussPointNumber(); ig++) {
      finiteElement.GetGaussQuantities(xv, ig, weight);
      phi = finiteElement.GetPhi(ig);

      double dist_g = 0.;
      for(unsigned j = 0; j < numberOfNodes; j++) {
        dist_g += dist[j] * phi[j];
      }
      double xig = 0.5 + C2 * atan(dist_g / C1);
      area += xig * weight;
    }
    if(printMesh) {
      fout.open("mesh.txt", std::ios::app);
      for(unsigned j = 0; j < refineElement.GetNumberOfLinearNodes(); j++) {
        fout << level << " " << j << " " << xv[0][j] << " " << xv[1][j] << std::endl;
      }
      fout << level << " " << 0 << " " << xv[0][0] << " " << xv[1][0] << std::endl;
      fout << std::endl;
      fout.close();
    }
  }

  return area;
}
