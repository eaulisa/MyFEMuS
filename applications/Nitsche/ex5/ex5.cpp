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



void GetChildernQuadElements(const std::vector < std::vector <double>> &xv,
                             const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &_PMatrix,
                             std::vector<std::vector < std::vector <double>>> &xvChildren);

double GetArea(const unsigned &level, const unsigned &LevelMax, const std::vector<std::vector<double>>&xv,
               const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &_PMatrix);

void BuildPMat(std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &_PMatrix);

double radius = 0.5001;
bool printMesh = true;
std::ofstream fout;
const elem_type *finiteElement;
std::vector <double > dist;
std::vector <double > phi;
std::vector <double > phi_x;
double weight;

int main(int argc, char** args) {
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  unsigned dim = 2;
  std::vector <double> xg(dim);

  xg[0] = 0.5;
  xg[1] = 0.25;

  std::vector < std::vector <double> > xv(dim);

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


  finiteElement = new const elem_type_2D("quad", "biquadratic", "seventh");
  std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > PMatrix;
  BuildPMat(PMatrix);
  if(printMesh) {
    fout.open("mesh.txt");
    fout.close();
  }

  double Area;
  double AnalyticArea =  M_PI * radius * radius ;
  Area = GetArea(0, 0, xv, PMatrix);
  std::cout << "Computed Area level 0 = " << Area << " Analytic Area = " << AnalyticArea << std::endl;
  double Error0 = fabs(Area - AnalyticArea) / AnalyticArea;
  std::cout << "Relative Error level 0 = " << Error0 << std::endl;

  for(unsigned l = 1; l < 12; l++) {
    Area = GetArea(0, l, xv, PMatrix);
    std::cout << "Computed Area level " << l << " = " << Area << " Analytic Area = " << AnalyticArea << std::endl;

    double Errorl = fabs(Area - AnalyticArea) / (M_PI * AnalyticArea);
    std::cout << "Order of Convergence = " << log(Error0 / Errorl) / log(pow(2., l)) << std::endl;
    std::cout << "Relative Error level " << l + 1 << " = " << Errorl << std::endl;

  }
  delete finiteElement;
  return 1;
}



double GetArea(const unsigned &level, const unsigned &levelMax, const std::vector<std::vector<double>>&xv,
               const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &PMatrix) {

  double area = 0.;

  if(printMesh) {
    fout.open("mesh.txt", std::ios::app);
    for(unsigned j = 0; j < 4; j++) {
      fout << level << " " << j << " " << xv[0][j] << " " << xv[1][j] << std::endl;
    }
    fout << level << " " << 0 << " " << xv[0][0] << " " << xv[1][0] << std::endl;
    fout << std::endl;
    fout.close();
  }
  const unsigned _numberOfNodes = 9;
  if(level < levelMax) {
    bool oneNodeIsInside = false;
    bool oneNodeIsOutside = false;
    for(unsigned j = 0; j < _numberOfNodes; j++) {
      double rj = sqrt(xv[0][j] * xv[0][j] + xv[1][j] * xv[1][j]);
      if(rj < radius) {
        oneNodeIsInside = true;
        if(!oneNodeIsOutside) { //loop on the left nodes
          for(unsigned jj = j + 1; jj < _numberOfNodes; jj++) {
            double rjj = sqrt(xv[0][jj] * xv[0][jj] + xv[1][jj] * xv[1][jj]);
            if(rjj >= radius) {
              oneNodeIsOutside = true;
              break;
            }
          }
        }
        break;
      }
      else {
        oneNodeIsOutside = true;
      }
    }

    if((oneNodeIsInside * oneNodeIsOutside)) {
      std::vector < std::vector < std::vector <double> > > xvChildren;
      GetChildernQuadElements(xv, PMatrix, xvChildren);
      for(unsigned i = 0; i < xvChildren.size(); i++) {
        area += GetArea(level + 1, levelMax, xvChildren[i], PMatrix);
      }
    }
    else if(!oneNodeIsOutside) {
      for(unsigned ig = 0; ig < finiteElement->GetGaussPointNumber(); ig++) {
        finiteElement->Jacobian(xv, ig, weight, phi, phi_x); 
        area += weight;
      }  
    }
  }
  else { // integration rule for interface elements
    dist.resize(_numberOfNodes, 0.);
    for(unsigned j = 0; j < _numberOfNodes; j++) {
      dist[j] = radius - sqrt(xv[0][j] * xv[0][j] + xv[1][j] * xv[1][j]);
    }
    double distMax = 2. * sqrt(2.) / pow(2., levelMax);
    double C1 = distMax / 50.;
    double C2 = 1. / (0.5 * M_PI + atan(distMax / C1));

  //  double areaE = 0.;

    for(unsigned ig = 0; ig < finiteElement->GetGaussPointNumber(); ig++) {
      finiteElement->Jacobian(xv, ig, weight, phi, phi_x);
      double dist_g = 0.;
      for(unsigned j = 0; j < _numberOfNodes; j++) {
        dist_g += dist[j] * phi[j];
      }
      double xig = 0.5 + C2 * atan(dist_g / C1);

//       std::cout << xig << " " << weight << " " << xig * weight;
//       std::cout << std::endl;
      area += xig * weight;

//      areaE += xig * weight;
    }

    //std::cout << areaE << std::endl;
  }

  return area;
}

const unsigned _dim = 2;
const unsigned _numberOfChildren = 4;
const unsigned _numberOfNodes = 9;

void GetChildernQuadElements(const std::vector < std::vector <double>> &xv,
                             const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &_PMatrix,
                             std::vector<std::vector < std::vector <double>>> &xvChildren) {

  xvChildren.resize(_numberOfChildren);
  for(unsigned i = 0; i < _numberOfChildren; i++) {
    xvChildren[i].resize(_dim);
    for(unsigned k = 0; k < _dim; k++) {
      xvChildren[i][k].resize(_numberOfNodes);
    }
  }

  for(unsigned i = 0; i < _numberOfChildren; i++) {
    for(unsigned j = 0; j < _numberOfNodes; j++) {
      xvChildren[i][0][j] = 0.;
      xvChildren[i][1][j] = 0.;
      for(unsigned jj = 0; jj < _PMatrix[i][j].size(); jj++) {
        xvChildren[i][0][j] += _PMatrix[i][j][jj].second * xv[0][_PMatrix[i][j][jj].first];
        xvChildren[i][1][j] += _PMatrix[i][j][jj].second * xv[1][_PMatrix[i][j][jj].first];
      }
    }
  }
}


void BuildPMat(std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &_PMatrix) {

  _PMatrix.resize(_numberOfChildren);
  for(unsigned i = 0; i < _numberOfChildren; i++) {
    _PMatrix[i].resize(_numberOfNodes);
  }

  const unsigned _ChildBottomLeftVertex[4] = {0, 4, 8, 7};

  const double _xi[9][2] = {{ -1., -1.}, {1., -1.}, {1., 1.}, { -1., 1.},
    {0., -1.}, {1., 0.}, {0., 1,}, { -1., 0.},  {0., 0.}
  };

  std::vector< double > xiChild(_dim);



  std::vector <double> phi;

  for(unsigned i = 0; i < _numberOfChildren; i++) {

    double xi0b = _xi[_ChildBottomLeftVertex[i]][0];
    double xi1b = _xi[_ChildBottomLeftVertex[i]][1];

    for(unsigned j = 0; j < _numberOfNodes; j++) {

      _PMatrix[i][j].resize(_numberOfNodes);

      xiChild[0] = xi0b + 0.5 * (_xi[j][0] + 1.);
      xiChild[1] = xi1b + 0.5 * (_xi[j][1] + 1.);

      unsigned cnt = 0;
      finiteElement->GetPhi(phi, xiChild);
      for(unsigned jj = 0; jj < _numberOfNodes; jj++) {
        if(fabs(phi[jj]) > 1.0e-10) {
          _PMatrix[i][j][cnt].first = jj;
          _PMatrix[i][j][cnt].second = phi[jj];
          cnt++;
        }
      }
      _PMatrix[i][j].resize(cnt);
    }
  }
}
