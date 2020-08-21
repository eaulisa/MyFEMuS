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

double GetIntegral2(const double &eps, const unsigned &level, const unsigned &levelMin, const unsigned &levelMax,
                    RefineElement &refineElement, const unsigned ii = 0);

double GetIntegral3(const double &eps, const double &tol, const unsigned &level, const unsigned &levelMax, OctTreeElement & element);

double GetIntegral4(const double &eps, const unsigned &level, const unsigned &levelMin, const unsigned &levelMax, 
                    OctTreeElement & element, const std::vector <unsigned> &igb);


double radius;
double xc = 0.;//+0.3473333;
double yc = 0.;//-0.782333;

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

/* Regularized Heaviside Function from
    * Efficient adaptive integration of functions with sharp gradients
    * and cusps in n-dimensional parallelepipeds, sec 4.1 (1)
    * https://arxiv.org/abs/1202.5341
    */
//std::cout << weight << " ";

double GetSmoothTestFunction(const double &dg1, const double &eps) {
  double U;
  if(dg1 < -eps)
    U = 0.;
  else if(dg1 > eps) {
    U = 1.;
  }
  else {
    double dg2 = dg1 * dg1;
    U = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
  }
  return U;
}

double GetDistance(const std::vector < double>  &x) {
  return radius - sqrt((x[0] - xc) * (x[0] - xc) + (x[1] - yc) * (x[1] - yc));
  //return -.5 - x[0] - x[1];
}

double GetIntegrand(const std::vector < double>  &x) {
  return 1;//(radius * radius) - ((x[0] - xc) * (x[0] - xc) + (x[1] - yc) * (x[1] - yc));
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
  fout << xv[0][0] << " " << xv[1][0] << " " << f << std::endl;
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
    radius = .235;
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

  double eps0 = (dMax < radius) ? dMax : radius;

  RefineElement refineElement = RefineElement(geometry, "quadratic", "first", "ninth", "fifteenth",  "lobatto");

  unsigned n1 = (15 + 3) / 2;
  std::vector <unsigned> igb(n1 * n1);
  unsigned cnt = 0;
  for(unsigned i = 0; i < n1; i++) {
    for(unsigned j = 0; j < n1; j += (i == 0 || i == 8) ? 1 : n1 - 1) {
      igb[cnt] = (i * n1) + j;
      cnt++;
    }
  }
  igb.resize(cnt);
  for(unsigned i = 0; i < igb.size(); i++) {
    std::cout << igb[i] << " ";
  }
  
  
  printMesh = false;

  for(unsigned k = 0; k < dim; k++) xv[k].resize(refineElement.GetNumberOfNodes());
  if(printMesh) {
    fout.open("mesh.txt");
    fout.close();
  }

  unsigned lmin = 0;
  unsigned lmax = 15;
  refineElement.InitElement(xv, lmax);

  OctTreeElement element;
  element.Init(xv, xv, refineElement.GetProlongationMatrix(), &refineElement.GetFEMMedium(), &refineElement.GetFEMFine());


  std::cout.precision(14);

  double integral;

  //for a given level max of refinement eps is the characteristic length really used for the unit step function: eps = eps0 * 0.5^lmax
  double eps = 0.01 * radius;
  double analyticIntegral =  M_PI * pow(radius - eps, 2.)
                             + 2. * M_PI * (-5. / 11. * pow(eps, 2) + eps * radius);

  SetConstants(eps);

  xc = +0.3473333;
  yc = -0.2333;

  std::clock_t c_start = std::clock();
  integral = GetIntegral3(eps, M_PI * radius * radius * 1.0e-7, 0, lmax, element);

  std::cout << "Computed Integral level " << lmin << " = " << integral << " Analytic Integral = " << analyticIntegral << " ";
  double Errorl0 = fabs(integral - analyticIntegral) / analyticIntegral;
  std::cout << "Relative Error level " << lmin << " = " << Errorl0 << std::endl;
  std::clock_t c_end = std::clock();
  long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

  element.PrintCounter();

  xc = +0;
  yc = -0;

  c_start = std::clock();
  integral = GetIntegral3(eps, M_PI * radius * radius * 1.0e-7, 0, lmax, element);

  std::cout << "Computed Integral level " << lmin << " = " << integral << " Analytic Integral = " << analyticIntegral << " ";
  Errorl0 = fabs(integral - analyticIntegral) / analyticIntegral;
  std::cout << "Relative Error level " << lmin << " = " << Errorl0 << std::endl;
  c_end = std::clock();
  time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

  element.PrintCounter();
  element.PrintElement("mesh.txt");

  xc = +0.3473333;
  yc = -0.2333;


  OctTreeElement element2;
  element2.Init(xv, xv, refineElement.GetProlongationMatrix(), &refineElement.GetFEMCoarse(), &refineElement.GetFEMFine());

  std::cout << std::endl << std::endl;

  printMesh = false;

  lmax = 10;

  c_start = std::clock();

  eps = std::min(0.25 * radius, 0.25 * dMax * pow(0.5, lmin));
  analyticIntegral =  M_PI * pow(radius - eps, 2.)
                      + 2. * M_PI * (-5. / 11. * pow(eps, 2) + eps * radius);
  lmin = 2;
  SetConstants(eps);
  integral = GetIntegral4(eps, 0, lmin, lmin, element2,igb);

  std::cout << "Computed Integral level " << lmin << " = " << integral << " Analytic Integral = " << analyticIntegral << " ";
  Errorl0 = fabs(integral - analyticIntegral) / analyticIntegral;
  std::cout << "Relative Error level " << lmin << " = " << Errorl0 << std::endl;

  double Errorlm1 = Errorl0;
  for(unsigned l = lmin + 1; l < lmax; l++) {
    eps = std::min(0.25 * radius, 0.25 * dMax * pow(0.5, l));
    analyticIntegral =  M_PI * pow(radius - eps, 2.)
                        + 2. * M_PI * (-5. / 11. * pow(eps, 2) + eps * radius);

    SetConstants(eps);
    integral = GetIntegral4(eps, 0, lmin, l, element2, igb);
    double Errorl = fabs(integral - analyticIntegral) / analyticIntegral;
    std::cout << "Order of Convergence1 = " << log(Errorlm1 / Errorl) / log(2) << " ";
    std::cout << "Order of Convergence2 = " << log(Errorl0 / Errorl) / log(pow(2, l - lmin)) << std::endl;
    std::cout << "Computed Integral level " << l << " = " << integral << " Analytic Integral = " << analyticIntegral << " ";
    std::cout << "Relative Error level " << l  << " = " << Errorl << std::endl;
    Errorlm1 = Errorl;
  }

  std::cout << ",\"+\" using (" << xc << "+" << radius << "*cos(2*pi*$0/400)):(" << yc << "+" << radius << "*sin(2*pi*$0/400)) w l";
  std::cout << ",\"+\" using (" << xc << "+" << radius + eps << "*cos(2*pi*$0/400)):(" << yc << "+" << radius + eps << "*sin(2*pi*$0/400)) w l";
  std::cout << ",\"+\" using (" << xc << "+" << radius - eps << "*cos(2*pi*$0/400)):(" << yc << "+" << radius - eps << "*sin(2*pi*$0/400)) w l" << std::endl;

  c_end = std::clock();
  time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

  element2.PrintCounter();
  element2.PrintElement("mesh2.txt");

  return 1;
}


double GetIntegral3(const double &eps, const double &tolerance, const unsigned &level, const unsigned &levelMax, OctTreeElement &element) {

  double thisIntegral = 0.;
  double integral[2] = {0., 0.};

  for(unsigned l = 0 ; l < 2; l++) {

    const std::vector < std::vector < double> > &xg = element.GetGaussCoordinates(l);
    const std::vector < std::vector < double > > &phi  = element.GetGaussShapeFunctions(l);;
    const std::vector < double> &weight = element.GetGaussWeights(l);

    for(unsigned ig = 0; ig < weight.size(); ig++) {
      integral[l] += GetSmoothTestFunction(GetDistance(xg[ig]) , eps) * GetIntegrand(xg[ig]) * weight[ig];
    }
  }

  if(fabs(integral[1] - integral[0]) > tolerance && level < levelMax - 1) {
    for(unsigned i = 0; i < element.GetNumberOfChildren(); i++) {
      thisIntegral += GetIntegral3(eps, tolerance / 1.4, level + 1, levelMax, *element.GetElement(std::vector<unsigned> {i}));
    }
  }
  else {
    thisIntegral += integral[1];
  }

  return thisIntegral;
}

double GetIntegral4(const double &eps, const unsigned &level,
                    const unsigned &levelMin, const unsigned &levelMax, OctTreeElement & element, const std::vector <unsigned> &igb) {

  double integral = 0.;

  const std::vector < std::vector < double> > &xg1 = element.GetGaussCoordinates(1);
  const unsigned &numberOfChildren = element.GetNumberOfChildren();

  bool oneNodeIsInside = true;
  bool oneNodeIsOutside = true;

  if(level < levelMin) {
  refine:
    for(unsigned i = 0; i < numberOfChildren; i++) {
      integral += GetIntegral4(eps, level + 1, levelMin, levelMax, *element.GetElement(std::vector<unsigned> {i}), igb);
    }
  }
  else if(level <= levelMax) {
    oneNodeIsInside = false;
    oneNodeIsOutside = false;
    double factor = (level == levelMax) ? 1. : 1.;
    double d;
    std::vector< double > x3(3, 0.);
    for(unsigned j = 0; j < igb.size(); j++) {
      d = GetDistance(xg1[igb[j]]);
      if(d > factor * eps) { // check if one node is inside thick interface
        oneNodeIsInside = true;
        if(oneNodeIsOutside)
          if(level == levelMax)  goto integrate;
          else goto refine;
      }
      else if(d < -factor * eps) { // check if one node is outside thick interface
        oneNodeIsOutside = true;
        if(oneNodeIsInside)
          if(level == levelMax)  goto integrate;
          else goto refine;

      }
      else { // node is inside layer
        oneNodeIsOutside = true;
        if(level == levelMax)  goto integrate;
        else goto refine;
      }
    }
    if(!oneNodeIsOutside) { // the entire element is inside the thick interface
      goto integrate;
    }
  }
  else { // integration rule for interface elements
  integrate:

    unsigned l = (oneNodeIsOutside) ? 1 : 0;

    const std::vector < std::vector < double> > &xg = element.GetGaussCoordinates(l);
    const std::vector < std::vector < double > > &phi  = element.GetGaussShapeFunctions(l);;
    const std::vector < double> &weight = element.GetGaussWeights(l);

    for(unsigned ig = 0; ig < weight.size(); ig++) {
      integral += GetSmoothTestFunction(GetDistance(xg[ig]) , eps) * GetIntegrand(xg[ig]) * weight[ig];
    }
  }
  return integral;
}


double GetIntegral2(const double &eps, const unsigned &level,
                    const unsigned &levelMin, const unsigned &levelMax,
                    RefineElement &refineElement, const unsigned ii) {

  double integral = 0.;
  const unsigned &numberOfNodes = refineElement.GetNumberOfNodes();
  const unsigned &numberOfChildren = refineElement.GetNumberOfChildren();
  const unsigned &dim = refineElement.GetDimension();

  const std::vector < std::vector <double> >  &xv = refineElement.GetNodeCoordinates(level, ii);

  bool oneNodeIsInside = true;
  bool oneNodeIsOutside = true;

  if(level < levelMin) {
  refine:
    refineElement.BuildElementProlongation(level, ii);
    for(unsigned i = 0; i < numberOfChildren; i++) {
      integral += GetIntegral2(eps, level + 1, levelMin, levelMax, refineElement, i);
    }
  }
  else if(level <= levelMax) {
    oneNodeIsInside = false;
    oneNodeIsOutside = false;
    double factor = (level == levelMax) ? 1. : 1.;
    double d;
    std::vector< double > x3(3, 0.);
    for(unsigned j = 0; j < numberOfNodes; j++) {
      for(unsigned k = 0; k < dim; k++) {
        x3[k] = xv[k][j];
      }
      d = GetDistance({x3[0], x3[1], x3[2]});
      if(d > factor * eps) { // check if one node is inside thick interface
        oneNodeIsInside = true;
        if(oneNodeIsOutside)
          if(level == levelMax)  goto integrate;
          else goto refine;
      }
      else if(d < -factor * eps) { // check if one node is outside thick interface
        oneNodeIsOutside = true;
        if(oneNodeIsInside)
          if(level == levelMax)  goto integrate;
          else goto refine;

      }
      else { // node is inside layer
        oneNodeIsOutside = true;
        if(level == levelMax)  goto integrate;
        else goto refine;
      }
    }
    if(!oneNodeIsOutside) { // the entire element is inside the thick interface
      goto integrate;
    }
  }
  else { // integration rule for interface elements
  integrate:

    const elem_type &finiteElement = (oneNodeIsOutside) ? refineElement.GetFEMFine() : refineElement.GetFEMCoarse();
    std::vector < double> xg(dim);
    std::vector < double> xiFg(dim);
    double f;
    double dg1;
    double dg2;
    double weight;
    std::vector<double> phiC;
    std::vector<double> phiCx;
    std::vector < double > phiF(numberOfNodes);
    double U;
    const std::vector < std::vector <double> >  &xiF = refineElement.GetNodeLocalCoordinates(level, ii);


//     for(unsigned i = 0; i < xv[0].size(); i++) {
//       std::cout << xv[0][i] << " " << xv[1][i] << std::endl;
//     }

    //std::cout << finiteElement.GetGaussPointNumber() <<std::endl;
    for(unsigned ig = 0; ig < finiteElement.GetGaussPointNumber(); ig++) {
      //finiteElement.GetGaussQuantities(xv, ig, weight, phiC);
      finiteElement.Jacobian(xv, ig, weight, phiC, phiCx);
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
      //std::cout << weight << " ";
      if(oneNodeIsOutside) {
        //if(level == levelMax) { // any element at level l = lmax
        dg1 = GetDistance(xg);
        //std::cout << dg1 <<" ";
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
        //if(U > 0.) std::cout <<"BBBB "<<ig<<" " << U << " " << integral << " "<<std::endl;
      }
      else { // interior element at level < lmax
        //std::cout << "AAAA";
        integral += f * weight;
      }
      //std::cout<<std::endl;
    }

    double xc = 0.5 * (xv[0][0] + xv[0][2]);
    double yc = 0.5 * (xv[1][0] + xv[1][2]);
//     if(xc > 0.12 && xc < 0.19 && yc > 0.21 && yc < 0.26) {
//       std::cout << level << " " << ii << " " << integral << std::endl;
//
//       std::cout << 0.5 * (xv[0][0] + xv[0][2]) << " " << 0.5 * (xv[1][0] + xv[1][2]) << std::endl;
//     }
    if(level == levelMax &&  printMesh) PrintElement(xv, refineElement);
  }
  //std::cout << level <<"\n";
  return integral;
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

    const elem_type &finiteElement = refineElement.GetFEMFine();
    std::vector < double> xg(dim);
    std::vector < double> xiFg(dim);
    double f;
    double dg1;
    double dg2;
    double weight;
    std::vector<double> phiC;
    std::vector<double> phiCx;
    std::vector < double > phiF(numberOfNodes);
    double U;
    const std::vector < std::vector <double> >  &xiF = refineElement.GetNodeLocalCoordinates(level, ii);


//     for(unsigned i = 0; i < xv[0].size(); i++) {
//       std::cout << xv[0][i] << " " << xv[1][i] << std::endl;
//     }


    for(unsigned ig = 0; ig < finiteElement.GetGaussPointNumber(); ig++) {
      //finiteElement.GetGaussQuantities(xv, ig, weight, phiC);
      finiteElement.Jacobian(xv, ig, weight, phiC, phiCx);
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
        //if(U > 0.) std::cout <<"BBBB "<<ig<<" " << U << " " << integral << " "<<std::endl;
      }
      else { // interior element at level < lmax
        std::cout << "AAAA";
        integral += f * weight;
      }
    }
    if(integral > 0) std::cout << level << " " << ii << " " << integral << std::endl;
    if(printMesh) PrintElement(xv, refineElement);
  }
  //std::cout << level <<"\n";
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
        const elem_type &finiteElement = refineElement.GetFEMFine();
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
    const elem_type &finiteElement = refineElement.GetFEMFine();
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



