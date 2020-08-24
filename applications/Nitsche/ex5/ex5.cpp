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

double GetIntegral2(const double &eps0, const double &eps, const unsigned &level, const unsigned &levelMin, const unsigned &levelMax,
                    RefineElement &refineElement, const unsigned ii = 0);

double GetIntegral3(const double &eps, const double &tol, const unsigned &level, const unsigned &levelMax, OctTreeElement & element);

double GetIntegral4(const double &eps0, const double &eps, const unsigned &level, const unsigned &levelMin, const unsigned &levelMax,
                    OctTreeElement & element);

void GetIntegral5(const double &eps0, const double &eps, const unsigned &level, const unsigned &levelMin, const unsigned &levelMax,
                  RefineElement &refineElement, const std::vector< std::vector<double>>&xc_i,
                  const std::vector <unsigned> &igrFather, std::vector<double> & integral_i, const unsigned ii = 0);


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


double GetSmoothTestFunction(const double &dg1, const double &eps, const bool &allNodesAreInside) {

  if(allNodesAreInside) return 1.;

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

double GetDistance(const std::vector < double>  &x, const std::vector < double>  &xc) {
  return radius - sqrt((x[0] - xc[0]) * (x[0] - xc[0]) + (x[1] - xc[1]) * (x[1] - xc[1]));
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

  unsigned lmin = 0;
  unsigned lmax = 15;

  RefineElement refineElement = RefineElement(geometry, "quadratic", "first", "seventh", "fifteenth",  "lobatto");
  refineElement.InitElement(xv, lmax);

  OctTreeElement element;
  element.Init(xv, xv, refineElement.GetProlongationMatrix(), &refineElement.GetFEMMedium(), &refineElement.GetFEMFine(), 0);


  std::cout.precision(14);

  double integral;

  //for a given level max of refinement eps is the characteristic length really used for the unit step function: eps = eps0 * 0.5^lmax
  double eps = 0.01 * radius;
  double analyticIntegral =  M_PI * pow(radius - eps, 2.)
                             + 2. * M_PI * (-5. / 11. * pow(eps, 2) + eps * radius);

  SetConstants(eps);

  xc = +0.;
  yc = -0.;

  std::clock_t c_start = std::clock();
  integral = GetIntegral3(eps, M_PI * radius * radius * 1.0e-7, 0, lmax, element);

  std::cout << "Computed Integral level " << lmin << " = " << integral << " Analytic Integral = " << analyticIntegral << " ";
  double Errorl0 = fabs(integral - analyticIntegral) / analyticIntegral;
  std::cout << "Relative Error level " << lmin << " = " << Errorl0 << std::endl;
  std::clock_t c_end = std::clock();
  long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

  element.PrintCounter();

  xc = +0.3473333;
  yc = -0.2333;

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
  element.Clear();


  lmax = 11;

  xc = +0.3473333;
  yc = -0.2333;

  RefineElement refineElement2 = RefineElement(geometry, "quadratic", "first", "ninth", "fifteenth",  "legendre");
  refineElement2.InitElement(xv, lmax);

  OctTreeElement element2;
  element2.Init(xv, xv, refineElement2.GetProlongationMatrix(), &refineElement2.GetFEMCoarse(), &refineElement2.GetFEMFine(), 0);

  std::cout << std::endl << std::endl;


  printMesh = false;

  for(unsigned k = 0; k < dim; k++) xv[k].resize(refineElement.GetNumberOfNodes());
  if(printMesh) {
    fout.open("mesh.txt");
    fout.close();
  }


  c_start = std::clock();




  xc = 0.2345 + 0.22344 ;
  yc = -0.12345;


  eps = std::min(0.1 * radius, 0.25 * dMax * pow(0.5, lmin));
  eps0 = 0.125 * dMax;
  analyticIntegral =  M_PI * pow(radius - eps, 2.)
                      + 2. * M_PI * (-5. / 11. * pow(eps, 2) + eps * radius);
  lmin = 2;
  lmax = 11;
  SetConstants(eps);

  integral = GetIntegral2(eps0, eps, 0, lmin, lmin, refineElement2);
  //integral = GetIntegral4(eps0, eps, 0, lmin, lmin, element2);

  std::cout << "Computed Integral level " << lmin << " = " << integral << " Analytic Integral = " << analyticIntegral << " ";
  Errorl0 = fabs(integral - analyticIntegral) / analyticIntegral;
  std::cout << "Relative Error level " << lmin << " = " << Errorl0 << std::endl;

  double Errorlm1 = Errorl0;
  for(unsigned l = lmin + 1; l < lmax; l++) {
    eps = std::min(0.1 * radius, 0.25 * dMax * pow(0.5, l));
    analyticIntegral =  M_PI * pow(radius - eps, 2.)
                        + 2. * M_PI * (-5. / 11. * pow(eps, 2) + eps * radius);

    SetConstants(eps);
    integral = GetIntegral2(eps0, eps, 0, lmin, l, refineElement2);
    //integral = GetIntegral4(eps0, eps, 0, lmin, l, element2);
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


  printMesh = false;
  if(printMesh) {
    fout.open("mesh.txt");
    fout.close();
  }




  std::cout << std::endl;
  c_start = std::clock();

  eps = std::min(0.1 * radius, 0.25 * dMax * pow(0.5, lmin));
  eps0 = 0.125 * dMax;
  analyticIntegral =  M_PI * pow(radius - eps, 2.)
                      + 2. * M_PI * (-5. / 11. * pow(eps, 2) + eps * radius);

  SetConstants(eps);


  unsigned ng = 10;
  std::vector<std::vector<double> > xc_i(ng);
  std::vector <unsigned> igr(ng);
  for(unsigned i = 0; i < ng; i++) {
    xc_i[i].resize(2);
    xc_i[i][0] = 0.2345 + 0.22344 * cos(2. * M_PI /  ng  * i);
    xc_i[i][1] = -0.12345 + 0.22344 * sin(2. * M_PI /  ng  * i);
    igr[i] = i;
  }
  std::vector<double> integral_i(ng, 0.);

  GetIntegral5(eps0, eps, 0, lmin, lmin, refineElement2, xc_i, igr, integral_i);

  integral = integral_i[0];

  std::vector<double> error0(ng);
  std::vector<double> errorlm1(ng);
  std::vector<double> errorl(ng);

  for(unsigned ig = 0; ig < ng; ig++) {
    std::cout << "Computed Integral level " << lmin << " = " << integral_i[ig] << " Analytic Integral = " << analyticIntegral << " ";
    error0[ig] = fabs(integral_i[ig] - analyticIntegral) / analyticIntegral;
    std::cout << "Relative Error level " << lmin << " = " << error0[ig] << std::endl;
    errorlm1[ig] = error0[ig];
  }
  Errorl0 = error0[0];
  Errorlm1 = Errorl0;
  for(unsigned l = lmin + 1; l < lmax; l++) {
    eps = std::min(0.1 * radius, 0.25 * dMax * pow(0.5, l));
    analyticIntegral =  M_PI * pow(radius - eps, 2.)
                        + 2. * M_PI * (-5. / 11. * pow(eps, 2) + eps * radius);
    SetConstants(eps);
    integral_i.assign(ng, 0.);

    for(unsigned i = 0; i < ng; i++) {
      igr[i] = i;
    }
    GetIntegral5(eps0, eps, 0, lmin, l, refineElement2, xc_i, igr, integral_i);

    for(unsigned ig = 0; ig < ng; ig++) {
      errorl[ig] = fabs(integral_i[ig] - analyticIntegral) / analyticIntegral;
      std::cout << "Order of Convergence1 = " << log(errorlm1[ig] / errorl[ig]) / log(2) << " ";
      std::cout << "Order of Convergence2 = " << log(error0[ig] / errorl[ig]) / log(pow(2, l - lmin)) << std::endl;
    }
    for(unsigned ig = 0; ig < ng; ig++) {
      std::cout << "Computed Integral level " << l << " = " << integral_i[ig] << " Analytic Integral = " << analyticIntegral << " ";
      std::cout << "Relative Error level " << l  << " = " << errorl[ig] << std::endl;
      errorlm1[ig] = errorl[ig];
    }
  }

  c_end = std::clock();
  time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";


  std::cout << std::endl;
  for(unsigned i = 0; i < ng; i++) {
    std::cout << ",\"+\" using (" << xc_i[i][0] << "+" << radius << "*cos(2*pi*$0/400)):(" << xc_i[i][1] << "+" << radius << "*sin(2*pi*$0/400)) w l title \"\"";
    std::cout << ",\"+\" using (" << xc_i[i][0] << "+" << radius + eps << "*cos(2*pi*$0/400)):(" << xc_i[i][1] << "+" << radius + eps << "*sin(2*pi*$0/400)) w l title \"\"";
    std::cout << ",\"+\" using (" << xc_i[i][0] << "+" << radius - eps << "*cos(2*pi*$0/400)):(" << xc_i[i][1] << "+" << radius - eps << "*sin(2*pi*$0/400)) w l title \"\"";
  }
  std::cout << std::endl;



  return 1;
}

void GetIntegral5(const double &eps0, const double &eps, const unsigned &level,
                  const unsigned &levelMin, const unsigned &levelMax,
                  RefineElement &refineElement, const std::vector< std::vector<double>>&xc,
                  const std::vector <unsigned> &igFather, std::vector<double> & integral, const unsigned ii) {

  const unsigned &numberOfNodes = refineElement.GetNumberOfNodes();
  const unsigned &numberOfChildren = refineElement.GetNumberOfChildren();
  const unsigned &dim = refineElement.GetDimension();

  const std::vector < std::vector <double> >  &xv = refineElement.GetNodeCoordinates(level, ii);

  std::vector <unsigned> igi(igFather.size()); //interface
  unsigned igiSize = 0;

  std::vector <unsigned> ig(igFather.size()); //refine or boundary integral
  unsigned igSize = 0;
  unsigned igbSize;
  unsigned igrSize;

  if(level < levelMin) {
    ig = igFather;
    goto refine;
  }
  else  {
    for(unsigned i = 0; i < igFather.size(); i++) {// loop only on the nodes the father asked for refinement
      bool oneNodeIsInside = false;
      bool oneNodeIsOutside = false;
      double d;
      std::vector< double > xv_j(dim, 0.);
      for(unsigned j = 0; j < numberOfNodes; j++) {
        for(unsigned k = 0; k < dim; k++) {
          xv_j[k] = xv[k][j];
        }
        d = GetDistance(xv_j, xc[igFather[i]]);
        if(d > std::max(eps0, eps)) { // check if the node is inside the thick interface
          if(oneNodeIsOutside) {
            ig[igSize] = igFather[i];
            igSize++;
            break;
          }
          oneNodeIsInside = true;
        }
        else if(d < -std::max(eps0, eps)) { // check if the node is outside the thick interface
          oneNodeIsOutside = true;
          if(oneNodeIsInside) {
            ig[igSize] = igFather[i];
            igSize++;
            break;
          }
        }
        else { // the node is within the thick interface
          oneNodeIsOutside = true;
          ig[igSize] = igFather[i];
          igSize++;
          break;
        }
      }
      if(!oneNodeIsOutside) { // the entire element is inside the thick interface
        igi[igiSize] = igFather[i];
        igiSize++;
      }
    }
    ig.resize(igSize);
  }

  if(level == levelMax) {
    igbSize = igSize;
    igrSize = 0;
  }
  else{
    igbSize = 0;
    igrSize = igSize;  
  }
  
  if(igbSize + igiSize) { // at least one of the integrals have to be computed
    const elem_type &finiteElement = (igbSize) ? refineElement.GetFEMFine() : refineElement.GetFEMCoarse();
    std::vector < double> xg(dim);
    double weight;

    std::vector < double> xiFg(dim);
    const double *phiC;
    std::vector < double > phiF(numberOfNodes);

    const std::vector < std::vector <double> >  &xiF = refineElement.GetNodeLocalCoordinates(level, ii);

    for(unsigned jg = 0; jg < finiteElement.GetGaussPointNumber(); jg++) {
      finiteElement.GetGaussQuantities(xv, jg, weight, phiC);
      xg.assign(dim, 0.);
      xiFg.assign(dim, 0.);
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned j = 0; j < numberOfNodes; j++) {
          xg[k] += xv[k][j] * phiC[j];
          xiFg[k] += xiF[k][j] * phiC[j];
        }
      }
      finiteElement.GetPhi(phiF, xiFg);

      for(unsigned i = 0; i < igbSize; i++) {
        integral[ig[i]] += GetSmoothTestFunction(GetDistance(xg , xc[ig[i]]) , eps) * GetIntegrand(xg) * weight;
      }
      for(unsigned i = 0; i < igiSize; i++) {
        integral[igi[i]] += GetIntegrand(xg) * weight;
      }
    }
    if(level == levelMax &&  printMesh) PrintElement(xv, refineElement);
  }

  if(igrSize) { // at least one of the integrals have to be refined
  refine:
    refineElement.BuildElementProlongation(level, ii);
    for(unsigned i = 0; i < numberOfChildren; i++) {
      GetIntegral5(eps0 / 2, eps, level + 1, levelMin, levelMax, refineElement, xc, ig, integral, i);
    }
  }

  return;
}

double GetIntegral2(const double & eps0, const double & eps, const unsigned & level,
                    const unsigned & levelMin, const unsigned & levelMax,
                    RefineElement & refineElement, const unsigned ii) {

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
      integral += GetIntegral2(eps0 / 2, eps, level + 1, levelMin, levelMax, refineElement, i);
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
      if(d > factor * std::max(eps0, eps)) { // check if one node is inside thick interface
        oneNodeIsInside = true;
        if(oneNodeIsOutside)
          if(level == levelMax)  goto integrate;
          else goto refine;
      }
      else if(d < -factor * std::max(eps0, eps)) { // check if one node is outside thick interface
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
    double weight;

    std::vector < double> xiFg(dim);
    const double *phiC;
    std::vector < double > phiF(numberOfNodes);

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
      integral += GetSmoothTestFunction(GetDistance(xg) , eps, !oneNodeIsOutside) * GetIntegrand(xg) * weight;
    }
    if(level == levelMax &&  printMesh) PrintElement(xv, refineElement);
  }
  return integral;
}



double GetIntegral3(const double & eps, const double & tolerance, const unsigned & level, const unsigned & levelMax, OctTreeElement & element) {

  double thisIntegral = 0.;
  double integral[2] = {0., 0.};

  for(unsigned l = 0 ; l < 2; l++) {

    const std::vector < std::vector < double> > &xg = element.GetGaussCoordinates(l);
    const std::vector < std::vector < double > > &phi  = element.GetGaussShapeFunctions(l);;
    const std::vector < double> &weight = element.GetGaussWeights(l);

    for(unsigned ig = 0; ig < weight.size(); ig++) {
      integral[l] += GetSmoothTestFunction(GetDistance(xg[ig]) , eps, false) * GetIntegrand(xg[ig]) * weight[ig];
    }
  }

  if(fabs(integral[1] - integral[0]) > tolerance && level < levelMax - 1) {
    for(unsigned i = 0; i < element.GetNumberOfChildren(); i++) {
      thisIntegral += GetIntegral3(eps, tolerance / 1., level + 1, levelMax, *element.GetElement(std::vector<unsigned> {i}));
    }
  }
  else {
    thisIntegral += integral[1];
  }

  return thisIntegral;
}

double GetIntegral4(const double & eps0, const double & eps, const unsigned & level,
                    const unsigned & levelMin, const unsigned & levelMax, OctTreeElement & element) {

  double integral = 0.;

//   const std::vector < std::vector < double> > &xg1 = element.GetGaussCoordinates(1);
//   const unsigned &numberOfChildren = element.GetNumberOfChildren();

  const unsigned &numberOfChildren = element.GetNumberOfChildren();
//   const unsigned &dim = element.GetDimension();

  const std::vector < std::vector <double> >  &xv = element.GetNodeCoordinates();


  bool oneNodeIsInside = true;
  bool oneNodeIsOutside = true;

  if(level < levelMin) {
  refine:
    for(unsigned i = 0; i < numberOfChildren; i++) {
      integral += GetIntegral4(eps0 / 2., eps, level + 1, levelMin, levelMax, *element.GetElement(std::vector<unsigned> {i}));
    }
  }
  else if(level <= levelMax) {
    oneNodeIsInside = false;
    oneNodeIsOutside = false;
    double factor = (level == levelMax) ? 1. : 1.;
    double d;
    std::vector< double > x3(3, 0.);
    for(unsigned j = 0; j < xv[0].size(); j++) {
      for(unsigned k = 0; k < xv.size(); k++) {
        x3[k] = xv[k][j];
      }
      d = GetDistance({x3[0], x3[1], x3[2]});
      if(d > factor * eps0) { // check if one node is inside thick interface
        oneNodeIsInside = true;
        if(oneNodeIsOutside)
          if(level == levelMax)  goto integrate;
          else goto refine;
      }
      else if(d < -factor * eps0) { // check if one node is outside thick interface
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

    unsigned l = (!oneNodeIsOutside) ? 0 : 1;

    const std::vector < std::vector < double> > &xg = element.GetGaussCoordinates(l);
    const std::vector < std::vector < double > > &phi  = element.GetGaussShapeFunctions(l);;
    const std::vector < double> &weight = element.GetGaussWeights(l);

    for(unsigned ig = 0; ig < weight.size(); ig++) {
      integral += GetSmoothTestFunction(GetDistance(xg[ig]) , eps, !oneNodeIsOutside) * GetIntegrand(xg[ig]) * weight[ig];
    }
  }
  return integral;
}




double GetIntegral(const double & eps, const unsigned & level,
                   const unsigned & levelMin, const unsigned & levelMax,
                   RefineElement & refineElement, const unsigned ii) {

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
      std::vector< double > xv_j(dim, 0.);
      for(unsigned j = 0; j < numberOfNodes; j++) {
        for(unsigned k = 0; k < dim; k++) {
          xv_j[k] = xv[k][j];
        }
        d = GetDistance(xv_j);
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






