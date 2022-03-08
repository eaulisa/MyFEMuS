#pragma once
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "MultiLevelSolution.hpp"


using namespace femus;

class Region {
  private:
    unsigned _size;

    std::vector<std::vector<std::vector<double>>> _x;
    std::vector<std::vector<std::vector<double>>> _minmax;
    std::vector<std::vector<unsigned>> _l2Gmap;
    std::vector<std::vector<double>> _sol;
    std::vector < const elem_type *> _fem;

    std::vector<std::vector<double>> _weight;
    std::vector<std::vector<std::vector<double>>> _xg;
    std::vector<std::vector<double>> _solug;


  public:

    void Reset() {
      _size = 0;
    }

    Region(const unsigned reserveSize = 0) {
      _size = 0;
      _x.reserve(reserveSize);
      _l2Gmap.reserve(reserveSize);
      _sol.reserve(reserveSize);
      _fem.reserve(reserveSize);
      _minmax.reserve(reserveSize);
      _weight.reserve(reserveSize);
      _xg.reserve(reserveSize);
      _solug.reserve(reserveSize);
    }

    void AddElement(const std::vector<std::vector<double>>&x,
                    const std::vector<unsigned> &l2GMap,
                    const std::vector<double> &sol,
                    const elem_type *fem,
                    const std::vector<std::vector<double>>&minmax) {
      _size++;

      _x.resize(_size);
      _l2Gmap.resize(_size);
      _sol.resize(_size);
      _fem.resize(_size);
      _minmax.resize(_size);
      _weight.resize(_size);
      _xg.resize(_size);
      _solug.resize(_size);

      unsigned jel = _size - 1;
      _x[jel] = x;
      _l2Gmap[jel] = l2GMap;
      _sol[jel] = sol;
      _fem[jel] = fem;
      _minmax[jel] = minmax;

      _weight[jel].resize(fem->GetGaussPointNumber());
      _xg[jel].resize(fem->GetGaussPointNumber());
      _solug[jel].resize(fem->GetGaussPointNumber());

      const double *phi;
      const double *phipt;
      std::vector<double>::const_iterator soluit;
      std::vector < std::vector<double> > ::const_iterator xvit;
      std::vector<double>::const_iterator xvkit;
      std::vector<double>::iterator xgit;

      for(unsigned jg = 0; jg < fem->GetGaussPointNumber(); jg++) {

        fem->GetGaussQuantities(x, jg, _weight[jel][jg], phi);

        _solug[jel][jg] = 0.;
        for(soluit = sol.begin(), phipt = phi; soluit != sol.end(); soluit++, phipt++) {
          _solug[jel][jg] += (*soluit) * (*phipt);
        }

        _xg[jel][jg].assign(x.size(), 0.);
        for(xgit = _xg[jel][jg].begin(), xvit = x.begin(); xgit != _xg[jel][jg].end(); xgit++, xvit++) {
          for(xvkit = (*xvit).begin(), phipt = phi;  xvkit != (*xvit).end(); phipt++, xvkit++) {
            *xgit += (*xvkit) * (*phipt);
          }
        }
      }
    }

    const std::vector<std::vector<double>>& GetGaussCoordinates(const unsigned &jel) const {
      return _xg[jel];
    }

    const std::vector<double>& GetGaussWeight(const unsigned &jel) const {
      return _weight[jel];
    }

    const std::vector<double>& GetGaussSolution(const unsigned &jel) const {
      return _solug[jel];
    }

    const std::vector<std::vector<double>>& GetCoordinates(const unsigned &jel) const {
      return _x[jel];
    }

    const std::vector<std::vector<double>>& GetMinMax(const unsigned &jel) const {
      return _minmax[jel];
    }

    const std::vector<unsigned>& GetMapping(const unsigned &jel) const {
      return _l2Gmap[jel];
    }

    const std::vector<double>& GetSolution(const unsigned &jel) const {
      return _sol[jel];
    }

    const elem_type *GetFem(const unsigned &jel) const {
      return _fem[jel];
    }

    const unsigned &size() const {
      return _size;
    }

    const unsigned GetDimension(const unsigned &jel) const {
      return _x[jel].size();
    }

    const unsigned GetDofNumber(const unsigned &jel) const {
      return _l2Gmap[jel].size();
    }

};


#include <boost/math/special_functions/sign.hpp>
#include "RefineElement.hpp"
#include "NonLocal.hpp"

//THIS IS THE 2D ASSEMBLY FOR THE NONLOCAL DIFFUSION PROBLEM with ADAPTIVE QUADRATURE RULE






using namespace femus;

// double GetExactSolutionValue(const std::vector < double >& x) {
//   double pi = acos(-1.);
//   return cos(pi * x[0]) * cos(pi * x[1]);
// };
//
// void GetExactSolutionGradient(const std::vector < double >& x, vector < double >& solGrad) {
//   double pi = acos(-1.);
//   solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
//   solGrad[1] = -pi * cos(pi * x[0]) * sin(pi * x[1]);
// };
//
// double GetExactSolutionLaplace ( const std::vector < double >& x )
// {
//     double pi = acos ( -1. );
//     return -pi * pi * cos ( pi * x[0] ) * cos ( pi * x[1] ) - pi * pi * cos ( pi * x[0] ) * cos ( pi * x[1] );
// };

bool nonLocalAssembly = true;

//DELTA sizes: martaTest1: 0.4, martaTest2: 0.01, martaTest3: 0.53, martaTest4: 0.2, maxTest1: both 0.4, maxTest2: both 0.01, maxTest3: both 0.53, maxTest4: both 0.2, maxTest5: both 0.1, maxTest6: both 0.8,  maxTest7: both 0.05, maxTest8: both 0.025, maxTest9: both 0.0125, maxTest10: both 0.00625

//double delta1 = 0.2; //cubic, quartic, consistency
double delta1 = 0.2; //parallel
double delta2 = 0.2;
// double epsilon = ( delta1 > delta2 ) ? delta1 : delta2;
double kappa1 = 1.;
double kappa2 = 1.;

double A1 = 1. / 16.;
double B1 = - 1. / 8.;
double A2 = 1. / 16.;
double B2 = - 1. / 24.;

// === New variables for adaptive assembly ===
double xc = -1;
double yc = -1;

const unsigned swap[4][9] = {
  {0, 1, 2, 3, 4, 5, 6, 7, 8},
  {3, 0, 1, 2, 7, 4, 5, 6, 8},
  {2, 3, 0, 1, 6, 7, 4, 5, 8},
  {1, 2, 3, 0, 5, 6, 7, 4, 8}
};

const unsigned swapI[4][9] = {
  {0, 1, 2, 3, 4, 5, 6, 7, 8},
  {1, 2, 3, 0, 5, 6, 7, 4, 8},
  {2, 3, 0, 1, 6, 7, 4, 5, 8},
  {3, 0, 1, 2, 7, 4, 5, 6, 8}
};

// === New variables for adaptive assembly ===

void GetBoundaryFunctionValue(double &value, const std::vector < double >& x) {

  //   double u1 = A1 + B1 * x[0] - 1. / (2. * kappa1) * x[0] * x[0] ;
  //   double u2 = A2 + B2 * x[0] - 1. / (2. * kappa2) * x[0] * x[0] ;

  double u1 = (A1 + B1 * x[0] - 1. / (2. * kappa1) * x[0] * x[0]) * (1. + x[0] * x[0]) * cos(x[1]) ;
  double u2 = (A2 + B2 * x[0] - 1. / (2. * kappa2) * x[0] * x[0]) * cos(x[0]) * cos(x[1]);

  value = (x[0] < 0.) ? u1 : u2;

//     value = 0.;
//     value = x[0];
//     value = x[0] * x[0];
//     value = ( x[0] < 0. ) ? x[0] * x[0] * x[0] : 3 * x[0] * x[0] * x[0];
//     value = x[0] * x[0] * x[0] + x[1] * x[1] * x[1];
//     value = x[0] * x[0] * x[0] * x[0] + 0.1 * x[0] * x[0];
//     value = x[0] * x[0] * x[0] * x[0];
//     value =  2 * x[0] + x[0] * x[0] * x[0] * x[0] * x[0]; //this is 2x + x^5


}

unsigned ReorderElement(std::vector < int > &dofs, std::vector < double > & sol, std::vector < std::vector < double > > & x);

void RectangleAndBallRelation(bool &theyIntersect, const std::vector<double> &ballCenter, const double &ballRadius, const std::vector < std::vector < double> > &elementCoordinates,  std::vector < std::vector < double> > &newCoordinates);

void RectangleAndBallRelation2(bool & theyIntersect, const std::vector<double> &ballCenter, const double & ballRadius, const std::vector < std::vector < double> > &elementCoordinates, std::vector < std::vector < double> > &newCoordinates);


//BEGIN New functions: GetIntegral on refined mesh (needs the RefineElement class)

void AssembleNonLocalRefined(MultiLevelProblem& ml_prob) {

  LinearImplicitSystem* mlPdeSys;

  mlPdeSys = &ml_prob.get_system<LinearImplicitSystem> ("NonLocal");

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*            KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  Region region2(10);

  const unsigned  dim = msh->GetDimension();

  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

  unsigned soluIndex;
  unsigned soluPdeIndex;

  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned cntIndex = mlSol->GetIndex("cnt");    // get the position of "u" in the ml_sol object

  std::vector < double >  solu1; // local solution for the nonlocal assembly
  std::vector < double >  solu2; // local solution for the nonlocal assembly

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  std::vector < vector < double > > x1(dim);
  std::vector < vector < double > > x2(dim);

  std::vector< unsigned > l2GMap1; // local to global mapping
  std::vector< unsigned > l2GMap2; // local to global mapping

  std::vector < std::vector < double > > xg1;
  std::vector <double> weight1;
  std::vector <const double *> phi1x;

  std::vector < std::pair<std::vector<double>::iterator, std::vector<double>::iterator> > x1MinMax(dim);
  std::vector < std::pair<std::vector<double>::iterator, std::vector<double>::iterator> > x2MinMax(dim);
  std::vector < std::vector <double> > minmax2;

  RES->zero();
  KK->zero(); // Set to zero all the entries of the Global Matrix

  //BEGIN setup for adaptive integration

  unsigned lmax1 = 4; // consistency form 3 -> 7
  //unsigned lmax1 = 2; // cubic or quartic
  unsigned lmin1 = 1;
  if(lmin1 > lmax1 - 1) lmin1 = lmax1 - 1;


  //consistency
  double dMax = 0.1;
  double eps = 0.125 * dMax *  pow(0.75, lmax1 - 3);

  //cubic
  //double dMax = 0.1 * pow(2./3., level - 1); //marta4, tri unstructured
  //double dMax = 0.1 * pow(2./3., level + 1); //marta4Fine
  //double eps = 0.125 * dMax;

  //quartic
  //double dMax = 0.1 * pow(2./3., level - 1); //marta4, tri unstructured
  //double dMax = 0.1 * pow(2./3., level + 1); //marta4Fine
  //double eps = 0.125 * dMax;

  //parallel
  //double dMax = 0.0125 * pow(1./2., level); //marta4finer
  //double eps = 0.5 * dMax;



  //consistency 3D
//   double dMax = 0.1;
//   double eps = 0.125 * dMax *  pow(0.75, lmax1 - 3);


//   //convergence 3D
//   double dMax = 0.1 * pow(2. / 3., level - 1); //marta4-3D
//   //double dMax = 0.1 * pow(2. / 3., level); //marta4-3D-fine
//   double eps = 0.125 * dMax;

  double areaEl = pow(0.1 * pow(1. / 2., level - 1), dim);

  std::cout << "level = " << level << " ";

  std::cout << "EPS = " << eps << " " << "delta1 = " << delta1 + eps << " " << " lmax1 = " << lmax1 << " lmin1 = " << lmin1 << std::endl;

  RefineElement *refineElement[6][3];
  RefineElement *refineElementCF[6][3];

  NonLocal *nonlocal;

  if(dim == 3) {
    refineElement[0][0] = new RefineElement(lmax1, "hex", "linear", "fifth", "fifth", "legendre");
    refineElement[0][1] = new RefineElement(lmax1, "hex", "quadratic", "fifth", "fifth", "legendre");
    refineElement[0][2] = new RefineElement(lmax1, "hex", "biquadratic", "fifth", "fifth", "legendre");

    refineElement[0][soluType]->SetConstants(eps);

    refineElementCF[0][0] = new RefineElement(lmax1, "hex", "linear", "seventh", "seventh", "legendre");
    refineElementCF[0][1] = new RefineElement(lmax1, "hex", "quadratic", "seventh", "seventh", "legendre");
    refineElementCF[0][2] = new RefineElement(lmax1, "hex", "biquadratic", "seventh", "seventh", "legendre");

    refineElementCF[0][soluType]->SetConstants(eps);


    nonlocal = new NonLocalBall3D();

  }
  else if(dim == 2) {
    refineElement[3][0] = new RefineElement(lmax1, "quad", "linear", "fifth", "fifth", "legendre");
    refineElement[3][1] = new RefineElement(lmax1, "quad", "quadratic", "fifth", "fifth", "legendre");
    refineElement[3][2] = new RefineElement(lmax1, "quad", "biquadratic", "fifth", "fifth", "legendre");

    refineElement[4][0] = new RefineElement(lmax1, "tri", "linear", "fifth", "fifth", "legendre");
    refineElement[4][1] = new RefineElement(lmax1, "tri", "quadratic", "fifth", "fifth", "legendre");
    refineElement[4][2] = new RefineElement(lmax1, "tri", "biquadratic", "fifth", "fifth", "legendre");

    refineElement[3][soluType]->SetConstants(eps);
    refineElement[4][soluType]->SetConstants(eps);


    refineElementCF[3][0] = new RefineElement(lmax1, "quad", "linear", "seventh", "seventh", "legendre");
    refineElementCF[3][1] = new RefineElement(lmax1, "quad", "quadratic", "seventh", "seventh", "legendre");
    refineElementCF[3][2] = new RefineElement(lmax1, "quad", "biquadratic", "seventh", "seventh", "legendre");

    refineElementCF[4][0] = new RefineElement(lmax1, "tri", "linear", "seventh", "seventh", "legendre");
    refineElementCF[4][1] = new RefineElement(lmax1, "tri", "quadratic", "seventh", "seventh", "legendre");
    refineElementCF[4][2] = new RefineElement(lmax1, "tri", "biquadratic", "seventh", "seventh", "legendre");

    refineElementCF[3][soluType]->SetConstants(eps);
    refineElementCF[4][soluType]->SetConstants(eps);




    nonlocal = new NonLocalBall();
    //nonlocal = new NonLocalBox();
  }

  nonlocal->SetKernel(kappa1, delta1, eps);

//   fout.open("mesh.txt");
//   fout.close();

  time_t sSearchTime = 0.;
  time_t sAssemblyTime = 0.;

  std::vector <double> kprocMinMax(nprocs * dim * 2);
  for(unsigned k = 0; k < dim; k++) {
    unsigned kk = iproc * (dim * 2) + k * 2;
    kprocMinMax[kk] = 1.0e10;
    kprocMinMax[kk + 1] = -1.0e10;
  }

  unsigned offset = msh->_elementOffset[iproc];
  unsigned offsetp1 = msh->_elementOffset[iproc + 1];

  for(unsigned iel = offset; iel < offsetp1; iel++) {
    unsigned nDof  = msh->GetElementDofNumber(iel, xType);
    for(unsigned i = 0; i < nDof; i++) {
      unsigned iDof  = msh->GetSolutionDof(i, iel, xType);
      for(unsigned k = 0; k < dim; k++) {
        unsigned kk = iproc * (dim * 2) + k * 2;
        double xk = (*msh->_topology->_Sol[k])(iDof);
        if(xk < kprocMinMax[kk]) kprocMinMax[kk] = xk;
        if(xk > kprocMinMax[kk + 1]) kprocMinMax[kk + 1] = xk;
      }
    }
  }

  for(unsigned kproc = 0; kproc < nprocs; kproc++) {
    MPI_Bcast(&kprocMinMax[kproc * dim * 2], dim * 2, MPI_DOUBLE, kproc, PETSC_COMM_WORLD);
  }
//   for(unsigned kproc = 0; kproc < nprocs; kproc++) {
//     for(unsigned k = 0; k < dim; k++) {
//       unsigned kk = kproc * (dim * 2) + k * dim;
//       //std::cout << kproc << " " << k << " " << kprocMinMax[kk] << " " << kprocMinMax[kk + 1] << "\n";
//     }
//   }

  std::vector < std::vector < unsigned > > orElements(nprocs);
  std::vector < unsigned > orCntSend(nprocs, 0);
  std::vector < unsigned > orCntRecv(nprocs, 0);
  std::vector < unsigned > orSizeSend(nprocs, 0);
  std::vector < unsigned > orSizeRecv(nprocs, 0);
  std::vector < std::vector < unsigned > > orGeomSend(nprocs);
  std::vector < std::vector < unsigned > > orGeomRecv(nprocs);

  std::vector < std::vector < unsigned > > orDofsSend(nprocs);
  std::vector < std::vector < unsigned > > orDofsRecv(nprocs);

  std::vector < std::vector < double > > orSolSend(nprocs);
  std::vector < std::vector < double > > orSolRecv(nprocs);

  std::vector < std::vector < std::vector < double > > > orXSend(nprocs);
  std::vector < std::vector < std::vector < double > > > orXRecv(nprocs);


  unsigned sizeAll = (offsetp1 - offset) * pow(3, dim);
  for(unsigned kproc = 0; kproc < nprocs; kproc++) {
    orElements[kproc].resize(offsetp1 - offset);
    orGeomSend[kproc].resize(offsetp1 - offset);
    orDofsSend[kproc].resize(sizeAll);
    orSolSend[kproc].resize(sizeAll);

    orXSend[kproc].resize(dim);
    orXRecv[kproc].resize(dim);

    for(unsigned k = 0; k < dim; k++) {
      orXSend[kproc][k].resize(sizeAll);
    }
  }

  sol->_Sol[cntIndex]->zero();




  //BEGIN nonlocal assembly
  for(unsigned iel = offset; iel < offsetp1; iel++) {

    std::cout << iel << " " << std::flush;

    unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDof1  = msh->GetElementDofNumber(iel, soluType);

    l2GMap1.resize(nDof1);
    solu1.resize(nDof1);
    for(unsigned k = 0; k < dim; k++) {
      x1[k].resize(nDof1);
    }

    for(unsigned i = 0; i < nDof1; i++) {

      unsigned uDof = msh->GetSolutionDof(i, iel, soluType);
      solu1[i] = (*sol->_Sol[soluIndex])(uDof);

      l2GMap1[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);

      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);
      for(unsigned k = 0; k < dim; k++) {
        x1[k][i] = (*msh->_topology->_Sol[k])(xDof);
      }
    }

    refineElement[ielGeom][soluType]->InitElement1(x1, lmax1);
    refineElementCF[ielGeom][soluType]->InitElement1(x1, lmax1);

    for(unsigned k = 0; k < dim; k++) {
      x1MinMax[k] = std::minmax_element(x1[k].begin(), x1[k].end());
    }


    for(unsigned kproc = 0; kproc < nprocs; kproc++) {
      if(kproc != iproc) {
        bool coarseIntersectionTest = true;
        for(unsigned k = 0; k < dim; k++) {
          unsigned kk = kproc * dim * 2 + k * 2;
          if((*x1MinMax[k].first  - kprocMinMax[kk + 1]) > delta1 + eps  || (kprocMinMax[kk]  - *x1MinMax[k].second) > delta1 + eps) {
            coarseIntersectionTest = false;
            break;
          }
        }
        if(coarseIntersectionTest) {
          orElements[kproc][orCntSend[kproc]] = iel;
          orGeomSend[kproc][orCntSend[kproc]] = ielGeom;
          orCntSend[kproc]++;
          for(unsigned i = 0; i < nDof1; i++) {
            orDofsSend[kproc][orSizeSend[kproc] + i] = l2GMap1[i];
            orSolSend[kproc][orSizeSend[kproc] + i] = solu1[i];
            for(unsigned k = 0; k < dim; k++) {
              orXSend[kproc][k][orSizeSend[kproc] + i] = x1[k][i];
            }
          }
          orSizeSend[kproc] += nDof1;
        }
      }
    }


    //assemble and store RHS
    std::vector <double> res1(nDof1, 0.);
    double weight1;
    const double *phi1;
    const elem_type* fem1 = refineElement[ielGeom][soluType]->GetFem1();
    for(unsigned ig = 0; ig < fem1->GetGaussPointNumber(); ig++) {
      fem1->GetGaussQuantities(x1, ig, weight1, phi1);
      std::vector< double > x1g(dim, 0);
      for(unsigned i = 0; i < nDof1; i++) {
        for(unsigned k = 0; k < dim; k++) {
          x1g[k] += x1[k][i] * phi1[i];
        }
      }
      for(unsigned i = 0; i < nDof1; i++) {

        for(unsigned k = 0; k < dim; k++) {
          res1[i] -= -2 * phi1[i] * weight1; // consistency
//          res1[i] -= -6.* x1g[k] * phi1[i] * weight1; //cubic
//          res1[i] -= ( -12.* x1g[k] * x1g[k] - delta1 * delta1 ) * phi1[i] * weight1; //quartic
        }
      }
    }
    RES->add_vector_blocked(res1, l2GMap1);


    time_t start = clock();
    region2.Reset();
    for(int jel = msh->_elementOffset[iproc]; jel < msh->_elementOffset[iproc + 1]; jel++) {

      short unsigned jelGeom = msh->GetElementType(jel);
      short unsigned jelGroup = msh->GetElementGroup(jel);
      unsigned nDof2  = msh->GetElementDofNumber(jel, soluType);

      for(int k = 0; k < dim; k++) {
        x2[k].resize(nDof2);
      }
      for(unsigned j = 0; j < nDof2; j++) {
        unsigned xDof  = msh->GetSolutionDof(j, jel, xType);
        for(unsigned k = 0; k < dim; k++) {
          x2[k][j] = (*msh->_topology->_Sol[k])(xDof);
        }
      }

      minmax2.resize(dim);
      for(unsigned k = 0; k < dim; k++) {
        minmax2[k].resize(2);
        x2MinMax[k] = std::minmax_element(x2[k].begin(), x2[k].end());
        minmax2[k][0] = *x2MinMax[k].first;
        minmax2[k][1] = *x2MinMax[k].second;
      }
      bool coarseIntersectionTest = true;
      for(unsigned k = 0; k < dim; k++) {
        if((*x1MinMax[k].first  - *x2MinMax[k].second) > delta1 + eps  || (*x2MinMax[k].first  - *x1MinMax[k].second) > delta1 + eps) {
          coarseIntersectionTest = false;
          break;
        }
      }

      if(coarseIntersectionTest) {
        sol->_Sol[cntIndex]->add(jel, 1);
        l2GMap2.resize(nDof2);
        solu2.resize(nDof2);
        for(unsigned j = 0; j < nDof2; j++) {
          unsigned uDof = msh->GetSolutionDof(j, jel, soluType);
          solu2[j] = (*sol->_Sol[soluIndex])(uDof);
          l2GMap2[j] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, j, jel);
        }
        region2.AddElement(x2, l2GMap2, solu2, refineElement[jelGeom][soluType]->GetFem2(), minmax2);
      }
    }

    sSearchTime += clock() - start;

    start = clock();

    nonlocal->ZeroLocalQuantities(nDof1, region2, lmax1);
    bool printMesh = false;

    std::vector<unsigned>jelIndex(region2.size());
    for(unsigned j = 0; j < jelIndex.size(); j++) {
      jelIndex[j] = j;
    }

//     nonlocal->Assembly1(0, lmin1, lmax1, 0, refineElement[ielGeom][soluType]->GetOctTreeElement1(),
//                         *refineElement[ielGeom][soluType], region2, jelIndex,
//                         solu1, kappa1, delta1, printMesh);

    nonlocal->AssemblyCutFem1(0, lmin1, lmax1, 0, refineElement[ielGeom][soluType]->GetOctTreeElement1(), refineElementCF[ielGeom][soluType]->GetOctTreeElement1(),
                              *refineElement[ielGeom][soluType], *refineElementCF[ielGeom][soluType], region2, jelIndex,
                              solu1, kappa1, delta1, printMesh);

    for(unsigned jel = 0; jel < region2.size(); jel++) {

      std::vector<double> & J21 = nonlocal->GetJac21(jel);
      for(unsigned ii = 0; ii < J21.size(); ii++) { // assembly only if one of the entries is different from zero
        if(fabs(J21[ii]) > 1.0e-12 * areaEl) {
          KK->add_matrix_blocked(J21, region2.GetMapping(jel), l2GMap1);
          break;
        }
      }

      KK->add_matrix_blocked(nonlocal->GetJac22(jel), region2.GetMapping(jel), region2.GetMapping(jel));
      RES->add_vector_blocked(nonlocal->GetRes2(jel), region2.GetMapping(jel));
    }
    sAssemblyTime += clock() - start;
  } //end iel loop

  KK->flush();

  std::cout << "[" << iproc << "]  ";
  std::cout << "serial Search Time = " << static_cast<double>(sSearchTime) / CLOCKS_PER_SEC << std::endl;
  std::cout << "[" << iproc << "]  ";
  std::cout << "serial Assembly Time = " << static_cast<double>(sAssemblyTime) / CLOCKS_PER_SEC << std::endl;
  std::cout << std::endl;

  //END serial nonlocal assembly


  //BEGIN parallel nonlocal assembly
  if(nprocs > 1) {

    time_t exchangeTime = clock();

    for(unsigned kproc = 0; kproc < nprocs; kproc++) {
      orElements[kproc].resize(orCntSend[kproc]);
      orGeomSend[kproc].resize(orCntSend[kproc]);
      orDofsSend[kproc].resize(orSizeSend[kproc]);
      orSolSend[kproc].resize(orSizeSend[kproc]);
      for(unsigned k = 0; k < dim; k++) {
        orXSend[kproc][k].resize(orSizeSend[kproc]);
      }
    }

    std::vector < std::vector < MPI_Request > >  reqsSend(nprocs) ;
    std::vector < std::vector < MPI_Request > >  reqsRecv(nprocs) ;
    for(unsigned kproc = 0; kproc < nprocs; kproc++) {
      reqsSend[kproc].resize(3 + dim);
      reqsRecv[kproc].resize(3 + dim);
    }

    for(unsigned kproc = 0; kproc < nprocs; kproc++) {
      MPI_Irecv(&orCntRecv[kproc], 1, MPI_UNSIGNED, kproc, 0, PETSC_COMM_WORLD, &reqsRecv[kproc][0]);
      MPI_Irecv(&orSizeRecv[kproc], 1, MPI_UNSIGNED, kproc, 1, PETSC_COMM_WORLD, &reqsRecv[kproc][1]);
    }

    for(unsigned kproc = 0; kproc < nprocs; kproc++) {
      MPI_Isend(&orCntSend[kproc], 1, MPI_UNSIGNED, kproc, 0, PETSC_COMM_WORLD, &reqsSend[kproc][0]);
      MPI_Isend(&orSizeSend[kproc], 1, MPI_UNSIGNED, kproc, 1, PETSC_COMM_WORLD, &reqsSend[kproc][1]);
    }

    //wait and check that all the sends and receives have been completed successfully
    MPI_Status status;
    for(unsigned kproc = 0; kproc < nprocs; kproc++) {
      for(unsigned m = 0; m < 2; m++) {
        int test = MPI_Wait(&reqsRecv[kproc][m], &status);
        if(test != MPI_SUCCESS) {
          abort();
        }
        test = MPI_Wait(&reqsSend[kproc][m], &status);
        if(test != MPI_SUCCESS) {
          abort();
        }
      }
    }

    for(unsigned kproc = 0; kproc < nprocs; kproc++) {
      orGeomRecv[kproc].resize(orCntRecv[kproc]);
      orDofsRecv[kproc].resize(orSizeRecv[kproc]);
      orSolRecv[kproc].resize(orSizeRecv[kproc]);
      for(unsigned k = 0; k < dim; k++) {
        orXRecv[kproc][k].resize(orSizeRecv[kproc]);
      }
    }

    for(unsigned kproc = 0; kproc < nprocs; kproc++) {
      MPI_Irecv(orGeomRecv[kproc].data(), orGeomRecv[kproc].size(), MPI_UNSIGNED, kproc, 0, PETSC_COMM_WORLD, &reqsRecv[kproc][0]);
      MPI_Irecv(orDofsRecv[kproc].data(), orDofsRecv[kproc].size(), MPI_UNSIGNED, kproc, 1, PETSC_COMM_WORLD, &reqsRecv[kproc][1]);
      MPI_Irecv(orSolRecv[kproc].data(), orSolRecv[kproc].size(), MPI_DOUBLE, kproc, 2, PETSC_COMM_WORLD, &reqsRecv[kproc][2]);
      for(unsigned k = 0; k < dim; k++) {
        MPI_Irecv(orXRecv[kproc][k].data(), orXRecv[kproc][k].size(), MPI_DOUBLE, kproc, 3 + k, PETSC_COMM_WORLD, &reqsRecv[kproc][3 + k]);
      }
    }

    for(unsigned kproc = 0; kproc < nprocs; kproc++) {
      MPI_Isend(orGeomSend[kproc].data(), orGeomSend[kproc].size(), MPI_UNSIGNED, kproc, 0, PETSC_COMM_WORLD, &reqsSend[kproc][0]);
      MPI_Isend(orDofsSend[kproc].data(), orDofsSend[kproc].size(), MPI_UNSIGNED, kproc, 1, PETSC_COMM_WORLD, &reqsSend[kproc][1]);
      MPI_Isend(orSolSend[kproc].data(), orSolSend[kproc].size(), MPI_DOUBLE, kproc, 2, PETSC_COMM_WORLD, &reqsSend[kproc][2]);
      for(unsigned k = 0; k < dim; k++) {
        MPI_Isend(orXSend[kproc][k].data(), orXSend[kproc][k].size(), MPI_DOUBLE, kproc, 3 + k, PETSC_COMM_WORLD, &reqsSend[kproc][3 + k]);
      }
    }

    //wait and check that all the sends and receives have been completed successfully
    for(unsigned kproc = 0; kproc < nprocs; kproc++) {
      for(unsigned m = 0; m < 3 + dim; m++) {
        int test = MPI_Wait(&reqsRecv[kproc][m], &status);
        if(test != MPI_SUCCESS) {
          abort();
        }
        test = MPI_Wait(&reqsSend[kproc][m], &status);
        if(test != MPI_SUCCESS) {
          abort();
        }
      }
    }
    std::cout << "[" << iproc << "]  ";
    std::cout << "Exchange Time = " << static_cast<double>(clock() - exchangeTime) / CLOCKS_PER_SEC << std::endl;
    std::cout << std::endl;

    time_t pSearchTime = 0.;
    time_t pAssemblyTime = 0.;

    std::vector<unsigned > procOrder(nprocs);
    for(unsigned i = 0; i < procOrder.size(); i++) {
      procOrder[i] = i;
    }
    for(unsigned i = 0; i < procOrder.size() - 1; i++) {
      for(unsigned j = i + 1; j < procOrder.size(); j++) {
        if(orGeomRecv[procOrder[i]].size() < orGeomRecv[procOrder[j]].size()) {
          unsigned procOrderi = procOrder[i];
          procOrder[i] = procOrder[j];
          procOrder[j] = procOrderi;
        }
      }
    }

    for(unsigned lproc = 0; lproc < nprocs; lproc++) {
      unsigned kproc = procOrder[lproc];
      if(kproc != iproc) {
        unsigned cnt1 = 0;
        for(unsigned iel = 0; iel < orGeomRecv[kproc].size(); iel++) { // these elements are not own by iproc

          short unsigned ielGeom = orGeomRecv[kproc][iel];
          unsigned nDof1  = el->GetNVE(ielGeom, soluType);

          l2GMap1.resize(nDof1);
          solu1.resize(nDof1);
          for(unsigned k = 0; k < dim; k++) {
            x1[k].resize(nDof1);
          }

          for(unsigned i = 0; i < nDof1; i++) {
            solu1[i] = orSolRecv[kproc][cnt1 + i];
            l2GMap1[i] = orDofsRecv[kproc][cnt1 + i];
            for(unsigned k = 0; k < dim; k++) {
              x1[k][i] =  orXRecv[kproc][k][cnt1 + i];
            }
          }

          refineElement[ielGeom][soluType]->InitElement1(x1, lmax1);
          refineElementCF[ielGeom][soluType]->InitElement1(x1, lmax1);
          for(unsigned k = 0; k < dim; k++) {
            x1MinMax[k] = std::minmax_element(x1[k].begin(), x1[k].end());
          }

          cnt1 += nDof1;
          time_t start = clock();

          region2.Reset();

          unsigned cnt2 = 0;
          for(unsigned jel = 0; jel < orGeomSend[kproc].size(); jel++) { // these elements are own by iproc

            short unsigned jelGeom = orGeomSend[kproc][jel];
            unsigned nDof2  = el->GetNVE(jelGeom, soluType);

            for(unsigned k = 0; k < dim; k++) {
              x2[k].assign(nDof2, 0.);
            }
            for(unsigned j = 0; j < nDof2; j++) {
              for(unsigned k = 0; k < dim; k++) {
                x2[k][j] = orXSend[kproc][k][cnt2 + j];
              }
            }
            minmax2.resize(dim);
            for(unsigned k = 0; k < dim; k++) {
              minmax2[k].resize(2);
              x2MinMax[k] = std::minmax_element(x2[k].begin(), x2[k].end());
              minmax2[k][0] = *x2MinMax[k].first;
              minmax2[k][1] = *x2MinMax[k].second;
            }
            bool coarseIntersectionTest = true;
            for(unsigned k = 0; k < dim; k++) {
              if((*x1MinMax[k].first  - *x2MinMax[k].second) > delta1 + eps  || (*x2MinMax[k].first  - *x1MinMax[k].second) > delta1 + eps) {
                coarseIntersectionTest = false;
                break;
              }
            }

            if(coarseIntersectionTest) {

              sol->_Sol[cntIndex]->add(orElements[kproc][jel], 1);

              l2GMap2.resize(nDof2);
              solu2.resize(nDof2);
              for(unsigned j = 0; j < nDof2; j++) {
                solu2[j] = orSolSend[kproc][cnt2 + j];
                l2GMap2[j] = orDofsSend[kproc][cnt2 + j];
              }
              region2.AddElement(x2, l2GMap2, solu2, refineElement[jelGeom][soluType]->GetFem2(), minmax2);
            }
            cnt2 += nDof2;
          }

          pSearchTime += clock() - start;
          start = clock();

          if(region2.size() > 0) {
            nonlocal->ZeroLocalQuantities(nDof1, region2, lmax1);
            bool printMesh = false;

            std::vector<unsigned>jelIndex(region2.size());
            for(unsigned j = 0; j < jelIndex.size(); j++) {
              jelIndex[j] = j;
            }

//             nonlocal->Assembly1(0, lmin1, lmax1, 0, refineElement[ielGeom][soluType]->GetOctTreeElement1(),
//                                 *refineElement[ielGeom][soluType], region2, jelIndex,
//                                 solu1, kappa1, delta1, printMesh);

            nonlocal->AssemblyCutFem1(0, lmin1, lmax1, 0, refineElement[ielGeom][soluType]->GetOctTreeElement1(), refineElementCF[ielGeom][soluType]->GetOctTreeElement1(),
                                      *refineElement[ielGeom][soluType], *refineElementCF[ielGeom][soluType], region2, jelIndex,
                                      solu1, kappa1, delta1, printMesh);

            for(unsigned jel = 0; jel < region2.size(); jel++) {
              /* The rows of J21, J22 and Res2 are mostly own by iproc, while the columns of J21 and J22 are mostly own by kproc
                 This is okay, since the rows of the global matrix KK and residual RES belong to iproc, and this should optimize
                 the bufferization and exchange of information when closing the KK matrix and the RES vector */

              std::vector<double> & J21 = nonlocal->GetJac21(jel);
              for(unsigned ii = 0; ii < J21.size(); ii++) { // assembly only if one of the entries is different from zero
                if(fabs(J21[ii]) > 1.0e-12 * areaEl) {
                  KK->add_matrix_blocked(J21, region2.GetMapping(jel), l2GMap1);
                  break;
                }
              }

              KK->add_matrix_blocked(nonlocal->GetJac22(jel), region2.GetMapping(jel), region2.GetMapping(jel));
              RES->add_vector_blocked(nonlocal->GetRes2(jel), region2.GetMapping(jel));
            }
          }
          pAssemblyTime += clock() - start;
        }//end iel loop
      }
      KK->flush();
    }

    std::cout << "[" << iproc << "]  ";
    std::cout << "parallel Search Time = " << static_cast<double>(pSearchTime) / CLOCKS_PER_SEC << std::endl;
    std::cout << "[" << iproc << "]  ";
    std::cout << "parallel Assembly Time = " << static_cast<double>(pAssemblyTime) / CLOCKS_PER_SEC << std::endl;
    std::cout << std::endl;

    time_t start = clock();
    RES->close();
    KK->close();
    std::cout << "[" << iproc << "]  ";
    std::cout << "parallel Closing Time = " << static_cast<double>(clock() - start) / CLOCKS_PER_SEC << std::endl;
    std::cout << std::endl;

    std::cout << "[" << iproc << "]  ";
    std::cout << "total Search Time = " << static_cast<double>(sSearchTime + pSearchTime) / CLOCKS_PER_SEC << std::endl;
    std::cout << "[" << iproc << "]  ";
    std::cout << "total Assembly Time = " << static_cast<double>(sAssemblyTime + pAssemblyTime) / CLOCKS_PER_SEC << std::endl;
    std::cout << std::endl;

  }
  else {
    RES->close();
    KK->close();
  }

  sol->_Sol[cntIndex]->close();

  double tolerance = 1.0e-12 * KK->linfty_norm();
  KK->RemoveZeroEntries(tolerance);
//
//   KK->draw();

  delete nonlocal;
  if(dim == 3) {
    delete refineElement[0][0];
    delete refineElement[0][1];
    delete refineElement[0][2];
  }
  else if(dim == 2) {
    delete refineElement[3][0];
    delete refineElement[3][1];
    delete refineElement[3][2];
    delete refineElement[4][0];
    delete refineElement[4][1];
    delete refineElement[4][2];
  }

// ***************** END ASSEMBLY *******************
}



void AssembleNonLocalSys(MultiLevelProblem& ml_prob) {
  adept::Stack& s = FemusInit::_adeptStack;

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("NonLocal");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

  unsigned soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu; // local solution for the local assembly (it uses adept)
  solu.reserve(maxSize);

  vector < double >  solu1; // local solution for the nonlocal assembly
  vector < double >  solu2; // local solution for the nonlocal assembly
  solu1.reserve(maxSize);
  solu2.reserve(maxSize);

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector < double > > x1(dim);
  vector < vector < double > > x2(dim);

  for(unsigned k = 0; k < dim; k++) {
    x1[k].reserve(maxSize);
    x2[k].reserve(maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  phi_xx.reserve(maxSize * dim2);

  vector< adept::adouble > aRes; // local redidual vector
  aRes.reserve(maxSize);

  vector< int > l2GMap1; // local to global mapping
  vector< int > l2GMap2; // local to global mapping
  l2GMap1.reserve(maxSize);
  l2GMap2.reserve(maxSize);

  vector< double > Res1; // local redidual vector
  Res1.reserve(maxSize);
  vector< double > Res2; // local redidual vector
  Res2.reserve(maxSize);

  vector < double > Jac11;
  Jac11.reserve(maxSize * maxSize);
  vector < double > Jac12;
  Jac12.reserve(maxSize * maxSize);

  vector < double > Jac21;
  Jac21.reserve(maxSize * maxSize);
  vector < double > Jac22;
  Jac22.reserve(maxSize * maxSize);


  KK->zero(); // Set to zero all the entries of the Global Matrix

  //BEGIN nonlocal assembly

  for(int kproc = 0; kproc < nprocs; kproc++) {
    for(int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++) {

      short unsigned jelGeom;
      short unsigned jelGroup;
      unsigned nDof2;
      //unsigned nDofx2;

      if(iproc == kproc) {
        jelGeom = msh->GetElementType(jel);
        jelGroup = msh->GetElementGroup(jel);
        nDof2  = msh->GetElementDofNumber(jel, soluType);
      }

      MPI_Bcast(&jelGeom, 1, MPI_UNSIGNED_SHORT, kproc, PETSC_COMM_WORLD);
      MPI_Bcast(&jelGroup, 1, MPI_UNSIGNED_SHORT, kproc, PETSC_COMM_WORLD);
      MPI_Bcast(&nDof2, 1, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);

      l2GMap2.resize(nDof2);
      solu2.resize(nDof2);

      for(int k = 0; k < dim; k++) {
        x2[k].resize(nDof2);
      }

      unsigned typej;
      if(iproc == kproc) {
        for(unsigned j = 0; j < nDof2; j++) {
          l2GMap2[j] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, j, jel);
          unsigned solDof = msh->GetSolutionDof(j, jel, soluType);
          solu2[j] = (*sol->_Sol[soluIndex])(solDof);
          unsigned xDof  = msh->GetSolutionDof(j, jel, xType);

          for(unsigned k = 0; k < dim; k++) {
            x2[k][j] = (*msh->_topology->_Sol[k])(xDof);
          }
        }

        typej = ReorderElement(l2GMap2, solu2, x2);
      }

      MPI_Bcast(&l2GMap2[0], nDof2, MPI_UNSIGNED, kproc, PETSC_COMM_WORLD);
      MPI_Bcast(&solu2[0], nDof2, MPI_DOUBLE, kproc, PETSC_COMM_WORLD);

      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& x2[k][0], nDof2, MPI_DOUBLE, kproc, PETSC_COMM_WORLD);
      }

      for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        bool midpointQuadrature = false;

        short unsigned ielGeom = msh->GetElementType(iel);
        short unsigned ielGroup = msh->GetElementGroup(iel);
        unsigned nDof1  = msh->GetElementDofNumber(iel, soluType);

        l2GMap1.resize(nDof1);
        solu1.resize(nDof1);

        Jac11.assign(nDof1 * nDof1, 0.);
        Jac12.assign(nDof1 * nDof2, 0.);
        Jac21.assign(nDof2 * nDof1, 0.);
        Jac22.assign(nDof2 * nDof2, 0.);
        Res1.assign(nDof1, 0.);
        Res2.assign(nDof2, 0.);

        for(int k = 0; k < dim; k++) {
          x1[k].resize(nDof1);
        }

        for(unsigned i = 0; i < nDof1; i++) {
          l2GMap1[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);
          unsigned solDof = msh->GetSolutionDof(i, iel, soluType);
          solu1[i] = (*sol->_Sol[soluIndex])(solDof);
          unsigned xDof  = msh->GetSolutionDof(i, iel, xType);

          for(unsigned k = 0; k < dim; k++) {
            x1[k][i] = (*msh->_topology->_Sol[k])(xDof);
          }
        }

        unsigned typei = ReorderElement(l2GMap1, solu1, x1);

        double sideLength = fabs(x1[0][0] - x1[0][1]);

        double leftBoundInterface = - sideLength;
        double rightBoundInterface = sideLength;

        unsigned igNumber = (midpointQuadrature) ? 4 : msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber();
        vector < vector < double > > xg1(igNumber);
        vector <double> weight1(igNumber);
        vector < vector <double> > phi1x(igNumber);

        if(midpointQuadrature) {

          for(unsigned ig = 0; ig < igNumber; ig++) {

            std::vector <double> xg1Local(dim);

            weight1[ig] = 0.25 * sideLength * sideLength;

            xg1[ig].assign(dim, 0.);

            unsigned midpointDof = ig + 4;
            unsigned xDof  = msh->GetSolutionDof(midpointDof, iel, xType);

            for(unsigned k = 0; k < dim; k++) {
              xg1[ig][k] = (*msh->_topology->_Sol[k])(xDof);
//                                 std::cout<< xg1[ig][k] << std::endl;
            }

            for(unsigned k = 0; k < dim; k++) {
              xg1Local[k] = - 1. + 2. * (xg1[ig][k] - x1[k][k]) / (x1[k][k + 1] - x1[k][k]);
            }

            double weightTemp;
            msh->_finiteElement[ielGeom][soluType]->Jacobian(x1, xg1Local, weightTemp, phi1x[ig], phi_x);
          }
        }

        else {

          for(unsigned ig = 0; ig < igNumber; ig++) {
            msh->_finiteElement[ielGeom][soluType]->Jacobian(x1, ig, weight1[ig], phi1x[ig], phi_x);

            xg1[ig].assign(dim, 0.);

            for(unsigned i = 0; i < nDof1; i++) {
              for(unsigned k = 0; k < dim; k++) {
                xg1[ig][k] += x1[k][i] * phi1x[ig][i];
              }
            }
          }

        }

        double kernel;
        double radius;

        if(ielGroup == 5 || ielGroup == 7) radius = delta1;       //if x is in Omega_1

        else if(ielGroup == 6 || ielGroup == 8) radius = delta2;       // if x is in Omega_2

        else if(ielGroup == 9 && (jelGroup == 5 || jelGroup == 7 || jelGroup == 9)) radius = delta1;       // use phi_11

        else if(ielGroup == 9 && (jelGroup == 6 || jelGroup == 8)) radius = delta2;    // use phi_22


        bool coarseIntersectionTest = true;

        for(unsigned k = 0; k < dim; k++) {
          double min = 1.0e10;
          min = (min < fabs(x1[k][k]   - x2[k][k])) ?    min :  fabs(x1[k][k]   - x2[k][k]);
          min = (min < fabs(x1[k][k]   - x2[k][k + 1])) ?  min :  fabs(x1[k][k]   - x2[k][k + 1]);
          min = (min < fabs(x1[k][k + 1] - x2[k][k])) ?    min :  fabs(x1[k][k + 1] - x2[k][k]);
          min = (min < fabs(x1[k][k + 1] - x2[k][k + 1])) ?  min :  fabs(x1[k][k + 1] - x2[k][k + 1]);

          if(min >= radius - 1.0e-10) {
            coarseIntersectionTest = false;
            break;
          }
        }

        if(coarseIntersectionTest) {

          bool ifAnyIntersection = false;

          for(unsigned ig = 0; ig < igNumber; ig++) {

            if(iel == jel) {
              for(unsigned i = 0; i < nDof1; i++) {
//                                 Res1[i] -= 0. * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = 0)
                Res1[i] -=  - 2. * weight1[ig]  * phi1x[ig][i]; //Ax - f (so f = - 2)
                //Res1[i] -=  - 6. * xg1[ig][0] * weight1[ig]  * phi1x[ig][i]; //Ax - f (so f = - 6 x)
                /*                if (xg1[ig][0] < 0.) {
                                  double resValue = cos (xg1[ig][1]) * (- 0.5 * xg1[ig][0] * xg1[ig][0] * xg1[ig][0] * xg1[ig][0] - kappa1 / 8. * xg1[ig][0] * xg1[ig][0] * xg1[ig][0] + 11. / 2. * xg1[ig][0] * xg1[ig][0] + kappa1 / 16. * xg1[ig][0] * xg1[ig][0] + kappa1 * 5. / 8. * xg1[ig][0] + 1. - 1. / 16. * kappa1);
                                  Res1[i] -=  resValue * weight1[ig]  * phi1x[ig][i]; //Ax - f (so f = cos(y) * ( - 0.5 * x^4 - kappa1 / 8 * x^3 + 11. / 2. * x^2 + kappa1 / 16. * x^2 + kappa1 * 5. / 8. * x + 1. - 1. / 16. * k1))
                                }
                                else {
                                  double resValue = cos (xg1[ig][1]) * (sin (xg1[ig][0]) * (-kappa2 / 12. - 2 * xg1[ig][0]) + cos (xg1[ig][0]) * (kappa2 / 8. + 1. - kappa2 / 12. * xg1[ig][0] - xg1[ig][0] * xg1[ig][0]));
                                  Res1[i] -=  resValue * weight1[ig]  * phi1x[ig][i];
                                }*///Ax - f (so f = cos(y) * (sin(x) * (-k2 / 12. - 2 * x) + cos(x) * (k2 / 8. + 1. - k2 / 12. * x - x^2)))
//                                 Res1[i] -=  - 6. * xg1[ig][0] * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 6 x)
                // Res1[i] -=  - 6. * ( xg1[ig][0] + xg1[ig][1] ) * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 6 (x + y))
//                                 Res1[i] -= ( - 12. * xg1[ig][0] * xg1[ig][0] - 6. / 5. * radius * radius - 2. * radius ) * weight1[ig] * phi1x[ig][i];  //Ax - f (so f = - 12x^2 - 6/5 * delta^2 - 2 delta)
//                                      Res1[i] -=  - 20. * ( xg1[ig][0] * xg1[ig][0] * xg1[ig][0] ) * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 20 x^3 )
//                                 Res1[i] -=  - 12. * ( xg1[ig][0] * xg1[ig][0] ) * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 12 x^2 )
              }
            }

            std::vector< std::vector < double > > x2New;
            bool theyIntersect;
            RectangleAndBallRelation2(theyIntersect, xg1[ig], radius, x2, x2New);

            if(theyIntersect) {

              ifAnyIntersection = true;

              unsigned jgNumber = msh->_finiteElement[jelGeom][soluType]->GetGaussPointNumber();
//                             unsigned jgNumber = fem->GetGaussPointNumber();

              for(unsigned jg = 0; jg < jgNumber; jg++) {

                vector <double>  phi2y;
                double weight2;

                msh->_finiteElement[jelGeom][soluType]->Jacobian(x2New, jg, weight2, phi2y, phi_x);
//                                 fem->Jacobian ( x2New, jg, weight2, phi2y, phi_x );

                std::vector< double > xg2(dim, 0.);

                for(unsigned j = 0; j < nDof2; j++) {
                  for(unsigned k = 0; k < dim; k++) {
                    xg2[k] += x2New[k][j] * phi2y[j];
                  }
                }

                std::vector <double> xg2Local(dim);

                for(unsigned k = 0; k < dim; k++) {
                  xg2Local[k] = - 1. + 2. * (xg2[k] - x2[k][k]) / (x2[k][k + 1] - x2[k][k]);
                }

                double weightTemp;
                msh->_finiteElement[jelGeom][soluType]->Jacobian(x2, xg2Local, weightTemp, phi2y, phi_x);
//                                 fem->Jacobian ( x2, xg2Local, weightTemp, phi2y, phi_x );

                if((ielGroup == 5 || ielGroup == 7) && (jelGroup == 5 || jelGroup == 7 || jelGroup == 9)) {         //both x and y are in Omega_1
                  kernel = 0.75 * kappa1 / (delta1 * delta1 * delta1 * delta1) ;
                }

                else if((ielGroup == 5 || ielGroup == 7) && (jelGroup == 6 || jelGroup == 8)) {        // x is in Omega_1 and y is in Omega_2
                  kernel = 0.75 * kappa2 / (delta1 * delta1 * delta1 * delta1) ;
                }

                else if((ielGroup == 6 || ielGroup == 8) && (jelGroup == 5 || jelGroup == 7)) {         // x is in Omega_2 and y is in Omega_1
                  kernel = 0.75 * kappa1 / (delta2 * delta2 * delta2 * delta2) ;
                }

                else if((ielGroup == 6 || ielGroup == 8) && (jelGroup == 6 || jelGroup == 8 || jelGroup == 9)) {        // both x and y are in Omega_2
                  kernel = 0.75 * kappa2 / (delta2 * delta2 * delta2 * delta2) ;
                }

                else if(ielGroup == 9 && (jelGroup == 5 || jelGroup == 7 || jelGroup == 9)) {        // use phi_11
                  kernel = 0.75 * kappa1 / (delta1 * delta1 * delta1 * delta1) ;
                }

                else if(ielGroup == 9 && (jelGroup == 6 || jelGroup == 8)) {        // use phi_22
                  kernel = 0.75 * kappa2 / (delta2 * delta2 * delta2 * delta2) ;
                }

                for(unsigned i = 0; i < nDof1; i++) {
                  for(unsigned j = 0; j < nDof1; j++) {
                    double jacValue11 = weight1[ig] * weight2 * kernel * (phi1x[ig][i]) * phi1x[ig][j];
                    Jac11[i * nDof1 + j] -= jacValue11;
                    Res1[i] +=  jacValue11 * solu1[j];
                  }

                  for(unsigned j = 0; j < nDof2; j++) {
                    double jacValue12 = - weight1[ig] * weight2 * kernel * (phi1x[ig][i]) * phi2y[j];
                    Jac12[i * nDof2 + j] -= jacValue12;
                    Res1[i] +=  jacValue12 * solu2[j];
                  }//endl j loop
                }

                for(unsigned i = 0; i < nDof2; i++) {
                  for(unsigned j = 0; j < nDof1; j++) {
                    double jacValue21 = weight1[ig] * weight2 * kernel * (- phi2y[i]) * phi1x[ig][j];
                    Jac21[i * nDof1 + j] -= jacValue21;
                    Res2[i] +=  jacValue21 * solu1[j];
                  }

                  for(unsigned j = 0; j < nDof2; j++) {
                    double jacValue22 = - weight1[ig] * weight2 * kernel * (- phi2y[i]) * phi2y[j];
                    Jac22[i * nDof2 + j] -= jacValue22;
                    Res2[i] +=  jacValue22 * solu2[j];
                  }//endl j loop
                } //endl i loop
              }//end jg loop
            }
          }//end ig loop

//           if(iel == 40 && jel == 28) {
//
//             std::cout.precision(14);
//
//             for(unsigned i = 0; i < nDof1;i++){
//               std::cout << Res1[i] <<" ";
//             }
//             std::cout<<std::endl;
//             for(unsigned i = 0; i < nDof1;i++){
//               std::cout << Res2[i] <<" ";
//             }
//             std::cout<<std::endl;
          /*
                      std::cout << std::endl;

                      for(unsigned i = 0; i < nDof1; i++) {
                        unsigned ii = swapI[typei][i];
                        for(unsigned j = 0; j < nDof1; j++) {
                          unsigned jj = swapI[typei][j];
                          std::cout << Jac11[ii * nDof1 + jj] << " ";
                        }
                        std::cout << std::endl;
                      }
                      std::cout << std::endl;

                      for(unsigned i = 0; i < nDof1; i++) {
                        unsigned ii = swapI[typei][i];
                        for(unsigned j = 0; j < nDof2; j++) {
                          unsigned jj = swapI[typej][j];
                          std::cout << Jac12[ii * nDof2 + jj] << " ";
                        }
                        std::cout << std::endl;
                      }
                      std::cout << std::endl;

                      for(unsigned i = 0; i < nDof2; i++) {
                        unsigned ii = swapI[typej][i];
                        for(unsigned j = 0; j < nDof1; j++) {
                          unsigned jj = swapI[typei][j];
                          std::cout << Jac21[ii * nDof1 + jj] << " ";
                        }
                        std::cout << std::endl;
                      }
                      std::cout << std::endl;

                      for(unsigned i = 0; i < nDof2; i++) {
                        unsigned ii = swapI[typej][i];
                        for(unsigned j = 0; j < nDof2; j++) {
                          unsigned jj = swapI[typej][j];
                          std::cout << Jac22[ii * nDof2 + jj] << " ";
                        }
                        std::cout << std::endl;
                      }
                      std::cout << std::endl;
                    }*/


          if(ifAnyIntersection) {
            KK->add_matrix_blocked(Jac11, l2GMap1, l2GMap1);
            KK->add_matrix_blocked(Jac12, l2GMap1, l2GMap2);
            RES->add_vector_blocked(Res1, l2GMap1);

            KK->add_matrix_blocked(Jac21, l2GMap2, l2GMap1);
            KK->add_matrix_blocked(Jac22, l2GMap2, l2GMap2);
            RES->add_vector_blocked(Res2, l2GMap2);
          }
        }
      } //end iel loop
    } // end jel loop
  } //end kproc loop

  RES->close();

  KK->close();

//     Mat A = ( static_cast<PetscMatrix*> ( KK ) )->mat();
//     MatAssemblyBegin ( A, MAT_FINAL_ASSEMBLY );
//     MatAssemblyEnd ( A, MAT_FINAL_ASSEMBLY );
//     PetscViewer viewer;
//     MatView ( A, viewer );

//     Vec v = ( static_cast< PetscVector* > ( RES ) )->vec();
//     VecView(v,PETSC_VIEWER_STDOUT_WORLD);

  // ***************** END ASSEMBLY *******************
}




void AssembleLocalSys(MultiLevelProblem& ml_prob) {
  adept::Stack& s = FemusInit::_adeptStack;

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Local");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

  unsigned soluIndex = mlSol->GetIndex("u_local");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u_local");    // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu; // local solution for the local assembly (it uses adept)
  solu.reserve(maxSize);

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector < double > > x1(dim);

  for(unsigned k = 0; k < dim; k++) {
    x1[k].reserve(maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  phi_xx.reserve(maxSize * dim2);

  vector< adept::adouble > aRes; // local redidual vector
  aRes.reserve(maxSize);

  vector< int > l2GMap1; // local to global mapping
  l2GMap1.reserve(maxSize);

  vector< double > Res1; // local redidual vector
  Res1.reserve(maxSize);

  vector < double > Jac11;
  Jac11.reserve(maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  //BEGIN local assembly

//   for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
//
//     short unsigned ielGroup = msh->GetElementGroup (iel);
//
//     if (ielGroup == 5 || ielGroup == 6) {   //5 and 6 are the boundary surfaces
//
//       unsigned nDofu  = msh->GetElementDofNumber (iel, soluType);
//       std::vector <double> dofCoordinates (dim);
//
//       for (unsigned i = 0; i < nDofu; i++) {
//         unsigned solDof = msh->GetSolutionDof (i, iel, soluType);
//         unsigned xDof = msh->GetSolutionDof (i, iel, xType);
//         sol->_Bdc[soluIndex]->set (solDof, 0.);
//
//         for (unsigned jdim = 0; jdim < dim; jdim++) {
//           dofCoordinates[jdim] = (*msh->_topology->_Sol[jdim]) (xDof);
//         }
//
//         double bdFunctionValue;
//         GetBoundaryFunctionValue (bdFunctionValue, dofCoordinates);
//         sol->_Sol[soluIndex]->set (solDof, bdFunctionValue);
//
//       }
//
//     }
//
//   }
//
//   sol->_Bdc[soluIndex]->close();
//   sol->_Sol[soluIndex]->close();


  // element loop: each process loops only on the elements that owns
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    l2GMap1.resize(nDofu);
    solu.resize(nDofu);

    for(int i = 0; i < dim; i++) {
      x1[i].resize(nDofx);
    }

    aRes.resize(nDofu);    //resize
    std::fill(aRes.begin(), aRes.end(), 0);    //set aRes to zero

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      l2GMap1[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for(unsigned jdim = 0; jdim < dim; jdim++) {
        x1[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x1, ig, weight, phi, phi_x, boost::none);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point

      vector < adept::adouble > gradSolu_gss(dim, 0.);
      vector < double > x_gss(dim, 0.);

      for(unsigned i = 0; i < nDofu; i++) {
        for(unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          x_gss[jdim] += x1[jdim][i] * phi[i];
        }
      }


//       double aCoeff = 1.;
      double aCoeff = (x_gss[0] < 0) ? kappa1 : kappa2;

      // *** phi_i loop ***
      for(unsigned i = 0; i < nDofu; i++) {

        adept::adouble laplace = 0.;

        for(unsigned jdim = 0; jdim < dim; jdim++) {
          laplace   +=  aCoeff * phi_x[i * dim + jdim] * gradSolu_gss[jdim];
        }


        double srcTerm = 0.;

        for(unsigned k = 0; k < dim; k++) {
          srcTerm +=  -2. ; // so f = - 2 //consistency
//           srcTerm +=  -6. * x_gss[k] ; // cubic
//          srcTerm +=  -12.* x_gss[k] * x_gss[k]; //quartic
        }
        aRes[i] += (-srcTerm * phi[i] + laplace) * weight;

      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    Res1.resize(nDofu);    //resize

    for(int i = 0; i < nDofu; i++) {
      Res1[i] = - aRes[i].value();
    }

    RES->add_vector_blocked(Res1, l2GMap1);

    // define the dependent variables
    s.dependent(&aRes[0], nDofu);

    // define the independent variables
    s.independent(&solu[0], nDofu);

    // get the jacobian matrix (ordered by row major )
    Jac11.resize(nDofu * nDofu);    //resize
    s.jacobian(&Jac11[0], true);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac11, l2GMap1, l2GMap1);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  //END local assembly

  RES->close();

  KK->close();

//     Mat A = ( static_cast<PetscMatrix*> ( KK ) )->mat();
//     MatAssemblyBegin ( A, MAT_FINAL_ASSEMBLY );
//     MatAssemblyEnd ( A, MAT_FINAL_ASSEMBLY );
//     PetscViewer viewer;
//     MatView ( A, viewer );

//     Vec v = ( static_cast< PetscVector* > ( RES ) )->vec();
//     VecView(v,PETSC_VIEWER_STDOUT_WORLD);

  // ***************** END ASSEMBLY *******************
}













void RectangleAndBallRelation(bool & theyIntersect, const std::vector<double> &ballCenter, const double & ballRadius, const std::vector < std::vector < double> > &elementCoordinates, std::vector < std::vector < double> > &newCoordinates) {

  //theyIntersect = true : element and ball intersect
  //theyIntersect = false : element and ball are disjoint

  //elementCoordinates are the coordinates of the vertices of the element

  theyIntersect = false; //by default we assume the two sets are disjoint

  unsigned dim = 2;
  unsigned nDofs = elementCoordinates[0].size();

  std::vector< std::vector < double > > ballVerticesCoordinates(dim);
  newCoordinates.resize(dim);


  for(unsigned n = 0; n < dim; n++) {
    newCoordinates[n].resize(nDofs);
    ballVerticesCoordinates[n].resize(4);

    for(unsigned i = 0; i < nDofs; i++) {
      newCoordinates[n][i] = elementCoordinates[n][i]; //this is just an initalization, it will be overwritten
    }
  }

  double xMinElem = elementCoordinates[0][0];
  double yMinElem = elementCoordinates[1][0];
  double xMaxElem = elementCoordinates[0][2];
  double yMaxElem = elementCoordinates[1][2];


  for(unsigned i = 0; i < 4; i++) {
    if(elementCoordinates[0][i] < xMinElem) xMinElem = elementCoordinates[0][i];

    if(elementCoordinates[0][i] > xMaxElem) xMaxElem = elementCoordinates[0][i];

    if(elementCoordinates[1][i] < yMinElem) yMinElem = elementCoordinates[1][i];

    if(elementCoordinates[1][i] > yMaxElem) yMaxElem = elementCoordinates[1][i];
  }

  //bottom left corner of ball (south west)
  ballVerticesCoordinates[0][0] =  ballCenter[0] - ballRadius;
  ballVerticesCoordinates[1][0] =  ballCenter[1] - ballRadius;

  //top right corner of ball (north east)
  ballVerticesCoordinates[0][2] = ballCenter[0] + ballRadius;
  ballVerticesCoordinates[1][2] = ballCenter[1] + ballRadius;

  newCoordinates[0][0] = (ballVerticesCoordinates[0][0] >= xMinElem) ? ballVerticesCoordinates[0][0] : xMinElem;
  newCoordinates[1][0] = (ballVerticesCoordinates[1][0] >= yMinElem) ? ballVerticesCoordinates[1][0] : yMinElem;

  newCoordinates[0][2] = (ballVerticesCoordinates[0][2] >= xMaxElem) ? xMaxElem : ballVerticesCoordinates[0][2];
  newCoordinates[1][2] = (ballVerticesCoordinates[1][2] >= yMaxElem) ? yMaxElem : ballVerticesCoordinates[1][2];

  if(newCoordinates[0][0] < newCoordinates[0][2] && newCoordinates[1][0] < newCoordinates[1][2]) {    //ball and rectangle intersect

    theyIntersect = true;

    newCoordinates[0][1] = newCoordinates[0][2];
    newCoordinates[1][1] = newCoordinates[1][0];

    newCoordinates[0][3] = newCoordinates[0][0];
    newCoordinates[1][3] = newCoordinates[1][2];

    if(nDofs > 4) {    //TODO the quadratic case has not yet been debugged

      newCoordinates[0][4] = 0.5 * (newCoordinates[0][0] + newCoordinates[0][1]);
      newCoordinates[1][4] = newCoordinates[1][0];

      newCoordinates[0][5] = newCoordinates[0][1];
      newCoordinates[1][5] = 0.5 * (newCoordinates[1][1] + newCoordinates[1][2]);

      newCoordinates[0][6] = newCoordinates[0][4];
      newCoordinates[1][6] = newCoordinates[1][2];

      newCoordinates[0][7] = newCoordinates[0][0];
      newCoordinates[1][7] = newCoordinates[1][5];

      if(nDofs > 8) {

        newCoordinates[0][8] = newCoordinates[0][4];
        newCoordinates[1][8] = newCoordinates[1][5];

      }

    }

  }

}



unsigned ReorderElement(std::vector < int > &dofs, std::vector < double > & sol, std::vector < std::vector < double > > & x) {

  unsigned type = 0;

  if(fabs(x[0][0] - x[0][1]) > 1.e-10) {
    if(x[0][0] - x[0][1] > 0) {
      type = 2;
    }
  }

  else {
    type = 1;

    if(x[1][0] - x[1][1] > 0) {
      type = 3;
    }
  }

  if(type != 0) {
    std::vector < int > dofsCopy = dofs;
    std::vector < double > solCopy = sol;
    std::vector < std::vector < double > > xCopy = x;

    for(unsigned i = 0; i < dofs.size(); i++) {
      dofs[i] = dofsCopy[swap[type][i]];
      sol[i] = solCopy[swap[type][i]];

      for(unsigned k = 0; k < x.size(); k++) {
        x[k][i] = xCopy[k][swap[type][i]];
      }
    }
  }

  return type;
}


void RectangleAndBallRelation2(bool & theyIntersect, const std::vector<double> &ballCenter, const double & ballRadius, const std::vector < std::vector < double> > &elementCoordinates, std::vector < std::vector < double> > &newCoordinates) {

  theyIntersect = false; //by default we assume the two sets are disjoint

  unsigned dim = 2;
  unsigned nDofs = elementCoordinates[0].size();

  newCoordinates.resize(dim);

  for(unsigned i = 0; i < dim; i++) {
    newCoordinates[i].resize(nDofs);
  }

  double xMin = elementCoordinates[0][0];
  double xMax = elementCoordinates[0][2];
  double yMin = elementCoordinates[1][0];
  double yMax = elementCoordinates[1][2];

  if(xMin > xMax || yMin > yMax) {
    std::cout << "error" << std::endl;

    for(unsigned i = 0; i < nDofs; i++) {
      std::cout <<  elementCoordinates[0][i] << " " << elementCoordinates[1][i] << std::endl;
    }

    exit(0);
  }


  double xMinBall = ballCenter[0] - ballRadius;
  double xMaxBall = ballCenter[0] + ballRadius;
  double yMinBall = ballCenter[1] - ballRadius;
  double yMaxBall = ballCenter[1] + ballRadius;


  xMin = (xMin > xMinBall) ? xMin : xMinBall;
  xMax = (xMax < xMaxBall) ? xMax : xMaxBall;
  yMin = (yMin > yMinBall) ? yMin : yMinBall;
  yMax = (yMax < yMaxBall) ? yMax : yMaxBall;

  if(xMin < xMax && yMin < yMax) {    //ball and rectangle intersect

    theyIntersect = true;

    //std::cout<< xMin <<" "<<xMax<<" "<<yMin<<" "<<yMax<<std::endl;

    newCoordinates[0][0] = xMin;
    newCoordinates[0][1] = xMax;
    newCoordinates[0][2] = xMax;
    newCoordinates[0][3] = xMin;

    newCoordinates[1][0] = yMin;
    newCoordinates[1][1] = yMin;
    newCoordinates[1][2] = yMax;
    newCoordinates[1][3] = yMax;

    if(nDofs > 4) {    //TODO the quadratic case has not yet been debugged

      double xMid = 0.5 * (xMin + xMax);
      double yMid = 0.5 * (yMin + yMax);

      newCoordinates[0][4] = xMid;
      newCoordinates[0][5] = xMax;
      newCoordinates[0][6] = xMid;
      newCoordinates[0][7] = xMin;

      newCoordinates[1][4] = yMin;
      newCoordinates[1][5] = yMid;
      newCoordinates[1][6] = yMax;
      newCoordinates[1][7] = yMid;

      if(nDofs > 8) {

        newCoordinates[0][8] = xMid;
        newCoordinates[1][8] = yMid;

      }

    }

  }

}















