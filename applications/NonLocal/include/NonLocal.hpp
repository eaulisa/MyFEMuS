#ifndef __femus_NonLocal_hpp__
#define __femus_NonLocal_hpp__

std::ofstream fout;

class NonLocal {
  public:
    NonLocal() {};
    ~NonLocal() {};
    virtual double GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &size) const = 0;
    virtual double GetKernel(const double  &kappa, const double &delta, const double &eps) const = 0;
    virtual double GetArea(const double &delta, const double &eps) const = 0;

    void ZeroLocalQuantities(const unsigned &nDof1, const unsigned &nDof2);

    void AddFineLevelLocalQuantities(const unsigned &level);

    void Assembly1(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax1, const unsigned &levelMax2,
                   const unsigned &iFather, RefineElement &refineElement1, RefineElement &refineElement2,
                   const vector < double >  &solu1, const vector < double > &solu2,
                   const double &kappa, const double &delta, const bool &ielEqualJel, const bool &printMesh);

    double Assembly2(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax1, const unsigned &levelMax2,
                     const unsigned &iFather, RefineElement &refineElement,
                     const unsigned &nDof1, const vector < double > &xg1, const double &weight1_ig, const vector < double > &phi1_ig,
                     const vector < double >  &solu1, const vector < double > &solu2,
                     const double &kappa, const double &delta, const bool &ielEqualJel, const bool &printMesh);



    double RefinedAssembly(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax, const unsigned &iFather,
                           RefineElement &refineElement,
                           const unsigned &nDof1, const vector < double > &xg1, const double &weight1_ig, const double *phi1_ig,
                           const vector < double >  &solu1, const vector < double > &solu2,
                           const double &kappa, const double &delta, const bool &printMesh, std::vector <double> &value);

    void RefinedAssembly5(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax, const unsigned &iFather,
                          const std::vector <unsigned> &igFather, std::vector <unsigned> &igi, RefineElement &refineElement,
                          const unsigned &nDof1, const std::vector< std::vector<double>>&xc, const std::vector<double> &weight1,
                          const std::vector <const double *> phi1, const vector < double >  &solu1, const vector < double > &solu2,
                          const double &kappa, const double &delta, const bool &printMesh);

    double GetSmoothTestFunction(const double &dg1, const double &eps);

    std::vector < double > & GetRes1() {
      return _res1;
    };
    std::vector < double > & GetRes2() {
      return _res2;
    };
    std::vector < double > & GetJac11() {
      return _jac11;
    };
    std::vector < double > & GetJac12() {
      return _jac12;
    };
    std::vector < double > & GetJac21() {
      return _jac21;
    };
    std::vector < double > & GetJac22() {
      return _jac22;
    };

  private:
    std::vector < double > _res1;
    std::vector < double > _res2;
    std::vector < double > _jac11;
    std::vector < double > _jac12;
    std::vector < double > _jac21;
    std::vector < double > _jac22;

    void PrintElement(const std::vector < std::vector < double> > &xv, const RefineElement &refineElement);

};

void NonLocal::ZeroLocalQuantities(const unsigned &nDof1, const unsigned &nDof2) {

  _jac11.assign(nDof1 * nDof1, 0.);
  _jac12.assign(nDof1 * nDof2, 0.);
  _jac21.assign(nDof2 * nDof1, 0.);
  _jac22.assign(nDof2 * nDof2, 0.);
  _res1.assign(nDof1, 0.);
  _res2.assign(nDof2, 0.);
}



double NonLocal::RefinedAssembly(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax, const unsigned &iFather,
                                 RefineElement &refineElement,
                                 const unsigned &nDof1, const vector < double > &xg1, const double &weight1_ig, const double *phi1_ig,
                                 const vector < double >  &solu1, const vector < double > &solu2,
                                 const double &kappa, const double &delta, const bool &printMesh, std::vector <double> &value) {

  double eps0l = std::max(refineElement.GetEps0() * pow(0.5, level) , refineElement.GetEps());

  double area = 0;
  std::vector <double> lvalue(4, 0.);

  const unsigned &nDof2 = refineElement.GetNumberOfNodes();
  const unsigned &numberOfChildren = refineElement.GetNumberOfChildren();
  const unsigned &dim = refineElement.GetDimension();

  const std::vector < std::vector <double> >  &xv2 = refineElement.GetElement2NodeCoordinates(level, iFather);

  bool oneNodeIsInside = true;
  bool oneNodeIsOutside = true;

  if(level < levelMax - 1) {
    if(level < levelMin) {
    refine:
      refineElement.BuildElement2Prolongation(level, iFather);
      //ZeroLocalQuantities(nDof1, nDof2, level + 1);
      for(unsigned i = 0; i < numberOfChildren; i++) {
        area += RefinedAssembly(level + 1, levelMin, levelMax, i, refineElement,
                                nDof1, xg1, weight1_ig, phi1_ig,
                                solu1, solu2, kappa, delta, printMesh, lvalue);
      }
      //AddFineLevelLocalQuantities(level);
    }
    else {
      oneNodeIsInside = false;
      oneNodeIsOutside = false;
      double d;
      std::vector< double > xv2j(dim);

      for(unsigned j = 0; j < nDof2; j++) {
        for(unsigned k = 0; k < dim; k++) {
          xv2j[k] = xv2[k][j];
        }
        d = this->GetDistance(xg1, xv2j, delta);
        if(d > eps0l) { // check if one node is inside thick interface
          if(oneNodeIsOutside) goto refine;
          oneNodeIsInside = true;
        }
        else if(d < - eps0l) { // check if one node is outside thick interface
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
    std::vector < double> xg2(dim);
    std::vector < double> xi2Fg(dim);
    double weight2;
    const double *phi2y;
    std::vector < double > phi2F(nDof2);
    double U;
    const std::vector < std::vector <double> >  &xi2F = refineElement.GetElement2LocalCoordinates(level, iFather);

    double kernel = this->GetKernel(kappa, delta, refineElement.GetEps());

    for(unsigned jg = 0; jg < finiteElement.GetGaussPointNumber(); jg++) {

      finiteElement.GetGaussQuantities(xv2, jg, weight2, phi2y);
      xg2.assign(dim, 0.);
      xi2Fg.assign(dim, 0.);
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned j = 0; j < nDof2; j++) {
          xg2[k] += xv2[k][j] * phi2y[j];
          xi2Fg[k] += xi2F[k][j] * phi2y[j];
        }
      }
      finiteElement.GetPhi(phi2F, xi2Fg);

      if(level == levelMax - 1) { // only for element at level l = lmax - 1
        U = refineElement.GetSmoothStepFunction(this->GetDistance(xg1, xg2, delta));
      }
      else {
        U = 1.;
      }

      if(U > 0.) {
        area += U * weight2;

        lvalue[0] += U * weight2;
        lvalue[1] += U * weight2 * (xg1[0] - xg2[0]) * (xg1[0] - xg2[0]);
        lvalue[2] += U * weight2 * (xg1[1] - xg2[1]) * (xg1[1] - xg2[1]);
        lvalue[3] += U * weight2 * (xg1[0] - xg2[0]) * (xg1[1] - xg2[1]);

        double C =  U * weight1_ig * weight2 * kernel;
        for(unsigned i = 0; i < nDof1; i++) {
          for(unsigned j = 0; j < nDof1; j++) {
            double jacValue11 =  1. * C * (phi1_ig[i]) * phi1_ig[j];
            _jac11[i * nDof1 + j] -= jacValue11;
            _res1[i] += jacValue11 * solu1[j];

            //_jacl11[level][i * nDof1 + j] -= jacValue11;
            //_resl1[level][i] +=  jacValue11 * solu1[j];

          }


          for(unsigned j = 0; j < nDof2; j++) {
            double jacValue12 = - 1. * C * (phi1_ig[i]) * phi2F[j];
            _jac12[i * nDof2 + j] -= jacValue12;
            _res1[i] +=  jacValue12 * solu2[j];

            //_jacl12[level][i * nDof2 + j] -= jacValue12;
            //_resl1[level][i] +=  jacValue12 * solu2[j];

          }//endl j loop

          //_res1[i] -=  2 * C * phi1_ig[i] * ( xg1[0] * xg1[0] -  xg2[0] * xg2[0]);

        }
        for(unsigned i = 0; i < nDof2; i++) {
          for(unsigned j = 0; j < nDof1; j++) {
            double jacValue21 = 1. * C * (- phi2F[i]) * phi1_ig[j];
            _jac21[i * nDof1 + j] -= jacValue21;
            _res2[i] +=  jacValue21 * solu1[j];

            //_jacl21[level][i * nDof1 + j] -= jacValue21;
            //_resl2[level][i] +=  jacValue21 * solu1[j];

          }


          for(unsigned j = 0; j < nDof2; j++) {
            double jacValue22 = - 1. * C * (- phi2F[i]) * phi2F[j];
            _jac22[i * nDof2 + j] -= jacValue22;
            _res2[i] += jacValue22 * solu2[j];

            //_jacl22[level][i * nDof2 + j] -= jacValue22;
            //_resl2[level][i] +=  jacValue22 * solu2[j];

          }//endl j loop

          //_res2[i] -=  C * phi2F[i] * (xg2[0] * xg2[0] - xg1[0] * xg1[0]);

        } //endl i loop
      }//end if U > 0.
    }//end jg loop

    if(printMesh) this->PrintElement(xv2, refineElement);
  }

  value[0] += lvalue[0];
  value[1] += lvalue[1];
  value[2] += lvalue[2];
  value[3] += lvalue[3];

  return area;
}


void NonLocal::Assembly1(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax1, const unsigned &levelMax2,
                         const unsigned &iFather, RefineElement &refineElement1, RefineElement &refineElement2,
                         const vector < double >  &solu1, const vector < double > &solu2,
                         const double &kappa, const double &delta, const bool &ielEqualJel, const bool &printMesh) {

  if(level < levelMax1 - 1) {
    refineElement1.BuildElement1Prolongation(level, iFather);
    for(unsigned i = 0; i < refineElement1.GetNumberOfChildren(); i++) {
      Assembly1(level + 1, levelMin, levelMax1, levelMax2, i, refineElement1, refineElement2,
                solu1, solu2, kappa, delta, ielEqualJel, printMesh);
    }
  }
  else {
    const unsigned &nDof1 = refineElement1.GetNumberOfNodes();
    const unsigned &dim = refineElement1.GetDimension();
    const std::vector < std::vector <double> >  &xv1 = refineElement1.GetElement1NodeCoordinates(level, iFather);
    const elem_type &finiteElement1 = refineElement1.GetFEMFine();

    std::vector < double> xg1(dim);
    std::vector < double> xi1Fg(dim);
    double weight1;
    const double *phi1;

    std::vector < double > phi1F(nDof1);
    const std::vector < std::vector <double> >  &xi1F = refineElement1.GetElement1LocalCoordinates(level, iFather);

    for(unsigned ig = 0; ig < finiteElement1.GetGaussPointNumber(); ig++) {

      finiteElement1.GetGaussQuantities(xv1, ig, weight1, phi1);
      xg1.assign(dim, 0.);
      xi1Fg.assign(dim, 0.);
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned i = 0; i < nDof1; i++) {
          xg1[k] += xv1[k][i] * phi1[i];
          xi1Fg[k] += xi1F[k][i] * phi1[i];
        }
      }
      finiteElement1.GetPhi(phi1F, xi1Fg);

      if(ielEqualJel) {
        for(unsigned i = 0; i < nDof1; i++) {
          _res1[i] -=  - 2. * weight1  * phi1F[i]; //Ax - f (so f = - 2)
        }
      }

      Assembly2(0, levelMin, levelMax1, levelMax2, 0, refineElement1,
                nDof1, xg1, weight1, phi1F,
                solu1, solu2, kappa, delta, ielEqualJel, printMesh);
    }

  }
}

double NonLocal::Assembly2(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax1, const unsigned &levelMax2, 
                           const unsigned &iFather, RefineElement &refineElement,
                           const unsigned &nDof1, const vector < double > &xg1, const double &weight1, const vector < double > &phi1,
                           const vector < double >  &solu1, const vector < double > &solu2,
                           const double &kappa, const double &delta, const bool &ielEqualJel, const bool &printMesh) {

  double eps0l = std::max(refineElement.GetEps0() * pow(0.5, level) , refineElement.GetEps());

  double area = 0;
  std::vector <double> lvalue(4, 0.);

  const unsigned &nDof2 = refineElement.GetNumberOfNodes();
  const unsigned &numberOfChildren = refineElement.GetNumberOfChildren();
  const unsigned &dim = refineElement.GetDimension();

  const std::vector < std::vector <double> >  &xv2 = refineElement.GetElement2NodeCoordinates(level, iFather);

  bool oneNodeIsInside = true;
  bool oneNodeIsOutside = true;

  if(level < levelMax2 - 1) {
    if(level < levelMin) {
    refine:
      refineElement.BuildElement2Prolongation(level, iFather);
      for(unsigned i = 0; i < numberOfChildren; i++) {
        area += Assembly2(level + 1, levelMin, levelMax1, levelMax2, i, refineElement,
                          nDof1, xg1, weight1, phi1,
                          solu1, solu2, kappa, delta, ielEqualJel, printMesh);
      }
    }
    else {
      oneNodeIsInside = false;
      oneNodeIsOutside = false;
      double d;
      std::vector< double > xv2j(dim);

      for(unsigned j = 0; j < nDof2; j++) {
        for(unsigned k = 0; k < dim; k++) {
          xv2j[k] = xv2[k][j];
        }
        d = this->GetDistance(xg1, xv2j, delta);
        if(d > eps0l) { // check if one node is inside thick interface
          if(oneNodeIsOutside) goto refine;
          oneNodeIsInside = true;
        }
        else if(d < - eps0l) { // check if one node is outside thick interface
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
    std::vector < double> xg2(dim);
    std::vector < double> xi2Fg(dim);
    double weight2;
    const double *phi2;
    std::vector < double > phi2F(nDof2);
    double U;
    const std::vector < std::vector <double> >  &xi2F = refineElement.GetElement2LocalCoordinates(level, iFather);

    double kernel = this->GetKernel(kappa, delta, refineElement.GetEps());

    for(unsigned jg = 0; jg < finiteElement.GetGaussPointNumber(); jg++) {

      finiteElement.GetGaussQuantities(xv2, jg, weight2, phi2);
      xg2.assign(dim, 0.);
      xi2Fg.assign(dim, 0.);
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned j = 0; j < nDof2; j++) {
          xg2[k] += xv2[k][j] * phi2[j];
          xi2Fg[k] += xi2F[k][j] * phi2[j];
        }
      }
      finiteElement.GetPhi(phi2F, xi2Fg);

      if(level == levelMax2 - 1) { // only for element at level l = lmax - 1
        U = refineElement.GetSmoothStepFunction(this->GetDistance(xg1, xg2, delta));
      }
      else {
        U = 1.;
      }

      if(U > 0.) {
        area += U * weight2;

        double C =  U * weight1 * weight2 * kernel;
        
        double C11, C12, C21, C22;
        
        if (levelMax1 == levelMax2){ // symmetric
          C11 = C12 = C21 = C22 = (1. + !ielEqualJel) * C;  
        }
        
        else if(levelMax1 < levelMax2) {  //coarse external - fine internal
          C11 = C12 = 2. * C;
          C21 = C22 = 0.;  
        }
        else { // fine external - coarse internal
          C11 = C12 = 0.;
          C21 = C22 = 2. * C;  
        }
        
        for(unsigned i = 0; i < nDof1; i++) {
          for(unsigned j = 0; j < nDof1; j++) {
            double jacValue11 =  C11 * phi1[i] * phi1[j];
            _jac11[i * nDof1 + j] -= jacValue11;
            _res1[i] += jacValue11 * solu1[j];
          }

          for(unsigned j = 0; j < nDof2; j++) {
            double jacValue12 = - C12 * phi1[i] * phi2F[j];
            _jac12[i * nDof2 + j] -= jacValue12;
            _res1[i] +=  jacValue12 * solu2[j];
          }//endl j loop
        }
        for(unsigned i = 0; i < nDof2; i++) {
          for(unsigned j = 0; j < nDof1; j++) {
            double jacValue21 = C21 * (- phi2F[i]) * phi1[j];
            _jac21[i * nDof1 + j] -= jacValue21;
            _res2[i] +=  jacValue21 * solu1[j];
          }

          for(unsigned j = 0; j < nDof2; j++) {
            double jacValue22 = - C22 * (- phi2F[i]) * phi2F[j];
            _jac22[i * nDof2 + j] -= jacValue22;
            _res2[i] += jacValue22 * solu2[j];
          }//endl j loop
        } //endl i loop
      }//end if U > 0.
    }//end jg loop

    if(printMesh) this->PrintElement(xv2, refineElement);
  }

  return area;
}



void NonLocal::RefinedAssembly5(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax, const unsigned &iFather,
                                const std::vector <unsigned> &igFather, std::vector <unsigned> &igi, RefineElement &refineElement,
                                const unsigned &nDof1, const std::vector< std::vector<double>>&x1, const std::vector<double> &weight1,
                                const std::vector <const double *> phi1, const vector < double >  &solu1, const vector < double > &solu2,
                                const double &kappa, const double &delta, const bool &printMesh) {


  double eps0l = std::max(refineElement.GetEps0() * pow(0.5, level) , refineElement.GetEps());

  const unsigned &nDof2 = refineElement.GetNumberOfNodes();
  const unsigned &numberOfChildren = refineElement.GetNumberOfChildren();
  const unsigned &dim = refineElement.GetDimension();

  const std::vector < std::vector <double> >  &xv = refineElement.GetElement2NodeCoordinates(level, iFather);

  igi.resize(igFather.size());//interface
  unsigned igiSize = 0;

  std::vector <unsigned> ig(igFather.size()); //refine or boundary integral
  unsigned igSize = 0;
  unsigned igbSize;
  unsigned igrSize;

  if(level < levelMin) {
    ig = igFather;
    goto refine;
  }
  else if(level < levelMax - 1) {
    for(unsigned i = 0; i < igFather.size(); i++) {// loop only on the nodes the father asked for refinement
      bool oneNodeIsInside = false;
      bool oneNodeIsOutside = false;
      double d;
      std::vector< double > xv_j(dim, 0.);
      for(unsigned j = 0; j < nDof2; j++) {
        for(unsigned k = 0; k < dim; k++) {
          xv_j[k] = xv[k][j];
        }
        d = this->GetDistance(xv_j, x1[igFather[i]], delta);
        if(d > eps0l) { // check if the node is inside the thick interface
          if(oneNodeIsOutside) {
            ig[igSize] = igFather[i];
            igSize++;
            break;
          }
          oneNodeIsInside = true;
        }
        else if(d < -eps0l) { // check if the node is outside the thick interface
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
    igbSize = 0;
    igrSize = igSize;
  }
  else { //if(level == levelMax - 1) {
    ig = igFather;
    igbSize = igFather.size();
    igrSize = 0;
    igiSize = 0;


//    igbSize = igSize;
//   igrSize = 0;
  }
//   else {
//     igbSize = 0;
//     igrSize = igSize;
//   }

  if(igbSize + igiSize) { // at least one of the integrals have to be computed
    const elem_type &finiteElement = (igbSize) ? refineElement.GetFEMFine() : refineElement.GetFEMCoarse();
    std::vector < double> xg2(dim);
    double weight2;

    std::vector < double> xi2Fg(dim);
    const double *phiC;
    std::vector < double > phi2(nDof2);

    const std::vector < std::vector <double> >  &xiF = refineElement.GetElement2LocalCoordinates(level, iFather);

    double kernel = this->GetKernel(kappa, delta, refineElement.GetEps());

    for(unsigned jg = 0; jg < finiteElement.GetGaussPointNumber(); jg++) {
      finiteElement.GetGaussQuantities(xv, jg, weight2, phiC);
      xg2.assign(dim, 0.);
      xi2Fg.assign(dim, 0.);
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned j = 0; j < nDof2; j++) {
          xg2[k] += xv[k][j] * phiC[j];
          xi2Fg[k] += xiF[k][j] * phiC[j];
        }
      }
      finiteElement.GetPhi(phi2, xi2Fg);

      for(unsigned gb = 0; gb < igbSize; gb++) {
        double C =  refineElement.GetSmoothStepFunction(this->GetDistance(xg2 , x1[ig[gb]], delta)) * weight1[ig[gb]] * weight2 * kernel;

        for(unsigned i = 0; i < nDof1; i++) {
          for(unsigned j = 0; j < nDof1; j++) {
            double jacValue11 =  C * (phi1[ig[gb]][i]) * phi1[ig[gb]][j];
            _jac11[i * nDof1 + j] -= jacValue11;
            _res1[i] +=  jacValue11 * solu1[j];
          }
          for(unsigned j = 0; j < nDof2; j++) {
            double jacValue12 = - C * (phi1[ig[gb]][i]) * phi2[j];
            _jac12[i * nDof2 + j] -= jacValue12;
            _res1[i] +=  jacValue12 * solu2[j];
          }//endl j loop
        }
        for(unsigned i = 0; i < nDof2; i++) {
          for(unsigned j = 0; j < nDof1; j++) {
            double jacValue21 = C * (- phi2[i]) * phi1[ig[gb]][j];
            _jac21[i * nDof1 + j] -= jacValue21;
            _res2[i] +=  jacValue21 * solu1[j];
          }
          for(unsigned j = 0; j < nDof2; j++) {
            double jacValue22 = - C * (- phi2[i]) * phi2[j];
            _jac22[i * nDof2 + j] -= jacValue22;
            _res2[i] +=  jacValue22 * solu2[j];
          }//endl j loop
        } //endl i loop

      }
      for(unsigned gi = 0; gi < igiSize; gi++) {

        double C = /*refineElement.GetSmoothStepFunction(this->GetDistance(xg2 , x1[igi[gi]], delta)) **/ weight1[igi[gi]] * weight2 * kernel;

        for(unsigned i = 0; i < nDof1; i++) {
          for(unsigned j = 0; j < nDof1; j++) {
            double jacValue11 =  C * (phi1[igi[gi]][i]) * phi1[igi[gi]][j];
            _jac11[i * nDof1 + j] -= jacValue11;
            _res1[i] +=  jacValue11 * solu1[j];
          }
          for(unsigned j = 0; j < nDof2; j++) {
            double jacValue12 = - C * (phi1[igi[gi]][i]) * phi2[j];
            _jac12[i * nDof2 + j] -= jacValue12;
            _res1[i] +=  jacValue12 * solu2[j];
          }//endl j loop
        }
        for(unsigned i = 0; i < nDof2; i++) {
          for(unsigned j = 0; j < nDof1; j++) {
            double jacValue21 = C * (- phi2[i]) * phi1[igi[gi]][j];
            _jac21[i * nDof1 + j] -= jacValue21;
            _res2[i] +=  jacValue21 * solu1[j];
          }
          for(unsigned j = 0; j < nDof2; j++) {
            double jacValue22 = - C * (- phi2[i]) * phi2[j];
            _jac22[i * nDof2 + j] -= jacValue22;
            _res2[i] +=  jacValue22 * solu2[j];
          }//endl j loop
        } //endl i loop
      }
    }
    if(/*level == levelMax - 1 &&*/  printMesh) PrintElement(xv, refineElement);
  }

  if(igrSize) { // at least one of the integrals have to be refined
  refine:
    refineElement.BuildElement2Prolongation(level, iFather);
    for(unsigned i = 0; i < numberOfChildren; i++) {
      RefinedAssembly5(level + 1, levelMin, levelMax, i, ig, igi, refineElement, nDof1, x1, weight1, phi1,
                       solu1, solu2, kappa, delta, printMesh);
    }
  }

  return;
}


void NonLocal::PrintElement(const std::vector < std::vector < double> > &xv, const RefineElement &refineElement) {
  fout.open("mesh.txt", std::ios::app);

  for(unsigned j = 0; j < refineElement.GetNumberOfLinearNodes(); j++) {
    fout << xv[0][j] << " " << xv[1][j] << " " << std::endl;
  }
  fout << xv[0][0] << " " << xv[1][0] << " " << std::endl;
  fout << std::endl;

  fout.close();
}


class NonLocalBall: public NonLocal {
  public:
    NonLocalBall(): NonLocal() {};
    ~NonLocalBall() {};
    double GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &radius) const;
    double GetKernel(const double &kappa, const double &delta, const double &eps) const;
    double GetArea(const double &delta, const double &eps) const {
      return M_PI * (delta * delta + eps * eps / 11.);
    };
};

class NonLocalBox: public NonLocal {
  public:
    NonLocalBox(): NonLocal() {};
    ~NonLocalBox() {};
    double GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &halfSide) const;
    double GetKernel(const double  &kappa, const double &delta, const double &eps) const;
    double GetArea(const double &delta, const double &eps) const {
      return delta * delta;
    };
};

double NonLocalBall::GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &radius) const {
  double distance  = 0.;
  for(unsigned k = 0; k < xc.size(); k++) {
    distance += (xp[k] - xc[k]) * (xp[k] - xc[k]);
  }
  distance = radius - sqrt(distance);
  return distance;
}

double NonLocalBall::GetKernel(const double  &kappa, const double &delta, const double &eps) const {
//   return 4. *  kappa / (M_PI * ((delta - eps) * (delta - eps)  +
//                                 + 2. * (-5. / 11. * eps * eps + eps * delta))
//                         * delta * delta) ;
//   std::cout.precision(14) ;
//   std::cout<<delta <<" "<<eps<<"\n";
//   std::cout<<(1. + 6./11. * pow( eps / delta, 2) + 3./143. * pow(eps / delta, 4.) ) <<"\n";
//
//   abort();

  return 4. * kappa / (M_PI  * delta * delta * delta * delta)
         / (1. + 6. / 11. * pow(eps / delta, 2) + 3. / 143. * pow(eps / delta, 4.))  ;

}

double NonLocalBox::GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &halfSide) const {

  double distance = 0.;
  unsigned dim = xc.size();
  std::vector < double > din(2 * dim); // used only if the point is inside
  std::vector < double > dout(dim, 0.); // used only if the point is outside

  bool inside = true;
  for(unsigned k = 0; k < dim; k++) {
    din[2 * k] = xp[k] - (xc[k] - halfSide); // point minus box left-side:  < 0 -> point is outside
    din[2 * k + 1] = (xc[k] + halfSide) - xp[k]; // box right-side minus point: < 0 -> point is outside
    if(din[2 * k] < 0.) {
      dout[k] = din[2 * k];
      inside = false;
    }
    else if(din[2 * k + 1] < 0.) {
      dout[k] = din[2 * k + 1];
      inside = false;
    }
  }

  if(inside) {
    distance = *std::min_element(din.begin(), din.end());
  }
  else {
    distance = 0.;
    for(unsigned k = 0; k < dim; k++) {
      distance += dout[k] * dout[k];
    }
    distance = -sqrt(distance);
  }
  return distance;
}

double NonLocalBox::GetKernel(const double  &kappa, const double &delta, const double &eps) const {
  return 0.75 * kappa / (delta * delta * delta * delta);
}

#endif
