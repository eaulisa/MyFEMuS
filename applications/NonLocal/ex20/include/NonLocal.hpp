#ifndef __femus_NonLocal_hpp__
#define __femus_NonLocal_hpp__

#include "GetNormal.hpp"

std::ofstream fout;

class NonLocal {
  public:
    NonLocal() {
      _ballAprx = new BallApproximation();
    };
    ~NonLocal() {
      delete _ballAprx;
    };
    double GetDistance(const std::vector < double>  &x1, const std::vector < double>  &x2) const {
      double distance  = 0.;
      for(unsigned k = 0; k < x1.size(); k++) {
        distance += (x2[k] - x1[k]) * (x2[k] - x1[k]);
      }
      return sqrt(distance);

    };
    virtual double GetInterfaceDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &size) const = 0;
    virtual void SetKernel(const double  &kappa, const double &delta, const double &eps) = 0;
    const double & GetKernel() const {
      return _kernel;
    };
    virtual double GetArea(const double &delta, const double &eps) const = 0;
    virtual double GetGamma(const double &d) const = 0;
    virtual double GetGamma(const std::vector < double>  &x1, const std::vector < double>  &x2) const = 0;


    void ZeroLocalQuantities(const unsigned &nDof1, const Region &region2, const unsigned &levelMax1);

    void Assembly1(const unsigned &level, const unsigned &levelMin1, const unsigned &levelMax1, const unsigned &iFather,
                   const OctTreeElement &octTreeElement1, RefineElement &element1,
                   const Region &region2, const std::vector <unsigned> &jelIndexF, const vector < double >  &solu1,
                   const double &kappa, const double &delta, const bool &printMesh);

    double Assembly2(const RefineElement &element1, const Region &region2, const std::vector<unsigned> & jelIndex,
                     const unsigned &nDof1, const vector < double > &xg1,
                     const double &twoWeigh1Kernel, const vector < double > &phi1, const vector < double >  &solu1,
                     const double &delta, const bool &printMesh);

    void AssemblyCutFem1(const unsigned &level, const unsigned &levelMin1, const unsigned &levelMax1, const unsigned &iFather,
                         const OctTreeElement &octTreeElement1, const OctTreeElement &octTreeElement1CF,
                         RefineElement &element1, Region &region2, const std::vector <unsigned> &jelIndexF, const vector < double >  &solu1,
                         const double &kappa, const double &delta, const bool &printMesh);

    void AssemblyCutFemI2(const unsigned &level, const unsigned &levelMin1, const unsigned &levelMax1, const unsigned &iFather,
                          const OctTreeElement &octTreeElement1, const OctTreeElement &octTreeElement1CF,
                          RefineElement &element1, Region &region2, const std::vector <unsigned> &jelIndexF, const vector < double >  &solu1,
                          const double &kappa, const double &delta, const bool &printMesh);


    void AssemblyCutFem2(const std::vector <double> &phi1W1,
                         const double &solu1W1W2,
                         const double &W1W2,
                         const unsigned &jel,
                         const unsigned &nDof2,
                         const double *phi2,
                         const double &solu2gW1W2,
                         const double &W2);

    double GetSmoothTestFunction(const double &dg1, const double &eps);

    std::vector < double > & GetRes2(const unsigned &jel) {
      return _res2[jel];
    };

    std::vector < double > & GetJac21(const unsigned &jel) {
      return _jac21[jel];
    };
    std::vector < double > & GetJac22(const unsigned &jel) {
      return _jac22[jel];
    };

  private:
    std::vector < std::vector < double > > _res2;
    std::vector < std::vector < double > > _jac21;
    std::vector < std::vector < double > > _jac22;

    std::vector <unsigned> _jelIndexI;
    std::vector < std::vector <unsigned> >_jelIndexR;

    BallApproximation *_ballAprx;
    std::vector<double> _a;

    std::vector< std::vector < double> > _xg1;
    std::vector< std::vector < double> > _xg1CF;

    std::vector < double> _weight1;
    std::vector < double> _weight1CF;
    std::vector<double> _eqPolyWeight;

    std::vector <double> _phi1W1;
    std::vector <double> _phi1W1CF;

    std::vector < double > _phi1W1W2;
    std::vector < double > _phi2W1W2;

    std::vector<double>::iterator _jac21It;
    std::vector<double>::iterator _jac21End;
    std::vector<double>::iterator _jac22It;
    std::vector<double>::iterator _res2It;
    std::vector<double>::const_iterator _phi1W1It;
    std::vector<double>::iterator _phi1W1W2It;
    std::vector<double>::iterator _phi1W1W2Begin;
    std::vector<double>::iterator _phi1W1W2End;
    std::vector<double>::iterator _phi2W1W2It;
    std::vector<double>::iterator _phi2W1W2Begin;
    std::vector<double>::iterator _phi2W1W2End;
    const double *_phi2pt;


    double _d;
    unsigned _cut;

    void PrintElement(const std::vector < std::vector < double> > &xv, const RefineElement &refineElement);

  protected:
    double _kernel;

};

void NonLocal::ZeroLocalQuantities(const unsigned &nDof1, const Region &region2, const unsigned &levelMax1) {
  _res2.resize(region2.size());
  _jac21.resize(region2.size());
  _jac22.resize(region2.size());

  for(unsigned jel = 0; jel < region2.size(); jel++) {
    unsigned nDof2 = region2.GetDofNumber(jel);
    _jac21[jel].assign(nDof2 * nDof1, 0.);
    _jac22[jel].assign(nDof2 * nDof2, 0.);
    _res2[jel].assign(nDof2, 0.);
  }

  _jelIndexI.reserve(region2.size());
  _jelIndexR.resize(levelMax1);
  for(unsigned level = 0; level < levelMax1; level++) {
    _jelIndexR[level].reserve(region2.size());
  }

}

void NonLocal::Assembly1(const unsigned &level, const unsigned &levelMin1, const unsigned &levelMax1, const unsigned &iFather,
                         const OctTreeElement &octTreeElement1, RefineElement &element1,
                         const Region &region2, const std::vector <unsigned> &jelIndexF, const vector < double >  &solu1,
                         const double &kappa, const double &delta, const bool &printMesh) {


  if(level < levelMin1) {
    element1.BuildElement1Prolongation(level, iFather);
    for(unsigned i = 0; i < element1.GetNumberOfChildren(); i++) {
      Assembly1(level + 1, levelMin1, levelMax1, i,
                *octTreeElement1.GetElement(std::vector<unsigned> {i}), element1, region2, jelIndexF,
                solu1, kappa, delta, printMesh);
    }
  }
  else if(level == levelMax1 - 1) {
    const unsigned &dim = element1.GetDimension();
    const std::vector < std::vector <double> >  &xv1 = element1.GetElement1NodeCoordinates(level, iFather);
    double eps = element1.GetEps();
    const unsigned &nDof1 = element1.GetNumberOfNodes();
    //double kernel = this->GetKernel(kappa, delta, eps);

    const elem_type *fem1 = element1.GetFem1();

    std::vector < double> xg1(dim);
    double weight1;
    const double *phi1;

    const std::vector < std::vector < double> > & phi1F = octTreeElement1.GetGaussShapeFunctions();

    for(unsigned ig = 0; ig < fem1->GetGaussPointNumber(); ig++) {
      fem1->GetGaussQuantities(xv1, ig, weight1, phi1);
      xg1.assign(dim, 0.);
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned i = 0; i < nDof1; i++) {
          xg1[k] += xv1[k][i] * phi1[i];
        }
      }

      Assembly2(element1, region2, jelIndexF, nDof1, xg1, 2. * weight1 * _kernel,
                phi1F[ig], solu1, delta, printMesh);
    }
  }
  else {
    const unsigned &dim = element1.GetDimension();
    const std::vector < std::vector <double> >  &xv1 = element1.GetElement1NodeCoordinates(level, iFather);
    double eps = element1.GetEps();

    _jelIndexR[level].resize(0);
    _jelIndexI.resize(0);

    std::vector < std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> > x1MinMax(dim);
    for(unsigned k = 0; k < dim; k++) {
      x1MinMax[k] = std::minmax_element(xv1[k].begin(), xv1[k].end());
    }

    std::vector < double > dmM2(dim);
    std::vector < double > dMm2(dim);
    std::vector < double > dist(pow(2, dim));

    for(unsigned j = 0; j < jelIndexF.size(); j++) {
      unsigned jel = jelIndexF[j];
      const std::vector<std::vector<double>>& x2MinMax = region2.GetMinMax(jel);

      for(unsigned k = 0; k < dim; k++) {
        dmM2[k] = (*(x1MinMax[k].first) - x2MinMax[k][1]) * (*(x1MinMax[k].first) - x2MinMax[k][1]);
        dMm2[k] = (*(x1MinMax[k].second) - x2MinMax[k][0]) * (*(x1MinMax[k].second) - x2MinMax[k][0]);
      }

      if(dim == 2) {
        dist[0] = sqrt(dmM2[0] + dmM2[1]);
        dist[1] = sqrt(dMm2[0] + dmM2[1]);
        dist[2] = sqrt(dmM2[0] + dMm2[1]);
        dist[3] = sqrt(dMm2[0] + dMm2[1]);
      }
      else if(dim == 3) {
        dist[0] = sqrt(dmM2[0] + dmM2[1] + dmM2[2]);
        dist[1] = sqrt(dMm2[0] + dmM2[1] + dmM2[2]);
        dist[2] = sqrt(dmM2[0] + dMm2[1] + dmM2[2]);
        dist[3] = sqrt(dMm2[0] + dMm2[1] + dmM2[2]);

        dist[4] = sqrt(dmM2[0] + dmM2[1] + dMm2[2]);
        dist[5] = sqrt(dMm2[0] + dmM2[1] + dMm2[2]);
        dist[6] = sqrt(dmM2[0] + dMm2[1] + dMm2[2]);
        dist[7] = sqrt(dMm2[0] + dMm2[1] + dMm2[2]);
      }

      if(*std::max_element(dist.begin(), dist.end()) < delta - eps) {
        _jelIndexI.resize(_jelIndexI.size() + 1, jel);
      }
      else {
        _jelIndexR[level].resize(_jelIndexR[level].size() + 1, jel);
      }
    }
    if(_jelIndexI.size() > 0) {
      const unsigned &nDof1 = element1.GetNumberOfNodes();
      const elem_type *fem1 = element1.GetFem1();

      std::vector < double> xg1(dim);
      double weight1;
      const double *phi1;

      const std::vector < std::vector < double> > & phi1F = octTreeElement1.GetGaussShapeFunctions();

      for(unsigned ig = 0; ig < fem1->GetGaussPointNumber(); ig++) {
        fem1->GetGaussQuantities(xv1, ig, weight1, phi1);
        xg1.assign(dim, 0.);
        for(unsigned k = 0; k < dim; k++) {
          for(unsigned i = 0; i < nDof1; i++) {
            xg1[k] += xv1[k][i] * phi1[i];
          }
        }
        Assembly2(element1, region2, _jelIndexI, nDof1, xg1, 2. * weight1 * _kernel,
                  phi1F[ig], solu1, delta, printMesh);
      }
    }
    if(_jelIndexR[level].size() > 0) {
      element1.BuildElement1Prolongation(level, iFather);
      for(unsigned i = 0; i < element1.GetNumberOfChildren(); i++) {
        Assembly1(level + 1, levelMin1, levelMax1, i,
                  *octTreeElement1.GetElement(std::vector<unsigned> {i}), element1, region2, _jelIndexR[level],
                  solu1, kappa, delta, printMesh);
      }
    }
  }
}






double NonLocal::Assembly2(const RefineElement & element1, const Region & region2, const std::vector<unsigned> &jelIndex,
                           const unsigned & nDof1, const vector < double > &xg1,
                           const double & twoWeigh1Kernel, const vector < double > &phi1, const vector < double >  &solu1,
                           const double & delta, const bool & printMesh) {

  double area = 0.;

  double solu1g = 0.;
  for(unsigned i = 0; i < nDof1; i++) {
    solu1g += solu1[i] * phi1[i];
  }

  const double *phi2;
  const double *phi2pt;
  double U;

  std::vector< double > mCphi2iSum;

  const double& eps = element1.GetEps();

  for(unsigned jj = 0; jj < jelIndex.size(); jj++) {

    unsigned jel = jelIndex[jj];

    const unsigned &dim = region2.GetDimension(jel);
    const std::vector<std::vector<double>>& x2MinMax = region2.GetMinMax(jel);

    bool coarseIntersectionTest = true;
    for(unsigned k = 0; k < dim; k++) {
      if((xg1[k]  - x2MinMax[k][1]) > delta + eps  || (x2MinMax[k][0] - xg1[k]) > delta + eps) {
        coarseIntersectionTest = false;
        break;
      }
    }

    if(coarseIntersectionTest) {

      const unsigned &nDof2 = region2.GetDofNumber(jel);
      const elem_type *fem = region2.GetFem(jel);
      const std::vector <double >  &solu2g = region2.GetGaussSolution(jel);
      const std::vector <double >  &weight2 = region2.GetGaussWeight(jel);
      const std::vector < std::vector <double> >  &xg2 = region2.GetGaussCoordinates(jel);

      mCphi2iSum.assign(nDof2, 0.);

      for(unsigned jg = 0; jg < fem->GetGaussPointNumber(); jg++) {
        phi2 = fem->GetPhi(jg);
        U = element1.GetSmoothStepFunction(this->GetInterfaceDistance(xg1, xg2[jg], delta));
        if(U > 0.) {
          double C =  U * GetGamma(xg1, xg2[jg]) *  weight2[jg] * twoWeigh1Kernel;
          double *jac22pt = &_jac22[jel][0];
          for(unsigned i = 0; i < nDof2; i++) {
            double cPhi2i = C * phi2[i];
            mCphi2iSum[i] -= cPhi2i;
            unsigned j = 0;
            for(phi2pt = phi2; j < nDof2; j++, phi2pt++, jac22pt++) {
              *jac22pt -= cPhi2i * (*phi2pt);
            }
            _res2[jel][i] += cPhi2i * solu2g[jg];
          }
        }//end if U > 0.
      }//end jg loop

      unsigned ijIndex = 0;
      for(unsigned i = 0; i < nDof2; i++) {
        for(unsigned j = 0; j < nDof1; j++, ijIndex++) {
          _jac21[jel][ijIndex] -= mCphi2iSum[i] * phi1[j];
        }
        _res2[jel][i] += mCphi2iSum[i] * solu1g;
      }
    }
  }
  return area;
}





void NonLocal::AssemblyCutFemI2(const unsigned &level, const unsigned &levelMin1, const unsigned &levelMax1, const unsigned &iFather,
                               const OctTreeElement &octTreeElement1, const OctTreeElement &octTreeElement1CF,
                               RefineElement &element1, Region &region2,
                               const std::vector <unsigned> &jelIndexF, const vector < double >  &solu1,
                               const double &kappa, const double &delta, const bool &printMesh) {


  if(level < levelMin1) {
    element1.BuildElement1Prolongation(level, iFather);
    for(unsigned i = 0; i < element1.GetNumberOfChildren(); i++) {
      AssemblyCutFemI2(level + 1, levelMin1, levelMax1, i,
                      *octTreeElement1.GetElement(std::vector<unsigned> {i}),
                      *octTreeElement1CF.GetElement(std::vector<unsigned> {i}),
                      element1, region2, jelIndexF,
                      solu1, kappa, delta, printMesh);
    }
  }
  else if(level == levelMax1 - 1) {
    const unsigned &dim = element1.GetDimension();
    std::vector < std::vector <double> >  xv1 = element1.GetElement1NodeCoordinates(level, iFather);

    const unsigned &nDof1 = element1.GetNumberOfNodes();
    const elem_type *fem1 = element1.GetFem1();
    const elem_type *fem1CF = element1.GetFem1CF();

    const unsigned &ng1 = fem1->GetGaussPointNumber();
    const unsigned &ng1CF = fem1CF->GetGaussPointNumber();
    _xg1.assign(ng1, std::vector<double>(dim, 0));
    _weight1.resize(ng1);
    for(unsigned ig = 0; ig < ng1; ig++) {
      const double *phi;
      fem1->GetGaussQuantities(xv1, ig, _weight1[ig], phi);
      for(unsigned i = 0; i < nDof1; i++) {
        for(unsigned k = 0; k < dim; k++) {
          _xg1[ig][k] += phi[i] * xv1[k][i];
        }
      }
    }

    _weight1CF.resize(ng1CF);
    _xg1CF.assign(ng1CF, std::vector<double>(dim, 0));
    for(unsigned ig = 0; ig < ng1CF; ig++) {
      const double *phi;
      fem1CF->GetGaussQuantities(xv1, ig, _weight1CF[ig], phi);
      for(unsigned i = 0; i < nDof1; i++) {
        for(unsigned k = 0; k < dim; k++) {
          _xg1CF[ig][k] += phi[i] * xv1[k][i];
        }
      }
    }

    //BEGIN NEW STUFF

    std::vector < std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> > x1MinMax(dim);
    for(unsigned k = 0; k < dim; k++) {
      x1MinMax[k] = std::minmax_element(xv1[k].begin(), xv1[k].end());
    }

    for(unsigned jj = 0; jj < jelIndexF.size(); jj++) {
      unsigned jel = jelIndexF[jj];
      const std::vector<std::vector<double>>& x2MinMax = region2.GetMinMax(jel);

      const elem_type *fem2 = region2.GetFem(jel);
      const std::vector < std::vector <double> >  &xg2 = region2.GetGaussCoordinates(jel);

      for(unsigned jg = 0; jg < fem2->GetGaussPointNumber(); jg++) {

        bool coarseIntersectionTest = true;
        for(unsigned k = 0; k < dim; k++) {
          if((xg2[jg][k]  - * (x1MinMax[k].second)) > delta || (*(x1MinMax[k].first) - xg2[jg][k]) > delta) {  // this can be improved with the l2 norm
            coarseIntersectionTest = false;
            break;
          }
        }

        if(coarseIntersectionTest) {
          _ballAprx->GetNormal(element1.GetElementType(), xv1, xg2[jg], delta, _a, _d, _cut);

          if(_cut == 0) { //interior element
            double d2W1 = 0.;
            for(unsigned ig = 0; ig < ng1; ig++) {
              double d2 = 0.;
              for(unsigned k = 0; k < dim; k++) {
                d2 += (xg2[jg][k] - _xg1[ig][k]) * (xg2[jg][k] - _xg1[ig][k]);
              }
              d2W1 += d2 * _weight1[ig];
              //d2W1 += _weight1[ig];
            }
            region2.AddI2(jel, jg, d2W1);
          }
          else if(_cut == 1) { //cut element
            element1.GetCutFem()->clear();
//             element1.GetCutFem()->GetWeightWithMap(0, _a, _d, _eqPolyWeight);
//             (*element1.GetCutFem())(0, _a, _d, _eqPolyWeight);
            // element1.GetCDweight()->GetWeight(_a, _d, _eqPolyWeight);

            // BEGIN Parabola integration
            bool twoInt = false;
            for (unsigned k = 0; k < xv1.size(); k++) xv1[k].resize(element1.GetNumberOfLinearNodes());

            std::vector<double> A(6, 0.);
            A[0] = -1;
            A[1] = 0;
            A[2] = -1;
            A[3] = + 2 * xg2[jg][0];
            A[4] = + 2 * xg2[jg][1];
            A[5] = - xg2[jg][0] * xg2[jg][0] - xg2[jg][1] * xg2[jg][1] + delta * delta;
            element1.GetCDWeightPar()->GetWeight(xv1,A,_eqPolyWeight,twoInt);
            if(!twoInt) element1.GetCDweight()->GetWeight(_a, _d, _eqPolyWeight);
            // END Parabola integration



            // //TODO TMP!!
            // xv1 = {{3,3,1,1,3,2,1,2},{-1,1,1,-1,0,1,0,-1}};
            // _weight1CF.resize(ng1CF);
            // _xg1CF.assign(ng1CF, std::vector<double>(dim, 0));
            // for(unsigned ig = 0; ig < ng1CF; ig++) {
            //   const double *phi;
            //   fem1CF->GetGaussQuantities(xv1, ig, _weight1CF[ig], phi);
            //   for(unsigned i = 0; i < nDof1; i++) {
            //     for(unsigned k = 0; k < dim; k++) {
            //       _xg1CF[ig][k] += phi[i] * xv1[k][i];
            //     }
            //   }
            // }
            // twoInt = false;
            // for (unsigned k = 0; k < xv1.size(); k++) xv1[k].resize(element1.GetNumberOfLinearNodes());
            //
            // A.resize(6, 0.);
            // A[0] = -1;
            // A[1] = 0;
            // A[2] = -1;
            // A[3] = 0;
            // A[4] = 0;
            // A[5] = 5;
            // element1.GetCDWeightPar()->GetWeight(xv1,A,_eqPolyWeight,twoInt);
            // // if(!twoInt) element1.GetCDweight()->GetWeight(_a, _d, _eqPolyWeight);
            // double Area = 0;
            // for(unsigned ig = 0; ig < ng1CF; ig++) {
            //   Area += _weight1CF[ig] * _eqPolyWeight[ig];
            // }
            // std::cout << "AAAAAAA " << Area << " correct: "<< 2.318238045 <<"\n";
            //
            // _ballAprx->GetNormal(element1.GetElementType(), xv1, {0,0}, sqrt(5), _a, _d, _cut);
            // element1.GetCDweight()->GetWeight(_a, _d, _eqPolyWeight);
            // Area = 0;
            // for(unsigned ig = 0; ig < ng1CF; ig++) {
            //   Area += _weight1CF[ig] * _eqPolyWeight[ig];
            // }
            // std::cout << "AAAAAAA " << Area << " correct: "<< 2.318238045 <<"\n";







            double d2W1CF = 0.;
            for(unsigned ig = 0; ig < ng1CF; ig++) {
              double d2 = 0.;
              for(unsigned k = 0; k < dim; k++) {
                d2 += (xg2[jg][k] - _xg1CF[ig][k]) * (xg2[jg][k] - _xg1CF[ig][k]);
              }
              d2W1CF += d2 * _weight1CF[ig] * _eqPolyWeight[ig];
              //d2W1CF += _weight1CF[ig] * _eqPolyWeight[ig];
            }
            region2.AddI2(jel, jg, d2W1CF);
          }
        }
      }
    }
  }
  else {
    const unsigned &dim = element1.GetDimension();
    const std::vector < std::vector <double> >  &xv1 = element1.GetElement1NodeCoordinates(level, iFather);

    _jelIndexR[level].resize(0);
    _jelIndexI.resize(0);

    std::vector < std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> > x1MinMax(dim);
    for(unsigned k = 0; k < dim; k++) {
      x1MinMax[k] = std::minmax_element(xv1[k].begin(), xv1[k].end());
    }

    std::vector < double > dmM2(dim);
    std::vector < double > dMm2(dim);
    std::vector < double > dist(pow(2, dim));

    for(unsigned j = 0; j < jelIndexF.size(); j++) {
      unsigned jel = jelIndexF[j];
      const std::vector<std::vector<double>>& x2MinMax = region2.GetMinMax(jel);

      for(unsigned k = 0; k < dim; k++) {
        dmM2[k] = (*(x1MinMax[k].first) - x2MinMax[k][1]) * (*(x1MinMax[k].first) - x2MinMax[k][1]);
        dMm2[k] = (*(x1MinMax[k].second) - x2MinMax[k][0]) * (*(x1MinMax[k].second) - x2MinMax[k][0]);
      }

      if(dim == 2) {
        dist[0] = sqrt(dmM2[0] + dmM2[1]);
        dist[1] = sqrt(dMm2[0] + dmM2[1]);
        dist[2] = sqrt(dmM2[0] + dMm2[1]);
        dist[3] = sqrt(dMm2[0] + dMm2[1]);
      }
      else if(dim == 3) {
        dist[0] = sqrt(dmM2[0] + dmM2[1] + dmM2[2]);
        dist[1] = sqrt(dMm2[0] + dmM2[1] + dmM2[2]);
        dist[2] = sqrt(dmM2[0] + dMm2[1] + dmM2[2]);
        dist[3] = sqrt(dMm2[0] + dMm2[1] + dmM2[2]);

        dist[4] = sqrt(dmM2[0] + dmM2[1] + dMm2[2]);
        dist[5] = sqrt(dMm2[0] + dmM2[1] + dMm2[2]);
        dist[6] = sqrt(dmM2[0] + dMm2[1] + dMm2[2]);
        dist[7] = sqrt(dMm2[0] + dMm2[1] + dMm2[2]);
      }

      if(*std::max_element(dist.begin(), dist.end()) < delta) {
        _jelIndexI.resize(_jelIndexI.size() + 1, jel);
      }
      else {
        _jelIndexR[level].resize(_jelIndexR[level].size() + 1, jel);
      }
    }
    if(_jelIndexI.size() > 0) {

      const unsigned &dim = element1.GetDimension();
      const std::vector < std::vector <double> >  &xv1 = element1.GetElement1NodeCoordinates(level, iFather);

      const unsigned &nDof1 = element1.GetNumberOfNodes();
      const elem_type *fem1 = element1.GetFem1();

      //these are the shape functions of iel evaluated in the gauss points of the refined elements
      const unsigned &ng1 = fem1->GetGaussPointNumber();

      _xg1.assign(ng1, std::vector<double>(dim, 0));
      _weight1.resize(ng1);
      for(unsigned ig = 0; ig < ng1; ig++) {
        const double *phi;
        fem1->GetGaussQuantities(xv1, ig, _weight1[ig], phi);
        for(unsigned i = 0; i < nDof1; i++) {
          for(unsigned k = 0; k < dim; k++) {
            _xg1[ig][k] += phi[i] * xv1[k][i];
          }
        }
      }

      for(unsigned jj = 0; jj < _jelIndexI.size(); jj++) {
        unsigned jel = _jelIndexI[jj];
        const elem_type *fem2 = region2.GetFem(jel);

        const std::vector < std::vector <double> >  &xg2 = region2.GetGaussCoordinates(jel);

        for(unsigned jg = 0; jg < fem2->GetGaussPointNumber(); jg++) {
          double d2W1 = 0.;
          for(unsigned ig = 0; ig < ng1; ig++) {
            double d2 = 0.;
            for(unsigned k = 0; k < dim; k++) {
              d2 += (xg2[jg][k] - _xg1[ig][k]) * (xg2[jg][k] - _xg1[ig][k]);
            }
            d2W1 += d2 * _weight1[ig];
            //d2W1 += _weight1[ig];
          }
          region2.AddI2(jel, jg, d2W1);
        }
      }
    }
    if(_jelIndexR[level].size() > 0) {
      element1.BuildElement1Prolongation(level, iFather);
      for(unsigned i = 0; i < element1.GetNumberOfChildren(); i++) {
        AssemblyCutFemI2(level + 1, levelMin1, levelMax1, i,
                        *octTreeElement1.GetElement(std::vector<unsigned> {i}),
                        *octTreeElement1CF.GetElement(std::vector<unsigned> {i}),
                        element1, region2, _jelIndexR[level],
                        solu1, kappa, delta, printMesh);
      }
    }
  }
}
























void NonLocal::AssemblyCutFem1(const unsigned &level, const unsigned &levelMin1, const unsigned &levelMax1, const unsigned &iFather,
                               const OctTreeElement &octTreeElement1, const OctTreeElement &octTreeElement1CF,
                               RefineElement &element1, Region &region2,
                               const std::vector <unsigned> &jelIndexF, const vector < double >  &solu1,
                               const double &kappa, const double &delta, const bool &printMesh) {


  if(level < levelMin1) {
    element1.BuildElement1Prolongation(level, iFather);
    for(unsigned i = 0; i < element1.GetNumberOfChildren(); i++) {
      AssemblyCutFem1(level + 1, levelMin1, levelMax1, i,
                      *octTreeElement1.GetElement(std::vector<unsigned> {i}),
                      *octTreeElement1CF.GetElement(std::vector<unsigned> {i}),
                      element1, region2, jelIndexF,
                      solu1, kappa, delta, printMesh);
    }
  }
  else if(level == levelMax1 - 1) {
    const unsigned &dim = element1.GetDimension();
    const std::vector < std::vector <double> >  &xv1 = element1.GetElement1NodeCoordinates(level, iFather);

    const unsigned &nDof1 = element1.GetNumberOfNodes();
    const elem_type *fem1 = element1.GetFem1();
    const elem_type *fem1CF = element1.GetFem1CF();

    //these are the shape functions of iel evaluated in the gauss points of the refined elements
    const std::vector < std::vector < double> > & phi1 = octTreeElement1.GetGaussShapeFunctions();
    const std::vector < std::vector < double> > & phi1CF = octTreeElement1CF.GetGaussShapeFunctions();

    //these are the shape functions of the refined element evaluated in the gauss points of the refined element

    const unsigned &ng1 = fem1->GetGaussPointNumber();
    const unsigned &ng1CF = fem1CF->GetGaussPointNumber();

    double W1 = 0.;
    _phi1W1.assign(nDof1, 0.);
    _xg1.assign(ng1, std::vector<double>(dim, 0));
    _weight1.resize(ng1);
    for(unsigned ig = 0; ig < ng1; ig++) {
      const double *phi;
      fem1->GetGaussQuantities(xv1, ig, _weight1[ig], phi);
      W1 += _weight1[ig];
      for(unsigned i = 0; i < nDof1; i++) {
        _phi1W1[i] += phi1[ig][i] * _weight1[ig];
        for(unsigned k = 0; k < dim; k++) {
          _xg1[ig][k] += phi[i] * xv1[k][i];
        }
      }
    }
    double solu1W1 = 0.;
    for(unsigned i = 0; i < nDof1; i++) {
      solu1W1 += solu1[i] * _phi1W1[i];
    }

    _weight1CF.resize(ng1CF);
    _xg1CF.assign(ng1CF, std::vector<double>(dim, 0));
    for(unsigned ig = 0; ig < ng1CF; ig++) {
      const double *phi;
      fem1CF->GetGaussQuantities(xv1, ig, _weight1CF[ig], phi);
      for(unsigned i = 0; i < nDof1; i++) {
        for(unsigned k = 0; k < dim; k++) {
          _xg1CF[ig][k] += phi[i] * xv1[k][i];
        }
      }
    }

    //BEGIN NEW STUFF

    std::vector < std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> > x1MinMax(dim);
    for(unsigned k = 0; k < dim; k++) {
      x1MinMax[k] = std::minmax_element(xv1[k].begin(), xv1[k].end());
    }

    for(unsigned jj = 0; jj < jelIndexF.size(); jj++) {
      unsigned jel = jelIndexF[jj];
      const std::vector<std::vector<double>>& x2MinMax = region2.GetMinMax(jel);

      const elem_type *fem2 = region2.GetFem(jel);

      const std::vector <double >  &solu2g = region2.GetGaussSolution(jel);
      const std::vector <double >  &weight2 = region2.GetGaussWeight(jel);
      const std::vector < std::vector <double> >  &xg2 = region2.GetGaussCoordinates(jel);
      const std::vector<double>& I2 = region2.GetI2(jel);

      for(unsigned jg = 0; jg < fem2->GetGaussPointNumber(); jg++) {

        bool coarseIntersectionTest = true;
        for(unsigned k = 0; k < dim; k++) {
          if((xg2[jg][k]  - * (x1MinMax[k].second)) > delta || (*(x1MinMax[k].first) - xg2[jg][k]) > delta) {  // this can be improved with the l2 norm
            coarseIntersectionTest = false;
            break;
          }
        }

        if(coarseIntersectionTest) {
          _ballAprx->GetNormal(element1.GetElementType(), xv1, xg2[jg], delta, _a, _d, _cut);

          if(_cut == 0) { //interior element
            double W2 = 2. * weight2[jg] * _kernel * I2[jg];
            double W1W2 = W1 * W2;
            AssemblyCutFem2(_phi1W1, solu1W1 * W2,  W1W2, jel, region2.GetDofNumber(jel), fem2->GetPhi(jg), solu2g[jg] * W1W2, W2);
          }
          else if(_cut == 1) { //cut element
            element1.GetCutFem()->clear();
            //       element1.GetCutFem()->GetWeightWithMap(0, _a, _d, _eqPolyWeight);
//             (*element1.GetCutFem())(0, _a, _d, _eqPolyWeight);
            element1.GetCDweight()->GetWeight(_a, _d, _eqPolyWeight);

            double W1CF = 0.;
            _phi1W1CF.assign(nDof1, 0.);
            for(unsigned ig = 0; ig != ng1CF; ++ig) {
              double weightigjg = _weight1CF[ig] * _eqPolyWeight[ig];
              W1CF += weightigjg;
              for(unsigned i = 0; i != nDof1; ++i) {
                _phi1W1CF[i] += phi1CF[ig][i] * weightigjg;
              }
            }
            double solu1W1CF = 0.;
            for(unsigned i = 0; i != nDof1; ++i) {
              solu1W1CF += solu1[i] * _phi1W1CF[i];
            }

            double W2 = 2. * weight2[jg] * _kernel * I2[jg];
            double W1CFW2 = W1CF * W2;
            AssemblyCutFem2(_phi1W1CF, solu1W1CF * W2,  W1CFW2, jel, region2.GetDofNumber(jel), fem2->GetPhi(jg), solu2g[jg] * W1CFW2, W2);
          }
        }
      }
    }
  }
  else {
    const unsigned &dim = element1.GetDimension();
    const std::vector < std::vector <double> >  &xv1 = element1.GetElement1NodeCoordinates(level, iFather);

    _jelIndexR[level].resize(0);
    _jelIndexI.resize(0);

    std::vector < std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> > x1MinMax(dim);
    for(unsigned k = 0; k < dim; k++) {
      x1MinMax[k] = std::minmax_element(xv1[k].begin(), xv1[k].end());
    }

    std::vector < double > dmM2(dim);
    std::vector < double > dMm2(dim);
    std::vector < double > dist(pow(2, dim));

    for(unsigned j = 0; j < jelIndexF.size(); j++) {
      unsigned jel = jelIndexF[j];
      const std::vector<std::vector<double>>& x2MinMax = region2.GetMinMax(jel);

      for(unsigned k = 0; k < dim; k++) {
        dmM2[k] = (*(x1MinMax[k].first) - x2MinMax[k][1]) * (*(x1MinMax[k].first) - x2MinMax[k][1]);
        dMm2[k] = (*(x1MinMax[k].second) - x2MinMax[k][0]) * (*(x1MinMax[k].second) - x2MinMax[k][0]);
      }

      if(dim == 2) {
        dist[0] = sqrt(dmM2[0] + dmM2[1]);
        dist[1] = sqrt(dMm2[0] + dmM2[1]);
        dist[2] = sqrt(dmM2[0] + dMm2[1]);
        dist[3] = sqrt(dMm2[0] + dMm2[1]);
      }
      else if(dim == 3) {
        dist[0] = sqrt(dmM2[0] + dmM2[1] + dmM2[2]);
        dist[1] = sqrt(dMm2[0] + dmM2[1] + dmM2[2]);
        dist[2] = sqrt(dmM2[0] + dMm2[1] + dmM2[2]);
        dist[3] = sqrt(dMm2[0] + dMm2[1] + dmM2[2]);

        dist[4] = sqrt(dmM2[0] + dmM2[1] + dMm2[2]);
        dist[5] = sqrt(dMm2[0] + dmM2[1] + dMm2[2]);
        dist[6] = sqrt(dmM2[0] + dMm2[1] + dMm2[2]);
        dist[7] = sqrt(dMm2[0] + dMm2[1] + dMm2[2]);
      }

      if(*std::max_element(dist.begin(), dist.end()) < delta) {
        _jelIndexI.resize(_jelIndexI.size() + 1, jel);
      }
      else {
        _jelIndexR[level].resize(_jelIndexR[level].size() + 1, jel);
      }
    }
    if(_jelIndexI.size() > 0) {

      const unsigned &dim = element1.GetDimension();
      const std::vector < std::vector <double> >  &xv1 = element1.GetElement1NodeCoordinates(level, iFather);

      const unsigned &nDof1 = element1.GetNumberOfNodes();
      const elem_type *fem1 = element1.GetFem1();

      //these are the shape functions of iel evaluated in the gauss points of the refined elements
      const std::vector < std::vector < double> > & phi1 = octTreeElement1.GetGaussShapeFunctions();
      const unsigned &ng1 = fem1->GetGaussPointNumber();


      double W1 = 0.;
      _phi1W1.assign(nDof1, 0.);
      _xg1.assign(ng1, std::vector<double>(dim, 0));
      _weight1.resize(ng1);
      for(unsigned ig = 0; ig < ng1; ig++) {
        const double *phi;
        fem1->GetGaussQuantities(xv1, ig, _weight1[ig], phi);
        W1 += _weight1[ig];
        for(unsigned i = 0; i < nDof1; i++) {
          _phi1W1[i] += phi1[ig][i] * _weight1[ig];
          for(unsigned k = 0; k < dim; k++) {
            _xg1[ig][k] += phi[i] * xv1[k][i];
          }
        }
      }

      double solu1W1 = 0.;
      for(unsigned i = 0; i < nDof1; i++) {
        solu1W1 += solu1[i] * _phi1W1[i];
      }

      for(unsigned jj = 0; jj < _jelIndexI.size(); jj++) {
        unsigned jel = _jelIndexI[jj];
        const elem_type *fem2 = region2.GetFem(jel);

        const std::vector <double >  &solu2g = region2.GetGaussSolution(jel);
        const std::vector <double >  &weight2 = region2.GetGaussWeight(jel);
        const std::vector < std::vector <double> >  &xg2 = region2.GetGaussCoordinates(jel);
        const std::vector<double>& I2 = region2.GetI2(jel);

        for(unsigned jg = 0; jg < fem2->GetGaussPointNumber(); jg++) {
          double W2 = 2. * weight2[jg] * _kernel * I2[jg];
          double W1W2 = W1 * W2;
          AssemblyCutFem2(_phi1W1, solu1W1 * W2,  W1W2, jel, region2.GetDofNumber(jel), fem2->GetPhi(jg), solu2g[jg] * W1W2, W2);
        }
      }
    }
    if(_jelIndexR[level].size() > 0) {
      element1.BuildElement1Prolongation(level, iFather);
      for(unsigned i = 0; i < element1.GetNumberOfChildren(); i++) {
        AssemblyCutFem1(level + 1, levelMin1, levelMax1, i,
                        *octTreeElement1.GetElement(std::vector<unsigned> {i}),
                        *octTreeElement1CF.GetElement(std::vector<unsigned> {i}),
                        element1, region2, _jelIndexR[level],
                        solu1, kappa, delta, printMesh);
      }
    }
  }
}












void NonLocal::AssemblyCutFem2(const std::vector <double> &phi1W1,
                               const double &solu1W1W2,
                               const double &W1W2,
                               const unsigned &jel,
                               const unsigned &nDof2,
                               const double *phi2,
                               const double &solu2W1W2,
                               const double &W2) {

  _phi1W1W2.resize(phi1W1.size());
  _phi1W1W2Begin = _phi1W1W2.begin();
  _phi1W1W2End = _phi1W1W2.end();

  for(_phi1W1W2It = _phi1W1W2Begin, _phi1W1It = phi1W1.begin(); _phi1W1W2It !=  _phi1W1W2End; ++_phi1W1W2It, ++_phi1W1It) {
    *_phi1W1W2It = *_phi1W1It * W2;
  }

  _phi2W1W2.resize(nDof2);
  _phi2W1W2Begin = _phi2W1W2.begin();
  _phi2W1W2End = _phi2W1W2.end();

  for( _phi2pt = phi2, _phi2W1W2It = _phi2W1W2Begin; _phi2W1W2It != _phi2W1W2End; ++_phi2pt, ++_phi2W1W2It) {
    *_phi2W1W2It = *_phi2pt * W1W2;
  }

  _jac21End = _jac21[jel].end();
  for(_jac21It = _jac21[jel].begin(), _jac22It = _jac22[jel].begin(), _res2It = _res2[jel].begin(), _phi2pt = phi2; _jac21It != _jac21End; ++_phi2pt, ++_res2It) {
    for(_phi1W1W2It = _phi1W1W2Begin; _phi1W1W2It != _phi1W1W2End; ++_phi1W1W2It, ++_jac21It) {
      *_jac21It += *_phi2pt * (*_phi1W1W2It);
    }
    for(_phi2W1W2It = _phi2W1W2Begin; _phi2W1W2It != _phi2W1W2End; ++_phi2W1W2It, ++_jac22It) {
      *_jac22It -= *_phi2pt * (*_phi2W1W2It);
    }
    *_res2It += *_phi2pt * (solu2W1W2 - solu1W1W2);
  }
}



void NonLocal::PrintElement(const std::vector < std::vector < double> > &xv, const RefineElement & refineElement) {
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

    double GetInterfaceDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &radius) const {
      double distance  = 0.;
      for(unsigned k = 0; k < xc.size(); k++) {
        distance += (xp[k] - xc[k]) * (xp[k] - xc[k]);
      }
      distance = radius - sqrt(distance);
      return distance;
    };

    void SetKernel(const double  &kappa, const double &delta, const double &eps) {
      _kernel = 4. * kappa / (M_PI  * delta * delta * delta * delta)
                / (1. + 6. / 11. * pow(eps / delta, 2) + 3. / 143. * pow(eps / delta, 4.));
    }

    double GetArea(const double &delta, const double &eps) const {
      return M_PI * (delta * delta + eps * eps / 11.);
    };

    double GetGamma(const double &d) const {
      return 1.;
    }

    double GetGamma(const std::vector < double>  &x1, const std::vector < double>  &x2) const {
      return 1.;
    }
};

class NonLocalBall3D: public NonLocal {
  public:
    NonLocalBall3D(): NonLocal() {};
    ~NonLocalBall3D() {};

    double GetInterfaceDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &radius) const {
      double distance  = 0.;
      for(unsigned k = 0; k < xc.size(); k++) {
        distance += (xp[k] - xc[k]) * (xp[k] - xc[k]);
      }
      distance = radius - sqrt(distance);
      return distance;
    };

    void SetKernel(const double  &kappa, const double &delta, const double &eps) {
      _kernel = 15. * kappa / (4. * M_PI  * delta * delta * delta * delta * delta)
                / (1. + 10. / 11. * pow(eps / delta, 2) + 15. / 143. * pow(eps / delta, 4.));
    }

    double GetArea(const double &delta, const double &eps) const {
      return 4. / 3. * M_PI * (delta * delta * delta) * (1. + 3. / 11. * pow(eps / delta, 2));
    };

    double GetGamma(const double &d) const {
      return 1.;
    }

    double GetGamma(const std::vector < double>  &x1, const std::vector < double>  &x2) const {
      return 1.;
    }
};



class NonLocalBall1: public NonLocal {
  public:
    NonLocalBall1(): NonLocal() {};
    ~NonLocalBall1() {};

    double GetInterfaceDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &radius) const {
      double distance  = 0.;
      for(unsigned k = 0; k < xc.size(); k++) {
        distance += (xp[k] - xc[k]) * (xp[k] - xc[k]);
      }
      distance = radius - sqrt(distance);
      return distance;
    }

    void SetKernel(const double  &kappa, const double &delta, const double &eps) {
      _kernel = 3. * kappa / (M_PI  * delta * delta * delta)
                / (1. + 3. / 11. * pow(eps / delta, 2.))  ;
    }

    double GetArea(const double &delta, const double &eps) const {
      return 2. * M_PI * delta;
    };

    double GetGamma(const double &d) const {
      return 1. / d;
    }

    double GetGamma(const std::vector < double>  &x1, const std::vector < double>  &x2) const {
      return 1. / GetDistance(x1, x2);
    }

};




class NonLocalBox: public NonLocal {
  public:
    NonLocalBox(): NonLocal() {};
    ~NonLocalBox() {};
    double GetInterfaceDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &halfSide) const;
    void SetKernel(const double  &kappa, const double &delta, const double &eps) {
      _kernel = 0.75 * kappa / (delta * delta * delta * delta);
    };
    double GetArea(const double &delta, const double &eps) const {
      return delta * delta;
    };

    double GetGamma(const double &d) const {
      return 1.;
    }
    double GetGamma(const std::vector < double>  &x1, const std::vector < double>  &x2) const {
      return 1.;
    }
};


double NonLocalBox::GetInterfaceDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double & halfSide) const {

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


#endif






