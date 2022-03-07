#ifndef __femus_NonLocal_hpp__
#define __femus_NonLocal_hpp__

std::ofstream fout;

class NonLocal {
  public:
    NonLocal() {};
    ~NonLocal() {};
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
                         const OctTreeElement &octTreeElement1, RefineElement &element1,
                         const Region &region2, const std::vector <unsigned> &jelIndexF, const vector < double >  &solu1,
                         const double &kappa, const double &delta, const bool &printMesh);

    double AssemblyCutFem2(const RefineElement & element1, const Region & region2, const std::vector<unsigned> &jelIndex,
                           const unsigned & nDof1, const vector < double > &xg1,
                           const double & twoWeigh1Kernel, const vector < double > &phi1, const vector < double >  &solu1,
                           const double & delta, const bool & printMesh);

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





void NonLocal::AssemblyCutFem1(const unsigned &level, const unsigned &levelMin1, const unsigned &levelMax1, const unsigned &iFather,
                               const OctTreeElement &octTreeElement1, RefineElement &element1,
                               const Region &region2, const std::vector <unsigned> &jelIndexF, const vector < double >  &solu1,
                               const double &kappa, const double &delta, const bool &printMesh) {


  if(level < levelMin1) {
    element1.BuildElement1Prolongation(level, iFather);
    for(unsigned i = 0; i < element1.GetNumberOfChildren(); i++) {
      AssemblyCutFem1(level + 1, levelMin1, levelMax1, i,
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

    //preloop to store commo quatiies in iel
    for(unsigned ig = 0; ig < fem1->GetGaussPointNumber(); ig++) {
      fem1->GetGaussQuantities(xv1, ig, weight1, phi1);
      xg1.assign(dim, 0.);
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned i = 0; i < nDof1; i++) {
          xg1[k] += xv1[k][i] * phi1[i];
        }
      }
      //close prelop

      std::vector < std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> > x1MinMax(dim);
      for(unsigned k = 0; k < dim; k++) {
        x1MinMax[k] = std::minmax_element(xv1[k].begin(), xv1[k].end());
      }

      //BEGIN NEW STUFF
      
      for(unsigned jj = 0; jj < jelIndexF.size(); jj++) {
        unsigned jel = jelIndexF[jj];
        const unsigned &dim = region2.GetDimension(jel);
        const std::vector<std::vector<double>>& x2MinMax = region2.GetMinMax(jel);

        const elem_type *fem = region2.GetFem(jel);

        const std::vector <double >  &solu2g = region2.GetGaussSolution(jel);
        const std::vector <double >  &weight2 = region2.GetGaussWeight(jel);
        const std::vector < std::vector <double> >  &xg2 = region2.GetGaussCoordinates(jel);

        for(unsigned jg = 0; jg < fem->GetGaussPointNumber(); jg++) {

          bool coarseIntersectionTest = true;
          for(unsigned k = 0; k < dim; k++) {
            if((xg2[jg][k]  - * (x1MinMax[k].second)) > delta + eps  || (*(x1MinMax[k].first) - xg2[jg][k]) > delta + eps) {
              coarseIntersectionTest = false;
              break;
            }
          }
        }
      }

      //END NEW STUFF
      AssemblyCutFem2(element1, region2, jelIndexF, nDof1, xg1, 2. * weight1 * _kernel,
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
        AssemblyCutFem2(element1, region2, _jelIndexI, nDof1, xg1, 2. * weight1 * _kernel,
                        phi1F[ig], solu1, delta, printMesh);
      }
    }
    if(_jelIndexR[level].size() > 0) {
      element1.BuildElement1Prolongation(level, iFather);
      for(unsigned i = 0; i < element1.GetNumberOfChildren(); i++) {
        AssemblyCutFem1(level + 1, levelMin1, levelMax1, i,
                        *octTreeElement1.GetElement(std::vector<unsigned> {i}), element1, region2, _jelIndexR[level],
                        solu1, kappa, delta, printMesh);
      }
    }
  }
}




double NonLocal::AssemblyCutFem2(const RefineElement & element1, const Region & region2, const std::vector<unsigned> &jelIndex,
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


