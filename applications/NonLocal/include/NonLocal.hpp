#ifndef __femus_NonLocal_hpp__
#define __femus_NonLocal_hpp__

std::ofstream fout;

class NonLocal {
  public:
    NonLocal() {};
    ~NonLocal() {};
    double GetRadius(const std::vector < double>  &xc, const std::vector < double>  &xp) {
      double distance  = 0.;
      for(unsigned k = 0; k < xc.size(); k++) {
        distance += (xp[k] - xc[k]) * (xp[k] - xc[k]);
      }
      return sqrt(distance);

    };
    virtual double GetInterfaceDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &size) const = 0;
    virtual double GetKernel(const double  &kappa, const double &delta, const double &eps) const = 0;
    virtual double GetArea(const double &delta, const double &eps) const = 0;
    virtual double GetGamma(const double &d) const = 0;


    void ZeroLocalQuantities(const unsigned &nDof1, const unsigned &nDof2);

    void AddFineLevelLocalQuantities(const unsigned &level);

    void Assembly1(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax1, const unsigned &levelMax2,
                   const unsigned &iFather, RefineElement &refineElement1, RefineElement &refineElement2,
                   const vector < double >  &solu1, const vector < double > &solu2,
                   const double &kappa, const double &delta, const bool &ielEqualJel, const bool &printMesh);

    void Assembly1(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax1, const unsigned &levelMax2,
                   const unsigned &iFather, const OctTreeElement & octTreeElement1,
                   const std::vector < std::pair<std::vector<double>::iterator, std::vector<double>::iterator> > &x2MinMax,
                   RefineElement &refineElement1, RefineElement &refineElement2,
                   const vector < double >  &solu1, const vector < double > &solu2,
                   const double &kappa, const double &delta, const bool &ielEqualJel, const bool &printMesh);


    double Assembly2(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax1, const unsigned &levelMax2,
                     const unsigned &iFather, RefineElement &refineElement,
                     const unsigned &nDof1, const vector < double > &xg1, const double &weight1_ig, const vector < double > &phi1_ig,
                     const vector < double >  &solu1, const vector < double > &solu2,
                     const double &kappa, const double &delta, const bool &ielEqualJel, const bool &printMesh);

    double GetSmoothTestFunction(const double &dg1, const double &eps);

    std::vector < double > & GetRes1() {
      return _res1;
    };
    std::vector < double > & GetRes2() {
      return _res2;
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
    std::vector < double > _jac21;
    std::vector < double > _jac22;

    void PrintElement(const std::vector < std::vector < double> > &xv, const RefineElement &refineElement);

};

void NonLocal::ZeroLocalQuantities(const unsigned &nDof1, const unsigned &nDof2) {

  _jac21.assign(nDof2 * nDof1, 0.);
  _jac22.assign(nDof2 * nDof2, 0.);
  _res1.assign(nDof1, 0.);
  _res2.assign(nDof2, 0.);
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
    const elem_type &finiteElement1 = refineElement1.GetFEM1();

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

      Assembly2(0, levelMin, levelMax1, levelMax2, 0, refineElement2,
                nDof1, xg1, weight1, phi1F,
                solu1, solu2, kappa, delta, ielEqualJel, printMesh);
    }

  }
}




void NonLocal::Assembly1(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax1, const unsigned &levelMax2,
                         const unsigned &iFather, const OctTreeElement &octTreeElement1,
                         const std::vector < std::pair<std::vector<double>::iterator, std::vector<double>::iterator> > &x2MinMax,
                         RefineElement &refineElement1, RefineElement &refineElement2,
                         const vector < double >  &solu1, const vector < double > &solu2,
                         const double &kappa, const double &delta, const bool &ielEqualJel, const bool &printMesh) {

  if(level < levelMax1 - 1) {
    refineElement1.BuildElement1Prolongation(level, iFather);
    for(unsigned i = 0; i < refineElement1.GetNumberOfChildren(); i++) {
      Assembly1(level + 1, levelMin, levelMax1, levelMax2, i,
                *octTreeElement1.GetElement(std::vector<unsigned> {i}), x2MinMax,
                refineElement1, refineElement2,
                solu1, solu2, kappa, delta, ielEqualJel, printMesh);
    }
  }
  else {
    const unsigned &nDof1 = refineElement1.GetNumberOfNodes();
    const unsigned &dim = refineElement1.GetDimension();
    const std::vector < std::vector <double> >  &xv1 = refineElement1.GetElement1NodeCoordinates(level, iFather);

    std::vector < std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> > x1MinMax(dim);
    for(unsigned k = 0; k < dim; k++) {
      x1MinMax[k] = std::minmax_element(xv1[k].begin(), xv1[k].end());
    }

    double eps = refineElement1.GetEps();

    bool coarseIntersectionTest = true;
    for(unsigned k = 0; k < dim; k++) {
      if((*x1MinMax[k].first  - *x2MinMax[k].second) > delta + eps  || (*x2MinMax[k].first  - *x1MinMax[k].second) > delta + eps) {
        coarseIntersectionTest = false;
        break;
      }
    }

    const elem_type &finiteElement1 = refineElement1.GetFEM1();

    std::vector < double> xg1(dim);
    double weight1 = 0.1;
    const double *phi1;

    const std::vector < std::vector < double> > & phi1F = octTreeElement1.GetGaussShapeFunctions();

    for(unsigned ig = 0; ig < finiteElement1.GetGaussPointNumber(); ig++) {

      finiteElement1.GetGaussQuantities(xv1, ig, weight1, phi1);
      //phi1 = finiteElement1.GetPhi(ig);
      xg1.assign(dim, 0.);
      for(unsigned k = 0; k < dim; k++) {
        for(unsigned i = 0; i < nDof1; i++) {
          xg1[k] += xv1[k][i] * phi1[i];
        }
      }

      if(ielEqualJel) {
        for(unsigned i = 0; i < nDof1; i++) {
          _res1[i] -=  - 2. * weight1  * phi1F[ig][i]; //Ax - f (so f = - 2)
        }
      }
      if(coarseIntersectionTest) {
        Assembly2(0, levelMin, levelMax1, levelMax2, 0, refineElement2,
                  nDof1, xg1, weight1, phi1F[ig],
                  solu1, solu2, kappa, delta, ielEqualJel, printMesh);
      }
    }

  }
}




double NonLocal::Assembly2(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax1, const unsigned &levelMax2,
                           const unsigned &iFather, RefineElement &refineElement,
                           const unsigned &nDof1, const vector < double > &xg1, const double &weight1, const vector < double > &phi1,
                           const vector < double >  &solu1, const vector < double > &solu2,
                           const double &kappa, const double &delta, const bool &ielEqualJel, const bool &printMesh) {


  double area = 0;
  const unsigned &nDof2 = refineElement.GetNumberOfNodes();
  const unsigned &dim = refineElement.GetDimension();

  const std::vector < std::vector <double> >  &xv2 = refineElement.GetElement2NodeCoordinates(level, iFather);

  const elem_type &finiteElement = refineElement.GetFEM2();
  std::vector < double> xg2(dim);
  double weight2;
  const double *phi2;
  double U;
  double twoWeigh1Kernel = 2. * weight1 * this->GetKernel(kappa, delta, refineElement.GetEps());
  const double *phi2pt;
  std::vector<double>::iterator xg2it;
  std::vector < std::vector<double> > ::const_iterator xv2it;
  std::vector<double>::const_iterator xv2kit;

  for(unsigned jg = 0; jg < finiteElement.GetGaussPointNumber(); jg++) {

    finiteElement.GetGaussQuantities(xv2, jg, weight2, phi2);
    xg2.assign(dim, 0.);

    for(xg2it = xg2.begin(), xv2it = xv2.begin() ; xg2it != xg2.end(); xg2it++, xv2it++) {
      for(xv2kit = (*xv2it).begin(), phi2pt = phi2;  xv2kit != (*xv2it).end(); phi2pt++, xv2kit++) {
        *xg2it += (*xv2kit) * (*phi2pt);
      }
    }

    U = refineElement.GetSmoothStepFunction(this->GetInterfaceDistance(xg1, xg2, delta)) * GetGamma(GetRadius(xg1, xg2));
    if(U > 0.) {
      area += U * weight2;
      double C =  U *  weight2 * twoWeigh1Kernel;
      for(unsigned i = 0; i < nDof2; i++) {
        for(unsigned j = 0; j < nDof1; j++) {
          double jacValue21 = C * (- phi2[i]) * phi1[j];
          _jac21[i * nDof1 + j] -= jacValue21;
          _res2[i] +=  jacValue21 * solu1[j];
        }

        for(unsigned j = 0; j < nDof2; j++) {
          double jacValue22 = C * phi2[i] * phi2[j];
          _jac22[i * nDof2 + j] -= jacValue22;
          _res2[i] += jacValue22 * solu2[j];
        }//endl j loop
      } //endl i loop
    }//end if U > 0.
  }//end jg loop

  if(printMesh) this->PrintElement(xv2, refineElement);

  return area;
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

    double GetInterfaceDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &radius) const {
      double distance  = 0.;
      for(unsigned k = 0; k < xc.size(); k++) {
        distance += (xp[k] - xc[k]) * (xp[k] - xc[k]);
      }
      distance = radius - sqrt(distance);
      return distance;
    }

    double GetKernel(const double  &kappa, const double &delta, const double &eps) const {
      return 4. * kappa / (M_PI  * delta * delta * delta * delta)
             / (1. + 6. / 11. * pow(eps / delta, 2) + 3. / 143. * pow(eps / delta, 4.))  ;
    }

    double GetArea(const double &delta, const double &eps) const {
      return M_PI * (delta * delta + eps * eps / 11.);
    };

    double GetGamma(const double &d) const {
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

    double GetKernel(const double  &kappa, const double &delta, const double &eps) const {
      return 3. * kappa / (M_PI  * delta * delta * delta)
             / (1. + 3. / 11. * pow(eps / delta, 2.))  ;
    }

    double GetArea(const double &delta, const double &eps) const {
      return 2. * M_PI * delta;
    };

    double GetGamma(const double &d) const {
      return 1. / d;
    }
};




class NonLocalBox: public NonLocal {
  public:
    NonLocalBox(): NonLocal() {};
    ~NonLocalBox() {};
    double GetInterfaceDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &halfSide) const;
    double GetKernel(const double  &kappa, const double &delta, const double &eps) const;
    double GetArea(const double &delta, const double &eps) const {
      return delta * delta;
    };

    double GetGamma(const double &d) const {
      return 1.;
    }
};


double NonLocalBox::GetInterfaceDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &halfSide) const {

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
