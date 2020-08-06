#ifndef __femus_NonLocal_hpp__
#define __femus_NonLocal_hpp__

std::ofstream fout;

class NonLocal {
  public:
    NonLocal() {};
    ~NonLocal() {};
    virtual double GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &size) const = 0;
    virtual double GetKernel(const double  &kappa, const double &delta) const = 0;

    void ZeroLocalQuantities(const unsigned &nDof1, const unsigned &nDof2);

    double RefinedAssembly(const unsigned &level, const unsigned &levelMin, const unsigned &levelMax, const unsigned &iFather,
                           RefineElement &refineElement,
                           const unsigned &nDof1, const vector < double > &xg1, const double &weight1_ig, const double *phi1_ig,
                           const vector < double >  &solu1, const vector < double > &solu2,
                           const double &kappa, const double &delta, const bool &printMesh);

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
                                 const double &kappa, const double &delta, const bool &printMesh) {

  double area = 0;

  const unsigned &nDof2 = refineElement.GetNumberOfNodes();
  const unsigned &numberOfChildren = refineElement.GetNumberOfChildren();
  const unsigned &dim = refineElement.GetDimension();

  const std::vector < std::vector <double> >  &xv2 = refineElement.GetNodeCoordinates(level, iFather);

  bool oneNodeIsInside = true;
  bool oneNodeIsOutside = true;

  if(level < levelMax - 1) {
    if(level < levelMin) {
    refine:
      refineElement.BuildElementProlongation(level, iFather);
      for(unsigned i = 0; i < numberOfChildren; i++) {
        area += RefinedAssembly(level + 1, levelMin, levelMax, i, refineElement,
                                nDof1, xg1, weight1_ig, phi1_ig,
                                solu1, solu2, kappa, delta, printMesh);
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
        if(d > refineElement.GetEps()) { // check if one node is inside thick interface
          if(oneNodeIsOutside) goto refine;
          oneNodeIsInside = true;
        }
        else if(d < -refineElement.GetEps()) { // check if one node is outside thick interface
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

    const elem_type &finiteElement = refineElement.GetFEM();
    std::vector < double> xg2(dim);
    std::vector < double> xi2Fg(dim);
    double weight2;
    const double *phi2y;
    std::vector < double > phi2F(nDof2);
    double U;
    const std::vector < std::vector <double> >  &xi2F = refineElement.GetNodeLocalCoordinates(level, iFather);

    double kernel = this->GetKernel(kappa, delta);

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
      
//       std::cout << xi2Fg[0] << " " << xi2Fg[0]<<std::endl;
//       for(unsigned j = 0; j < nDof2; j++){
//         std::cout << j<<" " << phi2F[j] << std::endl;    
//       }
      
//      exit(1);
      
      if(level == levelMax - 1) { // only for element at level l = lmax - 1
        U = refineElement.GetSmoothStepFunction(  this->GetDistance(xg1, xg2, delta) );
      }
      else{
        U = 1.;    
      }

      if(U > 0.) {
        area += weight2;
        double C =  U * weight1_ig * weight2 * kernel;
        for(unsigned i = 0; i < nDof1; i++) {
          for(unsigned j = 0; j < nDof1; j++) {
            double jacValue11 =  C * (phi1_ig[i]) * phi1_ig[j];
            _jac11[i * nDof1 + j] -= jacValue11;
            _res1[i] +=  jacValue11 * solu1[j];
          }
          for(unsigned j = 0; j < nDof2; j++) {
            double jacValue12 = - C * (phi1_ig[i]) * phi2F[j];
            _jac12[i * nDof2 + j] -= jacValue12;
            _res1[i] +=  jacValue12 * solu2[j];
          }//endl j loop
        }
        for(unsigned i = 0; i < nDof2; i++) {
          for(unsigned j = 0; j < nDof1; j++) {
            double jacValue21 = C * (- phi2F[i]) * phi1_ig[j];
            _jac21[i * nDof1 + j] -= jacValue21;
            _res2[i] +=  jacValue21 * solu1[j];
          }
          for(unsigned j = 0; j < nDof2; j++) {
            double jacValue22 = - C * (- phi2F[i]) * phi2F[j];
            _jac22[i * nDof2 + j] -= jacValue22;
            _res2[i] +=  jacValue22 * solu2[j];
          }//endl j loop
        } //endl i loop
      }//end if U > 0.
    }//end jg loop

    if(printMesh) PrintElement(xv2, refineElement);
  }

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
    double GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &radius) const;
    double GetKernel(const double &kappa, const double &delta) const;
};

class NonLocalBox: public NonLocal {
  public:
    NonLocalBox(): NonLocal() {};
    ~NonLocalBox() {};
    double GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &halfSide) const;
    double GetKernel(const double  &kappa, const double &delta) const;
};

double NonLocalBall::GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &radius) const {
  double distance  = 0.;
  for(unsigned k = 0; k < xc.size(); k++) {
    distance += (xp[k] - xc[k]) * (xp[k] - xc[k]);
  }
  distance = radius - sqrt(distance);
  return distance;
}

double NonLocalBall::GetKernel(const double  &kappa, const double &delta) const {
  return 4. / M_PI * kappa / (delta * delta * delta * delta) ;
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

double NonLocalBox::GetKernel(const double  &kappa, const double &delta) const {
  return 0.75 * kappa / (delta * delta * delta * delta);
}

#endif
