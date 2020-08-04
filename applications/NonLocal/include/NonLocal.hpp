#ifndef __femus_NonLocal_hpp__
#define __femus_NonLocal_hpp__

class NonLocal {
  public:
    NonLocal() {};
    ~NonLocal() {};
    virtual double GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &bSide) = 0;
    virtual double GetKernel(const double  &kappa, const double &delta) = 0;
};


class NonLocalBall: public NonLocal {
  public:
    NonLocalBall(): NonLocal() {};
    ~NonLocalBall() {};
    virtual double GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &bSide);
    virtual double GetKernel(const double  &kappa, const double &delta);
};

class NonLocalBox: public NonLocal {
  public:
    NonLocalBox(): NonLocal() {};
    ~NonLocalBox() {};
    virtual double GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &bSide);
    virtual double GetKernel(const double  &kappa, const double &delta);
};

double NonLocalBall::GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &bSide) {
  double distance  = 0.;
  for(unsigned k = 0; k < xc.size(); k++) {
    distance += (xp[k] - xc[k]) * (xp[k] - xc[k]);
  }
  distance = bSide - sqrt(distance);
  return distance;
}

double NonLocalBall::GetKernel(const double  &kappa, const double &delta) {
  return 4. / M_PI * kappa / (delta * delta * delta * delta) ;
}

double NonLocalBox::GetDistance(const std::vector < double>  &xc, const std::vector < double>  &xp, const double &bSide) {

  double distance = 0.;
  unsigned dim = xc.size();
  std::vector < double > din(2 * dim); // used only if the point is inside
  std::vector < double > dout(dim, 0.); // used only if the point is outside

  bool inside = true;
  for(unsigned k = 0; k < dim; k++) {
    din[2 * k] = xp[k] - (xc[k] - bSide); // point minus box left-side:  < 0 -> point is outside
    din[2 * k + 1] = (xc[k] + bSide) - xp[k]; // box right-side minus point: < 0 -> point is outside
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

double NonLocalBox::GetKernel(const double  &kappa, const double &delta) {
  return 0.75 * kappa / (delta * delta * delta * delta);
}


#endif
