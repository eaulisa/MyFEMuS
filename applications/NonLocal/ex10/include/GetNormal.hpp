#ifndef __femus_GetNormal_hpp__
#define __femus_GetNormal_hpp__



class BallApproximation {
  public:


    void GetNormalTet(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double & R,  std::vector<double> &b, double & db, unsigned & cut);
    void GetNormalQuad(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &b, double &db, unsigned &cut);
    void GetNormalTri(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &b, double &db, unsigned &cut);
    void GetNormal(const unsigned &elType, const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &b, double &db, unsigned &cut) {
      switch(elType)  {
        case 1:
          GetNormalTet(xv, xg, R, b, db, cut);
          break;
        case 3:
          GetNormalQuad(xv, xg, R, b, db, cut);
          break;
        case 4:
          GetNormalQuad(xv, xg, R, b, db, cut);
          break;
        default:
          std::cout << "Element type " << elType << " in GetNormal class not yet implemented\n";
          abort();
      }
    }


    double GetHeightPolyhedronSphereInt(const std::vector < std::vector <double> > &b, const std::vector <double> &a, const std::vector <double> &xg, const double &R);


  private:
    std::vector<double> _a;
    std::vector<double> _xm;
    std::vector<double> _xi;
    std::vector<double> _dist;
    std::vector<double> _dist0;
    std::vector <double> _theta;

};



void BallApproximation::GetNormalQuad(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &b, double &db, unsigned &cut) {

  const unsigned &dim =  xv.size();
  const unsigned nve =  4;

  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& x4 = xv[0][3];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];
  const double& y4 = xv[1][3];

  double hx = 0.5 * (fabs(x3 - x1) + fabs(x4 - x2));
  double hy = 0.5 * (fabs(y3 - y1) + fabs(y4 - y2));
  double h = sqrt(hx * hx + hy * hy);
  double eps = 1.0e-10 * h;

  _dist.assign(nve, 0);
  _dist0.resize(nve);
  unsigned cnt0 = 0;
  for(unsigned i = 0; i < nve; i++) {
    for(unsigned k = 0;  k < dim; k++) {
      _dist[i] += (xv[k][i] - xg[k]) * (xv[k][i] - xg[k]);
    }
    _dist[i] = sqrt(_dist[i]) - R;

    if(fabs(_dist[i]) < eps) {
      _dist0[i] = (_dist[i] < 0) ? -eps : eps;
      _dist[i] = 0.;
      cnt0++;
    }
    else {
      _dist0[i] = _dist[i];
    }
  }

  if(cnt0 > 0) {
    unsigned cntp = 0;
    for(unsigned i = 0; i < nve; i++) {
      if(_dist[i] > 0) cntp++;
      _dist[i] = _dist0[i];
    }
    if(cntp == 0) { // the element is inside the ball
      cut = 0;
      return;
    }
    else if(cntp == nve - cnt0) {  // the element in outside the ball
      cut = 2;
      return;
    }
  }

  _theta.resize(2);
  unsigned cnt = 0;
  for(unsigned e = 0; e < nve; e++) {
    unsigned ep1 = (e + 1) % nve;
    if(_dist[e] * _dist[ep1] < 0) {
      double s = 0.5  * (1 + (_dist[e] + _dist[ep1]) / (_dist[e] - _dist[ep1]));
      _theta[cnt] = atan2((1 - s) * xv[1][e] + s * xv[1][ep1]  - xg[1], (1 - s) * xv[0][e] + s * xv[0][ep1] - xg[0]) ;
      cnt++;
    }
  }

  if(cnt == 0) {
    if(_dist[0] < 0) cut = 0; // cell inside the ball
    else cut = 2; // cell outside the ball
    return;
  }
  else {
    cut = 1;
    if(_theta[0] > _theta[1]) {
      std::swap(_theta[0], _theta[1]);
    }
    double DT = _theta[1] - _theta[0];
    if(DT > M_PI) {
      std::swap(_theta[0], _theta[1]);
      _theta[1] += 2. * M_PI;
      DT = _theta[1] - _theta[0];
    }
    _xm.resize(dim);

    double d = R * sqrt(0.5 * DT / tan(0.5 * DT)) ;
    _a.resize(dim);
    _a[0] = -cos(_theta[0] + 0.5 * DT);
    _a[1] = -sin(_theta[0] + 0.5 * DT);

    for(unsigned k = 0; k < dim; k++) {
      _xm[k] = -_a[k] * d + xg[k];
    }
    //d += - _a[0] * xg[0] - _a[1] * xg[1]; //TODO

    _xi.resize(dim);
    double &u = _xi[0];
    double &v = _xi[1];

    double J[2][2];

    double dx12 = x1 - x2;
    double dx34 = x3 - x4;
    double dy12 = y1 - y2;
    double dy34 = y3 - y4;
    double hu = dx34 * dy12 - dx12 * dy34;

    double dx14 = (x1 - x4);
    double dy23 = (y2 - y3);
    double dx23 = (x2 - x3);
    double dy14 = (y1 - y4);
    double hv = dx14 * dy23 - dx23 * dy14;

    double eps2 = 1.0e-10 * h * h;

    if(fabs(hu) > eps2) {//edges 1 and 3 are not parallel
      double gu = -x4 * y1 + x3 * y2 - x2 * y3 + x1 * y4;
      double f = _xm[0] * (dy12 + dy34) - _xm[1] * (dx12 + dx34);
      double fpgu = f + gu;

      double det = sqrt(hu * (- 2. * _xm[0] * (dy14 + dy23)
                              + (2. * _xm[1] - y3 - y4) * (x1 + x2)
                              - (2. * _xm[1] - y1 - y2) * (x3 + x4))
                        + fpgu * fpgu);
      u = (fpgu + det) / hu;

      if(fabs(hv) > eps2) { //edges 2 and 4 are not parallel
        double gv = -x4 * y3 + x3 * y4 - x2 * y1 + x1 * y2;
        v = (f + gv - det) / hv;

        J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
        J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
        J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);
        J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);

      }
      else { //edges 2 and 4 are parallel
        //   std::cout << "2 and 4 are parallel\n";
        J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
        J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);

        v = (J[0][1] > eps) ?
            (0.25 * ((-1. + u) * (x1 + x4) - (1. + u) * (x3 + x2)) + _xm[0]) / J[0][1] :
            (0.25 * ((-1. + u) * (y1 + y4) - (1. + u) * (y3 + y2)) + _xm[1]) / J[1][1];

        J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
        J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);

      }
    }
    else if(fabs(hv) > eps2) {  //edges 1 and 3 are parallel, but edges 2 and 4 are not
      // std::cout << "1 and 3 are parallel\n";
      double f = _xm[0] * (dy12 + dy34) - _xm[1] * (dx12 + dx34);
      double gv = -x4 * y3 + x3 * y4 - x2 * y1 + x1 * y2;
      double fpgv = f + gv;

      double det = sqrt(hv * (- 2. * _xm[0] * (dy12 - dy34)
                              + (2. * _xm[1] - y2 - y3) * (x1 + x4)
                              - (2. * _xm[1] - y1 - y4) * (x2 + x3))
                        +  fpgv * fpgv);

      v = (fpgv - det) / hv;

      J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
      J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);

      u = (fabs(J[0][0]) > eps) ?
          (0.25 * ((-1. + v) * (x1 + x2) - (1. + v) * (x3 + x4)) + _xm[0]) / J[0][0] :
          (0.25 * ((-1. + v) * (y1 + y2) - (1. + v) * (y3 + y4)) + _xm[1]) / J[1][0];

      J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
      J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);
    }
    else { //edges 1 and 3, and  edges 2 and 4 are parallel
      //   std::cout << "Romboid\n";
      std::vector<std::vector<unsigned> > idx = {{3, 1}, {0, 2}};

      double A[2][2] = {{-dy14, dy23}, {dy12, dy34}};
      double B[2][2] = {{dx14, -dx23}, {-dx12, -dx34}};

      for(unsigned k = 0; k < 2; k++) {
        double d[2];
        for(unsigned j = 0 ; j < 2; j++) {
          double Ckj = - A[k][j] * xv[0][idx[k][j]] - B[k][j] * xv[1][idx[k][j]];
          d[j] = (A[k][j] * _xm[0] + B[k][j] * _xm[1] + Ckj) / sqrt(A[k][j] * A[k][j] + B[k][j] * B[k][j]);
        }
        _xi[k] = -1. + 2. * d[0] / (d[0] + d[1]);
      }

      J[0][0] = 0.25 * ((-1. + v) * dx12 + (1. + v) * dx34);
      J[0][1] = 0.25 * ((-1. + u) * dx14 - (1. + u) * dx23);
      J[1][0] = 0.25 * ((-1. + v) * dy12 + (1. + v) * dy34);
      J[1][1] = 0.25 * ((-1. + u) * dy14 - (1. + u) * dy23);
    }

    b.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        b[k] += J[j][k] * _a[j];
      }
    }
    double bNorm = sqrt(b[0] * b[0] + b[1] * b[1]);
    b[0] /= bNorm;
    b[1] /= bNorm;
    db = - b[0] * u - b[1] * v;
  }
}


void BallApproximation::GetNormalTri(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R, std::vector<double> &b, double &db, unsigned &cut) {

  const unsigned &dim =  xv.size();
  const unsigned nve =  3;

  //std::cout<<nve<<std::endl;

  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];

  double hx = (fabs(x2 - x1) + fabs(x3 - x2) + fabs(x3 - x1)) / 3.;
  double hy = (fabs(y2 - y1) + fabs(y3 - y2) + fabs(y3 - y1)) / 3.;

  double h = sqrt(hx * hx + hy * hy);
  double eps = 1.0e-10 * h;

  _dist.assign(nve, 0);
  _dist0.resize(nve);
  unsigned cnt0 = 0;
  for(unsigned i = 0; i < nve; i++) {
    for(unsigned k = 0;  k < dim; k++) {
      _dist[i] += (xv[k][i] - xg[k]) * (xv[k][i] - xg[k]);
    }
    _dist[i] = sqrt(_dist[i]) - R;

    if(fabs(_dist[i]) < eps) {
      _dist0[i] = (_dist[i] < 0) ? -eps : eps;
      _dist[i] = 0.;
      cnt0++;
    }
    else {
      _dist0[i] = _dist[i];
    }
  }

  if(cnt0 > 0) {
    unsigned cntp = 0;
    for(unsigned i = 0; i < nve; i++) {
      if(_dist[i] > 0) cntp++;
      _dist[i] = _dist0[i];
    }
    if(cntp == 0) { // the element is inside the ball
      cut = 0;
      return;
    }
    else if(cntp == nve - cnt0) {  // the element in outside the ball
      cut = 2;
      return;
    }
  }

  _theta.resize(2);
  unsigned cnt = 0;
  for(unsigned e = 0; e < nve; e++) {
    unsigned ep1 = (e + 1) % nve;
    if(_dist[e] * _dist[ep1] < 0) {
      double s = 0.5  * (1 + (_dist[e] + _dist[ep1]) / (_dist[e] - _dist[ep1]));
      _theta[cnt] = atan2((1 - s) * xv[1][e] + s * xv[1][ep1]  - xg[1], (1 - s) * xv[0][e] + s * xv[0][ep1] - xg[0]) ;
      cnt++;
    }
  }

  if(cnt == 0) {
    if(_dist[0] < 0) cut = 0; // cell inside the ball
    else cut = 2; // cell outside the ball
    return;
  }
  else {
    cut = 1;
    if(_theta[0] > _theta[1]) {
      std::swap(_theta[0], _theta[1]);
    }
    double DT = _theta[1] - _theta[0];
    if(DT > M_PI) {
      std::swap(_theta[0], _theta[1]);
      _theta[1] += 2. * M_PI;
      DT = _theta[1] - _theta[0];
    }
    _xm.resize(dim);

    double d = R * sqrt(0.5 * DT / tan(0.5 * DT)) ;
    _a.resize(dim);
    _a[0] = -cos(_theta[0] + 0.5 * DT);
    _a[1] = -sin(_theta[0] + 0.5 * DT);

    for(unsigned k = 0; k < dim; k++) {
      _xm[k] = -_a[k] * d + xg[k];
    }

    _xi.resize(dim);

    std::vector < std::vector < double > > J(2, std::vector<double>(2));
    J[0][0] = (-x1 + x2);
    J[0][1] = (-x1 + x3);

    J[1][0] = (-y1 + y2);
    J[1][1] = (-y1 + y3);

    double den = (x3 * y1 - x1 * y3 + x2 * J[1][1] - y2 * J[0][1]);

    _xi[0] = (x3 * y1 - x1 * y3 + _xm[0] * J[1][1] - _xm[1] * J[0][1]) / den;
    _xi[1] = (x1 * y2 - x2 * y1 - _xm[0] * J[1][0] + _xm[1] * J[0][0]) / den;

    b.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        b[k] += J[j][k] * _a[j];
      }
    }
    double bNorm = sqrt(b[0] * b[0] + b[1] * b[1]);
    b[0] /= bNorm;
    b[1] /= bNorm;
    db = - b[0] * _xi[0] - b[1] * _xi[1];
  }
}


void BallApproximation::GetNormalTet(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double & R,  std::vector<double> &a2, double & d2, unsigned & cut) {

  const unsigned dim =  3;
  const unsigned nve =  4;

  //std::cout<<nve<<std::endl;

  const double& x1 = xv[0][0];
  const double& x2 = xv[0][1];
  const double& x3 = xv[0][2];
  const double& x4 = xv[0][3];
  const double& y1 = xv[1][0];
  const double& y2 = xv[1][1];
  const double& y3 = xv[1][2];
  const double& y4 = xv[1][3];
  const double& z1 = xv[2][0];
  const double& z2 = xv[2][1];
  const double& z3 = xv[2][2];
  const double& z4 = xv[2][3];

  double hx = (fabs(x2 - x1) + fabs(x3 - x2) + fabs(x3 - x1) + fabs(x4 - x1) + fabs(x4 - x2) + fabs(x4 - x3)) / 6.;
  double hy = (fabs(y2 - y1) + fabs(y3 - y2) + fabs(y3 - y1) + fabs(y4 - y1) + fabs(y4 - y2) + fabs(y4 - y3)) / 6.;
  double hz = (fabs(z2 - z1) + fabs(z3 - z2) + fabs(z3 - z1) + fabs(z4 - z1) + fabs(z4 - z2) + fabs(z4 - z3)) / 6.;
  double h = sqrt(hx * hx + hy * hy + hz * hz);
  double eps = 1.0e-10 * h;

  _dist.assign(nve, 0);
  _dist0.resize(nve);
  unsigned cnt0 = 0;
  for(unsigned i = 0; i < nve; i++) {
    for(unsigned k = 0;  k < dim; k++) {
      _dist[i] += (xv[k][i] - xg[k]) * (xv[k][i] - xg[k]);
    }
    _dist[i] = sqrt(_dist[i]) - R;

    if(fabs(_dist[i]) < eps) {
      _dist0[i] = (_dist[i] < 0) ? -eps : eps;
      _dist[i] = 0.;
      cnt0++;
    }
    else {
      _dist0[i] = _dist[i];
    }
  }

  if(cnt0 > 0) {
    unsigned cntp = 0;
    for(unsigned i = 0; i < nve; i++) {
      if(_dist[i] > 0) cntp++;
      _dist[i] = _dist0[i];
    }
    if(cntp == 0) { // the element is inside the ball
      cut = 0;
      return;
    }
    else if(cntp == nve - cnt0) {  // the element in outside the ball
      cut = 2;
      return;
    }
  }

  std::vector < std::vector <double> > y(4, std::vector<double>(dim));
  std::vector < unsigned > i0(4);
  unsigned cnt = 0;
  for(unsigned i = 0; i < nve - 1; i++) {
    for(unsigned j = i + 1; j < nve; j++) {
      if(_dist[i] * _dist[j] < 0) {
        double s = _dist[i] / (_dist[i] - _dist[j]);
        for(unsigned k = 0; k < dim; k++) {
          y[cnt][k] = (1. - s) * xv[k][i] + s * xv[k][j];
        }
        i0[cnt] = (i + j) - (i == 0);
        cnt++;
      }
    }
  }

  if(cnt == 0) {
    if(_dist[0] < 0) cut = 0; // cell inside the ball
    else cut = 2; // cell outside the ball
    return;
  }
  else {
    cut = 1;

    if(cnt == 4) {

      if((i0[0] == 0 && i0[1] == 1) || (i0[0] == 1 && i0[1] == 2)) {
        std::swap(y[2], y[3]);
        std::swap(i0[2], i0[3]);
      }
      else {
        std::swap(y[1], y[3]);
        std::swap(y[1], y[2]);

        std::swap(i0[1], i0[3]);
        std::swap(i0[1], i0[2]);
      }

      //std::cout << i0[0] << " " << i0[1] << " " << i0[2] << " " << i0[3] << std::endl;
    }

    std::vector <double> yg(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned i = 0; i < cnt; i++) {
        yg[k] += y[i][k];
      }
      yg[k] /= cnt;
    }

    _a.resize(dim);

    std::vector < std::vector <double> > b(cnt, std::vector<double>(dim));
    for(unsigned k = 0; k < dim; k++) {
      _a[k] =  yg[k] - xg[k];
      for(unsigned i = 0; i < cnt; i++) {
        b[i][k] = y[i][k] - xg[k];
      }
    }
    double an = 0.;
    std::vector <double> bn(cnt, 0);
    for(unsigned k = 0; k < dim; k++) {
      an += _a[k] * _a[k];
      for(unsigned i = 0; i < cnt; i++) {
        bn[i] += b[i][k] * b[i][k];
      }
    }
    an = sqrt(an);
    for(unsigned i = 0; i < cnt; i++) {
      bn[i] = sqrt(bn[i]);
    }

    for(unsigned k = 0; k < dim; k++) {
      _a[k] /= an;
      for(unsigned i = 0; i < cnt; i++) {
        b[i][k] /= bn[i];
      }
    }

    double H = GetHeightPolyhedronSphereInt(b, _a, xg, R);

    _xm.resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      _xm[k] = xg[k] + _a[k] * H;
    }

    std::cout.precision(14);

    _xi.resize(dim);

    double J[3][3];
    
    J[0][0] = (-x1 + x2);
    J[0][1] = (-x1 + x3);
    J[0][2] = (-x1 + x4);

    J[1][0] = (-y1 + y2);
    J[1][1] = (-y1 + y3);
    J[1][2] = (-y1 + y4);

    J[2][0] = (-z1 + z2);
    J[2][1] = (-z1 + z3);
    J[2][2] = (-z1 + z4);

    double den =   J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
                   - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
                   + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

    double volume = den / 6.;

    _xi[0] = -(x3 * y4 * z1 - x3 * _xm[1] * z1 - x1 * y4 * z3 + x1 * _xm[1] * z3 - x3 * y1 * z4 + x1 * y3 * z4 - x1 * _xm[1] * z4 + x3 * _xm[1] * z4 +
               _xm[0] * (y3 * z1 - y4 * z1 - y1 * z3 + y4 * z3 + y1 * z4 - y3 * z4) +
               x3 * y1 * _xm[2] - x1 * y3 * _xm[2] + x1 * y4 * _xm[2] - x3 * y4 * _xm[2] +
               x4 * (_xm[1] * z1 + y1 * z3 - _xm[1] * z3 - y1 * _xm[2] + y3 * (-z1 + _xm[2]))) / den;

    _xi[1] = -(-(x2 * y4 * z1) + x2 * _xm[1] * z1 + x1 * y4 * z2 - x1 * _xm[1] * z2 + x2 * y1 * z4 - x1 * y2 * z4 + x1 * _xm[1] * z4 - x2 * _xm[1] * z4 +
               _xm[0] * (-(y2 * z1) + y4 * z1 + y1 * z2 - y4 * z2 - y1 * z4 + y2 * z4) +
               (-(x2 * y1) + x1 * y2 - x1 * y4 + x2 * y4) * _xm[2] +
               x4 * (-(_xm[1] * z1) - y1 * z2 + _xm[1] * z2 + y2 * (z1 - _xm[2]) + y1 * _xm[2])) / den;


    _xi[2] = -(x2 * y3 * z1 - x2 * _xm[1] * z1 - x1 * y3 * z2 + x1 * _xm[1] * z2 - x2 * y1 * z3 + x1 * y2 * z3 - x1 * _xm[1] * z3 + x2 * _xm[1] * z3 +
               _xm[0] * (y2 * z1 - y3 * z1 - y1 * z2 + y3 * z2 + y1 * z3 - y2 * z3) +
               x2 * y1 * _xm[2] - x1 * y2 * _xm[2] + x1 * y3 * _xm[2] - x2 * y3 * _xm[2] +
               x3 * (_xm[1] * z1 + y1 * z2 - _xm[1] * z2 - y1 * _xm[2] + y2 * (-z1 + _xm[2]))) / den;


    a2.assign(dim, 0);
    for(unsigned k = 0; k < dim; k++) {
      for(unsigned j = 0; j < dim; j++) {
        a2[k] -= J[j][k] * _a[j]; // this normal has to point toward the center of the ball, thus -=
      }
    }
    double bNorm = sqrt(a2[0] * a2[0] + a2[1] * a2[1] + a2[2] * a2[2]);
    a2[0] /= bNorm;
    a2[1] /= bNorm;
    a2[2] /= bNorm;
    d2 = - a2[0] * _xi[0] - a2[1] * _xi[1] - a2[2] * _xi[2];


    // std::cout << a2[0] << " " << a2[1] << " " << a2[2] << " " << d2 << " " << std::endl;
  }
}


double BallApproximation::GetHeightPolyhedronSphereInt(const std::vector < std::vector <double> > &b, const std::vector <double> &a, const std::vector <double> &xg, const double &R) {
  const unsigned& cnt = b.size();
  if(b.size() < 3) {
    abort();
  }
  const unsigned& dim = b[0].size();

  std::vector < std::vector <double> > v(cnt, std::vector<double>(dim, 0.));
  for(unsigned i = 0; i < cnt; i++) {
    unsigned ip1 = (i + 1) % cnt;
    for(unsigned k = 0; k < dim; k++) {
      v[i][k] += b[ip1][k] - b[i][k];
    }
  }


  double S = - M_PI * (cnt - 2u);
  for(unsigned i = 0; i < cnt; i++) {
    double dotf = 0.;
    double dotb = 0.;
    unsigned im1 = (cnt + i - 1u) % cnt;
    for(unsigned k = 0; k < dim; k++) {
      dotf += v[i][k] * b[i][k];
      dotb += v[im1][k] * b[i][k];
    }
    double PfdotPb = 0.;
    double normPf = 0.;
    double normPb = 0.;
    for(unsigned k = 0; k < dim; k++) {
      double pf = v[i][k] - dotf * b[i][k];
      double pb = - v[im1][k] + dotb * b[i][k];
      PfdotPb += pf * pb;
      normPf += pf * pf;
      normPb += pb * pb;
    }
    normPf = sqrt(normPf);
    normPb = sqrt(normPb);
    S += acos(PfdotPb / (normPf * normPb));

  }

  std::vector < std::vector <double> > x(cnt, xg);

  for(unsigned i = 0; i < cnt; i++) {
    double h = 0.;
    for(unsigned k = 0; k < dim; k++) {
      h += b[i][k] * a[k];
    }
    h = 1. / h;
    for(unsigned k = 0; k < dim; k++) {
      x[i][k] += h * b[i][k];
    }
  }
  for(unsigned i = 1; i < cnt; i++) {
    for(unsigned k = 0; k < dim; k++) {
      x[i][k] -= x[0][k];
    }
  }
  x[0] = {0., 0., 0.};

  double A = 0.;
  for(unsigned i = 1; i < cnt - 1; i++) {
    A += 0.5 * sqrt((x[i][1] * x[i + 1][2] - x[i][2] * x[i + 1][1]) * (x[i][1] * x[i + 1][2] - x[i][2] * x[i + 1][1]) +
                    (x[i][2] * x[i + 1][0] - x[i][0] * x[i + 1][2]) * (x[i][2] * x[i + 1][0] - x[i][0] * x[i + 1][2]) +
                    (x[i][0] * x[i + 1][1] - x[i][1] * x[i + 1][0]) * (x[i][0] * x[i + 1][1] - x[i][1] * x[i + 1][0]));
  }

  return R * pow(S / A, 1. / 3.);
}

#endif
