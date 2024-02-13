

bool GetPoint(const std::vector<double> &A,
              const std::vector<double> &N,
              std::vector<double> &PM);

void BuildMarkersOnConicArc(const double &dt1, const unsigned &nMin,
                            const std::vector<double> &A,
                            const std::vector <double>&Xg,
                            const std::vector<double> &Ps,
                            const std::vector<double> &Pe,
                            std::vector<std::vector<double>> &P);

void BuildMarkersOnConicArc(const double &dt1, const unsigned &nMin,
                            const std::vector<double> &A,
                            const std::vector <double>&Xg,
                            const std::vector<double> &P1,
                            const std::vector<double> &P2,
                            std::vector<std::vector<double>> &P) {


  std::vector<double> N1 = {A[3] + 2. * A[0] * P1[0] + A[1] * P1[1], A[4] + A[1] * P1[0] + 2. * A[2] * P1[1]};

  double den = sqrt(N1[0] * N1[0] + N1[1] * N1[1]);
  N1 = {N1[0] / den, N1[1] / den};

  double kappa1 = (2. * A[2] * N1[0] * N1[0] + A[1] * N1[0] * N1[1] + 2. * A[0] * N1[1] * N1[1]) / den;
  //std::cout<< "kappa1 = " << kappa1 << std::endl << std::endl;
  int signKappa = (kappa1 > 0) ? 1 : -1;

  std::vector<double> N2 = {A[3] + 2. * A[0] * P2[0] + A[1] * P2[1], A[4] + A[1] * P2[0] + 2. * A[2] * P2[1]};

  den = sqrt(N2[0] * N2[0] + N2[1] * N2[1]);
  N2 = {N2[0] / den, N2[1] / den};

  int sign = 1; // it takes points on the shorter arc
  bool pathological = false;
  if(A[1] * A[1] - 4. * A[0] *A[2] < 0.) { // only for ellipse
    std::vector<double> X1, X2;
    if(N1[0] * N2[0] + N1[1] * N2[1] > -1 + 1.0e-5) {
      GetPoint(A, { 0.5 * (N1[0] + N2[0]),  0.5 * (N1[1] + N2[1])}, X1);
      GetPoint(A, {-0.5 * (N1[0] + N2[0]), -0.5 * (N1[1] + N2[1])}, X2);
    }
    else { // pathological case for N1 =- N2
      pathological = true;
      GetPoint(A, { -N1[1], N1[0]}, X1);
      GetPoint(A, { -N2[1], N2[0]}, X2);
    }

    double dX1 = (Xg[0] - X1[0]) * (Xg[0] - X1[0]) + (Xg[1] - X1[1]) * (Xg[1] - X1[1]);
    double dX2 = (Xg[0] - X2[0]) * (Xg[0] - X2[0]) + (Xg[1] - X2[1]) * (Xg[1] - X2[1]);

    if(dX2  <  dX1) {
      sign = -1; // it takes points on the longer arc
    }
  }

  double t1 = atan2(N1[1], N1[0]);
  double t2 = atan2(N2[1], N2[0]);

  std::vector<double> N = N1;
  double theta = (t2 - t1 >= 0) ? t2 - t1 : 2. * M_PI + t2 - t1;
  bool reverse = false;
  if((sign * theta > sign * M_PI) || (pathological && sign < 0)) {
    N = N2;
    theta = 2. * M_PI - theta;
    reverse = true;
  }

  unsigned n = ceil(theta / dt1) + 1;
  if(n < nMin) n = nMin;
  double dt = theta / (n - 1);

  unsigned n0 = (signKappa > 0) ? 0 : n - 1;
  unsigned nk = (n - 1) - n0;

  P.assign(n, std::vector<double>(2));
  P[n0] = (!reverse) ? P1 : P2;
  P[nk] = (!reverse) ? P2 : P1;


  std::cout<<"aaaaaaaaaaa "<< n <<" "<<dt <<std::endl;

  std::vector<double> Nt = N;
  for(unsigned k = 1, nk = n0 + signKappa; k < n - 1; k++, nk += signKappa) {
    N = { cos(dt) * Nt[0] - sin(dt) * Nt[1], sin(dt) * Nt[0] + cos(dt) * Nt[1]};
    GetPoint(A, N, P[nk]);
    Nt = N;
  }

  for(unsigned k = 0; k < n; k++) {
    std::cout << k << " " << P[k][0] << " " << P[k][1] << std::endl;
  }
  std::cout << std::endl;
}

bool GetPoint(const std::vector<double> &A,
              const std::vector<double> &N,
              std::vector<double> &PM) {

  const double &a = A[0];
  const double &b = A[1];
  const double &c = A[2];
  const double &d = A[3];
  const double &e = A[4];
  const double &f = A[5];

  double n1 = N[0];
  double n2 = N[1];

  PM.resize(2);
  double &x = PM[0];
  double &y = PM[1];

  double deta0 = (b * n1 - a * n2);
  double deta1 = (deta0 - a * n2);
  double deta2 = deta1 * deta1;

  double detb0 =  b * n2 - c * n1;
  double detb1 = detb0 - c * n1;
  double detb2 = detb1 * detb1;

  if(deta2 > detb2) {
    double &det0 = deta0;
    double &det1 = deta1;
    double &det2 = deta2;
    double det3 = (-c * n1 * n1 + n2 * det0);
    double det4 = e * n1 - d * n2;
    double &det5 = detb1;

    double cc = f + (a * e * n1 - d * det0) * det4 / det2;
    double bb = 2. * (b * d - 2. * a * e) * det3 / det2;
    double aa = (b * b - 4. * a * c) * det3 / det2;

    if(fabs(aa) > 1.0e-5) {
      double dis  = bb * bb - 4. * aa * cc;
      if(dis >= 0.) {
        dis = sqrt(dis);
        y = (-bb + dis) / (2. * aa);
        x = (-det4 + y * det5) / det1;

        double m1 = d + 2. * a * x + b * y;
        double m2 = e + b * x + 2. * c * y;

        if(n1 * m1 + n2 * m2 < 0) {
          y = (-bb - dis) / (2. * aa);
          x = (-det4 + y * det5) / det1;
        }
      }
      else {
        return false;
      }
    }
    else {
      y = -cc / bb;
      x = (-det4 + y * det5) / det1;
    }
  }
  else {
    double &det0 = detb0;
    double &det1 = detb1;
    double &det2 = detb2;
    double det3 = -a * n2 * n2 + n1 * det0;
    double det4 = -e * n1 + d * n2;
    double &det5 = deta1;

    double cc = f + (c * d * n2 - e * det0) * det4 / det2;
    double bb = 2. * (b * e - 2. * c * d) * det3 / det2;
    double aa = (b * b - 4. * a * c) * det3 / det2;

    if(fabs(aa) > 1.0e-5) {
      double dis  = bb * bb - 4. * aa * cc;
      if(dis >= 0.) {
        dis = sqrt(dis);
        x = (-bb + dis) / (2. * aa);
        y = (-det4 + x * det5) / det1;

        double m1 = d + 2. * a * x + b * y;
        double m2 = e + b * x + 2. * c * y;

        if(n1 * m1 + n2 * m2 < 0) {
          x = (-bb - dis) / (2. * aa);
          y = (-det4 + x * det5) / det1;
        }
      }
      else {
        return false;
      }
    }
    else {
      x = -cc / bb;
      y = (-det4 + x * det5) / det1;
    }
  }
  return true;
}
