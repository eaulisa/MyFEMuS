
#include "FemusInit.hpp"

//2D NONLOCAL EX : nonlocal diffusion for a body with different material properties

#include <vector>
#include <cmath>

void GetNormalQuad(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R,  std::vector<double> &a, double &d);

int main() {

  std::vector < std::vector<double> > xva = {{1., 1.}, {2., 1.}, {2., 2.}, {1., 2.}};
  std::vector < std::vector<double> > xvb =  {{2., 1.}, {2., 2.}, {1., 2.}, {1., 1.}};
  std::vector < std::vector<double> > xvc = {{2., 2.}, {1., 2.}, {1., 1.}, {2., 1.}};
  std::vector < std::vector<double> > xvd = {{1., 2.}, {1., 1.}, {2., 1.}, {2., 2.}};
  std::vector < std::vector<double> > xve = {{-1., -0.5}, {-1., 0.5}, {-2., .5}, {-2., -0.5}};
  std::vector<double> xg(2, 0);
  double R = 2.;
  std::vector<double> a;
  double d;
  GetNormalQuad(xva, xg, R, a, d);
  GetNormalQuad(xvb, xg, R, a, d);
  GetNormalQuad(xvc, xg, R, a, d);
  GetNormalQuad(xvd, xg, R, a, d);
  GetNormalQuad(xve, xg, R, a, d);

  xg[0] -= 2;
  xg[1] -= 2;

  for(unsigned i = 0; i < xva.size(); i++) {
    for(unsigned k = 0; k < xva[i].size(); k++) {
      xva[i][k] -= 2;
      xvb[i][k] -= 2;
      xvc[i][k] -= 2;
      xvd[i][k] -= 2;
      xve[i][k] -= 2;
    }
  }
  GetNormalQuad(xva, xg, R, a, d);
  GetNormalQuad(xvb, xg, R, a, d);
  GetNormalQuad(xvc, xg, R, a, d);
  GetNormalQuad(xvd, xg, R, a, d);
  GetNormalQuad(xve, xg, R, a, d);


  return 1;
}


void GetNormalQuad(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R,  std::vector<double> &a, double &d) {

  unsigned nve =  xv.size();
  unsigned dim =  xv[0].size();

  std::vector<double> dist(nve, 0);
  for(unsigned i = 0; i < nve; i++) {
    for(unsigned k = 0;  k < dim; k++) {
      dist[i] += (xv[i][k] - xg[k]) * (xv[i][k] - xg[k]);
    }
    dist[i] = sqrt(dist[i]) - R;
    if(dist[i] == 0.) dist[i] = 1.0e-10;
  }

  std::vector <double> theta(2);
  unsigned cnt = 0;
  for(unsigned e = 0; e < nve; e++) {
    unsigned ep1 = (e + 1) % nve;
    if(dist[e] * dist[ep1] < 0) {

      double s = 0.5  * (1 + (dist[e] + dist[ep1]) / (dist[e] - dist[ep1]));
      theta[cnt] = atan2((1 - s) * xv[e][1] + s * xv[ep1 ][1]  - xg[1], (1 - s) * xv[e][0] + s * xv[ep1][0] - xg[0]) ;
      cnt++;

    }

  }

  if(theta[0] > theta[1]) {
    std::swap(theta[0], theta[1]);
  }

  double DT = theta[1] - theta[0];
  if(DT > M_PI) {
    std::swap(theta[0], theta[1]);
    theta[1] += 2. * M_PI;
    DT = theta[1] - theta[0];
  }

  std::cout << theta[0] / M_PI * 180 << " " << theta[1] / M_PI * 180 << " " << DT / M_PI * 180 << std::endl;

  std::vector < double > xm(dim);


  d = R * sqrt(0.5 * DT / tan(0.5 * DT)) ;
  a.resize(dim);
  a[0] = -cos(theta[0] + 0.5 * DT);
  a[1] = -sin(theta[0] + 0.5 * DT);


  xm[0] = (a[0] == 0) ? 0 : 0.5 * (-d / a[0]);

  for(unsigned k = 0; k < dim; k++) {
    xm[k] = -a[k] * d;
    xm[k] += xg[k];
  }
  d += - a[0] * xg[0] - a[1] * xg[1];

  std::cout << "xm = " << xm[0] << " " << xm[1] << std::endl;
  std::cout << "a = " << a[0] << " b = " << a[1] << " d = " << d << std::endl;

}




