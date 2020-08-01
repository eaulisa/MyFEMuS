#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <ctime>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include </usr/include/eigen3/Eigen/src/Core/util/DisableStupidWarnings.h>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <fstream>
#include <cmath>
#include "Marker.hpp"
#include "Line.hpp"
#include <algorithm>



void InitRectangleParticle(const unsigned &dim, const double &L, const double &H, const double &Lf,
                           const double &dL, const double &DB, const unsigned &nbl,
                           const std::vector < double> &xc,
                           std::vector < MarkerType > &markerType,
                           std::vector < std::vector <double> > &xp,
                           std::vector <double> &wp,
                           std::vector <double> &dist) {


  double L0 = L - DB;
  double H0 = H - 0.5 * DB;

  unsigned cols0 = ceil(L0 / dL);
  unsigned rows0 = ceil(H0 / dL);
  double dx0 = L0 / cols0;
  double dy0 = H0  / rows0;

  double dbl = DB / nbl;
  unsigned colsl = ceil(L0 / dbl);
  unsigned rowsl = ceil(H0 / dbl);
  double dxl = L0 / colsl;
  double dyl = H0 / rowsl;

  double L1 = L + DB;
  double H1 = H + 0.5 * DB;
  double DH = 0.5 * (Lf - L);

  unsigned cols1 = ceil(L1 / dL);    // number of cols above the inner bulk + interface
  unsigned rows1 = ceil(H1 / dL);    // number of rows upto tops
  unsigned nDH = ceil((DH - 0.5 * DB) / dL);       //number of side cols without tops

  double dx1 = L1 / cols1;
  double dy1 = H1 / rows1;
  double dH = (DH - 0.5 * DB) / nDH;

  unsigned size0 = rows0 * cols0;                                       // number of markers on the inner bulk
  unsigned sizel = 2 * (nbl * rowsl) + 2 * (nbl * nbl) + (nbl * colsl); // number of markers on the boundary layer
  unsigned size1 = 2 * (nDH * rows1) + 2 * (nDH * nDH) + (nDH * cols1); // number of markers on the outer shell
  unsigned sizeAll = size0 + sizel + size1;

  xp.resize(sizeAll);
  wp.resize(sizeAll);
  dist.resize(sizeAll);
  unsigned cnt = 0;
  std::vector<double> XP(dim);

//inner bulk
  std::vector<double> d(3); // d = {x-xc[0], L-x, H - y  }
  for(unsigned i = 0; i < rows0 ; i++) {
    for(unsigned j = 0; j < cols0; j++) {
      XP[0] = (xc[0] + 0.5 * DB + 0.5 * dx0) + j * dx0;
      XP[1] = (xc[1] + 0.5 * dy0) + i * dy0;
      xp[cnt] = XP;
      wp[cnt] = dx0 * dy0;
      d = {xp[cnt][0] - xc[0], (xc[0] + L) - xp[cnt][0], (xc[1] + H) - xp[cnt][1] };
      dist[cnt] = *std::min_element(d.begin(), d.end());
      cnt++;
    }
  }

//left and right boundary
  for(unsigned k = 0; k < nbl; k++) {
    for(unsigned j = 0; j < rowsl; j++) {
      XP[0] = (xc[0] - 0.5 * DB + 0.5 * dbl) + k * dbl;
      XP[1] = (xc[1] + 0.5 * dyl) + j * dyl;
      xp[cnt] = XP;
      wp[cnt] = dbl * dyl;
      dist[cnt] =  xp[cnt][0] - xc[0];
      cnt++;

      XP[0] = (xc[0] + L - 0.5 * DB + 0.5 * dbl) + k * dbl;
      XP[1] = (xc[1] + 0.5 * dyl) + j * dyl;
      xp[cnt] = XP;
      wp[cnt] = dbl * dyl;
      dist[cnt] += (xc[0] + L) - xp[cnt][0];
      cnt++;

    }
  }

  //top band without corners
  for(unsigned k = 0; k < nbl; k++) {
    for(unsigned j = 0; j < colsl; j++) {
      XP[0] = (xc[0] + 0.5 * DB + 0.5 * dxl) + j * dxl;
      XP[1] = (xc[1] + H0  + 0.5 * dbl) + k * dbl;
      xp[cnt] = XP;
      wp[cnt] = dbl * dxl;
      dist[cnt] = (xc[1] + H) - xp[cnt][1];
      cnt++;
    }
  }

//corner chuncks
  for(unsigned k = 0; k < nbl; k++) {
    for(unsigned j = 0; j < nbl; j++) {
      XP[0] = (xc[0] - 0.5 * DB + 0.5 * dbl) + k * dbl;
      XP[1] = (xc[1] + H0  + 0.5 * dbl) + j * dbl;
      xp[cnt] = XP;
      wp[cnt] = dbl * dbl;

      if(xp[cnt][0] > xc[0]) { //right
        if(xp[cnt][1] < xc[1] + H) { //bottom
          d = {xp[cnt][0] - xc[0], (xc[1] + H) - xp[cnt][1]};
          dist[cnt] = *std::min_element(d.begin(), d.end());
        }
        else { //top + interface
          dist[cnt] = (xc[1] + H) - xp[cnt][1];
        }
      }
      else { //left + interface
        if(xp[cnt][1] < H + xc[1]) { //bottom
          dist[cnt] = xc[0] - xp[cnt][0];
        }
        else { //top +interface
          dist[cnt] = -sqrt(pow(xp[cnt][0] - xc[0], 2) + pow(xp[cnt][1] - (xc[1] + H), 2));
        }
      }

      cnt++;

      XP[0] = (xc[0] + L - 0.5 * DB + 0.5 * dbl) + k * dbl;
      XP[1] = (xc[1] + H0  + 0.5 * dbl) + j * dbl;
      xp[cnt] = XP;
      wp[cnt] = dbl * dbl;


      if(xp[cnt][0] < xc[0] + L) { //left
        if(xp[cnt][1] < xc[1] + H) { //bottom
          d = {(xc[0] + L) - xp[cnt][0], (xc[1] + H) - xp[cnt][1]};
          dist[cnt] = *std::min_element(d.begin(), d.end());
        }
        else { //top + interface
          dist[cnt] = (xc[1] + H) - xp[cnt][1];
        }
      }
      else { //right + interface
        if(xp[cnt][1] < xc[1] + H) { //bottom
          dist[cnt] = (xc[0] + L) - xp[cnt][0];
        }
        else { //top +interface
          dist[cnt] = -sqrt(pow(xp[cnt][0] - (xc[0] + L), 2) + pow(xp[cnt][1] - (xc[1] + H), 2));
        }
      }

      cnt++;
    }

  }


  ////////////////OUTER SHELL

  //left and right outer
  for(unsigned k = 0; k < nDH; k++) {
    for(unsigned j = 0; j < rows1; j++) {
      XP[0] = (xc[0] - DH + 0.5 * dH) + k * dH;
      XP[1] = (xc[1] + 0.5 * dy1) + j * dy1;
      xp[cnt] = XP;
      wp[cnt] = dH * dy1;
      if(xp[cnt][1] < xc[1] + H) {
        dist[cnt] =  xp[cnt][0] - xc[0];
      }
      else {
        dist[cnt] = -sqrt(pow(xp[cnt][0] - xc[0], 2) + pow(xp[cnt][1] - (xc[1] + H), 2));
      }
      cnt++;

      XP[0] = (xc[0] + L + 0.5 * DB + 0.5 * dH) + k * dH;
      XP[1] = (xc[1] + 0.5 * dy1) + j * dy1;
      xp[cnt] = XP;
      wp[cnt] = dH * dy1;
      if(xp[cnt][1] < xc[1] + H) {
        dist[cnt] = (xc[0] + L) - xp[cnt][0];
      }
      else {
        dist[cnt] = -sqrt(pow(xp[cnt][0] - (xc[0] + L), 2) + pow(xp[cnt][1] - (xc[1] + H), 2));
      }
      cnt++;
    }
  }


  //top band without corners
  for(unsigned k = 0; k < nDH; k++) {
    for(unsigned j = 0; j < cols1; j++) {
      XP[0] = (xc[0] + 0.5 * dx1) + j * dx1;
      XP[1] = (xc[1] + H1  + 0.5 * dH) + k * dH;
      xp[cnt] = XP;
      wp[cnt] = dx1 * dH;
      dist[cnt] = (xc[1] + H) - xp[cnt][1];
      cnt++;
    }
  }


  //top two corners
  for(unsigned k = 0; k < nDH; k++) {
    for(unsigned j = 0; j < nDH; j++) {
      XP[0] = (xc[0] - DH + 0.5 * dH) + k * dH;
      XP[1] = (xc[1] + H1 + 0.5 * dH) + j * dH;
      xp[cnt] = XP;
      wp[cnt] = dH * dH;
      dist[cnt] = -sqrt(pow(xp[cnt][0] - xc[0], 2) + pow(xp[cnt][1] - (xc[1] + H), 2));
      cnt++;

      XP[0] = (xc[0] + L + 0.5 * DB + 0.5 * dH) + k * dH;
      XP[1] = (xc[1] + H1 + 0.5 * dH) + j * dH;
      xp[cnt] = XP;
      wp[cnt] = dH * dH;
      wp[cnt] = dH * dH;
      dist[cnt] = -sqrt(pow(xp[cnt][0] - (xc[0] + L), 2) + pow(xp[cnt][1] - (xc[1] + H), 2));
      cnt++;
    }
  }





  double sum = 0.;
  for(unsigned j = 0; j < xp.size(); j++) {
    sum += wp[j];
  }
  std::cout << "Volume difference = " << sum << " " << sum - (H + DH)*Lf << std::endl;
  //
  //   //could not fix
  //   double area;
  //   (dbl == 1) ? area = L * H : area = L * H + 2 * 0.5 * dL * H + (0.5 * dL) * (L + 2 * 0.5 * dL);

  //   std::setprecision(6);
  //   std::cout << " ExactArea = " << area << " ComputedArea = " << sum << std::endl;



  markerType.assign(cnt, VOLUME);


}






void InitRectangleInterface(const unsigned & dim, const double &L, const double &H, const double &DB, const unsigned nbl,
                            const unsigned & FI, const std::vector < double> &xc, std::vector < MarkerType > &markerType,
                            std::vector < std::vector <double> > &xp,
                            std::vector < std::vector < std::vector < double > > > &T) {


  double L0 = L - DB;
  double H0 = H - 0.5 * DB;



  double dbl = DB / nbl;
  unsigned colsl = ceil(L0 / dbl);
  unsigned rowsl = ceil(H0 / dbl);
  double dxl = L0 / colsl;
  double dyl = H0 / rowsl;

  //cols0 = FI * ceil(L / dx0);
  //rows0 = FI * ceil(H / dy0 );



  unsigned size = (2 * rowsl + colsl);// +  2 * nbl;
  xp.resize(size);


  for(unsigned i = 0; i < size; i++) {
    xp[i].assign(dim, 0.);
  }

  T.resize(size);
  for(unsigned i = 0; i < size; i++) {
    T[i].resize(dim - 1);
    for(unsigned k = 0; k < dim - 1; k++) {
      T[i][k].resize(dim, 0.);
    }
  }


  unsigned cnt = 0;

//left and right boundary
  std::vector<double> XP(dim, 0.);

  for(unsigned j = 0; j < rowsl; j++) {

    XP[0] = xc[0];
    XP[1] = (xc[1] + 0.5 * dyl) + j * dyl;
    xp[cnt] = XP;
    //std::cout << xp[cnt][0] << " " << xp[cnt][1] << std::endl;

    T[cnt][0][0] = 0.;
    T[cnt][0][1] = -dyl;
    //std::cout << T[cnt][0][0] << " " << T[cnt][0][1] << std::endl;
    cnt++;

    XP[0] = L + xc[0];
    XP[1] = (xc[1] + 0.5 * dyl) + j * dyl;
    xp[cnt] = XP;
    //std::cout << xp[cnt][0] << " " << xp[cnt][1] << std::endl;

    T[cnt][0][0] = 0.;
    T[cnt][0][1] = dyl;

    //std::cout << T[cnt][0][0] << " " << T[cnt][0][1] << std::endl;
    cnt++;
  }


  //top boundary
  for(unsigned j = 0; j < colsl; j++) {
    XP[0] = xc[0] + 0.5 * DB + 0.5 * dxl + j * dxl;
    XP[1] = H + xc[1];
    xp[cnt] = XP;
    //std::cout << xp[cnt][0] << " " << xp[cnt][1] << std::endl;

    T[cnt][0][0] = -dxl;
    T[cnt][0][1] = 0.;
    //std::cout << T[cnt][0][0] << " " << T[cnt][0][1] << std::endl;
    cnt++;
  }




//   // corner left
//   for(unsigned j = 0; j < nbl / 2; j++) {
//     XP[0] = xc[0] ;
//     XP[1] = (H0 + xc[1] + 0.5 * dbl) + j * dbl;
//     xp[cnt] = XP;
//     T[cnt][0][0] = 0.;
//     T[cnt][0][1] = -dbl;
//     cnt++;
//   }
//   XP[0] = xc[0] ;
//   XP[1] = H + xc[1];
//   xp[cnt] = XP;
//   T[cnt][0][0] = -dbl * acos(M_PI / 4.);
//   T[cnt][0][1] = -dbl * acos(M_PI / 4.);
//   cnt++;
//   for(unsigned j = 0; j < nbl / 2; j++) {
//     XP[0] = xc[0] + dbl + j * dbl ;
//     XP[1] = xc[1] + H ;
//     xp[cnt] = XP;
//     T[cnt][0][0] = -dbl;
//     T[cnt][0][1] = 0.;
//     cnt++;
//   }
// 
//   //corner right
//   for(unsigned j = 0; j < nbl / 2; j++) {
//     XP[0] = xc[0] + L ;
//     XP[1] = (H0 + xc[1] + 0.5 * dbl) + j * dbl;
//     xp[cnt] = XP;
//     T[cnt][0][0] = 0.;
//     T[cnt][0][1] = dbl;
//     cnt++;
//   }
// 
//   XP[0] = xc[0] + L;
//   XP[1] = H + xc[1];
//   xp[cnt] = XP;
//   T[cnt][0][0] = -dbl * acos(M_PI / 4.);
//   T[cnt][0][1] = dbl * acos(M_PI / 4.);
//   cnt++;
// 
// 
// 
//   for(unsigned j = 0; j < nbl / 2; j++) {
// 
//     XP[0] = xc[0] + L - dbl - j * dbl ;
//     XP[1] = xc[1] + H ;
//     xp[cnt] = XP;
// 
//     T[cnt][0][0] = -dbl;
//     T[cnt][0][1] = 0.;
//     cnt++;
// 
//   }



  std::cout << "size = " << size << " cnt = " << cnt << std::endl;

  markerType.assign(size, INTERFACE);

}




void InitBallVolumeParticles(const unsigned & dim, std::vector<double> &VxL, std::vector<double> &VxR,
                             const std::vector < double> &xc, std::vector < MarkerType > &markerType, const double & R, const double & Rmax, const double & DR, const unsigned & nbl, const unsigned & FI,
                             std::vector < std::vector <double> > &xp, std::vector <double> &wp, std::vector <double> &dist) {



  double theta0 = 0;
  double theta1 = 2 * M_PI;

  double phi1 = M_PI;
  double phi0 = 0;

  double R0 = 0. ;

  double dp = DR;


  //std::cout << "theta0 = " << theta0 << " theta1 = " << theta1 << " phi0 = " << phi0 << " phi1 = " << phi1 <<  " R0 = " << R0 <<  " R1 = " << R1 << std::endl;


  unsigned nr = ceil(((R - 0.5 * dp)) / dp);
  double dr = ((R - 0.5 * dp)) /  nr;
  double area = 0.;

  xp.reserve(pow(4 * nr, dim));
  wp.reserve(pow(4 * nr, dim));
  dist.reserve(pow(4 * nr, dim));
  unsigned cnt = 0;

  int i0 = floor(R0 / dr - 0.5);
  if(i0 < 0) i0 = 0;
  //int i1 = ceil(R1 / dr - 0.5);
  int i1 = ceil((R - 0.5 * dp) / dr - 0.5);
  if(i1 > nr - 1) i1 = nr - 1;

  std::vector<double> XP(dim);

  for(int i = i0; i <= i1; i++) {
    double ri = (i + 0.5) * dr;
    if(dim == 2) {

      unsigned nti = ceil((theta1 - theta0) * ri / dr);
      double dti = (theta1 - theta0) / nti;

      for(unsigned j = 0; j < nti; j++) {
        double tj = theta0 + (0.5 + j) * dti;
        XP[0] = xc[0] + ri * cos(tj);
        XP[1] = xc[1] + ri * sin(tj);
        unsigned flag  = 0;
        for(unsigned k = 0; k < dim; k++) {
          if(VxL[k] > XP[k] || XP[k] > VxR[k]) {
            flag++;
            break;
          }
        }
        if(!flag) {
          xp.resize(cnt + 1);
          xp[cnt] = XP;
          //std::cout << xp[cnt][0] << " " << xp[cnt][1] << std::endl;
          wp.resize(cnt + 1);
          dist.resize(cnt + 1);
          wp[cnt] = ri * dti * dr;
          dist[cnt] = (R - ri);

          area += ri * dti * dr;
          cnt++;
        }
      }
    }
    else {

      unsigned nphi = ceil((phi1 - phi0) * ri / dr);
      double dphi = (phi1 - phi0) / nphi;

      for(unsigned k = 0; k < nphi; k++) {
        double pk = phi0 + (0.5 + k) * dphi;

        unsigned nti = ceil((theta1 - theta0) * ri * sin(pk) / dr);
        double dti = (theta1 - theta0) / nti;

        for(unsigned j = 0; j < nti; j++) {
          double tj = theta0 + (0.5 + j) * dti;

          XP[0] = xc[0] + ri * sin(pk) * cos(tj);
          XP[1] = xc[1] + ri * sin(pk) * sin(tj);
          XP[2] = xc[2] + ri * cos(pk);

          unsigned flag  = 0;
          for(unsigned k = 0; k < dim; k++) {
            if(VxL[k] > XP[k] || XP[k] > VxR[k]) {
              flag++;
              break;
            }
          }

          if(!flag) {

            xp.resize(cnt + 1);
            xp[cnt] = XP;

            wp.resize(cnt + 1);
            dist.resize(cnt + 1);
            wp[cnt] = dr * (ri * dphi) * (ri * sin(pk) * dti);
            dist[cnt] = (R - ri);

            area += dr * (ri * dphi) * (ri * sin(pk) * dti);  // fix the volume
            cnt++;
          }
        }
      }
    }
  }

  //std::cout << xp[0].size() << " " << Rmax << "\n";


  double dbl = dp / nbl;
  for(unsigned i = 0; i < nbl; i++) {
    double ri = (R - 0.5 * dp) + (i + 0.5) * dbl;
    if(dim == 2) {
      unsigned nti = FI * ceil((theta1 - theta0) * R / dr);
      double dti = (theta1 - theta0) / nti;
      for(unsigned j = 0; j < nti; j++) {
        double tj = theta0 + (0.5 + j) * dti;
        XP[0] = xc[0] + ri * cos(tj);
        XP[1] = xc[1] + ri * sin(tj);
        unsigned flag  = 0;
        for(unsigned k = 0; k < dim; k++) {
          if(VxL[k] > XP[k]  || XP[k] > VxR[k]) {
            ++flag;
            break;
          }
        }
        if(!flag) {

          xp.resize(cnt + 1);
          xp[cnt] = XP;

          //std::cout << xp[cnt][0] << " " << xp[cnt][1] << std::endl;

          wp.resize(cnt + 1);
          dist.resize(cnt + 1);
          wp[cnt] = ri * dti * dbl;
          dist[cnt] = (R - ri);

          area += ri * dti * dbl;
          cnt++;
        }
      }
    }
    else {
      unsigned nphi = FI * ceil((phi1 - phi0) * ri / dr);
      double dphi = (phi1 - phi0) / nphi;

      for(unsigned k = 0; k < nphi; k++) {
        double pk = phi0 + (0.5 + k) * dphi;

        unsigned nti = FI * ceil((theta1 - theta0) * ri * sin(pk) / dr);
        double dti = (theta1 - theta0) / nti;

        for(unsigned j = 0; j < nti; j++) {
          double tj = theta0 + (0.5 + j) * dti;
          XP[0] = xc[0] + ri * sin(pk) * cos(tj);
          XP[1] = xc[1] + ri * sin(pk) * sin(tj);
          XP[2] = xc[2] + ri * cos(pk);

          unsigned flag  = 0;
          for(unsigned k = 0; k < dim; k++) {
            if(VxL[k] > XP[k] || XP[k] > VxR[k]) {
              flag++;
              break;
            }
          }

          if(!flag) {

            xp.resize(cnt + 1);
            xp[cnt] = XP;

            wp.resize(cnt + 1);
            dist.resize(cnt + 1);
            wp[cnt] = dbl * (ri * dphi) * (ri * sin(pk) * dti);
            dist[cnt] = (R - ri);

            area += dbl * (ri * dphi) * (ri * sin(pk) * dti);
            cnt++;
          }
        }

      }
    }
  }



  i0 = floor((R0 - (R + 0.5 * dp)) / dr - 0.5);
  if(i0 < 0) i0 = 0;
  i1 = floor((Rmax - (R + 0.5 * dp)) / dr - 0.5);

  for(int i = i0; i <= i1; i++) {

    double ri = (R + 0.5 * (dp + dr)) + i * dr;

    if(dim == 2) {
      unsigned nti = ceil((theta1 - theta0) * ri / dr);
      double dti = (theta1 - theta0) / nti;
      for(unsigned j = 0; j < nti; j++) {
        double tj = theta0 + (0.5 + j) * dti;
        XP[0] = xc[0] + ri * cos(tj);
        XP[1] = xc[1] + ri * sin(tj);
        //std::cout << XP[0] << " " << XP[1] << std::endl;
        unsigned flag  = 0;
        for(unsigned k = 0; k < dim; k++) {
          if(VxL[k] > XP[k]  || XP[k] > VxR[k]) {
            ++flag;
            break;
          }
        }
        if(!flag) {
          xp.resize(cnt + 1);
          xp[cnt] = XP;
          //std::cout << xp[cnt][0] << " " << xp[cnt][1] << std::endl;
          wp.resize(cnt + 1);
          dist.resize(cnt + 1);
          wp[cnt] = ri * dti * dr;
          dist[cnt] = (R - ri);

          area += ri * dti * dr;
          cnt++;

        }
      }
    }
    else {

      unsigned nphi = ceil((phi1 - phi0) * ri / dr);
      double dphi = (phi1 - phi0) / nphi;

      for(unsigned k = 0; k < nphi; k++) {
        double pk = phi0 + (0.5 + k) * dphi;


        unsigned nti = ceil((theta1 - theta0) * ri * sin(pk) / dr);
        double dti = (theta1 - theta0) / nti;

        for(unsigned j = 0; j < nti; j++) {
          double tj = theta0 + (0.5 + j) * dti;
          XP[0] = xc[0] + ri * sin(pk) * cos(tj);
          XP[1] = xc[1] + ri * sin(pk) * sin(tj);
          XP[2] = xc[2] + ri * cos(pk);

          unsigned flag  = 0;
          for(unsigned k = 0; k < dim; k++) {
            if(VxL[k] > XP[k] || XP[k] > VxR[k]) {
              flag++;
              break;
            }
          }

          if(!flag) {
            xp.resize(cnt + 1);
            xp[cnt] = XP;

            wp.resize(cnt + 1);
            dist.resize(cnt + 1);
            wp[cnt] = dr * (ri * dphi) * (ri * sin(pk) * dti);
            dist[cnt] = (R - ri);

            area += dr * (ri * dphi) * (ri * sin(pk) * dti);  // fix the volume
            cnt++;
          }
        }

      }
    }
  }
  if(dim == 2) {
    std::cout << "computed area = " << area << " " << " analytic area = " << M_PI * Rmax * Rmax << std::endl;
  }
  else {
    std::cout << "computed volume = " << area << " " << " analytic volume = " << 4. / 3. * M_PI * Rmax * Rmax * Rmax << std::endl;
  }

  markerType.assign(cnt, VOLUME);

}



void InitBallInterfaceParticles(const unsigned & dim, const double & R, const double & DR, const unsigned & FI, const std::vector < double> &xc, std::vector < MarkerType > &markerType, std::vector < std::vector <double> > &xp, std::vector < std::vector < std::vector < double > > > &T) {


  unsigned nr = ceil(((R - 0.5 * DR)) / DR);
  double dr = ((R - 0.5 * DR)) / nr;

  if(dim == 2) {

    unsigned Ntheta = FI * ceil(2. * M_PI * R / dr);

    xp.resize(Ntheta);
    markerType.assign(Ntheta, INTERFACE);

    for(unsigned i = 0; i < Ntheta; i++) {
      xp[i].assign(dim, 0.);
    }

    T.resize(Ntheta);
    for(unsigned i = 0; i < Ntheta; i++) {
      T[i].resize(dim - 1);
      for(unsigned k = 0; k < dim - 1; k++) {
        T[i][k].resize(dim, 0.);
      }
    }

    double arcLenght = 2. * M_PI * R / Ntheta;
    double dtheta = 2 * M_PI / Ntheta;


    for(unsigned i = 0; i < Ntheta; i++) {

      double ti = 0. + (FI * 0.5 + i) * dtheta;

      xp[i][0] = xc[0] + R * cos(ti);
      xp[i][1] = xc[1] + R * sin(ti);

      T[i][0][0] = -arcLenght * sin(ti);
      T[i][0][1] = arcLenght * cos(ti);

    }
    markerType.assign(Ntheta, INTERFACE);
  }
  else {

    std::vector<double> XP(dim);
    unsigned nphi = FI * ceil(M_PI * R / dr);
    double dphi = M_PI / nphi;

    xp.reserve(nphi * 2 * nphi);
    T.reserve(nphi * 2 * nphi);

    unsigned cnt = 0;

    for(unsigned k = 0; k < nphi; k++) {
      double pk = (0.5 + k) * dphi;

      unsigned nti = FI * ceil(2 * M_PI * R * sin(pk) / dr);
      double dti = 2 * M_PI / nti;

      for(unsigned j = 0; j < nti; j++) {
        double tj = (0.5 + j) * dti;

        XP[0] = xc[0] + R * sin(pk) * cos(tj);
        XP[1] = xc[1] + R * sin(pk) * sin(tj);
        XP[2] = xc[2] + R * cos(pk);

        xp.resize(cnt + 1);
        xp[cnt] = XP;

        T.resize(cnt + 1);
        T[cnt].resize(2);
        T[cnt][0].resize(3);
        T[cnt][1].resize(3);

        T[cnt][0][0] = R * cos(pk) * cos(tj) * dphi;
        T[cnt][0][1] = R * cos(pk) * sin(tj) * dphi;
        T[cnt][0][2] = -R * sin(pk) * dphi;

        T[cnt][1][0] = -R * sin(pk) * sin(tj) * dti;
        T[cnt][1][1] = R * sin(pk) * cos(tj) * dti;
        T[cnt][1][2] = 0.;

        cnt++;
      }
    }
    markerType.assign(cnt, INTERFACE);
  }

}




void InitParticlesDisk3D(const unsigned & dim, const unsigned & ng, std::vector<double> &VxL, std::vector<double> &VxR,
                         const std::vector < double> &xc, const double & R, std::vector < std::vector <double> > &xp,
                         std::vector <double> &wp, std::vector <double> &dist) {
//  VxL = {xL,yL,zL};
//  VxR = {xR,yR,zR};

  unsigned m1 = ceil(2 * ng - 1);

  double theta0 = atan2(VxL[1] - xc[1], VxR[0] - xc[0]);
  double theta1 = atan2(VxR[1] - xc[1], VxL[0] - xc[0]);
  if(theta0 < 0) theta0 += M_PI;
  if(theta1 < 0) theta1 += M_PI;


  double phi1 = M_PI / 2 - atan2(VxL[2] - xc[2], sqrt((VxR[0] - xc[0]) * (VxR[0] - xc[0]) + (VxR[1] - xc[1]) * (VxR[1] - xc[1])));
  double phi0 = M_PI / 2 - atan2(VxR[2] - xc[2], sqrt((VxL[0] - xc[0]) * (VxL[0] - xc[0]) + (VxL[1] - xc[1]) * (VxL[1] - xc[1])));
  if(phi0 < 0) phi0 += M_PI;
  if(phi1 < 0) phi1 += M_PI;

  double R0 = 0. ;
  double R1 = 0. ;
  double dp = 1.;


  for(unsigned k = 0; k < dim; k++) {
    R0 += (VxL[k] - xc[k]) * (VxL[k] - xc[k]);
    R1 += (VxR[k] - xc[k]) * (VxR[k] - xc[k]);
    dp *= VxR[k] - VxL[k];
  }

  R0 = sqrt(R0);
  R1 = sqrt(R1);

  dp = pow(dp, 1. / dim) / m1;


  //std::cout << "theta0 = " << theta0 << " theta1 = " << theta1 << " phi0 = " << phi0 << " phi1 = " << phi1 <<  " R0 = " << R0 <<  " R1 = " << R1 << std::endl;


  unsigned nr = ceil(((R - 0.5 * dp)) / dp);
  double dr = ((R - 0.5 * dp)) / nr;
  double area = 0.;

  xp.resize(dim);
  for(unsigned k = 0; k < dim; k++) {
    xp[k].reserve(2 * pow(m1, dim));
  }
  wp.reserve(2 * pow(m1, dim));
  dist.reserve(2 * pow(m1, dim));
  unsigned cnt = 0;

  int i0 = floor(R0 / dr - 0.5);
  if(i0 < 0) i0 = 0;
  //int i1 = ceil(R1 / dr - 0.5);
  int i1 = ceil((R - 0.5 * dp) / dr - 0.5);
  if(i1 > nr - 1) i1 = nr - 1;

  //std::cout << "i0 = " << i0 << " i1 = " << i1 <<std::endl;

  std::vector<double> XP(dim);

  for(int i = i0; i <= i1; i++) {
    double ri = (i + 0.5) * dr;
    if(dim == 2) {
      //unsigned nti = ceil(2 * M_PI * ri / dr);
      //double dti = 2 * M_PI / nti;
      double dti = dr / ri;
      int j0 = floor(theta0 / dti);
      int j1 = ceil(theta1 / dti);
      for(unsigned j = j0; j <= j1; j++) {
        double tj = j * dti;
        XP[0] = xc[0] + ri * cos(tj);
        XP[1] = xc[1] + ri * sin(tj);
        unsigned flag  = dim;
        for(unsigned k = 0; k < dim; k++) {
          if(VxL[k] < XP[k]  && XP[k] < VxR[k]) {
            --flag;
          }
        }
        if(!flag) {

          for(unsigned k = 0; k < dim; k++) {
            xp[k].resize(cnt + 1);
          }
          wp.resize(cnt + 1);
          dist.resize(cnt + 1);

          xp[0][cnt] = XP[0];
          xp[1][cnt] = XP[1];

          //std::cout << xp[0][cnt] << " " << xp[1][cnt] << std::endl;

          wp[cnt] = ri * dti * dr;
          dist[cnt] = (R - ri);
          area += ri * dti * dr;
          cnt++;
        }
      }
    }
    else {

      double dphi = dr / ri;
      int k0 = floor(phi0 / dphi);
      int k1 = ceil(phi1 / dphi);

      for(unsigned k = k0; k <= k1; k++) {
        double pk = (k + 0.5) * dphi;

        double dti = dr / (ri * sin(pk));
        int j0 = floor(theta0 / dti);
        int j1 = ceil(theta1 / dti);

        for(unsigned j = j0; j <= j1; j++) {
          double tj = j * dti;
          XP[0] = xc[0] + ri * sin(pk) * cos(tj);
          XP[1] = xc[1] + ri * sin(pk) * sin(tj);
          XP[2] = xc[2] + ri * cos(pk);

          unsigned flag  = dim;
          for(unsigned k = 0; k < dim; k++) {
            if(VxL[k] <= XP[k]  && XP[k] <= VxR[k]) {
              --flag;
            }
          }

          if(!flag) {
            for(unsigned k = 0; k < dim; k++) {
              xp[k].resize(cnt + 1);
            }
            wp.resize(cnt + 1);
            dist.resize(cnt + 1);

            xp[0][cnt] = XP[0];
            xp[1][cnt] = XP[1];
            xp[2][cnt] = XP[2];

            //std::cout << xp[0][cnt] << " " << xp[1][cnt] << " " << xp[2][cnt] << std::endl;

            wp[cnt] = dr * (ri * dphi) * (ri * sin(pk) * dti);
            dist[cnt] = (R - ri);

            area += dr * (ri * dphi) * (ri * sin(pk) * dti);  // fix the volume
            cnt++;
          }
        }

      }
    }
  }

  double ri = R;
  if(dim == 2) {
    double dti = dr / ri;
    int j0 = floor(theta0 / dti);
    int j1 = ceil(theta1 / dti);
    for(unsigned j = j0; j <= j1; j++) {
      double tj = j * dti;
      XP[0] = xc[0] + ri * cos(tj);
      XP[1] = xc[1] + ri * sin(tj);
      unsigned flag  = dim;
      for(unsigned k = 0; k < dim; k++) {
        if(VxL[k] < XP[k]  && XP[k] < VxR[k]) {
          --flag;
        }
      }
      if(!flag) {

        for(unsigned k = 0; k < dim; k++) {
          xp[k].resize(cnt + 1);
        }
        wp.resize(cnt + 1);
        dist.resize(cnt + 1);

        xp[0][cnt] = XP[0];
        xp[1][cnt] = XP[1];

        //std::cout << xp[0][cnt] << " " << xp[1][cnt] << std::endl;

        wp[cnt] = ri * dti * dr;
        dist[cnt] = (R - ri);
        area += ri * dti * dr;
        cnt++;
      }
    }
  }
  else {
    double dphi = dr / ri;
    int k0 = floor(phi0 / dphi);
    int k1 = ceil(phi1 / dphi);

    for(unsigned k = k0; k <= k1; k++) {
      double pk = (k + 0.5) * dphi;

      double dti = dr / (ri * sin(pk));
      int j0 = floor(theta0 / dti);
      int j1 = ceil(theta1 / dti);

      for(unsigned j = j0; j <= j1; j++) {
        double tj = j * dti;
        XP[0] = xc[0] + ri * sin(pk) * cos(tj);
        XP[1] = xc[1] + ri * sin(pk) * sin(tj);
        XP[2] = xc[2] + ri * cos(pk);

        unsigned flag  = dim;
        for(unsigned k = 0; k < dim; k++) {
          if(VxL[k] <= XP[k]  && XP[k] <= VxR[k]) {
            --flag;
          }
        }

        if(!flag) {
          for(unsigned k = 0; k < dim; k++) {
            xp[k].resize(cnt + 1);
          }
          wp.resize(cnt + 1);
          dist.resize(cnt + 1);

          xp[0][cnt] = XP[0];
          xp[1][cnt] = XP[1];
          xp[2][cnt] = XP[2];

          //std::cout << xp[0][cnt] << " " << xp[1][cnt] << " " << xp[2][cnt] << std::endl;

          wp[cnt] = dp * (ri * dphi) * (ri * sin(pk) * dti);
          dist[cnt] = (R - ri);

          area += dp * (ri * dphi) * (ri * sin(pk) * dti);  // fix the volume
          cnt++;
        }
      }

    }
  }


  i0 = floor((R0 - (R + 0.5 * dp)) / dr - 0.5);
  if(i0 < 0) i0 = 0;
  i1 = floor((R1 - (R + 0.5 * dp)) / dr - 0.5);

  for(int i = i0; i <= i1; i++) {
    double ri = (R + 0.5 * (dp + dr)) + i * dr;
    if(dim == 2) {
      double dti = dr / ri;
      int j0 = floor(theta0 / dti);
      int j1 = ceil(theta1 / dti);
      for(unsigned j = j0; j <= j1; j++) {
        double tj = j * dti;
        XP[0] = xc[0] + ri * cos(tj);
        XP[1] = xc[1] + ri * sin(tj);
        unsigned flag  = dim;
        for(unsigned k = 0; k < dim; k++) {
          if(VxL[k] < XP[k]  && XP[k] < VxR[k]) {
            --flag;
          }
        }
        if(!flag) {

          for(unsigned k = 0; k < dim; k++) {
            xp[k].resize(cnt + 1);
          }
          wp.resize(cnt + 1);
          dist.resize(cnt + 1);

          xp[0][cnt] = XP[0];
          xp[1][cnt] = XP[1];

          //std::cout << xp[0][cnt] << " " << xp[1][cnt] << std::endl;

          wp[cnt] = ri * dti * dr;
          dist[cnt] = (R - ri);

          area += ri * dti * dr;
          cnt++;
        }
      }
    }
    else {

      double dphi = dr / ri;
      int k0 = floor(phi0 / dphi);
      int k1 = ceil(phi1 / dphi);

      for(unsigned k = k0; k <= k1; k++) {
        double pk = (k + 0.5) * dphi;
        //unsigned nti = ceil(2 * M_PI * pk / dphi);
        //double dti =  M_PI / nti;
        double dti = dr / (ri * sin(pk));
        int j0 = floor(theta0 / dti);
        int j1 = ceil(theta1 / dti);

        for(unsigned j = j0; j <= j1; j++) {
          double tj = j * dti;
          XP[0] = xc[0] + ri * sin(pk) * cos(tj);
          XP[1] = xc[1] + ri * sin(pk) * sin(tj);
          XP[2] = xc[2] + ri * cos(pk);

          unsigned flag  = dim;
          for(unsigned k = 0; k < dim; k++) {
            if(VxL[k] <= XP[k]  && XP[k] <= VxR[k]) {
              --flag;
            }
          }

          if(!flag) {
            for(unsigned k = 0; k < dim; k++) {
              xp[k].resize(cnt + 1);
            }
            wp.resize(cnt + 1);
            dist.resize(cnt + 1);

            xp[0][cnt] = XP[0];
            xp[1][cnt] = XP[1];
            xp[2][cnt] = XP[2];

            //std::cout << xp[0][cnt] << " " << xp[1][cnt] << " " << xp[2][cnt] << std::endl;

            wp[cnt] = dr * (ri * dphi) * (ri * sin(pk) * dti);
            dist[cnt] = (R - ri);

            area += dr * (ri * dphi) * (ri * sin(pk) * dti);  // fix the volume
            cnt++;
          }
        }
      }
    }
  }

}



void Testing(double & a, double & b, const unsigned & m, const unsigned & dim, Eigen::MatrixXd & x,
             Eigen::VectorXd & w_new, std::vector<double> &dist, const double & eps, double & QuadSum, double & IntSum) {


  double deps = (b - a) * eps; // eps1

  double a0 = 0.5; // 128./256.;
  double a1 = pow(deps, -1.) * 1.23046875; // 315/256.;
  double a3 = -pow(deps, -3.) * 1.640625; //420./256.;
  double a5 = pow(deps, -5.) * 1.4765625; // 378./256.;
  double a7 = -pow(deps, -7.) * 0.703125; // 180./256.;
  double a9 = pow(deps, -9.) * 0.13671875; // 35./256.;

  QuadSum = 0.;
  IntSum = 0.;

  for(unsigned i = 0; i < w_new.size(); i++) {

    double dg1 = dist[i];
    double dg2 = dg1 * dg1;
    double xi;
    if(dg1 < -deps)
      xi = 0.;
    else if(dg1 > deps) {
      xi = 1.;
    }
    else {
      xi = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
    }

    double r = 0.75 * 0.75 ;
    for(unsigned k = 0; k < dim; k++) {
      r +=  - x(k, i) * x(k, i);
    }
    IntSum += xi * r * w_new(i);

    r = 1.;
    for(unsigned k = 0; k < dim; k++) {
      r *=  x(k, i);
    }
    QuadSum += pow(r, m) * w_new(i);
  }

}


void Cheb(const unsigned & m, Eigen::VectorXd & xg, Eigen::MatrixXd & C) {

  C.resize(xg.size(), m + 1);
  for(unsigned i = 0; i < xg.size(); i++) {
    C(i, 0) = 1;
    C(i, 1) = xg(i);
    for(unsigned j = 2; j <= m; j++) {
      C(i, j) =  2 * xg(i) * C(i, j - 1) - C(i, j - 2);
    }
  }
  C.transposeInPlace();

}



void AssembleMatEigen(std::vector<double>& VxL, std::vector<double> &VxR, const unsigned & m, const unsigned & dim, const unsigned & np, Eigen::Tensor<double, 3, Eigen::RowMajor>  &PmX, Eigen::MatrixXd & Pg,  Eigen::VectorXd & wg, Eigen::MatrixXd & A, Eigen::VectorXd & F) {


  A.resize(pow(m + 1, dim), np);
  F.resize(pow(m + 1, dim));
  Eigen::VectorXi I(dim);
  Eigen::VectorXi N(dim);


  for(unsigned k = 0; k < dim ; k++) {
    N(k) = pow(m + 1, dim - k - 1);
  }

  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    I(0) = t / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = t % N(k - 1);
      I(k) = pk / N(k); // dimensional index over on the space of polynomaials
    }
    for(unsigned j = 0; j < np; j++) {
      double r = 1;

      for(unsigned k = 0; k < dim; k++) {
        r *= PmX(k, I[k], j);
      }
      A(t, j) = r ;
    }

  }

  unsigned ng = Pg.row(0).size();
  Eigen::VectorXi J(dim);
  Eigen::VectorXi NG(dim);



  for(unsigned k = 0; k < dim ; k++) {
    NG(k) = pow(ng, dim - k - 1);
  }

  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    I(0) = t / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = t % N(k - 1);
      I(k) = pk / N(k); // dimensional index over on the space of polynomaials
    }
    F(t) = 0.;
    for(unsigned g = 0; g < pow(ng, dim) ; g++) { // multidimensional index on the space of polynomaials
      J(0) = g / NG(0);
      for(unsigned k = 1; k < dim ; k++) {
        unsigned pk = g % NG(k - 1);
        J(k) = pk / NG(k); // dimensional index over on the space of polynomaials
      }
      double value = 1.;

      for(unsigned k = 0; k < dim ; k++) {
        value *= 0.5 * (VxR[k] - VxL[k]) * Pg(I(k), J(k)) * wg(J(k)) ;
      }
      F(t) += value;
    }

  }

}




void GetChebGaussF(const unsigned & dim, const unsigned & m, std::vector<double> &VxL, std::vector<double> &VxU, Eigen::MatrixXd & Pg,  Eigen::VectorXd & wg, Eigen::VectorXd & F) {

  F.resize(pow(m + 1, dim));
  F.setZero();
  Eigen::VectorXi I(dim);
  Eigen::VectorXi N(dim);
  Eigen::VectorXi J(dim);
  Eigen::VectorXi NG(dim);

  unsigned ng = Pg.row(0).size();

  for(unsigned k = 0; k < dim ; k++) {
    N(k) = pow(m + 1, dim - k - 1);
  }

  for(unsigned k = 0; k < dim ; k++) {
    NG(k) = pow(ng, dim - k - 1);
  }

  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    I(0) = t / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = t % N(k - 1);
      I(k) = pk / N(k);
    }
    F(t) = 0.;
    for(unsigned g = 0; g < pow(ng, dim) ; g++) { // gauss loop
      J(0) = g / NG(0);
      for(unsigned k = 1; k < dim ; k++) {
        unsigned pk = g % NG(k - 1);
        J(k) = pk / NG(k);
      }
      double value = 1.;
      unsigned ig = 0;
      double jac = 1.;

      for(unsigned k = 0; k < dim ; k++) {
        value *= 0.5 * (VxU[k] - VxL[k]) * Pg(I(k), J(k)) * wg(J(k));
      }
      F(t) += value;
    }
  }

}


void GetChebGaussF(const unsigned & dim, const unsigned & m, const std::vector<double> &jac, Eigen::MatrixXd & Pg,  Eigen::VectorXd & wg, Eigen::VectorXd & F) {

  F.resize(pow(m + 1, dim));
  F.setZero();
  Eigen::VectorXi I(dim);
  Eigen::VectorXi N(dim);
  Eigen::VectorXi J(dim);
  Eigen::VectorXi NG(dim);

  unsigned ng = Pg.row(0).size();

  for(unsigned k = 0; k < dim ; k++) {
    N(k) = pow(m + 1, dim - k - 1);
  }

  for(unsigned k = 0; k < dim ; k++) {
    NG(k) = pow(ng, dim - k - 1);
  }

  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    I(0) = t / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = t % N(k - 1);
      I(k) = pk / N(k);
    }
    F(t) = 0.;
    for(unsigned g = 0; g < pow(ng, dim) ; g++) { // gauss loop
      J(0) = g / NG(0);
      for(unsigned k = 1; k < dim ; k++) {
        unsigned pk = g % NG(k - 1);
        J(k) = pk / NG(k);
      }

      double value = jac[g];
      for(unsigned k = 0; k < dim ; k++) {
        value *= Pg(I(k), J(k)) * wg(J(k));
      }

      F(t) += value;
    }
  }

}





void GetChebXInfo(const unsigned & m, const unsigned & dim, const unsigned & np, Eigen::MatrixXd & xL, Eigen::Tensor<double, 3, Eigen::RowMajor>& PmX) {
  // xL is taken in reference coordinate system
  PmX.resize(dim, m + 1, np);
  Eigen::MatrixXd Ptemp;
  Eigen::VectorXd xtemp;
  for(unsigned k = 0; k < dim; k++) {
    xtemp = xL.row(k);
    Cheb(m, xtemp, Ptemp);
    for(unsigned i = 0; i < m + 1; i++) {
      for(unsigned j = 0; j < np; j++) {
        PmX(k, i, j) = Ptemp(i, j);
      }
    }
  }
}



void GetMultiDimChebMatrix(const unsigned & dim, const unsigned & m, const unsigned & np, Eigen::Tensor<double, 3, Eigen::RowMajor>  &PmX, Eigen::MatrixXd & A) {


  A.resize(pow(m + 1, dim), np);
  Eigen::VectorXi I(dim);
  Eigen::VectorXi N(dim);


  for(unsigned k = 0; k < dim ; k++) {
    N(k) = pow(m + 1, dim - k - 1);
  }

  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    I(0) = t / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = t % N(k - 1);
      I(k) = pk / N(k); // dimensional index over on the space of polynomaials
    }
    for(unsigned j = 0; j < np; j++) {
      double r = 1;

      for(unsigned k = 0; k < dim; k++) {
        r *= PmX(k, I[k], j);
      }
      A(t, j) = r ;
    }
  }


}




void SolWeightEigen(Eigen::MatrixXd & A, Eigen::VectorXd & F, Eigen::VectorXd & wP, Eigen::VectorXd & w_new) {

  w_new.resize(wP.size());

  Eigen::VectorXd y = A.transpose() * (A * A.transpose()).partialPivLu().solve(F - A * wP);
  w_new = y + wP;


}

void PrintMarkers(const unsigned & dim, const Eigen::MatrixXd & xP, const std::vector <double> &dist,
                  const Eigen::VectorXd wP, const Eigen::VectorXd & w_new, const unsigned & l, const unsigned & t) {

  std::ofstream fout;

  char filename[100];
  sprintf(filename, "marker%d.txt", l);

  if(t == 0) {
    fout.open(filename);
  }
  else {
    fout.open(filename, std::ios::app);
  }

  for(unsigned i = 0; i < w_new.size(); i++) {

    for(unsigned k = 0; k < dim; k++) {
      fout << xP(k, i) << " ";
    }
    fout <<  dist[i] << " " << wP[i] << " " << w_new[i] << std::endl;
  }

  fout.close();

}


void  GetParticlesOnBox(const double & a, const double & b, const unsigned & n1, const unsigned & dim, Eigen::MatrixXd & x, Eigen::MatrixXd & xL) {

  double h = (b - a) / n1;
  x.resize(dim, pow(n1, dim));
  Eigen::VectorXi I(dim);
  Eigen::VectorXi N(dim);

  for(unsigned k = 0; k < dim ; k++) {
    N(k) = pow(n1, dim - k - 1);
  }

  for(unsigned p = 0; p < pow(n1, dim) ; p++) {
    I(0) = 1 + p / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = p % N(k - 1);
      I(k) = 1 + pk / N(k);
    }
    //std::cout << I(0) << " " << I(1) << std::endl;

    for(unsigned k = 0; k < dim ; k++) {
      std::srand(std::time(0));
      double r = 2 * ((double) rand() / (RAND_MAX)) - 1;
      x(k, p) = a + h / 2 + (I(k) - 1) * h; // + 0.1 * r;
    }
  }

  xL.resize(dim, pow(n1, dim));
  Eigen::MatrixXd ID;
  ID.resize(dim, pow(n1, dim));
  ID.fill(1.);
  xL = (2. / (b - a)) * x - ((b + a) / (b - a)) * ID;


}

// N-point gauss quadrature points and weights by finding the roots of Legendre polynomaial
void GetGaussPointsWeights(unsigned & N, Eigen::VectorXd & xg, Eigen::VectorXd & wg) {
  unsigned N1 = N ;
  unsigned N2 = N + 1;
  Eigen::VectorXd xu;
  xu.setLinSpaced(N1, -1, 1);
  xg.resize(N1);
  for(unsigned i = 0; i <= N - 1; i++) {
    xg(i) = cos((2 * i + 1) * M_PI / (2 * (N - 1) + 2)) + (0.27 / N1) * sin(M_PI * xu(i) * (N - 1) / N2) ;
  }
  Eigen::MatrixXd L(N1, N2);
  L.fill(0.);
  Eigen::VectorXd Lp(N1);
  Lp.fill(0.);
  Eigen::VectorXd y0(xg.size());
  y0.fill(2);
  double eps = 1e-15;
  Eigen::VectorXd d = xg - y0;

  double max = d.cwiseAbs().maxCoeff();

  //Newton step for finding the roots
  while(max > eps) {

    L.col(0).fill(1.);
    L.col(1) = xg;

    for(unsigned k = 2; k < N2; k++) {
      for(unsigned i = 0; i < N1; i++) {
        L(i, k) = ((2 * k - 1) * xg(i) * L(i, k - 1) - (k - 1) * L(i, k - 2)) / k;
      }
    }


    for(unsigned i = 0; i < N1; i++) {
      Lp(i) = N2 * (L(i, N1 - 1) - xg(i) * L(i, N2 - 1)) / (1 - xg(i) * xg(i));
    }

    y0 = xg;

    for(unsigned i = 0; i < N1; i++) {
      xg(i) =  y0(i) - L(i, N2 - 1) / Lp(i);
    }

    d = xg - y0;

    max = d.cwiseAbs().maxCoeff();
  }

  // compute the weights from the roots
  wg.resize(N1);
  for(unsigned i = 0; i < N1; i++) {
    double r = double(N2) / double(N1);
    wg(i) = (2) / ((1 - xg(i) * xg(i)) * Lp(i) * Lp(i)) * r * r;
  }
}


double GetDistance(const Eigen::VectorXd & x) {

  double radius = 0.75;
  Eigen::VectorXd xc(x.size());
  xc.fill(0.);

  double rx = 0;
  for(unsigned i = 0; i < x.size(); i++) {
    rx += (x[i] - xc[i]) * (x[i] - xc[i]);
  }
  return radius - sqrt(rx);

}


double get_g(const double & r, const double & T, const unsigned & n) {
  double rn = pow(r, n);
  return (-1. + r) * (-1 + r * rn + T - r * T) / (1. + (-1. + n * (-1 + r)) * rn);
}

double get_r(const double & T, const unsigned & n) {
  double r0 = 2.;
  double r = 0;
  while(fabs(r - r0) > 1.0e-10) {
    r0 = r;
    r = r0 - get_g(r0, T, n);
  }
  return r;
}

void PrintMat(std::vector< std::vector<double> >& M) {

  for(unsigned i = 0; i < M.size(); i++) {
    for(unsigned j = 0; j < M[i].size(); j++) {
      std::cout << M[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "\n" << std::endl;
}


void PrintVec(std::vector<double>& v) {
  for(unsigned i = 0; i < v.size(); i++) {

    std::cout << v[i] << " ";
  }
  std::cout << "\n" << std::endl;
}


