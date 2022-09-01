
#include "MyEigenFunctions.hpp"
#include <numeric>
#include <iostream>

namespace femus {

  void FindBestFit(const std::vector < std::vector < double > > &xp, boost::optional < const std::vector < double > & > w, const std::vector < double > &N, std::vector < double > &a, double &d) {

    const unsigned& dim = N.size();
    a.resize(dim);

    unsigned np = xp.size();
    Eigen::MatrixXd m(np, dim);

    std::vector < double > xg(dim, 0.);


    if(w) {
      //Calculate centroid
      double wSum = 0;
      for(unsigned i = 0; i < np; i++) {
        wSum += (*w)[i];
        for(unsigned j = 0; j < dim; j++) {
          xg[j] += (*w)[i] * xp[i][j];
        }
      }
      for(unsigned j = 0; j < dim; j++) {
        xg[j] /= wSum;
      }

      //Fill matrix to be passed to JacobiSVD
      for(unsigned i = 0; i < np; i++) {
        for(unsigned j = 0; j < dim; j++) {
          m(i, j) = sqrt((*w)[i]) * (xp[i][j] - xg[j]);
        }
      }
    }
    else {
      //Calculate centroid
      for(unsigned i = 0; i < np; i++) {
        for(unsigned j = 0; j < dim; j++) {
          xg[j] += xp[i][j];
        }
      }
      for(unsigned j = 0; j < dim; j++) {
        xg[j] /= np;
      }
      //Fill matrix to be passed to JacobiSVD
      for(unsigned i = 0; i < np; i++) {
        for(unsigned j = 0; j < dim; j++) {
          m(i, j) = (xp[i][j] - xg[j]);
        }
      }
    }


    Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Eigen::MatrixXd &v = svd.matrixV();

    // use singular vector associated with min singular vector
    double aDotN = 0.;
    for(unsigned i = 0; i < dim; i++) {
      a[i] = v(i, dim - 1);
      aDotN += a[i] * N[i];
    }

    //Rotate normal by pi if Normal dot coefficents is less than zero
    if(aDotN < 0) {
      for(unsigned i = 0; i < dim; i++) {
        a[i] *= -1.;
      }
    }
  
    //Calculate constant d in ax+by+d=0 or ax+by+cz+d=0
    d = 0.;
    for(unsigned i = 0; i < dim; i++) {
      d -= a[i] * xg[i];
    }

// // use singular vector associated with min singular vector
//     for(unsigned i = 0; i < dim; i++) {
//       a[i] = v(i, dim - 1);
//     }
// 
// 
//     double aDotN = std::inner_product(a.begin(), a.end(), N.begin(), 0);
// 
//     //Rotate normal by pi if Normal dot coefficients is less than zero
//     if(aDotN < 0) {
//       for(unsigned i = 0; i < dim; i++) {
//         a[i] = -a[i];
//       }
//     }
// 
//     //Calculate constant d in ax+by+d=0 or ax+by+cz+d=0
//     d = -std::inner_product(a.begin(), a.end(), xg.begin(), 0);
  }
  
  void FindQuadraticBestFit(const std::vector < std::vector < double > > &xp, boost::optional < const std::vector < double > & > w, const std::vector < double > &N, std::vector < double > &a) {
    const unsigned& dim = N.size();
    const unsigned nParam = 4 * dim - 2;
    std::vector < double > aT;
    aT.resize(nParam);
    a.resize(nParam);
    unsigned np = xp.size();
    Eigen::MatrixXd m(np, nParam);

    std::vector < double > xg(dim, 0.);
    std::vector < double > dx(np, 0.);
    std::vector < double > dy(np, 0.);
    
    double maxDX = 0.;
    double maxDY = 0.;
    
    if(w) {
      //Calculate centroid
      double wSum = 0;
      for(unsigned i = 0; i < np; i++) {
        wSum += (*w)[i];
        for(unsigned j = 0; j < dim; j++) {
          xg[j] += (*w)[i] * xp[i][j];
        }
      }
      for(unsigned j = 0; j < dim; j++) {
        xg[j] /= wSum;
      }
      
      for(unsigned i = 0; i < np; i++){
        dx[i] = (xp[i][0] - xg[0]);
        dy[i] = (xp[i][1] - xg[1]);
      }
      maxDX = *max_element(dx.begin(), dx.end());
      maxDY = *max_element(dy.begin(), dy.end());
      
      //Fill matrix to be passed to JacobiSVD
      for(unsigned i = 0; i < np; i++) {
        unsigned cnt = 0;  
        for(int o = 2; o >= 0; o--){
          for(int b = o; b >= 0; b--){
            m(i, cnt) = sqrt((*w)[i]) * pow(dx[i] / maxDX, b) * pow(dy[i] / maxDY, o - b);
            cnt ++;
          }
        }
        if(cnt != nParam) {std::cerr<<"3D best fit not yet implemented!"; abort();}
      }
    }
    else {
      //Calculate centroid
      for(unsigned i = 0; i < np; i++) {
        for(unsigned j = 0; j < dim; j++) {
          xg[j] += xp[i][j];
        }
      }
      for(unsigned j = 0; j < dim; j++) {
        xg[j] /= np;
      }
      
      for(unsigned i = 0; i < np; i++){
        dx[i] = (xp[i][0] - xg[0]);
        dy[i] = (xp[i][1] - xg[1]);
      }
      maxDX = *max_element(dx.begin(), dx.end());
      maxDY = *max_element(dy.begin(), dy.end());
      
      //Fill matrix to be passed to JacobiSVD for Ax2 + Bxy + Cy2+ Dx + Ey + F = 0
      for(unsigned i = 0; i < np; i++) {
        unsigned cnt = 0;  
        for(int o = 2; o >= 0; o--){
          for(int b = o; b >= 0; b--){
            m(i, cnt) = pow(dx[i] / maxDX, b) * pow(dy[i] / maxDY, o - b);
            cnt ++;
          }  
        }
        if(cnt != nParam) {std::cerr<<"3D best fit not yet implemented!"; abort();}  
      }
    }
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Eigen::MatrixXd &v = svd.matrixV();
    
    // use singular vector associated with min singular vector
    double aDotN = 0.;
    for(unsigned i = 0; i < nParam; i++) {
      aT[i] = v(i, nParam - 1);
//       aDotN += aT[i] * N[i]; //TODO 
    }
    
    a[0] = aT[0] / (maxDX * maxDX);
    a[1] = aT[1] / (maxDX * maxDY);
    a[2] = aT[2] / (maxDY * maxDY);
    a[3] = ( aT[3] / maxDX ) - ( 2 * aT[0] * xg[0] / (maxDX * maxDX) ) - ( aT[1] * xg[1] / (maxDX * maxDY) );
    a[4] = ( aT[4] / maxDY ) - ( 2 * aT[2] * xg[1] / (maxDY * maxDY) ) - ( aT[1] * xg[0] / (maxDX * maxDY) );
    a[5] = ( aT[0] * xg[0] * xg[0] / (maxDX * maxDX) ) + ( aT[1] * xg[0] * xg[1] / (maxDX * maxDY) ) + ( aT[2] * xg[1] * xg[1] / (maxDY * maxDY) )
           - ( aT[3] * xg[0] / maxDX ) - ( aT[4] * xg[1] / maxDY ) + aT[5]; 
    
    
//     //Rotate normal by pi if Normal dot coefficents is less than zero   TODO for quadric
//     if(aDotN < 0) {
//       for(unsigned i = 0; i < nParam; i++) {
//         aT[i] *= -1.;
//       }
//     }
  }

  
}
