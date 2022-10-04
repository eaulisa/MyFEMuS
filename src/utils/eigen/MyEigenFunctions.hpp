
#ifndef __femus_utils_myeigenfunctions_hpp__
#define __femus_utils_myeigenfunctions_hpp__

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/Eigenvalues> 
#include <vector>
#include <boost/optional.hpp>
namespace femus {
  void FindBestFit(const std::vector < std::vector < double > > &xp, boost::optional < const std::vector < double > & > w, const std::vector < double > &N, std::vector < double > &a, double &d);
  void FindQuadraticBestFit(const std::vector < std::vector < double > > &xp, boost::optional < const std::vector < double > & > w, const std::vector < double > &N, std::vector < double > &a);
  void FindParabolaBestFit(const std::vector < std::vector < double > > &xp, boost::optional < const std::vector < double > & > w, const std::vector < double > &N, std::vector < double > &a);
  void GetQuadricBestFit(const std::vector < std::vector < double > > &x, boost::optional < const std::vector < double > & > w, std::vector < double > &N, std::vector < double > &a);
}
#endif
