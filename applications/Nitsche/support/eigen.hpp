
#ifndef __femus_eigen_hpp__
#define __femus_eigen_hpp__

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include </usr/include/eigen3/Eigen/src/Core/util/DisableStupidWarnings.h>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

using namespace femus;

void Testing(double & a, double & b, const unsigned & m, const unsigned & dim, Eigen::MatrixXd & x,
             Eigen::VectorXd & w_new, std::vector<double> &dist, const double & eps, double & QuadSum, double & IntSum); 

void Cheb(const unsigned & m, Eigen::VectorXd & xg, Eigen::MatrixXd & C);

void AssembleMatEigen(std::vector<double>& VxL, std::vector<double> &VxR, const unsigned & m, const unsigned & dim, const unsigned & np, Eigen::Tensor<double, 3, Eigen::RowMajor>  &PmX, Eigen::MatrixXd & Pg,  Eigen::VectorXd & wg, Eigen::MatrixXd & A, Eigen::VectorXd & F);

void GetChebGaussF(const unsigned & dim, const unsigned & m, std::vector<double> &VxL, std::vector<double> &VxU, Eigen::MatrixXd & Pg,  Eigen::VectorXd & wg, Eigen::VectorXd & F);

void GetChebGaussF(const unsigned & dim, const unsigned & m, const std::vector<double> &jac, Eigen::MatrixXd & Pg,  Eigen::VectorXd & wg, Eigen::VectorXd & F);

void GetChebXInfo(const unsigned & m, const unsigned & dim, const unsigned & np, Eigen::MatrixXd & xL, Eigen::Tensor<double, 3, Eigen::RowMajor>& PmX);

void GetMultiDimChebMatrix(const unsigned & dim, const unsigned & m, const unsigned & np, Eigen::Tensor<double, 3, Eigen::RowMajor>  &PmX, Eigen::MatrixXd & A);

void SolWeightEigen(Eigen::MatrixXd & A, Eigen::VectorXd & F, Eigen::VectorXd & wP, Eigen::VectorXd & w_new);

void PrintMarkers(const unsigned & dim, const Eigen::MatrixXd & xP, const std::vector <double> &dist,
                  const Eigen::VectorXd wP, const Eigen::VectorXd & w_new, const unsigned & l, const unsigned & t);

void  GetParticlesOnBox(const double & a, const double & b, const unsigned & n1, const unsigned & dim, Eigen::MatrixXd & x, Eigen::MatrixXd & xL);

void GetGaussPointsWeights(unsigned & N, Eigen::VectorXd & xg, Eigen::VectorXd & wg);

double GetDistance(const Eigen::VectorXd & x);

double get_g(const double & r, const double & T, const unsigned & n);
double get_r(const double & T, const unsigned & n);
void PrintMat(std::vector< std::vector<double> >& M);

void PrintVec(std::vector<double>& v);
void PrintMatlabMatrix(Eigen::MatrixXd &A);

#endif

