#include "cutFemWeightParabola.hpp"
// #include "parabolaIntegration.hpp"
using namespace std;

int main() {
  typedef cpp_bin_float_oct Type;
  int s = 0;
  int table = 0;
  Type a(0);
  Type c(1);
  PointT <Type> p1, p2, p3;
  p1 = { static_cast<Type>(0), static_cast<Type>(0.125) };
  p2 = { static_cast<Type>(0.375), static_cast<Type>(1) };
  p3 = { static_cast<Type>((p1.x + p2.x) / 2.0), static_cast<Type>(0.125) };
  std::vector<double>weightCF;

  CutFemWeightParabola <double, Type> Pweights(QUAD, 4, "legendre");
  Pweights(s, a, c, table, p1, p2, p3, weightCF);




  double Area0 = 0;
  double Area = 0;
  double Ix = 0;
  double Iy = 0;

  double Ix3 = 0;
  double Ix2y = 0;
  double Ixy2 = 0;
  double Iy3 = 0;

  double Ix2y2 = 0;

  const double* gaussWeight =  Pweights.GetGaussWeightPointer();
  const double* xg = Pweights.GetGaussCoordinatePointer(0);
  const double* yg = Pweights.GetGaussCoordinatePointer(1);
  for(unsigned ig = 0; ig < weightCF.size(); ig++) {
    Area0 += gaussWeight[ig];
    Area += weightCF[ig] * gaussWeight[ig];
    Ix += xg[ig] * weightCF[ig] * gaussWeight[ig];
    Iy += yg[ig] * weightCF[ig] * gaussWeight[ig];


    Ix3 += xg[ig] * xg[ig] * xg[ig] * weightCF[ig] * gaussWeight[ig];
    Ix2y += xg[ig] * xg[ig] * yg[ig] * weightCF[ig] * gaussWeight[ig];
    Ixy2 += xg[ig] * yg[ig] * yg[ig] * weightCF[ig] * gaussWeight[ig];
    Iy3 += yg[ig] * yg[ig] * yg[ig] * weightCF[ig] * gaussWeight[ig];

    Ix2y2 += xg[ig] * xg[ig] * yg[ig] * yg[ig] * weightCF[ig] * gaussWeight[ig];


  }

  std::cout<<"Area0 = "<<Area0<<std::endl;
  std::cout<<"Area = "<<Area<<std::endl;
  std::cout<<"Ix = "<<Ix<<std::endl;
  std::cout<<"Iy = "<<Iy<<std::endl;

  std::cout<<"Ix3 = "<<Ix3<<std::endl;
  std::cout<<"Ix2y = "<<Ix2y<<std::endl;
  std::cout<<"Ixy2 = "<<Ixy2<<std::endl;
  std::cout<<"Iy3 = "<<Iy3<<std::endl;

  std::cout<<"Ix2y2 = "<<Ix2y2<<std::endl;



  /*
  void CutFemWeightParabola<TypeIO, TypeA>::operator()(const int &s, const int &a, const int &c,  const int &table, PointT <TypeIO> &p1,  PointT <TypeI0> &p2, const PointT <TypeIO> &p3,  std::vector <TypeIO> &weightCF)*/


//     typedef cpp_bin_float_oct Type;      //     typedef double Type;
// //     std::cout.precision(16);
//
//     srand(10); // Fixed seed for random number generation
//     std::vector<std::vector<double>> pointsVector;
//     unsigned int m = 0;
//     unsigned int n = 0;
//     int s = 0;
//     int table = 0;
//     Type a(0);
//     Type c (1) ;
//     int maxDepth = 5;
//     std::vector<OctreeNode<Type>> roots;
//
//     for (int ttable = 0; ttable < 1; ++ttable) {
//         OctreeNode<Type> root({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0},ttable, 0);
//         if(ttable == 0 || ttable == 1 || ttable == 2 || ttable == 4 || ttable == 6){
//           root.subdivideWithRelativeError(maxDepth, 0.01);
//         }
//         else {
//           root.subdivideWithRelativeError(3, 0.1);
//         }
//
//         std::cout << "Octree Structure:\n";
//         roots.push_back(root);
//     }
//
//     printOctreeStructure(&roots[0]);
//     return 0;


}
