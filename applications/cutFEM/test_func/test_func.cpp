#include "cutFemWeightParabola.hpp"
// #include "parabolaIntegration.hpp"
using namespace std;

int main() {
       typedef cpp_bin_float_oct Type;
       int s = 0;
       int table = 0;
       Type a(0);
       Type c (1);
       PointT <Type> p1,p2,p3;
       p1 = { static_cast<Type>(0), static_cast<Type>(0.125) };
       p2 = { static_cast<Type>( 0.375), static_cast<Type>(1) };
       p3 = { static_cast<Type>((p1.x+p2.x)/2.0), static_cast<Type>(0.125) };
      std::vector<double>weightCF;

      CutFemWeightParabola <double,Type> Pweights(QUAD, 3, "legendre");
      Pweights(s, a, c, table, p1, p2, p3, weightCF);

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
