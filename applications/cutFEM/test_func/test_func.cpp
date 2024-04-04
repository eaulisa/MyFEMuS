
#include "cutFemWeightParabola.hpp"

int main() {


      CutFemWeightParabola <double,cpp_bin_float_oct> a(QUAD, 3, "legendre");


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
