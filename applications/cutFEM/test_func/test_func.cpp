#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <algorithm>    // std::sort
#include <ctime>
#include <cstdlib>
#include <climits>
#include <typeinfo>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/filesystem.hpp>
#include <map>
#include <fstream>
#include <string>
// #include <filesystem>

namespace fs = boost::filesystem;

using namespace std;

#include "Rebuild.hpp"
#include "parabolaIntegration.hpp"
#include "PolynomialBases.hpp"

int main() {
  typedef cpp_bin_float_oct Type;
  int s = 0;
  int table = 0;
  Type a(0);
  Type c(1);
  PointT <Type> p1, p2, p3;
  p1 = { static_cast<Type>(0), static_cast<Type>(0.5) };
  p2 = { static_cast<Type>(0.5), static_cast<Type>(1) };
  p3 = { static_cast<Type>((p1.x + p2.x) / 2.0), static_cast<Type>(0.125) };
  double Area0 = 0, Area = 0, Ix = 0, Iy = 0, Ix3 = 0, Ix2y = 0, Ixy2 = 0, Iy3 = 0, Ix2y2 = 0;


  std::vector<double>weightCF;
  std::vector< double > interp_point_weights;
  CutFemWeightParabola <double, Type> Pweights(QUAD, 3, "legendre");
  Pweights(s, a, c, table, p1, p2, p3, weightCF);

  cout << " weightCF = " ;
  for(unsigned ig = 0; ig < weightCF.size(); ig++) {
   cout <<  weightCF[ig] << " " ;
  }
  cout << endl;


  const double* gaussWeight =  Pweights.GetGaussWeightPointer();
  const double* xg = Pweights.GetGaussCoordinatePointer(0);
  const double* yg = Pweights.GetGaussCoordinatePointer(1);





  int maxDepth = 5;
  int degree = 3;
  double percent = 0.001;
//   std::vector<OctreeNode<Type>> roots;
  std::vector<OctreeNode<Type>>loadedRoots;

  generateAndLoadOctrees<Type>(maxDepth, degree, percent, Pweights, /*roots,*/ loadedRoots);
// Example: Search for a point in the loaded Octree
  table = 0 ;
  Point3D searchPoint(0.5, 0.5, 0.125);
  OctreeNode<Type>* result = loadedRoots[table].search(searchPoint);
  if (result) {
      std::cout << "Found the smallest sub-cube containing the search point." << std::endl;
      std::cout << "\nSearch Point: (" << searchPoint.x << ", " << searchPoint.y << ", " << searchPoint.z << ")\n";
      std::cout << "Smallest Sub-cube Bounds: ";
      std::cout << "(" << result->minBounds.x << ", " << result->minBounds.y << ", " << result->minBounds.z << ") to ";
      std::cout << "(" << result->maxBounds.x << ", " << result->maxBounds.y << ", " << result->maxBounds.z << ")\n";
      std::cout << "depth : = " << result->depth << " \n";

      std::vector<double>interp_point = {searchPoint.x, searchPoint.y, searchPoint.z};
      std::cout << "\n interp Point: (" << interp_point[0] << ", " << interp_point[1] << ", " << interp_point[2] << ")\n";

      trilinier_interpolation_vector(result->corners, result->cornerWeights, interp_point, interp_point_weights);

      Area = GaussIntegral(0,0,xg,yg,interp_point_weights, gaussWeight);
      Ix  = GaussIntegral(1,0,xg,yg,interp_point_weights, gaussWeight);
      Iy  = GaussIntegral(0,1,xg,yg,interp_point_weights, gaussWeight);
      Ix3  = GaussIntegral(3,0,xg,yg,interp_point_weights, gaussWeight);
      Ix2y  = GaussIntegral(2,1,xg,yg,interp_point_weights, gaussWeight);
      Ixy2  = GaussIntegral(1,2,xg,yg,interp_point_weights, gaussWeight);
      Iy3 = GaussIntegral(0,3,xg,yg,interp_point_weights, gaussWeight);
      Ix2y2  = GaussIntegral(2,2,xg,yg,interp_point_weights, gaussWeight);

      std::cout << "Area0 = " << Area0 << std::endl;
      std::cout << "Area = " << Area << std::endl;
      std::cout << "Ix = " << Ix << std::endl;
      std::cout << "Iy = " << Iy << std::endl;
      std::cout << "Ix3 = " << Ix3 << std::endl;
      std::cout << "Ix2y = " << Ix2y << std::endl;
      std::cout << "Ixy2 = " << Ixy2 << std::endl;
      std::cout << "Iy3 = " << Iy3 << std::endl;
      std::cout << "Ix2y2 = " << Ix2y2 << std::endl;
  }
  else {
      std::cout << "Search point not found in the Octree." << std::endl;
  }

  //checking other side:



//   printOctreeStructure(&loadedRoot);








  unsigned nInt;
  std::vector<std::vector<double>> unitxv = {{0., 1., 1., 0.}, {0., 0., 1., 1.}};
//   std::vector<std::vector<double>> xv = {{1., 2., 2, 1.}, {1., 1., 2., 2.}};  // it was written like {{x1,x2,x3,x4},{y1,y2,y3,y4}}
      //   std::vector<std::vector<double>> xv = {{2., 2., 1., 1.}, {1., 2., 2., 1.}};
    //   std::vector<std::vector<double>> xv = {{2., 1., 1., 2.}, {2., 2., 1., 1.}};
    //   std::vector<std::vector<double>> xv = {{1., 1., 2., 2.}, {2., 1., 1., 2.}};
        std::vector<std::vector<double>> xv = {{0., 1., 1., 0.}, {0., 0., 1., 1.}};
  for ( unsigned i = 0; i < 4 ; i ++){
//        std::vector<double> A = {-10, 0, 0, 4, 1, -0.5}; // {a,b,c,d,e,f} means ax^2 + bxy + cy^2 + dx + ey + f = 0
       std::vector<double> A = {0, 0, -10, 1, 4, -0.5};
//       std::vector<double> A = {0, -1.8, 0, 1, 5, -4.5}; // horizotal prabola

      unsigned nPoints = 3;
      unsigned dim = 2;
      short unsigned ielType = 3; //quad
      unsigned femType = 0; //linear FEM


      cout << "..........."<<endl;
      cout <<"\n \n Box(" <<i<<"): xv = {"<<xv[0][0]<<" "<<xv[0][1]<<" "<<xv[0][2]<<" "<<xv[0][3]<<"},{"<<xv[1][0]<<" "<<xv[1][1]<<" "<<xv[1][2]<<" "<<xv[1][3]<<"}"  <<endl;

      std::pair<std::vector<std::vector<double>>, std::vector<double>> xp = GetCellPointsFromQuadric(xv, A, nPoints, nInt);     //This fins the points in physical space

      std::vector < std::vector < std::vector <double > > > aP(1);
      ProjectNodalToPolynomialCoefficients(aP[femType], xv, ielType, femType);   //TODO what does this do?


      std::vector<int> xvsign(4);
      std::vector<int> unitxvsign(4);

      std::vector<std::vector<double>> xi(nPoints, std::vector<double>(2, 0.));
      for(unsigned i = 0; i < nPoints; i++) {
        bool inverseMapping = GetInverseMapping(femType, ielType, aP, xp.first[i], xi[i], 100);        //This maps the phsical points to {(-1,-1),(1,1)} box
//         std::cout << " \nx[i] physical value " << i << " " << xp.first[i][0] << " " << xp.first[i][1] << std::endl;
//         std::cout << " x[i] value in (-1,1) " << i << " " << xi[i][0] << " " << xi[i][1] << std::endl;
        xi[i] = {0.5 * (xi[i][0] + 1.), 0.5 * (xi[i][1] + 1.)};                                        // //This maps the points to unit box
        std::cout <<"value in unit box"<< i << " " << xi[i][0] << " " << xi[i][1] << std::endl;
      }


      //finiding monotone
      if ((xi[0][0] < xi[1][0] && xi[1][0] < xi[2][0]) || (xi[0][0] > xi[1][0] && xi[1][0] > xi[2][0])) cout << " monotonous in x " << endl;
      else cout << " non-monotonous in x " << endl;

      if ((xi[0][1] < xi[1][1] && xi[1][1] < xi[2][1]) || (xi[0][1] > xi[1][1] && xi[1][1]> xi[2][1])) cout << " monotonous in y " << endl;
      else cout << " non-monotonous in y " << endl;

      cout << " xv sign = {" ;
      for (unsigned l = 0; l < 4 ; l++){
        xvsign[l] = ((A[0]*xv[0][l]*xv[0][l] + A[1]*xv[0][l]*xv[1][l] + A[2]*xv[1][l]*xv[1][l] + A[3]*xv[0][l] + A[4]*xv[1][l] + A[5]) >= 0)? 1:-1 ;
        cout << xvsign[l]<<", ";
      }
      cout << "} "<<endl;


      if ((xi[0][0] < xi[1][0] && xi[1][0] < xi[2][0]) || (xi[0][0] > xi[1][0] && xi[1][0] > xi[2][0])){  //vertical

          p1 = { static_cast<Type>(xi[0][0]), static_cast<Type>(xi[0][1]) };
          p2 = { static_cast<Type>(xi[2][0]), static_cast<Type>(xi[2][1]) };
          p3 = { static_cast<Type>(xi[1][0]), static_cast<Type>(xi[1][1]) };


          Parabola <Type> parabola = get_parabola_equation(p1,p2,p3);
          int normal;

          cout << " unit box sign = {" ;
          for (unsigned l = 0; l < 4 ; l++){
            unitxvsign[l] = ((static_cast<double>(parabola.k)*unitxv[0][l]*unitxv[0][l] + static_cast<double>(parabola.b)*unitxv[0][l] + static_cast<double>(parabola.d) + unitxv[1][l])>0)? 1:-1;
            cout << unitxvsign[l]<<", ";
          }
          cout << "} "<<endl;

          cout <<  "( "<<p1.x<<","<<p1.y<<" )"<< " , ( "<<p2.x<<","<<p2.y<<" )"<< " , ( "<<p3.x<<","<<p3.y<<" ) "<< endl;
          cout<< parabola.k << "x^2+ " << parabola.b << "x+ " << parabola.d << "+y =0 "<< endl;

          normal = checkVectorRelation(xvsign,unitxvsign);
          cout << " normal = " << normal <<endl;

          int intersect_number;
          unsigned table_number;
          std::vector <Type> intersection;
          std::vector <Type> interp_point;
          CheckIntersection<Type>(intersect_number, table_number, intersection, interp_point, parabola);
          cout << " table number = " << table_number <<endl;
          p3.x = (p1.x+p2.x)/2;
          p3.y = -parabola.k * p3.x * p3.x - parabola.b * p3.x - parabola.d ;
          cout <<  "( "<<p1.x<<","<<p1.y<<" )"<< " , ( "<<p2.x<<","<<p2.y<<" )"<< " , ( "<<p3.x<<","<<p3.y<<" ) "<< endl;

          if (interp_point.size() == 2){
            Point3D searchP(static_cast<double>(interp_point[0]),static_cast<double>(interp_point[1]), static_cast<double>(p3.y));
            OctreeNode<Type>* result = loadedRoots[table].search(searchP);
              if (result) {
                  std::cout << "Found the smallest sub-cube containing the search point." << std::endl;
                  std::cout << "\nSearch Point: (" << searchP.x << ", " << searchP.y << ", " << searchP.z << ")\n";
                  std::cout << "Smallest Sub-cube Bounds: ";
                  std::cout << "(" << result->minBounds.x << ", " << result->minBounds.y << ", " << result->minBounds.z << ") to ";
                  std::cout << "(" << result->maxBounds.x << ", " << result->maxBounds.y << ", " << result->maxBounds.z << ")\n";
                  std::cout << "depth : = " << result->depth << " \n";

                  std::vector<double>interp_point = {searchP.x, searchP.y, searchP.z};
                  std::cout << "\n interp Point: (" << interp_point[0] << ", " << interp_point[1] << ", " << interp_point[2] << ")\n";

                  trilinier_interpolation_vector(result->corners, result->cornerWeights, interp_point, interp_point_weights);

                  if (normal == -1){
                    for (unsigned aq=0;aq<interp_point_weights.size();aq++){
                        interp_point_weights[aq] = -1*interp_point_weights[aq];
                    }
                  }

                  Area = GaussIntegral(0,0,xg,yg,interp_point_weights, gaussWeight);
                  Ix  = GaussIntegral(1,0,xg,yg,interp_point_weights, gaussWeight);
                  Iy  = GaussIntegral(0,1,xg,yg,interp_point_weights, gaussWeight);
                  Ix3  = GaussIntegral(3,0,xg,yg,interp_point_weights, gaussWeight);
                  Ix2y  = GaussIntegral(2,1,xg,yg,interp_point_weights, gaussWeight);
                  Ixy2  = GaussIntegral(1,2,xg,yg,interp_point_weights, gaussWeight);
                  Iy3 = GaussIntegral(0,3,xg,yg,interp_point_weights, gaussWeight);
                  Ix2y2  = GaussIntegral(2,2,xg,yg,interp_point_weights, gaussWeight);

                  std::cout << "Area0 = " << Area0 << std::endl;
                  std::cout << "Area = " << Area << std::endl;
                  std::cout << "Ix = " << Ix << std::endl;
                  std::cout << "Iy = " << Iy << std::endl;
                  std::cout << "Ix3 = " << Ix3 << std::endl;
                  std::cout << "Ix2y = " << Ix2y << std::endl;
                  std::cout << "Ixy2 = " << Ixy2 << std::endl;
                  std::cout << "Iy3 = " << Iy3 << std::endl;
                  std::cout << "Ix2y2 = " << Ix2y2 << std::endl;
              }
              else {
                  std::cout << "Search point not found in the Octree." << std::endl;
              }



          }


      }

      else if ((xi[0][1] < xi[1][1] && xi[1][1] < xi[2][1]) || (xi[0][1] > xi[1][1] && xi[1][1]> xi[2][1])){  //horizontal

          p1 = { static_cast<Type>(xi[0][1]), static_cast<Type>(xi[0][0]) };
          p2 = { static_cast<Type>(xi[2][1]), static_cast<Type>(xi[2][0]) };
          p3 = { static_cast<Type>(xi[1][1]), static_cast<Type>(xi[1][0]) };



          Parabola <Type> parabola = get_parabola_equation(p1,p2,p3);
          int normal;


          //use horizotal parabola for the normal
          cout << " unit box sign = {" ;
          for (unsigned l = 0; l < 4 ; l++){
            unitxvsign[l] = ((static_cast<double>(parabola.k)*unitxv[1][l]*unitxv[1][l] + static_cast<double>(parabola.b)*unitxv[1][l] + static_cast<double>(parabola.d) + unitxv[0][l])>0)? 1:-1;
            cout << unitxvsign[l]<<", ";
          }
          cout << "} "<<endl;

          cout <<  "( "<<p1.y<<","<<p1.x<<" )"<< " , ( "<<p2.y<<","<<p2.x<<" )"<< " , ( "<<p3.y<<","<<p3.x<<" ) "<< endl;
          cout<< parabola.k << "y^2+ " << parabola.b << "y+ " << parabola.d << "+x =0 "<< endl;

          normal = checkVectorRelation(xvsign,unitxvsign);
          cout << " normal = " << normal<<endl ;


          int intersect_number;
          unsigned table_number;
          std::vector <Type> intersection;
          std::vector <Type> interp_point;
          CheckIntersection<Type>(intersect_number, table_number, intersection, interp_point, parabola);
          cout << " table number = " << table_number <<endl;

          p3.x = (p1.x+p2.x)/2;
          p3.y = -parabola.k * p3.x * p3.x - parabola.b * p3.x - parabola.d ;
          cout <<  "( "<<p1.x<<","<<p1.y<<" )"<< " , ( "<<p2.x<<","<<p2.y<<" )"<< " , ( "<<p3.x<<","<<p3.y<<" ) "<< endl;

          if (interp_point.size() == 2){
            Point3D searchP(static_cast<double>(interp_point[0]),static_cast<double>(interp_point[1]), static_cast<double>(p3.y));
            OctreeNode<Type>* result = loadedRoots[table].search(searchP);
              if (result) {
                  std::cout << "Found the smallest sub-cube containing the search point." << std::endl;
                  std::cout << "\nSearch Point: (" << searchP.x << ", " << searchP.y << ", " << searchP.z << ")\n";
                  std::cout << "Smallest Sub-cube Bounds: ";
                  std::cout << "(" << result->minBounds.x << ", " << result->minBounds.y << ", " << result->minBounds.z << ") to ";
                  std::cout << "(" << result->maxBounds.x << ", " << result->maxBounds.y << ", " << result->maxBounds.z << ")\n";
                  std::cout << "depth : = " << result->depth << " \n";

                  std::vector<double>interp_point = {searchP.x, searchP.y, searchP.z};
                  std::cout << "\n interp Point: (" << interp_point[0] << ", " << interp_point[1] << ", " << interp_point[2] << ")\n";

                  trilinier_interpolation_vector(result->corners, result->cornerWeights, interp_point, interp_point_weights);

                  if (normal == 1){
                    for (unsigned aq=0;aq<interp_point_weights.size();aq++){
                        interp_point_weights[aq] = 1-interp_point_weights[aq];
                    }
                  }
                  //interchange x and y power
                  Area = GaussIntegral(0,0,xg,yg,interp_point_weights, gaussWeight);
                  Ix  = GaussIntegral(0,1,xg,yg,interp_point_weights, gaussWeight);
                  Iy  = GaussIntegral(1,0,xg,yg,interp_point_weights, gaussWeight);
                  Ix3  = GaussIntegral(0,3,xg,yg,interp_point_weights, gaussWeight);
                  Ix2y  = GaussIntegral(1,2,xg,yg,interp_point_weights, gaussWeight);
                  Ixy2  = GaussIntegral(2,1,xg,yg,interp_point_weights, gaussWeight);
                  Iy3 = GaussIntegral(3,0,xg,yg,interp_point_weights, gaussWeight);
                  Ix2y2  = GaussIntegral(2,2,xg,yg,interp_point_weights, gaussWeight);

/*
                  Area = GaussIntegral(0,0,xg,yg,interp_point_weights, gaussWeight);
                  Ix  = GaussIntegral(0,1,xg,yg,interp_point_weights, gaussWeight);
                  Iy  = GaussIntegral(1,0,xg,yg,interp_point_weights, gaussWeight);
                  Ix3  = GaussIntegral(0,3,xg,yg,interp_point_weights, gaussWeight);
                  Ix2y  = GaussIntegral(1,2,xg,yg,interp_point_weights, gaussWeight);
                  Ixy2  = GaussIntegral(2,1,xg,yg,interp_point_weights, gaussWeight);
                  Iy3 = GaussIntegral(3,0,xg,yg,interp_point_weights, gaussWeight);
                  Ix2y2  = GaussIntegral(2,2,xg,yg,interp_point_weights, gaussWeight);*/


                  std::cout << "Area = " << Area << std::endl;
                  std::cout << "Ix = " << Ix << std::endl;
                  std::cout << "Iy = " << Iy << std::endl;
                  std::cout << "Ix3 = " << Ix3 << std::endl;
                  std::cout << "Ix2y = " << Ix2y << std::endl;
                  std::cout << "Ixy2 = " << Ixy2 << std::endl;
                  std::cout << "Iy3 = " << Iy3 << std::endl;
                  std::cout << "Ix2y2 = " << Ix2y2 << std::endl;
              }
              else {
                  std::cout << "Search point not found in the Octree." << std::endl;
              }



          }

      }






      double swap;
      //change the orientation of the xv
        swap = xv[0][0];
        xv[0][0] = xv[0][1];
        xv[0][1] = xv[0][2];
        xv[0][2] = xv[0][3];
        xv[0][3] = swap;

        swap = xv[1][0];
        xv[1][0] = xv[1][1];
        xv[1][1] = xv[1][2];
        xv[1][2] = xv[1][3];
        xv[1][3] = swap;
  }

  return 0;
}




/*
  Area0 = 0, Area = 0, Ix = 0, Iy = 0, Ix3 = 0, Ix2y = 0, Ixy2 = 0, Iy3 = 0, Ix2y2 = 0;
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
  std::cout << "Area0 = " << Area0 << std::endl;
  std::cout << "Area = " << Area << std::endl;
  std::cout << "Ix = " << Ix << std::endl;
  std::cout << "Iy = " << Iy << std::endl;
  std::cout << "Ix3 = " << Ix3 << std::endl;
  std::cout << "Ix2y = " << Ix2y << std::endl;
  std::cout << "Ixy2 = " << Ixy2 << std::endl;
  std::cout << "Iy3 = " << Iy3 << std::endl;
  std::cout << "Ix2y2 = " << Ix2y2 << std::endl;
  std::cout << std::endl;
*/









//   cout << " weight CF = " << weightCF.size() << endl;
//   for(unsigned ig = 0; ig < weightCF.size(); ig++) {
//     cout  << weightCF[ig] << " ";
//   }
//   cout << endl;





//     cout << " gaussWeight = " << *gaussWeight <<endl;
//     cout<< "gaussWeight = " << weightCF.size() << endl;
//     for(unsigned ig = 0; ig < weightCF.size(); ig++) {
//     cout  << gaussWeight[ig] << " ";
//     }
//     cout<<endl;






//   std::cout << "Interp weights:\n";
//   for(size_t k = 0; k < interp_point_weights.size(); ++k) {
//     const auto& entry = interp_point_weights[k];
//     std::cout << " " << entry;
//   }
//
//   std::cout << "Interp weights error:\n";
//   for(size_t k = 0; k < interp_point_weights.size(); ++k) {
//     const auto& entry = (interp_point_weights[k] - weightCF[k]) / weightCF[k];
//     std::cout << " " << entry;
//   }




/*
    std::cout <<" points on given lebel set";
      for (const auto& innerVec : xp.first) {
        for (const auto& val : innerVec) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
      }

      // Iterate over the second vector
      std::cout << "and:" << std::endl;
      for (const auto& val : xp.second) {
          std::cout << val << " ";
      }
      std::cout << std::endl;*/
