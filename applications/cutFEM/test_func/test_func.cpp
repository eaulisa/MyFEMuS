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
#include "Fem.hpp"


template <class Type>
std::vector<double> find_Weight_CF( std::vector<OctreeNode<Type>> &loadedRoots, const std::vector<std::vector<double>> &xv, const std::vector<double> &A);




int main() {
  typedef cpp_bin_float_oct Type;
  int s = 0;
  int table = 4;
  Type a(0);
  Type c(1);
  PointT <Type> p1, p2, p3;
  p1 = { static_cast<Type>(0.4471), static_cast<Type>(1) };
  p2 = { static_cast<Type>(1), static_cast<Type>(0.4471) };
  p3 = { static_cast<Type>((p1.x + p2.x) / 2.0), static_cast<Type>(0.8291) };
  double Area0 = 0, Area = 0, Ix = 0, Iy = 0,Ixy =0, Ix3 = 0, Ix2y = 0, Ixy2 = 0, Iy3 = 0, Ix2y2 = 0;


  std::vector<double>weightCF;
  std::vector< double > interp_point_weights;
  CutFemWeightParabola <double, Type> Pweights(QUAD, 3, "legendre");
  Pweights(s, a, c, table, p1, p2, p3, weightCF);

  Fem fem = Fem(3 * 2, 2);
  unsigned quad = 3;
  unsigned linear = 0;
  const elem_type *femQuad = fem.GetFiniteElement(quad, linear);



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
//   table = 0 ;
//   Point3D searchPoint(static_cast<double>(p1.x), static_cast<double>(p2.y), static_cast<double>(p3.y));
//   OctreeNode<Type>* result = loadedRoots[table].search(searchPoint);
//   if(result) {
//     std::cout << "Found the smallest sub-cube containing the search point." << std::endl;
//     std::cout << "\nSearch Point: (" << searchPoint.x << ", " << searchPoint.y << ", " << searchPoint.z << ")\n";
//     std::cout << "Smallest Sub-cube Bounds: ";
//     std::cout << "(" << result->minBounds.x << ", " << result->minBounds.y << ", " << result->minBounds.z << ") to ";
//     std::cout << "(" << result->maxBounds.x << ", " << result->maxBounds.y << ", " << result->maxBounds.z << ")\n";
//     std::cout << "depth : = " << result->depth << " \n";
//
//     std::vector<double>interp_point = {searchPoint.x, searchPoint.y, searchPoint.z};
//     std::cout << "\n interp Point: (" << interp_point[0] << ", " << interp_point[1] << ", " << interp_point[2] << ")\n";
//
//     trilinier_interpolation_vector(result->corners, result->cornerWeights, interp_point, interp_point_weights);
//
//     Area = GaussIntegral(0, 0, xg, yg, interp_point_weights, gaussWeight);
//     Ix  = GaussIntegral(1, 0, xg, yg, interp_point_weights, gaussWeight);
//     Iy  = GaussIntegral(0, 1, xg, yg, interp_point_weights, gaussWeight);
//     Ix3  = GaussIntegral(3, 0, xg, yg, interp_point_weights, gaussWeight);
//     Ix2y  = GaussIntegral(2, 1, xg, yg, interp_point_weights, gaussWeight);
//     Ixy2  = GaussIntegral(1, 2, xg, yg, interp_point_weights, gaussWeight);
//     Iy3 = GaussIntegral(0, 3, xg, yg, interp_point_weights, gaussWeight);
//     Ix2y2  = GaussIntegral(2, 2, xg, yg, interp_point_weights, gaussWeight);
//
//     std::cout << "Area0 = " << Area0 << std::endl;
//     std::cout << "Area = " << Area << std::endl;
//     std::cout << "Ix = " << Ix << std::endl;
//     std::cout << "Iy = " << Iy << std::endl;
//     std::cout << "Ix3 = " << Ix3 << std::endl;
//     std::cout << "Ix2y = " << Ix2y << std::endl;
//     std::cout << "Ixy2 = " << Ixy2 << std::endl;
//     std::cout << "Iy3 = " << Iy3 << std::endl;
//     std::cout << "Ix2y2 = " << Ix2y2 << std::endl;
//   }
//   else {
//     std::cout << "Search point not found in the Octree." << std::endl;
//   }
//
//   Type direct_area_00 = find_area_2intersection_formula<Type>(0, 0, 0, 0, 1, 4,  p1,  p2, p3);
//   Type direct_area_10 = find_area_2intersection_formula<Type>(1, 0, 0, 0, 1, 4,  p1,  p2, p3);
//   Type direct_area_01 = find_area_2intersection_formula<Type>(0, 1, 0, 0, 1, 4,  p1,  p2, p3);
//   cout << " area = " << direct_area_00 << endl;
//   cout << " area10 = " << direct_area_10 << endl;
//   cout << " area01 = " << direct_area_01 << endl;




  //checking other side:



//   printOctreeStructure(&loadedRoot);






  unsigned nInt;
  std::vector<std::vector<double>> unitxv = {{0., 1., 1., 0.}, {0., 0., 1., 1.}};
  //   std::vector<std::vector<double>> xv = {{1., 2., 2, 1.}, {1., 1., 2., 2.}};  // it was written like {{x1,x2,x3,x4},{y1,y2,y3,y4}}
  //   std::vector<std::vector<double>> xv = {{2., 2., 1., 1.}, {1., 2., 2., 1.}};
  //   std::vector<std::vector<double>> xv = {{2., 1., 1., 2.}, {2., 2., 1., 1.}};
  //   std::vector<std::vector<double>> xv = {{1., 1., 2., 2.}, {2., 1., 1., 2.}};
//   std::vector<std::vector<double>> xv = {{0., 1., 1., 0.}, {0., 0., 1., 1.}};
//   std::vector<std::vector<double>> xv = {{0., 2.2361, 2.2361, 0.}, {0., 0., 2.2361, 2.2361}};
//     std::vector<std::vector<double>> xv = {{0., 2.2361, 2.2361, 0.}, {-2.2361, -2.2361, 0., 0.}};
//     std::vector<std::vector<double>> xv = {{-1., 1., 1., -1}, {1., 1., 3., 3.}};
//     std::vector<double> A = {1., 0., 1., 0., 0., -6.};      // Table 1 or 5 gives trouble when rotated.

    // trouble case nonlocal 1 . jg=6, jj=4
      std::vector<std::vector<double>> xv = {{-0.4, -0.4, -0.5, -0.5}, {-0.4, -0.3, -0.3, -0.4}};
//       std::vector<double> A = {-1., 0., -1., -0.82254033307585184, -0.22254033307585197, -0.14152419984551109};


      std::vector<double> A = {-1., 0., -1., -1.1774596669241479, -0.49999999999999994, -0.36910281680828139};
//       std::vector<double> A = {-1., 0., -1., -1.1, -0.42254033307585198, -0.30713508326896299};
//       -1x^{2}-1.\ y^{2}-1.0225403330758518\ x-0.22254033307585197y\ -0.23377823315309629\ =0
//       std::vector<double> A = {-1., 0., -1., -1.0225403330758518, -0.22254033307585197, -0.23377823315309629};

//     std::vector<std::vector<double>> xv = {{0.55, 2.3, 2.3, 0.55}, {0.5, 0.5, 2.25, 2.25}};
//     std::cout <<  "     llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll" <<  std::endl;
//   std::vector<std::vector<double>> xv = {{0., 2., 2., 0.}, {0., 0., 2.2361, 2.2361}};
//   std::vector<std::vector<double>> xv = {{2., 3., 3., 2.}, {0., 0., 1., 1.}};
  //     std::vector<double> A = {-10, 0, 0, 4, 1, -0.5}; // {a,b,c,d,e,f} means ax^2 + bxy + cy^2 + dx + ey + f = 0
  //     std::vector<double> A = {0, 0, -10, 1, 4, -0.5};
  // horizotal prabola
//   std::vector<double> A = {1., 0., 1., 0., 0., -6.};
    PointT <double> q1, q2, q3;
    Point3D searchP(0., 0., 0.);



  for(unsigned i = 0; i < 4 ; i ++) {
    unsigned nPoints = 3;
    unsigned dim = 2;
    short unsigned ielType = 3; //quad
    unsigned femType = 0; //linear FEM


    cout << "\n....................................................................." << endl;
    cout << "\n \n Box(" << i << "): xv = {" << xv[0][0] << " " << xv[0][1] << " " << xv[0][2] << " " << xv[0][3] << "},{" << xv[1][0] << " " << xv[1][1] << " " << xv[1][2] << " " << xv[1][3] << "}"  << endl;

    std::pair<std::vector<std::vector<double>>, std::vector<double>> xp = GetCellPointsFromQuadric(xv, A, nPoints, nInt);     //This finds the points in physical space

   std::cout << "Intersection Point = (" << xp.first[0][0] << ", " << xp.first[0][1] << "), (" <<  xp.first[1][0] << ", " << xp.first[1][1] << "), (" <<  xp.first[2][0] << ", " << xp.first[2][1] << "), ("<<   std::endl;


//     for (size_t i = 0; i < xp.first.size(); ++i) {
//         std::cout << "(" << xp.first[i][0] << ", " << xp.first[i][1] << ")\n";
//     }

    vector<vector<double>> qvector = transformPoints(xv, unitxv, xp.first);
    cout << "affine transformation" << i+1;
    for (size_t i = 0; i < qvector.size(); ++i) {
        cout << ": (" << qvector[i][0] << ", " << qvector[i][1] << "), ";
    }

    cout <<  endl;


    std::vector < std::vector < std::vector <double > > > aP(1);
    ProjectNodalToPolynomialCoefficients(aP[femType], xv, ielType, femType);   //TODO what does this do?

//     cout << "======================================== >>>>>>>>>>>>>>> size of aP = "<<aP.size() << endl;
//         for (const auto& matrix : aP) {
//         std::cout << "{" << std::endl;
//         for (const auto& row : matrix) {
//             std::cout << "  {";
//             for (const auto& elem : row) {
//                 std::cout << elem << " ";
//               }
//             std::cout << "}" << std::endl;
//             }
//           std::cout << "}" << std::endl;
//         }


    std::vector<int> xvsign(4);
    std::vector<int> unitxvsign(4);
    std::vector<std::vector<double>> xi(nPoints, std::vector<double>(2, 0.));

    for(unsigned i = 0; i < nPoints; i++) {
      bool inverseMapping = GetInverseMapping(femType, ielType, aP, xp.first[i], xi[i], 100);        //This maps the phsical points to {(-1,-1),(1,1)} box
//         std::cout << " \nx[i] physical value " << i << " " << xp.first[i][0] << " " << xp.first[i][1] << std::endl;
//         std::cout << " x[i] value in (-1,1) " << i << " " << xi[i][0] << " " << xi[i][1] << std::endl;
      xi[i] = {0.5 * (xi[i][0] + 1.), 0.5 * (xi[i][1] + 1.)};                                        // //This maps the points to unit box
      std::cout << "value in unit box" << i << " " << xi[i][0] << " " << xi[i][1] << std::endl;
    }






    //finiding monotone
    if((xi[0][0] < xi[1][0] && xi[1][0] < xi[2][0]) || (xi[0][0] > xi[1][0] && xi[1][0] > xi[2][0])) cout << " monotonous in x " << endl;
    else cout << " non-monotonous in x " << endl;

    if((xi[0][1] < xi[1][1] && xi[1][1] < xi[2][1]) || (xi[0][1] > xi[1][1] && xi[1][1] > xi[2][1])) cout << " monotonous in y " << endl;
    else cout << " non-monotonous in y " << endl;

    cout << " xv sign = {" ;
    for(unsigned l = 0; l < 4 ; l++) {
      xvsign[l] = ((A[0] * xv[0][l] * xv[0][l] + A[1] * xv[0][l] * xv[1][l] + A[2] * xv[1][l] * xv[1][l] + A[3] * xv[0][l] + A[4] * xv[1][l] + A[5]) >= 0) ? 1 : -1 ;
      cout << xvsign[l] << ", ";
    }
    cout << "} " << endl;


    if((xi[0][0] < xi[1][0] && xi[1][0] < xi[2][0]) || (xi[0][0] > xi[1][0] && xi[1][0] > xi[2][0])) {  //vertical
      cout<<" it is a vertical parabola .......@.......@........@........"<<endl;

//       p1 = { static_cast<Type>(xi[0][0]), static_cast<Type>(xi[0][1]) };
//       p2 = { static_cast<Type>(xi[2][0]), static_cast<Type>(xi[2][1]) };
//       p3 = { static_cast<Type>(xi[1][0]), static_cast<Type>(xi[1][1]) };

      q1 = { xi[0][0], xi[0][1] };
      q2 = { xi[2][0], xi[2][1] };
      q3 = { xi[1][0], xi[1][1] };

      Parabola <double> parabola = get_parabola_equation(q1, q2, q3);
      int normal;

      cout << " unit box sign = {" ;
      for(unsigned l = 0; l < 4 ; l++) {
//         unitxvsign[l] = ((static_cast<double>(parabola.k) * unitxv[0][l] * unitxv[0][l] + static_cast<double>(parabola.b) * unitxv[0][l] + static_cast<double>(parabola.d) + unitxv[1][l]) > 0) ? 1 : -1;
        unitxvsign[l] = ((parabola.k * unitxv[0][l] * unitxv[0][l] + parabola.b * unitxv[0][l] + parabola.d + unitxv[1][l]) > 0) ? 1 : -1;
        cout << unitxvsign[l] << ", ";
      }
      cout << "} " << endl;

      cout <<  "( " << q1.x << "," << q1.y << " )" << " , ( " << q2.x << "," << q2.y << " )" << " , ( " << q3.x << "," << q3.y << " ) " << endl;
      cout << parabola.k << "x^2+ " << parabola.b << "x+ " << parabola.d << "+y =0 " << endl;

      normal = checkVectorRelation(xvsign, unitxvsign);
      cout << " normal = " << normal << endl;

      unsigned table_number;

//       CheckIntersection<Type>(intersect_number, table_number, intersection, interp_point, parabola);
//       cout << " table number = " << table_number << endl;
      q3.x = (q1.x + q2.x) / 2;
      q3.y = -parabola.k * q3.x * q3.x - parabola.b * q3.x - parabola.d ;
      cout <<  "( " << q1.x << "," << q1.y << " )" << " , ( " << q2.x << "," << q2.y << " )" << " , ( " << q3.x << "," << q3.y << " ) " << endl;


//       if(interp_point.size() == 2) {

      find_search_table(q1, q2, q3, table_number, searchP);

      cout << " table number = " << table_number << endl;


//         Point3D searchP(static_cast<double>(interp_point[0]), static_cast<double>(interp_point[1]), static_cast<double>(p3.y));
        OctreeNode<Type>* result = loadedRoots[table_number].search(searchP);
        if(result) {
          std::cout << "Found the smallest sub-cube containing the search point." << std::endl;
          std::cout << "\nSearch Point: (" << searchP.x << ", " << searchP.y << ", " << searchP.z << ")\n";
          std::cout << "Smallest Sub-cube Bounds: ";
          std::cout << "(" << result->minBounds.x << ", " << result->minBounds.y << ", " << result->minBounds.z << ") to ";
          std::cout << "(" << result->maxBounds.x << ", " << result->maxBounds.y << ", " << result->maxBounds.z << ")\n";
          std::cout << "depth : = " << result->depth << " \n";

          std::vector<double>interp_point = {searchP.x, searchP.y, searchP.z};
          std::cout << "\n interp Point: (" << interp_point[0] << ", " << interp_point[1] << ", " << interp_point[2] << ")\n";

          trilinier_interpolation_vector(result->corners, result->cornerWeights, interp_point, interp_point_weights);
          std::vector<double>modified_weights(interp_point_weights.size());
          if(normal == -1) {
            for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
//               modified_weights[aq] = 1 - interp_point_weights[interp_point_weights.size()-1-aq];
              modified_weights[aq] = 1 - interp_point_weights[aq];
            }
          }
          else modified_weights = interp_point_weights;
// modified_weights = interp_point_weights;

          std::cout << "AAAAA original weight\n";
          for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
            std::cout << interp_point_weights[aq] << " ";
          }
          std::cout << std::endl;
          std::cout << "AAAAA modified weight\n";
          for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
            std::cout << modified_weights[aq] << " ";
          }
          std::cout << std::endl;

          std::vector<double> phi, gradPhi;
          std::vector<double> Xg(femQuad->GetGaussPointNumber(),0);
          std::vector<double> Yg(femQuad->GetGaussPointNumber(),0);
          std::vector<double> Jg(femQuad->GetGaussPointNumber(),0);
          for(unsigned ig = 0; ig < femQuad->GetGaussPointNumber(); ig++) {
            // *** get gauss point weight, test function and test function partial derivatives ***
            femQuad->Jacobian(xv, ig, Jg[ig], phi, gradPhi);
            for(unsigned i =0;i<phi.size();i++){
              Xg[ig] += phi[i]*xv[0][i];
              Yg[ig] += phi[i]*xv[1][i];
            }
            //std::cout <<ig<<" "<< xg[ig] <<" "<<Xg[ig]<<" "<< yg[ig] <<" "<<Yg[ig]<<" "<<Jg[ig]<<std::endl;
          }

          // Area = GaussIntegral(0, 0, xg, yg, interp_point_weights, gaussWeight);

          Area = GaussIntegral(0, 0, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix  = GaussIntegral(1, 0, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Iy  = GaussIntegral(0, 1, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ixy  = GaussIntegral(1, 1, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix3  = GaussIntegral(3, 0, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix2y  = GaussIntegral(2, 1, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ixy2  = GaussIntegral(1, 2, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Iy3 = GaussIntegral(0, 3, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix2y2  = GaussIntegral(2, 2, Xg.data(), Yg.data(), modified_weights, Jg.data());


          std::cout << "Area = " << Area << std::endl;
          std::cout << "Ix = " << Ix << std::endl;
          std::cout << "Iy = " << Iy << std::endl;
          std::cout << "Ixy = " << Ixy << std::endl;
          std::cout << "Ix3 = " << Ix3 << std::endl;
          std::cout << "Ix2y = " << Ix2y << std::endl;
          std::cout << "Ixy2 = " << Ixy2 << std::endl;
          std::cout << "Iy3 = " << Iy3 << std::endl;
          std::cout << "Ix2y2 = " << Ix2y2 << std::endl;

/*
          Type AArea = find_area_2intersection_formula(0, 0, s, a, c, table_number, p1, p2, p3);

          Pweights(s, a, c, table_number, p1, p2, p3, weightCF);


//         std::cout << "corner points:\n";
//         for (unsigned ig = 0; ig < result->corners.size(); ig++) {
//             std::cout << "(" << result->corners[ig][0] << ", " << result->corners[ig][1] << ", " << result->corners[ig][2] << ") : ";
//               std::cout << result->cornerAreas[ig][0] << ", " << result->cornerAreas[ig][1] << ", " << result->cornerAreas[ig][2] << " ; "<<endl;
//               PointT <Type> p1, p2, p3 ;
//               get_p1_p2_p3(table, interp_point, p1, p2, p3);
//
//         }*/


          trilinier_interpolation_vector(result->corners, result->cornerAreas, interp_point, interp_point_weights);  // interpolating the integrals from corners.

          cout << "\n interpolated integrals " ;
          for(unsigned ig = 0; ig < interp_point_weights.size(); ig++) {
            cout <<  interp_point_weights[ig] << " " ;
          }
          cout << endl;
/*

          cout << " weightCF = " ;
          for(unsigned ig = 0; ig < weightCF.size(); ig++) {
            cout <<  weightCF[ig] << " " ;
          }
          cout << endl;

          cout << " Jacobian = " ;
          for(unsigned ig = 0; ig < Jg.size(); ig++) {
            cout <<  Jg[ig] << " " ;
          }
          cout << endl;

          cout << " Analytic area = " << AArea << endl ;*/
        }
        else {
          std::cout << "Search point not found in the Octree." << std::endl;
        }
//       }
    }



    if((xi[0][1] < xi[1][1] && xi[1][1] < xi[2][1]) || (xi[0][1] > xi[1][1] && xi[1][1] > xi[2][1])) { //horizontal

      cout<<" it is a Horizontal parabola .......@.......@........@........"<<endl;
//       p1 = { static_cast<Type>(xi[0][1]), static_cast<Type>(xi[0][0]) };
//       p2 = { static_cast<Type>(xi[2][1]), static_cast<Type>(xi[2][0]) };
//       p3 = { static_cast<Type>(xi[1][1]), static_cast<Type>(xi[1][0]) };

      q1 = { xi[0][1], xi[0][0] };
      q2 = { xi[2][1], xi[2][0] };
      q3 = { xi[1][1], xi[1][0] };

      Parabola <double> parabola = get_parabola_equation(q1, q2, q3);
      int normal;

      //use horizotal parabola for the normal
      cout << " unit box sign = {" ;
      for(unsigned l = 0; l < 4 ; l++) {
        unitxvsign[l] = ((static_cast<double>(parabola.k) * unitxv[1][l] * unitxv[1][l] + static_cast<double>(parabola.b) * unitxv[1][l] + static_cast<double>(parabola.d) + unitxv[0][l]) > 0) ? 1 : -1;
        cout << unitxvsign[l] << ", ";
      }
      cout << "} " << endl;

      cout <<  "( " << q1.x << "," << q1.y << " )" << " , ( " << q2.x << "," << q2.y << " )" << " , ( " << q3.x << "," << q3.y << " ) " << endl;
      cout << parabola.k << "y^2+ " << parabola.b << "y+ " << parabola.d << "+x =0 " << endl;

      normal = checkVectorRelation(xvsign, unitxvsign);
      cout << " normal = " << normal << endl ;

      unsigned table_number;
      q3.x = (q1.x + q2.x) / 2.;
      q3.y = -parabola.k * q3.x * q3.x - parabola.b * q3.x - parabola.d ;
      cout <<  "( " << q1.x << "," << q1.y << " )" << " , ( " << q2.x << "," << q2.y << " )" << " , ( " << q3.x << "," << q3.y << " ) " << endl;

      find_search_table(q1, q2, q3, table_number, searchP);

      cout << " table number = " << table_number << endl;

//       if(interp_point.size() == 2) {
//         Point3D searchP(static_cast<double>(interp_point[0]), static_cast<double>(interp_point[1]), static_cast<double>(p3.y));
        OctreeNode<Type>* result = loadedRoots[table_number].search(searchP);
        if(result) {
          std::cout << "Found the smallest sub-cube containing the search point." << std::endl;
          std::cout << "\nSearch Point: (" << searchP.x << ", " << searchP.y << ", " << searchP.z << ")\n";
          std::cout << "Smallest Sub-cube Bounds: ";
          std::cout << "(" << result->minBounds.x << ", " << result->minBounds.y << ", " << result->minBounds.z << ") to ";
          std::cout << "(" << result->maxBounds.x << ", " << result->maxBounds.y << ", " << result->maxBounds.z << ")\n";
          std::cout << "depth : = " << result->depth << " \n";

          std::vector<double>interp_point = {searchP.x, searchP.y, searchP.z};
          std::cout << "\n interp Point: (" << interp_point[0] << ", " << interp_point[1] << ", " << interp_point[2] << ")\n";

          trilinier_interpolation_vector(result->corners, result->cornerWeights, interp_point, interp_point_weights);

          std::vector<double>modified_weights(interp_point_weights.size());



          if (table_number == 2 || table_number == 4){
            if(normal == -1) {
              for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
//                 modified_weights[aq] = 1 - interp_point_weights[interp_point_weights.size()-1-aq];   // Originally I use this. I changed it to the bottom one.
                modified_weights[aq] = 1 - interp_point_weights[aq];
              }
            }
            else{
              modified_weights = interp_point_weights;
//               for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
//                 modified_weights[aq] = interp_point_weights[interp_point_weights.size()-1-aq];
//               }
            }
          }

          else if (table_number == 0 || table_number == 6){
            if(normal == -1) {
              for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
                modified_weights[aq] = 1 - interp_point_weights[interp_point_weights.size()-1-aq];
//                 modified_weights[aq] = 1 - interp_point_weights[aq];
              }
            }
            else{
//               modified_weights = interp_point_weights;
              for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
                modified_weights[aq] = interp_point_weights[interp_point_weights.size()-1-aq];
              }
            }
          }

          else if (table_number == 1) {

            for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
//                 modified_weights[aq] = 1 - interp_point_weights[interp_point_weights.size()-1-aq];
                modified_weights[aq] = interp_point_weights[interp_point_weights.size()-1-aq];
//                 modified_weights[aq] = 1 - interp_point_weights[aq];
              }
          }

          else{
            if(normal == -1) {
              for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
                modified_weights[aq] = 1 - interp_point_weights[aq];
              }
            }
            else{
              modified_weights = interp_point_weights;
            }
          }




          std::cout << "AAAAA original weight\n";
          for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
            std::cout << interp_point_weights[aq] << " ";
          }
          std::cout << std::endl;
          std::cout << "BBBBB modified weight\n";
          for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
            std::cout << modified_weights[aq] << " ";
          }
          std::cout << std::endl;


          std::vector<double> phi, gradPhi;
          std::vector<double> Xg(femQuad->GetGaussPointNumber(),0);
          std::vector<double> Yg(femQuad->GetGaussPointNumber(),0);
          std::vector<double> Jg(femQuad->GetGaussPointNumber(),0);
          for(unsigned ig = 0; ig < femQuad->GetGaussPointNumber(); ig++) {
            // *** get gauss point weight, test function and test function partial derivatives ***
            femQuad->Jacobian(xv, ig, Jg[ig], phi, gradPhi);
            for(unsigned i =0;i<phi.size();i++){
              Xg[ig] += phi[i]*xv[0][i];
              Yg[ig] += phi[i]*xv[1][i];
            }
            //std::cout <<ig<<" "<< xg[ig] <<" "<<Xg[ig]<<" "<< yg[ig] <<" "<<Yg[ig]<<" "<<Jg[ig]<<std::endl;
          }

          // Area = GaussIntegral(0, 0, xg, yg, interp_point_weights, gaussWeight);
/*
          Area = GaussIntegral(0, 0, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix  = GaussIntegral(0, 1, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Iy  = GaussIntegral(1, 0, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ixy  = GaussIntegral(1, 1, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix3  = GaussIntegral(0, 3, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix2y  = GaussIntegral(1, 2, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ixy2  = GaussIntegral(2, 1, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Iy3 = GaussIntegral(3, 0, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix2y2  = GaussIntegral(2, 2, Xg.data(), Yg.data(), modified_weights, Jg.data());*/

          Area = GaussIntegral(0, 0, Yg.data(), Xg.data(), modified_weights, Jg.data());
          Ix  = GaussIntegral(0, 1,  Yg.data(), Xg.data(), modified_weights, Jg.data());
          Iy  = GaussIntegral(1, 0,  Yg.data(), Xg.data(), modified_weights, Jg.data());
          Ixy  = GaussIntegral(1, 1, Yg.data(), Xg.data(), modified_weights, Jg.data());
          Ix3  = GaussIntegral(0, 3, Yg.data(), Xg.data(), modified_weights, Jg.data());
          Ix2y  = GaussIntegral(1, 2,Yg.data(), Xg.data(), modified_weights, Jg.data());
          Ixy2  = GaussIntegral(2, 1,Yg.data(), Xg.data(), modified_weights, Jg.data());
          Iy3 = GaussIntegral(3, 0, Yg.data(), Xg.data(), modified_weights, Jg.data());
          Ix2y2  = GaussIntegral(2, 2, Yg.data(), Xg.data(), modified_weights, Jg.data());

          std::cout << "Area = " << Area << std::endl;
          std::cout << "Ix = " << Ix << std::endl;
          std::cout << "Iy = " << Iy << std::endl;
          std::cout << "Ixy = " << Ixy << std::endl;
          std::cout << "Ix3 = " << Ix3 << std::endl;
          std::cout << "Ix2y = " << Ix2y << std::endl;
          std::cout << "Ixy2 = " << Ixy2 << std::endl;
          std::cout << "Iy3 = " << Iy3 << std::endl;
          std::cout << "Ix2y2 = " << Ix2y2 << std::endl;

//           Type AArea = find_area_2intersection_formula(0, 0, s, a, c, table_number, p1, p2, p3);
//           Pweights(s, a, c, table_number, p1, p2, p3, weightCF);
//
//
//           std::cout << "corner points:\n";
//           for (unsigned ig = 0; ig < result->corners.size(); ig++) {
//               std::cout << "(" << result->corners[ig][0] << ", " << result->corners[ig][1] << ", " << result->corners[ig][2] << ") : ";
//                 std::cout << result->cornerAreas[ig][0] << ", " << result->cornerAreas[ig][1] << ", " << result->cornerAreas[ig][2] << " ; "<<endl;
//                 PointT <Type> p1, p2, p3 ;
//                 get_p1_p2_p3(table_number, interp_point, p1, p2, p3);
//           }


          trilinier_interpolation_vector(result->corners, result->cornerAreas, interp_point, interp_point_weights);  // interpolating the integrals from corners.

          cout << "\n interpolated area " ;
          for(unsigned ig = 0; ig < interp_point_weights.size(); ig++) {
            cout <<  interp_point_weights[ig] << " " ;
          }
          cout << endl;
/*

          cout << " weightCF = " ;
          for(unsigned ig = 0; ig < weightCF.size(); ig++) {
            cout <<  weightCF[ig] << " " ;
          }
          cout << endl;

          cout << " 1 - weightCF = " ;
          for(unsigned ig = 0; ig < weightCF.size(); ig++) {
            cout <<  1 - weightCF[ig] << " " ;
          }
          cout << endl;

          cout << " Analytic area = " << AArea << endl ;*/
        }
        else {
          std::cout << "Search point not found in the Octree." << std::endl;
        }
//       }

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

  cout << "========================================>    " <<endl;

//   std::vector<double> interpolated_Weight_CF = find_Weight_CF<Type>(loadedRoots, xv, A);
//   cout << " interpolated weightCF = " <<endl;
//
//           for(unsigned ig = 0; ig < interpolated_Weight_CF.size(); ig++) {
//             cout <<  interpolated_Weight_CF[ig] << " " ;
//           }


  return 0;
}






















template <class Type>
std::vector<double> find_Weight_CF( std::vector<OctreeNode<Type>> &loadedRoots, const std::vector<std::vector<double>> &xv, const std::vector<double> &A){


    unsigned nInt;
    std::vector<std::vector<double>> unitxv = {{0., 1., 1., 0.}, {0., 0., 1., 1.}};
    unsigned nPoints = 3;
    unsigned dim = 2;
    short unsigned ielType = 3; //quad
    unsigned femType = 0; //linear FEM
    std::vector< double > interp_point_weights;
    PointT <Type> p1, p2, p3;
    std::vector<double>modified_weights;

    Fem fem = Fem(3 * 2, 2);
    unsigned quad = 3;
    unsigned linear = 0;
    const elem_type *femQuad = fem.GetFiniteElement(quad, linear);

    std::pair<std::vector<std::vector<double>>, std::vector<double>> xp = GetCellPointsFromQuadric(xv, A, nPoints, nInt);     //This fins the points in physical space

    std::vector < std::vector < std::vector <double > > > aP(1);
    ProjectNodalToPolynomialCoefficients(aP[femType], xv, ielType, femType);

    std::vector<int> xvsign(4);
    std::vector<int> unitxvsign(4);

    std::vector<std::vector<double>> xi(nPoints, std::vector<double>(2, 0.));
    for(unsigned i = 0; i < nPoints; i++) {
      bool inverseMapping = GetInverseMapping(femType, ielType, aP, xp.first[i], xi[i], 100);        //This maps the phsical points to {(-1,-1),(1,1)} box
//         std::cout << " \nx[i] physical value " << i << " " << xp.first[i][0] << " " << xp.first[i][1] << std::endl;
//         std::cout << " x[i] value in (-1,1) " << i << " " << xi[i][0] << " " << xi[i][1] << std::endl;
      xi[i] = {0.5 * (xi[i][0] + 1.), 0.5 * (xi[i][1] + 1.)};                                        // //This maps the points to unit box
//       std::cout << "value in unit box" << i << " " << xi[i][0] << " " << xi[i][1] << std::endl;
    }





        //finiding monotone
//     if((xi[0][0] < xi[1][0] && xi[1][0] < xi[2][0]) || (xi[0][0] > xi[1][0] && xi[1][0] > xi[2][0])) cout << " monotonous in x " << endl;
//     else cout << " non-monotonous in x " << endl;
//
//     if((xi[0][1] < xi[1][1] && xi[1][1] < xi[2][1]) || (xi[0][1] > xi[1][1] && xi[1][1] > xi[2][1])) cout << " monotonous in y " << endl;
//     else cout << " non-monotonous in y " << endl;

//     cout << " xv sign = {" ;
    for(unsigned l = 0; l < 4 ; l++) {
      xvsign[l] = ((A[0] * xv[0][l] * xv[0][l] + A[1] * xv[0][l] * xv[1][l] + A[2] * xv[1][l] * xv[1][l] + A[3] * xv[0][l] + A[4] * xv[1][l] + A[5]) >= 0) ? 1 : -1 ;
//       cout << xvsign[l] << ", ";
    }
//     cout << "} " << endl;


    if((xi[0][0] < xi[1][0] && xi[1][0] < xi[2][0]) || (xi[0][0] > xi[1][0] && xi[1][0] > xi[2][0])) {  //vertical


      p1 = { static_cast<Type>(xi[0][0]), static_cast<Type>(xi[0][1]) };
      p2 = { static_cast<Type>(xi[2][0]), static_cast<Type>(xi[2][1]) };
      p3 = { static_cast<Type>(xi[1][0]), static_cast<Type>(xi[1][1]) };

      Parabola <Type> parabola = get_parabola_equation(p1, p2, p3);
      int normal;

//       cout << " unit box sign = {" ;
      for(unsigned l = 0; l < 4 ; l++) {
        unitxvsign[l] = ((static_cast<double>(parabola.k) * unitxv[0][l] * unitxv[0][l] + static_cast<double>(parabola.b) * unitxv[0][l] + static_cast<double>(parabola.d) + unitxv[1][l]) > 0) ? 1 : -1;
//         cout << unitxvsign[l] << ", ";
      }
//       cout << "} " << endl;

//       cout <<  "( " << p1.x << "," << p1.y << " )" << " , ( " << p2.x << "," << p2.y << " )" << " , ( " << p3.x << "," << p3.y << " ) " << endl;
//       cout << parabola.k << "x^2+ " << parabola.b << "x+ " << parabola.d << "+y =0 " << endl;

      normal = checkVectorRelation(xvsign, unitxvsign);
//       cout << " normal = " << normal << endl;

      int intersect_number;
      unsigned table_number;
      std::vector <Type> intersection;
      std::vector <Type> interp_point;
      CheckIntersection<Type>(intersect_number, table_number, intersection, interp_point, parabola);
//       cout << " table number = " << table_number << endl;
      p3.x = (p1.x + p2.x) / 2;
      p3.y = -parabola.k * p3.x * p3.x - parabola.b * p3.x - parabola.d ;
//       cout <<  "( " << p1.x << "," << p1.y << " )" << " , ( " << p2.x << "," << p2.y << " )" << " , ( " << p3.x << "," << p3.y << " ) " << endl;

      if(interp_point.size() == 2) {
        Point3D searchP(static_cast<double>(interp_point[0]), static_cast<double>(interp_point[1]), static_cast<double>(p3.y));
        OctreeNode<Type>* result = loadedRoots[table_number].search(searchP);
        if(result) {
//           std::cout << "Found the smallest sub-cube containing the search point." << std::endl;
//           std::cout << "\nSearch Point: (" << searchP.x << ", " << searchP.y << ", " << searchP.z << ")\n";
//           std::cout << "Smallest Sub-cube Bounds: ";
//           std::cout << "(" << result->minBounds.x << ", " << result->minBounds.y << ", " << result->minBounds.z << ") to ";
//           std::cout << "(" << result->maxBounds.x << ", " << result->maxBounds.y << ", " << result->maxBounds.z << ")\n";
//           std::cout << "depth : = " << result->depth << " \n";

          std::vector<double>interp_point = {searchP.x, searchP.y, searchP.z};
//           std::cout << "\n interp Point: (" << interp_point[0] << ", " << interp_point[1] << ", " << interp_point[2] << ")\n";

          trilinier_interpolation_vector(result->corners, result->cornerWeights, interp_point, interp_point_weights);
          modified_weights.resize(interp_point_weights.size());
          if(normal == -1) {
            for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
//               modified_weights[aq] = 1 - interp_point_weights[interp_point_weights.size()-1-aq];
              modified_weights[aq] = 1 - interp_point_weights[aq];
            }
          }
          else modified_weights = interp_point_weights;
// modified_weights = interp_point_weights;

          std::cout << "AAAA\n";
          for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
            std::cout << modified_weights[aq] << " ";
          }
          std::cout << std::endl;

          std::vector<double> phi, gradPhi;
          std::vector<double> Xg(femQuad->GetGaussPointNumber(),0);
          std::vector<double> Yg(femQuad->GetGaussPointNumber(),0);
          std::vector<double> Jg(femQuad->GetGaussPointNumber(),0);
          for(unsigned ig = 0; ig < femQuad->GetGaussPointNumber(); ig++) {
            // *** get gauss point weight, test function and test function partial derivatives ***
            femQuad->Jacobian(xv, ig, Jg[ig], phi, gradPhi);
            for(unsigned i =0;i<phi.size();i++){
              Xg[ig] += phi[i]*xv[0][i];
              Yg[ig] += phi[i]*xv[1][i];
            }
            //std::cout <<ig<<" "<< xg[ig] <<" "<<Xg[ig]<<" "<< yg[ig] <<" "<<Yg[ig]<<" "<<Jg[ig]<<std::endl;
          }

          // Area = GaussIntegral(0, 0, xg, yg, interp_point_weights, gaussWeight);

/*
          double Area0 = 0, Area = 0, Ix = 0, Iy = 0, Ix3 = 0, Ix2y = 0, Ixy2 = 0, Iy3 = 0, Ix2y2 = 0;
          Area = GaussIntegral(0, 0, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix  = GaussIntegral(1, 0, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Iy  = GaussIntegral(0, 1, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix3  = GaussIntegral(3, 0, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix2y  = GaussIntegral(2, 1, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ixy2  = GaussIntegral(1, 2, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Iy3 = GaussIntegral(0, 3, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix2y2  = GaussIntegral(2, 2, Xg.data(), Yg.data(), modified_weights, Jg.data());


          std::cout << "Area = " << Area << std::endl;
          std::cout << "Ix = " << Ix << std::endl;
          std::cout << "Iy = " << Iy << std::endl;
          std::cout << "Ix3 = " << Ix3 << std::endl;
          std::cout << "Ix2y = " << Ix2y << std::endl;
          std::cout << "Ixy2 = " << Ixy2 << std::endl;
          std::cout << "Iy3 = " << Iy3 << std::endl;
          std::cout << "Ix2y2 = " << Ix2y2 << std::endl;*/
        }
        else {
          std::cout << "Search point not found in the Octree." << std::endl;
        }
      }
    }



    else if((xi[0][1] < xi[1][1] && xi[1][1] < xi[2][1]) || (xi[0][1] > xi[1][1] && xi[1][1] > xi[2][1])) { //horizontal

      p1 = { static_cast<Type>(xi[0][1]), static_cast<Type>(xi[0][0]) };
      p2 = { static_cast<Type>(xi[2][1]), static_cast<Type>(xi[2][0]) };
      p3 = { static_cast<Type>(xi[1][1]), static_cast<Type>(xi[1][0]) };

      Parabola <Type> parabola = get_parabola_equation(p1, p2, p3);
      int normal;

      //use horizotal parabola for the normal
//       cout << " unit box sign = {" ;
      for(unsigned l = 0; l < 4 ; l++) {
        unitxvsign[l] = ((static_cast<double>(parabola.k) * unitxv[1][l] * unitxv[1][l] + static_cast<double>(parabola.b) * unitxv[1][l] + static_cast<double>(parabola.d) + unitxv[0][l]) > 0) ? 1 : -1;
//         cout << unitxvsign[l] << ", ";
      }
//       cout << "} " << endl;

//       cout <<  "( " << p1.y << "," << p1.x << " )" << " , ( " << p2.y << "," << p2.x << " )" << " , ( " << p3.y << "," << p3.x << " ) " << endl;
//       cout << parabola.k << "x^2+ " << parabola.b << "x+ " << parabola.d << "+y =0 " << endl;

      normal = checkVectorRelation(xvsign, unitxvsign);
//       cout << " normal = " << normal << endl ;


      int intersect_number;
      unsigned table_number;
      std::vector <Type> intersection;
      std::vector <Type> interp_point;
      CheckIntersection<Type>(intersect_number, table_number, intersection, interp_point, parabola);
//       cout << " table number = " << table_number << endl;

      p3.x = (p1.x + p2.x) / 2;
      p3.y = -parabola.k * p3.x * p3.x - parabola.b * p3.x - parabola.d ;
//       cout <<  "( " << p1.x << "," << p1.y << " )" << " , ( " << p2.x << "," << p2.y << " )" << " , ( " << p3.x << "," << p3.y << " ) " << endl;

      if(interp_point.size() == 2) {

        Point3D searchP(static_cast<double>(interp_point[0]), static_cast<double>(interp_point[1]), static_cast<double>(p3.y));
        OctreeNode<Type>* result = loadedRoots[table_number].search(searchP);
        if(result) {
//           std::cout << "Found the smallest sub-cube containing the search point." << std::endl;
//           std::cout << "\nSearch Point: (" << searchP.x << ", " << searchP.y << ", " << searchP.z << ")\n";
//           std::cout << "Smallest Sub-cube Bounds: ";
//           std::cout << "(" << result->minBounds.x << ", " << result->minBounds.y << ", " << result->minBounds.z << ") to ";
//           std::cout << "(" << result->maxBounds.x << ", " << result->maxBounds.y << ", " << result->maxBounds.z << ")\n";
//           std::cout << "depth : = " << result->depth << " \n";

          std::vector<double>interp_point = {searchP.x, searchP.y, searchP.z};
//           std::cout << "\n interp Point: (" << interp_point[0] << ", " << interp_point[1] << ", " << interp_point[2] << ")\n";

          trilinier_interpolation_vector(result->corners, result->cornerWeights, interp_point, interp_point_weights);

          modified_weights.resize(interp_point_weights.size());


          if (table_number == 2 || table_number == 4){
            if(normal == -1) {
              for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
                modified_weights[aq] = 1 - interp_point_weights[interp_point_weights.size()-1-aq];
  //               modified_weights[aq] = 1 - interp_point_weights[aq];
              }
            }
            else{
  //             modified_weights = interp_point_weights;
              for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
                modified_weights[aq] = interp_point_weights[interp_point_weights.size()-1-aq];
              }
            }
          }

          else if (table_number == 1) {

            for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
//                 modified_weights[aq] = 1 - interp_point_weights[interp_point_weights.size()-1-aq];
                modified_weights[aq] = 1 - interp_point_weights[aq];
              }
          }

          else{
            if(normal == -1) {
              for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
                modified_weights[aq] = 1 - interp_point_weights[aq];
              }
            }
            else{
              modified_weights = interp_point_weights;
            }
          }




// //           std::cout << "BBBBB\n";
//           for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
// //             std::cout << modified_weights[aq] << " ";
//           }
// //           std::cout << std::endl;


          std::vector<double> phi, gradPhi;
          std::vector<double> Xg(femQuad->GetGaussPointNumber(),0);
          std::vector<double> Yg(femQuad->GetGaussPointNumber(),0);
          std::vector<double> Jg(femQuad->GetGaussPointNumber(),0);
          for(unsigned ig = 0; ig < femQuad->GetGaussPointNumber(); ig++) {
            // *** get gauss point weight, test function and test function partial derivatives ***
            femQuad->Jacobian(xv, ig, Jg[ig], phi, gradPhi);
            for(unsigned i =0;i<phi.size();i++){
              Xg[ig] += phi[i]*xv[0][i];
              Yg[ig] += phi[i]*xv[1][i];
            }
            //std::cout <<ig<<" "<< xg[ig] <<" "<<Xg[ig]<<" "<< yg[ig] <<" "<<Yg[ig]<<" "<<Jg[ig]<<std::endl;
          }
/*
          // Area = GaussIntegral(0, 0, xg, yg, interp_point_weights, gaussWeight);
          double Area0 = 0, Area = 0, Ix = 0, Iy = 0, Ix3 = 0, Ix2y = 0, Ixy2 = 0, Iy3 = 0, Ix2y2 = 0;
          Area = GaussIntegral(0, 0, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix  = GaussIntegral(0, 1, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Iy  = GaussIntegral(1, 0, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix3  = GaussIntegral(0, 3, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix2y  = GaussIntegral(1, 2, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ixy2  = GaussIntegral(2, 1, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Iy3 = GaussIntegral(3, 0, Xg.data(), Yg.data(), modified_weights, Jg.data());
          Ix2y2  = GaussIntegral(2, 2, Xg.data(), Yg.data(), modified_weights, Jg.data());

          std::cout << "Area = " << Area << std::endl;
          std::cout << "Ix = " << Ix << std::endl;
          std::cout << "Iy = " << Iy << std::endl;
          std::cout << "Ix3 = " << Ix3 << std::endl;
          std::cout << "Ix2y = " << Ix2y << std::endl;
          std::cout << "Ixy2 = " << Ixy2 << std::endl;
          std::cout << "Iy3 = " << Iy3 << std::endl;
          std::cout << "Ix2y2 = " << Ix2y2 << std::endl;*/

//           Type AArea = find_area_2intersection_formula(0, 0, s, a, c, table_number, p1, p2, p3);
//           Pweights(s, a, c, table_number, p1, p2, p3, weightCF);
//
//
//           std::cout << "corner points:\n";
//           for (unsigned ig = 0; ig < result->corners.size(); ig++) {
//               std::cout << "(" << result->corners[ig][0] << ", " << result->corners[ig][1] << ", " << result->corners[ig][2] << ") : ";
//                 std::cout << result->cornerAreas[ig][0] << ", " << result->cornerAreas[ig][1] << ", " << result->cornerAreas[ig][2] << " ; "<<endl;
//                 PointT <Type> p1, p2, p3 ;
//                 get_p1_p2_p3(table_number, interp_point, p1, p2, p3);
//           }

/*
          trilinier_interpolation_vector(result->corners, result->cornerAreas, interp_point, interp_point_weights);  // interpolating the integrals from corners.

          cout << "\n interpolated area " ;
          for(unsigned ig = 0; ig < interp_point_weights.size(); ig++) {
            cout <<  interp_point_weights[ig] << " " ;
          }
          cout << endl;*/
/*

          cout << " weightCF = " ;
          for(unsigned ig = 0; ig < weightCF.size(); ig++) {
            cout <<  weightCF[ig] << " " ;
          }
          cout << endl;

          cout << " 1 - weightCF = " ;
          for(unsigned ig = 0; ig < weightCF.size(); ig++) {
            cout <<  1 - weightCF[ig] << " " ;
          }
          cout << endl;

          cout << " Analytic area = " << AArea << endl ;*/
        }
        else {
          std::cout << "Search point not found in the Octree." << std::endl;
        }
      }

    }

    return modified_weights ;

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
