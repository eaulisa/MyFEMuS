/** tutorial/Ex1
 * This example shows how to:
 * initialize a femus application;
 * define the multilevel-mesh object mlMsh;
 * read from the file ./input/square.neu the coarse-level mesh and associate it to mlMsh;
 * add in mlMsh uniform refined level-meshes;
 * define the multilevel-solution object mlSol associated to mlMsh;
 * add in mlSol different types of finite element solution variables;
 * initialize the solution varables;
 * define vtk and gmv writer objects associated to mlSol;
 * print vtk and gmv binary-format files in ./output directory.
 **/

#include "GaussPoints.hpp"
#include "ElemType.hpp"
#include <cmath>
#include <iostream>


using namespace femus;

int main(int argc, char** args) {


  const double *Gauss[20] = { line_gauss::Gauss0[0],  line_gauss::Gauss1[0], line_gauss::Gauss2[0], line_gauss::Gauss3[0], line_gauss::Gauss4[0],
                              line_gauss::Gauss5[0],  line_gauss::Gauss6[0], line_gauss::Gauss7[0], line_gauss::Gauss8[0], line_gauss::Gauss9[0],
                              line_gauss::Gauss10[0], line_gauss::Gauss11[0], line_gauss::Gauss12[0], line_gauss::Gauss13[0], line_gauss::Gauss14[0],
                              line_gauss::Gauss15[0], line_gauss::Gauss16[0], line_gauss::Gauss17[0], line_gauss::Gauss18[0], line_gauss::Gauss19[0]
                            };

  //QUAD
  unsigned dim = 2;
  unsigned maxNG = 1;

//   std::cout << "const unsigned quad_gauss::GaussPoints[" << maxNG << "] = {";
//   for(unsigned ng = 1; ng <= maxNG; ng++) {
//     unsigned size = pow(ng, dim);
//     std::cout << size << ", ";
//   }
//   std::cout << "\b\b};\n\n";
//
//   std::cout << "const double * quad_gauss::Gauss[" << maxNG << "] = {";
//   for(unsigned ng = 1; ng <= maxNG; ng++) {
//     std::cout << "Gauss" << ng - 1 << "[0], ";
//   }
//   std::cout << "\b\b};\n\n";
//
//   for(unsigned ng = 1; ng <= maxNG; ng++) {
//
//     unsigned size = pow(ng, dim);
//     std::vector < double> weight(size);
//     std::vector < double> x(size);
//     std::vector < double> y(size);
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         weight[i * ng + j] = Gauss[ng - 1][i] * Gauss[ng - 1][j];
//         x[i * ng + j] = Gauss[ng - 1][ng + i];
//         y[i * ng + j] = Gauss[ng - 1][ng + j];
//       }
//     }
//     std::cout.precision(14);
//     std::cout << "const double quad_gauss::Gauss" << ng - 1 << "[" << dim + 1 << "][" <<  size << "] = {{";
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         std::cout << weight[i * ng + j] << ", ";
//       }
//     }
//     std::cout << "\b\b},\n{";
//
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         std::cout << x[i * ng + j] << ", ";
//       }
//     }
//     std::cout << "\b\b},\n{";
//
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         std::cout << y[i * ng + j] << ", ";
//       }
//     }
//     std::cout << "\b\b}};" << std::endl << std::endl;
//
//   }
//
//
//
//   //HEX
//   dim = 3;
//   maxNG = 1;
//
//   std::cout << "const unsigned hex_gauss::GaussPoints[" << maxNG << "] = {";
//   for(unsigned ng = 1; ng <= maxNG; ng++) {
//     unsigned size = pow(ng, dim);
//     std::cout << size << ", ";
//   }
//   std::cout << "\b\b};\n\n";
//
//   std::cout << "const double * hex_gauss::Gauss[" << maxNG << "] = {";
//   for(unsigned ng = 1; ng <= maxNG; ng++) {
//     std::cout << "Gauss" << ng - 1 << "[0], ";
//   }
//   std::cout << "\b\b};\n\n";
//
//   for(unsigned ng = 1; ng <= maxNG; ng++) {
//
//     unsigned size = pow(ng, dim);
//     std::vector < double> weight(size);
//     std::vector < double> x(size);
//     std::vector < double> y(size);
//     std::vector < double> z(size);
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           weight[i * ng * ng + j * ng + k] = Gauss[ng - 1][i] * Gauss[ng - 1][j] * Gauss[ng - 1][k];
//           x[i * ng * ng + j * ng + k] = Gauss[ng - 1][ng + i];
//           y[i * ng * ng + j * ng + k] = Gauss[ng - 1][ng + j];
//           z[i * ng * ng + j * ng + k] = Gauss[ng - 1][ng + k];
//         }
//       }
//     }
//     std::cout.precision(14);
//     std::cout << "const double hex_gauss::Gauss" << ng - 1 << "[" << dim + 1 << "][" <<  size << "] = {{";
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           std::cout << weight[i * ng * ng + j * ng + k]  << ", ";
//         }
//       }
//     }
//     std::cout << "\b\b},\n{";
//
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           std::cout << x[i * ng * ng + j * ng + k] << ", ";
//         }
//       }
//     }
//     std::cout << "\b\b},\n{";
//
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           std::cout << y[i * ng * ng + j * ng + k] << ", ";
//         }
//       }
//     }
//     std::cout << "\b\b},\n{";
//
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           std::cout << z[i * ng * ng + j * ng + k] << ", ";
//         }
//       }
//     }
//     std::cout << "\b\b}};" << std::endl << std::endl;
//
//   }
//
//
// //triangle 1
// 
//     std::vector < std::vector < double > > xv = {{0., 1., 0., 0.}, {0., 0., 1., 1.}}; //old
//     std::vector < std::vector < double > > xv = {{0., 1., 0.5, 0.}, {0., 0., 0.5, 1.}}; //symmetric
//     std::vector <double> phi;  // local test function
//     std::vector <double> phi_x; // local test function first order partial derivatives
//     double jac;
// 
//     const elem_type *fem = new const elem_type_2D("quad", "linear", "seventh");
// 
//     dim = 2;
//     maxNG = 3;
// 
//     std::cout << "const unsigned tri_gauss::GaussPoints[" << maxNG << "] = {";
//     for(unsigned ng = 1; ng <= maxNG; ng++) {
//       unsigned size = pow(ng, dim);
//       std::cout << size << ", ";
//     }
//     std::cout << "\b\b};\n\n";
// 
//     std::cout << "const double * tri_gauss::Gauss[" << maxNG << "] = {";
//     for(unsigned ng = 1; ng <= maxNG; ng++) {
//       std::cout << "Gauss" << ng - 1 << "[0], ";
//     }
//     std::cout << "\b\b};\n\n";
// 
//     for(unsigned ng = 1; ng <= maxNG; ng++) {
// 
//       unsigned size = pow(ng, dim);
//       std::vector < double> weight(size);
//       std::vector < double> x(size);
//       std::vector < double> y(size);
// 
//       std::vector < double > xi(dim);
//       for(unsigned i = 0; i < ng; i++) {
//         for(unsigned j = 0; j < ng; j++) {
//           weight[i * ng + j] = Gauss[ng - 1][i] * Gauss[ng - 1][j];
//           xi[0] = Gauss[ng - 1][ng + i];
//           xi[1] = Gauss[ng - 1][ng + j];
// 
//           fem->Jacobian(xv, xi, jac, phi, phi_x);
//           weight[i * ng + j] *= jac;
// 
//           x[i * ng + j] = 0.;
//           y[i * ng + j] = 0.;
//           for(unsigned ii = 0; ii < 4; ii++) {
//             x[i * ng + j] += phi[ii] * xv[0][ii];
//             y[i * ng + j] += phi[ii] * xv[1][ii];
//           }
// 
// 
//         }
//       }
//       std::cout.precision(14);
//       std::cout << "const double tri_gauss::Gauss" << ng - 1 << "[" << dim + 1 << "][" <<  size << "] = {{";
//       for(unsigned i = 0; i < ng; i++) {
//         for(unsigned j = 0; j < ng; j++) {
//           std::cout << weight[i * ng + j] << ", ";
//         }
//       }
//       std::cout << "\b\b},\n{";
// 
//       for(unsigned i = 0; i < ng; i++) {
//         for(unsigned j = 0; j < ng; j++) {
//           std::cout << x[i * ng + j] << ", ";
//         }
//       }
//       std::cout << "\b\b},\n{";
// 
//       for(unsigned i = 0; i < ng; i++) {
//         for(unsigned j = 0; j < ng; j++) {
//           std::cout << y[i * ng + j] << ", ";
//         }
//       }
//       std::cout << "\b\b}};" << std::endl << std::endl;
// 
//     }
// 
//     delete fem;
//
//

  //triangle 2

  std::vector < std::vector < double > > xv = {{0., 1., 1., 0.}, {0., 0., 0., 1.}};
  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives
  double jac;

  const elem_type *fem = new const elem_type_2D("quad", "linear", "seventh");

  dim = 2;
  maxNG = 20;

  std::cout << "const unsigned tri_gauss::GaussPoints[" << maxNG << "] = {";
  for(unsigned m = 2; m <= maxNG; m++) {
    unsigned size = (m * (m + 1)) / 2 - 1;  
    std::cout << size << ", ";
  }
  std::cout << "\b\b};\n\n";

  std::cout << "const double * tri_gauss::Gauss[" << maxNG << "] = {";
  for(unsigned m = 2; m <= maxNG; m++) {
    std::cout << "Gauss" << m - 1 << "[0], ";
  }
  std::cout << "\b\b};\n\n";

  for(unsigned m = 2; m <= maxNG; m++) {

    unsigned size = (m * (m + 1)) / 2 - 1;
      
    std::vector < double> weight(size);
    std::vector < double> x(size);
    std::vector < double> y(size);

    std::vector < double > xi(dim);

    unsigned cnt = 0;
    for(unsigned i = 0; i < m - 1; i++) {
      for(unsigned j = 0; j < m - i; j++) {
        std::cout << cnt << " ";  
        weight[cnt] = Gauss[m - 2][i] * Gauss[m - 1 - i ][j];
        xi[0] = Gauss[m - 2][m - 1 + i];
        xi[1] = Gauss[m - 1 - i][m - i  + j];

        fem->Jacobian(xv, xi, jac, phi, phi_x);
        weight[cnt] *= jac;

        x[cnt] = 0.;
        y[cnt] = 0.;
        for(unsigned ii = 0; ii < 4; ii++) {
          x[cnt] += phi[ii] * xv[0][ii];
          y[cnt] += phi[ii] * xv[1][ii];
        }
        cnt++;
      }
    }
    std::cout.precision(14);
    std::cout << "const double tri_gauss::Gauss" << m - 1 << "[" << dim + 1 << "][" <<  size << "] = {{";
    for(unsigned k = 0; k < size; k++) {
      std::cout << weight[k] << ", ";
    }
    std::cout << "\b\b},\n{";

    for(unsigned k = 0; k < size; k++) {
      std::cout << x[k] << ", ";

    }
    std::cout << "\b\b},\n{";

    for(unsigned k = 0; k < size; k++) {
      std::cout << y[k] << ", ";
    }
    std::cout << "\b\b}};" << std::endl << std::endl;

  }

  delete fem;
  
  
  
//tetrahedra

//   std::vector < std::vector < double > > xv = {
//       {0, 1, 0.5, 0., 0., 0.5, 1./3., 0}, 
//       {0, 0, 0.5, 1., 0., 0., 1./3., 0.5}, 
//       {0, 0, 0, 0, 1, 0.5, 1./3., 0.5}};
      
//   std::vector < std::vector < double > > xv = {
//       {0., 1., 0., 0., 0., 0., 0., 0.}, 
//       {0., 0., 1., 1., 0., 0., 0., 0.}, 
//       {0., 0., 0., 0., 1., 1., 1., 1.}};    
//       
//   std::vector <double> phi;  // local test function
//   std::vector <double> phi_x; // local test function first order partial derivatives
//   double jac;
// 
//   const elem_type *fem = new const elem_type_3D("hex", "linear", "seventh");
// 
//   dim = 3;
//   maxNG = 3;
// 
//   std::cout << "const unsigned tet_gauss::GaussPoints[" << maxNG << "] = {";
//   for(unsigned ng = 1; ng <= maxNG; ng++) {
//     unsigned size = pow(ng, dim);
//     std::cout << size << ", ";
//   }
//   std::cout << "\b\b};\n\n";
// 
//   std::cout << "const double * tet_gauss::Gauss[" << maxNG << "] = {";
//   for(unsigned ng = 1; ng <= maxNG; ng++) {
//     std::cout << "Gauss" << ng - 1 << "[0], ";
//   }
//   std::cout << "\b\b};\n\n";
// 
//   for(unsigned ng = 1; ng <= maxNG; ng++) {
// 
//     unsigned size = pow(ng, dim);
//     std::vector < double> weight(size);
//     std::vector < double> x(size);
//     std::vector < double> y(size);
//     std::vector < double> z(size);
//     std::vector < double > xi(dim);
// 
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           weight[i * ng * ng + j * ng + k] = Gauss[ng - 1][i] * Gauss[ng - 1][j] * Gauss[ng - 1][k];
//           xi[0] = Gauss[ng - 1][ng + i];
//           xi[1] = Gauss[ng - 1][ng + j];
//           xi[2] = Gauss[ng - 1][ng + k];
// 
// 
//           fem->Jacobian(xv, xi, jac, phi, phi_x);
// 
//           weight[i * ng * ng + j * ng + k] *= jac;
// 
//           x[i * ng * ng + j * ng + k] = 0.;
//           y[i * ng * ng + j * ng + k] = 0.;
//           z[i * ng * ng + j * ng + k] = 0.;
// 
//           for(unsigned ii = 0; ii < 8; ii++) {
//             x[i * ng * ng + j * ng + k] += phi[ii] * xv[0][ii];
//             y[i * ng * ng + j * ng + k] += phi[ii] * xv[1][ii];
//             z[i * ng * ng + j * ng + k] += phi[ii] * xv[2][ii];
//           }
// 
//         }
//       }
//     }
// 
//     std::cout.precision(14);
//     std::cout << "const double tet_gauss::Gauss" << ng - 1 << "[" << dim + 1 << "][" <<  size << "] = {{";
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           std::cout << weight[i * ng * ng + j * ng + k]  << ", ";
//         }
//       }
//     }
//     std::cout << "\b\b},\n{";
// 
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           std::cout << x[i * ng * ng + j * ng + k] << ", ";
//         }
//       }
//     }
//     std::cout << "\b\b},\n{";
// 
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           std::cout << y[i * ng * ng + j * ng + k] << ", ";
//         }
//       }
//     }
//     std::cout << "\b\b},\n{";
// 
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           std::cout << z[i * ng * ng + j * ng + k] << ", ";
//         }
//       }
//     }
//     std::cout << "\b\b}};" << std::endl << std::endl;
// 
//   }
// 
//   delete fem;


  
//   //wedge
// 
//   std::vector < std::vector < double > > xv = {{0, 1, 0, 0, 0, 1, 0, 0}, {0, 0, 1, 1, 0, 0, 1, 1}, {-1, -1, -1, -1, 1, 1, 1, 1}};
//   std::vector <double> phi;  // local test function
//   std::vector <double> phi_x; // local test function first order partial derivatives
//   double jac;
// 
//   const elem_type *fem = new const elem_type_3D("hex", "linear", "seventh");
// 
//   dim = 3;
//   maxNG = 3;
// 
//   std::cout << "const unsigned tet_gauss::GaussPoints[" << maxNG << "] = {";
//   for(unsigned ng = 1; ng <= maxNG; ng++) {
//     unsigned size = pow(ng, dim);
//     std::cout << size << ", ";
//   }
//   std::cout << "\b\b};\n\n";
// 
//   std::cout << "const double * tet_gauss::Gauss[" << maxNG << "] = {";
//   for(unsigned ng = 1; ng <= maxNG; ng++) {
//     std::cout << "Gauss" << ng - 1 << "[0], ";
//   }
//   std::cout << "\b\b};\n\n";
// 
//   for(unsigned ng = 1; ng <= maxNG; ng++) {
// 
//     unsigned size = pow(ng, dim);
//     std::vector < double> weight(size);
//     std::vector < double> x(size);
//     std::vector < double> y(size);
//     std::vector < double> z(size);
//     std::vector < double > xi(dim);
// 
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           weight[i * ng * ng + j * ng + k] = Gauss[ng - 1][i] * Gauss[ng - 1][j] * Gauss[ng - 1][k];
//           xi[0] = Gauss[ng - 1][ng + i];
//           xi[1] = Gauss[ng - 1][ng + j];
//           xi[2] = Gauss[ng - 1][ng + k];
// 
// 
//           fem->Jacobian(xv, xi, jac, phi, phi_x);
// 
//           weight[i * ng * ng + j * ng + k] *= jac;
// 
//           x[i * ng * ng + j * ng + k] = 0.;
//           y[i * ng * ng + j * ng + k] = 0.;
//           z[i * ng * ng + j * ng + k] = 0.;
// 
//           for(unsigned ii = 0; ii < 8; ii++) {
//             x[i * ng * ng + j * ng + k] += phi[ii] * xv[0][ii];
//             y[i * ng * ng + j * ng + k] += phi[ii] * xv[1][ii];
//             z[i * ng * ng + j * ng + k] += phi[ii] * xv[2][ii];
//           }
// 
//         }
//       }
//     }
// 
//     std::cout.precision(14);
//     std::cout << "const double tet_gauss::Gauss" << ng - 1 << "[" << dim + 1 << "][" <<  size << "] = {{";
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           std::cout << weight[i * ng * ng + j * ng + k]  << ", ";
//         }
//       }
//     }
//     std::cout << "\b\b},\n{";
// 
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           std::cout << x[i * ng * ng + j * ng + k] << ", ";
//         }
//       }
//     }
//     std::cout << "\b\b},\n{";
// 
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           std::cout << y[i * ng * ng + j * ng + k] << ", ";
//         }
//       }
//     }
//     std::cout << "\b\b},\n{";
// 
//     for(unsigned i = 0; i < ng; i++) {
//       for(unsigned j = 0; j < ng; j++) {
//         for(unsigned k = 0; k < ng; k++) {
//           std::cout << z[i * ng * ng + j * ng + k] << ", ";
//         }
//       }
//     }
//     std::cout << "\b\b}};" << std::endl << std::endl;
// 
//   }
// 
//   delete fem;








}

