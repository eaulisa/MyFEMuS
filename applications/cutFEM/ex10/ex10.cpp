
#include "CutFemIntegration.hpp"

int main(int, char**) {

  //typedef cpp_bin_float_oct TypeIO;  
  typedef double TypeIO;

  typedef cpp_bin_float_oct TypeA;
  //typedef double TypeA;

  
  std::cout.precision(20);

  {
    unsigned qM = 3;
    CutFemIntegral <TypeIO, TypeA>line   = CutFemIntegral<TypeIO, TypeA>(LINE, qM, "legendre");
    CutFemIntegral <TypeIO, TypeA> quad  = CutFemIntegral<TypeIO, TypeA >(QUAD, qM, "legendre");
    CutFemIntegral <TypeIO, TypeA> tri   = CutFemIntegral<TypeIO, TypeA >(TRI, qM, "legendre");
    CutFemIntegral <TypeIO, TypeA> hex   = CutFemIntegral<TypeIO, TypeA >(HEX, qM, "legendre");
    CutFemIntegral <TypeIO, TypeA> wedge = CutFemIntegral<TypeIO, TypeA >(WEDGE, qM, "legendre");
    CutFemIntegral <TypeIO, TypeA> tet   = CutFemIntegral<TypeIO, TypeA >(TET, qM, "legendre");

//    Line test
    
    std::vector <TypeIO> weightCF;
    line(qM, 0, {-1.},  0,  weightCF);

    const double* weight = line.GetGaussWeightPointer();
    const double* x = line.GetGaussCoordinatePointer(0);

    TypeIO sum = 0.;
    for(unsigned ig = 0; ig < weightCF.size(); ig++) sum += pow(x[ig], 3) * weight[ig] * weightCF[ig];

    std::cout << " sum line test = " << sum << std::endl;
    
//     Quad test
    std::vector <TypeIO> weightCFQuad;
    quad(qM, -1, {1./sqrt(2), 1./sqrt(2)}, 0.,  weightCFQuad);

    const double* weightQ = quad.GetGaussWeightPointer();
    const double* xQ = quad.GetGaussCoordinatePointer(0);
    const double* yQ = quad.GetGaussCoordinatePointer(1);

    sum = 0.;
    for(unsigned ig = 0; ig < weightCFQuad.size(); ig++) {
        sum += pow(xQ[ig], 0) * pow(yQ[ig], 0) * weightQ[ig] * weightCFQuad[ig];
    }

    std::cout << " sum quad test = " << sum << std::endl;
    
    //     Triangle test
    std::vector <TypeIO> weightCFTri;
    tri(qM, 0, {1., -1.}, 0.,  weightCFTri);

    const double* weightT = tri.GetGaussWeightPointer();
    const double* xT = tri.GetGaussCoordinatePointer(0);
    const double* yT = tri.GetGaussCoordinatePointer(1);

    sum = 0.;
    for(unsigned ig = 0; ig < weightCFTri.size(); ig++) {
        sum += pow(xT[ig], 1) * pow(yT[ig], 0) * weightT[ig] * weightCFTri[ig];
    }
    std::cout << " sum tri test = " << sum << std::endl;
    
    //     Hexahedron test
    std::vector <TypeIO> weightCFHex;
    hex(qM, 0, {-0.1, -0.1, +1.}, 0.05,  weightCFHex);

    const double* weightH = hex.GetGaussWeightPointer();
    const double* xH = hex.GetGaussCoordinatePointer(0);
    const double* yH = hex.GetGaussCoordinatePointer(1);
    const double* zH = hex.GetGaussCoordinatePointer(2);

    sum = 0.;
    for(unsigned ig = 0; ig < weightCFHex.size(); ig++) {
        sum += pow(xH[ig], 0) * pow(yH[ig], 1) * pow(zH[ig], 2) * weightH[ig] * weightCFHex[ig];
    }
    std::cout << " sum hex test = " << sum << std::endl;
    
    //     Wedge test
    std::vector <TypeIO> weightCFWed;
    wedge(qM, -1, {-0.1 / sqrt(1.02), 0.1 / sqrt(1.02), 1. / sqrt(1.02)}, 0.,  weightCFWed);

    const double* weightW = wedge.GetGaussWeightPointer();
    const double* xW = wedge.GetGaussCoordinatePointer(0);
    const double* yW = wedge.GetGaussCoordinatePointer(1);
    const double* zW = wedge.GetGaussCoordinatePointer(2);

    sum = 0.;
    for(unsigned ig = 0; ig < weightCFWed.size(); ig++) {
        sum += pow(xW[ig], 0) * pow(yW[ig], 0) * pow(zW[ig], 2) * weightW[ig] * weightCFWed[ig];
    }
    std::cout << " sum wedge test = " << sum << std::endl;
    
    //     Tet test
    std::vector <TypeIO> weightCFTet;
    tet(qM, 0, {-1, -1, -1}, +0.5,  weightCFTet);

    const double* weightTet = tet.GetGaussWeightPointer();
    const double* xTet = tet.GetGaussCoordinatePointer(0);
    const double* yTet = tet.GetGaussCoordinatePointer(1);
    const double* zTet = tet.GetGaussCoordinatePointer(2);

    sum = 0.;
    int m = 1;
    int n = 0;
    int o = 2;
    for(unsigned ig = 0; ig < weightCFTet.size(); ig++) {
        sum += pow(xTet[ig], m) * pow(yTet[ig], n) * pow(zTet[ig], o) * weightTet[ig] * weightCFTet[ig];
    }
    std::cout << " sum tet test = " << sum << std::endl;
    
    

  }
  return 1;

}






























