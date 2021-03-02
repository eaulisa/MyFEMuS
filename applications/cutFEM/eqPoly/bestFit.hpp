


#ifndef __femus_bestFit_hpp__
#define __femus_bestFit_hpp__

#include <iostream>
#include <iomanip>

#include "./LiSK/lisk.hpp"
#include <cmath>       /* exp */
#include "eqPoly.hpp"
#include </usr/include/eigen3/Eigen/Core>
#include </usr/include/eigen3/Eigen/SVD>


//using namespace femus;
using namespace std;
using namespace Eigen;



class BestFit {
  public:
    
      BestFit();


    std::vector < double > FindBestFit(const std::vector < double > &pts, const std::vector < double > &Npts, const unsigned &dim); 
    
    
        
};







#endif



