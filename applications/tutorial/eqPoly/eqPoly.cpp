//
//  main.cpp
//  example
//
//  Created by Sebastian on 18/05/16.
//  Copyright © 2016 Sebastian Kirchner. All rights reserved.
//
//  This file is part of LiSK.
//
//  LiSK is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  LiSK is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with LiSK.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <iomanip>
#include "./LiSK/lisk.hpp"

int main(int argc, const char * argv[]) {

  using namespace std;

  // Double precision example for Li_n and Li_22
  try {
    //Create a LiSK object. Constants are prepared for
    //computations of Li_n up to n=20
    LiSK::LiSK<complex<double>> lisk(20);

    const complex<double> x(1 / 4., 1 / 4.);
    const vector<int> weights = {1, 2, 3, 4, 5, 6, 10, 15, 20};

    // Compute Li_n(x) for n=1,2,3,4,5,6,10,15,20 at x=1/4+I*1/4
    for (auto n : weights) {
      cout << "Li(" << n << ",x) = " << setprecision(17) << lisk.Li(n, x) << endl;
    }

    // Compute Li_22(x,y) at x=y=1/4+I*1/4
    cout << "\nLi_22(x,y) = " << setprecision(17) << lisk.Li22(x, x) << endl << endl;


  }
  catch (runtime_error& e) {
    cout << e.what() << endl;
  }


  // Arbitrary precision example for Li_n and Li_22
  try {
    //Create a LiSK object. Constants are prepared for
    //computations of Li_n up to n=20 with 34 digit precision
    LiSK::LiSK<cln::cl_N> lisk(20, 34);

    const cln::cl_N x = cln::complex("1/4", "1/4");
    const vector<int> weights = {1, 2, 3, 4, 5, 6, 10, 15, 20};

    // Compute Li_n(x) for n=1,2,3,4,5,6,10,15,20 at x=1/4+I*1/4
    for (auto n : weights) {
      cout << "Li(" << n << ",x) = " << setprecision(17) << lisk.Li(n, x) << endl;
    }

    // Compute Li_22(x,y) at x=y=1/4+I*1/4
    cout << "\nLi_22(x,y) = " << setprecision(17) << lisk.Li22(x, x) << endl << endl;


  }
  catch (runtime_error& e) {
    cout << e.what() << endl;
  }

  return 0;
}


