#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <algorithm>    // std::sort
#include <ctime>
#include <cstdlib>
#include <climits>

#include <boost/math/special_functions/factorials.hpp>
//#include <boost/math/special_functions/pow.hpp>
using namespace std;

using boost::math::factorial;

struct Point {
    double x;
    double y;
};

struct Parabola {
    double k;
    double b;
    double d;
};

// Function to calculate the equation of a parabola
Parabola get_parabola_equation(Point p1, Point p2, Point p3) {
    double x1 = p1.x, x2 = p2.x, x3 = p3.x;
    double y1 = p1.y, y2 = p2.y, y3 = p3.y;
    double det = x1 * x1 * (x2 - x3) -x1* (x2*x2 - x3*x3)+ x2*x3*(x2 - x3) ;
//     double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
    double k = (y1 * (x2 - x3) + y2 * (x3 - x1) + y3 * (x1 - x2)) / det;
    double b = (y1 * (x3*x3 - x2*x2) + y2 * (x1*x1 - x3*x3)+ y3 * ((x2*x2 - x1*x1))) / det;
    double d = (y1 * x2 * x3 * (x2 - x3) + y2 * x3 * x1 * (x3 - x1) + y3 * x1 * x2 * (x1 - x2)) / det;
    return {k, b, d};
}

int main() {
  unsigned int count = 1;
  cout << " (x1,y1) , (x2,y2) , (x3,y3) , k , b , d , c = 1" <<endl;
  Point p1,p2,p3;
  for (int table = 1; table <=8; table++){
    for (double i1=0.;i1<=1.;i1+=0.1){
      for (double i2=0.;i2<=1.;i2+=0.1){
        for (double i3=0.;i3<=1.;i3+=0.1){

      switch (table) {
        case 1:
            p1 = {0, i1};
            p2 = {i2, 1};
            break;
        case 2:
            p1 = {0, i1};
            p2 = {1,i2};
            break;
        case 3:
            p1 = {0, i1};
            p2 = {i2, 0};
            break;
        case 4:
            p1 = {i1,1};
            p2 = {i2, 1};
            if(i2 >= i1){p1 = {0, 0};p2 = {0, 0};}
            break;
        case 5:
            p1 = {i1,1};
            p2 = {1, i2};
            break;
        case 6:
            p1 = {i1,1};
            p2 = {i2, 0};
            break;
        case 7:
            p1 = {1, i1};
            p2 = {i2, 0};
            break;
        case 8:
            p1 = {i1, 0};
            p2 = {i2, 0};
            if(i2 >= i1){p1 = {0, 0};p2 = {0, 0};}
            break;
    }
//           Point p1 = {0, i1};
//           Point p2 = {i2, 1};
           p3 = {0.5 * (p1.x + p2.x) , i3 };

          double det = p1.x * p1.x * (p2.x - p3.x) -p1.x* (p2.x*p2.x - p3.x*p3.x)+ p2.x*p3.x*(p2.x - p3.x) ;
          if (det !=0){
            Parabola parabola = get_parabola_equation(p1, p2, p3);
            cout << count << ". (" << p1.x << ", " << p1.y << ") ," << "(" << p2.x << ", " << p2.y << ") ," << "(" << p3.x << ", " << p3.y << ")  : "  <<  parabola.k << ", " << parabola.b << ", " << parabola.d << endl;
            count ++ ;
          }

        }
      }
    }
  }
    return 0;
}

    //           cout << "The equation of the parabola is: " << parabola.a << "x^2 + " << parabola.b << "x + " << parabola.c << endl;
    //         cout << "The input points are: " << endl;
    //     cout << " (" << p1.x << ", " << p1.y << ") ," << "(" << p2.x << ", " << p2.y << ") ," << "(" << p3.x << ", " << p3.y << ")" << endl;


// void get_parabola_equation(float x1, float y1, float x2, float y2, float x3, float y3, float& a, float& b, float& c) {
//     // Calculate the determinant of the 3x3 matrix
//     float detA = x1 * (y2 - y3) - x2 * (y1 - y3) + x3 * (y1 - y2);
//     if (detA == 0) {
//         cout << "Error: Points are collinear!" << endl;
//         return;
//     }
//
//     // Calculate the coefficients a, b, and c of the parabola equation
//     float A = y1 * (x2 - x3) - y2 * (x1 - x3) + y3 * (x1 - x2);
//     float B = pow(x1, 2) * (y2 - y3) - pow(x2, 2) * (y1 - y3) + pow(x3, 2) * (y1 - y2);
//     float C = pow(x1, 2) * (y3 - y2) + pow(x2, 2) * (y1 - y3) + pow(x3, 2) * (y2 - y1);
//     a = A / detA;
//     b = B / detA;
//     c = C / detA;
// }
//
// int main() {
//     // Example usage of the function
//     float x1 = -0.5, y1 = 3;
//     float x2 = 0, y2 = 4;
//     float x3 = -1, y3 = 3;
//     float a, b, c;
//     get_parabola_equation(x1, y1, x2, y2, x3, y3, a, b, c);
//     cout << "The equation of the parabola passing through (" << x1 << ", " << y1 << "), (" << x2 << ", " << y2 << "), and (" << x3 << ", " << y3 << ") is:" << endl;
//     cout << "y = " << a << "x^2 + " << b << "x + " << c << endl;
//     return 0;
// }
//

/*
int main() {
clock_t t = clock();
std::srand(std::time(NULL));
for(unsigned i=0;i<100000;i++){



}




   t = clock() - t;
   std::cout << "\nTime taken for predetermined cases: " <<(double)(t)/ CLOCKS_PER_SEC << std::endl;
   return 1;
}*/



/*
    Point p1 = {0, 0.5};
    Point p2 = {0.5, 0.75};
    Point p3 = {1, 0.5};*/
//     Point p1 = { -0.5, 3};
//     Point p2 = {0, 4};
//     Point p3 = {-1, 3};
