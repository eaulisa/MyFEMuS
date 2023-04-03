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
    double a;
    double b;
    double c;
};

// Function to calculate the equation of a parabola
Parabola get_parabola_equation(Point p1, Point p2, Point p3) {
    double x1 = p1.x, x2 = p2.x, x3 = p3.x;
    double y1 = p1.y, y2 = p2.y, y3 = p3.y;
    double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
    double a = (y1 * (x2 - x3) + y2 * (x3 - x1) + y3 * (x1 - x2)) / denom;
    double b = (y1 * pow(x2, 2) * x3 + y2 * pow(x3, 2) * x1 + y3 * pow(x1, 2) * x2) / denom;
    double c = (y1 * x2 * x3 + y2 * x3 * x1 + y3 * x1 * x2) / denom;
    return {a, b, c};
}

int main() {
    Point p1 = {1, 2};
    Point p2 = {3, 8};
    Point p3 = {5, 6};
    Parabola parabola = get_parabola_equation(p1, p2, p3);
    cout << "The equation of the parabola is: " << parabola.a << "x^2 + " << parabola.b << "x + " << parabola.c << endl;
    return 0;
}





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





