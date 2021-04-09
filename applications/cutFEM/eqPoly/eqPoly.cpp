#include "eqPoly.hpp"
#include <algorithm>


//TODO hexahedron = 0, tet = 1, wedge =2, quad = 3, tri = 4, line = 5, point = 6
//Flag = 1 for standard calculation of coeffiecients, flag = 2 for recursion in the degenerate dim cases
void EquivalentPolynomial::SetCoefficients(const unsigned &dim, const unsigned &degree, const double &p, const std::vector < double > &c, const unsigned &element) {


  _dim = dim;
  _degree = degree;



  if(dim == 1) {
    int ni = degree;
    unsigned n = std::max(ni, 3) + 1;
    _b_vector[0].resize(n);

    if(_flag) _lisk = new LiSK::LiSK< complex<double> > (4);


    double a = c[0];
    double d = c[1];
    double x1 = exp(d * p);
    double x2 = exp((a + d) * p);
    double x3 = exp((-a + d) * p);
    double ap = a * p;
    double a2p2 = ap * ap;
    double a3p3 = a2p2 * ap;
    double a4p4 = a3p3 * ap;

    double ln1px3 = log(1 + x3);
    double ln1px2 = log(1 + x2);

    complex <double> l2mx3 = _lisk->Li(2, -x3);
    complex <double> l2mx2 = _lisk->Li(2, -x2);
    complex <double> l3mx3 = _lisk->Li(3, -x3);
    complex <double> l3mx2 = _lisk->Li(3, -x2);



    if(degree <= 3) {

      std::cout << " IN HERE DIM 1 " << endl;
      
      std::cout << a << " a " << endl;
      
      std::cout << d << " d " << endl;


      // 1
      _b_vector[0][0] = (-2.0 * (log(exp(ap) + x1) - ln1px2)) / (ap);
      // x
      _b_vector[0][1] = (2. * (ap * (ln1px3 + ln1px2) - l2mx3 + l2mx2)) / a2p2;
      // xx
      _b_vector[0][2] = -0.6666666666666666 +
                        (2. * (a2p2 * (-ln1px3 + ln1px2) +
                               2. * ap * (l2mx3 + l2mx2) + 2. * l3mx3 - 2. * l3mx2)
                        ) / a3p3;
      //xxx
      _b_vector[0][3] = (2. * (a3p3 * (ln1px3 + ln1px2) -
                               3 * ap * (ap * (l2mx3 - l2mx2) + 2. * (l3mx3 + l3mx2)) -
                               6. * _lisk->Li(4, -x3) + 6. * _lisk->Li(4, -x2))
                        ) / a4p4;


      for(unsigned i = 0; i < degree+1; i++) {
        std::cout << _b_vector[0][i] <<  "  b vector ******************************" << endl;


      }

      if(degree > 3) {

        _b_vector[0][4] = -2. + (2 * pow(a, 4) * pow(p, 4) * (-log(1 + exp((-a + d) * p)) + log(1 + exp((a + d) * p))) +
                                 8 * a * p * (a * p * (a * p * (_lisk->Li(2, -exp((-a + d) * p)) + _lisk->Li(2, -exp((a + d) * p))) +
                                                       3. * (_lisk->Li(3, -exp((-a + d) * p)) - _lisk->Li(3, -exp((a + d) * p)))) +
                                              6. * (_lisk->Li(4, -exp((-a + d) * p)) + _lisk->Li(4, -exp((a + d) * p)))) + 48. * _lisk->Li(5, -exp((-a + d) * p)) -
                                 48. * _lisk->Li(5, -exp((a + d) * p))) / (pow(a, 5) * pow(p, 5));

for(unsigned i = 0; i < degree; i++) {
        std::cout << _b_vector[0][i] <<  "  b vector ******************************" << endl;


      }
      }

      // TODO more cleaning


      std::vector < complex < double > > bcolumn(_b_vector[0].begin(), _b_vector[0].end());

      if(_flag) {


        EquivalentPolynomial::MatrixVectorMultiply(_A13_inverse, bcolumn, _coefficients);

      }

    }





    else {
      std::cout << "Degree " << degree << " has not been implemented for dimension " << dim << std::endl;
      abort();
    }
  }
  else if(dim == 2) {



    double a = c[0];
    double b = c[1];
    double d = c[2];
    unsigned n = ((degree + 1u) * (degree + 2u)) / 2u;
    _lisk = new LiSK::LiSK< complex<double> > (7);
    _b_vector[1].resize(n);
    //std::vector <double> bvector(n);
    bool aorbiszero = false;

    if(element == 3) {


      if(abs(a) < 0.00001 || abs(b) < 0.00001) {

        _flag = false;
        bool aiszero = (abs(a) < 0.001) ? true : false;
        aorbiszero = true;

        EquivalentPolynomial::SetCoefficients(dim - 1, degree, p, {(aiszero) ? b : a, d}, 5);

        if(degree >= 0) {

          //std::vector < double > dim1coefficients(degree, 0.);


          // 1
          _b_vector[1][0] = 2 * _b_vector[0][0].real();
          // x
          _b_vector[1][(aiszero) ? 1 : 2] = 0.;
          // y
          _b_vector[1][(aiszero) ? 2 : 1] = 2 * _b_vector[0][1].real();
          // xx
          _b_vector[1][(aiszero) ? 3 : 4] = (2. / 3.) * _b_vector[0][0].real();
          // yy
          _b_vector[1][(aiszero) ? 4 : 3] = 2 * _b_vector[0][2].real();
          // xy
          _b_vector[1][5] = 0. ;

        }


        if(degree > 2) {
          // xxx
          _b_vector[1][(aiszero) ? 6 : 7] = 0.;
          // yyy
          _b_vector[1][(aiszero) ? 7 : 6] = 2 * _b_vector[0][3].real();
          // xx y
          _b_vector[1][(aiszero) ? 8 : 9] = (2. / 3.) * _b_vector[0][1].real();
          // yy x
          _b_vector[1][(aiszero) ? 9 : 8] = 0.;

        }


        if(degree > 3) {
          // xxxx
          _b_vector[1][(aiszero) ? 10 : 11] = (2. / 5.) * _b_vector[0][0].real();
          // yyyy
          _b_vector[1][(aiszero) ? 11 : 10] = 2 * _b_vector[0][4].real();
          // xx yy
          _b_vector[1][12] = (2. / 3.) * _b_vector[0][2].real();
          // xxx y
          _b_vector[1][13] = 0.;
          // yyy x
          _b_vector[1][14] = 0.;



        }

        else std::cout << "Degree " << degree << " has not been implemented for dimension " << dim << std::endl;


        std::vector < complex < double > > bcolumn(_b_vector[1].begin(), _b_vector[1].end());


        EquivalentPolynomial::MatrixVectorMultiply(_A22_inverse, bcolumn, _coefficients);

      }








      double a_squared = a * a;
      double p_squared = p * p;
      double x2 = exp((-a - b + d) * p);
      double x3 = exp((-a + b + d) * p);
      double x4 = exp((a - b + d) * p);
      double x5 = exp((a + b + d) * p);
      double b2 = b * b;
      double b3 = b * b2;
      double p4 = p * p * p * p;
      double ap = a * p;
      double ab = a * b;
      double bp = b * p;

      complex <double> l2mx2 = _lisk->Li(2, -x2);
      complex <double> l2mx3 = _lisk->Li(2, -x3);
      complex <double> l2mx4 = _lisk->Li(2, -x4);
      complex <double> l2mx5 = _lisk->Li(2, -x5);
      complex <double> l3mx2 = _lisk->Li(3, -x2);
      complex <double> l3mx4 = _lisk->Li(3, -x4);
      complex <double> l3mx3 = _lisk->Li(3, -x3);
      complex <double> l3mx5 = _lisk->Li(3, -x5);
      complex <double> l4mx2 = _lisk->Li(4, -x2);
      complex <double> l4mx4 = _lisk->Li(4, -x4);
      complex <double> l4mx3 = _lisk->Li(4, -x3);
      complex <double> l4mx5 = _lisk->Li(4, -x5);


      if(degree >= 0 && !aorbiszero) {




        // 1
        _b_vector[1][0] = (-2. * (2 * ab * p_squared + l2mx2 -
                                  l2mx4 - l2mx3 +
                                  l2mx5)) / (a * b * p_squared);
        // x
        _b_vector[1][1] = (2. * (ap * l2mx2 + ap * l2mx4 - ap * l2mx3 -
                                 ap * l2mx5 + l3mx2 - l3mx4 -
                                 l3mx3 + l3mx5)) / (a_squared * b * p_squared * p);
        // y
        _b_vector[1][2] = (2. * (b * p * l2mx2 - bp * l2mx4 + bp * l2mx3 -
                                 bp * l2mx5 + l3mx2 - l3mx4 -
                                 l3mx3 + l3mx5)) / (a * b2 * p_squared * p);
        // xx
        _b_vector[1][3] = (-2. * (2 * a_squared * ab * p4 + 3 * a_squared * p_squared * l2mx2 -
                                  3 * a_squared * p_squared * l2mx4 - 3 * a_squared * p_squared * l2mx3 +
                                  3 * a_squared * p_squared * l2mx5 + 6 * ap * l3mx2 +
                                  6 * ap * l3mx4 - 6 * ap * l3mx3 - 6 * ap * l3mx5 +
                                  6. * l4mx2 - 6. * l4mx4 - 6. * l4mx3 +
                                  6. * l4mx5)) / (3. * a_squared * ab * p4);
        // yy
        _b_vector[1][4] = (-2. * (2 * a * b3 * p4 + 3 * b2 * p_squared * l2mx2 -
                                  3 * b2 * p_squared * l2mx4 - 3 * b2 * p_squared * l2mx3 +
                                  3 * b2 * p_squared * l2mx5 + 6 * bp * l3mx2 -
                                  6 * bp * l3mx4 + 6 * bp * l3mx3 - 6 * bp * l3mx5 +
                                  6. * l4mx2 - 6. * l4mx4 - 6. * l4mx3 +
                                  6. * l4mx5)) / (3. * a * b3 * p4);
        // xy
        _b_vector[1][5] = (-2. * (a * b * p_squared * l2mx2 + ab * p_squared * l2mx4 +
                                  ab * p_squared * l2mx3 + ab * p_squared * l2mx5 +
                                  ap * l3mx2 + bp * l3mx2 + ap * l3mx4 -
                                  bp * l3mx4 - ap * l3mx3 + bp * l3mx3 -
                                  ap * l3mx5 - bp * l3mx5 + l4mx2 -
                                  l4mx4 - l4mx3 + l4mx5)) /
                          (a_squared * b2 * p4);



        for(unsigned i = 0; i < n; i++) {
          std::cout << _b_vector[1][i] <<  "  b vector ******************************" << endl;


        }


        if(degree > 2) {

          double x6 = exp(-((a + b - d) * p));
          double x7 = exp(-((-a + b + d) * p));
          double x8 = exp((b - d) * p);
          double x9 = exp((a + b - d) * p);
          double x10 = exp((-b + d) * p);
          double x11 = exp((b + d) * p);
          double x13 = exp((-a - b - d) * p);


          double bp = b * p;
          double a2p2 = ap * ap;
          double a3p3 = a2p2 * ap;
          double a4p4 = a3p3 * ap;
          double b2p2 = bp * bp;
          double b3p3 = b2p2 * bp;
          double b4p4 = b3p3 * bp;
          double a2 = a * a;
          double a4 = a2 * a2;
          double a3 = a * a2;
          double b4 = b2 * b2;
          double p5 = p4 * p;



          complex <double> l2mx6 = _lisk->Li(2, -x6);
          complex <double> l5mx4 = _lisk->Li(5, -x4);
          complex <double> l3mx6 = _lisk->Li(3, -x6);
          complex <double> l4mx6 = _lisk->Li(4, -x6);
          complex <double> l5mx6 = _lisk->Li(5, -x6);
          complex <double> l5mx3 = _lisk->Li(5, -x3);
          complex <double> l5mx5 = _lisk->Li(5, -x5);
          complex <double> l2mx9 = _lisk->Li(2, -x9);
          complex <double> l2mx7 = _lisk->Li(2, -x7);
          complex <double> l3mx9 = _lisk->Li(3, -x9);
          complex <double> l3mx7 = _lisk->Li(3, -x7);
          complex <double> l4mx9 = _lisk->Li(4, -x9);
          complex <double> l4mx7 = _lisk->Li(4, -x7);
          complex <double> l6mx6 = _lisk->Li(6, -x6);
          complex <double> l6mx4 = _lisk->Li(6, -x4);
          complex <double> l6mx3 = _lisk->Li(6, -x3);
          complex <double> l6mx5 = _lisk->Li(6, -x5);

          double ln1px9 = log(1 + x9);
          double ln1px4 = log(1 + x4);
          double ln1px7 = log(1 + x7);
          double ln1px5 = log(1 + x5);


          // xxx
          _b_vector[1][6] = (2. * (a3p3 * l2mx6 + a3p3 * l2mx4 -
                                   a3p3 * l2mx3 - a3p3 * l2mx5 +
                                   3 * a2p2 * l3mx6 - 3 * a2p2 * l3mx4 -
                                   3 * a2p2 * l3mx3 + 3 * a2p2 * l3mx5 +
                                   6. * ap * l4mx6 + 6. * ap * l4mx4 -
                                   6. * ap * l4mx3 - 6. * ap * l4mx5 +
                                   6. * l5mx6 - 6. * l5mx4 -
                                   6. * l5mx3 + 6. * l5mx5)) / (a4 * b * p5);
          // yyy
          _b_vector[1][7] = (2. * (b3p3 * l2mx6 -
                                   b3p3 * l2mx4 +
                                   b3p3 * l2mx3 -
                                   b3p3 * l2mx5 +
                                   3 * b2p2 * l3mx6 -
                                   3 * b2p2 * l3mx4 -
                                   3 * b2p2 * l3mx3 +
                                   3 * b2p2 * l3mx5 +
                                   6. * bp * l4mx6 -
                                   6. * bp * l4mx4 + 6. * bp * l4mx3 -
                                   6. * bp * l4mx5 + 6. * l5mx6 -
                                   6. * l5mx4 - 6. * l5mx3 +
                                   6. * l5mx5)) / (a * b4 * p5);

          // xxy
          _b_vector[1][8] = (6. * ap * (ap * (l3mx6 -
                                              l3mx4) + 2. * (l4mx6 +
                                                             l4mx4)) +
                             bp * (a4p4 - 2 * a3p3 * ln1px9 -
                                   2 * a3p3 * ln1px4 -
                                   2 * a3p3 * ln1px7 +
                                   2 * a3p3 * log(((1 + x9) * (1 + x3)) / x9) -
                                   2 * a3p3 * ln1px5 +
                                   2 * a3p3 * log((1 + x4) * (1 + x5)) -
                                   6 * a2p2 * l2mx9 -
                                   6 * a2p2 * l2mx4 -
                                   6 * a2p2 * l2mx7 -
                                   6 * a2p2 * l2mx5 +
                                   12. * ap * l3mx9 +
                                   12. * ap * l3mx4 +
                                   12. * ap * l3mx7 +
                                   12. * ap * l3mx5 +
                                   12. * _lisk->Li(4, -x8) -
                                   12. * l4mx9 +
                                   12. * _lisk->Li(4, -x10) -
                                   12. * l4mx4 +
                                   12. * _lisk->Li(4, -x13) +
                                   12. * _lisk->Li(4, -x11) -
                                   12. * l4mx7 -
                                   12. * l4mx5) -
                             6 * ap * (ap * (l3mx3 -
                                             l3mx5) + 2. * (l4mx3 +
                                                            l4mx5)) +
                             12. * (l5mx6 - l5mx4) -
                             12. * (l5mx3 - l5mx5)) /
                            (3. * a3 * b2 * p5);

          // yyx
          _b_vector[1][9] = (2. * (b2p2 * (ap * l2mx6 +
                                           ap * l2mx4 + l3mx6 -
                                           l3mx4) + b2p2 * (-(ap * (l2mx3 +
                                                                    l2mx5)) - l3mx3 +
                                                            l3mx5) + 2 * bp * (ap * (l3mx6 +
                                                                                     l3mx4) + l4mx6 -
                                                                               l4mx4) + 2 * bp * (ap * (l3mx3 +
                                                                                                        l3mx5) + l4mx3 -
                                                                                                  l4mx5) + 2. * (ap * (l4mx6 +
                                                                                                      l4mx4) + l5mx6 -
                                                                                                      l5mx4) - 2. * (ap * (l4mx3 +
                                                                                                          l4mx5) + l5mx3 -
                                                                                                          l5mx5))) / (a2 * b3 * p5);

          if(degree > 3) {

            double a5 = a * a4;
            double p6 = p * p5;
            double b5 = b * b4;


            // xxxx
            _b_vector[1][10] = (-2. * (2 * a5 * b * p6 +
                                       5 * a4p4 * l2mx6 -
                                       5 * a4p4 * l2mx4 -
                                       5 * a4p4 * l2mx3 +
                                       5 * a4p4 * l2mx5 +
                                       20. * (a3p3 * l3mx6 +
                                              a3p3 * l3mx4 +
                                              3 * a2p2 * l4mx6 -
                                              3 * a2p2 * l4mx4 +
                                              6. * ap * l5mx6 +
                                              6. * ap * l5mx4 +
                                              6. * l6mx6 -
                                              6. * l6mx4) -
                                       20. * (a3p3 * l3mx3 +
                                              a3p3 * l3mx5 +
                                              3 * a2p2 * l4mx3 -
                                              3 * a2p2 * l4mx5 +
                                              6. * ap * l5mx3 +
                                              6. * ap * l5mx5 +
                                              6. * l6mx3 -
                                              6. * l6mx5))) / (5. * a5 * b * p6);

            // yyyy
            _b_vector[1][11] = (-2. * (2 * a * b5 * p6 +
                                       5 * b4p4 * l2mx6 -
                                       5 * b4p4 * l2mx4 -
                                       5 * b4p4 * l2mx3 +
                                       5 * b4p4 * l2mx5 +
                                       20 * b3p3 * l3mx6 -
                                       20 * b3p3 * l3mx4 +
                                       20 * b3p3 * l3mx3 -
                                       20 * b3p3 * l3mx5 +
                                       60 * b2p2 * l4mx6 -
                                       60 * b2p2 * l4mx4 -
                                       60 * b2p2 * l4mx3 +
                                       60 * b2p2 * l4mx5 +
                                       120 * bp * l5mx6 -
                                       120 * bp * l5mx4 +
                                       120 * bp * l5mx3 -
                                       120 * bp * l5mx5 +
                                       120. * l6mx6 -
                                       120. * l6mx4 -
                                       120. * l6mx3 +
                                       120. * l6mx5)) / (5. * a * b5 * p6);

            // xxyy
            _b_vector[1][12] = (-2. * (2. * a3 * b3 * p6 +
                                       9 * b2p2 * (a2p2 * l2mx6 -
                                                   a2p2 * l2mx4 + 2. * (ap * l3mx6 +
                                                                        ap * l3mx4 + l4mx6 - l4mx4)) - 9 * b2p2 * (a2p2 * l2mx3 -
                                                                            a2p2 * l2mx5 + 2. * (ap * l3mx3 + ap * l3mx5 +
                                                                                                 l4mx3 - l4mx5)) + 18. * bp * (ap * (ap * (l3mx6 -
                                                                                                     l3mx4) + 2. * (l4mx6 + l4mx4)) + 2. * (l5mx6 - l5mx4)) +
                                       18. * bp * (ap * (ap * (l3mx3 - l3mx5) + 2. * (l4mx3 +
                                                                                      l4mx5)) + 2. * (l5mx3 - l5mx5)) + 18. * (ap * (ap * (l4mx6 -
                                                                                          l4mx4) + 2. * (l5mx6 + l5mx4)) +
                                                                                          2. * (l6mx6 - l6mx4)) - 18. * (ap * (ap * (l4mx3 -
                                                                                              l4mx5) + 2. * (l5mx3 + l5mx5)) +
                                                                                              2. * l6mx3 - 2. * l6mx5))) / (9. * a3 * b3 * p6);

            // xxxy
            _b_vector[1][13] = (b * p * (-2 * a5 * p5 +
                                         5 * a4p4 * ln1px9 -
                                         5 * a4p4 * ln1px4 +
                                         5 * a4p4 * ln1px7 -
                                         5 * a4p4 * log(((1 + x9) * (1 + x3)) / x9) -
                                         5 * a4p4 * ln1px5 +
                                         5 * a4p4 * log((1 + x4) * (1 + x5)) +
                                         20 * a3p3 * l2mx9 -
                                         20 * a3p3 * l2mx4 +
                                         20 * a3p3 * l2mx7 -
                                         20 * a3p3 * l2mx5 -
                                         60 * a2p2 * l3mx9 +
                                         60 * a2p2 * l3mx4 -
                                         60 * a2p2 * l3mx7 +
                                         60 * a2p2 * l3mx5 +
                                         120. * ap * l4mx9 -
                                         120. * ap * l4mx4 +
                                         120. * ap * l4mx7 -
                                         120. * ap * l4mx5 +
                                         120. * _lisk->Li(5, -x8) -
                                         120. * _lisk->Li(5, -x9) - 120. * _lisk->Li(5, -x10) +
                                         120. * l5mx4 + 120. * _lisk->Li(5, -x13) -
                                         120. * _lisk->Li(5, -x11) - 120. * _lisk->Li(5, -x7) +
                                         120. * l5mx5) -
                                20. * (a3p3 * l3mx6 +
                                       a3p3 * l3mx4 +
                                       3 * a2p2 * l4mx6 -
                                       3 * a2p2 * l4mx4 +
                                       6. * ap * l5mx6 +
                                       6. * ap * l5mx4 +
                                       6. * l6mx6 -
                                       6. * l6mx4) +
                                20. * (a3p3 * l3mx3 +
                                       a3p3 * l3mx5 +
                                       3 * a2p2 * l4mx3 -
                                       3 * a2p2 * l4mx5 +
                                       6. * ap * l5mx3 +
                                       6. * ap * l5mx5 + 6. * l6mx3 -
                                       6. * l6mx5)) / (10. * a4 * b2 * p6);

            // yyyx
            _b_vector[1][14] = (2. * (b3p3 * (-(ap * (l2mx6 + l2mx4)) - l3mx6 + l3mx4) +
                                      b3p3 * (-(ap * (l2mx3 + l2mx5)) - l3mx3 + l3mx5) -
                                      3 * b2p2 * (ap * (l3mx6 + l3mx4) + l4mx6 - l4mx4) +
                                      3 * b2p2 * (ap * (l3mx3 + l3mx5) + l4mx3 - l4mx5) -
                                      6 * bp * (ap * (l4mx6 + l4mx4) + l5mx6 - l5mx4) -
                                      6 * bp * (ap * (l4mx3 + l4mx5) + l5mx3 - l5mx5) -
                                      6. * (ap * (l5mx6 + l5mx4) + l6mx6 - l6mx4) +
                                      6. * (ap * (l5mx3 + l5mx5) + l6mx3 - l6mx5))) /
                               (a2 * b4 * p6);

            if(degree > 4) {
              //std::cout << "hope not here   " << degree << "  " << dim << std::endl;
              std::cout << "Degree " << degree << " has not been implemented for dimension " << dim << std::endl;
              abort();
            }
          }
        }
      }

      //std::vector < complex < double > > bcolumn(_b_vector[1].begin(), _b_vector[1].end());

      //std::cout << " _flag =  " << _flag   << endl;



      if(degree <= 2 && _flag) {
        //std::cout << _b_vector[0] << "   this is right before MatrixVectorMultiply" << std::endl;
        EquivalentPolynomial::MatrixVectorMultiply(_A22_inverse, _b_vector[1], _coefficients);

      }
      else if(degree == 3 && _flag) EquivalentPolynomial::MatrixVectorMultiply(_A23_inverse, _b_vector[1], _coefficients);
      else if(degree == 4 && _flag) EquivalentPolynomial::MatrixVectorMultiply(_A24_inverse, _b_vector[1], _coefficients);

    }


    else if(element == 4) {

      std::cout << " Triangle "   << endl;

      _b_vector[1][0] = -0.5 + (2. * ((-a + b) * _lisk->Li(2, -exp(d * p)) - b * _lisk->Li(2, -exp((a + d) * p)) + a * _lisk->Li(2, -exp((b + d) * p)))) / (a * (a - b) * b * pow(p, 2));

      _b_vector[1][1] = -0.16666666666666666 + (2. * (a * b * (-a + b) * p * _lisk->Li(2, -exp((a + d) * p)) +
                                                      pow(a - b, 2) * _lisk->Li(3, -exp(d * p)) + (2 * a - b) * b * _lisk->Li(3, -exp((a + d) * p))
                                                      - pow(a, 2) * _lisk->Li(3, -exp((b + d) * p)))) / (pow(a, 2) * pow(a - b, 2) * b * pow(p, 3));

      _b_vector[1][2] = -(pow(a, 3) * pow(b, 2) * pow(p, 3) - 2 * pow(a, 2) * pow(b, 3) * pow(p, 3) +
                          a * pow(b, 4) * pow(p, 3) - 12 * a * (a - b) * b * p * _lisk->Li(2, -exp((b + d) * p)) -
                          12 * pow(a - b, 2) * _lisk->Li(3, -exp(d * p)) + 12 * pow(b, 2) * _lisk->Li(3, -exp((a + d) * p)) +
                          12 * pow(a, 2) * _lisk->Li(3, -exp((b + d) * p)) - 24 * a * b * _lisk->Li(3, -exp((b + d) * p))) / (6. * a * pow(a - b, 2) * pow(b, 2) * pow(p, 3));

      _b_vector[1][3] = -0.08333333333333333 + (2. * (pow(a, 2) * pow(p, 2) * _lisk->Li(2, -exp((a + d) * p)) -
                                                      2. * (a * p * _lisk->Li(3, -exp((a + d) * p)) + _lisk->Li(4, -exp(d * p)) - _lisk->Li(4, -exp((a + d) * p))))) / (pow(a, 3) * b * pow(p, 4)) +
                        (2. * ((a - b) * p * ((-a + b) * p * _lisk->Li(2, -exp((a + d) * p)) + 2. * _lisk->Li(3, -exp((a + d) * p))) - 2. * _lisk->Li(4, -exp((a + d) * p)) + 2. * _lisk->Li(4, -exp((b + d) * p)))) /
                        (pow(a - b, 3) * b * pow(p, 4));

      _b_vector[1][4] = -0.08333333333333333 + (2 * a * pow(a - b, 2) * pow(b, 2) * pow(p, 2) * _lisk->Li(2, -exp((b + d) * p)) - 4 * a * (a - 2 * b) * (a - b) * b * p * _lisk->Li(3,                                 -exp((b + d) * p)) - 4 * pow(a - b, 3) * _lisk->Li(4, -exp(d * p)) - 4 * pow(b, 3) * _lisk->Li(4, -exp((a + d) * p)) +
                                                4 * a * (pow(a, 2) - 3 * a * b + 3 * pow(b, 2)) * _lisk->Li(4, -exp((b + d) * p))) / (a * pow(a - b, 3) * pow(b, 3) * pow(p, 4));

      _b_vector[1][5] = -(pow(a, 5) * pow(b, 2) * pow(p, 4) - 3 * pow(a, 4) * pow(b, 3) * pow(p, 4) + 3 * pow(a, 3) * pow(b, 4) * pow(p, 4) - pow(a, 2) * pow(b, 5) * pow(p, 4) +
                          48 * a * (a - b) * pow(b, 2) * p * _lisk->Li(3, -exp((a + d) * p)) + 48 * pow(a, 2) * (a - b) * b * p * _lisk->Li(3, -exp((b + d) * p)) + 48 * pow(a, 3) * _lisk->Li(4, -exp(d * p)) - 144 * pow(a, 2) * b * _lisk->Li(4, -exp(d * p)) + 144 * a * pow(b, 2) * _lisk->Li(4, -exp(d * p)) - 48 * pow(b, 3) * _lisk->Li(4, -exp(d * p)) - 144 * a * pow(b, 2) * _lisk->Li(4, -exp((a + d) * p)) + 48 * pow(b, 3) * _lisk->Li(4, -exp((a + d) * p)) -
                          48 * pow(a, 3) * _lisk->Li(4, -exp((b + d) * p)) +
                          144 * pow(a, 2) * b * _lisk->Li(4, -exp((b + d) * p))) / (24. * pow(a, 2) * pow(a - b, 3) * pow(b, 2) * pow(p, 4));



      std::vector < complex < double > > bcolumn(_b_vector[1].begin(), _b_vector[1].end());

      if(degree <= 2) EquivalentPolynomial::MatrixVectorMultiply(_A22T_inverse, bcolumn, _coefficients);
      //else if(degree == 3) EquivalentPolynomial::MatrixVectorMultiply(_A23_inverse, bcolumn, _coefficients);
      //else if(degree == 4) EquivalentPolynomial::MatrixVectorMultiply(_A24_inverse, bcolumn, _coefficients);

    }

  }



  else if(dim == 3) {


    double a = c[0];
    double b = c[1];
    double cz = c[2];
    double d = c[3];
    unsigned n = ((degree + 1u) * (degree + 2u) * (degree + 3u)) / 6u;
    _lisk = new LiSK::LiSK< complex<double> > (7);
    _b_vector[2].resize(n);
    std::vector <double> bvector(n);


    if(element == 3) {


      if(abs(a) < 0.0001) {

        _flag = false;


        if(abs(b) < 0.00001) {

          EquivalentPolynomial::SetCoefficients(dim - 2, degree, p, {cz, d}, 5);


          if(degree >= 2) {


            std::cout << " IN HERE right"   << endl;



            _b_vector[2][0] = 4. * _b_vector[0][0];

            _b_vector[2][1] = 0.;

            _b_vector[2][2] = 0.;

            _b_vector[2][3] = 4. * _b_vector[0][1];

            _b_vector[2][4] = (2. / 3.) * _b_vector[0][0];

            _b_vector[2][5] = (2. / 3.) * _b_vector[0][0];

            _b_vector[2][6] = 4. * _b_vector[0][2];

            _b_vector[2][7] = 0.;

            _b_vector[2][8] = 0.;

            _b_vector[2][9] = 0.;

            for(unsigned i = 0; i < n; i++) {
              std::cout << _b_vector[2][i] <<  "  b vector ******************************" << endl;


            }


            if(degree <= 2) EquivalentPolynomial::MatrixVectorMultiply(_A32_inverse, _b_vector[2], _coefficients);

          }



          if(degree > 2) {


          }

        }

        else if(abs(cz) < 0.00001) {

          EquivalentPolynomial::SetCoefficients(dim - 2, degree, p, {b, d}, 5);

          if(degree >= 2) {



            _b_vector[2][0] = 4. * _b_vector[0][0];

            _b_vector[2][1] = 0.;

            _b_vector[2][2] = 4. * _b_vector[0][1];

            _b_vector[2][3] = 0.;

            _b_vector[2][4] = (2. / 3.) * _b_vector[0][0];

            _b_vector[2][5] = 4. * _b_vector[0][2];

            _b_vector[2][6] = (2. / 3.) * _b_vector[0][0];

            _b_vector[2][7] = 0.;

            _b_vector[2][8] = 0.;

            _b_vector[2][9] = 0.;

            if(degree <= 2) EquivalentPolynomial::MatrixVectorMultiply(_A32_inverse, _b_vector[2], _coefficients);



          }



          if(degree > 2) {


          }




        }

        else {

          EquivalentPolynomial::SetCoefficients(dim - 1, degree, p, {b, cz, d}, 3);


          if(degree >= 2) {

            std::cout << " IN HERE wrong"   << endl;


            _b_vector[2][0] = 2. * _b_vector[1][0];

            _b_vector[2][1] = 0.;

            _b_vector[2][2] = 2. * _b_vector[1][1];

            _b_vector[2][3] = 2. * _b_vector[1][2];

            _b_vector[2][4] = (2. / 3.) * _b_vector[1][0];

            _b_vector[2][5] = 2. * _b_vector[1][3];

            _b_vector[2][6] = 2. * _b_vector[1][4];

            _b_vector[2][7] = 0.;

            _b_vector[2][8] = 0.;

            _b_vector[2][9] = 2. * _b_vector[1][5];

            if(degree <= 2) EquivalentPolynomial::MatrixVectorMultiply(_A32_inverse, _b_vector[2], _coefficients);


          }
        }
      }

      else if(abs(b) < 0.0001) {

        _flag = false;

        std::cout << " IN HERE wrong b = 0 spot"   << endl;


        if(abs(cz) < 0.00001) {

          EquivalentPolynomial::SetCoefficients(dim - 2, degree, p, {a, d}, 5);

          if(degree >= 2) {




            _b_vector[2][0] = 4. * _b_vector[0][0];

            _b_vector[2][1] = 4. * _b_vector[0][1];

            _b_vector[2][2] = 0.;

            _b_vector[2][3] = 0.;

            _b_vector[2][4] = 4. * _b_vector[0][2];

            _b_vector[2][5] = (2. / 3.) * _b_vector[0][0];

            _b_vector[2][6] = (2. / 3.) * _b_vector[0][0];

            _b_vector[2][7] = 0.;

            _b_vector[2][8] = 0.;

            _b_vector[2][9] = 0.;


            if(degree <= 2) EquivalentPolynomial::MatrixVectorMultiply(_A32_inverse, _b_vector[2], _coefficients);


          }



          if(degree > 2) {


          }




        }


        else {

          EquivalentPolynomial::SetCoefficients(dim - 1, degree, p, {a, cz, d}, 3);


          if(degree >= 2) {


            _b_vector[2][0] = 2. * _b_vector[1][0];

            _b_vector[2][1] = 2. * _b_vector[1][1];

            _b_vector[2][2] = 0.;

            _b_vector[2][3] = 2. * _b_vector[1][1];

            _b_vector[2][4] = 2. * _b_vector[1][2];

            _b_vector[2][5] = (2. / 3.) * _b_vector[1][0];

            _b_vector[2][6] = 2. * _b_vector[1][4];

            _b_vector[2][7] = 0.;

            _b_vector[2][8] = 2. * _b_vector[1][5];

            _b_vector[2][9] = 0.;

            if(degree <= 2) EquivalentPolynomial::MatrixVectorMultiply(_A32_inverse, _b_vector[2], _coefficients);

          }
        }
      }


      else if(abs(cz) < 0.0001) {

        _flag = false;

        EquivalentPolynomial::SetCoefficients(dim - 1, degree, p, {a, b, d}, 3);


        if(degree >= 2) {


          _b_vector[2][0] = 2. * _b_vector[1][0];

          _b_vector[2][1] = 2. * _b_vector[1][1];

          _b_vector[2][2] = 2. * _b_vector[1][2];

          _b_vector[2][3] = 0.;

          _b_vector[2][4] = 2. * _b_vector[1][3];

          _b_vector[2][5] = 2. * _b_vector[1][4];

          _b_vector[2][6] = (2. / 3.) * _b_vector[1][0];

          _b_vector[2][7] = 2. * _b_vector[1][5];

          _b_vector[2][8] = 0.;

          _b_vector[2][9] = 0.;

          if(degree <= 2) EquivalentPolynomial::MatrixVectorMultiply(_A32_inverse, _b_vector[2], _coefficients);


        }
      }




      if(abs(a) < 0.0001 && abs(b) < 0.0001 && _flag) {

        _flag = false;






      }


      if(abs(a) < 0.000001 && abs(cz) < 0.000001 && _flag) {

        _flag = false;


      }


      if(abs(b) < 0.000001 && abs(cz) < 0.000001 && _flag) {

        _flag = false;


      }

      if(_flag) {

        std::cout << " IN HERE really bad"   << endl;


        if(degree >= 2) {


          // 1
          _b_vector[2][0] = (-2. * (4 * a * b * cz * pow(p, 3) - _lisk->Li(3, -exp((-a - b - cz + d) * p)) + _lisk->Li(3, -exp((a - b - cz + d) * p)) + _lisk->Li(3, -exp((-a + b - cz + d) * p)) - _lisk->Li(3, -exp((a + b - cz + d) * p)) +
                                    _lisk->Li(3, -exp((-a - b + cz + d) * p)) - _lisk->Li(3, -exp((a - b + cz + d) * p)) - _lisk->Li(3, -exp((-a + b + cz + d) * p)) + _lisk->Li(3, -exp((a + b + cz + d) * p)))) / (a * b * cz * pow(p, 3));
          // x
          _b_vector[2][1] = (-2. * (a * p * _lisk->Li(3, -exp((-a - b - cz + d) * p)) + a * p * _lisk->Li(3, -exp((a - b - cz + d) * p)) - a * p * _lisk->Li(3, -exp((-a + b - cz + d) * p)) - a * p * _lisk->Li(3, -exp((a + b - cz + d) * p)) -
                                    a * p * _lisk->Li(3, -exp((-a - b + cz + d) * p)) - a * p * _lisk->Li(3, -exp((a - b + cz + d) * p)) + a * p * _lisk->Li(3, -exp((-a + b + cz + d) * p)) + a * p * _lisk->Li(3, -exp((a + b + cz + d) * p)) +
                                    _lisk->Li(4, -exp((-a - b - cz + d) * p)) - _lisk->Li(4, -exp((a - b - cz + d) * p)) - _lisk->Li(4, -exp((-a + b - cz + d) * p)) + _lisk->Li(4, -exp((a + b - cz + d) * p)) -
                                    _lisk->Li(4, -exp((-a - b + cz + d) * p)) + _lisk->Li(4, -exp((a - b + cz + d) * p)) + _lisk->Li(4, -exp((-a + b + cz + d) * p)) - _lisk->Li(4, -exp((a + b + cz + d) * p)))) /
                            (pow(a, 2) * b * cz * pow(p, 4));
          // y
          _b_vector[2][2] = (-2. * (b * p * _lisk->Li(3, -exp((-a - b - cz + d) * p)) - b * p * _lisk->Li(3, -exp((a - b - cz + d) * p)) + b * p * _lisk->Li(3, -exp((-a + b - cz + d) * p)) -
                                    b * p * _lisk->Li(3, -exp((a + b - cz + d) * p)) - b * p * _lisk->Li(3, -exp((-a - b + cz + d) * p)) + b * p * _lisk->Li(3, -exp((a - b + cz + d) * p)) - b * p * _lisk->Li(3, -exp((-a + b + cz + d) * p)) +
                                    b * p * _lisk->Li(3, -exp((a + b + cz + d) * p)) + _lisk->Li(4, -exp((-a - b - cz + d) * p)) - _lisk->Li(4, -exp((a - b - cz + d) * p)) - _lisk->Li(4, -exp((-a + b - cz + d) * p)) +
                                    _lisk->Li(4, -exp((a + b - cz + d) * p)) - _lisk->Li(4, -exp((-a - b + cz + d) * p)) + _lisk->Li(4, -exp((a - b + cz + d) * p)) + _lisk->Li(4, -exp((-a + b + cz + d) * p)) -
                                    _lisk->Li(4, -exp((a + b + cz + d) * p)))) / (a * pow(b, 2) * cz * pow(p, 4));
          // z
          _b_vector[2][3] = (-2. * (cz * p * _lisk->Li(3, -exp((-a - b - cz + d) * p)) - cz * p * _lisk->Li(3, -exp((a - b - cz + d) * p)) - cz * p * _lisk->Li(3, -exp((-a + b - cz + d) * p)) + cz * p * _lisk->Li(3, -exp((a + b - cz + d) * p)) +
                                    cz * p * _lisk->Li(3, -exp((-a - b + cz + d) * p)) - cz * p * _lisk->Li(3, -exp((a - b + cz + d) * p)) - cz * p * _lisk->Li(3, -exp((-a + b + cz + d) * p)) + cz * p * _lisk->Li(3, -exp((a + b + cz + d) * p)) +
                                    _lisk->Li(4, -exp((-a - b - cz + d) * p)) - _lisk->Li(4, -exp((a - b - cz + d) * p)) - _lisk->Li(4, -exp((-a + b - cz + d) * p)) + _lisk->Li(4, -exp((a + b - cz + d) * p)) -
                                    _lisk->Li(4, -exp((-a - b + cz + d) * p)) + _lisk->Li(4, -exp((a - b + cz + d) * p)) + _lisk->Li(4, -exp((-a + b + cz + d) * p)) - _lisk->Li(4, -exp((a + b + cz + d) * p)))) /
                            (a * b * pow(cz, 2) * pow(p, 4));
          // xy
          _b_vector[2][4] = (2. * (a * b * pow(p, 2) * _lisk->Li(3, -exp((-a - b - cz + d) * p)) + a * b * pow(p, 2) * _lisk->Li(3, -exp((a - b - cz + d) * p)) +
                                   a * b * pow(p, 2) * _lisk->Li(3, -exp((-a + b - cz + d) * p)) + a * b * pow(p, 2) * _lisk->Li(3, -exp((a + b - cz + d) * p)) - a * b * pow(p, 2) * _lisk->Li(3, -exp((-a - b + cz + d) * p)) -
                                   a * b * pow(p, 2) * _lisk->Li(3, -exp((a - b + cz + d) * p)) - a * b * pow(p, 2) * _lisk->Li(3, -exp((-a + b + cz + d) * p)) - a * b * pow(p, 2) * _lisk->Li(3, -exp((a + b + cz + d) * p)) +
                                   a * p * _lisk->Li(4, -exp((-a - b - cz + d) * p)) + b * p * _lisk->Li(4, -exp((-a - b - cz + d) * p)) + a * p * _lisk->Li(4, -exp((a - b - cz + d) * p)) - b * p * _lisk->Li(4, -exp((a - b - cz + d) * p)) -
                                   a * p * _lisk->Li(4, -exp((-a + b - cz + d) * p)) + b * p * _lisk->Li(4, -exp((-a + b - cz + d) * p)) - a * p * _lisk->Li(4, -exp((a + b - cz + d) * p)) - b * p * _lisk->Li(4, -exp((a + b - cz + d) * p)) -
                                   a * p * _lisk->Li(4, -exp((-a - b + cz + d) * p)) - b * p * _lisk->Li(4, -exp((-a - b + cz + d) * p)) - a * p * _lisk->Li(4, -exp((a - b + cz + d) * p)) + b * p * _lisk->Li(4, -exp((a - b + cz + d) * p)) +
                                   a * p * _lisk->Li(4, -exp((-a + b + cz + d) * p)) - b * p * _lisk->Li(4, -exp((-a + b + cz + d) * p)) + a * p * _lisk->Li(4, -exp((a + b + cz + d) * p)) + b * p * _lisk->Li(4, -exp((a + b + cz + d) * p)) +
                                   _lisk->Li(5, -exp((-a - b - cz + d) * p)) - _lisk->Li(5, -exp((a - b - cz + d) * p)) - _lisk->Li(5, -exp((-a + b - cz + d) * p)) + _lisk->Li(5, -exp((a + b - cz + d) * p)) -
                                   _lisk->Li(5, -exp((-a - b + cz + d) * p)) + _lisk->Li(5, -exp((a - b + cz + d) * p)) + _lisk->Li(5, -exp((-a + b + cz + d) * p)) - _lisk->Li(5, -exp((a + b + cz + d) * p)))) /
                            (pow(a, 2) * pow(b, 2) * cz * pow(p, 5));
          // yz
          _b_vector[2][5] = (2. * (a * cz * pow(p, 2) * _lisk->Li(3, -exp((-a - b - cz + d) * p)) + a * cz * pow(p, 2) * _lisk->Li(3, -exp((a - b - cz + d) * p)) -
                                   a * cz * pow(p, 2) * _lisk->Li(3, -exp((-a + b - cz + d) * p)) - a * cz * pow(p, 2) * _lisk->Li(3, -exp((a + b - cz + d) * p)) + a * cz * pow(p, 2) * _lisk->Li(3, -exp((-a - b + cz + d) * p)) +
                                   a * cz * pow(p, 2) * _lisk->Li(3, -exp((a - b + cz + d) * p)) - a * cz * pow(p, 2) * _lisk->Li(3, -exp((-a + b + cz + d) * p)) - a * cz * pow(p, 2) * _lisk->Li(3, -exp((a + b + cz + d) * p)) +
                                   a * p * _lisk->Li(4, -exp((-a - b - cz + d) * p)) + cz * p * _lisk->Li(4, -exp((-a - b - cz + d) * p)) + a * p * _lisk->Li(4, -exp((a - b - cz + d) * p)) - cz * p * _lisk->Li(4, -exp((a - b - cz + d) * p)) -
                                   a * p * _lisk->Li(4, -exp((-a + b - cz + d) * p)) - cz * p * _lisk->Li(4, -exp((-a + b - cz + d) * p)) - a * p * _lisk->Li(4, -exp((a + b - cz + d) * p)) + cz * p * _lisk->Li(4, -exp((a + b - cz + d) * p)) -
                                   a * p * _lisk->Li(4, -exp((-a - b + cz + d) * p)) + cz * p * _lisk->Li(4, -exp((-a - b + cz + d) * p)) - a * p * _lisk->Li(4, -exp((a - b + cz + d) * p)) - cz * p * _lisk->Li(4, -exp((a - b + cz + d) * p)) +
                                   a * p * _lisk->Li(4, -exp((-a + b + cz + d) * p)) - cz * p * _lisk->Li(4, -exp((-a + b + cz + d) * p)) + a * p * _lisk->Li(4, -exp((a + b + cz + d) * p)) + cz * p * _lisk->Li(4, -exp((a + b + cz + d) * p)) +
                                   _lisk->Li(5, -exp((-a - b - cz + d) * p)) - _lisk->Li(5, -exp((a - b - cz + d) * p)) - _lisk->Li(5, -exp((-a + b - cz + d) * p)) + _lisk->Li(5, -exp((a + b - cz + d) * p)) -
                                   _lisk->Li(5, -exp((-a - b + cz + d) * p)) + _lisk->Li(5, -exp((a - b + cz + d) * p)) + _lisk->Li(5, -exp((-a + b + cz + d) * p)) - _lisk->Li(5, -exp((a + b + cz + d) * p)))) /
                            (pow(a, 2) * b * pow(cz, 2) * pow(p, 5));
          // yz
          _b_vector[2][6] = (2. * (b * cz * pow(p, 2) * _lisk->Li(3, -exp((-a - b - cz + d) * p)) - b * cz * pow(p, 2) * _lisk->Li(3, -exp((a - b - cz + d) * p)) +
                                   b * cz * pow(p, 2) * _lisk->Li(3, -exp((-a + b - cz + d) * p)) - b * cz * pow(p, 2) * _lisk->Li(3, -exp((a + b - cz + d) * p)) + b * cz * pow(p, 2) * _lisk->Li(3, -exp((-a - b + cz + d) * p)) -
                                   b * cz * pow(p, 2) * _lisk->Li(3, -exp((a - b + cz + d) * p)) + b * cz * pow(p, 2) * _lisk->Li(3, -exp((-a + b + cz + d) * p)) - b * cz * pow(p, 2) * _lisk->Li(3, -exp((a + b + cz + d) * p)) +
                                   b * p * _lisk->Li(4, -exp((-a - b - cz + d) * p)) + cz * p * _lisk->Li(4, -exp((-a - b - cz + d) * p)) - b * p * _lisk->Li(4, -exp((a - b - cz + d) * p)) - cz * p * _lisk->Li(4, -exp((a - b - cz + d) * p)) +
                                   b * p * _lisk->Li(4, -exp((-a + b - cz + d) * p)) - cz * p * _lisk->Li(4, -exp((-a + b - cz + d) * p)) - b * p * _lisk->Li(4, -exp((a + b - cz + d) * p)) + cz * p * _lisk->Li(4, -exp((a + b - cz + d) * p)) -
                                   b * p * _lisk->Li(4, -exp((-a - b + cz + d) * p)) + cz * p * _lisk->Li(4, -exp((-a - b + cz + d) * p)) + b * p * _lisk->Li(4, -exp((a - b + cz + d) * p)) - cz * p * _lisk->Li(4, -exp((a - b + cz + d) * p)) -
                                   b * p * _lisk->Li(4, -exp((-a + b + cz + d) * p)) - cz * p * _lisk->Li(4, -exp((-a + b + cz + d) * p)) + b * p * _lisk->Li(4, -exp((a + b + cz + d) * p)) + cz * p * _lisk->Li(4, -exp((a + b + cz + d) * p)) +
                                   _lisk->Li(5, -exp((-a - b - cz + d) * p)) - _lisk->Li(5, -exp((a - b - cz + d) * p)) - _lisk->Li(5, -exp((-a + b - cz + d) * p)) + _lisk->Li(5, -exp((a + b - cz + d) * p)) -
                                   _lisk->Li(5, -exp((-a - b + cz + d) * p)) + _lisk->Li(5, -exp((a - b + cz + d) * p)) + _lisk->Li(5, -exp((-a + b + cz + d) * p)) - _lisk->Li(5, -exp((a + b + cz + d) * p)))) /
                            (a * pow(b, 2) * pow(cz, 2) * pow(p, 5));
          // xx
          _b_vector[2][7] = (-2. * (4 * pow(a, 3) * b * cz * pow(p, 5) - 3 * pow(a, 2) * pow(p, 2) * _lisk->Li(3, -exp((-a - b - cz + d) * p)) +
                                    3 * pow(a, 2) * pow(p, 2) * _lisk->Li(3, -exp((a - b - cz + d) * p)) + 3 * pow(a, 2) * pow(p, 2) * _lisk->Li(3, -exp((-a + b - cz + d) * p)) -
                                    3 * pow(a, 2) * pow(p, 2) * _lisk->Li(3, -exp((a + b - cz + d) * p)) + 3 * pow(a, 2) * pow(p, 2) * _lisk->Li(3, -exp((-a - b + cz + d) * p)) -
                                    3 * pow(a, 2) * pow(p, 2) * _lisk->Li(3, -exp((a - b + cz + d) * p)) - 3 * pow(a, 2) * pow(p, 2) * _lisk->Li(3, -exp((-a + b + cz + d) * p)) +
                                    3 * pow(a, 2) * pow(p, 2) * _lisk->Li(3, -exp((a + b + cz + d) * p)) - 6 * a * p * _lisk->Li(4, -exp((-a - b - cz + d) * p)) - 6 * a * p * _lisk->Li(4, -exp((a - b - cz + d) * p)) +
                                    6 * a * p * _lisk->Li(4, -exp((-a + b - cz + d) * p)) + 6 * a * p * _lisk->Li(4, -exp((a + b - cz + d) * p)) + 6 * a * p * _lisk->Li(4, -exp((-a - b + cz + d) * p)) +
                                    6 * a * p * _lisk->Li(4, -exp((a - b + cz + d) * p)) - 6 * a * p * _lisk->Li(4, -exp((-a + b + cz + d) * p)) - 6 * a * p * _lisk->Li(4, -exp((a + b + cz + d) * p)) - 6. * _lisk->Li(5, -exp((-a - b - cz + d) * p)) +
                                    6. * _lisk->Li(5, -exp((a - b - cz + d) * p)) + 6. * _lisk->Li(5, -exp((-a + b - cz + d) * p)) - 6. * _lisk->Li(5, -exp((a + b - cz + d) * p)) + 6. * _lisk->Li(5, -exp((-a - b + cz + d) * p)) -
                                    6. * _lisk->Li(5, -exp((a - b + cz + d) * p)) - 6. * _lisk->Li(5, -exp((-a + b + cz + d) * p)) + 6. * _lisk->Li(5, -exp((a + b + cz + d) * p)))) / (3. * pow(a, 3) * b * cz * pow(p, 5));
          // yy
          _b_vector[2][8] = (-2. * (4 * a * pow(b, 3) * cz * pow(p, 5) - 3 * pow(b, 2) * pow(p, 2) * _lisk->Li(3, -exp(-((a + b + cz - d) * p))) + 3 * pow(b, 2) * pow(p, 2) * _lisk->Li(3, -exp((a - b - cz + d) * p)) +
                                    3 * pow(b, 2) * pow(p, 2) * _lisk->Li(3, -exp((-a + b - cz + d) * p)) - 3 * pow(b, 2) * pow(p, 2) * _lisk->Li(3, -exp((a + b - cz + d) * p)) +
                                    3 * pow(b, 2) * pow(p, 2) * _lisk->Li(3, -exp((-a - b + cz + d) * p)) - 3 * pow(b, 2) * pow(p, 2) * _lisk->Li(3, -exp((a - b + cz + d) * p)) -
                                    3 * pow(b, 2) * pow(p, 2) * _lisk->Li(3, -exp((-a + b + cz + d) * p)) + 3 * pow(b, 2) * pow(p, 2) * _lisk->Li(3, -exp((a + b + cz + d) * p)) - 6. * b * p * _lisk->Li(4, -exp(-((a + b + cz - d) * p))) +
                                    6. * b * p * _lisk->Li(4, -exp((a - b - cz + d) * p)) - 6. * b * p * _lisk->Li(4, -exp((-a + b - cz + d) * p)) + 6. * b * p * _lisk->Li(4, -exp((a + b - cz + d) * p)) +
                                    6. * b * p * _lisk->Li(4, -exp((-a - b + cz + d) * p)) - 6. * b * p * _lisk->Li(4, -exp((a - b + cz + d) * p)) + 6. * b * p * _lisk->Li(4, -exp((-a + b + cz + d) * p)) -
                                    6. * b * p * _lisk->Li(4, -exp((a + b + cz + d) * p)) - 6. * _lisk->Li(5, -exp(-((a + b + cz - d) * p))) + 6. * _lisk->Li(5, -exp((a - b - cz + d) * p)) + 6. * _lisk->Li(5, -exp((-a + b - cz + d) * p)) -
                                    6. * _lisk->Li(5, -exp((a + b - cz + d) * p)) + 6. * _lisk->Li(5, -exp((-a - b + cz + d) * p)) - 6. * _lisk->Li(5, -exp((a - b + cz + d) * p)) - 6. * _lisk->Li(5, -exp((-a + b + cz + d) * p)) +
                                    6. * _lisk->Li(5, -exp((a + b + cz + d) * p)))) / (3. * a * pow(b, 3) * cz * pow(p, 5));
          // zz
          _b_vector[2][9] = (-2. * (4 * a * b * pow(cz, 3) * pow(p, 5) - 3 * pow(cz, 2) * pow(p, 2) * _lisk->Li(3, -exp(-((a + b + cz - d) * p))) + 3 * pow(cz, 2) * pow(p, 2) * _lisk->Li(3, -exp((a - b - cz + d) * p)) +
                                    3 * pow(cz, 2) * pow(p, 2) * _lisk->Li(3, -exp((-a + b - cz + d) * p)) - 3 * pow(cz, 2) * pow(p, 2) * _lisk->Li(3, -exp((a + b - cz + d) * p)) +
                                    3 * pow(cz, 2) * pow(p, 2) * _lisk->Li(3, -exp((-a - b + cz + d) * p)) - 3 * pow(cz, 2) * pow(p, 2) * _lisk->Li(3, -exp((a - b + cz + d) * p)) -
                                    3 * pow(cz, 2) * pow(p, 2) * _lisk->Li(3, -exp((-a + b + cz + d) * p)) + 3 * pow(cz, 2) * pow(p, 2) * _lisk->Li(3, -exp((a + b + cz + d) * p)) - 6. * cz * p * _lisk->Li(4, -exp(-((a + b + cz - d) * p))) +
                                    6. * cz * p * _lisk->Li(4, -exp((a - b - cz + d) * p)) + 6. * cz * p * _lisk->Li(4, -exp((-a + b - cz + d) * p)) - 6. * cz * p * _lisk->Li(4, -exp((a + b - cz + d) * p)) -
                                    6. * cz * p * _lisk->Li(4, -exp((-a - b + cz + d) * p)) + 6. * cz * p * _lisk->Li(4, -exp((a - b + cz + d) * p)) + 6. * cz * p * _lisk->Li(4, -exp((-a + b + cz + d) * p)) -
                                    6. * cz * p * _lisk->Li(4, -exp((a + b + cz + d) * p)) - 6. * _lisk->Li(5, -exp(-((a + b + cz - d) * p))) + 6. * _lisk->Li(5, -exp((a - b - cz + d) * p)) + 6. * _lisk->Li(5, -exp((-a + b - cz + d) * p)) -
                                    6. * _lisk->Li(5, -exp((a + b - cz + d) * p)) + 6. * _lisk->Li(5, -exp((-a - b + cz + d) * p)) - 6. * _lisk->Li(5, -exp((a - b + cz + d) * p)) - 6. * _lisk->Li(5, -exp((-a + b + cz + d) * p)) +
                                    6. * _lisk->Li(5, -exp((a + b + cz + d) * p)))) / (3. * a * b * pow(cz, 3) * pow(p, 5));
                                    
                                    
                                    
                                     for(unsigned i = 0; i < n; i++) {
              std::cout << _b_vector[2][i] <<  "  b vector ******************************" << endl;


            }



          std::vector < complex < double > > bcolumn(_b_vector[2].begin(), _b_vector[2].end());


          if(degree <= 2) EquivalentPolynomial::MatrixVectorMultiply(_A32_inverse, bcolumn, _coefficients);
          //else if(degree == 3) EquivalentPolynomial::MatrixVectorMultiply(_A23_inverse, bcolumn, _coefficients);
          //else if(degree == 4) EquivalentPolynomial::MatrixVectorMultiply(_A24_inverse, bcolumn, _coefficients);
        }

        if(degree > 2) {

        }
      }

    }



  }
  else {
    std::cout << "Wrong Dimension " << dim << std::endl;
    abort();
  }
}

void EquivalentPolynomial::MatrixVectorMultiply(const std::vector<std::vector <double>> &A, const std::vector < complex < double > > &bv, std::vector < complex < double > > &xv) {

  xv.resize(bv.size());
  for(unsigned i = 0; i < bv.size(); i++) {
    xv[i] = 0;
    for(unsigned j = 0; j < bv.size(); j++) {
      xv[i] += A[i][j] * bv[j].real();
      //std::cout << xv[i] <<  "  made it here in MatrixVectorMultiply" << endl;
    }
    //std::cout << "  MatrixVectorMultiply " << bv[i] <<  std::endl;
  }

  std::cout <<  "  made it here OUT OG MatrixVectorMultiply" <<  _dim << endl;
}

//This function takes a vector of points as inputs and calculates the best fit plane for those points

// void EquivalentPolynomial::FindBestFit(const std::vector < double > &pts, const std::vector < double > &Npts, const unsigned &dim) {
//
//
//   unsigned cnt = 0;
//   double normaldotcoefficients = 0.;
//   double d = 0;
//   unsigned numberofpoints = pts.size() / dim;
//   MatrixXd m(numberofpoints, dim);
//   _bestfit.resize(dim + 1);
//   std::vector < double > N(dim, 0.);
//   std::vector < double > centroid(dim, 0.);
//
// //Calculate average Normal and centroid from points
//   for(unsigned i = 0; i < numberofpoints; i++) {
//
//     for(unsigned j = 0; j < dim; j++, cnt++) {
//       centroid[j] += pts[cnt];
//       N[j] += Npts[cnt];
//     }
//
//   }
//   for(unsigned j = 0; j < dim; j++) {
//     N[j] /= numberofpoints;
//     centroid[j] /= numberofpoints;
//   }
//
//   cnt = 0;
//
//   //Fill matrix to be passed to JacobiSVD
//   for(unsigned i = 0; i < numberofpoints; i++) {
//
//     for(unsigned j = 0; j < dim; j++, cnt++) {
//       m(i, j) = pts[cnt] - centroid[j];
//     }
//
//   }
//
//   JacobiSVD<MatrixXd> svd(m, ComputeThinU | ComputeThinV);
//   MatrixXd v = svd.matrixV();
//
//
// //If dim = 2 and line of best fit is desired, use singualr vector associated with max singular vector
//
//   if(dim <= 2) {
//
//
//     for(unsigned i = 0; i < dim; i++) {
//
//       _bestfit[i] = v(i, 0);
//       normaldotcoefficients += _bestfit[i] * N[i];
//     }
//
//   }
//
//   //If dim = 3 and plane of best fit is desired, use singualr vector associated with min singular vector
//   if(dim == 3) {
//
//     for(unsigned i = 0; i < dim; i++) {
//
//       _bestfit[i] = v(i, dim - 1);
//       normaldotcoefficients += _bestfit[i] * N[i];
//     }
//
//   }
//
// //   for(unsigned i = 0; i < dim; i++) {
// //
// //     std::cout << " coefficent before dot product "  << _bestfit[i] << endl;
// //   }
//
//   //Rotate normal by pi if Normal dot coefficents is less than zero
//   if(normaldotcoefficients < 0) {
//
//     for(unsigned i = 0; i < dim; i++) {
//       _bestfit[i] *= -1.;
//
//     }
//
//   }
//
// //   double numerator = 0;
// //   double denominator = 0;
// //
// //   for(unsigned i = 0; i < numberofpoints - 1; i++) {
// //
// //     numerator += (pts[2*i] - centroid[0]) * (pts[2*i + 1] - centroid[1]);
// //     denominator += (pts[2*i + 1] - centroid[1]) * (pts[2*i + 1] - centroid[1]) - ((pts[2*i] - centroid[0]) * (pts[2*i] - centroid[0]));
// //
// //   }
// //
// //   if(denominator != 0.) {
// //
// //     _bestfit[0] = cos(0.5 * atan(2*numerator / denominator));
// //     _bestfit[1] = sin(0.5 * atan(2*numerator / denominator));
// //
// //   }
// //
// //Calculate constant d in ax+by+d=0 or ax+by+cz+d=0
//   for(unsigned i = 0; i < dim; i++) {
//
//     d -= _bestfit[i] * centroid[i];
//   }
//
//   _bestfit[dim] = d;
//
// //   for(unsigned i = 0; i < dim + 1; i++) {
// //
// //     std::cout << " coefficent "  << _bestfit[i] << endl;
// //   }
//
//
//
//
//   //std::cout << _bestfit[0] * _bestfit[0] + _bestfit[2] * _bestfit[2] + _bestfit[1] * _bestfit[1]  << "   = norm squared" << endl;
//   //std::cout << v << "   v matrix" << endl;
//   //std::cout << _bestfit[0] * 1. + _bestfit[1] * 2.  + _bestfit[2] << "  check equation" << std::endl;
//
//
//
// }


double EquivalentPolynomial::GetValue(std::vector <double> &x, unsigned & element) {

  unsigned numvalues = x.size();
  unsigned numcoefficients = _coefficients.size();
  //std::vector < double > values(numvalues);
  //std::fill(values.begin(), values.end(), 0.);
  double value = 0.;

  if(_dim == 1) {

    for(unsigned k = 0; k < numvalues; k++) {

      for(unsigned i = 0; i < numcoefficients; i++) {

        value += _coefficients[i].real() * pow(x[k], i);

      }

    }

  }

  if(_dim == 2) {

    if(element == 3) {
      if(_degree == 2) {

        value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] +
                _coefficients[3].real() * x[0] * x[0] + _coefficients[4].real() * x[1] * x[1] + _coefficients[5].real() * x[0] * x[1];

      }

      if(_degree == 3) {

        value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] +
                _coefficients[3].real() * x[0] * x[0] + _coefficients[4].real() * x[1] * x[1] + _coefficients[5].real() * x[0] * x[1]  +
                _coefficients[6].real() * x[0] * x[0] * x[0] + _coefficients[7].real() * x[1] * x[1] * x[1] +
                _coefficients[8].real() * x[0] * x[0] * x[1] + _coefficients[9].real() * x[1] * x[1] * x[0];

      }

      if(_degree == 4) {

        value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] +
                _coefficients[3].real() * x[0] * x[0] + _coefficients[4].real() * x[1] * x[1] + _coefficients[5].real() * x[0] * x[1] +
                _coefficients[6].real() * x[0] * x[0] * x[0] + _coefficients[7].real() * x[1] * x[1] * x[1] + _coefficients[8].real() * x[0] * x[0] * x[1] + _coefficients[9].real() * x[1] * x[1] * x[0]  + _coefficients[10].real() * x[0] * x[0] * x[0] * x[0] +
                _coefficients[11].real() * x[1] * x[1] * x[1] * x[1] + _coefficients[12].real() * x[0] * x[0] * x[1] * x[1] +
                _coefficients[13].real() * x[0] * x[0] * x[0] * x[1] + _coefficients[14].real() * x[1] * x[1] * x[1] * x[0];

      }

      if(_degree == 5) {

        value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] +
                _coefficients[3].real() * x[0] * x[0] + _coefficients[4].real() * x[1] * x[1] + _coefficients[5].real() * x[0] * x[0] * x[0] +
                _coefficients[6].real() * x[1] * x[1] * x[1] + _coefficients[7].real() * x[0] * x[0] * x[0] * x[0] + _coefficients[8].real() * x[1] * x[1] * x[1] * x[1] +
                _coefficients[9].real() * pow(x[0], 5) + _coefficients[10].real() * pow(x[1], 5) +
                _coefficients[11].real() * x[0] * x[0] * x[1] + _coefficients[12].real() * x[0] * x[0] * x[0] * x[1] +
                _coefficients[13].real() * x[0] * x[0] * x[0] * x[0] * x[1] + _coefficients[14].real() * x[1] * x[1] * x[0] +
                _coefficients[15].real() * x[1] * x[1] * x[1] * x[0] + _coefficients[16].real() * x[1] * x[1] * x[1] * x[1] * x[0] +
                _coefficients[17].real() * x[0] * x[0] * x[0] * x[1] * x[1] + _coefficients[18].real() * x[1] * x[1] * x[1] * x[0] * x[0] +
                _coefficients[19].real() * x[0] * x[1] + _coefficients[20].real() * x[0] * x[0] * x[1] * x[1];

      }

    }

    if(element == 4) {

      value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] +
              _coefficients[3].real() * x[0] * x[0] + _coefficients[4].real() * x[1] * x[1] + _coefficients[5].real() * x[0] * x[1];



    }
  }

  if(_dim == 3) {

    if(element == 0) {

      value = _coefficients[0].real() + _coefficients[1].real() * x[0] + _coefficients[2].real() * x[1] +
              _coefficients[3].real() * x[2] + _coefficients[4].real() * x[0] * x[1] + _coefficients[5].real() * x[0] * x[2] +
              _coefficients[6].real() * x[1] * x[2] + _coefficients[7].real() * x[0] * x[0] +  _coefficients[8].real() * x[1] * x[1] +
              _coefficients[9].real() * x[2] * x[2];
    }
  }

  return value;

}















std::vector<std::vector<double>> EquivalentPolynomial::_A13_inverse = {{9. / 8., 0, -(15. / 8), 0}, {0, 75. / 8, 0, -(105. / 8)}, {-(15. / 8), 0, 45. / 8, 0}, {0, -(105. / 8), 0, 175. / 8}};

std::vector<std::vector<double>> EquivalentPolynomial::_A22_inverse = {{7. / 8, 0, 0, -(15. / 16), -(15. / 16), 0}, {0, 3. / 4, 0, 0, 0, 0}, {0, 0, 3. / 4, 0, 0, 0},
  {-(15. / 16), 0, 0, 45. / 16, 0, 0}, {-(15. / 16), 0, 0, 0, 45. / 16, 0}, {0, 0, 0, 0, 0, 9. / 4}
};

std::vector<std::vector<double>> EquivalentPolynomial::_A23_inverse = {{7. / 8, 0, 0, -(15. / 16), -(15. / 16), 0, 0, 0, 0, 0}, {
    0, 45. / 8, 0, 0, 0,
    0, -(105. / 16), 0, 0, -(45. / 16)
  }, {
    0, 0, 45. / 8, 0, 0, 0,
    0, -(105. / 16), -(45. / 16), 0
  }, {
    -(15. / 16), 0, 0, 45. / 16, 0, 0, 0, 0, 0,
      0
    }, {-(15. / 16), 0, 0, 0, 45. / 16, 0, 0, 0, 0, 0}, {
    0, 0, 0, 0, 0, 9. / 4,
    0, 0, 0, 0
  }, {0, -(105. / 16), 0, 0, 0, 0, 175. / 16, 0, 0, 0}, {
    0,
    0, -(105. / 16), 0, 0, 0, 0, 175. / 16, 0, 0
  }, {
    0, 0, -(45. / 16), 0, 0, 0,
    0, 0, 135. / 16, 0
  }, {0, -(45. / 16), 0, 0, 0, 0, 0, 0, 0, 135. / 16}
};

std::vector<std::vector<double>> EquivalentPolynomial::_A24_inverse = {{
    243. / 128, 0, 0, -(675. / 128), -(675. / 128), 0, 0, 0, 0, 0, 945. / 256, 945. /
    256, 225. / 64, 0, 0
  }, {
    0, 45. / 8, 0, 0, 0, 0, -(105. / 16), 0, 0, -(45. / 16),
    0, 0, 0, 0, 0
  }, {
    0, 0, 45. / 8, 0, 0, 0, 0, -(105. / 16), -(45. / 16), 0, 0,
    0, 0, 0, 0
  }, {
    -(675. / 128), 0, 0, 1215. / 32, 225. / 64, 0, 0, 0, 0,
      0, -(4725. / 128), 0, -(675. / 64), 0, 0
    }, {
    -(675. / 128), 0, 0, 225. / 64,
      1215. / 32, 0, 0, 0, 0, 0, 0, -(4725. / 128), -(675. / 64), 0, 0
    }, {
    0, 0, 0,
    0, 0, 207. / 8, 0, 0, 0, 0, 0, 0,
    0, -(315. / 16), -(315. / 16)
  }, {
    0, -(105. / 16), 0, 0, 0, 0, 175. / 16, 0, 0,
    0, 0, 0, 0, 0, 0
  }, {
    0, 0, -(105. / 16), 0, 0, 0, 0, 175. / 16, 0, 0, 0, 0,
    0, 0, 0
  }, {
    0, 0, -(45. / 16), 0, 0, 0, 0, 0, 135. / 16, 0, 0, 0, 0, 0,
    0
  }, {0, -(45. / 16), 0, 0, 0, 0, 0, 0, 0, 135. / 16, 0, 0, 0, 0, 0}, {
    945. /
    256, 0, 0, -(4725. / 128), 0, 0, 0, 0, 0, 0, 11025. / 256, 0, 0, 0,
    0
  }, {
    945. / 256, 0, 0, 0, -(4725. / 128), 0, 0, 0, 0, 0, 0, 11025. / 256, 0,
    0, 0
  }, {
    225. / 64, 0, 0, -(675. / 64), -(675. / 64), 0, 0, 0, 0, 0, 0, 0,
    2025. / 64, 0, 0
  }, {
    0, 0, 0, 0, 0, -(315. / 16), 0, 0, 0, 0, 0, 0, 0, 525. /
    16, 0
  }, {0, 0, 0, 0, 0, -(315. / 16), 0, 0, 0, 0, 0, 0, 0, 0, 525. / 16}
};

std::vector<std::vector<double>> EquivalentPolynomial::_A22T_inverse = {
  {72, -240, -240, 180, 180, 360},
  {-240, 1200, 600, -1080, -360, -1440},
  {-240, 600, 1200, -360, -1080, -1440},
  {180, -1080, -360, 1080, 180, 1080},
  {180, -360, -1080, 180, 1080, 1080},
  {
    360, -1440, -1440, 1080, 1080, 2880
  }
};


std::vector<std::vector<double>> EquivalentPolynomial::_A32_inverse = {{19. / 32, 0, 0, 0, 0, 0, 0, -(15. / 32), -(15. / 32), -(15. / 32)}, {
    0, 3. / 8, 0,
    0, 0, 0, 0, 0, 0, 0
  }, {0, 0, 3. / 8, 0, 0, 0, 0, 0, 0, 0}, {
    0, 0, 0, 3. /
    8, 0, 0, 0, 0, 0, 0
  }, {0, 0, 0, 0, 9. / 8, 0, 0, 0, 0, 0}, {
    0, 0, 0, 0,
    0, 9. / 8, 0, 0, 0, 0
  }, {0, 0, 0, 0, 0, 0, 9. / 8, 0, 0, 0}, {
    -(15. / 32),
      0, 0, 0, 0, 0, 0, 45. / 32, 0, 0
    }, {
    -(15. / 32), 0, 0, 0, 0, 0, 0, 0, 45. /
      32, 0
    }, {-(15. / 32), 0, 0, 0, 0, 0, 0, 0, 0, 45. / 32}
};



std::vector<std::vector<double>> EquivalentPolynomial::_ATA_deg6_hex = {{3.56104, 0., 0., 0., -12.4292, 0., 0., -12.4292, 0., -12.4292, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 16.4575, 0., 0., 15.8203, 0., 

  15.8203, 0., 0., 0., 0., 16.4575, 0., 15.8203, 0., 16.4575, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -7.33154, 0., 0., -6.92139, 0., -6.92139, 0., 0., 0., 

  0., -6.92139, 0., -6.5918, 0., -6.92139, 0., 0., 0., 0., 0., 

  0., -7.33154, 0., -6.92139, 0., -6.92139, 0., -7.33154}, {0., 

  14.5723, 0., 0., 0., 0., 0., 0., 0., 0., -34.043, 0., 0., -15.293, 

  0., -15.293, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 20.3027, 0., 0., 12.3047, 0., 12.3047, 0., 0., 

  0., 0., 5.53711, 0., 5.27344, 0., 5.53711, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 14.5723, 0., 

  0., 0., 0., 0., 0., 0., 0., -15.293, 0., 0., 0., 0., -34.043, 

  0., -15.293, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 5.53711, 0., 0., 0., 0., 12.3047, 0., 5.27344, 0., 0., 

  0., 0., 0., 0., 20.3027, 0., 12.3047, 0., 5.53711, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 14.5723, 0., 0., 0., 

  0., 0., 0., 0., 0., -15.293, 0., 0., 0., 0., -15.293, 0., -34.043, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  5.53711, 0., 0., 0., 0., 5.27344, 0., 12.3047, 0., 0., 0., 0., 0., 

  0., 5.53711, 0., 12.3047, 0., 20.3027, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0.}, {-12.4292, 0., 0., 0., 136.033, 0., 0., 15.8203, 

  0., 15.8203, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -274.548, 0., 

  0., -88.9893, 0., -88.9893, 0., 0., 0., 0., -6.92139, 0., -6.5918, 

  0., -6.92139, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 153.962, 0., 0., 69.2139, 0., 

  69.2139, 0., 0., 0., 0., 20.7642, 0., 19.7754, 0., 20.7642, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 

  90.5449, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -141.504, 0., 0., 0., 0., -141.504, 0., -62.7539, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 60.9082, 0., 0., 0., 0., 86.1328, 

  0., 36.9141, 0., 0., 0., 0., 0., 0., 60.9082, 0., 36.9141, 0., 

  16.6113, 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 

  90.5449, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -141.504, 0., 0., 0., 0., -62.7539, 0., -141.504, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 60.9082, 0., 0., 0., 0., 36.9141, 

  0., 86.1328, 0., 0., 0., 0., 0., 0., 16.6113, 0., 36.9141, 0., 

  60.9082, 0., 0., 0., 0., 0., 0., 0.}, {-12.4292, 0., 0., 0., 

  15.8203, 0., 0., 136.033, 0., 15.8203, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -6.92139, 0., 0., -88.9893, 0., -6.5918, 0., 0., 0., 

  0., -274.548, 0., -88.9893, 0., -6.92139, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 20.7642, 0., 0., 0., 0., 0., 0., 69.2139, 0., 19.7754, 0., 0., 

  0., 0., 0., 0., 0., 0., 153.962, 0., 69.2139, 0., 20.7642, 0., 

  0.}, {0., 0., 0., 0., 0., 0., 0., 0., 90.5449, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -62.7539, 0., 0., 0., 0., 

  0., 0., -141.504, 0., -141.504, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  16.6113, 0., 0., 0., 0., 0., 0., 36.9141, 0., 36.9141, 0., 0., 0., 

  0., 0., 0., 0., 0., 60.9082, 0., 86.1328, 0., 60.9082, 

  0.}, {-12.4292, 0., 0., 0., 15.8203, 0., 0., 15.8203, 0., 136.033, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -6.92139, 0., 0., -6.5918, 

  0., -88.9893, 0., 0., 0., 0., -6.92139, 0., -88.9893, 0., -274.548, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 20.7642, 0., 0., 0., 0., 0., 0.,

   19.7754, 0., 69.2139, 0., 0., 0., 0., 0., 0., 0., 0., 20.7642, 0., 

  69.2139, 0., 153.962}, {0., -34.043, 0., 0., 0., 0., 0., 0., 0., 0.,

   124.414, 0., 0., 12.3047, 0., 12.3047, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -94.7461, 0., 

  0., -20.5078, 0., -20.5078, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., -15.293, 0., 0., 0., 0., 0., 0., 0., 0., 79.1016, 0., 

  0., 0., 0., 12.3047, 0., 5.27344, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -55.3711, 0., 0., 0., 

  0., -36.9141, 0., -15.8203, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 

  0., -15.293, 0., 0., 0., 0., 0., 0., 0., 0., 79.1016, 0., 0., 0., 

  0., 5.27344, 0., 12.3047, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -55.3711, 0., 0., 0., 0., -15.8203, 

  0., -36.9141, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., -15.293, 0., 0., 0., 0., 

  0., 0., 0., 0., 12.3047, 0., 0., 79.1016, 0., 5.27344, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -36.9141, 0., 0., 0., 0., 0., 0., -55.3711, 0., -15.8203, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  56.5313, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., -29.5312, 0., 0., 0., 0., 0., 

  0., -29.5312, 0., -29.5312, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0.}, {0., -15.293, 0., 0., 0., 0., 0., 

  0., 0., 0., 12.3047, 0., 0., 5.27344, 0., 79.1016, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -36.9141, 0., 0., 0., 0., 0., 0., -15.8203, 

  0., -55.3711, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0.}, {0., 0., -34.043, 0., 0., 0., 0., 0., 0., 0., 0., 

  12.3047, 0., 0., 0., 0., 124.414, 0., 12.3047, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -20.5078, 0., 0., 0., 0., 0., 0., 0., 0., -94.7461, 

  0., -20.5078, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., -15.293, 0., 0., 0., 0., 0., 0., 0., 0., 5.27344, 

  0., 0., 0., 0., 79.1016, 0., 12.3047, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -15.8203, 0., 0., 0., 0., 0., 0., 0., 0., -55.3711, 

  0., -36.9141, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., -15.293, 0., 0., 0., 0., 0., 0., 0., 0., 5.27344, 0., 

  0., 0., 0., 12.3047, 0., 79.1016, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -15.8203, 0., 0., 0., 0., 0., 0., 0., 0., -36.9141, 

  0., -55.3711, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., -34.043, 0., 0., 0., 0., 0., 0., 0., 0., 12.3047, 

  0., 0., 0., 0., 12.3047, 0., 124.414, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -20.5078, 0., 0., 0., 0., 0., 0., 0., 0., -20.5078, 

  0., -94.7461, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {16.4575, 0., 0., 0., -274.548, 0., 0., -6.92139, 0., -6.92139,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 705.212, 0., 0., 69.2139, 

  0., 69.2139, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -461.887, 0., 0., -80.7495, 0., -80.7495, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., -141.504, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 438.867, 0., 0., 0., 0., 86.1328, 0., 

  36.9141, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -284.238, 0., 

  0., 0., 0., -143.555, 0., -61.5234, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 

  0., -141.504, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 438.867, 0., 0., 0., 0., 36.9141, 0., 86.1328, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -284.238, 0., 0., 0., 0., -61.5234, 

  0., -143.555, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0.}, {15.8203, 0., 0., 0., -88.9893, 0., 

  0., -88.9893, 0., -6.5918, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  69.2139, 0., 0., 391.553, 0., 19.7754, 0., 0., 0., 0., 69.2139, 0., 

  19.7754, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -207.642, 0., 0., 0., 

  0., 0., 0., -207.642, 0., -59.3262, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 

  0., -62.7539, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 287.93, 0., 0., 0., 0., 0., 0., 36.9141, 0., 36.9141, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -166.113, 0., 0., 0., 0., 0., 

  0., -110.742, 0., -110.742, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0.}, {15.8203, 0., 0., 0., -88.9893, 0., 0., -6.5918, 

  0., -88.9893, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 69.2139, 0., 

  0., 19.7754, 0., 391.553, 0., 0., 0., 0., 0., 0., 19.7754, 0., 

  69.2139, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -207.642, 0., 0., 0., 

  0., 0., 0., -59.3262, 0., -207.642, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., -141.504, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 86.1328, 0., 0., 0., 0.,

   438.867, 0., 36.9141, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -143.555, 0., 0., 0., 0., 0., 0., 0., 

  0., -284.238, 0., -61.5234, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 0., -62.7539, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 36.9141, 0., 0., 0., 0., 287.93, 0.,

   36.9141, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -110.742, 0., 0., 0., 0., 0., 0., 0., 0., -166.113, 

  0., -110.742, 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 

  0., -62.7539, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 36.9141, 0., 0., 0., 0., 36.9141, 0., 287.93, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -110.742, 0., 0., 0., 0., 0., 0., 0., 0., -110.742, 

  0., -166.113, 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 

  0., -141.504, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 86.1328, 0., 0., 0., 0., 36.9141, 0., 438.867, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -143.555, 0., 0., 0., 0., 0., 0., 0., 0., -61.5234, 

  0., -284.238, 0., 0., 0., 0., 0., 0., 0.}, {16.4575, 0., 0., 

  0., -6.92139, 0., 0., -274.548, 0., -6.92139, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 69.2139, 0., 0., 0., 0., 0., 0., 

  705.212, 0., 69.2139, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -80.7495, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -461.887, 0., -80.7495, 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 

  0., 0., 0., -141.504, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 36.9141, 0., 0., 0., 0., 0., 0., 438.867, 0., 

  86.1328, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., -61.5234, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -284.238, 

  0., -143.555, 0., 0., 0.}, {15.8203, 0., 0., 0., -6.5918, 0., 

  0., -88.9893, 0., -88.9893, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 19.7754, 0., 19.7754, 0., 0., 0., 0., 69.2139, 0., 

  391.553, 0., 69.2139, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -59.3262, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -207.642, 0., -207.642, 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 

  0., -141.504, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 36.9141, 0., 0., 0., 0., 0., 0., 86.1328, 0., 438.867, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -61.5234, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -143.555, 

  0., -284.238, 0.}, {16.4575, 0., 0., 0., -6.92139, 0., 0., -6.92139,

   0., -274.548, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 69.2139, 0., 0., 0., 0., 0., 0., 69.2139, 0., 705.212, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -80.7495, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -80.7495, 

  0., -461.887}, {0., 20.3027, 0., 0., 0., 0., 0., 0., 0., 

  0., -94.7461, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 85.2715, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 5.53711, 0., 0., 

  0., 0., 0., 0., 0., 0., -55.3711, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  64.5996, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 

  0., 0., 5.53711, 0., 0., 0., 0., 0., 0., 0., 0., -55.3711, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 64.5996, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0.}, {0., 12.3047, 0., 0., 0., 0., 0., 0., 0., 0., -20.5078, 

  0., 0., -36.9141, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 61.5234, 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., -29.5312, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  49.2187, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 12.3047, 0., 

  0., 0., 0., 0., 0., 0., 0., -20.5078, 0., 0., 0., 0., -36.9141, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 61.5234, 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 12.3047, 0., 0., 0., 0., 0., 0., 0., 0., -36.9141, 0.,

   0., 0., 0., -20.5078, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 61.5234, 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0.}, {0., 0., 0., 5.27344, 0., 0., 0., 0., 0., 

  0., 0., 0., -15.8203, 0., 0., 0., 0., -15.8203, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 47.4609, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 5.27344, 0., 

  0., 0., 0., 0., 0., 0., 0., -15.8203, 0., 0., 0., 0., 0., 

  0., -15.8203, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 47.4609, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0.}, {0., 0., 0., 12.3047, 0., 0., 0., 0., 0., 0., 0., 

  0., -36.9141, 0., 0., 0., 0., 0., 0., -20.5078, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 61.5234, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 5.53711, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -55.3711, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 64.5996, 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -29.5312, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 49.2187, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 5.27344, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -15.8203, 0., -15.8203, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 47.4609, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., -29.5312, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 49.2187, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 5.53711, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -55.3711, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 64.5996, 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 

  0., 20.3027, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -94.7461, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 85.2715, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0.}, {0., 0., 0., 5.53711, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., -55.3711, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 64.5996, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 12.3047, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -20.5078, 0., -36.9141, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 61.5234, 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 

  12.3047, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -36.9141, 0., -20.5078, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 61.5234, 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0.}, {0., 0., 5.53711, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -55.3711, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 64.5996, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 20.3027, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -94.7461, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  85.2715, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {-7.33154, 0.,

   0., 0., 153.962, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -461.887, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 338.717, 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0.}, {0., 0., 0., 0., 0., 60.9082, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -284.238, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 255.814, 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 60.9082, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -284.238, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  255.814, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {-6.92139, 0., 0., 0., 

  69.2139, 0., 0., 20.7642, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -80.7495, 0., 0., -207.642, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 242.249, 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 16.6113, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -166.113, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 193.799,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0.}, {-6.92139, 0., 0., 0., 69.2139, 0., 0., 

  0., 0., 20.7642, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -80.7495, 

  0., 0., 0., 0., -207.642, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 242.249, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 

  0., 0., 0., 0., 86.1328, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., -143.555, 0., 0., 0., 0., -143.555, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 239.258, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 36.9141, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -61.5234, 0., 0., 0., 

  0., -110.742, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 184.57, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 

  36.9141, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -61.5234, 0., 0., 0., 0., 0., 0., -110.742, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 184.57, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 0., 86.1328, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -143.555, 0., 0., 0., 0., 0., 

  0., -143.555, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 239.258, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0.}, {-6.92139, 0., 0., 0., 20.7642, 

  0., 0., 69.2139, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., -207.642, 0., 0., 0., 0., 0., 0., -80.7495, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 242.249, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 0., 0., 0., 36.9141, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -110.742, 0., 0., 0., 0., 

  0., 0., -61.5234, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 184.57, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0.}, {-6.5918, 0., 0., 0., 19.7754, 0., 0., 

  19.7754, 0., 19.7754, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -59.3262, 0., -59.3262, 0., 0., 0., 0., 0., 0., -59.3262, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 177.979, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0.}, {0., 0., 0., 0., 0., 0., 0., 0., 36.9141, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -110.742, 0., 0., 0., 0., 

  0., 0., 0., 0., -61.5234, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 184.57, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0.}, {-6.92139, 0., 0., 0., 20.7642, 0., 0.,

   0., 0., 69.2139, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -207.642, 0., 0., 0., 0., 0., 0., 0., 0., -80.7495, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  242.249, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 

  0., 0., 0., 0., 60.9082, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., -284.238, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 255.814, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0.}, {0., 0., 0., 0., 0., 0., 16.6113, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -166.113, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 193.799, 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 36.9141, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -61.5234, 0., -110.742, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  184.57, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 

  0., 0., 36.9141, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., -110.742, 0., -61.5234, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 184.57, 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 16.6113, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -166.113, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 193.799, 0., 0., 0.,

   0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 60.9082, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -284.238, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  255.814, 0., 0., 0., 0., 0., 0., 0.}, {-7.33154, 0., 0., 0., 0., 0.,

   0., 153.962, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -461.887, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 338.717, 0., 0., 0., 0., 0., 0.}, {0., 

  0., 0., 0., 0., 0., 0., 0., 60.9082, 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -284.238, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 255.814, 0., 0., 0.,

   0., 0.}, {-6.92139, 0., 0., 0., 0., 0., 0., 69.2139, 0., 20.7642, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -80.7495, 0., -207.642, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 242.249, 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 

  0., 0., 86.1328, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., -143.555, 0., -143.555, 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 239.258, 0., 0., 

  0.}, {-6.92139, 0., 0., 0., 0., 0., 0., 20.7642, 0., 69.2139, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., -207.642, 0., -80.7495, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 242.249, 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 

  0., 60.9082, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -284.238, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 255.814, 0.}, {-7.33154, 0.,

   0., 0., 0., 0., 0., 0., 0., 153.962, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -461.887, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 338.717}
};











 
std::vector<std::vector<double>> EquivalentPolynomial::_ATA_deg8_hex = {{6.57012, 0., 0., 0., -35.63, 0., 0., -35.63, 0., -35.63, 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 80.4659, 0., 0., 77.8931, 0., 77.8931, 

  0., 0., 0., 0., 80.4659, 0., 77.8931, 0., 80.4659, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -80.1888, 0., 0., -76.4236, 0., -76.4236, 0., 0., 0., 

  0., -76.4236, 0., -73.3337, 0., -76.4236, 0., 0., 0., 0., 0., 

  0., -80.1888, 0., -76.4236, 0., -76.4236, 0., -80.1888, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  29.2116, 0., 0., 27.4933, 0., 27.4933, 0., 0., 0., 0., 27.253, 0., 

  25.9552, 0., 27.253, 0., 0., 0., 0., 0., 0., 27.4933, 0., 25.9552, 

  0., 25.9552, 0., 27.4933, 0., 0., 0., 0., 0., 0., 0., 0., 29.2116, 

  0., 27.4933, 0., 27.253, 0., 27.4933, 0., 29.2116}, {0., 45.8276, 

  0., 0., 0., 0., 0., 0., 0., 0., -184.622, 0., 0., -85.2759, 

  0., -85.2759, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 248.708, 0., 0., 153.809, 0., 153.809, 0., 

  0., 0., 0., 71.521, 0., 68.5547, 0., 71.521, 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -109.973, 0., 

  0., -76.1353, 0., -76.1353, 0., 0., 0., 0., -48.4497, 0., -46.1426, 

  0., -48.4497, 0., 0., 0., 0., 0., 0., -21.9946, 0., -20.7642, 

  0., -20.7642, 0., -21.9946, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 45.8276, 0., 0., 0., 

  0., 0., 0., 0., 0., -85.2759, 0., 0., 0., 0., -184.622, 

  0., -85.2759, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 71.521, 0., 0., 0., 0., 153.809, 0., 68.5547, 0., 

  0., 0., 0., 0., 0., 248.708, 0., 153.809, 0., 71.521, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -21.9946, 0., 0., 0., 

  0., -48.4497, 0., -20.7642, 0., 0., 0., 0., 0., 0., -76.1353, 

  0., -46.1426, 0., -20.7642, 0., 0., 0., 0., 0., 0., 0., 

  0., -109.973, 0., -76.1353, 0., -48.4497, 0., -21.9946, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 45.8276, 0., 0., 

  0., 0., 0., 0., 0., 0., -85.2759, 0., 0., 0., 0., -85.2759, 

  0., -184.622, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 71.521, 0., 0., 0., 0., 68.5547, 0., 153.809, 0., 

  0., 0., 0., 0., 0., 71.521, 0., 153.809, 0., 248.708, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -21.9946, 0., 0., 0., 

  0., -20.7642, 0., -48.4497, 0., 0., 0., 0., 0., 0., -20.7642, 

  0., -46.1426, 0., -76.1353, 0., 0., 0., 0., 0., 0., 0., 

  0., -21.9946, 0., -48.4497, 0., -76.1353, 0., -109.973, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0.}, {-35.63, 0., 0., 0., 589.686, 0., 

  0., 77.8931, 0., 77.8931, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -2007.49, 0., 0., -692.221, 0., -692.221, 0., 0., 0., 

  0., -76.4236, 0., -73.3337, 0., -76.4236, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2501.89,

   0., 0., 1176.64, 0., 1176.64, 0., 0., 0., 0., 392.789, 0., 375.732,

   0., 392.789, 0., 0., 0., 0., 0., 0., 27.4933, 0., 25.9552, 0., 

  25.9552, 0., 27.4933, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -1051.62, 0., 0., -577.359, 

  0., -577.359, 0., 0., 0., 0., -272.53, 0., -259.552, 0., -272.53, 

  0., 0., 0., 0., 0., 0., -82.4799, 0., -77.8656, 0., -77.8656, 

  0., -82.4799, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0.}, {0., 0., 0., 0., 0., 388.042, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1027.6, 0., 0., 0., 

  0., -1027.6, 0., -467.292, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 989.758, 0., 0., 0., 0., 1378.13, 0., 609.082, 0., 0., 0., 0., 

  0., 0., 989.758, 0., 609.082, 0., 281.008, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -329.919, 0., 0., 0., 0., -532.947, 

  0., -228.406, 0., 0., 0., 0., 0., 0., -532.947, 0., -322.998, 

  0., -145.349, 0., 0., 0., 0., 0., 0., 0., 0., -329.919, 

  0., -228.406, 0., -145.349, 0., -65.9839, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 388.042, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1027.6, 0., 0., 0., 

  0., -467.292, 0., -1027.6, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 989.758, 0., 0., 0., 0., 609.082, 0., 1378.13, 0., 0., 0., 0., 

  0., 0., 281.008, 0., 609.082, 0., 989.758, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -329.919, 0., 0., 0., 0., -228.406, 

  0., -532.947, 0., 0., 0., 0., 0., 0., -145.349, 0., -322.998, 

  0., -532.947, 0., 0., 0., 0., 0., 0., 0., 0., -65.9839, 

  0., -145.349, 0., -228.406, 0., -329.919, 0., 0., 0., 0., 0., 0., 

  0., 0., 0.}, {-35.63, 0., 0., 0., 77.8931, 0., 0., 589.686, 0., 

  77.8931, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -76.4236, 0., 

  0., -692.221, 0., -73.3337, 0., 0., 0., 0., -2007.49, 0., -692.221, 

  0., -76.4236, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 27.4933, 0., 0., 392.789, 0., 

  25.9552, 0., 0., 0., 0., 1176.64, 0., 375.732, 0., 25.9552, 0., 0., 

  0., 0., 0., 0., 2501.89, 0., 1176.64, 0., 392.789, 0., 27.4933, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., -82.4799, 0., 0., 0., 0., 0., 0., -272.53, 

  0., -77.8656, 0., 0., 0., 0., 0., 0., 0., 0., -577.359, 

  0., -259.552, 0., -77.8656, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -1051.62, 0., -577.359, 0., -272.53, 0., -82.4799, 0., 0.}, {0.,

   0., 0., 0., 0., 0., 0., 0., 388.042, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -467.292, 0., 0., 0., 0., 0., 

  0., -1027.6, 0., -1027.6, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  281.008, 0., 0., 0., 0., 0., 0., 609.082, 0., 609.082, 0., 0., 0., 

  0., 0., 0., 0., 0., 989.758, 0., 1378.13, 0., 989.758, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., -65.9839, 0., 0., 0., 0., 0., 0., -145.349, 

  0., -145.349, 0., 0., 0., 0., 0., 0., 0., 0., -228.406, 

  0., -322.998, 0., -228.406, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -329.919, 0., -532.947, 0., -532.947, 0., -329.919, 

  0.}, {-35.63, 0., 0., 0., 77.8931, 0., 0., 77.8931, 0., 589.686, 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., -76.4236, 0., 0., -73.3337, 

  0., -692.221, 0., 0., 0., 0., -76.4236, 0., -692.221, 0., -2007.49, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 27.4933, 0., 0., 25.9552, 0., 392.789, 0., 0., 0., 

  0., 25.9552, 0., 375.732, 0., 1176.64, 0., 0., 0., 0., 0., 0., 

  27.4933, 0., 392.789, 0., 1176.64, 0., 2501.89, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -82.4799, 0., 0., 0., 0., 0., 0., -77.8656, 0., -272.53, 0.,

   0., 0., 0., 0., 0., 0., 0., -77.8656, 0., -259.552, 0., -577.359, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -82.4799, 0., -272.53, 

  0., -577.359, 0., -1051.62}, {0., -184.622, 0., 0., 0., 0., 0., 0., 

  0., 0., 1136.73, 0., 0., 153.809, 0., 153.809, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1930.45, 

  0., 0., -510.132, 0., -510.132, 0., 0., 0., 0., -48.4497, 

  0., -46.1426, 0., -48.4497, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 989.758, 0., 0., 355.298, 0., 355.298, 

  0., 0., 0., 0., 80.7495, 0., 76.9043, 0., 80.7495, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., -85.2759, 0., 

  0., 0., 0., 0., 0., 0., 0., 684.954, 0., 0., 0., 0., 153.809, 0., 

  68.5547, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., -1045.13, 0., 0., 0., 0., -752.124, 0., -330.249, 0., 0., 

  0., 0., 0., 0., -76.1353, 0., -46.1426, 0., -20.7642, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 461.887, 0., 0., 0., 0., 

  484.497, 0., 207.642, 0., 0., 0., 0., 0., 0., 228.406, 0., 138.428, 

  0., 62.2925, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 

  0., -85.2759, 0., 0., 0., 0., 0., 0., 0., 0., 684.954, 0., 0., 0., 

  0., 68.5547, 0., 153.809, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -1045.13, 0., 0., 0., 0., -330.249, 

  0., -752.124, 0., 0., 0., 0., 0., 0., -20.7642, 0., -46.1426, 

  0., -76.1353, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  461.887, 0., 0., 0., 0., 207.642, 0., 484.497, 0., 0., 0., 0., 0., 

  0., 62.2925, 0., 138.428, 0., 228.406, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0.}, {0., -85.2759, 0., 0., 0., 0., 0., 0., 0., 0., 153.809, 0.,

   0., 684.954, 0., 68.5547, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -76.1353, 0., 0., -752.124, 

  0., -46.1426, 0., 0., 0., 0., -1045.13, 0., -330.249, 0., -20.7642, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 228.406, 0., 0., 0., 0., 0., 0., 484.497, 0., 138.428, 

  0., 0., 0., 0., 0., 0., 0., 0., 461.887, 0., 207.642, 0., 62.2925, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 466.119, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -542.637, 0., 0., 0., 

  0., 0., 0., -542.637, 0., -542.637, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 182.725, 0.,

   0., 0., 0., 0., 0., 258.398, 0., 258.398, 0., 0., 0., 0., 0., 0., 

  0., 0., 182.725, 0., 258.398, 0., 182.725, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., -85.2759, 0., 0., 0., 0., 0., 0., 0., 0., 153.809, 0., 0.,

   68.5547, 0., 684.954, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -76.1353, 0., 0., -46.1426, 

  0., -752.124, 0., 0., 0., 0., -20.7642, 0., -330.249, 0., -1045.13, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 228.406, 0., 0., 0., 0., 0., 0., 138.428, 0., 

  484.497, 0., 0., 0., 0., 0., 0., 0., 0., 62.2925, 0., 207.642, 0., 

  461.887, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0.}, {0., 0., -184.622, 0., 0., 0., 0., 0., 0., 0., 0.,

   153.809, 0., 0., 0., 0., 1136.73, 0., 153.809, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -48.4497, 0., 0., 

  0., 0., -510.132, 0., -46.1426, 0., 0., 0., 0., 0., 0., -1930.45, 

  0., -510.132, 0., -48.4497, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 80.7495, 0., 0., 0., 0., 0., 0., 0.,

   0., 355.298, 0., 76.9043, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  989.758, 0., 355.298, 0., 80.7495, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0.}, {0., 0., 0., -85.2759, 0., 0., 0., 0., 0., 

  0., 0., 0., 68.5547, 0., 0., 0., 0., 684.954, 0., 153.809, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -20.7642, 0., 0., 0., 0., -330.249, 0., -46.1426, 0., 0., 0., 

  0., 0., 0., -1045.13, 0., -752.124, 0., -76.1353, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 62.2925, 0.,

   0., 0., 0., 0., 0., 0., 0., 207.642, 0., 138.428, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 461.887, 0., 484.497, 0., 228.406, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., -85.2759, 0., 

  0., 0., 0., 0., 0., 0., 0., 68.5547, 0., 0., 0., 0., 153.809, 0., 

  684.954, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., -20.7642, 0., 0., 0., 0., -46.1426, 0., -330.249, 0., 0., 

  0., 0., 0., 0., -76.1353, 0., -752.124, 0., -1045.13, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  62.2925, 0., 0., 0., 0., 0., 0., 0., 0., 138.428, 0., 207.642, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 228.406, 0., 484.497, 0., 

  461.887, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 

  0., -184.622, 0., 0., 0., 0., 0., 0., 0., 0., 153.809, 0., 0., 0., 

  0., 153.809, 0., 1136.73, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -48.4497, 0., 0., 0., 0., -46.1426, 

  0., -510.132, 0., 0., 0., 0., 0., 0., -48.4497, 0., -510.132, 

  0., -1930.45, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 80.7495, 0., 0., 0., 0., 0., 0., 0., 0.,

   76.9043, 0., 355.298, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  80.7495, 0., 355.298, 0., 989.758, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0.}, {80.4659, 0., 0., 0., -2007.49, 0., 0., -76.4236, 

  0., -76.4236, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 8596.79, 0., 

  0., 1176.64, 0., 1176.64, 0., 0., 0., 0., 27.253, 0., 25.9552, 0., 

  27.253, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., -12413.2, 0., 0., -2816.14, 0., -2816.14, 

  0., 0., 0., 0., -272.53, 0., -259.552, 0., -272.53, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 5783.9, 0., 0., 1732.08,

   0., 1732.08, 0., 0., 0., 0., 317.951, 0., 302.811, 0., 317.951, 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 

  0., -1027.6, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 5011.85, 0., 0., 0., 0., 1378.13, 0., 609.082, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -6928.31, 0., 0., 0., 0., -4073.36, 

  0., -1776.49, 0., 0., 0., 0., 0., 0., -532.947, 0., -322.998, 

  0., -145.349, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2969.27,

   0., 0., 0., 0., 2487.08, 0., 1065.89, 0., 0., 0., 0., 0., 0., 

  888.245, 0., 538.33, 0., 242.249, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 0., -1027.6, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 5011.85, 0., 0., 0., 0., 609.082, 

  0., 1378.13, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -6928.31, 

  0., 0., 0., 0., -1776.49, 0., -4073.36, 0., 0., 0., 0., 0., 

  0., -145.349, 0., -322.998, 0., -532.947, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 2969.27, 0., 0., 0., 0., 1065.89, 0., 2487.08, 

  0., 0., 0., 0., 0., 0., 242.249, 0., 538.33, 0., 888.245, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0.}, {77.8931, 0., 0., 0., -692.221, 0., 

  0., -692.221, 0., -73.3337, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  1176.64, 0., 0., 4433.4, 0., 375.732, 0., 0., 0., 0., 1176.64, 0., 

  375.732, 0., 25.9552, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -577.359, 0., 0., -5165.08, 

  0., -259.552, 0., 0., 0., 0., -5165.08, 0., -1594.39, 0., -77.8656, 

  0., 0., 0., 0., 0., 0., -577.359, 0., -259.552, 0., -77.8656, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 1732.08, 0., 0., 0., 0., 0., 0., 2725.3, 0.,

   778.656, 0., 0., 0., 0., 0., 0., 0., 0., 1732.08, 0., 778.656, 0., 

  233.597, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., -467.292, 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 3087.93, 0., 

  0., 0., 0., 0., 0., 609.082, 0., 609.082, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -3799.84, 0., 0., 0., 0., 0., 0., -2699.34, 

  0., -2699.34, 0., 0., 0., 0., 0., 0., 0., 0., -228.406, 

  0., -322.998, 0., -228.406, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1385.66, 0.,

   0., 0., 0., 0., 0., 1453.49, 0., 1453.49, 0., 0., 0., 0., 0., 0., 

  0., 0., 685.217, 0., 968.994, 0., 685.217, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {77.8931, 0., 0., 

  0., -692.221, 0., 0., -73.3337, 0., -692.221, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 1176.64, 0., 0., 375.732, 0., 4433.4, 0., 0., 

  0., 0., 25.9552, 0., 375.732, 0., 1176.64, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -577.359, 0., 0., -259.552, 0., -5165.08, 0., 0., 0., 

  0., -77.8656, 0., -1594.39, 0., -5165.08, 0., 0., 0., 0., 0., 0., 

  0., 0., -77.8656, 0., -259.552, 0., -577.359, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 1732.08, 0., 0., 0., 0., 0., 0., 778.656, 0., 2725.3, 0., 

  0., 0., 0., 0., 0., 0., 0., 233.597, 0., 778.656, 0., 1732.08, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., -1027.6, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 1378.13, 0., 0., 0., 0., 5011.85, 0., 

  609.082, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -532.947, 0., 

  0., 0., 0., -4073.36, 0., -322.998, 0., 0., 0., 0., 0., 

  0., -6928.31, 0., -1776.49, 0., -145.349, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 888.245, 0., 0., 0., 0., 0.,

   0., 0., 0., 2487.08, 0., 538.33, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 2969.27, 0., 1065.89, 0., 242.249, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., -467.292, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 609.082, 0., 0.,

   0., 0., 3087.93, 0., 609.082, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -228.406, 0., 0., 0., 0., -2699.34, 0., -322.998, 0., 0., 

  0., 0., 0., 0., -3799.84, 0., -2699.34, 0., -228.406, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 685.217, 0., 0.,

   0., 0., 0., 0., 0., 0., 1453.49, 0., 968.994, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 1385.66, 0., 1453.49, 0., 685.217, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., -467.292, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 609.082, 0.,

   0., 0., 0., 609.082, 0., 3087.93, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -228.406, 0., 0., 0., 0., -322.998, 0., -2699.34, 0., 

  0., 0., 0., 0., 0., -228.406, 0., -2699.34, 0., -3799.84, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  685.217, 0., 0., 0., 0., 0., 0., 0., 0., 968.994, 0., 1453.49, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 685.217, 0., 1453.49, 0., 

  1385.66, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 

  0., 0., -1027.6, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 1378.13, 0., 0., 0., 0., 609.082, 0., 5011.85, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -532.947, 0., 0., 0., 0., -322.998, 

  0., -4073.36, 0., 0., 0., 0., 0., 0., -145.349, 0., -1776.49, 

  0., -6928.31, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 888.245, 0., 0., 0., 0., 0., 0., 0., 0., 538.33,

   0., 2487.08, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 242.249, 0., 

  1065.89, 0., 2969.27, 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {80.4659,

   0., 0., 0., -76.4236, 0., 0., -2007.49, 0., -76.4236, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 27.253, 0., 0., 1176.64, 0., 25.9552, 

  0., 0., 0., 0., 8596.79, 0., 1176.64, 0., 27.253, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -272.53, 0., 0., 0., 0., 0., 0., -2816.14, 0., -259.552,

   0., 0., 0., 0., 0., 0., 0., 0., -12413.2, 0., -2816.14, 

  0., -272.53, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   317.951, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1732.08, 0., 

  302.811, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 5783.9, 0.,

   1732.08, 0., 317.951, 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0.,

   0., -1027.6, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 609.082, 0., 0., 0., 0., 0., 0., 5011.85, 0., 1378.13, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -145.349, 0., 0., 0., 0., 0., 

  0., -1776.49, 0., -322.998, 0., 0., 0., 0., 0., 0., 0., 

  0., -6928.31, 0., -4073.36, 0., -532.947, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 242.249, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 1065.89, 0., 538.33, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 2969.27, 0., 2487.08, 0., 888.245, 0., 0., 

  0.}, {77.8931, 0., 0., 0., -73.3337, 0., 0., -692.221, 0., -692.221,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 25.9552, 0., 0., 375.732, 

  0., 375.732, 0., 0., 0., 0., 1176.64, 0., 4433.4, 0., 1176.64, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., -77.8656, 0., -77.8656, 0., 0., 0., 

  0., -259.552, 0., -1594.39, 0., -259.552, 0., 0., 0., 0., 0., 

  0., -577.359, 0., -5165.08, 0., -5165.08, 0., -577.359, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 233.597, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 778.656, 0., 778.656, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 1732.08, 0., 2725.3, 0., 1732.08, 0., 

  0.}, {0., 0., 0., 0., 0., 0., 0., 0., -1027.6, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 609.082, 0., 0., 0., 0., 0.,

   0., 1378.13, 0., 5011.85, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -145.349, 0., 0., 0., 0., 0., 0., -322.998, 0., -1776.49, 0., 

  0., 0., 0., 0., 0., 0., 0., -532.947, 0., -4073.36, 0., -6928.31, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  242.249, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 538.33, 0., 

  1065.89, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 888.245, 

  0., 2487.08, 0., 2969.27, 0.}, {80.4659, 0., 0., 0., -76.4236, 0., 

  0., -76.4236, 0., -2007.49, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  27.253, 0., 0., 25.9552, 0., 1176.64, 0., 0., 0., 0., 27.253, 0., 

  1176.64, 0., 8596.79, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -272.53,

   0., 0., 0., 0., 0., 0., -259.552, 0., -2816.14, 0., 0., 0., 0., 0.,

   0., 0., 0., -272.53, 0., -2816.14, 0., -12413.2, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 317.951, 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 302.811, 0., 1732.08, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 317.951, 0., 1732.08, 0., 

  5783.9}, {0., 248.708, 0., 0., 0., 0., 0., 0., 0., 0., -1930.45, 0.,

   0., -76.1353, 0., -76.1353, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 3815.9, 0., 0., 355.298, 0., 

  355.298, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -2177.47, 0., 

  0., -319.768, 0., -319.768, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 

  71.521, 0., 0., 0., 0., 0., 0., 0., 0., -1045.13, 0., 0., 0., 

  0., -48.4497, 0., -20.7642, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 2374.04, 0., 0., 0., 0., 484.497, 0., 

  207.642, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., -1385.66, 0., 0., 0., 

  0., -565.247, 0., -242.249, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 71.521, 

  0., 0., 0., 0., 0., 0., 0., 0., -1045.13, 0., 0., 0., 0., -20.7642, 

  0., -48.4497, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 2374.04, 0., 0., 0., 0., 207.642, 0., 484.497, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., -1385.66, 0., 0., 0., 0., -242.249, 

  0., -565.247, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0.}, {0., 153.809, 0., 0., 0., 0., 0., 0., 

  0., 0., -510.132, 0., 0., -752.124, 0., -46.1426, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 355.298,

   0., 0., 2014.89, 0., 76.9043, 0., 0., 0., 0., 484.497, 0., 138.428,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., -1065.89, 0., 0., 0., 0., 0., 0., -807.495, 

  0., -230.713, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -542.637, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1513.48, 0., 0., 0.,

   0., 0., 0., 258.398, 0., 258.398, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -852.715, 

  0., 0., 0., 0., 0., 0., -430.664, 0., -430.664, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 153.809, 0., 0., 

  0., 0., 0., 0., 0., 0., -510.132, 0., 0., -46.1426, 0., -752.124, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 355.298, 0., 0., 76.9043, 0., 2014.89, 0., 0., 0., 0., 0., 

  0., 138.428, 0., 484.497, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1065.89, 0., 0., 

  0., 0., 0., 0., -230.713, 0., -807.495, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 153.809, 0., 0., 0., 0., 

  0., 0., 0., 0., -752.124, 0., 0., 0., 0., -510.132, 0., -46.1426, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  484.497, 0., 0., 0., 0., 2014.89, 0., 138.428, 0., 0., 0., 0., 0., 

  0., 355.298, 0., 76.9043, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -807.495, 0., 0., 0., 0., 

  0., 0., 0., 0., -1065.89, 0., -230.713, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0.}, {0., 0., 0., 68.5547, 0., 0., 0., 0., 0., 0., 

  0., 0., -330.249, 0., 0., 0., 0., -330.249, 0., -46.1426, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 207.642,

   0., 0., 0., 0., 1364.5, 0., 138.428, 0., 0., 0., 0., 0., 0., 

  207.642, 0., 138.428, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -622.925, 0., 0., 0., 0., 0., 

  0., 0., 0., -622.925, 0., -415.283, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0.}, {0., 0., 68.5547, 0., 0., 0., 0., 0., 0., 0., 

  0., -330.249, 0., 0., 0., 0., -46.1426, 0., -330.249, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 207.642, 0.,

   0., 0., 0., 138.428, 0., 1364.5, 0., 0., 0., 0., 0., 0., 0., 0., 

  138.428, 0., 207.642, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -622.925, 0., 0., 0., 0., 0., 

  0., 0., 0., -415.283, 0., -622.925, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0.}, {0., 0., 0., 153.809, 0., 0., 0., 0., 0., 0., 0., 

  0., -752.124, 0., 0., 0., 0., -46.1426, 0., -510.132, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 484.497, 0.,

   0., 0., 0., 138.428, 0., 2014.89, 0., 0., 0., 0., 0., 0., 0., 0., 

  76.9043, 0., 355.298, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -807.495, 0., 0., 0., 0., 0., 

  0., 0., 0., -230.713, 0., -1065.89, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 71.521, 0., 0., 0., 0., 0., 0., 0., 0., -48.4497, 0., 

  0., -1045.13, 0., -20.7642, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 484.497, 0., 0., 0.,

   0., 0., 0., 2374.04, 0., 207.642, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., -565.247, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -1385.66, 0., -242.249, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -542.637, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 258.398, 0., 0., 0., 0., 0., 0., 1513.48, 0., 

  258.398, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -430.664, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -852.715, 0., -430.664, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0.}, {0., 68.5547, 0., 0., 0., 0., 0., 0., 0., 

  0., -46.1426, 0., 0., -330.249, 0., -330.249, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  138.428, 0., 138.428, 0., 0., 0., 0., 207.642, 0., 1364.5, 0., 

  207.642, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -415.283, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -622.925, 0., -622.925, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -542.637, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 258.398, 0., 0., 0., 0.,

   0., 0., 258.398, 0., 1513.48, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -430.664, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -430.664, 0., -852.715, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 71.521, 0., 0., 0., 

  0., 0., 0., 0., 0., -48.4497, 0., 0., -20.7642, 0., -1045.13, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 484.497, 0., 0., 0., 0., 0., 0., 207.642, 

  0., 2374.04, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -565.247, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -242.249, 

  0., -1385.66, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0.}, {0., 0., 248.708, 0., 0., 0., 0., 0., 0., 

  0., 0., -76.1353, 0., 0., 0., 0., -1930.45, 0., -76.1353, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 355.298, 0., 0., 0., 0., 0., 0., 0., 0., 3815.9, 0., 

  355.298, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -319.768, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -2177.47, 

  0., -319.768, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0.}, {0., 0., 0., 71.521, 0., 0., 0., 0., 0., 0., 0., 

  0., -20.7642, 0., 0., 0., 0., -1045.13, 0., -48.4497, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 207.642, 0., 0., 0., 0., 0., 0., 0., 0., 2374.04, 0., 

  484.497, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -242.249, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1385.66, 

  0., -565.247, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0.}, {0., 0., 153.809, 0., 0., 0., 0., 0., 0., 0., 0., -46.1426,

   0., 0., 0., 0., -510.132, 0., -752.124, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  76.9043, 0., 138.428, 0., 0., 0., 0., 0., 0., 355.298, 0., 2014.89, 

  0., 484.497, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -230.713, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -1065.89, 0., -807.495, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0.}, {0., 0., 0., 153.809, 0., 0., 0., 0., 0., 0., 0., 

  0., -46.1426, 0., 0., 0., 0., -752.124, 0., -510.132, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 138.428, 0., 76.9043, 0., 0., 0., 0., 0., 0., 484.497, 0., 

  2014.89, 0., 355.298, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -230.713, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -807.495, 0., -1065.89, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0.}, {0., 0., 71.521, 0., 0., 0., 0., 0., 0., 0., 

  0., -20.7642, 0., 0., 0., 0., -48.4497, 0., -1045.13, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 207.642, 0., 0., 0., 0., 0., 0., 0., 0., 484.497, 

  0., 2374.04, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -242.249, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -565.247, 0., -1385.66, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0.}, {0., 0., 0., 248.708, 0., 0., 0., 0., 0., 0., 0., 

  0., -76.1353, 0., 0., 0., 0., -76.1353, 0., -1930.45, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 355.298, 0., 0., 0., 0., 0., 0., 0., 0., 355.298, 

  0., 3815.9, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -319.768, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -319.768, 0., -2177.47, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {-80.1888, 0., 0., 0., 2501.89, 0., 0., 27.4933, 0., 27.4933, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -12413.2, 0., 0., -577.359, 

  0., -577.359, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  19899.6, 0., 0., 1732.08, 0., 1732.08, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -10025.4, 0., 0., -1270.19, 0., -1270.19, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 989.758, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -6928.31, 0., 0., 0., 0., -532.947, 

  0., -228.406, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 12471., 

  0., 0., 0., 0., 2487.08, 0., 1065.89, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -6532.4,

   0., 0., 0., 0., -2238.38, 0., -959.304, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 

  0., 0., 989.758, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., -6928.31, 0., 0., 0., 0., -228.406, 0., -532.947, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 12471., 0., 0., 0., 0., 1065.89,

   0., 2487.08, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -6532.4, 0., 0., 0., 

  0., -959.304, 0., -2238.38, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0.}, {-76.4236, 0., 0., 0., 1176.64, 0.,

   0., 392.789, 0., 25.9552, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -2816.14, 0., 0., -5165.08, 0., -259.552, 0., 0., 0., 

  0., -272.53, 0., -77.8656, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1732.08, 0., 0.,

   10356.1, 0., 302.811, 0., 0., 0., 0., 2725.3, 0., 778.656, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -5196.23, 0., 0., 0., 0., 0., 0., -3179.51, 0., -908.432, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 

  0., 0., 0., 0., 0., 0., 281.008, 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., -3799.84, 0., 0., 0., 0., 0., 0., -145.349,

   0., -145.349, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 7897.3, 0., 0., 

  0., 0., 0., 0., 1453.49, 0., 1453.49, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -4156.98, 0., 

  0., 0., 0., 0., 0., -1695.74, 0., -1695.74, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0.}, {-76.4236, 0., 0., 0., 1176.64, 0.,

   0., 25.9552, 0., 392.789, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -2816.14, 0., 0., -259.552, 0., -5165.08, 0., 0., 0., 0., 0., 

  0., -77.8656, 0., -272.53, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1732.08, 0., 0., 

  302.811, 0., 10356.1, 0., 0., 0., 0., 0., 0., 778.656, 0., 2725.3, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -5196.23, 0., 0., 0., 0., 0., 0., -908.432, 

  0., -3179.51, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 1378.13, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -4073.36, 0., 0., 0., 0., -4073.36, 

  0., -322.998, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2487.08,

   0., 0., 0., 0., 9749.76, 0., 538.33, 0., 0., 0., 0., 0., 0., 

  2487.08, 0., 538.33, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., -4145.14, 0., 0., 0., 0., 0., 0., 0., 

  0., -4145.14, 0., -897.217, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 0., 609.082, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -1776.49, 0., 0., 0., 0., -2699.34, 

  0., -322.998, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1065.89,

   0., 0., 0., 0., 6782.96, 0., 538.33, 0., 0., 0., 0., 0., 0., 

  1453.49, 0., 968.994, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -3197.68, 0., 0., 0., 0., 0., 0., 0., 

  0., -2422.49, 0., -1614.99, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 609.082, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -1776.49, 0., 0., 0., 0., -322.998, 

  0., -2699.34, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1065.89,

   0., 0., 0., 0., 538.33, 0., 6782.96, 0., 0., 0., 0., 0., 0., 0., 

  0., 968.994, 0., 1453.49, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -3197.68, 0., 0., 0., 0., 0., 0., 

  0., 0., -1614.99, 0., -2422.49, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 0., 1378.13, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -4073.36, 0., 0., 0., 0., -322.998, 

  0., -4073.36, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2487.08,

   0., 0., 0., 0., 538.33, 0., 9749.76, 0., 0., 0., 0., 0., 0., 0., 

  0., 538.33, 0., 2487.08, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., -4145.14, 0., 0., 0., 0., 0., 0., 0., 

  0., -897.217, 0., -4145.14, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {-76.4236, 

  0., 0., 0., 392.789, 0., 0., 1176.64, 0., 25.9552, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., -272.53, 0., 0., -5165.08, 0., -77.8656, 0.,

   0., 0., 0., -2816.14, 0., -259.552, 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 2725.3, 0., 0., 0., 0., 0., 0., 10356.1, 0., 778.656, 0., 0., 

  0., 0., 0., 0., 0., 0., 1732.08, 0., 302.811, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -3179.51, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -5196.23, 0., -908.432, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 0., 0., 0., 609.082, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -2699.34, 0., 0., 0., 0., 

  0., 0., -1776.49, 0., -322.998, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  1453.49, 0., 0., 0., 0., 0., 0., 6782.96, 0., 968.994, 0., 0., 0., 

  0., 0., 0., 0., 0., 1065.89, 0., 538.33, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., -2422.49, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -3197.68, 0., -1614.99, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {-73.3337, 

  0., 0., 0., 375.732, 0., 0., 375.732, 0., 375.732, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., -259.552, 0., 0., -1594.39, 0., -1594.39, 

  0., 0., 0., 0., -259.552, 0., -1594.39, 0., -259.552, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 778.656, 0., 778.656, 0., 0., 0., 0., 778.656, 0., 

  6184.75, 0., 778.656, 0., 0., 0., 0., 0., 0., 0., 0., 778.656, 0., 

  778.656, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., -2335.97, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -2335.97, 

  0., -2335.97, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 609.082, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -2699.34, 0., 0., 0., 0., 0., 0., -322.998, 0., -1776.49, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 1453.49, 0., 0., 0., 0., 0., 0., 

  968.994, 0., 6782.96, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  538.33, 0., 1065.89, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., -2422.49, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -1614.99, 0., -3197.68, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0.}, {-76.4236, 0., 0., 0., 392.789, 0.,

   0., 25.9552, 0., 1176.64, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -272.53, 0., 0., -77.8656, 0., -5165.08, 0., 0., 0., 0., 0., 

  0., -259.552, 0., -2816.14, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  2725.3, 0., 0., 0., 0., 0., 0., 778.656, 0., 10356.1, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 302.811, 0., 1732.08, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -3179.51, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -908.432, 0., -5196.23, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 

  0., 0., 0., 989.758, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., -532.947, 0., 0., 0., 0., -6928.31, 0., -228.406, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2487.08,

   0., 0., 0., 0., 0., 0., 0., 0., 12471., 0., 1065.89, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -2238.38, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -6532.4, 0., -959.304, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 281.008, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -145.349, 0., 0., 0., 0., -3799.84, 0., -145.349, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1453.49, 0., 0.,

   0., 0., 0., 0., 0., 0., 7897.3, 0., 1453.49, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -1695.74, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -4156.98, 0., -1695.74, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 609.082, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -322.998, 0., 0., 0., 

  0., -1776.49, 0., -2699.34, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 538.33, 0., 968.994, 0., 0., 0., 0., 0., 0.,

   1065.89, 0., 6782.96, 0., 1453.49, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -1614.99, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -3197.68, 0., -2422.49, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0.}, {0., 0., 0., 0., 0., 0., 609.082, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -322.998, 0., 0., 0., 

  0., -2699.34, 0., -1776.49, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 968.994, 0., 538.33, 0., 0., 0., 0., 0., 0.,

   1453.49, 0., 6782.96, 0., 1065.89, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -1614.99, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -2422.49, 0., -3197.68, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 281.008, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -145.349, 0., 0., 0., 0., -145.349, 

  0., -3799.84, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 1453.49, 0., 0., 0., 0., 0., 0., 0., 0., 

  1453.49, 0., 7897.3, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -1695.74, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -1695.74, 0., -4156.98, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 0., 989.758, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -532.947, 0., 0., 0., 0., -228.406, 

  0., -6928.31, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 2487.08, 0., 0., 0., 0., 0., 0., 0., 0., 

  1065.89, 0., 12471., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -2238.38, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -959.304, 0., -6532.4, 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {-80.1888, 0., 0., 0., 27.4933, 0., 0., 2501.89, 0., 27.4933, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -577.359, 0., 

  0., 0., 0., 0., 0., -12413.2, 0., -577.359, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1732.08, 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 19899.6, 0., 1732.08, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., -1270.19, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., -10025.4, 0., -1270.19, 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 0., 0., 0., 989.758, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -228.406, 0., 0., 0., 0., 

  0., 0., -6928.31, 0., -532.947, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 1065.89, 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 12471., 0., 2487.08, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -959.304, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -6532.4, 0., -2238.38, 0., 0., 0., 0., 0.}, {-76.4236, 0., 0., 

  0., 25.9552, 0., 0., 1176.64, 0., 392.789, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -259.552, 0., -77.8656, 0., 0., 0., 

  0., -2816.14, 0., -5165.08, 0., -272.53, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 302.811, 0., 778.656, 0., 0., 0., 0., 

  0., 0., 0., 0., 1732.08, 0., 10356.1, 0., 2725.3, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., -908.432, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., -5196.23, 0., -3179.51, 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., 0., 0., 0., 1378.13, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -322.998, 0., 0., 0., 0., 

  0., 0., -4073.36, 0., -4073.36, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 538.33, 0., 538.33, 0., 0., 0., 0., 0., 

  0., 0., 0., 2487.08, 0., 9749.76, 0., 2487.08, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -897.217, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -4145.14, 0., -4145.14, 0., 0., 0.}, {-76.4236, 

  0., 0., 0., 25.9552, 0., 0., 392.789, 0., 1176.64, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -77.8656, 0., -259.552, 0., 0., 

  0., 0., -272.53, 0., -5165.08, 0., -2816.14, 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 778.656, 0., 302.811, 0., 

  0., 0., 0., 0., 0., 0., 0., 2725.3, 0., 10356.1, 0., 1732.08, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -908.432, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -3179.51, 0., -5196.23, 0., 

  0.}, {0., 0., 0., 0., 0., 0., 0., 0., 989.758, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -228.406, 0., 0., 0., 0., 

  0., 0., -532.947, 0., -6928.31, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 1065.89, 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 2487.08, 0., 12471., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -959.304, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -2238.38, 0., -6532.4, 0.}, {-80.1888, 0., 0., 

  0., 27.4933, 0., 0., 27.4933, 0., 2501.89, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -577.359, 0., 0., 0., 0., 0., 

  0., -577.359, 0., -12413.2, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 1732.08, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 1732.08, 0., 19899.6, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., -1270.19, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., -1270.19, 0., -10025.4}, {0., -109.973, 0., 0., 0., 

  0., 0., 0., 0., 0., 989.758, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -2177.47, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1347.96, 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., -21.9946, 0.,

   0., 0., 0., 0., 0., 0., 0., 461.887, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -1385.66, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 1016.15, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 

  0., -21.9946, 0., 0., 0., 0., 0., 0., 0., 0., 461.887, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., -1385.66, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 1016.15, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., -76.1353, 0., 0., 0., 0., 0., 0., 0., 0., 355.298, 0., 0.,

   228.406, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -319.768, 0., 0., -1065.89, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 959.304, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 182.725, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -852.715, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 767.443, 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0.}, {0., -76.1353, 0., 0., 0., 0., 0., 0., 0., 0.,

   355.298, 0., 0., 0., 0., 228.406, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -319.768, 0., 0., 0., 

  0., -1065.89, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 959.304, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., -48.4497, 0., 0.,

   0., 0., 0., 0., 0., 0., 484.497, 0., 0., 0., 0., 80.7495, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -565.247, 0., 0., 0., 0., -807.495, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 942.078, 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 

  0., 0., -20.7642, 0., 0., 0., 0., 0., 0., 0., 0., 207.642, 0., 0., 

  0., 0., 62.2925, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., -242.249, 0., 0., 0., 0., -622.925, 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 726.746, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0.}, {0., 0., -20.7642, 0., 0., 0., 0., 0., 0., 0., 0., 

  207.642, 0., 0., 0., 0., 0., 0., 62.2925, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -242.249, 0., 0., 0., 

  0., 0., 0., -622.925, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 726.746, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., -48.4497, 0., 0.,

   0., 0., 0., 0., 0., 0., 484.497, 0., 0., 0., 0., 0., 0., 80.7495, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -565.247, 0., 0., 0., 0., 0., 0., -807.495, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 942.078, 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., -48.4497, 0., 0., 0., 0., 0., 0., 0., 0., 80.7495, 0., 0.,

   484.497, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -807.495, 0., 0., 0., 0., 

  0., 0., -565.247, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 942.078, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 258.398, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -430.664, 0., 

  0., 0., 0., 0., 0., -430.664, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 717.773, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0.}, {0., -46.1426, 0., 0., 0., 0., 0.,

   0., 0., 0., 76.9043, 0., 0., 138.428, 0., 138.428, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -230.713, 0., -230.713, 0., 0., 0., 0., 0., 0., -415.283, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 692.139, 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 258.398, 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., -430.664, 0., 0., 0., 0., 0., 0., 0., 

  0., -430.664, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 717.773, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0.}, {0., -48.4497, 0., 0., 0., 0., 0., 0., 0., 0., 80.7495, 

  0., 0., 0., 0., 484.497, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -807.495, 0., 

  0., 0., 0., 0., 0., 0., 0., -565.247, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 942.078, 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0.}, {0., 0., -76.1353, 0., 0., 0., 0., 0.,

   0., 0., 0., 228.406, 0., 0., 0., 0., 355.298, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -1065.89, 0., 0., 0., 0., 0., 0., 0., 0., -319.768, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 959.304, 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., -20.7642,

   0., 0., 0., 0., 0., 0., 0., 0., 62.2925, 0., 0., 0., 0., 207.642, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -622.925, 0., 0., 0., 0., 0., 0., 0., 

  0., -242.249, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 726.746, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., -46.1426, 0., 0., 0., 0., 0., 0., 0., 0., 138.428, 0.,

   0., 0., 0., 76.9043, 0., 138.428, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -230.713, 

  0., -415.283, 0., 0., 0., 0., 0., 0., 0., 0., -230.713, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 692.139, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., -46.1426, 0., 0., 0.,

   0., 0., 0., 0., 0., 138.428, 0., 0., 0., 0., 138.428, 0., 76.9043, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -415.283, 0., -230.713, 0., 0., 0., 0., 0., 0., 

  0., 0., -230.713, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 692.139, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 

  0., -20.7642, 0., 0., 0., 0., 0., 0., 0., 0., 62.2925, 0., 0., 0., 

  0., 0., 0., 207.642, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -622.925, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -242.249, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 726.746, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0.}, {0., 0., 0., -76.1353, 0., 0., 0., 0., 0., 0., 0., 0.,

   228.406, 0., 0., 0., 0., 0., 0., 355.298, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -1065.89, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -319.768, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 959.304, 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0.}, {0., -21.9946, 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 461.887, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -1385.66, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  1016.15, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 182.725, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., -852.715, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 767.443, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., -20.7642, 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 207.642, 0., 62.2925, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -242.249, 0., -622.925, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 726.746, 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 258.398, 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -430.664, 0., -430.664,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 717.773, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., -20.7642, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  62.2925, 0., 207.642, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -622.925, 0., -242.249, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  726.746, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 182.725, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., -852.715, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 767.443, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0.}, {0., -21.9946, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 461.887, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -1385.66, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 1016.15, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., -109.973, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 989.758, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -2177.47, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 1347.96, 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., -21.9946, 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 461.887, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1385.66, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1016.15, 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., -76.1353,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 355.298, 0., 

  228.406, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -319.768, 0., -1065.89, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  959.304, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0.}, {0., 0., 0., -48.4497, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 484.497, 0., 80.7495, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., -565.247, 0., -807.495, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 942.078, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0.}, {0., 0., -48.4497, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 80.7495, 0., 484.497, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -807.495, 0., -565.247, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 942.078, 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., -76.1353, 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 228.406, 0., 355.298, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -1065.89, 0., -319.768, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  959.304, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 

  0., -21.9946, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 461.887, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., -1385.66, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 1016.15, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., -109.973, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 989.758, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., -2177.47, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 1347.96, 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0.}, {29.2116, 0., 0., 0., -1051.62, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 5783.9, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -10025.4, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 5370.76, 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., -329.919, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2969.27, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -6532.4,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 4043.87, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., -329.919, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2969.27, 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -6532.4, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  4043.87, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {27.4933, 0., 0., 

  0., -577.359, 0., 0., -82.4799, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 1732.08, 0., 0., 1732.08, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -1270.19, 0., 0., -5196.23, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 3810.57, 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., -65.9839, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1385.66, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -4156.98, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  3048.46, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0.}, {27.4933, 0., 0., 0., -577.359, 

  0., 0., 0., 0., -82.4799, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  1732.08, 0., 0., 0., 0., 1732.08, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -1270.19, 0., 0., 0., 0., -5196.23, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 3810.57, 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., -532.947, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 2487.08, 0., 0., 0., 0., 888.245, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -2238.38, 0., 0., 

  0., 0., -4145.14, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  3730.63, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., -228.406, 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1065.89, 0., 

  0., 0., 0., 685.217, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., -959.304, 0., 0., 0., 0., -3197.68, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 2877.91, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 

  0., -228.406, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 1065.89, 0., 0., 0., 0., 0., 0., 685.217, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -959.304, 0., 0., 0., 0., 0., 

  0., -3197.68, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  2877.91, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., -532.947, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2487.08, 0., 0., 0., 

  0., 0., 0., 888.245, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -2238.38, 0., 0., 0., 0., 0., 0., -4145.14, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 3730.63, 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {27.253, 0., 0., 

  0., -272.53, 0., 0., -272.53, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 317.951, 0., 0., 2725.3, 0., 0., 0., 0., 0., 0., 

  317.951, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -3179.51, 0., 

  0., 0., 0., 0., 0., -3179.51, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 3709.43, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., -145.349, 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1453.49, 

  0., 0., 0., 0., 0., 0., 242.249, 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., -1695.74, 0., 0., 0., 0., 0., 0., -2422.49, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 2826.23, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {25.9552, 0., 0., 

  0., -259.552, 0., 0., -77.8656, 0., -77.8656, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 302.811, 0., 0., 778.656, 0., 778.656, 0., 0., 

  0., 0., 0., 0., 233.597, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -908.432, 0., -908.432, 0., 0., 0., 0., 0., 0., -2335.97, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2725.3, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 

  0., 0., 0., -145.349, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 1453.49, 0., 0., 0., 0., 0., 0., 0., 0., 242.249, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., -1695.74, 0., 0., 0., 0., 0., 

  0., 0., 0., -2422.49, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 2826.23, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {27.253, 0., 0., 0., -272.53, 0., 0., 0., 0., -272.53, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 317.951, 0., 0., 0., 0., 2725.3, 0.,

   0., 0., 0., 0., 0., 0., 0., 317.951, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -3179.51, 0., 0., 0., 0., 0., 0., 0., 0., -3179.51, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 3709.43, 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., -532.947, 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 888.245, 

  0., 0., 0., 0., 2487.08, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., -4145.14, 0., 0., 0., 0., 0., 0., 0., 

  0., -2238.38, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  3730.63, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 

  0., 0., 0., 0., -145.349, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 242.249, 0., 0., 0., 0., 1453.49, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -2422.49, 0., 0., 0., 0., 0., 0., 0., 0., -1695.74, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 2826.23, 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., -322.998, 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 538.33, 0., 0., 0., 0.,

   538.33, 0., 968.994, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -897.217, 0., -1614.99, 0., 0., 0., 0., 0., 0., 

  0., 0., -1614.99, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  2691.65, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 

  0., 0., -322.998, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 538.33, 0., 0., 0., 0., 968.994, 0., 538.33, 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1614.99, 

  0., -897.217, 0., 0., 0., 0., 0., 0., 0., 0., -1614.99, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 2691.65, 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0.}, {0., 0., 0., 0., 0., -145.349, 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 242.249, 0., 0., 0., 0., 0., 

  0., 1453.49, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., -2422.49, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -1695.74, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  2826.23, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 

  0., -532.947, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 888.245, 0., 0., 0., 0., 0., 0., 2487.08, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -4145.14, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -2238.38, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 3730.63, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {27.4933, 0., 0., 0., -82.4799, 0., 0., -577.359, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1732.08, 0., 0., 0.,

   0., 0., 0., 1732.08, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., -5196.23, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -1270.19, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 3810.57, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 

  0., -228.406, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 685.217, 0., 0., 0., 0., 0., 0., 1065.89, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -3197.68, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -959.304, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 2877.91, 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {25.9552, 0., 0., 0., -77.8656, 0., 0., -259.552, 0., -77.8656,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 778.656, 0., 

  233.597, 0., 0., 0., 0., 302.811, 0., 778.656, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -908.432, 0., -2335.97, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -908.432, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 2725.3, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 

  0., 0., 0., 0., 0., 0., -322.998, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 968.994, 0., 0., 0., 0., 0., 0., 538.33,

   0., 538.33, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., -1614.99, 0., -1614.99, 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., -897.217, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  2691.65, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0.}, {25.9552, 0., 0., 0., -77.8656, 0., 

  0., -77.8656, 0., -259.552, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 233.597, 0., 778.656, 0., 0., 0., 0., 0., 0., 778.656, 

  0., 302.811, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., -2335.97, 0., -908.432, 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., -908.432, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  2725.3, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., -228.406, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 685.217, 0., 0.,

   0., 0., 0., 0., 0., 0., 1065.89, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -3197.68, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -959.304, 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 2877.91, 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0.}, {27.4933, 0., 0., 

  0., -82.4799, 0., 0., 0., 0., -577.359, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 1732.08, 0., 0., 0., 0., 0., 0., 0.,

   0., 1732.08, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., -5196.23, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -1270.19, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 3810.57, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0.}, {0., 0., 0., 0., 0., -329.919, 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2969.27, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -6532.4, 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 4043.87, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., -65.9839, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  1385.66, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -4156.98, 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 3048.46, 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., -228.406, 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   1065.89, 0., 685.217, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -959.304, 

  0., -3197.68, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2877.91, 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 

  0., -145.349, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 1453.49, 0., 242.249, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., -1695.74, 0., -2422.49, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 2826.23, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., -145.349, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 242.249, 0., 

  1453.49, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -2422.49, 

  0., -1695.74, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2826.23, 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 

  0., -228.406, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 685.217, 0., 1065.89, 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., -3197.68, 0., -959.304, 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 2877.91, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0.}, {0., 0., 0., 0., 0., -65.9839, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1385.66, 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -4156.98, 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 3048.46, 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0.}, {0., 0., 0., 0., 0., 0., -329.919, 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  2969.27, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -6532.4, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 4043.87, 0., 0., 0., 0., 0.,

   0., 0., 0., 0.}, {29.2116, 0., 0., 0., 0., 0., 0., -1051.62, 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 5783.9, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., -10025.4, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 5370.76, 0.,

   0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 

  0., -329.919, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 2969.27, 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., -6532.4, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 4043.87, 0., 0., 0., 0., 0., 0., 0.}, {27.4933, 0., 0., 0., 

  0., 0., 0., -577.359, 0., -82.4799, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1732.08, 0., 

  1732.08, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1270.19, 

  0., -5196.23, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 3810.57, 0.,

   0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., -532.947, 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 2487.08, 0., 888.245, 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -2238.38, 0., -4145.14, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 3730.63, 0., 0., 0., 0., 0.}, {27.253, 0., 0., 0., 0., 0., 

  0., -272.53, 0., -272.53, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 317.951, 0., 2725.3, 0., 

  317.951, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -3179.51, 

  0., -3179.51, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 3709.43, 0.,

   0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., -532.947, 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 888.245, 0., 2487.08, 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., -4145.14, 0., -2238.38, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 3730.63, 0., 0., 0.}, {27.4933, 0., 0., 0., 0., 0., 

  0., -82.4799, 0., -577.359, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1732.08, 0., 

  1732.08, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -5196.23, 

  0., -1270.19, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 3810.57, 0.,

   0.}, {0., 0., 0., 0., 0., 0., 0., 0., -329.919, 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 2969.27, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., -6532.4, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  4043.87, 0.}, {29.2116, 0., 0., 0., 0., 0., 0., 0., 0., -1051.62, 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 5783.9, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., -10025.4, 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 5370.76}
    
};



 

 

 

 

 std::vector<std::vector<double>> EquivalentPolynomial::_ATA_deg3_prism = {{111.25, -615., 0., 0., 1050., 0., 0., 0., 0., -33.75, -560., 0., 0.,

   0., 0., 45., 0., 0., 0., 0.}, {-615., 5430., -2715., 0., -11340., 

  7560., 0., 0., 0., 45., 6720., -5040., 0., 0., 0., -90., 0., 0., 

  45., 0.}, {0., -2715., 5430., 0., 7560., -15120., 0., 0., 0., 

  0., -5040., 10080., 0., 0., 0., 45., 0., 0., -90., 0.}, {0., 0., 0.,

   123.75, 0., 0., -360., 0., 0., 0., 0., 0., 270., 0., 0., 0., 0., 

  0., 0., -26.25}, {1050., -11340., 7560., 0., 26460., -26460., 0., 

  4410., 0., 0., -16800., 20160., 0., -5040., 0., 0., 0., 0., 0., 

  0.}, {0., 7560., -15120., 0., -26460., 70560., 0., -26460., 0., 0., 

  20160., -60480., 0., 30240., 0., 0., 0., 0., 0., 0.}, {0., 0., 

  0., -360., 0., 0., 1800., 0., -900., 0., 0., 0., -1620., 0., 1080., 

  0., 0., 0., 0., 0.}, {0., 0., 0., 0., 4410., -26460., 0., 26460., 

  0., 0., -5040., 30240., 0., -30240., 0., 0., 0., 0., 0., 0.}, {0., 

  0., 0., 0., 0., 0., -900., 0., 1800., 0., 0., 0., 1080., 0., -2160.,

   0., 0., 0., 0., 0.}, {-33.75, 45., 0., 0., 0., 0., 0., 0., 0., 

  101.25, 0., 0., 0., 0., 0., -135., 0., 0., 0., 0.}, {-560., 

  6720., -5040., 0., -16800., 20160., 0., -5040., 0., 0., 

  11200., -16800., 0., 6720., 0., 0., -560., 0., 0., 0.}, {0., -5040.,

   10080., 0., 20160., -60480., 0., 30240., 0., 0., -16800., 60480., 

  0., -45360., 0., 0., 6720., 0., 0., 0.}, {0., 0., 0., 270., 0., 

  0., -1620., 0., 1080., 0., 0., 0., 1620., 0., -1620., 0., 0., 270., 

  0., 0.}, {0., 0., 0., 0., -5040., 30240., 0., -30240., 0., 0., 

  6720., -45360., 0., 60480., 0., 0., -16800., 0., 0., 0.}, {0., 0., 

  0., 0., 0., 0., 1080., 0., -2160., 0., 0., 0., -1620., 0., 4320., 

  0., 0., -1620., 0., 0.}, {45., -90., 45., 0., 0., 0., 0., 0., 

  0., -135., 0., 0., 0., 0., 0., 270., 0., 0., -135., 0.}, {0., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., -560., 6720., 0., -16800., 0., 0., 

  11200., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 270., 0., -1620., 0., 0., 1620., 0., 0.}, {0., 45., -90., 0., 

  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -135., 0., 0., 270., 

  0.}, {0., 0., 0., -26.25, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 

  0., 0., 0., 0., 0., 43.75}
     
};







 

 

 

 

 





/*







#include "eqPoly.hpp"
#include <time.h>
#include "FemusInit.hpp"
#include "Elem.hpp"
#include <stdlib.h>



using namespace femus;
using namespace std;
using namespace Eigen;
using std::cout;

int main(int argc, char** args) {

  //FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  //SetCoefficients takes ( dim, degree of equivalent polynomial, rho, vector < a, b, c > for the dicontinuity ( point, line, or plane ), and element (3=triangle/tet, 4=square/cube
  // a*x + b*y + c*z + d = 0, and d) as inputs
  EquivalentPolynomial eqP;

  //
   eqP.SetCoefficients(3, 2, 2, std::vector<double> {1., 2., 1.}, 0., 3);
   eqP.PrintCoefficients();
//   std::cout << eqP.GetValue(std::vector<double> {0.5, 0.5}) << " " << std::endl;


  std::vector < double >points {1.,2.,4.,5.,7.,8.,-5.,23.,12.,15.,14.,14.};
  //std::vector < double >points(2000);
  std::vector < double >normal {1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0.};
  unsigned dim = 3;
  eqP.FindBestFit(points, normal, dim);
  std::vector < double > onepoint {1.,2.,3.};
  unsigned element = 0;
  std::cout << eqP.GetValue(onepoint, element) << "  value" << endl;

//   clock_t t;
//   t = clock();
//
//   for(unsigned k = 0; k < 1000; k++) {
//
//     for(unsigned i = 0; i < 2000; i++) {
//
//       points[i] = rand() % 100;
//
//     }
//
//
//     eqP.FindBestFit(points, normal, dim);
//
//
//   }
//
//
//
//   t = clock() - t;
//   printf("It took me %d clicks (%f seconds).\n", t, ((float)t) / CLOCKS_PER_SEC);

*/
/*std::vector <double> tes(27);
MatrixXd C;
C.setRandom(10,3);
JacobiSVD<MatrixXd> svd( C, ComputeThinU | ComputeThinV);
MatrixXd Cp = svd.matrixV();
MatrixXd sigma = svd.singularValues().asDiagonal();
svd.matrixV().transpose();
MatrixXd diff = Cp;
tes[0] = Cp(1,2);

cout << "diff:\n" << Cp.col(2) << "\n";
cout << "diff:\n" << tes[0] << "\n";
*/

/*

  std::vector<double> phi;
  std::vector<double> gradPhi;
  double weight;
  {
    std::vector<std::vector<double>> xv = {{-1., 1.}};
    double integral = 0.;
    unsigned dim = 1;
    element = 6;
    eqP.SetCoefficients(1, 3, 50, std::vector<double> {1.}, -0.5, 4);
    const elem_type * fe = new const elem_type_1D("line", "linear", "ninth");
    for(unsigned ig = 0; ig < fe->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      fe->Jacobian(xv, ig, weight, phi, gradPhi);
      std::vector<double> xg(dim, 0.);
      for(unsigned i = 0; i < phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xg[k] += xv[k][i] * phi[i];
        }
      }
      integral += xg[0] * xg[0] * xg[0] * eqP.GetValue(xg, element) * weight;
    }
    std::cout << "Integral = " << integral << std::endl;
    delete fe;
  }

  element = 4;
  std::vector  <double> pt {0.5, 0.5};
  eqP.SetCoefficients(2, 2, 10, std::vector<double> {2., 1.}, -1., element);
  std::cout << eqP.GetValue(pt, element) << " triangle at 0.5, 0.5 " << std::endl;
     eqP.PrintCoefficients();

     //TODO Fix triangle integration

  {
    std::vector<std::vector<double>> xv = {{-1., 1., 1., -1.}, {-1., -1., 1., 1.}};
    double integral = 0.;
    unsigned dim = 2;
    element = 4;
    eqP.SetCoefficients(2, 2, 10, std::vector<double> {2., 1.}, -1., element);
    const elem_type * fe = new const elem_type_2D("tri", "linear", "ninth");
    for(unsigned ig = 0; ig < fe->GetGaussPointNumber(); ig++) {

      fe->Jacobian(xv, ig, weight, phi, gradPhi);
      std::vector<double> xg(dim, 0.);
      for(unsigned i = 0; i < phi.size(); i++) {
        for(unsigned k = 0; k < dim; k++) {
          xg[k] += xv[k][i] * phi[i];
        }
      }
      integral += xg[0] * xg[1] * eqP.GetValue(xg, element) * weight;
    }
    std::cout << "Integral = " << integral << std::endl;
    delete fe;
  }





  return 0;
}

*/

