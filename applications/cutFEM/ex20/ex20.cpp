
#include <algorithm>
#include <vector>
#include <cmath>
#include <algorithm>    // std::sort

int main() {
  return 1;
}


double b(const double &p) {
  return std::max(std::min(p, 1.), 0.);
}



void GetInterval4(const std::vector <double> &a1, const std::vector <double> &a2, std::vector< std::pair<double, double> > &I1, std::vector< std::pair<double, double> > &I2, std::vector<std::pair<double, double>> &I3) {

  std::vector <double> x(6);

  x[0] = 0 ;
  unsigned cnt = 1;

  for(unsigned k = 0; k < 2; k++) {
    const std::vector<double> &a = (k == 0) ? a1 : a2;
    double delta = a[1] * a[1] - 4 * a[0] * a[2] ;
    if(delta > 0) {
      for(unsigned i = 0; i < 2; i++) {
        double y = (-a[1] - pow(-1, i) * sqrt(delta)) / (2 * a[0]) ;
        if(y > 0 && y < 1) {
          x[cnt] = y ;
          cnt++ ;
        }
      }
    }
  }

  x[cnt] = 1 ;
  cnt++;

  x.resize(cnt);
  std::sort(x.begin(), x.end());
  I1.resize(0);
  I2.resize(0);
  I3.resize(0);

  for(unsigned i = 0 ; i < cnt - 1 ; i++) {
    double xm = (x[i] + x[i + 1]) * 0.5 ;
    double f1 = a1[0] * xm * xm + a1[1] * xm + a1[2] ;
    double f2 = a2[0] * xm * xm + a2[1] * xm + a2[2] ;
    if(f1 > 0) {
      if(f2 > 0) {
        I3.resize(I3.size() + 1, std::pair<double, double>(x[i], x[i + 1]));
      }
      else {
        I1.resize(I1.size() + 1, std::pair<double, double>(x[i], x[i + 1]));
      }
    }
    else if(f2 > 0) {
      I2.resize(I2.size() + 1, std::pair<double, double>(x[i], x[i + 1]));
    }
  }

  for(unsigned i = 0; i < I1.size() - 1; i++) {
    if(I1[i].second == I1[i + 1].first) {
      I1[i].second = I1[i + 1].second;
      I1.erase(I1.begin() + i + 1);
    }
  }
  
  for(unsigned i = 0; i < I2.size() - 1; i++) {
    if(I2[i].second == I2[i + 1].first) {
      I2[i].second = I2[i + 1].second;
      I2.erase(I2.begin() + i + 1);
    }
  }
  
  for(unsigned i = 0; i < I3.size() - 1; i++) {
    if(I3[i].second == I3[i + 1].first) {
      I3[i].second = I3[i + 1].second;
      I3.erase(I3.begin() + i + 1);
    }
  }


}



void GetInterval4(double p1, double p2, double q1, double q2, const double &k, std::vector< std::pair<double, double> > &I1, std::vector< std::pair<double, double> > &I2, std::vector<std::pair<double, double>> &I3) {

  bool swap = false;
  if(p1 > q1) {
    std::swap(p1, q1);
    std::swap(p2, q2);
    swap = true ;
  }


  if(k > 0) {
    if(p2 <= q1) { //case 1;
      I1.resize(1);
      I1[0].first = b(q1);
      I1[0].second = b(q2);

      I2.resize(1);
      I2[0].first = b(p1);
      I2[0].second = b(p2);


      I3.resize(3);
      I3[0].first = 0;
      I3[0].second = b(p1);
      I3[1].first = b(p2);
      I3[1].second = b(q1);
      I3[2].first = b(q2);
      I3[2].second = 1;

    }
    else if(p2 <= q2) {// case 2

      I1.resize(1);
      I1[0].first = b(p2);
      I1[0].second = b(q2);

      I2.resize(1);
      I2[0].first = b(p1);
      I2[0].second = b(q1);


      I3.resize(2);
      I3[0].first = 0;
      I3[0].second = b(p1);
      I3[1].first = b(q2);
      I3[1].second = 1;

    }
    else {// case 3

      I2.resize(2);
      I2[0].first = b(p1);
      I2[0].second = b(q1);
      I2[1].first = b(q2);
      I2[1].second = b(p2);


      I3.resize(2);
      I3[0].first = 0;
      I3[0].second = b(p1);
      I3[1].first = b(p1);
      I3[1].second = 1;


    }
  }
  else { // (k<0)
    if(p2 <= q1) { //case 1
      I1.resize(1);
      I1[0].first = b(p1);
      I1[0].second = b(p2);

      I2.resize(1);
      I2[0].first = b(q1);
      I2[0].second = b(q2);

      // I3 =empty
    }
    else if(p2 <= q2) {// case 2
      I1.resize(1);
      I1[0].first = b(p1);
      I1[0].second = b(q1);
      I2.resize(1);
      I2[0].first = b(p2);
      I2[0].second = b(q2);
      I3.resize(1);
      I3[0].first = b(q1);
      I3[0].second = b(p2);

    }
    else { // case 3
      I1.resize(2);
      I1[0].first = b(p1);
      I1[0].second = b(q1);
      I1[1].first = b(q2);
      I1[1].second = b(p2);
      // I2= empty
      I3.resize(1);
      I3[0].first = b(q1);
      I3[0].second = b(q2);
    }
  }

  if(swap) {
    I1.swap(I2);
  }

}
void GetInterval2(const double &r1, const double &r2, const bool &pIsComplex, const double &k, std::vector< std::pair<double, double> > &I1, std::vector< std::pair<double, double> > &I2, std::vector<std::pair<double, double>> &I3) {
  if(k > 0) {
    // I1 = empty
    I2.resize(1);
    I2[0].first = b(r1);
    I2[0].second = b(r2);
    I3.resize(2);
    I3[0].first = 0;
    I3[0].second = b(r1);
    I3[1].first = b(r2);
    I3[1].second = b(1);
    if(pIsComplex) I1.swap(I2); // swap I1 and I2 if p-roots are complex

  }

  else { //(k<0)
    I1.resize(1);
    I1[0].first = b(r1);
    I1[0].second = b(r2);
    // I2 = I3 = empty
    if(pIsComplex) I1.swap(I2); // swap I1 and I2 if p-roots are complex

  }

}

void GetInterval0(const double &k, std::vector< std::pair<double, double> > &I1, std::vector< std::pair<double, double> > &I2, std::vector<std::pair<double, double>> &I3) {
  if(k > 0) I3.assign(1, std::pair<double, double>(0., 1.));
}


void GetInterval4Old(double p1, double p2, double q1, double q2, const double &k, std::vector< std::pair<double, double> > &I1, std::vector< std::pair<double, double> > &I2, std::vector<std::pair<double, double>> &I3) {

  if(k > 0) {
    if(p2 <= q1) { //case 1;
      I1.resize(1);
      I1[0].first = b(q1);
      I1[0].second = b(q2);

      I2.resize(1);
      I2[0].first = b(p1);
      I2[0].second = b(p2);


      I3.resize(3);
      I3[0].first = 0;
      I3[0].second = b(p1);
      I3[1].first = b(p2);
      I3[1].second = b(q1);
      I3[2].first = b(q2);
      I3[2].second = 1;

    }
    else if(q2 <= p1) { // case 4
      I1.resize(1);
      I1[0].first = b(q1);
      I1[0].second = b(q2);

      I2.resize(1);
      I2[0].first = b(p1);
      I2[0].second = b(p2);


      I3.resize(3);
      I3[0].first = 0;
      I3[0].second = b(q1);
      I3[1].first = b(q2);
      I3[1].second = b(p1);
      I3[2].first = b(p2);
      I3[2].second = 1;


    }
    else if(p2 <= q2) {
      if(p1 <= q1) {// case 2

        I1.resize(1);
        I1[0].first = b(p2);
        I1[0].second = b(q2);

        I2.resize(1);
        I2[0].first = b(p1);
        I2[0].second = b(q1);


        I3.resize(2);
        I3[0].first = 0;
        I3[0].second = b(p1);
        I3[1].first = b(q2);
        I3[1].second = 1;
      }
      else { // case 6

        I1.resize(2);
        I1[0].first = b(q1);
        I1[0].second = b(p1);
        I1[1].first = b(p2);
        I1[1].second = b(q2);

        // I2= empty


        I3.resize(2);
        I3[0].first = 0;
        I3[0].second = b(q1);
        I3[1].first = b(q2);
        I3[1].second = 1;
      }
    }
    else {
      if(p1 <= q1) { // case 3


        I2.resize(2);
        I2[0].first = b(p1);
        I2[0].second = b(q1);
        I2[1].first = b(q2);
        I2[1].second = b(p2);


        I3.resize(2);
        I3[0].first = 0;
        I3[0].second = b(p1);
        I3[1].first = b(p1);
        I3[1].second = 1;


      }
      else { //case 5;
        I1.resize(1);
        I1[0].first = b(q1);
        I1[0].second = b(p1);

        I2.resize(1);
        I2[0].first = b(q2);
        I2[0].second = b(p2);


        I3.resize(2);
        I3[0].first = 0;
        I3[0].second = b(q1);
        I3[1].first = b(p2);
        I3[1].second = 1;


      }
    }
  }
  else { // (k<0)
    if(p2 <= q1 || q2 <= p1) { //case 1 and 4
      I1.resize(1);
      I1[0].first = b(p1);
      I1[0].second = b(p2);

      I2.resize(1);
      I2[0].first = b(q1);
      I2[0].second = b(q2);

      // I3 =empty
    }
    else if(p2 <= q2) {
      if(p1 <= q1) {// case 2
        I1.resize(1);
        I1[0].first = b(p1);
        I1[0].second = b(q1);
        I2.resize(1);
        I2[0].first = b(p2);
        I2[0].second = b(q2);
        I3.resize(1);
        I3[0].first = b(q1);
        I3[0].second = b(p2);
      }
      else { // case 6
        // I1 is empty

        I2.resize(2);
        I2[0].first = b(q1);
        I2[0].second = b(p1);
        I2[1].first = b(p2);
        I2[1].second = b(q2);
        // I2= empty
        I3.resize(1);
        I3[0].first = b(p1);
        I3[0].second = b(p2);

      }
    }
    else {
      if(p1 <= q1) { // case 3
        I1.resize(2);
        I1[0].first = b(p1);
        I1[0].second = b(q1);
        I1[1].first = b(q2);
        I1[1].second = b(p2);
        // I2= empty
        I3.resize(1);
        I3[0].first = b(q1);
        I3[0].second = b(q2);
      }
      else { //case 5;
        I1.resize(1);
        I1[0].first = b(q2);
        I1[0].second = b(p2);
        I2.resize(1);
        I2[0].first = b(q1);
        I2[0].second = b(p1);
        I3.resize(2);
        I3[0].first = b(p1);
        I3[0].second = b(q2);
      }
    }
  }
}
