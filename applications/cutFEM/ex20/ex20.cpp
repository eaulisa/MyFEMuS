
#include <algorithm>
#include <vector>

int main() {
  return 1;
}


double b(const double &p) {
  return std::max(std::min(p, 1.), 0.);
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
