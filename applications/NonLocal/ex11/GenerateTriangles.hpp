
void GenerateRandomTriangles(const unsigned seed, std::vector<std::vector<double> > &x) {


    srand(seed);
    x.assign(2, std::vector<double>(3));

restart:

    double t0 = 2. * M_PI * rand() / RAND_MAX;
    double h0 = 0.5 * rand() / RAND_MAX;
    x[0][0] = 1. + h0 * cos(t0);
    x[1][0] = h0 * sin(t0);

    double dt1 = M_PI / 8. + 3. / 4. * M_PI * rand() / RAND_MAX;
    double h1 = h0 * (0.5 + rand() / RAND_MAX);
    x[0][1] = 1. + h1 * cos(t0 + dt1);
    x[1][1] = h1 * sin(t0 + dt1);

    x[0][2] = 3. - x[0][0] - x[0][1];
    x[1][2] = 0. - x[1][0] - x[1][1];

    double det = (x[0][1] - x[0][0]) * (x[1][2] - x[1][0]) - (x[1][1] - x[1][0]) * (x[0][2] - x[0][0]);

    if(det > 0) {

      bool plus = false;
      bool minus = false;
      for(unsigned j = 0; j < x.size(); j++) {
        double d2 = 1 - x[0][j] * x[0][j] - x[1][j] * x[1][j];
        if(d2 < 0) minus = true;
        else if(d2 > 0) plus = true;
      }
      if(plus*minus == true){
        return;
      }
      else{
        goto restart;
      }
    }
    
    else {
      goto restart;
    }
  

}
