#include <iostream>
#include <vector>
#include <math.h>
#include <ctime>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include <eigen3/unsupported/Eigen//CXX11/Tensor>
#include <fstream>


// valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all ./executable
//

//first row-weights, second row: x-coordinates
const std::vector<std::vector < double > > Gauss1 = { {2}, {0} };

const std::vector<std::vector < double > > Gauss2 = {{1, 1},
  { -0.57735026918963, 0.57735026918963}
};

const std::vector<std::vector < double > > Gauss3 = {{0.55555555555556, 0.88888888888889, 0.55555555555556},
  { -0.77459666924148, 0, 0.77459666924148}
};

const std::vector<std::vector < double > > Gauss4 = {{0.34785484513745, 0.65214515486255, 0.65214515486255, 0.34785484513745},
  { -0.86113631159405, -0.33998104358486, 0.33998104358486, 0.86113631159405}
};

const std::vector<std::vector < double > > Gauss5 = {{0.23692688505619, 0.47862867049937, 0.56888888888889, 0.47862867049937, 0.23692688505619},
  { -0.90617984593866, -0.53846931010568, 0, 0.53846931010568, 0.90617984593866}
};

const std::vector<std::vector < double > > Gauss6 = {{0.3607615730481388, 0.3607615730481388, 0.4679139345726911, 0.4679139345726911, 0.17132449237917097, 0.17132449237917097},
  {0.6612093864662645, -0.6612093864662645, -0.23861918608319715, 0.23861918608319715, -0.9324695142031519, 0.9324695142031519}
};

const std::vector<std::vector<double>> Gauss7 = {{0.4179591836734694, 0.3818300505051188, 0.3818300505051188, 0.27970539148927687, 0.27970539148927703, 0.12948496616886915, 0.12948496616886915},
  {0., 0.4058451513773972, -0.4058451513773972, -0.7415311855993942, 0.7415311855993945, -0.9491079123427586, 0.9491079123427586}
};

const std::vector<std::vector < double > > Gauss8 = {{0.362683783378362, 0.313706645877887, 0.222381034453375, 0.101228536290376, 0.362683783378362, 0.313706645877887, 0.222381034453375, 0.101228536290376}, { -0.183434642495650, -0.525532409916329, -0.796666477413627, -0.960289856497536, 0.183434642495650, 0.525532409916329, 0.796666477413627, 0.960289856497536}
};

const std::vector<std::vector < double > > Gauss9 = {{0.330239355001260, 0.312347077040003, 0.260610696402935, 0.180648160694857, 0.081274388361574, 0.312347077040003, 0.260610696402935, 0.180648160694857, 0.081274388361574}, {0.000000000000000, -0.324253423403809, -0.613371432700590, -0.836031107326636, -0.968160239507626, 0.324253423403809, 0.613371432700590, 0.836031107326636, 0.968160239507626}
};

const std::vector<std::vector < double > > Gauss10 = {{0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688, 0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688}, { -0.148874338981631, -0.433395394129247, -0.679409568299024, -0.865063366688985, -0.973906528517172, 0.148874338981631, 0.433395394129247, 0.679409568299024, 0.865063366688985, 0.973906528517172}
};

const std::vector<std::vector < double > > Gauss11 = {{0.272925086777901, 0.262804544510247, 0.233193764591990, 0.186290210927734, 0.125580369464905, 0.055668567116174, 0.262804544510247, 0.233193764591990, 0.186290210927734, 0.125580369464905, 0.055668567116174}, {0.000000000000000, -0.269543155952345, -0.519096129206812, -0.730152005574049, -0.887062599768095, -0.978228658146057, 0.269543155952345, 0.519096129206812, 0.730152005574049, 0.887062599768095, 0.978228658146057}
};

const std::vector<std::vector < double > > Gauss12 = {{0.249147045813403, 0.233492536538355, 0.203167426723066, 0.160078328543346, 0.106939325995318, 0.047175336386512, 0.249147045813403, 0.233492536538355, 0.203167426723066, 0.160078328543346, 0.106939325995318, 0.047175336386512}, {0.125233408511469, 0.367831498998180, 0.587317954286617, 0.769902674194305, 0.904117256370475, 0.981560634246719, -0.125233408511469, -0.367831498998180, -0.587317954286617, -0.769902674194305, -0.904117256370475, -0.981560634246719}
};

const std::vector<std::vector < double > > Gauss13 = {{0.232551553230874, 0.226283180262897, 0.207816047536889, 0.178145980761946, 0.138873510219787, 0.092121499837728, 0.040484004765316, 0.226283180262897, 0.207816047536889, 0.178145980761946, 0.138873510219787, 0.092121499837728, 0.040484004765316}, {0.000000000000000, 0.230458315955135, 0.448492751036447, 0.642349339440340, 0.801578090733310, 0.917598399222978, 0.984183054718588, -0.230458315955135, -0.448492751036447, -0.642349339440340, -0.801578090733310, -0.917598399222978, -0.984183054718588}
};

const std::vector<std::vector < double > > Gauss14 = {{0.215263853463158, 0.205198463721296, 0.185538397477938, 0.157203167158194, 0.121518570687903, 0.080158087159760, 0.035119460331752, 0.215263853463158, 0.205198463721296, 0.185538397477938, 0.157203167158194, 0.121518570687903, 0.080158087159760, 0.035119460331752}, {0.108054948707344, 0.319112368927890, 0.515248636358154, 0.687292904811685, 0.827201315069765, 0.928434883663574, 0.986283808696812, -0.108054948707344, -0.319112368927890, -0.515248636358154, -0.687292904811685, -0.827201315069765, -0.928434883663574, -0.986283808696812}
};

const std::vector<std::vector < double > > Gauss15 = {{0.202578241925561, 0.198431485327112, 0.186161000015562, 0.166269205816994, 0.139570677926154, 0.107159220467172, 0.070366047488108, 0.030753241996117, 0.198431485327112, 0.186161000015562, 0.166269205816994, 0.139570677926154, 0.107159220467172, 0.070366047488108, 0.030753241996117}, {0.000000000000000, 0.201194093997435, 0.394151347077563, 0.570972172608539, 0.724417731360170, 0.848206583410427, 0.937273392400706, 0.987992518020485, -0.201194093997435, -0.394151347077563, -0.570972172608539, -0.724417731360170, -0.848206583410427, -0.937273392400706, -0.987992518020485}
};

const std::vector<std::vector < double > > Gauss16 = {{0.189450610455068, 0.182603415044922, 0.169156519395003, 0.149595988816577, 0.124628971255534, 0.095158511682492, 0.062253523938648, 0.027152459411754, 0.189450610455068, 0.182603415044922, 0.169156519395003, 0.149595988816577, 0.124628971255534, 0.095158511682492, 0.062253523938648, 0.027152459411754}, {0.095012509837637, 0.281603550779259, 0.458016777657227, 0.617876244402644, 0.755404408355003, 0.865631202387832, 0.944575023073233, 0.989400934991650, -0.095012509837637, -0.281603550779259, -0.458016777657227, -0.617876244402644, -0.755404408355003, -0.865631202387832, -0.944575023073233, -0.989400934991650}
};

const std::vector<std::vector < double > > Gauss17 = {{0.179446470356207, 0.176562705366993, 0.168004102156450, 0.154045761076810, 0.135136368468525, 0.111883847193404, 0.085036148317179, 0.055459529373987, 0.024148302868548, 0.176562705366993, 0.168004102156450, 0.154045761076810, 0.135136368468525, 0.111883847193404, 0.085036148317179, 0.055459529373987, 0.024148302868548}, {0.000000000000000, 0.178484181495848, 0.351231763453876, 0.512690537086477, 0.657671159216691, 0.781514003896801, 0.880239153726986, 0.950675521768768, 0.990575475314417, -0.178484181495848, -0.351231763453876, -0.512690537086477, -0.657671159216691, -0.781514003896801, -0.880239153726986, -0.950675521768768, -0.990575475314417}
};

const std::vector<std::vector < double > > Gauss18 = {{0.169142382963144, 0.164276483745833, 0.154684675126265, 0.140642914670651, 0.122555206711479, 0.100942044106287, 0.076425730254889, 0.049714548894970, 0.021616013526483, 0.169142382963144, 0.164276483745833, 0.154684675126265, 0.140642914670651, 0.122555206711479, 0.100942044106287, 0.076425730254889, 0.049714548894970, 0.021616013526483}, { -0.084775013041735, -0.251886225691505, -0.411751161462843, -0.559770831073948, -0.691687043060353, -0.803704958972523, -0.892602466497556, -0.955823949571398, -0.991565168420931, 0.084775013041735, 0.251886225691505, 0.411751161462843, 0.559770831073948, 0.691687043060353, 0.803704958972523, 0.892602466497556, 0.955823949571398, 0.991565168420931}
};

const std::vector<std::vector < double > > Gauss19 = {{0.161054449848784, 0.158968843393954, 0.152766042065860, 0.142606702173607, 0.128753962539336, 0.111566645547334, 0.091490021622450, 0.069044542737641, 0.044814226765700, 0.019461788229727, 0.158968843393954, 0.152766042065860, 0.142606702173607, 0.128753962539336, 0.111566645547334, 0.091490021622450, 0.069044542737641, 0.044814226765700, 0.019461788229727}, {0.000000000000000, -0.160358645640225, -0.316564099963630, -0.464570741375961, -0.600545304661681, -0.720966177335229, -0.822714656537143, -0.903155903614818, -0.960208152134830, -0.992406843843584, 0.160358645640225, 0.316564099963630, 0.464570741375961, 0.600545304661681, 0.720966177335229, 0.822714656537143, 0.903155903614818, 0.960208152134830, 0.992406843843584}
};

const std::vector<std::vector < double > > Gauss20 = {{0.152753387130726, 0.149172986472604, 0.142096109318382, 0.131688638449177, 0.118194531961518, 0.101930119817240, 0.083276741576705, 0.062672048334109, 0.040601429800387, 0.017614007139152, 0.152753387130726, 0.149172986472604, 0.142096109318382, 0.131688638449177, 0.118194531961518, 0.101930119817240, 0.083276741576705, 0.062672048334109, 0.040601429800387, 0.017614007139152}, {0.076526521133497, 0.227785851141645, 0.373706088715420, 0.510867001950827, 0.636053680726515, 0.746331906460151, 0.839116971822219, 0.912234428251326, 0.963971927277914, 0.993128599185095, -0.076526521133497, -0.227785851141645, -0.373706088715420, -0.510867001950827, -0.636053680726515, -0.746331906460151, -0.839116971822219, -0.912234428251326, -0.963971927277914, -0.993128599185095}
};


double a0, a1, a3, a5, a7, a9;
double exactInt = 0.8743332880;
void SetConstants(const double &eps);
double GetDistance(const Eigen::VectorXd &x);



// void PrintMat(std::vector< std::vector<double> >& M);
//void PrinVec(std::vector<double>& v);
void Cheb(const unsigned & m, Eigen::VectorXd &xg, Eigen::MatrixXd &C);
void GetParticle(const double & a, const double & b, const unsigned & n1, const unsigned& dim, Eigen::MatrixXd &x, Eigen::MatrixXd &xL);
void AssembleMatEigen(double& a, double& b, const unsigned& m, const unsigned& dim, const unsigned& np, Eigen::Tensor<double, 3, Eigen::RowMajor>  &PmX, Eigen::MatrixXd &Pg, Eigen::VectorXd &w, Eigen::MatrixXd &A, Eigen::VectorXd &F);
void SolWeightEigen(Eigen::MatrixXd &A, Eigen::VectorXd &F, Eigen::VectorXd &wp, Eigen::VectorXd &w_new);
void GetGaussInfo(const unsigned& m, const unsigned& dim, const std::vector<std::vector<double>> &GaussP, Eigen::MatrixXd &Pg, Eigen::VectorXd &w);
void GetChebXInfo(const unsigned& m, const unsigned& dim, const unsigned& np, Eigen::MatrixXd &xL, Eigen::Tensor<double, 3, Eigen::RowMajor>& PmX);
void Testing(double& a, double& b, const unsigned& m, const unsigned& dim, Eigen::MatrixXd &x, Eigen::VectorXd &w_new);

int main(int argc, char** args) {

//   std::vector<std::vector<std::vector<double>>> GR(13);
//   GR = {Gauss3, Gauss4, Gauss5, Gauss6, Gauss7, Gauss8, Gauss9, Gauss10, Gauss11, Gauss12, Gauss13, Gauss14, Gauss15};
//   Eigen::VectorXi N1(4);
//   N1 << 40, 80, 120, 150;
//
//   std::ofstream fout;
//   fout.open("PerformanceFullEigenLLT.txt", std::ofstream::app);
//   fout << "NG  " << "," << "NP     " << "," << "AssembleTime%  " << "," << "SolverTime%    " << "," << "ExecutionTime%" << std::endl;
//
//   for(unsigned i = 0; i <  N1.size(); i++) {
//     for(unsigned j = 0; j < GR.size(); j++) {

  clock_t begin1 = clock();
  unsigned dim = 2;
  double a = 0;
  double b = 1;
  unsigned n1 = 20;
  unsigned np = pow(n1, dim);
  unsigned m = 4;
  //unsigned ng = GR[j][0].size();


  Eigen::MatrixXd x;
  Eigen::MatrixXd xL;
  GetParticle(a, b, n1, dim, x, xL);


  Eigen::MatrixXd Pg;
  Eigen::VectorXd w;
  //GetGaussInfo(m, dim, GR[j], Pg, w);
  GetGaussInfo(m, dim, Gauss9, Pg, w);


  Eigen::Tensor<double, 3, Eigen::RowMajor> PmX;
  GetChebXInfo(m, dim, np, xL, PmX);



  Eigen::MatrixXd A;
  Eigen::VectorXd F;
  AssembleMatEigen(a, b, m, dim, np, PmX, Pg, w, A, F);


  Eigen::VectorXd wp(np);
  wp.fill(pow(b - a, dim) / np);

  double AssembleTime = (clock() - begin1) / ((float)CLOCKS_PER_SEC);


  clock_t begin2 = clock();

  Eigen::VectorXd w_new;
  SolWeightEigen(A, F, wp, w_new);

  double SolverTime = (clock() - begin2) / ((float)CLOCKS_PER_SEC);
  double TotalTime = AssembleTime + SolverTime;

  std::cout << "AssembleTime: " << (AssembleTime / TotalTime) * 100 << " SolverTime: " << (SolverTime / TotalTime) * 100 << "ExecutionTime: " << TotalTime << std::endl;
  //fout << ng << ",   " << np << ",   " << (AssembleTime / TotalTime) * 100  << ",        " << (SolverTime / TotalTime) * 100 << ",        " << TotalTime << std::endl;
  Testing(a, b, m, dim, x, w_new);


//   double integral = 0.;
//   double eps = 0.085;
//   SetConstants(eps);
// 
//   for(unsigned i = 0; i < w_new.size(); i++) {
//     double dg1 = GetDistance(x);
// 
//     double dg2 = dg1 * dg1;
//     double xi = 0;
//     if(dg1 < -eps){
//       xi = 0.;
//     }
//     else if(dg1 > eps) {
//       xi = 1.;
//     }
//     else {
//       xi = (a0 + dg1 * (a1 + dg2 * (a3 + dg2 * (a5 + dg2 * (a7 + dg2 * a9)))));
//     }
//     integral += xi * w_new(i);
//   }
//   
//   Eigen::VectorXd x1(2);
//   x1 << 1.,0.;
//   std::cout << GetDistance(x1) << std::endl;






//         }
//       }
//
//   fout.close();
  return 0;
}


void GetChebXInfo(const unsigned& m, const unsigned& dim, const unsigned& np, Eigen::MatrixXd &xL, Eigen::Tensor<double, 3, Eigen::RowMajor>& PmX) {

  PmX.resize(dim, m + 1, np);
  Eigen::MatrixXd Ptemp;
  Eigen::VectorXd xtemp;
  for(unsigned k = 0; k < dim; k++) {
    xtemp = xL.row(k);
    Cheb(m, xtemp, Ptemp);
    for(unsigned i = 0; i < m + 1; i++) {
      for(unsigned j = 0; j < np; j++) {
        PmX(k, i, j) = Ptemp(i, j);
      }
    }
  }
}

void GetGaussInfo(const unsigned& m, const unsigned& dim, const std::vector<std::vector<double>> &GaussP, Eigen::MatrixXd &Pg, Eigen::VectorXd &w) {
  Eigen::VectorXd xg = Eigen::VectorXd::Map(&GaussP[1][0], GaussP[1].size());
  Cheb(m, xg, Pg);
  Eigen::VectorXd wG = Eigen::VectorXd::Map(&GaussP[0][0], GaussP[0].size());
  if(dim == 1) {
    w = wG;
  }
  else {
    w = Eigen::kroneckerProduct(wG, wG).eval();
    if(dim == 3) {
      w = Eigen::kroneckerProduct(wG, w).eval();
    }
  }
}




void AssembleMatEigen(double& a, double& b, const unsigned& m, const unsigned& dim, const unsigned& np, Eigen::Tensor<double, 3, Eigen::RowMajor>  &PmX, Eigen::MatrixXd &Pg, Eigen::VectorXd &w, Eigen::MatrixXd &A, Eigen::VectorXd &F) {

  A.resize(pow(m + 1, dim), np);
  F.resize(pow(m + 1, dim));
  Eigen::VectorXi I(dim);
  Eigen::VectorXi N(dim);


  for(unsigned k = 0; k < dim ; k++) {
    N(k) = pow(m + 1, dim - k - 1);
  }

  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    I(0) = t / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = t % N(k - 1);
      I(k) = pk / N(k); // dimensional index over on the space of polynomaials
    }
    for(unsigned j = 0; j < np; j++) {
      double r = 1;
      for(unsigned k = 0; k < dim; k++) {
        r = r * PmX(k, I[k], j);
      }
      A(t, j) = r ;
    }

  }

  if(dim == 1) {
    F = (1 / pow(2, dim)) * pow(b - a, dim) * Pg * w;
  }
  else {
    Eigen::MatrixXd C ;
    C =  Eigen::kroneckerProduct(Pg, Pg);
    if(dim == 3) {
      C = Eigen::kroneckerProduct(Pg, C).eval();
    }
    F = (1 / pow(2, dim)) * pow(b - a, dim) * C * w;
  }
}


void SolWeightEigen(Eigen::MatrixXd &A, Eigen::VectorXd &F, Eigen::VectorXd &wp, Eigen::VectorXd &w_new) {

  Eigen::VectorXd y_temp = (A * A.transpose()).partialPivLu().solve(F - A * wp);
  Eigen::VectorXd w_new_temp = A.transpose() * y_temp;
  w_new = w_new_temp + wp;

}


void Testing(double &a, double &b, const unsigned &m, const unsigned &dim, Eigen::MatrixXd &x, Eigen::VectorXd &w_new) {

  double s = 0;
  for(unsigned i = 0; i < w_new.size(); i++) {
    double r = 1;
    for(unsigned k = 0; k < dim; k++) {
      r = r * x(k, i);
    }
    s = s + pow(r, m) * w_new(i);
  }
  double err = (s - pow((pow(b, m + 1) - pow(a, m + 1)) / (m + 1), dim)) / s;
  std::cout << "Error is: " << err << std::endl;

}


void Cheb(const unsigned &m, Eigen::VectorXd &xg, Eigen::MatrixXd &C) {

  C.resize(xg.size(), m + 1);
  for(unsigned i = 0; i < xg.size(); i++) {
    C(i, 0) = 1;
    C(i, 1) = xg(i);
    for(unsigned j = 2; j <= m; j++) {
      C(i, j) =  2 * xg(i) * C(i, j - 1) - C(i, j - 2);
    }
  }
  C.transposeInPlace();

}

void  GetParticle(const double &a, const double &b, const unsigned &n1, const unsigned &dim, Eigen::MatrixXd &x, Eigen::MatrixXd &xL) {
  double h = (b - a) / n1;
  x.resize(dim, pow(n1, dim));
  Eigen::VectorXi I(dim);
  Eigen::VectorXi N(dim);

  for(unsigned k = 0; k < dim ; k++) {
    N(k) = pow(n1, dim - k - 1);
  }

  for(unsigned p = 0; p < pow(n1, dim) ; p++) {
    I(0) = 1 + p / N(0);
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = p % N(k - 1);
      I(k) = 1 + pk / N(k);
    }
    //std::cout << I(0) << " " << I(1) << std::endl;

    for(unsigned k = 0; k < dim ; k++) {
      //std::srand(std::time(0));
      //double r = 2 * ((double) rand() / (RAND_MAX)) - 1;
      x(k, p) = a + h / 2 + (I(k) - 1) * h; // + 0.1 * r;

    }
  }


  xL.resize(dim, pow(n1, dim));
  Eigen::MatrixXd ID;
  ID.resize(dim, pow(n1, dim));
  ID.fill(1.);
  xL = (2. / (b - a)) * x - ((b + a) / (b - a)) * ID;


}




// void  GetParticleOnDisk(const double &a, const double &b, const unsigned &n1, const unsigned &dim, Eigen::MatrixXd &x, Eigen::MatrixXd &xL) {
//   double h = (b - a) / n1;
//   x.resize(dim, pow(n1, dim));
//   Eigen::VectorXi I(dim);
//   Eigen::VectorXi N(dim);
// 
//   for(unsigned k = 0; k < dim ; k++) {
//     N(k) = pow(n1, dim - k - 1);
//   }
// 
//   
//   double q=0;
//   
//   
//   
//   
//   for(unsigned p = 0; p < pow(n1, dim) ; p++) {
//     I(0) = 1 + p / N(0);
//     for(unsigned k = 1; k < dim ; k++) {
//       unsigned pk = p % N(k - 1);
//       I(k) = 1 + pk / N(k);
//     }
//     //std::cout << I(0) << " " << I(1) << std::endl;
// 
//     for(unsigned k = 0; k < dim ; k++) {
//       //std::srand(std::time(0));
//       //double r = 2 * ((double) rand() / (RAND_MAX)) - 1;
//       x(k, p) = a + h / 2 + (I(k) - 1) * h; // + 0.1 * r;
// 
//     }
//   }
// 
// 
//   xL.resize(dim, pow(n1, dim));
//   Eigen::MatrixXd ID;
//   ID.resize(dim, pow(n1, dim));
//   ID.fill(1.);
//   xL = (2. / (b - a)) * x - ((b + a) / (b - a)) * ID;
// 
// 
// }










void SetConstants(const double &eps) {
  a0 = 0.5; // 128./256.;
  a1 = pow(eps, -1.) * 1.23046875; // 315/256.;
  a3 = -pow(eps, -3.) * 1.640625; //420./256.;
  a5 = pow(eps, -5.) * 1.4765625; // 378./256.;
  a7 = -pow(eps, -7.) * 0.703125; // 180./256.;
  a9 = pow(eps, -9.) * 0.13671875; // 35./256.;
}

double GetDistance(const Eigen::VectorXd &x) {
  return (-x[0] + x[1] + 0.5) /sqrt(2.);

}







// void PrintMat(std::vector< std::vector<double> >& M) {
//
//   for(unsigned i = 0; i < M.size(); i++) {
//     for(unsigned j = 0; j < M[i].size(); j++) {
//       std::cout << M[i][j] << " ";
//     }
//     std::cout << std::endl;
//   }
//   std::cout << "\n" << std::endl;
// }
//
//
// void PrinVec(std::vector<double>& v) {
//   for(unsigned i = 0; i < v.size(); i++) {
//
//     std::cout << v[i] << " ";
//   }
//   std::cout << "\n" << std::endl;
// }
























