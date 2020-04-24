#include <iostream>
#include <vector>
#include <math.h>
#include <ctime>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include <fstream>


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

const std::vector<std::vector < double > > Gauss8 = {{0.362683783378362, 0.313706645877887, 0.222381034453375, 0.101228536290376, 0.362683783378362, 0.313706645877887, 0.222381034453375, 0.101228536290376}, {-0.183434642495650, -0.525532409916329, -0.796666477413627, -0.960289856497536, 0.183434642495650, 0.525532409916329, 0.796666477413627, 0.960289856497536}
};

const std::vector<std::vector < double > > Gauss9 = {{0.330239355001260, 0.312347077040003, 0.260610696402935, 0.180648160694857, 0.081274388361574, 0.312347077040003, 0.260610696402935, 0.180648160694857, 0.081274388361574}, {0.000000000000000, -0.324253423403809, -0.613371432700590, -0.836031107326636, -0.968160239507626, 0.324253423403809, 0.613371432700590, 0.836031107326636, 0.968160239507626}
};

const std::vector<std::vector < double > > Gauss10 = {{0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688, 0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688}, {-0.148874338981631, -0.433395394129247, -0.679409568299024, -0.865063366688985, -0.973906528517172, 0.148874338981631, 0.433395394129247, 0.679409568299024, 0.865063366688985, 0.973906528517172}
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

const std::vector<std::vector < double > > Gauss18 = {{0.169142382963144, 0.164276483745833, 0.154684675126265, 0.140642914670651, 0.122555206711479, 0.100942044106287, 0.076425730254889, 0.049714548894970, 0.021616013526483, 0.169142382963144, 0.164276483745833, 0.154684675126265, 0.140642914670651, 0.122555206711479, 0.100942044106287, 0.076425730254889, 0.049714548894970, 0.021616013526483}, {-0.084775013041735, -0.251886225691505, -0.411751161462843, -0.559770831073948, -0.691687043060353, -0.803704958972523, -0.892602466497556, -0.955823949571398, -0.991565168420931, 0.084775013041735, 0.251886225691505, 0.411751161462843, 0.559770831073948, 0.691687043060353, 0.803704958972523, 0.892602466497556, 0.955823949571398, 0.991565168420931}
};

const std::vector<std::vector < double > > Gauss19 = {{0.161054449848784, 0.158968843393954, 0.152766042065860, 0.142606702173607, 0.128753962539336, 0.111566645547334, 0.091490021622450, 0.069044542737641, 0.044814226765700, 0.019461788229727, 0.158968843393954, 0.152766042065860, 0.142606702173607, 0.128753962539336, 0.111566645547334, 0.091490021622450, 0.069044542737641, 0.044814226765700, 0.019461788229727}, {0.000000000000000, -0.160358645640225, -0.316564099963630, -0.464570741375961, -0.600545304661681, -0.720966177335229, -0.822714656537143, -0.903155903614818, -0.960208152134830, -0.992406843843584, 0.160358645640225, 0.316564099963630, 0.464570741375961, 0.600545304661681, 0.720966177335229, 0.822714656537143, 0.903155903614818, 0.960208152134830, 0.992406843843584}
};

const std::vector<std::vector < double > > Gauss20 = {{0.152753387130726, 0.149172986472604, 0.142096109318382, 0.131688638449177, 0.118194531961518, 0.101930119817240, 0.083276741576705, 0.062672048334109, 0.040601429800387, 0.017614007139152, 0.152753387130726, 0.149172986472604, 0.142096109318382, 0.131688638449177, 0.118194531961518, 0.101930119817240, 0.083276741576705, 0.062672048334109, 0.040601429800387, 0.017614007139152}, {0.076526521133497, 0.227785851141645, 0.373706088715420, 0.510867001950827, 0.636053680726515, 0.746331906460151, 0.839116971822219, 0.912234428251326, 0.963971927277914, 0.993128599185095, -0.076526521133497, -0.227785851141645, -0.373706088715420, -0.510867001950827, -0.636053680726515, -0.746331906460151, -0.839116971822219, -0.912234428251326, -0.963971927277914, -0.993128599185095}
};



void PrintMat(std::vector< std::vector<double> >& M);
void PrinVec(std::vector<double>& v);
void  GetParticle(const double & a, const double & b, const unsigned & n1, const unsigned& dim, std::vector < std::vector < double > >& x, std::vector < std::vector < double > >& xL);
void Cheb(const unsigned & m, const std::vector < double > &xg, std::vector< std::vector < double > > & C);
void GetGaussPointsWeights(const unsigned dim, const std::vector<std::vector<double>> GaussP, std::vector<std::vector<double>>& xg, std::vector<double>& w);
void Kron(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B, std::vector<std::vector<double>>& C);

// void AssembleMat(double& a, double& b, std::vector<std::vector< std::vector < double > >>&  PmX, std::vector<std::vector<double>>& Pg, std::vector<double>& w, const unsigned& dim, const unsigned& np, const unsigned& m, std::vector<std::vector<double>>& A, std::vector<double>& F);
// void SolWeight(std::vector<std::vector<double>>& A1, std::vector<double>& F1, std::vector<double>& wp, std::vector<double>& w_new);

void AssembleMatEigen(double& a, double& b, std::vector<std::vector< std::vector < double > >>&  PmX, std::vector<std::vector<double>>& Pg, std::vector<double>& w, const unsigned& dim, const unsigned& np, const unsigned& m, Eigen::MatrixXd &A, Eigen::VectorXd &F);
void SolWeightEigen(Eigen::MatrixXd &A, Eigen::VectorXd &F, std::vector<double>& wp, std::vector<double>& w_new);


int main(int argc, char** args) {


  std::vector<std::vector<std::vector<double>>> GR(15);
  GR = {Gauss3, Gauss4, Gauss5, Gauss6, Gauss7, Gauss8, Gauss9, Gauss10, Gauss11, Gauss12, Gauss13, Gauss14, Gauss15};

  std::vector<double> N1 = {10, 20, 40, 80};

  std::ofstream fout;
  fout.open("Performance.txt", std::ofstream::app);
  fout << "NG  " << "," << "NP     " << "," << "AssembleTime%  " << "," << "SolverTime%    " << "," << "ExecutionTime%" << std::endl;


  for(unsigned i = 0; i < N1.size(); i++) {
    for(unsigned j = 0; j < GR.size(); j++) {
      clock_t begin1 = clock();
      unsigned dim = 3;
      double a = 2;
      double b = 5;
      unsigned n1 = N1[i];
      unsigned np = pow(n1, dim);
      unsigned ng = GR[j][0].size();
      unsigned m = 4; // Note that m>=2*ng-1;
      
      std::vector<double> wp;
      wp.assign(np, pow(b - a, dim) / np);
      std::vector < std::vector < double > > x;
      std::vector<std::vector<double>> xL;
      GetParticle(a, b, n1, dim, x, xL);
      std::vector < std::vector < double > > xg;
      std::vector<double> w;
      GetGaussPointsWeights(dim, GR[j], xg, w);
      std::vector<std::vector<double>> Pg;
      Cheb(m, GR[j][1], Pg);
      std::vector<std::vector< std::vector < double > >>  PmX(dim);
      for(unsigned k = 0; k < dim; k++) {
        Cheb(m, xL[k], PmX[k]);
      }


      std::vector<std::vector<double>> A1;
      //std::vector<double> F1;
      //AssembleMat(a, b, PmX, Pg, w, dim, np, m, A1, F1);
      Eigen::MatrixXd A;
      Eigen::VectorXd F;
      AssembleMatEigen(a, b, PmX, Pg, w, dim, np, m, A, F);

      clock_t stop1 = clock();
      double AssembleTime = (stop1 - begin1) / ((float)CLOCKS_PER_SEC);

      clock_t begin2 = clock();
      std::vector<double> w_new;
      //SolWeight(A1, F1, wp, w_new);
      SolWeightEigen(A, F, wp, w_new);

      clock_t stop2 = clock();
      double SolverTime = (stop2 - begin2) / ((float)CLOCKS_PER_SEC);

      double TotalTime = AssembleTime + SolverTime;

      std::cout << "AssembleTime: " << (AssembleTime / TotalTime) * 100 << " SolverTime: " << (SolverTime / TotalTime) * 100 << std::endl;


      double s = 0;
      for(unsigned i = 0; i < w_new.size(); i++) {
        double r = 1;
        for(unsigned k = 0; k < dim; k++) {
          r = r * x[k][i];
        }
        s = s + pow(r, m) * w_new[i];
      }

      double err = (s - pow((pow(b, m + 1) - pow(a, m + 1)) / (m + 1), dim)) / s;
      std::cout << "Eror is: " << err << std::endl;

      fout << ng << ",   " << np << ",   " << (AssembleTime / TotalTime) * 100  << ",        " << (SolverTime / TotalTime) * 100 << ",        " << TotalTime << std::endl;

    }
  }
  fout.close();

  return 0;
}



// void SolWeight(std::vector<std::vector<double>>& A1, std::vector<double>& F1, std::vector<double>& wp, std::vector<double>& w_new) {
// 
//   Eigen::VectorXd F = Eigen::VectorXd::Map(&F1[0], F1.size());
//   Eigen::VectorXd wP = Eigen::VectorXd::Map(&wp[0], wp.size());
// 
//   Eigen::MatrixXd A(A1.size(), A1[0].size());
//   for(int i = 0; i < A1.size(); ++i) {
//     A.row(i) = Eigen::VectorXd::Map(&A1[i][0], A1[0].size());
//   }
// 
//   Eigen::VectorXd y_temp = (A * A.transpose()).fullPivLu().solve(F - A * wP);
//   Eigen::VectorXd w_new_temp = A.transpose() * y_temp + wP;
// 
//   w_new.resize(w_new_temp.size());
//   Eigen::VectorXd::Map(&w_new[0], w_new_temp.size()) = w_new_temp;
// 
// }


// void AssembleMat(double& a, double& b, std::vector<std::vector< std::vector < double > >>&  PmX, std::vector<std::vector<double>>& Pg, std::vector<double>& w, const unsigned& dim, const unsigned& np, const unsigned& m, std::vector<std::vector<double>>& A, std::vector<double>& F) {
// 
//   A.resize(pow(m + 1, dim));
//   F.resize(pow(m + 1, dim));
//   std::vector<unsigned> I(dim);
//   std::vector<unsigned> J(dim);
//   std::vector<unsigned> N(dim);
// 
// 
//   for(unsigned k = 0; k < dim ; k++) {
//     N[k] = pow(m + 1, dim - k - 1);
//   }
// 
//   unsigned ng = Pg[0].size();
//   for(unsigned k = 0; k < dim ; k++) {
//     J[k] = pow(ng + 1, dim - k - 1);
//   }
// 
//   for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
//     A[t].resize(np);
//     I[0] = t / N[0];
//     for(unsigned k = 1; k < dim ; k++) {
//       unsigned pk = t % N[k - 1];
//       I[k] = pk / N[k]; // dimensional index over on the space of polynomaials
//     }
//     for(unsigned j = 0; j < np; j++) {
//       double r = 1;
//       for(unsigned k = 0; k < dim; k++) {
//         r = r * PmX[k][I[k]][j];
//       }
//       A[t][j] = r ;
//     }
// 
//   }
// 
//   Eigen::VectorXd wE = Eigen::VectorXd::Map(&w[0], w.size());
//   Eigen::VectorXd F1(F.size());
//   Eigen::MatrixXd C ;
//   Eigen::MatrixXd PgE(Pg.size(), Pg[0].size());
// 
//   
//   for(int i = 0; i < Pg.size(); ++i) {
//     PgE.row(i) = Eigen::VectorXd::Map(&Pg[i][0], Pg[0].size());
//   }
// 
// 
//   if(dim == 1) {
//     F1 = (1 / pow(2, dim)) * pow(b - a, dim) * PgE * wE;
//   }
//   else {
//     C =  Eigen::kroneckerProduct(PgE, PgE);
//     if(dim == 3) {
//       C = Eigen::kroneckerProduct(PgE, C).eval();
//     }
//     F1 = (1 / pow(2, dim)) * pow(b - a, dim) * C * wE;
//   }
// 
//   Eigen::VectorXd::Map(&F[0], F1.size()) = F1;
// 
// }





void AssembleMatEigen(double& a, double& b, std::vector<std::vector< std::vector < double > >>&  PmX, std::vector<std::vector<double>>& Pg, std::vector<double>& w, const unsigned& dim, const unsigned& np, const unsigned& m, Eigen::MatrixXd &A, Eigen::VectorXd &F) {


  std::vector<std::vector<double>> A1;
  A1.resize(pow(m + 1, dim));
  F.resize(pow(m + 1, dim));
  std::vector<unsigned> I(dim);
  //std::vector<unsigned> J(dim);
  std::vector<unsigned> N(dim);


  for(unsigned k = 0; k < dim ; k++) {
    N[k] = pow(m + 1, dim - k - 1);
  }

 // unsigned ng = Pg[0].size();
//   for(unsigned k = 0; k < dim ; k++) {
//     J[k] = pow(ng + 1, dim - k - 1);
//   }

  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    A1[t].resize(np);
    I[0] = t / N[0];
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = t % N[k - 1];
      I[k] = pk / N[k]; // dimensional index over on the space of polynomaials
    }
    for(unsigned j = 0; j < np; j++) {
      double r = 1;
      for(unsigned k = 0; k < dim; k++) {
        r = r * PmX[k][I[k]][j];
      }
      A1[t][j] = r ;
    }

  }


  A.resize(A1.size(), A1[0].size());
  for(int i = 0; i < A1.size(); ++i) {
    A.row(i) = Eigen::VectorXd::Map(&A1[i][0], A1[0].size());
  }



  Eigen::VectorXd wE = Eigen::VectorXd::Map(&w[0], w.size());
  Eigen::MatrixXd C ;

  Eigen::MatrixXd PgE(Pg.size(), Pg[0].size());
  for(int i = 0; i < Pg.size(); ++i) {
    PgE.row(i) = Eigen::VectorXd::Map(&Pg[i][0], Pg[0].size());
  }


  if(dim == 1) {
    F = (1 / pow(2, dim)) * pow(b - a, dim) * PgE * wE;
  }
  else {
    C =  Eigen::kroneckerProduct(PgE, PgE);
    if(dim == 3) {
      C = Eigen::kroneckerProduct(PgE, C).eval();
    }
    F = (1 / pow(2, dim)) * pow(b - a, dim) * C * wE;
  }
}




void SolWeightEigen(Eigen::MatrixXd &A, Eigen::VectorXd &F, std::vector<double>& wp, std::vector<double>& w_new) {

  Eigen::VectorXd wP = Eigen::VectorXd::Map(&wp[0], wp.size());

  Eigen::VectorXd y_temp = (A * A.transpose()).fullPivLu().solve(F - A * wP);
  Eigen::VectorXd w_new_temp = A.transpose() * y_temp + wP;

  w_new.resize(w_new_temp.size());
  Eigen::VectorXd::Map(&w_new[0], w_new_temp.size()) = w_new_temp;

}




void Cheb(const unsigned & m, const std::vector < double > &xg, std::vector< std::vector < double > > & C) {
  // C1[i][j] = P_{j}(x_{i})
  std::vector< std::vector < double > >  C1;
  C1.resize(xg.size());
  for(unsigned i = 0; i < xg.size(); i++) {
    C1[i].resize(m + 1);
    C1[i][0] = 1;
    C1[i][1] = xg[i];
    for(unsigned j = 2; j <= m; j++) {
      C1[i][j] =  2 * xg[i] * C1[i][j - 1] - C1[i][j - 2];
    }
  }

  C.resize(C1[0].size());
  for(unsigned i = 0; i < m + 1; i++) {
    C[i].resize(C1.size());
    for(unsigned j = 0; j < C1.size(); j++) {
      C[i][j] = C1[j][i];
    }
  }
}

void  GetParticle(const double & a, const double & b, const unsigned & n1, const unsigned& dim, std::vector < std::vector < double > >& x, std::vector < std::vector < double > >& xL) {
  unsigned np = pow(n1, dim);
  double h = (b - a) / n1;
  x.resize(dim);
  std::vector < unsigned > I(dim);
  std::vector < unsigned > N(dim);

  for(unsigned k = 0; k < dim ; k++) {
    N[k] = pow(n1, dim - k - 1);
    x[k].resize(np);
  }

  for(unsigned p = 0; p < np ; p++) {
    I[0] = 1 + p / N[0];
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = p % N[k - 1];
      I[k] = 1 + pk / N[k];
    }
    //std::cout << I[0] << " " << I[1] << std::endl;

    for(unsigned k = 0; k < dim ; k++) {
      //std::srand(std::time(0));
      //double r = 2 * ((double) rand() / (RAND_MAX)) - 1;
      x[k][p] = a + h / 2 + (I[k] - 1) * h;// + 0.1 * r;

    }
  }

  xL.resize(dim);
  for(unsigned i = 0; i < dim; i++) {
    xL[i].resize(np);
    for(unsigned j = 0; j < np; j++) {
      xL[i][j] += (2 * x[i][j] - (b + a)) / (b - a);
    }
  }


}


void GetGaussPointsWeights(const unsigned dim, const std::vector<std::vector<double>> GaussP, std::vector<std::vector<double>>& xg, std::vector<double>& w) {

  unsigned ng = GaussP[0].size();
  xg.resize(dim);
  w.resize(pow(ng, dim));
  std::vector < unsigned > I(dim);
  std::vector < unsigned > N(dim);

  for(unsigned k = 0; k < dim ; k++) {
    xg[k].resize(pow(ng, dim));
    N[k] = pow(ng, dim - k - 1);
  }

  for(unsigned p = 0; p < pow(ng, dim) ; p++) {
    I[0] = p / N[0];
    for(unsigned k = 1; k < dim ; k++) {
      unsigned pk = p % N[k - 1];
      I[k] = pk / N[k];
    }
    //std::cout<< I[0] << " " << I[1] << " " <<I[2] << std::endl;
    for(unsigned i = 0; i < dim; i++) {
      xg[i][p] = GaussP[1][I[i]];
    }
    double r = 1;
    for(unsigned k = 0; k < dim; k++) {
      r = r * GaussP[0][I[k]];
      w[p] = r;
    }

  }
  //PrinVec(w);
  //PrintMat(xg);

}




// void Kron(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B, std::vector<std::vector<double>>& C) {
//   A:mxn, B:pxq, C--> pmxqn
//
//   unsigned rowA = A.size();
//   unsigned colA = A[0].size();
//   unsigned rowB = B.size();
//   unsigned colB = B[0].size();
//   C.resize(rowB * rowA);
//   for(unsigned k = 0; k < rowB * rowA; k++) {
//     C[k].resize(colA * colB);
//   }
//
//   for(unsigned i = 0; i < rowA; i++) {
//     for(unsigned j = 0; j < colA; j++) {
//       unsigned startRow = i * rowB;
//       unsigned startCol = j * colB;
//       for(unsigned k = 0; k < rowB; k++) {
//         for(unsigned l = 0; l < colB; l++) {
//           C[startRow + k][startCol + l] = A[i][j] * B[k][l];
//         }
//       }
//     }
//   }
//
// }



void PrintMat(std::vector< std::vector<double> >& M) {

  for(unsigned i = 0; i < M.size(); i++) {
    for(unsigned j = 0; j < M[i].size(); j++) {
      std::cout << M[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "\n" << std::endl;
}


void PrinVec(std::vector<double>& v) {
  for(unsigned i = 0; i < v.size(); i++) {

    std::cout << v[i] << " ";
  }
  std::cout << "\n" << std::endl;
}






// void Cheb1D(const unsigned & m, const std::vector < double > &xg, std::vector< std::vector < double > > & C1) {
//   // C1[i][j] = P_{j}(x_{i})
//   C1.resize(xg.size());
//   for(unsigned i = 0; i < xg.size(); i++) {
//     C1[i].resize(m + 1);
//     C1[i][0] = 1;
//     C1[i][1] = xg[i];
//     for(unsigned j = 2; j <= m; j++) {
//       C1[i][j] =  2 * xg[i] * C1[i][j - 1] - C1[i][j - 2];
//     }
//   }
//
// }




























