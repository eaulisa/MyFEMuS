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














void PrintMat(std::vector< std::vector<double> >& M);
void PrinVec(std::vector<double>& v);
void  GetParticle(const double & a, const double & b, const unsigned & n1, const unsigned& dim, std::vector < std::vector < double > >& x, std::vector < std::vector < double > >& xL);
void Cheb(const unsigned & m, const std::vector < double > &xg, std::vector< std::vector < double > > & C);
void GetGaussPointsWeights(const unsigned dim, const std::vector<std::vector<double>> GaussP, std::vector<std::vector<double>>& xg, std::vector<double>& w);
void Kron(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B, std::vector<std::vector<double>>& C);
void AssembleMat(double& a, double& b, const unsigned& ng, std::vector<std::vector< std::vector < double > >>&  PmX, std::vector<std::vector< std::vector < double > >>&  PmG, std::vector<double>& w, const unsigned& dim, const unsigned& np, const unsigned& m, std::vector<std::vector<double>>& A, std::vector<double>& F);
void SolWeight(std::vector<std::vector<double>>& A1, std::vector<double>& F1, std::vector<double>& wp, std::vector<double>& w_new);

int main(int argc, char** args) {


  std::vector<std::vector<std::vector<double>>> GR(7);
  GR = { Gauss1, Gauss2, Gauss3, Gauss4, Gauss5,Gauss6,Gauss7};

  std::vector<double> N1 = {40, 80 ,120, 150};
  

  std::ofstream fout;
  fout.open ("AssVsSolver.txt", std::ofstream::app);
  fout << "NG"<<","<<"NP"<<","<<"AS" <<","<< "Sol" <<std::endl;

  for(unsigned i = 0; i < N1.size(); i++) {
    for(unsigned j = 0; j < GR.size(); j++) {
      clock_t begin1 = clock();
      unsigned dim = 3;
      double a = 2;
      double b = 5;
      unsigned n1 = N1[i];
      unsigned np = pow(n1, dim);
      unsigned ng = GR[j][0].size();
      unsigned m = 2 * ng - 1;
      std::vector<double> wp;
      wp.assign(np, pow(b - a, dim) / np);
      std::vector < std::vector < double > > x;
      std::vector<std::vector<double>> xL;
      GetParticle(a, b, n1, dim, x, xL);
      std::vector < std::vector < double > > xg;
      std::vector<double> w;
      GetGaussPointsWeights(dim, GR[j], xg, w);
      std::vector<std::vector< std::vector < double > >>  PmG(dim);
      std::vector<std::vector< std::vector < double > >>  PmX(dim);
      for(unsigned k = 0; k < dim; k++) {
        Cheb(m, xg[k], PmG[k]);
        Cheb(m, xL[k], PmX[k]);
      }


      std::vector<std::vector<double>> A1;
      std::vector<double> F1;
      AssembleMat(a, b, ng, PmX, PmG, w, dim, np, m, A1, F1);
      clock_t stop1 = clock();
      double AssembleTime = (stop1-begin1) /  ((float)CLOCKS_PER_SEC);
          
      clock_t begin2 = clock();
      std::vector<double> w_new;
      SolWeight(A1, F1, wp, w_new);
      clock_t stop2 = clock();
      double SolverTime = (stop2-begin2) /  ((float)CLOCKS_PER_SEC);
      
      double TotalTime = AssembleTime + SolverTime;
      
      std::cout << "AssembleTime: " << (AssembleTime / TotalTime)*100 << " SolverTime: " << (SolverTime / TotalTime)*100 << std::endl;
      

      double s = 0;
      for(unsigned i = 0; i < w_new.size(); i++) {
        double r = 1;
        for(unsigned k = 0; k < dim; k++) {
          r = r * x[k][i];
        }
        s = s + pow(r, m) * w_new[i];
      }

      double err = (s - pow((pow(b, m + 1) - pow(a, m + 1)) / (m + 1), dim)) / s;

      std::cout << "Error is: " << err << std::endl;

      fout << ng <<","<< np <<","<< (AssembleTime / TotalTime)*100  << "," << (SolverTime / TotalTime)*100 <<std::endl;

    }
  }
  fout.close();
  return 0;
}



void SolWeight(std::vector<std::vector<double>>& A1, std::vector<double>& F1, std::vector<double>& wp, std::vector<double>& w_new) {

  Eigen::VectorXd F = Eigen::VectorXd::Map(&F1[0], F1.size());
  Eigen::VectorXd wP = Eigen::VectorXd::Map(&wp[0], wp.size());

  Eigen::MatrixXd A(A1.size(), A1[0].size());
  for(int i = 0; i < A1.size(); ++i) {
    A.row(i) = Eigen::VectorXd::Map(&A1[i][0], A1[0].size());
  }

  Eigen::VectorXd y_temp = (A * A.transpose()).fullPivLu().solve(F - A * wP);
  Eigen::VectorXd w_new_temp = A.transpose() * y_temp + wP;

  w_new.resize(w_new_temp.size());
  Eigen::VectorXd::Map(&w_new[0], w_new_temp.size()) = w_new_temp;

}


void AssembleMat(double& a, double& b, const unsigned& ng, std::vector<std::vector< std::vector < double > >>&  PmX, std::vector<std::vector< std::vector < double > >>&  PmG, std::vector<double>& w, const unsigned& dim, const unsigned& np, const unsigned& m, std::vector<std::vector<double>>& A, std::vector<double>& F) {

  A.resize(pow(m + 1, dim));
  F.resize(pow(m + 1, dim));
  std::vector<unsigned> I(dim);
  std::vector<unsigned> J(dim);
  std::vector<unsigned> N(dim);


  for(unsigned k = 0; k < dim ; k++) {
    N[k] = pow(m + 1, dim - k - 1);
  }




  for(unsigned t = 0; t < pow(m + 1, dim) ; t++) { // multidimensional index on the space of polynomaials
    A[t].resize(np);
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
      A[t][j] = r ;
    }


    F[t] = 0.;
    for(unsigned j = 0; j < pow(ng, dim); j++) { // multidimensional index on the gauss points
      double r = 1;
      for(unsigned k = 0; k < dim; k++) {
        r = r * PmG[k][I[k]][j] ;
      }
      F[t] += (1 / pow(2, dim)) * pow(b - a, dim) * r * w[j] ;
    }
  }


  // unsigned ng = Pg[0].size();



//   Eigen::VectorXd wE = Eigen::VectorXd::Map(&w[0], w.size());
//
//   Eigen::VectorXd F1(F.size());
//   Eigen::MatrixXd C ;
//
//   Eigen::MatrixXd PgE(Pg.size(), Pg[0].size());
//   for(int i = 0; i < Pg.size(); ++i) {
//     PgE.row(i) = Eigen::VectorXd::Map(&Pg[i][0], Pg[0].size());
//   }
//
//
//   if(dim == 1) {
//
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



//     if(dim == 1) {
//       for(unsigned t = 0; t < pow(m + 1, dim) ; t++) {
//         for(unsigned j = 0; j < ng; j++) {
//           F[t] += (1 / pow(2, dim)) * pow(b - a, dim) * Pg[t][j] * w[j];
//         }
//       }
//     }
//
//     else {
//       std::vector<std::vector<double>> C2;
//
//       Kron(Pg, Pg, C2);
//       if(dim == 3) {
//         Kron(Pg, C2, C2);
//       }
//
//       for(unsigned t = 0; t < pow(m + 1, dim) ; t++) {
//         for(unsigned j = 0; j < pow(ng, dim); j++) {
//           F[t] += (1 / pow(2, dim)) * pow(b - a, dim) * C2[t][j] * w[j];
//         }
//       }
//
//     }











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




void Kron(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B, std::vector<std::vector<double>>& C) {
  //A:mxn, B:pxq, C--> pmxqn

  unsigned rowA = A.size();
  unsigned colA = A[0].size();
  unsigned rowB = B.size();
  unsigned colB = B[0].size();
  C.resize(rowB * rowA);
  for(unsigned k = 0; k < rowB * rowA; k++) {
    C[k].resize(colA * colB);
  }

  for(unsigned i = 0; i < rowA; i++) {
    for(unsigned j = 0; j < colA; j++) {
      unsigned startRow = i * rowB;
      unsigned startCol = j * colB;
      for(unsigned k = 0; k < rowB; k++) {
        for(unsigned l = 0; l < colB; l++) {
          C[startRow + k][startCol + l] = A[i][j] * B[k][l];
        }
      }
    }
  }

}



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



























