#include "FemusInit.hpp"


#include "ConicAdaptiveRefinement.hpp"

using namespace femus;


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  //example inputs
  // std::vector<std::vector<double>> y = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  // std::vector<std::vector<double>> y = {{1, -1}, {1, 1}, {-1, 1},{-1, -1}};
  // std::vector<std::vector<double>> yi = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};

  unsigned VType = 2, PType = 3;
  unsigned dofsV = 7;
  unsigned dofsP = 1;
  unsigned dofsAll = 2 * dofsV + 2 * dofsP;


  std::vector<std::vector<double>> V(2,std::vector<double>(dofsV, 2.));
  std::vector<double> P1(dofsP, 1.);
  std::vector<double> P2(dofsP, 2.);

  std::vector<double> res(dofsAll, 0.);
  std::vector<double> jac(dofsAll * dofsAll, 0.);

  // unsigned elType = 3;
  //std::vector<std::vector<double>> xv = {{-1., 1., 1., -1., 0., 1., 0., -1., 0.}, {-1., -1., 1., 1., -1., 0., 1., 0., 0.}};
  //std::vector<std::vector<double>> xv = {{1., 1., -1., -1.,  1., 0., -1., 0., 0.}, { -1., 1., 1.,-1., 0., 1., 0.,-1., 0.}};
  unsigned elType = 4;
  std::vector<std::vector<double>> xv = {{-1., 3., -1., 1., 1., -1., 1. / 3.}, {-1., -1., 3., -1., 1., 1., 1. / 3.}};
  //TODO

  double rho1 = 1., rho2 = 2., mu1 = .2, mu2 = 0.4, sigma = 1., dt = 0.01;

  //Matrix to store calculated coefficients
  //Coeficients of conics in physical system
  std::vector <double> A = {1, 0, 1., 0, 0, -0.5};
  std::vector <double> Ap;

  ConicAdaptiveRefinement cad;

  Data *data = new Data(VType, PType, V, P1, P2, res, jac, xv, elType, rho1, rho2, mu1, mu2, sigma, dt, A);

  cad.SetDataPointer(data);

  std::vector<std::vector<double>> y;
  cad.GetXInParentElement(xv, y);
  std::vector<std::vector<double>> yi = cad.GetXiInParentElement();

  cad.GetConicsInTargetElement(y, A, Ap);
  tuple <double, double, double> a = cad.AdaptiveRefinement(8, y, yi, Ap);

  delete data;

  double exact;



  exact = 0.39269908169872414;

  //exact = M_PI * 0.5;
  std::cout.precision(14);
  std::cout << " Analytic Area1 = " << exact << std::endl;
  std::cout << " Computed Area1 = " << std::get<0>(a) << std::endl;
  std::cout << " Relative Error = " << fabs((exact - std::get<0>(a)) / exact) << std::endl;

  exact = 8 - M_PI * 0.5;
  std::cout << " Analytic Area2 = " << exact << std::endl;
  std::cout << " Computed Area2 = " << std::get<1>(a) << std::endl;
  std::cout << " Relative Error = " << fabs((exact - std::get<1>(a)) / exact) << std::endl;

  exact = 2 * M_PI * sqrt(2.) / 2.;
  std::cout << " Analytic perimeter = " << exact << std::endl;
  std::cout << " Computed perimeter = " << std::get<2>(a) << std::endl;
  std::cout << " Relative Error = " << fabs((exact - std::get<2>(a)) / exact) << std::endl;
  //

  return 0;
}


void AssembleNavierStokes(Data *data, const std::vector <double> &phiV, const std::vector <double> &phiV_x, const std::vector <double> &phiP,
                          const double &C, const double &weight, const double &weight1, const double &weight2, const double &weightI,
                          const std::vector <double> &N, const double &kappa, const double &dsN, const double &eps) {


  const unsigned &dim = data->_V.size();
  const unsigned &nDofsV = phiV.size();
  const unsigned &nDofsP = phiP.size();
  unsigned nDofsVP = dim * nDofsV + 2 * nDofsP;

  std::vector < double > solVg(dim, 0);
  std::vector < std::vector < double > > SolVg_x(dim, std::vector<double> (dim, 0));
  for(unsigned i = 0; i < nDofsV; i++) {
    for(unsigned  K = 0; K < dim; K++) {
      solVg[K] += data->_V[K][i] * phiV[i];
      for(unsigned J = 0; J < dim; J++) {
        SolVg_x[K][J] += data->_V[K][i] * phiV_x[i * dim + J];
      }
    }
  }

  double solP1g = 0;
  double solP2g = 0;

  for(unsigned i = 0; i < nDofsP; i++) {
    solP1g += phiP[i] * data->_P1[i];
    solP2g += phiP[i] * data->_P2[i];
  }

  double rho = data->_rho1 * weight1 + data->_rho2 * weight2;
  double mu = data->_mu1 * weight1 + data->_mu2 * weight2;
  double rhoC = rho;

  //double rho = rho1 * C + rho2 * (1. - C);
  //double mu = mu1 * C + mu2 * (1. - C);
  //double rhoC = rho1 * C + rho2 * (1. - C);

  // *** phiV_i loop ***
  for(unsigned i = 0; i < nDofsV; i++) {
    for(unsigned  I = 0; I < dim; I++) {  //momentum equation in I
      double NSV = 0.;
      for(unsigned J = 0; J < dim; J++) {  // second index J in each equation
        NSV   +=  mu * phiV_x[i * dim + J] * (SolVg_x[I][J] + SolVg_x[J][I]); // diffusion
        //NSV   +=  rho * phiV[i] * solVg[J] * SolVg_x[I][J]; // nonlinear term
      }
      NSV += - phiV_x[i * dim + I] * (solP1g * weight1 + solP2g * weight2);  // pressure gradient
      NSV += rho * phiV[i] * solVg[I] / data->_dt ;
      //NSV += - rhoC * phiV[i] * g[I]; // gravity term
      data->_res[I * nDofsV + i] -=  NSV * weight;
      if(weightI != 0.) {
        data->_res[I * nDofsV + i] += -data->_sigma * phiV[i] * N[I] * weight * weightI * kappa * dsN;
      }
    }
  } // end phiV_i loop

  // *** phiP_i loop ***
  for(unsigned i = 0; i < nDofsP; i++) {
    for(int I = 0; I < dim; I++) {
      data->_res[dim * nDofsV + i] += - SolVg_x[I][I] * phiP[i]  * weight * weight1; //continuity
      data->_res[dim * nDofsV + nDofsP + i] += - SolVg_x[I][I] * phiP[i]  * weight * weight2; //continuity
    }
    if(C == 0.)
      data->_res[dim * nDofsV + i] += - solP1g * phiP[i]  * weight * eps; //penalty
    if(C == 1.)
      data->_res[dim * nDofsV + nDofsP + i] += - solP2g * phiP[i]  * weight * eps; //penalty
  } // end phiP_i loop


  for(unsigned i = 0; i < nDofsV; i++) {
    for(unsigned I = 0; I < dim; I++) { //row velocity blocks or dimension
      unsigned VIrow = I * nDofsV + i;
      for(unsigned j = 0; j < nDofsV; j++) {
        unsigned VIcolumn = I * nDofsV + j;
        data->_jac[ VIrow * nDofsVP + VIcolumn] += rho * phiV[i] * phiV[j] * weight / data->_dt; // inertia

        for(unsigned J = 0; J < dim ; J++) { //column velocity blocks or dimension
          unsigned VJcolumn = J * nDofsV + j;
          data->_jac[ VIrow * nDofsVP + VIcolumn ] += mu * phiV_x[i * dim + J] * phiV_x[j * dim + J] * weight; //diagonal diffusion
          data->_jac[ VIrow * nDofsVP + VJcolumn ] += mu * phiV_x[i * dim + J] * phiV_x[j * dim + I] * weight; //off-diagonal diffusion

          //data->_jac[ VIrow * nDofsVP + VIcolumn ] += rho * phiV[i] * solVg[J] * phiV_x[j * dim + J] * weight; //diagonal nonlinear
          //data->_jac[ VIrow * nDofsVP + VJcolumn ] += rho * phiV[i] * phiV[j] * SolVg_x[I][J] * weight; //off-diagonal nonlinear
        }
      }

      for(unsigned j = 0; j < nDofsP; j++) {
        unsigned P1column = dim * nDofsV + j;
        unsigned P2column = dim * nDofsV + nDofsP + j;
        data->_jac[VIrow * nDofsVP + P1column] += - phiV_x[i * dim + I] * phiP[j] * weight * weight1; //pressure gradient
        data->_jac[P1column * nDofsVP + VIrow] -= - phiV_x[i * dim + I] * phiP[j] * weight * weight1; //continuity
        data->_jac[VIrow * nDofsVP + P2column] += - phiV_x[i * dim + I] * phiP[j] * weight * weight2; //pressure gradient
        data->_jac[P2column * nDofsVP + VIrow] -= - phiV_x[i * dim + I] * phiP[j] * weight * weight2; //continuity
      }
    }
  }
  for(unsigned i = 0; i < nDofsP; i++) {
    unsigned P1row = dim * nDofsV + i;
    unsigned P2row = dim * nDofsV + nDofsP + i;
    for(unsigned j = 0; j < nDofsP; j++) {
      unsigned P1column = dim * nDofsV + j;
      unsigned P2column = dim * nDofsV + nDofsP + j;
      if(C == 0.)
        data->_jac[P1row * nDofsVP + P1column] += phiP[i] * phiP[j] * weight * eps; // continuity
      if(C == 1.)
        data->_jac[P2row * nDofsVP + P2column] += phiP[i] * phiP[j] * weight * eps; //continuity
    }
  }
}

