
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"

//2D NONLOCAL EX : nonlocal diffusion for a body with different material properties

#include <vector>
#include <cmath>

using namespace femus;

#define N_UNIFORM_LEVELS  1
#define N_ERASED_LEVELS   0

#define EX_1       -1.
#define EX_2        1.
#define EY_1       -1.
#define EY_2        1.

void GetNormalQuad(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R,  std::vector<double> &a, double &d, bool &cut);

int main(int argc, char** argv) {

  const std::string fe_quad_rule_1 = "seventh";
  const std::string fe_quad_rule_2 = "eighth";

  // ======= Init ========================
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  unsigned numberOfUniformLevels = N_UNIFORM_LEVELS;


  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;

  mlMsh.GenerateCoarseBoxMesh(50, 50, 0, EX_1, EX_2, EY_1, EY_2, 0., 0., QUAD9, fe_quad_rule_1.c_str());
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);

  // erase all the coarse mesh levels
  const unsigned erased_levels = N_ERASED_LEVELS;
  mlMsh.EraseCoarseLevels(erased_levels);

  const unsigned level = N_UNIFORM_LEVELS - N_ERASED_LEVELS - 1;

  Mesh*                    msh = mlMsh.GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)

  unsigned xType = 2;

  const unsigned  dim = msh->GetDimension();

  std::vector < std::vector < double > > x1;

  FILE * fp;

  fp = fopen ("lines.dat", "w");


  for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    unsigned nDof = 4; //TODO msh->GetElementDofNumber(iel, xType);  // number of coordinate element dofs
    x1.resize(nDof);
    std::vector < std::vector < double > > x1(nDof);
    for(unsigned i = 0; i < nDof; i++) {
      x1[i].resize(dim);
    }

    for(unsigned i = 0; i < nDof; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        x1[i][k] = (*msh->_topology->_Sol[k])(xDof);  // global extraction and local storage for the element coordinates
      }
    }

    std::vector<double> xg(2, 0);
    xg[0] -= 0.2;
    xg[1] -= 0.2;
    double R = 0.5;
    std::vector<double> a;
    double d;
    bool cut;

    GetNormalQuad(x1, xg, R, a, d, cut);

    /* trivial print for xmgrace */
    if(cut) {

      double xx = 0.;
      double yy = 0.;
      for( unsigned j = 0; j < 50; j++ ) {
        xx = x1[0][0] + ( j * (x1[1][0] - x1[0][0]) / 50 );
        yy = - ( a[0] / a[1] ) * xx - ( d / a[1] );
        fprintf(fp, "%f %f \n", xx, yy );
      }

      fprintf(fp, "\n \n");
    }
  }

  /* Basic numerical tests for the circle - no mesh involved*/
  std::vector < std::vector<double> > xva = {{1., 1.}, {2., 1.}, {2., 2.}, {1., 2.}};
  std::vector < std::vector<double> > xvb =  {{2., 1.}, {2., 2.}, {1., 2.}, {1., 1.}};
  std::vector < std::vector<double> > xvc = {{2., 2.}, {1., 2.}, {1., 1.}, {2., 1.}};
  std::vector < std::vector<double> > xvd = {{1., 2.}, {1., 1.}, {2., 1.}, {2., 2.}};
  std::vector < std::vector<double> > xve = {{-1., -0.5}, {-1., 0.5}, {-2., .5}, {-2., -0.5}};
  std::vector<double> xg(2, 0);
  double R = 2.;
  std::vector<double> a;
  double d;
  bool cut;
  GetNormalQuad(xva, xg, R, a, d, cut);
  GetNormalQuad(xvb, xg, R, a, d, cut);
  GetNormalQuad(xvc, xg, R, a, d, cut);
  GetNormalQuad(xvd, xg, R, a, d, cut);
  GetNormalQuad(xve, xg, R, a, d, cut);

  xg[0] -= 2;
  xg[1] -= 2;

  for(unsigned i = 0; i < xva.size(); i++) {
    for(unsigned k = 0; k < xva[i].size(); k++) {
      xva[i][k] -= 2;
      xvb[i][k] -= 2;
      xvc[i][k] -= 2;
      xvd[i][k] -= 2;
      xve[i][k] -= 2;
    }
  }
  GetNormalQuad(xva, xg, R, a, d, cut);
  GetNormalQuad(xvb, xg, R, a, d, cut);
  GetNormalQuad(xvc, xg, R, a, d, cut);
  GetNormalQuad(xvd, xg, R, a, d, cut);
  GetNormalQuad(xve, xg, R, a, d, cut);

  fclose(fp);

  return 1;
}


void GetNormalQuad(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R,  std::vector<double> &a, double &d, bool &cut) {

  unsigned nve =  xv.size();
  unsigned dim =  xv[0].size();

  std::vector<double> dist(nve, 0);
  for(unsigned i = 0; i < nve; i++) {
    for(unsigned k = 0;  k < dim; k++) {
      dist[i] += (xv[i][k] - xg[k]) * (xv[i][k] - xg[k]);
    }
    dist[i] = sqrt(dist[i]) - R;
    if(dist[i] == 0.) dist[i] = 1.0e-10;
  }

  std::vector <double> theta(2);
  unsigned cnt = 0;
  for(unsigned e = 0; e < nve; e++) {
    unsigned ep1 = (e + 1) % nve;
    if(dist[e] * dist[ep1] < 0) {

      double s = 0.5  * (1 + (dist[e] + dist[ep1]) / (dist[e] - dist[ep1]));
      theta[cnt] = atan2((1 - s) * xv[e][1] + s * xv[ep1 ][1]  - xg[1], (1 - s) * xv[e][0] + s * xv[ep1][0] - xg[0]) ;
      cnt++;

    }

  }
  if(cnt == 0 ) {
    cut = 0;
//     std::cout << "Empty/Full cell\n";
    return;
  }
  else {

    cut = 1;

    if(theta[0] > theta[1]) {
      std::swap(theta[0], theta[1]);
    }

    double DT = theta[1] - theta[0];
    if(DT > M_PI) {
      std::swap(theta[0], theta[1]);
      theta[1] += 2. * M_PI;
      DT = theta[1] - theta[0];
    }

// std::cout << theta[0] / M_PI * 180 << " " << theta[1] / M_PI * 180 << " " << DT / M_PI * 180 << std::endl;

    std::vector < double > xm(dim);


    d = R * sqrt(0.5 * DT / tan(0.5 * DT)) ;
    a.resize(dim);
    a[0] = -cos(theta[0] + 0.5 * DT);
    a[1] = -sin(theta[0] + 0.5 * DT);

    xm[0] = (a[0] == 0) ? 0 : 0.5 * (-d / a[0]);

    for(unsigned k = 0; k < dim; k++) {
      xm[k] = -a[k] * d;
      xm[k] += xg[k];
    }
    d += - a[0] * xg[0] - a[1] * xg[1];

    std::cout << "xm = " << xm[0] << " " << xm[1] << std::endl;
    std::cout << "a = " << a[0] << " b = " << a[1] << " d = " << d << std::endl;

  }

}



