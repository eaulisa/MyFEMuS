
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "CurrentElem.hpp"
#include "LinearImplicitSystem.hpp"

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

double InitialValueU(const std::vector < double >& x)
{
  return 0. * x[0] * x[0];
}
bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time)
{
  bool dirichlet = true; //dirichlet
  value = 0.;

//   if(facename == 1) {
//     dirichlet = true; //dirichlet
//     value = 0.;
//   }
//   else if(facename == 2) {
//     dirichlet = true; //dirichlet
//     value = 0.;
//   }

  return dirichlet;
}

void GetNormalQuad(const std::vector < std::vector<double> > &xv, const std::vector<double> &xg, const double &R,  std::vector<double> &a, double &d, bool &cut);
void SimpleNonlocalAssembly (MultiLevelProblem& ml_prob);

int main(int argc, char** argv) {

  const std::string fe_quad_rule_1 = "seventh";
  const std::string fe_quad_rule_2 = "eighth";

  // ======= Init ========================
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  unsigned numberOfUniformLevels = N_UNIFORM_LEVELS;


  MultiLevelMesh mlMsh;
    double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.GenerateCoarseBoxMesh(20, 20, 0, EX_1, EX_2, EY_1, EY_2, 0., 0., QUAD9, fe_quad_rule_1.c_str());
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);

  // erase all the coarse mesh levels
  const unsigned erased_levels = N_ERASED_LEVELS;
  mlMsh.EraseCoarseLevels(erased_levels);

  const unsigned level = N_UNIFORM_LEVELS - N_ERASED_LEVELS - 1;
  
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("u", LAGRANGE, SECOND, 2);


  mlSol.Initialize("All");
  mlSol.Initialize("u", InitialValueU);

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("All");

  MultiLevelProblem ml_prob(&mlSol);
  
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("FracProblem");


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
  
  /*Testing the function GetNormalQuad inside a nonlocal assembly-like function*/
  SimpleNonlocalAssembly(ml_prob);

//   /* Basic numerical tests for the circle - no mesh involved*/
//   std::vector < std::vector<double> > xva = {{1., 1.}, {2., 1.}, {2., 2.}, {1., 2.}};
//   std::vector < std::vector<double> > xvb =  {{2., 1.}, {2., 2.}, {1., 2.}, {1., 1.}};
//   std::vector < std::vector<double> > xvc = {{2., 2.}, {1., 2.}, {1., 1.}, {2., 1.}};
//   std::vector < std::vector<double> > xvd = {{1., 2.}, {1., 1.}, {2., 1.}, {2., 2.}};
//   std::vector < std::vector<double> > xve = {{-1., -0.5}, {-1., 0.5}, {-2., .5}, {-2., -0.5}};
//   std::vector<double> xg(2, 0);
//   double R = 2.;
//   std::vector<double> a;
//   double d;
//   bool cut;
//   GetNormalQuad(xva, xg, R, a, d, cut);
//   GetNormalQuad(xvb, xg, R, a, d, cut);
//   GetNormalQuad(xvc, xg, R, a, d, cut);
//   GetNormalQuad(xvd, xg, R, a, d, cut);
//   GetNormalQuad(xve, xg, R, a, d, cut);
// 
//   xg[0] -= 2;
//   xg[1] -= 2;
// 
//   for(unsigned i = 0; i < xva.size(); i++) {
//     for(unsigned k = 0; k < xva[i].size(); k++) {
//       xva[i][k] -= 2;
//       xvb[i][k] -= 2;
//       xvc[i][k] -= 2;
//       xvd[i][k] -= 2;
//       xve[i][k] -= 2;
//     }
//   }
//   GetNormalQuad(xva, xg, R, a, d, cut);
//   GetNormalQuad(xvb, xg, R, a, d, cut);
//   GetNormalQuad(xvc, xg, R, a, d, cut);
//   GetNormalQuad(xvd, xg, R, a, d, cut);
//   GetNormalQuad(xve, xg, R, a, d, cut);

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


void SimpleNonlocalAssembly (MultiLevelProblem& ml_prob) {

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("FracProblem");

  const unsigned level = N_UNIFORM_LEVELS - N_ERASED_LEVELS - 1;


  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension();
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  constexpr unsigned int space_dim = 3;

  unsigned    iproc = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  vector < vector < double > > x1(dim);    // local coordinates
  vector < vector < double > > x2(dim);    // local coordinates
  for(unsigned k = 0; k < dim; k++) {
    x1[k].reserve(maxSize);
    x2[k].reserve(maxSize);
  }

  std::vector < std::vector < double > >  Jac_qp(dim);
  std::vector < std::vector < double > >  JacI_qp(space_dim);
  double detJac_qp;

  CurrentElem < double > geom_element1(dim, msh);            // must be adept if the domain is moving, otherwise double
  CurrentElem < double > geom_element2(dim, msh);

  unsigned soluIndex = mlSol->GetIndex("u");
  unsigned solType   = mlSol->GetSolutionType(soluIndex);

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < std::vector < /*const*/ elem_type_templ_base< double, double > *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
  
  FILE * fp1;

  fp1 = fopen ("lines_qp.dat", "w");

  for(int kproc = 0; kproc < nprocs; kproc++) {
    for(int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++) {
      short unsigned ielGeom2;
      unsigned nDof2;
      unsigned nDofx2;
      unsigned nDofLin;
      unsigned nDofu2;
      //unsigned n_face;

      if(iproc == kproc) {
        ielGeom2 = msh->GetElementType(jel);
        nDof2  = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs
        nDofx2 = msh->GetElementDofNumber(jel, xType);    // number of coordinate element dofs
        nDofLin = 4; //TODO

      }

      MPI_Bcast(&ielGeom2, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDof2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDofx2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      //MPI_Bcast(&n_face, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);

      for(int k = 0; k < dim; k++) {
        x2[k].resize(nDofx2);
      }
      vector < vector < double > > x2i(nDofLin);
      for(unsigned j = 0; j < nDofLin; j++) {
        x2i[j].resize(dim);
      }

      vector < double > phi;
      vector < double > phi_x;

      phi.reserve(maxSize);
      phi_x.reserve(maxSize * dim);


      // local storage of coordinates  #######################################
      if(iproc == kproc) {
        for(unsigned j = 0; j < nDofx2; j++) {
          unsigned xDof  = msh->GetSolutionDof(j, jel, xType);  // global to global mapping between coordinates node and coordinate dof
          for(unsigned k = 0; k < dim; k++) {
            x2[k][j] = (*msh->_topology->_Sol[k])(xDof);  // global extraction and local storage for the element coordinates
          }
        }
        for(unsigned j = 0; j < nDofLin; j++) {
           for(unsigned k = 0; k < dim; k++) {
              x2i[j][k] = x2[k][j];
        } 
        }
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& x2[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }

      if(iproc == kproc) {
        geom_element2.set_coords_at_dofs_and_geom_type(jel, xType);
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& geom_element2.get_coords_at_dofs()[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
      for(unsigned k = 0; k < space_dim; k++) {
        MPI_Bcast(& geom_element2.get_coords_at_dofs_3d()[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }

      const unsigned jgNumber = msh->_finiteElement[ielGeom2][solType]->GetGaussPointNumber();
//       const unsigned jgNumber = ml_prob.GetQuadratureRule(ielGeom2).GetGaussPointsNumber();

      vector < vector < double > > xg2(jgNumber);
      vector <double> weight2(jgNumber);
      vector < vector <double> > phi2(jgNumber);  // local test function

      for(unsigned jg = 0; jg < jgNumber; jg++) {

        msh->_finiteElement[ielGeom2][solType]->Jacobian(x2, jg, weight2[jg], phi2[jg], phi_x);

//         elem_all[ielGeom2][xType]->JacJacInv(/*x2*/geom_element2.get_coords_at_dofs_3d(), jg, Jac_qp, JacI_qp, detJac_qp, space_dim);
//         weight2[jg] = detJac_qp * ml_prob.GetQuadratureRule(ielGeom2).GetGaussWeightsPointer()[jg];
//         elem_all[ielGeom2][solType]->shape_funcs_current_elem(jg, JacI_qp, phi2[jg], phi_x /*boost::none*/, boost::none /*phi_u_xx*/, space_dim);




        xg2[jg].assign(dim, 0.);

        for(unsigned j = 0; j < nDof2; j++) {
          for(unsigned k = 0; k < dim; k++) {
            xg2[jg][k] += x2[k][j] * phi2[jg][j];
          }
        }
      }

      std::vector <int> bd_face(0);
      unsigned nFaces;

      for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        short unsigned ielGeom1 = msh->GetElementType(iel);
        unsigned nDof1  = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
        unsigned nDofx1 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

        for(int k = 0; k < dim; k++) {
          x1[k].resize(nDofx1);
        }

        // local storage of coordinates
        for(unsigned i = 0; i < nDofx1; i++) {
          unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
          for(unsigned k = 0; k < dim; k++) {
            x1[k][i] = (*msh->_topology->_Sol[k])(xDof);
          }
        }

        const unsigned igNumber = msh->_finiteElement[ielGeom1][solType]->GetGaussPointNumber();
        double weight1;
        vector < double > phi1;  // local test function


        for(unsigned ig = 0; ig < igNumber; ig++) {
          msh->_finiteElement[ielGeom1][solType]->Jacobian(x1, ig, weight1, phi1, phi_x);
          vector < double > xg1(dim, 0.);
          for(unsigned i = 0; i < nDof1; i++) {
            for(unsigned k = 0; k < dim; k++) {
              xg1[k] += x1[k][i] * phi1[i];
            }
          }
          double R = 0.5;
          std::vector<double> a;
          double d;
          bool cut;
          GetNormalQuad(x2i, xg1, R, a, d, cut);

          if(cut && iel == 110 && ig == 0) { //TODO
            double xx = 0.;
            double yy = 0.;
            for( unsigned j = 0; j < 50; j++ ) {
              xx = x2i[0][0] + ( j * (x2i[1][0] - x2i[0][0]) / 50 );
              yy = - ( a[0] / a[1] ) * xx - ( d / a[1] );
              fprintf(fp1, "%f %f \n", xx, yy );
            }

            fprintf(fp1, "\n \n");
          }

        }


      }
    }
  }
  fclose(fp1);
}



