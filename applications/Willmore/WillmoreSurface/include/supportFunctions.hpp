#include <boost/math/quaternion.hpp>

double q_dot(const boost::math::quaternion <double> &a, const boost::math::quaternion <double> &b) {
  return a.R_component_1() * b.R_component_1() +
         a.R_component_2() * b.R_component_2() +
         a.R_component_3() * b.R_component_3() +
         a.R_component_4() * b.R_component_4();
}

double operator % (const boost::math::quaternion <double> &a, const boost::math::quaternion <double> &b) { //dot Product overloading
  return a.R_component_1() * b.R_component_1() +
         a.R_component_2() * b.R_component_2() +
         a.R_component_3() * b.R_component_3() +
         a.R_component_4() * b.R_component_4();
}

void  q_set_component_1(boost::math::quaternion <double> &a, const double &w) {
  a = boost::math::quaternion <double> (w, a.R_component_2(), a.R_component_3(), a.R_component_4());
}

void  q_set_component_2(boost::math::quaternion <double> &a, const double &x) {
  a = boost::math::quaternion <double> (a.R_component_1(), x, a.R_component_3(), a.R_component_4());
}

void  q_set_component_3(boost::math::quaternion <double> &a, const double &y) {
  a = boost::math::quaternion <double> (a.R_component_1(), a.R_component_2(), y, a.R_component_4());
}

void  q_set_component_4(boost::math::quaternion <double> &a, const double &z) {
  a = boost::math::quaternion <double> (a.R_component_1(), a.R_component_2(), a.R_component_3(), z);
}

void GetElementNearVertexNumber(MultiLevelSolution &mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* sol = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);

  unsigned envIndex = mlSol.GetIndex("env");
  unsigned angleIndex = mlSol.GetIndex("vAngle");

  unsigned solType = mlSol.GetSolutionType(envIndex);

  sol->_Sol[envIndex]->zero();

  unsigned iproc = msh->processor_id();

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Number of solution element dofs.
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, solType);

    for(unsigned i = 0; i < nDofs; i++) {

      unsigned iDof = msh->GetSolutionDof(i, iel, solType);

      sol->_Sol[envIndex]->add(iDof, 1);

    }
  }

  sol->_Sol[envIndex]->close();

  for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {

    unsigned env = (*sol->_Sol[envIndex])(i);
    double angle = (*sol->_Sol[angleIndex])(i);

    sol->_Sol[angleIndex]->set(i, (2. * M_PI - angle) / env);

  }
  sol->_Sol[angleIndex]->close();

}

// Eugenio's standard FEMuS function.
void CopyDisplacement(MultiLevelSolution &mlSol,  const bool &forward) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);

  unsigned DIM = 3u;
  vector < unsigned > solDxIndex(DIM);
  solDxIndex[0] = mlSol.GetIndex("Dx1");
  solDxIndex[1] = mlSol.GetIndex("Dx2");
  solDxIndex[2] = mlSol.GetIndex("Dx3");

  vector < unsigned > solNDxIndex(DIM);
  solNDxIndex[0] = mlSol.GetIndex("nDx1");
  solNDxIndex[1] = mlSol.GetIndex("nDx2");
  solNDxIndex[2] = mlSol.GetIndex("nDx3");

  if(forward) {
    for(unsigned i = 0; i < DIM; i++) {
      * (solution->_Sol[solNDxIndex[i]]) = * (solution->_Sol[solDxIndex[i]]);
    }
  }
  else {
    for(unsigned i = 0; i < DIM; i++) {
      * (solution->_Sol[solDxIndex[i]]) = * (solution->_Sol[solNDxIndex[i]]);
    }
  }
}

double max(const double &a , const double &b) {
  return (a > b) ? a : b;
}






void GetConformalStructure(std::vector <double> &angle, std::vector < std::vector <double> > &xC) {

  unsigned n = angle.size();

  unsigned imax = n;
  for(unsigned i = 0; i < n; i++) {
    unsigned cnt = 0;
    for(unsigned j = 1; j < n; j++) {
      unsigned ip = (i + j) % n;
      if(angle[i] > angle[ip]) {
        cnt++;
      }
      else {
        break;
      }
    }
    if(cnt == n - 1) {
      imax = i;
      break;
    }
  }
  if(imax == n) {
    double angleSum = 0.;
    for(unsigned i = 0; i < n; i++) {
      angleSum += angle[i];
    }
    double scale = M_PI * (n - 2) / angleSum;
    for(unsigned i = 0; i < n; i++) {
      angle[i] *= scale;
    }
  }
  else {
    double angleSum = 0.;
    for(unsigned j = 1; j < n; j++) {
      unsigned ip = (imax + j) % n;
      angleSum += angle[ip];
    }
    double scale = (M_PI * (n - 2) - angle[imax]) / angleSum;
    for(unsigned j = 1; j < n; j++) {
      unsigned ip = (imax + j) % n;
      angle[ip] *= scale;
    }
  }

  if(n == 3) {
    xC[0][0] = -0.5;
    xC[1][0] = 0.;

    double b = 1;
    double h = b * sin(angle[0]) * sin(angle[1]) / sin(angle[0] + angle[1]);
    double scale = sqrt((sqrt(3.) / 2.) / (b * h));
    b *= scale;
    h *= scale;

    xC[0][1] = xC[0][0] + b;
    xC[1][1] = 0.;

    xC[0][2] = xC[0][0] + h / tan(angle[0]);
    xC[1][2] = h;

    if(xC[0].size() > 3) {
      xC[0][3] = 0.5 * (xC[0][0] + xC[0][1]);
      xC[1][3] = 0.5 * (xC[1][0] + xC[1][1]);

      xC[0][4] = 0.5 * (xC[0][1] + xC[0][2]);
      xC[1][4] = 0.5 * (xC[1][1] + xC[1][2]);

      xC[0][5] = 0.5 * (xC[0][2] + xC[0][0]);
      xC[1][5] = 0.5 * (xC[1][2] + xC[1][0]);
      if(xC[0].size() > 6) {
        xC[0][6] = (xC[0][0] + xC[0][1] + xC[0][2]) / 3.;
        xC[1][6] = (xC[1][0] + xC[1][1] + xC[1][2]) / 3.;
      }
    }
  }
  else if(n == 4) {
    xC[0][0] = 0.;
    xC[1][0] = 0.;

    double a0 = angle[0];
    double a1 = angle[1];
    double a2 = angle[2];
    double a3 = angle[3];

    double a01 = a0 + a1;
    double a12 = a1 + a2;

    double a = 1;
    double b = a * sin(0.5 * a0) * sin(0.5 * a12) / (sin(0.5 * a2) * sin(0.5 * a01));
    double csca3 = 1. / sin(a3);
    double d = csca3 * (b * sin(a2) - a * sin(a12));
    double c = csca3 * (a * sin(a0) - b * sin(a01));

    double scale = sqrt(2. / (a * d * sin(a0) + b * c * sin(a2)));

    a *= scale;
    b *= scale;
    c *= scale;
    d *= scale;

    xC[0][1] = a;
    xC[1][1] = 0.;

    xC[0][2] = (a - b * cos(a1));
    xC[1][2] = (b * sin(a1));

    xC[0][3] = (d * cos(a0));
    xC[1][3] = (d * sin(a0));

    if(xC[0].size() > 4) {
      xC[0][4] = 0.5 * (xC[0][0] + xC[0][1]);
      xC[1][4] = 0.5 * (xC[1][0] + xC[1][1]);

      xC[0][5] = 0.5 * (xC[0][1] + xC[0][2]);
      xC[1][5] = 0.5 * (xC[1][1] + xC[1][2]);

      xC[0][6] = 0.5 * (xC[0][2] + xC[0][3]);
      xC[1][6] = 0.5 * (xC[1][2] + xC[1][3]);

      xC[0][7] = 0.5 * (xC[0][3] + xC[0][0]);
      xC[1][7] = 0.5 * (xC[1][3] + xC[1][0]);

      if(xC[0].size() > 8) {
        xC[0][8] = (xC[0][0] + xC[0][1] + xC[0][2] + xC[0][3]) * 0.25;
        xC[1][8] = (xC[1][0] + xC[1][1] + xC[1][2] + xC[1][3]) * 0.25;
      }
    }
  }
}


void ProjectSolution(MultiLevelSolution& mlSol) {
  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  elem* el = msh->el;

  //unsigned  dim = msh->GetDimension();
  unsigned dim = 2;
  unsigned DIM = 3;
  std::vector < unsigned > indexDx(DIM);
  indexDx[0] = mlSol.GetIndex("Dx1");
  indexDx[1] = mlSol.GetIndex("Dx2");
  indexDx[2] = mlSol.GetIndex("Dx3");
  unsigned solType = mlSol.GetSolutionType(indexDx[0]);

  unsigned iproc = msh->processor_id();

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned nDofs  = msh->GetElementDofNumber(iel, solType);

    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solType);
      unsigned xdof = msh->GetSolutionDof(i, iel, 2);

      double Yhat = (*msh->_topology->_Sol[1])(xdof);
      double Zhat = (*msh->_topology->_Sol[2])(xdof);

      double DY  = (*sol->_Sol[indexDx[1]])(idof);
      double DZ  = (*sol->_Sol[indexDx[2]])(idof);

      double theta = atan2(Zhat + DZ, Yhat + DY);

      sol->_Sol[indexDx[1]]->set(idof, 0.5 * cos(theta) - Yhat);
      sol->_Sol[indexDx[2]]->set(idof, 0.5 * sin(theta) - Zhat);
    }
  }

  sol->_Sol[indexDx[1]]->close();
  sol->_Sol[indexDx[2]]->close();

}


double GetPWillmoreEnergy(MultiLevelSolution &mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* sol  = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  elem* el =  msh->el;

  // Convenience variables to encode the dimension.
  const unsigned dim = 2;
  const unsigned DIM = 3;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Extract the solution vector; get solDx positions in the ml_sol object.
  unsigned solDxIndex[DIM];
  solDxIndex[0] = mlSol.GetIndex("Dx1");
  solDxIndex[1] = mlSol.GetIndex("Dx2");
  solDxIndex[2] = mlSol.GetIndex("Dx3");

  // Extract the finite element type for solx.
  unsigned solxType;
  solxType = mlSol.GetSolutionType(solDxIndex[0]);

  // Define solx and solxOld.
  std::vector < double > solx[DIM];
  std::vector < double > solxOld[DIM];

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get positions of Y in the ml_sol object.
  unsigned solYIndex[DIM];
  solYIndex[0] = mlSol.GetIndex("Y1");
  solYIndex[1] = mlSol.GetIndex("Y2");
  solYIndex[2] = mlSol.GetIndex("Y3");

  // Extract the finite element type for Y.
  unsigned solYType;
  solYType = mlSol.GetSolutionType(solYIndex[0]);


  // Define solY and solYOld.
  std::vector < double > solY[DIM];
  std::vector < double > solYOld[DIM];

  // Initialize area, volume, P-Willmore energy.

  double energy = 0.;

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Number of solution element dofs.
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nxDofs  = msh->GetElementDofNumber(iel, solxType);
    unsigned nYDofs  = msh->GetElementDofNumber(iel, solYType);

    // Resize solution vectors.
    for(unsigned K = 0; K < DIM; K++) {
      solx[K].resize(nxDofs);
      solxOld[K].resize(nxDofs);
      solY[K].resize(nYDofs);
      solYOld[K].resize(nYDofs);
    }


    // Loop which handles local storage of global mapping and solution X.
    for(unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof(i, iel, solxType);
      unsigned iXDof  = msh->GetSolutionDof(i, iel, xType);

      for(unsigned K = 0; K < DIM; K++) {
        solxOld[K][i] = (*msh->_topology->_Sol[K])(iXDof) + (*sol->_SolOld[solDxIndex[K]])(iDDof);
        solx[K][i] = (*msh->_topology->_Sol[K])(iXDof) + (*sol->_Sol[solDxIndex[K]])(iDDof);
      }
    }

    // Loop which handles local storage of global mapping and solution Y.
    for(unsigned i = 0; i < nYDofs; i++) {

      // Global-to-local mapping between solution node and solution dof.
      unsigned iYDof = msh->GetSolutionDof(i, iel, solYType);
      for(unsigned K = 0; K < DIM; K++) {

        // Global-to-local solutions.
        solYOld[K][i] = (*sol->_SolOld[solYIndex[K]])(iYDof);
        solY[K][i] = (*sol->_Sol[solYIndex[K]])(iYDof);

      }
    }

    // begin GAUSS POINT LOOP
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solxType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function
      const double *phix_uv[dim]; // first order derivatives in (u,v)

      const double *phiY;  // local test function


      double weight; // gauss point weight

      //Extract Gauss point weight, test functions, and their partial derivatives.
      // "0" is derivative in u, "1" is derivative in v.
      weight = msh->_finiteElement[ielGeom][solxType]->GetGaussWeight(ig);

      phix = msh->_finiteElement[ielGeom][solxType]->GetPhi(ig);
      phix_uv[0] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDXi(ig);
      phix_uv[1] = msh->_finiteElement[ielGeom][solxType]->GetDPhiDEta(ig);

      phiY = msh->_finiteElement[ielGeom][solYType]->GetPhi(ig);


      double solYg[3] = {0., 0., 0.};
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nYDofs; i++) {
          solYg[K] += phiY[i] * 0.5 * (0.* solYOld[K][i] + 2.* solY[K][i]);
        }
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * 0.5 * ( 0. * solxOld[K][i] +  2. * solx[K][i]);
          }
        }
      }

      // Computing the metric, metric determinant, and area element.
      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      double Area = weight * sqrt(detg);

      // Computing the unit normal vector N.
      double normal[DIM];
      normal[0] = normalSign * (solx_uv[1][0] * solx_uv[2][1]
                                - solx_uv[2][0] * solx_uv[1][1]) / sqrt(detg);
      normal[1] = normalSign * (solx_uv[2][0] * solx_uv[0][1]
                                - solx_uv[0][0] * solx_uv[2][1]) / sqrt(detg);
      normal[2] = normalSign * (solx_uv[0][0] * solx_uv[1][1]
                                - solx_uv[1][0] * solx_uv[0][1]) / sqrt(detg);

      // Computing Y.N and |Y|^2, which are essentially 2H and 4H^2.
      double YdotN = 0.;
      double YdotY = 0.;
      for(unsigned K = 0; K < DIM; K++) {
        YdotN += solYg[K] * normal[K];
        YdotY += solYg[K] * solYg[K];
      }
      // double signYdotN = (YdotN.value() >= 0.) ? 1. : -1.;
      double signYdotN = 1.;

      // Some necessary quantities when working with polynomials.
      double sumP3 = 0.;
      for(unsigned p = 0; p < 3; p++) {
        double signP = (P[p] % 2u == 0) ? 1. : signYdotN;
        sumP3 += signP * ap[p] * pow(YdotY, P[p] / 2.);
      }

      energy += sumP3 * Area;

    } // end GAUSS POINT LOOP.


  } // End ELEMENT LOOP for each process.


  double energyAll;
  MPI_Reduce(&energy, &energyAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return energyAll;

}
