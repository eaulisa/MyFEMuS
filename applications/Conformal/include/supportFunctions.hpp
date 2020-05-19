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

  unsigned DIM = 3u;
  unsigned solIndex = mlSol.GetIndex("ENVN");
  unsigned angleIndex = mlSol.GetIndex("bAngle");

  unsigned solType = mlSol.GetSolutionType(solIndex);

  sol->_Sol[solIndex]->zero();

  unsigned iproc = msh->processor_id();

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // Number of solution element dofs.
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, solType);

    for(unsigned i = 0; i < nDofs; i++) {

      unsigned iDof = msh->GetSolutionDof(i, iel, solType);

      sol->_Sol[solIndex]->add(iDof, 1);

    }
  }

  sol->_Sol[solIndex]->close();

  for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {

    unsigned nve = (*sol->_Sol[solIndex])(i);
    double angle = (*sol->_Sol[angleIndex])(i);

    sol->_Sol[angleIndex]->set(i, (2. * M_PI - angle) / nve);

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


void ChangeTriangleConfiguration1(const std::vector<unsigned> & ENVN, std::vector <double> &angle) {
  double scale;
  if(ENVN[0] < ENVN[1] && ENVN[0] < ENVN[2]) {
    scale = (M_PI - angle[0]) / (angle[1] + angle [2]);
    angle[1] *= scale;
    angle[2] *= scale;
  }
  else if(ENVN[0] < ENVN[1] && ENVN[0] == ENVN[2]) {
    angle[1] = M_PI - 2. * angle[0];
  }
  else if(ENVN[0] <= ENVN[1]  && ENVN[0] > ENVN[2]) {
    scale = (M_PI - angle[2]) / (angle[1] + angle [0]);
    angle[1] *= scale;
    angle[0] *= scale;
  }
  else if(ENVN[0] == ENVN[1] && ENVN[0] < ENVN[2]) {
    angle[2] = M_PI - 2. * angle[0];
  }
  else if(ENVN[0] == ENVN[1] && ENVN[0] == ENVN[2]) {
    angle[0] = angle[1] = angle[2] =  M_PI / 3.;
  }
  else if(ENVN[0] > ENVN[1] && ENVN[0] <= ENVN[2]) {
    scale = (M_PI - angle[1]) / (angle[0] + angle [2]);
    angle[0] *= scale;
    angle[2] *= scale;
  }
  else if(ENVN[0] > ENVN[1] && ENVN[0] > ENVN[2]) {
    if(ENVN[1] < ENVN[2]) {
      scale = (M_PI - angle[1]) / (angle[0] + angle [2]);
      angle[0] *= scale;
      angle[2] *= scale;
    }
    else if(ENVN[1] == ENVN[2]) {
      angle[0] = M_PI - 2. * angle[1];
    }
    else if(ENVN[1] > ENVN[2]) {
      scale = (M_PI - angle[2]) / (angle[0] + angle [1]);
      angle[0] *= scale;
      angle[1] *= scale;
    }
  }
}


void ChangeTriangleConfiguration2(const std::vector<unsigned> & ENVN, std::vector <double> &angle) {
  unsigned type = 3; // there are 2 or 3 leading angles
  if(ENVN[0] < ENVN[1]) {  // 0 leads on 1
    if(ENVN[0] < ENVN[2]) type = 0;  // 0 is leading angle
    else if(ENVN[0] > ENVN[2]) type = 2;  // 2 is leading angle
  }
  else if(ENVN[0] > ENVN[1]) {  // 1 leads on 0
    if(ENVN[1] < ENVN[2]) type = 1;  // 1 is leading angle
    else if(ENVN[1] > ENVN[2]) type = 2;  // 2 is leading angle
  }
  else { // 0 equals 1
    if(ENVN[0] > ENVN[2]) type = 2;  // 2 is leading angle
  }

  double scale;
  if(type == 0) {
    scale = (M_PI - angle[0]) / (angle[1] + angle [2]);
    angle[1] *= scale;
    angle[2] *= scale;
  }
  else if(type == 1) {
    scale = (M_PI - angle[1]) / (angle[0] + angle [2]);
    angle[0] *= scale;
    angle[2] *= scale;
  }
  else if(type == 2) {
    scale = (M_PI - angle[2]) / (angle[1] + angle [0]);
    angle[1] *= scale;
    angle[0] *= scale;
  }
  else {
    scale = M_PI / (angle[0] + angle[1] + angle[2]);
    angle[0] *= scale;
    angle[1] *= scale;
    angle[2] *= scale;
  }
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
    xC[1][1] = 0.;

    double l = 1;
    double d = l * sin(angle[0]) * sin(angle[1]) / sin(angle[0] + angle[1]);
    double scale = sqrt((sqrt(3.) / 2.) / (l * d));
    l = l * scale;
    d = d * scale;

    xC[0][1] = xC[0][0] + l;
    xC[0][2] = xC[0][0] + d / tan(angle[0]);
    xC[1][2] = d;

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
    xC[1][1] = 0.;

    double a0 = angle[0];
    double a1 = angle[1];
    double a2 = angle[2];

    double a01 = a0 + a1;
    double a12 = a1 + a2;

    double a = 1;
    double b;
    if(fabs(a01 - M_PI) < 1.0e-12) {
      b = 1. / sin(a0) + 0.5 / sin(a2) * sin(a12);
    }
    else if(fabs(2 * angle[2] - M_PI) < 1.0e-12) {

      b = 1. / tan(a01) *
          (1. / cos(a01) * sin(a0) +
           sqrt((2.* cos(a1) * sin(a0) + (2. * cos(a0) - sin(a0)) * sin(a1)) / cos(a01)));
    }
    double csca3 = 1. / sin(angle[3]);
    double d = csca3 * (b * sin(angle[2]) - a * sin(angle[1] + angle[2]));

//     double a = 1.;
//     double csc3 = 1. / sin(angle[3]);
//     double den = sin(angle[1]) + csc3 * sin(angle[0]) * sin(angle[2]);
//     double b = (2. + csc3 * sin(angle[0]) * sin(angle[1] + angle[2])) / den;
//     double d = csc3 * (2. * sin(angle[2]) - sin(angle[1]) * sin(angle[1] + angle[2])) / den;
//     double c1 = (1. - b * cos(angle[1]) - d * cos(angle[0]));
//     double c2 = (b * sin(angle[1]) - d * sin(angle[0]));
//     double c = sqrt(c1 * c1 + c2 * c2);
//     double scale = 2. / ( a * d * sin(angle[0]) +  b * c *sin(angle[2]));
//
//
//     xC[0][1] =  a * scale;
//     xC[0][2] = (a - b * cos(angle[1])) * scale;
//     xC[1][2] = (b * sin(angle[1])) * scale;
//     xC[0][3] = (d * cos(angle[0])) * scale;
//     xC[1][3] = (d * sin(angle[0])) * scale;

    if(xC[0].size() > 4) {
      xC[0][4] = 0.5 * (xC[0][0] + xC[0][1]);
      xC[1][4] = 0.5 * (xC[1][0] + xC[1][1]);

      xC[0][5] = 0.5 * (xC[0][1] + xC[0][2]);
      xC[1][5] = 0.5 * (xC[1][1] + xC[1][2]);

      xC[0][6] = 0.5 * (xC[0][2] + xC[0][3]);
      xC[1][6] = 0.5 * (xC[1][2] + xC[1][3]);

      xC[0][7] = 0.5 * (xC[0][3] + xC[0][0]);
      xC[1][7] = 0.5 * (xC[1][3] + xC[1][0]);

      if(xC[0].size() > 7) {
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

