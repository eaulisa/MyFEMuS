#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "MultiLevelSolution.hpp"

#include "MeshRefinement.hpp"

const double WEIGHT[6][27] = {
  {
    1. / 64., 1. / 64., 1. / 64., 1. / 64.,
    1. / 64., 1. / 64., 1. / 64., 1. / 64.,
    1. / 32., 1. / 32., 1. / 32., 1. / 32.,
    1. / 32., 1. / 32., 1. / 32., 1. / 32.,
    1. / 32., 1. / 32., 1. / 32., 1. / 32.,
    1. / 16., 1. / 16., 1. / 16., 1. / 16., 1. / 16., 1. / 16., 1. / 8
  }, //hex
  {
    1. / 32., 1. / 32., 1. / 32., 1. / 32.,
    7. / 48., 7. / 48., 7. / 48., 7. / 48., 7. / 48., 7. / 48.
  }, //tet
  {
    1. / 48., 1. / 48., 1. / 48., 1. / 48., 1. / 48., 1. / 48.,
    1. / 16., 1. / 16., 1. / 16., 1. / 16., 1. / 16., 1. / 16.,
    1. / 24., 1. / 24., 1. / 24., 1. / 8., 1. / 8., 1. / 8.
  }, //wedge
  {
    0.0625, 0.0625, 0.0625, 0.0625,
    0.125, 0.125, 0.125, 0.125, 0.25
  },
  {
    1. / 12., 1. / 12., 1. / 12.,
    0.25, 0.25, 0.25, 0.
  },
  {0.25, 0.25, .5}
} ;


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int facename, const double t) {
  bool test = 1; //dirichlet
  value = 0.;
  return test;
}

void BuildMarkers(MultiLevelMesh& mlMesh, const double &dminCut, const double &dmaxCut);
void FlagElements(MultiLevelMesh& mlMesh, const unsigned &layers);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  //mlMsh.ReadCoarseMesh("../input/beam.neu", "fifth", scalingFactor);
  mlMsh.ReadCoarseMesh("../input/3dbeam.neu", "fifth", scalingFactor);
  //mlMsh.ReadCoarseMesh("../input/blades.neu", "fifth", scalingFactor);
  //mlMsh.ReadCoarseMesh("../input/mindcraft_valve.neu", "fifth", scalingFactor);


  //mlMsh.RefineMesh(2, 2, NULL);

  unsigned numberOfRefinement = 2;

  for(unsigned i = 0; i < numberOfRefinement; i++) {
    FlagElements(mlMsh, 1);
    mlMsh.AddAMRMeshLevel();
  }

  mlMsh.EraseCoarseLevels(numberOfRefinement);

  unsigned dim = mlMsh.GetDimension();

  FEOrder femOrder = FIRST;

  MultiLevelSolution mlSol(&mlMsh);

  //add variables to mlSol
  mlSol.AddSolution("u", LAGRANGE, femOrder, 2);

  mlSol.Initialize("All");
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("u", "Steady");

  BuildMarkers(mlMsh, -0.6, 1.0E10);

  //******* Print solution *******
  mlSol.SetWriter(VTK);
  std::vector<std::string> mov_vars;
  std::vector<std::string> print_vars;
  print_vars.push_back("All");
  mlSol.GetWriter()->SetDebugOutput(false);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);

  return 0;

} //end main




void BuildMarkers(MultiLevelMesh& mlMesh, const double &dminCut, const double &dmaxCut) {

  unsigned level = mlMesh.GetNumberOfLevels() - 1;
  Mesh *msh   = mlMesh.GetLevel(level);
  unsigned iproc  = msh->processor_id();
  const unsigned dim = msh->GetDimension();

  unsigned solType = 2;
  unsigned xpSize = msh->el->GetNodeNumber();



  //BEGIN bulk markers
  std::vector < std::vector < double > > xp(xpSize);
  std::vector < double > dist(xpSize);
  
  for(unsigned i = 0; i < xpSize; i++) {
    xp[i].resize(dim);
    dist[i] = (msh->GetSolidMark(i) > 0.5) ? 1. : -1.;
    for(unsigned k = 0; k < dim; k++) {
      xp[i][k] = (*msh->_topology->_Sol[k])(i);
    }
  }
    
  std::vector < double > weight(xpSize, 0.);
 
  std::vector < double > phi;
  std::vector < double > gradPhi;
  std::vector <std::vector < double> > vx(dim);

  double totalArea = 0.;
  for(unsigned iel = 0; iel < msh->el->GetElementNumber(); iel++) {
    short unsigned ielt = msh->GetElementType(iel);
    unsigned nDofs = msh->GetElementDofNumber(iel, solType);
    for(unsigned  k = 0; k < dim; k++) {
      vx[k].resize(nDofs);
    }

    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idofX = msh->GetSolutionDof(i, iel, solType);
      for(unsigned  k = 0; k < dim; k++) {
        vx[k][i] = (*msh->_topology->_Sol[k])(idofX);
      }
    }
    double ielArea = 0.;
    for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {
      double jac;
      msh->_finiteElement[ielt][solType]->Jacobian(vx, ig, jac, phi, gradPhi);
      ielArea += jac;
    }
    for(unsigned i = 0; i < nDofs; i++) {
      weight[msh->GetSolutionDof(i, iel, solType)] += ielArea * WEIGHT[ielt][i];
      totalArea += ielArea * WEIGHT[ielt][i];
    }
  }
  //END bulk markers

  std::cout << totalArea << " " << 2 * 4.99 + 0.5 * M_PI * 1. * 1. << std::endl;

  //BEGIN interface markers
  std::vector < std::vector < double > > xpI(0);
  std::vector < std::vector < std::vector < double > > > xpT(0);

  xpI.reserve(xpSize);
  xpT.reserve(xpSize);

  totalArea = 0.;

  for(unsigned iel = 0; iel < msh->el->GetElementNumber(); iel++) {
    unsigned ielMat = msh->GetElementMaterial(iel);

    if(ielMat == 4) {
      unsigned nDofs = msh->GetElementDofNumber(iel, solType);
      for(unsigned  k = 0; k < dim; k++) {
        vx[k].resize(nDofs);
      }
      for(unsigned i = 0; i < nDofs; i++) {
        unsigned idofX = msh->GetSolutionDof(i, iel, solType);
        for(unsigned  k = 0; k < dim; k++) {
          vx[k][i] = (*msh->_topology->_Sol[k])(idofX);
        }
      }

      for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
        int jel = msh->el->GetFaceElementIndex(iel, iface) - 1;
        if(jel >= 0) { // iface is not a boundary of the domain
          unsigned jelMat = msh->GetElementMaterial(jel);
          if(ielMat != jelMat) { //iel and jel are on the FSI interface

            const unsigned faceGeom = msh->GetElementFaceType(iel, iface);
            unsigned faceDofs = msh->GetElementFaceDofNumber(iel, iface, solType);
            std::vector  < std::vector  <  double > > faceVx(dim);    // A matrix holding the face coordinates rowwise.
            for(int k = 0; k < dim; k++) {
              faceVx[k].resize(faceDofs);
            }
            for(unsigned i = 0; i < faceDofs; i++) {
              unsigned inode = msh->GetLocalFaceVertexIndex(iel, iface, i);    // face-to-element local node mapping.
              for(unsigned k = 0; k < dim; k++) {
                faceVx[k][i] =  vx[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
              }
            }

            double faceArea = 0.;
            for(unsigned ig = 0; ig < msh->_finiteElement[faceGeom][solType]->GetGaussPointNumber(); ig++) {
              double jac;
              std::vector < double> normal;
              msh->_finiteElement[faceGeom][solType]->JacobianSur(faceVx, ig, jac, phi, gradPhi, normal);
              faceArea += jac;
            }

            for(unsigned i = 0; i < faceDofs; i++) {

              unsigned ip =  xpI.size();
              xpI.resize(ip + 1);
              xpI[ip].resize(dim, 0.);

              for(unsigned k = 0; k < dim; k++) {
                xpI[ip][k] = faceVx[k][i];
              }

              std::vector<double> xi(dim - 1); //local coordinate of the i-faceDof
              for(int k = 0; k < dim - 1; k++) {
                xi[k] = *(msh->_finiteElement[faceGeom][solType]->GetBasis()->GetXcoarse(i) + k);
              }

              std::vector < std::vector < double > > dphidu(dim - 1); // first order derivatives in (u,v)

              msh->_finiteElement[faceGeom][solType]->GetDPhiDXi(dphidu[0], xi);

              std::vector <double> a(dim, 0.);
              for(unsigned k = 0; k < dim; k++) {
                for(unsigned i = 0; i < faceDofs; i++) {
                  a[k] += dphidu[0][i] * faceVx[k][i];
                }
              }

              double area = faceArea * WEIGHT[faceGeom][i];

              totalArea += area;

              xpT.resize(ip + 1);
              xpT[ip].resize(dim - 1);
              for(unsigned j = 0; j < dim - 1; j++) {
                xpT[ip][j].resize(dim, 0.);
              }

              if(dim == 2) {
                double norm = sqrt(a[0] * a[0] + a[1] * a[1]);
                for(unsigned k = 0; k < dim; k++) {
                  xpT[ip][0][k] = area / norm * a[k];
                }
              }
              else {
                msh->_finiteElement[faceGeom][solType]->GetDPhiDEta(dphidu[1], xi);
                std::vector <double> b(3, 0.);

                for(unsigned k = 0; k < dim; k++) {
                  for(unsigned i = 0; i < faceDofs; i++) {
                    b[k] += dphidu[1][i] * faceVx[k][i];
                  }
                }

                double norm = sqrt(pow(a[1] * b[2] - a[2] * b[1], 2) +
                                   pow(a[2] * b[0] - a[0] * b[2], 2) +
                                   pow(a[0] * b[1] - a[1] * b[0], 2));

                for(unsigned k = 0; k < dim; k++) {
                  xpT[ip][0][k] = sqrt(area / norm) * a[k];
                  xpT[ip][1][k] = sqrt(area / norm) * b[k];
                }
              }
            }
          }
        }
      }
    }
  }
  //END interface markers

  std::cout << xpI.size() << std::endl;
  std::cout << totalArea << " " << 2 * 4.99 + M_PI * 0.25 << std::endl;

  //BEGIN remove interface marker duplicates
  totalArea = 0.;
  for(unsigned i = 0; i < xpI.size(); i++) {
    unsigned cnt = 1;
    std::vector < double > N(dim);
    double area;
    if(dim == 2) {
      N[0] = xpT[i][0][1];
      N[1] = -xpT[i][0][0];
      area = sqrt(N[0] * N[0] + N[1] * N[1]);
    }
    else if(dim == 3) {
      N[0] = xpT[i][0][1] * xpT[i][1][2] - xpT[i][0][2] * xpT[i][1][1];
      N[1] = xpT[i][0][2] * xpT[i][1][0] - xpT[i][0][0] * xpT[i][1][2];
      N[2] = xpT[i][0][0] * xpT[i][1][1] - xpT[i][0][1] * xpT[i][1][0];
      area = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
    }

    for(unsigned j = i + 1; j < xpI.size(); j++) {
      bool same = true;
      for(unsigned k = 0; k < dim; k++) {
        if(fabs(xpI[i][k] - xpI[j][k]) > 1.0e-12) {
          same = false;
          break;
        }
      }
      if(same) {
        cnt++;
        std::vector < double > Nj(dim);
        if(dim == 2) {
          Nj[0] = xpT[j][0][1];
          Nj[1] = -xpT[j][0][0];
          area += sqrt(Nj[0] * Nj[0] + Nj[1] * Nj[1]);
        }
        else if(dim == 3) {
          Nj[0] = xpT[j][0][1] * xpT[j][1][2] - xpT[j][0][2] * xpT[j][1][1];
          Nj[1] = xpT[j][0][2] * xpT[j][1][0] - xpT[j][0][0] * xpT[j][1][2];
          Nj[2] = xpT[j][0][0] * xpT[j][1][1] - xpT[j][0][1] * xpT[j][1][0];
          area += sqrt(Nj[0] * Nj[0] + Nj[1] * Nj[1] + Nj[2] * Nj[2]);
        }

        for(unsigned k = 0; k < dim; k++) {
          N[k] += Nj[k];
        }
        xpI.erase(xpI.begin() + j);
        xpT.erase(xpT.begin() + j);
        j--;
      }
    }

    totalArea += area;

    if(cnt > 1) {
      double normN2 = 0;
      for(unsigned k = 0; k < dim; k++) {
        normN2 += N[k] * N[k];
      }
      for(unsigned l = 0; l < dim - 1; l++) {
        double prj = 0.;
        for(unsigned k = 0; k < dim; k++) {
          prj += xpT[i][l][k] * N[k];
        }
        prj /= normN2;
        double normT2 = 0.;
        for(unsigned k = 0; k < dim; k++) {
          xpT[i][l][k] -= prj * N[k];
          normT2 += xpT[i][l][k] * xpT[i][l][k];
        }
        if(dim == 2) {
          double scale = area / sqrt(normT2);
          for(unsigned k = 0; k < dim; k++) {
            xpT[i][l][k] *= scale;
          }
        }
      }
      if(dim == 3) {
        double area2 = 0.;
        N[0] = xpT[i][0][1] * xpT[i][1][2] - xpT[i][0][2] * xpT[i][1][1];
        N[1] = xpT[i][0][2] * xpT[i][1][0] - xpT[i][0][0] * xpT[i][1][2];
        N[2] = xpT[i][0][0] * xpT[i][1][1] - xpT[i][0][1] * xpT[i][1][0];
        area2 += sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
        for(unsigned l = 0; l < dim - 1; l++) {
          for(unsigned k = 0; k < dim; k++) {
            xpT[i][l][k] *= sqrt(area / area2);
          }
        }
      }
    }
  }
  //END remove interface marker duplicates

  std::cout << xpI.size() << std::endl;
  std::cout << totalArea << " " << 2 * 4.99 + M_PI * .25 << std::endl;

  //BEGIN get the distance of bulk particles from the interface
  for(unsigned ip = 0; ip < xp.size(); ip++) {
    double dmin = 1.0e100;
    for(unsigned jp = 0; jp < xpI.size(); jp++) {
      double d = 0.;
      for(unsigned k = 0; k < dim; k++) {
        d += (xp[ip][k] - xpI[jp][k]) * (xp[ip][k] - xpI[jp][k]);
      }
      d = sqrt(d);
      if(d < dmin) dmin = d;
    }
    dist[ip] *= dmin;
  }
  //END get the distance of bulk particles from the interface


  //BEGIN remove bulk markers that are too far from the interface
  for(unsigned ip = 0; ip < xp.size();) {
    unsigned dip = 1;
    if(dist[ip] < dminCut || dist[ip] > dmaxCut) {
      xp.erase(xp.begin() + ip);
      weight.erase(weight.begin() + ip);
      dist.erase(dist.begin() + ip);
      dip = 0;
    }
    ip += dip;
  }
  //END remove bulk markers that are too far from the interface
 
  std::ofstream fout;

  //BEGIN bulk printing
  fout.open("../input/bulk.txt");
  std::cout << dim << " " << xp.size() << std::endl;
  for(unsigned ip = 0; ip < xp.size(); ip++) {
    for(unsigned k = 0; k < dim; k++) {
      fout << xp[ip][k] << " ";
    }
    fout << weight[ip] << " ";
    fout << dist[ip] << std::endl;
  }
  fout.close();
  //END bulk printing

  //BEGIN bulk printing
  fout.open("../input/interface.txt");
  std::cout << dim << " " << xpI.size() << std::endl;
  for(unsigned ip = 0; ip < xpI.size(); ip++) {
    for(unsigned k = 0; k < dim; k++) {
      fout << xpI[ip][k] << " ";
    }
    for(unsigned l = 0; l < dim - 1; l++) {
      for(unsigned k = 0; k < dim; k++) {
        fout << xpT[ip][l][k] << " ";
      }
    }
    fout << std::endl;
  }
  fout.close();
  //END bulk printing

  std::vector < std::vector < std::vector < double > > >  bulkPoints(1);
  bulkPoints[0] = xp;
  PrintLine(DEFAULT_OUTPUTDIR, "bulk", bulkPoints, 0);

  std::vector < std::vector < std::vector < double > > >  interfacePoints(1);
  interfacePoints[0] = xpI;
  PrintLine(DEFAULT_OUTPUTDIR, "interface", interfacePoints, 0);

  std::vector < std::vector < double > > xpI2;
  xpI2 = xpI;

  totalArea = 0.;
  for(unsigned i = 0; i < xpI.size(); i++) {
    std::vector < double > N(dim);
    if(dim == 2) {
      N[0] = xpT[i][0][1];
      N[1] = -xpT[i][0][0];
      totalArea += sqrt(N[0] * N[0] + N[1] * N[1]);
    }
    else if(dim == 3) {
      N[0] = xpT[i][0][1] * xpT[i][1][2] - xpT[i][0][2] * xpT[i][1][1];
      N[1] = xpT[i][0][2] * xpT[i][1][0] - xpT[i][0][0] * xpT[i][1][2];
      N[2] = xpT[i][0][0] * xpT[i][1][1] - xpT[i][0][1] * xpT[i][1][0];
      totalArea += sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
    }

    for(unsigned k = 0; k < dim; k++) {
      xpI2[i][k] += N[k];
    }
  }

  std::cout << totalArea << " " << 2 * 4.99 + M_PI * .25 << std::endl;
  interfacePoints[0] = xpI2;
  PrintLine(DEFAULT_OUTPUTDIR, "interface2", interfacePoints, 0);

}

void FlagElements(MultiLevelMesh& mlMesh, const unsigned &layers) {

  unsigned level = mlMesh.GetNumberOfLevels() - 1;
  Mesh *msh   = mlMesh.GetLevel(level);
  unsigned iproc  = msh->processor_id();
  const unsigned dim = msh->GetDimension();

  for(unsigned iel = 0; iel < msh->el->GetElementNumber(); iel++) {
    if(msh->el->GetIfElementCanBeRefined(iel)) {
      bool refine = 0;
      unsigned ielMat = msh->GetElementMaterial(iel);

      
      for(unsigned j = 0; j < msh->el->GetElementNearElementSize(iel,1);j++){
	   int jel = msh->el->GetElementNearElement(iel,j);
      
      
      
      //for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
       // int jel = msh->el->GetFaceElementIndex(iel, iface) - 1;
        if(jel >= 0) { // iface is not a boundary of the domain
          unsigned jelMat = msh->GetElementMaterial(jel);
          if(ielMat != jelMat) { //iel and jel are on the FSI interface
            refine = true;
            break;
          }
        }
      }
      if(refine)  {
        msh->_topology->_Sol[msh->GetAmrIndex()]->set(iel, 1.);
      }
    }
  }

  msh->_topology->_Sol[msh->GetAmrIndex()]->close();

  for(unsigned k = 1; k < layers; k++) {
    for(unsigned iel = 0; iel < msh->el->GetElementNumber(); iel++) {
      if(msh->el->GetIfElementCanBeRefined(iel) && (*msh->_topology->_Sol[msh->GetAmrIndex()])(iel) == 0.) {
        bool refine = 0;
        
        for(unsigned j = 0; j < msh->el->GetElementNearElementSize(iel,1);j++){
	      int jel = msh->el->GetElementNearElement(iel,j);
        //for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          //int jel = msh->el->GetFaceElementIndex(iel, iface) - 1;
          if(jel >= 0 && (*msh->_topology->_Sol[msh->GetAmrIndex()])(jel) == 1.) { // iface is not a boundary of the domain
            refine = true;
            break;
          }
        }
        if(refine)  {
          msh->_topology->_Sol[msh->GetAmrIndex()]->set(iel, 0.5);
        }
      }
    }
    for(unsigned iel = 0; iel < msh->el->GetElementNumber(); iel++) {
      if((*msh->_topology->_Sol[msh->GetAmrIndex()])(iel) == 0.5) {
        msh->_topology->_Sol[msh->GetAmrIndex()]->set(iel, 1.);
      }
    }
    msh->_topology->_Sol[msh->GetAmrIndex()]->close();
  }
}
