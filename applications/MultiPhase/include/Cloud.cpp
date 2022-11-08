void Cloud::ComputeQuadraticBestFit() {
  _A.clear();

  map<unsigned, bool> pSearch;

  Mesh *msh = _sol->GetMesh();

  unsigned dim = _sol->GetMesh()->GetDimension();
  std::vector<std::vector<double>> coord;
  std::vector<double> weight;
  coord.reserve(_elem.size());
  std::vector<double> norm;
  std::vector<double> xn;

  unsigned iproc = _sol->processor_id();
  unsigned nprocs = _sol->n_processors();

  for(_itElMrkIdx = _elMrkIdx.begin(); _itElMrkIdx != _elMrkIdx.end(); _itElMrkIdx++) {
    unsigned iel = _itElMrkIdx->first;
    unsigned i0 = _itElMrkIdx->second[0];
    unsigned i1 = _itElMrkIdx->second[1];

    unsigned cnt0 = i1 - i0;
    coord.resize(i1 - i0, std::vector<double> (dim));
    weight.resize(i1 - i0);
    norm.assign(dim, 0);
    unsigned cnt = 0;
    double weightSum = 0;

    for(unsigned i = i0; i < i1; i++, cnt++) {
      for(unsigned k = 0; k < dim; k++) {
        coord[cnt][k] = _yp[_map[i]][k];
        norm[k] += _N[_map[i]][k];
      }
      weight[cnt] = _ds[_map[i]];
      weightSum += weight[cnt];
    }

    double det = 0;
    for(unsigned k = 0; k < dim; k++) {
      det += norm[k] * norm[k];
    }
    det = sqrt(det);
    for(unsigned k = 0; k < dim; k++) {
      norm[k] /= det;
    }

    xn.assign(dim, 0);

    for(unsigned i = 0; i < cnt; i++) {
      for(unsigned k = 0; k < dim; k++) {
        xn[k] += weight[i] * coord[i][k];
      }
    }
    for(unsigned k = 0; k < dim; k++) {
      xn[k] /= weightSum;
    }

    unsigned nFaces = msh->GetElementFaceNumber(iel);
    std::vector<unsigned> jelFace(nFaces);
    unsigned cntF = 0;
    for(unsigned iface = 0; iface < nFaces; iface++) {
      int jel = msh->el->GetFaceElementIndex(iel, iface) - 1;
      if(jel >= 0) {
        jelFace[cntF] = jel;
        cntF++;
      }
    }
    jelFace.resize(cntF);

    bool parallelSearch = false;

    for(unsigned i = 1; i < msh->el->GetElementNearElementSize(iel, 1); i++) {
      int jel = msh->el->GetElementNearElement(iel, i);
      unsigned jproc = msh->IsdomBisectionSearch(jel, 3);
      if(jproc != iproc) {
        pSearch[iel] = true;
        parallelSearch = true;
        break;
      }

      bool isFaceElement = false;

      for(unsigned iFace = 0; iFace < jelFace.size(); iFace++) {
        if(jel == jelFace[iFace]) {
          isFaceElement = true;
          break;
        }
      }

      if(_elMrkIdx.find(jel) != _elMrkIdx.end()) { //jel is a cut fem inside iproc
        unsigned j0 = _elMrkIdx[jel][0];
        unsigned j1 = _elMrkIdx[jel][1];
        coord.resize(coord.size() + (j1 - j0), std::vector<double> (dim));
        weight.resize(coord.size() + (j1 - j0), 0.250 * !isFaceElement + isFaceElement);
        for(unsigned j = j0; j < j1; j++) {
          double dotProduct = 0.;
          for(unsigned k = 0; k < dim; k++) {
            coord[cnt][k] = _yp[_map[j]][k];
          }
          weight[cnt] = _ds[_map[j]] * (0.250 * !isFaceElement + isFaceElement);
          cnt++;
        }
        coord.resize(cnt);
        weight.resize(cnt);
      }
    }
    if(parallelSearch == false) {
      double sigma2 = 0.;
      double sigma = 0.;
      if(cnt > 1) {
        for(unsigned i = 0; i < cnt; i++) {
          for(unsigned k = 0; k < dim; k++) {
            sigma2 += (coord[i][k] - xn[k]) * (coord[i][k] - xn[k]);
          }
        }
        sigma2 /= cnt;
        sigma2 /= 2;
        sigma = sqrt(sigma2);
        for(unsigned i = 0; i < cnt; i++) {
          double a = 0;
          for(unsigned k = 0; k < dim; k++) {
            a += -0.5 / sigma2 * (coord[i][k] - xn[k]) * (coord[i][k] - xn[k]);
          }
          weight[i] *= exp(a) / (sigma * sqrt(2. * M_PI));
        }
      }
      else {
        std::cerr << "Abbiamo solo un marker!!!!!!!!!!!!!!!!!!!!!!!!!\n" ;
        abort();
      }

      //This return an Ellipse or an Hperbola
      FindQuadraticBestFit(coord, weight, norm, _A[iel]);
      std::vector<double> dotProduct(i1 - i0, 0.);
      double n1Dotn = 0;
      for(unsigned i = i0; i < i1; i++) {
        std::vector <double> n1 = GetNormal(iel,  _yp[_map[i]]);
        for(unsigned k = 0; k < dim; k++) {
          dotProduct[i - i0] += _N[_map[i]][k] * n1[k];
        }
        n1Dotn += dotProduct[i - i0] * weight[i - i0];
      }
      if(n1Dotn < 0) {
        for(unsigned  i = 0; i < _A[iel].size(); i++) {
          _A[iel][i] *= -1.;
        }
        for(unsigned i = 0; i < cnt0; i++) {
          dotProduct[i] *= -1.;
        }
      }
      double cost1 = GetCost(coord, dotProduct, weight, iel, cnt0);

//       unsigned cnt1 = coord.size() - cnt0;
//       std::vector<unsigned> j(cnt1);
//       std::vector<double*> weightP(cnt1);
//       for(unsigned i = 0; i < cnt1; i++) {
//         weightP[i] = &weight[cnt0 + i];
//       }
//       std::sort(weightP.begin(), weightP.end(), [](const double * a, const double * b) {
//         return *a < *b;
//       });
//       for(unsigned i = 0; i < cnt1; i++) {
//         j[i] =  cnt0 + static_cast<unsigned>(weightP[i] - &weight[cnt0]);
//       }
//       for(unsigned k = cnt0 + cnt1; k < coord.size(); k++) {
//         bool keepWhile = true;
//         unsigned kn = cnt1 - 1;
//         while(keepWhile) {
//           keepWhile = false;
//           if(weight[k] > weight[j[kn]]) {
//             std::swap(k, j[kn]);
//             if(kn > 0) {
//               keepWhile = true;
//               kn--;
//             }
//           }
//         }
//       }

//       std::vector<std::vector<double>> coordNew = coord;
//       std::vector<double> weightNew = weight;
//       for(unsigned k = 0; k < cnt1; k++) {
//         std::swap(weight[j[k]], weightNew[cnt0 + k]);
//         coord[j[k]].swap(coordNew[cnt0 + k]);
//       }

      //This return a parabola





      std::vector<double> Acon = _A[iel];
      if(cnt0 <= 2) {
        unsigned nEP = 3;
        nEP = (nEP > coord.size() - cnt0) ? (coord.size() - cnt0) : nEP;
        if(coord.size() > 6) {
          for(unsigned j = cnt0; j < cnt0 + nEP; j++) {
            unsigned k = max_element(weight.begin() + j, weight.end()) - weight.begin();
            swap(weight[j], weight[k]);
            coord[j].swap(coord[k]);
          }
          if(iel == 298 || iel == 208) {
            std::cout << xn[0] << " " << xn[1] << std::endl;
            for(unsigned j = 0; j < coord.size(); j++) {
              std::cout << cnt0 << " " << iel << " " << j << " " << weight[j]  << " " << coord[j][0]  << " " << coord[j][1] << std::endl;
            }
          }
          //cnt0+=nEP;
          femus::GetQuadricBestFit(coord, weight, norm, _A[iel], cnt0 + nEP); //parabola
        }
        else {
          femus::GetQuadricBestFit(coord, weight, norm, _A[iel], coord.size()); //parabola

        }
      }
      else {
        femus::GetQuadricBestFit(coord, weight, norm, _A[iel], cnt0); //parabola
      }

      dotProduct.assign(i1 - i0, 0);
      n1Dotn = 0;
      for(unsigned i = i0; i < i1; i++) {
        std::vector <double> n1 = GetNormal(iel, coord[i - i0]);
        for(unsigned k = 0; k < dim; k++) {
          dotProduct[i - i0] += _N[_map[i]][k] * n1[k];
        }
        n1Dotn += dotProduct[i - i0] * weight[i - i0]; //TODO is this ok if we swap the weights above?
      }

      if(n1Dotn < 0) {
        for(unsigned  i = 0; i < _A[iel].size(); i++) {
          _A[iel][i] *= -1.;
        }
        for(unsigned i = 0; i < cnt0; i++) {
          dotProduct[i] *= -1.;
        }
      }
      double cost2 = GetCost(coord, dotProduct, weight, iel, cnt0);
      std::vector<double> Apar = _A[iel];

      if(false && iel == 4784) {
        std::cout << "a = " << Acon[0] << "; b=" << Acon[1] << "; c=" << Acon[2] << "; d=" << Acon[3] << "; e=" << Acon[4] << "; f= " << Acon[5] << ";\n";
        std::cout << "a = " << Apar[0] << "; b=" << Apar[1] << "; c=" << Apar[2] << "; d=" << Apar[3] << "; e=" << Apar[4] << "; f= " << Apar[5] << ";\n";
      }

      if(cnt0 <= 2)  cnt0 += 3;

      if((cost2 < 0.0001 && cost1 / cost2 > 0.00001) && cost2 / cost1 > 0.00001) {

        double counter = 0.;
        std::vector <double> xp(dim);
        unsigned nDofs =  msh->GetElementDofNumber(iel, 2);
        for(unsigned i = 0; i < nDofs; i++) {
          unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
          for(unsigned k = 0; k < dim; k++) {
            xp[k] = (*msh->_topology->_Sol[k])(xDof);
          }
          std::vector<double> n1 = GetNormal(Acon, xp);
          std::vector<double> n2 = GetNormal(Apar, xp);
          double dotProd = 0;
          for(unsigned k = 0; k < dim; k++) {
            dotProd += n1[k] * n2[k];
          }
          if(dotProd < 0 && GetValue(Apar, xp) * GetValue(Acon, xp) < 0) {
            double distMin = 1.0E10;
            unsigned jMin = 0;
            for(unsigned j = 1; j < cnt0; j++) {
              double distj = 0.;
              for(unsigned k = 0; k < dim; k++) {
                distj += (coord[j][k] - xp[k]) * (coord[j][k] - xp[k]);
              }
              if(distj > distMin) {
                jMin = j;
                distMin = distj;
              }
            }

            double dotProd1 = 0;
            for(unsigned k = 0; k < dim; k++) {
              dotProd1 += n1[k] * _N[_map[i0 + jMin]][k];
            }

            double dotProd2 = 0;
            for(unsigned k = 0; k < dim; k++) {
              dotProd2 += n2[k] * _N[_map[i0 + jMin]][k];
            }
            counter += dotProd1 - dotProd2;
          }
        }
        if(counter < 0) {
          Acon = Apar;
        }
        else if(counter > 0) {
          Apar = Acon;
        }
      }

      _A[iel] = (cost1 < cost2) ? Acon : Apar;


      if(false && iel == 4784) {
        double xMin[2] = { 1.0e10, 1.0e10};
        double xMax[2] = { -1.0e10, -1.0e10};

        unsigned nDofs =  msh->GetElementDofNumber(iel, 2);
        for(unsigned i = 0; i < nDofs; i++) {
          unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
          for(unsigned k = 0; k < dim; k++) {
            double x = (*msh->_topology->_Sol[k])(xDof);
            if(x < xMin[k]) xMin[k] = x;
            if(x > xMax[k]) xMax[k] = x;
          }
        }

        std::cout << "BBBBBBB\n" ;//<< treshold <<  " ";
        std::cout << cost1 << " " << cost2 << std::endl;
        std::cout << "a = " << _A[iel][0] << "; b=" << _A[iel][1] << "; c=" << _A[iel][2] << "; d=" << _A[iel][3] << "; e=" << _A[iel][4] << "; f= " << _A[iel][5] << ";\n";
        std::cout << " {x," << xMin[0] << ", " << xMax[0] << "},{" << "{y," << xMin[1] << ", " << xMax[1] << "}" << std::endl;
      }
    }
  }

  map<unsigned, bool>::iterator it;


  if(nprocs > 1) {

    for(unsigned kp = 0; kp < nprocs; kp++) {

      unsigned nel;
      if(iproc == kp) {
        nel = pSearch.size();
      }
      MPI_Bcast(&nel, 1, MPI_UNSIGNED, kp, MPI_COMM_WORLD);

      if(nel > 0) {
        if(iproc == kp) {
          it =  pSearch.begin();
        }
        for(unsigned cntEl = 0; cntEl < nel; cntEl++) {
          unsigned kel;
          unsigned nNgbElms;
          unsigned i0;
          unsigned i1;

          std::vector<unsigned> jelFace;
          unsigned cnt0 = 0;

          if(iproc == kp) {
            kel = it->first;
            i0 = _elMrkIdx[kel][0];
            i1 = _elMrkIdx[kel][1];
            coord.resize(i1 - i0, std::vector<double> (dim));
            weight.resize(i1 - i0);
            cnt0 = i1 - i0;
            norm.assign(dim, 0);
            unsigned cnt = 0;
            double weightSum = 0;
            for(unsigned i = i0; i < i1; i++, cnt++) {
              for(unsigned k = 0; k < dim; k++) {
                coord[cnt][k] = _yp[_map[i]][k];
                norm[k] += _N[_map[i]][k];
              }
              weight[cnt] = _ds[_map[i]];
              weightSum += weight[cnt];
            }

            xn.assign(dim, 0);

            for(unsigned i = 0; i < cnt; i++) {
              for(unsigned k = 0; k < dim; k++) {
                xn[k] += weight[i] * coord[i][k];
              }
            }
            for(unsigned k = 0; k < dim; k++) {
              xn[k] /= weightSum;
            }

            unsigned nFaces = msh->GetElementFaceNumber(kel);
            jelFace.resize(nFaces);
            unsigned cntF = 0;
            for(unsigned iface = 0; iface < nFaces; iface++) {
              int jel = msh->el->GetFaceElementIndex(kel, iface) - 1;
              if(jel >= 0) {
                jelFace[cntF] = jel;
                cntF++;
              }
            }
            jelFace.resize(cntF);
            nNgbElms = msh->el->GetElementNearElementSize(kel, 1);
          }
          MPI_Bcast(&nNgbElms, 1, MPI_UNSIGNED, kp, PETSC_COMM_WORLD);
          jelFace.resize(nNgbElms);
          MPI_Bcast(jelFace.data(), jelFace.size(), MPI_UNSIGNED, kp, PETSC_COMM_WORLD);
//             MPI_Bcast(&i1, 1, MPI_UNSIGNED, kp, PETSC_COMM_WORLD);
//             MPI_Bcast(&i0, 1, MPI_UNSIGNED, kp, PETSC_COMM_WORLD);

          for(unsigned i = 1; i < nNgbElms; i++) {

            int jel;
            if(iproc == kp) {
              jel = msh->el->GetElementNearElement(kel, i);
            }
            MPI_Bcast(&jel, 1, MPI_INT, kp, PETSC_COMM_WORLD);

            unsigned jp = msh->IsdomBisectionSearch(jel, 3);  // return  jproc for piece-wise constant discontinuous type (3)

            bool isFaceElement = false;
            for(unsigned iFace = 0; iFace < jelFace.size(); iFace++) {
              if(jel == jelFace[iFace]) {
                isFaceElement = true;
                break;
              }
            }

            std::vector<std::vector<double>> coordJel;
            std::vector<double> weightJel;
            unsigned cntJel = 0;
            if(iproc == jp) {
              if(_elMrkIdx.find(jel) != _elMrkIdx.end()) {   // if jel is cut cell
                unsigned j0 = _elMrkIdx[jel][0];
                unsigned j1 = _elMrkIdx[jel][1];
                coordJel.resize(dim, std::vector<double> (j1 - j0));
                weightJel.resize(j1 - j0);
                for(unsigned j = j0; j < j1; j++, cntJel++) {
                  for(unsigned k = 0; k < dim; k++) {
                    coordJel[k][cntJel] = _yp[_map[j]][k];
                  }
                  weightJel[cntJel] = _ds[_map[j]] * (0.250 * !isFaceElement + isFaceElement);
                }
              }

              if(jp != kp) {



                MPI_Send(&cntJel, 1, MPI_UNSIGNED, kp, 0, MPI_COMM_WORLD);
                if(cntJel != 0) {
                  for(unsigned k = 0; k < dim; k++) {
                    MPI_Send(coordJel[k].data(), coordJel[k].size(), MPI_DOUBLE, kp, 1 + k, MPI_COMM_WORLD);
                  }
                  MPI_Send(weightJel.data(), weightJel.size(), MPI_DOUBLE, kp, dim + 1, MPI_COMM_WORLD);
                }
              }
            }

            if(iproc == kp) {
              if(kp != jp) {

                MPI_Recv(&cntJel, 1, MPI_UNSIGNED, jp, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(cntJel != 0) {
                  coordJel.resize(dim, std::vector<double> (cntJel));
                  weightJel.resize(cntJel);
                  for(unsigned k = 0; k < dim; k++) {
                    MPI_Recv(coordJel[k].data(), coordJel[k].size(), MPI_DOUBLE, jp, 1 + k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                  }
                  MPI_Recv(weightJel.data(), weightJel.size(), MPI_DOUBLE, jp, dim + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
              }
              if(cntJel != 0) {

                unsigned size0 = coord.size();
                coord.resize(size0 + cntJel, std::vector<double> (dim));
                weight.resize(size0 + cntJel);
                for(unsigned j = 0; j < cntJel; j++) {
                  for(unsigned k = 0; k < dim; k++) {
                    coord[size0 + j][k] = coordJel[k][j];
                  }
                  weight[size0 + j] = weightJel[j];
                }
              }
            }
          }//face loop

          if(iproc == kp) {

            unsigned cnt = coord.size();
            double sigma2 = 0.;
            double sigma = 0.;
            if(cnt > 1) {
              for(unsigned i = 0; i < cnt; i++) {
                for(unsigned k = 0; k < dim; k++) {
                  sigma2 += (coord[i][k] - xn[k]) * (coord[i][k] - xn[k]);
                }
              }
              sigma2 /= cnt;
              sigma2 /= 2;
              sigma = sqrt(sigma2);
              for(unsigned i = 0; i < cnt; i++) {
                double a = 0;
                for(unsigned k = 0; k < dim; k++) {
                  a += -0.5 / sigma2 * (coord[i][k] - xn[k]) * (coord[i][k] - xn[k]);
                }
                weight[i] *= exp(a) / (sigma * sqrt(2. * M_PI));
              }
            }
            else {
              std::cerr << "Abbiamo solo un marker!!!!!!!!!!!!!!!!!!!!!!!!!\n" ;
              abort();
            }

            FindQuadraticBestFit(coord, weight, norm, _A[kel]); //conica

            std::vector<double> dotProduct(i1 - i0, 0.);
            double n1Dotn = 0;
            for(unsigned i = i0; i < i1; i++) {
              std::vector <double> n1 = GetNormal(kel,  _yp[_map[i]]);
              for(unsigned k = 0; k < dim; k++) {
                dotProduct[i - i0] += _N[_map[i]][k] * n1[k];
              }
              n1Dotn += dotProduct[i - i0] * weight[i - i0];
            }
            if(n1Dotn < 0) {
              for(unsigned  i = 0; i < _A[kel].size(); i++) {
                _A[kel][i] *= -1.;
              }
              for(unsigned i = 0; i < cnt0; i++) {
                dotProduct[i] *= -1.;
              }
            }



            double cost1 = GetCost(coord, dotProduct, weight, kel, cnt0);
            std::vector<double> Acon = _A[kel];
            if(cnt0 <= 2) {
              unsigned nEP = 3;
              nEP = (nEP > coord.size() - cnt0) ? (coord.size() - cnt0) : nEP;
              if(coord.size() > 6) {

                std::vector<unsigned> j(nEP);

                std::vector<double*> weightP(nEP);
                for(unsigned i = 0; i < nEP; i++) {
                  weightP[i] = &weight[cnt0 + i];
                }
                std::sort(weightP.begin(), weightP.end(), [](const double * a, const double * b) {
                  return *a < *b;
                });
                for(unsigned i = 0; i < nEP; i++) {
                  j[i] =  cnt0 + static_cast<unsigned>(weightP[i] - &weight[cnt0]);
                }

                for(unsigned k = cnt0 + nEP; k < coord.size(); k++) {
                  bool keepWhile = true;
                  unsigned kn = nEP - 1;
                  while(keepWhile) {
                    keepWhile = false;
                    if(weight[k] > weight[j[kn]]) {
                      std::swap(k, j[kn]);
                      if(kn > 0) {
                        keepWhile = true;
                        kn--;
                      }
                    }
                  }
                }
                for(unsigned k = 0; k < nEP; k++) {
                  std::swap(weight[j[k]], weight[cnt0 + k]);
                  coord[j[k]].swap(coord[cnt0 + k]);
                }

                femus::GetQuadricBestFit(coord, weight, norm, _A[kel], cnt0 + nEP); //parabola

              }
              else {
                femus::GetQuadricBestFit(coord, weight, norm, _A[kel], coord.size()); //parabola
              }
            }
            else {
              femus::GetQuadricBestFit(coord, weight, norm, _A[kel], cnt0); //parabola
            }

            dotProduct.assign(i1 - i0, 0);
            n1Dotn = 0;
            for(unsigned i = i0; i < i1; i++) {
              std::vector <double> n1 = GetNormal(kel,  _yp[_map[i]]);
              for(unsigned k = 0; k < dim; k++) {
                dotProduct[i - i0] += _N[_map[i]][k] * n1[k];
              }
              n1Dotn += dotProduct[i - i0] * weight[i - i0]; //TODO is this ok if we swap the weights above?
            }

            if(n1Dotn < 0) {
              for(unsigned  i = 0; i < _A[kel].size(); i++) {
                _A[kel][i] *= -1.;
              }
              for(unsigned i = 0; i < cnt0; i++) {
                dotProduct[i] *= -1.;
              }
            }

            double cost2 = GetCost(coord, dotProduct, weight, kel, cnt0);
            std::vector<double> Apar = _A[kel];

            if(false && kel == 4016) {
              std::cerr << "a = " << Acon[0] << "; b=" << Acon[1] << "; c=" << Acon[2] << "; d=" << Acon[3] << "; e=" << Acon[4] << "; f= " << Acon[5] << ";\n";
              std::cerr << "a = " << Apar[0] << "; b=" << Apar[1] << "; c=" << Apar[2] << "; d=" << Apar[3] << "; e=" << Apar[4] << "; f= " << Apar[5] << ";\n";
            }


            if((cost2 < 0.0001 && cost1 / cost2 > 0.00001) && cost2 / cost1 > 0.00001) {

              double counter = 0.;
              std::vector <double> xp(dim);
              unsigned nDofs =  msh->GetElementDofNumber(kel, 2);
              for(unsigned i = 0; i < nDofs; i++) {
                unsigned xDof  = msh->GetSolutionDof(i, kel, 2);
                for(unsigned k = 0; k < dim; k++) {
                  xp[k] = (*msh->_topology->_Sol[k])(xDof);
                }
                std::vector<double> n1 = GetNormal(Acon, xp);
                std::vector<double> n2 = GetNormal(Apar, xp);
                double dotProd = 0;
                for(unsigned k = 0; k < dim; k++) {
                  dotProd += n1[k] * n2[k];
                }

                if(dotProd < 0 && GetValue(Apar, xp) * GetValue(Acon, xp) < 0) {
                  double distMin = 1.0E10;
                  unsigned jMin = 0;
                  for(unsigned j = 1; j < cnt0; j++) {
                    double distj = 0.;
                    for(unsigned k = 0; k < dim; k++) {
                      distj += (coord[j][k] - xp[k]) * (coord[j][k] - xp[k]);
                    }
                    if(distj > distMin) {
                      jMin = j;
                      distMin = distj;
                    }
                  }

                  double dotProd1 = 0;
                  for(unsigned k = 0; k < dim; k++) {
                    dotProd1 += n1[k] * _N[_map[i0 + jMin]][k];
                  }
                  double dotProd2 = 0;
                  for(unsigned k = 0; k < dim; k++) {
                    dotProd2 += n2[k] * _N[_map[i0 + jMin]][k];
                  }
                  counter += dotProd1 - dotProd2;
                }
              }
              if(counter < 0) {
                Acon = Apar;
              }
              else if(counter > 0) {
                Apar = Acon;
              }
            }

            _A[kel] = (cost1 < cost2) ? Acon : Apar;

            if(false && kel == 4016) {
              double xMin[2] = { 1.0e10, 1.0e10};
              double xMax[2] = { -1.0e10, -1.0e10};

              unsigned nDofs =  msh->GetElementDofNumber(kel, 2);
              for(unsigned i = 0; i < nDofs; i++) {
                unsigned xDof  = msh->GetSolutionDof(i, kel, 2);
                for(unsigned k = 0; k < dim; k++) {
                  double x = (*msh->_topology->_Sol[k])(xDof);
                  if(x < xMin[k]) xMin[k] = x;
                  if(x > xMax[k]) xMax[k] = x;
                }
              }
              std::cerr << "BBBBBBB\n" ;//<< treshold <<  " ";
              std::cerr << cost1 << " " << cost2 << std::endl;
              std::cerr << "a = " << _A[kel][0] << "; b=" << _A[kel][1] << "; c=" << _A[kel][2] << "; d=" << _A[kel][3] << "; e=" << _A[kel][4] << "; f= " << _A[kel][5] << ";\n";
              std::cerr << " {x," << xMin[0] << ", " << xMax[0] << "},{" << "{y," << xMin[1] << ", " << xMax[1] << "}" << std::endl;
            }
            it++;
          }
        }//element loop
      }
    }


  }


}





