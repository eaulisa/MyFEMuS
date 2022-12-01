void Cloud::ComputeQuadraticBestFit() {
  _A.clear();



  map<unsigned, bool> pSearch;

  Mesh *msh = _sol->GetMesh();

  unsigned SolQIndex = _sol->GetIndex("Q");
  _sol->_Sol[SolQIndex]->zero();

  unsigned dim = _sol->GetMesh()->GetDimension();
  std::vector<std::vector<double>> normOld;
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
    normOld.resize(i1 - i0, std::vector<double> (dim));
    norm.assign(dim, 0);
    unsigned cnt = 0;
    double weightSum = 0;

    for(unsigned i = i0; i < i1; i++, cnt++) {
      for(unsigned k = 0; k < dim; k++) {
        coord[cnt][k] = _yp[_map[i]][k];
        norm[k] += _N[_map[i]][k];
        normOld[cnt][k] = _N[_map[i]][k];
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
        weight.resize(coord.size() + (j1 - j0));
        normOld.resize(normOld.size() + (j1 - j0), std::vector<double> (dim));
        for(unsigned j = j0; j < j1; j++) {
          double dotProduct = 0.;
          for(unsigned k = 0; k < dim; k++) {
            coord[cnt][k] = _yp[_map[j]][k];
            normOld[cnt][k] =  _N[_map[j]][k];
          }
          weight[cnt] = _ds[_map[j]];
          cnt++;
        }
        coord.resize(cnt);
        weight.resize(cnt);
      }
    }
    if(parallelSearch == false) {
      double sigma2 = 0.;
      std::vector<double> d2(cnt, 0.);
      if(cnt > 1) {
        for(unsigned i = 0; i < cnt; i++) {
          for(unsigned k = 0; k < dim; k++) {
            d2[i] += (coord[i][k] - xn[k]) * (coord[i][k] - xn[k]);
          }
          sigma2 += d2[i];
        }
        sigma2 /= cnt * 7.;
        for(unsigned i = 0; i < cnt; i++) {
          weight[i] *= exp(-d2[i] / sigma2);
        }
      }
      else {
        std::cerr << "Abbiamo solo un marker!!!!!!!!!!!!!!!!!!!!!!!!!\n" ;
        abort();
      }

      std::vector<double*> d2P(d2.size());
      for(unsigned i = 0; i < d2P.size(); i++) {
        d2P[i] = &d2[i];
      }
      std::sort(d2P.begin() + cnt0, d2P.end(), [](const double * a, const double * b) {
        return *a < *b;
      });
      std::vector<unsigned> mapj(d2.size());
      for(unsigned i = 0; i < d2.size(); i++) {
        mapj[i] =  static_cast<unsigned>(d2P[i] - &d2[0]);
      }

      std::vector<std::vector<double>> coord1(coord.size(), std::vector<double>(dim));
      std::vector<std::vector<double>> normOld1(normOld.size(), std::vector<double>(dim));
      std::vector<double> weight1(weight.size());

      for(unsigned i = 0; i < coord.size(); i++) {
        coord1[i] = coord[mapj[i]];
        normOld1[i] = normOld[mapj[i]];
        weight1[i] = weight[mapj[i]];
      }
      coord1.swap(coord);
      normOld1.swap(normOld);
      weight1.swap(weight);

      if(cnt0 <= 4) {
        cnt0 = (coord.size() >= 9) ? cnt0 + 5 : coord.size();
      }

      unsigned newSize = (weight.size() - cnt0 > 30) ? (cnt0 + 30) : weight.size();
      weight.resize(newSize);
      coord.resize(newSize);
      normOld.resize(newSize);

      //This return an Ellipse or an Hyperbola
      FindQuadraticBestFit(coord, weight, norm, _A[iel]);

      std::vector<double> dotProduct(cnt0, 0.);

      double n1Dotn = 0;
      for(unsigned i = 0; i < cnt0; i++) {
        std::vector <double> n1 = GetNormal(iel, coord[i]);
        for(unsigned k = 0; k < dim; k++) {
          dotProduct[i] += normOld[i][k] * n1[k];
        }
        n1Dotn += dotProduct[i] * weight[i]; //TODO is this ok if we swap the weights above?
      }

      if(n1Dotn < 0) {
        for(unsigned  i = 0; i < _A[iel].size(); i++) {
          _A[iel][i] *= -1.;
        }
        for(unsigned i = 0; i < dotProduct.size(); i++) {
          dotProduct[i] *= -1.;
        }
      }
      double cost1 = GetCost(coord, dotProduct, weight, iel, cnt0);
      std::vector<double> Acon = _A[iel];

      //This return a parabola
      double minDP = 1;
      for(unsigned i = 0; i < cnt0 - 1; i++) {
        for(unsigned j = i + 1; j < cnt0; j++) {
          double dp = 0.;
          for(unsigned k = 0; k < dim; k++) {
            dp += normOld[j][k] * normOld[i][k];
          }
          if(dp < minDP) minDP = dp;
        }
      }

      femus::GetQuadricBestFit(coord, weight, norm, _A[iel], cnt0, minDP); //parabola
      dotProduct.assign(cnt0, 0);
      n1Dotn = 0;
      for(unsigned i = 0; i < cnt0; i++) {
        std::vector <double> n1 = GetNormal(iel, coord[i]);
        for(unsigned k = 0; k < dim; k++) {
          dotProduct[i] += normOld[i][k] * n1[k];
        }
        n1Dotn += dotProduct[i] * weight[i]; //TODO is this ok if we swap the weights above?
      }

      if(n1Dotn < 0) {
        for(unsigned  i = 0; i < _A[iel].size(); i++) {
          _A[iel][i] *= -1.;
        }
        for(unsigned i = 0; i < dotProduct.size(); i++) {
          dotProduct[i] *= -1.;
        }
      }
      double cost2 = GetCost(coord, dotProduct, weight, iel, cnt0);
      std::vector<double> Apar = _A[iel];

//       if(true && iel == 165) {
//         std::cout << "a = " << Acon[0] << "; b=" << Acon[1] << "; c=" << Acon[2] << "; d=" << Acon[3] << "; e=" << Acon[4] << "; f= " << Acon[5] << ";\n";
//         std::cout << "a = " << Apar[0] << "; b=" << Apar[1] << "; c=" << Apar[2] << "; d=" << Apar[3] << "; e=" << Apar[4] << "; f= " << Apar[5] << ";\n";
//       }

      _A[iel] = (cost1 < cost2) ? Acon : Apar;
      double useOldPoints = 1.;

      if(iel == 419) std::cout << cost1 << " " << cost2 << std::endl;


      if(cost2 < 1.e-10) {
        _A[iel] = Apar;
      }
      else if(cost1 > 1.0e-6 && cost2 > 1.0e-6) {
        useOldPoints = -1.;
      }
      else if(cost1 / cost2 > 1.0e-5 && cost2 / cost1 > 1.0e-5 && (cost1 < 1.0e-4 && cost2 < 1.0e-4))   {
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
            for(unsigned k = 0; k < dim; k++) {
              counter += (n1[k] - n2[k]) * normOld[jMin][k];
            }
          }
        }
        if(counter < 0) {
          _A[iel] = Apar;
        }
        else if(counter > 0) {
          _A[iel] = Acon;
        }
      }

      double delta = _A[iel][1] * _A[iel][1] - 4. * _A[iel][0] * _A[iel][2];
      if(fabs(delta) < 1.0e-5) _sol->_Sol[SolQIndex]->set(iel, useOldPoints * 1); // parabola;
      else if(delta < 0) _sol->_Sol[SolQIndex]->set(iel, useOldPoints * 2); // ellipse;
      else _sol->_Sol[SolQIndex]->set(iel, useOldPoints * 3); // hyperpola;


      if(true && iel == 165) {
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
            normOld.resize(i1 - i0, std::vector<double> (dim));
            cnt0 = i1 - i0;
            norm.assign(dim, 0);
            unsigned cnt = 0;
            double weightSum = 0;
            for(unsigned i = i0; i < i1; i++, cnt++) {
              for(unsigned k = 0; k < dim; k++) {
                coord[cnt][k] = _yp[_map[i]][k];
                norm[k] += _N[_map[i]][k];
                normOld[cnt][k] = _N[_map[i]][k];
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

          for(unsigned i = 1; i < nNgbElms; i++) {

            int jel;
            if(iproc == kp) {
              jel = msh->el->GetElementNearElement(kel, i);
            }
            MPI_Bcast(&jel, 1, MPI_INT, kp, PETSC_COMM_WORLD);

            unsigned jp = msh->IsdomBisectionSearch(jel, 3);  // return  jproc for piece-wise constant discontinuous type (3)

            std::vector<std::vector<double>> coordJel;
            std::vector<std::vector<double>> normOldJel;
            std::vector<double> weightJel;
            unsigned cntJel = 0;
            if(iproc == jp) {
              if(_elMrkIdx.find(jel) != _elMrkIdx.end()) {   // if jel is cut cell
                unsigned j0 = _elMrkIdx[jel][0];
                unsigned j1 = _elMrkIdx[jel][1];
                coordJel.resize(dim, std::vector<double> (j1 - j0));
                normOldJel.resize(dim, std::vector<double> (j1 - j0));
                weightJel.resize(j1 - j0);
                for(unsigned j = j0; j < j1; j++, cntJel++) {
                  for(unsigned k = 0; k < dim; k++) {
                    coordJel[k][cntJel] = _yp[_map[j]][k];
                    normOldJel[k][cntJel] =  _N[_map[j]][k];
                  }
                  weightJel[cntJel] = _ds[_map[j]];
                }
              }

              if(jp != kp) {


                MPI_Send(&cntJel, 1, MPI_UNSIGNED, kp, 0, MPI_COMM_WORLD);
                if(cntJel != 0) {
                  for(unsigned k = 0; k < dim; k++) {
                    MPI_Send(coordJel[k].data(), coordJel[k].size(), MPI_DOUBLE, kp, 1 + k, MPI_COMM_WORLD);
                    MPI_Send(normOldJel[k].data(), normOldJel[k].size(), MPI_DOUBLE, kp, dim + 1 + k, MPI_COMM_WORLD);
                  }
                  MPI_Send(weightJel.data(), weightJel.size(), MPI_DOUBLE, kp, 2 * dim + 1, MPI_COMM_WORLD);
                }
              }
            }

            if(iproc == kp) {
              if(kp != jp) {

                MPI_Recv(&cntJel, 1, MPI_UNSIGNED, jp, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(cntJel != 0) {
                  coordJel.resize(dim, std::vector<double> (cntJel));
                  normOldJel.resize(dim, std::vector<double> (cntJel));
                  weightJel.resize(cntJel);
                  for(unsigned k = 0; k < dim; k++) {
                    MPI_Recv(coordJel[k].data(), coordJel[k].size(), MPI_DOUBLE, jp, 1 + k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(normOldJel[k].data(), normOldJel[k].size(), MPI_DOUBLE, jp, dim + 1 + k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                  }
                  MPI_Recv(weightJel.data(), weightJel.size(), MPI_DOUBLE, jp, 2 * dim + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
              }
              if(cntJel != 0) {

                unsigned size0 = coord.size();
                coord.resize(size0 + cntJel, std::vector<double> (dim));
                normOld.resize(size0 + cntJel, std::vector<double> (dim));
                weight.resize(size0 + cntJel);
                for(unsigned j = 0; j < cntJel; j++) {
                  for(unsigned k = 0; k < dim; k++) {
                    coord[size0 + j][k] = coordJel[k][j];
                    normOld[size0 + j][k] = normOldJel[k][j];
                  }
                  weight[size0 + j] = weightJel[j];
                }
              }
            }
          }//face loop

          if(iproc == kp) {

            unsigned cnt = coord.size();

            double sigma2 = 0.;
            std::vector<double> d2(cnt, 0.);
            if(cnt > 1) {
              for(unsigned i = 0; i < cnt; i++) {
                for(unsigned k = 0; k < dim; k++) {
                  d2[i] += (coord[i][k] - xn[k]) * (coord[i][k] - xn[k]);
                }
                sigma2 += d2[i];
              }
              sigma2 /= cnt * 2.;
              for(unsigned i = 0; i < cnt; i++) {
                weight[i] *= exp(-0.5 / sigma2 * d2[i]);
              }
            }
            else {
              std::cerr << "Abbiamo solo un marker!!!!!!!!!!!!!!!!!!!!!!!!!\n" ;
              abort();
            }

            std::vector<double*> d2P(d2.size());
            for(unsigned i = 0; i < d2P.size(); i++) {
              d2P[i] = &d2[i];
            }
            std::sort(d2P.begin() + cnt0, d2P.end(), [](const double * a, const double * b) {
              return *a < *b;
            });
            std::vector<unsigned> mapj(d2.size());
            for(unsigned i = 0; i < d2.size(); i++) {
              mapj[i] =  static_cast<unsigned>(d2P[i] - &d2[0]);
            }

            std::vector<std::vector<double>> coord1(coord.size(), std::vector<double>(dim));
            std::vector<std::vector<double>> normOld1(normOld.size(), std::vector<double>(dim));
            std::vector<double> weight1(weight.size());

            for(unsigned i = 0; i < coord.size(); i++) {
              coord1[i] = coord[mapj[i]];
              normOld1[i] = normOld[mapj[i]];
              weight1[i] = weight[mapj[i]];
            }
            coord1.swap(coord);
            normOld1.swap(normOld);
            weight1.swap(weight);

            if(cnt0 <= 4) {
              cnt0 = (coord.size() >= 9) ? cnt0 + 5 : coord.size();
            }

            unsigned newSize = (weight.size() - cnt0 > 30) ? (cnt0 + 30) : weight.size();
            weight.resize(newSize);
            coord.resize(newSize);
            normOld.resize(newSize);

            //This return an Ellipse or an Hyperbola
            FindQuadraticBestFit(coord, weight, norm, _A[kel]);

            std::vector<double> dotProduct(cnt0, 0.);

            double n1Dotn = 0;
            for(unsigned i = 0; i < cnt0; i++) {
              std::vector <double> n1 = GetNormal(kel, coord[i]);
              for(unsigned k = 0; k < dim; k++) {
                dotProduct[i] += normOld[i][k] * n1[k];
              }
              n1Dotn += dotProduct[i] * weight[i]; //TODO is this ok if we swap the weights above?
            }

            if(n1Dotn < 0) {
              for(unsigned  i = 0; i < _A[kel].size(); i++) {
                _A[kel][i] *= -1.;
              }
              for(unsigned i = 0; i < dotProduct.size(); i++) {
                dotProduct[i] *= -1.;
              }
            }
            double cost1 = GetCost(coord, dotProduct, weight, kel, cnt0);
            std::vector<double> Acon = _A[kel];

            //This return a parabola
            double minDP = 1;
            for(unsigned i = 0; i < cnt0 - 1; i++) {
              for(unsigned j = i + 1; j < cnt0; j++) {
                double dp = 0.;
                for(unsigned k = 0; k < dim; k++) {
                  dp += normOld[j][k] * normOld[i][k];
                }
                if(dp < minDP) minDP = dp;
              }
            }

            femus::GetQuadricBestFit(coord, weight, norm, _A[kel], cnt0, minDP); //parabola
            dotProduct.assign(cnt0, 0);
            n1Dotn = 0;
            for(unsigned i = 0; i < cnt0; i++) {
              std::vector <double> n1 = GetNormal(kel, coord[i]);
              for(unsigned k = 0; k < dim; k++) {
                dotProduct[i] += normOld[i][k] * n1[k];
              }
              n1Dotn += dotProduct[i] * weight[i]; //TODO is this ok if we swap the weights above?
            }

            if(n1Dotn < 0) {
              for(unsigned  i = 0; i < _A[kel].size(); i++) {
                _A[kel][i] *= -1.;
              }
              for(unsigned i = 0; i < dotProduct.size(); i++) {
                dotProduct[i] *= -1.;
              }
            }
            double cost2 = GetCost(coord, dotProduct, weight, kel, cnt0);
            std::vector<double> Apar = _A[kel];

            if(true && kel == 165) {
              std::cout << "a = " << Acon[0] << "; b=" << Acon[1] << "; c=" << Acon[2] << "; d=" << Acon[3] << "; e=" << Acon[4] << "; f= " << Acon[5] << ";\n";
              std::cout << "a = " << Apar[0] << "; b=" << Apar[1] << "; c=" << Apar[2] << "; d=" << Apar[3] << "; e=" << Apar[4] << "; f= " << Apar[5] << ";\n";
            }

            _A[kel] = (cost1 < cost2) ? Acon : Apar;
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
                  for(unsigned k = 0; k < dim; k++) {
                    counter += (n1[k] - n2[k]) * normOld[jMin][k];
                  }
                }
              }
              if(counter < 0) {
                _A[kel] = Apar;
              }
              else if(counter > 0) {
                _A[kel] = Acon;
              }
            }

            double delta = _A[kel][1] * _A[kel][1] - 4. * _A[kel][0] * _A[kel][2];
            if(fabs(delta) < 1.0e-5) _sol->_Sol[SolQIndex]->set(kel, 1); // parabola;
            else if(delta < 0) _sol->_Sol[SolQIndex]->set(kel, 2); // ellipse;
            else _sol->_Sol[SolQIndex]->set(kel, 3); // hyperpola;


            if(true && kel == 165) {
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

              std::cout << "BBBBBBB\n" ;//<< treshold <<  " ";
              std::cout << cost1 << " " << cost2 << std::endl;
              std::cout << "a = " << _A[kel][0] << "; b=" << _A[kel][1] << "; c=" << _A[kel][2] << "; d=" << _A[kel][3] << "; e=" << _A[kel][4] << "; f= " << _A[kel][5] << ";\n";
              std::cout << " {x," << xMin[0] << ", " << xMax[0] << "},{" << "{y," << xMin[1] << ", " << xMax[1] << "}" << std::endl;
            }





            it++;
          }
        }//element loop
      }
    }


  }
  _sol->_Sol[SolQIndex]->close();

}





