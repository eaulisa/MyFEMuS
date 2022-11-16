

#ifndef __femus_AdaptiveSplit_hpp__
#define __femus_AdaptiveSplit_hpp__

namespace femus {

  class AdaptiveSplit {
    public:
      void Split(const std::vector<std::vector<double>> &xv, const std::vector<std::vector<double>> &Ni, const unsigned ielType, const unsigned &level, const unsigned &father);

      void SetFather(const std::vector<std::vector<double>> &xp, const std::vector<std::vector<double>> &xi, const std::vector<std::vector<double>> &Np, const std::vector<std::vector<double>> Ni) {
        _xp0 = xp;
      };

      std::vector<std::vector<double>> &GetXpFather() {
        return _xp0;
      }
      std::vector<std::vector<double>> &GetNpFather() {
        return _Np0;
      }
      std::vector<std::vector<double>> &GetNiFather() {
        return _Ni0;
      }
      std::vector<std::vector<double>> &GetXiFather() {
        if(_xi.size() == 0) _xi.resize(1);
        if(_xi[0].size() == 0) _xi[0].resize(1);
        return _xi[0][0];
      }


    private:

      std::vector<std::vector<double>> _xp0;
      std::vector<std::vector<double>> _Np0;
      std::vector<std::vector<double>> _Ni0;
      std::vector<std::vector<std::vector<std::vector<double>>>> _xi;
      std::vector<std::vector<std::vector<unsigned>>> _i;
      std::vector<std::vector<std::vector<double>>> _signNi; //?

  };

  void AdaptiveSplit::Split(const std::vector<std::vector<double>> &xv, const std::vector<std::vector<double>> &Ni, const unsigned ielType, const unsigned &level, const unsigned &father) {

    if(level == 0) {
      if(_i.size() < 1) _i.resize(1);
      _i[0].assign(1, std::vector<unsigned> (_xp0.size()));
      for(unsigned j = 0; j < _xp0.size(); j++) _i[0][0][j] = j;
    }

    const unsigned &np = _i[level][father].size();
    if(np > 1) {
     
      const unsigned &dim = xv.size();

      //std::cout << level << " " << child << std::endl;

      unsigned nChilds = (dim == 2) ? 4 : 8;
      std::vector<double> j(nChilds, 0);
      //std::vector<std::vector<std::vector<double>>> xpj(nChilds, std::vector<std::vector<double>>(xp.size()));
      std::vector<std::vector<std::vector<double>>> Nij(nChilds, std::vector<std::vector<double>>(np));

      if(_xi.size() < level + 2) _xi.resize(level + 2);
      _xi[level + 1].assign(nChilds, std::vector<std::vector<double>>(np, std::vector<double>(dim)));
      std::vector<std::vector<double>> &xi = _xi[level][father];


      if(_i.size() < level + 2) _i.resize(level + 2);
      _i[level + 1].assign(nChilds, std::vector<unsigned> (_i[level][father].size()));

      std::vector<unsigned> &iFather = _i[level][father];


      //std::vector<std::vector<std::vector<double>>> xij(nChilds, std::vector<std::vector<double>>(xp.size(), std::vector<double>(dim)));
      for(unsigned i = 0; i < np; i++) {
        if(ielType == 3) {
          if(xi[i][0] < 0) {
            if(xi[i][1] < 0) {
              //xpj[0][j[0]] = xp[i];
              Nij[0][j[0]] = Ni[i];
              _xi[level + 1][0][j[0]][0] = -1 + 2 * (xi[i][0] + 1);
              _xi[level + 1][0][j[0]][1] = -1 + 2 * (xi[i][1] + 1);
              _i[level + 1][0][j[0]] = iFather[i];

              j[0]++;
            }
            else {
              //xpj[3][j[3]] = xp[i];
              Nij[3][j[3]] = Ni[i];
              _xi[level + 1][3][j[3]][0] = -1 + 2 * (xi[i][0] + 1);
              _xi[level + 1][3][j[3]][1] = -1 + 2 * xi[i][1];
              _i[level + 1][3][j[3]] = iFather[i];

              j[3]++;
            }
          }
          else {
            if(xi[i][1] < 0) {
             // xpj[1][j[1]] = xp[i];
              Nij[1][j[1]] = Ni[i];
              _xi[level + 1][1][j[1]][0] = -1 + 2 * xi[i][0];
              _xi[level + 1][1][j[1]][1] = -1 + 2 * (xi[i][1] + 1);
              _i[level + 1][1][j[1]] = iFather[i];

              j[1]++;
            }
            else {
           //   xpj[2][j[2]] = xp[i];
              Nij[2][j[2]] = Ni[i];
              _xi[level + 1][2][j[2]][0] = -1 + 2 * xi[i][0];
              _xi[level + 1][2][j[2]][1] = -1 + 2 * xi[i][1];
              _i[level + 1][2][j[2]] = iFather[i];

              j[2]++;
            }

          }
        }

        else if(ielType == 4) {
          if(xi[i][0] > 0.5) {
          //  xpj[1][j[1]] = xp[i];
            Nij[1][j[1]] = Ni[i];
            _xi[level + 1][1][j[1]][0] = 2 * (xi[i][0] - 0.5);
            _xi[level + 1][1][j[1]][1] = 2 * xi[i][1];
            _i[level + 1][1][j[1]] = iFather[i];
            j[1]++;
          }
          else if(xi[i][1] > 0.5) {
           // xpj[2][j[2]] = xp[i];
            Nij[2][j[2]] = Ni[i];
            _xi[level + 1][2][j[2]][0] = 2 * xi[i][0];
            _xi[level + 1][2][j[2]][1] = 2 * (xi[i][1] - 0.5);
            _i[level + 1][2][j[2]] = iFather[i];
            j[2]++;
          }
          else if(xi[i][0] + xi[i][1] < 0.5) {
           // xpj[0][j[0]] = xp[i];
            Nij[0][j[0]] = Ni[i];
            _xi[level + 1][0][j[0]][0] = 2 * xi[i][0];
            _xi[level + 1][0][j[0]][1] = 2 * xi[i][1];
            _i[level + 1][0][j[0]] = iFather[i];
            j[0]++;
          }
          else {
          //  xpj[3][j[3]] = xp[i];

            Nij[3][j[3]].resize(2);
            Nij[3][j[3]][0] = -Ni[i][0];
            Nij[3][j[3]][1] = -Ni[i][1];

            _xi[level + 1][3][j[3]][0] = 1 - 2 * xi[i][0];
            _xi[level + 1][3][j[3]][1] = 1 - 2 * xi[i][1];
            _i[level + 1][3][j[3]] = iFather[i];
            j[3]++;
          }

        }
      }

      for(unsigned l = 0; l < nChilds; l++) {

        //xpj[l].resize(j[l]);
        Nij[l].resize(j[l]);
        _xi[level + 1][l].resize(j[l]);
        _i[level + 1][l].resize(j[l]);
        unsigned nve = (ielType == 3) ? 4 : 3;
        std::vector<std::vector<double> > xvj(dim, std::vector<double>(nve));


        xvj.assign(dim, std::vector<double>(nve, 0.));
        for(unsigned k = 0; k < dim; k++) {
          for(unsigned I = 0; I < nve; I++) {
            for(unsigned J = 0 ; J < nve; J++) {
              xvj[k][I] += PJ[ielType][l][I][J] * xv[k][J];
            }
          }
        }

        this->Split(xvj, Nij[l], ielType, level + 1, l);
      }
    }
    else if(np == 1) {
      unsigned nve = (ielType == 3) ? 4 : 3;

      double xmin = 100;
      double xmax =-100;
      double ymin = 100;
      double ymax =-100;
      
      
      for(unsigned i = 0; i < nve; i++) {
       // std::cout << xv[0][i] << " " << xv[1][i] << std::endl;
        if(xmin > xv[0][i]) xmin = xv[0][i];
        if(xmax < xv[0][i]) xmax = xv[0][i];
        if(ymin > xv[1][i]) ymin = xv[1][i];
        if(ymax < xv[1][i]) ymax = xv[1][i];
      }
    //  std::cout << xv[0][0] << " " << xv[1][0] << std::endl;
    //  std::cout << std::endl;
      
      
      unsigned i0 =_i[level][father][0];
      double x1 = _xp0[i0][0];
      double y1 = _xp0[i0][1];
      
      if(x1 < xmin || x1 > xmax ) std::cout<<"AAAAAAAAAAAAAA "<<i0 <<" "<<father <<" "<< xmin <<  " " << x1 << " " << xmax <<std::endl;
      if(y1 < ymin || y1 > ymax ) std::cout<<"BBBBBBBBBBBBBB "<<i0 <<" "<<father <<" "<< ymin <<  " " << y1 << " " << ymax <<std::endl;
      
      //std::cout << level << " " << child << " " << xp[0][0] << " " << xp[0][1] << " " << xi[0][0] << " " << xi[0][1] << std::endl;
    }
    else {
      unsigned nve = (ielType == 3) ? 4 : 3;

      for(unsigned i = 0; i < nve; i++) {
   //     std::cout << xv[0][i] << " " << xv[1][i] << std::endl;
      }
    //  std::cout << xv[0][0] << " " << xv[1][0] << std::endl;
     // std::cout << std::endl;
      //std::cout << level << " " << child << " empty or full\n";
    }
  }

}
#endif
