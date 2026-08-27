#pragma once

class LevelMarkers {
  public:
    LevelMarkers() = default;

    std::vector<MyVector<double>>& GetFields() {
      return _field;
    }

    const std::vector<MyVector<double>>& GetFields() const {
      return _field;
    }

    std::vector<MyVector<double>>& GetLocalCoordinates() {
      return _Xi;
    }

    const std::vector<MyVector<double>>& GetLocalCoordinates() const {
      return _Xi;
    }

    MyVector<unsigned>& GetElements() {
      return _XIel;
    }

    const MyVector<unsigned>& GetElements() const {
      return _XIel;
    }

    std::vector<bool>& GetPointInsideDomain() {
      return _isInsideDomain;
    }

    const std::vector<bool>& GetPointInsideDomain() const {
      return _isInsideDomain;
    }


    unsigned& GetLevel() {
      return _level;
    }

    unsigned GetLevel() const {
      return _level;
    }

    std::vector<std::vector<unsigned>> & GetMap_s() {
      return _map_s;
    }

    const std::vector<std::vector<unsigned>> & GetMap_s() const {
      return _map_s;
    }

    std::vector<std::vector<unsigned>> & GetMap_r() {
      return _map_r;
    }

    const std::vector<std::vector<unsigned>> & GetMap_r() const {
      return _map_r;
    }



    void SetLevel(unsigned level) {
      _level = level;
    }

    void SetFieldLocalSize(const unsigned localSize) {
      _fieldLocalSize = localSize;
    }

    unsigned& GetFieldLocalSize() {
      return _fieldLocalSize;
    }

    const unsigned& GetFieldLocalSize() const {
      return _fieldLocalSize;
    }


    void RebuildLocalFromField(std::vector<std::vector<double>> &field_r, const unsigned nFields, const bool backward = true) {

      const std::vector<std::vector<unsigned>> &map = (backward) ? _map_r : _map_s;

      const unsigned nprocs = map.size();

      field_r.resize(nprocs);

      assert(_fieldLocalSize != uint_max);
      assert(nFields  == _field.size());
      if(nFields > 0) {
        for (unsigned j = 0; j < nFields; ++j) {
          assert(_field[j].begin() == _field[0].begin());
          assert(_field[j].end()   == _field[0].end());
        }
        const unsigned offset1 = _field[0].begin();
        for (unsigned kproc = 0; kproc < nprocs; ++kproc) {
          field_r[kproc].assign(map[kproc].size() * nFields, 0.);
          for (unsigned i = 0; i < map[kproc].size(); ++i) {
            if(map[kproc][i] != uint_max) {
              assert(map[kproc][i] < _fieldLocalSize);
              for (unsigned j = 0; j < nFields; ++j) {
                field_r[kproc][i * nFields + j] = _field[j][offset1 + map[kproc][i]];
              }
            }
          }
        }
      }
    }

    void SendLocalField(const std::vector<std::vector<double>>& field_s,
                        std::vector<std::vector<double>>& field_r) {
      int iproc_i, nprocs_i;
      MPI_Comm_rank(MPI_COMM_WORLD, &iproc_i);
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs_i);

      const unsigned iproc  = static_cast<unsigned>(iproc_i);
      const unsigned nprocs = static_cast<unsigned>(nprocs_i);

      if (field_s.size() != nprocs) {
        std::cerr << "Error: field_s.size() = " << field_s.size()
                  << " but MPI_Comm_size = " << nprocs << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      std::vector<unsigned> size_s(nprocs, 0);
      std::vector<unsigned> size_r(nprocs, 0);

      for (unsigned kproc = 0; kproc < nprocs; ++kproc) {
        size_s[kproc] = static_cast<unsigned>(field_s[kproc].size());
      }

      MPI_Alltoall(size_s.data(), 1, MPI_UNSIGNED,
                   size_r.data(), 1, MPI_UNSIGNED,
                   MPI_COMM_WORLD);

      field_r.resize(nprocs);
      for (unsigned kproc = 0; kproc < nprocs; ++kproc) {
        field_r[kproc].resize(size_r[kproc]);
      }

      std::vector<MPI_Request> requests;
      requests.reserve(2 * nprocs);

      // self-copy
      if (!field_s[iproc].empty()) {
        field_r[iproc] = field_s[iproc];
      }

      // receives
      for (unsigned kproc = 0; kproc < nprocs; ++kproc) {
        if (kproc == iproc) continue;
        if (!field_r[kproc].empty()) {
          MPI_Request req;
          MPI_Irecv(field_r[kproc].data(),
                    static_cast<int>(field_r[kproc].size()),
                    MPI_DOUBLE,
                    static_cast<int>(kproc),
                    100,
                    MPI_COMM_WORLD,
                    &req);
          requests.push_back(req);
        }
      }

      // sends
      for (unsigned kproc = 0; kproc < nprocs; ++kproc) {
        if (kproc == iproc) continue;
        if (!field_s[kproc].empty()) {
          MPI_Request req;
          MPI_Isend(field_s[kproc].data(),
                    static_cast<int>(field_s[kproc].size()),
                    MPI_DOUBLE,
                    static_cast<int>(kproc),
                    100,
                    MPI_COMM_WORLD,
                    &req);
          requests.push_back(req);
        }
      }

      if (!requests.empty()) {
        MPI_Waitall(static_cast<int>(requests.size()),
                    requests.data(),
                    MPI_STATUSES_IGNORE);
      }
    }

    void RebuildFieldFromLocal(const std::vector<std::vector<double>> &field_s, const unsigned nFields, const bool backward) {

      if (nFields == 0) {
        _field.clear();
        return;
      }

      const std::vector<std::vector<unsigned>> &map = (backward) ? _map_s : _map_r;

      assert(map.size() == field_s.size());
      unsigned nprocs = field_s.size();

      for (unsigned p = 0; p < nprocs; ++p) {
        assert(field_s[p].size() % nFields == 0);
        assert(map[p].size() * nFields == field_s[p].size());
      }
      assert(_fieldLocalSize != uint_max);

      std::vector<std::vector<double>> locField(nFields);
      for(unsigned k = 0; k < nFields; k++) locField[k].assign(_fieldLocalSize, 0.);

      for(unsigned p = 0; p < nprocs; p++) {
        for(unsigned i = 0; i < map[p].size(); i++ ) {
          if(map[p][i] != uint_max) {
            assert(map[p][i] < _fieldLocalSize);
            for(unsigned j = 0; j < nFields; j++ ) {
              locField[j][map[p][i]] = field_s[p][ i * nFields + j];
            }
          }
        }
      }
      _field.resize(nFields);
      for(unsigned k = 0; k < nFields; k++) _field[k].buildFromLocal(locField[k]);
    }


  private:
    std::vector<MyVector<double>> _field;
    std::vector<MyVector<double>> _Xi;
    MyVector<unsigned> _XIel;

    std::vector<std::vector<unsigned>> _map_s;
    std::vector<std::vector<unsigned>> _map_r;

    std::vector<bool> _isInsideDomain;

    unsigned _level = 0;
    static constexpr unsigned uint_max = std::numeric_limits<unsigned>::max();
    unsigned _fieldLocalSize = uint_max;

};
