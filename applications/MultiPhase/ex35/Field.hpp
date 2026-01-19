#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include "Mesh.hpp"

#include "FemProjection.hpp"


class Field {
  public:
    using Vec = std::vector<double>;

    enum class Location : unsigned {
      Nodal     = 0,
      Elemental = 1
    };

    Field() = delete;                 // must be bound
    explicit Field(Mesh& m) noexcept
      : _mesh(m) {}

    ~Field() = default;

    // Copy construction still not allowed (no implicit binding target)
    //Field(const Field&) = delete;

    // Deep-copy into my already-bound mesh + containers
    Field& operator=(const Field& other) {
      if (this == &other) return *this;

      _mesh = other._mesh;            // copies mesh CONTENTS (does not rebind)

      _field      = other._field;
      _fieldName  = other._fieldName;
      _fieldIndex = other._fieldIndex;
      _fieldLoc   = other._fieldLoc;
      _nameToId   = other._nameToId;
      _idToPos    = other._idToPos;
      _nextId     = other._nextId;

      return *this;
    }

    // Swap everything, keeping bindings (swap mesh contents, not attachments)
    friend void swap(Field& a, Field& b) {
      using std::swap;
      swap(a._mesh,      b._mesh);       // swaps mesh CONTENTS
      swap(a._field,     b._field);
      swap(a._fieldName, b._fieldName);
      swap(a._fieldIndex, b._fieldIndex);
      swap(a._fieldLoc,  b._fieldLoc);
      swap(a._nameToId,  b._nameToId);
      swap(a._idToPos,   b._idToPos);
      swap(a._nextId,    b._nextId);
    }

    const Mesh& mesh() const {
      // if (_mesh == nullptr) throw std::runtime_error("Field::mesh: _mesh is null");
      return _mesh;
    }

    // ------------------------------------------------------------
    // Add / Erase / Clear
    // ------------------------------------------------------------
    unsigned addField(const std::string& name, Location loc, double initVal = 0.0) {
      if (_nameToId.find(name) != _nameToId.end()) {
        throw std::runtime_error("Field::addField: field already exists: " + name);
      }

      const std::size_t n = sizeFromMesh(loc);

      const unsigned id_ = _nextId++;
      const std::size_t pos = _field.size();

      _field.emplace_back(n, initVal);
      _fieldName.emplace_back(name);
      _fieldIndex.emplace_back(id_);
      _fieldLoc.emplace_back(loc);

      _nameToId[name] = id_;
      _idToPos[id_]   = pos;

      return id_;
    }

    void eraseByName(const std::string& name) {
      eraseById(id(name));
    }

    void eraseById(unsigned id_) {
      const std::size_t pos  = idToPos(id_);
      const std::size_t last = _field.size() - 1;

      const std::string nameDel = _fieldName[pos];

      if (pos != last) {
        std::swap(_field[pos],      _field[last]);
        std::swap(_fieldName[pos],  _fieldName[last]);
        std::swap(_fieldIndex[pos], _fieldIndex[last]);
        std::swap(_fieldLoc[pos],   _fieldLoc[last]);

        _idToPos[_fieldIndex[pos]] = pos;
      }

      _field.pop_back();
      _fieldName.pop_back();
      _fieldIndex.pop_back();
      _fieldLoc.pop_back();

      _idToPos.erase(id_);
      _nameToId.erase(nameDel);
    }

    void clear() noexcept {
      _field.clear();
      _fieldName.clear();
      _fieldIndex.clear();
      _fieldLoc.clear();

      _nameToId.clear();
      _idToPos.clear();

      _nextId = 0;
    }



    void evalNodalAtLocatedPointsById(
      unsigned id_,
      const std::vector<PointLocatorResult>& results,
      const std::array<std::unique_ptr<FemProjection>, 6>& elProj,
      std::vector<double>& out,
      double outsideVal = -1.0) const {
      requireLocation(id_, Location::Nodal, "Field::evalNodalAtLocatedPointsById");

      const Mesh& m = mesh();
      const std::size_t nEl = m.numElements();
      const std::size_t nN  = m.numNodes();

      const auto& elTplgy = m.elTplgy();
      const auto& elType  = m.elType();
      const Vec&  U       = getNodalById(id_);

      if (elTplgy.size() != nEl) throw std::runtime_error("Field::eval...: elTplgy size mismatch");
      if (elType.size()  != nEl) throw std::runtime_error("Field::eval...: elType size mismatch");
      if (U.size()       != nN)  throw std::runtime_error("Field::eval...: nodal field size mismatch");

      out.assign(results.size(), outsideVal);

      std::vector<double> phi;

      for (std::size_t i = 0; i < results.size(); ++i) {
        const PointLocatorResult& r = results[i];

        // Outside / invalid -> keep outsideVal, but DO NOT skip copying logic because out is already set
        if (!r.ok || r.elem == UMAX) {
          continue;
        }

        const std::size_t e = static_cast<std::size_t>(r.elem);
        if (e >= nEl) {
          throw std::runtime_error("Field::eval...: results contains elem out of range");
        }

        const unsigned et = elType[e];
        if (et >= elProj.size()) throw std::runtime_error("Field::eval...: element type out of [0,5]");
        if (!elProj[et])         throw std::runtime_error("Field::eval...: elProj[et] is null");

        // Get basis at xi (fem() returns reference)
        phi.clear();
        elProj[et]->fem().GetPhi(phi, r.xi);


        const auto& conn = elTplgy[e];
        if (conn.empty()) {
          throw std::runtime_error("Field::eval...: empty element connectivity");
        }

        if (phi.size() != conn.size()) {
          throw std::runtime_error("Field::eval...: phi.size() != conn.size()");
        }

        double val = 0.0;
        for (std::size_t j = 0; j < conn.size(); ++j) {
          const unsigned node = conn[j];
          if (node >= nN) {
            throw std::runtime_error("Field::eval...: connectivity node out of range");
          }
          val += U[node] * phi[j];
        }

        out[i] = val;
      }
    }

    void evalNodalAtLocatedPointsById(
      unsigned id_,
      const std::vector<PointLocatorResult>& results,
      const std::array<std::unique_ptr<FemProjection>, 6>& elProj,
      std::vector<double>& out_val,
      std::vector<std::vector<double>>& out_grad,
      double outsideVal) const {
        requireLocation(id_, Location::Nodal, "Field::evalNodalAtLocatedPointsById");

        const Mesh& m = mesh();
        const std::size_t nEl = m.numElements();
        const std::size_t nN  = m.numNodes();
        const std::vector<std::vector<double> > &X = m.X(); // global coordinates

        const auto& elTplgy = m.elTplgy();
        const auto& elType  = m.elType();
        const Vec&  U       = getNodalById(id_);

        if (elTplgy.size() != nEl) throw std::runtime_error("Field::eval...: elTplgy size mismatch");
        if (elType.size()  != nEl) throw std::runtime_error("Field::eval...: elType size mismatch");
        if (U.size()       != nN)  throw std::runtime_error("Field::eval...: nodal field size mismatch");

        out_val.assign(results.size(), outsideVal);
        for(unsigned d = 0; d < m.dim(); ++d)
          out_grad[d].assign(results.size(), outsideVal);

        std::vector<double> phi;

        for (std::size_t i = 0; i < results.size(); ++i) {
          const PointLocatorResult& r = results[i];

          // Outside / invalid -> keep outsideVal, but DO NOT skip copying logic because out is already set
          if (!r.ok || r.elem == UMAX) {
            continue;
          }

          const std::size_t e = static_cast<std::size_t>(r.elem);
          if (e >= nEl) {
            throw std::runtime_error("Field::eval...: results contains elem out of range");
          }

          const unsigned et = elType[e];
          if (et >= elProj.size()) throw std::runtime_error("Field::eval...: element type out of [0,5]");
          if (!elProj[et])         throw std::runtime_error("Field::eval...: elProj[et] is null");

          const auto& conn = elTplgy[e];
          
          std::vector<std::vector<double> > xv(m.dim()); // local coordinates
          for(unsigned d = 0; d < m.dim(); d++) xv[d].resize(conn.size());
          
          for (std::size_t j = 0; j < conn.size(); ++j) {
            const unsigned node = conn[j];
            for(unsigned d = 0; d<m.dim();++d){
              xv[d][j] =  X[d][node];
            }
          }

          std::vector < double > gradphi;
          double weight;
          
          elProj[et]->fem().Jacobian(xv, r.xi, weight, phi, gradphi);

          if (conn.empty()) {
            throw std::runtime_error("Field::eval...: empty element connectivity");
          }

          if (phi.size() != conn.size()) {
            throw std::runtime_error("Field::eval...: phi.size() != conn.size()");
          }

          if (gradphi.size() != conn.size() * m.dim()) {
            throw std::runtime_error("Field::eval...: gradphi.size() != conn.size() * mesh.dim()");
          }

          double val = 0.0;
          std::vector<double> grad(m.dim());

          for (std::size_t j = 0; j < conn.size(); ++j) {
            const unsigned node = conn[j];
            if (node >= nN) {
              throw std::runtime_error("Field::eval...: connectivity node out of range");
            }
            for(unsigned d = 0; d < m.dim(); d++)
              grad[d] += U[node] * gradphi[m.dim() * j + d];
            val += U[node] * phi[j];
          }

          for(unsigned d = 0; d < m.dim(); d++)
            out_grad[d][i] = grad[d];
          out_val[i] = val;
        }
      }





      // ------------------------------------------------------------
      // Queries / Iteration support
      // ------------------------------------------------------------
      std::size_t numFields() const noexcept {
        return _field.size();
      }
      bool empty() const noexcept {
        return _field.empty();
      }

      bool hasName(const std::string & name) const {
        return _nameToId.find(name) != _nameToId.end();
      }
      bool hasId(unsigned id_) const {
        return _idToPos.find(id_) != _idToPos.end();
      }

      unsigned idByPos(std::size_t pos) const {
        if (pos >= _fieldIndex.size()) throw std::runtime_error("Field::idByPos: pos out of range");
        return _fieldIndex[pos];
      }

      const std::string& nameByPos(std::size_t pos) const {
        if (pos >= _fieldName.size()) throw std::runtime_error("Field::nameByPos: pos out of range");
        return _fieldName[pos];
      }

      Location locationByPos(std::size_t pos) const {
        if (pos >= _fieldLoc.size()) throw std::runtime_error("Field::locationByPos: pos out of range");
        return _fieldLoc[pos];
      }

      // ------------------------------------------------------------
      // Name <-> id
      // ------------------------------------------------------------
      unsigned id(const std::string & name) const {
        auto it = _nameToId.find(name);
        if (it == _nameToId.end()) throw std::runtime_error("Field::id: name not found: " + name);
        return it->second;
      }

      const std::string& name(unsigned id_) const {
        return _fieldName[idToPos(id_)];
      }

      Location location(unsigned id_) const {
        return _fieldLoc[idToPos(id_)];
      }

      // ------------------------------------------------------------
      // Get field vectors
      // ------------------------------------------------------------
      Vec& getById(unsigned id_) {
        return _field[idToPos(id_)];
      }
      const Vec& getById(unsigned id_) const {
        return _field[idToPos(id_)];
      }

      Vec& getByName(const std::string & name_) {
        return getById(id(name_));
      }
      const Vec& getByName(const std::string & name_) const {
        return getById(id(name_));
      }

      Vec& getNodalById(unsigned id_) {
        requireLocation(id_, Location::Nodal, "Field::getNodalById");
        return getById(id_);
      }
      const Vec& getNodalById(unsigned id_) const {
        requireLocation(id_, Location::Nodal, "Field::getNodalById");
        return getById(id_);
      }

      Vec& getElementalById(unsigned id_) {
        requireLocation(id_, Location::Elemental, "Field::getElementalById");
        return getById(id_);
      }
      const Vec& getElementalById(unsigned id_) const {
        requireLocation(id_, Location::Elemental, "Field::getElementalById");
        return getById(id_);
      }

      // ------------------------------------------------------------
      // Scalar access
      // ------------------------------------------------------------
      double valueById(unsigned id_, std::size_t k) const {
        const Vec& v = getById(id_);
        if (k >= v.size()) throw std::runtime_error("Field::valueById: entry out of range");
        return v[k];
      }

      void setValueById(unsigned id_, std::size_t k, double val) {
        Vec& v = getById(id_);
        if (k >= v.size()) throw std::runtime_error("Field::setValueById: entry out of range");
        v[k] = val;
      }

      double valueByName(const std::string & name_, std::size_t k) const {
        return valueById(id(name_), k);
      }

      void setValueByName(const std::string & name_, std::size_t k, double val) {
        setValueById(id(name_), k, val);
      }

      // ------------------------------------------------------------
      // Replace whole vector (size must match mesh size for that location)
      // ------------------------------------------------------------
      void setFieldById(unsigned id_, const Vec & data) {
        const Location loc = location(id_);
        const std::size_t target = sizeFromMesh(loc);
        if (data.size() != target) throw std::runtime_error("Field::setFieldById: size mismatch vs mesh");
        getById(id_) = data;
      }

      void setFieldByName(const std::string & name_, const Vec & data) {
        setFieldById(id(name_), data);
      }

      // ------------------------------------------------------------
      // Mesh sync (after refine/coarsen)
      // ------------------------------------------------------------
      void resizeToMesh(double fillVal = 0.0) {
        for (std::size_t pos = 0; pos < _field.size(); ++pos) {
          const std::size_t target = sizeFromMesh(_fieldLoc[pos]);
          _field[pos].assign(target, fillVal);
        }
      }

      // void rebindMeshAndResize(Mesh& newMesh, double fillVal = 0.0) {
      //   _mesh = &newMesh;
      //   resizeToMesh(fillVal);
      // }



      // ------------------------------------------------------------
      // DOF coordinate (by id / by name)
      // - Nodal: returns mesh node coordinate at index k
      // - Elemental: returns coordinate of the "center node" = last connectivity entry of element k
      // ------------------------------------------------------------
      std::vector<double> dofCoordById(unsigned id_, std::size_t k) const {
        //if (_mesh == nullptr) throw std::runtime_error("Field::dofCoordById: _mesh is null");

        const Location loc = location(id_);
        const std::size_t d = _mesh.dim();
        if (d == 0) throw std::runtime_error("Field::dofCoordById: mesh.dim()==0");

        std::vector<double> x(d, 0.0);

        if (loc == Location::Nodal) {
          if (k >= _mesh.numNodes()) throw std::runtime_error("Field::dofCoordById: nodal k out of range");
          const auto& X = _mesh.X();
          if (X.size() != d) throw std::runtime_error("Field::dofCoordById: mesh.X().size() != dim()");
          for (std::size_t a = 0; a < d; ++a) {
            if (X[a].size() <= k) throw std::runtime_error("Field::dofCoordById: mesh.X()[a] size mismatch");
            x[a] = X[a][k];
          }
          return x;
        }

        if (loc == Location::Elemental) {
          if (k >= _mesh.numElements()) throw std::runtime_error("Field::dofCoordById: elemental k out of range");
          const auto& elTplgy = _mesh.elTplgy();
          if (elTplgy.size() <= k) throw std::runtime_error("Field::dofCoordById: mesh.elTplgy() size mismatch");
          const auto& conn = elTplgy[k];
          if (conn.empty()) throw std::runtime_error("Field::dofCoordById: empty element connectivity");
          const unsigned centerNode = conn.back();
          if (centerNode >= _mesh.numNodes()) throw std::runtime_error("Field::dofCoordById: center node out of range");

          const auto& X = _mesh.X();
          if (X.size() != d) throw std::runtime_error("Field::dofCoordById: mesh.X().size() != dim()");
          for (std::size_t a = 0; a < d; ++a) {
            if (X[a].size() <= centerNode) throw std::runtime_error("Field::dofCoordById: mesh.X()[a] size mismatch");
            x[a] = X[a][centerNode];
          }
          return x;
        }

        throw std::runtime_error("Field::dofCoordById: unknown Location");
      }

      std::vector<double> dofCoordByName(const std::string & name_, std::size_t k) const {
        return dofCoordById(id(name_), k);
      }


      std::size_t extractInterfaceVerticesAndCentersById(
        unsigned id_,
        std::vector<std::vector<double>>& Xuniq,
        unsigned min_level = 0u,
        unsigned max_level = std::numeric_limits<unsigned>::max()) const {
        requireLocation(id_, Location::Nodal, "Field::extractInterfaceVerticesAndCentersById");

        const Mesh& m = mesh();
        const std::size_t d   = m.dim();
        const std::size_t nEl = m.numElements();
        const std::size_t nN  = m.numNodes();

        if (d == 0) throw std::runtime_error("Field::extractInterfaceVerticesAndCentersById: mesh.dim()==0");

        // Prepare output (reuse capacity if caller keeps Xuniq around)
        Xuniq.clear();
        Xuniq.resize(d);
        for (std::size_t a = 0; a < d; ++a) {
          Xuniq[a].clear();
          Xuniq[a].reserve(nN); // conservative; reusing capacity if already larger
        }

        if (min_level >= max_level) return 0;

        const auto& elTplgy = m.elTplgy();
        const auto& elType  = m.elType();
        const auto& elLevel = m.elLevel();
        const auto& X       = m.X();
        const Vec&  U       = getNodalById(id_);

        auto numVerticesFromType = [&](unsigned t) -> unsigned {
          switch (t) {
            case 0:
              return 8; // Hex27
            case 1:
              return 4; // Tet15
            case 2:
              return 6; // Wedge21
            case 3:
              return 4; // Quad9
            case 4:
              return 3; // Tri7
            case 5:
              return 2; // Line3
            default:
              throw std::runtime_error("Field::extractInterfaceVerticesAndCentersById: unknown element type");
          }
        };

        std::vector<char> markedNode(nN, 0);
        std::size_t nAdded = 0;

        auto pushNodeIfNew = [&](unsigned node) {
          if (node >= nN) throw std::runtime_error("Field::extractInterfaceVerticesAndCentersById: node out of range");
          if (!markedNode[node]) {
            markedNode[node] = 1;
            ++nAdded;
            for (std::size_t a = 0; a < d; ++a) {
              Xuniq[a].push_back(X[a][node]);
            }
          }
        };

        for (std::size_t e = 0; e < nEl; ++e) {
          const unsigned L = elLevel[e];
          if (L < min_level || L >= max_level) continue;

          const auto& conn = elTplgy[e];
          if (conn.empty()) continue;

          const unsigned Nv = numVerticesFromType(elType[e]);
          if (conn.size() < static_cast<std::size_t>(Nv + 1u)) {
            throw std::runtime_error("Field::extractInterfaceVerticesAndCentersById: connectivity too small for element type");
          }

          // sign change on vertices only
          double umin = 0.0, umax = 0.0;
          bool first = true;
          for (unsigned i = 0; i < Nv; ++i) {
            const unsigned node = conn[i];
            if (node >= nN) throw std::runtime_error("Field::extractInterfaceVerticesAndCentersById: vertex node out of range");
            const double u = U[node];
            if (first) {
              umin = umax = u;
              first = false;
            }
            else {
              if (u < umin) umin = u;
              if (u > umax) umax = u;
            }
          }

          if (!(umin <= 0.0 && umax >= 0.0)) continue;

          for (unsigned i = 0; i < Nv; ++i) pushNodeIfNew(conn[i]); // vertices only
          pushNodeIfNew(conn.back());                               // center node (last)
        }

        return nAdded;
      }

      std::size_t extractInterfaceVerticesAndCentersByName(
        const std::string & name_,
        std::vector<std::vector<double>>& Xuniq,
        unsigned min_level = 0u,
        unsigned max_level = std::numeric_limits<unsigned>::max()) const {
        return extractInterfaceVerticesAndCentersById(id(name_), Xuniq, min_level, max_level);
      }


    private:
      //  Mesh* _mesh = nullptr;
      Mesh& _mesh; // bound, non-owning

      std::vector<Vec>         _field;
      std::vector<std::string> _fieldName;
      std::vector<unsigned>    _fieldIndex; // stable ids, parallel to _field
      std::vector<Location>    _fieldLoc;

      std::unordered_map<std::string, unsigned> _nameToId; // name -> stable id
      std::unordered_map<unsigned, std::size_t> _idToPos;  // id -> packed position

      unsigned _nextId = 0;

      std::size_t sizeFromMesh(Location loc) const {
        //if (_mesh == nullptr) throw std::runtime_error("Field::sizeFromMesh: _mesh is null");
        if (loc == Location::Nodal)     return _mesh.numNodes();
        if (loc == Location::Elemental) return _mesh.numElements();
        throw std::runtime_error("Field::sizeFromMesh: unknown Location");
      }

      std::size_t idToPos(unsigned id_) const {
        auto it = _idToPos.find(id_);
        if (it == _idToPos.end()) throw std::runtime_error("Field::idToPos: id not found");
        return it->second;
      }

      void requireLocation(unsigned id_, Location want, const char* where) const {
        const Location have = _fieldLoc[idToPos(id_)];
        if (have != want) throw std::runtime_error(std::string(where) + ": wrong field location");
      }
    };

