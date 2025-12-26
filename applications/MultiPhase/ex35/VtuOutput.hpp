#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <array>
#include <cstdint>

#include "Mesh.hpp"
#include "Field.hpp"
#include "Encoder.hpp"   // must provide write_binary_array(os, type, name, ncomp, data)

// VTK cell types (UInt8)
static constexpr unsigned char VTK_TRIQUADRATIC_HEXAHEDRON     = 29; // Hex27
static constexpr unsigned char VTK_QUADRATIC_TETRA             = 24; // Tet10 (we export Tet15 as Tet10)
static constexpr unsigned char VTK_BIQUADRATIC_QUADRATIC_WEDGE = 32; // Wedge18 (we export Wedge21 as Wedge18)
static constexpr unsigned char VTK_BIQUADRATIC_QUAD            = 28; // Quad9
static constexpr unsigned char VTK_QUADRATIC_TRIANGLE          = 22; // Tri6 (we export Tri7 as Tri6)
static constexpr unsigned char VTK_QUADRATIC_EDGE              = 21; // Line3

inline void writeMeshFieldVTU(const std::string& filename,
                              const Field& field) {
  const Mesh& mesh = field.mesh();

  const auto& X       = mesh.X();
  const auto& elType  = mesh.elType();
  const auto& elTplgy = mesh.elTplgy();
  const auto& level   = mesh.elLevel();
  const std::size_t dim = mesh.dim();

  if (dim != 2 && dim != 3) throw std::runtime_error("writeMeshFieldVTU: dim must be 2 or 3");
  if (X.size() != dim)      throw std::runtime_error("writeMeshFieldVTU: X.size()!=dim");

  const std::size_t numPoints = X[0].size();
  const std::size_t numCells  = elTplgy.size();

  for (std::size_t d = 1; d < dim; ++d)
    if (X[d].size() != numPoints) throw std::runtime_error("writeMeshFieldVTU: X[d] size mismatch");

  if (elType.size() != numCells) throw std::runtime_error("writeMeshFieldVTU: elType size mismatch");
  if (level.size()  != numCells) throw std::runtime_error("writeMeshFieldVTU: elLevel size mismatch");

  // ---- flatten points (avoid push_back) ----
  std::vector<double> flatPoints;
  flatPoints.resize(numPoints * 3);
  for (std::size_t i = 0; i < numPoints; ++i) {
    flatPoints[3 * i + 0] = X[0][i];
    flatPoints[3 * i + 1] = (dim >= 2) ? X[1][i] : 0.0;
    flatPoints[3 * i + 2] = (dim == 3) ? X[2][i] : 0.0;
  }

  // ---- precompute total connectivity size ----
  auto cellNConn = [](unsigned et) -> std::size_t {
    switch (et) {
    case 0: return 27; // Hex27
    case 1: return 10; // Tet15 -> Tet10
    case 2: return 18; // Wedge21 -> Wedge18
    case 3: return 9;  // Quad9
    case 4: return 6;  // Tri7 -> Tri6
    case 5: return 3;  // Line3
    default: throw std::runtime_error("writeMeshFieldVTU: unknown elType");
    }
  };

  std::size_t totalConn = 0;
  for (std::size_t e = 0; e < numCells; ++e) totalConn += cellNConn(elType[e]);

  std::vector<int32_t> connectivity(totalConn);
  std::vector<int32_t> offsets(numCells);
  std::vector<uint8_t> types(numCells);

  const std::array<int, 7> hex27_tail_map{{3, 1, 0, 2, 4, 5, 6}};

  // ---- fill cells arrays by index (no push_back) ----
  std::size_t w = 0;      // write cursor in connectivity
  int32_t off = 0;

  for (std::size_t e = 0; e < numCells; ++e) {
    const unsigned et = elType[e];
    const auto& connU = elTplgy[e];

    const std::size_t n = cellNConn(et);
    if (connU.size() < n) throw std::runtime_error("writeMeshFieldVTU: conn too small");

    if (et == 0) { // Hex27 special reorder
      for (int j = 0; j < 20; ++j) connectivity[w++] = static_cast<int32_t>(connU[j]);
      for (int i = 0; i < 7;  ++i) connectivity[w++] = static_cast<int32_t>(connU[20 + hex27_tail_map[i]]);
      types[e] = VTK_TRIQUADRATIC_HEXAHEDRON;
      off += 27;
    }
    else {
      for (std::size_t j = 0; j < n; ++j) connectivity[w++] = static_cast<int32_t>(connU[j]);
      switch (et) {
      case 1: types[e] = VTK_QUADRATIC_TETRA; break;
      case 2: types[e] = VTK_BIQUADRATIC_QUADRATIC_WEDGE; break;
      case 3: types[e] = VTK_BIQUADRATIC_QUAD; break;
      case 4: types[e] = VTK_QUADRATIC_TRIANGLE; break;
      case 5: types[e] = VTK_QUADRATIC_EDGE; break;
      }
      off += static_cast<int32_t>(n);
    }

    offsets[e] = off;
  }

  // ---- per-cell level ----
  std::vector<int32_t> cellLevel(numCells);
  for (std::size_t e = 0; e < numCells; ++e) cellLevel[e] = static_cast<int32_t>(level[e]);

  // ---- output: binary mode + big buffer ----
  std::ofstream os(filename, std::ios::binary);
  if (!os) throw std::runtime_error("writeMeshFieldVTU: cannot open file");

  std::vector<char> buf(8 * 1024 * 1024);
  os.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));

  os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
     << "  <UnstructuredGrid>\n"
     << "    <Piece NumberOfPoints=\"" << numPoints << "\" NumberOfCells=\"" << numCells << "\">\n";

  os << "      <Points>\n";
  write_binary_array(os, "Float64", "", 3, flatPoints);
  os << "      </Points>\n";

  os << "      <Cells>\n";
  write_binary_array(os, "Int32", "connectivity", 1, connectivity);
  write_binary_array(os, "Int32", "offsets",      1, offsets);
  write_binary_array(os, "UInt8", "types",        1, types);
  os << "      </Cells>\n";

  os << "      <CellData Scalars=\"level\">\n";
  write_binary_array(os, "Int32", "level", 1, cellLevel);
  for (std::size_t p = 0; p < field.numFields(); ++p) {
    const unsigned fid = field.idByPos(p);
    if (field.location(fid) != Field::Location::Elemental) continue;
    const auto& v = field.getById(fid);
    if (v.size() == numCells) write_binary_array(os, "Float64", field.name(fid), 1, v);
  }
  os << "      </CellData>\n";

  os << "      <PointData>\n";
  for (std::size_t p = 0; p < field.numFields(); ++p) {
    const unsigned fid = field.idByPos(p);
    if (field.location(fid) != Field::Location::Nodal) continue;
    const auto& v = field.getById(fid);
    if (v.size() == numPoints) write_binary_array(os, "Float64", field.name(fid), 1, v);
  }
  os << "      </PointData>\n";

  os << "    </Piece>\n"
     << "  </UnstructuredGrid>\n"
     << "</VTKFile>\n";

  std::cout << "Mesh+Field written to " << filename << " (VTU UNSTRUCTURED_GRID)\n";
}


