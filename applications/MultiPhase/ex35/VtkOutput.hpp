#pragma once
#include <fstream>
#include <iostream>
#include <array>
#include <vector>
#include <string>


// elType convention (adjust if you use different codes):
//   0 -> Hex27
//   1 -> Tet15
//   2 -> Wedge21
//   3 -> Quad9
//   4 -> Tri7

//   VTK cell types:
//   VTK_TRIQUADRATIC_HEXAHEDRON     = 29
//   VTK_QUADRATIC_TETRA             = 24   (10-node tetra)
//   VTK_BIQUADRATIC_QUADRATIC_WEDGE = 32   (18-node wedge)
//   VTK_BIQUADRATIC_QUAD            = 28
//   VTK_QUADRATIC_TRIANGLE          = 22   (6-node triangle)

// Field is linked to a Mesh internally (Field stores Mesh& _mesh)

void writeMeshFieldVTK(const std::string& filename,
                       const Field& field) {
  const Mesh& mesh = field.mesh();

  const std::vector<std::vector<double>>& X         = mesh.X();
  const std::vector<unsigned>&            elType    = mesh.elType();
  const std::vector<std::vector<unsigned>>& elTplgy = mesh.elTplgy();
  const std::vector<unsigned>&            level     = mesh.elLevel();
  const std::size_t                       dim       = mesh.dim();

  if (dim != 2 && dim != 3) {
    std::cerr << "writeMeshFieldVTK: mesh.dim() must be 2 or 3, got " << dim << "\n";
    return;
  }
  if (X.size() != dim) {
    std::cerr << "writeMeshFieldVTK: X.size() != dim (" << X.size() << " != " << dim << ")\n";
    return;
  }
  if (X[0].empty()) {
    std::cerr << "writeMeshFieldVTK: X[0] is empty\n";
    return;
  }

  const std::size_t numPoints = X[0].size();
  const std::size_t numCells  = elTplgy.size();

  for (std::size_t d = 1; d < dim; ++d) {
    if (X[d].size() != numPoints) {
      std::cerr << "writeMeshFieldVTK: X[" << d << "] size mismatch ("
                << X[d].size() << " != " << numPoints << ")\n";
      return;
    }
  }

  if (elType.size() != numCells) {
    std::cerr << "writeMeshFieldVTK: elType.size() != elTplgy.size() ("
              << elType.size() << " != " << numCells << ")\n";
    return;
  }

  std::ofstream out(filename.c_str());
  if (!out.is_open()) {
    std::cerr << "writeMeshFieldVTK: cannot open file: " << filename << "\n";
    return;
  }

  // ------------------------------------------------------
  // VTK legacy header
  // ------------------------------------------------------
  out << "# vtk DataFile Version 3.0\n";
  out << "FEM mesh + fields\n";
  out << "ASCII\n";
  out << "DATASET UNSTRUCTURED_GRID\n";

  // ------------------------------------------------------
  // POINTS
  // ------------------------------------------------------
  out << "POINTS " << numPoints << " double\n";
  if (dim == 2) {
    for (std::size_t i = 0; i < numPoints; ++i) {
      out << X[0][i] << " " << X[1][i] << " 0.0\n";
    }
  }
  else {   // dim == 3
    for (std::size_t i = 0; i < numPoints; ++i) {
      out << X[0][i] << " " << X[1][i] << " " << X[2][i] << "\n";
    }
  }

  // ------------------------------------------------------
  // CELLS: compute total ints
  // ------------------------------------------------------
  std::size_t totalCellInts = 0;
  for (std::size_t e = 0; e < numCells; ++e) {
    switch (elType[e]) {
    case 0: totalCellInts += 1 + 27; break; // Hex27
    case 1: totalCellInts += 1 + 10; break; // Tet15 -> export as 10-node quadratic tetra
    case 2: totalCellInts += 1 + 18; break; // Wedge21 -> export as 18-node wedge
    case 3: totalCellInts += 1 + 9;  break; // Quad9
    case 4: totalCellInts += 1 + 6;  break; // Tri7 -> export as 6-node quadratic tri
    default:
      std::cerr << "writeMeshFieldVTK: unknown elType[" << e << "] = " << elType[e] << "\n";
      out.close();
      return;
    }
  }

  out << "CELLS " << numCells << " " << totalCellInts << "\n";

  for (std::size_t e = 0; e < numCells; ++e) {
    const auto& conn = elTplgy[e];

    switch (elType[e]) {
    case 0: {
      // HEX27: last 7 DOFs remapped with tail_map
      static const std::array<unsigned, 7> tail_map{{3u, 1u, 0u, 2u, 4u, 5u, 6u}};
      if (conn.size() != 27) {
        std::cerr << "writeMeshFieldVTK: Hex element with conn.size() != 27 (" << conn.size() << ")\n";
      }

      out << 27;
      for (unsigned j = 0; j < 27 && j < conn.size(); ++j) {
        unsigned src = j;
        if (j >= 20) {
          const unsigned i = j - 20; // i in [0,6]
          src = 20u + tail_map[i];
        }
        out << " " << conn[src];
      }
      out << "\n";
      break;
    }

    case 1: {
      // TET15 exported as 10-node VTK_QUADRATIC_TETRA (0..9)
      if (conn.size() < 10) {
        std::cerr << "writeMeshFieldVTK: Tet element with conn.size() < 10 (" << conn.size() << ")\n";
      }
      out << 10;
      for (unsigned j = 0; j < 10 && j < conn.size(); ++j) out << " " << conn[j];
      out << "\n";
      break;
    }

    case 2: {
      // WEDGE21 exported as 18-node vtkBiQuadraticQuadraticWedge (0..17)
      if (conn.size() < 18) {
        std::cerr << "writeMeshFieldVTK: Wedge element with conn.size() < 18 (" << conn.size() << ")\n";
      }
      out << 18;
      for (unsigned j = 0; j < 18 && j < conn.size(); ++j) out << " " << conn[j];
      out << "\n";
      break;
    }

    case 3: {
      // QUAD9 exported as 9-node VTK_BIQUADRATIC_QUAD (0..8)
      if (conn.size() != 9) {
        std::cerr << "writeMeshFieldVTK: Quad element with conn.size() != 9 (" << conn.size() << ")\n";
      }
      out << 9;
      for (unsigned j = 0; j < 9 && j < conn.size(); ++j) out << " " << conn[j];
      out << "\n";
      break;
    }

    case 4: {
      // TRI7 exported as 6-node VTK_QUADRATIC_TRIANGLE (0..5)
      if (conn.size() < 6) {
        std::cerr << "writeMeshFieldVTK: Tri element with conn.size() < 6 (" << conn.size() << ")\n";
      }
      out << 6;
      for (unsigned j = 0; j < 6 && j < conn.size(); ++j) out << " " << conn[j];
      out << "\n";
      break;
    }

    default:
      break;
    }
  }

  // ------------------------------------------------------
  // CELL_TYPES
  // ------------------------------------------------------
  out << "CELL_TYPES " << numCells << "\n";
  for (std::size_t e = 0; e < numCells; ++e) {
    switch (elType[e]) {
    case 0: out << 29 << "\n"; break; // VTK_TRIQUADRATIC_HEXAHEDRON
    case 1: out << 24 << "\n"; break; // VTK_QUADRATIC_TETRA
    case 2: out << 32 << "\n"; break; // VTK_BIQUADRATIC_QUADRATIC_WEDGE
    case 3: out << 28 << "\n"; break; // VTK_BIQUADRATIC_QUAD
    case 4: out << 22 << "\n"; break; // VTK_QUADRATIC_TRIANGLE
    default:
      std::cerr << "writeMeshFieldVTK: unknown elType in CELL_TYPES, e=" << e
                << ", type=" << elType[e] << "\n";
      out << 0 << "\n";
      break;
    }
  }

  // ------------------------------------------------------
  // CELL_DATA: level + elemental fields
  // ------------------------------------------------------
  bool hasElementalField = false;
  for (std::size_t p = 0; p < field.numFields(); ++p) {
    const unsigned fid = field.idByPos(p);
    if (field.location(fid) == Field::Location::Elemental) {
      hasElementalField = true;
      break;
    }
  }

  const bool canWriteLevel = (level.size() == numCells);

  if (canWriteLevel || hasElementalField) {
    out << "CELL_DATA " << numCells << "\n";

    if (canWriteLevel) {
      out << "SCALARS level int 1\n";
      out << "LOOKUP_TABLE default\n";
      for (std::size_t e = 0; e < numCells; ++e) out << level[e] << "\n";
    }
    else {
      std::cerr << "writeMeshFieldVTK: level.size() != numCells ("
                << level.size() << " != " << numCells << "), skipping level\n";
    }

    for (std::size_t p = 0; p < field.numFields(); ++p) {
      const unsigned fid = field.idByPos(p);
      if (field.location(fid) != Field::Location::Elemental) continue;

      const auto& v = field.getById(fid);
      if (v.size() != numCells) {
        std::cerr << "writeMeshFieldVTK: elemental field '" << field.name(fid)
                  << "' size mismatch (" << v.size() << " != " << numCells << "), skipping\n";
        continue;
      }

      out << "SCALARS " << field.name(fid) << " double 1\n";
      out << "LOOKUP_TABLE default\n";
      for (std::size_t e = 0; e < numCells; ++e) out << v[e] << "\n";
    }
  }

  // ------------------------------------------------------
  // POINT_DATA: nodal fields
  // ------------------------------------------------------
  bool hasNodalField = false;
  for (std::size_t p = 0; p < field.numFields(); ++p) {
    const unsigned fid = field.idByPos(p);
    if (field.location(fid) == Field::Location::Nodal) {
      hasNodalField = true;
      break;
    }
  }

  if (hasNodalField) {
    out << "POINT_DATA " << numPoints << "\n";

    for (std::size_t p = 0; p < field.numFields(); ++p) {
      const unsigned fid = field.idByPos(p);
      if (field.location(fid) != Field::Location::Nodal) continue;

      const auto& v = field.getById(fid);
      if (v.size() != numPoints) {
        std::cerr << "writeMeshFieldVTK: nodal field '" << field.name(fid)
                  << "' size mismatch (" << v.size() << " != " << numPoints << "), skipping\n";
        continue;
      }

      out << "SCALARS " << field.name(fid) << " double 1\n";
      out << "LOOKUP_TABLE default\n";
      for (std::size_t i = 0; i < numPoints; ++i) out << v[i] << "\n";
    }
  }

  out.close();
  if (hasElementalField || hasNodalField) std::cout << "Mesh+Field written to " << filename << " (VTK UNSTRUCTURED_GRID)\n";
  else std::cout << "Mesh written to " << filename << " (VTK UNSTRUCTURED_GRID)\n";
}

void writeMeshVTK(const std::string& filename, Mesh& mesh) {
  Field field(mesh);
  writeMeshFieldVTK(filename, field);
}




static void writePointsVTK(const std::string& filename,
                           const std::vector<std::vector<double>>& Xp)
{
  const std::size_t dim = Xp.size();
  if (dim == 0) throw std::runtime_error("writePointsVTK: Xp.size()==0");

  const std::size_t nPts = Xp[0].size();
  for (std::size_t d = 1; d < dim; ++d) {
    if (Xp[d].size() != nPts) throw std::runtime_error("writePointsVTK: inconsistent Xp[d].size()");
  }
  if (dim > 3) throw std::runtime_error("writePointsVTK: dim > 3 not supported");

  std::ofstream out(filename);
  if (!out) throw std::runtime_error("writePointsVTK: cannot open file");

  out << "# vtk DataFile Version 3.0\n";
  out << "Point cloud\n";
  out << "ASCII\n";
  out << "DATASET POLYDATA\n";
  out << "POINTS " << nPts << " double\n";

  for (std::size_t i = 0; i < nPts; ++i) {
    const double x = (dim >= 1) ? Xp[0][i] : 0.0;
    const double y = (dim >= 2) ? Xp[1][i] : 0.0;
    const double z = (dim >= 3) ? Xp[2][i] : 0.0;
    out << x << " " << y << " " << z << "\n";
  }

  // Create one vertex cell per point (so ParaView can render them).
  out << "VERTICES " << nPts << " " << (2 * nPts) << "\n";
  for (std::size_t i = 0; i < nPts; ++i) {
    out << "1 " << i << "\n";
  }
}


