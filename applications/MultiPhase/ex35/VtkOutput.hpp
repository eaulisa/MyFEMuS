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

//
//   VTK cell types:
//   VTK_TRIQUADRATIC_HEXAHEDRON     = 29
//   VTK_QUADRATIC_TETRA             = 24   (10-node tetra)
//   VTK_BIQUADRATIC_QUADRATIC_WEDGE = 32   (18-node wedge)
//   VTK_BIQUADRATIC_QUAD            = 28
//   VTK_QUADRATIC_TRIANGLE          = 22   (6-node triangle)

void writeMeshVTK(const std::string& filename,
                  const std::vector<std::vector<double>>& X,
                  const std::vector<unsigned>& elType,
                  const std::vector<std::vector<unsigned>>& elTplgy,
                  const std::vector<unsigned>& level) {

  const std::size_t dim = X.size();
  if (dim != 2 && dim != 3) {
    std::cerr << "writeMeshVTK: X must have size 2 (2D) or 3 (3D), got "
              << dim << "\n";
    return;
  }

  if (X[0].empty()) {
    std::cerr << "writeMeshVTK: X[0] is empty\n";
    return;
  }

  const std::size_t numPoints = X[0].size();
  const std::size_t numCells  = elTplgy.size();

  // check all coordinate components have same number of points
  for (std::size_t d = 1; d < dim; ++d) {
    if (X[d].size() != numPoints) {
      std::cerr << "writeMeshVTK: X[" << d << "] size mismatch ("
                << X[d].size() << " != " << numPoints << ")\n";
      return;
    }
  }

  if (elType.size() != numCells) {
    std::cerr << "writeMeshVTK: elType.size() != elTplgy.size() ("
              << elType.size() << " != " << numCells << ")\n";
    return;
  }

  std::ofstream out(filename.c_str());
  if (!out.is_open()) {
    std::cerr << "Cannot open VTK file: " << filename << "\n";
    return;
  }

  // ------------------------------------------------------
  // VTK legacy header
  // ------------------------------------------------------
  out << "# vtk DataFile Version 3.0\n";
  out << "FEM refined mesh\n";
  out << "ASCII\n";
  out << "DATASET UNSTRUCTURED_GRID\n";

  // ------------------------------------------------------
  // POINTS: 2D or 3D (z = 0 if dim == 2)
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
  // CELLS
  // ------------------------------------------------------
  std::size_t totalCellInts = 0;
  for (std::size_t e = 0; e < numCells; ++e) {
    switch (elType[e]) {

    case 0: // Hex27 -> 27-node tri-quadratic hex
      totalCellInts += 1 + 27;
      break;
    case 1: // Tet15 -> export as 10-node quadratic tetra
      totalCellInts += 1 + 10;
      break;
    case 2: // Wedge21 -> export as 18-node wedge
      totalCellInts += 1 + 18;
      break;
    case 3: // Quad9
      totalCellInts += 1 + 9;
      break;
    case 4: // Tri7 -> 6-node quadratic triangle
      totalCellInts += 1 + 6;
      break;
    default:
      std::cerr << "writeMeshVTK: unknown elType[" << e
                << "] = " << elType[e] << "\n";
      return;
    }
  }

  out << "CELLS " << numCells << " " << totalCellInts << "\n";

  for (std::size_t e = 0; e < numCells; ++e) {
    const auto &conn = elTplgy[e];

    switch (elType[e]) {


    case 0: {
      // HEX27: last 7 DOFs remapped with tail_map
      static const std::array<unsigned, 7> tail_map{{3u, 1u, 0u, 2u, 4u, 5u, 6u}};

      if (conn.size() != 27) {
        std::cerr << "writeMeshVTK: Hex element with conn.size() != 27 ("
                  << conn.size() << ")\n";
      }

      out << 27;
      for (unsigned j = 0; j < 27 && j < conn.size(); ++j) {
        unsigned src = j;
        if (j >= 20) {
          const unsigned i = j - 20;          // i in [0,6]
          src = 20u + static_cast<unsigned>(tail_map[i]);
        }
        out << " " << conn[src];
      }
      out << "\n";
      break;
    }

    case 1: {
      // TET15 exported as 10-node VTK_QUADRATIC_TETRA:
      // use nodes 0..9 (4 vertices + 6 edge mids)
      if (conn.size() < 10) {
        std::cerr << "writeMeshVTK: Tet element with conn.size() < 10 ("
                  << conn.size() << ")\n";
      }

      out << 10;
      for (unsigned j = 0; j < 10 && j < conn.size(); ++j) {
        out << " " << conn[j];
      }
      out << "\n";
      break;
    }

    case 2: {
      // WEDGE21 exported as 18-node vtkBiQuadraticQuadraticWedge:
      // nodes 0..17 only
      if (conn.size() < 18) {
        std::cerr << "writeMeshVTK: Wedge element with conn.size() < 18 ("
                  << conn.size() << ")\n";
      }

      out << 18;
      for (unsigned j = 0; j < 18 && j < conn.size(); ++j) {
        out << " " << conn[j];
      }
      out << "\n";
      break;
    }

    case 3: {
      // QUAD9: local ordering 0..8
      if (conn.size() != 9) {
        std::cerr << "writeMeshVTK: Quad element with conn.size() != 9 ("
                  << conn.size() << ")\n";
      }
      out << 9;
      for (unsigned j = 0; j < 9 && j < conn.size(); ++j) {
        out << " " << conn[j];
      }
      out << "\n";
      break;
    }

    case 4: {
      // TRI7: VTK_QUADRATIC_TRIANGLE uses only 6 nodes: 0..5
      if (conn.size() < 6) {
        std::cerr << "writeMeshVTK: Tri element with conn.size() < 6 ("
                  << conn.size() << ")\n";
      }
      out << 6;
      for (unsigned j = 0; j < 6 && j < conn.size(); ++j) {
        out << " " << conn[j];
      }
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
    case 0:
      out << 29 << "\n"; // VTK_TRIQUADRATIC_HEXAHEDRON
      break;
    case 1:
      out << 24 << "\n"; // VTK_QUADRATIC_TETRA
      break;
    case 2:
      out << 32 << "\n"; // VTK_BIQUADRATIC_QUADRATIC_WEDGE
      break;
    case 3:
      out << 28 << "\n"; // VTK_BIQUADRATIC_QUAD
      break;
    case 4:
      out << 22 << "\n"; // VTK_QUADRATIC_TRIANGLE
      break;
    default:
      std::cerr << "writeMeshVTK: unknown elType in CELL_TYPES, e = "
                << e << ", type = " << elType[e] << "\n";
      out << 0 << "\n"; // fallback
      break;
    }
  }

  // ------------------------------------------------------
  // CELL DATA: LEVEL
  // ------------------------------------------------------
  if (level.size() != numCells) {
    std::cerr << "writeMeshVTK: level.size() != numCells ("
              << level.size() << " != " << numCells << ")\n";
  }
  else {
    out << "CELL_DATA " << numCells << "\n";
    out << "SCALARS level int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (std::size_t e = 0; e < numCells; ++e) {
      out << level[e] << "\n";
    }
  }

  out.close();
  std::cout << "Mesh written to " << filename << " (VTK UNSTRUCTURED_GRID)\n";
/*

  // ------------------------------------------------------
  // CELL_TYPES
  // ------------------------------------------------------
  out << "CELL_TYPES " << numCells << "\n";
  for (std::size_t e = 0; e < numCells; ++e) {
    switch (elType[e]) {
    case 0:
      out << 29 << "\n"; // VTK_TRIQUADRATIC_HEXAHEDRON
      break;
    case 1:
      out << 24 << "\n"; // VTK_QUADRATIC_TETRA
      break;
    case 2:
      out << 32 << "\n"; // VTK_BIQUADRATIC_QUADRATIC_WEDGE
      break;
    case 3:
      out << 28 << "\n"; // VTK_BIQUADRATIC_QUAD
      break;
    case 4:
      out << 22 << "\n"; // VTK_QUADRATIC_TRIANGLE
      break;
    default:
      std::cerr << "writeMeshVTK: unknown elType in CELL_TYPES, e = "
                << e << ", type = " << elType[e] << "\n";
      out << 0 << "\n"; // VTK_EMPTY_CELL as fallback
      break;
    }
  }*/

  out.close();
  std::cout << "Mesh written to " << filename << " (VTK UNSTRUCTURED_GRID)\n";
}

void writeMeshVTK(const std::string& filename,
                  const Mesh &mesh){


  const std::vector<std::vector<double>>& X = mesh.X();
  const std::vector<unsigned>& elType = mesh.elType();
  const std::vector<std::vector<unsigned>>& elTplgy = mesh.elTplgy();
  const std::vector<unsigned>& level = mesh.elLevel();
  const std::size_t dim = mesh.dim();

  if (dim != 2 && dim != 3) {
    std::cerr << "writeMeshVTK: X must have size 2 (2D) or 3 (3D), got "
              << dim << "\n";
    return;
  }

  if (X[0].empty()) {
    std::cerr << "writeMeshVTK: X[0] is empty\n";
    return;
  }

  const std::size_t numPoints = X[0].size();
  const std::size_t numCells  = elTplgy.size();

  // check all coordinate components have same number of points
  for (std::size_t d = 1; d < dim; ++d) {
    if (X[d].size() != numPoints) {
      std::cerr << "writeMeshVTK: X[" << d << "] size mismatch ("
                << X[d].size() << " != " << numPoints << ")\n";
      return;
    }
  }

  if (elType.size() != numCells) {
    std::cerr << "writeMeshVTK: elType.size() != elTplgy.size() ("
              << elType.size() << " != " << numCells << ")\n";
    return;
  }

  std::ofstream out(filename.c_str());
  if (!out.is_open()) {
    std::cerr << "Cannot open VTK file: " << filename << "\n";
    return;
  }

  // ------------------------------------------------------
  // VTK legacy header
  // ------------------------------------------------------
  out << "# vtk DataFile Version 3.0\n";
  out << "FEM refined mesh\n";
  out << "ASCII\n";
  out << "DATASET UNSTRUCTURED_GRID\n";

  // ------------------------------------------------------
  // POINTS: 2D or 3D (z = 0 if dim == 2)
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
  // CELLS
  // ------------------------------------------------------
  std::size_t totalCellInts = 0;
  for (std::size_t e = 0; e < numCells; ++e) {
    switch (elType[e]) {

    case 0: // Hex27 -> 27-node tri-quadratic hex
      totalCellInts += 1 + 27;
      break;
    case 1: // Tet15 -> export as 10-node quadratic tetra
      totalCellInts += 1 + 10;
      break;
    case 2: // Wedge21 -> export as 18-node wedge
      totalCellInts += 1 + 18;
      break;
    case 3: // Quad9
      totalCellInts += 1 + 9;
      break;
    case 4: // Tri7 -> 6-node quadratic triangle
      totalCellInts += 1 + 6;
      break;
    default:
      std::cerr << "writeMeshVTK: unknown elType[" << e
                << "] = " << elType[e] << "\n";
      return;
    }
  }

  out << "CELLS " << numCells << " " << totalCellInts << "\n";

  for (std::size_t e = 0; e < numCells; ++e) {
    const auto &conn = elTplgy[e];

    switch (elType[e]) {


    case 0: {
      // HEX27: last 7 DOFs remapped with tail_map
      static const std::array<unsigned, 7> tail_map{{3u, 1u, 0u, 2u, 4u, 5u, 6u}};

      if (conn.size() != 27) {
        std::cerr << "writeMeshVTK: Hex element with conn.size() != 27 ("
                  << conn.size() << ")\n";
      }

      out << 27;
      for (unsigned j = 0; j < 27 && j < conn.size(); ++j) {
        unsigned src = j;
        if (j >= 20) {
          const unsigned i = j - 20;          // i in [0,6]
          src = 20u + static_cast<unsigned>(tail_map[i]);
        }
        out << " " << conn[src];
      }
      out << "\n";
      break;
    }

    case 1: {
      // TET15 exported as 10-node VTK_QUADRATIC_TETRA:
      // use nodes 0..9 (4 vertices + 6 edge mids)
      if (conn.size() < 10) {
        std::cerr << "writeMeshVTK: Tet element with conn.size() < 10 ("
                  << conn.size() << ")\n";
      }

      out << 10;
      for (unsigned j = 0; j < 10 && j < conn.size(); ++j) {
        out << " " << conn[j];
      }
      out << "\n";
      break;
    }

    case 2: {
      // WEDGE21 exported as 18-node vtkBiQuadraticQuadraticWedge:
      // nodes 0..17 only
      if (conn.size() < 18) {
        std::cerr << "writeMeshVTK: Wedge element with conn.size() < 18 ("
                  << conn.size() << ")\n";
      }

      out << 18;
      for (unsigned j = 0; j < 18 && j < conn.size(); ++j) {
        out << " " << conn[j];
      }
      out << "\n";
      break;
    }

    case 3: {
      // QUAD9: local ordering 0..8
      if (conn.size() != 9) {
        std::cerr << "writeMeshVTK: Quad element with conn.size() != 9 ("
                  << conn.size() << ")\n";
      }
      out << 9;
      for (unsigned j = 0; j < 9 && j < conn.size(); ++j) {
        out << " " << conn[j];
      }
      out << "\n";
      break;
    }

    case 4: {
      // TRI7: VTK_QUADRATIC_TRIANGLE uses only 6 nodes: 0..5
      if (conn.size() < 6) {
        std::cerr << "writeMeshVTK: Tri element with conn.size() < 6 ("
                  << conn.size() << ")\n";
      }
      out << 6;
      for (unsigned j = 0; j < 6 && j < conn.size(); ++j) {
        out << " " << conn[j];
      }
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
    case 0:
      out << 29 << "\n"; // VTK_TRIQUADRATIC_HEXAHEDRON
      break;
    case 1:
      out << 24 << "\n"; // VTK_QUADRATIC_TETRA
      break;
    case 2:
      out << 32 << "\n"; // VTK_BIQUADRATIC_QUADRATIC_WEDGE
      break;
    case 3:
      out << 28 << "\n"; // VTK_BIQUADRATIC_QUAD
      break;
    case 4:
      out << 22 << "\n"; // VTK_QUADRATIC_TRIANGLE
      break;
    default:
      std::cerr << "writeMeshVTK: unknown elType in CELL_TYPES, e = "
                << e << ", type = " << elType[e] << "\n";
      out << 0 << "\n"; // fallback
      break;
    }
  }

  // ------------------------------------------------------
  // CELL DATA: LEVEL
  // ------------------------------------------------------
  if (level.size() != numCells) {
    std::cerr << "writeMeshVTK: level.size() != numCells ("
              << level.size() << " != " << numCells << ")\n";
  }
  else {
    out << "CELL_DATA " << numCells << "\n";
    out << "SCALARS level int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (std::size_t e = 0; e < numCells; ++e) {
      out << level[e] << "\n";
    }
  }

  out.close();
  std::cout << "Mesh written to " << filename << " (VTK UNSTRUCTURED_GRID)\n";
/*

  // ------------------------------------------------------
  // CELL_TYPES
  // ------------------------------------------------------
  out << "CELL_TYPES " << numCells << "\n";
  for (std::size_t e = 0; e < numCells; ++e) {
    switch (elType[e]) {
    case 0:
      out << 29 << "\n"; // VTK_TRIQUADRATIC_HEXAHEDRON
      break;
    case 1:
      out << 24 << "\n"; // VTK_QUADRATIC_TETRA
      break;
    case 2:
      out << 32 << "\n"; // VTK_BIQUADRATIC_QUADRATIC_WEDGE
      break;
    case 3:
      out << 28 << "\n"; // VTK_BIQUADRATIC_QUAD
      break;
    case 4:
      out << 22 << "\n"; // VTK_QUADRATIC_TRIANGLE
      break;
    default:
      std::cerr << "writeMeshVTK: unknown elType in CELL_TYPES, e = "
                << e << ", type = " << elType[e] << "\n";
      out << 0 << "\n"; // VTK_EMPTY_CELL as fallback
      break;
    }
  }*/

  out.close();
  std::cout << "Mesh written to " << filename << " (VTK UNSTRUCTURED_GRID)\n";
}

