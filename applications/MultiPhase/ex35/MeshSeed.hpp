
#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>

struct MeshSeed {
  std::vector<unsigned>               elLevel;
  std::vector<unsigned>               elType;
  std::vector<std::vector<unsigned>>  elTplgy;
  std::vector<std::vector<double>>    X;   // dim x nNodes (dim=2 or 3)
};

class MeshSeedOps {
  public:
    // In-place shift: works for 2D or 3D (dims inferred from X.size()).
    static void shift(MeshSeed& s, double dx, double dy, double dz = 0.0) {
      const std::size_t dim = s.X.size();
      if (dim < 2 || dim > 3) throw std::runtime_error("MeshSeedOps::shift: X must have dim 2 or 3");
      for (double& x : s.X[0]) x += dx;
      for (double& y : s.X[1]) y += dy;
      if (dim == 3) for (double& z : s.X[2]) z += dz;
    }

    // In-place scale about origin: works for 2D or 3D
    static void scale(MeshSeed& s, double sx, double sy, double sz = 1.0) {
      const std::size_t dim = s.X.size();
      if (dim < 2 || dim > 3) throw std::runtime_error("MeshSeedOps::scale: X must have dim 2 or 3");
      for (double& x : s.X[0]) x *= sx;
      for (double& y : s.X[1]) y *= sy;
      if (dim == 3) for (double& z : s.X[2]) z *= sz;
    }

    // 2D rotation about origin by theta (radians): (x,y) -> (c x - s y, s x + c y)
    static void rotate2D(MeshSeed& s, double theta) {
      if (s.X.size() != 2) throw std::runtime_error("MeshSeedOps::rotate2D: X must have dim 2");
      const double c = std::cos(theta);
      const double sn = std::sin(theta);
      auto& X = s.X[0];
      auto& Y = s.X[1];
      for (std::size_t i = 0; i < X.size(); ++i) {
        const double x = X[i], y = Y[i];
        X[i] = c * x - sn * y;
        Y[i] = sn * x + c * y;
      }
    }

    // 3D rotation about axis through origin with direction <a,b,c>, angle theta (radians).
    // Uses Rodrigues' rotation formula. Axis is normalized internally.
    static void rotate3D(MeshSeed& s, double a, double b, double c, double theta) {
      if (s.X.size() != 3) throw std::runtime_error("MeshSeedOps::rotate3D: X must have dim 3");

      const double n2 = a * a + b * b + c * c;
      if (n2 == 0.0) throw std::runtime_error("MeshSeedOps::rotate3D: axis must be nonzero");
      const double invn = 1.0 / std::sqrt(n2);
      const double ux = a * invn, uy = b * invn, uz = c * invn;

      const double ct = std::cos(theta);
      const double st = std::sin(theta);
      const double one = 1.0 - ct;

      auto& X = s.X[0];
      auto& Y = s.X[1];
      auto& Z = s.X[2];

      for (std::size_t i = 0; i < X.size(); ++i) {
        const double x = X[i], y = Y[i], z = Z[i];

        // u x v
        const double cx = uy * z - uz * y;
        const double cy = uz * x - ux * z;
        const double cz = ux * y - uy * x;

        // u · v
        const double dot = ux * x + uy * y + uz * z;

        X[i] = x * ct + cx * st + ux * dot * one;
        Y[i] = y * ct + cy * st + uy * dot * one;
        Z[i] = z * ct + cz * st + uz * dot * one;
      }
    }
};


class MeshSeedFactory {
  public:
    enum class Type {
      SquareQuad9,
      SquareTri7,
      CubeHex27,
      CubeWedge21
    };

    static MeshSeed make(Type t) {
      switch (t) {
      case Type::SquareQuad9:  return square_quad9();
      case Type::SquareTri7: return square_tri7();
      case Type::CubeHex27: return cube_hex27();
      case Type::CubeWedge21: return cube_wedge21();
      }
      return {}; // unreachable, quiet compilers
    }

  private:
    //unit square centered at 0,0
    static MeshSeed square_quad9() {
      MeshSeed s;

      s.elLevel = {0};
      s.elType  = {3};
      s.elTplgy = { {0, 1, 2, 3, 4, 5, 6, 7, 8} };
      s.X = {
        { -0.5,  0.5,  0.5, -0.5,  0.0,  0.5,  0.0, -0.5,  0.0},
        { -0.5, -0.5,  0.5,  0.5, -0.5,  0.0,  0.5,  0.0,  0.0}
      };

      return s;
    }

    //unit square centered at 0,0
    static MeshSeed square_tri7() {
      MeshSeed s;

      s.elLevel = {0, 0};
      s.elType  = {4, 4};
      s.elTplgy = {
        {0, 1, 3, 4, 8, 7, 9},
        {1, 2, 3, 5, 6, 8, 10}
      };
      s.X = {
        { -0.5,  0.5,  0.5, -0.5,  0.0,  0.5,  0.0, -0.5,  0.0,  -1.0 / 6.0,  1.0 / 6.0 },
        { -0.5, -0.5,  0.5,  0.5, -0.5,  0.0,  0.5,  0.0,  0.0,  -1.0 / 6.0,  1.0 / 6.0 }
      };

      return s;
    }

    //unit cube centered at 0,0,0
    static MeshSeed cube_hex27() {
      MeshSeed s;

      s.elType  = {0};
      s.elLevel = {0};
      s.elTplgy = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26}};

      s.X = {
        {
          // x (shifted by -0.5 already)
          -0.5,  0.5,  0.5, -0.5, -0.5,  0.5,  0.5, -0.5,
            0.0,  0.5,  0.0, -0.5,  0.0,  0.5,  0.0, -0.5, -0.5,  0.5,  0.5, -0.5,
            0.0,  0.5,  0.0, -0.5,  0.0,  0.0,
            0.0
          },
        {
          // y (shifted by -0.5 already)
          -0.5, -0.5,  0.5,  0.5, -0.5, -0.5,  0.5,  0.5,
            -0.5,  0.0,  0.5,  0.0, -0.5,  0.0,  0.5,  0.0, -0.5, -0.5,  0.5,  0.5,
            -0.5,  0.0,  0.5,  0.0,  0.0,  0.0,
            0.0
          },
        {
          // z (shifted by -0.5 already)
          -0.5, -0.5, -0.5, -0.5,  0.5,  0.5,  0.5,  0.5,
            -0.5, -0.5, -0.5, -0.5,  0.5,  0.5,  0.5,  0.5,  0.0,  0.0,  0.0,  0.0,
            0.0,  0.0,  0.0,  0.0, -0.5,  0.5,
            0.0
          }
      };

      return s;
    }


    static MeshSeed cube_wedge21() {
      MeshSeed s;

      s.elType  = {2, 2};
      s.elLevel = {0, 0};
      s.elTplgy = {
        {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20},
        {21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41}
      };

      s.X = {
        {
          // X
          // Wedge0
          -0.5,  0.5, -0.5, -0.5,  0.5, -0.5,
            0.0,  0.0, -0.5,  0.0,  0.0, -0.5, -0.5,  0.5, -0.5,
            0.0,  0.0, -0.5, -1.0 / 6.0, -1.0 / 6.0,
            -1.0 / 6.0,
            // Wedge1
            -0.5,  0.5, -0.5, -0.5,  0.5, -0.5,
            0.0,  0.0, -0.5,  0.0,  0.0, -0.5, -0.5,  0.5, -0.5,
            0.0,  0.0, -0.5, -1.0 / 6.0, -1.0 / 6.0,
            -1.0 / 6.0
          },
        {
          // Y
          // Wedge0
          -0.5, -0.5,  0.5, -0.5, -0.5,  0.5,
            -0.5,  0.0,  0.0, -0.5,  0.0,  0.0, -0.5, -0.5,  0.5,
            -0.5,  0.0,  0.0, -1.0 / 6.0, -1.0 / 6.0,
            -1.0 / 6.0,
            // Wedge1
            -0.5, -0.5,  0.5, -0.5, -0.5,  0.5,
            -0.5,  0.0,  0.0, -0.5,  0.0,  0.0, -0.5, -0.5,  0.5,
            -0.5,  0.0,  0.0, -1.0 / 6.0, -1.0 / 6.0,
            -1.0 / 6.0
          },
        {
          // Z
          // Wedge0
          -0.5, -0.5, -0.5,  0.5,  0.5,  0.5,
            -0.5, -0.5, -0.5,  0.5,  0.5,  0.5,  0.0,  0.0,  0.0,
            0.0,  0.0,  0.0, -0.5,  0.5,
            0.0,
            // Wedge1
            -0.5, -0.5, -0.5,  0.5,  0.5,  0.5,
            -0.5, -0.5, -0.5,  0.5,  0.5,  0.5,  0.0,  0.0,  0.0,
            0.0,  0.0,  0.0, -0.5,  0.5,
            0.0
          }
      };

      // Rotate ONLY 2nd element nodes (21..41): (x,y)->(-x,-y), z unchanged
      for (unsigned i = 21; i < 42; ++i) {
        s.X[0][i] = -s.X[0][i];
        s.X[1][i] = -s.X[1][i];
      }

      return s;
    }
};







