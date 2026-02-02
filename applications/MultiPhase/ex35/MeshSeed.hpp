
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

    // ------------------------------------------------------------
    // Reshape: map unit square/cube centered at origin to disk/ball.
    // Input is expected to be in [-0.5,0.5]^dim.
    // Output fits inside radius 0.5 (disk/ball).
    // ------------------------------------------------------------
    static void reshapeBall(MeshSeed& s) {
      const std::size_t dim = s.X.size();
      if (dim != 2 && dim != 3) throw std::runtime_error("MeshSeedOps::reshapeBall: X must have dim 2 or 3");

      auto& X = s.X[0];
      auto& Y = s.X[1];
      const std::size_t N = X.size();
      if (Y.size() != N) throw std::runtime_error("MeshSeedOps::reshapeBall: size mismatch");
      auto* Zp = (dim == 3 ? &s.X[2] : nullptr);
      if (dim == 3 && Zp->size() != N) throw std::runtime_error("MeshSeedOps::reshapeBall: size mismatch (Z)");

      for (std::size_t i = 0; i < N; ++i) {
        // Normalize to [-1,1]
        const double u = 2.0 * X[i];
        const double v = 2.0 * Y[i];

        if (dim == 2) {
          const double m  = std::max(std::fabs(u), std::fabs(v)); // L_inf radius in [0,1]
          const double r2 = u * u + v * v;
          if (r2 == 0.0) continue;
          const double invr = 1.0 / std::sqrt(r2);
          const double u1 = m * u * invr;
          const double v1 = m * v * invr;
          X[i] = 0.5 * u1;
          Y[i] = 0.5 * v1;
        }
        else {
          const double w = 2.0 * (*Zp)[i];
          const double m  = std::max(std::fabs(u), std::max(std::fabs(v), std::fabs(w)));
          const double r2 = u * u + v * v + w * w;
          if (r2 == 0.0) continue;
          const double invr = 1.0 / std::sqrt(r2);
          const double u1 = m * u * invr;
          const double v1 = m * v * invr;
          const double w1 = m * w * invr;
          X[i] = 0.5 * u1;
          Y[i] = 0.5 * v1;
          (*Zp)[i] = 0.5 * w1;
        }
      }
    }

    // ------------------------------------------------------------
    // Reshape: funnel profile with axis in y (2D) or z (3D).
    // Constraints (same as you asked):
    //   at axisCoord = ±0.5: scale = r_end  (default 1.0)
    //   at axisCoord =  0  : scale = r_mid  (default 0.8)
    //
    // 3D: scales (x,y) based on z.
    // 2D: scales x based on y (a "hyperbola-like" funnel).
    //
    // profile="hyperbolic" uses rational curve (hyperbola-like);
    // profile="quadratic" uses simple t^2.
    // ------------------------------------------------------------
    enum class FunnelProfile { Quadratic, Hyperbolic };

    static void reshapeFunnel(MeshSeed& s,
                              double r_mid = 0.8,
                              double r_end = 1.0,
                              FunnelProfile profile = FunnelProfile::Hyperbolic) {
      const std::size_t dim = s.X.size();
      if (dim != 2 && dim != 3) throw std::runtime_error("MeshSeedOps::reshapeFunnel: X must have dim 2 or 3");
      if (!(r_mid > 0.0) || !(r_end > 0.0)) throw std::runtime_error("MeshSeedOps::reshapeFunnel: radii must be > 0");

      auto& X = s.X[0];
      auto& Y = s.X[1];
      const std::size_t N = X.size();
      if (Y.size() != N) throw std::runtime_error("MeshSeedOps::reshapeFunnel: size mismatch");

      auto* Zp = (dim == 3 ? &s.X[2] : nullptr);
      if (dim == 3 && Zp->size() != N) throw std::runtime_error("MeshSeedOps::reshapeFunnel: size mismatch (Z)");

      // axis coordinate and half-range are always 0.5 in your seeds
      constexpr double half = 0.5;

      // Precompute parameter for hyperbolic profile so that:
      // s(t) = r_mid + (r_end-r_mid) * t^2 / (t^2 + a)
      // gives s(0)=r_mid and s(1)=r_end exactly for any a>0, and shape depends on a.
      // Pick a so that mid growth is moderate; a=0.25 is a good default.
      const double a = 0.25;

      for (std::size_t i = 0; i < N; ++i) {
        const double axis = (dim == 2 ? Y[i] : (*Zp)[i]);

        double t = std::fabs(axis) / half;  // nominally in [0,1]
        if (t > 1.0) t = 1.0;               // safety clamp

        double sxy;
        if (profile == FunnelProfile::Quadratic) {
          sxy = r_mid + (r_end - r_mid) * (t * t);
        }
        else {   // Hyperbolic (rational)
          const double t2 = t * t;
          sxy = r_mid + (r_end - r_mid) * (t2 / (t2 + a));
        }

        if (dim == 2) {
          // hyperbola-like funnel: scale x based on y
          X[i] *= sxy;
        }
        else {
          // circular funnel: scale (x,y) based on z
          X[i] *= sxy;
          Y[i] *= sxy;
        }
      }
    }
};


class MeshSeedFactory {
  public:
    enum class Type {
      SquareQuad9,
      SquareTri7,
      CubeHex27,
      CubeWedge21,
      CubeTet15,
      SphereTet15
    };

    enum class Reshape {
      Identity,
      Ball,
      Funnel
    };

    static MeshSeed make(Type t, Reshape reshape = Reshape::Identity, std::vector<double> shift = {0, 0, 0}) {
      switch (t) {
      case Type::SquareQuad9:  return square_quad9(reshape, shift);
      case Type::SquareTri7: return square_tri7(reshape, shift);
      case Type::CubeHex27: return cube_hex27(reshape, shift);
      case Type::CubeWedge21: return cube_wedge21(reshape, shift);
      case Type::CubeTet15: return cube_tet15(reshape, shift);
      case Type::SphereTet15: return sphere_tet15(reshape, shift);
      }
      return {}; // unreachable, quiet compilers
    }

  private:
    //unit square centered at 0,0
    static MeshSeed square_quad9(Reshape reshape = Reshape::Identity, std::vector<double> shift = {0, 0, 0}) {
      MeshSeed s;

      s.elLevel = {0};
      s.elType  = {3};
      s.elTplgy = { {0, 1, 2, 3, 4, 5, 6, 7, 8} };
      s.X = {
        { -0.5,  0.5,  0.5, -0.5,  0.0,  0.5,  0.0, -0.5,  0.0},
        { -0.5, -0.5,  0.5,  0.5, -0.5,  0.0,  0.5,  0.0,  0.0}
      };

      if (reshape == Reshape::Ball)   MeshSeedOps::reshapeBall(s);
      if (reshape == Reshape::Funnel) MeshSeedOps::reshapeFunnel(s, 0.8, 1.0);

      MeshSeedOps::shift(s, shift[0], shift[1], 0);

      return s;
    }

    //unit square centered at 0,0
    static MeshSeed square_tri7(Reshape reshape = Reshape::Identity, std::vector<double> shift = {0, 0, 0}) {
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

      if (reshape == Reshape::Ball)   MeshSeedOps::reshapeBall(s);
      if (reshape == Reshape::Funnel) MeshSeedOps::reshapeFunnel(s, 0.8, 1.0);

      MeshSeedOps::shift(s, shift[0], shift[1], 0);

      return s;
    }

    //unit cube centered at 0,0,0
    static MeshSeed cube_hex27(Reshape reshape = Reshape::Identity, std::vector<double> shift = {0, 0, 0}) {
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

      if (reshape == Reshape::Ball)   MeshSeedOps::reshapeBall(s);
      if (reshape == Reshape::Funnel) MeshSeedOps::reshapeFunnel(s, 0.8, 1.0);

      MeshSeedOps::shift(s, shift[0], shift[1], shift[2]);

      return s;
    }


    static MeshSeed cube_wedge21(Reshape reshape = Reshape::Identity, std::vector<double> shift = {0, 0, 0}) {
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

      if (reshape == Reshape::Ball)   MeshSeedOps::reshapeBall(s);
      if (reshape == Reshape::Funnel) MeshSeedOps::reshapeFunnel(s, 0.8, 1.0);

      MeshSeedOps::shift(s, shift[0], shift[1], shift[2]);

      return s;
    }



    static MeshSeed cube_tet15(Reshape reshape = Reshape::Identity, std::vector<double> shift = {0, 0, 0}) {
      MeshSeed s;

      s.elType  = {1, 1, 1, 1, 1};
      s.elLevel = {0, 0, 0, 0, 0};

      s.elTplgy = { {0, 1, 2, 3, 8, 9, 10, 11, 12, 13, 26, 27, 28, 29, 30},
        {2, 1, 0, 4, 9, 8, 10, 14, 15, 16, 26, 31, 32, 33, 34},
        {0, 2, 4, 5, 10, 14, 16, 17, 18, 19, 33, 35, 36, 37, 38},
        {1, 0, 4, 6, 8, 16, 15, 20, 21, 22, 32, 39, 40, 41, 42},
        {2, 1, 4, 7, 9, 15, 14, 23, 24, 25, 31, 43, 44, 45, 46},
      };

      s.X = { {-0.5, 0.5, -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, 0, 0, -0.5, -0.5, 0, -0.5, 0, 0.5, 0, -0.5, -0.5, 0, 0.5, 0, 0.5, 0, 0.5, 0.5, -0.166667, -0.166667, -0.166667, -0.5, -0.25, 0.166667, 0.166667, -0.166667, 0, -0.5, -0.166667, -0.166667, -0.25, 0.166667, 0.166667, 0.5, 0.25, 0.166667, 0.5, 0.166667, 0.25},
        {0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, -0.5, 0, -0.5, 0, 0, -0.5, -0.5, 0, 0, 0.5, 0.5, 0, 0.5, 0, 0.5, 0.5, -0.5, -0.5, 0, -0.166667, -0.166667, -0.5, -0.166667, -0.25, -0.166667, 0.166667, 0.166667, 0, 0.166667, 0.166667, 0.5, 0.25, 0.166667, 0.5, 0.166667, 0.25, -0.5, -0.166667, -0.166667, -0.25},
        {-0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, -0.5, 0, 0, -0.5, -0.5, 0, 0.5, 0, 0, 0, 0.5, 0.5, -0.5, -0.5, 0, 0.5, 0, 0.5, -0.166667, -0.5, -0.166667, -0.166667, -0.25, 0.166667, -0.166667, 0.166667, 0, 0.166667, 0.5, 0.166667, 0.25, -0.5, -0.166667, -0.166667, -0.25, 0.166667, 0.166667, 0.5, 0.25},
      };


      if (reshape == Reshape::Ball)   MeshSeedOps::reshapeBall(s);
      if (reshape == Reshape::Funnel) MeshSeedOps::reshapeFunnel(s, 0.8, 1.0);


      MeshSeedOps::shift(s, shift[0], shift[1], shift[2]);

      return s;
    }


    static MeshSeed sphere_tet15(Reshape reshape = Reshape::Identity, std::vector<double> shift = {0, 0, 0}) {
      MeshSeed s;

      s.elType  = {1, 1, 1, 1, 1, 1, 1, 1};
      s.elLevel = {0, 0, 0, 0, 0, 0, 0, 0};


      s.X = { {-4.63511e-34, -0.5, 5.20334e-50, 3.06162e-17, 6.12323e-17, 0.5, 5.20334e-50, -0.25, -0.353553, -2.31755e-34, 1.53081e-17, -0.353553, 2.16489e-17, 3.06162e-17, -2.16489e-17, -0.353553, 0.353553, 0.353553, 0.25, -0.353553, 2.16489e-17, -2.31755e-34, -2.16489e-17, 0.353553, 0.353553, -0.21269, -0.21269, -0.258714, 1.30235e-17, -0.176777, -2.81814e-18, -0.258714, -0.21269, -0.176777, 0.258714, 0.21269, 0.21269, 0.176777, -0.258714, 1.30235e-17, -0.21269, -0.176777, -0.258714, -2.81814e-18, -0.176777, 0.21269, 0.21269, 0.258714, 0.176777, 0.258714, 0.176777, 0.258714, 0.176777},
        {-4.63511e-34, -2.36515e-50, 0.5, 3.06162e-17, -3.06162e-17, -2.36515e-50, -0.5, -2.31755e-34, 0.353553, 0.25, 1.53081e-17, 2.16489e-17, 0.353553, -1.53081e-17, 0.353553, -2.16489e-17, 0.353553, -2.16489e-17, -2.31755e-34, -0.353553, -0.353553, -0.25, -0.353553, -0.353553, 2.16489e-17, 0.21269, 1.30235e-17, 0.258714, 0.21269, 0.176777, 0.21269, 0.258714, -1.30235e-17, 0.176777, 0.258714, -1.30235e-17, 0.21269, 0.176777, -0.258714, -0.21269, -0.21269, -0.176777, -0.258714, -0.21269, -0.176777, -0.21269, 1.30235e-17, -0.258714, -0.176777, -0.258714, -0.176777, 0.258714, 0.176777},
        {0, 0, 0, -0.5, 0.5, 0, 0, 0, 0, 0, -0.25, -0.353553, -0.353553, 0.25, 0.353553, 0.353553, 0, 0.353553, 0, 0, -0.353553, 0, 0.353553, 0, -0.353553, 0, -0.21269, -0.258714, -0.21269, -0.176777, 0.21269, 0.258714, 0.21269, 0.176777, 0.258714, 0.21269, 0, 0.176777, -0.258714, -0.21269, 0, -0.176777, 0.258714, 0.21269, 0.176777, 0, -0.21269, -0.258714, -0.176777, 0.258714, 0.176777, -0.258714, -0.176777},
      };

      // s.X = { { -4.63511e-34, -0.5, 5.20334e-50, 3.06162e-17, 6.12323e-17, 0.5, 5.20334e-50, -0.25, -0.353553390593, -2.31755e-34, 1.53081e-17, -0.353553390593, 2.1648923917e-17, 3.06162e-17, -2.1648923917e-17, -0.353553390593, 0.353553390593, 0.353553390593, 0.25, -0.353553390593, 2.1648923917e-17, -2.31755e-34, -2.1648923917e-17, 0.353553390593, 0.353553390593, -0.21269, -0.21269, -0.288675134595, 1.30235e-17, -0.176777, -2.81814e-18, -0.288675134595, -0.21269, -0.176777, 0.288675134595, 0.21269, 0.21269, 0.176777, -0.288675134595, 1.30235e-17, -0.21269, -0.176777, -0.288675134595, -2.81814e-18, -0.176777, 0.21269, 0.21269, 0.288675134595, 0.176777, 0.288675134595, 0.176777, 0.288675134595, 0.176777 },
      //   { -4.63511e-34, -2.36515e-50, 0.5, 3.06162e-17, -3.06162e-17, -2.36515e-50, -0.5, -2.31755e-34, 0.353553390593, 0.25, 1.53081e-17, 2.1648923917e-17, 0.353553390593, -1.53081e-17, 0.353553390593, -2.1648923917e-17, 0.353553390593, -2.1648923917e-17, -2.31755e-34, -0.353553390593, -0.353553390593, -0.25, -0.353553390593, -0.353553390593, 2.1648923917e-17, 0.21269, 1.30235e-17, 0.288675134595, 0.21269, 0.176777, 0.21269, 0.288675134595, -1.30235e-17, 0.176777, 0.288675134595, -1.30235e-17, 0.21269, 0.176777, -0.288675134595, -0.21269, -0.21269, -0.176777, -0.288675134595, -0.21269, -0.176777, -0.21269, 1.30235e-17, -0.288675134595, -0.176777, -0.288675134595, -0.176777, 0.288675134595, 0.176777 },
      //   { 0, 0, 0, -0.5, 0.5, 0, 0, 0, 0, 0, -0.25, -0.353553390593, -0.353553390593, 0.25, 0.353553390593, 0.353553390593, 0, 0.353553390593, 0, 0, -0.353553390593, 0, 0.353553390593, 0, -0.353553390593, 0, -0.21269, -0.288675134595, -0.21269, -0.176777, 0.21269, 0.288675134595, 0.21269, 0.176777, 0.288675134595, 0.21269, 0, 0.176777, -0.288675134595, -0.21269, 0, -0.176777, 0.288675134595, 0.21269, 0.176777, 0, -0.21269, -0.288675134595, -0.176777, 0.288675134595, 0.176777, -0.288675134595, -0.176777 },
      // };

      s.elTplgy = { {0, 1, 2, 3, 7, 8, 9, 10, 11, 12, 25, 26, 27, 28, 29},
        {0, 2, 1, 4, 9, 8, 7, 13, 14, 15, 25, 30, 31, 32, 33},
        {2, 4, 0, 5, 14, 13, 9, 16, 17, 18, 30, 34, 35, 36, 37},
        {1, 3, 0, 6, 11, 10, 7, 19, 20, 21, 26, 38, 39, 40, 41},
        {1, 4, 6, 0, 15, 22, 19, 7, 13, 21, 42, 32, 43, 40, 44},
        {0, 5, 6, 3, 18, 23, 21, 10, 24, 20, 45, 46, 47, 39, 48},
        {0, 6, 5, 4, 21, 23, 18, 13, 22, 17, 45, 43, 49, 35, 50},
        {2, 3, 5, 0, 12, 24, 16, 9, 10, 18, 51, 28, 46, 36, 52},
      };

      //if (reshape == Reshape::Ball)   MeshSeedOps::reshapeBall(s);
      //if (reshape == Reshape::Funnel) MeshSeedOps::reshapeFunnel(s, 0.8, 1.0);

      MeshSeedOps::shift(s, shift[0], shift[1], shift[2]);

      return s;
    }

};







// femus::FemusInit mpinit(argc, args, MPI_COMM_WORLD);
//
// // define multilevel mesh
// femus::MultiLevelMesh mlMsh;
// double scalingFactor = 1.;
// //mlMsh.ReadCoarseMesh("./input/cube_hex.neu", "fifth", scalingFactor);
// mlMsh.ReadCoarseMesh("./input/cube.neu", "fifth", scalingFactor);
// unsigned dim = mlMsh.GetDimension();
//
// unsigned numberOfUniformLevels = 1;
// unsigned numberOfSelectiveLevels = 0;
// mlMsh.RefineMesh(numberOfUniformLevels, numberOfUniformLevels + numberOfSelectiveLevels, NULL);
//
// femus::MultiLevelSolution mlSol(&mlMsh);
// mlSol.Initialize("All");
//
//
// femus::Mesh* msh = mlMsh.GetLevel(0);
//
// std::cout << "s.X = {";
// for (unsigned d = 0; d < 3; ++d) {
//   std::cout << " {";
//   for (unsigned i = msh->_dofOffset[2][0]; i < msh->_dofOffset[2][1]; i++) {
//     std::cout << (*msh->_topology->_Sol[d])(i) << ", ";
//   }
//   std::cout << "\b\b},\n";
// }
// std::cout << "\b\b};\n\n";
//
// std::cout << "s.elTplgy = {";
// for (unsigned iel = msh->_elementOffset[0]; iel < msh->_elementOffset[1]; iel++) {
//   short unsigned ielGeom = msh->GetElementType(iel);
//   unsigned nDofs = msh->GetElementDofNumber(iel, 2);
//   std::cout << " {";
//   for (unsigned i = 0; i < nDofs; i++) {
//     std::cout << msh->GetSolutionDof(i, iel, 2) << ", ";
//   }
//   std::cout << "\b\b},\n";
// }
// std::cout << "\b\b};\n\n";
//
//
//
// std::vector < std::string > variablesToBePrinted;
// variablesToBePrinted.push_back("All");
// femus::VTKWriter vtkIO(&mlSol);
// vtkIO.SetDebugOutput(true);
// vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);
//
//
// return 0;







// elType0[0] = {0, 2, 1};
// elLevel0[0] = {0, 0, 0};
// elTplgy0[0] = {
//   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26},
//   {27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47},
//   {48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62}
// };
//
// X0[0] = {{
//     //hex
//     0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
//     0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.0, 1.0, 1.0, 0.0,
//     0.5, 1.0, 0.5, 0.0, 0.5, 0.5,
//     0.5,
//     //wedge
//     1.0, 2.0, 1.0, 1.0, 2.0, 1.0,
//     1.5, 1.5, 1.0, 1.5, 1.5, 1.0, 1.0, 2.0, 1.0,
//     1.5, 1.5, 1.0, 1.3333333333333333, 1.3333333333333333,
//     1.3333333333333333,
//     //tet
//     1.0, 2.0, 1.0, 1.0,
//     1.5, 1.5, 1.0, 1.0, 1.5, 1.0,
//     1.3333333333333333, 1.3333333333333333,
//     1.3333333333333333, 1.0,
//     1.25
//   }, {
//     //hex
//     0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
//     0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.0, 1.0, 1.0,
//     0.0, 0.5, 1.0, 0.5, 0.5, 0.5,
//     0.5,
//     //wedge
//     0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
//     0.0, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.0, 1.0,
//     0.0, 0.5, 0.5, 1.0 / 3.0, 1.0 / 3.0,
//     1.0 / 3.0,
//     //tet
//     0.0, 0.0, 1.0, 0.0,
//     0.0, 0.5, 0.5, 0.0, 0.0, 0.5,
//     0.3333333333333333, 0.0,
//     0.3333333333333333, 0.3333333333333333,
//     0.25
//   }, {
//     //hex
//     0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
//     0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5,
//     0.5, 0.5, 0.5, 0.5, 0.0, 1.0,
//     0.5,
//     //wedge
//     0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
//     0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5,
//     0.5, 0.5, 0.5, 0.0, 1.0,
//     0.5,
//     //tet
//     1.0, 1.0, 1.0, 2.0,
//     1.0, 1.0, 1.0, 1.5, 1.5, 1.5,
//     1.0, 1.3333333333333333,
//     1.3333333333333333, 1.3333333333333333,
//     1.25
//   }
// };




