// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2023 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: sameeragarwal@google.com (Sameer Agarwal)
//
// This example fits the curve f(x;m,c) = e^(m * x + c) to data, minimizing the
// sum squared loss.
#include "ceres/ceres.h"
#include "glog/logging.h"
// Data generated using the following octave code.
//   randn('seed', 23497);
//   m = 0.3;
//   c = 0.1;
//   x=[0:0.075:5];
//   y = exp(m * x + c);
//   noise = randn(size(x)) * 0.2;
//   y_observed = y + noise;
//   data = [x', y_observed'];
const int kNumObservations = 7;
// clang-format off
const double data[] = {
  -3., 10., 2 * (-3.) / sqrt (4. * 9. + 1.), -1. / sqrt (4. * 9. + 1.),
    -2., 5.,  2 * (-2.) / sqrt (4. * 4. + 1.), -1. / sqrt (4. * 4. + 1.),
    -1., 2.,  2 * (-1.) / sqrt (4. * 1. + 1.), -1. / sqrt (4. * 1. + 1.),
    0., 1.,  2 * (-0.) / sqrt (4. * 0. + 1.), -1. / sqrt (4. * 0. + 1),
    1., 2.,  2 * (1.) / sqrt (4. * 1. + 1.), -1. / sqrt (4. * 1. + 1),
    2., 5.,  2 * (2.) / sqrt (4. * 4. + 1.), -1. / sqrt (4. * 4. + 1),
    3., 10., 2 * (3.) / sqrt (4. * 9. + 1.), -1. / sqrt (4. * 9. + 1)
  };


// const double data[] = {
//   -3., 10., -2*(-3.)/sqrt(4. * 9. + 1.), 1./sqrt(4. * 9. + 1.),
//   -2., 5.,  -2*(-2.)/sqrt(4. * 4. + 1.), 1./sqrt(4. * 4. + 1.),
//   -1., 2.,  -2*(-1.)/sqrt(4. * 1. + 1.), 1./sqrt(4. * 1. + 1.),
//    0., 1.,  -2*(-0.)/sqrt(4. * 0. + 1.), 1./sqrt(4. * 0. + 1),
//    1., 2.,  -2*( 1.)/sqrt(4. * 1. + 1.), 1./sqrt(4. * 1. + 1),
//    2., 5.,  -2*( 2.)/sqrt(4. * 4. + 1.), 1./sqrt(4. * 4. + 1),
//    3., 10., -2*( 3.)/sqrt(4. * 9. + 1.), 1./sqrt(4. * 9. + 1)
// };

int main1 (int argc, char** argv);
int main2 (int argc, char** argv);

int main3 (int argc, char** argv) {
  main2 (argc, argv);
  return 0;
}

class ConicCostFunction : public ceres::SizedCostFunction<1, 6> {
  public:
    ConicCostFunction (const double &x, const double &y, const double &w) : _x (x), _y (y), _x2 (x * x), _xy (x * y), _y2 (y * y), _sqrtw (sqrt (w)) {}
    ~ConicCostFunction() {}
    bool Evaluate (double const* const* p, double* residuals, double** jacobians) const override {

      residuals[0] = _sqrtw * (p[0][0] * _x2 + p[0][1] * _xy + p[0][2] * _y2 + p[0][3] * _x + p[0][4] * _y + p[0][5]);
      if (jacobians != nullptr && jacobians[0] != nullptr) {
        jacobians[0][0] = _sqrtw * _x2;
        jacobians[0][1] = _sqrtw * _xy;
        jacobians[0][2] = _sqrtw * _y2;
        jacobians[0][3] = _sqrtw * _x;
        jacobians[0][4] = _sqrtw * _y;
        jacobians[0][5] = _sqrtw;
      }
      return true;
    }
  private:
    const double &_x;
    const double &_y;
    const double _x2;
    const double _xy;
    const double _y2;
    const double _sqrtw;
};

class N1CostFunction : public ceres::SizedCostFunction<1, 6> {
  public:
    N1CostFunction (const double &x, const double &y, const double &n1, const double &n2, const double &w) :
      _x (x), _y (y), _n1 (n1), _n2 (n2), _sqrtw (sqrt (w)) {}
    ~N1CostFunction() {}
    bool Evaluate (double const* const* p, double* residuals, double** jacobians) const override {
      double N1 = 2. * p[0][0] * _x + p[0][1] * _y + p[0][3];
      double N2 = p[0][1] * _x + 2. * p[0][2] * _y + p[0][4];
      double DET2 = N1 * N1 + N2 * N2;
      double DET = sqrt (DET2);
      residuals[0] = _sqrtw * (N1 / DET - _n1);

      if (jacobians != nullptr && jacobians[0] != nullptr) {
        double DET3 = DET2 * DET;
        double C1 = N2 / DET3;
        jacobians[0][0] = 2. * _sqrtw * _x * N2 * C1;
        jacobians[0][1] = _sqrtw * (N2 * _y - N1 * _x)  * C1;
        jacobians[0][2] = -2. * _sqrtw * _y * N1 * C1;
        jacobians[0][3] = _sqrtw * N2 * C1;
        jacobians[0][4] = -_sqrtw * N1 * C1;
        jacobians[0][5] = 0.;
      }
      return true;
    }
  private:
    const double &_x;
    const double &_y;
    const double &_n1;
    const double &_n2;
    const double &_sqrtw;
};


class N2CostFunction : public ceres::SizedCostFunction<1, 6> {

  public:

    N2CostFunction (const double &x, const double &y, const double &n1, const double &n2, const double &w) :
      _x (x), _y (y), _n1 (n1), _n2 (n2), _sqrtw (sqrt (w)) {}
    ~N2CostFunction() {}
    bool Evaluate (double const* const* p, double* residuals, double** jacobians) const override {

      double N1 = 2. * p[0][0] * _x + p[0][1] * _y + p[0][3];
      double N2 = p[0][1] * _x + 2. * p[0][2] * _y + p[0][4];
      double DET2 = N1 * N1 + N2 * N2;
      double DET = sqrt (DET2);
      residuals[0] = _sqrtw * (N2 / DET - _n2);

      if (jacobians != nullptr && jacobians[0] != nullptr) {
        double DET3 = DET2 * DET;
        double C2 = N1 / DET3;
        jacobians[0][0] = -2. * _sqrtw * _x * N2 * C2;
        jacobians[0][1] = _sqrtw * (N1 * _x - N2 * _y)  * C2;
        jacobians[0][2] = 2. * _sqrtw * _y * N1 * C2;
        jacobians[0][3] = -N2 * _sqrtw * C2;
        jacobians[0][4] = _sqrtw * N1 * C2;
        jacobians[0][5] = 0.;
      }
      return true;
    }

  private:
    const double &_x;
    const double &_y;
    const double &_n1;
    const double &_n2;
    const double &_sqrtw;
};




int main1 (int argc, char** argv) {
  google::InitGoogleLogging (argv[0]);
  // The variable to solve for with its initial value. It will be
  // mutated in place by the solver.
  const double initial_a = 1.0;
  const double initial_b = 0.0;
  const double initial_c = 1.0;
  const double initial_d = 0.0;
  const double initial_e = 0.0;
  const double initial_f = 1.0;
  std::vector<double> p = {initial_a, initial_b, initial_c, initial_d, initial_e, initial_f};
  double &a = p[0];
  double &b = p[1];
  double &c = p[2];
  double &d = p[3];
  double &e = p[4];
  double &f = p[5];

  // Build the problem.
  ceres::Problem problem;
  // Set up the only cost function (also known as residual).

  unsigned dim = 2;
  std::vector<double> xg = {0., 0.};
  for (unsigned i = 0; i < kNumObservations; i++) {
    xg[0] += data[4 * i];
    xg[1] += data[4 * i + 1];
  }
  for (unsigned k = 0; k < dim; k++) {
    xg[k] /= kNumObservations;
  }
  std::vector < std::vector<double> >  xp (kNumObservations, std::vector<double> (dim));
  for (unsigned i = 0; i < kNumObservations; i++) {
    xp[i][0] = data[4 * i] - xg[0];
    xp[i][1] = data[4 * i + 1] - xg[1];
  }


  for (int i = 0; i < kNumObservations; ++i) {
    // ceres::CostFunction* ccf = new ConicCostFunction (data[4 * i], data[4 * i + 1], 0.8);
    // problem.AddResidualBlock (ccf, nullptr, p);
    // ceres::CostFunction* n1cf = new N1CostFunction (data[4 * i], data[4 * i + 1], data[4 * i + 2], data[4 * i + 3], .5);
    // problem.AddResidualBlock (n1cf, nullptr, p);
    // ceres::CostFunction* n2cf = new N2CostFunction (data[4 * i], data[4 * i + 1], data[4 * i + 2], data[4 * i + 3], .6);
    // problem.AddResidualBlock (n2cf, nullptr, p);

    ceres::CostFunction* ccf = new ConicCostFunction (xp[i][0], xp[i][1], 0.8);
    problem.AddResidualBlock (ccf, nullptr, &p[0]);
    ceres::CostFunction* n1cf = new N1CostFunction (xp[i][0], xp[i][1], data[4 * i + 2], data[4 * i + 3], .5);
    problem.AddResidualBlock (n1cf, nullptr, &p[0]);
    ceres::CostFunction* n2cf = new N2CostFunction (xp[i][0], xp[i][1], data[4 * i + 2], data[4 * i + 3], .6);
    problem.AddResidualBlock (n2cf, nullptr, &p[0]);

  }
  // Run the solver!
  ceres::Solver::Options options;
  options.max_num_iterations = 25;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = true;
  ceres::Solver::Summary summary;
  ceres::Solve (options, &problem, &summary);
  std::cout << summary.BriefReport() << "\n";


  std::vector<double> q(6);
  double maxD2 = 1;
  double maxD = 1;

  q[0] = p[0] / maxD2;
  q[1] = p[1] / maxD2;
  q[2] = p[2] / maxD2;
  q[3] = p[3] / maxD;
  q[4] = p[4] / maxD;
  q[5] = p[5] + q[0] * xg[0] * xg[0] +  q[1] * xg[0] * xg[1] + q[2] * xg[1] * xg[1] - q[3] * xg[0] - q[4] * xg[1];

  q[3] -= 2. * q[0] * xg[0] + q[1] * xg[1];
  q[4] -= 2. * q[2] * xg[1] + q[1] * xg[0];

  p = q;


  double det = sqrt (a * a + b * b + c * c + d * d + e * e + f * f);
  a /= det;
  b /= det;
  c /= det;
  d /= det;
  e /= det;
  f /= det;

  std::cout << "Initial a: " << initial_a << " Final a: " << a << "\n";
  std::cout << "Initial b: " << initial_b << " Final b: " << b << "\n";
  std::cout << "Initial c: " << initial_c << " Final c: " << c << "\n";
  std::cout << "Initial d: " << initial_d << " Final d: " << d << "\n";
  std::cout << "Initial e: " << initial_e << " Final e: " << e << "\n";
  std::cout << "Initial f: " << initial_f << " Final f: " << f << "\n";

  return 0;
}































// clang-format on
struct ConicResidual {
    ConicResidual (double x, double y) : x_ (x), y_ (y) {}
    template <typename T>
    bool operator() (const T* const a, const T* const b, const T* const c,
                     const T* const d, const T* const e, const T* const f, T*
                     residual) const {

      const T one (-1);
      residual[0] = a[0] * x_ * x_ + b[0] * x_ * y_ + c[0] * y_ * y_ + d[0] *
                    x_ + e[0] * y_ + f[0];
      return true;
    }
  private:
    const double x_;
    const double y_;
};

struct N1Residual {
    N1Residual (double x, double y, double n1, double n2)
      : x_ (x), y_ (y), n1_ (n1), n2_ (n2) {}
    template <typename T>
    bool operator() (const T *const a, const T *const b, const T *const c,
                     const T *const d, const T *const e, const T *const f,
                     T *residual) const {

      const T N1 (2. * a[0] * x_ + b[0] * y_ + d[0]);
      const T N2 (b[0] * x_ + 2. * c[0] * y_ + e[0]);
      residual[0] = N1 / sqrt (N1 * N1 + N2 * N2) - n1_;
      return true;
    }

  private:
    const double x_;
    const double y_;
    const double n1_;
    const double n2_;
};

struct N2Residual {
    N2Residual (double x, double y, double n1, double n2)
      : x_ (x), y_ (y), n1_ (n1), n2_ (n2) {}
    template <typename T>
    bool operator() (const T *const a, const T *const b, const T *const c,
                     const T *const d, const T *const e, const T *const f,
                     T *residual) const {

      const T N1 (2. * a[0] * x_ + b[0] * y_ + d[0]);
      const T N2 (b[0] * x_ + 2. * c[0] * y_ + e[0]);
      residual[0] = N2 / sqrt (N1 * N1 + N2 * N2) - n2_;
      return true;
    }

  private:
    const double x_;
    const double y_;
    const double n1_;
    const double n2_;
};





int main2 (int argc, char **argv) {
  google::InitGoogleLogging (argv[0]);
  const double initial_a = 1.0;
  const double initial_b = 0.0;
  const double initial_c = 1.0;
  const double initial_d = 0.0;
  const double initial_e = 0.0;
  const double initial_f = 1.0;
  double a = initial_a;
  double b = initial_b;
  double c = initial_c;
  double d = initial_d;
  double e = initial_e;
  double f = initial_f;
  ceres::Problem problem;
  for (int i = 0; i < kNumObservations; ++i) {
    problem.AddResidualBlock (
      new ceres::AutoDiffCostFunction<ConicResidual, 1, 1, 1, 1, 1, 1, 1> (
        new ConicResidual (data[4 * i], data[4 * i + 1])),
      nullptr, &a, &b, &c, &d, &e, &f);

    problem.AddResidualBlock (
      new ceres::AutoDiffCostFunction<N1Residual, 1, 1, 1, 1, 1, 1, 1> (
        new N1Residual (data[4 * i], data[4 * i + 1], data[4 * i + 2],
                        data[4 * i + 3])),
      nullptr, &a, &b, &c, &d, &e, &f);

    problem.AddResidualBlock (
      new ceres::AutoDiffCostFunction<N2Residual, 1, 1, 1, 1, 1, 1, 1> (
        new N2Residual (data[4 * i], data[4 * i + 1], data[4 * i + 2],
                        data[4 * i + 3])),
      nullptr, &a, &b, &c, &d, &e, &f);
  }
  ceres::Solver::Options options;
  options.max_num_iterations = 25;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = true;
  ceres::Solver::Summary summary;
  ceres::Solve (options, &problem, &summary);
  std::cout << summary.BriefReport() << "\n";

  double det = sqrt (a * a + b * b + c * c + d * d + e * e + f * f);
  a /= det;
  b /= det;
  c /= det;
  d /= det;
  e /= det;
  f /= det;

  std::cout << "Initial a: " << initial_a << " Final a: " << a << "\n";
  std::cout << "Initial b: " << initial_b << " Final b: " << b << "\n";
  std::cout << "Initial c: " << initial_c << " Final c: " << c << "\n";
  std::cout << "Initial d: " << initial_d << " Final d: " << d << "\n";
  std::cout << "Initial e: " << initial_e << " Final e: " << e << "\n";
  std::cout << "Initial f: " << initial_f << " Final f: " << f << "\n";

  return 0;
}
