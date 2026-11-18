#pragma once

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

// Diagnostic-only I/O: no PETSc/FEMuS dependency and no alteration of KK/RES.
// Every rank writes element contributions using global FEMuS system DOF numbers.
// Repeated rows and triplets MUST be summed by post-processing.
class SurfacePressureRangeDiagnostic {
public:
  void Begin(const std::string& prefix, int rank, unsigned dim,
             double pressureJump, unsigned surfaceWeights) {
    if (capturing_ || complete_) {
      throw std::logic_error("Surface range diagnostic already started");
    }
    if (prefix.empty() || !std::isfinite(pressureJump) || pressureJump <= 0) {
      throw std::invalid_argument("Invalid surface range diagnostic metadata");
    }
    pressureJump_ = pressureJump;
    const std::string suffix = "_rank" + std::to_string(rank) + ".csv";
    open(bt_,     prefix + "_bt"     + suffix, "row,col,value\n");
    open(force_,  prefix + "_forces" + suffix, "row,value\n");
    open(action_, prefix + "_action" + suffix, "row,value\n");
    open(dofs_,   prefix + "_dofs"   + suffix, "sys_dof,field,bdc,fixed_value\n");
    if (rank == 0) {
      std::ofstream meta(prefix + "_meta.csv", std::ios::trunc);
      if (!meta) throw std::runtime_error("Cannot write surface range metadata");
      meta << "key,value\n"
           << "dim," << dim << '\n'
           << "pressure_jump," << std::setprecision(17) << pressureJump << '\n'
           << "surface_weights," << surfaceWeights << '\n';
    }
    capturing_ = true;
  }

  bool IsCapturing() const noexcept { return capturing_; }
  bool IsComplete() const noexcept { return complete_; }
  double PressureJump() const noexcept { return pressureJump_; }

  void AddDof(unsigned sysDof, const std::string& field,
              double bdc, double fixedValue) {
    requireCapturing();
    if (!std::isfinite(bdc) || !std::isfinite(fixedValue)) {
      throw std::runtime_error("Non-finite DOF metadata");
    }
    dofs_ << sysDof << ',' << field << ',' << bdc << ',' << fixedValue << '\n';
  }

  // The operator is G = d(RES_velocity)/d(pressure) = -Jac_velocity,pressure.
  // The surface vector and the analytical pressure action are in RES sign.
  void AddElement(const std::vector<unsigned>& sysDof, unsigned nV,
                  unsigned nP, const std::vector<double>& surfaceResidual,
                  const std::vector<double>& pressureOperator,
                  const std::vector<double>& pressureActionYL) {
    requireCapturing();
    if (sysDof.size() != nV + nP || surfaceResidual.size() != nV ||
        pressureOperator.size() != static_cast<size_t>(nV)*nP ||
        pressureActionYL.size() != nV) {
      throw std::runtime_error("Surface range diagnostic: invalid local dimensions");
    }
    for (unsigned i = 0; i < nV; ++i) {
      if (!std::isfinite(surfaceResidual[i]) || !std::isfinite(pressureActionYL[i])) {
        throw std::runtime_error("Surface range diagnostic: non-finite force/action");
      }
      if (surfaceResidual[i] != 0.) {
        force_ << sysDof[i] << ',' << surfaceResidual[i] << '\n';
      }
      if (pressureActionYL[i] != 0.) {
        action_ << sysDof[i] << ',' << pressureActionYL[i] << '\n';
      }
      for (unsigned j = 0; j < nP; ++j) {
        const double x = pressureOperator[static_cast<size_t>(i)*nP+j];
        if (!std::isfinite(x)) {
          throw std::runtime_error("Surface range diagnostic: non-finite operator");
        }
        if (x != 0.) {
          bt_ << sysDof[i] << ',' << sysDof[nV+j] << ',' << x << '\n';
        }
      }
    }
  }

  void FinalizeAssembly() {
    requireCapturing();
    bt_.close(); force_.close(); action_.close(); dofs_.close();
    capturing_ = false;
    complete_ = true;
  }

private:
  static void open(std::ofstream& stream, const std::string& filename,
                   const char* header) {
    stream.open(filename, std::ios::trunc);
    if (!stream) throw std::runtime_error("Cannot open diagnostic file: " + filename);
    stream << std::setprecision(17) << header;
  }
  void requireCapturing() const {
    if (!capturing_ || complete_) {
      throw std::logic_error("Surface range diagnostic is not capturing");
    }
  }
  std::ofstream bt_, force_, action_, dofs_;
  double pressureJump_ = 0.;
  bool capturing_ = false;
  bool complete_ = false;
};
