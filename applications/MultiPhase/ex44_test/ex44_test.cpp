#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "Fem.hpp"
#include "GeomElTypeEnum.hpp"
#include "VTKWriter.hpp"
#include "../includeLS/Psi.hpp"
#include "IntegrateCut.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace femus;

template <class PsiFunction>
void FlagFinestMeshLevel(MultiLevelMesh& mlMsh, const PsiFunction& psi);

template <class PsiFunction>
void InitLevelSet(MultiLevelSolution& mlSol, const std::string& name, const PsiFunction& psi);

namespace {

// Physical-coordinate evaluator of the manufactured BALL level set.
// Geometry: x = q + epsilon*(radius^2 - |q-center|^2)*e_1.
// The inverse recovers q from any physical point x, including on refined meshes.
struct ManufacturedBallPsi {
  std::vector<double> center;
  double radius;
  double epsilon;

  double operator()(const std::vector<double>& x) const {
    if(x.size() != center.size() || x.empty()) {
      throw std::runtime_error("ManufacturedBallPsi: incompatible point dimension");
    }

    double A = radius * radius;
    for(unsigned d = 1; d < x.size(); ++d) {
      const double v = x[d] - center[d];
      A -= v * v;
    }

    // Let u=q[0]-center[0], t=x[0]-center[0].
    // Then t - epsilon*A = u - epsilon*u*u.
    const double y = (x[0] - center[0]) - epsilon * A;
    const double discriminant = 1. - 4. * epsilon * y;

    if(discriminant <= 0. || !std::isfinite(discriminant)) {
      throw std::runtime_error("ManufacturedBallPsi: inverse mapping is not regular");
    }

    // Numerically stable branch continuous with the identity as epsilon -> 0.
    const double u = 2. * y / (1. + std::sqrt(discriminant));
    return A - u * u;
  }
};

struct DriverOptions {
  std::vector<GeomElType> elements = {TET};
  std::vector<PsiType> psis = {PsiType::BALL};
  bool curved = false;
  unsigned levels = 3; // Number of AMR refinements; also evaluate level 0.
  double curvature = 0.25;
  double epsilon = 0.20; // Manufactured BALL mapping strength.
  bool vtk = true;
  IntegrationOptions integration;
};

std::string Lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

std::vector<std::string> Split(const std::string& s) {
  std::vector<std::string> out;
  std::istringstream in(s);
  std::string token;
  while(std::getline(in, token, ',')) {
    token.erase(std::remove_if(token.begin(), token.end(), [](unsigned char c) { return std::isspace(c); }), token.end());
    if(!token.empty()) out.push_back(Lower(token));
  }
  if(out.empty()) throw std::invalid_argument("Empty comma-separated list");
  return out;
}

unsigned Unsigned(const std::string& value, const std::string& option) {
  std::size_t used = 0;
  if(value.empty() || value[0] == '-') throw std::invalid_argument("Invalid " + option + ": " + value);
  const unsigned long x = std::stoul(value, &used);
  if(used != value.size() || x > std::numeric_limits<unsigned>::max()) throw std::invalid_argument("Invalid " + option + ": " + value);
  return static_cast<unsigned>(x);
}

double Real(const std::string& value, const std::string& option) {
  std::size_t used = 0;
  const double x = std::stod(value, &used);
  if(used != value.size() || !std::isfinite(x)) throw std::invalid_argument("Invalid " + option + ": " + value);
  return x;
}

bool Boolean(const std::string& value) {
  const std::string s = Lower(value);
  if(s == "1" || s == "true" || s == "yes" || s == "on") return true;
  if(s == "0" || s == "false" || s == "no" || s == "off") return false;
  throw std::invalid_argument("Expected true/false or 1/0, got: " + value);
}

GeomElType Element(const std::string& value) {
  if(value == "tri") return TRI;
  if(value == "quad") return QUAD;
  if(value == "tet") return TET;
  if(value == "hex") return HEX;
  throw std::invalid_argument("Unknown element: " + value + " (expected tri,quad,tet,hex)");
}

const char* ElementName(GeomElType value) {
  switch(value) {
    case TRI: return "tri";
    case QUAD: return "quad";
    case TET: return "tet";
    case HEX: return "hex";
    default: throw std::invalid_argument("Unsupported element type");
  }
}

PsiType PsiKind(const std::string& value) {
  if(value == "ball") return PsiType::BALL;
  if(value == "saye") return PsiType::SAYE;
  throw std::invalid_argument("Unknown psi type: " + value + " (supported: ball,saye)");
}

const char* PsiName(PsiType value) {
  switch(value) {
    case PsiType::BALL: return "ball";
    case PsiType::SAYE: return "saye";
    default: throw std::invalid_argument("Unsupported psi type");
  }
}

void Usage() {
  std::cout << "Usage: MultiPhase_ex44_test [options]\n"
            << "  --elements tri,quad,tet,hex  Element types (comma-separated)\n"
            << "  --psi ball,saye              Level-set types (comma-separated)\n"
            << "  --curved 0|1                 Use curved physical elements (default 0)\n"
            << "  --levels N                   Number of AMR refinements, plus level 0 (default 3)\n"
            << "  --curvature X                Simplex curvature for curved SAYE (default 0.25)\n"
            << "  --epsilon X                  Manufactured BALL curvature (default 0.20)\n"
            << "  --mode moment|algoim|compare Integration method (default moment)\n"
            << "  --degree N                   Moment-fitting degree (default 7)\n"
            << "  --algoim-order N             Algoim quadrature order (default 5)\n"
            << "  --algoim-subdivisions N      Algoim subdivisions (default 0)\n"
            << "  --vtk 0|1                    Write VTK at each level (default 1)\n"
            << "  --help                       Show this help\n";
}

DriverOptions Parse(int argc, char** argv) {
  DriverOptions options;
  for(int i = 1; i < argc; ++i) {
    const std::string key = argv[i];
    if(key == "--help" || key == "-h") { Usage(); std::exit(0); }
    if(i + 1 >= argc) throw std::invalid_argument("Missing value for " + key);
    const std::string value = argv[++i];
    if(key == "--elements") {
      options.elements.clear();
      for(const auto& token : Split(value)) {
        const GeomElType e = Element(token);
        if(std::find(options.elements.begin(), options.elements.end(), e) == options.elements.end()) options.elements.push_back(e);
      }
    }
    else if(key == "--psi") {
      options.psis.clear();
      for(const auto& token : Split(value)) {
        const PsiType p = PsiKind(token);
        if(std::find(options.psis.begin(), options.psis.end(), p) == options.psis.end()) options.psis.push_back(p);
      }
    }
    else if(key == "--curved") options.curved = Boolean(value);
    else if(key == "--levels") options.levels = Unsigned(value, key);
    else if(key == "--curvature") options.curvature = Real(value, key);
    else if(key == "--epsilon") options.epsilon = Real(value, key);
    else if(key == "--degree") options.integration.momentDegree = Unsigned(value, key);
    else if(key == "--algoim-order") options.integration.algoimOrder = Unsigned(value, key);
    else if(key == "--algoim-subdivisions") options.integration.algoimSubdivisions = Unsigned(value, key);
    else if(key == "--vtk") options.vtk = Boolean(value);
    else if(key == "--mode") {
      const auto s = Lower(value);
      if(s == "moment") options.integration.mode = IntegrationMode::MomentFitting;
      else if(s == "algoim") options.integration.mode = IntegrationMode::Algoim;
      else if(s == "compare") options.integration.mode = IntegrationMode::Compare;
      else throw std::invalid_argument("Unknown integration mode: " + value);
    }
    else throw std::invalid_argument("Unknown option: " + key);
  }
  if(options.levels > 12) throw std::invalid_argument("--levels must be <= 12");
  if(options.integration.momentDegree == 0) throw std::invalid_argument("--degree must be positive");
  if(options.curvature < 0.) throw std::invalid_argument("--curvature must be nonnegative");
  return options;
}

std::vector<IntegrationTest> MakeTests(const PsiType kind, const unsigned dim,
                                       const std::vector<double>& center, const double radius) {
  const double pi = std::acos(-1.);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const bool ball = kind == PsiType::BALL;
  return {
    {"1", [](const std::vector<double>&) { return 1.; },
      ball ? (dim == 2 ? pi * radius * radius : 4. * pi * std::pow(radius, 3) / 3.) : nan,
      ball ? (dim == 2 ? 2. * pi * radius : 4. * pi * radius * radius) : nan},
    {"r^4", [center](const std::vector<double>& x) {
        double r2 = 0.;
        for(unsigned d = 0; d < x.size(); ++d) r2 += (x[d] - center[d]) * (x[d] - center[d]);
        return r2 * r2;
      },
      ball ? (dim == 2 ? pi * std::pow(radius, 6) / 3. : 4. * pi * std::pow(radius, 7) / 7.) : nan,
      ball ? (dim == 2 ? 2. * pi * std::pow(radius, 5) : 4. * pi * std::pow(radius, 6)) : nan}
  };
}

void PrintHeader() {
  std::cout << std::string(90, '-') << '\n';
  std::cout << std::left << std::setw(8) << "Test" << std::setw(9) << "Method"
            << std::setw(9) << "Domain" << std::right << std::setw(16) << "Numerical"
            << std::setw(16) << "Exact" << std::setw(16) << "Abs. error"
            << std::setw(16) << "Rel. error" << '\n';
  std::cout << std::string(90, '-') << '\n';
}

void PrintRow(const std::string& test, const char* method, const char* domain,
              double numerical, double exact) {
  std::cout << std::left << std::setw(8) << test << std::setw(9) << method
            << std::setw(9) << domain << std::right << std::scientific << std::setprecision(8)
            << std::setw(16) << numerical;
  if(std::isfinite(exact)) {
    const double absoluteError = std::abs(numerical - exact);
    std::cout << std::setw(16) << exact << std::setw(16) << absoluteError;
    if(exact != 0.) std::cout << std::setw(16) << absoluteError / std::abs(exact);
    else std::cout << std::setw(16) << "N/A";
  } else {
    std::cout << std::setw(16) << "N/A" << std::setw(16) << "N/A"
              << std::setw(16) << "N/A";
  }
  std::cout << '\n';
}

void SaveResults(std::ofstream& out, unsigned level, const IntegrationTest& test,
                 const IntegrationResult& result) {
  const auto write = [&](const char* method, const char* domain, double numerical, double exact) {
    out << level << ',' << test.name << ',' << method << ',' << domain << ','
        << std::setprecision(17) << numerical << ',';
    if(std::isfinite(exact)) {
      const double absoluteError = std::abs(numerical - exact);
      out << exact << ',' << absoluteError << ',';
      if(exact != 0.) out << absoluteError / std::abs(exact);
      else out << "NA";
    } else {
      out << "NA,NA,NA";
    }
    out << ',' << result.cutElements << ',' << result.algoimVolumePoints
        << ',' << result.algoimSurfacePoints << '\n';
  };
  if(result.momentComputed) {
    write("moment", "omega", result.omega1, test.exactOmega1);
    write("moment", "gamma", result.gamma, test.exactGamma);
  }
  if(result.algoimComputed) {
    write("algoim", "omega", result.omega1Algoim, test.exactOmega1);
    write("algoim", "gamma", result.gammaAlgoim, test.exactGamma);
  }
}

void PrintResults(const IntegrationTest& test, const IntegrationResult& result) {
  if(result.momentComputed) {
    PrintRow(test.name, "Moment", "Omega+", result.omega1, test.exactOmega1);
    PrintRow(test.name, "Moment", "Gamma", result.gamma, test.exactGamma);
  }
  if(result.algoimComputed) {
    PrintRow(test.name, "Algoim", "Omega+", result.omega1Algoim, test.exactOmega1);
    PrintRow(test.name, "Algoim", "Gamma", result.gammaAlgoim, test.exactGamma);
  }
}

} // namespace

int main(int argc, char** argv) {
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);
  int iproc = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

  DriverOptions options;
  try {
    options = Parse(argc, argv);
  } catch(const std::exception& error) {
    if(iproc == 0) std::cerr << "Argument error: " << error.what() << "\nUse --help for options.\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  const std::string psiName = "Psi";
  for(const GeomElType element : options.elements) {
    const unsigned dim = (element == TET || element == HEX) ? 3u : 2u;
    for(const PsiType kind : options.psis) {
      std::vector<double> center(dim, 0.);
      double radius = 0.25;
      if(element == TRI || element == TET) {
        const double rin = 1. / (dim + std::sqrt(static_cast<double>(dim)));
        center.assign(dim, -0.5 + rin);
        radius = 0.5 * rin;
      }

      PsiConfig cfg;
      cfg.type = kind;
      cfg.center = center;
      cfg.radius = radius;
      cfg.amplitude = 0.1;
      cfg.waveNumber = 2. * std::acos(-1.);
      Psi psi(cfg);
      const bool manufactured = options.curved && kind == PsiType::BALL;
      const ManufacturedBallPsi manufacturedPsi = {center, radius, options.epsilon};
      const auto tests = MakeTests(kind, dim, center, radius);

      MultiLevelMesh mlMsh;
      const std::array<double, 3> min = {{-0.5, -0.5, -0.5}};
      const std::array<double, 3> max = {{ 0.5,  0.5,  0.5}};
      const std::array<double, 3> ballCenter = {{center[0], center[1], dim == 3 ? center[2] : 0.}};

      if(iproc == 0) {
        std::cout<<std::endl;
        std::cout << "\n" << std::string(90, '=') << "\n"
                  << "Element: " << ElementName(element) << " | Psi: " << PsiName(kind)
                  << " | Mesh: " << (manufactured ? "manufactured curved" : (options.curved ? "curved" : "straight"));
        if(manufactured) std::cout << " | epsilon: " << options.epsilon;
        std::cout << "\n" << std::string(90, '=') << '\n';
        std::cout<<std::endl;
        std::cout<<std::endl;
      }

      if(!options.curved) {
        mlMsh.GenerateSingleElementMesh(element, min, max, "fifth");
      } else if(manufactured) {
        // All four element types; F(q)=q+epsilon*psi_ref(q)*e_1.
        mlMsh.GenerateSingleManufacturedBallElementMesh(
          element, min, max, ballCenter, radius, options.epsilon, "fifth");
      } else {
        throw std::runtime_error("Not implemented");
      }

      std::string meshTag = options.curved ? "curved" : "straight";
      if(manufactured) {
        std::ostringstream tag;
        tag << "manufactured_eps" << std::fixed << std::setprecision(3) << options.epsilon;
        meshTag = tag.str();
        std::replace(meshTag.begin(), meshTag.end(), '.', 'p');
      }

      std::ofstream csv;
      if(iproc == 0) {
        const std::string file = "results_" + std::string(ElementName(element)) + "_" + PsiName(kind)
                               + "_" + meshTag + ".csv";
        csv.open(file.c_str(), std::ios::trunc);
        if(!csv) throw std::runtime_error("Cannot create CSV: " + file);
        csv << "level,test,method,domain,numerical,exact,absolute_error,relative_error,cut_elements,algoim_volume_points,algoim_surface_points\n";
      }

      // Initialize a new solution on every level: Integrate() always reads the finest level.
      for(unsigned level = 0; level <= options.levels; ++level) {
        if (iproc == 0) {
            // std::cout<<std::endl;
            std::cout  << std::string(90, '-') << std::endl;
            std::cout << "\nLevel " << level << " (" << mlMsh.GetNumberOfLevels() << " mesh level(s))\n";
            // std::cout  << std::string(90, '-') << std::endl;
            std::cout<<std::endl;
        }
        if(level != 0) {
          // AMR marking must use the same level set used for initialization.
          if(manufactured) FlagFinestMeshLevel(mlMsh, manufacturedPsi);
          else FlagFinestMeshLevel(mlMsh, psi);
          mlMsh.AddAMRMeshLevel(false);
        }
        {
          MultiLevelSolution mlSol(&mlMsh);
          mlSol.AddSolution(psiName.c_str(), LAGRANGE, SECOND, 0, false);
          mlSol.Initialize("All");
          // This also works on refined levels because manufacturedPsi
          // inverts the global quadratic geometry map at each physical DoF.
          if(manufactured) InitLevelSet(mlSol, psiName, manufacturedPsi);
          else InitLevelSet(mlSol, psiName, psi);

          if(options.vtk) {
            const std::string outputdir = "output_" + std::string(ElementName(element));
            VTKWriter writer(&mlSol);
            writer.SetDebugOutput(true);
            writer.Write(outputdir, "biquadratic", std::vector<std::string>{psiName}, 0);
          }

          if(iproc == 0) {
            std::cout<<std::endl;
            PrintHeader();
          }
          for(const auto& test : tests) {
            const IntegrationResult result = Integrate(mlSol, psiName, test.f, options.integration);
            if(iproc == 0) {
              PrintResults(test, result);
              SaveResults(csv, level, test, result);
            }
          }
        }
      }
    }
  }
  return 0;
}

template <class PsiType>
void FlagFinestMeshLevel(MultiLevelMesh& mlMsh, const PsiType& psi) {

  const unsigned level = mlMsh.GetNumberOfLevels() - 1;
  Mesh& msh = *mlMsh.GetLevel(level);

  const unsigned iproc = msh.processor_id();
  const unsigned dim = msh.GetDimension();
  const unsigned xType = 2;
  const unsigned amrIndex = msh.GetAmrIndex();

  auto& xv = msh._topology->_Sol;
  auto* amrFlag = msh._topology->_Sol[amrIndex];

  amrFlag->zero();

  const unsigned offset = msh._elementOffset[iproc];
  const unsigned offsetp1 = msh._elementOffset[iproc + 1];

  std::vector<double> x(dim);

  for(unsigned iel = offset; iel < offsetp1; ++iel) {

    const unsigned nDof = msh.GetElementDofNumber(iel, xType);

    const unsigned xDof0 = msh.GetSolutionDof(0, iel, xType);

    for(unsigned k = 0; k < dim; ++k)
      x[k] = (*xv[k])(xDof0);

    const double psi0 = psi(x);
    const int sign0 = (psi0 > 0.) ? 1 : -1;

    bool signChange = false;

    for(unsigned i = 1; i < nDof; ++i) {

      const unsigned xDof = msh.GetSolutionDof(i, iel, xType);

      for(unsigned k = 0; k < dim; ++k)
        x[k] = (*xv[k])(xDof);

      const double psii = psi(x);
      const int signi = (psii > 0.) ? 1 : -1;

      if(signi != sign0) {
        signChange = true;
        break;
      }
    }

    if(signChange)
      amrFlag->set(iel, 2);
  }

  amrFlag->close();

  for(unsigned iel = offset; iel < offsetp1; ++iel) {

    if((*amrFlag)(iel) == 2) {

      for(unsigned j = 0; j < msh.el->GetElementNearElementSize(iel, 1); ++j) {

        const unsigned jel = msh.el->GetElementNearElement(iel, j);

        if(offset <= jel && jel < offsetp1) {

          if((*amrFlag)(jel) < 0.5)
            amrFlag->set(jel, 1);
        }
        else {
          amrFlag->add(jel, 1);
        }
      }
    }
  }

  amrFlag->close();

  unsigned localRefined = 0;

  for(unsigned iel = offset; iel < offsetp1; ++iel) {

    if((*amrFlag)(iel) > 0.5) {

      amrFlag->set(iel, 1);

      ++localRefined;
    }
  }

  amrFlag->close();

  unsigned globalRefined = 0;

  MPI_Allreduce(
    &localRefined,
    &globalRefined,
    1,
    MPI_UNSIGNED,
    MPI_SUM,
    MPI_COMM_WORLD);

  msh.el->SetRefinedElementNumber(globalRefined);
}

template <class PsiType>
void InitLevelSet(MultiLevelSolution & mlSol, const std::string & name,
                  const PsiType & psi2D) {

  MultiLevelMesh &mlMsh = *mlSol.GetMultilevelMesh();
  const unsigned level = mlMsh.GetNumberOfLevels() - 1u;
  Mesh &msh = *mlMsh.GetLevel(level);
  Solution &sol = *mlSol.GetLevel(level);
  const unsigned iproc = msh.processor_id();
  const unsigned dim = msh.GetDimension();

  unsigned solIndex = mlSol.GetIndex(name.c_str());
  unsigned solType = mlSol.GetSolutionType(name.c_str());

  const unsigned xType = 2u; // coordinate field type

  auto &xv = msh._topology->_Sol;
  std::vector<double> x(dim);

  const unsigned offset = msh._elementOffset[iproc];
  const unsigned offsetp1 = msh._elementOffset[iproc + 1];

  auto &solVec = sol._Sol[solIndex];

  solVec->zero();

  // Loop over local elements and interpolate psi2D at solution DoFs
  for (unsigned iel = offset; iel < offsetp1; ++iel) {

    const unsigned nDof = msh.GetElementDofNumber(iel, solType);

    for (unsigned i = 0; i < nDof; ++i) {

      // Get physical coordinates of the current DoF
      const unsigned xDof = msh.GetSolutionDof(i, iel, xType);
      for (unsigned k = 0; k < dim; ++k) {
        x[k] = (*xv[k])(xDof);
      }

      // Evaluate and assign field value
      const unsigned solDof = msh.GetSolutionDof(i, iel, solType);
      solVec->set(solDof, psi2D(x));
    }
  }

  solVec->close();
}