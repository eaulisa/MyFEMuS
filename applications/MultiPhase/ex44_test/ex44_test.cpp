#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "Fem.hpp"
#include "GeomElTypeEnum.hpp"
#include "VTKWriter.hpp"


using namespace femus;
#include "../includeLS/Psi.hpp"
#include "../includeLS/Psi.hpp"
#include "IntegrateCut.hpp"

template <class PsiType>
void FlagFinestMeshLevel(MultiLevelMesh& mlMsh, const PsiType& psi);

template <class PsiType>
void InitLevelSet(MultiLevelSolution & mlSol, const std::string & name, const PsiType & psi2D);

int main(int argc, char **argv) {

    FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int iproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

    std::vector<GeomElType> el_type_list = {TET};
    std::vector<PsiType> psi_type_list = {PsiType::BALL};
    std::vector<std::vector<Psi>> psi_list(el_type_list.size());

    std::string psiName = "Psi";

    for (unsigned iel_type = 0; iel_type < el_type_list.size(); iel_type++) {

            unsigned eltype = el_type_list[iel_type];

            unsigned dim = (eltype == 0 || eltype == 1) ? 3 : 2;

            std::vector<double> psi_center(dim, 0.);
            double psi_radius = 0.25;

            if(eltype == TRI || eltype == TET) {
                const double L = 1.;
                const double rin = L / (dim + std::sqrt(double(dim)));
                psi_center.assign(dim, -0.5 + rin);
                psi_radius = 0.5 * rin;
            }

            for (unsigned ipsi_type = 0; ipsi_type < psi_type_list.size(); ipsi_type++) {

                PsiConfig psi_cfg;

                psi_cfg.type = psi_type_list[ipsi_type];
                psi_cfg.center = psi_center;
                psi_cfg.radius = psi_radius;

                Psi psi(psi_cfg);

                MultiLevelMesh mlMsh;
                femus::Mesh* msh = mlMsh.GenerateSingleElementMesh(el_type_list[iel_type],{{-0.5, -0.5, -0.5}},{{ 0.5,  0.5,  0.5}},"fifth");

                unsigned adaptive_levels = 4;

                for (unsigned k = 0; k < adaptive_levels; ++k) {
                    FlagFinestMeshLevel(mlMsh, psi);
                    mlMsh.AddAMRMeshLevel(false);
                } 

                MultiLevelSolution mlSol(&mlMsh);

                mlSol.AddSolution(psiName.c_str(), LAGRANGE, SECOND, 0, false);
                mlSol.Initialize("All");

                InitLevelSet(mlSol, psiName, psi);

                const std::string outputdir = "output_"+std::to_string(eltype);

                VTKWriter vtkIO(&mlSol);
                vtkIO.SetDebugOutput(true);

                vtkIO.Write(outputdir, "biquadratic", std::vector<std::string>{"Psi"}, 0);

                
                const double pi = M_PI;
                const double r = psi_radius;
                const unsigned dim = psi_center.size();

                std::vector<IntegrationTest> tests = {

                {
                    "f(x) = 1",
                    [](const std::vector<double>& x) -> double { return 1.; },
                    (dim == 2) ? pi * r * r : 4. * pi * std::pow(r, 3) / 3.,
                    (dim == 2) ? 2. * pi * r : 4. * pi * r * r
                },

                {
                    "f(x) = r^4",
                    [psi_center](const std::vector<double>& x) -> double {
                    double r2 = 0.;
                    for(unsigned d = 0; d < x.size(); ++d) r2 += (x[d] - psi_center[d]) * (x[d] - psi_center[d]);
                    return r2 * r2;
                    },
                    (dim == 2) ? pi * std::pow(r, 6) / 3. : 4. * pi * std::pow(r, 7) / 7.,
                    (dim == 2) ? 2. * pi * std::pow(r, 5) : 4. * pi * std::pow(r, 6)
                }

                };

                for(const auto& test : tests) RunIntegrationTest(mlSol, psiName, test);


            }
    }


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