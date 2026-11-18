# Native level-set quadrature

`CutFemWeight_LS.hpp` is the global native weight class; `CutFemLevelSetData.hpp` prepares the enriched nodal polynomial on unit coordinates. The implementation was promoted from ex47 at baseline `e7f9724713afcb6bab24d5ec388d84b244e14b3b`.

The active headers contain the ex47 data, graph/endpoint integration, fitting, quadrature rules, and their required polynomial and geometry operations. `FELevelSetMoments::LineIntegration` is the single LINE kernel class. The former empty derived wrapper is removed. The coefficient graph qualification and exact univariate validation used by ex47 are retained directly in the active planar integrator.

The 28 former native headers and their TODO are preserved byte-for-byte in `old/`. That directory is absent from the include paths and build sources. Application-only tangential factors belong to `applications/MultiPhase/includeLS/TangentialSurfaceFactors.hpp`.

For the public interface and migration validation, see `applications/MultiPhase/ex47/INTEGRATION.md` and `MIGRATION.md`.
