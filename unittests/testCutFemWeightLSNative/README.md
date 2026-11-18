# Native CutFEM tests

The active tests execute the library through `applications/MultiPhase/ex47`.
They cover all six enriched elements, prepared data, fitted weights, independent
analytic references, and degree-13 spatial moments.

Build this directory independently or enable `BUILD_TESTING` in the main build,
then run `ctest --output-on-failure`.

The previous primitive-interface suite is retained in `old/` for reference.
It belongs to the implementation archived under `src/quadrature/cutfem_ls/old/`
and is excluded from the active build.
