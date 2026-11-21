# Native CutFEM tests

The active tests execute the library through `applications/MultiPhase/ex47`.
They cover all six enriched elements, prepared data, fitted weights, independent
analytic references, and degree-13 spatial moments.

Build this directory independently or enable `BUILD_TESTING` in the main build,
then run `ctest --output-on-failure`.

The `cutfem_hex27_elements` test replays the 28 reported HEX27 cells through the combined physical-weight API at order 5 and tolerance `1e-8`. Each fixture record contains the element number, its 27 original rows of `psi x y z`, and 56 rows of independent reference moments in the order `volume tangential nx ny nz`. The references use a separate graph integrator with 64-point Gauss–Legendre quadrature. Unit-cube monomials are ordered by increasing total degree, then decreasing x and y exponents. The checks use absolute error for volume, relative error for the tangential weight, and error divided by `1 + abs(reference)` for signed normals.

The previous primitive-interface suite is retained in `old/` for reference.
It belongs to the implementation archived under `src/quadrature/cutfem_ls/old/`
and is excluded from the active build.
