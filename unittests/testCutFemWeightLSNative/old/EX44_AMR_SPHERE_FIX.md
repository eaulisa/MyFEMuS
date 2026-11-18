# ex44_test adaptive sphere fix

Baseline: `2b11e34d8bf8a14990895503e26fcc1658322c8c` on `rifat-cutfem-endpoint-only` (the standalone ex44_test.cpp push).

## Failure and correction

The unchanged `run_ex44_test.sh` reproduced the reported abort at HEX refinement level 2, element 6, with 64 mesh elements. The exact 27 nodal values supplied in the report are preserved in `hex_sphere_boundary_tangencies.cpp`.

At the sphere contact with a cell face, the face trace has zero tangential gradient and a definite Hessian, while the derivative of the original parent polynomial transverse to that face is nonzero. For example, on the face x=-1/4 the sphere polynomial G=1/16-x*x-y*y-z*z reduces to -y*y-z*z, although the full gradient is (1/2,0,0) at the contact. The old classifier treated this isolated face tangency as an unsupported non-transverse surface projection.

`EndpointProjection.hpp` now carries the exact parent derivative restricted to each face and checks that derivative and the face Hessian locally. Only a definite, isolated face tangency with nonzero transverse parent derivative receives this certificate. Fiber discriminants, singular parent gradients, and unsupported saddle contacts retain their explicit rejection. The original level set and p=2 endpoint map are preserved.

After the HEX fix, the same runner exposed a separate TET failure at refinement level 2, element 38. Nearly coincident projected events create endpoint strips whose GK15 nodes can round to the endpoint. `SurfaceEndpointIntegration.hpp` now chooses GK7 or GK3 within the same transformed-integrand quadrature when those nodes are representable. The embedded estimate also includes an adjacent-representable-slice comparison and its inner error. Every strip contribution is integrated and all selective and accumulated error checks remain active. If representability or tolerance cannot be met, the explicit failure remains.

The existing debug output now prints actual coordinates, rather than three copies of Psi, and prints nodal values at double precision.

## Application validation

Built target `MultiPhase_ex44_test` and ran the checked-in `run_ex44_test.sh` unchanged: HEX, TET, QUAD, TRI; ball; straight geometry; levels 0 through 4; degree 8; moment mode; VTK enabled. One MPI rank. Exit status: 0. All 80 reported volume/surface rows are finite and satisfy the default absolute-plus-relative analytical error threshold (1e-8 each).

Final-level errors for the constant integrand:

| Element | Elements at level 4 | Cut elements | Volume/area absolute error | Surface/curve absolute error |
|---|---:|---:|---:|---:|
| HEX | 2024 | 296 | 4.49640325e-15 | 1.22124533e-15 |
| TET | 2339 | 248 | 2.24039537e-15 | 1.58774438e-10 |
| QUAD | 172 | 36 | 3.60822483e-16 | 1.88737914e-14 |
| TRI | 184 | 32 | 2.77555756e-17 | 2.88657986e-15 |

The four CSV files in `ex44_amr_sphere_results/` preserve every level, both integrands (1 and r^4), numerical values, exact references, and separate errors.

## Maintained regressions

- 20/20 standalone native CTests pass.
- `hex_sphere_boundary_tangencies`: 48 cases and 27,459 checks, including the exact element-6 nodal data, all eight octants, monomials through total degree 10, raw delta/reference surface/physical surface metrics, and tolerances 1e-10 and 1e-12. Maximum physical volume error: 1.3894649e-18; physical surface error: 7.5035208e-17.
- `endpoint_coordinate_resolution`: 336 checks covering both strip orientations, powers 2/4/6, GK7/GK3 representability, exact transformation Jacobians, selective masks, and explicit rejection below coordinate resolution.
- AddressSanitizer and UndefinedBehaviorSanitizer pass both new regression executables. Leak detection was disabled because the execution environment does not support LeakSanitizer process inspection.
- The retained LINE and 2D implementations are unchanged; their maintained tests pass.

## Build environment

GCC 14.2.0; C++17; release -O3 -DNDEBUG; OpenMP enabled; AMD EPYC 9V74. Local development dependencies were installed outside the repository. PETSc 3.22.5 needs the compiler-only name mappings `KSPOrthogonalizationSetCGSRefinementType=KSPGMRESSetCGSRefinementType` and `KSP_ORTHOGONALIZATION_CGS_REFINE_IFNEEDED=KSP_GMRES_CGS_REFINE_IFNEEDED` for the pre-existing solver sources. ParMETIS and METIS were built from their upstream sources; no solver or partitioning source changes are part of this fix.

The complete script took 169.34 seconds wall time, including node printing, VTK, and concurrent validation.

## Scope

No geometric subdivision, polynomial refitting, discarded strip contributions, Taylor path, or alternative surface integrator was introduced. The supported non-transverse extension is the certified isolated face tangency described above. Other unsupported projected degeneracies and unresolved numerical error estimates still raise errors.
