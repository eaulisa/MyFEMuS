# ex47 to native CutFEM migration

Baseline: `origin/barbi-ex47-levelset-data`, commit `e7f9724713afcb6bab24d5ec388d84b244e14b3b`. Validation date: November 18, 2026.

## Source layout

| Previous location or name | Current location or name |
| --- | --- |
| 28 native headers and TODO in `src/quadrature/cutfem_ls` | Preserved byte-for-byte under `src/quadrature/cutfem_ls/old` |
| `ex47/include` and `ex47/include/prepared` | Native `src/quadrature/cutfem_ls`; the local copies are removed |
| `Ex47::CutFemWeight_LS` | Global `CutFemWeight_LS`, matching the previous native class name |
| `FELevelSetMoments::Ex47Prepared` | `FELevelSetMoments` |
| `LineDomainIntegration` plus its empty derived wrapper | One `FELevelSetMoments::LineIntegration` class |
| Native application-only tangential factors | `applications/MultiPhase/includeLS/TangentialSurfaceFactors.hpp` |
| Primitive tests for the former native implementation | `unittests/testCutFemWeightLSNative/old` |
| Side-by-side comparison harness for the two former implementations | `ex47/old/native_tests` |

There are 23 active native headers, totaling 7,055 lines. They contain the ex47 implementation and its required algebra, fitting, quadrature, and geometry operations. The former native tensor/simplex integrators, nodal classifier, rule-owner wrapper, and unused helpers are archived. The two operations that ex47 previously obtained through the old planar class—coefficient graph qualification and exact univariate validation—are retained directly as required math in the active planar integrator. No active code includes or calls the archived implementation. A compiler dependency audit found zero archived or ex47-local includes in 167 dependency files.

## Interface and application changes

The preferred interface constructs `CutFemLevelSetData` once for a cut element and passes it to `GetWeightsFromLevelSetData`, `GetPhysicalSurfaceWeightsFromLevelSetData`, `GetTangentialSurfaceWeightsFromLevelSetData`, or `GetNormalSurfaceWeightsFromLevelSetData`. Set the GK profile with `SetOuterQuadratureRules` before these calls. The nodal-input methods already required by ex47 remain. The former native physical-surface overload's extra Boolean argument is removed.

| Application | Migration |
| --- | --- |
| ex47 | Includes the native headers; checks independent references and nodal/prepared interface equivalence; timing labels distinguish those two interfaces |
| ex44 | Uses prepared element data and the selected GK profile for physical/reference surface weights |
| ex44_test | Uses prepared data for native moments; Algoim adapters consume unit coefficients and exact simplex Duffy substitutions |
| ex45 | Already uses its own direct Algoim implementation and needs no source adaptation; it builds and runs independently of the archived code |
| ex46 | Uses application-owned per-geometry rule objects and prepared data for tangential, physical, and signed-normal measures |
| ex46_par | Uses the same migrated rule/data interface; its MPI assembly and reductions are preserved |

The shared AMR/color classifier uses the new data class with the previous certification depth of 12. ex47's default depth and integration budgets are unchanged. The launcher setting `SURFACE_EXTENDED_BUDGET` and its `FEMUS_CUTFEM_UNLIMITED` export are removed because they applied to the archived integrator. No obsolete budget interface is added to the replacement.

Two benchmark-only coordinates, `-0.23` and `-0.13`, were corrected to valid unit coordinates `0.23` and `0.13`. The before/after kernel runs use the same corrected inputs. One unrelated existing compiler error required qualifying `isnan` as `std::isnan` in `LinearImplicitSystem.cpp`.

## Build and accuracy

All six targets built: `MultiPhase_ex44`, `MultiPhase_ex44_test`, `MultiPhase_ex45`, `MultiPhase_ex46`, `MultiPhase_ex46_par`, and `MultiPhase_ex47`. The full application build used GCC 14.2, PETSc 3.26.1, Open MPI 5.0.7, METIS 5.1, ParMETIS 4.0.3, Eigen 3.4, Boost 1.83, LAPACKE, HDF5, and MUMPS. SLEPc was disabled. The branch uses a PETSc API unavailable in the initially available PETSc 3.22 package; a matching PETSc was built without changing the solver API. Matching METIS/ParMETIS dependencies were also required for the parallel checks.

| Check | Result |
| --- | --- |
| Four standalone CTest groups | All pass: data, LINE, planar weights, spatial weights |
| Data checks, 256 cases and 64 points | 5,021,366 pass; maximum scaled error `1.4422057e-14`, unchanged |
| 16 fitted/raw/reference-floor accuracy rows | Every printed `prepared_error`, `raw_error`, and `fitting_floor` matches the pre-migration run |
| Raw degree-13 moments | HEX `1.7896944e-17`; TET `2.0305602e-16`; WEDGE `7.1717491e-17` maximum reported errors |
| Random TRI/TET Duffy polynomial checks against FEM bases | 20,000 pass; maximum absolute error `9.71011e-16` |
| ex44_test straight balls, TRI/QUAD/TET/HEX, degree 5 | Maximum absolute native error over volume/surface integrals of `1`, `r^2`, and `r^4`: `3.99680289e-14` |
| ex44_test curved QUAD/HEX, degree 10 | All twelve volume/surface comparisons of `1`, `r^2`, and `r^4` agree within `2.19824159e-14` absolute |
| ex44_test Algoim moment backend | Curved QUAD/HEX smoke runs finish; its order-7 quadrature retains its own finite quadrature error |
| ex44 and ex46 | Symmetric 2D PAR, two uniform levels, no adaptive levels, zero level offsets, one requested time interval; runs finish and exercise cut elements |
| ex45 | Nonsymmetric 2D PAR smoke run finishes |

The supplied finite-precision Gauss tables impose a fitting floor; for example, LINE degree 10 has a reported floor near `1.18e-11`. Raw integration accuracy and fitted-weight accuracy remain separate in the tests. Degree 5 underintegrates some curved `r^2`/`r^4` tests; the QUAD/HEX accuracy check above uses degree 10. ex44's printed area/perimeter diagnostics use a hard-coded radius of 0.25, whereas PAR uses radius 0.002; those printed relative errors are not PAR accuracy criteria.

## MPI check

The fixed-interface 3D PAR case used a generated `2x2x2` HEX mesh, degree 5, tolerance `1e-8`, one solve, `--nsteps 1000`, Young–Laplace initial pressure, and pressure/velocity ghost penalties plus stabilization enabled. Both surface-weight modes completed with one and two MPI processes.

| Surface weights | L2 velocity, 1 rank | L2 velocity, 2 ranks | Largest relative difference across L2, maximum Gauss, and maximum nodal velocity |
| --- | ---: | ---: | ---: |
| 1 | `3.051484113561887e-7` | `3.051484113561903e-7` | `4.70e-13` |
| 3 | `8.458102405006516e-7` | `8.458102405006506e-7` | `1.65e-13` |

These are comparisons of the same spatial/time-discrete problem across MPI process counts, not exact zero-current assertions for the coarse discrete sphere.

## Before/after timing

The complete weight benchmark uses the same harness before and after migration, GCC `-O2 -DNDEBUG`, one pinned CPU, and no concurrent build or simulation. It covers all six elements, degrees 2 and 5, delta/reference/physical/tangential modes, and the three-normal mode in 3D: 54 cases. Constructors are outside the timed region. Every integration clears the final fitted-weight cache and alternates two fields. Each process reports the median of five timings; three before/after process pairs are compared.

All 54 checksums are identical at the printed 20 significant digits. The median ratio `before/after` is **1.0022**, effectively unchanged. The largest individual median slowdown is 9.8%, with overlapping before/after repeat ranges. No case shows a slowdown above 10% with separated repeat ranges. The earlier 133 data/algebra/quadrature kernel measurements have median `before/after = 1.0042`. These measurements show no systematic migration regression; individual short timings retain measurement noise.

Raw paired timings are in [weight_timing.csv](migration_validation/weight_timing.csv) and [kernel_timing.csv](migration_validation/kernel_timing.csv). [summary.json](migration_validation/summary.json), [weight_accuracy.json](migration_validation/weight_accuracy.json), and [application_runs.json](migration_validation/application_runs.json) retain the numerical results and exact application arguments. Application wall times include startup and one or two steps and are not scaling benchmarks.

## Repeating the checks

```sh
cmake -S unittests/testCutFemWeightLSNative -B build-cutfem -DCMAKE_BUILD_TYPE=Release
cmake --build build-cutfem --parallel 2
ctest --test-dir build-cutfem --output-on-failure
cmake --build build --target b64 jsoncpp adept
cmake --build build --target MultiPhase_ex44 MultiPhase_ex44_test MultiPhase_ex45 MultiPhase_ex46 MultiPhase_ex46_par MultiPhase_ex47
```

Configure the main build with the dependencies above. Run application commands from their build directories, preserving their input folders. `application_runs.json` uses `./migration_validation_run` for output and solver-option paths; create that directory and copy the supplied `par-solver.options` there. Use `mpiexec -n 2` for each two-process entry.

Long production transients, the default degree-13/17 application profiles, and large-process MPI scaling were not run. This migration preserves the existing ex47 numerical behavior, including the extreme-small-cut limits described in the archived review; it is not a numerical-method redesign.
