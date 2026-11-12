# Native CutFEM endpoint integration

This branch starts at `8bad259c5a51060e4c2674db49a7ae34b75a0367`. `rifat` remains the benchmark at `0ff0669494e603be248d02c44c0a2e9b81f72b4d`.

The hierarchy is HEX -> QUAD -> LINE and TET -> full TRI -> LINE. General curves and surfaces use actual-integrand endpoint quadrature. Regular and endpoint-mapped outer intervals have independent 2D and 3D Gauss-Kronrod choices, defaulting to G7/K15, followed by selective subdivision. Every G_n/K_(2n+1) pair for n=1,...,15 is supported. A single tangent uses the entire interval; two tangents use two half maps. Coordinate-collapse emergencies halve the selected Gauss order down to G1/K3. The exact quadratic LINE3 algorithm is unchanged. Taylor, Full GK dispatch, and geometric subdivision have been removed from the native implementation. The planar and independently used old CutFEM components are separate.

See `ENDPOINT_ONLY_MATHEMATICS.md` for the precise definition of p, the transverse-contact assumptions, TET15 enrichment, and limitations. See `ENDPOINT_ONLY_VALIDATION.md` for the executed results and remaining requirements. Historical files in `review_results` describe earlier commits and are not acceptance evidence for this branch.

## Build and run

The standalone CMake target needs a C++14 compiler, Boost headers, and Eigen headers. Set `BOOST_INCLUDEDIR` and `EIGEN3_INCLUDE_DIR` if needed.

```sh
cmake -S unittests/testCutFemWeightLSNative -B build-cutfem -DCMAKE_BUILD_TYPE=Release
cmake --build build-cutfem --parallel 2
ctest --test-dir build-cutfem --output-on-failure
```

The wider standalone suite uses the actual FEMuS FE and Gauss sources:

```sh
python3 unittests/testCutFemWeightLSNative/build.py --jobs 2
python3 unittests/testCutFemWeightLSNative/run_focused.py --extended
python3 unittests/testCutFemWeightLSNative/run.py --selection original --jobs 2
python3 unittests/testCutFemWeightLSNative/run_outer_interval_benchmark.py --baseline /path/to/baseline-benchmark --candidate build-cutfem/cutfem_outer_interval_benchmark --output unittests/testCutFemWeightLSNative/unified_outer_results
```

All moment tests are limited to total degree ten. `endpoint_surface_regressions` preserves planar, face, scaling, sphere, geometry, and metric cases formerly exercised by the subdivision test. `tri_endpoint_regressions` preserves the former closed-curve and scale cases using full-triangle endpoint integration. `endpoint_exponents` tests analytical transformation models and exact projected constraints. `endpoint_nested_contacts` probes more difficult projected contacts and reports unsupported integrations explicitly; its results must not be interpreted as convergence when it reports UNSUPPORTED.

`endpoint_comparison.cpp` is a common numerical/timing driver for this branch, the accepted start, and `rifat`. Define `FEMUS_ENDPOINT_ONLY` only when compiling this branch, because the obsolete Taylor-selection argument was removed from its public wrappers. Numerical exceptions and timing failures must be reported separately from successful comparisons.

See `UNIFIED_OUTER_QUADRATURE.md` for the f5f81a8 baseline comparison, full-interval mapping rules, coordinate-resolution tests, and final measured results.

See `QUADRATURE_RULE_SELECTION.md` for the current method arguments, application options, stored rules, downward failsafe, and validation against `8a2fa3f` with the default settings.
