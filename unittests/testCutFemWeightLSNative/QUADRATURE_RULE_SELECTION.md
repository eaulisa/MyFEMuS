# Configurable outer Gauss-Kronrod rules

This change extends `8a2fa3fef5f4d8a7a9e327ffddad5ef9e1242535` on `rifat-cutfem-endpoint-only`. All four defaults remain G7/K15.

| Setting | Interval |
| --- | --- |
| `regular2D` / `--gk-regular-2d` | Regular outer intervals in QUAD/TRI, including slices inside HEX/TET |
| `tangent2D` / `--gk-tangent-2d` | Endpoint-mapped intervals in QUAD/TRI, including slices inside HEX/TET |
| `regular3D` / `--gk-regular-3d` | Regular outer intervals in HEX/TET |
| `tangent3D` / `--gk-tangent-3d` | Endpoint-mapped outer intervals in HEX/TET |

`src/quadrature/cutfem_ls/GaussKronrodRules.hpp` stores immutable shared descriptors for every G_n/K_(2n+1) pair with n=1,...,15: 1/3, 2/5, 3/7, 4/9, 5/11, 6/13, 7/15, 8/17, 9/19, 10/21, 11/23, 12/25, 13/27, 14/29, 15/31. The descriptors reuse Boost's persistent long-double nodes and embedded weights. Rules are initialized once; integration panels perform direct table access. Even and odd Gauss orders use their respective embedded-node parity.

## C++ interface

The four `OuterQuadratureRules` fields store the **Kronrod point count**, an odd integer from 3 through 31. `GaussKronrodRules::Parse` converts either slash or dash pair notation. Invalid pairs and unsupported orders are rejected before changing a weight object's configuration.

```cpp
FELevelSetMoments::OuterQuadratureRules rules;
rules.regular2D = FELevelSetMoments::GaussKronrodRules::Parse("3/7");
rules.tangent2D = FELevelSetMoments::GaussKronrodRules::Parse("10/21");
rules.regular3D = FELevelSetMoments::GaussKronrodRules::Parse("7/15");
rules.tangent3D = FELevelSetMoments::GaussKronrodRules::Parse("15/31");
cfw.GetWeightsFromNodalValues(psi, weightOmega, weightDelta, 1e-8, 1e-8, false, rules);
cfw.GetPhysicalSurfaceWeightsFromNodalGeometry(psi, coordinates, weightOmega, weightGamma, 1e-8, 1e-8, rules);
```

Both methods accept the configuration as their final argument. These per-call overrides do not change the object's persistent defaults. Alternatively, call `cfw.SetOuterQuadratureRules(rules)` and use the existing method signatures. `GetOuterQuadratureRules()` returns those persistent settings. The nodal-weight cache includes all four rule choices, so the same nodal and geometry data cannot return weights computed with a different configuration.

The native QUAD, TRI, HEX, TET and unit-element entry points also accept a final configuration argument with the same defaults. Configuration is passed explicitly through the hierarchy; no global mutable integration policy is used by the library. The exact LINE3 integration is unchanged.

## Applications and run scripts

Both `ex44` and `ex44_test` forward the four choices into the actual weight-method calls. For example:

```sh
./MultiPhase_ex44_test --elements quad,tri,hex,tet --psi ball --levels 1 --degree 8 --vtk 0 --gk-regular-2d 7/15 --gk-tangent-2d 10/21 --gk-regular-3d 7/15 --gk-tangent-3d 15/31
```

The same four flags are accepted by `MultiPhase_ex44`. `--help` lists their defaults and supported pair notation. `run_ex44.sh` and `run_ex44_test.sh` expose `GK_REGULAR_2D`, `GK_TANGENT_2D`, `GK_REGULAR_3D`, `GK_TANGENT_3D` in their configuration sections. Each defaults to `7/15` and can also be supplied through the environment:

```sh
GK_TANGENT_2D=10/21 GK_TANGENT_3D=15/31 ./run_ex44_test.sh
```

These settings control the adaptive outer GK integration, separately from the moment-fitting degree and background FEM Gauss quadrature.

## Endpoint maps and precision failsafe

The interval strategy remains unchanged: no vertical tangent uses a regular panel; one tangent uses a full-length map with delta=L; two tangents use two maps with delta=L/2. Higher endpoint exponents and selective moment refinement are preserved. Ordinary accuracy failure subdivides using the chosen rule and retains the same endpoint map.

When the nearest mapped point rounds onto its endpoint, the failsafe starts below the **selected** rule and repeatedly halves its Gauss order using integer division. Thus 15/31 -> 7/15 -> 3/7 -> 1/3, while 10/21 -> 5/11 -> 2/5 -> 1/3. The default path remains 7/15 -> 3/7 -> 1/3. It never raises the user's selected order. Rules whose nearest point still collapses are skipped. A resolvable candidate must satisfy the embedded estimate, propagated inner errors and adjacent-coordinate sensitivity check for every remaining moment. Failure at the available precision is reported explicitly; it does not discard a contribution or repeatedly subdivide a collapsed interval.

An inner precision failure also reaches the enclosing panel as `OuterQuadratureResolutionError`. The parent retries its own selected regular or tangent rule downward, discarding the incomplete panel and checking the complete replacement panel before accepting its moments. This matters when a G15/K31 outer node creates a 2D slice whose own endpoint interval is too small even for G1/K3: reducing the outer rule changes the slice locations. Budget, topology and ordinary convergence errors do not enter this retry. Work spent on failed inner slices still counts against the parent evaluation budget. Regular retries probe adjacent coordinates directly; mapped retries retain their map and Jacobian.

Low rules can require substantially more subdivisions at tight tolerances, especially in nested 3D integration. Evaluation and depth limits are unchanged. Selecting G1/K3 does not guarantee that an arbitrary high-degree case will finish within those limits.

## Validation

The reproducible native test is `cutfem_quadrature_rule_selection`. It verifies cached tables and embedded polynomial exactness through degree ten, all supported primary rules on regular/mapped intervals in both dimensions, independently mixed HEX/TET slice settings, analytic QUAD/TRI/HEX/TET volume and delta references, per-call and persistent cache behavior, unchanged exact LINE3 moments, and downward coordinate-collapse transitions. Existing endpoint, geometry and refinement tests remain enabled.

Executed results are stored in `rule_selection_results`:

- All 26 native CTests pass, including 18,522 configurable-rule checks. The latter also pass AddressSanitizer and UndefinedBehaviorSanitizer; LeakSanitizer is disabled because this environment uses ptrace.
- With all four defaults, 180/180 degree-ten cases complete and meet the independent-reference criterion. All 63,360 moment values match the `8a2fa3f` baseline exactly when parsed as decimal numbers, and root, slice and outer sample counts are unchanged. Maximum absolute errors against the independent references are 6.39e-16 for H and 8.89e-15 for the second moment family, covering raw delta, reference surface and curved physical surface metrics. Timing fields are diagnostic; no performance conclusion is drawn from these runs.
- Both full FEMuS targets, `MultiPhase_ex44` and `MultiPhase_ex44_test`, compile. Six application parser checks cover help, intermediate/even-order settings, dash notation, invalid pairs and missing arguments. Both run scripts pass `bash -n`.
- The test preserves the actual curved HEX27 at refinement level two which exposed the nested precision failure, using all 27 level-set values and geometry nodes. Default, G15/K31, mixed G10/K21/G15/K31, and independent inner/outer choices are checked against analytic volume and physical surface references; all degree-eight moments are also compared with the default result. Separate synthetic cases verify parent retries on regular and mapped intervals, prevent accumulation of a failed panel, and ensure unrelated errors propagate.

The checked-in `run_ex44_test.sh` completes for all four element types through refinement level four at degree eight, both with defaults and with `GK_REGULAR_2D=10/21 GK_TANGENT_2D=15/31 GK_REGULAR_3D=10/21 GK_TANGENT_3D=15/31`:

| Full curved-ball run | Reference comparisons | Maximum volume error | Maximum surface error |
| --- | ---: | ---: | ---: |
| Defaults | 80 | 2.15e-11 | 1.27e-8 |
| Mixed regular/tangent rules | 80 | 5.00e-14 | 1.27e-8 |

Eleven completed application runs produce 280 successful reference comparisons in total. Additional curved QUAD/TRI cases exercise 1/3, 2/5, 3/7, 7/15, 10/21 and 15/31 at degree eight; straight HEX/TET cases exercise mixed, even-Gauss and highest-order choices at degree four. An all-G1/K3 3D ball run at degree four reaches the unchanged evaluation budget and is recorded as a known limit, separately from the successful runs. G1/K3 still passes the native analytic 3D tests at their stated looser tolerance. No rule is raised automatically to conceal budget exhaustion.

To repeat the default comparison, retain the benchmark executable built from `8a2fa3f` and run:

```sh
python3 unittests/testCutFemWeightLSNative/run_outer_interval_benchmark.py --baseline /path/to/8a2fa3f-benchmark --baseline-sha 8a2fa3fef5f4d8a7a9e327ffddad5ef9e1242535 --candidate build-cutfem/cutfem_outer_interval_benchmark --output unittests/testCutFemWeightLSNative/rule_selection_results/default_comparison
```

`moment_records.jsonl.gz` retains every compared H/D moment; `comparison.csv` and `default_compatibility.json` summarize accuracy and unchanged defaults. Application commands and reference rows are recorded separately under `applications`. The earlier `unified_outer_results` and `UNIFIED_OUTER_QUADRATURE.md` describe the accepted fixed-G7/K15 implementation and remain historical evidence for that revision.
