# Unit-coordinate correction

Historical record from before promotion into the native library; see [MIGRATION.md](MIGRATION.md) for the current layout and validation.

Recovered from the interrupted working session on top of commit `4314f8d46`. Production files under `src/` remain unchanged.

| Component | Coordinate contract | Implementation |
|---|---|---|
| Nodal inputs | Existing FEMuS ordering | Direct unit Lagrange expansion; no node-order change |
| Tensor polynomials | LINE, QUAD, HEX on `[0,1]` | Unit power coefficients and derivatives |
| Simplex polynomials | Standard TRI/TET simplex | Existing enriched basis expansion |
| WEDGE polynomials | Unit triangle × `[0,1]` | Unit height expansion |
| Bernstein data | Same chart as power coefficients | Removed redundant coordinate transforms |
| Restrictions | Unit children; simplex width `1-s` | Corrected scalar, component, Bernstein, and ancestor maps |
| Roots and intervals | Unit coordinates | Updated endpoint pinning; boundary tags remain -1/+1 |
| Reflection | `u -> 1-u` | Exact coefficient tests; removed native odd-moment symmetry |
| Projections | Current unit-slice polynomial | Consistent local topology, discriminants, and endpoint slopes |
| Sensitive algebra | Retained nodal/parent context | Compensated fallbacks preserve tiny nodal information |
| Surface measures | Unit derivatives and geometry cofactors | Existing production density factors applied at the public boundary |
| Quadrature setup | Fixed rule objects | Constructor/shared caches; no rule reinitialization per element |
| Inner quadrature | Direct local fiber/metric evaluation | No LINE construction or per-node metric restriction |
| Primitive test adapters | Explicit legacy input/output conversion | Conversion occurs at the test boundary only |

## Recovery validation

Freshly rerun after reconstruction:

- 5,901,177 data checks passed; maximum scaled error 2.3647e-14.
- 32,513 fitted-weight/moment checks passed across all six elements, including degree-13 HEX/TET/WEDGE moments.
- All 386 public native cases pass 36,860 assertions per backend; all repeated timing calls also pass.
- All 32 CTest commands pass for both backends at matched detection depth 12; the requested default depth 1 retains its conservative classification difference.
- The curved HEX27 regression passes 970 checks across 48 cases; warning contracts pass 316 checks.
- Graph and degree-13 surface benchmarks pass 59,520 paired moment comparisons.
- Address/undefined-behavior sanitizers pass 155,239 main checks with leak detection disabled.
- All 50 native sources were compared over 385 extended requests: prepared passes 383, original 382; only the requested default-depth detection difference is prepared-only.
- Retained-element benchmarks cover 58 element/degree/measure rows. Full accuracy and timing results, including slower cases, are recorded separately in `native_tests/results/unit_coordinates/`.

The last-ULP WEDGE discriminant comparison uses 4e-16 relative tolerance after conversion (previously 3e-16); the earlier observed error was 3.0385e-16, with retained-slope density error 1.4626e-16. Public integration tolerances are unchanged.

Raw degree-13 scaled moment errors are 1.790e-17 for HEX, 2.031e-16 for TET, and 7.172e-17 for WEDGE. The smallest degree-10 LINE cuts retain a fitted-weight error floor of 1.185e-11 despite raw moment error of 2.587e-19; raw integration accuracy and fitted-weight conditioning are reported separately.
