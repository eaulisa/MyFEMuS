# Unit-coordinate correction

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
- Public native cases, curved-HEX regressions, warning contracts, and timing comparisons are being rerun before the final report.

The last-ULP WEDGE discriminant comparison uses 4e-16 relative tolerance after conversion (previously 3e-16); the earlier observed error was 3.0385e-16, with retained-slope density error 1.4626e-16. Public integration tolerances are unchanged.
