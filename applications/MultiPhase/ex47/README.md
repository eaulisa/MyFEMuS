# ex47: prepared element level-set data

`include/CutFemLevelSetData.hpp` is a local, header-only prototype for LINE3, QUAD9, TRI7, HEX27, TET15, and WEDGE21. All checks are functions called by `ex47.cpp` main. The class needs only standard C++14 and the FEMuS element enum; the test main additionally uses the actual FEMuS bases and Boost multiprecision/Gauss-Kronrod headers. Existing CutFEM integration is unchanged.

Build independently of MPI/PETSc:

```sh
cmake -S applications/MultiPhase/ex47 -B build-ex47 -DCMAKE_BUILD_TYPE=Release
cmake --build build-ex47 --parallel 2
./build-ex47/MultiPhase_ex47
./build-ex47/MultiPhase_ex47 --cases 256 --points 64 --benchmark
./build-ex47/MultiPhase_ex47 --benchmark-quadrature
./build-ex47/MultiPhase_ex47 --benchmark-surface
```

If Boost is outside the compiler's search path, add `-DEX47_BOOST_INCLUDE_DIR=/path/to/boost`. The normal FEMuS build also exposes `MultiPhase_ex47`. `run_ex47.sh` offers the standalone build and accepts the main's arguments; `EX47_BUILD_DIR`, `EX47_JOBS`, `EX47_BOOST_INCLUDE_DIR`, and `CMAKE` override its defaults.

For sanitizers, configure another build with `-DEX47_SANITIZE=ON -DCMAKE_BUILD_TYPE=Debug` and run its main. The main defaults to 64 random cases and 32 points per case; `--benchmark` additionally prints six arithmetic timing/reference-error rows. `--benchmark-quadrature` benchmarks prepared quadrature kernels; `--benchmark-surface` compares each requested measure with computing all measures and compares slicing with/without original-gradient context. Configure another Release build with `-DEX47_TRACK_ALLOCATIONS=ON` to assert zero heap allocations in both warmed benchmark paths.

## Core interface

```cpp
CutFemLevelSetData parent(HEX, nodalValues, 1);
parent.SetContext(meshLevel, elementId, nodalVersion);
int phase = parent.GetPhase();
const std::vector<long double>& power = parent.GetPowerCoefficients();
const std::vector<long double>& bernstein = parent.GetBernsteinCoefficients();
CutFemLevelSetData quad;
CutFemLevelSetData line;
parent.Restrict(2, 0.3L, quad);
quad.Restrict(1, -0.2L, line);
long double value = line.Evaluate({0.1L});
```

| Operation | Interface |
| --- | --- |
| Replace nodal state | `Reset(geometry, nodalValues, maximumDepth)` |
| Change classification depth | `SetMaximumDepth(depth)`; default 1, allowed 0 through 32 |
| Classification | `GetPhase()`, `GetClassification()`, `IsCutCandidate()` |
| Preparation state | `HasPowerCoefficients()`, `HasBernsteinCoefficients()` |
| Polynomial data | `GetPowerCoefficients()`, `GetNormalizedPowerCoefficients()`, `GetNormalizationScale()` |
| Bernstein data | `GetBernsteinCoefficients()`, `GetBernsteinIndices()`, `GetBernsteinOrders()`, `GetBounds()` |
| Evaluation | `Evaluate(point)`, `EvaluateBernstein(point)` |
| Scalar line evaluation | `EvaluateLine(coordinate)` |
| Derivatives | `GetGradient(point)`, `GetDerivativeCoefficients(axis)`, `EvaluateDerivative(axis, point)`, `EvaluateSecondDerivative(first, second, point)` |
| Graph certification | `GetDerivativeBounds(axis)`, `GetGraphDirection()`; direction is an active axis or -1 if none is certified |
| Slice metadata | `GetDimension()`, `GetRootDimension()`, `GetGeometry()`, `IsDegenerate()`, `GetRevision()` |
| Coordinate maps | `MapToParent(point)`, `MapToRoot(point)`, `GetParentOrigin()`, `GetParentMap()`, `GetParentMeasureScale()` |
| Original gradient components | `GetRootGradient(point)` |
| Optional original-gradient context | `Restrict(axis, coordinate, child, retainRootGradient)`, `HasRootGradientContext()` |
| Allocation-free gradients | `GetGradient(point, std::array<long double,3>&)`, `GetRootGradient(point, std::array<long double,3>&)` |
| External geometry polynomials | `RestrictComponents(axis, coordinate, degree, components, output)` |
| One-dimensional roots | `GetLineRoots()`; repeated roots appear once, endpoint roots are retained |
| Allocation-free roots | `GetLineRoots(std::array<long double,2>&)` returns the root count |
| Shared regular-root data | `GetLineSurfaceData(coordinate, rootGradient)` returns `1/abs(dG/dt)` and all original reference-gradient components |
| Delta only | `GetLineDeltaFactor(coordinate)` |
| Reference surface only | `GetLineReferenceSurfaceFactor(coordinate)` |
| Physical surface only | `GetLinePhysicalSurfaceFactor(coordinate, cofactorGradient)` |
| Inverse-metric tangential only | `GetLineTangentialSurfaceFactor(coordinate, cofactorGradient)` |
| One normal component | `GetLineNormalSurfaceFactor(coordinate, cofactorGradient, component)` |
| Shared normal family | `GetLineNormalSurfaceFactors(coordinate, cofactorGradient, normalFactors)` |
| Prepared cofactor matrix | Physical/tangential/normal methods also have `FromCofactor` variants taking a prepared matrix |
| Explicit all-measures batch | `GetLineAllSurfaceFactors(...)`, `GetLineAllSurfaceFactorsFromCofactor(...)` |
| Ownership metadata | `SetContext(level, elementId, nodalVersion)` and the corresponding getters |
| Preparation counters | `GetWorkCounts()` returns nodal expansions, power-to-Bernstein conversions, local derivatives, and restrictions |

## Coordinates and coefficient layout

Nodal input follows the actual FEMuS enriched bases and node ordering. Public evaluation coordinates are the FEMuS reference coordinates: LINE/QUAD/HEX use `[-1,1]`, TRI/TET use the unit simplex, and WEDGE uses the unit triangle times `[-1,1]`.

Power coefficients use three slots per active coordinate. In 3D the index is `9*i+3*j+k`; in 2D it is `3*i+j`; in 1D it is `i`. These are native-reference coefficients, with structural zeros retained. `GetNormalizedPowerCoefficients()` divides them by the reported positive power-of-two normalization scale. Slices inherit the original scale.

Bernstein indices are explicit. Tensor coordinates and the wedge height are mapped to `[0,1]` for Bernstein evaluation. Simplex indices store the free barycentric exponents; the zeroth exponent is the degree minus their sum. Root orders are 2 per tensor coordinate, 3 for TRI7, 4 for TET15, and `(3,3,2)` for WEDGE21. Exact simplex restriction retains a degree-elevated control net when needed: a tetrahedral TRI slice can have order 4 even though its power polynomial is lower degree. No coefficient is discarded to reduce that order.

`Restrict()` fixes an active reference coordinate and fills an existing child. HEX produces QUAD, TET produces TRI, and WEDGE produces TRI when fixing height or QUAD when fixing a triangle coordinate. QUAD/TRI produce LINE. The remaining coordinates are rescaled to the child's standard reference domain. Endpoint slices with zero measure remain valid data objects marked degenerate.

`GetParentMeasureScale()` is the reference-coordinate Jacobian for the unfixed variables, not a physical surface metric. `GetRootGradient()` retains all original reference derivative components, including derivatives along fixed coordinates. Physical geometry and `cof(J)*grad(G)` remain separate; `RestrictComponents()` preserves every supplied component and supports coordinate degrees up to 32.

## Classification and ownership

Nodal sign changes or zero nodal values are handled before coefficient construction. Otherwise the class builds and retains power/Bernstein coefficients. Strict sign certification uses a scaled long-double roundoff margin; slices inherit their parent's margin and include their own contraction error allowance. Derivative bounds include the corresponding control-net uncertainty. Ambiguous cases can receive the requested number of Bernstein bisections. Phase +1/-1 means certified positive/negative; phase 0 includes crossings, contacts, and uncertified elements. A single zero Bernstein coefficient is not treated as proof of a crossing.

The class owns its nodal snapshot. Coefficients, normalization, derivatives, and Hessians are prepared once per state. Changing depth reuses polynomial data. A reset with identical nodal data preserves preparation; changed nodal data invalidates it. Context metadata does not alter the polynomial. Slices own their restricted state and remain valid if the parent is reset or destroyed. Reusing a child retains its coefficient/derivative storage and invalidates slice-dependent caches.

The three-argument `Restrict()` inherits the parent's original-gradient availability; a nodal parent starts with it available. Passing `false` as the fourth argument skips preparing and contracting every original-gradient polynomial. Descendants inherit this choice. Values, roots, classification, maps, local derivatives, and graph certification remain available. Original-gradient/reference-surface/matrix-input requests on such a slice throw a clear error; restore that context by restricting again from an ancestor that retained it. `HasRootGradientContext()` makes the dependency explicit.

No quadrature degree, GK rule, tolerance, moment matrix, or final quadrature weight array is owned by this class. It is mutable preparation/workspace data, so concurrent mutation of the same instance requires caller synchronization.

## Quadrature reuse and geometry

Reference index tables, Bernstein restriction plans, and simplex multipliers are shared immutable data initialized once. Restrictions use strided Horner contraction, retained child storage, sparse coordinate-map composition, and fixed-size subdivision buffers. Bounds and graph directions are cached. Derivative bounds come directly from the Bernstein control net, without reconstructing derivative polynomials. The vector-returning interfaces remain available; use the fixed-output overloads in quadrature loops.

Select the requested measure and its data dependencies before the quadrature loops, then call its dedicated evaluator at regular roots. There is no per-point dispatch to other measures. Delta/reference integration never prepares or restricts geometry polynomials. Physical modes require geometry: retain a constant cofactor matrix for affine geometry and use the selected `FromCofactor` method, or prepare `h = cof(J)*grad_ref(G)` once, retain its components separately, and restrict them alongside `G` with `RestrictComponents()`. Geometry-dependent preparation must be refreshed when either geometry or level-set data change.

For a regular root in the final LINE chart, let `a = 1/abs(dG/dt)`, `g = grad_ref(G)`, and `m = abs(h)`. Use these individual integration densities:

| Requested measure | Density | Retain original gradients during restriction? |
| --- | --- | --- |
| Delta | `a` | No |
| Reference surface | `a*abs(g)` | Yes |
| Physical surface, prepared `h` | `a*m` | No |
| Inverse-metric tangential, prepared `h` | `a/m` | No |
| Normal component/family, prepared `h` | `a*h[i]/m` | No |
| Physical/tangential/normal, cofactor matrix | Corresponding density above | Yes, to form `h` |
| Explicit all-measures batch | All four scalars and all normal components | Yes |

For example, a delta-only traversal can omit original-gradient context at its first restriction and inherit that choice at every lower level:

```cpp
parent.Restrict(2, s, quad, false);
quad.Restrict(1, t, line);
long double density = line.GetLineDeltaFactor(root);
```

Prepared-`h` physical/tangential/normal traversals can use the same restriction choice. Reference-surface and matrix-input traversals retain original gradients. Computing one normal component still needs every component of `h` to obtain its norm; request `GetLineNormalSurfaceFactors()` when the whole normal family is needed so that the derivative and norm are shared.

Only `GetLineAllSurfaceFactors()` and its matrix variant compute all measures. Their scalar array is `[a, a*abs(g), a*m, a/m]`, and their normal array is `a*h/m`. These explicit batch names replace the prototype's earlier `GetLineSurfaceFactors` names.

Multiply these densities by the outer chart measure scales and the selected quadrature weights. No extra physical determinant belongs in these surface factors. They are integrand densities, not final fitted CutFEM Gauss weights: moment integration and the reusable moment-to-weight operator remain the integration layer's responsibility. Norms and inverse norms are evaluated directly, without polynomial approximation. Zero line derivatives require endpoint handling and zero physical metrics are rejected.

`PrepareSurfaceTestPolynomials()`, `TestFastQuadratureInterfaces()`, `TestSelectiveSurfacePreparation()`, and `TestSelectiveSurfaceScaling()` in the main cover every element family, affine/non-affine geometry, bi/triquadratic QUAD/HEX geometry, omitted dependencies, and small cuts over widely scaled inputs. The production integration interfaces have not yet been migrated.

The agreed tasks are tracked in [the shared TODO](../../../src/quadrature/cutfem_ls/TODO.md), with [local status](TODO.md), [validation results](VALIDATION.md), and [kernel timings](PERFORMANCE.md).
