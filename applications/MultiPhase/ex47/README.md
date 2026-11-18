# ex47: prepared element level-set data

`include/CutFemLevelSetData.hpp` is a local, header-only prototype for LINE3, QUAD9, TRI7, HEX27, TET15, and WEDGE21. All checks are functions called by `ex47.cpp` main. The class needs only standard C++14 and the FEMuS element enum; the test main additionally uses the actual FEMuS bases, Boost multiprecision/Gauss-Kronrod headers, and Eigen. The local `include/CutFemWeight_LS.hpp` duplicate is `Ex47::CutFemWeight_LS`; the production `::CutFemWeight_LS` remains available for comparisons. Production integration is unchanged.

Build independently of MPI/PETSc:

```sh
cmake -S applications/MultiPhase/ex47 -B build-ex47 -DCMAKE_BUILD_TYPE=Release
cmake --build build-ex47 --parallel 2
./build-ex47/MultiPhase_ex47
./build-ex47/MultiPhase_ex47 --cases 256 --points 64 --benchmark
./build-ex47/MultiPhase_ex47 --benchmark-quadrature
./build-ex47/MultiPhase_ex47 --benchmark-surface
./build-ex47/MultiPhase_ex47 --benchmark-data
./build-ex47/MultiPhase_ex47 --benchmark-algebra
./build-ex47/MultiPhase_ex47 --test-weights --weight-dimension 1
./build-ex47/MultiPhase_ex47 --test-weights --weight-dimension 2
./build-ex47/MultiPhase_ex47 --test-weights --weight-dimension 3
./build-ex47/MultiPhase_ex47 --benchmark-weights
```

Weight integration details and comparison methodology are in [INTEGRATION.md](INTEGRATION.md).

If Eigen is outside the compiler search path, set `-DEX47_EIGEN_INCLUDE_DIR=/path/to/eigen3`. If Boost is outside the compiler's search path, add `-DEX47_BOOST_INCLUDE_DIR=/path/to/boost`. The normal FEMuS build also exposes `MultiPhase_ex47`. `run_ex47.sh` offers the standalone build and accepts the main's arguments; `EX47_BUILD_DIR`, `EX47_JOBS`, `EX47_BOOST_INCLUDE_DIR`, and `CMAKE` override its defaults.

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
| Original gradient components | `GetRootGradient(point)`, `GetRootDerivativeCoefficients(axis)` |
| Optional original-gradient context | `Restrict(axis, coordinate, child, retainRootGradient)`, `HasRootGradientContext()` |
| Allocation-free gradients | `GetGradient(point, std::array<long double,3>&)`, `GetRootGradient(point, std::array<long double,3>&)` |
| External geometry polynomials | `RestrictComponents(axis, coordinate, degree, components, output)` |
| Nodal polynomial expansion without classification | `static ExpandNodalPolynomial(geometry, longDoubleNodalValues, coefficients)` |
| One-dimensional roots | `GetLineRoots()`; repeated roots appear once, endpoint roots are retained |
| Allocation-free roots | `GetLineRoots(std::array<long double,2>&)` returns the root count |
| Shared regular-root data | `GetLineSurfaceData(coordinate, rootGradient)` returns `1/abs(dG/dt)` and all original reference-gradient components |
| Delta only | `GetLineDeltaFactor(coordinate)` |
| Known regular-root slope | `static GetRegularRootDeltaFactor(root)` |
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

The class owns its nodal snapshot. Coefficients, normalization, derivatives, and Hessians are prepared once per state. Changing depth reuses polynomial data. A reset with identical nodal data preserves preparation; changed nodal data invalidates it. Context metadata does not alter the polynomial. Slice sign certification is evaluated lazily on the first classification query; value/root-only traversals avoid unused certification. Slices own their restricted state and remain valid if the parent is reset or destroyed. Reusing a child retains its coefficient/derivative storage and invalidates slice-dependent caches.

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

## Prepared topology, boundaries, and charts

| Operation | Interface and contract |
| --- | --- |
| Complete LINE analysis | `GetLineAnalysis()` returns a retained record with degree, zero-fiber status, isolated roots, signed derivatives, multiplicities, endpoint tags, and up to two positive intervals |
| Root densities | Every selected/batch surface evaluator also accepts `const LineRoot&`; use records from the current LINE analysis |
| Coordinate compatibility | `GetLineRoots()` uses the retained analysis; coordinate-based density methods reuse its stable slope when their coordinate matches a retained root |
| Exact zero polynomial | `IsIdenticallyZero()` distinguishes an entire zero fiber from isolated roots |
| Boundary information | `GetBoundaryCount()`, `IsZeroOnBoundary(boundary)`, `GetZeroBoundaryMask()` |
| Boundary restriction | `RestrictBoundary(boundary, child)` returns the same data type with original-gradient context and coordinate maps |
| Axis chart | `PermuteAxes(axes, chart)` uses `axes[newAxis] = oldAxis`; WEDGE preserves its height axis |
| Barycentric chart | `PermuteSimplexVertices(vertices, chart)` uses `vertices[newVertex] = oldVertex` for TRI/TET and the WEDGE triangle |
| External chart fields | `chart.TransformParentComponents(degree, components, outputDegree, output)` composes every scalar component with the chart's parent map |
| Graph information | `GetGraphCertificate(axis)` reports derivative bounds, sign, strictness, and zero-fiber status |
| Weak directions | `GetGraphDirection(true)` allows a weak monotone direction only when zero fibers have been excluded; the default accepts strict directions |
| Zero-fiber information | `GetZeroFiberStatus(axis)` is +1 for a detected zero fiber, -1 when excluded, and 0 when unresolved |
| Structural information | `DependsOn(axis)` and `GetReflectionParity(axis)` inspect actual retained coefficients without tolerance-zeroing monomials |
| Projection algebra | `GetFiberProjection(axis)` retains the three normalized fiber-coefficient polynomials and their full degree-four discriminant |
| Stable discriminant | `EvaluateFiberDiscriminant(axis, outerPoint)` evaluates the fiber coefficients with compensation; planar data also accept a scalar coordinate |
| Boundary/event algebra | `GetFiberBoundaryPolynomial(axis, constant, slope, output)` and `GetShiftedFiberDiscriminant(axis, endpoint, output)` expose planar projection polynomials |
| Unit-coordinate adapter | `GetUnitPowerCoefficients()` substitutes tensor coordinates `xi = 2*u-1`; `GetUnitToReferenceMeasureScale()` returns the corresponding reference-volume factor |

Boundary numbering is local to this algebraic interface: tensor boundaries are `2*axis` at -1 and `2*axis+1` at +1; simplex boundary 0 is `lambda0 = 0` and boundary `axis+1` is `xi[axis] = 0`. WEDGE uses its three triangle boundaries followed by height -1 and +1. Map these indices explicitly when connecting FEMuS face-connectivity numbering. LINE boundaries are points and support zero queries, while positive-dimensional boundary restriction starts at dimension two.

Boundary restriction includes sloping simplex faces. Its parent measure is the Euclidean reference embedding Jacobian obtained from the returned map; this agrees with the existing coordinate-slice factor on coordinate faces. The root map and all original gradient components retain their original reference ordering. Supply cofactors in that original reference basis, or compose prepared original-component fields with the chart. Vertex charts use closure of the supported enriched FEM spaces; external fields use full polynomial composition and retain any increase in coordinate degree. Degrees exceeding the existing limit of 32 are rejected explicitly.

For a simple quadratic root the retained derivative magnitude is `sqrt(b*b-4*a*c)`, including the normalization scale. The root record avoids subtraction in `2*a*rounded_root+b`. Root-surface gradients also enforce the directional chain rule with that slope, preserving the gradient/denominator cancellation in the reference and matrix-based measures. Near cancellation, LINE analysis recovers its coefficients from an immutable retained parent snapshot using compensated restriction, including simplex rescaling. This handles both nearly repeated roots and a quadratic coefficient rounded away during ordinary contraction. Positive intervals follow root topology instead of sign sampling between very close roots. Repeated roots and identically-zero fibers remain explicit requests for special integration handling; no half/full boundary contribution convention is imposed here.

The parent snapshot is allocated once when first needed by restriction and shared immutably by its descendants. LINE, graph, and projection caches are optional; copied records share prepared results and detach before a reused workspace modifies them. No unrequested surface density or geometry field is prepared. Prepare a chosen chart or boundary outside point loops, then retain it and reuse one child per lower dimension. The sloping-boundary convenience method constructs a temporary chart; a repeated traversal can prepare that chart once explicitly and restrict it directly.

Default graph selection requests strict bounds only. It uses Bernstein bounds and a QUAD corner-sign rejection before requesting tighter quadratic extrema; weak metadata and zero-fiber searches remain lazy. Weak directions retain their zero-fiber qualification. Unresolved three-dimensional zero-fiber searches remain unresolved and cannot enable a weak graph path. Reflection parity is +1 for even, -1 for odd, and 0 otherwise; coordinate reflections that do not preserve a simplex domain return 0. The integration layer must also check geometry/density parity before using symmetry or independence.

Projection coefficients are normalized by the parent's positive scale, and the discriminant therefore by its square. Auxiliary discriminants and boundary polynomials have separate degree-four storage; they are not quadratic level-set records. Compensated discriminant evaluation is the sensitive-event path; on a 3D-derived planar slice it also recovers the retained ancestor information. Expanded-polynomial evaluation can lose relative accuracy near a zero. Endpoint placement, endpoint transformations, resultant-root searches, integration budgets, moment restoration, and final fitted weights remain integration responsibilities.

An endpoint-mapped sample can retain a nonzero offset even when its absolute coordinate rounds to the endpoint. Preserve that integration context and pass its known root and signed native-chart derivative in a `LineRoot` to the root-aware evaluators; ordinary analysis of the rounded coordinate cannot reconstruct the lost offset. This is an explicit requirement for the future adapter to the existing endpoint kernels.

`--benchmark-data` measures fresh slicing, fresh slicing plus roots/physical density, and retained roots/densities. `--benchmark-algebra` measures cached metadata, supported chart construction, compensated discriminants, and external component transforms. Both modes assert zero warmed C++ allocations in allocation-tracking builds; correctness checks run in the normal main. See the performance and validation reports for measured results and limits.
