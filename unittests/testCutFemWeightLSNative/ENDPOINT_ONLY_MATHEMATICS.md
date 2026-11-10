# Endpoint transformations and the supported local models

Starting implementation: `8bad259c5a51060e4c2674db49a7ae34b75a0367`, the verified tip of `rifat-cutfem-2d-endpoint-quadrature` and direct child of `aea2c13efd1ac883103461d3218790cfbe70e4f0`. Benchmark: `0ff0669494e603be248d02c44c0a2e9b81f72b4d` on `rifat`.

## Definition of p

The left strip uses `s = a + delta*t^p`, `0 <= t <= 1`. Its positive integration Jacobian is `p*delta*t^(p-1)`. The right strip uses `s = b - delta*t^p`, with the same positive measure after reversing the limits. Neither Jacobian is a geometric surface metric. The implementation multiplies the strip Jacobian, the slice-coordinate Jacobian, and the original surface metric once each.

If an integrand has leading behavior `s^alpha`, its transformed leading behavior is `t^(p*(alpha+1)-1)`. For an integrable endpoint (`alpha > -1`), boundedness requires `p >= 1/(alpha+1)`. Smoothness requires more than boundedness: fractional exponents in the expansion must also be cleared. There is no largest admissible p; a larger exponent can work but crowd samples near the endpoint and worsen floating-point conditioning. Coefficient steepness affects constants, strip widths, and conditioning, but the required powers depend on contact multiplicities.

## An ordinary inner fold

For `G = A*z^2+B*z+C`, a simple fold has `A != 0`, `D = B^2-4*A*C = 0`, a simple transverse zero of D, and a nonzero full level-set gradient. At its two roots, `|G_z| = sqrt(D)`. Surface/coarea density therefore contains `s^(-1/2)`. Its minimum integer exponent for boundedness is p=2; that substitution also removes the ordinary square-root branch. The volume fiber integral is nonsingular, but it benefits from the same map.

## A projected 3D endpoint is a different issue

For HEX27, A, B, and C have degree at most two **in each** remaining coordinate. D can have degree four in each remaining coordinate and total degree eight. It must remain a full polynomial in the projection calculation; fitting it at nine points changes the problem.

Near a regular projected critical point, consider the transverse unfolding `D(u,s) = c*s+d*u^m+...`, `c*d != 0`, `A != 0`. On a contracting branch, `u = O(s^(1/m))`, so integrating the inner fold density over u gives leading behavior `s^(1/m-1/2)`. For this local model:

| Contact m | Leading exponent alpha | Minimum integer p for boundedness | Power clearing both the m-th-root and square-root branches |
|---|---:|---:|---:|
| 2 | 0 | 1 | 2 |
| 3 | -1/6 | 2 | 6 |
| 4 | -1/4 | 2 | 4 |

The last column is a sufficient common power `lcm(2,m)`, not an assertion that the leading scalar term alone requires that power. The code takes the least common multiple when several supported contacts meet an event; a common power for m=2,3,4 is 12. These statements assume the specified transverse local model. They do not prove that every arrangement of intersections of projected constraints has that model.

Actual triquadratic examples include `G=z^2+(2*y+y^2)*z+y^2+s`, with `D=4*y^3+y^4-4*s`, and `G=z^2+y^2*z+s`, with `D=y^4-4*s`. Both are regular at the origin because `G_s=1`, despite the cubic or quartic projected contact. A plain p=2 map makes the corresponding model slice integrands bounded but leaves fractional powers.

## FEMuS TET15

The actual enriched basis spans P2 plus the four cubic barycentric face bubbles and the quartic interior bubble `lambda0*lambda1*lambda2*lambda3`. With `lambda0=1-x-y-z`, its original polynomial has total degree at most four and degree at most two in each reference coordinate. It is not a plain P2 tetrahedron. Its coordinate-direction discriminant can again be quartic in a remaining coordinate. Coordinate restrictions retain these enrichment terms. A restriction to a tetrahedron face is cubic because the interior bubble vanishes there.

The following examples belong to this enriched space, with `lambda0=1-s-y-z`:

- `G=z^2+lambda0*y*z+y^2/4+s` has `D(y,0)=-y^3*(1-y)`.
- `G=z^2+lambda0*y*z-y*z+s` has `D(y,0)=y^4`.

Both have `G_s=1` at the origin. Thus enrichment does not justify assuming that every projected endpoint is an ordinary quadratic fold. The tests also compare expansion and gradients against the actual FEMuS basis rather than a substitute polynomial basis.

## Cases requiring a separate conclusion

- **Vanishing quadratic coefficient:** A=0 and B nonzero is a regular linear fiber and is evaluated as such. A=B=0 and C nonzero contributes no finite root. A=B=C=0 is an entire zero fiber; a valid transverse direction may handle it, otherwise integration fails explicitly. The quadratic discriminant formula alone does not classify this case.
- **Repeated factors:** a zero Sylvester resultant can indicate a persistent repeated projection factor. The general projection route rejects it rather than treating the resultant as an empty event list. A repeated quadratic trace on the boundary of an already certified regular 3D graph is treated separately: the bounded parent graph density is retained and the finite boundary branches use quadratic endpoint maps. Repeated factors of G, such as `G=z^2`, can make the raw delta moment divergent; an endpoint substitution cannot repair that integral.
- **Non-transverse critical endpoints:** where the derivative with respect to s also vanishes, the general surface implementation does not claim the transverse local model. The multiplicity classifier rejects the detected case, except for a certified isolated face tangency: the face Hessian must be definite and the original parent derivative transverse to that face must be nonzero. For example, `G=2*u-u*u-v*v-w*w` is regular at the origin (`G_u=2`), although its face trace `G(0,v,w)=-v*v-w*w` has zero gradient. The implicit function theorem gives the regular parent graph `u=1-sqrt(1-v*v-w*w)`; this isolated face event needs no higher exponent than the retained p=2 map. This certificate is applied only to face traces, never to the fiber discriminant. Boundary traces of a separately certified regular graph remain distinguished from singular parent fields as before.
- **Regular saddle projections:** even a regular level set can produce a logarithmic slice singularity. For `G=z^2-y^2-s`, a slice contains `integral dy/sqrt(y^2+s)`, whose leading behavior is `-log(s)`. A power map with p>1 makes the transformed density bounded, but no finite p makes the logarithm analytic. Quadrature may converge, but an analytic-convergence claim would be false.
- **Floating-point resolution:** narrow strips retain their events and contributions. When GK15 nodes would round to an endpoint, the same transformed-integrand quadrature uses GK7 or GK3 only if those nodes remain representable. Its embedded error estimate includes the change at an adjacent representable slice and that slice's inner error. The usual selective and accumulated moment checks remain mandatory. If no rule resolves the strip, or its error does not pass, integration fails explicitly; no interval is discarded or interpolated. Exhausted quadrature budgets also remain explicit failures.

There is consequently no finite universal exponent that makes every HEX27/TET15 input regular and analytic, or even makes every singular input integrable. There is a finite common power (12) for the stated transverse algebraic contact models of orders 2,3,4. The implementation and numerical results must be assessed within those assumptions. Classification tests and analytical model quadrature do not, by themselves, establish successful full 3D integration of every cubic/quartic projected geometry. The separate `endpoint_nested_contacts` test compares forced non-graph projections of representative cubic/quartic contacts in both HEX and the enriched TET space against an independently integrated regular graph direction, for raw delta and reference surface moments.

## Numerical checks and architecture

Projection events are refined against the original slice polynomial using compensated long-double sums and products. This keeps a boundary crossing and a nearly coincident discriminant zero consistent without replacing the level set or increasing an event-merging tolerance. It also avoids creating a spurious integration interval on the opposite side of a fold.

A simple fold just outside a cell can leave a small positive endpoint discriminant D(a), producing a sharp but finite surface density. The same quadratic endpoint substitution is used when D(a) is small compared with the discriminant's inward variation across the strip. Its constant is retained: for h=+/-delta*t^2, `D(a+h)=D(a)+h*R(h)` and the transformed surface ratio is `2*delta/sqrt(D(a)/t^2 +/- delta*R(h))`. The quadrature still evaluates the actual transformed integrand, including D(a); no zero endpoint value is imposed on these near-fold strips.

`endpoint_exponents.cpp` checks p=2,6,4 for transverse contacts, verifies transformed model integrals analytically, checks that quartic triquadratic discriminants remain quartic, and requires rejection of a non-transverse model. The integration hierarchy is HEX -> QUAD -> LINE and TET -> full TRI -> LINE. LINE graph quadrature and exact quadratic two-root integration are retained; general curves use actual-integrand endpoint quadrature. General 3D surfaces use full parent slices, projection events, endpoint strips, and adaptive quadrature panels. There is no geometric subdivision dispatcher.

The derivations above are local calculations for this implementation. For broader background on polynomial projection and dimension reduction, see [Saye, JCP 2022](https://doi.org/10.1016/j.jcp.2021.110720); that reference is not a proof of this branch's supported input class.
