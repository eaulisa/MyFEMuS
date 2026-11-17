# CutFemWeight_LS tasks

- [ ] Reverse explicit cut-detection checks: test nodal level-set signs first, then compute Bernstein coefficients only when all nodal values have the same strict sign.
- [ ] Expose the maximum Bernstein subdivision depth as a configurable parameter instead of the hard-coded depth of 12; allow depth 0 for a single Bernstein test without subdivision.
- [ ] Set a substantially lower default subdivision depth; the exact default remains to be chosen.
- [ ] Retain and share the parent power coefficients and parent Bernstein coefficients between cut detection and integration, including the different surface-weight computations; rebuild them only when the nodal level set or reference-element context changes, independently of moment order, integration tolerances, and Gauss-Kronrod settings.

Preserve the roundoff margins and classify opposite nodal signs, zero nodal values, and unresolved Bernstein cases as possible cuts.
