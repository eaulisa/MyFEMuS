#!/usr/bin/env python3
"""Analyze rank-local surface/pressure contributions exported by ex46.

The operator G is the pressure contribution to the FEMuS RES velocity rows;
force F is the physical capillary contribution to those same rows.
The measured constrained range distance is min_{p_free} ||F + Gp|| / ||F||.
"""
import argparse
import csv
import glob
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsmr


def read_rows(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        raise RuntimeError(f"No files matching {pattern}")
    for filename in files:
        with open(filename, newline="") as stream:
            for row in csv.DictReader(stream):
                yield row


def ldnorm(x):
    z = np.asarray(x, dtype=np.longdouble)
    return float(np.sqrt(np.sum(z*z, dtype=np.longdouble)))


def matvec_ld(A, x):
    # Recompute residual after LSMR with extended-precision row accumulation.
    x = np.asarray(x, dtype=np.longdouble)
    out = np.zeros(A.shape[0], dtype=np.longdouble)
    for i in range(A.shape[0]):
        k, l = A.indptr[i:i+2]
        out[i] = np.sum(A.data[k:l].astype(np.longdouble) *
                        x[A.indices[k:l]], dtype=np.longdouble)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix", help="Common file prefix passed to --surface-pressure-prefix")
    ap.add_argument("--maxiter", type=int, default=6000)
    ap.add_argument("--max-relative", type=float, default=None)
    args = ap.parse_args()

    metadata = {}
    with open(args.prefix + "_meta.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            metadata[row["key"]] = float(row["value"])
    jump = metadata["pressure_jump"]

    # Owner-derived metadata; duplicate ownership or disagreements are errors.
    dofs = {}
    for record in read_rows(args.prefix + "_dofs_rank*.csv"):
        sid = int(record["sys_dof"])
        if sid in dofs:
            raise RuntimeError(f"Duplicate owner metadata for system DOF {sid}")
        field = record["field"]
        bdc = float(record["bdc"])
        val = float(record["fixed_value"])
        if field not in ("U", "V", "W", "P1", "P2"):
            raise RuntimeError(f"Unexpected field {field}")
        dofs[sid] = (field, bdc, val)

    velocity_rows = sorted(k for k, (field, bdc, _) in dofs.items()
                           if field in ("U", "V", "W") and bdc > 1.5)
    pressure_cols = sorted(k for k, (field, _, _) in dofs.items()
                           if field in ("P1", "P2"))
    if not velocity_rows or not pressure_cols:
        raise RuntimeError("No velocity rows or pressure columns")
    ri = {key: i for i, key in enumerate(velocity_rows)}
    ci = {key: i for i, key in enumerate(pressure_cols)}

    # Duplicate entries from distinct elements/ranks must be ADDED.
    ii, jj, vv = [], [], []
    for rec in read_rows(args.prefix + "_bt_rank*.csv"):
        row, col = int(rec["row"]), int(rec["col"])
        if row not in dofs or col not in dofs:
            raise RuntimeError(f"Matrix references DOF without owner metadata: {row}, {col}")
        if row in ri:
            if col not in ci:
                raise RuntimeError(f"Matrix column {col} is not a pressure DOF")
            ii.append(ri[row]); jj.append(ci[col]); vv.append(float(rec["value"]))

    G = coo_matrix((vv, (ii, jj)), shape=(len(ri), len(ci))).tocsr()
    G.sum_duplicates()
    G.eliminate_zeros()
    if not np.all(np.isfinite(G.data)):
        raise RuntimeError("Non-finite pressure operator")

    F = np.zeros(len(ri), dtype=np.longdouble)
    direct = np.zeros(len(ri), dtype=np.longdouble)
    for rec in read_rows(args.prefix + "_forces_rank*.csv"):
        row = int(rec["row"])
        if row not in dofs:
            raise RuntimeError(f"Force references unknown DOF {row}")
        if row in ri:
            F[ri[row]] += np.longdouble(rec["value"])
    for rec in read_rows(args.prefix + "_action_rank*.csv"):
        row = int(rec["row"])
        if row not in dofs:
            raise RuntimeError(f"Action references unknown DOF {row}")
        if row in ri:
            direct[ri[row]] += np.longdouble(rec["value"])

    f_norm = ldnorm(F)
    if f_norm == 0.0 or not np.isfinite(f_norm):
        raise RuntimeError("Zero or non-finite capillary force")

    p_nominal = np.array([jump if dofs[col][0] == "P1" else 0.
                          for col in pressure_cols], dtype=np.longdouble)
    matrix_action = matvec_ld(G, p_nominal)
    pressure_mismatch = ldnorm(matrix_action - direct)
    action_norm = ldnorm(direct)
    action_relative = pressure_mismatch / max(action_norm, f_norm)
    if action_relative > 1.e-10:
        raise RuntimeError(f"Pressure action check failed: relative mismatch {action_relative}")

    # The YL reference is a physical constant jump before pressure pinning.
    # The range calculation below enforces the production fixed-pressure values.
    yl_abs = ldnorm(F + direct)
    yl_rel = yl_abs / f_norm

    free = np.array([dofs[col][1] > 1.5 for col in pressure_cols], dtype=bool)
    fixed_val = np.array([0. if free[j] else dofs[col][2]
                          for j, col in enumerate(pressure_cols)], dtype=float)
    G_free = G[:, free].tocsr()
    rhs = -F.astype(float) - G @ fixed_val
    colnorm = np.sqrt(np.asarray(G_free.power(2).sum(axis=0)).ravel())
    active = np.flatnonzero(colnorm > 0)
    if not len(active):
        raise RuntimeError("No active free pressure columns")
    A = G_free[:, active].multiply(1. / colnorm[active]).tocsr()
    A.eliminate_zeros()

    fit = lsmr(A, rhs, atol=0, btol=0, conlim=1e15, maxiter=args.maxiter)
    if fit[1] in (3, 6, 7):
        raise RuntimeError(f"LSMR did not converge: stop={fit[1]}, iterations={fit[2]}")
    p_free = np.zeros(G_free.shape[1], dtype=float)
    p_free[active] = fit[0] / colnorm[active]
    p_opt = fixed_val.copy()
    p_opt[free] = p_free
    residual = F + matvec_ld(G, p_opt.astype(np.longdouble))
    range_abs = ldnorm(residual)
    range_rel = range_abs / f_norm

    # Optimality on the scaled least-squares problem.
    normal = ldnorm(A.T @ np.asarray(residual, dtype=float))
    rzero = np.random.default_rng(12345).standard_normal(len(ri))
    neg = lsmr(A, -rzero, atol=0, btol=0, conlim=1e15, maxiter=args.maxiter)
    if neg[1] in (3, 6, 7):
        raise RuntimeError("LSMR negative control failed to converge")
    neg_rel = ldnorm(rzero + A @ neg[0]) / ldnorm(rzero)
    if len(ri) > len(active) and neg_rel < 0.1:
        raise RuntimeError(f"Negative control failed: {neg_rel}")
    if args.max_relative is not None and range_rel > args.max_relative:
        raise RuntimeError(f"Range distance {range_rel} exceeds {args.max_relative}")

    result = {
        "status": "converged", "free_velocity_dofs": len(ri),
        "pressure_dofs_total": len(ci), "active_pressure_dofs": len(active),
        "surface_norm": f_norm,
        "prescribed_absolute": yl_abs, "prescribed_relative": yl_rel,
        "range_absolute": range_abs, "range_relative": range_rel,
        "pressure_action_mismatch": pressure_mismatch,
        "pressure_action_mismatch_relative": action_relative,
        "lsmr_iterations": int(fit[2]), "lsmr_stop": int(fit[1]),
        "normal_residual": normal,
        "negative_control_relative": neg_rel,
        "metadata": metadata,
    }
    with open(args.prefix + "_result.json", "w") as fh:
        json.dump(result, fh, indent=2)
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
