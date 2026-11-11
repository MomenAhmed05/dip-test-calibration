"""Step 7 (DT5) analysis: intervals and the frozen decision rules R1-R5 of work/DT5-PREREG.md.

Reads deposits/dt5[_quick]/power.json and the per-job parts written by dt5_power.py. Writes analysis.json and
DT5-TABLES.md next to them. No random numbers.
Run: python code/dt5_analysis.py [--quick]
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import beta, binom

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dt5_power as P  # noqa: E402  (same QUICK flag through sys.argv)

TOL = {0.05: 0.005, 0.01: 0.0025}  # R1 practical tolerance Delta_alpha
DIP_METHODS = ("uniform", "LF", "kuiper2", "dkw")
SSE_METHODS = ("sse_asymptotic", "sse_uniformcal")
PIN = ("pin(1.4,.1)", "pin(2.5,1)")


def cp(k, n, conf=0.95):
    """Two-sided Clopper-Pearson interval."""
    a = 1 - conf
    lo = 0.0 if k == 0 else float(beta.ppf(a / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a / 2, k + 1, n - k))
    return lo, hi


def load_law(li, n):
    return np.concatenate([np.load(P.part_path(("law", li, n, c))) for c in range(P.LAW_R // P.CHUNK)])


def holm(pvals):
    order = np.argsort(pvals)
    m = len(pvals)
    adj = np.empty(m)
    run = 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (m - r) * pvals[i]))
        adj[i] = run
    return adj


def main():
    try:
        d = json.loads((P.OUT / "power.json").read_text())
    except (OSError, json.JSONDecodeError) as e:
        print(f"cannot read {P.OUT / 'power.json'}: {e} (run dt5_power.py first)")
        return
    crit = {int(n): {float(a): v for a, v in c.items()} for n, c in d["crit"].items()}
    man = d["manifest"]
    null_idx = [m["index"] for m in man if m["pop_modes"] == 1]
    alt_idx = [m["index"] for m in man if m["pop_modes"] >= 2]
    cells, nest_viol, order_ok = {}, {}, {}
    for m in man:
        li = m["index"]
        for n in P.NS:
            arr = load_law(li, n)
            dip, th = arr[:, 0], arr[:, 1]
            N = len(dip)
            for a in P.ALPHAS:
                c = crit[n][a]
                rej = {k: dip > c[k] for k in DIP_METHODS}
                failed = int(np.sum(np.isnan(th)))
                rej["sse_asymptotic"] = th > c["sse_asymptotic"]
                rej["sse_uniformcal"] = th > c["sse_uniformcal"]
                for k, r in rej.items():
                    kk = int(np.sum(r))
                    cells[(li, n, a, k)] = {"k": kk, "N": N, "rate": kk / N, "ci": cp(kk, N),
                                            "failed_rows": failed if k in SSE_METHODS else 0}
                # R4: realised threshold order and dataset-level nesting
                order_ok[(n, a)] = bool(c["LF"] < c["kuiper2"] < c["dkw"])
                nest_viol[(li, n, a)] = int(np.sum(rej["kuiper2"] & ~rej["LF"]) + np.sum(rej["dkw"] & ~rej["kuiper2"]))
                lf_extra = int(np.sum(rej["LF"] & ~rej["kuiper2"]))
                cells[(li, n, a, "LF-kuiper2")] = {"k": lf_extra, "N": N, "rate": lf_extra / N, "ci": cp(lf_extra, N)}
                un_extra = int(np.sum(rej["uniform"] & ~rej["LF"]))
                cells[(li, n, a, "uniform-LF")] = {"k": un_extra, "N": N, "rate": un_extra / N, "ci": cp(un_extra, N)}

    # R1 per method, n, alpha over all null/control laws
    r1 = {}
    for k in DIP_METHODS + SSE_METHODS:
        for n in P.NS:
            for a in P.ALPHAS:
                ups = [cells[(li, n, a, k)]["ci"][1] for li in null_idx]
                los = [cells[(li, n, a, k)]["ci"][0] for li in null_idx]
                fails = sum(cells[(li, n, a, k)]["failed_rows"] for li in null_idx)
                if k in SSE_METHODS and fails:
                    status = "INVALID (failed rows)"
                elif max(ups) <= a + TOL[a]:
                    status = "certified within tolerance"
                elif max(los) > a + TOL[a]:
                    status = "demonstrated above tolerance"
                else:
                    status = "inconclusive"
                worst = int(np.argmax([cells[(li, n, a, k)]["rate"] for li in null_idx]))
                r1[(k, n, a)] = {"status": status, "max_upper": max(ups), "max_rate": max(cells[(li, n, a, k)]["rate"] for li in null_idx),
                                 "worst_law": man[null_idx[worst]]["law"],
                                 "laws_lower_above_alpha": [man[li]["law"] for li in null_idx if cells[(li, n, a, k)]["ci"][0] > a]}

    # R2: persistent finite-sample overshoot of LF on the pinning sequences (Holm over 4 claims)
    claims, pv = [], []
    for name in PIN:
        li = next(m["index"] for m in man if m["law"] == ["null", name])
        for a in P.ALPHAS:
            ps = []
            for n in (6400, 12800):
                c = cells[(li, n, a, "LF")]
                ps.append(float(binom.sf(c["k"] - 1, c["N"], a)))  # P(K >= k | size = alpha)
            claims.append((name, a))
            pv.append(max(ps))
    adj = holm(np.array(pv))
    r2 = [{"law": cl[0], "alpha": cl[1], "p_max": p, "p_holm": float(q), "alarm": bool(q < 0.05)}
          for cl, p, q in zip(claims, pv, adj, strict=True)]

    # power crossings (no monotonicity assumed: first registered n that reaches the level)
    cross = {}
    for li in alt_idx:
        for a in P.ALPHAS:
            for k in DIP_METHODS + SSE_METHODS:
                out = {}
                for lev in (0.5, 0.8):
                    pt = next((n for n in P.NS if cells[(li, n, a, k)]["rate"] >= lev), None)
                    lb = next((n for n in P.NS if cells[(li, n, a, k)]["ci"][0] >= lev), None)
                    out[str(lev)] = {"point": pt, "lower_bound": lb}
                cross[(li, a, k)] = out

    res = {"tier": d["tier"], "tolerance": {str(a): t for a, t in TOL.items()},
           "R1": [{"method": k, "n": n, "alpha": a, **v} for (k, n, a), v in r1.items()],
           "R2": r2,
           "R4": {"threshold_order_ok": {f"{n}/{a}": v for (n, a), v in order_ok.items()},
                  "nesting_violations_total": int(sum(nest_viol.values()))},
           "cells": [{"law": man[li]["law"], "n": n, "alpha": a, "method": k, **v}
                     for (li, n, a, k), v in cells.items()],
           "crossings": [{"law": man[li]["law"], "alpha": a, "method": k, **v} for (li, a, k), v in cross.items()]}
    (P.OUT / "analysis.json").write_text(json.dumps(res, indent=1))

    lines = [f"# DT5 tables ({d['tier']}; generated by code/dt5_analysis.py)", "",
             "## R1 size status per method (all null/control laws; two-sided CP-95% bounds; tolerance 0.5 pt at 5%, 0.25 pt at 1%)", "",
             "| method | alpha | " + " | ".join(str(n) for n in P.NS) + " |", "|---|---|" + "---|" * len(P.NS)]
    short = {"certified within tolerance": "C", "demonstrated above tolerance": "ABOVE", "inconclusive": "?",
             "INVALID (failed rows)": "INV"}
    for k in DIP_METHODS + SSE_METHODS:
        for a in P.ALPHAS:
            lines.append(f"| {k} | {a} | " + " | ".join(f"{short[r1[(k, n, a)]['status']]} {100 * r1[(k, n, a)]['max_rate']:.2f}"
                                                         for n in P.NS) + " |")
    lines += ["", "Cell = status, max null rate (%). C = certified within tolerance; ABOVE = demonstrated above tolerance; ? = inconclusive.", "",
              "## R2 LF on pinning sequences (Holm over 4 claims)", ""]
    lines += [f"- {r['law']} alpha={r['alpha']}: p_max={r['p_max']:.3g}, Holm={r['p_holm']:.3g}, alarm={r['alarm']}" for r in r2]
    lines += ["", f"## R4 sanity: threshold order LF < kuiper2 < dkw holds in {sum(order_ok.values())}/{len(order_ok)} (n, alpha); "
              f"dataset-level nesting violations: {sum(nest_viol.values())}", "",
              "## Power at alpha = 5% (%; bimodal laws)", "",
              "| law | method | " + " | ".join(str(n) for n in P.NS) + " |", "|---|---|" + "---|" * len(P.NS)]
    for li in alt_idx:
        for k in DIP_METHODS + SSE_METHODS:
            lines.append(f"| {man[li]['law']} | {k} | " + " | ".join(f"{100 * cells[(li, n, 0.05, k)]['rate']:.1f}" for n in P.NS) + " |")
    lines += ["", "## First registered n with power >= 80% at 5% (point / CP lower bound)", "",
              "| law | " + " | ".join(DIP_METHODS) + " |", "|---|" + "---|" * len(DIP_METHODS)]
    for li in alt_idx:
        lines.append(f"| {man[li]['law']} | " + " | ".join(
            f"{cross[(li, 0.05, k)]['0.8']['point']} / {cross[(li, 0.05, k)]['0.8']['lower_bound']}" for k in DIP_METHODS) + " |")
    (P.OUT / "DT5-TABLES.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:4 + 2 * len(DIP_METHODS + SSE_METHODS) + 2]))
    print("\n".join(f"R2 {r}" for r in r2))
    print(f"R4 order ok {sum(order_ok.values())}/{len(order_ok)}; nesting violations {sum(nest_viol.values())}")
    print("saved", P.OUT / "analysis.json", P.OUT / "DT5-TABLES.md")


if __name__ == "__main__":
    main()
