"""DT5 R5 + P6 report (post-hoc presentation of frozen rule R5; no new draws, no RNG).

The frozen analyser (dt5_analysis.py, sha b4fbe8ad) produced R1/R2/R4 and the crossings,
but not the R5 paired-gain intervals or the P6 sample-size ratios. This file computes them
from deposits/dt5/analysis.json only. R4 found 0 nesting violations, so on every dataset
rejection by kuiper2 (or dkw) implies rejection by LF; the paired gain LF - other is then
(k_LF - k_other) / N with a Clopper-Pearson interval on the discordant count.

Output: deposits/dt5/r5_posthoc.json
"""
import json
import math
from pathlib import Path

from scipy.stats import beta

DEP = Path(__file__).resolve().parents[1] / "deposits" / "dt5"
NS = [50, 100, 200, 400, 800, 1600, 3200, 6400, 12800]


def cp(k, n):
    lo = float(beta.ppf(0.025, k, n - k + 1)) if k > 0 else 0.0
    hi = float(beta.ppf(0.975, k + 1, n - k)) if k < n else 1.0
    return lo, hi


def n_at(cells, law, method, p=0.8, alpha=0.05):
    """n where power first reaches p: logit(power) linear in log n between grid points."""
    pts = [(n, cells[(law, n, alpha, method)]["rate"]) for n in NS]

    def lg(r):
        r = min(max(r, 1e-6), 1 - 1e-6)
        return math.log(r / (1 - r))

    for (n0, r0), (n1, r1) in zip(pts, pts[1:], strict=False):
        if r0 < p <= r1:
            t = (lg(p) - lg(r0)) / (lg(r1) - lg(r0))
            return math.exp(math.log(n0) + t * (math.log(n1) - math.log(n0)))
    return None


def main():
    try:
        a = json.loads((DEP / "analysis.json").read_text())
    except (OSError, json.JSONDecodeError) as e:
        raise SystemExit(f"cannot read analysis.json: {e}") from e
    cells = {(json.dumps(c["law"]), c["n"], c["alpha"], c["method"]): c for c in a["cells"]}
    alts = sorted({k[0] for k in cells if not k[0].startswith('["null')})
    gains = []
    for law in alts:
        for n in NS:
            for alpha in (0.05, 0.01):
                lf = cells[(law, n, alpha, "LF")]
                row = {"law": law, "n": n, "alpha": alpha, "LF": lf["rate"]}
                for other in ("kuiper2", "dkw"):
                    o = cells[(law, n, alpha, other)]
                    d = lf["k"] - o["k"]
                    row[f"gain_vs_{other}"] = d / lf["N"]
                    row[f"gain_vs_{other}_ci"] = cp(max(d, 0), lf["N"])
                u = cells[(law, n, alpha, "uniform")]
                row["uniform_advantage"] = (u["k"] - lf["k"]) / lf["N"]
                gains.append(row)
    n80 = {}
    for law in alts:
        v = {m: n_at(cells, law, m) for m in ("uniform", "LF", "kuiper2", "dkw")}
        n80[law] = {
            "n80": v,
            "ratio_kuiper2_LF": v["kuiper2"] / v["LF"] if v["kuiper2"] and v["LF"] else None,
            "ratio_dkw_LF": v["dkw"] / v["LF"] if v["dkw"] and v["LF"] else None,
            "ratio_uniform_LF": v["uniform"] / v["LF"] if v["uniform"] and v["LF"] else None,
        }
    neg = [g for g in gains if g["gain_vs_kuiper2"] < 0 or g["gain_vs_dkw"] < 0]
    out = {"note": "post-hoc presentation of frozen R5/P6; no new draws", "gains": gains,
           "n80": n80, "negative_gain_cells": len(neg)}
    (DEP / "r5_posthoc.json").write_text(json.dumps(out, indent=1))
    rk = [r["ratio_kuiper2_LF"] for r in n80.values() if r["ratio_kuiper2_LF"]]
    rd = [r["ratio_dkw_LF"] for r in n80.values() if r["ratio_dkw_LF"]]
    ru = [r["ratio_uniform_LF"] for r in n80.values() if r["ratio_uniform_LF"]]
    print(f"n80 ratio kuiper2/LF {min(rk):.2f}-{max(rk):.2f} ({len(rk)} laws); dkw/LF {min(rd):.2f}-{max(rd):.2f}; "
          f"uniform/LF {min(ru):.2f}-{max(ru):.2f}; negative-gain cells {len(neg)}")
    best = max(gains, key=lambda g: g["gain_vs_kuiper2"])
    print("max gain vs kuiper2", round(100 * best["gain_vs_kuiper2"], 1), "pt", best["law"], best["n"], best["alpha"])
    print("saved", DEP / "r5_posthoc.json")


if __name__ == "__main__":
    main()
