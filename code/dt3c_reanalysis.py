"""POST-HOC re-analysis of DT3c (Step 5 CONFIRM). Read-only: no new random draws.

Why: the frozen DT3c verdicts REFUTED L2/L3 only at the 0.50 tail point. There the formula needs
F_M(x) at x < 0.90, below the transcribed Balabdaoui-Filali table, and dt2_farm.FM used a guessed
extension 0.393*(x/0.9)**8. This script recomputes the SAME predictions with F_M from the independent
Step 6 check (deposits/dt4/fm_check.json, N = 1e7, grid 0.60..2.00 step 0.02) and re-applies the
frozen decision rule. The frozen verdicts stay the record; this is a labelled post-hoc readout.

Run: python code/dt3c_reanalysis.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import dt3c_confirm as dc  # noqa: E402
from dt1b_confirm import bf_table  # noqa: E402
from fm_table import FM_pub  # noqa: E402

W = HERE.parents[0]
FMJ = json.loads((W / "deposits" / "dt4" / "fm_check.json").read_text())
FX = np.array([r["x"] for r in FMJ["rows"]])
FV = np.array([r["F_M_est"] for r in FMJ["rows"]])


def fm_legacy(x):
    """Verbatim copy of dt2_farm.FM as used in the frozen DT3c run (before 2026-11-06 03:15Z)."""
    xs, fv = bf_table()
    if x < xs[0]:
        return max(0.0, fv[0] * (x / xs[0]) ** 8)
    if x <= xs[-1]:
        return float(np.interp(x, xs, fv))
    return 1 - (1 - fv[-1]) * np.exp(-8 * (x - xs[-1]))


def fm_verified(x):
    """F_M from the independent MC check; flags arguments outside its grid."""
    if x < FX[0]:
        raise ValueError(f"F_M needed at x={x:.3f} < {FX[0]} (outside the verified grid)")
    if x > FX[-1]:
        return float(1 - (1 - FV[-1]) * np.exp(-8 * (x - FX[-1])))  # > 2.0: F_M > 0.9987, tail negligible here
    return float(np.interp(x, FX, FV))


def pred(t, n, li, ref_sorted, fm):
    _, masses, mi = dc.spec(li, n)
    prod = 1.0
    for i, p in enumerate(masses):
        if i == mi:
            prod *= np.searchsorted(ref_sorted, t / np.sqrt(p), side="right") / len(ref_sorted)
        else:
            prod *= fm(2 * t / np.sqrt(p))
    return 1.0 - prod


def main():
    d = json.loads((dc.OUT / "confirm.json").read_text())
    ref = np.concatenate([np.load(dc.part_path(("ref", None, c))) for c in range(dc.REF_R // dc.CHUNK)])
    refs = {n: np.sort(np.sqrt(n) * ref[:, j]) for j, n in enumerate(dc.NS)}
    cells = {(c["law"], c["n"]): c for c in d["cells"]}
    out = {"note": "POST-HOC: same draws, F_M from fm_check.json instead of the x<0.9 extension", "laws": {}}
    for li in range(len(dc.BF_LAWS)):
        rows = {}
        for n in dc.NS:
            c = cells[(li, n)]
            r = {}
            for p in dc.TAILS:
                tp = c["tails"][str(p)]
                old = fm_legacy  # the function used in the frozen run
                pv = pred(tp["t"], n, li, refs[n], fm_verified)
                pp = pred(tp["t"], n, li, refs[n], lambda x: FM_pub(x, clip_tails=True))
                po = pred(tp["t"], n, li, refs[n], old)
                args = [2 * tp["t"] / np.sqrt(m) for i, m in enumerate(dc.spec(li, n)[1]) if i != dc.spec(li, n)[2]]
                r[str(p)] = {"t": tp["t"], "obs": tp["obs"], "pred_frozen": po, "pred_verified": pv,
                             "pred_published": pp, "diff_frozen_pts": 100 * (tp["obs"] - po),
                             "diff_verified_pts": 100 * (tp["obs"] - pv), "diff_published_pts": 100 * (tp["obs"] - pp),
                             "min_FM_arg": float(min(args)), "below_table": bool(min(args) < 0.90)}
            rows[n] = r
        c5, c51 = rows[dc.NS[-1]], rows[dc.NS[-2]]
        maxtail = max(abs(v["diff_verified_pts"]) for v in c5.values())
        refuted = any(abs(c51[k]["diff_verified_pts"]) > 1.5 and abs(c5[k]["diff_verified_pts"]) > 1.5
                      and np.sign(c51[k]["diff_verified_pts"]) == np.sign(c5[k]["diff_verified_pts"])
                      and abs(c5[k]["diff_verified_pts"]) >= abs(c51[k]["diff_verified_pts"]) - 0.3 for k in c5)
        fr = d["readings"][dc.NAMES[li]]
        q99a, q99b = fr["q99absR_by_n"][0], fr["q99absR_by_n"][-1]
        conf = (abs(fr["diff5_by_n"][-1]) <= 1.0 and abs(fr["diff1_by_n"][-1]) <= 0.7 and maxtail <= 1.5
                and (q99b <= 0.5 * q99a or q99b <= 0.002))
        verdict = "REFUTED" if refuted else ("CONFIRMED" if conf else "INCONCLUSIVE (slow)")
        out["laws"][dc.NAMES[li]] = {"frozen_verdict": fr["verdict"], "posthoc_verdict": verdict,
                                     "max_tail_diff_204800_verified": round(maxtail, 3),
                                     "rows": {str(n): rows[n] for n in dc.NS}}
        print(f"{dc.NAMES[li]:28s} frozen={fr['verdict']:12s} post-hoc={verdict:12s} maxtail={maxtail:.2f}")
        for n in dc.NS:
            print("   n", n, " ".join(f"p{k}:{v['diff_frozen_pts']:+.2f}->{v['diff_verified_pts']:+.2f}"
                                     f"{'*' if v['below_table'] else ''}" for k, v in rows[n].items()))
    pubmax = max(abs(v["diff_published_pts"] - v["diff_verified_pts"]) for law in out["laws"].values()
                 for r in law["rows"].values() for v in r.values())
    out["max_abs_published_minus_mc_pts"] = pubmax
    print(f"published-table vs MC evaluator: max difference in any prediction = {pubmax:.4f} pt")
    (dc.OUT / "reanalysis_posthoc.json").write_text(json.dumps(out, indent=1))
    print("saved", dc.OUT / "reanalysis_posthoc.json", "(* = some F_M argument below the 0.90 table start)")


if __name__ == "__main__":
    main()
