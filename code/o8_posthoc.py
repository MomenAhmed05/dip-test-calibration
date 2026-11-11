"""O8 audit, POST-HOC readout (no new draws): the frozen audit null "ramp" is not unimodal.

o8_audit.s_ramp draws 0.9 * U on (0, 0.9) w.p. 0.9 and 0.9 + 0.1 * sqrt(U) w.p. 0.1, i.e. density 1 on (0, 0.9) and
density 20 (x - 0.9) on (0.9, 1): it drops from 1 to 0 at x = 0.9 and rises again to 2 at x = 1 (two modes).
It was meant to be the DT1 flat-plus-ramp null (ramp starting at the flat level). The frozen labels
(audit.json 'worst') are kept as computed; this file applies the same O8-R1/O8-R2 rules with that cell removed
(Bonferroni m = 6 methods x 33 null cells), and the O8-R2 baseline comparison.
Output: deposits/o8/audit_posthoc.json
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta

SRC = Path(__file__).resolve().parents[1] / "deposits" / "o8"
ALPHA, TOL = 0.05, 0.055
CONST = ("LF", "lam0.75", "lam0.5", "lam0.25")  # most to least conservative; LF = lambda 1


def main():
    try:
        a = json.loads((SRC / "audit.json").read_text())
    except (OSError, json.JSONDecodeError) as e:
        raise SystemExit(f"cannot read audit.json: {e}") from e
    R = a["R"]
    nulls = [r for r in a["rows"] if r["role"] == "null" and r["law"] != "ramp"]
    methods = list(nulls[0]["rej"])
    m = len(nulls) * len(methods)
    labels = {}
    for k in methods:
        ub = [float(beta.ppf(1 - 0.05 / m, r["rej"][k] + 1, R - r["rej"][k])) for r in nulls]
        lb = [float(beta.ppf(0.05 / m, r["rej"][k], R - r["rej"][k] + 1)) if r["rej"][k] > 0 else 0.0 for r in nulls]
        w = int(np.argmax([r["rej"][k] for r in nulls]))
        lab = ("certified within tolerance" if max(ub) <= TOL else
               "demonstrated above tolerance" if max(lb) > TOL else "inconclusive")
        labels[k] = {"label": lab, "worst_law": nulls[w]["law"], "n": nulls[w]["n"],
                     "rate": nulls[w]["rej"][k] / R, "ub": ub[w], "max_ub": max(ub)}
    certified = [k for k in CONST if labels[k]["label"] == "certified within tolerance"]
    best = certified[-1] if certified else None  # most powerful certified constant baseline
    gains = []
    if best:
        for r in a["rows"]:
            if r["role"] == "alt" and r["n"] == 800:
                lo, hi = r["disc_learned_vs"][best]
                gains.append({"law": r["law"], "gain": (lo - hi) / R, "learned_only": lo, "baseline_only": hi})
    mean_gain = float(np.mean([g["gain"] for g in gains])) if gains else None
    out = {"note": "post-hoc: ramp removed (not unimodal); same rules", "bonferroni_m": m, "labels": labels,
           "best_certified_constant": best, "gains_n800_vs_best": gains, "mean_gain_n800": mean_gain,
           "learning_helped_(>=2pt)": (mean_gain is not None and labels["learned"]["label"] ==
                                       "certified within tolerance" and mean_gain >= 0.02)}
    (SRC / "audit_posthoc.json").write_text(json.dumps(out, indent=1))
    for k, v in labels.items():
        print(f"{k:8s} {v['label']:30s} worst {v['worst_law']} n={v['n']} {100 * v['rate']:.2f}% UB {100 * v['ub']:.2f}%")
    print("best certified constant:", best, "| mean paired gain of learned at n=800:",
          None if mean_gain is None else round(100 * mean_gain, 2), "pt")
    for g in gains:
        print(f"  {g['law']:12s} gain {100 * g['gain']:+.2f} pt ({g['learned_only']}/{g['baseline_only']})")


if __name__ == "__main__":
    main()
