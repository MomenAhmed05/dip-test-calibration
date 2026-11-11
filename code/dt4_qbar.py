"""Step 6 (DT4): q-bar (uniform limit quantile of sqrt(n)*dip) = EVIDENTIAL, plus the size numbers.

Inputs (read-only): DT3c uniform references (deposits/dt3c/parts, 1e6 per n, nested CRN, root 2026110630, used),
the published F_M table (code/fm_table.py). No new random draws except the bootstrap below (root 2026110641, DEV).

Per n: point quantile + exact distribution-free order-statistic 95% CI (binomial).
Extrapolation models fitted to the 4 point estimates (weights 1/SE^2, SE from a bootstrap over the 1e6 rows that keeps
the nesting across n, so the correlation across n is respected):
  M1: q_n = qbar + a n^-1/2
  M2: q_n = qbar + a n^-1/2 + b log(n) n^-1/2
  M3: q_n = qbar + a n^-1/3
Report qbar per model, the spread across models, and the size numbers implied by every value in the range.

Run: python code/dt4_qbar.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import binom

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import dt3c_confirm as dc  # noqa: E402
from fm_table import FM_pub  # noqa: E402

W = HERE.parents[0]
OUT = W / "deposits" / "dt4" / "qbar.json"
ROOT_BOOT = 2026110641
B = 400
ALPHAS = (0.05, 0.01)


def os_ci(sorted_x, p, conf=0.95):
    """Exact order-statistic CI for the p-quantile (distribution-free, continuous data)."""
    n = len(sorted_x)
    lo_k = int(binom.ppf((1 - conf) / 2, n, p))
    hi_k = int(binom.ppf(1 - (1 - conf) / 2, n, p)) + 1
    return float(sorted_x[max(lo_k - 1, 0)]), float(sorted_x[min(hi_k - 1, n - 1)])


def fit(ns, q, se, model):
    n = np.asarray(ns, float)
    cols = [np.ones_like(n), n ** -0.5]
    if model == "M2":
        cols.append(np.log(n) * n ** -0.5)
    if model == "M3":
        cols = [np.ones_like(n), n ** (-1 / 3)]
    X = np.stack(cols, axis=1)
    w = 1 / np.asarray(se) ** 2
    beta, *_ = np.linalg.lstsq(X * np.sqrt(w)[:, None], np.asarray(q) * np.sqrt(w), rcond=None)
    resid = (np.asarray(q) - X @ beta) / np.asarray(se)
    return beta, float(np.sum(resid ** 2))


def main():
    ref = np.concatenate([np.load(dc.part_path(("ref", None, c))) for c in range(dc.REF_R // dc.CHUNK)])
    ns = list(dc.NS)
    scaled = np.sqrt(np.asarray(ns, float))[None, :] * ref  # rows = datasets, cols = n (nested)
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_BOOT)))
    out = {"ns": ns, "reps": int(ref.shape[0]), "bootstrap_B": B, "root_boot": ROOT_BOOT, "alphas": {}}
    for a in ALPHAS:
        p = 1 - a
        pts, cis = [], []
        for j in range(len(ns)):
            s = np.sort(scaled[:, j])
            pts.append(float(np.quantile(s, p, method="inverted_cdf")))
            cis.append(os_ci(s, p))
        boots = np.empty((B, len(ns)))
        for b in range(B):
            idx = rng.integers(0, ref.shape[0], ref.shape[0])
            boots[b] = np.quantile(scaled[idx], p, axis=0, method="inverted_cdf")
        se = boots.std(axis=0, ddof=1)
        models = {}
        for m in ("M1", "M2", "M3"):
            beta, chi2 = fit(ns, pts, se, m)
            qb = [fit(ns, boots[b], se, m)[0][0] for b in range(B)]
            models[m] = {"qbar": float(beta[0]), "coef": [float(v) for v in beta[1:]], "chi2": chi2,
                         "dof": len(ns) - len(beta), "qbar_boot_sd": float(np.std(qb, ddof=1))}
        qbs = [v["qbar"] for v in models.values()]
        lo = min(v["qbar"] - 2 * v["qbar_boot_sd"] for v in models.values())
        hi = max(v["qbar"] + 2 * v["qbar_boot_sd"] for v in models.values())
        lo = min(lo, pts[-1])  # the limit cannot be judged below the largest-n estimate without a monotonicity proof
        size = {}
        for name, qv in (("lo", lo), ("M1", models["M1"]["qbar"]), ("hi", hi), ("q_204800", pts[-1])):
            ecdf = np.sort(scaled[:, -1])
            h_top = np.searchsorted(ecdf, qv / np.sqrt(0.1), side="right") / len(ecdf)
            size[name] = {"q": qv, "sup_step_laws_pct": 100 * (1 - FM_pub(2 * qv)),
                          "W0.1_limit_pct": 100 * (1 - FM_pub(2 * qv / np.sqrt(0.9)) * h_top),
                          "W0.1_lower_bound_pct": 100 * (1 - FM_pub(2 * qv / np.sqrt(0.9)))}
        out["alphas"][str(a)] = {"q_point": pts, "q_os_ci95": cis, "q_boot_se": se.tolist(), "models": models,
                                 "qbar_models_range": [min(qbs), max(qbs)], "qbar_evidential_band": [lo, hi],
                                 "size": size}
        print(f"alpha={a}: q_n = " + ", ".join(f"{n}:{q:.4f}[{c[0]:.4f},{c[1]:.4f}]" for n, q, c in zip(ns, pts, cis, strict=True)))
        print("   boot SE", np.round(se, 5))
        for m, v in models.items():
            print(f"   {m}: qbar={v['qbar']:.4f} (boot sd {v['qbar_boot_sd']:.4f}) coef={np.round(v['coef'], 4)} chi2={v['chi2']:.2f} dof={v['dof']}")
        print(f"   evidential band [{lo:.4f}, {hi:.4f}]")
        for k, v in size.items():
            print(f"   {k:9s} q={v['q']:.4f}: sup={v['sup_step_laws_pct']:.2f}%  W(0.1) limit={v['W0.1_limit_pct']:.2f}%  (lower bd {v['W0.1_lower_bound_pct']:.2f}%)")
    OUT.write_text(json.dumps(out, indent=1))
    print("saved", OUT)


if __name__ == "__main__":
    main()
