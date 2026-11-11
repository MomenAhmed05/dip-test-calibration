"""O5-ML DipMeans diagnostic: is the distance law H_x = law of ||X - x|| unimodal for the d = 5 null laws?
X = (law coordinate, 4 independent N(0,1)); viewers x = (law quantile 0.1 / 0.5 / 0.9, 0, 0, 0, 0).
H_x is judged from a LARGE independent reference sample (N = 1e6, not the tested data): sqrt(N) * dip and the
LF-asymptotic p-value; p_LF < 1e-6 is read as "H_x not unimodal" (diagnostic, no claim tier).
Run with the venv python: python code/o5_hx.py      DEV root 2026110664.
"""
import json
import sys
from pathlib import Path

import diptest
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from fm_table import FM_pub  # noqa: E402
from o5_ml import LAWS, sample_1d  # noqa: E402

N = 1_000_000


def main():
    out = []
    for li, (name, role, (kind, par)) in enumerate(LAWS):
        if role == "price":
            continue
        g = np.random.default_rng(np.random.SeedSequence([2026110664, li]))
        n_law = 500  # pinning laws depend on n: use the DipMeans n
        x1 = sample_1d(kind, par, N, g) if kind != "pin" else _pin(g, N, par / np.sqrt(n_law))
        Z = g.standard_normal((N, 4))
        for qv in (0.1, 0.5, 0.9):
            v = float(np.quantile(x1, qv))
            d = np.sqrt((x1 - v) ** 2 + np.sum(Z ** 2, axis=1))
            dip = float(diptest.dipstat(np.sort(d), sort_x=False))  # type: ignore[arg-type]
            t = 2 * np.sqrt(N) * dip
            p = float(1 - FM_pub(t, clip_tails=True))
            out.append({"law": name, "viewer_q": qv, "sqrtN_dip": np.sqrt(N) * dip, "p_LF": p,
                        "H_x_unimodal": "no" if p < 1e-6 else "no evidence against"})
            print(f"{name:12s} viewer q{qv}: sqrt(N)*dip = {np.sqrt(N) * dip:.3f}  p_LF = {p:.2g}")
    (HERE.parent / "deposits" / "o5ml_hx.json").write_text(json.dumps(out, indent=1))


def _pin(g, n, w):
    pile = g.random(n) < w
    return np.where(pile, 0.995 + 0.005 * g.random(n), 0.995 * g.random(n))


if __name__ == "__main__":
    main()
