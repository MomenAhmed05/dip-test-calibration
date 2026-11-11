"""Tool smoke test for O5-ML (DEV, not evidence): ClustPy UniDip / SkinnyDip / DipMeans with the uniform
('table') calibration vs the LF calibration, swapped by monkeypatching clustpy's dip_pval inside each module.
Run with the isolated venv:  python code/o5_smoke.py
Seeds: DEV root 2026110662 (DT5-ML-DESIGN §1).
"""
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fm_table import FM_pub  # noqa: E402

import clustpy.centroid.dipmeans as cdm  # noqa: E402
import clustpy.density.skinnydip as csd  # noqa: E402
from clustpy.centroid import DipMeans  # noqa: E402
from clustpy.density import SkinnyDip, UniDip  # noqa: E402
from clustpy.utils import dip_pval as dip_pval_uniform  # noqa: E402


def dip_pval_lf(dip_value, n_points, *args, **kwargs):
    """LF p-value 1 - F_M(2 sqrt(n) dip); F_M only from the published table (clip outside 0.33-2.54)."""
    return float(1.0 - FM_pub(2.0 * np.sqrt(n_points) * dip_value, clip_tails=True))


def set_calibration(name):
    f = dip_pval_uniform if name == "uniform" else dip_pval_lf
    cdm.dip_pval = f
    csd.dip_pval = f


def w_law(rng, n, w=0.1, eps=0.005):
    pile = rng.random(n) < w
    return np.where(pile, 1 - eps + eps * rng.random(n), (1 - eps) * rng.random(n))


def main():
    rng = np.random.default_rng(np.random.SeedSequence(2026110662))
    R, n = 20, 500
    for cal in ("uniform", "LF"):
        set_calibration(cal)
        t0 = time.time()
        k1 = []
        for r in range(R):
            x = w_law(np.random.default_rng([2026110662, r]), n)
            k1.append(UniDip().fit(x.reshape(-1, 1)).n_clusters_)
        t1 = time.time()
        k5s, k5m = [], []
        for r in range(R // 4):
            g = np.random.default_rng([2026110662, 1000 + r])
            X = np.column_stack([w_law(g, n)] + [g.standard_normal(n) for _ in range(4)])
            k5s.append(SkinnyDip(random_state=r).fit(X).n_clusters_)
            k5m.append(DipMeans(random_state=r).fit(X).n_clusters_)
        t2 = time.time()
        print(f"{cal:8s} UniDip 1-D W(0.1) n={n}: P(k>1) = {np.mean(np.array(k1) > 1):.2f} ({R} reps, {t1 - t0:.1f}s)"
              f" | d=5 SkinnyDip P(k>1) = {np.mean(np.array(k5s) > 1):.2f}, DipMeans P(k>1) = {np.mean(np.array(k5m) > 1):.2f}"
              f" ({R // 4} reps, {t2 - t1:.1f}s)")
    _ = rng


if __name__ == "__main__":
    main()
