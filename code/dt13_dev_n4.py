"""Step 14 (DT13) DEV checks, root 2026110740 (declared DEV in CONVENTIONS CV-6 before use).

1. H&H 1985 n = 4 reduction: T = max(1, s2 / max(s1, s3)), s_i = spacings of the 4 order statistics.
   Check: the dip D_alg is a function of T alone, nondecreasing; fit its closed form.
2. P_U(T >= r) = 2 / ((r + 1)(r + 2)) (derived by S from iid exponential spacings): MC check.
3. Decoupled pinned form: D = min over u of max(V_[0,u], C_[u,1]) ?  (pinned_dip.Z_dec), n = 5..8, 20.
Usage: python dt13_dev_n4.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt0c_screen import dips_of  # noqa: E402
from pinned_dip import Z_dec  # noqa: E402

ROOT = 2026110740


def dip(x):
    return float(dips_of(np.sort(x)[None, :])[0])


def main():
    rng = np.random.default_rng(np.random.SeedSequence([ROOT, 1]))
    # 1. D as a function of T at n = 4
    X = np.sort(rng.random((20000, 4)), axis=1)
    s = np.diff(X, axis=1)
    T = np.maximum(1.0, s[:, 1] / np.maximum(s[:, 0], s[:, 2]))
    D = np.array([dip(x) for x in X])
    o = np.argsort(T)
    dmono = float(np.min(np.diff(D[o])))
    print(f"n=4: min successive diff of D sorted by T (>= -1e-12 if nondecreasing): {dmono:.3e}")
    cand = {"T/(4(T+1))": T / (4 * (T + 1)), "(T+1)/(4(T+2))": (T + 1) / (4 * (T + 2)),
            "T/(2(2T+2))": T / (2 * (2 * T + 2)), "(2T+1)/(8(T+1))": (2 * T + 1) / (8 * (T + 1))}
    for k, v in cand.items():
        print(f"  max|D - {k}| = {np.max(np.abs(D - v)):.3e}")
    print(f"  P(D = 1/8 within 1e-12) = {np.mean(np.abs(D - 0.125) < 1e-12):.4f} (pred 2/3); P(T>1)={np.mean(T > 1):.4f}")
    # 2. tail law of T
    X = np.sort(rng.random((2_000_000, 4)), axis=1)
    s = np.diff(X, axis=1)
    T = np.maximum(1.0, s[:, 1] / np.maximum(s[:, 0], s[:, 2]))
    for r in (1.0, 1.5, 3.0, (np.sqrt(161) - 3) / 2, (np.sqrt(801) - 3) / 2):
        p = np.mean(T > 1) if r <= 1 else np.mean(T >= r)
        pr = 2 / ((r + 1) * (r + 2))
        se = np.sqrt(pr * (1 - pr) / len(T))
        print(f"  r={r:.4f}: MC {p:.5f} formula {pr:.5f} z={(p - pr) / se:+.2f}")
    # 3. decoupled pinned form vs diptest
    for n in (5, 6, 7, 8, 20):
        worst, worst_signed = 0.0, 0.0
        for _ in range(300):
            x = np.sort(rng.random(n))
            us = np.unique(np.concatenate([np.linspace(0, 1, 201), x, (x[:-1] + x[1:]) / 2,
                                           np.clip(x - 1e-12, 0, 1), np.clip(x + 1e-12, 0, 1)]))
            z = min(Z_dec(x, u) for u in us)
            d = dip(x)
            worst = max(worst, abs(z - d))
            worst_signed = min(worst_signed, z - d)
        print(f"n={n}: max|min_u Z_dec - D| over 300 samples = {worst:.3e} (most negative z-d {worst_signed:.3e})")


if __name__ == "__main__":
    main()
