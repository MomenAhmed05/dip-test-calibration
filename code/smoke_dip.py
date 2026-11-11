"""Step 2 smoke: D_alg (diptest 0.12) vs D_U / D_cont (code/lp_dip.py), deterministic fixtures,
invariances, ties, and timing. Writes deposits/dt0b/smoke_dip.json. Seeds: CV-6 root SMOKE.
Run: python code/smoke_dip.py   (plain python: diptest is a --user install)
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import diptest

sys.path.insert(0, str(Path(__file__).parent))
from lp_dip import dip_lp  # noqa: E402

ROOT = 2026110530  # CV-6: Step 2 smoke root (DEV, banned for confirmation)
out = {"diptest_version": diptest.__version__, "root": ROOT}
def d_alg(x):
    return float(diptest.dipstat(np.asarray(x, float)))  # type: ignore[arg-type]


# 1. deterministic fixtures
fx = {}
fx["n4_equispaced"] = [0, 1, 2, 3]
fx["all_identical_n10"] = [1.0] * 10
fx["two_values_5_5"] = [0.0] * 5 + [1.0] * 5
fx["two_values_2_8"] = [0.0] * 2 + [1.0] * 8
fx["n4_atom_end_HH"] = [0.0, 1 / 3, 2 / 3, 1.0, 1.0]  # uniform-ish + duplicate at the top
rng = np.random.default_rng(np.random.SeedSequence(ROOT, spawn_key=(1,)))
base = np.sort(rng.random(30))
fx["cont30_plus_endblock8"] = list(base) + [1.0] * 8
fx["cont30_plus_interiorblock8"] = list(base) + [0.5] * 8
fx["cont30"] = list(base)
fixtures = {}
for name, x in fx.items():
    fixtures[name] = {"n": len(x), "D_alg": d_alg(x), "D_U": dip_lp(x, True), "D_cont": dip_lp(x, False)}
out["fixtures"] = fixtures

# 2. invariances (on cont30): reflection, positive affine
x = np.array(fx["cont30"])
out["invariance"] = {
    "D_alg": d_alg(x), "D_alg_reflect": d_alg(-x), "D_alg_affine": d_alg(3.7 * x - 11.0),
    "D_U": dip_lp(x), "D_U_reflect": dip_lp(-x),
}

# 3. near-ties vs exact ties
xt = np.array(list(base) + [1.0] * 8)
xn = np.array(list(base) + list(1.0 + 1e-9 * np.arange(8)))
out["near_vs_exact_tie"] = {"exact_D_alg": d_alg(xt), "near_D_alg": d_alg(xn),
                            "exact_D_U": dip_lp(xt), "near_D_U": dip_lp(xn)}

# 4. random continuous samples: LP vs diptest, 1,000 total over n in {10, 30, 50}
agree = {}
for n in (10, 30, 50):
    g = np.random.default_rng(np.random.SeedSequence(ROOT, spawn_key=(2, n)))
    m = 400 if n == 10 else 300
    diffs = []
    for _ in range(m):
        # mixture of shapes so that both small and large dips occur
        u = g.random(n)
        sel = g.integers(0, 3)
        y = u if sel == 0 else (np.where(g.random(n) < 0.5, 0.3 * u, 0.7 + 0.3 * u) if sel == 1 else g.normal(size=n))
        diffs.append(d_alg(y) - dip_lp(y))
    diffs = np.array(diffs)
    agree[n] = {"m": m, "max_abs_diff": float(np.abs(diffs).max()), "n_over_1e-9": int((np.abs(diffs) > 1e-9).sum())}
out["lp_vs_diptest_continuous"] = agree

# 5. tied samples: atom-at-top laws (the atom rows) and rounded data
tied = {}
g = np.random.default_rng(np.random.SeedSequence(ROOT, spawn_key=(3,)))
for kind in ("top_atom", "interior_atom", "rounded"):
    diffs, dcont = [], []
    for _ in range(100):
        n = 30
        u = g.random(n)
        if kind == "top_atom":
            y = np.where(g.random(n) < 0.3, 1.0, u)
        elif kind == "interior_atom":
            y = np.where(g.random(n) < 0.3, 0.5, u)
        else:
            y = np.round(g.normal(size=n) * 2) / 2
        a, b, c = d_alg(y), dip_lp(y, True), dip_lp(y, False)
        diffs.append(a - b)
        dcont.append(c - b)
    diffs = np.array(diffs)
    tied[kind] = {"m": 100, "max_abs(D_alg-D_U)": float(np.abs(diffs).max()),
                  "n_over_1e-9": int((np.abs(diffs) > 1e-9).sum()),
                  "mean(D_cont-D_U)": float(np.mean(dcont))}
out["lp_vs_diptest_tied"] = tied

# 6. timing of diptest
tm = {}
for n in (100, 200, 1000, 10000):
    g = np.random.default_rng(np.random.SeedSequence(ROOT, spawn_key=(4, n)))
    X = g.random((2000, n))
    t0 = time.perf_counter()
    for row in X:
        diptest.dipstat(np.sort(row), sort_x=False)
    tm[n] = (time.perf_counter() - t0) / 2000 * 1e6
out["diptest_us_per_dip_presorted"] = tm

Path("deposits/dt0b").mkdir(parents=True, exist_ok=True)
Path("deposits/dt0b/smoke_dip.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
