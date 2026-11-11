"""Fidelity checks for code/ch98.py before the DT5b freeze. DEV, root 2026110684. No hypothesis test.

1. KDE: binned d_hat vs unbinned d_hat_direct on the same grid (relative difference), incl. the two-point case.
2. Inversion: reference(d) round trip on a fine d grid, incl. the band around 2 pi.
3. Table vs direct: dip quantile at a d between two nodes, computed directly (M samples) vs interpolated from the
   two neighbouring node quantiles; compared with the Monte Carlo standard error.
4. Curved laws vs the paper's Table 1 (alpha = 0.05; n = 50, 200): normal, t6, beta(4, 8), beta(12, 4); plus beta(2, 2).
Output: deposits/dt5b_fidelity.json
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch98  # noqa: E402
from dt0c_screen import dips_of  # noqa: E402

ROOT = 2026110684
W = Path(__file__).resolve().parents[1]
PAPER_T1 = {"normal": (0.036, 0.048), "t6": (0.032, 0.024), "beta(4,8)": (0.056, 0.042), "beta(12,4)": (0.062, 0.056)}


def rng(*key):
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT, spawn_key=key)))


def draw(g, law, size):
    if law == "normal":
        return g.standard_normal(size)
    if law == "t6":
        return g.standard_t(6, size)
    if law == "uniform":
        return g.random(size)
    if law.startswith("beta("):
        a, b = (float(v) for v in law[5:-1].split(","))
        return g.beta(a, b, size)
    if law == "W(0.1)":
        pile = g.random(size) < 0.1
        return np.where(pile, 0.995 + 0.005 * g.random(size), 0.995 * g.random(size))
    raise ValueError(law)


def node_quantiles(n, d_values, M, a, key):
    out = []
    for j, d in enumerate(d_values):
        g = rng(key, n, j)
        x = np.sort(ch98.sample_reference(g, ch98.reference(float(d)), (M, n)), axis=1)
        out.append(float(np.quantile(dips_of(x), 1 - a)))
    return np.array(out)


def main():
    t0 = time.time()
    res = {}
    # 1. KDE
    kde = {}
    for li, law in enumerate(["normal", "beta(4,8)", "uniform", "W(0.1)", "t6"]):
        for n in (50, 800):
            g = rng(1, li, n)
            rel = []
            for _ in range(100):
                x = draw(g, law, n)
                a, b = ch98.d_hat(x), ch98.d_hat_direct(x)
                rel.append(abs(a - b) / b)
            kde[f"{law}/n{n}"] = float(np.max(rel))
    x = np.tile([-1.0, 1.0], 25)
    kde["two-point"] = abs(ch98.d_hat(x) - ch98.d_hat_direct(x)) / ch98.d_hat_direct(x)
    res["kde_max_rel_diff"] = kde
    print("1 KDE max rel diff:", {k: round(v, 4) for k, v in kde.items()}, flush=True)
    # 2. inversion
    bad = []
    for d in np.concatenate([np.geomspace(1e-4, 1e4, 400), ch98.TWO_PI * (1 + np.linspace(-1e-2, 1e-2, 201))]):
        try:
            r = ch98.reference(float(d))
            if r[0] == "beta":
                err = abs(ch98.d_beta(r[1]) - d) / d
            elif r[0] == "t":
                err = abs(ch98.d_t(r[1]) - d) / d
            else:
                err = abs(ch98.TWO_PI - d) / d
            if err > 2 * ch98.NEAR_NORMAL:
                bad.append((float(d), r[0], err))
        except (ValueError, RuntimeError) as e:
            bad.append((float(d), "error", str(e)))
    res["inversion_bad"] = bad
    print("2 inversion failures:", len(bad), flush=True)
    # 3. table vs direct at between-node d
    n, M, a = 200, 20000, 0.05
    tv = []
    for k, d in enumerate((0.0005, 0.5, 3.0, 6.0, 9.0, 40.0)):
        j = 0 if d < ch98.D_GRID[0] else int(np.searchsorted(ch98.D_GRID, d))  # node index in D_NODES of the lower node
        nodes = ch98.D_NODES[j:j + 2]
        qn = node_quantiles(n, nodes, M, a, 30 + k)
        table = np.full(len(ch98.D_NODES), np.nan)
        table[j], table[j + 1] = qn
        q_interp = float(ch98.crit(table, np.array([d]))[0])
        q_direct = float(node_quantiles(n, [d], M, a, 40 + k)[0])
        tv.append({"d": d, "nodes": nodes.tolist(), "q_interp": q_interp, "q_direct": q_direct,
                   "rel_diff": abs(q_interp - q_direct) / q_direct})
    res["table_vs_direct_n200"] = tv
    print("3 table vs direct:", [(t["d"], round(t["rel_diff"], 4)) for t in tv], flush=True)
    # 4. curved laws, small table around the observed d_hat values
    rep = {}
    for n in (50, 200):
        tab = {a: node_quantiles(n, ch98.D_NODES, 5000, a, 50) for a in (0.05,)}
        for li, law in enumerate(["normal", "t6", "beta(4,8)", "beta(12,4)", "beta(2,2)", "uniform"]):
            g = rng(60, li, n)
            x = np.sort(draw(g, law, (2000, n)), axis=1)
            dip = dips_of(x)
            dh = np.array([ch98.d_hat(r) for r in x])
            rate = float(np.mean(dip > ch98.crit(tab[0.05], dh)))
            rep[f"{law}/n{n}"] = {"rate": rate, "paper": PAPER_T1.get(law, (None, None))[0 if n == 50 else 1],
                                  "dhat_median": float(np.median(dh)),
                                  "frac_below_dmin": float(np.mean(dh < ch98.D_GRID[0])),
                                  "frac_above_dmax": float(np.mean(dh > ch98.D_GRID[-1]))}
    res["curved_laws_5pct"] = rep
    for k, v in rep.items():
        print(f"4 {k:14s} CH98 {100 * v['rate']:5.2f}%  paper {v['paper']}  dhat~{v['dhat_median']:.2f}", flush=True)
    res["seconds"] = time.time() - t0
    (W / "deposits" / "dt5b_fidelity.json").write_text(json.dumps(res, indent=1))
    print(f"saved; {res['seconds']:.0f}s")


if __name__ == "__main__":
    main()
