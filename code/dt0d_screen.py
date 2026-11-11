"""Step 2c (DT0-d) follow-up screen. DEV tier only. Frozen design: work/DT0d-PREREG.md.

(i) local slope scaling beta_n = c / sqrt(n) on the witness W(0.1, 0.005);
(ii) fixed Gaussian noise on W up to n = 12800;
(iii) laws with a structurally flat shoulder: p-value mixtures, z-test p-values, biased-forecast PIT.
Seeds per CV-6: root 2026110550, reference spawn (0, n, chunk), laws spawn (1, f, n).

Run:  python code/dt0d_screen.py [--quick]
"""

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.special import ndtr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import bracket, cp_lower, cp_upper, dips_of, g_uniform, g_witness  # noqa: E402, I001

ROOT = 2026110550
NS = (200, 800, 3200, 12800)
ALPHAS = (0.05, 0.01)
R_REF = 1_000_000
REF_CHUNK = 50_000
REPS = 20_000
GAMMA = 0.001
OUT = Path(__file__).resolve().parents[1] / "deposits" / "dt0d"


def g_pmix(pi0, a):
    def g(rng, size):
        null = rng.random(size) < pi0
        return np.where(null, rng.random(size), rng.random(size) ** (1 / a))
    return g


def g_zp(pi0, mu):
    def g(rng, size):
        null = rng.random(size) < pi0
        return np.where(null, rng.random(size), ndtr(-(rng.standard_normal(size) + mu)))
    return g


def g_pit(delta):
    def g(rng, size):
        return ndtr(rng.standard_normal(size) - delta)
    return g


C_LOCAL = (0, 2, 4, 8)
LAWS = {  # f: (name, block)
    0: ("W c=0", "i"), 1: ("W c=2", "i"), 2: ("W c=4", "i"), 3: ("W c=8", "i"),
    4: ("W beta=0.1 fixed", "i"), 5: ("W beta=0.25 fixed", "i"),
    6: ("W+N(0,.005)", "ii"), 7: ("W+N(0,.025)", "ii"), 8: ("W+N(0,.05)", "ii"),
    9: ("pmix pi0=.9 a=.1", "iii"), 10: ("pmix pi0=.9 a=.25", "iii"),
    11: ("pmix pi0=.8 a=.1", "iii"), 12: ("pmix pi0=.8 a=.25", "iii"),
    13: ("z-p pi0=.9 mu=2", "iii"), 14: ("z-p pi0=.9 mu=3", "iii"),
    15: ("PIT delta=.25", "iii"), 16: ("PIT delta=.5", "iii"),
    17: ("U(0,1)", "anchor"),
}


def law(f, n):
    if f <= 3:
        return g_witness(beta=C_LOCAL[f] / np.sqrt(n))
    return {
        4: g_witness(beta=0.1), 5: g_witness(beta=0.25),
        6: g_witness(noise_sd=0.005), 7: g_witness(noise_sd=0.025), 8: g_witness(noise_sd=0.05),
        9: g_pmix(0.9, 0.1), 10: g_pmix(0.9, 0.25), 11: g_pmix(0.8, 0.1), 12: g_pmix(0.8, 0.25),
        13: g_zp(0.9, 2.0), 14: g_zp(0.9, 3.0), 15: g_pit(0.25), 16: g_pit(0.5), 17: g_uniform,
    }[f]


def run_cell(args):
    key, f, n, reps = args
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT, spawn_key=key)))
    gen = g_uniform if f == "ref" else law(f, n)
    batch = max(100, min(2000, 2_000_000 // n))
    out = np.empty(reps)
    for b0 in range(0, reps, batch):
        m = min(batch, reps - b0)
        out[b0:b0 + m] = dips_of(gen(rng, (m, n)))
    return key, out


def main():
    quick = "--quick" in sys.argv
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    r_ref, reps, chunk = (10_000, 200, 10_000) if quick else (R_REF, REPS, REF_CHUNK)
    ns = (200, 800) if quick else NS
    jobs = [((0, n, c), "ref", n, chunk) for n in ns for c in range(r_ref // chunk)]
    jobs += [((1, f, n), f, n, reps) for f in LAWS for n in ns]
    jobs.sort(key=lambda j: -j[2] * j[3])
    res = {}
    with ProcessPoolExecutor(max_workers=11) as ex:
        for key, d in ex.map(run_cell, jobs, chunksize=1):
            res[key] = d
    print(f"simulated {len(jobs)} jobs in {time.time() - t0:.0f}s", flush=True)

    g_tail = GAMMA / (len(ns) * len(ALPHAS)) / 2
    summary = {"tier": "DEV", "root": ROOT, "R_ref": r_ref, "reps": reps, "quick": quick,
               "reference": {}, "cells": []}
    crit = {}
    for n in ns:
        ref = np.sort(np.concatenate([res[(0, n, c)] for c in range(r_ref // chunk)]))
        if not quick:
            np.save(OUT / f"ref_n{n}.npy", ref)
        for a in ALPHAS:
            k = int(np.ceil((1 - a) * r_ref))
            lo, hi = bracket(r_ref, 1 - a, g_tail)
            crit[(n, a)] = (ref[k - 1], ref[lo - 1], ref[hi - 1])
            summary["reference"][f"n{n}_a{a}"] = {
                "q_hat": float(ref[k - 1]), "X_l": float(ref[lo - 1]), "X_u": float(ref[hi - 1]),
                "sqrt_n_q": float(np.sqrt(n) * ref[k - 1]), "ties_at_q": int(np.sum(ref == ref[k - 1])) - 1}
    for f, (name, block) in LAWS.items():
        for n in ns:
            d = res[(1, f, n)]
            if not quick:
                np.save(OUT / f"dips_f{f}_n{n}.npy", d)
            sq = np.sqrt(n) * d
            for a in ALPHAS:
                q, xl, xu = crit[(n, a)]
                kq, kl, ku = int(np.sum(d >= q)), int(np.sum(d >= xl)), int(np.sum(d >= xu))
                summary["cells"].append({
                    "f": f, "name": name, "block": block, "n": n, "alpha": a,
                    "rej_point": kq / reps, "k_at_q": kq, "k_at_Xl": kl, "k_at_Xu": ku,
                    "cp95_lower_at_Xu": cp_lower(ku, reps, 0.975), "cp95_upper_at_Xl": cp_upper(kl, reps, 0.975),
                    "bonf_lower_at_Xu": cp_lower(ku, reps, 1 - 0.05 / 32),
                    "bonf_upper_at_Xl": cp_upper(kl, reps, 1 - 0.05 / 32),
                    "sqrt_n_dip_q50_q90_q95": [float(x) for x in np.quantile(sq, (0.5, 0.9, 0.95))],
                    "beta_n": (C_LOCAL[f] / np.sqrt(n)) if f <= 3 else None})
    cells = {(c["f"], c["n"], c["alpha"]): c for c in summary["cells"]}

    # readings (PREREG section 4), alpha = 0.05
    def r(f, n):
        return cells[(f, n, 0.05)]["rej_point"]
    uni = [(n, a, cells[(17, n, a)]["rej_point"],
            bool(cp_lower(cells[(17, n, a)]["k_at_q"], reps, 1 - 0.05 / 16) <= a
                 <= cp_upper(cells[(17, n, a)]["k_at_q"], reps, 1 - 0.05 / 16))) for n in ns for a in ALPHAS]
    valid = all(u[3] for u in uni) and r(0, 200) >= 0.15
    rd = {"anchors": {"uniform": uni, "W_n200": r(0, 200), "valid": valid}}
    if not quick:
        a_ok = all(abs(r(f, 3200) - r(f, 12800)) <= 0.025 for f in range(4))
        b_ok = all(r(f, 12800) >= r(f + 1, 12800) for f in range(3))
        c_ok = r(4, 800) - r(4, 12800) > 0.025
        rd["i_local"] = {"a_stable": a_ok, "b_monotone_c": b_ok, "c_fixed_beta_falls": c_ok,
                         "verdict": "SUPPORTED" if (a_ok and b_ok and c_ok) else "NOT SUPPORTED"}
        rd["ii_noise"] = {}
        for f in (6, 7, 8):
            up = cells[(f, 12800, 0.05)]["cp95_upper_at_Xl"]
            lo = cells[(f, 3200, 0.05)]["cp95_lower_at_Xu"]
            rd["ii_noise"][LAWS[f][0]] = "TURNOVER" if up < lo else "no turnover to n=12800"
    rd["iii"] = {}
    for f in range(9, 17):
        cs = [cells[(f, n, 0.05)] for n in ns]
        st = ("MATERIAL" if any(c["bonf_lower_at_Xu"] >= 0.075 for c in cs) else
              "CONTROLLED" if all(c["bonf_upper_at_Xl"] < 0.075 for c in cs) else "INCONCLUSIVE")
        rd["iii"][LAWS[f][0]] = st
    summary["readings"] = rd
    summary["seconds"] = time.time() - t0
    name = "screen_quick.json" if quick else "screen.json"
    (OUT / name).write_text(json.dumps(summary, indent=1))
    print(json.dumps(rd, indent=1))
    for f, (nm, b) in LAWS.items():
        print(f, nm.ljust(20), b.ljust(6), " ".join(f"{r(f, n):.4f}" for n in ns))
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
