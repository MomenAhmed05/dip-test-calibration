"""Step 2b (DT0-c) early breadth screen. DEV tier only. Frozen design: work/DT0c-PREREG.md.

Statistic D_alg (CONVENTIONS CV-1), our own uniform reference H_n (CV-3), rejection
iff D >= c (CV-4), seeds per CV-6: root 2026110540, reference spawn (0, n), anchors and
families spawn (1, f, n). Writes every dip to deposits/dt0c/ and a summary JSON.

Run:  python -I code/dt0c_screen.py            (full screen)
      python -I code/dt0c_screen.py --quick    (tiny smoke of every generator)
"""

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import diptest
import numpy as np
from scipy import stats

ROOT = 2026110540
NS = (100, 200, 500, 1000)
ALPHAS = (0.05, 0.01)
R_REF = 1_000_000
REPS = 20_000
BATCH = 2_000
GAMMA = 0.001  # bracket failure, spread over the 8 (n, alpha) pairs
OUT = Path(__file__).resolve().parents[1] / "deposits" / "dt0c"

IQR_Z = 1.3489795
LN3 = float(np.log(3.0))


def iqr_beta(a):
    return 0.75 ** (1 / a) - 0.25 ** (1 / a)


# ---------------------------------------------------------------- generators
# Each generator takes (rng, size) with size = (rows, n) and returns float64 draws.


def g_uniform(rng, size):
    return rng.random(size)


def g_cens_normal(c, noise=0.0):
    def g(rng, size):
        y = np.minimum(rng.standard_normal(size), c)
        if noise:
            y = y + noise * IQR_Z * rng.standard_normal(size)
        return y
    return g


def g_topbeta(a, noise=0.0):
    q = 0.9 ** (1 / a)

    def g(rng, size):
        y = np.minimum(rng.random(size) ** (1 / a), q)
        if noise:
            y = y + noise * iqr_beta(a) * rng.standard_normal(size)
        return y
    return g


def g_zeroexp(noise=0.0):
    def g(rng, size):
        atom = rng.random(size) < 0.1
        y = np.where(atom, 0.0, rng.standard_exponential(size))
        if noise:
            y = y + noise * LN3 * rng.standard_normal(size)
        return y
    return g


def g_spread_normal(delta):
    def g(rng, size):
        z = rng.standard_normal(size)
        u = rng.random(size)
        return np.where(z < 0, z, u * delta)
    return g


def g_spread_beta(a, delta):
    q = 0.9 ** (1 / a)

    def g(rng, size):
        b = rng.random(size) ** (1 / a)
        u = rng.random(size)
        return np.where(b < q, b, q + u * delta)
    return g


def g_spread_exp(delta):
    def g(rng, size):
        atom = rng.random(size) < 0.1
        e = rng.standard_exponential(size)
        u = rng.random(size)
        return np.where(atom, u * delta, e)
    return g


def g_unif_cap(cap, noise_sd=0.0):
    def g(rng, size):
        y = np.minimum(rng.random(size), cap)
        if noise_sd:
            y = y + noise_sd * rng.standard_normal(size)
        return y
    return g


def g_relu(b, noise=0.0):
    def g(rng, size):
        y = np.maximum(0.0, rng.standard_normal(size) + b)
        if noise:
            y = y + noise * IQR_Z * rng.standard_normal(size)
        return y
    return g


def g_witness(w=0.1, eps=0.005, beta=0.0, noise_sd=0.0):
    L = 1 - eps

    def g(rng, size):
        top = rng.random(size) < w
        u = rng.random(size)
        if beta == 0:
            t = u
        else:  # shoulder density prop. to 1 + beta (t - 1/2) on t = x / L in (0, 1)
            b1 = 1 - beta / 2
            t = (-b1 + np.sqrt(b1 * b1 + 2 * beta * u)) / beta
        y = np.where(top, L + eps * rng.random(size), L * t)
        if noise_sd:
            y = y + noise_sd * rng.standard_normal(size)
        return y
    return g


def _softplus(x):
    return np.logaddexp(0.0, x)


def g_logistic_top(tau, w=0.1, eps=0.005):
    L = 1 - eps
    s0 = float(_softplus(-L / tau))
    ctot = tau * (float(_softplus((1 - L) / tau)) - s0)

    def g(rng, size):
        top = rng.random(size) < w
        v = rng.random(size) * ctot
        s = v / tau + s0  # softplus((x - L) / tau) = s  ->  (x - L) / tau = s + log(1 - e^{-s})
        xt = L + tau * (s + np.log(-np.expm1(-s)))
        xt = np.clip(xt, 0.0, 1.0)
        return np.where(top, xt, rng.random(size))
    return g


def g_ramp(w=0.1, eps_r=0.01):
    Lp = 1 - eps_r
    d0 = (1 - w) / Lp
    ctop = 2 * w / eps_r - d0
    k = (ctop - d0) / eps_r  # density slope on the ramp

    def g(rng, size):
        top = rng.random(size) < w
        v = rng.random(size) * w  # mass inside ramp: d0 s + k s^2 / 2 = v
        s = (-d0 + np.sqrt(d0 * d0 + 2 * k * v)) / k
        return np.where(top, Lp + s, Lp * rng.random(size))
    return g


def g_offedge(m0=0.8, w=0.1, eps=0.005):
    def g(rng, size):
        top = rng.random(size) < w
        y = rng.random(size) * (1 - eps)
        y = np.where(y >= m0 - eps, y + eps, y)
        return np.where(top, m0 - eps + eps * rng.random(size), y)
    return g


def g_t5_clip(rng, size):
    return np.clip(rng.standard_t(5, size), -2.0, 2.0)


def g_round_normal(rng, size):
    return np.round(rng.standard_normal(size) * 2) / 2


FAMILIES = {
    0: ("U(0,1)", "anchor", g_uniform),
    1: ("min(Z,0)+1%", "decision F1", g_cens_normal(0.0, 0.01)),
    2: ("min(Z,0)+5%", "decision F1", g_cens_normal(0.0, 0.05)),
    3: ("min(B,q) a1.2 +1%", "decision F2", g_topbeta(1.2, 0.01)),
    4: ("min(B,q) a1.2 +5%", "decision F2", g_topbeta(1.2, 0.05)),
    5: ("0.1d0+0.9Exp +1%", "decision F3", g_zeroexp(0.01)),
    6: ("0.1d0+0.9Exp +5%", "decision F3", g_zeroexp(0.05)),
    7: ("min(Z,0)", "atom row", g_cens_normal(0.0)),
    8: ("min(B,q) a1.2", "atom row", g_topbeta(1.2)),
    9: ("0.1d0+0.9Exp", "atom row", g_zeroexp()),
    10: ("min(Z,0)+10%", "descriptive", g_cens_normal(0.0, 0.10)),
    11: ("min(B,q) a1.2 +10%", "descriptive", g_topbeta(1.2, 0.10)),
    12: ("0.1d0+0.9Exp +10%", "descriptive", g_zeroexp(0.10)),
    13: ("F1 spread U[0,1%]", "descriptive twin", g_spread_normal(0.01 * IQR_Z)),
    14: ("F1 spread U[0,5%]", "descriptive twin", g_spread_normal(0.05 * IQR_Z)),
    15: ("F2 spread 1%", "descriptive twin", g_spread_beta(1.2, 0.01 * iqr_beta(1.2))),
    16: ("F2 spread 5%", "descriptive twin", g_spread_beta(1.2, 0.05 * iqr_beta(1.2))),
    17: ("F3 spread 1%", "descriptive twin", g_spread_exp(0.01 * LN3)),
    18: ("F3 spread 5%", "descriptive twin", g_spread_exp(0.05 * LN3)),
    19: ("min(Z,-0.5)", "atom row", g_cens_normal(-0.5)),
    20: ("min(Z,-0.5)+5%", "descriptive", g_cens_normal(-0.5, 0.05)),
    21: ("min(Z,-1)", "atom row", g_cens_normal(-1.0)),
    22: ("min(Z,-1)+5%", "descriptive", g_cens_normal(-1.0, 0.05)),
    23: ("min(U,0.9)", "atom row", g_unif_cap(0.9)),
    24: ("min(U,0.9)+5%", "descriptive", g_unif_cap(0.9, 0.05 * 0.5)),
    25: ("min(B,q) a1.5", "atom row", g_topbeta(1.5)),
    26: ("min(B,q) a1.5 +5%", "descriptive", g_topbeta(1.5, 0.05)),
    27: ("max(0,Z-0.5)", "atom row", g_relu(-0.5)),
    28: ("max(0,Z-0.5)+5%", "descriptive", g_relu(-0.5, 0.05)),
    29: ("W(0.1,0.005)", "witness anchor", g_witness()),
    30: ("W beta0.1", "K10", g_witness(beta=0.1)),
    31: ("W beta0.25", "K10", g_witness(beta=0.25)),
    32: ("W beta0.5", "K10", g_witness(beta=0.5)),
    33: ("W beta1", "K10", g_witness(beta=1.0)),
    34: ("W beta2", "K10", g_witness(beta=2.0)),
    35: ("W logistic tau.003", "K10", g_logistic_top(0.003)),
    36: ("W logistic tau.01", "K10", g_logistic_top(0.01)),
    37: ("W logistic tau.03", "K10", g_logistic_top(0.03)),
    38: ("ramp w.1 eps.01", "K10", g_ramp()),
    39: ("off-edge m0.8", "K10", g_offedge()),
    40: ("W w.02", "K10", g_witness(w=0.02)),
    41: ("W w.05", "K10", g_witness(w=0.05)),
    42: ("W w.2", "K10", g_witness(w=0.2)),
    43: ("W w.4", "K10", g_witness(w=0.4)),
    44: ("t5 clip[-2,2]", "non-null law", g_t5_clip),
    45: ("max(0,Z+0.5)", "non-null law", g_relu(0.5)),
    46: ("N(0,1) round .5", "non-null law", g_round_normal),
    47: ("W+N(0,.005)", "descriptive", g_witness(noise_sd=0.005)),
    48: ("W+N(0,.025)", "descriptive", g_witness(noise_sd=0.025)),
}


# ---------------------------------------------------------------- dips


def dips_of(x):
    x = np.sort(x, axis=1)
    out = np.empty(x.shape[0])
    for i in range(x.shape[0]):
        row = x[i]
        if row[0] == row[-1]:
            out[i] = 0.0  # CV-1: all-identical sample
        else:
            out[i] = float(diptest.dipstat(row, allow_zero=False, sort_x=False))  # type: ignore[arg-type]
    return out


def run_cell(args):
    key, gen_id, n, reps = args
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT, spawn_key=key)))
    gen = g_uniform if gen_id == "ref" else FAMILIES[gen_id][2]
    out = np.empty(reps)
    for b0 in range(0, reps, BATCH):
        m = min(BATCH, reps - b0)
        out[b0:b0 + m] = dips_of(gen(rng, (m, n)))
    return key, out


# ---------------------------------------------------------------- statistics


def bracket(R, p, g_tail):
    """Order-stat indices (1-based) l <= u with P(X_(l) > q) <= g_tail and P(X_(u) < q) <= g_tail."""
    lo = int(stats.binom.ppf(g_tail, R, p))  # P(Bin < lo) <= g_tail
    while stats.binom.cdf(lo - 1, R, p) > g_tail:
        lo -= 1
    u = int(stats.binom.isf(g_tail, R, p)) + 1  # P(Bin >= u) <= g_tail
    while stats.binom.sf(u - 1, R, p) > g_tail:
        u += 1
    return max(lo, 1), min(u, R)


def cp_lower(k, N, conf):
    return 0.0 if k == 0 else float(stats.beta.ppf(1 - conf, k, N - k + 1))


def cp_upper(k, N, conf):
    return 1.0 if k == N else float(stats.beta.ppf(conf, k + 1, N - k))


def main():
    quick = "--quick" in sys.argv
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    r_ref, reps = (20_000, 200) if quick else (R_REF, REPS)
    ref_chunk = 50_000 if not quick else r_ref
    # reference: 20 child streams of SeedSequence(ROOT, (0, n)) for parallel speed, concatenated in order
    jobs = []
    for n in NS:
        for c in range(r_ref // ref_chunk):
            jobs.append(((0, n, c), "ref", n, ref_chunk))
    for f in FAMILIES:
        for n in NS:
            jobs.append(((1, f, n), f, n, reps))
    jobs.sort(key=lambda j: -j[2] * j[3])  # big first
    res = {}
    with ProcessPoolExecutor(max_workers=11) as ex:
        for key, d in ex.map(run_cell, jobs, chunksize=1):
            res[key] = d
    print(f"simulated {len(jobs)} jobs in {time.time() - t0:.0f}s", flush=True)

    summary = {"tier": "DEV", "root": ROOT, "R_ref": r_ref, "reps": reps, "gamma": GAMMA,
               "quick": quick, "reference": {}, "cells": []}
    crit = {}
    g_tail = GAMMA / (len(NS) * len(ALPHAS)) / 2
    for n in NS:
        ref = np.concatenate([res[(0, n, c)] for c in range(r_ref // ref_chunk)])
        ref.sort()
        if not quick:
            np.save(OUT / f"ref_n{n}.npy", ref)
        for a in ALPHAS:
            k = int(np.ceil((1 - a) * r_ref))
            q = ref[k - 1]
            lo, hi = bracket(r_ref, 1 - a, g_tail)
            crit[(n, a)] = (q, ref[lo - 1], ref[hi - 1])
            summary["reference"][f"n{n}_a{a}"] = {
                "q_hat": float(q), "X_l": float(ref[lo - 1]), "X_u": float(ref[hi - 1]),
                "l": lo, "u": hi, "ties_at_q": int(np.sum(ref == q)) - 1,
                "sqrt_n_q": float(np.sqrt(n) * q)}

    for f, (name, role, _) in FAMILIES.items():
        for n in NS:
            d = res[(1, f, n)]
            if not quick:
                np.save(OUT / f"dips_f{f}_n{n}.npy", d)
            for a in ALPHAS:
                q, xl, xu = crit[(n, a)]
                kq, kl, ku = int(np.sum(d >= q)), int(np.sum(d >= xl)), int(np.sum(d >= xu))
                summary["cells"].append({
                    "f": f, "name": name, "role": role, "n": n, "alpha": a,
                    "rej_point": kq / reps, "k_at_q": kq, "k_at_Xl": kl, "k_at_Xu": ku,
                    "cp95_lower_at_Xu": cp_lower(ku, reps, 0.975),
                    "cp95_upper_at_Xl": cp_upper(kl, reps, 0.975),
                    "mean_sqrt_n_dip": float(np.sqrt(n) * d.mean())})

    # anchors (PREREG 4); uniform check simultaneous over 8 (n, alpha) pairs
    cells = {(c["f"], c["n"], c["alpha"]): c for c in summary["cells"]}
    anc = {"uniform": [], "witness_n200": None}
    for n in NS:
        for a in ALPHAS:
            c = cells[(0, n, a)]
            lo = cp_lower(c["k_at_q"], reps, 1 - 0.05 / 16)
            up = cp_upper(c["k_at_q"], reps, 1 - 0.05 / 16)
            se3 = 3 * np.sqrt(a * (1 - a) / reps)
            anc["uniform"].append({"n": n, "alpha": a, "rej": c["rej_point"], "cp_lo": lo, "cp_up": up,
                                   "cp_ok": bool(lo <= a <= up),
                                   "se3_ok": bool(abs(c["rej_point"] - a) <= se3)})
    w = cells[(29, 200, 0.05)]["rej_point"]
    anc["witness_n200"] = {"rej": w, "ok": bool(w >= 0.15)}
    anc["valid"] = all(x["cp_ok"] for x in anc["uniform"]) and anc["witness_n200"]["ok"]
    summary["anchors"] = anc

    # decision (PREREG 5): 12 cells, one-sided CP at 1 - 0.05/12
    conf = 1 - 0.05 / 12
    dec = {}
    for fam, fs in (("F1", (1, 2)), ("F2", (3, 4)), ("F3", (5, 6))):
        rows = []
        for f in fs:
            for n in (200, 1000):
                c = cells[(f, n, 0.05)]
                rows.append({"f": f, "n": n, "rej": c["rej_point"],
                             "sim_lower": cp_lower(c["k_at_Xu"], reps, conf),
                             "sim_upper": cp_upper(c["k_at_Xl"], reps, conf)})
        material = any(r["sim_lower"] >= 0.075 for r in rows)
        controlled = all(r["sim_upper"] < 0.075 for r in rows)
        dec[fam] = {"cells": rows, "status": "MATERIAL" if material else
                    ("CONTROLLED" if controlled else "INCONCLUSIVE")}
    nm = sum(v["status"] == "MATERIAL" for v in dec.values())
    nc = sum(v["status"] == "CONTROLLED" for v in dec.values())
    verdict = ("VOID (anchor failed)" if not anc["valid"] else
               "PRACTICALLY BROAD ON THE PANEL" if nm >= 2 else
               "FALSIFIED ON THE PANEL" if nc == 3 else "NOT ESTABLISHED")
    summary["decision"] = {"families": dec, "verdict": verdict}
    summary["seconds"] = time.time() - t0
    name = "screen_quick.json" if quick else "screen.json"
    (OUT / name).write_text(json.dumps(summary, indent=1))
    print("verdict:", verdict, "| anchors valid:", anc["valid"], "| witness n200:", round(w, 4))
    for fam, v in dec.items():
        print(fam, v["status"], [(r["f"], r["n"], round(r["rej"], 4), round(r["sim_lower"], 4),
                                  round(r["sim_upper"], 4)) for r in v["cells"]])
    print(f"done in {time.time() - t0:.0f}s -> {OUT / name}")


if __name__ == "__main__":
    main()
