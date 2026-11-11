"""Step 3b (DT1-b) confirmation. CONFIRM tier. Frozen design: work/DT1b-PREREG.md (receipt RC-005).

Roots (CV-6): 2026110601 uniform reference (0, n, chunk); 2026110602 confirmation draws (1, f, n, chunk);
2026110603 LP-audit selection (2, f, n).
Run: python code/dt1b_confirm.py
"""

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import cp_lower, cp_upper, dips_of, g_ramp, g_uniform, g_witness  # noqa: E402, I001
from lp_dip import dip_lp  # noqa: E402
from path_smoke import vn_sample  # noqa: E402

QUICK = "--quick" in sys.argv  # tiny dry run of every code path (spawned workers see the same argv)
R_REF, REF_CHUNKS = (20_000, 20) if QUICK else (1_000_000, 20)
REPS, CHUNKS = (2_000, 10) if QUICK else (100_000, 10)
ROOT_REF, ROOT_CONF, ROOT_AUDIT = 2026110601, 2026110602, 2026110603
NS_ALL = (10, 20, 30, 50, 100, 200, 500, 1000, 2000, 5000, 12800)
NS_MAIN = (50, 100, 200, 500, 1000, 2000, 5000, 12800)
NS_PROBE = (10, 20, 30, 50, 100, 200)
ALPHAS = (0.05, 0.01)
GAMMA, DELTA = 0.001, 0.001
XS = (0.8, 1.0, 1.2, 1.4, 1.6)
FM = (0.220644, 0.560962, 0.806055, 0.928968, 0.978091)
OUT = Path(__file__).resolve().parents[1] / "deposits" / ("dt1b_quick" if QUICK else "dt1b")


def g_probe_step(rng, size):
    a, b = 0.033 * 0.994, 1.0 * 0.006
    w = b / (a + b)
    u = rng.random(size)
    return np.where(u < 1 - w, 0.994 * u / (1 - w), 0.994 + 0.006 * (u - (1 - w)) / w)


def g_probe_spike(rng, size):
    x = rng.random(size)
    hit = rng.random(size) < 0.1
    return np.where(hit, 1.0 - np.abs(rng.normal(0, 0.003, size)) % 1.0, x)


def g_probe_atom(rng, size):
    x = rng.random(size)
    return np.where(rng.random(size) < 0.1, 1.0, x)


LAWS = {0: ("U(0,1)", None, None, g_uniform, NS_ALL),
        1: ("W(0.1,0.005)", 0.1, 0.995, g_witness(w=0.1), NS_MAIN),
        2: ("ramp(w=.1,eps=.01)", 0.1, 0.99, g_ramp(), NS_MAIN),
        3: ("W(0.4,0.005)", 0.4, 0.995, g_witness(w=0.4), NS_MAIN),
        4: ("probe step#0", None, None, g_probe_step, NS_PROBE),
        5: ("probe edge spike s=.003", None, None, g_probe_spike, NS_PROBE),
        6: ("probe atom", None, None, g_probe_atom, NS_PROBE)}


def draw(f, n, c, m):
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_CONF, spawn_key=(1, f, n, c))))
    return np.sort(LAWS[f][3](rng, (m, n)), axis=1)


def run_job(job):
    kind, f, n, c = job
    if kind == "ref":
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_REF, spawn_key=(0, n, c))))
        out = np.empty(R_REF // REF_CHUNKS)
        b = max(100, min(2000, 2_000_000 // n))
        for b0 in range(0, len(out), b):
            m = min(b, len(out) - b0)
            out[b0:b0 + m] = dips_of(rng.random((m, n)))
        return job, out
    m = REPS // CHUNKS
    sub = max(100, min(m, 2_000_000 // n))
    d = np.empty(m)
    extra = None
    if f in (1, 2, 3) and n == 12800:
        extra = np.empty((m, 2))
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_CONF, spawn_key=(1, f, n, c))))
    for b0 in range(0, m, sub):
        mm = min(sub, m - b0)
        x = np.sort(LAWS[f][3](rng, (mm, n)), axis=1)
        d[b0:b0 + mm] = dips_of(x)
        if extra is not None:
            L = LAWS[f][2]
            for i in range(mm):
                k = int(np.searchsorted(x[i], L))
                vsh = vn_sample(x[i, :k].copy(), n, 0)
                top = x[i, k:]
                tt = (len(top) / n) * float(dips_of(top[None, :])[0]) if len(top) >= 2 else 0.0
                extra[b0 + i] = (vsh, tt)
    return job, (d, extra)


def audit_job(args):
    """Regenerate the chunks (same seeds and batch order) and LP-audit the selected samples."""
    f, n, idx_global, q05 = args
    m = REPS // CHUNKS
    sub = max(100, min(m, 2_000_000 // n))
    want = {int(i) for i in idx_global}
    rows = []
    for c in range(CHUNKS):
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_CONF, spawn_key=(1, f, n, c))))
        for b0 in range(0, m, sub):
            mm = min(sub, m - b0)
            x = np.sort(LAWS[f][3](rng, (mm, n)), axis=1)
            d = dips_of(x)
            for i in range(mm):
                g = c * m + b0 + i
                near = abs(d[i] - q05) < 1e-6
                if g in want or near:
                    dl = dip_lp(x[i])
                    rows.append((g, float(d[i]), float(dl), abs(float(d[i]) - dl), bool(near)))
    return (f, n), rows


def bracket(R, p, g):
    lo = int(stats.binom.ppf(g, R, p))
    while stats.binom.cdf(lo - 1, R, p) > g:
        lo -= 1
    u = int(stats.binom.isf(g, R, p)) + 1
    while stats.binom.sf(u - 1, R, p) > g:
        u += 1
    return max(lo, 1), min(u, R)


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [("ref", None, n, c) for n in NS_ALL for c in range(REF_CHUNKS)]
    jobs += [("law", f, n, c) for f, v in LAWS.items() for n in v[4] for c in range(CHUNKS)]
    jobs.sort(key=lambda j: -j[2])
    res = {}
    with ProcessPoolExecutor(max_workers=11) as ex:
        for job, val in ex.map(run_job, jobs, chunksize=1):
            res[job] = val
    print(f"simulated {len(jobs)} jobs in {time.time() - t0:.0f}s", flush=True)

    g_tail = GAMMA / (len(NS_ALL) * len(ALPHAS)) / 2
    ref_info, crit = {}, {}
    for n in NS_ALL:
        ref = np.sort(np.concatenate([res[("ref", None, n, c)] for c in range(REF_CHUNKS)]))
        np.save(OUT / f"ref_n{n}.npy", ref)
        for a in ALPHAS:
            k = int(np.ceil((1 - a) * R_REF))
            lo, hi = bracket(R_REF, 1 - a, g_tail)
            q, xl, xu = ref[k - 1], ref[lo - 1], ref[hi - 1]
            crit[(n, a)] = (q, xl, xu)
            ref_info[f"n{n}_a{a}"] = {"q_hat": float(q), "X_l": float(xl), "X_u": float(xu),
                                      "sqrt_n_q": float(np.sqrt(n) * q), "ties_at_q": int(np.sum(ref == q)) - 1}
    cells, raw = {}, {}
    for f, (name, _w, _L, _, ns) in LAWS.items():
        for n in ns:
            d = np.concatenate([res[("law", f, n, c)][0] for c in range(CHUNKS)])
            raw[(f, n)] = d
            np.save(OUT / f"dips_f{f}_n{n}.npy", d)
            for a in ALPHAS:
                q, xl, xu = crit[(n, a)]
                kq, kl, ku = int(np.sum(d >= q)), int(np.sum(d >= xl)), int(np.sum(d >= xu))
                cells[(f, n, a)] = {"f": f, "law": name, "n": n, "alpha": a, "rej_point": kq / REPS,
                                    "k_at_q": kq, "k_at_Xl": kl, "k_at_Xu": ku}
    K = len(cells)
    for c in cells.values():
        c["bonf_lower_at_Xu"] = cp_lower(c["k_at_Xu"], REPS, 1 - 0.05 / (2 * K))
        c["bonf_upper_at_Xl"] = cp_upper(c["k_at_Xl"], REPS, 1 - 0.05 / (2 * K))

    # primary: IUT over 12 cells
    prim = []
    for f in (1, 2):
        for n in (200, 500, 1000, 2000, 5000, 12800):
            c = cells[(f, n, 0.05)]
            lb = cp_lower(c["k_at_Xu"], REPS, 1 - DELTA)
            prim.append({"f": f, "n": n, "rej_point": c["rej_point"], "k_at_Xu": c["k_at_Xu"],
                         "cp_lower_1-delta": lb, "pass": bool(lb > 0.15)})
    primary = {"cells": prim, "claim_E_3b": "CONFIRMED" if all(p["pass"] for p in prim) else "NOT CONFIRMED",
               "error_budget": DELTA + GAMMA}

    # S2 anchor
    s2 = []
    for n in NS_ALL:
        for a in ALPHAS:
            c = cells[(0, n, a)]
            lo = cp_lower(c["k_at_q"], REPS, 1 - 0.05 / 44)
            up = cp_upper(c["k_at_q"], REPS, 1 - 0.05 / 44)
            s2.append({"n": n, "alpha": a, "rej": c["rej_point"], "ok": bool(lo <= a <= up)})

    # S3 probe re-derivation
    probe = {(4, 0.05): (9.1, 11.4, 14.8, 19.2, 24.3, 24.3), (4, 0.01): (1.9, 3.5, 4.4, 6.5, 9.8, 10.5),
             (5, 0.05): (7.4, 8.7, 10.8, 14.2, 21.2, 26.1)}
    s3 = []
    for (f, a), vals in probe.items():
        for n, pv in zip((10, 20, 30, 50, 100, 200), vals, strict=True):
            nv = 100 * cells[(f, n, a)]["rej_point"]
            s3.append({"law": LAWS[f][0], "alpha": a, "n": n, "probe": pv, "new": round(nv, 2),
                       "diff": round(nv - pv, 2), "within_1pt": bool(abs(nv - pv) <= 1.0)})
    for a, pv in ((0.05, 26.5), (0.01, 11.1)):
        nv = 100 * cells[(6, 200, a)]["rej_point"]
        s3.append({"law": "probe atom", "alpha": a, "n": 200, "probe": pv, "new": round(nv, 2),
                   "diff": round(nv - pv, 2), "within_1pt": bool(abs(nv - pv) <= 1.0)})

    # S4 / S5 at n = 12800
    s4, s5 = {}, {}
    sref = np.sqrt(12800) * np.load(OUT / "ref_n12800.npy")
    for f in (1, 2, 3):
        name, w, L, _, _ = LAWS[f]
        d = raw[(f, 12800)]
        ex = np.concatenate([res[("law", f, 12800, c)][1] for c in range(CHUNKS)])
        z = 2 * np.sqrt(12800) * d / np.sqrt(1 - w)
        target = []
        for x, fm in zip(XS, FM, strict=True):
            hu = float(np.mean(sref <= (x * np.sqrt(1 - w) / 2) / np.sqrt(w))) if f == 3 else 1.0
            target.append(fm * hu)
        pts = []
        zc = stats.norm.ppf(1 - 0.05 / 30)
        for x, tv in zip(XS, target, strict=True):
            e = float(np.mean(z <= x))
            se = np.sqrt(e * (1 - e) / REPS)
            pts.append({"x": x, "ecdf": e, "target": tv, "diff": e - tv,
                        "ci": [e - tv - zc * se, e - tv + zc * se],
                        "inside_0.01": bool(abs(e - tv) + zc * se < 0.01)})
        qh = {a: crit[(12800, a)][0] * np.sqrt(12800) for a in ALPHAS}
        pred = {}
        for a in ALPHAS:
            t2 = 2 * qh[a] / np.sqrt(1 - w)
            fmv = float(np.interp(t2, *bf_table()))
            hu = float(np.mean(sref <= qh[a] / np.sqrt(w)))
            pred[a] = {"predicted": 1 - fmv * hu, "observed": cells[(f, 12800, a)]["rej_point"]}
        s4[name] = {"points": pts, "five_point_equivalence": all(p["inside_0.01"] for p in pts), "rejection": pred}
        vsh, tt = ex[:, 0], ex[:, 1]
        sq = np.sqrt(12800)
        s5[name] = {"S_n_q01_q50_q99_max": [float(v) for v in np.quantile(sq * (d - np.maximum(vsh, tt)), (0.01, 0.5, 0.99, 1))],
                    "sqrtn(D-Vsh)_q01_q50_q99_max": [float(v) for v in np.quantile(sq * (d - vsh), (0.01, 0.5, 0.99, 1))],
                    "P(D<Vsh-1e-12)": float(np.mean(d < vsh - 1e-12)),
                    "P(D=max)": float(np.mean(np.abs(d - np.maximum(vsh, tt)) <= 1e-12))}

    # S6 LP audit
    ajobs = []
    for f in (0, 1, 2, 3):
        for n in (50, 100, 200):
            rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_AUDIT, spawn_key=(2, f, n))))
            idx = rng.choice(REPS, size=5 if QUICK else 100, replace=False)
            ajobs.append((f, n, idx, crit[(n, 0.05)][0]))
    audit = {}
    with ProcessPoolExecutor(max_workers=11) as ex:
        for key, rows in ex.map(audit_job, ajobs, chunksize=1):
            audit[f"f{key[0]}_n{key[1]}"] = {"audited": len(rows), "near_threshold": sum(r[4] for r in rows),
                                             "max_abs_diff": max(r[3] for r in rows), "pass": all(r[3] <= 1e-9 for r in rows)}
    s6 = {"cells": audit, "all_pass": all(v["pass"] for v in audit.values())}

    summary = {"tier": "CONFIRM", "roots": [ROOT_REF, ROOT_CONF, ROOT_AUDIT], "R_ref": R_REF, "reps": REPS,
               "reference": ref_info, "primary": primary, "S2_anchor": {"cells": s2, "all_ok": all(x["ok"] for x in s2)},
               "S3_probe": s3, "S4_limit": s4, "S5_reduction": s5, "S6_lp_audit": s6,
               "cells": [cells[k] for k in sorted(cells)], "seconds": time.time() - t0}
    (OUT / "confirm.json").write_text(json.dumps(summary, indent=1))
    print("PRIMARY:", primary["claim_E_3b"], [(p["f"], p["n"], round(p["rej_point"], 4), round(p["cp_lower_1-delta"], 4)) for p in prim])
    print("S2 anchor all ok:", summary["S2_anchor"]["all_ok"], "| S6 LP audit all pass:", s6["all_pass"])
    print(json.dumps(s3, indent=0)[:3000])
    print(json.dumps(s4, indent=0)[:4000])
    print(json.dumps(s5, indent=0))
    print(json.dumps(audit, indent=0))
    print(f"done in {time.time() - t0:.0f}s")


def bf_table():
    """F_M on a grid from Balabdaoui & Filali 2011 Table 2 (transcribed values used in this project)."""
    xs = np.round(np.arange(0.90, 1.80, 0.01), 2)
    fv = [0.393227821925, 0.410762393335, 0.428198510241, 0.445505815829, 0.462656184022, 0.479623704507,
          0.496384649825, 0.512917427370, 0.529202518946, 0.545222410267, 0.560961512572, 0.576406078288,
          0.591544112452, 0.606365281379, 0.620860819881, 0.635023438140, 0.648847229167, 0.662327577632,
          0.675461070708, 0.688245411417, 0.700679334912, 0.712762527962, 0.724495551883, 0.735879769041,
          0.746917273011, 0.757610822418, 0.767963778457, 0.777980046011, 0.787664018322, 0.797020525083,
          0.806054783852, 0.814772354647, 0.823179097587, 0.831281133430, 0.839084806872, 0.846596652448,
          0.853823362907, 0.860771759907, 0.867448766898, 0.873861384073, 0.880016665251, 0.885921696569,
          0.891583576893, 0.897009399815, 0.902206237159, 0.907181123890, 0.911941044337, 0.916492919660,
          0.920843596472, 0.924999836546, 0.928968307546, 0.932755574715, 0.936368093452, 0.939812202742,
          0.943094119365, 0.946219932852, 0.949195601129, 0.952026946811, 0.954719654107, 0.957279266289,
          0.959711183695, 0.962020662224, 0.964212812296, 0.966292598240, 0.968264838078, 0.970134203687,
          0.971905221291, 0.973582272276, 0.975169594296, 0.976671282637, 0.978091291833, 0.979433437494,
          0.980701398339, 0.981898718407, 0.983028809430, 0.984094953345, 0.985100304937, 0.986047894590,
          0.986940631128, 0.987781304739, 0.988572589969, 0.989317048762, 0.990017133547, 0.990675190351,
          0.991293461936, 0.991874090944, 0.992419123041, 0.992930510053, 0.993410113095, 0.993859705663]
    return xs, np.array(fv)


if __name__ == "__main__":
    main()
