"""Step 4 (DT2) conjecture farm, part 1: the block formula (BF). DEV tier. Frozen design: work/DT2-PREREG.md.

Root 2026110610: reference spawn (0, n, chunk); laws spawn (1, f, n, chunk).
Run: python code/dt2_farm.py [--quick]
"""

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import dips_of  # noqa: E402, I001
from dt1b_confirm import bf_table  # noqa: E402
from path_smoke import vn_sample  # noqa: E402

QUICK = "--quick" in sys.argv
ROOT = 2026110610
REPS, CHUNKS = (400, 4) if QUICK else (20_000, 10)
REF_R = {800: 1_000_000, 3200: 1_000_000, 12800: 1_000_000, 51200: 500_000}
if QUICK:
    REF_R = dict.fromkeys(REF_R, 20_000)
REF_CHUNK = 20_000 if QUICK else 50_000
NS_A = (800, 3200, 12800)
NS_C = (3200, 12800, 51200)
OUT = Path(__file__).resolve().parents[1] / "deposits" / ("dt2_quick" if QUICK else "dt2")


def law_spec(f, n):
    """Return (edges, masses, modal_index)."""
    if f == 0:
        return [0, .995, 1], [.9, .1], 1
    if f == 1:
        return [0, .5, .995, 1], [.3, .6, .1], 2
    if f == 2:
        return [0, .5, .995, 1], [.43, .47, .1], 2
    if f == 3:
        return [0, .33, .66, .995, 1], [.2, .3, .4, .1], 3
    if f == 4:
        return [0, .6, .605, 1], [.6, .1, .3], 1
    if f == 5:
        return [0, .4975, .5025, 1], [.45, .1, .45], 1
    if f == 6:
        return [0, .9995, 1], [.9, .1], 1
    if f == 7:
        return [0, .95, 1], [.9, .1], 1
    if f == 8:
        return [0, .005, 1], [.1, .9], 0
    if f == 9:
        d = 4 / np.sqrt(n)
        c = 0.9 / (0.5 + 0.495 * (1 + d))
        return [0, .5, .995, 1], [c * 0.5, c * (1 + d) * 0.495, .1], 2
    if 10 <= f <= 13:
        w = (.05, .02, .01, .005)[f - 10]
        eps = 0.25 * w
        return [0, 1 - eps, 1], [1 - w, w], 1
    if f == 14:
        return [0, .995, 1], [.98, .02], 1
    raise ValueError(f)


NAMES = {0: "W(0.1,.005)", 1: "2 shoulders big jump", 2: "2 shoulders small jump", 3: "3-block staircase",
         4: "interior plateau unequal", 5: "interior plateau symmetric", 6: "W narrow top", 7: "W wide top",
         8: "W reflected", 9: "shrinking jump 4/sqrt(n)", 10: "K4 w=.05 eps=.25w", 11: "K4 w=.02 eps=.25w",
         12: "K4 w=.01 eps=.25w", 13: "K4 w=.005 eps=.25w", 14: "W(0.02,.005)"}
NS_OF = {f: (NS_A if f <= 9 else NS_C) for f in NAMES}


def sample(rng, size, edges, masses):
    e, p = np.asarray(edges, float), np.asarray(masses, float)
    cm = np.concatenate([[0.0], np.cumsum(p)])
    u = rng.random(size)
    k = np.clip(np.searchsorted(cm, u, side="right") - 1, 0, len(p) - 1)
    return e[k] + (e[k + 1] - e[k]) * (u - cm[k]) / p[k]


def block_terms(x, n, edges, m):
    """Block terms for one sorted sample x."""
    idx = np.searchsorted(x, edges[1:-1])
    bounds = np.concatenate([[0], idx, [n]])
    terms = []
    for i in range(len(edges) - 1):
        a, b = int(bounds[i]), int(bounds[i + 1])
        pts = x[a:b]
        if b - a < 2:
            terms.append(0.0)
        elif i < m:
            terms.append(vn_sample(pts.copy(), n, a))
        elif i > m:
            terms.append(vn_sample((-pts[::-1]).copy(), n, n - b))
        else:
            terms.append((b - a) / n * float(dips_of(pts[None, :])[0]))
    return terms


def run_job(job):
    kind, f, n, c, m_reps = job
    if kind == "ref":
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT, spawn_key=(0, n, c))))
        out = np.empty(m_reps)
        b = max(50, min(2000, 2_000_000 // n))
        for b0 in range(0, m_reps, b):
            mm = min(b, m_reps - b0)
            out[b0:b0 + mm] = dips_of(rng.random((mm, n)))
        return job, out
    edges, masses, mi = law_spec(f, n)
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT, spawn_key=(1, f, n, c))))
    sub = max(50, min(m_reps, 2_000_000 // n))
    d = np.empty(m_reps)
    bf = np.empty(m_reps)
    for b0 in range(0, m_reps, sub):
        mm = min(sub, m_reps - b0)
        x = np.sort(sample(rng, (mm, n), edges, masses), axis=1)
        d[b0:b0 + mm] = dips_of(x)
        for i in range(mm):
            bf[b0 + i] = max(block_terms(x[i], n, edges, mi))
    return job, (d, bf)


XS_BF, FV_BF = bf_table()


def FM(x):
    # 2026-11-06 (DT3c-RESULT §2): was a guessed extension below 0.90 and above 1.79.
    # Now the full published table (0.33..2.54); identical to the old values on 0.90..1.79.
    # Old function kept verbatim as dt3c_reanalysis.fm_legacy for reproducing frozen predictions.
    from fm_table import FM_pub
    return FM_pub(float(x), clip_tails=True)


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = []
    for n, R in REF_R.items():
        for c in range(R // REF_CHUNK):
            jobs.append(("ref", None, n, c, REF_CHUNK))
    for f in NAMES:
        for n in NS_OF[f]:
            for c in range(CHUNKS):
                jobs.append(("law", f, n, c, REPS // CHUNKS))
    jobs.sort(key=lambda j: -j[2])
    res = {}
    with ProcessPoolExecutor(max_workers=11) as ex:
        for job, val in ex.map(run_job, jobs, chunksize=1):
            res[job[:4]] = val
    print(f"simulated {len(jobs)} jobs in {time.time() - t0:.0f}s", flush=True)
    ref = {}
    for n, R in REF_R.items():
        r = np.sort(np.concatenate([res[("ref", None, n, c)] for c in range(R // REF_CHUNK)]))
        ref[n] = np.sqrt(n) * r
        if not QUICK:
            np.save(OUT / f"ref_n{n}.npy", r)
    out = {"tier": "DEV", "root": ROOT, "reps": REPS, "cells": []}
    for f in NAMES:
        for n in NS_OF[f]:
            d = np.concatenate([res[("law", f, n, c)][0] for c in range(CHUNKS)])
            bf = np.concatenate([res[("law", f, n, c)][1] for c in range(CHUNKS)])
            if not QUICK:
                np.save(OUT / f"farm_f{f}_n{n}.npy", np.stack([d, bf], axis=1))
            edges, masses, mi = law_spec(f, n)
            row = {"f": f, "name": NAMES[f], "n": n}
            for a in (0.05, 0.01):
                q = float(np.quantile(ref[n], 1 - a, method="inverted_cdf"))
                obs = float(np.mean(np.sqrt(n) * d >= q))
                prod = 1.0
                for i, p in enumerate(masses):
                    if i == mi:
                        prod *= float(np.mean(ref[n] <= q / np.sqrt(p)))
                    else:
                        prod *= FM(2 * q / np.sqrt(p))
                row[f"a{a}"] = {"obs": obs, "r_BF": 1 - prod, "diff_pts": 100 * (obs - (1 - prod)),
                                "obs_BF_n": float(np.mean(np.sqrt(n) * bf >= q))}
            rn = np.sqrt(n) * (d - bf)
            row["R_n_q01_q50_q99_max"] = [float(v) for v in np.quantile(rn, (0.01, 0.5, 0.99, 1.0))]
            row["P(|R_n|<=1e-12)"] = float(np.mean(np.abs(d - bf) <= 1e-12))
            row["q99_absR"] = float(np.quantile(np.abs(rn), 0.99))
            out["cells"].append(row)
    cells = {(c["f"], c["n"]): c for c in out["cells"]}
    rd = {}
    p1 = {f: bool(abs(cells[(f, 12800)]["a0.05"]["diff_pts"]) <= 1.5) for f in (0, 1, 3, 4, 5, 6, 7, 8)}
    rd["P1_BF_consistent_n12800"] = {NAMES[f]: v for f, v in p1.items()}
    rd["P1"] = "BF SURVIVES" if all(p1.values()) else "BF KILLED for " + str([NAMES[f] for f, v in p1.items() if not v])
    o2 = [cells[(2, n)]["a0.05"]["obs"] for n in NS_A]
    rd["f2_slow_merger"] = bool(o2[0] > o2[1] > o2[2])
    o9 = [cells[(9, n)]["a0.05"]["obs"] for n in NS_A]
    rd["f9_local_transition"] = bool(abs(o9[1] - o9[2]) <= 0.02)
    k4 = {}
    for f in (10, 11, 12, 13):
        dd = [abs(cells[(f, n)]["a0.05"]["diff_pts"]) for n in NS_C]
        k4[NAMES[f]] = {"absdiff_pts": dd, "approach": bool(dd[0] >= dd[1] >= dd[2] - 0.5 and dd[2] <= 2.0)}
    rd["K4"] = k4
    o14 = [cells[(14, n)]["a0.05"]["obs"] for n in NS_C]
    rd["f14_pinning_transition"] = bool(o14[0] < o14[1] < o14[2])
    rd["Rn_q99_decreasing"] = {NAMES[f]: [cells[(f, n)]["q99_absR"] for n in NS_A] for f in range(9)}
    out["readings"] = rd
    out["seconds"] = time.time() - t0
    (OUT / "farm.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(rd, indent=1))
    for c in out["cells"]:
        a5, a1 = c["a0.05"], c["a0.01"]
        print(f"{c['f']:2d} {c['name'][:26]:26s} n={c['n']:5d} obs5={100*a5['obs']:6.2f} BF5={100*a5['r_BF']:6.2f} "
              f"d={a5['diff_pts']:+6.2f} | obs1={100*a1['obs']:5.2f} BF1={100*a1['r_BF']:5.2f} | "
              f"P(R=0)={c['P(|R_n|<=1e-12)']:.3f} q99|R|={c['q99_absR']:.4f}")
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
