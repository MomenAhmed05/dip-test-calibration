"""Step 4 farm part 2: the pinning curve g(kappa, lambda) and the shrinking jump. DEV. Design: work/DT2b-PREREG.md.

Root 2026110611. Run: python code/dt2b_pinning.py [--quick]
"""

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import dips_of  # noqa: E402, I001
from dt2_farm import FM, sample  # noqa: E402

QUICK = "--quick" in sys.argv
ROOT = 2026110611
NS = (3200, 12800)
REPS, CHUNKS = (400, 4) if QUICK else (20_000, 10)
REF_R, REF_CHUNK = (20_000, 20_000) if QUICK else (1_000_000, 50_000)
KAPPAS = (0, .2, .4, .6, .8, 1, 1.2, 1.4, 1.7, 2, 2.5, 3, 4)
LAMBDAS = (0.1, 1, 4)
CS = (0, 1, 2, 4, 8, 16)
OUT = Path(__file__).resolve().parents[1] / "deposits" / ("dt2b_quick" if QUICK else "dt2b")

LAWS = [("P", k, lam) for k in KAPPAS for lam in LAMBDAS] + [("S", c, None) for c in CS]


def spec(law, n):
    kind, a, b = law
    if kind == "P":
        eps, e = b / np.sqrt(n), a / np.sqrt(n)
        w = eps + e * (1 - eps)
        return [0, 1 - eps, 1], [1 - w, w], 1
    d = a / np.sqrt(n)
    h1 = 0.9 / (0.5 + 0.495 * (1 + d))
    return [0, .5, .995, 1], [h1 * .5, h1 * (1 + d) * .495, .1], 2


def run_job(job):
    kind, li, n, c, m = job
    if kind == "ref":
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT, spawn_key=(0, n, c))))
        out = np.empty(m)
        b = max(50, min(2000, 2_000_000 // n))
        for b0 in range(0, m, b):
            mm = min(b, m - b0)
            out[b0:b0 + mm] = dips_of(rng.random((mm, n)))
        return job, out
    edges, masses, _ = spec(LAWS[li], n)
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT, spawn_key=(1, li, n, c))))
    sub = max(50, min(m, 2_000_000 // n))
    d = np.empty(m)
    for b0 in range(0, m, sub):
        mm = min(sub, m - b0)
        d[b0:b0 + mm] = dips_of(np.sort(sample(rng, (mm, n), edges, masses), axis=1))
    return job, d


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [("ref", None, n, c, REF_CHUNK) for n in NS for c in range(REF_R // REF_CHUNK)]
    jobs += [("law", li, n, c, REPS // CHUNKS) for li in range(len(LAWS)) for n in NS for c in range(CHUNKS)]
    jobs.sort(key=lambda j: -j[2])
    res = {}
    with ProcessPoolExecutor(max_workers=11) as ex:
        for job, val in ex.map(run_job, jobs, chunksize=1):
            res[job[:4]] = val
    ref = {n: np.sqrt(n) * np.sort(np.concatenate([res[("ref", None, n, c)] for c in range(REF_R // REF_CHUNK)]))
           for n in NS}
    q = {n: float(np.quantile(ref[n], 0.95, method="inverted_cdf")) for n in NS}
    q1 = {n: float(np.quantile(ref[n], 0.99, method="inverted_cdf")) for n in NS}
    rsup = {n: 1 - FM(2 * q[n]) for n in NS}
    rows = []
    g = {}
    for li, law in enumerate(LAWS):
        for n in NS:
            d = np.sqrt(n) * np.concatenate([res[("law", li, n, c)] for c in range(CHUNKS)])
            if not QUICK:
                np.save(OUT / f"dips_l{li}_n{n}.npy", d / np.sqrt(n))
            r5, r1 = float(np.mean(d >= q[n])), float(np.mean(d >= q1[n]))
            g[(li, n)] = r5
            row = {"law": law, "n": n, "rej5": r5, "rej1": r1, "q50_q95": [float(v) for v in np.quantile(d, (0.5, 0.95))]}
            if law[0] == "S":
                _, masses, _ = spec(law, n)
                sep = 1 - FM(2 * q[n] / np.sqrt(masses[0])) * FM(2 * q[n] / np.sqrt(masses[1])) * float(np.mean(ref[n] <= q[n] / np.sqrt(.1)))
                mer = 1 - FM(2 * q[n] / np.sqrt(.9)) * float(np.mean(ref[n] <= q[n] / np.sqrt(.1)))
                row.update({"separated": sep, "merged": mer})
            rows.append(row)
    idx = {law: i for i, law in enumerate(LAWS)}
    rd = {}
    spread = {k: 100 * (max(g[(idx[("P", k, lam)], 12800)] for lam in LAMBDAS)
                        - min(g[(idx[("P", k, lam)], 12800)] for lam in LAMBDAS)) for k in KAPPAS}
    rd["a_lambda_spread_pts_n12800"] = spread
    rd["a"] = "COLLAPSE" if all(v <= 2.0 for v in spread.values()) else \
        "LAMBDA MATTERS at kappa " + str([k for k, v in spread.items() if v > 2.0])
    stab = {k: 100 * abs(g[(idx[("P", k, 1)], 3200)] - g[(idx[("P", k, 1)], 12800)]) for k in KAPPAS}
    rd["b_nstab_pts"] = stab
    rd["b"] = "STABLE" if all(v <= 2.0 for v in stab.values()) else "NOT STABLE at " + str([k for k, v in stab.items() if v > 2.0])
    mx = {n: max(g[(idx[("P", k, lam)], n)] for k in KAPPAS for lam in LAMBDAS) for n in NS}
    rd["c_max_vs_rsup"] = {n: (100 * mx[n], 100 * rsup[n]) for n in NS}
    rd["c"] = "NO OVERSHOOT" if all(mx[n] <= rsup[n] + 0.01 for n in NS) else "OVERSHOOT"
    srows = [r for r in rows if r["law"][0] == "S"]
    inside = all(r["separated"] - 0.01 <= r["rej5"] <= r["merged"] + 0.01 for r in srows)
    mono = all(g[(idx[("S", CS[i + 1], None)], n)] <= g[(idx[("S", CS[i], None)], n)] + 0.01 for n in NS for i in range(len(CS) - 1))
    sst = all(abs(g[(idx[("S", c, None)], 3200)] - g[(idx[("S", c, None)], 12800)]) <= 0.02 for c in CS)
    rd["d"] = {"inside_bounds": bool(inside), "nonincreasing": bool(mono), "n_stable": bool(sst),
               "verdict": "SANDWICH + LOCAL LIMIT" if (inside and mono and sst) else "NOT ALL"}
    out = {"tier": "DEV", "root": ROOT, "reps": REPS, "q95_sqrtn": q, "q99_sqrtn": q1, "r_sup": rsup,
           "rows": rows, "readings": rd, "seconds": time.time() - t0}
    (OUT / "pinning.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(rd, indent=1))
    print("g(kappa) at 5%, n=12800, by lambda:")
    for k in KAPPAS:
        print(f"  kappa={k:4}: " + "  ".join(f"lam={lam}: {100 * g[(idx[('P', k, lam)], 12800)]:6.2f}" for lam in LAMBDAS)
              + f"   | lam=1 n=3200: {100 * g[(idx[('P', k, 1)], 3200)]:6.2f}")
    for r in srows:
        print(f"  S c={r['law'][1]:2} n={r['n']:5}: {100 * r['rej5']:6.2f}  [sep {100 * r['separated']:5.2f}, merged {100 * r['merged']:5.2f}]")
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
