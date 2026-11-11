"""Step 5 CONFIRM (DT3c): block formula at large n + pinning-overshoot check. Frozen design: work/DT3c-PREREG.md.

Roots: 2026110630 reference spawn (0, chunk); 2026110631 laws spawn (1, law_index, chunk).
Nested common random numbers: one stream of NMAX uniforms per replicate; each n uses the prefix u[:n].
Chunks are saved to deposits/dt3c/parts/ as they finish; a rerun skips finished chunks (resume).
Run: python code/dt3c_confirm.py [--quick]
"""

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import dips_of  # noqa: E402, I001
from dt2_farm import FM, block_terms  # noqa: E402

QUICK = "--quick" in sys.argv
# quick dry runs use DEV roots 2026110632/33 so the CONFIRM streams are never seen before the full run
ROOT_REF, ROOT_LAW = (2026110632, 2026110633) if QUICK else (2026110630, 2026110631)
NS = (3200, 12800, 51200, 204800)
NMAX = NS[-1]
REF_R, LAW_R, CHUNK = (440, 110, 110) if QUICK else (1_000_000, 100_000, 5_000)
BATCH = 10
OUT = Path(__file__).resolve().parents[1] / "deposits" / ("dt3c_quick" if QUICK else "dt3c")
PARTS = OUT / "parts"
ALPHAS = (0.05, 0.01)
TAILS = (0.50, 0.75, 0.90, 0.95, 0.99)

BF_LAWS = [
    ("L1 interior modal plateau", [0, .30, .45, .55, 1], [.24, .15, .20, .41], 2),
    ("L2 W(0.1,.005)", [0, .995, 1], [.90, .10], 1),
    ("L3 unbalanced masses", [0, .98, .99, 1], [.95, .02, .03], 2),
    ("L4 close heights", [0, .50, .995, 1], [.43, .47, .10], 2),
    ("L5 5-step staircase", [0, .2, .4, .6, .8, 1], [.12, .16, .20, .24, .28], 4),
]
PIN = [("P1", 1.4, 0.1), ("P2", 1.7, 0.1), ("P3", 1.2, 0.1), ("P4", 2.5, 1.0)]
NAMES = [b[0] for b in BF_LAWS] + [f"{p[0]} pin k={p[1]} lam={p[2]}" for p in PIN]


def spec(li, n):
    """Return (edges, masses, modal_index) for law index li at sample size n."""
    if li < len(BF_LAWS):
        _, e, m, mi = BF_LAWS[li]
        return e, m, mi
    _, k, lam = PIN[li - len(BF_LAWS)]
    eps, e = lam / np.sqrt(n), k / np.sqrt(n)
    w = eps + e * (1 - eps)
    return [0, 1 - eps, 1], [1 - w, w], 1


def qmap(u, edges, masses):
    """Quantile map of a step density applied to uniforms u (same logic as dt2_farm.sample)."""
    e, p = np.asarray(edges, float), np.asarray(masses, float)
    cm = np.concatenate([[0.0], np.cumsum(p)])
    k = np.clip(np.searchsorted(cm, u, side="right") - 1, 0, len(p) - 1)
    return e[k] + (e[k + 1] - e[k]) * (u - cm[k]) / p[k]


def part_path(job):
    kind, li, c = job
    return PARTS / (f"ref_c{c}.npy" if kind == "ref" else f"law{li}_c{c}.npy")


def run_job(job):
    kind, li, c = job
    path = part_path(job)
    if path.exists():
        return job, "skip"
    root, key = (ROOT_REF, (0, c)) if kind == "ref" else (ROOT_LAW, (1, li, c))
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(root, spawn_key=key)))
    nr = CHUNK
    # columns: per n -> D (ref) or (D, BF_n) (law)
    width = len(NS) if kind == "ref" else 2 * len(NS)
    out = np.empty((nr, width))
    for b0 in range(0, nr, BATCH):
        mm = min(BATCH, nr - b0)
        u = rng.random((mm, NMAX))
        for j, n in enumerate(NS):
            if kind == "ref":
                x = np.sort(u[:, :n], axis=1)
                out[b0:b0 + mm, j] = dips_of(x)
            else:
                edges, masses, mi = spec(li, n)
                x = np.sort(qmap(u[:, :n], edges, masses), axis=1)
                out[b0:b0 + mm, 2 * j] = dips_of(x)
                for i in range(mm):
                    out[b0 + i, 2 * j + 1] = max(block_terms(x[i], n, edges, mi))
    tmp = path.with_suffix(".tmp.npy")
    np.save(tmp, out)
    tmp.replace(path)
    return job, "done"


def bf_pred(t, n, li, ref_sorted):
    """BF prediction P(limit > t) with own-n reference for the modal block."""
    _, masses, mi = spec(li, n)
    prod = 1.0
    for i, p in enumerate(masses):
        if i == mi:
            prod *= np.searchsorted(ref_sorted, t / np.sqrt(p), side="right") / len(ref_sorted)
        else:
            prod *= FM(2 * t / np.sqrt(p))
    return 1.0 - prod


def main():
    t0 = time.time()
    PARTS.mkdir(parents=True, exist_ok=True)
    jobs = [("law", li, c) for li in range(len(NAMES)) for c in range(LAW_R // CHUNK)]
    jobs += [("ref", None, c) for c in range(REF_R // CHUNK)]
    done = 0
    with ProcessPoolExecutor(max_workers=11) as ex:
        for job, status in ex.map(run_job, jobs, chunksize=1):
            done += 1
            if done % 20 == 0 or done == len(jobs):
                print(f"[{time.time() - t0:7.0f}s] {done}/{len(jobs)} chunks ({status} {job})", flush=True)
    print(f"simulated in {time.time() - t0:.0f}s", flush=True)

    ref = np.concatenate([np.load(part_path(("ref", None, c))) for c in range(REF_R // CHUNK)])
    refs = {n: np.sort(np.sqrt(n) * ref[:, j]) for j, n in enumerate(NS)}
    qs = {n: {a: float(np.quantile(refs[n], 1 - a, method="inverted_cdf")) for a in ALPHAS} for n in NS}
    tails = {n: {p: float(np.quantile(refs[n], p, method="inverted_cdf")) for p in TAILS} for n in NS}
    out = {"tier": "QUICK" if QUICK else "CONFIRM", "roots": [ROOT_REF, ROOT_LAW], "ref_reps": REF_R,
           "law_reps": LAW_R, "q_sqrtn": {str(n): {str(a): qs[n][a] for a in ALPHAS} for n in NS},
           "sup_bound": {str(n): {str(a): 1 - FM(2 * qs[n][a]) for a in ALPHAS} for n in NS}, "cells": []}
    for li, name in enumerate(NAMES):
        arr = np.concatenate([np.load(part_path(("law", li, c))) for c in range(LAW_R // CHUNK)])
        if not QUICK:
            np.save(OUT / f"law{li}.npy", arr)
        for j, n in enumerate(NS):
            d, bf = arr[:, 2 * j], arr[:, 2 * j + 1]
            sd = np.sqrt(n) * d
            row = {"law": li, "name": name, "n": n}
            for a in ALPHAS:
                q = qs[n][a]
                obs = float(np.mean(sd > q))
                row[f"a{a}"] = {"obs": obs, "se": float(np.sqrt(obs * (1 - obs) / len(d))),
                                "r_BF": bf_pred(q, n, li, refs[n]), "bench": 1 - FM(2 * q),
                                "obs_BFn": float(np.mean(np.sqrt(n) * bf > q))}
                row[f"a{a}"]["diff_pts"] = 100 * (obs - row[f"a{a}"]["r_BF"])
            row["tails"] = {str(p): {"t": tails[n][p], "obs": float(np.mean(sd > tails[n][p])),
                                     "pred": bf_pred(tails[n][p], n, li, refs[n])} for p in TAILS}
            for p in TAILS:
                tp = row["tails"][str(p)]
                tp["diff_pts"] = 100 * (tp["obs"] - tp["pred"])
            rn = np.sqrt(n) * (d - bf)
            row["Rn_q50_q99_q999_max"] = [float(v) for v in np.quantile(rn, (0.5, 0.99, 0.999, 1.0))]
            row["q99_absR"] = float(np.quantile(np.abs(rn), 0.99))
            row["P(|R|>0.01)"] = float(np.mean(np.abs(rn) > 0.01))
            row["P(|R|>0.001)"] = float(np.mean(np.abs(rn) > 0.001))
            row["P(R<-1e-9)"] = float(np.mean(d - bf < -1e-9))
            out["cells"].append(row)

    cells = {(c["law"], c["n"]): c for c in out["cells"]}
    rd = {}
    for li in range(len(BF_LAWS)):
        c5, c51 = cells[(li, NS[-1])], cells[(li, NS[-2])]
        maxtail = max(abs(v["diff_pts"]) for v in c5["tails"].values())
        q99a, q99b = cells[(li, NS[0])]["q99_absR"], c5["q99_absR"]
        conf = (abs(c5["a0.05"]["diff_pts"]) <= 1.0 and abs(c5["a0.01"]["diff_pts"]) <= 0.7 and maxtail <= 1.5
                and (q99b <= 0.5 * q99a or q99b <= 0.002))
        refuted = False
        for key in ["a0.05"] + [None]:
            if key:
                d51, d5 = c51[key]["diff_pts"], c5[key]["diff_pts"]
                pairs = [(d51, d5)]
            else:
                pairs = [(c51["tails"][str(p)]["diff_pts"], c5["tails"][str(p)]["diff_pts"]) for p in TAILS]
            for d51, d5 in pairs:
                if abs(d51) > 1.5 and abs(d5) > 1.5 and np.sign(d51) == np.sign(d5) and abs(d5) >= abs(d51) - 0.3:
                    refuted = True
        rd[NAMES[li]] = {"verdict": "REFUTED" if refuted else ("CONFIRMED" if conf else "INCONCLUSIVE (slow)"),
                         "diff5_by_n": [round(cells[(li, n)]["a0.05"]["diff_pts"], 3) for n in NS],
                         "diff1_by_n": [round(cells[(li, n)]["a0.01"]["diff_pts"], 3) for n in NS],
                         "max_tail_diff_204800": round(maxtail, 3),
                         "q99absR_by_n": [round(cells[(li, n)]["q99_absR"], 5) for n in NS],
                         "P|R|>0.01_by_n": [cells[(li, n)]["P(|R|>0.01)"] for n in NS],
                         "P(R<-1e-9)_by_n": [cells[(li, n)]["P(R<-1e-9)"] for n in NS]}
    for k in range(len(PIN)):
        li = len(BF_LAWS) + k
        ex = [(cells[(li, n)]["a0.05"]["obs"] - cells[(li, n)]["a0.05"]["bench"]) / cells[(li, n)]["a0.05"]["se"]
              for n in NS]
        over = ex[2] > 3 and ex[3] > 3
        none = all(e <= 1 for e in ex)
        rd[NAMES[li]] = {"verdict": "OVERSHOOT" if over else ("NO-OVERSHOOT" if none else "UNCLEAR"),
                         "obs5_by_n": [round(100 * cells[(li, n)]["a0.05"]["obs"], 3) for n in NS],
                         "bench5_by_n": [round(100 * cells[(li, n)]["a0.05"]["bench"], 3) for n in NS],
                         "z_excess_by_n": [round(float(e), 2) for e in ex],
                         "obs1_by_n": [round(100 * cells[(li, n)]["a0.01"]["obs"], 3) for n in NS]}
    out["readings"] = rd
    out["seconds"] = time.time() - t0
    (OUT / "confirm.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(rd, indent=1))
    print("q_sqrtn:", {n: {a: round(qs[n][a], 4) for a in ALPHAS} for n in NS})
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
