"""Step 4b (DT2-ML): learned adversary for the uniform-calibrated dip test. DEV search + fresh re-estimates.

Frozen design: work/DT2ML-PREREG.md. Generator: unimodal 64-bin step density with learned monotone log-heights
(mode bin fixed per run). Signal: score-function (REINFORCE) gradient of a smoothed rejection indicator.
Roots: 2026110620 training, 2026110621 re-estimation, 2026110622 baseline.
Run: python code/ml_adversary.py [--quick]
"""

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import dips_of  # noqa: E402, I001
from dt2_farm import sample  # noqa: E402

QUICK = "--quick" in sys.argv
K = 64
N_TRAIN = 200
ITERS, BATCH = (20, 256) if QUICK else (300, 2048)
REEST = 2_000 if QUICK else 100_000
GRID = 6 if QUICK else 30
GRID_REPS = 400 if QUICK else 4_000
TAU = 0.05
W = Path(__file__).resolve().parents[1]
OUT = W / "deposits" / ("dt2ml_quick" if QUICK else "dt2ml")
RUNS = [(track, m, init) for track in ("unrestricted", "curved") for m in (63, 48, 32) for init in ("uniform", "peaked")]


def qhat(n, a=0.05):
    ref = np.load(W / "deposits" / "dt1b" / f"ref_n{n}.npy")
    return float(ref[int(np.ceil((1 - a) * len(ref))) - 1])


def logheights(torch, th_s, th_t, m, delta):
    s = torch.nn.functional.softplus(th_s) + delta  # m steps rising into the mode
    t = torch.nn.functional.softplus(th_t) + delta  # K-1-m steps falling after the mode
    left = -torch.flip(torch.cumsum(torch.flip(s, [0]), 0), [0])  # l_j = -sum_{i=j}^{m-1} s_i
    right = -torch.cumsum(t, 0)
    return torch.cat([left, torch.zeros(1), right])


def train(run):
    import torch
    torch.set_num_threads(1)
    track, m, init = run
    ri = RUNS.index(run)
    delta = 0.0 if track == "unrestricted" else 0.05
    torch.manual_seed(2026110620 + ri)
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(2026110620, spawn_key=(ri,))))
    v0 = -6.0 if init == "uniform" else 0.0
    th_s = torch.full((m,), v0, requires_grad=True)
    th_t = torch.full((K - 1 - m,), v0, requires_grad=True)
    opt = torch.optim.Adam([th_s, th_t], lr=0.05)
    q = qhat(N_TRAIN)
    sq = np.sqrt(N_TRAIN)
    hist = []
    for _ in range(ITERS):
        lh = logheights(torch, th_s, th_t, m, delta)
        logZ = torch.logsumexp(lh, 0) - np.log(K)
        p = torch.softmax(lh, 0).detach().numpy().astype(float)
        p /= p.sum()
        bins = rng.choice(K, size=(BATCH, N_TRAIN), p=p)
        x = np.sort((bins + rng.random((BATCH, N_TRAIN))) / K, axis=1)
        d = dips_of(x)
        ind = (d >= q).astype(float)
        s = 1 / (1 + np.exp(-(sq * d - sq * q) / TAU))
        counts = np.stack([np.bincount(b, minlength=K) for b in bins]).astype(np.float32)
        c = torch.from_numpy(counts)
        loglik = c @ lh - N_TRAIN * logZ
        adv = torch.from_numpy((s - s.mean()).astype(np.float32))
        loss = -(adv * loglik).mean() / N_TRAIN
        opt.zero_grad()
        loss.backward()
        opt.step()
        hist.append(float(ind.mean()))
    with torch.no_grad():
        lh = logheights(torch, th_s, th_t, m, delta).numpy().astype(float)
    p = np.exp(lh - lh.max())
    p /= p.sum()
    return run, {"probs": p.tolist(), "logh": lh.tolist(), "history": hist}


def describe(p):
    p = np.asarray(p)
    lh = np.log(p)
    m = int(np.argmax(p))
    flat = np.abs(np.diff(lh)) < 0.01
    best, cur, start, bs = 0, 0, 0, 0
    for i, f in enumerate(flat):
        if f:
            if cur == 0:
                start = i
            cur += 1
            if cur > best:
                best, bs = cur, start
        else:
            cur = 0
    run_mass = float(p[bs:bs + best + 1].sum()) if best else 0.0
    return {"mode_bin": m, "mode_bin_mass": float(p[m]), "longest_flat_run_bins": int(best + 1 if best else 0),
            "longest_flat_run_mass": run_mass, "max_over_min_density": float(p.max() / p.min())}


def reest(job):
    tag, p, n, key = job
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(2026110621, spawn_key=key)))
    q = qhat(n)
    hits, done = 0, 0
    b = max(100, min(2000, 2_000_000 // n))
    while done < REEST:
        mm = min(b, REEST - done)
        bins = rng.choice(K, size=(mm, n), p=np.asarray(p))
        x = np.sort((bins + rng.random((mm, n))) / K, axis=1)
        hits += int(np.sum(dips_of(x) >= q))
        done += mm
    return tag, n, hits / REEST


def grid_cell(args):
    i, j, w, eps = args
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(2026110622, spawn_key=(i, j))))
    x = np.sort(sample(rng, (GRID_REPS, N_TRAIN), [0, 1 - eps, 1], [1 - w, w]), axis=1)
    return (w, eps), float(np.mean(dips_of(x) >= qhat(N_TRAIN)))


def reest_twoblock(job):
    w, eps, n, key = job
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(2026110621, spawn_key=key)))
    q = qhat(n)
    hits, done = 0, 0
    b = max(100, min(2000, 2_000_000 // n))
    while done < REEST:
        mm = min(b, REEST - done)
        hits += int(np.sum(dips_of(np.sort(sample(rng, (mm, n), [0, 1 - eps, 1], [1 - w, w]), axis=1)) >= q))
        done += mm
    return n, hits / REEST


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=11) as ex:
        trained = dict(ex.map(train, RUNS))
    print(f"trained {len(RUNS)} runs in {time.time() - t0:.0f}s", flush=True)
    ws = np.geomspace(0.01, 0.5, GRID)
    es = np.geomspace(0.001, 0.2, GRID)
    cells = [(i, j, w, e) for i, w in enumerate(ws) for j, e in enumerate(es) if e < w]
    with ProcessPoolExecutor(max_workers=11) as ex:
        grid = dict(ex.map(grid_cell, cells))
    (wb, eb), rb = max(grid.items(), key=lambda kv: kv[1])
    print(f"baseline grid best w={wb:.4f} eps={eb:.4f} rej={rb:.4f} ({time.time() - t0:.0f}s)", flush=True)
    jobs = []
    for ri, run in enumerate(RUNS):
        for n in (200, 1000):
            jobs.append((run, trained[run]["probs"], n, (ri, n)))
    with ProcessPoolExecutor(max_workers=11) as ex:
        re = list(ex.map(reest, jobs))
        rb2 = list(ex.map(reest_twoblock, [(wb, eb, n, (99, n)) for n in (200, 1000)]))
    out = {"tier": "DEV search; fresh-seed re-estimates", "K": K, "n_train": N_TRAIN, "iters": ITERS, "batch": BATCH,
           "runs": [], "baseline": {"best_w": wb, "best_eps": eb, "grid_rej": rb,
                                    "reest": {str(n): r for n, r in rb2}, "grid_top5": sorted(
                                        [{"w": k[0], "eps": k[1], "rej": v} for k, v in grid.items()],
                                        key=lambda d: -d["rej"])[:5]}}
    for run in RUNS:
        tr = trained[run]
        rr = {str(n): r for tag, n, r in re if tag == run}
        h = tr["history"]
        out["runs"].append({"track": run[0], "mode_bin": run[1], "init": run[2],
                            "train_rej_first10": float(np.mean(h[:10])), "train_rej_last20": float(np.mean(h[-20:])),
                            "reest": rr, "shape": describe(tr["probs"]), "probs": tr["probs"]})
        np.save(OUT / f"law_{run[0]}_m{run[1]}_{run[2]}.npy", np.asarray(tr["probs"]))
    out["seconds"] = time.time() - t0
    (OUT / "adversary.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out["baseline"], indent=1))
    for r in out["runs"]:
        print(f"{r['track']:12s} m={r['mode_bin']:2d} {r['init']:8s} train {r['train_rej_first10']:.3f}->{r['train_rej_last20']:.3f} "
              f"| fresh n200 {r['reest'].get('200', float('nan')):.4f} n1000 {r['reest'].get('1000', float('nan')):.4f} | {r['shape']}")
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
