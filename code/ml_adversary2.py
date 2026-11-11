"""Step 4b part 2: longer, finer adversary training (K = 256, 3000 steps, annealed tau). DEV + fresh re-estimates.

Frozen design: work/DT2ML2-PREREG.md. Roots: 2026110623 training, 2026110624 re-estimation.
Run: python code/ml_adversary2.py [--quick]
"""

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import dips_of  # noqa: E402, I001
from ml_adversary import qhat  # noqa: E402

QUICK = "--quick" in sys.argv
K = 256
N_TRAIN = 200
ITERS, BATCH = (30, 256) if QUICK else (3000, 2048)
REEST = 2_000 if QUICK else 100_000
TAU0, TAU1 = 0.2, 0.02
W = Path(__file__).resolve().parents[1]
OUT = W / "deposits" / ("dt2ml2_quick" if QUICK else "dt2ml2")
RUNS = [("unrestricted", m, "uniform") for m in (255, 224, 192, 160, 128)] + \
       [("unrestricted", 128, "asym"), ("curved", 255, "uniform")]


def logheights(torch, th_s, th_t, delta):
    s = torch.nn.functional.softplus(th_s) + delta
    t = torch.nn.functional.softplus(th_t) + delta
    left = -torch.flip(torch.cumsum(torch.flip(s, [0]), 0), [0])
    right = -torch.cumsum(t, 0)
    return torch.cat([left, torch.zeros(1), right])


def train(run):
    import torch
    torch.set_num_threads(1)
    track, m, init = run
    ri = RUNS.index(run)
    delta = 0.0 if track == "unrestricted" else 0.0125
    torch.manual_seed(2026110623 + ri)
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(2026110623, spawn_key=(ri,))))
    th_s = torch.full((m,), -6.0)
    th_t = torch.full((K - 1 - m,), -6.0 if init != "asym" else -1.0)
    th_s.requires_grad_(True)
    th_t.requires_grad_(True)
    opt = torch.optim.Adam([th_s, th_t], lr=0.03)
    q = qhat(N_TRAIN)
    sq = np.sqrt(N_TRAIN)
    hist = []
    for it in range(ITERS):
        tau = TAU0 * (TAU1 / TAU0) ** (it / max(1, ITERS - 1))
        lh = logheights(torch, th_s, th_t, delta)
        logZ = torch.logsumexp(lh, 0) - np.log(K)
        p = torch.softmax(lh, 0).detach().numpy().astype(float)
        p /= p.sum()
        bins = rng.choice(K, size=(BATCH, N_TRAIN), p=p)
        x = np.sort((bins + rng.random((BATCH, N_TRAIN))) / K, axis=1)
        d = dips_of(x)
        s = 1 / (1 + np.exp(-(sq * d - sq * q) / tau))
        counts = np.stack([np.bincount(b, minlength=K) for b in bins]).astype(np.float32)
        loglik = torch.from_numpy(counts) @ lh - N_TRAIN * logZ
        adv = torch.from_numpy((s - s.mean()).astype(np.float32))
        loss = -(adv * loglik).mean() / N_TRAIN
        opt.zero_grad()
        loss.backward()
        opt.step()
        hist.append(float(np.mean(d >= q)))
    with torch.no_grad():
        lh = logheights(torch, th_s, th_t, delta).numpy().astype(float)
    p = np.exp(lh - lh.max())
    p /= p.sum()
    return run, {"probs": p.tolist(), "history": hist}


def reest(job):
    run, p, n, key = job
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(2026110624, spawn_key=key)))
    q = qhat(n)
    hits, done = 0, 0
    b = max(100, min(2000, 2_000_000 // n))
    while done < REEST:
        mm = min(b, REEST - done)
        bins = rng.choice(K, size=(mm, n), p=np.asarray(p))
        hits += int(np.sum(dips_of(np.sort((bins + rng.random((mm, n))) / K, axis=1)) >= q))
        done += mm
    return run, n, hits / REEST


def shape(p):
    p = np.asarray(p)
    m = int(np.argmax(p))
    lh = np.log(p)
    flat = np.abs(np.diff(lh)) < 0.005
    runs, cur, st = [], 0, 0
    for i, f in enumerate(flat):
        if f:
            if cur == 0:
                st = i
            cur += 1
        elif cur:
            runs.append((cur + 1, st))
            cur = 0
    if cur:
        runs.append((cur + 1, st))
    big = max(runs, default=(0, 0))
    top = np.sort(p)[::-1]
    return {"mode_bin": m, "mode_x": (m + 0.5) / K, "longest_flat_bins": int(big[0]),
            "longest_flat_mass": float(p[big[1]:big[1] + big[0]].sum()) if big[0] else 0.0,
            "top3_bin_mass": float(top[:3].sum()), "mass_left_of_mode": float(p[:m].sum()),
            "mass_right_of_mode": float(p[m + 1:].sum()), "max_over_min": float(p.max() / p.min())}


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=7) as ex:
        trained = dict(ex.map(train, RUNS))
    print(f"trained {len(RUNS)} runs in {time.time() - t0:.0f}s", flush=True)
    jobs = [(run, trained[run]["probs"], n, (ri, n)) for ri, run in enumerate(RUNS) for n in (200, 1000)]
    with ProcessPoolExecutor(max_workers=11) as ex:
        re = list(ex.map(reest, jobs))
    base = 0.25749
    out = {"tier": "DEV search; fresh re-estimates", "K": K, "iters": ITERS, "batch": BATCH, "runs": []}
    for run in RUNS:
        h = trained[run]["history"]
        rr = {str(n): r for (rn, n, r) in re if rn == run}
        np.save(OUT / f"law_{run[0]}_m{run[1]}_{run[2]}.npy", np.asarray(trained[run]["probs"]))
        out["runs"].append({"track": run[0], "mode_bin": run[1], "init": run[2],
                            "train_first50": float(np.mean(h[:50])), "train_last200": float(np.mean(h[-200:])),
                            "history_every100": [float(np.mean(h[i:i + 100])) for i in range(0, len(h), 100)],
                            "reest": rr, "shape": shape(trained[run]["probs"])})
    best = max(out["runs"], key=lambda r: r["reest"]["200"])
    se = np.sqrt(best["reest"]["200"] * (1 - best["reest"]["200"]) / REEST)
    out["Q_b2"] = "ADVERSARY BEATS 2-BLOCK" if best["reest"]["200"] - base > 2 * se else "2-BLOCK NOT BEATEN"
    asym = [r for r in out["runs"] if r["init"] == "asym"][0]
    out["Q_c2"] = "CENTRAL PEAK CAN FAIL" if asym["reest"]["200"] > 0.075 else "CENTRAL PEAK STAYS <= 7.5%"
    cur = [r for r in out["runs"] if r["track"] == "curved"][0]
    out["Q_d2"] = "FLATS NEEDED" if cur["reest"]["200"] <= 0.055 else "CURVED EXCEEDS 5.5%"
    out["seconds"] = time.time() - t0
    (OUT / "adversary2.json").write_text(json.dumps(out, indent=1))
    print(out["Q_b2"], "|", out["Q_c2"], "|", out["Q_d2"])
    for r in out["runs"]:
        print(f"{r['track']:12s} m={r['mode_bin']:3d} {r['init']:7s} train {r['train_first50']:.3f}->{r['train_last200']:.3f} "
              f"| fresh n200 {r['reest']['200']:.4f} n1000 {r['reest']['1000']:.4f} | {r['shape']}")
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
