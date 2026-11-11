"""O8 red team (Step 7, after DT5ML-RESULT): search for UNIMODAL laws on which the FROZEN learned calibrator
(deposits/o8/calibrator.pt, sha 67c5727d) rejects more than 5% of the time.

Two parts (DEV, root 2026110690; quick 2026110692):
  --families  fixed families, R_FAM samples per (law, n), n in NS:
              spike-on-flat at 5 locations (left end, 0.25, 0.5, 0.75, right end) x pile mass w x pile width eps;
              two-block step-up laws (level 1 on [0, a), level h on [a, 1]); pinning sequences (mass k/sqrt(n),
              width lam/sqrt(n)) at both ends; the corrected flat-plus-ramp null (ramp rising FROM the flat level).
  --shard i/k REINFORCE adversaries over unimodal step densities with K = 256 bins (training used K = 64),
              N_ADV adversaries with spread mode positions and 3 initialisations, n cycling over NS; reward =
              smoothed rejection of the frozen calibrator; then a fresh re-estimate (R_RE samples per n) of each
              adversary's final law.
Summary: python o8_redteam.py [--quick]   -> deposits/o8_redteam[_quick]/redteam.json
Every law here is unimodal by construction (step densities with one maximal block or plateau; spike on flat).
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import dips_of  # noqa: E402, I001
from o8_calibrator import PROBS, Q_LF, features, logheights, q_uniform  # noqa: E402

QUICK = "--quick" in sys.argv
W = Path(__file__).resolve().parents[1]
OUT = W / "deposits" / ("o8_redteam_quick" if QUICK else "o8_redteam")
ROOT = 2026110692 if QUICK else 2026110690
NS = (200, 800, 3200)
K = 256
N_ADV = 6 if QUICK else 24
ITERS, B = (12, 64) if QUICK else (600, 256)
R_RE, R_FAM = (2000, 2000) if QUICK else (20000, 20000)
TAU = 0.05


def load_net(torch):
    last = torch.nn.Linear(16, 1)
    net = torch.nn.Sequential(torch.nn.Linear(len(PROBS), 16), torch.nn.Tanh(), torch.nn.Linear(16, 16),
                              torch.nn.Tanh(), last)
    net.load_state_dict(torch.load(W / "deposits" / "o8" / "calibrator.pt"))
    net.eval()
    return net


def decide(torch, net, x, n):
    """Rejection decisions (learned, LF, uniform) for sorted rows x."""
    qU = q_uniform(n)
    T = np.sqrt(n) * dips_of(x)
    with torch.no_grad():
        s = torch.sigmoid(net(torch.from_numpy(features(x))).squeeze(1)).numpy()
    q = qU + (Q_LF - qU) * s
    return T > q, T > Q_LF, T > qU, T, q  # noqa: SIM300


# ------------------------------------------------------------------ fixed families
def step_sampler(edges, dens):
    """Sampler for a step density with breakpoints `edges` (increasing, from 0 to 1) and block densities `dens`."""
    edges, dens = np.asarray(edges, float), np.asarray(dens, float)
    mass = dens * np.diff(edges)
    mass = mass / mass.sum()
    cm = np.concatenate([[0.0], np.cumsum(mass)])

    def f(g, b, n):
        u = g.random((b, n))
        j = np.clip(np.searchsorted(cm, u, side="right") - 1, 0, len(dens) - 1)
        v = (u - cm[j]) / mass[j]
        return edges[j] + v * (edges[j + 1] - edges[j])
    return f


def spike(loc, w, eps):
    """Flat on [0, 1] plus extra mass w spread uniformly on a window of width eps at location loc."""
    a = min(max(loc - eps / 2, 0.0), 1.0 - eps)
    base = 1.0 - w
    return step_sampler([0.0, a, a + eps, 1.0], [base, base + w / eps, base]) if 0 < a < 1 - eps else (
        step_sampler([0.0, eps, 1.0], [base + w / eps, base]) if a <= 0 else
        step_sampler([0.0, 1.0 - eps, 1.0], [base, base + w / eps]))


def families():
    fam = {}
    for loc in (0.0, 0.25, 0.5, 0.75, 1.0):
        for w in (0.02, 0.05, 0.1, 0.2, 0.4):
            for eps in (0.0005, 0.005, 0.05):
                fam[f"spike(loc={loc},w={w},eps={eps})"] = ("fixed", spike(loc, w, eps))
    for a in (0.5, 0.8, 0.9, 0.95):
        for h in (1.2, 1.5, 2.0, 4.0):
            fam[f"twoblock(a={a},h={h})"] = ("fixed", step_sampler([0.0, a, 1.0], [1.0, h]))
    for side in ("right", "left"):
        for k in (0.85, 1.4, 2.5):
            for lam in (0.1, 1.0):
                fam[f"pin({side},k={k},lam={lam})"] = ("pin", (side, k, lam))
    ramp_ok = ("fixed", lambda g, b, n: _ramp_ok(g, b, n))
    fam["ramp_from_flat(0.9,top2)"] = ramp_ok
    return fam


def _ramp_ok(g, b, n):
    """Density 1 on (0, 0.9), then linear from 1 up to 2 on (0.9, 1) (unimodal; the intended DT1-type ramp)."""
    # mass: 0.9 on the flat, 0.15 on the ramp (integral of 1 + 10 t over t in (0, 0.1)); normalise
    pf = 0.9 / 1.05
    u, z = g.random((b, n)), g.random((b, n))
    # inverse CDF on the ramp for density proportional to 1 + 10 t, t in (0, 0.1): F(t) = (t + 5 t^2) / 0.15
    t = (-1.0 + np.sqrt(1.0 + 20.0 * 0.15 * u)) / 10.0
    return np.where(z < pf, 0.9 * u, 0.9 + t)


def sample_family(g, spec, b, n):
    kind, p = spec
    if kind == "fixed":
        return p(g, b, n)
    side, k, lam = p
    eps, e = lam / np.sqrt(n), k / np.sqrt(n)
    w = eps + e * (1 - eps)
    pile = g.random((b, n)) < w
    x = np.where(pile, 1 - eps + eps * g.random((b, n)), (1 - eps) * g.random((b, n)))
    return 1.0 - x if side == "left" else x


def run_families(torch, net):
    OUT.mkdir(parents=True, exist_ok=True)
    rows, t0 = [], time.time()
    for li, (name, spec) in enumerate(families().items()):
        for n in NS:
            g = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT, spawn_key=(1, li, n))))
            k = np.zeros(3, int)
            for _ in range(R_FAM // 2000 if not QUICK else 1):
                b = 2000
                x = np.sort(sample_family(g, spec, b, n), axis=1)
                dl, dlf, du, _, _ = decide(torch, net, x, n)
                k += [int(dl.sum()), int(dlf.sum()), int(du.sum())]
            R = (R_FAM // 2000) * 2000 if not QUICK else 2000
            rows.append({"law": name, "n": n, "R": R, "learned": int(k[0]), "LF": int(k[1]), "uniform": int(k[2])})
        top = max(rows[-len(NS):], key=lambda r: r["learned"])
        print(f"[{time.time() - t0:6.0f}s] {name:34s} worst learned {100 * top['learned'] / top['R']:5.2f}% "
              f"(n={top['n']}) LF {100 * top['LF'] / top['R']:5.2f}% U {100 * top['uniform'] / top['R']:5.2f}%",
              flush=True)
    (OUT / "families.json").write_text(json.dumps({"root": ROOT, "rows": rows}, indent=1))


# ------------------------------------------------------------------ REINFORCE adversaries
def run_adversary(torch, net, i):
    g = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT, spawn_key=(2, i))))
    torch.manual_seed(ROOT + i)
    mode = int(np.linspace(0, K - 1, N_ADV)[i])
    init = i % 3  # 0 flat, 1 random, 2 pile at the mode
    if init == 1:
        a0 = torch.tensor(g.normal(-3, 1, mode), dtype=torch.float32)
        b0 = torch.tensor(g.normal(-3, 1, K - 1 - mode), dtype=torch.float32)
    else:
        a0, b0 = torch.full((mode,), -6.0), torch.full((K - 1 - mode,), -6.0)
        if init == 2:
            if mode > 0:
                a0[-1] = 2.0
            if K - 1 - mode > 0:
                b0[0] = 2.0
    ts, tt = a0.clone().requires_grad_(True), b0.clone().requires_grad_(True)
    opt = torch.optim.Adam([ts, tt], lr=0.05)
    hist = []
    for it in range(ITERS):
        n = NS[it % len(NS)]
        lh = logheights(torch, ts, tt)
        p = torch.softmax(lh, 0).detach().numpy().astype(float)
        p /= p.sum()
        bins = g.choice(K, size=(B, n), p=p)
        x = np.sort((bins + g.random((B, n))) / K, axis=1)
        T = np.sqrt(n) * dips_of(x)
        qU = q_uniform(n)
        with torch.no_grad():
            q = qU + (Q_LF - qU) * torch.sigmoid(net(torch.from_numpy(features(x))).squeeze(1)).numpy()
        rew = 1.0 / (1.0 + np.exp(-(T - q) / TAU))
        counts = torch.from_numpy(np.stack([np.bincount(r, minlength=K) for r in bins]).astype(np.float32))
        logZ = torch.logsumexp(lh, 0) - np.log(K)
        loglik = counts @ lh - n * logZ
        adv = torch.from_numpy((rew - rew.mean()).astype(np.float32))
        loss = -(adv * loglik).mean() / n
        opt.zero_grad()
        loss.backward()
        opt.step()
        if it % 10 == 0 or it == ITERS - 1:
            hist.append({"it": it, "n": n, "hard": float(np.mean(T > q)), "soft": float(rew.mean())})  # noqa: SIM300
    with torch.no_grad():
        p = torch.softmax(logheights(torch, ts, tt), 0).numpy().astype(float)
    p /= p.sum()
    re = []
    for n in NS:  # fresh re-estimate of the final law
        g2 = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT, spawn_key=(3, i, n))))
        k = np.zeros(3, int)
        for _ in range(R_RE // 1000):
            bins = g2.choice(K, size=(1000, n), p=p)
            x = np.sort((bins + g2.random((1000, n))) / K, axis=1)
            dl, dlf, du, _, _ = decide(torch, net, x, n)
            k += [int(dl.sum()), int(dlf.sum()), int(du.sum())]
        re.append({"n": n, "R": (R_RE // 1000) * 1000, "learned": int(k[0]), "LF": int(k[1]), "uniform": int(k[2])})
    return {"adv": i, "mode_bin": mode, "init": init, "probs": p.tolist(), "hist": hist, "reestimate": re}


def run_shard(s, k):
    import torch
    torch.set_num_threads(1)
    net = load_net(torch)
    (OUT / "adv").mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    for i in range(s, N_ADV, k):
        path = OUT / "adv" / f"adv{i:02d}.json"
        if path.exists():
            continue
        rec = run_adversary(torch, net, i)
        path.write_text(json.dumps(rec))
        worst = max(rec["reestimate"], key=lambda r: r["learned"])
        print(f"[{time.time() - t0:7.0f}s] adv {i} mode {rec['mode_bin']} init {rec['init']}: re-est learned "
              f"{100 * worst['learned'] / worst['R']:.2f}% at n={worst['n']} (LF {100 * worst['LF'] / worst['R']:.2f}%,"
              f" U {100 * worst['uniform'] / worst['R']:.2f}%)", flush=True)


def summarize():
    advs = sorted(OUT.glob("adv/adv*.json"))
    rows = []
    for pth in advs:
        try:
            r = json.loads(pth.read_text())
        except (OSError, json.JSONDecodeError) as e:
            raise SystemExit(f"cannot read {pth}: {e}") from e
        for x in r["reestimate"]:
            rows.append({"adv": r["adv"], "mode_bin": r["mode_bin"], "init": r["init"], **x})
    fam = []
    fp = OUT / "families.json"
    if fp.exists():
        try:
            fam = json.loads(fp.read_text())["rows"]
        except (OSError, json.JSONDecodeError, KeyError) as e:
            raise SystemExit(f"cannot read {fp}: {e}") from e
    allrows: list[dict] = [{"src": "adversary", **r} for r in rows] + [{"src": "family", **r} for r in fam]
    worst = sorted(allrows, key=lambda r: -float(r["learned"]) / float(r["R"]))[:15]
    (OUT / "redteam.json").write_text(json.dumps({"root": ROOT, "n_adv_files": len(advs), "worst15": worst,
                                                  "rows": allrows}, indent=1))
    for r in worst:
        lab = r.get("law", f"adv{r.get('adv')}")
        rr = float(r["R"])
        print(f"{r['src']:9s} {lab:36s} n={r['n']:5d} learned {100 * float(r['learned']) / rr:5.2f}% "
              f"LF {100 * float(r['LF']) / rr:5.2f}% U {100 * float(r['uniform']) / rr:5.2f}%")


def main():
    if "--families" in sys.argv:
        import torch
        torch.set_num_threads(2)
        run_families(torch, load_net(torch))
    elif "--shard" in sys.argv:
        s, k = (int(v) for v in sys.argv[sys.argv.index("--shard") + 1].split("/"))
        run_shard(s, k)
    else:
        summarize()


if __name__ == "__main__":
    main()
