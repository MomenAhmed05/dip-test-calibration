"""O8 (Step 7, DT5-ML-DESIGN v0.2 §2): a learned critical value for T = sqrt(n) * dip, bounded between the
uniform own-n value q_U(n) and the LF value q_LF = 0.73139526 (alpha = 0.05):
    q_theta(x) = q_U(n) + (q_LF - q_U(n)) * sigmoid(g_theta(phi(x))).
The form gives power >= LF; it gives NO level guarantee. The level is only checked by the audit (o8_audit, later).

Training (DEV): a calibrator MLP maximises smoothed power on fixed bimodal alternatives subject to smoothed size
<= ALPHA_T on (a) a frozen TRAIN null grid and (b) a population of step-density adversaries (REINFORCE, as in
ml_adversary2.py) that maximise the learned test's rejection rate. Dual ascent on one Lagrange weight per null.
Run: python code/o8_calibrator.py [--quick]    Roots: DEV training 2026110670 (quick 2026110672).
Output: deposits/o8[_quick]/calibrator.pt + train.json
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import dips_of  # noqa: E402

QUICK = "--quick" in sys.argv
W = Path(__file__).resolve().parents[1]
OUT = W / "deposits" / ("o8_quick" if QUICK else "o8")
ROOT = 2026110672 if QUICK else 2026110670
ALPHA, ALPHA_T, TAU = 0.05, 0.045, 0.05
Q_LF = 0.73139526
NS_TRAIN = (200, 800)
ITERS, B = (20, 64) if QUICK else (1500, 256)
K, N_ADV = 64, 6
PROBS = np.array([0, .01, .025, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .975, .99, 1.0])


def q_uniform(n, a=ALPHA):
    """sqrt(n) x the uniform own-n critical value from the DT5 reference (full run if present, else quick)."""
    for sub in ("dt5", "dt5_quick"):
        p = W / "deposits" / sub / "power.json"
        if p.exists():
            try:
                c = json.loads(p.read_text())["crit"]
            except (OSError, json.JSONDecodeError, KeyError):
                continue
            if str(n) in c:
                return float(np.sqrt(n) * c[str(n)][str(a)]["uniform"])
    raise RuntimeError(f"no uniform critical value for n={n}")


def features(xs):
    """Scale-free fingerprint of each sorted row: log normalised quantile-block widths per unit probability, plus log n."""
    n = xs.shape[1]
    idx = np.clip(np.round(PROBS * (n - 1)).astype(int), 0, n - 1)
    qv = xs[:, idx]
    width = np.diff(qv, axis=1) / np.diff(PROBS)
    width = np.maximum(width, 1e-12)
    lw = np.log(width)
    lw -= np.median(lw, axis=1, keepdims=True)
    return np.column_stack([np.clip(lw, -8, 8), np.full(xs.shape[0], np.log(n) / np.log(12800))]).astype(np.float32)


# ------------------------------------------------------------------ fixed laws (samplers of sorted rows)
def s_W(g, b, n, w):
    pile = g.random((b, n)) < w
    return np.where(pile, 0.995 + 0.005 * g.random((b, n)), 0.995 * g.random((b, n)))


def s_mix(g, b, n, sep, p=0.5):
    return g.standard_normal((b, n)) + sep * (g.random((b, n)) < p)


def s_twobump(g, b, n):
    u = g.random((b, n))
    z = g.random((b, n))
    return np.where(z < 0.8, u, np.where(z < 0.9, 0.25, 0.75) + 0.03 * g.standard_normal((b, n)))


TRAIN_NULLS = {"uniform": lambda g, b, n: g.random((b, n)),
               "W0.05": lambda g, b, n: s_W(g, b, n, 0.05),
               "W0.2": lambda g, b, n: s_W(g, b, n, 0.2),
               "W0.4": lambda g, b, n: s_W(g, b, n, 0.4)}
ALTS = {"mix2.5": lambda g, b, n: s_mix(g, b, n, 2.5), "mix3": lambda g, b, n: s_mix(g, b, n, 3.0),
        "mix4": lambda g, b, n: s_mix(g, b, n, 4.0), "mix3.5_0.7": lambda g, b, n: s_mix(g, b, n, 3.5, 0.3),
        "twobump": s_twobump}


def logheights(torch, th_s, th_t):
    s = torch.nn.functional.softplus(th_s)
    t = torch.nn.functional.softplus(th_t)
    left = -torch.flip(torch.cumsum(torch.flip(s, [0]), 0), [0])
    return torch.cat([left, torch.zeros(1), -torch.cumsum(t, 0)])


def main():
    import torch
    torch.set_num_threads(4)
    torch.manual_seed(ROOT)
    g = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT)))
    OUT.mkdir(parents=True, exist_ok=True)
    qU = {n: q_uniform(n) for n in NS_TRAIN}
    last = torch.nn.Linear(16, 1)
    net = torch.nn.Sequential(torch.nn.Linear(len(PROBS), 16), torch.nn.Tanh(), torch.nn.Linear(16, 16),
                              torch.nn.Tanh(), last)
    with torch.no_grad():
        last.bias.fill_(3.0)  # start near LF (conservative)
    opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    modes = np.linspace(0, K - 1, N_ADV).astype(int)
    adv = [(torch.full((m,), -4.0, requires_grad=True), torch.full((K - 1 - m,), -4.0, requires_grad=True)) for m in modes]
    adv_opt = [torch.optim.Adam(list(a), lr=0.05) for a in adv]
    nulls = list(TRAIN_NULLS) + [f"adv{i}" for i in range(N_ADV)]
    lam = dict.fromkeys(nulls, 1.0)
    hist = []
    t0 = time.time()
    for it in range(ITERS):
        n = NS_TRAIN[it % len(NS_TRAIN)]
        sq = np.sqrt(n)
        batches, advinfo = {}, {}
        for k, f in TRAIN_NULLS.items():
            batches[k] = np.sort(f(g, B, n), axis=1)
        for i, (ts, tt) in enumerate(adv):
            lh = logheights(torch, ts, tt)
            p = torch.softmax(lh, 0).detach().numpy().astype(float)
            p /= p.sum()
            bins = g.choice(K, size=(B, n), p=p)
            batches[f"adv{i}"] = np.sort((bins + g.random((B, n))) / K, axis=1)
            advinfo[i] = (lh, np.stack([np.bincount(r, minlength=K) for r in bins]).astype(np.float32))
        for k, f in ALTS.items():
            batches[k] = np.sort(f(g, B, n), axis=1)
        T = {k: torch.from_numpy((sq * dips_of(x)).astype(np.float32)) for k, x in batches.items()}
        phi = {k: torch.from_numpy(features(x)) for k, x in batches.items()}

        def qth(k, n=n, phi=phi):
            return qU[n] + (Q_LF - qU[n]) * torch.sigmoid(net(phi[k]).squeeze(1))

        size = {k: torch.sigmoid((T[k] - qth(k)) / TAU).mean() for k in nulls}
        power = torch.stack([torch.sigmoid((T[k] - qth(k)) / TAU).mean() for k in ALTS]).mean()
        loss = -power + sum(lam[k] * size[k] for k in nulls)
        opt.zero_grad()
        loss.backward()
        opt.step()
        for k in nulls:  # dual ascent on the hard (indicator) size
            hard = float(torch.mean((T[k] > qth(k)).float()))
            lam[k] = max(0.0, lam[k] + 2.0 * (hard - ALPHA_T))
        with torch.no_grad():  # adversaries: REINFORCE on the soft rejection of the current calibrator
            rew = {i: torch.sigmoid((T[f"adv{i}"] - qth(f"adv{i}")) / TAU) for i in range(N_ADV)}
        for i in range(N_ADV):
            lh, counts = advinfo[i]
            logZ = torch.logsumexp(lh, 0) - np.log(K)
            loglik = torch.from_numpy(counts) @ lh - n * logZ
            a = rew[i] - rew[i].mean()
            al = -(a * loglik).mean() / n
            adv_opt[i].zero_grad()
            al.backward()
            adv_opt[i].step()
        if it % 10 == 0 or it == ITERS - 1:
            row = {"it": it, "n": n, "power": float(power.detach()), "sec": round(time.time() - t0, 1),
                   **{f"hard_{k}": float(torch.mean((T[k] > qth(k)).float())) for k in nulls},
                   **{f"lam_{k}": lam[k] for k in nulls}}
            hist.append(row)
            print(f"[{row['sec']:6.0f}s] it {it} n {n} power {row['power']:.3f} | worst null "
                  f"{max(nulls, key=lambda k: row[f'hard_{k}'])} {max(row[f'hard_{k}'] for k in nulls):.3f}", flush=True)
    torch.save(net.state_dict(), OUT / "calibrator.pt")
    advp = []
    for ts, tt in adv:
        with torch.no_grad():
            advp.append(torch.softmax(logheights(torch, ts, tt), 0).numpy().tolist())
    (OUT / "train.json").write_text(json.dumps({"root": ROOT, "iters": ITERS, "batch": B, "ns": NS_TRAIN,
                                                "qU": qU, "q_LF": Q_LF, "alpha_t": ALPHA_T, "history": hist,
                                                "adversary_probs": advp, "probs": PROBS.tolist()}, indent=1))
    print(f"saved {OUT} in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
