"""O8 audit (Step 7): level and power of the FROZEN learned calibrator from o8_calibrator.py, on AUDIT laws that
were not used in training (held out by mechanism and parameter region; DT5-ML-DESIGN v0.2 §2).

Reports, per (law, n): rejection rate of T = sqrt(n) * dip against
  uniform q_U(n) | LF 0.73139526 | constant interpolations lambda in {0.25, 0.5, 0.75} | learned q_theta(x),
all on the same samples, with one-sided Clopper-Pearson upper bounds (Bonferroni over the null cells).
No selection is made on audit data: the whole lambda curve is reported.
Run: python code/o8_audit.py [--quick]     Root: AUDIT 2026110671 (quick: DEV 2026110673).
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import beta

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import dips_of  # noqa: E402, I001
from o8_calibrator import ALPHA, ALTS, PROBS, Q_LF, features, q_uniform, s_mix, s_W  # noqa: E402

QUICK = "--quick" in sys.argv
W = Path(__file__).resolve().parents[1]
SRC = W / "deposits" / ("o8_quick" if QUICK else "o8")
ROOT = 2026110673 if QUICK else 2026110671
NS = (200, 800) if QUICK else (200, 800, 3200)
R, CH = (2000, 1000) if QUICK else (50000, 2000)  # 20,000 -> 50,000 (pass probability at 5%)
LAMBDAS = (0.25, 0.5, 0.75)


def s_pin(g, b, n, kappa):
    return s_W(g, b, n, min(kappa / np.sqrt(n), 1.0))


def s_ramp(g, b, n):
    """Flat on (0, 0.9) then a linear ramp up to the right end (the DT1 'ramp' shape, mass 0.1 on the ramp)."""
    u = g.random((b, n))
    z = g.random((b, n))
    return np.where(z < 0.9, 0.9 * u, 0.9 + 0.1 * np.sqrt(u))


AUDIT_NULLS = {"W0.02": lambda g, b, n: s_W(g, b, n, 0.02), "W0.1": lambda g, b, n: s_W(g, b, n, 0.1),
               "pin0.85": lambda g, b, n: s_pin(g, b, n, 0.85), "pin1.4": lambda g, b, n: s_pin(g, b, n, 1.4),
               "pin2.5": lambda g, b, n: s_pin(g, b, n, 2.5), "ramp": lambda g, b, n: s_ramp(g, b, n),
               "clipjit0.002": lambda g, b, n: np.minimum(1.2 * g.random((b, n)), 1) + 0.002 * g.standard_normal((b, n)),
               "clipjit0.05": lambda g, b, n: np.minimum(1.2 * g.random((b, n)), 1) + 0.05 * g.standard_normal((b, n)),
               "relu0": lambda g, b, n: np.maximum(g.standard_normal((b, n)), 0) + 0.002 * g.standard_normal((b, n)),
               "normal": lambda g, b, n: g.standard_normal((b, n)),
               "beta22": lambda g, b, n: g.beta(2, 2, (b, n)),
               "uniform": lambda g, b, n: g.random((b, n))}
AUDIT_ALTS = dict(ALTS, **{"mix2.2": lambda g, b, n: s_mix(g, b, n, 2.2)})


def main():
    import torch
    torch.set_num_threads(2)
    last = torch.nn.Linear(16, 1)
    net = torch.nn.Sequential(torch.nn.Linear(len(PROBS), 16), torch.nn.Tanh(), torch.nn.Linear(16, 16),
                              torch.nn.Tanh(), last)
    net.load_state_dict(torch.load(SRC / "calibrator.pt"))
    net.eval()
    rows, t0 = [], time.time()
    laws = [(k, f, "null") for k, f in AUDIT_NULLS.items()] + [(k, f, "alt") for k, f in AUDIT_ALTS.items()]
    for li, (name, f, role) in enumerate(laws):
        for n in NS:
            qU = q_uniform(n)
            g = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT, spawn_key=(li, n))))
            T, Q = [], []
            for _ in range(R // CH):
                x = np.sort(f(g, CH, n), axis=1)
                T.append(np.sqrt(n) * dips_of(x))
                with torch.no_grad():
                    s = torch.sigmoid(net(torch.from_numpy(features(x))).squeeze(1)).numpy()
                Q.append(qU + (Q_LF - qU) * s)
            T, Q = np.concatenate(T), np.concatenate(Q)
            dec = {"uniform": T > qU, "LF": T > Q_LF, "learned": T > Q}  # noqa: SIM300
            dec.update({f"lam{lam}": T > qU + lam * (Q_LF - qU) for lam in LAMBDAS})  # noqa: SIM300
            rej = {k: int(v.sum()) for k, v in dec.items()}
            # paired discordant counts learned vs each constant baseline (gain over the simplest passing one)
            disc = {k: [int(np.sum(dec["learned"] & ~dec[k])), int(np.sum(dec[k] & ~dec["learned"]))]
                    for k in dec if k not in ("learned", "uniform")}
            rows.append({"law": name, "role": role, "n": n, "R": R, "rej": rej, "disc_learned_vs": disc,
                         "learned_q_mean": float(Q.mean()), "learned_q_p05": float(np.quantile(Q, 0.05))})
            print(f"[{time.time() - t0:6.0f}s] {role:4s} {name:13s} n={n:5d} " +
                  " ".join(f"{k} {100 * v / R:5.2f}%" for k, v in rej.items()) + f" | q_mean {Q.mean():.3f}", flush=True)
    nulls = [r for r in rows if r["role"] == "null"]
    m = len(nulls) * len(nulls[0]["rej"])  # Bonferroni over method x null cell (6 x 36 = 216)
    for r in nulls:  # simultaneous one-sided 95% CP bounds
        r["ub_simul"] = {k: float(beta.ppf(1 - 0.05 / m, v + 1, R - v)) if v < R else 1.0 for k, v in r["rej"].items()}
        r["lb_simul"] = {k: float(beta.ppf(0.05 / m, v, R - v + 1)) if v > 0 else 0.0 for k, v in r["rej"].items()}
    tol = ALPHA + 0.005
    summary = {}
    for k in nulls[0]["rej"]:
        w = max(nulls, key=lambda r, k=k: r["rej"][k])
        if all(r["ub_simul"][k] <= tol for r in nulls):
            label = "certified within tolerance"
        elif any(r["lb_simul"][k] > tol for r in nulls):
            label = "demonstrated above tolerance"
        else:
            label = "inconclusive"
        summary[k] = {"worst_law": w["law"], "n": w["n"], "rate": w["rej"][k] / R, "ub_simul": w["ub_simul"][k],
                      "label": label, "bonferroni_m": m}
    (SRC / "audit.json").write_text(json.dumps({"root": ROOT, "R": R, "ns": NS, "rows": rows, "worst": summary}, indent=1))
    for k, s in summary.items():
        print(f"WORST {k:8s}: {s['worst_law']} n={s['n']} {100 * s['rate']:.2f}% (simul UB {100 * s['ub_simul']:.2f}%)"
              f" -> {s['label']}")


if __name__ == "__main__":
    main()
