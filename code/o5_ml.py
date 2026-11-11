# pyright: reportMissingImports=false
# (clustpy is installed only in the isolated venv work/venv_o5)
"""O5-ML (Step 7, DT5-ML-DESIGN v0.2): calibration-substitution sensitivity of ClustPy dip-based clustering.

Same data, same algorithm seeds, two arms: the stock uniform 'table' p-value (clustpy.utils.dip_pval) and the
LF-asymptotic p-value 1 - F_M(2 sqrt(n) D_n) (F_M only from the published table, code/fm_table.py).
The swap is a monkeypatch of `dip_pval` inside the method modules; every call is counted and logged.

Run with the isolated venv:
  python code/o5_ml.py [--quick]
Outputs: deposits/o5ml[_quick]/cells/*.npz (per cell) and summary.json. Resumable (existing cells skipped).
Roots: CONFIRM data 2026110660 / init 2026110661; DEV (--quick) 2026110662 / 2026110663.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
QUICK = "--quick" in sys.argv
ROOT_DATA, ROOT_INIT = (2026110662, 2026110663) if QUICK else (2026110660, 2026110661)
OUT = HERE.parent / "deposits" / ("o5ml_quick" if QUICK else "o5ml")
NS = (500, 2000)
R_MAIN = 60 if QUICK else 2000
R_DIPMEANS = 30 if QUICK else 1000

# law id -> (kind, parameter); NULL = certified unimodal; PRICE = not unimodal (power the LF arm gives up)
LAWS = [
    ("U", "null", ("uniform", None)),
    ("G", "null", ("normal", None)),
    ("N1", "null", ("W", 0.1)),
    ("N2_s0.0005", "null", ("clipjit", 0.0005)),
    ("N2_s0.002", "null", ("clipjit", 0.002)),
    ("N2_s0.01", "null", ("clipjit", 0.01)),
    ("N2_s0.05", "null", ("clipjit", 0.05)),
    ("N2pin_k0.85", "null", ("pin", 0.85)),
    ("N2pin_k1.4", "null", ("pin", 1.4)),
    ("N2pin_k2.5", "null", ("pin", 2.5)),
    ("R0_s0.002", "null", ("relujit", 0.002)),
    ("R0_s0.05", "null", ("relujit", 0.05)),
    ("A_N2", "atom", ("clipjit", 0.0)),
    ("A_R0", "atom", ("relujit", 0.0)),
    ("P_relu0.5", "price", ("relu_shift", 0.5)),
    ("P_mix2.5", "price", ("mix", 2.5)),  # equal-weight N(0,1)/N(s,1) is bimodal iff s > 2 (s = 2 was unimodal: removed)
    ("P_mix3", "price", ("mix", 3.0)),
    ("P_mix4", "price", ("mix", 4.0)),
    ("P_tanh", "price", ("tanh", 3.0)),
]


def sample_1d(kind, par, n, g):
    if kind == "uniform":
        return g.random(n)
    if kind == "normal":
        return g.standard_normal(n)
    if kind == "W":  # mass par on (0.995, 1), the rest uniform on (0, 0.995)
        pile = g.random(n) < par
        return np.where(pile, 0.995 + 0.005 * g.random(n), 0.995 * g.random(n))
    if kind == "pin":  # local sequence: pile-up mass kappa / sqrt(n)
        return _w(g, n, min(par / np.sqrt(n), 1.0))
    if kind == "clipjit":  # clip(U(0, 1.2), max 1) + N(0, s^2)
        return np.minimum(1.2 * g.random(n), 1.0) + par * g.standard_normal(n)
    if kind == "relujit":  # ReLU(N(0, 1)) + N(0, s^2)
        return np.maximum(g.standard_normal(n), 0.0) + par * g.standard_normal(n)
    if kind == "relu_shift":  # ReLU(N(0.5, 1)) + N(0, 0.002^2): NOT unimodal (atom at 0, density rising to 0.5)
        return np.maximum(g.standard_normal(n) + par, 0.0) + 0.002 * g.standard_normal(n)
    if kind == "mix":
        return g.standard_normal(n) + par * (g.random(n) < 0.5)
    if kind == "tanh":  # tanh of a flat input on (-par, par): density rises towards both ends
        return np.tanh(par * (2 * g.random(n) - 1))
    raise ValueError(kind)


def _w(g, n, w):
    pile = g.random(n) < w
    return np.where(pile, 0.995 + 0.005 * g.random(n), 0.995 * g.random(n))


# ---------------------------------------------------------------- calibration swap (worker side)
CALLS = []


def install(arm):
    import clustpy.centroid.dipmeans as cdm
    import clustpy.density.skinnydip as csd
    from clustpy.utils import dip_pval as p_uniform
    from fm_table import FM_pub

    def p_lf(dip_value, n_points, *a, **k):
        return float(1.0 - FM_pub(2.0 * np.sqrt(n_points) * dip_value, clip_tails=True))

    base = p_uniform if arm == "uniform" else p_lf

    def logged(dip_value, n_points, *a, **k):
        p = base(dip_value, n_points, *a, **k)
        CALLS.append((int(n_points), float(dip_value), float(p)))
        return p

    cdm.dip_pval = logged
    csd.dip_pval = logged


def n_clusters(labels):
    lab = np.asarray(labels)
    return int(np.unique(lab[lab >= 0]).size), float(np.mean(lab < 0))


def run_cell(cell):
    """cell = (method, law_index, n, alpha). Returns per-rep arrays for both arms."""
    from clustpy.centroid import DipMeans
    from clustpy.density import SkinnyDip, UniDip
    from clustpy.utils import dip_pval as p_uniform
    from clustpy.utils import dip_test
    from fm_table import FM_pub

    method, li, n, alpha = cell
    name, role, (kind, par) = LAWS[li]
    R = R_DIPMEANS if method == "dipmeans" else R_MAIN
    mi = ("unidip", "skinnydip", "dipmeans").index(method)
    rec = {k: [] for k in ("k_U", "k_LF", "noise_U", "noise_LF", "root_pU", "root_pLF",
                           "calls_U", "calls_LF", "minn_U", "minn_LF")}
    for r in range(R):
        g = np.random.default_rng(np.random.SeedSequence([ROOT_DATA, li, n, r]))
        x1 = sample_1d(kind, par, n, g)
        X = x1.reshape(-1, 1) if method == "unidip" else np.column_stack([x1] + [g.standard_normal(n) for _ in range(4)])
        seed = int(np.random.SeedSequence([ROOT_INIT, mi, li, n, r]).generate_state(1)[0] % (2 ** 31 - 1))
        d = dip_test(x1)
        rec["root_pU"].append(float(p_uniform(d, n)))
        rec["root_pLF"].append(float(1.0 - FM_pub(2.0 * np.sqrt(n) * d, clip_tails=True)))
        for arm, tag in (("uniform", "U"), ("lf", "LF")):
            install(arm)
            CALLS.clear()
            if method == "unidip":
                m = UniDip(significance=alpha, random_state=seed).fit(X)
            elif method == "skinnydip":
                m = SkinnyDip(significance=alpha, random_state=seed).fit(X)
            else:
                m = DipMeans(significance=alpha, random_state=seed).fit(X)
            k, noise = n_clusters(m.labels_)
            rec[f"k_{tag}"].append(k)
            rec[f"noise_{tag}"].append(noise)
            rec[f"calls_{tag}"].append(len(CALLS))
            rec[f"minn_{tag}"].append(min((c[0] for c in CALLS), default=0))
    return cell, {k: np.asarray(v) for k, v in rec.items()}


def cells():
    out = []
    for li in range(len(LAWS)):
        for n in NS:
            out.append(("unidip", li, n, 0.05))
            out.append(("unidip", li, n, 0.01))
            out.append(("skinnydip", li, n, 0.05))
        out.append(("dipmeans", li, 500, 0.001))
    return out


def cell_path(c):
    return OUT / "cells" / f"{c[0]}_L{c[1]:02d}_n{c[2]}_a{c[3]}.npz"


def arg_value(flag):
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else None


def run_shard(i, k):
    """No process pool (Windows pool start-up is flaky in this shell: BrokenProcessPool / WinError 6).
    Cells with index % k == i run sequentially in this process; start k processes from the shell."""
    (OUT / "cells").mkdir(parents=True, exist_ok=True)
    todo = [c for j, c in enumerate(cells()) if j % k == i and not cell_path(c).exists()]
    t0 = time.time()
    print(f"shard {i}/{k}: {len(todo)} cells, root {ROOT_DATA}/{ROOT_INIT}, quick={QUICK}", flush=True)
    for m, c in enumerate(todo, 1):
        c, rec = run_cell(c)
        tmp = cell_path(c).with_suffix(".tmp.npz")
        np.savez(tmp, allow_pickle=False, **rec)
        tmp.replace(cell_path(c))
        print(f"[{time.time() - t0:7.0f}s] {m}/{len(todo)} {c}", flush=True)


def main():
    shard = arg_value("--shard")
    if shard is not None:
        i, k = (int(v) for v in shard.split("/"))
        run_shard(i, k)
        return
    missing = [c for c in cells() if not cell_path(c).exists()]
    if missing:
        print(f"{len(missing)} cells missing; run the shards first (--shard i/k)")
        return
    t0 = time.time()
    summary = []
    for c in cells():
        z = np.load(cell_path(c))
        kU, kL = z["k_U"] > 1, z["k_LF"] > 1
        summary.append({"method": c[0], "law": LAWS[c[1]][0], "role": LAWS[c[1]][1], "n": c[2], "alpha": c[3],
                        "R": int(kU.size), "split_U": float(kU.mean()), "split_LF": float(kL.mean()),
                        "U_only": float(np.mean(kU & ~kL)), "LF_only": float(np.mean(kL & ~kU)),
                        "root_rej_U": float(np.mean(z["root_pU"] < c[3])),
                        "root_rej_LF": float(np.mean(z["root_pLF"] < c[3])),
                        "noise_U": float(z["noise_U"].mean()), "noise_LF": float(z["noise_LF"].mean()),
                        "calls_U": float(z["calls_U"].mean()), "calls_LF": float(z["calls_LF"].mean()),
                        "minn_U_median": float(np.median(z["minn_U"])), "minn_LF_median": float(np.median(z["minn_LF"]))})
    (OUT / "summary.json").write_text(json.dumps({"roots": [ROOT_DATA, ROOT_INIT], "quick": QUICK,
                                                  "laws": LAWS, "rows": summary}, indent=1))
    for s in summary:
        print(f"{s['method']:9s} {s['law']:12s} {s['role']:5s} n={s['n']:5d} a={s['alpha']:<5} "
              f"split U {100 * s['split_U']:5.1f}% LF {100 * s['split_LF']:5.1f}% | root U {100 * s['root_rej_U']:5.1f}% "
              f"LF {100 * s['root_rej_LF']:5.1f}% | calls {s['calls_U']:.1f}/{s['calls_LF']:.1f} minn {s['minn_U_median']:.0f}")
    print(f"done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
