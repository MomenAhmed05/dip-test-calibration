"""Step 3a (DT1-a) mechanism + convention audit. DEV-mechanism tier. Frozen design: work/DT1a-PREREG.md.

Root 2026110560; each (f, n) cell is drawn in 10 chunks with spawn key (f, n, chunk).
Run: python code/dt1a_mechanism.py
"""

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import diptest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import g_ramp, g_witness  # noqa: E402, I001
from path_smoke import vn_sample  # noqa: E402

ROOT = 2026110560
NS = (800, 3200, 12800)
REPS, CHUNKS = 10_000, 10
LAWS = {0: ("W(0.1,0.005)", 0.1, 0.995, g_witness(w=0.1)),
        1: ("W(0.4,0.005)", 0.4, 0.995, g_witness(w=0.4)),
        2: ("ramp(w=.1,eps=.01)", 0.1, 0.99, g_ramp())}
XS = (0.8, 1.0, 1.2, 1.4, 1.6)
FM = (0.220644, 0.560962, 0.806055, 0.928968, 0.978091)  # Balabdaoui & Filali 2011, Table 2
OUT = Path(__file__).resolve().parents[1] / "deposits" / "dt1a"


def dalg(x):
    return 0.0 if x[0] == x[-1] else float(diptest.dipstat(x, allow_zero=False, sort_x=False))  # type: ignore[arg-type]


def run_chunk(args):
    f, n, c = args
    _, w, L, gen = LAWS[f]
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT, spawn_key=(f, n, c))))
    m = REPS // CHUNKS
    x = np.sort(gen(rng, (m, n)), axis=1)
    rows = np.empty((m, 6))
    for i in range(m):
        xi = x[i]
        d, info = diptest.dipstat(xi, full_output=True, allow_zero=False, sort_x=False)  # type: ignore[misc]
        k = int(np.searchsorted(xi, L))
        vsh = vn_sample(xi[:k].copy(), n, 0)
        top = xi[k:]
        ttop = (len(top) / n) * dalg(top) if len(top) >= 2 else 0.0
        rows[i] = (d, vsh, ttop, n - k, info["xl"], info["xu"])
    return (f, n, c), rows


def fixtures():
    out = {}
    out["all_identical_n5"] = dalg(np.full(5, 2.0))
    out["n2"] = float(diptest.dipstat(np.array([0.0, 1.0]), allow_zero=False))  # type: ignore[arg-type]
    out["n3"] = float(diptest.dipstat(np.array([0.0, 0.4, 1.0]), allow_zero=False))  # type: ignore[arg-type]
    out["ties_2_2_1"] = float(diptest.dipstat(np.array([0, 0, 1, 1, 2.0]), allow_zero=False))  # type: ignore[arg-type]
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT, spawn_key=(99,))))
    u = np.concatenate([rng.random(40), np.ones(10)])
    out["mode_atom_uniform40+atom10"] = float(diptest.dipstat(u, allow_zero=False))  # type: ignore[arg-type]
    errs = []
    for _ in range(200):
        z = rng.random(300) ** 2
        d0 = float(diptest.dipstat(z, allow_zero=False))  # type: ignore[arg-type]
        d1 = float(diptest.dipstat(-z, allow_zero=False))  # type: ignore[arg-type]
        d2 = float(diptest.dipstat(3 * z + 7, allow_zero=False))  # type: ignore[arg-type]
        errs.append(max(abs(d0 - d1), abs(d0 - d2)))
    out["reflect_affine_max_diff"] = float(max(errs))
    ties = {}
    dt0d = Path(__file__).resolve().parents[1] / "deposits" / "dt0d"
    for n in (200, 800, 3200, 12800):
        ref = np.load(dt0d / f"ref_n{n}.npy")
        for a in (0.05, 0.01):
            q = ref[int(np.ceil((1 - a) * len(ref))) - 1]
            ties[f"n{n}_a{a}"] = int(np.sum(ref == q)) - 1
    out["ref_ties_at_q"] = ties
    return out


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(f, n, c) for f in LAWS for n in NS for c in range(CHUNKS)]
    jobs.sort(key=lambda j: -j[1])
    res = {}
    with ProcessPoolExecutor(max_workers=11) as ex:
        for key, rows in ex.map(run_chunk, jobs, chunksize=1):
            res[key] = rows
    summary = {"tier": "DEV-mechanism", "root": ROOT, "reps": REPS, "F_M": dict(zip(XS, FM, strict=True)), "cells": {}}
    for f, (name, w, _, _) in LAWS.items():
        for n in NS:
            r = np.concatenate([res[(f, n, c)] for c in range(CHUNKS)])
            np.save(OUT / f"mech_f{f}_n{n}.npy", r)
            d, vsh, ttop = r[:, 0], r[:, 1], r[:, 2]
            s = np.sqrt(n)
            z = 2 * s * d / np.sqrt(1 - w)
            y = 2 * s * vsh / np.sqrt(1 - w)
            ez = [float(np.mean(z <= x)) for x in XS]
            ey = [float(np.mean(y <= x)) for x in XS]
            dz = [a - b for a, b in zip(ez, FM, strict=True)]
            dy = [a - b for a, b in zip(ey, FM, strict=True)]
            sgn = s * (d - vsh)
            mx = np.maximum(vsh, ttop)
            summary["cells"][f"f{f}_n{n}"] = {
                "law": name, "n": n, "w": w,
                "ECDF_Z": ez, "Z_minus_FM": dz, "Z_max_abs": float(max(map(abs, dz))),
                "Z_consistent": bool(max(map(abs, dz)) <= 0.015),
                "ECDF_Y": ey, "Y_minus_FM": dy, "Y_max_abs": float(max(map(abs, dy))),
                "Y_consistent": bool(max(map(abs, dy)) <= 0.015),
                "P(|D-Vsh|<=1e-12)": float(np.mean(np.abs(d - vsh) <= 1e-12)),
                "P(D>Vsh+1e-12)": float(np.mean(d > vsh + 1e-12)),
                "P(D<Vsh-1e-12)": float(np.mean(d < vsh - 1e-12)),
                "signed_sqrtn(D-Vsh)_q01_q50_q99_max": [float(v) for v in np.quantile(sgn, (0.01, 0.5, 0.99, 1.0))],
                "P(|D-max(Vsh,Ttop)|<=1e-12)": float(np.mean(np.abs(d - mx) <= 1e-12)),
                "sqrtn|D-max(Vsh,Ttop)|_q99_max": [float(v) for v in np.quantile(s * np.abs(d - mx), (0.99, 1.0))],
                "P(Ttop>Vsh)": float(np.mean(ttop > vsh)),
                "mean_k_top": float(r[:, 3].mean()),
                "modal_xl_q05_q50": [float(v) for v in np.quantile(r[:, 4], (0.05, 0.5))],
            }
    readings = {}
    for f, (name, *_rest) in LAWS.items():
        c800, c12 = summary["cells"][f"f{f}_n800"], summary["cells"][f"f{f}_n12800"]
        holds = c12["Z_consistent"] and c12["Z_max_abs"] <= c800["Z_max_abs"] + 0.005
        if holds:
            readings[name] = "FORMULA HOLDS"
        elif c12["Y_consistent"]:
            readings[name] = "SHOULDER-ONLY FAILS (Y consistent, Z not)"
        else:
            readings[name] = "BOTH FAIL (audit V_sh / conventions)"
    summary["readings"] = readings
    summary["fixtures"] = fixtures()
    summary["seconds"] = time.time() - t0
    (OUT / "mechanism.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary["readings"], indent=1))
    print(json.dumps(summary["fixtures"], indent=1))
    for k, v in summary["cells"].items():
        print(k, "Zmax", round(v["Z_max_abs"], 4), "Ymax", round(v["Y_max_abs"], 4),
              "eq", round(v["P(|D-Vsh|<=1e-12)"], 4), "D>V", round(v["P(D>Vsh+1e-12)"], 4),
              "D<V", round(v["P(D<Vsh-1e-12)"], 4), "sgn", [round(t, 4) for t in v["signed_sqrtn(D-Vsh)_q01_q50_q99_max"]],
              "eqMax", round(v["P(|D-max(Vsh,Ttop)|<=1e-12)"], 4), "T>V", round(v["P(Ttop>Vsh)"], 4))
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
