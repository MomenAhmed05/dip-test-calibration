"""DT5b (post-DT5 SUPPLEMENT): size and power of the Cheng & Hall (1998) calibrated dip test
next to the uniform and LF calibrations, on the same samples, on the DT5 null laws and 5 DT5 alternatives.

Phase "table": per n and per d on ch98.D_GRID, M dips of n-samples from the Cheng-Hall calibration law
               -> upper quantiles at alpha = 0.05, 0.01                    (root ROOT_TAB, spawn (n, d_index, chunk))
Phase "laws":  per (law, n), R samples -> dip, d_hat                         (root ROOT_LAW, spawn (law_index, n, chunk))
Summary:       rejection counts for uniform (DT5 own-n reference cutoff from deposits/dt5/power.json),
               LF (0.73139526 / 0.85987427 over sqrt(n)), and CH98 (table quantile interpolated at d_hat).
Runner:        python dt5b_ch98.py [--quick] --shard i/k   (k processes started from the shell; no process pool)
               python dt5b_ch98.py [--quick]                (summary; requires all parts)
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import beta as beta_dist

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch98  # noqa: E402
from dt0c_screen import dips_of  # noqa: E402
from dt5_power import LAWS as DT5_LAWS  # noqa: E402
from dt5_power import sample  # noqa: E402

QUICK = "--quick" in sys.argv
ROOT_TAB, ROOT_LAW = (2026110682, 2026110683) if QUICK else (2026110680, 2026110681)
W = Path(__file__).resolve().parents[1]
OUT = W / "deposits" / ("dt5b_quick" if QUICK else "dt5b")
NS = (50, 200, 800, 3200, 12800)
ALPHAS = (0.05, 0.01)
C_LF = {0.05: 0.73139526, 0.01: 0.85987427}
M, R, CH = (2000, 1000, 500) if QUICK else (20000, 10000, 2000)
ALT_KEEP = [("sym", 1.5), ("sym", 2.0), ("asym", (0.25, 2.5)), ("valley", (0.5, 1.0)),
            ("sat", (0.9, 0.0, 1.0, 0.1, 3.0, 0.25))]
LAWS = [lw for lw in DT5_LAWS if lw[0] == "null"] + ALT_KEEP


def jobs():
    out = [("tab", n, j, c) for n in NS for j in range(len(ch98.D_NODES)) for c in range(M // CH)]
    out += [("law", n, li, c) for n in NS for li in range(len(LAWS)) for c in range(R // CH)]
    return sorted(out, key=lambda t: -t[1])


def part(job):
    kind, n, i, c = job
    return OUT / "parts" / f"{kind}_n{n}_{i:02d}_c{c:02d}.npz"


def run_job(job):
    kind, n, i, c = job
    if kind == "tab":
        g = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_TAB, spawn_key=(n, i, c))))
        ref = ch98.reference(float(ch98.D_NODES[i]))
        x = np.sort(ch98.sample_reference(g, ref, (CH, n)), axis=1)
        return {"dip": dips_of(x)}
    g = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_LAW, spawn_key=(i, n, c))))
    x = np.sort(sample(g, LAWS[i], (CH, n)), axis=1)
    return {"dip": dips_of(x), "dhat": np.array([ch98.d_hat(r) for r in x])}


def run_shard(s, k):
    (OUT / "parts").mkdir(parents=True, exist_ok=True)
    todo = [j for m, j in enumerate(jobs()) if m % k == s and not part(j).exists()]
    t0 = time.time()
    print(f"shard {s}/{k}: {len(todo)} jobs, roots {ROOT_TAB}/{ROOT_LAW}, quick={QUICK}", flush=True)
    for m, j in enumerate(todo, 1):
        rec = run_job(j)
        tmp = part(j).with_suffix(".tmp.npz")
        np.savez(tmp, allow_pickle=False, **rec)
        tmp.replace(part(j))
        if m % 20 == 0 or m == len(todo):
            print(f"[{time.time() - t0:7.0f}s] {m}/{len(todo)} {j}", flush=True)


def cp(k, n):
    lo = float(beta_dist.ppf(0.025, k, n - k + 1)) if k > 0 else 0.0
    hi = float(beta_dist.ppf(0.975, k + 1, n - k)) if k < n else 1.0
    return lo, hi


def summarize():
    missing = [j for j in jobs() if not part(j).exists()]
    if missing:
        raise SystemExit(f"{len(missing)} parts missing; run the shards first")
    try:
        dt5 = json.loads((W / "deposits" / "dt5" / "power.json").read_text())
    except (OSError, json.JSONDecodeError) as e:
        raise SystemExit(f"cannot read DT5 power.json: {e}") from e
    table = {}
    for n in NS:
        for a in ALPHAS:
            q = []
            for j in range(len(ch98.D_NODES)):
                d = np.concatenate([np.load(part(("tab", n, j, c)))["dip"] for c in range(M // CH)])
                q.append(float(np.quantile(d, 1 - a)))
            table[(n, a)] = np.array(q)
    rows = []
    for li, law in enumerate(LAWS):
        for n in NS:
            z = [np.load(part(("law", n, li, c))) for c in range(R // CH)]
            dip = np.concatenate([p["dip"] for p in z])
            dh = np.concatenate([p["dhat"] for p in z])
            row = {"law": law, "n": n, "R": int(dip.size), "dhat_median": float(np.nanmedian(dh)),
                   "dhat_nan": int(np.sum(~np.isfinite(dh))),
                   "frac_above_dmax": float(np.mean(dh > ch98.D_GRID[-1])),
                   "frac_below_dmin": float(np.mean(dh < ch98.D_GRID[0]))}
            for a in ALPHAS:
                cut = {"uniform": dt5["crit"][str(n)][str(a)]["uniform"], "LF": C_LF[a] / np.sqrt(n)}
                ch = ch98.crit(table[(n, a)], np.where(np.isfinite(dh), dh, 0.0))  # NaN d_hat (zero spread) -> uniform node
                k = {m: int(np.sum(dip > c)) for m, c in cut.items()}
                k["CH98"] = int(np.sum(dip > ch))
                row[str(a)] = {m: {"k": v, "rate": v / dip.size, "ci": cp(v, dip.size)} for m, v in k.items()}
                row[str(a)]["CH98_crit_median"] = float(np.median(ch))
            rows.append(row)
    out = {"tier": "DEV-quick" if QUICK else "SUPPLEMENT-DT5b", "roots": [ROOT_TAB, ROOT_LAW], "M": M, "R": R,
           "NS": NS, "d_nodes": ch98.D_NODES.tolist(),
           "table": {f"{n}/{a}": v.tolist() for (n, a), v in table.items()}, "rows": rows}
    (OUT / "dt5b.json").write_text(json.dumps(out, indent=1))
    for r in rows:
        s = " | ".join(f"{m} {100 * r['0.05'][m]['rate']:5.2f}%" for m in ("uniform", "LF", "CH98"))
        print(f"{str(r['law']):45s} n={r['n']:6d} dhat~{r['dhat_median']:8.2f} | 5%: {s}")
    print("saved", OUT / "dt5b.json")


def main():
    if "--shard" in sys.argv:
        s, k = (int(v) for v in sys.argv[sys.argv.index("--shard") + 1].split("/"))
        run_shard(s, k)
    else:
        summarize()


if __name__ == "__main__":
    main()
