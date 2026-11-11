"""Step 7 (DT5) size and power study for dip calibrations (simulator). Design: work/DT5-PREREG.md.

One statistic (the dip D_n, diptest, allow_zero=False) at several critical values, on the SAME datasets:
  uniform  : own-n uniform reference quantile (the procedure being repaired; not level-valid; never ranked)
  LF       : c_alpha / sqrt(n), c_alpha = F_M^{-1}(1 - alpha) / 2 (sharp pointwise asymptotic over step laws)
  kuiper2  : conservative order-statistic upper bound (one-sided 1 - GAMMA) on the own-n Kuiper quantile, / 2
             (finite-n valid for continuous unimodal F, Davies & Kovac 2004 s3.6; GAMMA = calibration failure budget)
  dkw      : sqrt(log(2/alpha) / (2n)) (DKW-Massart)
Secondary diagnostic (appendix only): a dip-selected single-shoulder excursion statistic inspired by
Hartigan (2000), code/hartigan2000.py, with its asymptotic Kennedy reference (sse_asymptotic) and with own-n
uniform critical values (sse_uniformcal).

Outputs: deposits/dt5[_quick]/parts/*.npy (one per job; a restart skips finished jobs), ref_n*.npy, law*_n*.npy,
power.json (critical values, manifest, versions, point rates). Intervals and decisions: code/dt5_analysis.py.
Run: python code/dt5_power.py [--quick]
"""

import contextlib
import json
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parent))
import diptest  # noqa: E402
from dt0c_screen import dips_of  # noqa: E402, I001
from hartigan2000 import KENNEDY_Q, stats_full_of  # noqa: E402

QUICK = "--quick" in sys.argv
ROOT_REF, ROOT_LAW = (2026110652, 2026110653) if QUICK else (2026110650, 2026110651)
NS = (50, 100, 200, 400, 800, 1600, 3200, 6400, 12800)
ALPHAS = (0.05, 0.01)
REF_R, LAW_R, CHUNK = (4_000, 400, 200) if QUICK else (1_000_000, 50_000, 2_000)
H_REF_R = 2_000 if QUICK else 200_000  # reference samples that also get the shoulder statistic
GAMMA = 1e-3  # per (n, alpha) failure probability of the conservative Kuiper cutoff (18 cutoffs: total <= 0.018)
C_LF = {0.05: 0.73139526, 0.01: 0.85987427}  # Balabdaoui-Filali Table 5 halves; checked in DT4-CERT
OUT = Path(__file__).resolve().parents[1] / "deposits" / ("dt5_quick" if QUICK else "dt5")
PARTS = OUT / "parts"

DELTAS = (0.8, 1.0, 1.05, 1.1, 1.25, 1.5, 2.0, 2.5, 3.0)
ASYM = [(w, d) for w in (0.1, 0.25) for d in (1.25, 1.5, 2.0, 2.5, 3.0)]
VALLEY = [(rho, s) for rho in (0.9, 0.75, 0.5, 0.25) for s in (1.0, 0.5)]
NULLS = ["uniform", "W(0.1,.005)", "normal", "beta(2,2)", "pmix(.9U+.1Beta(.5,1))",
         "W(0.02,.001)", "pin(1.4,.1)", "pin(2.5,1)",  # added 2026-11-06 ~03:20Z
         "smoothW(0.02,.001,1e-4)"]  # added (Gaussian-smoothed near-LF step null)
SATELLITE = [(0.9, 0.0, 1.0, 0.1, 3.0, 0.25)]  # 0.9 N(0,1) + 0.1 N(3, 0.25^2)


def two_step(rng, size, w, eps):
    """W(w, eps): flat shoulder [0, 1-eps] with mass 1-w, top block (1-eps, 1] with mass w."""
    u = rng.random(size)
    return np.where(u < 1 - w, u / (1 - w) * (1 - eps), 1 - eps + (u - (1 - w)) / w * eps)


def laws():
    out = [("null", k) for k in NULLS]
    out += [("sym", d) for d in DELTAS]
    out += [("asym", wd) for wd in ASYM]
    out += [("valley", rs) for rs in VALLEY]
    out += [("sat", s) for s in SATELLITE]
    return out


LAWS = laws()


def mixture_modes(ws, mus, sds):
    """Number of local maxima of sum w_i N(mu_i, sd_i^2), by sign changes of f' on a fine grid, plus a
    near-tangency flag (a local minimum of |f'| below 1e-6 * max|f'| that is not a sign change). Numerical, labelled."""
    ws, mus, sds = (np.asarray(v, float) for v in (ws, mus, sds))
    x = np.linspace(mus.min() - 8 * sds.max(), mus.max() + 8 * sds.max(), 2_000_001)
    z = (x[:, None] - mus) / sds
    fp = np.sum(ws * (-z) * np.exp(-0.5 * z ** 2) / sds ** 2, axis=1)
    s = np.sign(fp)
    maxima = int(np.sum((s[:-1] > 0) & (s[1:] <= 0)))
    a = np.abs(fp)
    locmin = (a[1:-1] < a[:-2]) & (a[1:-1] < a[2:]) & (s[:-2] == s[2:])
    near = bool(np.any(a[1:-1][locmin] < 1e-6 * a.max()))
    return maxima, near


def pop_modes(law):
    """(mode count, how it was obtained)."""
    kind, p = law
    if kind == "asym":
        w, d = p
        m, near = mixture_modes((w, 1 - w), (-d, d), (1, 1))
        return m, "grid f' sign changes" + (" NEAR-TANGENT" if near else "")
    if kind == "sat":
        w1, m1, s1, w2, m2, s2 = p
        m, near = mixture_modes((w1, w2), (m1, m2), (s1, s2))
        return m, "grid f' sign changes" + (" NEAR-TANGENT" if near else "")
    if kind == "sym":
        return (2 if float(p) > 1 else 1), "closed form (delta > 1)"
    if kind == "valley":
        rho, s = p
        return (2 if float(rho) < min(1.0, float(s)) else 1), "closed form (rho < min(1, s))"
    if p.startswith("smoothW"):
        return 1, "proof: f' changes sign once"
    return 1, "null by construction"


def sample(rng, law, size):
    kind, p = law
    if kind == "null":
        if p == "uniform":
            return rng.random(size)
        if p == "W(0.1,.005)":
            return two_step(rng, size, 0.1, 0.005)
        if p == "W(0.02,.001)":
            return two_step(rng, size, 0.02, 0.001)
        if p == "smoothW(0.02,.001,1e-4)":
            return two_step(rng, size, 0.02, 0.001) + 1e-4 * rng.standard_normal(size)
        if p.startswith("pin("):
            # pinning sequence (DT2b / DT3c): top width lam/sqrt(n), excess k/sqrt(n); the law depends on n
            k, lam = (float(v) for v in p[4:-1].split(","))
            n = size[1]
            eps, e = lam / np.sqrt(n), k / np.sqrt(n)
            return two_step(rng, size, eps + e * (1 - eps), eps)
        if p == "normal":
            return rng.standard_normal(size)
        if p == "beta(2,2)":
            return rng.beta(2, 2, size)
        if p == "pmix(.9U+.1Beta(.5,1))":
            return np.where(rng.random(size) < 0.9, rng.random(size), rng.beta(0.5, 1.0, size))
    if kind == "sym":
        sign = np.where(rng.random(size) < 0.5, -1.0, 1.0)
        return sign * p + rng.standard_normal(size)
    if kind == "asym":
        w, d = p
        sign = np.where(rng.random(size) < w, -1.0, 1.0)
        return sign * d + rng.standard_normal(size)
    if kind == "sat":
        w1, m1, s1, _, m2, s2 = p
        first = rng.random(size) < w1
        return np.where(first, m1 + s1 * rng.standard_normal(size), m2 + s2 * rng.standard_normal(size))
    if kind == "valley":
        rho, s = p
        widths = np.array([0.4, 0.2, 0.4])
        heights = np.array([1.0, rho, s])
        masses = widths * heights / np.sum(widths * heights)
        edges = np.array([0.0, 0.4, 0.6, 1.0])
        cm = np.concatenate([[0.0], np.cumsum(masses)])
        u = rng.random(size)
        k = np.clip(np.searchsorted(cm, u, side="right") - 1, 0, 2)
        return edges[k] + (edges[k + 1] - edges[k]) * (u - cm[k]) / masses[k]
    raise ValueError(law)


def kuiper_stat(u_sorted):
    n = u_sorted.shape[1]
    i = np.arange(1, n + 1)
    return np.max(i / n - u_sorted, axis=1) + np.max(u_sorted - (i - 1) / n, axis=1)


def safe_stats_full(xs):
    """(T, m, R, eligible) per row; a failing row (ties, index inside the modal interval) is NaN, never fixed."""
    try:
        return stats_full_of(xs)
    except ValueError:
        out = np.full((xs.shape[0], 4), np.nan)
        for i in range(xs.shape[0]):
            with contextlib.suppress(ValueError):
                out[i] = stats_full_of(xs[i:i + 1])[0]
        return out


def part_path(job):
    kind, li, n, c = job
    return PARTS / f"{kind}_{'x' if li is None else li}_{n}_{c}.npy"


def run_job(job):
    kind, li, n, c = job
    p = part_path(job)
    if p.exists():
        return job, "skip"
    if kind == "ref":
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_REF, spawn_key=(0, n, c))))
        u = np.sort(rng.random((CHUNK, n)), axis=1)
        th = safe_stats_full(u)[:, 0] if c < H_REF_R // CHUNK else np.full(CHUNK, np.nan)
        arr = np.stack([dips_of(u), kuiper_stat(u), th], axis=1)
    else:
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_LAW, spawn_key=(1, li, n, c))))
        xs = np.sort(sample(rng, LAWS[li], (CHUNK, n)), axis=1)
        arr = np.concatenate([dips_of(xs)[:, None], safe_stats_full(xs)], axis=1)
    tmp = p.with_suffix(".tmp.npy")
    np.save(tmp, arr)
    tmp.replace(p)
    return job, "done"


def upper_os(x_sorted, p, gamma):
    """Order statistic X_(k) with P(X_(k) >= true p-quantile) >= 1 - gamma (distribution-free)."""
    n = len(x_sorted)
    k = int(binom.ppf(1 - gamma, n, p)) + 1
    return float(x_sorted[min(k, n) - 1])


def main():
    t0 = time.time()
    PARTS.mkdir(parents=True, exist_ok=True)
    jobs = [("ref", None, n, c) for n in NS for c in range(REF_R // CHUNK)]
    jobs += [("law", li, n, c) for li in range(len(LAWS)) for n in NS for c in range(LAW_R // CHUNK)]
    jobs.sort(key=lambda j: -j[2])
    workers = 3 if QUICK else 11
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for job, status in ex.map(run_job, jobs, chunksize=2):
            done += 1
            if done % 200 == 0 or done == len(jobs):
                print(f"[{time.time() - t0:7.0f}s] {done}/{len(jobs)} jobs ({status} {job})", flush=True)
    print(f"simulated in {time.time() - t0:.0f}s", flush=True)

    crit = {}
    for n in NS:
        r = np.concatenate([np.load(part_path(("ref", None, n, c))) for c in range(REF_R // CHUNK)])
        if not QUICK:
            np.save(OUT / f"ref_n{n}.npy", r)
        dsort, ksort = np.sort(r[:, 0]), np.sort(r[:, 1])
        th = np.sort(r[:, 2][~np.isnan(r[:, 2])])
        crit[n] = {}
        for a in ALPHAS:
            crit[n][a] = {"uniform": float(np.quantile(dsort, 1 - a, method="inverted_cdf")),
                          "LF": C_LF[a] / np.sqrt(n),
                          "kuiper2": upper_os(ksort, 1 - a, GAMMA) / 2,
                          "kuiper2_point": float(np.quantile(ksort, 1 - a, method="inverted_cdf")) / 2,
                          "dkw": float(np.sqrt(np.log(2 / a) / (2 * n))),
                          "sse_uniformcal": float(np.quantile(th, 1 - a, method="inverted_cdf")),
                          "sse_asymptotic": KENNEDY_Q[a],
                          "ref_shoulder_failed": int(np.sum(np.isnan(r[:H_REF_R, 2])))}
    rows, manifest = [], []
    for li, law in enumerate(LAWS):
        modes, how = pop_modes(law)
        manifest.append({"index": li, "law": [law[0], law[1]], "pop_modes": modes, "modes_from": how})
        for n in NS:
            arr = np.concatenate([np.load(part_path(("law", li, n, c))) for c in range(LAW_R // CHUNK)])
            if not QUICK:
                np.save(OUT / f"law{li}_n{n}.npy", arr)
            d, th = arr[:, 0], arr[:, 1]
            row = {"law": [law[0], law[1]], "n": n, "pop_modes": modes,
                   "sse_failed_rows": int(np.sum(np.isnan(th))),
                   "sse_eligible_frac": float(np.nanmean(arr[:, 4])),
                   "sse_median_m": float(np.nanmedian(arr[:, 2])),
                   "sse_median_R": float(np.nanmedian(arr[:, 3]))}
            for a in ALPHAS:
                c = crit[n][a]
                row[str(a)] = {m: float(np.mean(d > c[m])) for m in ("uniform", "LF", "kuiper2", "kuiper2_point", "dkw")}
                row[str(a)]["sse_asymptotic"] = float(np.mean(th > c["sse_asymptotic"]))
                row[str(a)]["sse_uniformcal"] = float(np.mean(th > c["sse_uniformcal"]))
            rows.append(row)
    out = {"tier": "QUICK" if QUICK else "CONFIRM-DT5", "roots": [ROOT_REF, ROOT_LAW],
           "reps": {"ref": REF_R, "law": LAW_R, "ref_shoulder": H_REF_R}, "gamma_kuiper": GAMMA,
           "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                        "diptest": getattr(diptest, "__version__", "unknown")},
           "manifest": manifest,
           "crit": {str(n): {str(a): crit[n][a] for a in ALPHAS} for n in NS}, "rows": rows,
           "seconds": time.time() - t0}
    (OUT / "power.json").write_text(json.dumps(out, indent=1))
    print("laws:", len(LAWS), "| modes:", [(m["index"], m["pop_modes"], m["modes_from"]) for m in manifest
                                         if "grid" in m["modes_from"]])
    for n in (NS[0], NS[-1]):
        print(n, {a: {m: round(np.sqrt(n) * v, 4) for m, v in crit[n][a].items() if m in ("uniform", "LF", "kuiper2", "dkw")}
                  for a in ALPHAS})
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
