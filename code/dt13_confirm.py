"""Step 14 (DT13) CONFIRM driver, Parts W (fixed witnesses, every n) and L (per-n selected laws, grid n).
Frozen design: work/DT13-PREREG.md. Conventions: CONVENTIONS.md CV-1..CV-11 (CV-4 amendment v1.1).

Phase 1 (ref): for every n in NS_ALL, R_REF dips of iid U(0,1), root ROOT_REF, spawn (0, n), one stream per n in
  batch order. Per n and alpha: q_hat = X_(ceil((1-a)R)), bracket X_(l) <= q <= X_(u) with binomial tails <= GT
  each (GT = GAMMA / (2 * #pairs)). Saved: deposits/dt13/ref_brackets.json (+ a +-2000 order-statistic window).
Phase 2 (laws): Part W: f = 0 uniform anchor, 1 W(0.1, 0.005), 2 ramp(w=.1, eps_r=.01), 3 W(0.4, 0.005),
  4 path witness W(w_n, 1e-4), w_n = min(0.4, 1.5 / sqrt(n)); N_LAW dips each, root ROOT_W, spawn (1, f, n).
  Part L: per (n in NS_L, alpha index ia) the frozen DEV-selected W(w, eps) (SEL below), root ROOT_L,
  spawn (2, n, ia). Counts at q_hat, X_(l), X_(u) for all three alpha.
Phase 3 (audit): LP audit (lp_dip, independent implementation): for f in AUDIT_F and n in AUDIT_NS, 100 random
  indices (root ROOT_AUDIT, spawn (3, f, n)) + every sample with |D - X| < 1e-9 for any bracket value X; the
  stream is regenerated with the same seeds. Rule: |D_alg - D_LP| <= 1e-9.
Usage: python dt13_confirm.py [--quick]   (quick = tiny R/N, DEV root offsets, never a confirmation)
"""
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import cp_lower, cp_upper, dips_of, g_ramp, g_uniform, g_witness  # noqa: E402, I001
from lp_dip import dip_lp  # noqa: E402

QUICK = "--quick" in sys.argv
ROOT_REF, ROOT_W, ROOT_L, ROOT_AUDIT = (2026110750, 2026110751, 2026110752, 2026110756)
if QUICK:   # quick dry runs use DEV roots only (CONVENTIONS: 2026110746)
    ROOT_REF = ROOT_W = ROOT_L = ROOT_AUDIT = 2026110746
NS_ALL = tuple(range(4, 1001)) if not QUICK else (4, 5, 8, 20, 50, 200)
R_REF = 20_000 if QUICK else 1_000_000
N_LAW = 2_000 if QUICK else 100_000
ALPHAS = (0.01, 0.05, 0.10)
GAMMA, DELTA_W, DELTA_L = 0.002, 0.002, 0.002
NS_L = (4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500, 750, 1000)
# n = 4: the exact two-step maximiser w* = 1 - p* of f(p) = a(p^4 + 4p^3 q + q^4) + 6 p^2 q^2 / (1 + r_a)
# (DT13 Part X; sympy, 12 digits), with a thin top eps = 1e-6 (the supremum is approached as eps -> 0).
W4 = {0: 0.473789162301, 1: 0.433348270040, 2: 0.396743906214}
# Frozen Part L selections (DEV root 2026110743, deposits/dt13/dev_search.json, best per (n, alpha)):
SEL = {}           # filled from the DEV file at freeze time and written into the prereg; see load_sel()
AUDIT_F = (0, 1, 2, 3, 4)          # every Part W law (ramp added); Part L laws audited too (see main)
AUDIT_NS = (4, 5, 6, 7, 8, 10, 20, 50, 100, 200) if not QUICK else (4, 8, 20)
OUT = Path(__file__).resolve().parents[1] / "deposits" / ("dt13_quick" if QUICK else "dt13")


def w_path(n):
    return min(0.4, 1.5 / math.sqrt(n))


def law_gen(f, n):
    if f == 0:
        return g_uniform
    if f == 1:
        return g_witness(w=0.1, eps=0.005)
    if f == 2:
        return g_ramp()
    if f == 3:
        return g_witness(w=0.4, eps=0.005)
    if f == 4:
        return g_witness(w=w_path(n), eps=1e-4)
    raise ValueError(f)


LAW_NAMES = {0: "U(0,1)", 1: "W(0.1,0.005)", 2: "ramp(w=.1,eps_r=.01)", 3: "W(0.4,0.005)",
             4: "path W(min(.4,1.5/sqrt n),1e-4)"}


def load_sel():
    p = Path(__file__).resolve().parents[1] / "deposits" / "dt13" / "dev_search.json"
    try:
        dev = json.loads(p.read_text())
    except (OSError, ValueError) as err:
        raise SystemExit(f"cannot read {p}: {err}") from err
    sel = {(4, ia): (W4[ia], 1e-6) for ia in range(len(ALPHAS))}
    for n in NS_L:
        b = dev["by_n"].get(str(n))
        if b is None:
            continue
        for ia, a in enumerate(ALPHAS):
            r = b["best"][str(a)]
            sel[(n, ia)] = (float(r["w"]), float(r["eps"]))
    return sel


def stream_dips(rng, gen, m, n):
    b = max(200, min(20000, 2_000_000 // n))
    out = np.empty(m)
    for s in range(0, m, b):
        k = min(b, m - s)
        out[s:s + k] = dips_of(gen(rng, (k, n)))
    return out


def bracket_idx(R, p, g):
    """1-based order-statistic indices (l, u) with P(Bin(R,p) <= l-1) <= g and P(Bin(R,p) >= u) <= g."""
    lo = int(stats.binom.ppf(g, R, p))
    while stats.binom.cdf(lo - 1, R, p) > g:
        lo -= 1
    u = int(stats.binom.isf(g, R, p)) + 1
    while stats.binom.sf(u - 1, R, p) > g:
        u += 1
    return max(lo, 1), min(u, R)


def ref_job(args):
    n, idx = args
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_REF, spawn_key=(0, n))))
    ref = np.sort(stream_dips(rng, lambda r, s: r.random(s), R_REF, n))
    out = {}
    for a in ALPHAS:
        k = int(np.ceil((1 - a) * R_REF))
        lo, hi = idx[a]
        w0, w1 = max(0, k - 2001), min(R_REF, k + 2000)
        out[str(a)] = {"k": k, "l": lo, "u": hi, "q_hat": float(ref[k - 1]), "X_l": float(ref[lo - 1]),
                       "X_u": float(ref[hi - 1]), "ties_at_q": int(np.sum(ref == ref[k - 1])) - 1,
                       "window_start": w0 + 1, "window": ref[w0:w1].tolist() if n in NS_L or n <= 8 else None}
    out["min"] = float(ref[0])
    out["p_at_min"] = float(np.mean(ref <= ref[0] + 1e-12))
    return n, out


def law_job(args):
    kind, f, n, ia, thr, w, eps = args
    if kind == "W":
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_W, spawn_key=(1, f, n))))
        gen = law_gen(f, n)
    else:
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_L, spawn_key=(2, n, ia))))
        gen = g_witness(w=w, eps=eps)
    d = stream_dips(rng, gen, N_LAW, n)
    cnt = {a: {key: int(np.sum(d >= thr[a][key])) for key in ("q_hat", "X_l", "X_u")} for a in ALPHAS}
    keep = d if (n in NS_L or n <= 8) else None
    return (kind, f, n, ia), cnt, keep


def audit_job(args):
    """Regenerate one law stream (same seeds, same batches as law_job) and LP-audit the selected samples.
    Row: (index, D_alg, D_LP, |diff|, near, decision_mismatch) where decision_mismatch = the two implementations
    disagree on D >= X for some bracket value X (q_hat, X_l, X_u at any alpha): counts as an audit FAILURE."""
    kind, f, n, ia, w, eps, idx, thr = args
    if kind == "W":
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_W, spawn_key=(1, f, n))))
        gen = law_gen(f, n)
    else:
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_L, spawn_key=(2, n, ia))))
        gen = g_witness(w=w, eps=eps)
    b = max(200, min(20000, 2_000_000 // n))
    want = {int(i) for i in idx}
    xs = [v for a in ALPHAS for v in thr[a].values()]
    rows = []
    g0 = 0
    for s in range(0, N_LAW, b):
        k = min(b, N_LAW - s)
        x = gen(rng, (k, n))
        d = dips_of(x)
        for i in range(k):
            near = any(abs(d[i] - v) < 1e-9 for v in xs)
            if g0 + i in want or near:
                dl = float(dip_lp(np.sort(x[i])))
                mism = any((d[i] >= v) != (dl >= v) for v in xs)
                rows.append((g0 + i, float(d[i]), dl, abs(float(d[i]) - dl), bool(near), bool(mism)))
        g0 += k
    return (kind, f, n, ia), rows


def main():
    t0 = time.time()
    try:
        OUT.mkdir(parents=True, exist_ok=True)
    except OSError as err:
        sys.exit(f"cannot create {OUT}: {err}")
    pairs = len(NS_ALL) * len(ALPHAS)
    gt = GAMMA / (2 * pairs)
    idx = {a: bracket_idx(R_REF, 1 - a, gt) for a in ALPHAS}
    ref = {}
    with ProcessPoolExecutor(max_workers=11) as ex:
        for n, r in ex.map(ref_job, [(n, idx) for n in sorted(NS_ALL, key=lambda v: -v)], chunksize=1):
            ref[n] = r
    print(f"phase 1 (reference) done in {time.time() - t0:.0f}s", flush=True)
    thr = {n: {a: {"q_hat": ref[n][str(a)]["q_hat"], "X_l": ref[n][str(a)]["X_l"], "X_u": ref[n][str(a)]["X_u"]}
               for a in ALPHAS} for n in NS_ALL}
    _dump({"roots": [ROOT_REF], "R_ref": R_REF, "gamma": GAMMA, "gamma_tail": gt,
           "idx": {str(a): idx[a] for a in ALPHAS}, "by_n": {str(n): ref[n] for n in NS_ALL}}, "ref_brackets.json")
    sel = load_sel()
    jobs = [("W", f, n, -1, thr[n], 0.0, 0.0) for f in LAW_NAMES for n in NS_ALL]
    jobs += [("L", -1, n, ia, thr[n], w, e) for (n, ia), (w, e) in sel.items() if n in thr]
    jobs.sort(key=lambda j: -j[2])
    cells, keep = {}, {}
    with ProcessPoolExecutor(max_workers=11) as ex:
        for key, cnt, d in ex.map(law_job, jobs, chunksize=4):
            cells[key] = cnt
            if d is not None:
                keep[key] = d
    print(f"phase 2 (laws) done in {time.time() - t0:.0f}s", flush=True)
    for (kind, f, n, ia), d in keep.items():
        np.save(OUT / f"dips_{kind}_f{f}_n{n}_ia{ia}.npy", d)
    KW = len(LAW_NAMES) * len(NS_ALL) * len(ALPHAS)
    KL = len(sel) * len(ALPHAS)
    table = []
    for (kind, f, n, ia), cnt in sorted(cells.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], kv[0][3])):
        K = KW if kind == "W" else KL
        dl = (DELTA_W if kind == "W" else DELTA_L) / (2 * K)
        for a in ALPHAS:
            c = cnt[a]
            row = {"part": kind, "f": f, "law": LAW_NAMES[f] if kind == "W" else f"W{sel[(n, ia)]}",
                   "n": n, "ia": ia, "alpha": a, "N": N_LAW, "k_q": c["q_hat"], "k_Xl": c["X_l"],
                   "k_Xu": c["X_u"], "rej_point": c["q_hat"] / N_LAW,
                   "lower": cp_lower(c["X_u"], N_LAW, 1 - dl), "upper": cp_upper(c["X_l"], N_LAW, 1 - dl)}
            row["class"] = "ABOVE" if row["lower"] > a else ("BELOW" if row["upper"] < a else "UNRESOLVED")
            table.append(row)
    # anchor (P-3): the uniform enclosure must contain alpha in every cell
    anchor_fail = [r for r in table if r["part"] == "W" and r["f"] == 0 and not (r["lower"] <= r["alpha"] <= r["upper"])]
    # audit
    ajobs = []
    m_audit = 10 if QUICK else 100
    for f in AUDIT_F:
        for n in AUDIT_NS:
            rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_AUDIT, spawn_key=(3, f, n))))
            ajobs.append(("W", f, n, -1, 0.0, 0.0, rng.choice(N_LAW, size=m_audit, replace=False), thr[n]))
    for (n, ia), (w, e) in sel.items():
        if n in AUDIT_NS:
            rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT_AUDIT, spawn_key=(4, n, ia))))
            ajobs.append(("L", -1, n, ia, w, e, rng.choice(N_LAW, size=m_audit, replace=False), thr[n]))
    audit = {}
    with ProcessPoolExecutor(max_workers=11) as ex:
        for (kind, f, n, ia), rows in ex.map(audit_job, ajobs, chunksize=1):
            audit[f"{kind}_f{f}_n{n}_ia{ia}"] = {
                "audited": len(rows), "near": sum(r[4] for r in rows),
                "max_abs_diff": max((r[3] for r in rows), default=0.0),
                "decision_mismatches": sum(r[5] for r in rows),
                "pass": all(r[3] <= 1e-9 and not r[5] for r in rows)}
    summary = {"tier": "DEV-QUICK" if QUICK else "CONFIRM", "roots": [ROOT_REF, ROOT_W, ROOT_L, ROOT_AUDIT],
               "R_ref": R_REF, "N_law": N_LAW, "deltas": [GAMMA, DELTA_W, DELTA_L], "KW": KW, "KL": KL,
               "sel": {f"{n}_{ia}": v for (n, ia), v in sel.items()}, "table": table,
               "anchor_fail": anchor_fail, "anchor_ok": not anchor_fail,
               "audit": audit, "audit_ok": all(v["pass"] for v in audit.values()), "seconds": time.time() - t0}
    _dump(summary, "confirm_WL.json")
    print(f"anchor ok: {summary['anchor_ok']} ({len(anchor_fail)} failures) | audit ok: {summary['audit_ok']} | "
          f"{time.time() - t0:.0f}s", flush=True)


def _dump(obj, name):
    try:
        (OUT / name).write_text(json.dumps(obj, indent=1))
    except OSError as err:
        sys.exit(f"cannot write {name}: {err}")


if __name__ == "__main__":
    main()
