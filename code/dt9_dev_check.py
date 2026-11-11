"""DT9 DEV checks (root 2026110720; not a CONFIRM run).
(1) decision B&B == brute force over all pieces (same certified piece_ub), at cutoffs at the true value +- tiny;
(2) size of the repair term (certified bound minus band error) vs the DT8 heuristic;
(3) speed of the decision driver at n = 800, 1600 (per sample, one core);
(4) one-sidedness under injected LP failures (a true exceedance must never be reported as none).
Usage: python dt9_dev_check.py [quick]
"""
import json
import os
import random
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dt9_block_cert as D  # noqa: E402
from sup_pinned import piece_dist, pieces  # noqa: E402

ROOT = 2026110720


def brute(x, lo_u, hi_u):
    return max(D.piece_ub(x, a, b) for a, b in pieces(x, lo_u, hi_u))


def repair_exact(rng, trials):
    """(5) Perturb LP solutions, then check in exact rationals that g + h is in C_[a,b] and that its band
    error is <= the certified bound returned by bound_of."""
    bad_shape, bad_bound, nonzero = 0, 0, 0
    for _ in range(trials):
        n = int(rng.integers(8, 40))
        x = np.sort(rng.uniform(size=n))
        P = pieces(x, 0.5, 1.0)
        a, b = P[int(rng.integers(len(P)))]
        ts, lo, hi, ia, ib, g = D.solve(x, a, b)
        g = g + rng.normal(scale=10.0 ** rng.uniform(-9, -4), size=len(g))
        ub = D.bound_of(ts, lo, hi, ia, ib, g)
        T = [Fraction(float(v)) for v in ts]
        Gq = [Fraction(float(v)) for v in g]
        K = len(T)
        sl = [(Gq[k + 1] - Gq[k]) / (T[k + 1] - T[k]) for k in range(K - 1)]
        J = {k: sl[k] - sl[k - 1] for k in range(1, K - 1)}
        conv = [k for k in range(1, K - 1) if k + 1 <= ib]
        conc = [k for k in range(1, K - 1) if k - 1 >= ia]
        dl = {k: max(-J[k], Fraction(0)) for k in conv}
        dr = {k: max(J[k], Fraction(0)) for k in conc}
        nonzero += int(any(dl.values()) or any(dr.values()))
        h = [sum((dl[k] * max(T[k] - T[i], Fraction(0)) for k in conv), Fraction(0))
             - sum((dr[k] * max(T[i] - T[k], Fraction(0)) for k in conc), Fraction(0)) for i in range(K)]
        g2 = [Gq[i] + h[i] for i in range(K)]
        s2 = [(g2[k + 1] - g2[k]) / (T[k + 1] - T[k]) for k in range(K - 1)]
        for k in conv:
            bad_shape += int(s2[k] - s2[k - 1] < 0)
        for k in conc:
            bad_shape += int(s2[k] - s2[k - 1] > 0)
        err = max(max(Fraction(float(hi[i])) - g2[i], g2[i] - Fraction(float(lo[i]))) for i in range(K))
        bad_bound += int(err > Fraction(ub))
    return {"trials": trials, "with_violations": nonzero, "bad_shape": bad_shape, "bad_bound": bad_bound}


def adversarial(rng):
    """(6) failure rates 0.5 and 1.0; equality t = v; near-coincident points; point on a block boundary."""
    real = D.piece_ub
    res = {}
    for rate in (0.5, 1.0):
        random.seed(ROOT + int(rate * 10))

        def flaky(x, a, b, rate=rate):
            if random.random() < rate:
                raise RuntimeError("injected")
            return real(x, a, b)

        bad = 0
        for _ in range(15):
            x = np.sort(rng.uniform(size=60))
            for lo_u, hi_u in D.BLOCKS:
                P = pieces(x, lo_u, hi_u)
                v = brute(x, lo_u, hi_u)
                for t in (v * (1 - 1e-9), 0.7315 / np.sqrt(60)):
                    D.piece_ub = flaky
                    e, _ = D.exceeds(x, P, t)
                    D.piece_ub = real
                    bad += int((v > t) and not e)
        res[f"fail_rate_{rate}"] = bad
    eq_bad, odd = 0, 0
    for kind in range(20):
        x = np.sort(rng.uniform(size=50))
        if kind % 3 == 1:
            x[10] = x[9] + 1e-13            # near-coincident points
        elif kind % 3 == 2:
            x[30] = 0.9                     # a sample point exactly on a block boundary
        x = np.sort(x)
        for lo_u, hi_u in D.BLOCKS:
            P = pieces(x, lo_u, hi_u)
            v = brute(x, lo_u, hi_u)
            e, _ = D.exceeds(x, P, v)       # equality: the brute bound itself is not "> v"
            odd += int(e)
            e2, _ = D.exceeds(x, P, v * (1 - 1e-9))
            eq_bad += int(not e2)        # safe by construction (a node bound <= t exists); information only
    res["equality_reported_exceed"] = odd
    res["decision_no_brute_yes_at_v(1-1e-9)"] = eq_bad
    return res


def main():
    quick = len(sys.argv) > 1
    out = {}
    rng5 = np.random.default_rng(np.random.SeedSequence([ROOT, 5]))
    out["repair_exact"] = repair_exact(rng5, 20 if quick else 80)
    out["adversarial"] = adversarial(rng5)
    rng = np.random.default_rng(np.random.SeedSequence([ROOT, 1]))
    mism, checks, reps, gaps, detail = 0, 0, [], [], []
    for n, R in ((50, 20 if quick else 120), (100, 10 if quick else 60), (200, 4 if quick else 25)):
        for _ in range(R):
            x = np.sort(rng.uniform(size=n))
            for lo_u, hi_u in D.BLOCKS:
                P = pieces(x, lo_u, hi_u)
                v = brute(x, lo_u, hi_u)
                for t in (v * (1 - 1e-9), v * (1 + 1e-9), 0.7315 / np.sqrt(n), 0.86 / np.sqrt(n), v - 1e-6):
                    e, _ = D.exceeds(x, P, t)
                    checks += 1
                    if e != (v > t):
                        mism += 1
                        # parent bound below t while a child's own (looser) bound is above t?
                        root_ub = D.piece_ub(x, P[0][0], P[-1][1])
                        detail.append({"n": n, "block": [lo_u, hi_u], "decision": bool(e), "brute_v": v,
                                       "t": t, "v_minus_t": v - t, "root_ub": root_ub, "npieces": len(P)})
                for a, b in P[:: max(1, len(P) // 3)]:
                    u1 = D.piece_ub(x, a, b)
                    u0 = piece_dist(x, a, b)
                    gaps.append(u1 - u0)
            reps.append(D.STATS["max_repair"])
    out["decision_vs_brute"] = {"checks": checks, "mismatches": mism, "detail": detail}
    out["repair_max"] = max(reps)
    out["ub_minus_dt8_value"] = [float(min(gaps)), float(max(gaps))]
    # (4) injected failures
    real = D.piece_ub
    random.seed(ROOT)

    def flaky(x, a, b):
        if random.random() < 0.1:
            raise RuntimeError("injected")
        return real(x, a, b)

    bad, tot = 0, 0
    rng2 = np.random.default_rng(np.random.SeedSequence([ROOT, 2]))
    for _ in range(10 if quick else 60):
        x = np.sort(rng2.uniform(size=100))
        for lo_u, hi_u in D.BLOCKS:
            P = pieces(x, lo_u, hi_u)
            v = brute(x, lo_u, hi_u)
            for t in (v * (1 - 1e-9), v * (1 + 1e-9), 0.7315 / 10):
                D.piece_ub = flaky
                e, _ = D.exceeds(x, P, t)
                D.piece_ub = real
                tot += 1
                bad += int((v > t) and not e)
    out["injected_failures"] = {"checks": tot, "true_exceed_reported_none": bad}
    # (3) speed
    sp = {}
    for n, R in ((800, 6 if quick else 40), (1600, 3 if quick else 20)):
        rng3 = np.random.default_rng(np.random.SeedSequence([ROOT, 3, n]))
        D.STATS["lps"] = 0
        t0 = time.time()
        ex = []
        for _ in range(R):
            r, _ = D.sample_row(np.sort(rng3.uniform(size=n)))
            ex.append([b[0] for b in r])
        sp[n] = {"sec_per_sample": (time.time() - t0) / R, "lps_per_sample": D.STATS["lps"] / R,
                 "edge_exceed_rate": float(np.mean([e[3] for e in ex]))}
    out["speed"] = sp
    out["max_repair_overall"] = D.STATS["max_repair"]
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits", "dt9",
                        "dev_check_quick.json" if quick else "dev_check.json")
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            json.dump(out, fh, indent=1)
    except OSError as err:
        sys.exit(f"cannot write {path}: {err}")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
