"""Step 16 (DT15) A2: paired two-bridge Monte Carlo certificate of C-edge at the standard cutoffs, with cell envelopes
and adaptive refinement.

Target per level (alpha, t'): for every u = 1 - eps^2 with eps in [EPS_LO(level), 1/sqrt 2]:
    P(d_u >= t') - P(d_1 >= t') < 0,    d_u = d(H_u, K_1/2) (DT10 Lemma S), d_1 = G(b1)/2.
Per path and cell C = [eps_lo, eps_hi] (evaluation point u_c = 1 - eps_lo^2, a_c = sqrt(u_c)):
    Y_C := hi_C - lo_1 >= sup_{u in C} 1{d_u >= t'} - 1{d_1 >= t'},
    hi_C = 1{ the PL split path at u_c is NOT certainly feasible at e = t' - E_c - w_C },
    lo_1 = 1{ G(b1_PL)/2 - E1 >= t' + PAD },
    E_c = max(a_c E1, eps_lo E2) (1-Lipschitz, Lemma S), w_C = max((a_c - a_hi)(|b1_PL|+E1), (eps_hi - eps_lo)(|b2_PL|+E2))
    (DT10 Corollary S': disjoint supports, fixed seam; |b| <= |b_PL| + E).
Feasibility (DT10 Lemma F(d)): d <= e iff 2e >= max(a G1, eps G2) and a Lam1(e/a) + eps Lam2(e/eps) <= 2e.
MONO (v2, replaces the Lipschitz cell modulus w): as eps grows in a cell, a G1 and a Lam1(e/a) decrease (Lam
non-increasing, DT10 Lam0(ii)), eps G2 and eps Lam2(e/eps) increase. So feasibility of the PL path at the MIXED corner
(a_lo = sqrt(1 - eps_lo^2) on side 1, eps_hi on side 2) at e = t - max(a_lo E1, eps_hi E2) implies d_u < t for every u in
the cell (1-Lipschitz residual: d_u <= d_PL(u) + max(a E1, eps E2) <= d_PL(u) + max(a_lo E1, eps_hi E2)).
"Certainly feasible" uses: (quick) Lemma G + Lambda <= G: 2e >= a (G1 - Q1) + eps G2 with Q1 = G1 - Lam1(G1/2);
(full) exact-predicate Lambda (int64, unit 2^-32) evaluated at y' < y (Lambda non-increasing => upper bound), PAD margins.
Refinement: decided from SKELETON values only with the deterministic window W(m) = (1 - 1e-6)-quantile of the
exact law of E at depth m; residual sups E1, E2 drawn ONCE at the final depth from P(E <= x) = K(x sqrt m)^m.
Refinement doubles m by conditional bridge midpoints (N(0, h/4) given the end values).
Usage: python dt16_pair_cert.py M0 MMAX R shard nshards root tag  -> deposits/dt16/pair/<tag>_r<root>_s<shard>.json
"""
import json
import math
import os
import sys
import time

import numpy as np
from numba import njit
from scipy.stats import kstwobign

ROOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SCALE = 2.0 ** 32
BMAX = 2 ** 35
PAD = 1e-9
ROUND = 2.0 ** -32
E_HI = 0.70710679                        # > 1/sqrt(2) (outward)
LEVELS = ((0.10, 0.4935, 0.125), (0.05, 0.5378, 0.130), (0.01, 0.6161, 0.150))   # (alpha, t' = q_- - 1e-4, EPS_LO)


def cell_edges(eps_lo):
    """Cells with width max(3e-4, min(0.03, 0.08 eps^2)) from eps_lo to E_HI (frozen rule)."""
    e = [eps_lo]
    while e[-1] < E_HI:
        w = max(3e-4, min(0.03, 0.08 * e[-1] ** 2))
        e.append(min(E_HI, e[-1] + w))
    return e


@njit(cache=True)
def gap_int(B):
    n = len(B)
    hx = np.empty(n, np.int64)
    hy = np.empty(n, np.int64)
    k = 0
    for i in range(n):
        while k >= 2 and (hy[k - 1] - hy[k - 2]) * (i - hx[k - 2]) >= (B[i] - hy[k - 2]) * (hx[k - 1] - hx[k - 2]):
            k -= 1
        hx[k] = i
        hy[k] = B[i]
        k += 1
    g = 0.0
    j = 0
    for i in range(n):
        while j < k - 2 and hx[j + 1] <= i:
            j += 1
        dx = hx[j + 1] - hx[j]
        num = (B[i] - hy[j]) * dx - (hy[j + 1] - hy[j]) * (i - hx[j])
        v = float(num) / float(dx) / 4294967296.0
        if v > g:
            g = v
    return g


@njit(cache=True)
def lam_int(B, Y2):
    """Lambda_b(y') with 2y' = Y2 / 2^32: max(0, sup_{s<x<1}[b(x) + (1-x)(b(x)-b(s)-2y')/(x-s)]) (DT10 Lemma Lambda0(i))."""
    n = len(B)
    m = n - 1
    hx = np.empty(n, np.int64)
    hy = np.empty(n, np.int64)
    k = 0
    best = 0.0
    for j in range(n - 1):
        if k >= 1 and j >= 1:
            qy = B[j] - Y2
            lo = 0
            hi = k - 1
            while lo < hi:
                mid = (lo + hi) // 2
                if (qy - hy[mid]) * (j - hx[mid + 1]) >= (qy - hy[mid + 1]) * (j - hx[mid]):
                    hi = mid
                else:
                    lo = mid + 1
            N = qy - hy[lo]
            D = j - hx[lo]
            v = float(B[j]) / 4294967296.0 + float(m - j) * float(N) / (float(D) * 4294967296.0)
            if v > best:
                best = v
        while k >= 2 and (hy[k - 1] - hy[k - 2]) * (j - hx[k - 2]) >= (B[j] - hy[k - 2]) * (hx[k - 1] - hx[k - 2]):
            k -= 1
        hx[k] = j
        hy[k] = B[j]
        k += 1
    return best


def lam_upper(B, y):
    """Upper bound of Lambda_{b_PL}(y) (Lambda non-increasing in y: evaluate at y' slightly below y)."""
    Y2 = int(math.floor(2.0 * y * SCALE)) - 2
    if Y2 <= 0 or Y2 > 2 ** 36:
        return math.inf
    return lam_int(B, Y2) + PAD


def lam_lower(B, y):
    """Lower bound of Lambda_{b_PL}(y) (evaluate at y' slightly above y)."""
    Y2 = int(math.ceil(2.0 * y * SCALE)) + 2
    if Y2 > 2 ** 36:
        return 0.0
    return max(0.0, lam_int(B, Y2) - PAD)


def bridge(rng, m):
    z = rng.standard_normal(m)
    w = np.concatenate([[0.0], np.cumsum(z)]) / math.sqrt(m)
    return w - (np.arange(m + 1) / m) * w[-1]


def refine(rng, b):
    """Double the grid: conditional Brownian-bridge midpoints, N(mean of ends, h/4) with h = 1/m."""
    m = len(b) - 1
    mid = 0.5 * (b[:-1] + b[1:]) + rng.standard_normal(m) * math.sqrt(0.25 / m)
    out = np.empty(2 * m + 1)
    out[0::2] = b
    out[1::2] = mid
    return out


EPS_REFINE = 0.2                       # refinement is decided only from cells with eps_lo < EPS_REFINE (frozen rule)


def window(m):
    """(1 - 1e-4)-quantile of E_m: P(E_m <= x) = K(x sqrt m)^m (deterministic; used only to DECIDE refinement)."""
    return float(kstwobign.ppf((1.0 - 1e-4) ** (1.0 / m))) / math.sqrt(m)


def draw_E(rng, m):
    v = rng.uniform()
    return float(kstwobign.isf(-np.expm1(np.log(v) / m))) / math.sqrt(m)


class Path:
    def __init__(self, b1, b2):
        self.m = len(b1) - 1           # depth of the long side b1 (the only side that is refined)
        self.m2 = len(b2) - 1
        self.B1 = np.rint(b1 * SCALE).astype(np.int64)
        self.B2 = np.rint(b2 * SCALE).astype(np.int64)
        self.ok = bool(np.max(np.abs(self.B1)) <= BMAX and np.max(np.abs(self.B2)) <= BMAX)
        self.G1 = gap_int(self.B1) if self.ok else math.inf
        self.G2 = gap_int(self.B2) if self.ok else math.inf
        self.n1 = float(np.max(np.abs(b1)))
        self.n2 = float(np.max(np.abs(b2)))
        self.Q1 = None

    def q1(self):
        if self.Q1 is None:
            lam = lam_upper(self.B1, self.G1 / 2.0 + 2 * PAD)
            self.Q1 = max(0.0, self.G1 - lam - 4 * PAD)        # Q1 = S(G1/2) = G1 - Lam1(G1/2), lower bound
        return self.Q1

    def certainly_feasible(self, e, a, eps):
        """True only if the PL split path at (a, eps) has distance <= e (with margins)."""
        if not self.ok or e <= 0:
            return False
        if 2 * e < max(a * self.G1, eps * self.G2) + 4 * PAD:
            return False
        if 2 * e >= a * (self.G1 - self.q1()) + eps * self.G2 + 8 * PAD:          # Lemma G + Lambda <= G
            return True
        l1 = lam_upper(self.B1, e / a)
        l2 = lam_upper(self.B2, e / eps)
        return a * l1 + eps * l2 <= 2 * e - 4 * PAD

    def certainly_infeasible(self, e, a, eps):
        """True only if distance > e (used ONLY for skeleton-only refinement decisions, not for Y)."""
        if not self.ok:
            return True
        if 2 * e < max(a * self.G1, eps * self.G2) - 4 * PAD:
            return True
        if e <= 0:
            return True
        return a * lam_lower(self.B1, e / a) + eps * lam_lower(self.B2, e / eps) > 2 * e + 4 * PAD


def path_Y(P, E1, E2, cells_by_level):
    out = []
    for (_alpha, t, _), edges in zip(LEVELS, cells_by_level, strict=True):
        lo1 = 1 if P.ok and t + PAD <= P.G1 / 2.0 - E1 else 0
        row = []
        for i in range(len(edges) - 1):
            el, eh = edges[i], edges[i + 1]
            a = math.sqrt(1 - el * el)
            Ec = max(a * E1, eh * E2)
            hi = 0 if P.certainly_feasible(t - Ec, a, eh) else 1          # MIXED corner (a_lo, eps_hi), see MONO
            row.append(hi - lo1)
        out.append(row)
    return out


def ambiguous(P, cells_by_level):
    """Skeleton-only refinement rule (no residual draw is inspected): window W = window(m)."""
    W, W2 = window(P.m), window(P.m2)
    for (_alpha, t, _), edges in zip(LEVELS, cells_by_level, strict=True):
        if t - W <= P.G1 / 2.0 < t + W:
            return True
        for i in range(len(edges) - 1):
            el, eh = edges[i], edges[i + 1]
            if el >= EPS_REFINE:
                break
            a = math.sqrt(1 - el * el)
            Wc = max(a * W, eh * W2)
            if P.certainly_feasible(t - Wc, a, eh):
                continue                       # hi = 0 already with the window: refinement cannot help this cell
            if not P.certainly_infeasible(t, a, el):
                return True        # PL distance may lie in (t - Wc - w, t]: refinement can turn hi = 1 into 0
    return False


def main():
    m0, mmax, R, shard, ns, root = (int(x) for x in sys.argv[1:7])
    tag = sys.argv[7]
    cells = [cell_edges(lv[2]) for lv in LEVELS]
    rng_paths = np.random.Generator(np.random.PCG64(np.random.SeedSequence([root, m0, shard, 0])))
    rng_E = np.random.Generator(np.random.PCG64(np.random.SeedSequence([root, m0, shard, 1])))
    reps = R // ns + (1 if shard < R % ns else 0)
    S1 = [np.zeros(len(c) - 1) for c in cells]
    S2 = [np.zeros(len(c) - 1) for c in cells]
    NP = [np.zeros(len(c) - 1) for c in cells]
    NM = [np.zeros(len(c) - 1) for c in cells]
    depth_hist, bad, t0 = {}, 0, time.time()
    for _ in range(reps):
        b1, b2 = bridge(rng_paths, m0), bridge(rng_paths, m0)
        P = Path(b1, b2)
        while P.m < mmax and ambiguous(P, cells):
            b1 = refine(rng_paths, b1)
            P = Path(b1, b2)
        depth_hist[P.m] = depth_hist.get(P.m, 0) + 1
        bad += 0 if P.ok else 1
        E1, E2 = draw_E(rng_E, P.m), draw_E(rng_E, P.m2)       # drawn once each, at the final depths
        E1, E2 = E1 * (1 + 1e-9) + ROUND, E2 * (1 + 1e-9) + ROUND   # isf inversion slack; rounding |b - B/2^32| <= 2^-33
        Y = path_Y(P, E1, E2, cells)
        for li in range(len(LEVELS)):
            y = np.array(Y[li], float)
            S1[li] += y
            S2[li] += y * y
            NP[li] += y > 0
            NM[li] += y < 0
    res = {"m0": m0, "mmax": mmax, "R": R, "shard": shard, "nshards": ns, "root": root, "tag": tag, "reps": reps,
           "levels": LEVELS, "edges": cells, "S1": [s.tolist() for s in S1], "S2": [s.tolist() for s in S2],
           "n_plus": [s.tolist() for s in NP], "n_minus": [s.tolist() for s in NM], "depths": depth_hist,
           "unresolved_paths": bad, "sec": time.time() - t0}
    out = os.path.join(ROOT_DIR, "deposits", "dt16", "pair", f"{tag}_r{root}_s{shard}.json")
    try:
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out + ".tmp", "w") as fh:
            json.dump(res, fh)
        os.replace(out + ".tmp", out)
    except OSError as err:
        sys.exit(f"cannot write {out}: {err}")
    print(f"{tag} shard {shard}: {reps} paths {res['sec']:.0f}s depths {depth_hist} unresolved {bad}", flush=True)


if __name__ == "__main__":
    main()
