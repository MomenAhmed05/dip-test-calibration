"""One-shoulder excursion test, Hartigan (2000)-inspired — S's implementation for the Step 7 benchmark.

Source: J. A. Hartigan (2000), Testing for Antimodes, in Gaul, Opitz & Schader (eds.), Data Analysis, Springer,
pp. 169-181 (our copy: prior/priorart/step1/user_book_hartigan2000/hartigan2000_ch14.txt, sections 1, 4, 6).

Statistic (Hartigan §4, frozen 2026-11-06 BEFORE any Step 7 size/power number): for a shoulder segment with
m interior points, T = N*D_s/sqrt(m), D_s = the maximum excursion of the ECDF above (left shoulder) or below (right
shoulder) the chord of the segment. Reference: Kennedy's (1976) maximum-Brownian-excursion law = Kuiper law,
G(x) = 1 - 2*sum_j (4 j^2 x^2 - 1) exp(-2 j^2 x^2); exact quantiles 1.74725995 (95%), 2.00091812 (99%) (checked by S).
Normalisation choice: §4's raw excursion N*D_s/sqrt(m) is used because §4 states its reference law for exactly this
statistic; §6's scaled dip N*DIP/sqrt(m) is recorded through R = D_s / dip for validation only.

SIMPLIFICATION (label every use): Hartigan builds a hierarchical tree of shoulders from m-modal fits.
We take ONE shoulder: the hull segment of the dip's own unimodal fit (GCM of the ECDF left of the modal interval,
LCM right of it) that contains the point of maximal deviation. Name in all reports: "one-shoulder excursion test,
Hartigan (2000)-inspired". Its finite reference G_{m+1} (§4, eq. (9)) is NOT implemented: we use the asymptotic
Kennedy law, and separately an own-n uniform calibration of the whole procedure (two calibrations, labelled).
Eligibility: a selected shoulder with m < 5 interior points => no rejection (T = 0), counted.
Ties: continuous samples only; a tie or an out-of-range dip index raises an error (no fallback).

Self-test: python code/hartigan2000.py --test
"""

import sys

import diptest
import numpy as np

KENNEDY_Q = {0.05: 1.74725995, 0.01: 2.00091812}
M_MIN = 5


def _lower_hull(px, py):
    """Indices of the lower convex hull of points sorted by x (monotone chain)."""
    h = []
    for i in range(len(px)):
        while len(h) >= 2:
            a, b = h[-2], h[-1]
            if (py[b] - py[a]) * (px[i] - px[a]) >= (py[i] - py[a]) * (px[b] - px[a]):
                h.pop()
            else:
                break
        h.append(i)
    return h


def shoulder_full(x):
    """For a sorted continuous sample x return (T, m, side, ratio R = D_s/dip, eligible)."""
    n = len(x)
    if np.any(np.diff(x) <= 0):
        raise ValueError("ties or unsorted input: the benchmark is defined for continuous samples only")
    res = diptest.dipstat(x, full_output=True, allow_zero=False, sort_x=False)
    if not isinstance(res, tuple):
        raise TypeError("diptest.dipstat(full_output=True) did not return a tuple")
    dip, o = res
    xl, xu, di = o["xl"], o["xu"], int(o["dipidx"])
    if not 0 <= di < n:
        raise ValueError(f"dip index {di} out of range (n={n})")
    xi = x[di]
    if xl < xi < xu:
        # explicit policy: a maximal-deviation point inside the modal interval has no shoulder; recorded failure
        raise ValueError("dip index inside the modal interval")
    j = np.arange(1, n + 1)
    if xi <= xl:
        idx = np.nonzero(x <= xl)[0]
        px = x[idx]
        hull = _lower_hull(px, (j[idx] - 1) / n)  # GCM of the left limits (lower jaw)
        side = "left"
    else:
        idx = np.nonzero(x >= xu)[0]
        px = x[idx]
        hull = _lower_hull(px, -(j[idx] / n))  # LCM of the right limits = -(lower hull of the negated values)
        side = "right"
    if len(hull) < 2:
        return 0.0, 0, side, 0.0, False
    hx = px[hull]
    k = int(np.clip(np.searchsorted(hx, xi, side="right") - 1, 0, len(hull) - 2))
    a, b = hull[k], hull[k + 1]
    m = b - a - 1  # interior points of the shoulder segment
    if m < 1:
        return 0.0, max(m, 0), side, 0.0, False
    seg = slice(a, b + 1)
    xs = px[seg]
    if side == "left":
        ya, yb = (j[idx][a] - 1) / n, (j[idx][b] - 1) / n
        line = ya + (xs - px[a]) * (yb - ya) / (px[b] - px[a])
        ds = float(np.max(j[idx][seg] / n - line))
    else:
        ya, yb = j[idx][a] / n, j[idx][b] / n
        line = ya + (xs - px[a]) * (yb - ya) / (px[b] - px[a])
        ds = float(np.max(line - (j[idx][seg] - 1) / n))
    eligible = m >= M_MIN
    t = n * ds / np.sqrt(m) if eligible else 0.0
    return t, m, side, (ds / dip if dip > 0 else np.nan), eligible


def shoulder_stat(x):
    """Back-compatible: (T, m, side)."""
    t, m, side, _, _ = shoulder_full(x)
    return t, m, side


def stats_of(xs):
    """Row-wise T for a 2-D array of sorted samples (ineligible shoulders give T = 0 = no rejection)."""
    out = np.empty(xs.shape[0])
    for i in range(xs.shape[0]):
        out[i] = shoulder_full(xs[i])[0]
    return out


def stats_full_of(xs):
    """Row-wise (T, m, R, eligible) for validation readouts."""
    out = np.empty((xs.shape[0], 4))
    for i in range(xs.shape[0]):
        t, m, _, r, e = shoulder_full(xs[i])
        out[i] = (t, m, r, float(e))
    return out


def _self_test():
    rng = np.random.default_rng(12345)  # scratch seed for the self-test only (not a study root)
    bad = 0
    for _ in range(200):
        x = np.sort(rng.random(300))
        t1, m1, s1, _, _ = shoulder_full(x)
        t2, m2, s2, _, _ = shoulder_full(np.sort(1 - x))  # reflection
        t3, m3, s3, _, _ = shoulder_full(5 + 3 * x)  # affine
        if not (np.isclose(t1, t3) and m1 == m3 and s1 == s3):
            bad += 1
        if not (np.isclose(t1, t2, rtol=1e-9, atol=1e-9) and m1 == m2):
            bad += 1000  # reflection mismatch counted separately
    try:
        shoulder_full(np.array([0.1, 0.2, 0.2, 0.5]))
        raise AssertionError("ties not rejected")
    except ValueError:
        pass
    print(f"self-test: affine mismatches {bad % 1000}, reflection mismatches {bad // 1000} (of 200); ties rejected OK")
    if bad:
        raise AssertionError("hartigan2000 self-test FAILED")
    return bad


if __name__ == "__main__":
    if "--test" in sys.argv:
        _self_test()
