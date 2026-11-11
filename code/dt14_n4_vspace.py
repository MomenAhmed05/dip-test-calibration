"""Step 15 / DT14 WP-B2: exact n = 4 rejection probability for a STEP quantile density (v-space engine).

Object: law F with quantile function Q = F^{-1}, quantile density q = Q' piecewise constant:
q = c_j on cell [b_j, b_{j+1}), 0 = b_0 < ... < b_K = 1, c_j > 0. F is unimodal iff (c_j) is
valley-shaped (nonincreasing then nondecreasing). With V_(1) < ... < V_(4) uniform order statistics and
I_k = int_{V_(k)}^{V_(k+1)} q (k = 1, 2, 3):
    P_F(T >= r) = P(I_2 >= r I_1, I_2 >= r I_3)   (H&H n = 4 reduction, B0).
For each assignment of the 4 ordered points to cells (multiset j_1 <= j_2 <= j_3 <= j_4) the I_k are
linear in v, so the event is a polytope in v-space: exact volume (exact_poly, Fractions) times 4!.
Independent of the Khintchine engine (dt14_n4_kernel.py); cross-checked in _validate().
"""
import sys
from fractions import Fraction
from itertools import combinations_with_replacement
from math import factorial

import numpy as np

sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
from dt14_n4_kernel import _vol_float  # noqa: E402
from exact_poly import volume as exact_volume  # noqa: E402


def _seg_row(i, j, ki, kj, b, c):
    """Coefficients (row over v_1..v_4, const) of int_{v_ki}^{v_kj} q when v_ki in cell i, v_kj in cell j>=i."""
    row = [Fraction(0)] * 4
    if i == j:
        row[kj] += c[i]
        row[ki] -= c[i]
        return row, Fraction(0)
    const = Fraction(0)
    row[ki] -= c[i]
    const += c[i] * b[i + 1]
    for m in range(i + 1, j):
        const += c[m] * (b[m + 1] - b[m])
    row[kj] += c[j]
    const -= c[j] * b[j]
    return row, const


def prob_exact(b, c, r):
    b = [Fraction(x) for x in b]
    c = [Fraction(x) for x in c]
    r = Fraction(r)
    K = len(c)
    tot = Fraction(0)
    for js in combinations_with_replacement(range(K), 4):
        A, rhs = [], []
        for k in range(4):  # cell boxes: b_j <= v_k <= b_{j+1}
            e = [Fraction(0)] * 4
            e[k] = Fraction(-1)
            A.append(e)
            rhs.append(-b[js[k]])
            e = [Fraction(0)] * 4
            e[k] = Fraction(1)
            A.append(e)
            rhs.append(b[js[k] + 1])
        for k in range(3):  # order v_k <= v_{k+1}
            e = [Fraction(0)] * 4
            e[k] = Fraction(1)
            e[k + 1] = Fraction(-1)
            A.append(e)
            rhs.append(Fraction(0))
        segs = [_seg_row(js[k], js[k + 1], k, k + 1, b, c) for k in range(3)]
        (r1, c1), (r2, c2), (r3, c3) = segs
        # I2 - r I1 >= 0  <=>  r*I1 - I2 <= 0 :  (r*r1 - r2) . v <= c2 - r*c1
        A.append([r * x - y for x, y in zip(r1, r2, strict=True)])
        rhs.append(c2 - r * c1)
        A.append([r * x - y for x, y in zip(r3, r2, strict=True)])
        rhs.append(c2 - r * c3)
        tot += exact_volume(A, rhs)
    return tot * factorial(4)


def prob_float(b, c, r):
    """Float twin of prob_exact (DEV searches only). Each event row is normalised by its max |coef|."""
    bf = [Fraction(x) for x in b]
    cf = [Fraction(x) for x in c]
    K = len(cf)
    tot = 0.0
    for js in combinations_with_replacement(range(K), 4):
        A, rhs = [], []
        for k in range(4):
            e = [0.0] * 4
            e[k] = -1.0
            A.append(e)
            rhs.append(-float(bf[js[k]]))
            e = [0.0] * 4
            e[k] = 1.0
            A.append(e)
            rhs.append(float(bf[js[k] + 1]))
        for k in range(3):
            e = [0.0] * 4
            e[k] = 1.0
            e[k + 1] = -1.0
            A.append(e)
            rhs.append(0.0)
        segs = [_seg_row(js[k], js[k + 1], k, k + 1, bf, cf) for k in range(3)]
        (r1, c1), (r2, c2), (r3, c3) = segs
        for ra, ca in ((r1, c1), (r3, c3)):
            row = np.array([float(r * x - y) for x, y in zip(ra, r2, strict=True)])
            cst = float(c2 - r * ca)
            s = max(np.abs(row).max(), abs(cst), 1e-300)
            A.append(list(row / s))
            rhs.append(cst / s)
        tot += _vol_float(A, rhs)
    return tot * 24


def _validate():
    ok = True
    for r in (Fraction(3), Fraction(5)):
        got = prob_exact([0, Fraction(1, 3), Fraction(7, 10), 1], [2, 2, 2], r)
        want = Fraction(2) / ((r + 1) * (r + 2))
        good = got == want
        ok &= good
        print(f"{'PASS' if good else 'FAIL'} uniform (3 equal cells) r={r}: {got} vs {want}")
    # cross-check vs Khintchine engine: one-sided two-step, mode at x = 0 (left), q nondecreasing
    from dt14_n4_kernel import phi_exact
    for (tau, c0, c1, r) in ((Fraction(2, 5), Fraction(1, 3), Fraction(2), Fraction(3)),
                            (Fraction(1, 2), Fraction(1, 10), Fraction(1), Fraction(5)),
                            (Fraction(3, 10), Fraction(1), Fraction(7, 2), Fraction(3))):
        got = prob_exact([0, tau, 1], [c0, c1], r)
        # x-space: density 1/c0 on [0, c0*tau), 1/c1 on [c0*tau, c0*tau + c1*(1-tau)]
        x1 = c0 * tau
        x2 = x1 + c1 * (1 - tau)
        w_out = x2 / c1
        w_in = (1 / c0 - 1 / c1) * x1
        want = phi_exact([w_out, w_in], [x2, x1], r)
        good = got == want
        ok &= good
        print(f"{'PASS' if good else 'FAIL'} v-space vs Khintchine two-step tau={tau} c=({c0},{c1}) r={r}: {float(got):.12f} vs {float(want):.12f}")
    # two-sided valley: mode inside
    got = prob_exact([0, Fraction(1, 4), Fraction(1, 2), 1], [3, 1, 2], Fraction(3))
    # density: 1/3 on [0, 3/4), 1 on [3/4, 1), 1/2 on [1, 2]; mode cell [3/4, 1): take m = 3/4
    # Khintchine about m: right side layers: 1/2 on [m, 2] -> radius 5/4 weight (1/2)(5/4);
    #   (1 - 1/2) on [m, 1] -> radius 1/4 weight (1/2)(1/4); left side: 1/3 on [0, m] -> radius -3/4 weight (1/3)(3/4)
    want = phi_exact([Fraction(5, 8), Fraction(1, 8), Fraction(1, 4)], [Fraction(5, 4), Fraction(1, 4), Fraction(-3, 4)], Fraction(3))
    good = got == want
    ok &= good
    print(f"{'PASS' if good else 'FAIL'} v-space vs Khintchine two-sided valley: {float(got):.12f} vs {float(want):.12f}")
    print("ALL VALIDATORS PASS" if ok else "VALIDATOR FAILURE")
    return ok


if __name__ == "__main__":
    sys.exit(0 if _validate() else 1)
