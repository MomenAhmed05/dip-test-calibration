"""DT10 exact-predicate kernels (blocker: one-sided numerics for the certificate).

The bridge skeleton is rounded to integers B_k = rint(b_k * 2^40) (|B_k| <= 2^43 required, else the caller sets
V = 1). The piecewise-linear bridge with values B_k / 2^40 differs from the true bridge by at most
E' = E + 2^-40 (E = exact interpolation error; rounding <= 2^-41 plus float generation error << 2^-41).
All hull and tangent decisions use exact int64 arithmetic:
  hull cross products  |dB| * |dk| <= 2^44 * 2^17 = 2^61  (difference <= 2^62);
  tangent comparisons  |qy - hy| * |dk| <= 2^45.5 * 2^17 < 2^63  (requires 2y <= 16; callers have y <= 7.2).
Values are then evaluated from exact integers with <= ~1e-14 absolute float error; callers add 1e-12.
"""
import numpy as np
from numba import njit

SCALE = 2.0 ** 40
BMAX = 2 ** 43


@njit(cache=True)
def gcm_gap_int(B):
    """G(b) = max_k (b_k - GCM(b)(k)) for the PL bridge with integer values B (unit 2^-40), up to ~1e-14."""
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
        num = (B[i] - hy[j]) * dx - (hy[j + 1] - hy[j]) * (i - hx[j])     # exact; = dx * (b_i - hull_i) * 2^40
        v = float(num) / float(dx) / SCALE
        if v > g:
            g = v
    return g


@njit(cache=True)
def lam_int(B, Y2):
    """Lambda_{b_int}(y') with 2y' = Y2 / 2^40 (Y2 integer): max(0, sup_{s<x<1} [b(x) + (1-x)(b(x)-b(s)-2y')/(x-s)]),
    sup over vertices (Lemma F(e)); exact hull/tangent predicates; value error <= ~1e-14."""
    n = len(B)
    m = n - 1
    hx = np.empty(n, np.int64)
    hy = np.empty(n, np.int64)
    k = 0
    best = 0.0
    for j in range(n - 1):              # x = j/m < 1
        if k >= 1 and j >= 1:
            qy = B[j] - Y2
            lo = 0
            hi = k - 1
            while lo < hi:                  # max slope from hull vertex to Q = (j, qy): first index where it stops rising
                mid = (lo + hi) // 2
                if (qy - hy[mid]) * (j - hx[mid + 1]) >= (qy - hy[mid + 1]) * (j - hx[mid]):
                    hi = mid
                else:
                    lo = mid + 1
            N = qy - hy[lo]
            D = j - hx[lo]
            # (1-x)/(x-s) = (m-j)/D, so b(x) + (1-x)(b(x)-b(s)-2y')/(x-s) = B_j/2^40 + (m-j) N / (D 2^40)
            v = float(B[j]) / SCALE + float((m - j) * N) / (float(D) * SCALE)
            if v > best:
                best = v
        while k >= 2 and (hy[k - 1] - hy[k - 2]) * (j - hx[k - 2]) >= (B[j] - hy[k - 2]) * (hx[k - 1] - hx[k - 2]):
            k -= 1
        hx[k] = j
        hy[k] = B[j]
        k += 1
    return best


def to_int(b):
    B = np.rint(np.asarray(b) * SCALE).astype(np.int64)
    return B, bool(np.max(np.abs(B)) <= BMAX)
