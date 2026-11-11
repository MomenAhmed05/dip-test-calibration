"""Step 15 / DT14 WP-B1: the n = 4 kernel engine.

Object (A-13): Hartigan dip at n = 4 via the H&H (1985) reduction D_4 = T/(4(T+1)),
T = s2 / max(s1, s3) (s_i spacings of the ordered sample); the dip test rejects at the uniform
alpha-quantile iff T >= r_alpha = (sqrt(1 + 8/alpha) - 3)/2 (T >= r with r > 1; ties have
probability 0 for continuous laws).  [B0 re-proves the reduction separately.]

Khintchine: a unimodal law with mode 0 is the law of U*Z (U ~ U(0,1) independent of Z ~ G).
For one sample X_i = U_i z_i (i = 1..4, U_i iid U(0,1), fixed nonzero z_i):
    K_r(z_1..z_4) = P(T >= r)
                  = sum over the 24 orders sigma of vol{u in [0,1]^4 : x_s1 <= x_s2 <= x_s3 <= x_s4,
                                                      s2 >= r s1, s2 >= r s3}
with x_i = z_i u_i, linear in u.  Then for G = sum_j w_j delta_{z_j}:
    Phi_r(G) = P_F(T >= r) = sum_{i,j,k,l} w_i w_j w_k w_l K_r(z_i, z_j, z_k, z_l).

Two engines:
  kernel_exact(z, r): rational z and r (Fractions), exact volumes via exact_poly.volume (no floats).
  kernel_float(z, r): float vertex enumeration + qhull volume (DEV searches only; A-14: never a CERT).
Validators (run: python dt14_n4_kernel.py): see _validate().
"""
import sys
from fractions import Fraction
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations
from math import factorial

import numpy as np
from scipy.spatial import ConvexHull

sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
from exact_poly import volume as exact_volume  # noqa: E402

ORDERS = list(permutations(range(4)))


def _constraints(z, r, order):
    """Rows (A, b) of A u <= b for one order, u in [0,1]^4, x_i = z_i u_i."""
    A, b = [], []
    zero = z[0] * 0
    one = zero + 1
    for i in range(4):  # 0 <= u_i <= 1
        e = [zero] * 4
        e[i] = -one
        A.append(e)
        b.append(zero)
        e = [zero] * 4
        e[i] = one
        A.append(e)
        b.append(one)

    def x(i):  # coefficient row of x_i
        e = [zero] * 4
        e[i] = z[i]
        return e

    def sub(p, q):
        return [pi - qi for pi, qi in zip(p, q, strict=True)]

    s = [sub(x(order[k + 1]), x(order[k])) for k in range(3)]  # spacings as rows
    for k in range(3):  # s_k >= 0  <=>  -s_k <= 0
        A.append([-v for v in s[k]])
        b.append(zero)
    # s2 >= r s1  <=>  r s1 - s2 <= 0 ; s2 >= r s3 likewise (s indices 0,1,2 = s1,s2,s3)
    A.append([r * a - c for a, c in zip(s[0], s[1], strict=True)])
    b.append(zero)
    A.append([r * a - c for a, c in zip(s[2], s[1], strict=True)])
    b.append(zero)
    return A, b


def kernel_exact(z, r):
    z = [Fraction(v) for v in z]
    r = Fraction(r)
    if any(v == 0 for v in z):
        raise ValueError("z must be nonzero (continuous laws)")
    tot = Fraction(0)
    for o in ORDERS:
        A, b = _constraints(z, r, o)
        tot += exact_volume(A, b)
    return tot


_COMBOS = list(combinations(range(13), 4))


def _vol_float(A, b):
    A = np.asarray(A, float)
    b = np.asarray(b, float)
    idx = np.array(_COMBOS)
    M = A[idx]  # (C,4,4)
    y = b[idx]
    det = np.linalg.det(M)
    ok = np.abs(det) > 1e-13
    if not ok.any():
        return 0.0
    sol = np.linalg.solve(M[ok], y[ok][..., None])[..., 0]
    feas = np.all(sol @ A.T <= b + 1e-11, axis=1)
    V = sol[feas]
    if len(V) < 5:
        return 0.0
    V = np.unique(np.round(V, 12), axis=0)
    if len(V) < 5:
        return 0.0
    try:
        return float(ConvexHull(V).volume)
    except Exception:  # degenerate (lower-dimensional) polytope
        return 0.0


def kernel_float(z, r):
    z = [float(v) for v in z]
    tot = 0.0
    for o in ORDERS:
        A, b = _constraints(z, float(r), o)
        tot += _vol_float(A, b)
    return tot


@cache
def _kf_sorted(zs, r):
    return kernel_float(zs, r)


def kf(z, r):
    """Cached float kernel (K is symmetric and scale invariant: normalise by max |z|, then sort).
    Float engine guard: refuse within-tuple ratios below 1e-7 (use kernel_exact there)."""
    m = max(abs(float(v)) for v in z)
    zn = tuple(sorted(round(float(v) / m, 15) for v in z))
    if min(abs(v) for v in zn) < 1e-7:
        raise ValueError("float kernel: within-tuple scale ratio < 1e-7; use kernel_exact")
    return _kf_sorted(zn, float(r))


def phi_float(weights, zs, r):
    """Phi_r(G) for G = sum_j w_j delta_{z_j} (weights normalised here)."""
    w = np.asarray(weights, float)
    w = w / w.sum()
    m = len(zs)
    tot = 0.0
    for c in combinations_with_replacement(range(m), 4):
        cnt = np.bincount(c, minlength=m)
        mult = factorial(4) / np.prod([factorial(k) for k in cnt])
        tot += mult * np.prod(w[list(c)]) * kf([zs[i] for i in c], r)
    return tot


def phi_exact(weights, zs, r):
    w = [Fraction(v) for v in weights]
    s = sum(w)
    w = [v / s for v in w]
    m = len(zs)
    tot = Fraction(0)
    for c in combinations_with_replacement(range(m), 4):
        cnt = [c.count(i) for i in range(m)]
        mult = factorial(4)
        for k in cnt:
            mult //= factorial(k)
        p = Fraction(mult)
        for i in c:
            p *= w[i]
        tot += p * kernel_exact([zs[i] for i in c], r)
    return tot


# ------------------------------------------------------------------ closed forms used by validators
def a_of(r):
    return 2 / ((r + 1) * (r + 2))


def b_of(r):
    return 1 / (r + 1)


def gamma_of(r):  # (input checked here)
    return (2 * r + 1) / (3 * r * (r + 1) ** 2)


def delta_of(r):
    return (3 * r + 2) / (4 * r * (r + 1) ** 2 * (r + 2))


def eta_of(r):
    return (14 * r * r + 11 * r - 2) / (3 * r * (r + 1) ** 2 * (r + 2))


def _validate():
    rows = []
    r = Fraction(3)  # alpha = 0.10 exactly
    for rr in (Fraction(3), Fraction(5), Fraction(13, 2)):
        rows.append(("uniform K(1,1,1,1) = a(r)", rr, kernel_exact([1, 1, 1, 1], rr), a_of(rr)))
        rows.append(("4 coarse +++- = delta(r)", rr, kernel_exact([1, 1, 1, -1], rr), delta_of(rr)))
        rows.append(("4 coarse ++-- = eta(r)", rr, kernel_exact([1, 1, -1, -1], rr), eta_of(rr)))
        rows.append(("4delta+3eta = 7a", rr, 4 * delta_of(rr) + 3 * eta_of(rr), 7 * a_of(rr)))
        # rescaled copies: scale invariance + reflection
        rows.append(("scale/reflect K(-2,-2,-2,-2) = a", rr, kernel_exact([-2, -2, -2, -2], rr), a_of(rr)))
    ok = True
    for name, rr, got, want in rows:
        good = got == want
        ok &= good
        print(f"{'PASS' if good else 'FAIL'} {name} r={rr}: got {got} want {want}")
    # separated limits (float, t -> 0): must approach the limit values
    for rr in (3.0, 5.0):
        for t in (1e-3, 1e-5):
            cases = [("m=3 one finer K(1,1,1,t) -> a", [1, 1, 1, t], a_of(rr)),
                     ("m=2 K(1,1,t,t) -> b", [1, 1, t, t], b_of(rr)),
                     ("m=1 K(1,t,t,t) -> 0", [1, t, t, t], 0.0),
                     ("signed m=3 ++- + finer -> gamma", [1, 1, -1, t], gamma_of(rr)),
                     ("signed m=2 +- + finer finer -> 0", [1, -1, t, t], 0.0)]
            for name, z, want in cases:
                got = kernel_float(z, rr)
                good = abs(got - want) < 50 * t + 1e-9
                ok &= good
                print(f"{'PASS' if good else 'FAIL'} {name} r={rr} t={t}: got {got:.8f} want {want:.8f}")
    # float vs exact agreement on random rational tuples
    rng = np.random.default_rng(2026110801)  # DT14 DEV root 0801 (validator only)
    worst = 0.0
    for _ in range(12):
        z = [Fraction(int(rng.integers(1, 40)), int(rng.integers(1, 40))) * (1 if rng.random() < 0.6 else -1) for _ in range(4)]
        e = kernel_exact(z, r)
        f = kernel_float(z, float(r))
        worst = max(worst, abs(float(e) - f))
    good = worst < 1e-9
    ok &= good
    print(f"{'PASS' if good else 'FAIL'} float vs exact on 12 random rational tuples: max |diff| = {worst:.2e}")
    print("ALL VALIDATORS PASS" if ok else "VALIDATOR FAILURE")
    return ok


if __name__ == "__main__":
    sys.exit(0 if _validate() else 1)
