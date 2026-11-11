"""Pinned-mode dip of a sample on [0, 1] (Step 10 / DT8, DEV).

Objects (DT7-UNIF.md section 1). For a sample x in [0, 1] with empirical CDF G_n and a fixed u in [0, 1]:
  Y(u) = d(G_n, K_u),  K_u = {g convex on [0, u], concave on [u, 1], g(u-) <= g(u+)}  (JOINED; an
         upward jump = atom at the mode is allowed), constant outside [0, 1];
  Z(u) = max(V_[0,u](G_n), C_[u,1](G_n))  (DECOUPLED branches; Z <= Y);
  V1   = V_[0,1](G_n) = Y(1) (edge; limit M/2 after scaling by sqrt(n)).
Distances are sup-norm over both one-sided values of G_n (G_n(x-) and G_n(x)).
Reduction: G_n is constant between knots, so g may be taken piecewise linear on the knots
{0, x_1, ..., x_n, 1, u}; the error on an open gap is attained at its ends.
Y(u) is an LP (HiGHS). V and C use the half-gap identity on the knots (Lemma 3 of DT3-THEOREM):
  V = 1/2 max_k (G_hi - GCM(G_lo))_k,  C = 1/2 max_k (LCM(G_hi) - G_lo)_k.
"""
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


def _hull(t, y, lower=True):
    """Indices of the lower (GCM) or upper (LCM) convex hull of points (t, y), t increasing."""
    idx = []
    for i in range(len(t)):
        while len(idx) >= 2:
            a, b = idx[-2], idx[-1]
            cross = (t[b] - t[a]) * (y[i] - y[a]) - (y[b] - y[a]) * (t[i] - t[a])
            if (cross <= 0) if lower else (cross >= 0):
                idx.pop()
            else:
                break
        idx.append(i)
    return np.array(idx)


def _envelope(t, y, lower=True):
    h = _hull(t, y, lower)
    return np.interp(t, t[h], y[h])


def _side_knots(x, u, n, left=True):
    """Knot times and (lo, hi) values of G_n on [0, u] (left) or [u, 1] (right)."""
    k = int(np.searchsorted(x, u, side="right"))
    c = k / n
    if left:
        xs = x[:k]
        t = np.concatenate([[0.0], xs, [u]])
        hi = np.concatenate([[0.0], np.arange(1, k + 1) / n, [c]])
        lo = np.concatenate([[0.0], np.arange(0, k) / n, [c]])
    else:
        xs = x[k:]
        m = len(xs)
        t = np.concatenate([[u], xs, [1.0]])
        hi = np.concatenate([[c], (k + np.arange(1, m + 1)) / n, [1.0]])
        lo = np.concatenate([[c], (k + np.arange(0, m)) / n, [1.0]])
    return t, lo, hi


def V_left(x, u):
    n = len(x)
    if u <= 0:
        return 0.0
    t, lo, hi = _side_knots(x, u, n, True)
    return 0.5 * float(np.max(hi - _envelope(t, lo, True)))


def C_right(x, u):
    n = len(x)
    if u >= 1:
        return 0.0
    t, lo, hi = _side_knots(x, u, n, False)
    return 0.5 * float(np.max(_envelope(t, hi, False) - lo))


def Z_dec(x, u):
    return max(V_left(x, u), C_right(x, u))


def Y_joined(x, u):
    """LP: minimise e over g on knots; left convex, right concave, a = g(u-) <= b = g(u+)."""
    x = np.sort(np.asarray(x, float))
    n = len(x)
    if u >= 1:
        return V_left(x, 1.0)
    if u <= 0:
        return C_right(x, 0.0)
    tL, loL, hiL = _side_knots(x, u, n, True)
    tR, loR, hiR = _side_knots(x, u, n, False)
    nL, nR = len(tL), len(tR)
    # variables: gL[0..nL-1] (last = a), gR[0..nR-1] (first = b), e
    E = nL + nR
    rows, cols, vals, b = [], [], [], []
    r = 0

    def add(coefs, rhs):
        nonlocal r
        for c_, v_ in coefs:
            rows.append(r); cols.append(c_); vals.append(v_)
        b.append(rhs); r += 1

    for j in range(nL):          # hi - e <= g <= lo + e
        add([(j, -1.0), (E, -1.0)], -hiL[j])
        add([(j, 1.0), (E, -1.0)], loL[j])
    for j in range(nR):
        add([(nL + j, -1.0), (E, -1.0)], -hiR[j])
        add([(nL + j, 1.0), (E, -1.0)], loR[j])
    for j in range(nL - 2):      # convex: s_j <= s_{j+1}
        h1, h2 = tL[j + 1] - tL[j], tL[j + 2] - tL[j + 1]
        if h1 <= 0 or h2 <= 0:
            continue
        add([(j, -1.0 / h1), (j + 1, 1.0 / h1 + 1.0 / h2), (j + 2, -1.0 / h2)], 0.0)
    for j in range(nR - 2):      # concave: s_j >= s_{j+1}
        h1, h2 = tR[j + 1] - tR[j], tR[j + 2] - tR[j + 1]
        if h1 <= 0 or h2 <= 0:
            continue
        o = nL
        add([(o + j, 1.0 / h1), (o + j + 1, -1.0 / h1 - 1.0 / h2), (o + j + 2, 1.0 / h2)], 0.0)
    add([(nL - 1, 1.0), (nL, -1.0)], 0.0)   # a <= b
    A = coo_matrix((vals, (rows, cols)), shape=(r, E + 1)).tocsr()
    cvec = np.zeros(E + 1); cvec[E] = 1.0
    res = linprog(cvec, A_ub=A, b_ub=np.array(b), bounds=[(None, None)] * E + [(0, None)],
                  method="highs")
    if res.status != 0:
        raise RuntimeError(f"LP failed: {res.message}")
    return float(res.x[E])


if __name__ == "__main__":
    import sys
    sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
    rng = np.random.default_rng(2026110700)   # DEV root (CONVENTIONS: 2026110700 = DT8 smoke)
    # 1. polygon check (deterministic): joined 2a/3 vs decoupled a/2 is for continuous h;
    # here check basic identities on samples instead.
    import diptest
    worst = 0.0
    for rep in range(20):
        x = np.sort(rng.uniform(size=60))
        us = np.unique(np.concatenate([np.linspace(0, 1, 41), (x[:-1] + x[1:]) / 2]))
        ys = [Y_joined(x, u) for u in us]
        d = float(diptest.dipstat(x))
        worst = max(worst, abs(min(ys) - d))
        assert abs(Y_joined(x, 1.0) - V_left(x, 1.0)) < 1e-12
        for u in (0.3, 0.7):
            assert Y_joined(x, u) >= Z_dec(x, u) - 1e-9
    print("min_u Y vs diptest max abs diff over 20 samples:", worst)
