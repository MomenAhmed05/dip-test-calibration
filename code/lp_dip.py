"""Independent dip by linear programming, written from work/CONVENTIONS.md (CV-1, CV-2, CV-7).

Object: D_U(F_n) = inf over G in U of sup_x |F_n(x) - G(x)|, where U is H&H's null: G is a CDF that
is convex on (-inf, m) and concave on (m, inf) for some m, with a jump allowed only at m.

Reduction (checked in the Step 2 smoke, not assumed elsewhere). Let v_1 < ... < v_k be the distinct
sample values, F_j = F_n(v_j) and F_j^- = F_n(v_j-). Only G's left limits L_j = G(v_j-) and values
R_j = G(v_j) matter, because F_n is constant between knots and G may be taken linear there. Segment
slopes s_j = (L_{j+1} - R_j) / (v_{j+1} - v_j) must be >= 0 and unimodal. Two mode types:
  * peak segment p (mode inside (v_p, v_{p+1})): G continuous at every knot; s nondecreasing up to p,
    nonincreasing from p;
  * mode at knot j: G may jump at v_j only (R_j - L_j = atom >= 0); slopes left of v_j nondecreasing,
    slopes right of v_j nonincreasing. (atom = 0 also covers a slope that jumps up at the mode.)
Each type is an LP in (L, R, t): minimise t subject to |L_j - F_j^-| <= t, |R_j - F_j| <= t.
D_U = min over all mode types. D_cont (continuous G only, a deliberately different comparator) uses
the same LPs with every atom forced to 0. Ties are handled through the multiplicities.
Slow: about 2k LPs per sample. For audits at small n only.
"""
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


def _knots(x):
    v, c = np.unique(np.asarray(x, float), return_counts=True)
    n = c.sum()
    Fr = np.cumsum(c) / n           # F_n(v_j)
    Fl = Fr - c / n                 # F_n(v_j-)
    return v, Fl, Fr, n


def _solve(v, Fl, Fr, mode_kind, idx, allow_atom):
    k = len(v)
    h = np.diff(v)
    # variable layout: L_0..L_{k-1}, R_0..R_{k-1}, t   (0-based)
    nv = 2 * k + 1
    T = 2 * k
    rows, cols, vals, b = [], [], [], []
    r = 0

    def add(entries, rhs):
        nonlocal r
        for cidx, val in entries:
            rows.append(r)
            cols.append(cidx)
            vals.append(val)
        b.append(rhs)
        r += 1

    for j in range(k):
        L, R = j, k + j
        add([(L, 1), (T, -1)], Fl[j])          # L_j - t <= F_j^-
        add([(L, -1), (T, -1)], -Fl[j])        # -L_j - t <= -F_j^-
        add([(R, 1), (T, -1)], Fr[j])          # R_j - t <= F_j
        add([(R, -1), (T, -1)], -Fr[j])        # -R_j - t <= -F_j
        add([(L, 1), (R, -1)], 0.0)            # L_j <= R_j
    # slopes s_j = (L_{j+1} - R_j)/h_j >= 0  ->  R_j - L_{j+1} <= 0
    for j in range(k - 1):
        add([(k + j, 1), (j + 1, -1)], 0.0)

    def slope(j):  # coefficients of s_j
        return [(j + 1, 1.0 / h[j]), (k + j, -1.0 / h[j])]

    def le(sa, sb):  # s_a <= s_b
        add([(c, val) for c, val in slope(sa)] + [(c, -val) for c, val in slope(sb)], 0.0)

    if mode_kind == "seg":
        p = idx
        for j in range(k - 2):
            if j < p:
                le(j, j + 1)
            else:
                le(j + 1, j)
    else:  # mode at knot idx: segments 0..idx-1 are left of v_idx, idx..k-2 right of it
        m = idx
        for j in range(k - 2):
            if j + 1 <= m - 1:
                le(j, j + 1)
            elif j >= m:
                le(j + 1, j)
    # equality constraints: continuity at every knot except the mode knot (if atoms allowed)
    Ae_r, Ae_c, Ae_v, be = [], [], [], []
    er = 0
    for j in range(k):
        if mode_kind == "knot" and j == idx and allow_atom:
            continue
        Ae_r += [er, er]
        Ae_c += [j, k + j]
        Ae_v += [1.0, -1.0]
        be.append(0.0)
        er += 1
    A = coo_matrix((vals, (rows, cols)), shape=(r, nv)).tocsr()
    Ae = coo_matrix((Ae_v, (Ae_r, Ae_c)), shape=(er, nv)).tocsr() if er else None
    cobj = np.zeros(nv)
    cobj[T] = 1.0
    bounds = [(0.0, 1.0)] * (2 * k) + [(0.0, 1.0)]
    res = linprog(cobj, A_ub=A, b_ub=np.array(b), A_eq=Ae, b_eq=np.array(be) if er else None,
                  bounds=bounds, method="highs")
    if res.status != 0:
        return np.inf
    return float(res.fun)


def dip_lp(x, allow_atom=True):
    """D_U (allow_atom=True) or D_cont (allow_atom=False) of the sample x."""
    v, Fl, Fr, n = _knots(x)
    k = len(v)
    if k == 1:
        return 0.0 if allow_atom else 0.5
    best = np.inf
    for p in range(k - 1):
        best = min(best, _solve(v, Fl, Fr, "seg", p, allow_atom))
    for m in range(k):
        best = min(best, _solve(v, Fl, Fr, "knot", m, allow_atom))
    return best
