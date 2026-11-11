"""Step 16 (DT15-NEPRIME): rigorous ball-arithmetic enclosure of the density f = F_M' on cells of [Z0, Z1].

Renewal identity (P1 = Balabdaoui & Pitman 2011 Thm 1 + the GEM(1) stick-breaking property): with V ~ U(0,1) the first
stick, M = max(sqrt(V) E, sqrt(1-V) M'), E an excursion maximum (CDF F_E, density f_E), M' ~ M, all independent. So
    F(x) = E[F_E(x/sqrt V) F(x/sqrt(1-V))],
    f(x) = E[V^{-1/2} f_E(x/sqrt V) F(x/sqrt(1-V))] + int_x^inf F_E(x y / sqrt(y^2 - x^2)) f(y) 2x/y^2 dy.   (T)
One application of (T) to the provable sandwich of DT11 Lemma 1 (valid for every x > 0):
    f_lo0(y) = I1(y) (1 - J(y)) <= f(y) <= I1(y) = f_hi0(y),  I1 = 2S/y, J(y) = int_y^inf 2S(v)/v dv, S = 1 - F_E,
gives f_lo1 = T(f_lo0, F_lo) <= f <= T(f_hi0, F_hi) = f_hi1 (T is monotone: all kernels >= 0).
F_lo / F_hi: the certified Arb enclosure table `deposits/dt10/fm_arb_low.npz` (grid 0.5..8.5, floor / ceiling lookup).
Rigour: lower (upper) Riemann sums with Arb interval bounds of every integrand on every piece, outward rounding;
integrands' monotonicity is NOT assumed except where stated (F_E increasing, F increasing, xy/sqrt(y^2-x^2) decreasing
in y and increasing in x, I1 decreasing on [0.9, inf), J decreasing).
Output: deposits/dt16/fm_density.json  {cells: [[z_a, z_b]], lo: [...], hi: [...]} with lo <= f(z) <= hi on each cell.
Usage: python dt16_fm_density.py [Z0 Z1 H]   (default 0.95 1.40 0.0025)
"""
import json
import os
import sys
import time

import numpy as np
from flint import arb, ctx

ctx.prec = 80
K = 16
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
YTOP = 8.0          # integrals over y are cut at YTOP; the rest is bounded analytically (see tails below)


def S_ball(z):
    """Arb enclosure of S(z) = 1 - F_E(z) = 2 sum_k (4k^2z^2 - 1) e^{-2k^2z^2}, z >= 0.9 (alternating-free for z > 1/2)."""
    s = arb(0)
    for k in range(1, K + 1):
        kk = k * k
        s += (4 * kk * z * z - 1) * (-2 * kk * z * z).exp()
    k = K + 1
    rem = 2 * (4 * k * k * z * z) * (-2 * k * k * z * z).exp()
    v = 2 * s
    return arb(v.mid(), v.rad() + 2 * rem.upper())


def fE_ball(y):
    """Arb enclosure of f_E(y) = -S'(y) = -2 sum_k e^{-2k^2y^2} [8k^2 y - 4k^2 y (4k^2 y^2 - 1)], y >= 0.9."""
    s = arb(0)
    for k in range(1, K + 1):
        kk = k * k
        s += (4 * kk * y * (4 * kk * y * y - 1) - 8 * kk * y) * (-2 * kk * y * y).exp()
    k = K + 1
    rem = 2 * (16 * k ** 4 * y ** 3) * (-2 * k * k * y * y).exp()
    v = 2 * s
    return arb(v.mid(), v.rad() + 2 * rem.upper())


def interval(a, b):
    """Arb ball containing the real interval [a, b] (a, b floats or arbs)."""
    a, b = arb(a), arb(b)
    lo, hi = a.lower(), b.upper()
    return arb((lo + hi) / 2, (hi - lo) / 2)


def load_table():
    d = np.load(os.path.join(ROOT, "deposits", "dt10", "fm_arb_low.npz"))
    ys, lo, hi = np.asarray(d["ys"], float), np.asarray(d["lo"], float), np.asarray(d["hi"], float)
    if not np.all(np.diff(ys) > 0):
        sys.exit("fm table grid not increasing")
    return ys, np.maximum.accumulate(lo), np.minimum.accumulate(hi[::-1])[::-1]


def F_lo(tab, x):
    ys, lo, _ = tab
    k = int(np.searchsorted(ys, x, side="right")) - 1
    return 0.0 if k < 0 else max(0.0, float(lo[k]) - 1e-12)


def F_hi(tab, x):
    ys, _, hi = tab
    k = int(np.searchsorted(ys, x, side="left"))
    return 1.0 if k >= len(ys) else min(1.0, float(hi[k]) + 1e-12)


def build_J(z_lo, z_hi, h):
    """Upper and lower bounds of J(y) = int_y^inf 2S(v)/v dv on the grid y_i = z_lo + i h (integrand decreasing)."""
    n = int(round((YTOP - z_lo) / h))
    g = [arb(z_lo) + arb(h) * i for i in range(n + 1)]
    vals = [2 * S_ball(v) / v for v in g]
    # tail beyond YTOP: 2S(v)/v <= 2 * 2 * 4 v^2 e^{-2v^2} * 1.02 / v = 16.32 v e^{-2v^2}; integral = 4.08 e^{-2 Y^2}
    tail = 4.08 * float((-2 * arb(YTOP) ** 2).exp().upper())
    J_hi, J_lo = [0.0] * (n + 1), [0.0] * (n + 1)
    acc_hi, acc_lo = tail, 0.0
    for i in range(n - 1, -1, -1):
        acc_hi += float((vals[i] * arb(h)).upper())          # left endpoint = max of a decreasing integrand
        acc_lo += float((vals[i + 1] * arb(h)).lower())
        J_hi[i], J_lo[i] = acc_hi, acc_lo
    return [float(x.mid()) for x in g], J_lo, J_hi


def term2(x0, x1, yg, fy_lo, fy_hi):
    """Bounds of int_x^inf F_E(x y/sqrt(y^2-x^2)) f(y) 2x/y^2 dy for x in [x0, x1], using f in [fy_lo, fy_hi] per y-piece.
    Lower: drop y in (x, x1); kernel inf at (x0, right end of piece). Upper: kernel <= 1 * 2 x1 / y_a^2, plus (x0, x1) part."""
    lo, hi = arb(0), arb(0)
    for j in range(len(yg) - 1):
        ya, yb = yg[j], yg[j + 1]
        if yb <= x0:
            continue
        w = arb(yb) - arb(ya)
        if ya >= x1:
            arg = arb(x0) * arb(yb) / (arb(yb) ** 2 - arb(x0) ** 2).sqrt()
            ker_lo = (1 - S_ball(arg)) * 2 * arb(x0) / arb(yb) ** 2
            lo += ker_lo * arb(fy_lo[j]) * w
        hi += 2 * arb(x1) / arb(max(ya, x0)) ** 2 * arb(fy_hi[j]) * w
    hi += 2 * arb(x1) / arb(yg[-1]) * 1.0       # y > YTOP: f <= I1 <= 1 there, int 2x/y^2 = 2x/Y
    return lo, hi


def term1(x0, x1, tab, nv=600):
    """Bounds of E[V^{-1/2} f_E(x/sqrt V) F(x/sqrt(1-V))] for x in [x0, x1]: v-pieces, interval inf/sup of the integrand."""
    lo, hi = arb(0), arb(0)
    vs = np.concatenate([np.geomspace(1e-6, 0.05, nv // 3), np.linspace(0.05, 1.0, 2 * nv // 3 + 1)[1:]])
    vs = np.concatenate([[0.0], vs])
    X = interval(x0, x1)
    for j in range(len(vs) - 1):
        va, vb = vs[j], vs[j + 1]
        if va <= 0.0:
            # v in (0, 1e-6]: v^{-1/2} f_E(x/sqrt v) <= v^{-1/2} 64 (x^2/v)^{3/2} e^{-2x^2/v} < 1e-100 for x >= 0.9;
            # contributes [0, 1e-100 * 1e-6] (lower bound: drop it)
            hi += arb("1e-100")
            continue
        V = interval(va, vb)
        g = fE_ball(X / V.sqrt()) / V.sqrt()
        glo, ghi = max(0.0, float(g.lower())), float(g.upper())
        flo = F_lo(tab, x0 / np.sqrt(1 - va))
        fhi = 1.0 if vb >= 1.0 else F_hi(tab, x1 / np.sqrt(1 - vb))
        w = arb(vb) - arb(va)
        lo += arb(glo) * arb(flo) * w
        hi += arb(ghi) * arb(fhi) * w
    # v in (0, 1e-6): integrand <= v^{-1/2} f_E(x/sqrt v) <= v^{-1/2} * 64 (x^2/v)^{3/2} e^{-2 x^2/v}: negligible; bound by 1e-300
    return lo, hi


def main():
    # usage: python dt16_fm_density.py [H NIT ZF]: cells of width H on [0.9, ZF), NIT renewal steps (default .0025 2 3.0)
    h = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0025
    nit = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    zf = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
    t0 = time.time()
    tab = load_table()
    nf = int(round((zf - 0.9) / h))
    yg = [0.9 + i * h for i in range(nf + 1)] + list(np.arange(zf + 0.05, YTOP + 1e-9, 0.05))
    _, _, J_hi = build_J(0.9, 0.9, h)

    def jhi(y):                                    # J decreasing: J(y) <= J_hi at the grid point <= y
        return J_hi[min(len(J_hi) - 1, max(0, int(np.floor((y - 0.9) / h + 1e-9))))]

    I1_lo = [float((2 * S_ball(arb(y)) / arb(y)).lower()) for y in yg]
    I1_hi = [float((2 * S_ball(arb(y)) / arb(y)).upper()) for y in yg]
    # step 0 (DT11 Lemma 1, any x > 0): lo0 = I1(y_b)(1 - J(y_a)) on each piece (I1 decreasing, 1 - J increasing)
    fy_lo = [max(0.0, I1_lo[j + 1] * (1 - jhi(yg[j]))) for j in range(len(yg) - 1)]
    fy_hi = [I1_hi[j] for j in range(len(yg) - 1)]
    hist = []
    LO, HI = [], []
    for it in range(nit):
        LO, HI = [], []
        for i in range(nf):
            a, b = yg[i], yg[i + 1]
            l1, h1 = term1(a, b, tab)
            l2, h2 = term2(a, b, yg, fy_lo, fy_hi)
            LO.append(max(fy_lo[i], float((l1 + l2).lower())))
            HI.append(min(fy_hi[i], float((h1 + h2).upper())))
        fy_lo[:nf], fy_hi[:nf] = LO, HI                   # T is monotone: improved bounds feed the next step
        k1 = int(round((1.0 - 0.9) / h))
        hist.append({"iter": it + 1, "ratio_at_1.0": LO[k1] / HI[k1], "sec": time.time() - t0})
        print(f"iteration {it + 1}: lo/hi at z = 1.0: {LO[k1] / HI[k1]:.4f}  ({time.time() - t0:.0f}s)", flush=True)
    cells = [[yg[i], yg[i + 1]] for i in range(nf)]
    n = nf
    out = {"h": h, "iterations": nit, "zf": zf, "cells": cells, "lo": LO, "hi": HI, "prec": ctx.prec, "K": K,
           "history": hist, "table": "deposits/dt10/fm_arb_low.npz", "sec": time.time() - t0}
    path = os.path.join(ROOT, "deposits", "dt16", "fm_density.json")
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            json.dump(out, fh)
    except OSError as err:
        sys.exit(f"cannot write {path}: {err}")
    for i in range(0, n, max(1, n // 9)):
        print(f"z=[{cells[i][0]:.4f},{cells[i][1]:.4f}] lo={LO[i]:.5f} hi={HI[i]:.5f} ratio={LO[i] / HI[i]:.4f}")
    print(f"{n} cells, {out['sec']:.1f}s -> {path}")


if __name__ == "__main__":
    main()
