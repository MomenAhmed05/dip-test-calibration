"""Step 16 (DT15) A1': validated (Arb) check of Theorem NE' at the single target t' per level.

Bound (DT15-NEPRIME §4-§7, DT11 §5):  g(1 - eps^2, t) - Q(2t) <= H(eps, t) :=
    2 S(lam t/eps) + E psi(delta X) - [F(x) - F(2t)],  x = 2t/s, delta = eps/s, s = sqrt(1 - eps^2),
    P(X > r) <= T(r) := exp(-2 r (r + kappa)), kappa = (2 - lam) t/eps,
    psi(w) <= min{ F(x) - F(x - w),  w [d(x - w) Hg(w) + w d(x - w) d(x - 2w)] }  (2nd branch only if x - 2w >= ZD),
    d(z) = z^-3 sup_{v >= z} v^3 f_E(v), Hg = certified upper step bound of H_geom (dt16_hgeom_cert.py, MC, delta_H).
Inputs: F' enclosure step cells (deposits/dt16/fm_density.json, code/dt16_fm_density.py), widened outward by REL.
(1) Boxes [ea, eb] in [EPS1, EPS_MAX]: H <= 2S(lam t/eb) + sum_j T_b(r_j)[g(r_{j+1}) - g(r_j)] + T_b(r_J)
    - [CL(x_a) - CL(2t)],  g(r) = psibar(delta_b r) non-decreasing (Stieltjes / layer cake, X dominated by T_b),
    psibar on x in [x_a, x_b]: min{CU(x_b) - CU(x_a - w), w[D(x_a - w) Hg(w) + w D(x_a - w) D(x_a - 2w)]}.
(2) (0, EPS1] analytically: Z = X/eps has P(Z > z) <= exp(-2 a z), a = (2 - lam)t; W = delta X = eps^2 Z/s;
    H/eps^2 <= 16.01 h^4 e^{-2h^2}/(lam t)^2 + L Hg(eta)/(2 a s1) + L^2 EPS1^2/(2 a^2 s1^2) + e^{-B/EPS1^2}/EPS1^2 - t m,
    h = lam t/EPS1 >= 2, L = D(2t - 2 eta), B = 2 a eta s1 >= EPS1^2, m = min lo-cell on [2t, 2t/s1].
All arithmetic in Arb balls (python-flint); decisions use certified upper bounds < 0.
Usage: python dt16_neprime_arb.py HGEOM_JSON tag   -> deposits/dt16/neprime_arb_<tag>.json
"""
import json
import math
import os
import sys
import time

from flint import arb, ctx

ctx.prec = 120
ROOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
LEVELS = ((0.10, "0.4935"), (0.05, "0.5378"), (0.01, "0.6161"))
REL = 1e-9                       # outward widening of the stored float cell bounds (float conversions upstream)
ZD = 0.5                         # D(z) is computed for z >= ZD only
EPS1 = 0.02
EPS_MAX = 0.30
RATIO = 1.004
LAMS = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)


def load_density():
    with open(os.path.join(ROOT_DIR, "deposits", "dt16", "fm_density.json")) as fh:
        d = json.load(fh)
    edges = [c[0] for c in d["cells"]] + [d["cells"][-1][1]]
    lo = [max(0.0, v) * (1 - REL) for v in d["lo"]]
    hi = [v * (1 + REL) for v in d["hi"]]
    PL, PU = [arb(0)], [arb(0)]
    for k in range(len(lo)):
        w = arb(edges[k + 1]) - arb(edges[k])
        PL.append(PL[-1] + arb(lo[k]) * w)
        PU.append(PU[-1] + arb(hi[k]) * w)
    return edges, lo, hi, PL, PU


EDGES, LO, HI, PL, PU = load_density()
Z0, ZF, HCELL = EDGES[0], EDGES[-1], EDGES[1] - EDGES[0]


def cell_of(z):
    k = int(math.floor((z - Z0) / HCELL))
    k = min(max(k, 0), len(LO) - 1)
    while k > 0 and z < EDGES[k]:
        k -= 1
    while k < len(LO) - 1 and z >= EDGES[k + 1]:
        k += 1
    return k


def CU(z, up):
    """Upper bound of int_{Z0}^{z} f for a float z (cumulative of the hi step); up=False for the lower-endpoint use."""
    if z < Z0 or z > ZF:
        raise ValueError("outside enclosure")
    k = cell_of(z)
    v = PU[k] + arb(HI[k]) * (arb(z) - arb(EDGES[k]))
    return arb(v.upper()) if up else arb(v.lower())


def CL(z, up):
    if z < Z0 or z > ZF:
        raise ValueError("outside enclosure")
    k = cell_of(z)
    v = PL[k] + arb(LO[k]) * (arb(z) - arb(EDGES[k]))
    return arb(v.upper()) if up else arb(v.lower())


def fE_ball(v):
    """f_E(v) = 2 sum_k 4k^2 v (4k^2 v^2 - 3) exp(-2 k^2 v^2), v >= 0.5 ball; tail k > 40 bounded by 1e-100."""
    s = arb(0)
    for k in range(1, 41):
        kk = arb(k * k)
        s += 4 * kk * v * (4 * kk * v * v - 3) * (-2 * kk * v * v).exp()
    return 2 * s + arb(0, 1e-100)


def build_D(hw=0.002, wmax=8.0):
    """Suffix maxima of upper bounds of v^3 f_E(v) on cells [w_j, w_j + hw], w_0 = ZD; beyond wmax: <= 4e-49 (§ doc)."""
    n = int(round((wmax - ZD) / hw))
    ups = []
    for j in range(n):
        a, b = ZD + j * hw, ZD + (j + 1) * hw
        v = arb((a + b) / 2, (b - a) / 2 + 1e-15)
        ups.append(float((v ** 3 * fE_ball(v)).upper()) * (1 + REL) + 1e-300)
    suf = [4e-49] * (n + 1)
    for j in range(n - 1, -1, -1):
        suf[j] = max(ups[j], suf[j + 1])
    return hw, suf


DHW, DSUF = build_D()


def Dz(z):
    """Upper bound of d(z) = z^-3 sup_{v >= z} v^3 f_E(v) (float z >= ZD)."""
    if z < ZD:
        return math.inf
    j = min(int((z - ZD) / DHW), len(DSUF) - 1)
    return DSUF[j] / (z * (1 - 1e-15)) ** 3 * (1 + 1e-12)


def S_up(h):
    """S(h) = 2 sum (4k^2h^2 - 1) e^{-2k^2h^2} (arb h >= 1), k <= 40 + tail."""
    s = arb(0)
    for k in range(1, 41):
        kk = arb(k * k)
        s += (4 * kk * h * h - 1) * (-2 * kk * h * h).exp()
    return 2 * s + arb(0, 1e-100)


def load_hg(path):
    with open(path) as fh:
        d = json.load(fh)
    return list(zip(d["rs"], d["bound"], strict=True)), d


def Hg_up(w, HG):
    for r, b in HG:
        if w <= r:
            return b
    return 1.0


def box_H(t, ea, eb, lam, HG):
    """Certified upper bound of H(eps, t) over eps in [ea, eb] (arb)."""
    ta, ea_, eb_ = arb(t), arb(ea), arb(eb)
    sa, sb = (1 - ea_ * ea_).sqrt(), (1 - eb_ * eb_).sqrt()
    xa, xb = 2 * ta / sa, 2 * ta / sb
    xa_lo, xb_hi = float(xa.lower()), float(xb.upper())
    db = float((eb_ / sb).upper())
    kap = float(((2 - arb(lam)) * ta / eb_).lower())
    h = arb(lam) * ta / eb_
    if float(h.lower()) < 1.0:
        return None
    term_s = 2 * S_up(arb(h.lower()))
    rmax = (xa_lo - Z0) / db                               # beyond: x_a - w < Z0, cap psi <= 1
    nr = 400
    rtop = min(rmax, 8.0 / math.sqrt(kap) + 12.0 / kap)
    rs = [rtop * (i / nr) ** 1.5 for i in range(nr + 1)]
    cu_xb = CU(min(xb_hi, ZF), True)

    def g(r):
        w = db * r
        if w <= 0:
            return arb(0)
        b1 = cu_xb - CU(xa_lo - w, False) if xa_lo - w >= Z0 else arb(1)
        if xa_lo - 2 * w >= ZD:
            d1, d2 = Dz(xa_lo - w), Dz(xa_lo - 2 * w)
            b2 = arb(w) * (arb(d1) * arb(Hg_up(w, HG)) + arb(w) * arb(d1) * arb(d2))
            if float(b2.upper()) < float(b1.upper()):
                return arb(b2.upper())
        return arb(b1.upper())

    def T(r):
        return (-2 * arb(r) * (arb(r) + arb(kap))).exp()

    gain = arb(0)
    gprev = arb(0)
    for j in range(nr):
        gnext = g(rs[j + 1])
        if float(gnext.upper()) < float(gprev.upper()):
            gnext = gprev                          # running max: a monotone upper bound (Stieltjes needs monotone g)
        gain += T(rs[j]) * (gnext - gprev)
        gprev = gnext
    gain += T(rs[-1]) * (1 - gprev)
    loss = CL(xa_lo, False) - CL(2 * t, True)
    return term_s + gain - loss


def analytic(t, HG):
    """Certified upper bound of sup_{0 < eps <= EPS1} H/eps^2 over a few (lam, eta); returns best (value, lam, eta)."""
    best = None
    ta = arb(t)
    s1 = (1 - arb(EPS1) ** 2).sqrt()
    x1 = float((2 * ta / s1).upper())
    m = min(LO[cell_of(2 * t):cell_of(x1) + 1])
    for lam in (0.5, 1.0, 1.5):
        a = (2 - arb(lam)) * ta
        h = arb(lam) * ta / arb(EPS1)
        if float(h.lower()) < 2:
            continue
        ms = arb(16.01) * h ** 4 * (-2 * h * h).exp() / (arb(lam) * ta) ** 2
        for eta in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2):
            B = 2 * a * arb(eta) * s1
            if float(B.lower()) < EPS1 ** 2:
                continue
            L = arb(Dz(2 * t - 2 * eta))
            hg = arb(Hg_up(eta, HG))
            phi = (ms + L * hg / (2 * a * s1) + L * L * arb(EPS1) ** 2 / (2 * a * a * s1 * s1)
                   + (-B / arb(EPS1) ** 2).exp() / arb(EPS1) ** 2 - ta * arb(m))
            v = float(phi.upper())
            if best is None or v < best[0]:
                best = (v, lam, eta)
    return best


def main():
    hg_path, tag = sys.argv[1], sys.argv[2]
    HG, hgd = load_hg(hg_path)
    t0 = time.time()
    out = {"tag": tag, "hgeom": os.path.basename(hg_path), "hgeom_root": hgd["root"], "hgeom_N": hgd["N"],
           "hgeom_delta": hgd["delta"], "EPS1": EPS1, "RATIO": RATIO, "REL": REL, "levels": []}
    for alpha, ts in LEVELS:
        t = float(ts)
        an = analytic(t, HG)
        boxes, ea, reach = [], EPS1, EPS1 if an and an[0] < 0 else 0.0
        while reach > 0 and ea < EPS_MAX:
            eb = ea * RATIO
            vals = [(box_H(t, ea, eb, lam, HG), lam) for lam in LAMS]
            vals = [(float(v.upper()), lam) for v, lam in vals if v is not None]
            if not vals:
                break
            v, lam = min(vals)
            if v >= 0:
                break
            boxes.append((ea, eb, lam, v))
            reach, ea = eb, eb
        out["levels"].append({"alpha": alpha, "t": ts, "analytic": an, "reach": reach, "nboxes": len(boxes),
                              "boxes_sample": boxes[:: max(1, len(boxes) // 12)]})
        print(f"alpha={alpha} t'={ts}: analytic (0,{EPS1}] sup H/eps^2 <= {an[0] if an else None} "
              f"(lam,eta)={an[1:] if an else None}; boxes {len(boxes)} -> certified (0, {reach:.4f}]", flush=True)
    out["sec"] = time.time() - t0
    path = os.path.join(ROOT_DIR, "deposits", "dt16", f"neprime_arb_{tag}.json")
    try:
        with open(path + ".tmp", "w") as fh:
            json.dump(out, fh)
        os.replace(path + ".tmp", path)
    except OSError as err:
        sys.exit(f"cannot write {path}: {err}")
    print(f"{out['sec']:.0f}s -> {path}")


if __name__ == "__main__":
    main()
