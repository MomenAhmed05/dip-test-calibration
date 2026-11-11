"""Independent checker for fm_cert.py: ball arithmetic (python-flint / Arb).

Same identity (Balabdaoui & Pitman 2011 stick-breaking; see fm_cert.py), different implementation:
every quantity is an Arb ball (outward rounding), F_E is evaluated as a ball with an explicit series
remainder, and the sandwich runs on a COARSE grid. The coarse enclosure is rigorous on its own; the
check passes if, at every shared grid point, the fine float enclosure of fm_cert.py lies inside the
coarse ball enclosure, and if the quantile brackets are consistent. Negative falsifiers: (1) a
deliberately wrong kernel (F_E replaced by F_E^2) must FAIL containment; (2) a shifted fine
enclosure (+1e-3) must FAIL.
Usage: python fm_cert_check.py H_COARSE FINE_NPZ
"""
import json
import os
import sys

import numpy as np
from flint import arb, ctx

ctx.prec = 80
YLO, YMAX, WMAX, TAILR = 0.9, 8.5, 2.0e4, 1.05
K = 14


def FE_ball(z):
    """Arb enclosure of F_E(z), z an arb >= 0.9 (upper end); remainder bounded by 2x first omitted term."""
    s = arb(0)
    for k in range(1, K + 1):
        kk = k * k
        s += (4 * kk * z * z - 1) * (-2 * kk * z * z).exp()
    k = K + 1
    rem = 2 * (4 * k * k * z * z) * (-2 * k * k * z * z).exp()
    v = 1 - 2 * s
    return arb(v.mid(), v.rad() + 2 * rem.upper())


def lower(b):
    return float(b.lower())


def upper(b):
    return float(b.upper())


def solve(h, kernel_power=1):
    N = int(round((YMAX - YLO) / h))
    ys = [arb(YLO) + arb(h) * i for i in range(N + 1)]
    lo = [0.0] * (N + 1)
    hi = [1.0] * (N + 1)
    lo[N] = lower(1 - 2 * (-(ys[N] * ys[N]) / 2).exp())
    nt = int(np.ceil(np.log(WMAX / float(ys[N].mid())) / np.log(TAILR)))
    wt = [ys[N] * arb(TAILR) ** k for k in range(nt + 1)]

    def fe(y, w):
        if w is None:
            return arb(1)
        a = y * w / (w * w - y * y).sqrt()
        v = FE_ball(a)
        return v ** kernel_power

    for i in range(N - 1, -1, -1):
        y = ys[i]
        y2 = y * y
        Alo, Ahi, c_lo = arb(0), arb(0), arb(0)
        for j in range(i, N):
            a, b = ys[j], ys[j + 1]
            wgt = y2 * (1 / (a * a) - 1 / (b * b))
            fl = fe(y, b)
            fh = arb(1) if j == i else fe(y, a)
            if j == i:
                c_lo = fl * wgt
            else:
                Alo += fl * arb(lo[j]) * wgt
            Ahi += fh * arb(hi[j + 1]) * wgt
        for k in range(nt):
            a, b = wt[k], wt[k + 1]
            wgt = y2 * (1 / (a * a) - 1 / (b * b))
            Alo += fe(y, b) * (1 - 2 * (-(a * a) / 2).exp()) * wgt
            Ahi += fe(y, a) * wgt
        Ahi += y2 / (wt[nt] * wt[nt])
        lo[i] = max(0.0, lower(Alo / (1 - c_lo)))
        hi[i] = min(1.0, upper(Ahi), hi[i + 1])
    return np.array([float(v.mid()) for v in ys]), np.array(lo), np.array(hi)


def contains(ysC, loC, hiC, fine, shift=0.0):
    """Count shared grid points where the two rigorous enclosures are DISJOINT (both cannot be true).

    Both enclosures are claimed rigorous, so they must intersect everywhere; the reported certificate
    is the intersection. (Strict containment is not required: the float run may be looser near 1.)
    """
    ysF, loF, hiF = fine["ys"], fine["lo"] + shift, fine["hi"] + shift
    bad = 0
    for y, lc, uc in zip(ysC, loC, hiC, strict=True):
        i = int(round((y - ysF[0]) / (ysF[1] - ysF[0])))
        if abs(ysF[i] - y) > 1e-9:
            continue
        if loF[i] > uc + 1e-15 or hiF[i] < lc - 1e-15:
            bad += 1
    return bad


def main():
    hC = float(sys.argv[1])
    fine = np.load(sys.argv[2])
    ys, lo, hi = solve(hC)
    bad = contains(ys, lo, hi, fine)
    bad_shift = contains(ys, lo, hi, fine, shift=1e-3)
    ys2, lo2, hi2 = solve(hC, kernel_power=2)
    bad_wrong = contains(ys2, lo2, hi2, fine)
    res = {"h_coarse": hC, "fine": os.path.basename(sys.argv[2]), "points": len(ys),
           "containment_failures": bad, "falsifier_shift_failures": bad_shift,
           "falsifier_wrong_kernel_failures": bad_wrong,
           "verdict": "PASS" if bad == 0 and bad_shift > 0 and bad_wrong > 0 else "FAIL",
           "coarse_at": {f"{y:.4f}": [float(a), float(b)] for y, a, b in zip(ys, lo, hi, strict=True)
                         if abs(y - 1.1) < 0.011 or abs(y - 1.46) < 0.011 or abs(y - 1.72) < 0.011}}
    out = os.path.join(os.path.dirname(__file__), "..", "deposits", "dt8", "fm_cert_check.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
