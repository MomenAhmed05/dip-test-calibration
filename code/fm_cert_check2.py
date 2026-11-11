"""Rigorous ball-arithmetic enclosure of F_M at the two cutoff quantiles (Step 10, gap G1; F-8 checker).

Same sandwich and Arb arithmetic as fm_cert_check.py (Balabdaoui & Pitman 2011 stick-breaking renewal),
but on a NON-uniform grid: fine cells near 2*c_.05 ~ 1.463 and 2*c_.01 ~ 1.72, coarse cells elsewhere.
The sandwich is valid for any increasing grid (each cell is bounded from both sides), so only the width
changes. Certifies, with outward rounding:
  F_M(1.4626) < 0.95 <= F_M(1.4630)   =>  q_.95(M) in (1.4626, 1.4630]  =>  c_.05 = q/2 in (0.7313, 0.7315]
  F_M(1.7196) < 0.99 <= F_M(1.7198)   =>  c_.01 in (0.8598, 0.8599]
Output: deposits/dt8/fm_cert_check2.json
"""
import json
import os
import sys

import numpy as np
from flint import arb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fm_cert_check as C  # noqa: E402  (FE_ball, constants, ctx.prec = 80)


def grid(mode):
    if mode == "wide":     # rigorous table on [1.0, 2.2] for quantile brackets and 1 - F_M(2 qbar)
        segs = [(1.0, 2.2, 0.001), (2.2, 3.0, 0.01), (3.0, C.YMAX, 0.05)]
    elif mode == "low":    # DT10 (Step 11): the wide table extended down to 0.5. FE_ball's remainder bound needs
        # consecutive-term ratio <= 1/2 for k >= 15; for z >= 0.3 it is <= 1.2 * exp(-2*31*0.09) < 0.005.
        segs = [(0.5, 1.0, 0.002), (1.0, 2.2, 0.001), (2.2, 3.0, 0.01), (3.0, C.YMAX, 0.05)]
    else:
        segs = [(1.40, 1.46, 0.002), (1.46, 1.466, 0.0002), (1.466, 1.716, 0.002), (1.716, 1.722, 0.0002),
                (1.722, 2.2, 0.002), (2.2, 3.0, 0.01), (3.0, C.YMAX, 0.05)]
    pts = []
    for a, b, h in segs:
        k = int(round((b - a) / h))
        pts += [(a, i, h) for i in range(k)]
    pts.append((C.YMAX, 0, 0.0))
    # exact decimal construction of each grid point as an arb (a + i*h)
    return [arb(f"{a:.6f}") + arb(f"{h:.6f}") * i for a, i, h in pts]


def solve(ys, kernel_power=1):
    N = len(ys) - 1
    lo = [0.0] * (N + 1)
    hi = [1.0] * (N + 1)
    lo[N] = C.lower(1 - 2 * (-(ys[N] * ys[N]) / 2).exp())
    nt = int(np.ceil(np.log(C.WMAX / float(ys[N].mid())) / np.log(C.TAILR)))
    wt = [ys[N] * arb(C.TAILR) ** k for k in range(nt + 1)]

    def fe(y, w):
        a = y * w / (w * w - y * y).sqrt()
        return C.FE_ball(a) ** kernel_power

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
        lo[i] = max(0.0, C.lower(Alo / (1 - c_lo)))
        hi[i] = min(1.0, C.upper(Ahi), hi[i + 1])
    return lo, hi


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "local"
    ys = grid(mode)
    lo, hi = solve(ys)
    mids = [float(y.mid()) for y in ys]
    if mode == "low":
        if mids[0] < 0.3:
            raise SystemExit("FE_ball remainder bound not justified below 0.3")
        out = os.path.join(os.path.dirname(__file__), "..", "deposits", "dt10", "fm_arb_low.npz")
        w, mm = np.array(hi) - np.array(lo), np.array(mids)
        try:
            os.makedirs(os.path.dirname(out), exist_ok=True)
            np.savez(out, ys=mm, lo=np.array(lo), hi=np.array(hi))
        except OSError as exc:
            raise SystemExit(f"cannot write {out}: {exc}") from exc
        print(json.dumps({"points": len(mids), "max_width_0.5_1.0": float(np.max(w[mm < 1.0])),
                          "max_width_1.0_2.2": float(np.max(w[(mm >= 1.0) & (mm <= 2.2)]))}))
        return
    if mode == "wide":
        np.savez(os.path.join(os.path.dirname(__file__), "..", "deposits", "dt8", "fm_arb_wide.npz"),
                 ys=np.array(mids), lo=np.array(lo), hi=np.array(hi))
        m = np.array(mids)
        rep = {}
        for p in (0.90, 0.95, 0.99):
            up = m[np.argmax(np.array(lo) >= p)]                 # F_M(up) >= p  =>  q_p <= up
            dn = m[np.nonzero(np.array(hi) < p)[0][-1]]          # F_M(dn) < p   =>  q_p > dn
            rep[str(p)] = {"q_bracket": [float(dn), float(up)], "c_bracket": [float(dn) / 2, float(up) / 2]}
        rep["max_width_1.0_2.2"] = float(np.max((np.array(hi) - np.array(lo))[m <= 2.2]))
        print(json.dumps(rep, indent=1))
        out = os.path.join(os.path.dirname(__file__), "..", "deposits", "dt8", "fm_arb_wide.json")
        try:
            with open(out, "w") as f:
                json.dump(rep, f, indent=1)
        except OSError as exc:
            raise SystemExit(f"cannot write {out}: {exc}") from exc
        return

    def at(x):
        i = int(np.argmin(np.abs(np.array(mids) - x)))
        if abs(mids[i] - x) > 1e-9:
            raise SystemExit(f"{x} is not a grid point")
        return lo[i], hi[i]

    pts = {x: at(x) for x in (1.4626, 1.4628, 1.4630, 1.7196, 1.7198)}
    c05 = pts[1.4626][1] < 0.95 <= pts[1.4630][0]
    c01 = pts[1.7196][1] < 0.99 <= pts[1.7198][0]
    res = {"grid_points": len(ys), "enclosures": {f"{x:.4f}": [a, b] for x, (a, b) in pts.items()},
           "certified_c05_in_(0.7313,0.7315]": bool(c05), "certified_c01_in_(0.8598,0.8599]": bool(c01)}
    out = os.path.join(os.path.dirname(__file__), "..", "deposits", "dt8", "fm_cert_check2.json")
    try:
        with open(out, "w") as f:
            json.dump(res, f, indent=1)
    except OSError as exc:
        raise SystemExit(f"cannot write {out}: {exc}") from exc
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
