"""Step 14 (DT13) Part Q analysis: certified-MC bracket of qbar_alpha (the (1 - alpha) quantile of D(B), the limit
of sqrt(n) D_n under the uniform) for alpha in {0.10, 0.05, 0.01}, and the asymptotic worst-case size of the
uniform-calibrated test, 1 - F_M(2 qbar) (C4 / DT3 Thm 1', Thm 4), bracketed with the certified F_M enclosure.
Same rule as dt8_cert_analysis.part_a (order statistics hi_(k), lo_(l)); eps per side per alpha = DELTA_Q / 6.
F_M enclosure: deposits/dt8/fm_cert_h0.0002.npz on [0.9, 8.5]; deposits/dt10/fm_arb_low.npz below 0.9.
Usage: python dt13_qbar_analysis.py root tag  (reads deposits/qbar/{tag}_s*.json)
"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt8_cert_analysis import _largest, _smallest, binom_cdf  # noqa: E402

DEP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits")
DELTA_Q = 0.002


def fm_bracket(y_lo, y_hi):
    """Certified [lower, upper] of F_M on [y_lo, y_hi] endpoints, outward grid points (F_M nondecreasing)."""
    out = []
    for f in (os.path.join(DEP, "dt8", "fm_cert_h0.0002.npz"), os.path.join(DEP, "dt10", "fm_arb_low.npz")):
        z = np.load(f)
        ys, lo, hi = z["ys"], z["lo"], z["hi"]
        if ys[0] <= y_lo and y_hi <= ys[-1]:
            i_lo = int(np.searchsorted(ys, y_lo, side="right")) - 1     # ys[i_lo] <= y_lo
            i_hi = int(np.searchsorted(ys, y_hi, side="left"))          # ys[i_hi] >= y_hi
            out.append((float(lo[i_lo]), float(hi[i_hi]), os.path.basename(f)))
    if not out:
        raise SystemExit(f"no F_M enclosure covers [{y_lo}, {y_hi}]")
    return min(out, key=lambda t: t[1] - t[0])


def main():
    root, tag = int(sys.argv[1]), sys.argv[2]
    lo, hi = [], []
    for f in sorted(glob.glob(os.path.join(DEP, "qbar", f"{tag}_s*.json"))):
        try:
            with open(f) as fh:
                r = json.load(fh)
        except (OSError, ValueError) as err:
            raise SystemExit(f"cannot read {f}: {err}") from err
        if r["root"] != root:
            raise SystemExit(f"root mismatch in {f}")
        lo += r["lo"]
        hi += r["hi"]
    lo, hi = np.sort(lo), np.sort(hi)
    R = len(lo)
    eps = DELTA_Q / 6
    res = {"root": root, "tag": tag, "R": R, "eps_per_side": eps, "alphas": {}}
    for a in (0.10, 0.05, 0.01):
        p = 1 - a
        k = _smallest(lambda kk, p=p: 1 - binom_cdf(kk - 1, R, p) <= eps, R)
        ell = _largest(lambda ll, p=p: binom_cdf(ll - 1, R, p) <= eps, R)
        if k is None or ell is None:
            print(f"alpha={a}: R too small")
            continue
        q_lo, q_up = float(lo[ell - 1]), float(hi[k - 1])
        f_lo, f_hi, src = fm_bracket(2 * q_lo, 2 * q_up)
        res["alphas"][str(a)] = {"qbar_bracket": [q_lo, q_up], "worst_size_bracket": [1 - f_hi, 1 - f_lo],
                                 "fm_source": src, "k": k, "l": ell}
        print(f"alpha={a}: qbar in [{q_lo:.5f}, {q_up:.5f}]; asymptotic worst-case size in [{1 - f_hi:.4f}, {1 - f_lo:.4f}] ({src})")
    try:
        with open(os.path.join(DEP, "qbar", f"{tag}_dt13_analysis.json"), "w") as fh:
            json.dump(res, fh, indent=1)
    except OSError as err:
        raise SystemExit(f"cannot write: {err}") from err


if __name__ == "__main__":
    main()
