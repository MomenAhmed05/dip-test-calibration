"""DT11 (Step 12) independent interval check of the near-edge box constants.

Re-computes, with Arb balls (python-flint, 200 bits), every scalar step of DT11-NEAREDGE.md §5
and of lean/NearEdge.lean `scalar_box`, on the box t in [0.72, 0.88], eps in (0, 0.15]:
  - a = sqrt(1 - eps^2) >= 0.9886..., x = 2t/a <= 1.7802, 2x^2 <= 6.34, e^{-2x^2} >= 0.0016
  - J-bar(1.44) <= 4/(e^{4.1472} - 1) <= 0.065;  r = (16/3) e^{-6}/(1 - e^{-5}) <= 0.014
  - loss coefficient 4(8t^2 - 1/2)(1 - 0.065) >= 13.6 at t = 0.72
  - gain coefficient 4(1.014)(8 - a^2/(2t^2))/(2.8a^2 - 8eps^2) <= 11.7 (worst corner)
  - range tail 2S(0.6t/eps)/eps^2 <= 0.001, truncation tail e^{-2R(R+kappa)}/eps^2 <= 0.0001
  - final margin 13.6*0.0016 - 11.7*0.0016 - 0.0011 > 0
Also a dense float scan of the exact (un-coarsened) bound H/eps^2 over the box.
Deterministic: no random numbers. Usage: python -I dt11_scalar_check.py  -> deposits/dt11/scalar_check.json
"""
import json
import math
import os
import sys

from flint import arb, ctx

ctx.prec = 200
OUT = os.path.join(os.path.dirname(__file__), "..", "deposits", "dt11", "scalar_check.json")


def S_arb(z, K=12):
    return 2 * sum((4 * k * k * z * z - 1) * (-2 * k * k * z * z).exp() for k in range(1, K))


def check(name, ok, value, rec):
    rec[name] = {"pass": bool(ok), "value": str(value)}
    return ok


def main():
    rec = {}
    t1, t2, e1 = arb("0.72"), arb("0.88"), arb("0.15")
    a2 = 1 - e1 * e1
    a = a2.sqrt()
    xmax = 2 * t2 / a
    ok = True
    ok &= check("a_min >= 0.9886", a >= arb("0.9886"), a, rec)
    ok &= check("x_max <= 1.8", xmax <= arb("1.8"), xmax, rec)
    ok &= check("x_max^2 <= 3.17", xmax * xmax <= arb("3.17"), xmax * xmax, rec)
    e2 = (-arb("6.34")).exp()
    ok &= check("exp(-6.34) >= 0.0016", e2 >= arb("0.0016"), e2, rec)
    jb = 4 / (arb("4.1472").exp() - 1)
    ok &= check("Jbar(1.44) <= 0.065", jb <= arb("0.065"), jb, rec)
    r = arb(16) / 3 * arb(-6).exp() / (1 - arb(-5).exp())
    ok &= check("r <= 0.014", r <= arb("0.014"), r, rec)
    loss = 4 * (8 * t1 * t1 - arb("0.5")) * (1 - arb("0.065"))
    ok &= check("loss coef >= 13.6", loss >= arb("13.6"), loss, rec)
    gain = 4 * arb("1.014") * (8 - a2 / (2 * t2 * t2)) / (arb("2.8") * a2 - 8 * e1 * e1)
    ok &= check("gain coef <= 11.7", gain <= arb("11.7"), gain, rec)
    s0 = (arb("0.6") * t1 / e1) ** 2
    tail1 = 16 * arb("1.014") * s0 * s0 * (-2 * s0).exp() / (arb("0.36") * t1 * t1)
    ok &= check("range tail/eps^2 <= 0.001", tail1 <= arb("0.001"), tail1, rec)
    sig = 1 / (e1 * e1)
    tail2 = sig * (-arb("1.27") * sig).exp()
    ok &= check("trunc tail/eps^2 <= 0.0001", tail2 <= arb("0.0001"), tail2, rec)
    margin = (arb("13.6") - arb("11.7")) * arb("0.0016") - arb("0.0011")
    ok &= check("final margin > 0", margin > 0, margin, rec)
    sharp = (loss - gain) * (-2 * xmax * xmax).exp() - tail1 - tail2
    rec["sharp margin (B0*e^{-2x^2} - tails)"] = str(sharp)

    # dense float scan of the exact bound H/eps^2 (DT11 §5 before coarsening)
    def S(z):
        return 2 * sum((4 * k * k * z * z - 1) * math.exp(-2 * k * k * z * z) for k in range(1, 10))

    def Jb(z):
        return 4 * sum(math.exp(-2 * k * k * z * z) for k in range(1, 10))

    worst = -1e9
    for i in range(41):
        t = 0.72 + 0.16 * i / 40
        for j in range(1, 601):
            e = 0.15 * j / 600
            aa = math.sqrt(1 - e * e)
            x = 2 * t / aa
            loss_f = 4 * (8 * t * t - 0.5) * (1 - Jb(2 * t)) * math.exp(-2 * x * x)
            gain_f = 4 * 1.014 * (8 - aa * aa / (2 * t * t)) / (2.8 * aa * aa - 8 * e * e) * math.exp(-2 * x * x)
            R = (2 * t - aa) / e
            kap = 1.4 * t / e
            h = -loss_f + gain_f + (2 * S(0.6 * t / e) + math.exp(-2 * R * (R + kap))) / (e * e)
            worst = max(worst, h)
    rec["float scan max H/eps^2 (coarse bound, 41x600 grid)"] = worst
    ok &= check("float scan negative", worst < 0, worst, rec)
    rec["ALL_PASS"] = bool(ok)
    try:
        os.makedirs(os.path.dirname(OUT), exist_ok=True)
        with open(OUT, "w") as fh:
            json.dump(rec, fh, indent=1)
    except OSError as err:
        sys.exit(f"cannot write {OUT}: {err}")
    print(json.dumps({k: (v["pass"] if isinstance(v, dict) else v) for k, v in rec.items()}, indent=1))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
