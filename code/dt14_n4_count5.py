"""Step 15 / DT14 WP-B (O1-a at n = 4, alpha = 10%): certify the COUNTING LEMMA for 5 radii.

Claim to certify: for every z = (+1, s2 a2, .., s5 a5) (1 >= a2 >= .. >= a5 > 0, s in {+1, -1}^4), at least two of the five
4-subsets I have K_r(z_I) <= THETA. With sup K_r <= KAPPA (dt14_n4_supk.py) and the urn bound (DT14-N4UB §1):
  S_4 <= M_5 <= (3 KAPPA + 2 THETA) / 5.
Per box and subset: upper bound = min over (B1) monotone volumes and (B2) collapse of every deep suffix of the subset's
sorted coordinates (S2 sign-aware constants, DT14-N4UB §4b). Prune the box if >= 2 subset bounds are <= THETA.
Kernel: exact rational (fmpq, CERT) or --float (DEV pilot).
Usage: python dt14_n4_count5.py r_num r_den theta_num theta_den [--float] [--signs i0-i1] [--log out.json]
"""
import json
import sys
import time
from fractions import Fraction
from itertools import combinations, product
from math import log

sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
from dt14_n4_supk import Engine, collapse_err_signed, geo_mid  # noqa: E402

MIN_LOGW = 1e-4
DEEP_ETA = 1e-3          # split heuristic only (does not affect validity)
NF = 4
SUBSETS = list(combinations(range(5), 4))


def tighten(box):
    lo = [b[0] for b in box]
    hi = [b[1] for b in box]
    for j in range(1, NF):
        hi[j] = min(hi[j], hi[j - 1])
    for j in range(NF - 2, -1, -1):
        lo[j] = max(lo[j], lo[j + 1])
    if any(lo[j] > hi[j] or hi[j] <= 0 for j in range(NF)):
        return None
    return list(zip(lo, hi, strict=True))


def subset_bound(eng, lo, hi, sg):
    """Upper bound on K over the box for one 4-tuple; lo/hi/sg in DESCENDING magnitude order (lo[0] may be 1 = z1)."""
    r = eng.r
    best = Fraction(1)
    if all(v > 0 for v in lo):
        fac = Fraction(1)
        for i in range(4):
            fac *= hi[i] / lo[i]
        best = min(best, eng.kernel([sg[i] * hi[i] for i in range(4)]) * fac)
    for j0 in range(1, 4):                    # deep suffix j0..3, retained 0..j0-1 (at least one retained)
        rho = min(lo[:j0])
        if rho == 0:
            continue
        eta = hi[j0] / rho
        if eta >= 1:
            continue
        fac = Fraction(1)
        for i in range(j0):
            fac *= hi[i] / lo[i]
        kp = eng.kernel([sg[i] * hi[i] for i in range(j0)] + [Fraction(0)] * (4 - j0))
        val = kp * fac + collapse_err_signed(eta, r, sg[:j0], sg[j0:])
        best = min(best, val)
    return best


def box_bounds(eng, s, box):
    lo = [Fraction(1)] + [b[0] for b in box]
    hi = [Fraction(1)] + [b[1] for b in box]
    sg = [1] + list(s)
    return [subset_bound(eng, [lo[i] for i in sub], [hi[i] for i in sub], [sg[i] for i in sub]) for sub in SUBSETS]


MODE = "count"      # "count": >= 2 subset bounds <= THETA; "mean": mean of the 5 subset bounds <= THETA (= M_5 <= THETA)


def resolved(b, theta):
    if MODE == "mean":
        return sum(b) / 5 <= theta
    return sum(1 for v in b if v <= theta) >= 2


def run(r, theta, exact, sign_range):
    eng = Engine(r, exact)
    st = {"boxes": 0, "pruned": 0, "failed": []}
    t0 = time.time()
    for s in list(product((1, -1), repeat=NF))[sign_range[0]:sign_range[1]]:
        stack = [[(Fraction(0), Fraction(1))] * NF]
        while stack:
            box = tighten(stack.pop())
            if box is None:
                continue
            st["boxes"] += 1
            if st["boxes"] % 5000 == 0:
                print(f"  .. {st['boxes']} boxes, stack {len(stack)}, calls {eng.calls}, box "
                      f"{[(float(a), float(b)) for a, b in box]} {time.time() - t0:.0f}s", flush=True)
            b = box_bounds(eng, s, box)
            if resolved(b, theta):
                st["pruned"] += 1
                continue
            # effective log-width; a deep interval (0, h] counts as [DEEP_ETA * rho, h] (rho = lower end of the previous
            # coordinate, 1 for a2): below that the collapse bounds are already tight, so deepening it further is useless
            lw = []
            for j, x in enumerate(box):
                if x[0] > 0:
                    lw.append(log(float(x[1] / x[0])))
                elif x[1] > 0:
                    rho = Fraction(1) if j == 0 else box[j - 1][0]
                    lw.append(log(float(x[1]) / (DEEP_ETA * float(rho))) if rho > 0 else 1e9)
                else:
                    lw.append(-1.0)
            cand = [j for j in range(NF) if lw[j] >= MIN_LOGW]
            if not cand:
                st["failed"].append({"s": s, "box": [[str(a), str(c)] for a, c in box],
                                     "bounds": [float(v) for v in b]})
                continue
            j = max(cand, key=lambda i: lw[i])
            pt = geo_mid(*box[j])
            a, c = box[j]
            left, right = list(box), list(box)
            left[j], right[j] = (a, pt), (pt, c)
            stack.extend([left, right])
        print(f"signs {s}: boxes {st['boxes']} pruned {st['pruned']} failed {len(st['failed'])} calls {eng.calls} "
              f"{time.time() - t0:.0f}s", flush=True)
    st.update({"mode": MODE, "r": str(r), "theta": str(theta), "exact": exact, "kernel_calls": eng.calls,
               "seconds": time.time() - t0, "sign_range": sign_range, "certified": exact and not st["failed"]})
    return st


def main():
    r = Fraction(int(sys.argv[1]), int(sys.argv[2]))
    theta = Fraction(int(sys.argv[3]), int(sys.argv[4]))
    exact = "--float" not in sys.argv
    sr = (0, 2 ** NF)
    if "--signs" in sys.argv:
        a, b = sys.argv[sys.argv.index("--signs") + 1].split("-")
        sr = (int(a), int(b))
    out = sys.argv[sys.argv.index("--log") + 1] if "--log" in sys.argv else None
    global MODE
    if "--mean" in sys.argv:
        MODE = "mean"
    st = run(r, theta, exact, sr)
    print(json.dumps({k: v for k, v in st.items() if k != "failed"} | {"n_failed": len(st["failed"])}), flush=True)
    if out:
        try:
            with open(out, "w") as fh:
                json.dump(st, fh, default=str)
        except OSError as err:
            sys.exit(f"cannot write {out}: {err}")


if __name__ == "__main__":
    main()
