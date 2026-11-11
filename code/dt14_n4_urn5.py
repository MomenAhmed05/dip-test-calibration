"""Step 15 / DT14 WP-B (O1-a at n = 4): certified upper bound on the urn value M_5 = sup_z Ubar_5(z).

Urn bound: S_4 = sup_G Phi_r(G) <= M_5, Ubar_5(z) = mean of K_r over the five 4-subsets of z.
Domain: z = (+1, s2 a2, .., s5 a5), 1 >= a2 >= .. >= a5 > 0, s in {+1, -1}^4 (scale, global reflection, permutation).
Box: magnitude intervals [l_j, h_j] (rational; l_j = 0 = "includes deep values"); j0 = first index with l = 0.
Per 4-subset I (bound of K over the box, each valid for every z in the box):
  - no deep member: K(s.h_I) * prod_{i in I} h_i / l_i                       (monotone volumes)
  - deep members D = I cap {j0..}, retained R = I minus D, |R| >= 1:
      Kpin(s.h_R, zeros) * prod_{i in R} h_i / l_i + E(|D|, eta), eta = h_{j0} / min_{i in R} l_i   (collapse)
  - R empty (I = {2,3,4,5} all deep): bound 1.
Ubar bound = mean of the five subset bounds. Prune if <= TAU. Split: if a deep collapse is too coarse (eta > ETA_MAX), split
coordinate j0 at ~ETA_MAX * l_{j0-1}; else split the retained coordinate with the largest log-width.
Kernel: exact rational (CERT) or --float (DEV pilot only). Usage: python dt14_n4_urn5.py r_num r_den tau_num tau_den
  [--float] [--signs i0-i1] [--log out.json]
"""
import json
import sys
import time
from fractions import Fraction
from itertools import combinations, product
from math import log

sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
from dt14_n4_supk import Engine, collapse_err, geo_mid  # noqa: E402

ETA_MAX = Fraction(1, 400)
MIN_LOGW = 1e-4
NF = 4                      # free coordinates a2..a5
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


def ubar_bound(eng, s, box, r):
    """Returns (bound, split_hint). Indices 0..4 = z1..z5; free coordinate j <-> index j + 1."""
    lo = [Fraction(1)] + [b[0] for b in box]
    hi = [Fraction(1)] + [b[1] for b in box]
    sg = [1] + list(s)
    deep0 = next((i for i in range(5) if lo[i] == 0), None)        # first deep index (>= 1)
    total = Fraction(0)
    need_eta_split = False
    for sub in SUBSETS:
        D = [i for i in sub if deep0 is not None and i >= deep0]
        R = [i for i in sub if i not in D]
        if not R:
            total += 1
            continue
        fac = Fraction(1)
        for i in R:
            fac *= hi[i] / lo[i]
        if not D or deep0 is None:
            total += eng.kernel([sg[i] * hi[i] for i in sub]) * fac
            continue
        rho = min(lo[i] for i in R)
        eta = hi[deep0] / rho
        if eta > ETA_MAX:
            need_eta_split = True
        same = len(R) == 2 and sg[R[0]] == sg[R[1]]
        kp = eng.kernel([sg[i] * hi[i] for i in R] + [Fraction(0)] * len(D))
        total += kp * fac + collapse_err(len(D), eta, r, same)
    return total / 5, (deep0 if need_eta_split else None)


def run(r, tau, exact, sign_range):
    eng = Engine(r, exact)
    st = {"boxes": 0, "pruned": 0, "failed": [], "max_bound_pruned": 0.0}
    t0 = time.time()
    patterns = list(product((1, -1), repeat=NF))[sign_range[0]:sign_range[1]]
    for s in patterns:
        stack = [[(Fraction(0), Fraction(1))] * NF]
        while stack:
            box = tighten(stack.pop())
            if box is None:
                continue
            st["boxes"] += 1
            if st["boxes"] % 5000 == 0:
                print(f"  .. {st['boxes']} boxes, stack {len(stack)}, calls {eng.calls} {time.time() - t0:.0f}s", flush=True)
            val, deep0 = ubar_bound(eng, s, box, r)
            if val <= tau:
                st["pruned"] += 1
                st["max_bound_pruned"] = max(st["max_bound_pruned"], float(val))
                continue
            j, pt = None, None
            if deep0 is not None:                       # deep collapse too coarse: cut deep coordinate
                jf = deep0 - 1                          # free-coordinate index
                rho = Fraction(1) if jf == 0 else box[jf - 1][0]
                if rho > 0:
                    pt = (ETA_MAX * rho).limit_denominator(10**12)
                    if 0 < pt < box[jf][1]:
                        j = jf
                if j is None:
                    pt = None
            if j is None:
                widths = [log(float(b[1] / b[0])) if b[0] > 0 else -1.0 for b in box]
                j = max(range(NF), key=lambda i: widths[i])
                if widths[j] < MIN_LOGW:
                    if widths[j] < 0:          # only deep coordinates left and still above tau: cut the first deep one
                        j = next(i for i in range(NF) if box[i][0] == 0)
                        pt = box[j][1] / 2
                    else:
                        st["failed"].append({"s": s, "box": [[str(a), str(b)] for a, b in box], "bound": float(val)})
                        continue
                else:
                    pt = geo_mid(*box[j])
            if pt is None:
                sys.exit("internal: no split point")
            a, b = box[j]
            left, right = list(box), list(box)
            left[j], right[j] = (a, pt), (pt, b)
            stack.extend([left, right])
        print(f"signs {s}: boxes {st['boxes']} pruned {st['pruned']} failed {len(st['failed'])} calls {eng.calls} "
              f"{time.time() - t0:.0f}s", flush=True)
    st.update({"r": str(r), "tau": str(tau), "exact": exact, "kernel_calls": eng.calls, "seconds": time.time() - t0,
               "ETA_MAX": str(ETA_MAX), "sign_range": sign_range, "certified": exact and not st["failed"]})
    return st


def main():
    r = Fraction(int(sys.argv[1]), int(sys.argv[2]))
    tau = Fraction(int(sys.argv[3]), int(sys.argv[4]))
    exact = "--float" not in sys.argv
    sr = (0, 2 ** NF)
    if "--signs" in sys.argv:
        a, b = sys.argv[sys.argv.index("--signs") + 1].split("-")
        sr = (int(a), int(b))
    out = sys.argv[sys.argv.index("--log") + 1] if "--log" in sys.argv else None
    st = run(r, tau, exact, sr)
    print(json.dumps({k: v for k, v in st.items() if k != "failed"} | {"n_failed": len(st["failed"])}), flush=True)
    if out:
        try:
            with open(out, "w") as fh:
                json.dump(st, fh, default=str)
        except OSError as err:
            sys.exit(f"cannot write {out}: {err}")


if __name__ == "__main__":
    main()
