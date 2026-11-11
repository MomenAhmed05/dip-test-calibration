"""Step 15 / DT14 WP-B (O1-a at n = 4): certified upper bound on sup_z K_r(z) (the 4-point Khintchine kernel).

Why: S_4 = sup_G Phi_r(G) <= M_4 = sup_z K_r(z) (urn bound with k = 4). At alpha = 5% the DEV value
M_4 ~ 0.1752 is below the width-10 target 0.18596, so a certificate sup K <= KAPPA closes the n = 4, 5% width.

Domain (scale + global reflection + permutation invariance): z = (+1, s2 a2, s3 a3, s4 a4), 1 >= a2 >= a3 >= a4 > 0,
signs s in {+1, -1}^3. Boxes [l_j, h_j] in magnitude (rational endpoints; l_j = 0 means "includes deep values").
Bounds (each valid for every z in the box with the sorting constraint):
  (B1) all l > 0 (monotone volumes): K(z) <= K(s.h) * prod_j h_j / l_j.
  (B2) first index j0 with l = 0: deep coordinates j0..4 have |x| <= eps = h_{j0}; retained radii >= rho = l_{j0-1}
       (rho = 1 if j0 = 2); eta = eps / rho. If eta <= ETA_MAX: K(z) <= Kpin(s.h_ret) * prod_ret h/l + E(m, eta),
       Kpin = kernel with the m = 5 - j0 deep points at 0 (exact: order sum / m!); E = collapse constants (
       to be S-verified; see DT14-N4UB.md). Otherwise split coordinate j0 at a rational point ~ ETA_MAX * rho.
Prune if bound <= KAPPA; else split the coordinate with the largest log-width (deep coordinates first if unresolved).
Kernel: --exact (rational, CERT) or --float (DEV pilot only).
Usage: python dt14_n4_supk.py r_num r_den kappa_num kappa_den [--float] [--log out.json]
"""
import json
import sys
import time
from fractions import Fraction
from itertools import product
from math import factorial, log

sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
import dt14_n4_kernel as KN  # noqa: E402
import exact_poly_fq as EQ  # noqa: E402
from flint import fmpq  # noqa: E402

ETA_MAX = Fraction(1, 4000)
MIN_LOGW = 1e-4          # stop splitting below this log-width (report as FAILED box)
CONSTANTS = "S2"         # "S2" = sign-aware constants (DT14-N4UB §4b, CERT); "S" = crude two-sided (§4);
                         # "P062" = reviewer's (DEV only)


def collapse_err_signed(eta, r, ret_signs, deep_signs):
    """DT14-N4UB §4b (sign-aware; the deep points have FIXED signs, so each shift is one-sided)."""
    m = len(deep_signs)
    one = Fraction(1)
    if m == 1:
        nu = sum(1 for s in ret_signs if s < 0)
        lin = 3 * r / (1 + r) if nu in (0, 3) else (1 + r) / r + 2 / (1 + r)
        return min(one, (3 + lin) * eta)
    if m == 2:
        if ret_signs[0] == ret_signs[1]:
            rel = [d * ret_signs[0] for d in deep_signs]          # deep signs relative to the retained sign
            beta = (1 + r) if rel == [1, 1] else (r if rel == [-1, -1] else 1 + 2 * r)
            B = (1 + 1 / r) * beta
            return min(one, 2 * eta / (1 + r) + (2 * (1 + 2 / r) + B * B) * eta * eta)
        return min(one, (1 + 2 / r) ** 2 * eta * eta)
    if m == 3:
        return min(one, (1 + 2 / r) * eta)
    raise ValueError(m)


def collapse_err(m, eta, r, same_sign_retained):
    """Collapse error for moving m deep points (|x| <= eps) to 0; eta = eps / (min retained radius)."""
    if CONSTANTS in ("S", "S2"):
        if m == 1:
            return min(Fraction(1), (3 + 12 * (3 + r + 1 / r)) * eta)
        if m == 2:
            if same_sign_retained:
                B = (1 + 1 / r) * (1 + 2 * r)
                return min(Fraction(1), (2 + 4 / (1 + r)) * eta + B * B * eta * eta)
            return min(Fraction(1), (1 + 2 / r) ** 2 * eta * eta)
        if m == 3:
            return min(Fraction(1), (1 + 2 / r) * eta)
        raise ValueError(m)
    if m == 1:
        return min(Fraction(1), (3 + 6 * r / (r + 1)) * eta)
    if m == 2:
        if same_sign_retained:
            B = (1 + 1 / r) * (1 + 2 * r)
            return min(Fraction(1), (2 + 2 / (r + 1)) * eta + B * B * eta * eta)
        return min(Fraction(1), (1 + 2 / r) ** 2 * eta * eta)
    if m == 3:
        return min(Fraction(1), (1 + 2 / r) * eta)
    raise ValueError(m)


class Engine:
    def __init__(self, r, exact):
        self.r, self.exact, self.cache, self.calls = r, exact, {}, 0

    def kernel(self, z):
        """K or pinned K (zeros allowed). Cached on the exact tuple."""
        key = tuple(z)
        if key in self.cache:
            return self.cache[key]
        self.calls += 1
        m = sum(1 for v in z if v == 0)
        if self.exact:
            # Only sign-consistent orders have positive volume (negatives < zeros < positives); the m! orders inside
            # the zero block give identical polytopes, so one canonical order per zero block replaces the /m!.
            zq = [fmpq(v.numerator, v.denominator) for v in z]
            rq = fmpq(self.r.numerator, self.r.denominator)
            sgn = [(v > 0) - (v < 0) for v in z]
            tot = fmpq(0)
            for o in KN.ORDERS:
                ks = [sgn[i] for i in o]
                if any(ks[k] > ks[k + 1] for k in range(3)):
                    continue
                zs = [i for i in o if sgn[i] == 0]
                if zs != sorted(zs):
                    continue
                A, b = KN._constraints(zq, rq, o)
                tot += EQ.volume(A, b)
            val = Fraction(int(tot.p), int(tot.q))
        else:
            if m == 0:
                val = Fraction(KN.kernel_float([float(v) for v in z], float(self.r)))
            else:
                tot = 0.0
                for o in KN.ORDERS:
                    A, b = KN._constraints([float(v) for v in z], float(self.r), o)
                    tot += KN._vol_float(A, b)
                val = Fraction(tot / factorial(m))
        self.cache[key] = val
        return val


def tighten(box):
    """Sorting 1 >= a2 >= a3 >= a4: h_j <= h_{j-1}, l_j >= l_{j+1}. Returns None if empty."""
    lo = [b[0] for b in box]
    hi = [b[1] for b in box]
    for j in range(1, 3):
        hi[j] = min(hi[j], hi[j - 1])
    for j in range(1, -1, -1):
        lo[j] = max(lo[j], lo[j + 1])
    if any(lo[j] > hi[j] or hi[j] <= 0 for j in range(3)):
        return None
    return list(zip(lo, hi, strict=True))


def bound(eng, s, box):
    """Upper bound on K over the box: returns (value, main, err, j0) with value = main + err; j0 = first deep free
    coordinate (None if no deep coordinate). (B1) if no deep coordinate; (B2) collapse otherwise."""
    r = eng.r
    lo = [b[0] for b in box]
    hi = [b[1] for b in box]
    best = None
    if all(v > 0 for v in lo):
        fac = Fraction(1)
        for j in range(3):
            fac *= hi[j] / lo[j]
        main = eng.kernel([Fraction(1)] + [s[j] * hi[j] for j in range(3)]) * fac
        best = (main, main, Fraction(0), None)
    # (B2) for every suffix j0..2 of deep coordinates (valid for any box: deep |x| <= h_{j0}, retained radii >= rho)
    for j0 in range(3):
        rho = Fraction(1) if j0 == 0 else lo[j0 - 1]
        if rho == 0 or any(lo[j] == 0 for j in range(j0)):
            continue
        eta = hi[j0] / rho
        if eta >= 1:
            continue
        m = 3 - j0
        ret = [Fraction(1)] + [s[j] * hi[j] for j in range(j0)]
        fac = Fraction(1)
        for j in range(j0):
            fac *= hi[j] / lo[j]
        main = eng.kernel(ret + [Fraction(0)] * m) * fac
        if CONSTANTS == "S2":
            err = collapse_err_signed(eta, r, [1] + [s[j] for j in range(j0)], [s[j] for j in range(j0, 3)])
        else:
            same = len(ret) == 2 and (ret[0] > 0) == (ret[1] > 0)
            err = collapse_err(m, eta, r, same)
        if best is None or main + err < best[0]:
            best = (main + err, main, err, j0)
    if best is None:                                      # e.g. coordinate 0 includes 0 with eta >= 1
        return Fraction(10**9), Fraction(10**9), Fraction(0), next(j for j in range(3) if lo[j] == 0)
    return best


def rat_between(a, b, x):
    """A simple rational strictly between a and b, close to the float x."""
    D = 1000
    while True:
        m = Fraction(x).limit_denominator(D)
        if a < m < b:
            return m
        D *= 10
        if D > 10**15:
            return (a + b) / 2


def geo_mid(a, b):
    if a == 0:
        return rat_between(Fraction(0), b, float(b) / 16)
    return rat_between(a, b, (float(a) * float(b)) ** 0.5)


def choose_split(box, main, err, j0, kappa):
    """Coordinate to split (never a zero-width one), or None (box fails).
    Collapse-dominated (deep variant with main < kappa but err too big): cut the first coordinate with l = 0 (at h/16),
    else the first suffix coordinate j >= j0 of positive width. Otherwise: the widest retained coordinate of log-width
    >= MIN_LOGW; failing that, any coordinate with l = 0; else None."""
    big = 1e9
    lw = [log(float(b[1] / b[0])) if b[0] > 0 else (big if b[1] > 0 else -1.0) for b in box]   # l = 0: "infinite"
    if j0 is None:                                   # (B1) was best: split the widest coordinate
        cand = [j for j in range(3) if lw[j] >= MIN_LOGW]
        return max(cand, key=lambda j: lw[j]) if cand else None
    ret = [j for j in range(j0) if lw[j] >= MIN_LOGW]
    err_dominated = main < kappa and 2 * err > kappa - main
    cand = [j for j in range(3) if lw[j] >= MIN_LOGW]
    if err_dominated and lw[j0] >= 0.05:
        return j0          # eta = h_{j0} / l_{j0-1} drives the error and h_{j0} can still shrink
    if err_dominated or not ret:     # otherwise resolve by (B1): widest overall
        return max(cand, key=lambda j: lw[j]) if cand else None
    return max(ret, key=lambda j: lw[j])


def run(r, kappa, exact, sign_range=(0, 8), a2_range=(Fraction(0), Fraction(1))):
    eng = Engine(r, exact)
    stats = {"boxes": 0, "pruned": 0, "failed": [], "max_bound_pruned": 0.0}
    t0 = time.time()
    for s in list(product((1, -1), repeat=3))[sign_range[0]:sign_range[1]]:
        stack = [[a2_range, (Fraction(0), Fraction(1)), (Fraction(0), Fraction(1))]]
        while stack:
            box = tighten(stack.pop())
            if box is None:
                continue
            stats["boxes"] += 1
            if stats["boxes"] % 5000 == 0:
                print(f"  .. {stats['boxes']} boxes, stack {len(stack)}, calls {eng.calls}, box "
                      f"{[(float(a), float(b)) for a, b in box]} {time.time() - t0:.0f}s", flush=True)
            val, main, err, j0 = bound(eng, s, box)
            if val <= kappa:
                stats["pruned"] += 1
                stats["max_bound_pruned"] = max(stats["max_bound_pruned"], float(val))
                continue
            j = choose_split(box, main, err, j0, kappa)
            if j is None:
                stats["failed"].append({"s": s, "box": [[str(a), str(b)] for a, b in box], "bound": float(val)})
                continue
            pt = geo_mid(*box[j])
            a, b = box[j]
            left, right = list(box), list(box)
            left[j], right[j] = (a, pt), (pt, b)
            stack.extend([left, right])
        print(f"signs {s}: boxes {stats['boxes']} pruned {stats['pruned']} failed {len(stats['failed'])} "
              f"kernel calls {eng.calls} {time.time() - t0:.0f}s", flush=True)
    stats.update({"r": str(r), "kappa": str(kappa), "exact": exact, "kernel_calls": eng.calls,
                  "seconds": time.time() - t0, "ETA_MAX": str(ETA_MAX),
                  "certified": exact and not stats["failed"]})
    return stats


def main():
    r = Fraction(int(sys.argv[1]), int(sys.argv[2]))
    kappa = Fraction(int(sys.argv[3]), int(sys.argv[4]))
    exact = "--float" not in sys.argv
    out = sys.argv[sys.argv.index("--log") + 1] if "--log" in sys.argv else None
    sr = (0, 8)
    if "--signs" in sys.argv:
        a, b = sys.argv[sys.argv.index("--signs") + 1].split("-")
        sr = (int(a), int(b))
    ar = (Fraction(0), Fraction(1))
    if "--a2" in sys.argv:
        a, b = sys.argv[sys.argv.index("--a2") + 1].split(":")
        ar = (Fraction(a), Fraction(b))
    st = run(r, kappa, exact, sr, ar)
    st["a2_range"] = [str(ar[0]), str(ar[1])]
    print(json.dumps({k: v for k, v in st.items() if k != "failed"} | {"n_failed": len(st["failed"])}), flush=True)
    if out:
        try:
            with open(out, "w") as fh:
                json.dump(st, fh, default=str)
        except OSError as err:
            sys.exit(f"cannot write {out}: {err}")


if __name__ == "__main__":
    main()
