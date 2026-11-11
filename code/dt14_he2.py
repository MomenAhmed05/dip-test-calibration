"""Step 15 / DT14 WP-B2: exact certificate of HE_2 (lemma A, m = 2) at rational r.

HE_2: for every two-radius Khintchine mixing law G = p delta_1 + (1-p) delta_{s t}
(s = +1 same side, s = -1 opposite side, 0 < t <= 1 after scaling and reflection),
    Phi_r(G) = sum_k C(4,k) p^k q^(4-k) K_k(t, s) <= M_2(r) = max_p { g_r(p) + q^4 a(r) },
where K_k(t, s) = K_r(1 x k, s t x (4-k)), K_4 = K_0 = a(r) (scale + reflection invariance),
g_r(p) = a (p^4 + 4 p^3 q) + 6 b p^2 q^2, a = 2/((r+1)(r+2)), b = 1/(r+1).

Method (all exact, Fractions):
 1. t^(4-k) K_k(t, s) is the x-space volume of a parametric polytope whose right-hand sides are
    linear in t, so it is a polynomial of degree <= 4 in t between breakpoints. Candidate
    breakpoints B(r) = {1/(r+2), 1/(r+1), 1/r, r/(r+1), (r+1)/(r+2)} (found by kink scans at r = 3).
    On each piece: interpolate from 5 exact kernel values, then VERIFY at 4 further rational points
    (exact equality required; any mismatch -> FAIL, the breakpoint list is incomplete).
 2. Bernstein subdivision on (p, t) boxes, exact:
      region N (near t = 0, first piece, p in [pL, pH]):  (H(p) - Phi(p, t)) / t >= 0,
        with H(p) = g(p) + q^4 a <= M_2 (so Phi <= M_2 there);
      everywhere else: t^3 (M_lo - Phi(p, t)) >= 0 with M_lo = H(p0) <= M_2 (p0 rational).
    A box is accepted when all Bernstein coefficients are > 0 (or >= 0 for region N on t = 0 edge
    handled by the division by t). Max depth bounded; failure is reported, never hidden.
Run: python dt14_he2.py [r_num r_den]   (default r = 3, alpha = 10%)
"""
import json
import sys
import time
from fractions import Fraction as F
from math import comb

import sympy as sp

sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
from dt14_n4_kernel import kernel_exact  # noqa: E402

TS = sp.symbols("t")
PS = sp.symbols("p")
CRIT_LOG = {}


def R(x):
    return sp.Rational(x.numerator, x.denominator)


def breakpoints(r):
    c = [F(1) / (r + 2), F(1) / (r + 1), F(1) / r, r / (r + 1), (r + 1) / (r + 2)]
    return sorted({x for x in c if 0 < x < 1})


def _det(m):
    """Exact determinant (Fraction Gaussian elimination)."""
    m = [row[:] for row in m]
    n, d = len(m), F(1)
    for c in range(n):
        piv = next((i for i in range(c, n) if m[i][c] != 0), None)
        if piv is None:
            return F(0)
        if piv != c:
            m[c], m[piv] = m[piv], m[c]
            d = -d
        d *= m[c][c]
        for i in range(c + 1, n):
            f = m[i][c] / m[c][c]
            if f:
                m[i] = [x - f * y for x, y in zip(m[i], m[c], strict=True)]
    return d


def critical_ts(k, s, r):
    """All t in (0,1) where 5 hyperplanes of the x-space arrangement (box faces x_i = 0, x_i = z_i(t);
    order planes x_i = x_j; spacing planes s2 = r s1, s2 = r s3 for every order) have a common point.
    Between consecutive values the combinatorial type of every order polytope is constant, its
    vertices are affine in t, hence t^(4-k) K_k is a polynomial of degree <= 4 there (rigour step)."""
    from itertools import combinations, permutations
    zt = [(F(1), F(0))] * k + [(F(0), F(s))] * (4 - k)  # z_i = c0 + c1 t
    box = []
    for i in range(4):  # x_i = 0 and x_i = z_i
        e = tuple(F(int(j == i)) for j in range(4))
        box.append((e, F(0), F(0)))
        box.append((e, zt[i][0], zt[i][1]))

    def sp_row(a_, b_):  # x_b - x_a
        return [F(int(m == b_)) - F(int(m == a_)) for m in range(4)]

    combos = set()
    for o in permutations(range(4)):  # each order polytope has its own 13-plane arrangement
        s1, s2, s3 = sp_row(o[0], o[1]), sp_row(o[1], o[2]), sp_row(o[2], o[3])
        planes = list(box) + [(tuple(x), F(0), F(0)) for x in (s1, s2, s3)]
        planes.append((tuple(x - r * y for x, y in zip(s2, s1, strict=True)), F(0), F(0)))
        planes.append((tuple(x - r * y for x, y in zip(s2, s3, strict=True)), F(0), F(0)))
        combos.update(tuple(sorted(c)) for c in combinations(planes, 5))
    out = set()
    for combo in combos:
        # det [A | b0 + b1 t] = D0 + D1 t (linear in t: only the last column depends on t)
        A = [list(pl[0]) for pl in combo]
        D0 = _det([row + [pl[1]] for row, pl in zip(A, combo, strict=True)])
        D1 = _det([row + [pl[2]] for row, pl in zip(A, combo, strict=True)])
        if D1 != 0:
            t = -D0 / D1
            if 0 < t < 1:
                out.add(t)
    return sorted(out)


def piece_polys(k, s, r):
    """Exact pieces [(lo, hi, poly in t)] of t^(4-k) K_k(t, s) on (0, 1]; verified."""
    crit = critical_ts(k, s, r)
    CRIT_LOG[(k, s)] = [str(x) for x in crit]
    bps = [F(0)] + sorted(set(breakpoints(r)) | set(crit)) + [F(1)]
    out = []
    for lo, hi in zip(bps[:-1], bps[1:], strict=True):
        w = hi - lo
        fit_t = [lo + w * F(j, 10) for j in (1, 3, 5, 7, 9)]
        chk_t = [lo + w * F(j, 13) for j in (2, 5, 8, 11)]

        def v(t):
            return kernel_exact([F(1)] * k + [s * t] * (4 - k), r) * t ** (4 - k)

        poly = sp.expand(sp.interpolate([(R(t), R(v(t))) for t in fit_t], TS))
        for t in chk_t:
            if poly.subs(TS, R(t)) != R(v(t)):
                raise RuntimeError(f"piece check FAIL k={k} s={s} [{lo},{hi}] t={t}")
        out.append((lo, hi, poly))
    return out


def bern_coeffs(poly2, pr, tr, dp, dt):
    """Bernstein coefficients of poly2(p, t) on box pr x tr (degrees dp, dt), exact."""
    (p0, p1), (t0, t1) = pr, tr
    u, v = sp.symbols("u v")
    sub = sp.Poly(sp.expand(poly2.subs({PS: R(p0) + (R(p1) - R(p0)) * u, TS: R(t0) + (R(t1) - R(t0)) * v})), u, v)
    a = [[sp.Integer(0)] * (dt + 1) for _ in range(dp + 1)]
    for (i, j), c in sub.terms():
        a[i][j] = c
    # power -> Bernstein: b_{IJ} = sum_{i<=I, j<=J} C(I,i)/C(dp,i) C(J,j)/C(dt,j) a_ij
    out = []
    for bi in range(dp + 1):
        for bj in range(dt + 1):
            sv = 0
            for i in range(bi + 1):
                for j in range(bj + 1):
                    if a[i][j] != 0:
                        sv += sp.Rational(comb(bi, i), comb(dp, i)) * sp.Rational(comb(bj, j), comb(dt, j)) * a[i][j]
            out.append(sv)
    return out


def certify_box(poly2, pr, tr, dp, dt, stats, depth=0, maxdepth=14):
    """True if poly2 > 0 on the box (Bernstein subdivision). stats counts boxes."""
    stats["boxes"] += 1
    bc = bern_coeffs(poly2, pr, tr, dp, dt)
    if min(bc) > 0:
        return True
    if depth >= maxdepth:
        stats["fail"].append((str(pr), str(tr), float(min(bc))))
        return False
    (p0, p1), (t0, t1) = pr, tr
    if (p1 - p0) >= (t1 - t0):
        pm = (p0 + p1) / 2
        return (certify_box(poly2, (p0, pm), tr, dp, dt, stats, depth + 1, maxdepth)
                and certify_box(poly2, (pm, p1), tr, dp, dt, stats, depth + 1, maxdepth))
    tm = (t0 + t1) / 2
    return (certify_box(poly2, pr, (t0, tm), dp, dt, stats, depth + 1, maxdepth)
            and certify_box(poly2, pr, (tm, t1), dp, dt, stats, depth + 1, maxdepth))


def main(r):
    t_start = time.time()
    a = F(2) / ((r + 1) * (r + 2))
    b = F(1) / (r + 1)
    q = 1 - PS
    H = sp.expand(R(a) * (PS**4 + 4 * PS**3 * q) + 6 * R(b) * PS**2 * q**2 + q**4 * R(a))
    # p* : root of H' in (0,1); p0 rational close to it; M_lo = H(p0) <= M_2
    crit = [c for c in sp.Poly(sp.diff(H, PS), PS).nroots(n=40) if c.is_real and 0 < c < 1]
    pstar = max(crit, key=lambda c: H.subs(PS, c))
    p0 = F(int(pstar * 10**12), 10**12)
    M_lo = H.subs(PS, R(p0))
    report = {"r": str(r), "a": str(a), "b": str(b), "p_star": float(pstar), "M2_float": float(H.subs(PS, pstar)),
              "M_lo": str(M_lo), "M_lo_float": float(M_lo), "breakpoints": [str(x) for x in breakpoints(r)],
              "pieces": {}, "boxes": {}, "result": None}
    allok = True
    for s in (1, -1):
        pcs = {k: piece_polys(k, s, r) for k in (1, 2, 3)}
        report["pieces"][str(s)] = {str(k): [(str(lo), str(hi), str(pl)) for lo, hi, pl in v] for k, v in pcs.items()}
        nP = len(pcs[1])
        for j in range(nP):
            lo, hi = pcs[1][j][0], pcs[1][j][1]
            # t^3 * Phi as a polynomial in (p, t) on this piece
            t3phi = (q**4 + PS**4) * R(a) * TS**3
            for k in (1, 2, 3):
                t3phi += comb(4, k) * PS**k * q ** (4 - k) * pcs[k][j][2] * TS ** (k - 1)
            t3phi = sp.expand(t3phi)
            if j == 0:
                # region N: (H - Phi)/t = (t^3 H - t3phi)/t^4 must be > 0 for p in [pL,pH], t in (0, hi]
                num = sp.expand(TS**3 * H - t3phi)
                quo, rem = sp.div(sp.Poly(num, TS), sp.Poly(TS**4, TS))
                if rem.as_expr() != 0:
                    raise RuntimeError("first piece: t^4 does not divide t^3(H-Phi)")
                gN = sp.expand(quo.as_expr())
                # band where d/dt Phi(p, 0+) < 0 for this sign (r = 3 run 1: +-0.08 failed at p ~ 0.523,
                # s = -1, because there d/dt Phi > 0; the band is a CLI parameter, default r = 3 values)
                pL, pH = BAND
                stats = {"boxes": 0, "fail": []}
                okN = certify_box(gN, (pL, pH), (F(0), hi), 4, sp.degree(gN, TS) or 1, stats)
                report["boxes"][f"s{s}_N"] = {"ok": okN, "p": [str(pL), str(pH)], "t": [0, str(hi)], **stats}
                allok &= okN
                # rest of the first piece: p outside [pL, pH]
                gM = sp.expand(TS**3 * M_lo - t3phi)
                for pr in ((F(0), pL), (pH, F(1))):
                    if pr[1] <= pr[0]:
                        continue
                    stats = {"boxes": 0, "fail": []}
                    # t^3(M_lo - Phi) vanishes at t = 0 like t^3 M_lo - t3phi(t=0)... divide by t^? not needed:
                    # use (M_lo - Phi) = gM / t^3 with gM/t^3 polynomial-in-t only after division; divide.
                    quo2, rem2 = sp.div(sp.Poly(gM, TS), sp.Poly(TS**3, TS))
                    if rem2.as_expr() != 0:
                        raise RuntimeError("first piece: t^3 does not divide t^3(M_lo-Phi)")
                    g2 = sp.expand(quo2.as_expr())
                    ok = certify_box(g2, pr, (F(0), hi), 4, sp.degree(g2, TS) or 1, stats)
                    report["boxes"][f"s{s}_P0_{pr[0]}"] = {"ok": ok, **stats}
                    allok &= ok
            else:
                gM = sp.expand(TS**3 * M_lo - t3phi)
                stats = {"boxes": 0, "fail": []}
                ok = certify_box(gM, (F(0), F(1)), (lo, hi), 4, sp.degree(gM, TS) or 1, stats)
                report["boxes"][f"s{s}_piece{j}"] = {"ok": ok, "t": [str(lo), str(hi)], **stats}
                allok &= ok
            print(f"s={s} piece {j} [{lo},{hi}] done {time.time() - t_start:.0f}s", flush=True)
    report["critical_ts"] = {f"k{k}_s{s}": v for (k, s), v in CRIT_LOG.items()}
    report["result"] = "HE2 CERTIFIED" if allok else "NOT CERTIFIED (see fails)"
    report["seconds"] = round(time.time() - t_start, 1)
    print(json.dumps({k: v for k, v in report.items() if k != "pieces"}, indent=1, default=str))
    with open(__file__.rsplit("\\", 1)[0].rsplit("/", 1)[0] + f"/../deposits/dt14/he2_r{r.numerator}_{r.denominator}.json",
              "w") as fh:
        json.dump(report, fh, indent=1, default=str)


BAND = (F(545, 1000), F(655, 1000))

if __name__ == "__main__":
    rr = F(int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else F(3)
    if len(sys.argv) > 4:
        BAND = (F(sys.argv[3]), F(sys.argv[4]))
    main(rr)
