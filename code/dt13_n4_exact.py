"""Step 14 (DT13) Part X, n = 4: EXACT rejection probabilities of the uniform-calibrated dip test.

Facts used (all checked rows 5-6):
  (a) H&H 1985: at n = 4, T = max{1, s2 / max(s1, s3)} (s_i = spacings of the order statistics) and the dip
      is D = T / (4 (T + 1))  (S check: DEV root 2026110740, max error 5.6e-17 on 20,000 samples).
  (b) Under U(0,1): P(T >= r) = 2 / ((r + 1)(r + 2)) for r >= 1 (S derivation from iid exponential spacings;
      reproduced EXACTLY below by the polytope code as a check).
  So q_4(1 - alpha) = r_a / (4 (r_a + 1)) with (r_a + 1)(r_a + 2) = 2 / alpha, the uniform size is exactly alpha
  (H_4 is continuous above 1/8), and for any law: P(reject) = P(T >= r_a).
Two-step law W(w, eps): density (1 - w)/L on (0, L), w/eps on (L, 1), L = 1 - eps (CONVENTIONS CV-9).
  P_W(T >= r) = sum_k 24 c1^k c2^(4-k) vol{0 <= x1 <= ... <= x4 <= 1, x_k <= L <= x_(k+1), s2 >= r s1, s2 >= r s3}
  computed exactly (exact_poly.volume, Fractions). r_a irrational for alpha = .05, .01: bracketed by rationals
  r_lo < r_a < r_hi (width < 1e-30); P is nonincreasing in r, so P(T >= r_hi) <= P <= P(T >= r_lo).
Atom family F_p = p U + (1 - p) delta_1 (H&H; atom row, CV-7): P = [p^4 + 4 p^3 (1 - p)] alpha + 6 p^2 (1 - p)^2 / (1 + r_a)
  (S derivation: m = #atom points ~ Bin(4, 1-p); m = 0, 1 give the uniform spacing law; m = 2 gives
  P(E3 >= r E2) = 1/(1 + r); m >= 3 gives T = 1; checked by MC below with diptest on tied samples).
Output: deposits/dt13/n4_exact.json (exact rationals as strings + floats). MC cross-check: DEV root 2026110741.
"""
import json
import os
import sys
from fractions import Fraction as F

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt0c_screen import dips_of  # noqa: E402
from exact_poly import volume  # noqa: E402

ALPHAS = (F(1, 100), F(1, 20), F(1, 10))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits", "dt13")


def r_bracket(a, width=F(1, 10 ** 30)):
    """Rationals r_lo < r_a <= r_hi with (r+1)(r+2) = 2/a; exact when r_a is rational."""
    g = lambda r: (r + 1) * (r + 2) - 2 / a  # noqa: E731  (increasing for r >= 1)
    lo, hi = F(1), F(1000)
    if g(lo) == 0:
        return lo, lo
    while hi - lo > width:
        mid = (lo + hi) / 2
        if g(mid) == 0:
            return mid, mid
        if g(mid) < 0:
            lo = mid
        else:
            hi = mid
    return lo, hi


def region(k, r, L):
    """Constraints A x <= b for: order, x_k <= L <= x_(k+1) (k points below L), s2 >= r s1, s2 >= r s3."""
    A = [[-1, 0, 0, 0], [1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1], [0, 0, 0, 1]]
    b = [0, 0, 0, 0, 1]
    if k >= 1:
        row = [0, 0, 0, 0]
        row[k - 1] = 1
        A.append(row)
        b.append(L)
    if k <= 3:
        row = [0, 0, 0, 0]
        row[k] = -1
        A.append(row)
        b.append(-L)
    A += [[-r, r + 1, -1, 0], [0, 1, -(1 + r), r]]
    b += [0, 0]
    return A, b


def p_two_step(w, eps, r):
    L = 1 - eps
    c1, c2 = (1 - w) / L, w / eps
    tot = F(0)
    for k in range(5):
        A, b = region(k, r, L)
        tot += 24 * c1 ** k * c2 ** (4 - k) * volume(A, b)
    return tot


def p_atom(p, a, r):
    return (p ** 4 + 4 * p ** 3 * (1 - p)) * a + 6 * p ** 2 * (1 - p) ** 2 / (1 + r)


def mc_check(rows, brs):
    """DEV MC cross-check (root 2026110741): rows = (name, w, eps, alpha, P_lo, P_hi)."""
    rng = np.random.default_rng(np.random.SeedSequence([2026110741, 4]))
    N = 400_000
    out = []
    for name, w, eps, a, plo, phi in rows:
        L = 1 - eps
        u = rng.random((N, 4))
        top = rng.random((N, 4)) < float(w)
        x = np.where(top, float(L) + float(eps) * u, float(L) * u)
        d = dips_of(x)
        rr = float(brs[a][0])
        q = rr / (4 * (rr + 1))
        mc = float(np.mean(d >= q - 1e-12))
        pe = float(plo)
        se = (pe * (1 - pe) / N) ** 0.5
        out.append({"law": name, "alpha": float(a), "exact": [float(plo), float(phi)], "mc": mc,
                    "z": (mc - pe) / se})
    # atom family, p = 0.9, alpha = 0.10 (P's example 0.10692) and p = 0.62
    for p, a in ((F(9, 10), F(1, 10)), (F(62, 100), F(1, 10)), (F(62, 100), F(1, 20))):
        rlo, _ = r_bracket(a)
        ex = p_atom(p, a, rlo)
        x = rng.random((N, 4))
        x = np.where(rng.random((N, 4)) < 1 - float(p), 1.0, x)
        d = dips_of(x)
        q = float(rlo) / (4 * (float(rlo) + 1))
        mc = float(np.mean(d >= q - 1e-12))
        se = (float(ex) * (1 - float(ex)) / N) ** 0.5
        out.append({"law": f"atom p={float(p)}", "alpha": float(a), "exact": float(ex), "mc": mc,
                    "z": (mc - float(ex)) / se})
    return out


def main():
    try:
        os.makedirs(OUT, exist_ok=True)
    except OSError as err:
        sys.exit(f"cannot create {OUT}: {err}")
    res = {"n": 4, "alphas": {}, "laws": [], "w_scan": [], "atom_scan": []}
    brs = {}
    for a in ALPHAS:
        lo, hi = r_bracket(a)
        brs[a] = (lo, hi)
        res["alphas"][str(float(a))] = {"r_lo": str(lo), "r_hi": str(hi), "r_float": float(lo),
                                         "q_dip": float(lo / (4 * (lo + 1))), "sqrt_n_q": 2 * float(lo / (4 * (lo + 1))),
                                         "uniform_size_exact": str(a)}
    laws = [("W(0.1,0.005)", F(1, 10), F(1, 200)), ("W(0.4,0.005)", F(2, 5), F(1, 200)),
            ("uniform as W(1/2,1/2) (polytope check)", F(1, 2), F(1, 2))]
    for name, w, eps in laws:
        for a in ALPHAS:
            lo, hi = brs[a]
            plo, phi = p_two_step(w, eps, hi), p_two_step(w, eps, lo)
            res["laws"].append({"law": name, "alpha": float(a), "P_lo": str(plo), "P_hi": str(phi),
                                "P": [float(plo), float(phi)], "width": float(phi - plo)})
            print(f"{name:22s} alpha={float(a):.2f}: P(reject) in [{float(plo):.10f}, {float(phi):.10f}]", flush=True)
    # scan w at eps = 0.005 and eps -> small, alpha all levels (exact; descriptive, not a search on noise)
    for eps in (F(1, 200), F(1, 10 ** 6)):
        for wi in range(5, 100, 5):
            w = F(wi, 100)
            row: dict[str, object] = {"eps": float(eps), "w": float(w)}
            for a in ALPHAS:
                lo, hi = brs[a]
                row[str(float(a))] = [float(p_two_step(w, eps, hi)), float(p_two_step(w, eps, lo))]
            res["w_scan"].append(row)
    # uniform end-to-end check: at both rational endpoints the polytope value equals 2/((r+1)(r+2)) EXACTLY,
    # and the resulting probability interval contains alpha
    res["uniform_check"] = []
    for a in ALPHAS:
        lo, hi = brs[a]
        v_hi, v_lo = p_two_step(F(1, 2), F(1, 2), hi), p_two_step(F(1, 2), F(1, 2), lo)
        ok = (v_hi == 2 / ((hi + 1) * (hi + 2)) and v_lo == 2 / ((lo + 1) * (lo + 2)) and v_hi <= a <= v_lo)
        res["uniform_check"].append({"alpha": float(a), "interval": [str(v_hi), str(v_lo)],
                                     "width": float(v_lo - v_hi), "identity_and_contains_alpha": bool(ok)})
        print(f"uniform check alpha={float(a)}: identity + contains alpha = {ok}, width {float(v_lo - v_hi):.2e}", flush=True)
    # exact values for every law that the CONFIRM run draws at n = 4 (validator V2: [P_lo, P_hi] inside the MC enclosure)
    w4 = {0: "0.473789162301", 1: "0.433348270040", 2: "0.396743906214"}   # dt13_confirm.W4 (Part L, n = 4)
    claws = [("W(0.1,0.005)", F(1, 10), F(1, 200)), ("W(0.4,0.005)", F(2, 5), F(1, 200)),
             ("path W(min(.4,1.5/sqrt n),1e-4)", F(2, 5), F(1, 10 ** 4))]
    claws += [(f"W{(float(w4[i]), 1e-06)}", F(w4[i]), F(1, 10 ** 6)) for i in range(3)]
    res["confirm_laws"] = []
    for name, w, eps in claws:
        for a in ALPHAS:
            lo, hi = brs[a]
            plo, phi = p_two_step(w, eps, hi), p_two_step(w, eps, lo)
            res["confirm_laws"].append({"law": name, "alpha": float(a), "P": [float(plo), float(phi)],
                                        "P_lo": str(plo), "P_hi": str(phi), "width": float(phi - plo)})
    # supremum of the two-step family (eps -> 0 limit, maximised over the top mass): closed form, see DT13 Part X
    res["two_step_sup"] = {"0.01": 0.03135178825414085, "0.05": 0.0845928225214854, "0.1": 0.1364852075326348,
                           "formula": "max_p a(p^4+4p^3q+q^4)+6p^2q^2/(1+r_a), q=1-p (sympy nroots, 25 digits)"}
    for pi in range(5, 100, 5):
        p = F(pi, 100)
        res["atom_scan"].append({"p": float(p), **{str(float(a)): float(p_atom(p, a, brs[a][1])) for a in ALPHAS}})
    mcrows = [(nm, w, e, a, p_two_step(w, e, brs[a][1]), p_two_step(w, e, brs[a][0]))
              for nm, w, e in laws[:2] for a in ALPHAS]
    res["mc_check"] = mc_check(mcrows, brs)
    for m in res["mc_check"]:
        print("MC check:", m, flush=True)
    try:
        with open(os.path.join(OUT, "n4_exact.json"), "w") as fh:
            json.dump(res, fh, indent=1)
    except OSError as err:
        sys.exit(f"cannot write n4_exact.json: {err}")
    best = {}
    for a in ALPHAS:
        k = str(float(a))
        bw = max(res["w_scan"], key=lambda r, k=k: r[k][0])
        ba = max(res["atom_scan"], key=lambda r, k=k: r[k])
        best[k] = (bw["eps"], bw["w"], bw[k][0], ba["p"], ba[k])
        print(f"alpha={k}: best two-step on scan eps={bw['eps']} w={bw['w']}: {bw[k][0]:.6f}; "
              f"best atom p={ba['p']}: {ba[k]:.6f}")


if __name__ == "__main__":
    main()
