"""Step 14 (DT13) DEV/exact exploration at n = 4: the uniform-calibrated rejection probability P_F(T >= r_a) for
piecewise-constant unimodal densities with K pieces, computed EXACTLY (exact_poly volumes, Fractions), then
maximised by a coarse-to-fine search over the parameters (the search itself is DEV: a maximum over a searched
family is a LOWER bound on the size, never the size; A-3).

Density: breakpoints 0 = z0 < z1 < ... < zK = 1, levels c_j > 0, unimodal (nondecreasing then nonincreasing).
P(T >= r) = sum over placements of the 4 order statistics into pieces (k_1 <= ... <= k_4) of
  4! * prod_i c_{k_i} * vol{order, x_i in piece k_i, s2 >= r s1, s2 >= r s3}.
r: alpha = .10 -> r = 3 exactly; .05, .01 -> rational r just ABOVE r_a (so the value is a lower bound).
Usage: python dt13_n4_kstep.py alpha_index K iterations -> deposits/dt13/n4_kstep_a{i}_K{K}.json
"""
import json
import os
import random
import sys
from fractions import Fraction as F
from itertools import combinations_with_replacement

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt13_n4_exact import r_bracket  # noqa: E402
from exact_poly import volume  # noqa: E402

ALPHAS = (F(1, 100), F(1, 20), F(1, 10))


def prob(z, c, r):
    """z: breakpoints [0, z1, ..., 1] (Fractions); c: levels (normalised inside). Exact P(T >= r)."""
    K = len(c)
    mass = sum(c[j] * (z[j + 1] - z[j]) for j in range(K))
    c = [cj / mass for cj in c]
    tot = F(0)
    for ks in combinations_with_replacement(range(K), 4):
        A = [[-1, 0, 0, 0], [1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1], [0, 0, 0, 1]]
        b = [0, 0, 0, 0, 1]
        for i, k in enumerate(ks):
            row = [0, 0, 0, 0]
            row[i] = 1
            A.append(row)
            b.append(z[k + 1])
            row = [0, 0, 0, 0]
            row[i] = -1
            A.append(row)
            b.append(-z[k])
        A += [[-r, r + 1, -1, 0], [0, 1, -(1 + r), r]]
        b += [0, 0]
        v = volume(A, b)
        if v:
            w = F(24)
            for k in ks:
                w *= c[k]
            tot += w * v
    return tot


def unimodal(c):
    i = max(range(len(c)), key=lambda j: c[j])
    return all(c[j] <= c[j + 1] for j in range(i)) and all(c[j] >= c[j + 1] for j in range(i, len(c) - 1))


def rnd(x, d=10 ** 6):
    return F(round(x * d), d)


def search(K, r, iters, seed):
    rng = random.Random(seed)
    best_p, best_z, best_c = F(-1), [F(0), F(1, 2), F(1)], [F(1), F(1)]
    for _ in range(iters):
        zs = sorted(rng.random() for _ in range(K - 1))
        z = [F(0)] + [rnd(v) for v in zs] + [F(1)]
        if any(z[j + 1] <= z[j] for j in range(K)):
            continue
        c = [F(1)] + [rnd(10 ** rng.uniform(-1, 4), 10 ** 3) for _ in range(K - 1)]
        c = [max(v, F(1, 1000)) for v in c]
        if not unimodal(c):
            c = sorted(c)          # make it monotone (mode at the right end)
        p = prob(z, c, r)
        if p > best_p:
            best_p, best_z, best_c = p, z, c
    # local refinement around the best
    p0, z, c = best_p, best_z, best_c
    for _step in range(iters):
        z2 = [z[0]] + [min(max(rnd(float(v) + rng.gauss(0, 0.02)), F(1, 10 ** 6)), F(1) - F(1, 10 ** 6)) for v in z[1:-1]] + [z[-1]]
        if any(z2[j + 1] <= z2[j] for j in range(K)):
            continue
        c2 = [max(rnd(float(v) * 10 ** rng.gauss(0, 0.15), 10 ** 3), F(1, 1000)) for v in c]
        if not unimodal(c2):
            continue
        p = prob(z2, c2, r)
        if p > p0:
            p0, z, c = p, z2, c2
    return p0, z, c


def main():
    """argv: alpha index (0, 1, 2), K, iterations. Seeds: DEV root 2026110745 + 10 K + alpha index."""
    ia, K, it = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    a = ALPHAS[ia]
    _lo, hi = r_bracket(a, F(1, 10 ** 12))
    r = hi                       # P is nonincreasing in r: the value at r_hi is a valid lower bound
    p, z, c = search(K, r, it, seed=2026110745 + 10 * K + ia)
    row = {"alpha": float(a), "K": K, "iters": it, "P_lower": float(p), "P_lower_exact": str(p),
           "z": [str(v) for v in z], "c": [str(v) for v in c]}
    print(f"alpha={float(a):.2f} K={K}: best P = {float(p):.6f}  z={[round(float(v), 5) for v in z]} "
          f"c={[round(float(v), 3) for v in c]}", flush=True)
    od = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits", "dt13")
    try:
        os.makedirs(od, exist_ok=True)
        with open(os.path.join(od, f"n4_kstep_a{ia}_K{K}.json"), "w") as fh:
            json.dump(row, fh, indent=1)
    except OSError as err:
        sys.exit(f"cannot write: {err}")


if __name__ == "__main__":
    main()
