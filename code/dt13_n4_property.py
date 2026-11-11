"""Step 14 (DT13) property test of Lemma N4. DEV root 2026110748. No claim rests on this test;
it tries to FALSIFY the lemma's sample-wise implication and the block containment.

Implication tested: for H = F^{-1} concave on (0, u] and convex on [u, 1) (F in U_c with u = F(m)), and r > 1,
  T(H(U)) >= r  ==>  phi(u) >= r,
with T = max{1, s2 / max(s1, s3)}, s_i = spacings of H(U_(i)), phi as in dt13_n4_bound.
Random H: piecewise linear with random knots; slopes on (0, u] nonincreasing (concave), on [u, 1) nondecreasing
(convex), drawn log-uniform over 12 decades (spikes, plateaus of the density); u includes 0, 1, and u near the
sample points; r includes 1 + 1e-9 and the three r_alpha.
Also: block containment E_v subset E_L u E_R for v in [L, R] (random L < v < R), and the r = 1 counterexample.
"""
import numpy as np
from dt13_n4_bound import phi  # noqa: I001

ROOT = 2026110748


def rand_H(rng, u):
    """Random increasing piecewise-linear H on [0, 1]: concave on [0, u], convex on [u, 1]. Returns callable."""
    kl = rng.integers(0, 6)
    kr = rng.integers(0, 6)
    xl = np.sort(rng.uniform(0, u, kl)) if u > 0 else np.array([])
    xr = np.sort(rng.uniform(u, 1, kr)) if u < 1 else np.array([])
    knots = np.concatenate([[0.0], xl, [u] if 0 < u < 1 else [], xr, [1.0]])
    knots = np.unique(knots)
    m = len(knots) - 1
    sl = 10 ** rng.uniform(-6, 6, m)
    mid = (knots[:-1] + knots[1:]) / 2
    left = mid < u
    sl_left = np.sort(sl[left])[::-1]        # nonincreasing slopes left of u (H concave there)
    sl_right = np.sort(sl[~left])            # nondecreasing slopes right of u (H convex there)
    slopes = np.concatenate([sl_left, sl_right])
    vals = np.concatenate([[0.0], np.cumsum(slopes * np.diff(knots))])
    return lambda v: np.interp(v, knots, vals)


def main():
    rng = np.random.default_rng(np.random.SeedSequence([ROOT, 1]))
    rs = [1 + 1e-9, 1.5, 3.0, (np.sqrt(161) - 3) / 2, (np.sqrt(801) - 3) / 2]
    viol, tested, hits = 0, 0, 0
    worst = None
    for _ in range(4000):
        u = float(rng.choice([0.0, 1.0, rng.uniform(), rng.uniform(0.45, 0.55)]))
        H = rand_H(rng, u)
        U = np.sort(rng.random((500, 4)), axis=1)
        if rng.random() < 0.3:              # put u very close to a sample point
            U[:, 1] = np.clip(u - 10 ** rng.uniform(-12, -3), U[:, 0] + 1e-15, None)
            U = np.sort(U, axis=1)
        X = H(U)
        s = np.diff(X, axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            T = np.maximum(1.0, s[:, 1] / np.maximum(s[:, 0], s[:, 2]))
        ph = phi(U, u)
        for r in rs:
            rej = np.greater_equal(T, r)
            ok = ~rej | (ph >= r * (1 - 1e-9))
            tested += len(U)
            hits += int(rej.sum())
            if not ok.all():
                viol += int((~ok).sum())
                worst = (u, r)
    print(f"implication T>=r => phi(u)>=r: {tested} checks, {hits} rejections, violations {viol} {worst or ''}")
    # block containment
    bad = 0
    for _ in range(2000):
        U = np.sort(rng.random((200, 4)), axis=1)
        L, R = np.sort(rng.random(2))
        v = rng.uniform(L, R)
        for r in rs[1:]:
            ev = phi(U, v) >= r
            eI = (phi(U, L) >= r) | (phi(U, R) >= r)
            bad += int((ev & ~eI).sum())
    print(f"block containment E_v subset E_L u E_R: violations {bad}")
    # r = 1 counterexample: uniform, u = 1: P(T >= 1) = 1 but P(phi(1) >= 1) = P(g2 >= g1) = 1/2
    U = np.sort(rng.random((200000, 4)), axis=1)
    print(f"r = 1 counterexample: P(phi(1) >= 1) = {np.mean(phi(U, 1.0) >= 1):.4f} (P(T >= 1) = 1)")


if __name__ == "__main__":
    main()
