"""Cheng & Hall (1998) calibration of the dip test ("Calibrating the excess mass and dip tests of modality",
JRSS-B 60, 579-589), as described in their Sections 2-3 (read from prior/priorart/step1/dl_cheng_hall_1998).

Their recipe:
  1. Normal-kernel estimates f_hat (bandwidth h0) and f2_hat = estimate of f'' (bandwidth h2), each with its
     asymptotically optimal global (AMISE) bandwidth under a normal N(0, S^2) reference, S^2 = sample variance;
     values on an equally spaced grid of 512 points.
  2. x0_hat = argmax f_hat (the largest mode); d_hat = |f2_hat(x0_hat)| / f_hat(x0_hat)^3.
  3. Calibration law g with the same d: symmetric Beta(b, b) for d in (0, 2 pi); normal for d = 2 pi;
     Student t_nu for d > 2 pi (d is scale invariant).
  4. Critical point z_alpha = (1 - alpha) quantile of the dip of an n-sample from g; reject if dip > z_alpha.

Implementation choices (ours, stated): the bootstrap of step 4 depends on the data only through d_hat, so we
replace the per-dataset resample by a pre-computed table of dip quantiles on a log-spaced d grid (per n), with
linear interpolation in log d (B -> large limit of their procedure). d below the grid uses the d_min row; above it
the d_max row. Kernel estimates use linear binning onto the 512-point grid and a discrete convolution.
Normal-reference AMISE bandwidths with a Gaussian kernel: h0 = (4/3)^(1/5) S n^(-1/5) = 1.0592 S n^(-1/5);
h2 = (4/7)^(1/9) S n^(-1/9) = 0.9397 S n^(-1/9)  [r = 2: R(K'') = 3/(8 sqrt(pi)), R(phi_S^(4)) = 105/(32 sqrt(pi) S^9)].
"""
import math

import numpy as np
from scipy.optimize import root_scalar
from scipy.special import betaln, gammaln

GRID = 512
TWO_PI = 2.0 * math.pi
D_GRID = np.geomspace(1e-3, 1000.0, 49)  # log-spaced calibration grid for d (8 points per decade)
D_NODES = np.concatenate([[0.0], D_GRID])  # node 0 = d = 0 = Beta(1, 1) = uniform (floor must reach it)
NEAR_NORMAL = 1e-3  # |d - 2 pi| < NEAR_NORMAL * 2 pi -> normal reference (inversion unstable there)
H0_C = (4.0 / 3.0) ** 0.2
H2_C = (4.0 / 7.0) ** (1.0 / 9.0)


def d_beta(b):
    """d = |g''(1/2)| / g(1/2)^3 for symmetric Beta(b, b): 8 (b - 1) / g(1/2)^2, g(1/2) = 4^(1-b) / B(b, b)."""
    log_g = (1.0 - b) * math.log(4.0) - betaln(b, b)
    return 8.0 * (b - 1.0) * math.exp(-2.0 * log_g)


def d_t(nu):
    """d for Student t_nu at 0: ((nu + 1) / nu) / g(0)^2."""
    log_g = gammaln((nu + 1) / 2) - gammaln(nu / 2) - 0.5 * math.log(nu * math.pi)
    return (nu + 1.0) / nu * math.exp(-2.0 * log_g)


def reference(d):
    """('beta', b) / ('normal', None) / ('t', nu) with the given d."""
    if d <= 0.0:
        return ("uniform", None)
    if abs(d - TWO_PI) < NEAR_NORMAL * TWO_PI:
        return ("normal", None)
    if d < TWO_PI:
        return ("beta", _root(lambda b: d_beta(b) - d, 1.0 + 1e-12, 1e7))
    return ("t", _root(lambda nu: d_t(nu) - d, 1e-4, 1e9))


def _root(fun, a, b):
    return float(root_scalar(fun, bracket=(a, b), method="brentq", xtol=1e-12).root)


def sample_reference(rng, ref, size):
    kind, p = ref
    if kind == "beta":
        return rng.beta(p, p, size)
    if kind == "normal":
        return rng.standard_normal(size)
    if kind == "uniform":
        return rng.random(size)
    return rng.standard_t(p, size)


def _binned(x, lo, hi):
    """Linear binning of the sample x onto GRID points spanning [lo, hi]. Returns (grid, counts, delta)."""
    g = np.linspace(lo, hi, GRID)
    delta = g[1] - g[0]
    pos = (x - lo) / delta
    i = np.clip(np.floor(pos).astype(int), 0, GRID - 2)
    w = pos - i
    c = np.bincount(i, weights=1.0 - w, minlength=GRID) + np.bincount(i + 1, weights=w, minlength=GRID)
    return g, c[:GRID], delta


def _kde(counts, delta, h, n, deriv):
    L = int(min(GRID - 1, math.ceil(5.0 * h / delta)))
    u = np.arange(-L, L + 1) * delta / h
    phi = np.exp(-0.5 * u * u) / math.sqrt(2.0 * math.pi)
    k = phi / (n * h) if deriv == 0 else (u * u - 1.0) * phi / (n * h ** 3)
    full = np.convolve(counts, k, mode="full")  # bug fix: mode="same" changes length when 2L+1 > GRID
    out = full[L:L + counts.size]
    if out.size != counts.size:
        raise RuntimeError("kde grid length mismatch")
    return out


def d_hat(x):
    """Cheng-Hall d_hat for one sample (1-D array)."""
    n = x.size
    s = float(np.std(x, ddof=1))
    if not np.isfinite(s) or s <= 0:
        return float("nan")
    h0, h2 = H0_C * s * n ** -0.2, H2_C * s * n ** (-1.0 / 9.0)
    pad = 3.0 * max(h0, h2)
    _, c, delta = _binned(x, float(x.min()) - pad, float(x.max()) + pad)
    f = _kde(c, delta, h0, n, 0)
    f2 = _kde(c, delta, h2, n, 2)
    j = int(np.argmax(f))
    return abs(float(f2[j])) / float(f[j]) ** 3


def crit(table_q, dval):
    """Interpolated dip critical value. table_q: array (len(D_NODES),) of quantiles for one (n, alpha), node 0 = d 0;
    dval: d_hat values (array). Linear in d on [0, D_GRID[0]], linear in log d above, clipped at D_GRID[-1]."""
    dval = np.asarray(dval, float)
    lo = dval < D_GRID[0]
    ld = np.log(np.clip(dval, D_GRID[0], D_GRID[-1]))
    hi_part = np.interp(ld, np.log(D_GRID), table_q[1:])
    lo_part = table_q[0] + (table_q[1] - table_q[0]) * np.clip(dval, 0.0, D_GRID[0]) / D_GRID[0]
    return np.where(lo, lo_part, hi_part)


def d_hat_direct(x):
    """Unbinned reference computation of d_hat on the same 512-point grid (fidelity check)."""
    n = x.size
    s = float(np.std(x, ddof=1))
    h0, h2 = H0_C * s * n ** -0.2, H2_C * s * n ** (-1.0 / 9.0)
    pad = 3.0 * max(h0, h2)
    g = np.linspace(float(x.min()) - pad, float(x.max()) + pad, GRID)
    u0 = (g[:, None] - x[None, :]) / h0
    u2 = (g[:, None] - x[None, :]) / h2
    f = np.exp(-0.5 * u0 * u0).sum(1) / (n * h0 * math.sqrt(2 * math.pi))
    f2 = ((u2 * u2 - 1) * np.exp(-0.5 * u2 * u2)).sum(1) / (n * h2 ** 3 * math.sqrt(2 * math.pi))
    j = int(np.argmax(f))
    return abs(float(f2[j])) / float(f[j]) ** 3


def _self_test():
    for b in (1.5, 3.0, 20.0, 2000.0):
        d = d_beta(b)
        r = reference(d)
        assert r[0] == "beta" and abs(r[1] - b) / b < 1e-6, (b, r)
    for nu in (1.0, 5.0, 50.0):
        d = d_t(nu)
        r = reference(d)
        assert r[0] == "t" and abs(r[1] - nu) / nu < 1e-6, (nu, r)
    assert abs(d_beta(1e6) - TWO_PI) / TWO_PI < 1e-4 and abs(d_t(1e8) - TWO_PI) / TWO_PI < 1e-4
    rng = np.random.default_rng(1)  # DEV check only: d_hat of large normal samples near 2 pi
    vals = [d_hat(rng.standard_normal(20000)) for _ in range(20)]
    print("beta/t inversion OK; normal n=20000 d_hat median", round(float(np.median(vals)), 3), "(2 pi = 6.283)")
    vals = [d_hat(rng.beta(4, 8, 20000)) for _ in range(20)]
    print("beta(4,8) n=20000 d_hat median", round(float(np.median(vals)), 3))


if __name__ == "__main__":
    _self_test()
