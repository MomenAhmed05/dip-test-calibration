"""Figures for the dip-test paper.

Every plotted number is read from a filed deposit under deposits/ (paths below);
nothing is re-estimated here. The only random draw is the single illustrative
Brownian-bridge path in Figure 1 (fixed seed, illustration only).

Run:  python -I paper/figs/make_figs.py   (from the repository root)
"""
import json
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DEP = os.path.join(ROOT, "deposits")
OUT = os.path.dirname(os.path.abspath(__file__))

# Okabe-Ito palette (colour-blind safe)
C_UNIF = "#D55E00"   # vermillion: uniform calibration
C_LF = "#0072B2"     # blue: least-favourable constant
C_KUI = "#009E73"    # green: Kuiper/2
C_DKW = "#CC79A7"    # purple: DKW
C_GREY = "#6b6b6b"
C_GOLD = "#E69F00"
C_SKY = "#56B4E9"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "Times New Roman"],
    "mathtext.fontset": "stix",
    "font.size": 9.5,
    "axes.titlesize": 10,
    "axes.labelsize": 9.5,
    "legend.fontsize": 8.3,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.5,
    "savefig.dpi": 300,
})

# Certified brackets for the limiting uniform-law quantile qbar (DT8 Part A, DT13 section 6)
QBAR = {0.10: (0.4936, 0.5108), 0.05: (0.5379, 0.5558), 0.01: (0.6162, 0.6438)}
C_ALPHA = {0.10: 0.6676, 0.05: 0.73139526, 0.01: 0.85987427}  # 0.5 F_M^{-1}(1-a) (DT4-CERT, DT5-PREREG)


def load(rel):
    path = os.path.join(DEP, rel)
    try:
        with open(path, encoding="utf8") as fh:
            return json.load(fh)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"cannot read deposit {path}: {exc}") from exc


BF = load("dt4/bf_table_full.json")  # Balabdaoui and Filali (2011), Table 2
BF_X = np.array(BF["x"])
BF_F = np.array(BF["F_M"])


def F_M(x):
    x = np.asarray(x, dtype=float)
    return np.interp(x, BF_X, BF_F, left=0.0, right=1.0)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 1
def lower_hull(x, y):
    """Greatest convex minorant of the points (x_i, y_i), x increasing."""
    hull = []
    for i in range(len(x)):
        while len(hull) >= 2:
            i1, i2 = hull[-2], hull[-1]
            if (y[i2] - y[i1]) * (x[i] - x[i1]) >= (y[i] - y[i1]) * (x[i2] - x[i1]):
                hull.pop()
            else:
                break
        hull.append(i)
    return np.interp(x, x[hull], y[hull])


def fig1():
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.95), gridspec_kw={"width_ratios": [1, 1.25]})
    # (a) the two-step density W(w, eps); eps drawn wider than 0.005 so the top is visible
    w, eps = 0.10, 0.04
    lo, hi = (1 - w) / (1 - eps), w / eps
    a.fill_between([0, 1 - eps], 0, lo, color=C_SKY, alpha=0.35, lw=0)
    a.fill_between([1 - eps, 1], 0, hi, color=C_GOLD, alpha=0.55, lw=0)
    a.plot([0, 0, 1 - eps, 1 - eps, 1, 1], [0, lo, lo, hi, hi, 0], color="k", lw=1.1)
    a.text(0.45, lo * 0.45, "flat shoulder\nmass $p = 1-w$", ha="center", va="center", fontsize=9)
    a.annotate("thin top: mass $w$,\nwidth $\\varepsilon$", xy=(1 - eps / 2, hi * 0.85), xytext=(0.38, hi * 1.05),
               fontsize=9, ha="center", va="center",
               arrowprops={"arrowstyle": "->", "lw": 0.8, "color": "k", "shrinkB": 3})
    a.annotate("the upward jump\npins the mode here", xy=(1 - eps - 0.01, (lo + hi) / 2), xytext=(0.40, lo + 0.62),
               fontsize=8.3, ha="center", va="center", color=C_GREY,
               arrowprops={"arrowstyle": "->", "lw": 0.7, "color": C_GREY, "shrinkB": 4})
    a.set_xlim(-0.03, 1.05)
    a.set_ylim(0, hi * 1.22)
    a.set_xlabel("$x$")
    a.set_ylabel("density")
    a.set_title("(a) A null law $W(w,\\varepsilon)$ (not to scale)", loc="left")

    # (b) one Brownian-bridge path on the shoulder, its greatest convex minorant, and the best convex fit
    rng = np.random.default_rng(20261110)
    m = 2000
    t = np.linspace(0, 1, m + 1)
    z = np.concatenate([[0.0], np.cumsum(rng.standard_normal(m)) / np.sqrt(m)])
    br = z - t * z[-1]
    gcm = lower_hull(t, br)
    gap = br - gcm
    k = int(np.argmax(gap))
    G = gap[k]
    b.plot(t, br, color="k", lw=0.8, label="bridge $b$ on the shoulder")
    b.plot(t, gcm, color=C_LF, lw=1.4, label="greatest convex minorant")
    b.fill_between(t, gcm, gcm + G, color=C_LF, alpha=0.13, lw=0, label="minorant to minorant $+\\,M$")
    b.plot(t, gcm + G / 2, color=C_LF, lw=1.0, ls="--", label="best convex fit (minorant $+\\,M/2$)")
    b.annotate("", xy=(t[k], br[k]), xytext=(t[k], gcm[k]),
               arrowprops={"arrowstyle": "<->", "lw": 1.0, "color": C_UNIF})
    b.text(t[k] + 0.02, (br[k] + gcm[k]) / 2, "$M$", color=C_UNIF, fontsize=11, va="center")
    b.set_xlabel("position in the shoulder (probability scale)")
    b.set_title("(b) The shoulder seen at scale $n^{-1/2}$", loc="left")
    handles, labels = b.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=4, frameon=False,
               fontsize=7.4, handlelength=1.5, columnspacing=1.0)
    ymin, ymax = br.min(), br.max()
    top = max(ymax, float((gcm + G).max()))
    b.set_ylim(ymin - 0.08 * (top - ymin), top + 0.04 * (top - ymin))
    b.set_yticks([])
    b.spines["left"].set_visible(False)
    fig.tight_layout(w_pad=1.6, rect=(0, 0.09, 1, 1))
    save(fig, "fig1_mechanism")


# ---------------------------------------------------------------- Figure 2
def fig2():
    d = np.load(os.path.join(DEP, "dt1b", "dips_f0_n12800.npy")) * np.sqrt(12800)
    xs = np.linspace(0.15, 1.15, 600)
    H = np.searchsorted(np.sort(d), xs, side="right") / d.size
    fig, ax = plt.subplots(figsize=(4.9, 3.0))
    ax.plot(xs, H, color=C_UNIF, label="uniform law: $\\sqrt{n} D_n$, $n=12{,}800$")
    ax.plot(xs, F_M(2 * xs), color=C_LF, label="$M/2$ (worst shoulder)")
    ax.plot(xs, F_M(2 * xs / np.sqrt(0.9)), color=C_LF, ls=":", lw=1.3,
            label="$\\sqrt{0.9} M/2$ (shoulder mass 0.9)")
    q = 0.5465
    hq = np.interp(q, xs, H)
    fq = float(F_M(2 * q))
    ax.axvline(q, color=C_GREY, lw=0.7, ls="--")
    ax.annotate("", xy=(q + 0.012, fq), xytext=(q + 0.012, hq),
                arrowprops={"arrowstyle": "<->", "lw": 1.1, "color": "k"})
    ax.annotate("excess rejection $\\approx 0.25$\nat the uniform 5% cutoff", xy=(q + 0.016, (fq + hq) / 2),
                xytext=(0.78, 0.80), fontsize=8.3, va="center",
                arrowprops={"arrowstyle": "-", "lw": 0.6, "color": C_GREY})
    c = C_ALPHA[0.05]
    ax.axvline(c, color=C_LF, lw=0.7, ls="--")
    ax.text(c + 0.012, 0.42, "$c_{0.05}=0.7314$", color=C_LF, fontsize=8.3)
    ax.text(q - 0.012, 0.03, "$\\bar q_{0.05}\\approx 0.546$", color=C_GREY, fontsize=8.3, ha="right")
    ax.set_xlim(0.15, 1.15)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel("$x$ (scale $\\sqrt{n}\\times$ dip)")
    ax.set_ylabel("$P( \\cdot\\leq x)$")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.04), frameon=False)
    save(fig, "fig2_limit_laws")


# ---------------------------------------------------------------- Figure 3
def fig3():
    a = load("dt5/analysis.json")
    rows = [x for x in a["cells"] if x["law"][1] == "W(0.1,.005)" and x["method"] in ("uniform", "LF")]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.7), sharex=True)
    for ax, alpha in zip(axes, (0.05, 0.01), strict=True):
        for meth, col, lab in (("uniform", C_UNIF, "uniform calibration"), ("LF", C_LF, "corrected cutoff $c_\\alpha/\\sqrt{n}$")):
            r = sorted((x for x in rows if x["alpha"] == alpha and x["method"] == meth), key=lambda x: x["n"])
            n = np.array([x["n"] for x in r])
            p = np.array([x["rate"] for x in r]) * 100
            lo = np.array([x["ci"][0] for x in r]) * 100
            hi = np.array([x["ci"][1] for x in r]) * 100
            ax.fill_between(n, lo, hi, color=col, alpha=0.18, lw=0)
            ax.plot(n, p, "o-", color=col, ms=3.2, label=lab)
        qlo, qhi = QBAR[alpha]
        blo = (1 - F_M(2 * qhi / np.sqrt(0.9))) * 100
        bhi = (1 - F_M(2 * qlo / np.sqrt(0.9))) * 100
        ax.axhspan(blo, bhi, color=C_UNIF, alpha=0.12, lw=0)
        ax.text(60, bhi + (1.0 if alpha == 0.05 else 0.4), "asymptotic band (Corollary 1)", color=C_UNIF, fontsize=8)
        ax.axhline(alpha * 100, color="k", lw=0.7, ls="--")
        ax.text(13000, alpha * 100, " nominal", fontsize=8, va="center")
        ax.set_xscale("log")
        ax.set_xlim(40, 13000)
        ax.set_ylim(0, 30 if alpha == 0.05 else 13)
        ax.set_xlabel("sample size $n$")
        ax.set_title(f"$\\alpha = {int(alpha * 100)}\\%$", loc="left")
    axes[0].set_ylabel("rejection probability (%)")
    axes[0].legend(loc="center right", frameon=False)
    fig.tight_layout(w_pad=1.4)
    save(fig, "fig3_fixed_law")


# ---------------------------------------------------------------- Figure 4
def fig4():
    sm = load("dt13/size_map.json")["size_map"]
    band = {0.10: (40.45, 46.22), 0.05: (28.40, 33.20), 0.01: (12.67, 16.79)}
    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.75), sharey=False)
    for ax, alpha in zip(axes, (0.10, 0.05, 0.01), strict=True):
        r = sorted((x for x in sm if abs(x["alpha"] - alpha) < 1e-12 and x["U"] is not None), key=lambda x: x["n"])
        n = np.array([x["n"] for x in r], dtype=float)
        L = np.array([x["L"] for x in r]) * 100
        U = np.array([x["U"] for x in r]) * 100
        ax.vlines(n, L, U, color=C_UNIF, lw=2.2, alpha=0.75)
        ax.plot(n, L, "_", color=C_UNIF, ms=5)
        ax.plot(n, U, "_", color=C_UNIF, ms=5)
        xinf = 3200
        ax.vlines([xinf], band[alpha][0], band[alpha][1], color="k", lw=3.0)
        ax.text(xinf * 1.25, sum(band[alpha]) / 2, "$n\\to\\infty$", fontsize=8, va="center")
        ax.axhline(alpha * 100, color="k", lw=0.7, ls="--")
        ax.set_xscale("log")
        ax.set_xlim(3.2, 9000)
        ax.set_xticks([4, 10, 100, 1000])
        ax.set_xticklabels(["4", "10", "100", "1000"])
        ax.set_title(f"$\\alpha = {int(round(alpha * 100))}\\%$", loc="left")
        ax.set_xlabel("$n$")
        ax.set_ylim(0, {0.10: 80, 0.05: 70, 0.01: 55}[alpha])
    axes[0].set_ylabel("worst-case size $S_n(\\alpha)$ (%)")
    fig.tight_layout(w_pad=1.0)
    save(fig, "fig4_worst_size")


# ---------------------------------------------------------------- Figure 5
def fig5():
    ref = load("dt13/ref_brackets.json")["by_n"]
    small = np.genfromtxt(os.path.join(DEP, "dt14", "r1small_cutoffs.csv"), delimiter=",", names=True)
    n4 = {0.10: 0.2335, 0.05: 0.2420, 0.01: 0.2488}
    big = {0.05: [(n, 0.7315) for n in (50, 100, 200, 400, 800, 1600)],
           0.01: [(n, 0.86) for n in (50, 100, 200, 400, 800)] + [(1600, 0.87)]}
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.8), sharex=True)
    for ax, alpha in zip(axes, (0.05, 0.01), strict=True):
        ns = sorted(int(k) for k in ref)
        qn = [ref[str(k)][str(alpha)]["q_hat"] * np.sqrt(k) for k in ns]
        ax.plot(ns, qn, color=C_UNIF, lw=1.3, label="uniform-law quantile $\\sqrt{n} q_n$")
        qlo, qhi = QBAR[alpha]
        ax.axhspan(qlo, qhi, xmin=0.0, xmax=1.0, color=C_UNIF, alpha=0.10, lw=0)
        ax.axhline(C_ALPHA[alpha], color=C_LF, lw=1.0, ls="--", label="corrected constant $c_\\alpha$")
        sel = small[np.isclose(small["alpha"], alpha)]
        cert_n = [4] + list(sel["n"].astype(int)) + [x[0] for x in big[alpha]]
        cert_v = [n4[alpha] * 2] + list(sel["d_cert"] * np.sqrt(sel["n"])) + [x[1] for x in big[alpha]]
        ax.plot(cert_n, cert_v, "o", color=C_LF, ms=3.0, label="certified finite-$n$ cutoffs")
        ax.set_xscale("log")
        ax.set_xlim(3.4, 2200)
        ax.set_ylim(0.33, 1.0)
        ax.set_xlabel("sample size $n$")
        ax.set_title(f"$\\alpha = {int(alpha * 100)}\\%$", loc="left")
        ax.text(2100, (qlo + qhi) / 2, "$\\bar q_\\alpha$", color=C_UNIF, fontsize=8.5, va="center", ha="right")
    axes[0].set_ylabel("cutoff on the $\\sqrt{n}\\times$dip scale")
    axes[1].legend(loc="lower right", frameon=False, fontsize=7.8)
    fig.tight_layout(w_pad=1.4)
    save(fig, "fig5_cutoffs")


# ---------------------------------------------------------------- Figure 6
def fig6():
    a = load("dt5/analysis.json")
    meths = (("LF", C_LF, "corrected, 0.731"), ("kuiper2", C_KUI, "Kuiper/2, 0.874"),
             ("dkw", C_DKW, "DKW, 1.358"), ("uniform", C_GREY, "uniform calibration (size not controlled)"))
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.7), sharey=True)
    for ax, delta in zip(axes, ("1.5", "2.0"), strict=True):
        for meth, col, lab in meths:
            r = sorted((x for x in a["cells"] if x["law"][0] == "sym" and str(x["law"][1]) == delta
                        and x["alpha"] == 0.05 and x["method"] == meth), key=lambda x: x["n"])
            n = [x["n"] for x in r]
            p = [100 * x["rate"] for x in r]
            if meth == "uniform":
                ax.plot(n, p, ls=":", color=col, lw=1.1, label=lab)
            else:
                ax.plot(n, p, "o-", color=col, ms=3.0, label=lab)
        ax.axhline(80, color="k", lw=0.6, ls="--")
        ax.set_xscale("log")
        ax.set_xlabel("sample size $n$")
        ax.set_title(f"$\\frac{{1}}{{2}}N(-{delta},1)+\\frac{{1}}{{2}}N({delta},1)$, $\\alpha=5\\%$", loc="left")
    axes[0].set_ylabel("power (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=7.8,
               bbox_to_anchor=(0.5, -0.02), columnspacing=1.4, handlelength=1.8)
    fig.tight_layout(w_pad=1.2, rect=(0, 0.08, 1, 1))
    save(fig, "fig6_power")


# ---------------------------------------------------------------- Figure 7
def fig7():
    """Geometry of Section 5: probability scale, seam split, and the two ranges of turning points."""
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.0, 2.7), gridspec_kw={"width_ratios": [1, 1.15, 1.1]})
    # (a) a unimodal CDF: convex before the mode, concave after; u = F(m)
    x = np.linspace(0, 1, 801)
    dens = x * (1 - x) ** 3
    cdf = np.concatenate([[0.0], np.cumsum((dens[1:] + dens[:-1]) / 2 * np.diff(x))])
    cdf /= cdf[-1]
    m = 0.25
    um = float(np.interp(m, x, cdf))
    a.fill_between(x[x <= m], 0, cdf[x <= m], color=C_SKY, alpha=0.30, lw=0)
    a.fill_between(x[x >= m], 0, cdf[x >= m], color=C_GOLD, alpha=0.30, lw=0)
    a.plot(x, cdf, color="k", lw=1.3)
    a.plot([m, m], [0, um], color=C_GREY, lw=0.8, ls="--")
    a.plot([0, m], [um, um], color=C_GREY, lw=0.8, ls="--")
    a.text(0.03, um + 0.05, "$u=F(m)$", fontsize=8.5)
    a.text(m + 0.02, 0.04, "mode $m$", fontsize=8.5)
    a.text(0.07, 0.62, "convex\nbefore $m$", fontsize=8, color=C_LF)
    a.text(0.55, 0.45, "concave\nafter $m$", fontsize=8, color=C_UNIF)
    a.set_xlim(0, 1)
    a.set_ylim(0, 1.05)
    a.set_xlabel("$x$")
    a.set_ylabel("$F(x)$")
    a.set_title("(a) The mode on the probability scale", loc="left", fontsize=9)

    # (b) one bridge cut at u: removing the broken line leaves two independent bridges
    rng = np.random.default_rng(20261111)
    k = 1500
    t = np.linspace(0, 1, k + 1)
    zz = np.concatenate([[0.0], np.cumsum(rng.standard_normal(k)) / np.sqrt(k)])
    br = zz - t * zz[-1]
    u = 0.7
    bu = float(np.interp(u, t, br))
    b.axvspan(0, u, color=C_SKY, alpha=0.18, lw=0)
    b.axvspan(u, 1, color=C_GOLD, alpha=0.22, lw=0)
    b.plot(t, br, color="k", lw=0.7)
    b.plot([0, u, 1], [0, bu, 0], color=C_LF, lw=1.2, ls="--")
    b.axvline(u, color=C_GREY, lw=0.8)
    span = br.max() - br.min()
    ytop = br.max() + 0.45 * span
    b.text(0.28, ytop - 0.05 * span, "$\\sqrt{u}\\,b_1$\nconvex side", ha="center", va="top", fontsize=8)
    b.text((1 + u) / 2 + 0.01, ytop - 0.05 * span, "$\\sqrt{1-u}\\,r_2$\nconcave\nside", ha="center", va="top",
           fontsize=7.6)
    b.set_ylim(br.min() - 0.08 * span, ytop)
    b.set_xlim(0, 1)
    b.set_yticks([])
    b.spines["left"].set_visible(False)
    b.set_xlabel("probability scale; cut at $u=0.7$")
    b.set_title("(b) A bridge split at the turning point", loc="left", fontsize=9)

    # (c) the two ranges of turning points at the 5 percent level, on the scale eps = sqrt(1 - u)
    near_hi, mid_lo, eps_max = 0.169, 0.130, 1 / np.sqrt(2)
    c.barh(1.0, near_hi, left=0, height=0.34, color=C_LF, alpha=0.55)
    c.barh(0.45, eps_max - mid_lo, left=mid_lo, height=0.34, color=C_GOLD, alpha=0.75)
    c.text(near_hi + 0.02, 1.0, "near the edge:\nwritten bound", va="center", fontsize=7.8)
    c.text(mid_lo + 0.03, 0.45, "middle: paired simulation", va="center", fontsize=7.8)
    c.axvline(0, color="k", lw=1.6)
    c.text(0.02, 1.62, "edge $u=1$: distance has the law of $M/2$", va="center", fontsize=7.6)
    c.text(0.36, 0.13, "$u<\\frac{1}{2}$: by reflection $u\\mapsto1-u$", va="center", ha="center", fontsize=7.4,
           color=C_GREY)
    c.set_xlim(0, 0.75)
    c.set_ylim(-0.05, 1.85)
    c.set_yticks([])
    c.spines["left"].set_visible(False)
    c.set_xlabel("distance from the edge, $\\varepsilon=\\sqrt{1-u}$", labelpad=9)
    c.set_xticks([0, 0.2, 0.4, eps_max])
    c.set_xticklabels(["0", "0.2", "0.4", "$1/\\sqrt{2}$"])
    c.set_title("(c) Turning points covered (5% level)", loc="left", fontsize=9)
    fig.tight_layout(w_pad=1.0)
    save(fig, "fig7_dependency")


def _old_fig7_unused():
    fig, ax = plt.subplots(figsize=(6.8, 3.9))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 62)
    ax.axis("off")

    tags = {"#eeeeee": "P", "#e8f1fa": "W", "#fdf0e0": "C", "#fff7cc": "W+C"}

    def box(x, y, w, h, text, fc, ec="k", fs=8.2):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                    fc=fc, ec=ec, lw=0.8))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, linespacing=1.25)
        ax.text(x + 0.5, y + h + 0.9, tags[fc], ha="left", va="bottom", fontsize=6.5, fontweight="bold",
                color=C_GREY)
        return (x, y, w, h)

    def arrow(p, q):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=9, lw=0.8, color="k",
                                     shrinkA=2, shrinkB=2))

    proven = "#e8f1fa"
    cert = "#fdf0e0"
    pub = "#eeeeee"
    b1 = box(2, 47, 26, 11, "Published inputs\nlaw of $M$; bridge, excursion\nand range laws; H&H Thm 3", pub)
    b2 = box(37, 47, 26, 11, "Domination (Theorem 5)\nworst law $\\to$ uniform\nsample with a fixed\nturning point", proven, fs=7.8)
    b3 = box(72, 47, 26, 11, "Seam decoupling\n(Lemma 5)", proven)
    b4 = box(37, 27, 26, 11, "Brownian reduction\n(Theorem 6)\n$\\limsup_n \\; \\sup_F P \\; \\leq \\; \\sup_u g(u,t)$", proven)
    b5 = box(2, 27, 26, 11, "Near the edge\nfour written lemmas\n(Lean-checked cores)\n+ interval arithmetic", proven, fs=7.8)
    b6 = box(72, 27, 26, 11, "Middle range\npaired Monte Carlo,\n33,000 paths, sign test\n+ reflection $u\\mapsto 1-u$", cert, fs=7.8)
    b7 = box(8, 6, 84, 11, "The edge is worst at $t'$: $g(u,t') < g(1,t') = 1-F_M(2t')$ for all $0<u<1$\n"
                           "Theorem 7: $\\liminf_n S_n$ and $\\limsup_n S_n$ both lie in $[1-F_M(2q_+), 1-F_M(2t')]$",
             "#fff7cc", fs=8.4)
    arrow((b2[0] + b2[2] / 2, b2[1] - 0.6), (b4[0] + b4[2] / 2, b4[1] + b4[3] + 0.6))
    arrow((b3[0] + b3[2] / 2, b3[1] - 0.6), (b4[0] + b4[2] - 3, b4[1] + b4[3] + 0.6))
    arrow((b1[0] + b1[2] / 2, b1[1] - 0.6), (b5[0] + b5[2] / 2, b5[1] + b5[3] + 0.6))
    arrow((b1[0] + b1[2] - 2, b1[1] - 0.6), (b4[0] + 8, b4[1] + b4[3] + 0.6))
    arrow((b5[0] + b5[2] / 2, b5[1] - 0.6), (b7[0] + 10, b7[1] + b7[3] + 0.6))
    arrow((b4[0] + b4[2] / 2, b4[1] - 0.6), (b7[0] + b7[2] / 2, b7[1] + b7[3] + 0.6))
    arrow((b6[0] + b6[2] / 2, b6[1] - 0.6), (b7[0] + b7[2] - 10, b7[1] + b7[3] + 0.6))
    # legend
    for i, (fc, lab) in enumerate(((pub, "P: published, cited"), (proven, "W: written proof in this paper"),
                                    (cert, "C: certified computation"))):
        ax.add_patch(FancyBboxPatch((2 + 34 * i, 0.3), 3.2, 2.2, boxstyle="round,pad=0.2", fc=fc, ec="k", lw=0.6))
        ax.text(6.5 + 34 * i, 1.4, lab, fontsize=7.8, va="center")
    save(fig, "fig7_dependency")


if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3", "4", "5", "6", "7"]
    for w in which:
        globals()["fig" + w]()
    print("done", which)
