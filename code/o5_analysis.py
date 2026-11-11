"""O5-ML frozen analyser (DT5ML-PREREG §1, rules M-R1..M-R3; written before the CONFIRM run; no RNG).

Reads deposits/o5ml[_quick]/cells/*.npz from code/o5_ml.py and writes analysis.json + O5ML-TABLES.md.
Per cell: b = #(uniform splits, LF does not), c = #(LF splits, uniform does not),
Delta = (b - c) / R with the Newcombe (1998, method 10) paired score interval,
exact one-sided McNemar p = P(Bin(b + c, 1/2) >= b) (p = 1 when b + c = 0),
Clopper-Pearson 95% per arm. M-R1 family: the 10 non-control continuous null laws
(N1, N2 x 4, N2-pin x 3, R0 x 2) x (UniDip a=.05/.01 x 2 n, SkinnyDip a=.05 x 2 n, DipMeans n=500) = 70 cells;
declared when p < 0.05 / 70 and Delta > 0. U, G, atom laws and price laws are reported only.
Run: python code/o5_analysis.py [--quick]   (any Python with numpy/scipy; imports o5_ml for the law list)
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import beta, binom, norm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import o5_ml  # noqa: E402

FAMILY_LAWS = ("N1", "N2_s0.0005", "N2_s0.002", "N2_s0.01", "N2_s0.05",
               "N2pin_k0.85", "N2pin_k1.4", "N2pin_k2.5", "R0_s0.002", "R0_s0.05")


def cp(k, n):
    lo = float(beta.ppf(0.025, k, n - k + 1)) if k > 0 else 0.0
    hi = float(beta.ppf(0.975, k + 1, n - k)) if k < n else 1.0
    return lo, hi


def wilson(k, n, z):
    p = k / n
    c = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    m = p + z * z / (2 * n)
    d = 1 + z * z / n
    return (m - c) / d, (m + c) / d


def newcombe_paired(a, b, c, d, z=norm.ppf(0.975)):
    """Newcombe (1998) method 10 for p1 - p2 with paired 2x2 counts
    a = both, b = arm1 only, c = arm2 only, d = neither."""
    n = a + b + c + d
    p1, p2 = (a + b) / n, (a + c) / n
    l1, u1 = wilson(a + b, n, z)
    l2, u2 = wilson(a + c, n, z)
    denom = (a + b) * (c + d) * (a + c) * (b + d)
    if denom == 0:
        phi = 0.0
    else:
        num = a * d - b * c
        num = max(num - n / 2, 0) if num > 0 else num
        phi = num / math.sqrt(denom)
    dd = p1 - p2
    lo = dd - math.sqrt(max((p1 - l1) ** 2 - 2 * phi * (p1 - l1) * (u2 - p2) + (u2 - p2) ** 2, 0))
    hi = dd + math.sqrt(max((u1 - p1) ** 2 - 2 * phi * (u1 - p1) * (p2 - l2) + (p2 - l2) ** 2, 0))
    return lo, hi


def main():
    quick = "--quick" in sys.argv
    out = Path(__file__).resolve().parents[1] / "deposits" / ("o5ml_quick" if quick else "o5ml")
    cells = o5_ml.cells()
    missing = [c for c in cells if not (out / "cells" / o5_ml.cell_path(c).name).exists()]
    if missing:
        raise SystemExit(f"{len(missing)} cells missing in {out}")
    fam = [c for c in cells if o5_ml.LAWS[c[1]][0] in FAMILY_LAWS]
    K = len(fam)
    rows = []
    for c in cells:
        try:
            z = np.load(out / "cells" / o5_ml.cell_path(c).name)
        except OSError as e:
            raise SystemExit(f"cannot read cell {c}: {e}") from e
        u, lf = z["k_U"] > 1, z["k_LF"] > 1
        R = int(u.size)
        a, b, cc = int(np.sum(u & lf)), int(np.sum(u & ~lf)), int(np.sum(~u & lf))
        d = R - a - b - cc
        p = float(binom.sf(b - 1, b + cc, 0.5)) if b + cc > 0 else 1.0
        name, role = o5_ml.LAWS[c[1]][0], o5_ml.LAWS[c[1]][1]
        in_fam = name in FAMILY_LAWS
        rows.append({"method": c[0], "law": name, "role": role, "n": c[2], "alpha": c[3], "R": R,
                     "split_U": (a + b) / R, "split_U_ci": cp(a + b, R), "split_LF": (a + cc) / R,
                     "split_LF_ci": cp(a + cc, R), "U_only": b / R, "LF_only": cc / R, "delta": (b - cc) / R,
                     "delta_ci": newcombe_paired(a, b, cc, d), "mcnemar_p": p, "in_family": in_fam,
                     "M_R1_declared": bool(in_fam and p < 0.05 / K and b > cc),
                     "root_rej_U": float(np.mean(z["root_pU"] < c[3])),
                     "root_rej_LF": float(np.mean(z["root_pLF"] < c[3])),
                     "noise_U": float(z["noise_U"].mean()), "noise_LF": float(z["noise_LF"].mean()),
                     "minn_U_median": float(np.median(z["minn_U"]))})
    (out / "analysis.json").write_text(json.dumps({"quick": quick, "K": K, "rows": rows}, indent=1))
    lines = [f"# O5-ML tables ({'DEV quick' if quick else 'CONFIRM'}; code/o5_analysis.py; M-R1 family K = {K})", "",
             "| method | law | role | n | alpha | split U % [CP] | split LF % [CP] | Delta pt [Newcombe] | McNemar p | M-R1 |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['method']} | {r['law']} | {r['role']} | {r['n']} | {r['alpha']} | "
                     f"{100 * r['split_U']:.1f} [{100 * r['split_U_ci'][0]:.1f}, {100 * r['split_U_ci'][1]:.1f}] | "
                     f"{100 * r['split_LF']:.1f} [{100 * r['split_LF_ci'][0]:.1f}, {100 * r['split_LF_ci'][1]:.1f}] | "
                     f"{100 * r['delta']:.1f} [{100 * r['delta_ci'][0]:.1f}, {100 * r['delta_ci'][1]:.1f}] | "
                     f"{r['mcnemar_p']:.2g} | {'YES' if r['M_R1_declared'] else ('no' if r['in_family'] else '-')} |")
    (out / "O5ML-TABLES.md").write_text("\n".join(lines) + "\n")
    print(f"K = {K}; declared {sum(r['M_R1_declared'] for r in rows)}; saved {out / 'analysis.json'}")


if __name__ == "__main__":
    main()
