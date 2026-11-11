"""Step 14 (DT13) analysis (frozen with DT13-PREREG). No random numbers are drawn here.

Inputs: deposits/dt13/confirm_WL.json (Parts W, L), deposits/dt13/u1/{tag}_n*_r{rootS}_*.json (Part S),
deposits/dt13/n4_bound_r{rootN4}.json (Part N4), deposits/dt13/n4_exact.json (Part X), qbar analysis (Part Q).
Size enclosure per (n, alpha):
  L = max over confirmed laws at that n of the simultaneous CP lower limit at X_(u) (Parts W, L) and, at n = 4,
      the EXACT rejection probability of the named n = 4 laws (Part X, CERT);
  U = (grid n only) max over the 4 blocks of the one-sided CP upper limit of the Part S count at X_(l),
      confidence 1 - DELTA_S / (#grid n * 3 * 4); at n = 4 also Lemma N4 (Part N4); U = min of the two.
Validators (P-3 / section 12): V1 uniform anchor contains alpha everywhere; V2 at n = 4 every MC enclosure of a
law with an exact value contains it; V6 L <= U on the grid; V8 k_Xu <= k_q <= k_Xl in every cell.
Usage: python dt13_analysis.py rootS tag rootN4 -> deposits/dt13/size_map.json and a printed table.
"""
import glob
import json
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DEP = os.path.join(HERE, "..", "deposits", "dt13")
ALPHAS = (0.01, 0.05, 0.10)
GRID_S = (4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500, 750, 1000)   # = G_L
DELTA_S = 0.002


def _load(p):
    try:
        with open(p) as fh:
            return json.load(fh)
    except (OSError, ValueError) as err:
        raise SystemExit(f"cannot read {p}: {err}") from err


def part_s(root, tag):
    out = {}
    nb = 4
    conf = 1 - DELTA_S / (len(GRID_S) * len(ALPHAS) * nb)   # one-sided, 21 x 3 x 4 = 252 bounds
    for n in GRID_S:
        files = sorted(glob.glob(os.path.join(DEP, "u1", f"{tag}_n{n}_r{root}_s*_c*.json")))
        if not files:
            continue
        ex = np.concatenate([np.array(_load(f)["ex"], bool) for f in files], axis=0)   # (N, blocks, alphas .10,.05,.01)
        N = ex.shape[0]
        cnt = ex.sum(axis=0)
        row = {}
        for ja, a in enumerate((0.10, 0.05, 0.01)):
            ups = [1.0 if k == N else float(stats.beta.ppf(conf, k + 1, N - k)) for k in cnt[:, ja]]
            row[str(a)] = {"N": N, "counts": cnt[:, ja].tolist(), "U": max(ups), "argmax_block": int(np.argmax(ups))}
        out[n] = row
    return out


def main():
    rootS, tag, rootN4 = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
    wl = _load(os.path.join(DEP, "confirm_WL.json"))
    ex4 = _load(os.path.join(DEP, "n4_exact.json"))
    n4b = _load(os.path.join(DEP, f"n4_bound_r{rootN4}.json"))
    S = part_s(rootS, tag)
    table = wl["table"]
    # V8
    v8 = all(r["k_Xu"] <= r["k_q"] <= r["k_Xl"] for r in table)
    # V1
    v1 = all(r["lower"] <= r["alpha"] <= r["upper"] for r in table if r["part"] == "W" and r["f"] == 0)
    # V2: n = 4 exact values (Part X) inside the n = 4 MC enclosures of the same law
    exact4 = {(e["law"], e["alpha"]): e for e in ex4["confirm_laws"]}
    v2_rows = []
    for r in table:
        if r["n"] != 4:
            continue
        key = (r["law"], r["alpha"])
        if key in exact4:
            lo, hi = exact4[key]["P"]
            v2_rows.append({"law": r["law"], "alpha": r["alpha"], "exact": [lo, hi], "mc": [r["lower"], r["upper"]],
                            "ok": bool(r["lower"] <= lo and hi <= r["upper"])})
    v2 = bool(v2_rows) and all(x["ok"] for x in v2_rows)
    smap = {}
    for r in table:
        if r["f"] == 0:
            continue
        key = (r["n"], r["alpha"])
        cur = smap.get(key)
        if cur is None or r["lower"] > cur["L"]:
            smap[key] = {"n": r["n"], "alpha": r["alpha"], "L": r["lower"], "L_law": r["law"], "L_tier": "MC (CP, simultaneous)"}
    for e in ex4["confirm_laws"]:
        key = (4, e["alpha"])
        if e["P"][0] > smap[key]["L"]:
            smap[key].update({"L": e["P"][0], "L_law": e["law"], "L_tier": "EXACT (Part X)"})
    for (n, a), v in smap.items():
        v["U"], v["U_src"] = None, None
        if n in S:
            v["U"], v["U_src"] = S[n][str(a)]["U"], f"U1' block {S[n][str(a)]['argmax_block']}"
        if n == 4:
            u4 = n4b["alphas"][str(a)]["U"]
            if v["U"] is None or u4 < v["U"]:
                v["U"], v["U_src"] = u4, "Lemma N4"
    v6 = all(v["U"] is None or v["L"] <= v["U"] for v in smap.values())
    res = {"validators": {"V1_anchor": v1, "V2_exact_n4": v2, "V6_L_le_U": v6, "V8_counts": v8},
           "v2_rows": v2_rows, "S": {str(k): v for k, v in S.items()},
           "size_map": [smap[k] for k in sorted(smap)]}
    try:
        with open(os.path.join(DEP, "size_map.json"), "w") as fh:
            json.dump(res, fh, indent=1)
    except OSError as err:
        sys.exit(f"cannot write size_map.json: {err}")
    print("validators:", res["validators"])
    for v in res["size_map"]:
        if v["U"] is not None:
            print(f"n={v['n']:5d} a={v['alpha']:.2f}: size in [{v['L']:.4f}, {v['U']:.4f}]  L by {v['L_law']} ({v['L_tier']}); U by {v['U_src']}")


if __name__ == "__main__":
    main()
