"""Step 14 (DT13) Part S: upper bound on the finite-n SIZE of the uniform-calibrated dip test via Corollary U1'
(DT7-UNIF): sup_{F in U_c} P_F(D_n >= q) <= max_j P(Sbar_{n,I_j} >= q) for the 4 blocks I_j of [1/2, 1].

Cutoffs: the LOWER bracket end X_(l) of the uniform quantile for each alpha (deposits/dt13/ref_brackets.json,
DT13 Part W reference; X_(l) <= q_n(1 - alpha) on the bracket event), raw dip scale.
Conservative indicator (dt9_block_cert.exceeds certified primal bound): a sample counts as
non-exceeding for block j only if a certified upper bound shows Sbar < X_(l) strictly (threshold X_(l) - 1e-12);
LP failure at a leaf counts as an exceedance. Nested cutoffs: the smallest cutoff (alpha = 0.10) is tested first.
Seeds: chunk k of shard s uses SeedSequence([root, n, s, k]); CH samples per chunk; one JSON per chunk;
reruns skip existing chunks (resume, same root).
Usage: python dt13_u1_cert.py n chunks shard nshards root tag [CH]
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt9_block_cert import STATS, exceeds  # noqa: E402
from sup_pinned import pieces  # noqa: E402

BLOCKS = [(0.5, 0.9), (0.9, 0.97), (0.97, 0.995), (0.995, 1.0)]
ALPHAS = (0.10, 0.05, 0.01)          # increasing cutoffs


def cutoffs(n, ref_dir):
    p = os.path.join(ref_dir, "ref_brackets.json")
    try:
        with open(p) as fh:
            ref = json.load(fh)
    except (OSError, ValueError) as err:
        sys.exit(f"cannot read {p}: {err}")
    return [float(ref["by_n"][str(n)][str(a)]["X_l"]) for a in ALPHAS]


def sample_row(x, cuts):
    row = []
    for lo_u, hi_u in BLOCKS:
        P = pieces(x, lo_u, hi_u)
        ex, prev = [], True
        for c in cuts:
            if not prev:
                ex.append(False)
                continue
            e, _f = exceeds(x, P, c - 1e-12)
            ex.append(bool(e))
            prev = e
        row.append(ex)
    return row


def main():
    n, chunks, shard, ns = (int(a) for a in sys.argv[1:5])
    root, tag = int(sys.argv[5]), sys.argv[6]
    CH = int(sys.argv[7]) if len(sys.argv) > 7 else 200
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits")
    ref_dir = os.path.join(base, "dt13_quick" if tag.startswith("quick") else "dt13")
    out = os.path.join(base, "dt13", "u1")
    try:
        os.makedirs(out, exist_ok=True)
    except OSError as err:
        sys.exit(f"cannot create {out}: {err}")
    cuts = cutoffs(n, ref_dir)
    for k in range(chunks):
        path = os.path.join(out, f"{tag}_n{n}_r{root}_s{shard}_c{k:04d}.json")
        if os.path.exists(path):
            continue
        rng = np.random.default_rng(np.random.SeedSequence([root, n, shard, k]))
        t0 = time.time()
        for key in ("lps", "fallbacks"):
            STATS[key] = 0
        STATS["max_repair"] = 0.0
        R = [sample_row(np.sort(rng.uniform(size=n)), cuts) for _ in range(CH)]
        res = {"n": n, "root": root, "shard": shard, "chunk": k, "CH": CH, "blocks": BLOCKS, "alphas": ALPHAS,
               "cuts": cuts, "ex": R, "stats": dict(STATS), "seconds": time.time() - t0}
        tmp = path + ".tmp"
        try:
            with open(tmp, "w") as fh:
                json.dump(res, fh)
            os.replace(tmp, path)
        except OSError as err:
            sys.exit(f"cannot write {path}: {err}")
        print(f"n={n} shard={shard} chunk={k} {time.time() - t0:.1f}s lps={STATS['lps']} fb={STATS['fallbacks']}",
              flush=True)


if __name__ == "__main__":
    main()
