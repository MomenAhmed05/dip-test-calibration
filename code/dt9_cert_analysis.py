"""DT9 frozen analysis (DT9-PREREG). Group-sequential, intersection-union certificate.

For each n and look l (chunk prefix 0..L_l-1 of EVERY shard; N_l = NSH * CH * L_l samples):
  K_j(c) = number of samples with block-j exceedance at cutoff c (LP failures count as exceedances).
  b = accept_count(N_l, alpha, D) = largest x with P(Bin(N_l, alpha) <= x) <= D   (exact; dt8_cert_analysis).
  alpha = 0.05: cutoff 0.7315 is certified at look l iff max_j K_j <= b.
  alpha = 0.01: fixed sequence 0.87 then 0.86; a cutoff is certified at look l iff it and every earlier one in
  the sequence pass at that look.
D = 0.01 / 12 (2 n x 2 alpha x 3 looks; intersection-union: a false all-blocks certificate needs the true
failing block to pass, so no block factor). Family error of all DT9 claims <= 0.01.
Stop rule (frozen): an n stops after the first look at which 0.7315 is certified (alpha = 0.05).
Usage: python dt9_cert_analysis.py root [look]   -> prints and writes deposits/dt9/cert_analysis_r<root>.json
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt8_cert_analysis import accept_count  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DEP = os.path.join(HERE, "..", "deposits", "dt9")
NSH, CH = 11, 500
LOOKS = {800: [36, 72, 108], 1600: [24, 48, 72]}      # chunks per shard at looks 1, 2, 3
CUTS = (0.7315, 0.86, 0.87)
SEQ = {0.05: [0.7315], 0.01: [0.87, 0.86]}
D = 0.01 / 12


def counts(root, n, L):
    """Exceedance counts K[j][ci] over chunk prefix 0..L-1 of all shards; None if any chunk is missing."""
    K = [[0] * len(CUTS) for _ in range(4)]
    N, fails, fb, rep = 0, 0, 0, 0.0
    for s in range(NSH):
        for k in range(L):
            p = os.path.join(DEP, f"cert_n{n}_r{root}_s{s}_c{k:04d}.json")
            if not os.path.exists(p):
                return None
            try:
                with open(p) as fh:
                    d = json.load(fh)
            except (OSError, ValueError) as err:
                sys.exit(f"cannot read {p}: {err}")
            for row, fl in zip(d["ex"], d["failed"], strict=True):
                N += 1
                fails += sum(fl)
                for j in range(4):
                    for ci in range(len(CUTS)):
                        K[j][ci] += int(row[j][ci])
            fb += d["stats"]["fallbacks"]
            rep = max(rep, d["stats"]["max_repair"])
    return {"N": N, "K": K, "failed_leaf_or_node": fails, "fallbacks": fb, "max_repair_raw": rep}


def look(root, n, li):
    L = LOOKS[n][li]
    c = counts(root, n, L)
    if c is None:
        return None
    out = {"n": n, "look": li + 1, "chunks_per_shard": L, **c, "D": D, "alpha": {}}
    for a, seq in SEQ.items():
        b = int(accept_count(c["N"], a, D))
        trail, cert = [], None
        for cut in seq:
            ci = CUTS.index(cut)
            Kj = [c["K"][j][ci] for j in range(4)]
            ok = max(Kj) <= b
            trail.append({"c": cut, "K": Kj, "pass": ok, "rate_edge": Kj[3] / c["N"]})
            if not ok:
                break
            cert = cut
        out["alpha"][str(a)] = {"b": b, "certified": cert, "trail": trail}
    return out


def main():
    root = int(sys.argv[1])
    res = {"root": root, "D": D, "rows": []}
    for n in LOOKS:
        for li in range(3):
            r = look(root, n, li)
            if r is None:
                break
            res["rows"].append(r)
            a5, a1 = r["alpha"]["0.05"], r["alpha"]["0.01"]
            print(f"n={n} look={li + 1} N={r['N']} 5%: b={a5['b']} cert={a5['certified']} "
                  f"K={a5['trail'][0]['K']} | 1%: b={a1['b']} cert={a1['certified']} "
                  f"trail={[(t['c'], t['K'], t['pass']) for t in a1['trail']]} fails={r['failed_leaf_or_node']} "
                  f"fb={r['fallbacks']} rep={r['max_repair_raw']:.1e}")
            if a5["certified"] == 0.7315:
                break
    p = os.path.join(DEP, f"cert_analysis_r{root}.json")
    try:
        with open(p, "w") as fh:
            json.dump(res, fh, indent=1)
    except OSError as err:
        sys.exit(f"cannot write {p}: {err}")


if __name__ == "__main__":
    main()
