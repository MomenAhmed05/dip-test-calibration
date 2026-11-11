"""DT9 Part R analysis: DT8 Part B rule applied to the replayed (certified-numerics) indicators.

Rule (unchanged from DT8-PREREG): d = 0.01/48; for each (n, alpha) test the DT8 grid from the largest cutoff
down, stop at the first failure; a cutoff passes iff all 4 block counts <= accept_count(N, alpha, d).
Identity check: EVERY replayed edge value V must equal DT8's stored V of the same sample (|diff| < 1e-9),
and the sample counts must match per shard.
A DT8 claim is KEPT iff the replay certifies a cutoff <= the DT8 certified cutoff; otherwise it is WITHDRAWN
(replaced by the replay result). No claim can be strengthened by the replay.
Usage: python dt9_replay_analysis.py   -> deposits/dt9/replay_analysis.json
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt8_cert_analysis import accept_count  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DEP = os.path.join(HERE, "..", "deposits")
NS = {50: 100000, 100: 100000, 200: 100000, 400: 60000, 800: 30000, 1600: 15000}
GRID = {0.05: [0.75, 0.74, 0.735, 0.7315], 0.01: [0.90, 0.88, 0.87, 0.86]}
CUTS = (0.7315, 0.735, 0.74, 0.75, 0.86, 0.87, 0.88, 0.90)
D8 = 0.01 / 48
NSH = 11


def _load(p):
    try:
        with open(p) as fh:
            return json.load(fh)
    except (OSError, ValueError) as err:
        sys.exit(f"cannot read {p}: {err}")


def main():
    dt8 = _load(os.path.join(DEP, "dt8blk", "cert_analysis_r2026110711.json"))
    dt8c = {(r["n"], r["alpha"]): r["certified_cutoff"] for r in dt8["rows"]}
    rows, ident_bad = [], 0
    for n, R in NS.items():
        K = [[0] * len(CUTS) for _ in range(4)]
        N, fails, rep = 0, 0, 0.0
        for s in range(NSH):
            d = _load(os.path.join(DEP, "dt9", f"replay_n{n}_s{s}.json"))
            old = _load(os.path.join(DEP, "dt8blk", f"cert_n{n}_r2026110711_s{s}.json"))
            if len(d["V"]) != len(old["V"]):
                ident_bad += 10 ** 9
            ident_bad += sum(abs(a - b) > 1e-9 for a, b in zip(d["V"], old["V"], strict=False))
            for row, fl in zip(d["ex"], d["failed"], strict=True):
                N += 1
                fails += sum(fl)
                for j in range(4):
                    for ci in range(len(CUTS)):
                        K[j][ci] += int(row[j][ci])
            rep = max(rep, d["stats"]["max_repair"])
        assert N == R, (n, N, R)
        for a, grid in GRID.items():
            b = int(accept_count(N, a, D8))
            cert, trail = None, []
            for c in grid:
                Kj = [K[j][CUTS.index(c)] for j in range(4)]
                ok = max(Kj) <= b
                trail.append({"c": c, "K": Kj, "pass": ok})
                if not ok:
                    break
                cert = c
            old_c = dt8c[(n, a)]
            keep = cert is not None and old_c is not None and cert <= old_c
            rows.append({"n": n, "alpha": a, "N": N, "b": b, "replay_certified": cert, "dt8_certified": old_c,
                         "verdict": "KEPT" if keep else "WITHDRAWN", "trail": trail, "failed": fails,
                         "max_repair_raw": rep})
            print(f"n={n} a={a} b={b} replay={cert} dt8={old_c} {'KEPT' if keep else 'WITHDRAWN'} "
                  f"trail={[(t['c'], max(t['K']), t['pass']) for t in trail]} fails={fails} rep={rep:.1e}")
    print(f"identity mismatches: {ident_bad}")
    p = os.path.join(DEP, "dt9", "replay_analysis.json")
    try:
        with open(p, "w") as fh:
            json.dump({"d": D8, "identity_mismatches": ident_bad, "rows": rows}, fh, indent=1)
    except OSError as err:
        sys.exit(f"cannot write {p}: {err}")


if __name__ == "__main__":
    main()
