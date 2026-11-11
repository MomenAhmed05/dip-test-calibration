"""Summarise DT8 block-envelope runs: per n and block, P(sqrt(n) Sbar_I > c) on a cutoff grid, the edge
law P(sqrt(n) V > c), and the smallest grid cutoff whose worst block is <= alpha (point estimate).
Usage: python dt8_block_summary.py root [n ...]
"""
import glob
import json
import os
import sys

import numpy as np

CS = [0.70, 0.72, 0.7315, 0.735, 0.74, 0.745, 0.75, 0.76, 0.78, 0.80, 0.84, 0.86, 0.87, 0.88, 0.90]


def load(root, n):
    d = os.path.join(os.path.dirname(__file__), "..", "deposits", "dt8blk")
    S, V, fb, blocks = [], [], 0, None
    for f in sorted(glob.glob(os.path.join(d, f"*_n{n}_r{root}_s*.json"))):
        with open(f) as fh:
            r = json.load(fh)
        S += r["S"]
        V += r["V"]
        fb += r.get("fallbacks", 0)
        blocks = r["blocks"]
    return np.array(S), np.array(V), fb, blocks


def main():
    root = int(sys.argv[1])
    for n in sys.argv[2:]:
        S, V, fb, blocks = load(root, int(n))
        if len(S) == 0:
            continue
        print(f"n={n} N={len(S)} fallbacks={fb} blocks={blocks}")
        print("   c      edgeV  " + "  ".join(f"B{j}" + " " * 4 for j in range(len(blocks))) + " worst")
        for c in CS:
            pv = (V > c).mean()
            pb = (S > c).mean(axis=0)
            print(f"  {c:.4f}  {pv:.4f}  " + "  ".join(f"{p:.4f}" for p in pb) + f"  {pb.max():.4f}")
        for a in (0.05, 0.01):
            ok = [c for c in CS if (S > c).mean(axis=0).max() <= a]
            qv = float(np.quantile(V, 1 - a))
            qs = [float(np.quantile(S[:, j], 1 - a)) for j in range(S.shape[1])]
            print(f"  alpha={a}: edge quantile {qv:.4f}; block quantiles " +
                  ", ".join(f"{q:.4f}" for q in qs) + f"; smallest grid c with worst<=alpha: {ok[0] if ok else None}")


if __name__ == "__main__":
    main()
