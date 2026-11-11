"""DT14-ODDMAX numerical check of Lemma M: D(F_n) <= floor(n/2)/(2n).

DEV only (root 2026110835, CONVENTIONS); never a claim by itself.
Usage: python dt14_oddmax_check.py [root]
Prints, per n: M_n, the dip of two tight clusters, and a hill-climb maximum;
then the number of violations of D <= M_n (expected 0).
"""
import sys

import diptest
import numpy as np

ROOT = int(sys.argv[1]) if len(sys.argv) > 1 else 2026110835
NS = [5, 7, 8, 9, 10, 11, 13, 15, 20, 25, 49]


def main():
    rng = np.random.default_rng(ROOT)
    viol = 0
    print("n, M_n, two-cluster dip, hill-climb max")
    for n in NS:
        k = n // 2
        m_n = k / (2 * n)
        clus = 0.0
        for delta in (1e-2, 1e-4, 1e-6):
            x = np.concatenate([rng.uniform(0, delta, k), 1 + rng.uniform(0, delta, n - k)])
            d = diptest.dipstat(x)
            clus = max(clus, d)
            viol += d > m_n + 1e-12
        best = 0.0
        for _ in range(300):
            x = np.sort(rng.standard_normal(n) * rng.uniform(0.1, 3)
                        + rng.choice([0, rng.uniform(0, 6)], n))
            d = diptest.dipstat(x)
            for _ in range(60):
                y = x + rng.standard_normal(n) * 0.2 * rng.uniform(0, 1) ** 3
                dy = diptest.dipstat(y)
                if dy >= d:
                    x, d = y, dy
            best = max(best, d)
            viol += d > m_n + 1e-12
        print(n, round(m_n, 6), round(clus, 6), round(best, 6))
    print("violations", viol)


if __name__ == "__main__":
    main()
