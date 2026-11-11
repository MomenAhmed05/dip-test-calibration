"""Independent exact polytope volume (F-8 checker for exact_poly.py; Step 14 / DT13 Part X).

Different algorithm: Lasserre's (1983) facet recursion, in exact rationals, with no vertex enumeration and
no triangulation:
  vol_d(P) = (1/d) * sum_i  b_i / |a_ij| * vol_{d-1}( proj_j F_i ),   P = {x : A x <= b},
where F_i = P cap {a_i x = b_i}, j is any coordinate with a_ij != 0, and proj_j F_i is the (d-1)-dimensional
polytope obtained by eliminating x_j through the equality. (Euler's identity for homogeneous functions; the
projection divides the facet area by |a_ij| / ||a_i||, which cancels the ||a_i|| of the original formula.)
Duplicate constraints (positive multiples) are merged first; constraints that become 0 <= c after elimination
are kept as feasibility checks. d = 1: length of an interval.
Usage: python exact_poly_check.py  -> compares with exact_poly.volume on the n = 4 regions and on seeded faults.
"""
from fractions import Fraction


def _normalise(A, b):
    """Merge positive multiples; drop 0 x <= c rows if c >= 0; return None if some 0 x <= c with c < 0."""
    seen = {}
    for row, bi in zip(A, b, strict=True):
        piv = next((abs(v) for v in row if v != 0), None)
        if piv is None:
            if bi < 0:
                return None
            continue
        key = tuple(v / piv for v in row)
        val = bi / piv
        if key not in seen or val < seen[key]:
            seen[key] = val
    return [list(k) for k in seen], list(seen.values())


def volume_lasserre(A, b):
    A = [[Fraction(v) for v in row] for row in A]
    b = [Fraction(v) for v in b]
    nb = _normalise(A, b)
    if nb is None:
        return Fraction(0)
    A, b = nb
    if not A:
        raise ValueError("unbounded (no constraints)")
    d = len(A[0])
    if d == 1:
        lo, hi = None, None
        for (a,), bi in zip(A, b, strict=True):
            if a > 0:
                hi = bi / a if hi is None else min(hi, bi / a)
            else:
                lo = bi / a if lo is None else max(lo, bi / a)
        if lo is None or hi is None:
            raise ValueError("unbounded interval")
        return max(Fraction(0), hi - lo)
    tot = Fraction(0)
    for i, (ai, bi) in enumerate(zip(A, b, strict=True)):
        if bi == 0:
            continue                       # the term vanishes
        j = next(k for k in range(d) if ai[k] != 0)
        # x_j = (bi - sum_{k != j} ai_k x_k) / ai_j ; substitute into every other row
        A2, b2 = [], []
        for t, (at, bt) in enumerate(zip(A, b, strict=True)):
            if t == i:
                continue
            f = at[j] / ai[j]
            A2.append([at[k] - f * ai[k] for k in range(d) if k != j])
            b2.append(bt - f * bi)
        sub = volume_lasserre(A2, b2) if A2 else Fraction(0)
        tot += bi / abs(ai[j]) * sub
    return tot / d


def main():
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from dt13_n4_exact import region
    from exact_poly import volume
    F = Fraction
    bad = 0
    checked = 0
    for r in (F(1), F(3), F(4844, 1000), F(12651, 1000)):
        for L in (F(199, 200), F(1, 2), F(999999, 1000000)):
            for k in range(5):
                A, b = region(k, r, L)
                v1, v2 = volume(A, b), volume_lasserre(A, b)
                checked += 1
                if v1 != v2:
                    bad += 1
                    print("MISMATCH", r, L, k, v1, v2)
    print(f"exact_poly vs Lasserre: {checked} regions, mismatches {bad}")
    # negative falsifiers: a perturbed input must change the value (the checker is not a constant)
    A, b = region(2, F(3), F(199, 200))
    base = volume_lasserre(A, b)
    A2, b2 = region(2, F(3) + F(1, 10 ** 6), F(199, 200))
    A3, b3 = region(2, F(3), F(199, 200) - F(1, 10 ** 6))
    print("falsifier r+1e-6 changes value:", volume_lasserre(A2, b2) != base,
          "| falsifier L-1e-6 changes value:", volume_lasserre(A3, b3) != base)
    print("checker PASS" if bad == 0 else "checker FAIL")


if __name__ == "__main__":
    main()
