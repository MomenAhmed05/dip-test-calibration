"""AUTO-GENERATED from exact_poly.py (Step 15 / DT14): identical algorithm with python-flint fmpq rationals in place of
fractions.Fraction (2.5x faster; exact). Regenerate: sed "s/from fractions import Fraction/from flint import fmpq as Fraction/" exact_poly.py"""
"""Exact volume of a bounded polytope {x in R^d : A x <= b} in rational arithmetic (Step 14 / DT13, Part X).

Method: exact vertex enumeration (every d-subset of constraints, solved with Fractions, kept if feasible),
then a recursive fan triangulation over the face lattice given by tight-vertex sets:
  vol(P) = sum over facets F not containing v0 of vol(conv(v0, F)), with each facet triangulated the same way.
Simplex volume = |det(v_1 - v_0, ..., v_d - v_0)| / d!.
All arithmetic is exact (fractions.Fraction); no floating point enters a reported value.
Self-test: python exact_poly.py  (unit cube, standard simplex, order simplex, a random polytope vs MC).
"""
from flint import fmpq as Fraction
from itertools import combinations
from math import factorial


def _solve(M, y):
    """Exact Gaussian elimination; returns the solution or None if singular."""
    d = len(M)
    A = [list(M[i]) + [y[i]] for i in range(d)]
    for c in range(d):
        p = next((i for i in range(c, d) if A[i][c] != 0), None)
        if p is None:
            return None
        A[c], A[p] = A[p], A[c]
        piv = A[c][c]
        A[c] = [v / piv for v in A[c]]
        for i in range(d):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[c], strict=True)]
    return [A[i][d] for i in range(d)]


def _det(M):
    d = len(M)
    A = [list(r) for r in M]
    det = Fraction(1)
    for c in range(d):
        p = next((i for i in range(c, d) if A[i][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            det = -det
        det *= A[c][c]
        for i in range(c + 1, d):
            if A[i][c] != 0:
                f = A[i][c] / A[c][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[c], strict=True)]
    return det


def _rank(vs):
    """Affine rank (dimension of the affine hull) of a list of points."""
    if len(vs) <= 1:
        return 0
    rows = [[a - b for a, b in zip(v, vs[0], strict=True)] for v in vs[1:]]
    r, m = 0, [list(x) for x in rows]
    ncol = len(vs[0])
    for c in range(ncol):
        p = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if p is None:
            continue
        m[r], m[p] = m[p], m[r]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c] / m[r][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r], strict=True)]
        r += 1
    return r


def vertices(A, b):
    """Exact vertices of {A x <= b} (A, b lists of Fractions). Returns (verts, tight) with tight[i] = set of
    constraint indices tight at vertex i."""
    m, d = len(A), len(A[0])
    found = {}
    for S in combinations(range(m), d):
        x = _solve([A[i] for i in S], [b[i] for i in S])
        if x is None:
            continue
        if all(sum(a * v for a, v in zip(A[i], x, strict=True)) <= b[i] for i in range(m)):
            found[tuple(x)] = None
    verts = list(found)
    tight = [frozenset(i for i in range(m) if sum(a * v for a, v in zip(A[i], x, strict=True)) == b[i])
             for x in verts]
    return verts, tight


def _faces_simplices(vidx, verts, tight, dim):
    """Triangulate the face with vertex index set vidx (affine dimension dim) into simplices (lists of
    vertex indices of length dim + 1), by a fan from its first vertex."""
    vidx = sorted(vidx)
    if dim == 0:
        return [[vidx[0]]]
    v0 = vidx[0]
    # facets of this face: maximal tight-subsets of affine dimension dim - 1
    cons = set().union(*(tight[i] for i in vidx))
    facets = set()
    for c in cons:
        sub = frozenset(i for i in vidx if c in tight[i])
        if len(sub) >= dim and sub != frozenset(vidx) and _rank([verts[i] for i in sub]) == dim - 1:
            facets.add(sub)
    out = []
    for F in facets:
        if v0 in F:
            continue
        for s in _faces_simplices(F, verts, tight, dim - 1):
            out.append([v0] + s)
    return out


def volume(A, b):
    """Exact d-dimensional volume of the bounded polytope {A x <= b}; 0 if lower-dimensional or empty."""
    A = [[Fraction(v) for v in row] for row in A]
    b = [Fraction(v) for v in b]
    d = len(A[0])
    verts, tight = vertices(A, b)
    if len(verts) <= d or _rank(verts) < d:
        return Fraction(0)
    tot = Fraction(0)
    for s in _faces_simplices(range(len(verts)), verts, tight, d):
        p0 = verts[s[0]]
        M = [[a - c for a, c in zip(verts[j], p0, strict=True)] for j in s[1:]]
        tot += abs(_det(M))
    return tot / factorial(d)


def _selftest():
    import random
    F = Fraction
    # unit cube in 3D
    A = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1]]
    b = [1, 1, 1, 0, 0, 0]
    assert volume(A, b) == 1, volume(A, b)
    # standard simplex in 4D: x_i >= 0, sum <= 1 -> 1/24
    A = [[-1 if j == i else 0 for j in range(4)] for i in range(4)] + [[1, 1, 1, 1]]
    assert volume(A, [0, 0, 0, 0, 1]) == F(1, 24)
    # order simplex 0 <= x1 <= ... <= x4 <= 1 -> 1/24
    A = [[-1, 0, 0, 0], [1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1], [0, 0, 0, 1]]
    assert volume(A, [0, 0, 0, 0, 1]) == F(1, 24)
    # n = 4 uniform tail: 24 vol{order, s2 >= r s1, s2 >= r s3} = 2/((r+1)(r+2))
    for r in (F(1), F(3, 2), F(3), F(7, 2)):
        # s2 >= r s1:  -(x3 - x2) + r (x2 - x1) <= 0  ->  -r x1 + (r + 1) x2 - x3 <= 0
        # s2 >= r s3:  -(x3 - x2) + r (x4 - x3) <= 0  ->  x2 - (1 + r) x3 + r x4 <= 0
        Ar = A + [[-r, r + 1, -1, 0], [0, 1, -(1 + r), r]]
        v = 24 * volume(Ar, [0, 0, 0, 0, 1, 0, 0])
        assert v == 2 / ((r + 1) * (r + 2)), (r, v)
    # random 3D polytope vs Monte Carlo (sanity of the triangulation on irregular shapes)
    rng = random.Random(7)
    A2 = [[F(v) for v in row] for row in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1])]
    b2 = [F(v) for v in (1, 1, 1, 0, 0, 0)]
    for _ in range(4):
        A2.append([F(rng.randint(-5, 5)) for _ in range(3)])
        b2.append(F(rng.randint(1, 6)))
    v = volume(A2, b2)
    hits = 0
    N = 200000
    for _ in range(N):
        x = [rng.random() for _ in range(3)]
        hits += all(sum(float(a) * xi for a, xi in zip(row, x, strict=True)) <= float(bb) for row, bb in zip(A2, b2, strict=True))
    print(f"random polytope: exact {float(v):.5f} vs MC {hits / N:.5f} (se {((hits / N) * (1 - hits / N) / N) ** .5:.5f})")
    print("exact_poly self-test PASS")


if __name__ == "__main__":
    _selftest()
