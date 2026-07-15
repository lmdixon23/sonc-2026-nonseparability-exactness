#!/usr/bin/env python3
"""Exact verification of the canonical planar minimal obstruction."""

from fractions import Fraction
import sympy as sp

x, y = sp.symbols("x y", positive=True)

# Canonical parameters a=1, b=2.
f_square = (x - 1)**2 * (x - 2)**2
f_transverse = 2 * (y + 1/y - 2)
f = sp.expand(f_square + f_transverse)

expected = x**4 + 13*x**2 + 2*y + 2/y - 12*x - 6*x**3
assert sp.simplify(f - expected) == 0

# Exact coefficient dictionary in exponent coordinates.
coeffs = {
    (4, 0): Fraction(1),
    (2, 0): Fraction(13),
    (0, 1): Fraction(2),
    (0, -1): Fraction(2),
    (1, 0): Fraction(-12),
    (3, 0): Fraction(-6),
}
assert all(c > 0 for a, c in coeffs.items() if a in {(4,0),(2,0),(0,1),(0,-1)})
assert all(c < 0 for a, c in coeffs.items() if a in {(1,0),(3,0)})

# Two exact zeros and singularity.
for point in [(1, 1), (2, 1)]:
    subs = {x: point[0], y: point[1]}
    assert sp.simplify(f.subs(subs)) == 0
    assert sp.simplify(sp.diff(f, x).subs(subs)) == 0
    assert sp.simplify(sp.diff(f, y).subs(subs)) == 0

# Interior barycentric coordinates in the outer triangle
# vertices W+ = (0,1), W- = (0,-1), U4 = (4,0).
def barycentric_axis(t):
    # (t,0) = l1*(0,1) + l2*(0,-1) + l3*(4,0)
    return (Fraction(4-t, 8), Fraction(4-t, 8), Fraction(t, 4))

for t in (1, 3):
    lambdas = barycentric_axis(t)
    assert sum(lambdas) == 1
    assert all(v > 0 for v in lambdas)

# Lower-facet check for the regular lifting:
# h(2,0)=-1 and h=0 at the three outer vertices.
points = [(0,1), (0,-1), (4,0), (2,0)]
heights = {(0,1):0, (0,-1):0, (4,0):0, (2,0):-1}
cells = [
    [(0,1), (0,-1), (2,0)],
    [(0,1), (2,0), (4,0)],
    [(0,-1), (2,0), (4,0)],
]

X, Y = sp.symbols("X Y")
for cell in cells:
    aa, bb, cc = sp.symbols("aa bb cc")
    equations = [
        sp.Eq(aa*p[0] + bb*p[1] + cc, heights[p])
        for p in cell
    ]
    plane = sp.solve(equations, (aa, bb, cc), dict=True)[0]
    for p in points:
        predicted = plane[aa]*p[0] + plane[bb]*p[1] + plane[cc]
        # Lower hull: every lifted point is on or above the supporting plane.
        assert sp.Rational(heights[p]) >= sp.simplify(predicted)

# Membership separation:
# (1,0) lies in the left/base triangle; (3,0) lies on the right shared edge.
def point_in_triangle(pt, tri):
    M = sp.Matrix([
        [tri[0][0], tri[1][0], tri[2][0]],
        [tri[0][1], tri[1][1], tri[2][1]],
        [1, 1, 1],
    ])
    rhs = sp.Matrix([pt[0], pt[1], 1])
    lam = M.LUsolve(rhs)
    return all(v >= 0 for v in lam), tuple(sp.simplify(v) for v in lam)

membership = {}
for pt in [(1,0), (3,0)]:
    membership[str(pt)] = []
    for idx, cell in enumerate(cells):
        inside, lam = point_in_triangle(pt, cell)
        if inside:
            membership[str(pt)].append({"cell": idx, "barycentric": [str(v) for v in lam]})

assert {entry["cell"] for entry in membership[str((1,0))]} == {0}
assert {entry["cell"] for entry in membership[str((3,0))]} == {1, 2}

print("VERDICT: VERIFIED")
print("expanded_signomial =", expected)
print("zeros = [(1, 1), (2, 1)]")
print("negative_interior_barycentrics =", {
    "u": [str(v) for v in barycentric_axis(1)],
    "3u": [str(v) for v in barycentric_axis(3)],
})
print("regular_subdivision_membership =", membership)
