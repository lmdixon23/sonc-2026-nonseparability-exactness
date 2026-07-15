#!/usr/bin/env python3
"""Symbolic checks for the planar square converse formulas."""

import sympy as sp

m1, m2, q, k = sp.symbols("m1 m2 q k", positive=True)
p00 = 1 - m1 - m2 + q
p10 = m1 - q
p01 = m2 - q
p11 = q

odds = sp.simplify(p00 * p11 / (p10 * p01))
assert sp.simplify(odds - q * (1 - m1 - m2 + q) / ((m1 - q) * (m2 - q))) == 0

Sigma = sp.Matrix([
    [m1 * (1 - m1), q - m1 * m2],
    [q - m1 * m2, m2 * (1 - m2)],
])

t = sp.symbols("t", positive=True)

main_limit = sp.simplify(Sigma.subs({m1: t, m2: t, q: t}))
expected_main = t * (1 - t) * sp.Matrix([[1, 1], [1, 1]])
assert sp.simplify(main_limit - expected_main) == sp.zeros(2)

anti_limit = sp.simplify(Sigma.subs({m1: t, m2: 1 - t, q: 0}))
expected_anti = t * (1 - t) * sp.Matrix([[1, -1], [-1, 1]])
assert sp.simplify(anti_limit - expected_anti) == sp.zeros(2)

assert expected_main * sp.Matrix([1, -1]) == sp.zeros(2, 1)
assert expected_anti * sp.Matrix([1, 1]) == sp.zeros(2, 1)

lam = sp.symbols("lam", positive=True)
D1, D2 = sp.symbols("D1 D2", real=True)
D = sp.Matrix([D1, D2])
Sinv = sp.MatrixSymbol("Sinv", 2, 2)

# Rank-one determinant identity:
# det(Sigma-a DD^T)=det(Sigma)*(1-a D^T Sigma^{-1}D).
a = sp.symbols("a", positive=True)
s11, s12, s22 = sp.symbols("s11 s12 s22", real=True)
S = sp.Matrix([[s11, s12], [s12, s22]])
lhs = sp.factor((S - a * D * D.T).det())
rhs = sp.factor(S.det() * (1 - a * (D.T * S.inv() * D)[0]))
assert sp.simplify(lhs - rhs) == 0

print("VERDICT: VERIFIED")
print("main_diagonal_covariance_limit =", expected_main)
print("main_null_direction = (1, -1)")
print("anti_diagonal_covariance_limit =", expected_anti)
print("anti_null_direction = (1, 1)")
print("rank_one_hessian_identity = VERIFIED")
