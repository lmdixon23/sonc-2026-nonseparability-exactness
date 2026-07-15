#!/usr/bin/env python3
"""Exact verification of the convex-position diagonal obstruction."""

import sympy as sp

X, Y = sp.symbols("X Y", positive=True)
U = sp.sqrt(X / Y)
V = sp.sqrt(X * Y)

H = U - 6 * sp.sqrt(U) + 13 - 12 / sp.sqrt(U) + 4 / U
factor_H = sp.simplify((sp.sqrt(U) - 1)**2 * (sp.sqrt(U) - 2)**2 / U)
assert sp.simplify(H - factor_H) == 0

F_factor = sp.Rational(13, 2) * (V - 1)**2 + V * H
F_expanded = (
    sp.Rational(13, 2) * (1 + X * Y)
    + X + 4 * Y
    - 12 * X**sp.Rational(1, 4) * Y**sp.Rational(3, 4)
    - 6 * X**sp.Rational(3, 4) * Y**sp.Rational(1, 4)
)
assert sp.simplify(F_factor - F_expanded) == 0

zeros = [(sp.Integer(1), sp.Integer(1)), (sp.Integer(4), sp.Rational(1, 4))]
for x0, y0 in zeros:
    assert sp.simplify(F_expanded.subs({X: x0, Y: y0})) == 0

# Hessian in logarithmic coordinates s=log V, t=log U.
s, t = sp.symbols("s t", real=True)
Ht = sp.exp(-t) * (sp.exp(t / 2) - 1)**2 * (sp.exp(t / 2) - 2)**2
G = sp.Rational(13, 2) * (sp.exp(s) - 1)**2 + sp.exp(s) * Ht
for t0 in [sp.Integer(0), sp.log(4)]:
    Hess = sp.hessian(G, (s, t)).subs({s: 0, t: t0})
    Hess = sp.simplify(Hess)
    assert Hess == sp.diag(13, sp.Rational(1, 2))

# Strict interior barycentric coordinates for b,c in the square are immediate.
b = sp.Matrix([sp.Rational(1,4), sp.Rational(3,4)])
c = sp.Matrix([sp.Rational(3,4), sp.Rational(1,4)])
y1 = sp.Matrix([0, 0])
y2 = sp.Matrix([sp.log(4), -sp.log(4)])
transversality = sp.simplify((b-c).dot(y1-y2))
assert transversality == sp.log(4)

print("VERDICT: VERIFIED")
print("signomial =", F_expanded)
print("factorization =", F_factor)
print("zeros =", zeros)
print("log_hessians = diag(13, 1/2)")
print("bitangency_transversality =", transversality)
