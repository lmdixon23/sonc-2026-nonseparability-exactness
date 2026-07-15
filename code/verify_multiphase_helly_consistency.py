#!/usr/bin/env python3
"""Exact consistency checks for the symmetric no-pair family."""
import sympy as sp
n,k=sp.symbols('n k', positive=True)
coefQ=sp.Rational(1,4)
coefP=(n+1)/(n+1+k)
expr=sp.factor(coefQ-coefP)
expected=sp.factor((k-3*(n+1))/(4*(n+1+k)))
assert sp.simplify(expr-expected)==0
# Positive denominator means the sign is exactly the sign of k-3(n+1).
for nv in range(1,8):
    assert sp.simplify(expr.subs({n:nv,k:3*(nv+1)+1}))>0
    assert sp.simplify(expr.subs({n:nv,k:3*(nv+1)-1}))<0
print('VERDICT: VERIFIED')
print('hessian_factor =',expr)
print('threshold = kappa > 3(n+1)')
