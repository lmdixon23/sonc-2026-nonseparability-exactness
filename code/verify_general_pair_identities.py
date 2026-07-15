#!/usr/bin/env python3
"""Symbolic identities used by the general two-negative theorem."""

import sympy as sp

lam = sp.symbols('lam', positive=True)
D1, D2 = sp.symbols('D1 D2', real=True)
s11, s12, s22 = sp.symbols('s11 s12 s22', real=True)
D = sp.Matrix([D1, D2])
S = sp.Matrix([[s11, s12], [s12, s22]])
a = lam * (1-lam)

lhs = sp.factor((S - a*D*D.T).det())
rhs = sp.factor(S.det() * (1 - a*(D.T*S.inv()*D)[0]))
assert sp.simplify(lhs-rhs) == 0

# Cauchy-Schwarz inverse covariance inequality in a diagonal symbolic model.
x1,x2,n1,n2=sp.symbols('x1 x2 n1 n2', positive=True)
Sd=sp.diag(x1,x2)
Dv=sp.Matrix([D1,D2])
nv=sp.Matrix([n1,n2])
expr=(Dv.T*Sd.inv()*Dv)[0]*(nv.T*Sd*nv)[0]-(Dv.dot(nv))**2
assert sp.simplify(expr-(D1*n2*x2-D2*n1*x1)**2/(x1*x2)) == 0

print('VERDICT: VERIFIED')
print('rank_one_hessian_identity = VERIFIED')
print('inverse_covariance_cauchy_schwarz = VERIFIED')
