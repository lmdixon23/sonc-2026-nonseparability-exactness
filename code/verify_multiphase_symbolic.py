#!/usr/bin/env python3
"""Symbolic and algebraic checks for the multiphase theorem."""
import sympy as sp

# Generic differential identity for an affine moment map m=B lambda.
l1,l2,l3=sp.symbols('l1 l2 l3', positive=True)
# Work in free simplex coordinates l1,l2, with l3=1-l1-l2.
x,y=sp.symbols('x y', real=True)
l3e=1-l1-l2
E=l1*sp.log(l1)+l2*sp.log(l2)+l3e*sp.log(l3e)
HE=sp.hessian(E,(l1,l2))
expected=sp.Matrix([[1/l1+1/l3e,1/l3e],[1/l3e,1/l2+1/l3e]])
assert sp.simplify(HE-expected)==sp.zeros(2)

# Projective invariance: rho -> rho+c1 changes rho.lambda by c.
c=sp.symbols('c')
r1,r2,r3=sp.symbols('r1 r2 r3')
expr=(r1+c)*l1+(r2+c)*l2+(r3+c)*l3e-(r1*l1+r2*l2+r3*l3e)
assert sp.simplify(expr-c)==0

# k=2 reduction. Tangent vector v=(-1,1) up to scaling.
lam=sp.symbols('lam', positive=True)
Dq=sp.symbols('Dq', positive=True)
scalar=1/lam+1/(1-lam)-Dq
assert sp.simplify(scalar-(1/(lam*(1-lam))-Dq))==0

# Rank-one local maximum condition used by the two-phase specialization.
a,s11,s12,s22,d1,d2=sp.symbols('a s11 s12 s22 d1 d2')
S=sp.Matrix([[s11,s12],[s12,s22]])
D=sp.Matrix([d1,d2])
lhs=sp.factor((S-a*D*D.T).det())
rhs=sp.factor(S.det()*(1-a*(D.T*S.inv()*D)[0]))
assert sp.simplify(lhs-rhs)==0

print('VERDICT: VERIFIED')
print('entropy_hessian =', HE)
print('k2_reduction = VERIFIED')
print('projective_rho_shift = VERIFIED')
