#!/usr/bin/env python3
"""Separately written optimizer-based replay of the multiphase Hessian sign."""
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

A=np.array([[0.,0.],[1.,0.],[0.,1.],[1.,1.]])
h=np.array([0.,0.,0.,-1.])
tau=8.0
Bneg=np.array([[.2,.2],[.8,.2],[.2,.8]])
lam=np.array([1/3,1/3,1/3])
v=np.array([0.,5/3,-5/3])


def phi(z):
    return logsumexp(A@z-tau*h)


def phistar(m):
    sol=minimize(lambda z: phi(z)-m@z, np.array([-4.,-4.]), method='BFGS', options={'gtol':1e-12,'maxiter':2000})
    if not sol.success and np.linalg.norm(sol.jac)>1e-7:
        raise RuntimeError(sol.message)
    return m@sol.x-phi(sol.x), sol.x


def covariance(z):
    logits=A@z-tau*h
    p=np.exp(logits-logsumexp(logits))
    mean=p@A
    X=A-mean
    return (X.T*p)@X


def G(x):
    m=x@Bneg
    ps,_=phistar(m)
    return np.sum(x*np.log(x))-ps

m=lam@Bneg
ps,z=phistar(m)
S=covariance(z)
D=Bneg.T@v
analytic=np.sum(v*v/lam)-D@np.linalg.inv(S)@D
eps=2e-4
fd=(G(lam+eps*v)-2*G(lam)+G(lam-eps*v))/eps**2
assert analytic<0 and fd<0, (analytic,fd,z,S)
assert abs((analytic-fd)/analytic)<5e-3, (analytic,fd)
print('VERDICT: VERIFIED')
print('analytic_directional_hessian',analytic)
print('finite_difference_directional_hessian',fd)
print('inverse_mean_z',z.tolist())
