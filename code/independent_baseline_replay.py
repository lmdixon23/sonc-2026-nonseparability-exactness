#!/usr/bin/env python3
"""Independent covariance-collapse replay using direct softmax covariance."""
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


def invert_mean(A,h,tau,m):
    def obj(z): return logsumexp(A@z-tau*h)-m@z
    def jac(z):
        logits=A@z-tau*h
        p=np.exp(logits-logsumexp(logits))
        return p@A-m
    sol=minimize(obj,np.zeros(A.shape[1]),jac=jac,method='BFGS',options={'gtol':1e-11,'maxiter':1000})
    if not sol.success and np.linalg.norm(sol.jac)>1e-7:
        raise RuntimeError(sol.message)
    logits=A@sol.x-tau*h
    p=np.exp(logits-logsumexp(logits))
    mean=p@A
    X=A-mean
    S=(X.T*p)@X
    return sol.x,S

checks=[
    (np.array([[0,0],[1,0],[0,1],[1,1.]],float),np.array([0,0,0,-1.]),np.array([.4,.4]),np.array([1.,-1.])),
    (np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[.2,.2,.2]],float),np.array([0,0,0,0,-1.]),np.array([.3,.1,.1]),np.array([0.,1.,-1.])),
]
for A,h,m,norm in checks:
    values=[]
    for tau in [2.,4.,6.,8.]:
        z,S=invert_mean(A,h,tau,m)
        values.append(float(norm@S@norm))
    assert values[-1] < values[0]*0.08, values
    print('collapse_values',values)
print('VERDICT: VERIFIED')
