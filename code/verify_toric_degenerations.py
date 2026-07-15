#!/usr/bin/env python3
"""Numerical corroboration of covariance collapse in four representative geometries.

This is non-load-bearing evidence. The theorem is analytic.
"""

import json
from pathlib import Path
import numpy as np
from scipy.optimize import root, minimize

OUT=Path(__file__).resolve().parents[1]/'results'/'toric_degenerations.json'

def softmax(v):
    m=np.max(v); e=np.exp(v-m); return e/e.sum()

def data(A,h,tau,z):
    logits=A@z-tau*h
    p=softmax(logits)
    mean=p@A
    C=A-mean
    S=(C.T*p)@C
    phi=np.max(logits)+np.log(np.exp(logits-np.max(logits)).sum())
    return mean,S,p,phi

def inverse_mean(A,h,tau,m):
    def fun(z): return data(A,h,tau,z)[0]-m
    sol=root(fun,np.zeros(A.shape[1]),tol=1e-11)
    if not sol.success or np.linalg.norm(fun(sol.x))>1e-8:
        def obj(z):
            md=data(A,h,tau,z)
            return md[3]-m@z
        sol2=minimize(obj,np.zeros(A.shape[1]),method='BFGS',tol=1e-11)
        z=sol2.x
    else: z=sol.x
    return z,data(A,h,tau,z)

def rprime(A,h,tau,b,c,lam):
    D=b-c; m=c+lam*D
    z,(mm,S,p,phi)=inverse_mean(A,h,tau,m)
    return 1/(lam*(1-lam))-D@np.linalg.inv(S)@D, np.linalg.eigvalsh(S).tolist()

examples=[
    dict(name='square_boundary',
         A=[[0,0],[1,0],[0,1],[1,1]],h=[0,0,0,1],tau=8,
         b=[.25,.25],c=[.8,.4],lam=2/7),
    dict(name='convex_quadrilateral',
         A=[[0,0],[2,0],[3,1],[0,2]],
         h=[-2,3,-2,1],tau=1.5,
         b=[1.5,.4],c=[.7,1.1],lam=26/29),
    dict(name='triangle_plus_interior',
         A=[[0,0],[3,0],[0,3],[.8,.9]],h=[0,0,0,-1],tau=4,
         b=[1.0,.3],c=[.25,1.0],lam=0.46558704453441296),
    dict(name='three_dimensional_circuit',
         A=[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[.2,.2,.2]],
         h=[0,0,0,0,-1],tau=4,
         b=[.3,.15,.05],c=[.3,.05,.15],lam=.5),
]
records=[]
for e in examples:
    A=np.array(e['A'],float); h=np.array(e['h'],float)
    b=np.array(e['b'],float); c=np.array(e['c'],float)
    rp,eigs=rprime(A,h,e['tau'],b,c,e['lam'])
    assert rp < 0, (e['name'],rp)
    records.append(dict(name=e['name'],R_prime=float(rp),covariance_eigenvalues=eigs))

report={'status':'non-load-bearing numerical corroboration','records':records}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(report,indent=2)+'\n')
print('VERDICT: VERIFIED')
for r in records: print(r['name'],r['R_prime'])
