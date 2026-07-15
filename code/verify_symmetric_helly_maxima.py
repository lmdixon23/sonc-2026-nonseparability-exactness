#!/usr/bin/env python3
"""Orbit-polished numerical corroboration for the symmetric Helly family."""

from itertools import permutations
import mpmath as mp

mp.mp.dps = 45


def logratio(zfree, n, kappa):
    t = [mp.mpf(1)] + [mp.e**z for z in zfree]
    g = mp.e**(sum(mp.log(x) for x in t)/(n+1))
    P = sum(t) + kappa*g
    Q = sum(mp.sqrt(x*g) for x in t)
    return mp.log(Q)-mp.log(P)


def polish(seed, n, kappa):
    z = [mp.mpf(x) for x in seed]
    for _ in range(100):
        h = mp.mpf('1e-7')
        f0 = logratio(z,n,kappa)
        grad=[]; H=mp.matrix(n,n)
        for i in range(n):
            zp=z.copy(); zm=z.copy(); zp[i]+=h; zm[i]-=h
            grad.append((logratio(zp,n,kappa)-logratio(zm,n,kappa))/(2*h))
        for i in range(n):
            for j in range(i,n):
                zpp=z.copy(); zpm=z.copy(); zmp=z.copy(); zmm=z.copy()
                zpp[i]+=h; zpp[j]+=h
                zpm[i]+=h; zpm[j]-=h
                zmp[i]-=h; zmp[j]+=h
                zmm[i]-=h; zmm[j]-=h
                H[i,j]=H[j,i]=(logratio(zpp,n,kappa)-logratio(zpm,n,kappa)-logratio(zmp,n,kappa)+logratio(zmm,n,kappa))/(4*h*h)
        gvec=mp.matrix(grad)
        if mp.norm(gvec)<mp.mpf('1e-20'): break
        try: step=mp.lu_solve(H,-gvec)
        except Exception: step=-gvec
        a=mp.mpf(1)
        while a>mp.mpf('1e-12'):
            trial=[z[k]+a*step[k] for k in range(n)]
            if logratio(trial,n,kappa)>=f0-mp.mpf('1e-28'):
                z=trial; break
            a/=2
        else: break
    return z, logratio(z,n,kappa)


def distance(a,b):
    return mp.sqrt(sum((a[i]-b[i])**2 for i in range(len(a))))


for n in (2,3):
    kappa=mp.mpf(4*(n+1))
    seeds=[]
    for i in range(n):
        v=[mp.mpf('-1.5')]*n; v[i]=mp.mpf(3); seeds.append(v)
    seeds += [[mp.mpf(2)]*n, [mp.mpf(-2)]*n]
    candidates=[polish(s,n,kappa) for s in seeds]
    best=max(candidates,key=lambda x:x[1])
    tbest=[mp.mpf(1)]+[mp.e**z for z in best[0]]
    maxima=[]
    for perm in permutations(range(n+1)):
        tp=[tbest[perm[i]] for i in range(n+1)]
        seed=[mp.log(tp[i+1]/tp[0]) for i in range(n)]
        z,val=polish(seed,n,kappa)
        if all(distance(z,w)>mp.mpf('1e-8') for w,_ in maxima):
            maxima.append((z,val))
    vals=[v for _,v in maxima]
    spread=max(vals)-min(vals)
    dmin=min(distance(maxima[i][0],maxima[j][0]) for i in range(len(maxima)) for j in range(i+1,len(maxima)))
    assert len(maxima)>=2
    assert spread<mp.mpf('1e-14')
    assert dmin>mp.mpf('0.1')
    print('dimension',n,'orbit_maxima',len(maxima),'value_spread',mp.nstr(spread,6),'min_distance',mp.nstr(dmin,8))
print('VERDICT: VERIFIED')
