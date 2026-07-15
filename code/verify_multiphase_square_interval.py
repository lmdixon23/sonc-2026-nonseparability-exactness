#!/usr/bin/env python3
"""High-precision and interval checks for a three-negative square phase example."""
import mpmath as mp
mp.mp.dps=80

# Positive square family 1+e^x+e^y+kappa e^(x+y).
kappa=mp.mpf('1000000')
t=mp.mpf(2)/5
# Negative exponents and positive phase lambda0=(1/3,1/3,1/3).
B=[(mp.mpf(1)/5,mp.mpf(1)/5),(mp.mpf(4)/5,mp.mpf(1)/5),(mp.mpf(1)/5,mp.mpf(4)/5)]
lam=[mp.mpf(1)/3]*3
v=[mp.mpf(0),mp.mpf(5)/3,-mp.mpf(5)/3]
D=[sum(B[j][i]*v[j] for j in range(3)) for i in range(2)]
assert abs(D[0]-1)<mp.mpf('1e-70') and abs(D[1]+1)<mp.mpf('1e-70')

# Solve the square mean-fiber odds ratio for q=p11.
A=kappa-1
Bb=-(kappa*(2*t)+1-2*t)
C=kappa*t*t
disc=Bb*Bb-4*A*C
roots=[(-Bb+mp.sqrt(disc))/(2*A),(-Bb-mp.sqrt(disc))/(2*A)]
q=[r for r in roots if max(0,2*t-1)<r<t][0]
Sigma=mp.matrix([[t*(1-t),q-t*t],[q-t*t,t*(1-t)]])
Sinv=Sigma**-1
Dvec=mp.matrix(D)
metric=(Dvec.T*Sinv*Dvec)[0]
entropy=sum(v[j]*v[j]/lam[j] for j in range(3))
hess_dir=entropy-metric
assert hess_dir<mp.mpf('-1000')

# Legendre conjugate for the weighted square family.
def q_from_mean(m1,m2):
    A=kappa-1
    Bq=-(kappa*(m1+m2)+1-m1-m2)
    Cq=kappa*m1*m2
    di=Bq*Bq-4*A*Cq
    rs=[(-Bq+mp.sqrt(di))/(2*A),(-Bq-mp.sqrt(di))/(2*A)]
    lo=max(mp.mpf(0),m1+m2-1); hi=min(m1,m2)
    return [r for r in rs if lo<r<hi][0]

def phi_star(m1,m2):
    qq=q_from_mean(m1,m2)
    p=[1-m1-m2+qq,m1-qq,m2-qq,qq]
    return sum(x*mp.log(x) for x in p)-qq*mp.log(kappa)

def G(L):
    m1=sum(L[j]*B[j][0] for j in range(3)); m2=sum(L[j]*B[j][1] for j in range(3))
    return sum(x*mp.log(x) for x in L)-phi_star(m1,m2)

eps=mp.mpf('0.002')
Lp=[lam[j]+eps*v[j] for j in range(3)]
Lm=[lam[j]-eps*v[j] for j in range(3)]
gap=G(lam)-(G(Lp)+G(Lm))/2
assert min(Lp+Lm)>0
assert gap>mp.mpf('1e-6')

# Independent interval enclosure of the directional Hessian using mpmath.iv.
iv=mp.iv
ki=iv.mpf([str(kappa),str(kappa)])
ti=iv.mpf([str(t),str(t)])
Ai=ki-1
Bi=-(ki*(2*ti)+1-2*ti)
Ci=ki*ti*ti
disci=Bi*Bi-4*Ai*Ci
# identify the smaller feasible root, which approaches t from below
qi=(-Bi-iv.sqrt(disci))/(2*Ai)
a=ti*(1-ti); off=qi-ti*ti
# D=(1,-1) is an eigenvector with eigenvalue a-off.
metric_i=2/(a-off)
entropy_i=iv.mpf([str(entropy),str(entropy)])
hess_i=entropy_i-metric_i
assert float(hess_i.b)<0

print('VERDICT: VERIFIED')
print('q =',mp.nstr(q,30))
print('inverse_metric =',mp.nstr(metric,30))
print('entropy_curvature =',mp.nstr(entropy,30))
print('directional_hessian =',mp.nstr(hess_dir,30))
print('certified_chord_gap =',mp.nstr(gap,30))
print('interval_hessian_upper =',hess_i.b)
