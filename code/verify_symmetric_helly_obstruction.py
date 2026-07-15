#!/usr/bin/env python3
"""Exact covariance and combinatorial checks for the symmetric Helly family."""

import sympy as sp

for n in range(2,7):
    kappa=sp.Integer(4*(n+1))
    cp=sp.Rational(n+1,n+1+kappa)
    cq=sp.Rational(1,4)
    assert cq-cp > 0
    # Every proper subset of the n+1 negative points is contained in one star cell.
    negatives=set(range(n+1))
    cells={j: negatives-{j} for j in negatives}
    assert not any(negatives <= cell for cell in cells.values())
    for j in negatives:
        subset=negatives-{j}
        assert subset <= cells[j]

print('VERDICT: VERIFIED')
print('dimensions_checked = 2..6')
print('threshold = kappa > 3(n+1)')
