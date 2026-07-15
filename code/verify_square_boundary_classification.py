#!/usr/bin/env python3
"""Exact rational checks for square common-refinement boundary cases."""

from fractions import Fraction as F


def L1(p): return p[0]-p[1]
def L2(p): return p[0]+p[1]-1

def shared_cell(b,c):
    return L1(b)*L1(c) >= 0 and L2(b)*L2(c) >= 0

cases = [
    ((F(1,4),F(1,4)), (F(4,5),F(2,5)), False), # main diagonal point, anti separated
    ((F(1,4),F(1,4)), (F(2,5),F(4,5)), False),
    ((F(1,2),F(1,2)), (F(4,5),F(1,5)), True),  # center belongs to all cells
    ((F(1,4),F(3,4)), (F(1,5),F(4,5)), True),  # same wall/cell closure
    ((F(1,5),F(1,10)), (F(4,5),F(9,10)), False),
    ((F(1,5),F(1,10)), (F(2,5),F(1,5)), True),
]
for b,c,expected in cases:
    assert shared_cell(b,c) is expected, (b,c)

print('VERDICT: VERIFIED')
print('boundary_cases =', len(cases))
