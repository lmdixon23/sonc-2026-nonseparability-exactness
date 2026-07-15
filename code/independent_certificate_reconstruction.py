#!/usr/bin/env python3
"""Independent implementation of the exact witness reconstruction.

The script begins from the exact benchmark polynomials P and Q, does not import
the primary certificate module, and checks the committed certificate as a
regression test. It is not an external review.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
from hashlib import sha256
import json
import gzip
import sys

import sympy as sp
from flint import arb, fmpq, ctx

ctx.prec = 400
try:
    sys.set_int_max_str_digits(0)
except AttributeError:
    pass

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BUNDLE = ROOT / "results" / "constructed_witness_certificate"
OUT = ROOT / "results" / "independent_certificate_reconstruction.json"

X, Y, s, p = sp.symbols("X Y s p")
P = 1 + X**5 + Y**5 + 4096 * X**5 * Y**5
Q = X * Y * (1 + X**3 + Y**3)
DEN = 1 + X**3 + Y**3

# Derive the logarithmic critical equations directly.
E1 = sp.expand(5 * DEN * X**5 * (1 + 4096 * Y**5) - P * (1 + 4 * X**3 + Y**3))
E2 = sp.expand(5 * DEN * Y**5 * (1 + 4096 * X**5) - P * (1 + X**3 + 4 * Y**3))
assert sp.expand(X * sp.diff(P, X) * Q - X * sp.diff(Q, X) * P - X * Y * E1) == 0
assert sp.expand(Y * sp.diff(P, Y) * Q - Y * sp.diff(Q, Y) * P - X * Y * E2) == 0


def to_symmetric(expr: sp.Expr) -> sp.Expr:
    sym, rem, mapping = sp.symmetrize(sp.expand(expr), [X, Y], formal=True)
    assert rem == 0
    return sp.expand(sym.subs({mapping[0][0]: s, mapping[1][0]: p}))


D = sp.cancel((E1 - E2) / (X - Y))
assert sp.expand((X - Y) * D - (E1 - E2)) == 0
S = to_symmetric(E1 + E2)
J = to_symmetric(D)
assert sp.expand(S.subs({s: X + Y, p: X * Y}) - (E1 + E2)) == 0
assert sp.expand((X - Y) * J.subs({s: X + Y, p: X * Y}) - (E1 - E2)) == 0

# Diagonal equation and complete original-system elimination.
diagonal = sp.Poly(sp.factor(E1.subs(Y, X)), X, domain=sp.ZZ)
assert diagonal.degree() == 13
resultant_expr = sp.resultant(S, J, p)
constant, raw_factors = sp.factor_list(resultant_expr)
factors = [(sp.Poly(f, s, domain=sp.ZZ), int(e)) for f, e in raw_factors]
R12 = next(poly for poly, exp in factors if poly.degree() == 12)
R35 = next(poly for poly, exp in factors if poly.degree() == 35)
linear = next((poly, exp) for poly, exp in factors if poly.degree() == 1)
assert linear[0].as_expr() == s + 1 and linear[1] == 6
assert sp.expand(sp.Integer(constant) * (s + 1)**6 * R12.as_expr() * R35.as_expr() - resultant_expr) == 0
assert R12.count_roots(0, sp.oo) == 0
assert R35.count_roots(0, sp.oo) == 3
assert diagonal.count_roots(0, sp.oo) == 1

# Independently recover the linear subresultant.
subresultants = sp.subresultants(S, J, p)
linear_sub = next(q for q in subresultants if sp.degree(q, p) == 1)
lin_const, lin_factors = sp.factor_list(linear_sub)
splus_exp = 0
core = sp.Integer(1)
for factor, exponent in lin_factors:
    if sp.degree(factor, p) == 0 and sp.expand(factor - (s + 1)) == 0:
        splus_exp += int(exponent)
    else:
        core *= factor**exponent
core = sp.expand(core)
A = sp.Poly(sp.Poly(core, p).coeff_monomial(p), s, domain=sp.ZZ)
B = sp.Poly(sp.Poly(core, p).coeff_monomial(1), s, domain=sp.ZZ)
assert sp.expand(sp.Integer(lin_const) * (s + 1)**splus_exp * (A.as_expr()*p + B.as_expr()) - linear_sub) == 0


def qarb(value: Fraction | sp.Rational | int) -> arb:
    if isinstance(value, sp.Rational):
        value = Fraction(int(value.p), int(value.q))
    elif isinstance(value, int):
        value = Fraction(value, 1)
    return arb(fmpq(value.numerator, value.denominator))


def interval(lo: sp.Rational, hi: sp.Rational) -> arb:
    mid = (lo + hi) / 2
    rad = (hi - lo) / 2
    return arb(fmpq(int(mid.p), int(mid.q)), fmpq(int(rad.p), int(rad.q)))


def peval(expr: sp.Expr, values: dict[sp.Symbol, arb], gens: tuple[sp.Symbol, ...]) -> arb:
    poly = sp.Poly(expr, *gens)
    total = arb(0)
    for monomial, coeff in poly.terms():
        coeff = sp.Rational(coeff)
        term = qarb(coeff)
        for var, degree in zip(gens, monomial):
            term *= values[var] ** degree
        total += term
    return total


def bounds(x: arb, digits: int = 80) -> list[str]:
    return [x.lower().str(digits), x.upper().str(digits)]


# Derive candidate boxes from exact R35 isolation, not from committed boxes.
eps = sp.Rational(1, 10**80)
r35_intervals = [I for I, mult in R35.intervals(eps=eps) if I[1] > 0]
assert len(r35_intervals) == 3
candidates = []
physical_input = None
for lo, hi in r35_intervals:
    sb = interval(lo, hi)
    av = peval(A.as_expr(), {s: sb}, (s,))
    bv = peval(B.as_expr(), {s: sb}, (s,))
    assert not av.contains(0)
    pb = -bv / av
    disc = sb*sb - 4*pb
    physical = pb.lower() > 0 and disc.lower() > 0
    candidates.append({
        "s_interval": [str(lo), str(hi)],
        "p_interval": bounds(pb),
        "discriminant_interval": bounds(disc),
        "physical": bool(physical),
    })
    if physical:
        physical_input = [sb, pb]
assert sum(int(c["physical"]) for c in candidates) == 1
assert physical_input is not None

# Independently implemented two-dimensional Krawczyk operator for S,J.
funcs = [S, J]
derivs = [[sp.diff(f, v) for v in (s, p)] for f in funcs]
box = physical_input
mid = [box[0].mid(), box[1].mid()]
fmid = [peval(f, {s: mid[0], p: mid[1]}, (s, p)) for f in funcs]
Jmid = [[peval(derivs[i][j], {s: mid[0], p: mid[1]}, (s, p)) for j in range(2)] for i in range(2)]
Jbox = [[peval(derivs[i][j], {s: box[0], p: box[1]}, (s, p)) for j in range(2)] for i in range(2)]
det = Jmid[0][0]*Jmid[1][1] - Jmid[0][1]*Jmid[1][0]
assert not det.contains(0)
C = [[Jmid[1][1]/det, -Jmid[0][1]/det], [-Jmid[1][0]/det, Jmid[0][0]/det]]
y = [mid[i] - sum(C[i][j]*fmid[j] for j in range(2)) for i in range(2)]
M = [[arb(int(i == j)) - sum(C[i][k]*Jbox[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
delta = [box[i] - mid[i] for i in range(2)]
image = [y[i] + sum(M[i][j]*delta[j] for j in range(2)) for i in range(2)]
assert box[0].contains_interior(image[0]) and box[1].contains_interior(image[1])

# Reconstruct X,Y, alpha and the diagonal comparison.
sb, pb = image
disc = sb*sb - 4*pb
assert disc.lower() > 0
root_disc = disc.sqrt()
xb = (sb-root_disc)/2
yb = (sb+root_disc)/2
assert xb.lower() > 0 and xb.upper() < yb.lower()

def ratio(x: arb, y: arb) -> arb:
    return (1+x**5+y**5+arb(4096)*x**5*y**5)/(x*y*(1+x**3+y**3))

alpha = ratio(xb, yb)
diag_interval = [I for I, mult in diagonal.intervals(eps=eps) if I[1] > 0][0]
xd = interval(diag_interval[0], diag_interval[1])
diag_ratio = ratio(xd, xd)
assert alpha.upper() < diag_ratio.lower()

# Compare exact algebraic data and final intervals to delivered artifacts.
resultant_artifact = json.loads((BUNDLE / "resultant_factorization.json").read_text())
rf = {entry["name"]: entry for entry in resultant_artifact["factorization"]}
assert resultant_artifact["factorization_constant"] == str(constant)
assert [str(R12.nth(i)) for i in range(13)] == rf["R12"]["coefficients_ascending"]
assert [str(R35.nth(i)) for i in range(36)] == rf["R35"]["coefficients_ascending"]
assert resultant_artifact["sha256_of_expanded_resultant_decimal"] == sha256(str(sp.expand(resultant_expr)).encode()).hexdigest()

sturm_artifact = json.loads((BUNDLE / "sturm_certificates.json").read_text())
assert sturm_artifact["R12"]["positive_root_count"] == 0
assert sturm_artifact["R35"]["positive_root_count"] == 3
assert sturm_artifact["diagonal"]["positive_root_count"] == 1
with gzip.open(BUNDLE / "sturm_sequences_full.json.gz", "rt", encoding="utf-8") as stream:
    full_sturm = json.load(stream)
for name in ("R12", "R35", "diagonal"):
    assert len(full_sturm[name]) == sturm_artifact[name]["sturm_sequence_length"]

sub_artifact = json.loads((BUNDLE / "subresultant_linear_relation.json").read_text())
assert [str(A.nth(i)) for i in range(A.degree()+1)] == sub_artifact["A_coefficients_ascending"]
assert [str(B.nth(i)) for i in range(B.degree()+1)] == sub_artifact["B_coefficients_ascending"]
assert sum(int(item["physical"]) for item in sub_artifact["candidate_classification"]) == 1

final_artifact = json.loads((BUNDLE / "final_witness_intervals.json").read_text())
def encloses(outer: list[str], inner: arb) -> bool:
    return arb(outer[0]).lower() <= inner.lower() and arb(outer[1]).upper() >= inner.upper()
assert encloses(final_artifact["negative_coefficient_alpha_interval"], alpha)
assert encloses(final_artifact["contact_XY_intervals"][0], xb)
assert encloses(final_artifact["contact_XY_intervals"][1], yb)
assert final_artifact["exact_zero_count"] == 2

report = {
    "status": "INDEPENDENT CERTIFICATE RECONSTRUCTION PASS",
    "independence_limit": "Separate implementation within the repository; not an external review.",
    "input": {"P": str(P), "Q": str(Q)},
    "derived": {
        "diagonal_degree": diagonal.degree(),
        "resultant_degree": int(sp.degree(resultant_expr, s)),
        "resultant_constant": str(constant),
        "R12_positive_roots": int(R12.count_roots(0, sp.oo)),
        "R35_positive_roots": int(R35.count_roots(0, sp.oo)),
        "diagonal_positive_roots": int(diagonal.count_roots(0, sp.oo)),
        "subresultant_degrees": [int(sp.degree(q, p)) for q in subresultants],
        "physical_candidate_count": sum(int(c["physical"]) for c in candidates),
    },
    "krawczyk": {
        "determinant": bounds(det),
        "input_box": [bounds(box[0]), bounds(box[1])],
        "image": [bounds(image[0]), bounds(image[1])],
        "strict_inclusion": True,
    },
    "global_classification": {
        "alpha_interval": bounds(alpha),
        "diagonal_ratio_interval": bounds(diag_ratio),
        "strict_value_comparison": True,
        "exact_global_minimizer_count_by_symmetry_and_complete_critical_enumeration": 2,
    },
    "artifact_comparison": {
        "resultant_coefficients": "MATCH",
        "sturm_counts_and_lengths": "MATCH",
        "linear_subresultant_coefficients": "MATCH",
        "final_interval_enclosures": "MATCH",
    },
}
OUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print("VERDICT: INDEPENDENT CERTIFICATE RECONSTRUCTION PASSED")
print("R12 positive roots =", report["derived"]["R12_positive_roots"])
print("R35 positive roots =", report["derived"]["R35_positive_roots"])
print("physical candidates =", report["derived"]["physical_candidate_count"])
print("output =", OUT.relative_to(ROOT))
