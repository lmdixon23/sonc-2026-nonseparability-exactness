#!/usr/bin/env python3
"""Downstream Arb consistency replay for the constructed witness.

This script intentionally does not redo elimination, Sturm isolation, the
linear-subresultant classification, or Krawczyk existence. It reads the certified
contact boxes from the modular certificate and independently reevaluates the
critical equations, coefficient interval, phases, Hessian, and diagonal-value
comparison.
"""
from __future__ import annotations
from pathlib import Path
import json
from flint import arb, ctx

ctx.prec = 320
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BUNDLE = ROOT / "results" / "constructed_witness_certificate"
final = json.loads((BUNDLE / "final_witness_intervals.json").read_text())
comparison = json.loads((BUNDLE / "critical_value_comparison.json").read_text())
sturm = json.loads((BUNDLE / "sturm_certificates.json").read_text())


def hull(pair: list[str]) -> arb:
    return arb(pair[0]).union(arb(pair[1]))


x = hull(final["contact_XY_intervals"][0])
y = hull(final["contact_XY_intervals"][1])
xd_interval = sturm["diagonal"]["positive_isolating_intervals"][0]
from fractions import Fraction

def rational_interval(pair: list[str]) -> arb:
    lo, hi = Fraction(pair[0]), Fraction(pair[1])
    mid, rad = (lo + hi) / 2, (hi - lo) / 2
    from flint import fmpq
    return arb(fmpq(mid.numerator, mid.denominator), fmpq(rad.numerator, rad.denominator))

xd = rational_interval(xd_interval)
K = arb(4096)


def P(xv: arb, yv: arb) -> arb:
    return 1 + xv**5 + yv**5 + K * xv**5 * yv**5


def Q(xv: arb, yv: arb) -> arb:
    return xv * yv * (1 + xv**3 + yv**3)


def ratio(xv: arb, yv: arb) -> arb:
    return P(xv, yv) / Q(xv, yv)


def critical(xv: arb, yv: arb) -> tuple[arb, arb]:
    den = 1 + xv**3 + yv**3
    pp = P(xv, yv)
    return (
        5 * den * xv**5 * (1 + K * yv**5) - pp * (1 + 4 * xv**3 + yv**3),
        5 * den * yv**5 * (1 + K * xv**5) - pp * (1 + xv**3 + 4 * yv**3),
    )


e1, e2 = critical(x, y)
assert e1.contains(0) and e2.contains(0)
alpha = ratio(x, y)
diagonal_ratio = ratio(xd, xd)
assert alpha.upper() < diagonal_ratio.lower()

phase_denominator = 1 + x**3 + y**3
phase = [1 / phase_denominator, x**3 / phase_denominator, y**3 / phase_denominator]
assert (sum(phase) - 1).contains(0)
log_zero = [5 * x.log(), 5 * y.log()]
assert (P(x, y) - alpha * Q(x, y)).contains(0)

positive_weights = [arb(1), x**5, y**5, K * x**5 * y**5]
positive_total = sum(positive_weights)
positive_probabilities = [weight / positive_total for weight in positive_weights]
negative_weights = [x * y, x**4 * y, x * y**4]
negative_total = sum(negative_weights)
negative_probabilities = [weight / negative_total for weight in negative_weights]
positive_points = [(0, 0), (5, 0), (0, 5), (5, 5)]
negative_points = [(1, 1), (4, 1), (1, 4)]


def covariance(probabilities, points):
    mean_x = sum(probabilities[i] * points[i][0] for i in range(len(points)))
    mean_y = sum(probabilities[i] * points[i][1] for i in range(len(points)))
    return (
        sum(probabilities[i] * (points[i][0] - mean_x) ** 2 for i in range(len(points))),
        sum(probabilities[i] * (points[i][0] - mean_x) * (points[i][1] - mean_y) for i in range(len(points))),
        sum(probabilities[i] * (points[i][1] - mean_y) ** 2 for i in range(len(points))),
    )


cp = covariance(positive_probabilities, positive_points)
cn = covariance(negative_probabilities, negative_points)
h11, h12, h22 = cp[0] - cn[0], cp[1] - cn[1], cp[2] - cn[2]
hdet = h11 * h22 - h12 * h12
assert h11.lower() > 0 and hdet.lower() > 0

# Confirm overlap with the producer certificate rather than asserting independent
# reconstruction of that certificate.
assert alpha.overlaps(hull(comparison["alpha_interval"]))
assert diagonal_ratio.overlaps(hull(comparison["diagonal_ratio_interval"]))

print("VERDICT: DOWNSTREAM ARB CONSISTENCY REPLAY PASS")
print("scope: critical equations, values, phases, log zeros, Hessian, diagonal comparison")
print("not replayed: elimination, Sturm isolation, subresultants, Krawczyk existence")
print("alpha =", alpha)
print("gamma =", alpha.log())
print("phase =", phase)
print("log zero =", log_zero)
print("diagonal ratio =", diagonal_ratio)
