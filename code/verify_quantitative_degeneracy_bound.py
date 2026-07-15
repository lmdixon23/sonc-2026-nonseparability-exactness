#!/usr/bin/env python3
"""Directed-interval verification of the explicit degeneration bound.

The analytic proposition is proved in the manuscript. This script certifies the
three-phase square benchmark using Arb rather than binary floating point.
"""
from __future__ import annotations
from pathlib import Path
import json
from flint import arb, fmpq, ctx

ctx.prec = 256
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT_JSON = ROOT / "results" / "quantitative_degeneracy_bound.json"
OUT_TXT = ROOT / "results" / "quantitative_degeneracy_bound.txt"


def q(num: int, den: int = 1) -> arb:
    return arb(fmpq(num, den))


def bounds(x: arb, digits: int = 70) -> list[str]:
    return [x.lower().str(digits), x.upper().str(digits)]

N = arb(4)
Delta = q(1, 2)
Rnu = arb(1)
D_dot_nu = arb(2)
lambda_coordinate = q(1, 3)
v2_sum = q(50, 9)  # 0^2+(5/3)^2+(-5/3)^2
entropy_curvature = v2_sum / lambda_coordinate  # 50/3
Hq = -(q(3, 5) * q(3, 5).log() + q(2, 5) * q(2, 5).log())
log_gap = N.log() - Hq
tau_bound = Rnu**2 * log_gap / (Delta * D_dot_nu**2) * entropy_curvature
assert tau_bound.lower() > arb("5.9440") and tau_bound.upper() < arb("5.9441")

tau = arb(6)
lower_metric = tau * Delta * D_dot_nu**2 / (Rnu**2 * log_gap)
assert lower_metric.lower() > entropy_curvature.upper()

# Exact mean-fiber covariance for P_tau=1+e^z1+e^z2+e^tau e^(z1+z2)
# at m=(2/5,2/5). The feasible q11 solves
# (kappa-1)q^2-(4kappa/5+1/5)q+4kappa/25=0.
kappa = tau.exp()
A = kappa - 1
B = -(q(4, 5) * kappa + q(1, 5))
C = q(4, 25) * kappa
discriminant = B * B - 4 * A * C
root_plus = (-B + discriminant.sqrt()) / (2 * A)
root_minus = (-B - discriminant.sqrt()) / (2 * A)
roots = [root_plus, root_minus]
feasible = [root for root in roots if root.lower() > 0 and root.upper() < q(2, 5).lower()]
assert len(feasible) == 1
q11 = feasible[0]
anti_diagonal_eigenvalue = q(2, 5) - q11
actual_metric = 2 / anti_diagonal_eigenvalue
actual_hessian = entropy_curvature - actual_metric
assert actual_metric.lower() >= lower_metric.lower()
assert actual_hessian.upper() < 0

report = {
    "status": "DIRECTED INTERVAL VERIFICATION PASS",
    "arb_precision_bits": ctx.prec,
    "carrier_entropy_interval": bounds(Hq),
    "entropy_curvature_exact": "50/3",
    "entropy_curvature_interval": bounds(entropy_curvature),
    "explicit_tau_bound_interval": bounds(tau_bound),
    "chosen_tau_exact": "6",
    "certified_metric_lower_bound_interval": bounds(lower_metric),
    "actual_inverse_metric_interval": bounds(actual_metric),
    "actual_directional_hessian_interval": bounds(actual_hessian),
    "mean_fiber_p11_interval": bounds(q11),
    "strict_checks": [
        "tau=6 is strictly above the explicit sufficient bound",
        "the certified metric lower endpoint exceeds the entropy-curvature upper endpoint",
        "the actual directional Hessian upper endpoint is negative",
    ],
}
OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
text = "\n".join([
    "VERDICT: DIRECTED INTERVAL QUANTITATIVE CHECK PASSED",
    f"carrier_entropy = {Hq}",
    f"entropy_curvature = {entropy_curvature}",
    f"explicit_tau_bound = {tau_bound}",
    "chosen_tau = 6",
    f"certified_metric_lower_bound = {lower_metric}",
    f"actual_inverse_metric = {actual_metric}",
    f"actual_directional_hessian = {actual_hessian}",
    f"mean_fiber_p11 = {q11}",
]) + "\n"
OUT_TXT.write_text(text, encoding="utf-8")
print(text, end="")
