#!/usr/bin/env python3
"""Exact reconstruction core for the certified three-negative SONC witness.

This module starts only from the exact polynomials P and Q. It derives the
logarithmic critical equations, the diagonal equation, the original symmetric
off-diagonal system, its full resultant factorization, Sturm root counts, the
linear subresultant, the physical-root classification, and an Arb Krawczyk
certificate. It writes a modular machine-readable certificate bundle.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from typing import Any
import json
import platform
import sys
import gzip

import sympy as sp
from flint import arb, fmpq, ctx
import flint

ctx.prec = 400
try:
    sys.set_int_max_str_digits(0)
except AttributeError:
    pass

X, Y, s, p = sp.symbols("X Y s p")
K = sp.Integer(4096)
P = 1 + X**5 + Y**5 + K * X**5 * Y**5
Q = X * Y * (1 + X**3 + Y**3)
DEN = 1 + X**3 + Y**3

# Cleared logarithmic critical equations X*P_X*Q - X*Q_X*P = 0,
# Y*P_Y*Q - Y*Q_Y*P = 0, divided by the common positive factor XY.
E1 = sp.expand(5 * DEN * X**5 * (1 + K * Y**5) - P * (1 + 4 * X**3 + Y**3))
E2 = sp.expand(5 * DEN * Y**5 * (1 + K * X**5) - P * (1 + X**3 + 4 * Y**3))


def _to_sp_symmetric(expr: sp.Expr) -> sp.Expr:
    """Express a symmetric polynomial in s=X+Y and p=XY."""
    sym, rem, mapping = sp.symmetrize(sp.expand(expr), [X, Y], formal=True)
    if rem != 0:
        raise AssertionError(f"expression is not symmetric; remainder={rem}")
    replacements = {mapping[0][0]: s, mapping[1][0]: p}
    return sp.expand(sym.subs(replacements))


# Original off-diagonal system. For X != Y,
# E1=E2=0 iff S(s,p)=J(s,p)=0 with s=X+Y, p=XY.
DIFF_QUOTIENT = sp.cancel((E1 - E2) / (X - Y))
assert sp.expand((X - Y) * DIFF_QUOTIENT - (E1 - E2)) == 0
S_EXPR = _to_sp_symmetric(E1 + E2)
J_EXPR = _to_sp_symmetric(DIFF_QUOTIENT)
assert sp.expand(S_EXPR.subs({s: X + Y, p: X * Y}) - (E1 + E2)) == 0
assert sp.expand((X - Y) * J_EXPR.subs({s: X + Y, p: X * Y}) - (E1 - E2)) == 0
S_POLY = sp.Poly(S_EXPR, p, s, domain=sp.ZZ)
J_POLY = sp.Poly(J_EXPR, p, s, domain=sp.ZZ)

DIAG_POLY = sp.Poly(sp.factor(E1.subs(Y, X)), X, domain=sp.ZZ)
assert DIAG_POLY.degree() == 13

RESULTANT_EXPR = sp.resultant(S_EXPR, J_EXPR, p)
RESULTANT_FACTOR = sp.factor_list(RESULTANT_EXPR)
RESULTANT_CONST = sp.Integer(RESULTANT_FACTOR[0])
RESULTANT_FACTORS = [(sp.Poly(f, s, domain=sp.ZZ), int(e)) for f, e in RESULTANT_FACTOR[1]]
R12 = next(f for f, e in RESULTANT_FACTORS if f.degree() == 12)
R35 = next(f for f, e in RESULTANT_FACTORS if f.degree() == 35)
LINEAR_FACTOR = next((f, e) for f, e in RESULTANT_FACTORS if f.degree() == 1)
assert LINEAR_FACTOR[0].as_expr() == s + 1 and LINEAR_FACTOR[1] == 6
assert sp.expand(RESULTANT_CONST * (s + 1) ** 6 * R12.as_expr() * R35.as_expr() - RESULTANT_EXPR) == 0

SUBRESULTANTS = sp.subresultants(S_EXPR, J_EXPR, p)
SUBRESULTANT_DEGREES = [int(sp.degree(q, p)) for q in SUBRESULTANTS]
LINEAR_SUBRESULTANT_EXPR = next(q for q in SUBRESULTANTS if sp.degree(q, p) == 1)
lin_factor = sp.factor_list(LINEAR_SUBRESULTANT_EXPR)
lin_const = sp.Integer(lin_factor[0])
lin_splus_exp = 0
lin_core = sp.Integer(1)
for factor, exponent in lin_factor[1]:
    factor_poly = sp.Poly(factor, p, s)
    if factor_poly.degree(p) == 0 and sp.expand(factor - (s + 1)) == 0:
        lin_splus_exp += int(exponent)
    else:
        lin_core *= factor ** exponent
lin_core = sp.expand(lin_core)
assert sp.degree(lin_core, p) == 1
A_POLY = sp.Poly(sp.Poly(lin_core, p).coeff_monomial(p), s, domain=sp.ZZ)
B_POLY = sp.Poly(sp.Poly(lin_core, p).coeff_monomial(1), s, domain=sp.ZZ)
assert sp.expand(lin_const * (s + 1) ** lin_splus_exp * (A_POLY.as_expr() * p + B_POLY.as_expr()) - LINEAR_SUBRESULTANT_EXPR) == 0


def fraction_arb(value: Fraction | int | str) -> arb:
    if isinstance(value, str):
        value = Fraction(value)
    elif isinstance(value, int):
        value = Fraction(value, 1)
    return arb(fmpq(value.numerator, value.denominator))


def interval_arb(lo: Fraction | str, hi: Fraction | str) -> arb:
    lo_f, hi_f = Fraction(lo), Fraction(hi)
    mid, rad = (lo_f + hi_f) / 2, (hi_f - lo_f) / 2
    return arb(fmpq(mid.numerator, mid.denominator), fmpq(rad.numerator, rad.denominator))


def arb_bounds(value: arb, digits: int = 90) -> list[str]:
    return [value.lower().str(digits), value.upper().str(digits)]


def arb_mid(value: arb, digits: int = 90) -> str:
    return value.mid().str(digits)


def poly_eval(poly: sp.Poly, values: dict[sp.Symbol, arb]) -> arb:
    total = arb(0)
    for monom, coeff in poly.terms():
        coeff_q = sp.Rational(coeff)
        term = fraction_arb(Fraction(int(coeff_q.p), int(coeff_q.q)))
        for generator, exponent in zip(poly.gens, monom):
            term *= values[generator] ** exponent
        total += term
    return total


def polynomial_coefficients(poly: sp.Poly, variable: sp.Symbol) -> list[str]:
    p1 = sp.Poly(poly.as_expr(), variable, domain=sp.ZZ)
    return [str(int(p1.nth(i))) for i in range(p1.degree() + 1)]


def bivariate_terms(poly: sp.Poly) -> list[dict[str, Any]]:
    return [
        {"monomial": list(map(int, monomial)), "coefficient": str(coeff)}
        for monomial, coeff in poly.terms()
    ]


def rational_interval(interval: tuple[sp.Rational, sp.Rational]) -> list[str]:
    return [str(interval[0]), str(interval[1])]


def _sign_of_rational(value: sp.Rational) -> int:
    return 1 if value > 0 else (-1 if value < 0 else 0)


def _variations(signs: list[int]) -> int:
    nz = [sign for sign in signs if sign]
    return sum(int(a != b) for a, b in zip(nz, nz[1:]))


def _sturm_signs_at(seq: list[sp.Poly], variable: sp.Symbol, point: sp.Rational | str) -> list[int]:
    if point == "+infinity":
        return [_sign_of_rational(sp.Rational(q.LC())) for q in seq]
    return [_sign_of_rational(sp.Rational(q.eval(point))) for q in seq]


def sturm_record(poly: sp.Poly, variable: sp.Symbol, eps: sp.Rational) -> tuple[dict[str, Any], list[list[str]]]:
    seq = [sp.Poly(q, variable, domain=sp.QQ) for q in sp.sturm(poly.as_expr(), variable)]
    full_coefficients = [[str(c) for c in reversed(q.all_coeffs())] for q in seq]
    canonical = json.dumps(full_coefficients, separators=(",", ":"), ensure_ascii=True)
    positive_intervals = [I for I, multiplicity in poly.intervals(eps=eps) if I[1] > 0]
    interval_records = []
    for left, right in positive_intervals:
        left_signs = _sturm_signs_at(seq, variable, left)
        right_signs = _sturm_signs_at(seq, variable, right)
        interval_records.append({
            "interval": rational_interval((left, right)),
            "left_signs": left_signs,
            "right_signs": right_signs,
            "left_variations": _variations(left_signs),
            "right_variations": _variations(right_signs),
            "root_count": _variations(left_signs) - _variations(right_signs),
        })
    zero_signs = _sturm_signs_at(seq, variable, sp.Rational(0))
    infinity_signs = _sturm_signs_at(seq, variable, "+infinity")
    record = {
        "degree": poly.degree(),
        "coefficients_ascending": polynomial_coefficients(poly, variable),
        "positive_root_count": int(poly.count_roots(0, sp.oo)),
        "positive_isolating_intervals": [rational_interval(I) for I in positive_intervals],
        "sturm_sequence_length": len(seq),
        "sturm_sequence_sha256": sha256(canonical.encode()).hexdigest(),
        "signs_at_zero": zero_signs,
        "signs_at_positive_infinity": infinity_signs,
        "variations_at_zero": _variations(zero_signs),
        "variations_at_positive_infinity": _variations(infinity_signs),
        "interval_variation_certificates": interval_records,
        "full_sequence_artifact": "sturm_sequences_full.json.gz",
    }
    return record, full_coefficients


@dataclass
class KrawczykData:
    input_box: list[arb]
    midpoint: list[arb]
    f_mid: list[arb]
    jacobian_mid: list[list[arb]]
    preconditioner: list[list[arb]]
    jacobian_box: list[list[arb]]
    image: list[arb]
    margins: list[dict[str, arb]]
    determinant_mid: arb


def krawczyk2(expr1: sp.Expr, expr2: sp.Expr, box: list[arb]) -> KrawczykData:
    polys = [sp.Poly(expr1, s, p), sp.Poly(expr2, s, p)]
    derivatives = [[sp.Poly(sp.diff(expr, var), s, p) for var in (s, p)] for expr in (expr1, expr2)]
    midpoint = [box[0].mid(), box[1].mid()]
    f_mid = [poly_eval(poly, {s: midpoint[0], p: midpoint[1]}) for poly in polys]
    jacobian_box = [[poly_eval(derivatives[i][j], {s: box[0], p: box[1]}) for j in range(2)] for i in range(2)]
    jacobian_mid = [[poly_eval(derivatives[i][j], {s: midpoint[0], p: midpoint[1]}) for j in range(2)] for i in range(2)]
    determinant_mid = jacobian_mid[0][0] * jacobian_mid[1][1] - jacobian_mid[0][1] * jacobian_mid[1][0]
    if determinant_mid.contains(0):
        raise AssertionError("midpoint Jacobian determinant contains zero")
    preconditioner = [
        [jacobian_mid[1][1] / determinant_mid, -jacobian_mid[0][1] / determinant_mid],
        [-jacobian_mid[1][0] / determinant_mid, jacobian_mid[0][0] / determinant_mid],
    ]
    y = [
        midpoint[i] - sum(preconditioner[i][j] * f_mid[j] for j in range(2))
        for i in range(2)
    ]
    matrix = [
        [
            arb(int(i == j)) - sum(preconditioner[i][k] * jacobian_box[k][j] for k in range(2))
            for j in range(2)
        ]
        for i in range(2)
    ]
    displacement = [box[i] - midpoint[i] for i in range(2)]
    image = [y[i] + sum(matrix[i][j] * displacement[j] for j in range(2)) for i in range(2)]
    if not (box[0].contains_interior(image[0]) and box[1].contains_interior(image[1])):
        raise AssertionError("Krawczyk image is not strictly inside input box")
    margins = [
        {
            "lower": image[i].lower() - box[i].lower(),
            "upper": box[i].upper() - image[i].upper(),
        }
        for i in range(2)
    ]
    if not all(m["lower"].lower() > 0 and m["upper"].lower() > 0 for m in margins):
        raise AssertionError("Krawczyk interior margins are not strictly positive")
    return KrawczykData(
        input_box=box,
        midpoint=midpoint,
        f_mid=f_mid,
        jacobian_mid=jacobian_mid,
        preconditioner=preconditioner,
        jacobian_box=jacobian_box,
        image=image,
        margins=margins,
        determinant_mid=determinant_mid,
    )


def _serialize_matrix(matrix: list[list[arb]]) -> list[list[list[str]]]:
    return [[arb_bounds(entry) for entry in row] for row in matrix]


def _write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_certificate_bundle(root: Path) -> dict[str, Any]:
    root.mkdir(parents=True, exist_ok=True)
    eps = sp.Rational(1, 10) ** 80

    diagonal_sturm, diagonal_sturm_full = sturm_record(DIAG_POLY, X, eps)
    r12_sturm, r12_sturm_full = sturm_record(R12, s, eps)
    r35_sturm, r35_sturm_full = sturm_record(R35, s, eps)
    assert diagonal_sturm["positive_root_count"] == 1
    assert r12_sturm["positive_root_count"] == 0
    assert r35_sturm["positive_root_count"] == 3

    r35_intervals = [(sp.Rational(a), sp.Rational(b)) for a, b in r35_sturm["positive_isolating_intervals"]]
    candidates: list[dict[str, Any]] = []
    physical_box: list[arb] | None = None
    for interval in r35_intervals:
        s_box = interval_arb(interval[0], interval[1])
        a_value = poly_eval(A_POLY, {s: s_box})
        b_value = poly_eval(B_POLY, {s: s_box})
        if a_value.contains(0):
            raise AssertionError("linear subresultant leading coefficient contains zero")
        p_box = -b_value / a_value
        discriminant = s_box * s_box - 4 * p_box
        physical = p_box.lower() > 0 and discriminant.lower() > 0
        candidate = {
            "s_interval": rational_interval(interval),
            "A_interval": arb_bounds(a_value),
            "B_interval": arb_bounds(b_value),
            "p_interval": arb_bounds(p_box),
            "xy_discriminant_interval": arb_bounds(discriminant),
            "physical": bool(physical),
        }
        candidates.append(candidate)
        if physical:
            physical_box = [s_box, p_box]
    assert sum(int(c["physical"]) for c in candidates) == 1
    assert physical_box is not None

    # The original-system Krawczyk certificate begins from the exact Sturm and
    # linear-subresultant boxes. No root midpoint or contact box is hard-coded.
    krawczyk = krawczyk2(S_EXPR, J_EXPR, physical_box)
    s_box, p_box = krawczyk.image
    xy_discriminant = s_box * s_box - 4 * p_box
    if xy_discriminant.lower() <= 0:
        raise AssertionError("physical Krawczyk image does not have positive discriminant")
    sqrt_discriminant = xy_discriminant.sqrt()
    x_small = (s_box - sqrt_discriminant) / 2
    y_large = (s_box + sqrt_discriminant) / 2
    if not (x_small.lower() > 0 and y_large.lower() > 0 and x_small.upper() < y_large.lower()):
        raise AssertionError("contact boxes are not positive and disjoint")

    diagonal_interval = [(sp.Rational(a), sp.Rational(b)) for a, b in diagonal_sturm["positive_isolating_intervals"]][0]
    x_diagonal = interval_arb(diagonal_interval[0], diagonal_interval[1])

    def ratio(xv: arb, yv: arb) -> arb:
        return (1 + xv**5 + yv**5 + arb(4096) * xv**5 * yv**5) / (xv * yv * (1 + xv**3 + yv**3))

    alpha = ratio(x_small, y_large)
    diagonal_ratio = ratio(x_diagonal, x_diagonal)
    if not alpha.upper() < diagonal_ratio.lower():
        raise AssertionError("off-diagonal critical value is not strictly below diagonal critical value")
    gamma = alpha.log()

    phase_denominator = 1 + x_small**3 + y_large**3
    phase = [1 / phase_denominator, x_small**3 / phase_denominator, y_large**3 / phase_denominator]
    log_zero = [5 * x_small.log(), 5 * y_large.log()]

    # Direct original critical-equation evaluation.
    e1_value = poly_eval(sp.Poly(E1, X, Y), {X: x_small, Y: y_large})
    e2_value = poly_eval(sp.Poly(E2, X, Y), {X: x_small, Y: y_large})
    if not (e1_value.contains(0) and e2_value.contains(0)):
        raise AssertionError("contact boxes fail the original critical equations")

    # Hessian of log(P/Q) in log X, log Y coordinates: covariance difference.
    positive_weights = [arb(1), x_small**5, y_large**5, arb(4096) * x_small**5 * y_large**5]
    positive_total = sum(positive_weights)
    positive_probabilities = [weight / positive_total for weight in positive_weights]
    negative_weights = [x_small * y_large, x_small**4 * y_large, x_small * y_large**4]
    negative_total = sum(negative_weights)
    negative_probabilities = [weight / negative_total for weight in negative_weights]
    positive_points = [(0, 0), (5, 0), (0, 5), (5, 5)]
    negative_points = [(1, 1), (4, 1), (1, 4)]

    def covariance(probabilities: list[arb], points: list[tuple[int, int]]) -> tuple[arb, arb, arb]:
        mean_x = sum(probabilities[i] * points[i][0] for i in range(len(points)))
        mean_y = sum(probabilities[i] * points[i][1] for i in range(len(points)))
        return (
            sum(probabilities[i] * (points[i][0] - mean_x) ** 2 for i in range(len(points))),
            sum(probabilities[i] * (points[i][0] - mean_x) * (points[i][1] - mean_y) for i in range(len(points))),
            sum(probabilities[i] * (points[i][1] - mean_y) ** 2 for i in range(len(points))),
        )

    cov_positive = covariance(positive_probabilities, positive_points)
    cov_negative = covariance(negative_probabilities, negative_points)
    h11 = cov_positive[0] - cov_negative[0]
    h12 = cov_positive[1] - cov_negative[1]
    h22 = cov_positive[2] - cov_negative[2]
    hdet = h11 * h22 - h12 * h12
    if not (h11.lower() > 0 and hdet.lower() > 0):
        raise AssertionError("contact Hessian is not certified positive definite")

    original_equations = {
        "variables": ["X", "Y", "s", "p"],
        "P": str(P),
        "Q": str(Q),
        "E1": str(E1),
        "E2": str(E2),
        "diagonal_polynomial": str(DIAG_POLY.as_expr()),
        "S_original_symmetric_sum": str(S_EXPR),
        "J_original_symmetric_difference_quotient": str(J_EXPR),
        "identities": [
            "E1(X,Y)+E2(X,Y)=S(X+Y,XY)",
            "E1(X,Y)-E2(X,Y)=(X-Y)J(X+Y,XY)",
            "For X != Y, E1=E2=0 iff S=J=0.",
        ],
        "E1_terms": bivariate_terms(sp.Poly(E1, X, Y)),
        "E2_terms": bivariate_terms(sp.Poly(E2, X, Y)),
        "S_terms": bivariate_terms(S_POLY),
        "J_terms": bivariate_terms(J_POLY),
    }
    _write_json(root / "original_critical_equations.json", original_equations)

    resultant_data = {
        "elimination_variable": "p",
        "resultant_degree_in_s": int(sp.degree(RESULTANT_EXPR, s)),
        "factorization_constant": str(RESULTANT_CONST),
        "factorization": [
            {
                "name": "s_plus_1",
                "exponent": 6,
                "coefficients_ascending": ["1", "1"],
            },
            {
                "name": "R12",
                "exponent": 1,
                "coefficients_ascending": polynomial_coefficients(R12, s),
            },
            {
                "name": "R35",
                "exponent": 1,
                "coefficients_ascending": polynomial_coefficients(R35, s),
            },
        ],
        "identity": "Res_p(S,J)=constant*(s+1)^6*R12(s)*R35(s)",
        "sha256_of_expanded_resultant_decimal": sha256(str(sp.expand(RESULTANT_EXPR)).encode()).hexdigest(),
    }
    _write_json(root / "resultant_factorization.json", resultant_data)

    sturm_data = {
        "epsilon": str(eps),
        "diagonal": diagonal_sturm,
        "R12": r12_sturm,
        "R35": r35_sturm,
    }
    _write_json(root / "sturm_certificates.json", sturm_data)
    with gzip.open(root / "sturm_sequences_full.json.gz", "wt", encoding="utf-8", compresslevel=9) as stream:
        json.dump({
            "diagonal": diagonal_sturm_full,
            "R12": r12_sturm_full,
            "R35": r35_sturm_full,
        }, stream, separators=(",", ":"), sort_keys=True)
        stream.write("\n")

    subresultant_data = {
        "subresultant_degrees_in_p": SUBRESULTANT_DEGREES,
        "linear_subresultant_factorization_constant": str(lin_const),
        "linear_subresultant_s_plus_1_exponent": lin_splus_exp,
        "A_degree": A_POLY.degree(),
        "B_degree": B_POLY.degree(),
        "A_coefficients_ascending": polynomial_coefficients(A_POLY, s),
        "B_coefficients_ascending": polynomial_coefficients(B_POLY, s),
        "identity": "L(s,p)=constant*(s+1)^2*(A(s)*p+B(s))",
        "candidate_classification": candidates,
    }
    _write_json(root / "subresultant_linear_relation.json", subresultant_data)

    krawczyk_data = {
        "system": ["S(s,p)", "J(s,p)"],
        "input_box": [arb_bounds(entry) for entry in krawczyk.input_box],
        "midpoint": [arb_mid(entry) for entry in krawczyk.midpoint],
        "f_mid": [arb_bounds(entry) for entry in krawczyk.f_mid],
        "jacobian_mid": _serialize_matrix(krawczyk.jacobian_mid),
        "midpoint_jacobian_determinant": arb_bounds(krawczyk.determinant_mid),
        "preconditioner": _serialize_matrix(krawczyk.preconditioner),
        "jacobian_box": _serialize_matrix(krawczyk.jacobian_box),
        "krawczyk_image": [arb_bounds(entry) for entry in krawczyk.image],
        "strict_interior_margins": [
            {"lower": arb_bounds(item["lower"]), "upper": arb_bounds(item["upper"])}
            for item in krawczyk.margins
        ],
        "verdict": "Krawczyk image lies strictly inside the input box; one and only one S=J root lies in the box.",
    }
    _write_json(root / "krawczyk_certificate.json", krawczyk_data)

    comparison_data = {
        "alpha_interval": arb_bounds(alpha),
        "diagonal_ratio_interval": arb_bounds(diagonal_ratio),
        "strict_comparison": "alpha.upper < diagonal_ratio.lower",
        "supporting_plane_height_gamma_interval": arb_bounds(gamma),
        "contact_hessian_log_coordinates": {
            "h11": arb_bounds(h11),
            "h12": arb_bounds(h12),
            "h22": arb_bounds(h22),
            "determinant": arb_bounds(hdet),
        },
        "original_critical_equation_residuals": {
            "E1": arb_bounds(e1_value),
            "E2": arb_bounds(e2_value),
        },
    }
    _write_json(root / "critical_value_comparison.json", comparison_data)

    final_data = {
        "positive_coefficients": [1, 1, 1, 4096],
        "positive_support": [[0, 0], [1, 0], [0, 1], [1, 1]],
        "negative_support": [["1/5", "1/5"], ["4/5", "1/5"], ["1/5", "4/5"]],
        "negative_coefficient_alpha_interval": arb_bounds(alpha),
        "contact_XY_intervals": [arb_bounds(x_small), arb_bounds(y_large)],
        "swapped_contact_XY_intervals": [arb_bounds(y_large), arb_bounds(x_small)],
        "log_zero_intervals": [arb_bounds(log_zero[0]), arb_bounds(log_zero[1])],
        "swapped_log_zero_intervals": [arb_bounds(log_zero[1]), arb_bounds(log_zero[0])],
        "phase_intervals": [arb_bounds(entry) for entry in phase],
        "swapped_phase_intervals": [arb_bounds(phase[0]), arb_bounds(phase[2]), arb_bounds(phase[1])],
        "exact_zero_count": 2,
        "global_nonnegativity_basis": [
            "R=P/Q is coercive in logarithmic coordinates.",
            "The diagonal and original off-diagonal critical systems are exhaustively enumerated.",
            "R12 has zero positive roots; R35 has exactly three positive roots.",
            "Only one R35 candidate has p>0 and s^2-4p>0.",
            "The original-system Krawczyk operator certifies the unique physical off-diagonal root.",
            "The off-diagonal critical value is strictly below the only diagonal critical value.",
            "Symmetry gives exactly two global minimizers.",
        ],
    }
    _write_json(root / "final_witness_intervals.json", final_data)

    crosswalk = {
        "manuscript_claims": [
            {"claim": "original critical equations", "artifact": "original_critical_equations.json"},
            {"claim": "complete original-system resultant factorization", "artifact": "resultant_factorization.json"},
            {"claim": "R12 has zero and R35 has three positive roots", "artifact": "sturm_certificates.json"},
            {"claim": "unique p recovery and physical candidate classification", "artifact": "subresultant_linear_relation.json"},
            {"claim": "existence and uniqueness of the physical root", "artifact": "krawczyk_certificate.json"},
            {"claim": "strict off-diagonal versus diagonal value comparison", "artifact": "critical_value_comparison.json"},
            {"claim": "coefficient, phases, and two zero boxes", "artifact": "final_witness_intervals.json"},
        ]
    }
    _write_json(root / "certificate_claim_crosswalk.json", crosswalk)

    readme = """# Constructed-witness certificate bundle

This directory is generated from the exact benchmark polynomials

P = 1 + X^5 + Y^5 + 4096 X^5 Y^5,
Q = X Y (1 + X^3 + Y^3).

Regenerate it from the repository root with:

    python code/reconstruct_certificate_from_exact_data.py

The bundle contains the original critical equations, complete original-system
resultant factorization, exact Sturm certificates and compressed full sequences,
the linear subresultant and physical-root classification, an original-system Arb
Krawczyk inclusion, the strict critical-value comparison, final witness intervals,
and a manuscript claim crosswalk. certificate_manifest.json hashes every bundle
artifact and every load-bearing reconstruction source file.

replay_downstream_interval_arithmetic.py is intentionally narrower: it checks the
downstream interval arithmetic from certified boxes and does not redo elimination,
Sturm isolation, subresultants, or Krawczyk existence.
"""
    (root / "README.md").write_text(readme, encoding="utf-8")

    file_hashes = {}
    for path in sorted(root.iterdir()):
        if not path.is_file() or path.name == "certificate_manifest.json":
            continue
        file_hashes[path.name] = sha256(path.read_bytes()).hexdigest()
    source_paths = [
        Path(__file__).resolve(),
        Path(__file__).with_name("reconstruct_certificate_from_exact_data.py"),
        Path(__file__).with_name("replay_downstream_interval_arithmetic.py"),
        Path(__file__).with_name("verify_certificate_schema.py"),
        Path(__file__).with_name("generate_computational_appendix.py"),
    ]
    source_hashes = {
        path.name: sha256(path.read_bytes()).hexdigest()
        for path in source_paths if path.exists()
    }
    manifest = {
        "schema": "sonc-constructed-witness-certificate-v1",
        "generator_command": "python code/reconstruct_certificate_from_exact_data.py",
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "sympy": sp.__version__,
            "python_flint": getattr(flint, "__version__", "unknown"),
            "arb_precision_bits": ctx.prec,
        },
        "source_hashes": source_hashes,
        "artifact_hashes": file_hashes,
        "verdict": "RIGOROUS COMPLETE CONSTRUCTED-WITNESS CERTIFICATE PASS",
    }
    _write_json(root / "certificate_manifest.json", manifest)

    summary = {
        "status": manifest["verdict"],
        "bundle_directory": root.name,
        "bundle_manifest_sha256": sha256((root / "certificate_manifest.json").read_bytes()).hexdigest(),
        "alpha_interval": arb_bounds(alpha),
        "supporting_plane_height_gamma_interval": arb_bounds(gamma),
        "contact_XY_intervals": [arb_bounds(x_small), arb_bounds(y_large)],
        "log_zero_intervals": [arb_bounds(log_zero[0]), arb_bounds(log_zero[1])],
        "phase_intervals": [arb_bounds(entry) for entry in phase],
        "diagonal_ratio_interval": arb_bounds(diagonal_ratio),
        "R12_positive_root_count": r12_sturm["positive_root_count"],
        "R35_positive_root_count": r35_sturm["positive_root_count"],
        "offdiagonal_candidates": candidates,
        "contact_hessian_log_coordinates": comparison_data["contact_hessian_log_coordinates"],
        "proof_summary": final_data["global_nonnegativity_basis"],
    }
    return summary
