#!/usr/bin/env python3
"""Validate the modular constructed-witness certificate, schema, and hashes."""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import gzip
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BUNDLE = ROOT / "results" / "constructed_witness_certificate"
MANIFEST = BUNDLE / "certificate_manifest.json"
manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
assert manifest["schema"] == "sonc-constructed-witness-certificate-v1"
assert manifest["generator_command"] == "python code/reconstruct_certificate_from_exact_data.py"
required = {
    "original_critical_equations.json",
    "resultant_factorization.json",
    "sturm_certificates.json",
    "sturm_sequences_full.json.gz",
    "subresultant_linear_relation.json",
    "krawczyk_certificate.json",
    "critical_value_comparison.json",
    "final_witness_intervals.json",
    "certificate_claim_crosswalk.json",
    "README.md",
}
assert required <= set(manifest["artifact_hashes"])
for name, expected in manifest["artifact_hashes"].items():
    path = BUNDLE / name
    assert path.is_file(), name
    assert sha256(path.read_bytes()).hexdigest() == expected, name

required_sources = {
    "constructed_witness_exact_core.py",
    "reconstruct_certificate_from_exact_data.py",
    "replay_downstream_interval_arithmetic.py",
    "verify_certificate_schema.py",
    "verify_certificate_portability.py",
    "generate_computational_appendix.py",
}
assert required_sources <= set(manifest["source_hashes"])
for name, expected in manifest["source_hashes"].items():
    path = HERE / name
    assert path.is_file(), name
    assert sha256(path.read_bytes()).hexdigest() == expected, name

original = json.loads((BUNDLE / "original_critical_equations.json").read_text())
assert "E1(X,Y)+E2(X,Y)=S(X+Y,XY)" in original["identities"]
assert "E1(X,Y)-E2(X,Y)=(X-Y)J(X+Y,XY)" in original["identities"]

resultant = json.loads((BUNDLE / "resultant_factorization.json").read_text())
factors = {entry["name"]: entry for entry in resultant["factorization"]}
assert resultant["resultant_degree_in_s"] == 53
assert factors["s_plus_1"]["exponent"] == 6
assert factors["R12"]["exponent"] == 1 and len(factors["R12"]["coefficients_ascending"]) == 13
assert factors["R35"]["exponent"] == 1 and len(factors["R35"]["coefficients_ascending"]) == 36

sturm = json.loads((BUNDLE / "sturm_certificates.json").read_text())
assert sturm["R12"]["positive_root_count"] == 0
assert sturm["R35"]["positive_root_count"] == 3
assert sturm["diagonal"]["positive_root_count"] == 1
for key in ("diagonal", "R12", "R35"):
    assert sturm[key]["sturm_sequence_sha256"]
    for item in sturm[key]["interval_variation_certificates"]:
        assert item["root_count"] == 1
with gzip.open(BUNDLE / "sturm_sequences_full.json.gz", "rt", encoding="utf-8") as stream:
    full_sturm = json.load(stream)
for key in ("diagonal", "R12", "R35"):
    assert len(full_sturm[key]) == sturm[key]["sturm_sequence_length"]

sub = json.loads((BUNDLE / "subresultant_linear_relation.json").read_text())
assert sub["linear_subresultant_s_plus_1_exponent"] == 2
assert len(sub["A_coefficients_ascending"]) == sub["A_degree"] + 1
assert len(sub["B_coefficients_ascending"]) == sub["B_degree"] + 1
assert sum(int(item["physical"]) for item in sub["candidate_classification"]) == 1

kraw = json.loads((BUNDLE / "krawczyk_certificate.json").read_text())
for key in (
    "input_box", "midpoint", "f_mid", "jacobian_mid", "preconditioner",
    "jacobian_box", "krawczyk_image", "strict_interior_margins",
):
    assert key in kraw
assert "strictly inside" in kraw["verdict"]

comparison = json.loads((BUNDLE / "critical_value_comparison.json").read_text())
assert comparison["strict_comparison"] == "alpha.upper < diagonal_ratio.lower"
assert set(comparison["contact_hessian_log_coordinates"]) == {"h11", "h12", "h22", "determinant"}

final = json.loads((BUNDLE / "final_witness_intervals.json").read_text())
assert final["exact_zero_count"] == 2
assert len(final["negative_support"]) == 3
assert len(final["global_nonnegativity_basis"]) >= 7

print("VERDICT: CERTIFICATE SCHEMA, SOURCE HASHES, AND ARTIFACT HASHES VERIFIED")
print("artifacts =", len(manifest["artifact_hashes"]))
print("sources =", len(manifest["source_hashes"]))
