#!/usr/bin/env python3
"""Cross-check manuscript certificate claims against delivered artifacts."""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
TEX = (ROOT / "paper" / "main.tex").read_text(encoding="utf-8")
APPENDIX = (ROOT / "paper" / "computational_appendix.tex").read_text(encoding="utf-8")
RUN_ALL = (ROOT / "run_all.sh").read_text(encoding="utf-8")
BUNDLE = ROOT / "results" / "constructed_witness_certificate"
MANIFEST = json.loads((BUNDLE / "certificate_manifest.json").read_text(encoding="utf-8"))

required_manuscript_phrases = [
    "standalone scalar proof",
    "original symmetric equations",
    "R_{12}",
    "downstream Arb replay",
    "it is not presented as an independent elimination",
    "Appendix~\\ref{app:certificate}",
    "membership in this global cone",
    "relative to the affine hull of $\\DeltaK$",
    "Carath\\'eodory's theorem",
    "$Q=\\overline{\\relint Q}$",
]
for phrase in required_manuscript_phrases:
    assert phrase in TEX, phrase

prohibited_overclaims = [
    "the JSON records the complete resultant",
    "independently replayable certificate",
    "independently proves the two-negative theorem",
    "first complete characterization",
    "For a fixed signed support $(A_+,A_-)$, the SONC cone is the conic hull",
]
for phrase in prohibited_overclaims:
    assert phrase.lower() not in TEX.lower(), phrase

required_run_labels = [
    "== Corroborative structural and symbolic checks ==",
    "== Multiphase consistency checks ==",
    "== Exact and validated certificate checks ==",
]
for phrase in required_run_labels:
    assert phrase in RUN_ALL, phrase
assert "== Analytic and structural checks ==" not in RUN_ALL

required_appendix_phrases = [
    "E_1+E_2=S(X+Y,XY)",
    "E_1-E_2=(X-Y)J(X+Y,XY)",
    "R_{12}",
    "N_{(0,\\infty)}(R_{12})=0",
    "Krawczyk image",
    "linear subresultant",
]
for phrase in required_appendix_phrases:
    assert phrase in APPENDIX, phrase

required_bundle = {
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
assert required_bundle <= set(MANIFEST["artifact_hashes"])
for name, expected in MANIFEST["artifact_hashes"].items():
    path = BUNDLE / name
    assert path.is_file(), name
    assert sha256(path.read_bytes()).hexdigest() == expected, name

crosswalk = json.loads((BUNDLE / "certificate_claim_crosswalk.json").read_text())
artifacts_named = {entry["artifact"] for entry in crosswalk["manuscript_claims"]}
assert artifacts_named <= required_bundle
assert len(artifacts_named) >= 7

print("VERDICT: MANUSCRIPT-TO-CERTIFICATE CLAIMS MATCH DELIVERED ARTIFACTS")
print("crosswalk claims =", len(crosswalk["manuscript_claims"]))
