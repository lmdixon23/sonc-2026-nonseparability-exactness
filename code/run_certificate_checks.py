#!/usr/bin/env python3
"""Run the maintained proof-bearing certificate and manuscript-consistency checks."""
from pathlib import Path
from runner_utils import run_scripts

HERE = Path(__file__).resolve().parent
SCRIPTS = [
    "verify_environment_lock.py",
    "verify_quantitative_degeneracy_bound.py",
    "reconstruct_certificate_from_exact_data.py",
    "verify_certificate_schema.py",
    "replay_downstream_interval_arithmetic.py",
    "generate_computational_appendix.py",
    "verify_manuscript_certificate_claims.py",
    "independent_certificate_reconstruction.py",
]
run_scripts(
    HERE,
    SCRIPTS,
    timeout_seconds=900,
    output_json=HERE.parent / "results" / "certificate_checks.json",
)
print("VERDICT: ALL CERTIFICATE CHECKS PASSED")
print("scripts =", len(SCRIPTS))
