#!/usr/bin/env python3
"""Run the maintained analytic and corroborative baseline checks."""
from pathlib import Path
from runner_utils import run_scripts

HERE = Path(__file__).resolve().parent
SCRIPTS = [
    "verify_minimal_obstruction.py",
    "verify_diagonal_parallelogram.py",
    "verify_square_converse_formulas.py",
    "verify_general_pair_identities.py",
    "verify_square_boundary_classification.py",
    "verify_symmetric_helly_obstruction.py",
    "verify_symmetric_helly_maxima.py",
    "verify_toric_degenerations.py",
    "independent_baseline_replay.py",
]
run_scripts(
    HERE,
    SCRIPTS,
    timeout_seconds=180,
    output_json=HERE.parent / "results" / "baseline_checks.json",
)
print("VERDICT: ALL BASELINE CHECKS PASSED")
print("scripts =", len(SCRIPTS))
