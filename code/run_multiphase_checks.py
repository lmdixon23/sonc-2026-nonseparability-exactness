#!/usr/bin/env python3
"""Run the maintained multiphase symbolic and numerical consistency checks."""
from pathlib import Path
from runner_utils import run_scripts

HERE = Path(__file__).resolve().parent
SCRIPTS = [
    "verify_multiphase_symbolic.py",
    "verify_multiphase_square_interval.py",
    "verify_multiphase_helly_consistency.py",
    "independent_multiphase_replay.py",
]
run_scripts(
    HERE,
    SCRIPTS,
    timeout_seconds=180,
    output_json=HERE.parent / "results" / "multiphase_checks.json",
)
print("VERDICT: ALL MULTIPHASE CHECKS PASSED")
print("scripts =", len(SCRIPTS))
