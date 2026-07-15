#!/usr/bin/env bash
# Reproduces the maintained mathematical and certificate claims in claim order.
# No network access is required after dependency installation.
set -euo pipefail
export PYTHONHASHSEED=0

PY=python3
"$PY" -c "" 2>/dev/null || PY=python
"$PY" -c "" 2>/dev/null || { echo "No working Python interpreter found"; exit 1; }

echo "== Analytic and structural checks =="
"$PY" code/run_baseline_checks.py

echo "== Multiphase consistency checks =="
"$PY" code/run_multiphase_checks.py

echo "== Exact and validated certificate checks =="
"$PY" code/run_certificate_checks.py

echo "VERDICT: ALL MAINTAINED CHECKS PASSED"
