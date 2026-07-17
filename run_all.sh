#!/usr/bin/env bash
# Reproduces the maintained mathematical and certificate claims in claim order.
# No network access is required after dependency installation.
set -euo pipefail
export PYTHONHASHSEED=0
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

PY=python3
"$PY" -c "" 2>/dev/null || PY=python
"$PY" -c "" 2>/dev/null || { echo "No working Python interpreter found"; exit 1; }

echo "== Canonical-file integrity =="
"$PY" code/verify_sha256_manifest.py

echo "== Corroborative structural and symbolic checks =="
"$PY" code/run_baseline_checks.py

echo "== Multiphase consistency checks =="
"$PY" code/run_multiphase_checks.py

echo "== Exact and validated certificate checks =="
"$PY" code/run_certificate_checks.py

echo "== Post-run canonical-file integrity =="
"$PY" code/verify_sha256_manifest.py

echo "VERDICT: ALL MAINTAINED CHECKS PASSED"
