#!/usr/bin/env python3
"""Reconstruct the complete exact/interval witness certificate from P and Q."""
from pathlib import Path
import json
from constructed_witness_exact_core import build_certificate_bundle

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BUNDLE = ROOT / "results" / "constructed_witness_certificate"
SUMMARY = ROOT / "results" / "constructed_witness_certificate.json"

summary = build_certificate_bundle(BUNDLE)
SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print("VERDICT: COMPLETE CERTIFICATE RECONSTRUCTED FROM EXACT DATA")
print("bundle =", BUNDLE)
print("summary =", SUMMARY)
print("alpha =", summary["alpha_interval"])
