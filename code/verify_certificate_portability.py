#!/usr/bin/env python3
"""Verify checkout-stable line endings and deterministic compressed certificate data."""
from __future__ import annotations

from pathlib import Path
import json
import struct

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BUNDLE = ROOT / "results" / "constructed_witness_certificate"

text_artifacts = sorted(BUNDLE.glob("*.json")) + [BUNDLE / "README.md"]
for path in text_artifacts:
    data = path.read_bytes()
    assert b"\r\n" not in data, f"noncanonical CRLF bytes: {path.relative_to(ROOT)}"
    assert data.endswith(b"\n"), f"missing final newline: {path.relative_to(ROOT)}"

compressed = (BUNDLE / "sturm_sequences_full.json.gz").read_bytes()
assert compressed[:3] == b"\x1f\x8b\x08", "invalid gzip header"
mtime = struct.unpack("<I", compressed[4:8])[0]
assert mtime == 0, f"gzip mtime must be zero, found {mtime}"

manifest = json.loads(
    (BUNDLE / "certificate_manifest.json").read_text(
        encoding="utf-8"
    )
)
environment = manifest["environment"]
assert environment["python_major_minor"] == "3.13"
assert "platform" not in environment
assert environment["serialization"] == (
    "platform-independent UTF-8 LF; deterministic gzip mtime=0"
)

attributes = (ROOT / ".gitattributes").read_text(encoding="utf-8")
for rule in (
    "*.py text eol=lf",
    "*.tex text eol=lf",
    "*.json text eol=lf",
    "*.md text eol=lf",
    "*.txt text eol=lf",
):
    assert rule in attributes, rule

core = (HERE / "constructed_witness_exact_core.py").read_text(encoding="utf-8")
reconstruction = (HERE / "reconstruct_certificate_from_exact_data.py").read_text(encoding="utf-8")
assert 'mtime=0' in core
assert 'newline="\\n"' in core
assert 'newline="\\n"' in reconstruction

print("VERDICT: CERTIFICATE PORTABILITY AND DETERMINISTIC-BYTE POLICY VERIFIED")
print("text_artifacts =", len(text_artifacts))
print("gzip_mtime =", mtime)
