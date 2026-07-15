#!/usr/bin/env python3
"""Verify the Python and external-tool release environment."""
from __future__ import annotations

import shutil
import subprocess
import sys
import sympy
import numpy
import scipy
import mpmath
import flint

expected = {
    "sympy": "1.14.0",
    "python-flint": "0.9.0",
    "numpy": "2.4.4",
    "scipy": "1.17.1",
    "mpmath": "1.3.0",
}
actual = {
    "sympy": sympy.__version__,
    "python-flint": flint.__version__,
    "numpy": numpy.__version__,
    "scipy": scipy.__version__,
    "mpmath": mpmath.__version__,
}
assert actual == expected, {"expected": expected, "actual": actual}
assert sys.version_info[:2] == (3, 13), sys.version
for tool in ("pdflatex", "pdfinfo", "pdffonts", "pdftotext", "zip"):
    assert shutil.which(tool), tool
version = subprocess.run(["pdflatex", "--version"], capture_output=True, text=True, timeout=30, check=True).stdout.splitlines()[0]
print("VERDICT: ENVIRONMENT LOCK VERIFIED")
print("python =", sys.version.split()[0])
for key in expected:
    print(f"{key} = {actual[key]}")
print("pdflatex =", version)
