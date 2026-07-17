# Status ledger

Version: v0.9-rc2, 2026-07-17

This file records the status and primary evidence for the public mathematical and computational claims.
The manuscript remains subject to ordinary scholarly review.

| ID | Claim | Status | Primary evidence | Scope or limitation |
|---|---|---|---|---|
| C001 | SONC exactness is equivalent to nonseparability for finite interior negative support. | PROVED | `paper/main.tex`, Theorem 1.1 | Assumes full-dimensional positive support and strict interiority of every negative exponent. |
| C002 | The two-negative converse has the standalone scalar proof stated in the manuscript. | PROVED | `paper/main.tex`; `code/verify_general_pair_identities.py` | The script checks the displayed identities; the proof is analytic. |
| C003 | The explicit degeneration bound has the stated numerical benchmark. | CERTIFIED | `code/verify_quantitative_degeneracy_bound.py`; `results/quantitative_degeneracy_bound.json` | Sufficient and conservative, not optimal. |
| C004 | The minimum-cardinality, square, and quadrilateral consequences follow as stated. | PROVED | `paper/main.tex`; baseline check suite | The scripts are corroborative. |
| C005 | The explicit three-negative signomial is globally nonnegative and has exactly two zeros. | CERTIFIED | `results/constructed_witness_certificate/`; reconstruction and replay scripts | Depends on the arithmetic semantics in `verification/ENVIRONMENT-LOCK.md`. |
| C006 | A second implementation reproduces the exact elimination and interval classification without importing the primary certificate module. | REPRODUCED | `code/independent_certificate_reconstruction.py` | Corroborates the certificate; it is not a separate peer review. |
| C007 | The committed certificate bundle has checkout-stable canonical text bytes and deterministic compressed Sturm data. | VERIFIED | `code/verify_certificate_portability.py`; `.gitattributes` | Environment provenance in the manifest is expected to change when the bundle is intentionally regenerated elsewhere. |
| C008 | The canonical public repository files match the release manifest. | VERIFIED | `SHA256SUMS.txt`; `code/verify_sha256_manifest.py` | Generated ignored run summaries and private workspaces are outside the manifest. |

Current main-result status: proved in the manuscript. The preprint has not undergone formal journal peer review.
