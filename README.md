# Nonseparability characterizes SONC exactness for interior signed supports

![verify](https://github.com/lmdixon23/sonc-2026-nonseparability-exactness/actions/workflows/verify.yml/badge.svg)

Code, exact certificate data, and reproducibility material for the accompanying preprint.

## Main result

Let $A_+,A_-\subset\mathbb R^n$ be finite and disjoint, with $A_+$ full dimensional and
with $A_-$ contained in the interior of the convex hull of $A_+$. The paper proves that the signed support
$(A_+,A_-)$ is SONC-exact if and only if it is nonseparable in the sense of Feliu, Ferrer, and Telek.

The universal theorem is analytic. Computation is proof-bearing only for the explicit three-negative
witness and the stated quantitative benchmark.

## Reproduce the maintained claims

Use Python 3.13 and the locked package versions:

```bash
python -m pip install -r requirements-lock.txt
bash run_all.sh
```

The final expected line is:

```text
VERDICT: ALL MAINTAINED CHECKS PASSED
```

Run scripts from the repository root. No network access is required after dependency installation.
The suite first verifies `SHA256SUMS.txt`, then runs the corroborative structural checks, multiphase checks, and proof-bearing certificate checks.

## Claim-to-script map

| Claim | Primary implementation | Evidence |
|---|---|---|
| Structural consequences and the two-negative identities | `code/run_baseline_checks.py` | exact symbolic checks and independent replay |
| Multiphase formulas and consistency checks | `code/run_multiphase_checks.py` | symbolic and interval checks with an independent replay |
| Quantitative degeneration benchmark | `code/verify_quantitative_degeneracy_bound.py` | directed Arb interval evaluation |
| Explicit three-negative witness | `code/reconstruct_certificate_from_exact_data.py` | exact elimination, Sturm isolation, and Krawczyk certification |
| Downstream contact, coefficient, phase, Hessian, and value checks | `code/replay_downstream_interval_arithmetic.py` | separately written outward-rounded Arb replay |
| Independent certificate reconstruction | `code/independent_certificate_reconstruction.py` | second implementation that does not import the primary certificate module |
| Manuscript-to-certificate agreement | `code/verify_manuscript_certificate_claims.py` | exact source and certificate consistency checks |
| Certificate portability | `code/verify_certificate_portability.py` | canonical LF text and deterministic gzip bytes |
| Canonical public repository integrity | `code/verify_sha256_manifest.py` | SHA-256 binding of public source, code, data, and verification records |

The proof-bearing certificate data are under `results/constructed_witness_certificate/`. The detailed
claim crosswalk is `verification/CLAIM-TO-ARTIFACT-MAP.md`.

## Independent reconstruction and integrity

`code/independent_certificate_reconstruction.py` follows a separate implementation path and does not import the primary certificate module. `SHA256SUMS.txt` binds the canonical public source, code, data, and verification records; `code/verify_sha256_manifest.py` fails on missing, altered, duplicated, malformed, or unsafe manifest entries. Generated run summaries are ignored and are not part of the citation record.

## Repository map

- `paper/`: canonical manuscript source and computational appendix
- `code/`: maintained mathematical and certificate checks
- `results/constructed_witness_certificate/`: exact and validated witness certificate
- `verification/`: public status, evidence, environment, and claim-to-artifact documentation

## Scope

The theorem assumes that every negative exponent lies in the interior of the positive Newton polytope.
The quantitative estimate is sufficient rather than optimal. The explicit certificate supports the
constructed example; the universal characterization is proved analytically in the manuscript.

## Cite

See `CITATION.cff`. Versioned PDFs and arXiv source archives are distributed through GitHub Releases. Add the arXiv identifier and archival DOI after posting.
