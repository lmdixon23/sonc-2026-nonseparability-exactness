# Claim-to-artifact map

| Manuscript claim | Evidence type | Analytic source | Computational artifact |
|---|---|---|---|
| Universal SONC exactness if and only if nonseparability for arbitrary finite interior negative support | analytic proof | Theorem 1.1 and the carrier, degeneration, phase, and duality sections of `paper/main.tex` | baseline and multiphase checks are corroborative |
| Two-negative characterization | analytic proof | standalone two-negative section | `code/verify_general_pair_identities.py` |
| Explicit toric degeneration bound | analytic proof plus certified evaluation | quantitative proposition and corollary | `code/verify_quantitative_degeneracy_bound.py`; `results/quantitative_degeneracy_bound.json` |
| Minimum-cardinality, square, and quadrilateral consequences | analytic proof | consequences section | baseline suite called by `code/run_baseline_checks.py` |
| Supporting discriminant self-intersection | analytic proof | discriminant subsection | no computation is load-bearing |
| Certified three-negative witness | exact and validated computation | witness proposition and computational appendix | `code/reconstruct_certificate_from_exact_data.py`; `results/constructed_witness_certificate/` |
| Complete critical-point classification | exact elimination and Sturm theory | computational appendix | `resultant_factorization.json`, `sturm_certificates.json`, and `sturm_sequences_full.json.gz` |
| Physical root existence and uniqueness | validated interval proof | computational appendix | `krawczyk_certificate.json` |
| Downstream coefficient, phase, Hessian, and value checks | independent interval replay | verification section | `code/replay_downstream_interval_arithmetic.py` |
| Manuscript-to-certificate agreement | source and schema checks | verification section | `code/verify_manuscript_certificate_claims.py`; `code/verify_certificate_schema.py` |
| Checkout-stable certificate bytes | artifact-integrity check | verification section | `code/verify_certificate_portability.py`; `.gitattributes` |
| Canonical public repository integrity | artifact-integrity check | public release boundary | `SHA256SUMS.txt`; `code/verify_sha256_manifest.py` |
