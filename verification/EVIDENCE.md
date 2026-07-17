# Public evidence model

The manuscript distinguishes four kinds of support.

1. **Analytic proof.** The universal characterization, the two-negative specialization, and the stated consequences are established in `paper/main.tex`. Symbolic scripts that check displayed identities are corroborative rather than substitutes for those proofs.
2. **Exact computation.** The explicit witness certificate uses exact polynomial arithmetic, resultant factorization, and Sturm root counts. The canonical records are under `results/constructed_witness_certificate/`.
3. **Validated numerical computation.** Krawczyk inclusion and downstream inequalities use outward-rounded Arb arithmetic through `python-flint`. The required arithmetic environment is recorded in `ENVIRONMENT-LOCK.md`.
4. **Independent reconstruction.** `code/independent_certificate_reconstruction.py` follows a separate implementation path and does not import the primary certificate module. It provides corroboration, not peer review.
5. **Repository integrity.** `SHA256SUMS.txt` binds the canonical public source, code, data, and verification records. `code/verify_sha256_manifest.py` checks the committed checkout before regeneration.

Generated run summaries are convenience outputs. They are ignored by Git and are not independent evidence.

The baseline and multiphase stages invoked by `run_all.sh` are corroborative structural and symbolic checks. They do not validate the universal analytic proof.
The proof-bearing text artifacts use canonical LF bytes, and the compressed Sturm data use a zero gzip timestamp. `code/verify_certificate_portability.py` checks these properties so the committed manifest does not depend on the checkout platform.
