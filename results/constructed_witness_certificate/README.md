# Constructed-witness certificate bundle

This directory is generated from the exact benchmark polynomials

P = 1 + X^5 + Y^5 + 4096 X^5 Y^5,
Q = X Y (1 + X^3 + Y^3).

Regenerate it from the repository root with:

    python code/reconstruct_certificate_from_exact_data.py

The bundle contains the original critical equations, complete original-system
resultant factorization, exact Sturm certificates and compressed full sequences,
the linear subresultant and physical-root classification, an original-system Arb
Krawczyk inclusion, the strict critical-value comparison, final witness intervals,
and a manuscript claim crosswalk. certificate_manifest.json hashes every bundle
artifact and every load-bearing reconstruction source file.

replay_downstream_interval_arithmetic.py is intentionally narrower: it checks the
downstream interval arithmetic from certified boxes and does not redo elimination,
Sturm isolation, subresultants, or Krawczyk existence.
