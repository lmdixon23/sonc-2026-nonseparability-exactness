# Author independent check

Run this check before opening `results/constructed_witness_certificate/` after any change to the witness code.

1. Starting only from
   $P(X,Y)=1+X^5+Y^5+4096X^5Y^5$ and
   $Q(X,Y)=XY(1+X^3+Y^3)$, derive the two logarithmic critical equations.
2. Symmetrize the off-diagonal system with $s=X+Y$ and $p=XY$.
3. Recompute the resultant factorization and count positive roots of the diagonal factor, $R_{12}$, and $R_{35}$.
4. Recover $p$ from the linear subresultant, classify the three positive $s$ candidates by $p>0$ and $s^2-4p>0$, and certify the surviving box with a separately written Krawczyk operator.
5. Compare the off-diagonal ratio enclosure with the unique positive diagonal critical value.
6. Only after recording the independently obtained intervals, open the committed certificate and compare hashes and enclosures.

The maintained implementation is `code/independent_certificate_reconstruction.py`. It is a separate implementation that does not import the primary certificate module. It provides corroborating reconstruction evidence and does not replace the analytic proof or the proof-bearing certificate.
