# Run notes

The protocol was written before the new bound and sample calculations.
Before evaluation, independent review supplied the active-column unitary
remainder, recorded as a protocol amendment. Model/horizons/orders/budget
were unchanged.

The primary calculation completed successfully. Independent review verified
the same certificate using Gaussian rational matrices. Its first symbolic
toy run had a false structural-equality mismatch; simplifying the polynomial
difference resolved it. No failed physical or approximation hypothesis was
removed.

Following review, the encoded-matrix reader was hardened to reject floating
and Boolean JSON integers. Parent-model reconstruction was added as a
diagnostic gate. The final primary certificate was regenerated to capture
the final code hashes and reproduced the same bounds. The review's earlier
integer-intake suggestion is resolved in this version.

Historical matrix arithmetic uses NumPy2.3.5/SciPy1.17.0, one OpenBLAS thread.
The certificate target is the stored finite dyadic matrix, so BLAS differences
in rebuilding BQ do not change the exact replay target. Numerical reconstruction
is checked within1e-12; it is exactly equal on the recorded runtime.

No prior study payload or MTFT source file was changed by this experiment.
