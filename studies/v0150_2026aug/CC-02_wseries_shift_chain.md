# CC-02 (append-only) — the Dirichlet series of w_n

**Filed:** 2026-08-12, v0.15.0 integration. **Class:** [Pr], E2-certified (three routes, no shared steps).
**Source:** `W1_weil_compression_study.md` §1 (byte-preserved in this directory).

## Corrected statement

    sum_{n>=1} w_n n^{-s} = F(s+1) = -zeta(s) * zeta'(s+1)

where w_n = sum_{d|n} (log d)/d and F(s) = -zeta(s-1) zeta'(s) is the settled
Arithmetica Generale series (AG Pr 4.1.3).

## What the corpus printed (both wrong, both preserved)

- Paper 1, Prop. 1.5, eq. (4):  sum w_n n^{-s} = -zeta'(s+1).
- Arithmetica Generale, Pr 4.1.4:  sum w_n n^{-s} = -zeta'(s)/zeta(s-1),
  "equivalently F(s+1)/zeta(s)".

## Adjudication (three routes sharing no steps — E2)

1. Direct sieve of w_n = sum_{d|n} (log d)/d, n <= 3e5, versus the convolution
   sieve w = (Lambda/id) * sigma_{-1} (equivalently w = Lambda_1 * 1,
   Lambda_1(n) = (log n)/n): max pointwise difference **7.1e-15**.
2. Numeric series at s = 3: partial sum 0.0828352629; -zeta(3) zeta'(4) =
   0.0828352629, |diff| = **5.2e-12**, exactly the size of the estimated tail.
   The printed alternatives miss by 1.39e-2 and 3.76e-2 respectively — excluded.
3. Consistency with the settled F: the same partial sum equals the F(4) partial
   sum identically (w_n n^{-3} = f(n) n^{-4}), and F(4) matches -zeta(3) zeta'(4)
   to the same precision.

## Diagnosis

-zeta'(s+1) is the Dirichlet series of the *summand* Lambda_1(n) = (log n)/n.
Since w_n = (Lambda_1 * 1)(n), the divisor sum multiplies by zeta(s); Paper 1's
own proof line contains the zeta(s) factor and the final equality drops it (the
same pathology as CC-01: a printed conclusion contradicting its own intermediate
line). AG Pr 4.1.4's "-zeta'(s)/zeta(s-1)" is additionally inconsistent with its
companion clause "F(s+1)/zeta(s)" (which itself is wrong — no division). Per
protocol both statements receive append-only annotations; nothing is rewritten.

## Downstream exposure

Expected nil. Paper 18's explicit-formula pipeline and the Laplace-ensemble
kernel use F(s) = -zeta(s-1) zeta'(s), which is correct and unaffected. The
corpus grep for uses of "-zeta'(s+1)" as the w-series found no downstream
dependence.

## Bridge

The corrected series places zeta' at the shifted argument s+1 inside the
w-generating function — the object of the Speiser-Hadamard lab, and of the
source paper's xi' extension (its Remark 7.3). See W1 study §5(b) (W3 proposal).
