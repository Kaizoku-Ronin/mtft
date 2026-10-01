# Matter-02 mode lab

Open `MTFT_Matter02_Mode_Lab.html` in a modern browser. It is a self-contained offline file; the primary-reference link is optional. No dependency installation or network access is needed to explore the lab.

The lab shows four relationships:

- Exact h0, h1, h2 and charged-Phi dimensions for the specified cusp representative, including the distinction between a net index and all available modes.
- The fat-point quotient for the q=1 line. Five local classes survive; the global image removes z. Three selected surviving classes enter the coupling witness.
- An invertible three-by-three holomorphic residue minor. Its determinant is -1 when the leading coordinate coefficient is fixed to1, even as the next coefficient changes. This proves rank at least3 in one18-by44 holomorphic pairing; it does not provide physical masses or select light families.
- Hypothetical axion-shift rows. Their input kernel conditions illustrate why exact flux projections alone cannot prove photon masslessness.

The frozen dimensions, jet determinant, canonical coefficients, matrix dimensions and source hashes are embedded in the file. They match `../cohomology/cohomology_exact.json` and `../yukawa/yukawa_exact.json`. The parent condition and open physical gates are visible in the lab.

Run the standard-library JavaScript check:

```bash
node verify_lab.cjs MTFT_Matter02_Mode_Lab.html lab_runtime_checks.json
```

The verifier checks138 numerical, coefficient-extraction, rank, Riemann-Roch, duality, state-update, selection and theme gates. A minimal mock DOM exercises the actual runtime script. **No actual browser layout test was performed**; successful mock-DOM execution is not a layout or cross-browser certificate.

The light and dark themes, responsive breakpoints, semantic controls, SVG text labels and scope statements are authored in the file. Root integration should retain the explicit QA limitation.
