# Reproducing SC7-01

This is a separate investigation of a proposed 7D extension. It does not modify MTFT.

The main result is a negative physical gate: the smooth commuting twisted 7D extension of the existing degree-three bundles on the spin-circle space has no massless bulk modes in the required M1 sectors. The exact calculations succeed; the proposed spectrum does not.

Files:

- MTFT_SC7_01_Spin_Circle_Investigation.md — derivations, interpretation, limits and references.
- PREREGISTRATION.md — model and success criteria fixed before running the code.
- investigate_spin_circle.py — reproducible exact calculation.
- results.json — machine-readable results and input hashes.
- run.log — final execution output.
- requirements.txt — exact dependency versions used.
- input/mtft-0.31.4.tar.gz — unchanged supplied source archive.
- manifest.json — hashes of the bundled files.

Use Python 3.12 and an environment with the dependencies in requirements.txt. Extract the supplied source archive to a directory, then run:

~~~bash
python investigate_spin_circle.py --source /path/to/mtft-0.31.4 --archive input/mtft-0.31.4.tar.gz --output reproduced_results.json
~~~

The program imports the supplied MTFT source directly. It does not download a different release, call PARI/GP, or require network access. Source-path and runtime-version metadata may differ between machines; the algebraic results should agree.

The 36 primary local-system configurations are exact calculations, with inverse characters also computed for Poincare duality. The spectral-sequence proof in the report covers arbitrary base characters.

The gravitational test assumes its own unwarped Einstein-action ansatz. It does not certify a supersymmetric gravitational embedding.
