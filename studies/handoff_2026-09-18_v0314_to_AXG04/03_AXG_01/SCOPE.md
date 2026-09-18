# AXG-01 scope and conventions

Date: 2026-09-17. Additive study of the axion/radius/Higgs proposal in the
preceding conversation, using the supplied MTFT 0.31.4 candidate.

The bounded questions are:
1. Does the full four-dimensional anomaly polynomial factor through the two
   proposed Stückelberg charge rows, including mixed anomalies involving the
   surviving Abelian generators?
2. What mass-rank and charge-lattice conclusions follow without inventing an
   axion scale or kinetic metric?
3. What Higgs and matter operators are allowed by these axion charges?
4. Does the minimal classical extension alter the inherited radius runaway
   or select Higgs flavor directions?

The candidate has stack order (c,L,a,b,d), fundamental central charges +1,
two dimensionless 2π-periodic axions with a -> a-K lambda, and
K=[[0,-2,1,1,0],[3,-4,0,0,1]]. Gauge connections transform A -> A+d lambda.
Both axion and gauge kinetic matrices are positive definite. A common stack
coupling is an explicit conditional assumption, as in the existing model.

Four-dimensional local anomaly cancellation is not six-dimensional or full
global consistency. No absent CW-01 parent, potential, axion kinetic matrix,
absolute scale, or counterterm matching condition will be supplied by a fit.
Optional equal-kinetic controls and negative-vacuum-energy controls are
marked as additional assumptions rather than MTFT predictions.

Scripts have separate responsibilities: anomaly_audit.py, vacuum_audit.py,
and operator_audit.py. No package code or earlier correction register is
modified. This scope records the current investigation, not a claim of a
blind preregistered discovery of the previously displayed charge matrix.
