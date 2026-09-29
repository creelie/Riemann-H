# Riemann-H

Research on the nontrivial zeros of the Riemann zeta function.

## Paper

`paper/zeta-zeros-unconditional.tex`: *Unconditional bounds for zeros of the Riemann zeta
function and limits of recent approaches to the Riemann Hypothesis*, by Deep Bhattacharjee
and Priyabrata Mandal (corresponding author). Draft, math.NT, not submitted.

It does **not** prove the Riemann Hypothesis. It proves, without assuming RH:

1. `liminf N_d(T)/N(T) >= 62359683640669/74504434380000 = 0.83699291404...` for distinct
   zeros, by re-choosing the parameters of Knausgard (arXiv:2609.33043) with the same
   inputs, and shows the block length used is the largest the argument admits.
2. An obstruction for trace-moment methods on the Gabor-compressed Weil matrix: an upper
   bound on an even trace moment of order b forces every extremal, weakly isolated zero to
   satisfy `beta <= 1/2 + 1/(b lambda) + o(1)`. The fourth- and sixth-moment inputs of the
   recent 0.7962 and density-one preprints are of this form, so those claims are unproven
   (the density-one preprint has been retracted by its authors).
3. The invalid step in X.-J. Li's claimed proof of RH (arXiv:0807.0090v10), with a
   numerical counterexample to the combined positivity statement it asserts.

A proportion of 1 (a density-one statement) would not by itself prove RH, since a
density-zero set of zeros could still lie off the critical line.

## Code

- `code/constants.py`: exact rational check of the distinct-zeros constant and optimality.
- `code/offline_pair.py`: complex Gabor-Poisson identity, the off-line pair block, and the
  thresholds of the obstruction (output in `code/offline_pair.txt`, about 5 minutes).
- `code/li_audit/`: the Li audit (`AUDIT.md`), `weil_counterexample.py`, `padic_check.py`,
  `li_coefficients.py`, `liouville_independent.py`, and the zero cache.

Requires Python 3 with numpy, mpmath and sympy.

## Build

    cd paper && pdflatex zeta-zeros-unconditional.tex && pdflatex zeta-zeros-unconditional.tex
