# Riemann-H

Research notes on proportions of nontrivial zeros of the Riemann zeta function.

## Contents

- `paper/distinct-zeros-retuning.tex`: *A parameter refinement of the lower bound for
  distinct zeros of the Riemann zeta function* (draft, math.NT, not submitted).
  It shows `liminf N_d(T)/N(T) >= 62359683640669/74504434380000 = 0.83699291404...`
  by re-choosing the free parameters in Knausgard, arXiv:2609.33043, with the same
  imported inputs, and shows that the block length used is the largest admissible one.
- `code/constants.py`: exact rational verification of every constant in the note.

## Build

    python3 code/constants.py
    cd paper && pdflatex distinct-zeros-retuning.tex && pdflatex distinct-zeros-retuning.tex

## Scope

These results concern proportions of zeros. A proportion of 1 (a density-one
statement) would not by itself prove the Riemann Hypothesis, since a density-zero
set of zeros could still lie off the critical line.
