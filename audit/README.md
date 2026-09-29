# Audit of X.-J. Li, "A proof of the Riemann hypothesis" (arXiv:0807.0090v10)

Authors of this audit: Deep Bhattacharjee and Priyabrata Mandal (corresponding author).

## Verdict

The paper does not prove the Riemann Hypothesis. Its two key positivity steps,
Theorems 1.3 and 1.4, rely on the same invalid change of variables. That step
cannot be repaired with computation, because the missing statement is about
every zero and every test function, and code can only check finitely many.
The code in `code/` shows that the argument proves too much: applied to an
ordinary smooth test function it yields an inequality that is false by a
factor of about 26,000.

No unconditional proof of RH is given here. None is known.

## What the paper claims

Li uses Connes' semilocal adelic setting. S is a finite set of places, C_S =
J_S / O_S^* and T_h is a convolution operator with two cutoffs. The chain is:

1. **Theorem 1.1** (cited from Meyer and Li): trace(T_h) = Δ(h) − ĥ(0) − ĥ(1),
   where Δ(h) = Σ_ρ ĥ(ρ) and h(x) = ∫ g(xt) g(t) dt, so ĥ(s) = ĝ(s) ĝ(1−s).
2. **Theorem 1.2**: special functions g_{n,ε} have ĥ(0) = ĥ(1) = 0 and Δ(h_{n,ε}) → 2λ_n, where λ_n is the n-th Li coefficient.
3. **Theorem 1.3**: the trace of T_h on E_S(Q_Λ^⊥) is 0.
4. **Theorem 1.4**: the trace of T_h on E_S(Q_Λ) is ≥ 0.
5. Therefore λ_n ≥ 0 for all n, and RH follows from Li's criterion.

## Claim-by-claim status

| Item | Status | How checked |
|---|---|---|
| Lemma 2.1 (fundamental domain I_S) | Correct | Standard (Tate's thesis) |
| Lemma 2.2 (local integrals over Z_p^*) | Correct | `padic_check.py`, exact finite sums for p = 2, 3, 5, 7 |
| Lemma 2.3, Plancherel | Plausible, not the issue | Read |
| Theorem 1.1 (trace = Δ(h) − ĥ(0) − ĥ(1)) | Cited, taken as given | The right side equals minus the prime and archimedean terms of Weil's explicit formula. `weil_counterexample.py` checks the explicit formula to about 2·10⁻¹² |
| Theorem 1.2 (Δ(h_{n,ε}) → 2λ_n) | Plausible (follows Bombieri 2000) | The target values λ_n are computed in `li_coefficients.py` |
| Theorem 1.3 | **Proof invalid** at (4.13)–(4.15) | See below, `padic_check.py` |
| Theorem 1.4 | **Proof invalid** at (5.8)–(5.10) | Same error |
| Theorems 1.1 + 1.3 + 1.4 together | **False** for a smooth bump g | `weil_counterexample.py` |
| Theorem 1.5 (RH) | Not established | Depends on 1.3 and 1.4 |

## The broken step

Section 4 ends with an absolutely convergent sum over the infinite group
O_S^* (4.13). Li substitutes x → ξ⁻¹x, then u → uξ, and says every term of
(4.14) is "the same number". An infinite sum of one repeated number can only be
finite if that number is 0, so the paper concludes (4.15) and trace = 0.
Section 5 makes the same move at (5.8) → (5.9) → (5.10).

The substitution x → ξ⁻¹x inside ∫_{C_S} … d^×x is valid only when the
integrand is a function on C_S, meaning it is unchanged when x is multiplied
by any η ∈ O_S^*. The integrand contains Ψ_S(xξv), and that is not invariant.
In Lemmas 4.3 and 4.4 the x-integral was made concrete by using the
fundamental domain I_S. Replacing x by ξ⁻¹x moves it to the different domain
ξ⁻¹I_S, so the ξ-terms are not equal. Their sum is just the original integral
over all of J_S.

`padic_check.py` shows this concretely. For S = {∞, 2}, 1/2 is in O_S^*.
Multiplying v by 1/2 flips the normalised 2-adic factor ∫_{Z_2^*} ψ_2(−x v) d^×x
from +1 to −1.

## Why the argument must be wrong: it proves too much

Sections 4 and 5 use only three facts about g: it is real, smooth, and
supported in [1−ε, μ_ε]. The special properties of g_{n,ε}, including
ĝ(0) = 0, never appear there. So if Theorems 1.3 and 1.4 were right, then with
Theorem 1.1 they would give, for every such g,

  Σ_ρ ĝ(ρ) ĝ(1−ρ)  ≥  ĥ(0) + ĥ(1) = 2 ĝ(0) ĝ(1).   (*)

Take g(x) = G(log x), where G is the C^∞ bump exp(−1/(1−y²)) with
y = 2u/3 − 1 on [0, 3]. The support is [1, e³]; take S to contain all primes
up to e³ and ε small. `weil_counterexample.py` gives:

```
test function: g(x) = G(log x), C^oo bump on [1, e^3.0] = [1, 20.0855]
g^(0) = 0.665990724252   g^(1) = 3.550155164339
h^(0) + h^(1)                   = 4.728740818213
Delta(h), first 50 zeros (T = 143.11) = 0.000182557559
  bound on the rest (|Im rho| > T)  <= 1.242e-07
Delta(h), explicit formula         = 0.000182557561
  |zeros - explicit|               = 1.739e-12
claimed trace(T_h) = Delta(h) - h^(0) - h^(1) = -4.728558260654
Li's Theorems 1.1+1.3+1.4 would force this to be >= 0: VIOLATED
```

The zero sum and the explicit formula are computed independently, one from the
zeros and one from the primes, and they agree. (*) fails by more than four
orders of magnitude. So the combined statement trace(T_h) ≥ 0 from
Theorems 1.3 and 1.4 is false for this g. Since their proofs cannot tell this g
apart from g_{n,ε}, they are not valid proofs for g_{n,ε} either.

The trace is not always nonnegative, which is exactly what Theorems 1.3 and 1.4
assert. When ĥ(0) = ĥ(1) = 0, trace ≥ 0 is Weil's positivity criterion, which
is equivalent to RH. Weil positivity for h supported in a fixed interval is a
statement about finitely many primes, but it still involves all zeros.

## History

Version 1 appeared in July 2008. It drew public criticism, notably from
A. Connes, and has been revised many times since. Version 10 (October 2025) is
listed under math.GM and has not been accepted by the community. This audit concerns v10
only.

## The other uploads

- **Eswaran, "The Final and Exhaustive Proof of the Riemann Hypothesis from First Principles" (2018).**
  This is not a proof. It reduces RH to Littlewood's criterion L(N) = O(N^{1/2+ε}).
  It then argues that λ(n) "behaves like coin tosses", using equal densities of
  λ = ±1, non-periodicity, and matching finite-range statistics. Equal density
  is equivalent to the prime number theorem, which gives only L(N) = o(N). A
  deterministic ±1 sequence can have equal densities, never repeat, and match
  the paper's statistical tests on any finite range, yet still have partial sums as large as
  N^{0.99}. Random-walk laws give almost-sure statements about random
  sequences. They say nothing about one fixed sequence. The argument needs a
  bound on L(N) that is equivalent to RH, and it assumes that bound.
- **Knausgård, arXiv:2609.33043.** This is about distinct zeros (83.69%), not RH.
  It is covered in the "Critical line proportion paper" thread.
- **Bhattacharjee, "Liouville Correlations and Zero-Free Half-Planes" (preprint 232162) and its code.**
  The preprint does not claim RH. It states equivalences and conditional
  results. We reran the supplementary `liouville_data.py` for N ≤ 10⁶, and it
  reproduces Table 1 and all three internal checks. `liouville_independent.py`
  recomputes L(N) with a separate sieve for all six rows up to 2·10⁷, and every
  value matches: −14, −94, −288, −530, −842, −4510.

## Code

All scripts are in `code/` and need Python 3 with numpy, mpmath and sympy.

| Script | What it does | Runtime |
|---|---|---|
| `weil_counterexample.py [L] [N]` | The counterexample to (*). Computes Δ(h) from zeros (with a tail bound valid for zeros off the line) and from Weil's explicit formula | seconds, once zeros are cached |
| `zeros_cache.py N` | Caches the first N zeros via `mpmath.zetazero` | minutes for 2000 |
| `padic_check.py` | Exact checks of Lemma 2.2, and the non-invariance behind (4.13)→(4.14) | < 1 s |
| `li_coefficients.py N` | λ_1..λ_N from derivatives of log ξ at s = 1, with no zeros used. All positive, consistent with RH | ~1 min for N = 30 |
| `liouville_independent.py` | Independent L(N) for Table 1 of the Liouville preprint | ~5 s |

Saved outputs are the `*.txt` files next to each script.

The Li-coefficient output is a sanity check, not evidence toward a proof.
Positivity of λ_n for finitely many n is already known far beyond n = 30, and
no finite range implies RH.
