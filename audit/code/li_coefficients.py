"""Li coefficients lambda_n = sum_rho [1 - (1 - 1/rho)^n], n = 1..N.

Computed without zeros, from
    lambda_n = 1/(n-1)! * d^n/ds^n [ s^(n-1) log xi(s) ] at s = 1,
xi(s) = s(s-1)/2 * pi^(-s/2) * Gamma(s/2) * zeta(s), via Cauchy integrals on
a circle around s = 1 (mpmath.diff, method='quad'), so zeta is never evaluated
at its pole.

Li's criterion: RH holds iff lambda_n >= 0 for every n >= 1.  A computation
like this checks finitely many n only.  It is consistent with RH and cannot
prove it; lambda_n > 0 is already known numerically far beyond the range here.
"""
import sys

import mpmath as mp

N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
mp.mp.dps = 60


def log_xi(s):
    return mp.log(s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s))


# Taylor coefficients of log xi at s = 1: log xi(1 + t) = sum a_k t^k
coeffs = mp.taylor(log_xi, 1, N, method="quad", radius=mp.mpf("0.5"))
# s^(n-1) = (1+t)^(n-1); the n-th derivative / (n-1)! of the product at t = 0
# equals n * [t^n] ( (1+t)^(n-1) * log xi(1+t) ).
print(" n   lambda_n")
for n in range(1, N + 1):
    c = sum(mp.binomial(n - 1, j) * coeffs[n - j] for j in range(0, n))
    lam = mp.re(n * c)  # imaginary part is quadrature noise (< 1e-50)
    print(f"{n:2d}   {mp.nstr(lam, 20)}")
# check: lambda_1 = 1 + gamma/2 - log(4 pi)/2
print("lambda_1 closed form:", mp.nstr(1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2, 20))
