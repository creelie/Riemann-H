"""Exact p-adic checks for arXiv:0807.0090v10.

1. Lemma 2.2's local integral, with d lambda_p the additive Haar measure
   giving Z_p mass 1:
       int_{Z_p^*} psi_p(-lambda gamma) d lambda =  1 - 1/p  if |gamma|_p <= 1
                                                  -1/p      if |gamma|_p = p
                                                   0        if |gamma|_p > p
   This is a finite exact sum over residues, so it can be checked exactly.
   It holds, so the error is not in this computation.

2. The step (4.13) -> (4.14) (and (5.8) -> (5.9)) replaces x by x/xi inside
   an integral over C_S = J_S / O_S^*, and claims every xi-term then equals
   the xi = 1 term.  That is legitimate only when the integrand is invariant
   under x -> eta x for eta in O_S^*.  The integrand contains Psi_S(x xi v),
   which is not invariant.  For S = {oo, 2}, the element 2 (and 1/2) lies in
   O_S^*.  Multiplying v by 1/2 moves v_2 from Z_2^* to 2^-1 Z_2^*, and the
   normalised 2-adic factor of the inner integral flips from +1 to -1 (the
   first and second rows of the table above).  The script prints both values.
"""
from fractions import Fraction
import cmath
import math


def frac_part_padic(a, b, p):
    """{a/b}_p: the p-adic fractional part of the rational a/b, as a Fraction."""
    q = Fraction(a, b)
    k = 0
    den = q.denominator
    while den % p == 0:
        den //= p
        k += 1
    if k == 0:
        return Fraction(0)
    # write q = m / (p^k * den) with gcd(den, p) = 1; the fractional part is
    # r / p^k where r = m * den^{-1} mod p^k
    pk = p**k
    m = q.numerator
    r = (m * pow(den, -1, pk)) % pk
    return Fraction(r, pk)


def psi(a, b, p):
    return cmath.exp(2j * math.pi * float(frac_part_padic(a, b, p)))


def unit_integral(gamma_num, gamma_den, p):
    """int_{Z_p^*} psi_p(-lambda gamma) d lambda, gamma = gamma_num/gamma_den.

    psi_p(-lambda gamma) depends only on lambda mod p^k when |gamma|_p <= p^k,
    so the integral is an exact average over residues mod p^k (each residue
    class has measure p^-k)."""
    k = 0
    d = Fraction(gamma_num, gamma_den).denominator
    while d % p == 0:
        d //= p
        k += 1
    k = max(k, 1)
    pk = p**k
    total = 0
    for lam in range(pk):
        if lam % p:
            total += psi(-lam * gamma_num, gamma_den, p)
    return total / pk


def expected(gamma_num, gamma_den, p):
    q = Fraction(gamma_num, gamma_den)
    v = 0
    d = q.denominator
    while d % p == 0:
        d //= p
        v += 1
    if v == 0:
        return 1 - 1 / p
    if v == 1:
        return -1 / p
    return 0.0


ok = True
for p in (2, 3, 5, 7):
    for num in (1, 2, 3, 5, 7, 11, -4, 9):
        for den in (1, 3, p, p * 5, p**2, p**3 * 7):
            got = unit_integral(num, den, p)
            exp = expected(num, den, p)
            if abs(got - exp) > 1e-12:
                ok = False
                print("MISMATCH", p, num, den, got, exp)
print("Lemma 2.2 local integrals, p in {2,3,5,7}: all match" if ok else "Lemma 2.2: mismatch")

# ---- non-invariance of the integrand used in (4.13) -> (4.14)
# Normalised by the measure of Z_p^* (which is 1 - 1/p) the table becomes
# +1 / (-1/(p-1)) / 0.  For p = 2 that is +1 / -1 / 0.
p = 2
for label, (num, den) in (("v_2 = 1", (1, 1)), ("v_2 = 1/2 (= v_2 times 2^-1, 2 in O_S^*)", (1, 2))):
    val = unit_integral(num, den, p) / (1 - 1 / p)
    print(f"normalised int_(Z_2^*) psi_2(-x_2 v_2) d*x_2 at {label}: {val.real:+.6f}")
print("The x-integrand changes when v is moved by an element of O_S^*, so the")
print("xi-terms of (4.13) are not all equal, and (4.14)-(4.15) do not follow.")
