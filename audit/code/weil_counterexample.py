"""Numerical counterexample to the chain Theorem 1.1 + 1.3 + 1.4 of
X.-J. Li, "A proof of the Riemann hypothesis", arXiv:0807.0090v10.

The paper states (Theorem 1.1) trace(T_h) = Delta(h) - h^(0) - h^(1), with
Delta(h) = sum over zeros rho of h^(rho), and h(x) = int_0^oo g(xt) g(t) dt,
so that h^(s) = g^(s) g^(1-s).  Theorems 1.3 and 1.4 then give
trace(T_h) >= 0.  Their proofs (Sections 4 and 5) use only that g is real,
smooth and supported in [1-eps, mu_eps]; they never use g^(0) = 0 or any
other property of the special functions g_{n,eps}.  So if the argument were
valid it would prove, for every such g,

    Delta(h) = sum_rho g^(rho) g^(1-rho)  >=  h^(0) + h^(1) = 2 g^(0) g^(1).   (*)

This script takes a nonnegative C^oo bump g and shows that (*) is false by a
wide margin.  Delta(h) is computed two independent ways:
  (a) directly from the zeros (first N zeros plus a bound for the rest);
  (b) from Weil's explicit formula (primes + archimedean term), in the form
      of Bombieri, Rend. Mat. Acc. Lincei 11 (2000), which uses no zeros.
The two agree to many digits, which checks the whole computation.
"""
import math
import sys

import numpy as np
import sympy as sp

from zeros_cache import zeros

# ---------------------------------------------------------------- test function
# g(x) = G(log x), G(u) = exp(-1/(1 - y^2)), y = 2u/L - 1, supported on [0, L].
L = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
N_ZEROS = int(sys.argv[2]) if len(sys.argv) > 2 else 2000

_u = sp.symbols("u")
_y = 2 * _u / sp.Float(L) - 1
_G = sp.exp(-1 / (1 - _y**2))
G_derivs = [sp.lambdify(_u, sp.diff(_G, _u, j), "numpy") for j in range(5)]

# Gauss-Legendre nodes on [0, L]; the bump is C^oo so this converges fast.
NODES = 6000
xg, wg = np.polynomial.legendre.leggauss(NODES)
U = 0.5 * L * (xg + 1)
W = 0.5 * L * wg
GU = G_derivs[0](U)


def gmellin(s):
    """g^(s) = int_0^oo g(x) x^(s-1) dx = int_0^L G(u) e^(s u) du."""
    return np.sum(W * GU * np.exp(s * U))


def h_of_log(a):
    """h(e^a) = int g(e^a t) g(t) dt = int G(a+b) G(b) e^b db (b in log t)."""
    a = np.atleast_1d(a)
    out = np.empty(a.shape)
    for i, ai in enumerate(a):
        v = U + ai
        m = (v > 0) & (v < L)
        out[i] = np.sum(W[m] * G_derivs[0](v[m]) * GU[m] * np.exp(U[m]))
    return out


# ------------------------------------------------------------ (a) zero side
gammas = np.array([float(g) for g in zeros(N_ZEROS)])
T = gammas[-1]
vals = np.array([abs(gmellin(0.5 + 1j * g)) ** 2 for g in gammas])
delta_zeros = 2.0 * vals.sum()  # rho = 1/2 +- i*gamma

# Tail bound for |Im rho| > T, valid even for zeros off the line.  For
# 0 <= sigma <= 1, integrating by parts k times gives
#   |g^(sigma+it)| <= B_k / |t|^k,  B_k = sum_j C(k,j) int |G^(j)(u)| e^u du,
# and |h^(rho)| <= |g^(rho)| |g^(1-rho)| <= B_k^2 / |t|^(2k).
# With N(t) <= t log t for t >= 10 (far weaker than Riemann-von Mangoldt),
#   sum_{|gamma|>T} |t|^(-2k) <= 2 * 2k * int_T^oo t^(-2k) log t dt.
k = 4
Bk = sum(math.comb(k, j) * np.sum(W * np.abs(G_derivs[j](U)) * np.exp(U)) for j in range(k + 1))
a_ = 2 * k
tail_int = T ** (1 - a_) * (math.log(T) / (a_ - 1) + 1 / (a_ - 1) ** 2)
tail = Bk**2 * 2 * (2 * k) * tail_int

# ------------------------------------------------- (b) explicit-formula side
# Bombieri (2000), for h in C_c^oo(0, oo):
# sum_rho h^(rho) = h^(1) + h^(0) - sum_n Lambda(n) [h(n) + h(1/n)/n]
#                   - (log 4pi + gamma) h(1)
#                   - int_1^oo {h(x) + h(1/x)/x - 2h(1)/x} x dx / (x^2 - 1)
g0, g1 = gmellin(0.0).real, gmellin(1.0).real
h0 = h1 = g0 * g1  # h^(0) = g^(0) g^(1) = h^(1)
h_at_1 = h_of_log(0.0)[0]

prime_sum = 0.0
for n in range(2, int(math.exp(L)) + 1):
    f = sp.factorint(n)
    if len(f) == 1:
        p = next(iter(f))
        ln = math.log(n)
        prime_sum += math.log(p) * (h_of_log(ln)[0] + h_of_log(-ln)[0] / n)

# archimedean integral in the variable a = log x, x in [1, e^L]
xa, wa = np.polynomial.legendre.leggauss(400)
A = 0.5 * L * (xa + 1)
WA = 0.5 * L * wa
X = np.exp(A)
integrand = (h_of_log(A) + h_of_log(-A) / X - 2 * h_at_1 / X) * X / (X**2 - 1) * X
# beyond x = e^L, h(x) = h(1/x) = 0 and the integrand is -2h(1)/(x^2-1)
arch_tail = -2 * h_at_1 * 0.5 * math.log((math.exp(L) + 1) / (math.exp(L) - 1))
arch = np.sum(WA * integrand) + arch_tail
euler_gamma = float(sp.EulerGamma.evalf(30))
delta_explicit = h1 + h0 - prime_sum - (math.log(4 * math.pi) + euler_gamma) * h_at_1 - arch

rhs = h0 + h1
print(f"test function: g(x) = G(log x), C^oo bump on [1, e^{L}] = [1, {math.exp(L):.4f}]")
print(f"g^(0) = {g0:.12f}   g^(1) = {g1:.12f}")
print(f"h^(0) + h^(1)                   = {rhs:.12f}")
print(f"Delta(h), first {N_ZEROS} zeros (T = {T:.2f}) = {delta_zeros:.12f}")
print(f"  bound on the rest (|Im rho| > T)  <= {tail:.3e}")
print(f"Delta(h), explicit formula         = {delta_explicit:.12f}")
print(f"  |zeros - explicit|               = {abs(delta_zeros - delta_explicit):.3e}")
print(f"claimed trace(T_h) = Delta(h) - h^(0) - h^(1) = {delta_zeros - rhs:.12f}")
viol = delta_zeros + tail < rhs
print("Li's Theorems 1.1+1.3+1.4 would force this to be >= 0:",
      "VIOLATED" if viol else "not violated")
