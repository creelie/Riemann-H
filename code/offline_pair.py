"""Numerical companion to the off-line pair obstruction (Section 3 of the paper).

Setting (Yang--Yang, "More than 0.7962 ...", Section 2, following [C26]):
  phi(u) = rho((L/2 - |u|)/w)   (C^3 taper, ramp width w, support [-L/2, L/2]),
  phihat(z) = int phi(u) e^{i z u} du,   Phi = (phi^2)^,
  tau_k = T0 + 2 pi k / L,   u_rho = (phihat(gamma_rho - tau_k))_k,
  G = sum_rho m_rho u_rho u_rho^T,   Gtilde = G / L.
For a zero rho = 1/2 + eta + i gamma we write gamma_rho = gamma - i eta.

The script checks, with a concrete taper:
  1. the Gabor--Poisson identity sum_k phihat(a - tau_k) phihat(b - tau_k) = L Phi(a - b)
     for COMPLEX a, b (analytic continuation of [C26, Lem. 2.2]);
  2. the eigenvalues of the block u u^T + conj(u) conj(u)^T of one off-line pair,
     against the closed form Re s +- sqrt(n^2 - (Im s)^2), n = L J(2 eta), s ~ a L^2;
  3. the resulting thresholds: an upper bound tr Gtilde^b <= K d l1^b (b even) is
     incompatible with an isolated extremal off-line zero once eta > 1/(b lambda) + o(1).
"""
from mpmath import mp, mpf, quad, exp, log, pi, sqrt, mpc
import numpy as np

mp.dps = 30


def ramp(t):
    """C^3 smoothstep: 0 for t <= 0, 1 for t >= 1, nondecreasing."""
    if t <= 0:
        return mpf(0)
    if t >= 1:
        return mpf(1)
    return 35 * t**4 - 84 * t**5 + 70 * t**6 - 20 * t**7


def make_phi(L, w):
    L, w = mpf(L), mpf(w)
    return lambda u: ramp((L / 2 - abs(u)) / w)


def phihat(phi, L, z):
    """int_{-L/2}^{L/2} phi(u) e^{i z u} du (phi even)."""
    L = mpf(L)
    f = lambda u: phi(u) * exp(1j * z * u)
    return quad(f, [-L / 2, -L / 2 + 1, 0, L / 2 - 1, L / 2])


def Phi(phi, L, z):
    L = mpf(L)
    f = lambda u: phi(u) ** 2 * exp(1j * z * u)
    return quad(f, [-L / 2, -L / 2 + 1, 0, L / 2 - 1, L / 2])


def J(phi, L, y):
    """J(y) = int phi^2 e^{y u} du = Phi(-i y)."""
    return Phi(phi, L, mpc(0, -1) * y).real


def check_poisson(L=12, w=1, K=100):
    phi = make_phi(L, w)
    T0 = mpf(0)
    a = mpc("0.37", "-0.30")   # a zero at height 0.37 with eta = 0.30
    b = mpc("-0.81", "0.15")
    tau = [T0 + 2 * pi * k / L for k in range(-K, K + 1)]
    lhs = sum(phihat(phi, L, a - t) * phihat(phi, L, b - t) for t in tau)
    rhs = L * Phi(phi, L, a - b)
    return lhs, rhs


def pair_block(L=12, w=1, eta=mpf("0.3"), K=150):
    """Eigenvalues of u u^T + conj(u) conj(u)^T for one off-line pair."""
    phi = make_phi(L, w)
    z0 = mpc(0, -eta)
    u = np.array([complex(phihat(phi, L, z0 - 2 * pi * k / L)) for k in range(-K, K + 1)])
    B = np.real(np.outer(u, u) + np.outer(u.conj(), u.conj()))
    ev = np.linalg.eigvalsh(B)
    n = float(np.vdot(u, u).real)
    s = complex(np.sum(u * u))
    closed = (s.real + np.sqrt(n**2 - s.imag**2), s.real - np.sqrt(n**2 - s.imag**2))
    aL2 = float(L * quad(lambda t: phi(t) ** 2, [-L / 2, L / 2]))
    return ev[-1], ev[0], closed, n, float(L * J(phi, L, 2 * eta)), s, aL2


def thresholds():
    """Smallest eta at which J(2 eta)/2 exceeds (K d)^{1/b} l1 (b = 2, 4, 6)."""
    rows = []
    K = mpf(13) / 4 + mpf("0.04")
    for logT in (20, 40, 80, 160):
        T = exp(logT)
        for lam in (mpf(1), mpf("0.9")):
            l = log(T / (2 * pi))
            L = lam * l
            l1 = l + 2 * log(2) - 1
            d = L * T / (2 * pi)
            for b in (2, 4, 6):
                bound = (K * d) ** (mpf(1) / b) * l1
                # lower bound J(2eta) >= (e^{2 eta (L/2 - 1)} - 1)/(2 eta), w = 1
                lo, hi = mpf("1e-6"), mpf("0.5")
                f = lambda e: (exp(2 * e * (L / 2 - 1)) - 1) / (2 * e) / 2 - bound
                if f(hi) < 0:
                    rows.append((logT, float(lam), b, None, 1 / (b * float(lam))))
                    continue
                for _ in range(80):
                    mid = (lo + hi) / 2
                    if f(mid) > 0:
                        hi = mid
                    else:
                        lo = mid
                rows.append((logT, float(lam), b, float(hi), 1 / (b * float(lam))))
    return rows


if __name__ == "__main__":
    lhs, rhs = check_poisson()
    print("1. Complex Gabor-Poisson identity (L = 12, w = 1, 201 grid points):")
    print("   sum_k phihat(a - tau_k) phihat(b - tau_k) =", mp.nstr(lhs, 12))
    print("   L Phi(a - b)                            =", mp.nstr(rhs, 12))
    print("   |difference|                             =", mp.nstr(abs(lhs - rhs), 3))

    top, bot, closed, n, nJ, s, aL2 = pair_block()
    print("\n2. Off-line pair block, L = 12, w = 1, eta = 0.3 (301 grid points):")
    print(f"   eigenvalues (numerical)   : {top:.10g}, {bot:.10g}")
    print(f"   Re s +- sqrt(n^2 - Im^2 s): {closed[0]:.10g}, {closed[1]:.10g}")
    print(f"   n = |u|^2 = {n:.10g};  L J(2 eta) = {nJ:.10g};  s = {s.real:.8g}{s.imag:+.2g}i;  a L^2 = {aL2:.8g}")

    print("\n3. Threshold eta* above which one isolated extremal off-line pair alone forces")
    print("   tr Gtilde^b > K d l1^b  (K = 13/4 + 0.04, w = 1); compare 1/(b lambda):")
    print("   log T  lambda  b    eta*       1/(b lambda)")
    for logT, lam, b, eta, ref in thresholds():
        es = "  > 1/2 " if eta is None else f"{eta:9.5f}"
        print(f"   {logT:5d}  {lam:5.2f}  {b:2d}  {es}   {ref:9.5f}")
