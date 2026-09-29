"""Independent recomputation of L(N) = sum_{n<=N} lambda(n) (Liouville), to
cross-check the L(N) column of Table 1 in "Liouville Correlations and
Zero-Free Half-Planes".  Omega(n) comes from a plain Eratosthenes sieve over
prime powers; none of the preprint's code is reused."""
import numpy as np

NMAX = 20_000_000
is_p = np.ones(NMAX + 1, dtype=bool)
is_p[:2] = False
for p in range(2, int(NMAX**0.5) + 1):
    if is_p[p]:
        is_p[p * p :: p] = False
omega = np.zeros(NMAX + 1, dtype=np.int16)
for p in np.nonzero(is_p)[0]:
    q = int(p)
    while q <= NMAX:
        omega[q::q] += 1
        q *= int(p)
lam = np.where(omega % 2 == 0, 1, -1).astype(np.int64)
lam[0] = 0
L = np.cumsum(lam)
print("N          L(N)")
for N in (10**3, 10**4, 10**5, 10**6, 10**7, 2 * 10**7):
    print(f"{N:<10d} {L[N]}")
