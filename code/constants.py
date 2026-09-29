"""Exact-arithmetic check of the constants in the note.

Reproduces Knausgard's constant q (arXiv:2609.33043, Theorem 1) from
(tau, c, m) = (12/5, 17/5, 1298) and verifies the admissibility conditions
and the value q' for the retuned choice (tau, c, m) = (2411/1000, 3411/1000, 1310).
It also shows that m = 1310 is the largest admissible block length and
prints the limiting value of the same bound as m -> infinity.
"""
from fractions import Fraction as F

H = F(3362285207, 5000000000)   # C(f) <= 2 - H   (imported, [trmdy], Prop. 5)
DELTA = F(891, 200000)          # local constant delta (imported, Prop. 6)
PRESSURE = F(1, 2736)           # spread weight in Prop. 6 (imported)
WINDOW = 7                      # seven-point windows


def data(m):
    Delta_m = (m - (WINDOW - 1)) * DELTA
    a = Delta_m / m
    beta = (WINDOW - 1) * (m - (WINDOW - 1)) * PRESSURE / m
    return Delta_m, a, beta


def admissible(tau, c, m):
    Delta_m, a, beta = data(m)
    rh = 6 * c - 7 - c * c
    rk = 4 * c - 2 - c * c
    checks = {
        "tau >= 0": tau >= 0,
        "c >= max(1+tau, 2+tau/2)": c >= max(1 + tau, 2 + tau / 2),
        "Delta_m <= tau^2": Delta_m <= tau * tau,
        "r_h - a > 0": rh - a > 0,
        "r_k - 2a > 0": rk - 2 * a > 0,
    }
    return checks, (rh - a, rk - 2 * a)


def q(m):
    _, a, beta = data(m)
    return (1 + H - beta) / (2 - a)


def best_m():
    """Largest m for which some (tau, c) is admissible.

    For tau >= 2 the binding constraints are c >= 1 + tau, tau^2 >= Delta_m
    and r_k = 4c - 2 - c^2 > 2a. Since r_k decreases for c > 2, the optimal
    choice is c = 1 + tau with tau^2 = Delta_m, i.e. r_k = 1 + 2 tau - tau^2.
    We test this with exact rationals by squaring (tau = sqrt(Delta_m)).
    """
    m = 7
    last = None
    while True:
        Delta_m, a, _ = data(m)
        x = Delta_m  # tau^2
        # condition 1 + 2 tau - x - 2a > 0  <=>  2 tau > x + 2a - 1
        rhs = x + 2 * a - 1
        ok = rhs < 0 or 4 * x > rhs * rhs
        if ok:
            last = m
        elif last is not None and m > last + 50:
            return last
        m += 1


if __name__ == "__main__":
    q_old = q(1298)
    assert q_old == F(16260119298029, 19426831050000), q_old
    chk, _ = admissible(F(12, 5), F(17, 5), 1298)
    assert all(chk.values())
    print("Knausgard (tau,c,m)=(12/5,17/5,1298): q =", q_old, "=", f"{float(q_old):.14f}")

    tau, c, m = F(2411, 1000), F(3411, 1000), 1310
    chk, slack = admissible(tau, c, m)
    for k, v in chk.items():
        print(f"  {k:28s} {v}")
    assert all(chk.values())
    Delta_m, a, beta = data(m)
    q_new = q(m)
    print("Retuned (tau,c,m)=(2411/1000,3411/1000,1310):")
    print("  Delta_m =", Delta_m, " a =", a, " beta =", beta)
    print("  slacks r_h - a =", slack[0], ", r_k - 2a =", slack[1])
    print("  q' =", q_new, "=", f"{float(q_new):.14f}")
    print("  q' - q =", f"{float(q_new - q_old):.4e}")
    print("Largest admissible m:", best_m())
    lim = (1 + H - (WINDOW - 1) * PRESSURE) / (2 - DELTA)
    print("Limit m -> infinity with no (tau,c) constraint:", f"{float(lim):.14f}")
