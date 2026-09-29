"""Imaginary parts of the first N nontrivial zeros of zeta, cached to disk.

mpmath.zetazero locates each zero on the critical line.  All zeros with
0 < Im(rho) < 3*10^12 are known to lie on the line and be simple
(Platt and Trudgian, Bull. LMS 53 (2021)), so for the heights used here the
list is complete and exact up to the working precision.
"""
import json
import os
import sys

import mpmath as mp

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "zeta_zeros.json")


def zeros(n, dps=30):
    have = []
    if os.path.exists(CACHE):
        with open(CACHE) as f:
            have = json.load(f)
    if len(have) < n:
        mp.mp.dps = dps
        for k in range(len(have) + 1, n + 1):
            have.append(mp.nstr(mp.zetazero(k).imag, dps))
        with open(CACHE, "w") as f:
            json.dump(have, f)
    return [mp.mpf(x) for x in have[:n]]


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    z = zeros(n)
    print(f"{len(z)} zeros, last gamma = {mp.nstr(z[-1], 15)}")
