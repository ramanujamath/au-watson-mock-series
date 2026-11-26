"""Section 6: Campbell's summand versus Au's Example III summand (exact)."""
import sys
from fractions import Fraction as Fr

from awlib import Report, poch


def p_n(q, n):
    return (1 + 3 * q ** (2 * n + 1) + 3 * q ** (4 * n + 2) - 2 * q ** (6 * n + 1) + q ** (6 * n + 3)
            - 3 * q ** (8 * n + 2) - 3 * q ** (10 * n + 3))


def u_n(q, n):
    return (q ** (6 * n * n) * p_n(q, n) * poch(q, q * q, n) ** 6
            / ((1 + q ** (2 * n + 1)) ** 3 * poch(q * q, q * q, 2 * n) ** 3))


def c_n(q, n):
    r = -(1 - q ** (2 * n - 1)) * p_n(q, n)
    first = ((-1) ** n * q ** ((n - 1) * (3 * n + 5)) * poch(q, q * q, n) ** 4 * poch(q ** 3, q * q, n)
             / (poch(-q * q, q * q, n) ** 2 * poch(q * q, q * q, n) ** 3
                * poch(q ** (3 - 4 * n), q * q, n) * poch(q ** (2 * n + 2), q * q, n)))
    second = (poch(q ** -1, q ** 4, n) * poch(q, q ** 4, n)
              / (poch(q * q, q ** 4, n) * poch(q ** 6, q ** 4, n)))
    return first * second * r / (1 + q ** (2 * n + 1)) ** 2


def main():
    R = Report("check_campbell")
    QS = (Fr(1, 3), Fr(2, 5), Fr(-3, 7), Fr(5, 9), Fr(7, 3), Fr(-11, 4))

    R.section("Proposition 6.1: c_n = (1-q^2) u_n / q^6")
    ok = all(c_n(q, n) == (1 - q * q) * u_n(q, n) / q ** 6 for q in QS for n in range(0, 12))
    R.check(ok, "c_n = (1-q^2) u_n/q^6 for 0 <= n < 12 and six rational q (exact)")

    R.section("The product identities used in the proof")
    ok = True
    for q in QS:
        Q = q * q
        for n in range(0, 10):
            ok &= poch(q ** (3 - 4 * n), Q, n) == (-1) ** n * q ** (-3 * n * n + 2 * n) * poch(q ** (2 * n - 1), Q, n)
            ok &= poch(q ** -1, q ** 4, n) * poch(q, q ** 4, n) == poch(q ** -1, Q, 2 * n)
            ok &= ((1 - q ** (2 * n - 1)) * poch(q ** -1, Q, 2 * n) / poch(q ** (2 * n - 1), Q, n)
                   == (1 - q ** -1) * poch(q, Q, n))
            ok &= poch(-Q, Q, n) * poch(Q, Q, n) == poch(q ** 4, q ** 4, n)
            ok &= poch(Q, q ** 4, n) * poch(q ** 4, q ** 4, n) == poch(Q, Q, 2 * n)
            ok &= poch(q ** 3, Q, n) / poch(q, Q, n) == (1 - q ** (2 * n + 1)) / (1 - q)
            ok &= poch(q ** 6, q ** 4, n) / poch(Q, q ** 4, n) == (1 - q ** (4 * n + 2)) / (1 - Q)
    R.check(ok, "the seven displayed product identities hold for 0 <= n < 10, six q (exact)")

    R.section("Classical limit of Au's Example III summand")
    import math
    ok = True
    for n in range(0, 12):
        terms = [(1, 0), (3, 2 * n + 1), (3, 4 * n + 2), (-2, 6 * n + 1), (1, 6 * n + 3),
                 (-3, 8 * n + 2), (-3, 10 * n + 3)]
        ok &= sum(c for c, _ in terms) == 0                       # p_n(1) = 0
        ok &= -sum(c * ex for c, ex in terms) == 42 * n + 5       # -p_n'(1) = 42n+5
        lim = Fr(1)                                               # (q;q^2)_n^6/(q^2;q^2)_{2n}^3 at q=1
        for j in range(n):
            lim *= Fr(2 * j + 1) ** 6
        for j in range(2 * n):
            lim /= Fr(2 * j + 2) ** 3
        half = Fr(1)
        for j in range(n):
            half *= Fr(1, 2) + j
        target = half ** 3 / (math.factorial(n) ** 3 * Fr(64) ** n)
        ok &= Fr(8, 8) * (42 * n + 5) * lim == target * (42 * n + 5)   # 8 u_n/(1-q) -> 8 (42n+5) lim / (1+1)^3
    R.check(ok, "8 u_n/(1-q) -> (1/2)_n^3 (42n+5)/(n!^3 64^n) as q -> 1 (exact termwise limit, n < 12)")
    return R.done()


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
