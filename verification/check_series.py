"""Section 1 and the end of Section 3: the summand, Theorem 1.1, the bilateral
rearrangements and Corollary 3.2, checked as exact power-series identities
(and numerically at complex q)."""
import cmath
import sys
from fractions import Fraction as Fr

from awlib import (Report, P_TERMS, B, a_term, poch, ser_one, ser_mul_binom,
                   ser_div_binom, ser_mul, ser_poch_mul, ser_J)

M = 220  # precision: identities are checked modulo q^M


def S_series(M):
    """sum_n a_n as a power series, each a_n built from binomial factors."""
    tot = [0] * M
    for n in range(0, M // 2 + 1):
        s = ser_one(M)
        # B_n = (q;Q)_n^4 (Q;Q^2)_n^2 (-Q;Q)_{4n} / (Q^2;Q^2)_{2n}^3
        ser_poch_mul(s, 1, 1, 2, n, power=4)
        ser_poch_mul(s, 1, 2, 4, n, power=2)
        ser_poch_mul(s, -1, 2, 2, 4 * n)
        ser_poch_mul(s, 1, 4, 4, 2 * n, power=3, divide=True)
        # q^{2n} P(q, q^{2n})
        p = [0] * M
        for c, u, v in P_TERMS:
            e = 2 * n + u + 2 * n * v
            if e < M:
                p[e] += c
        s = ser_mul(s, p)
        # denominators (1+q^{8n})(1+q^{2n+1})(1+q^{4n+2})^3
        if n == 0:
            assert all(c % 2 == 0 for c in s)
            s = [c // 2 for c in s]
        else:
            ser_div_binom(s, 1, 8 * n)
        ser_div_binom(s, 1, 2 * n + 1)
        for _ in range(3):
            ser_div_binom(s, 1, 4 * n + 2)
        tot = [x + y for x, y in zip(tot, s)]
    return tot


def geom(M, sign, coef_exp, k):
    """sign * q^coef_exp / (1 - q^k) as a series (k >= 1)."""
    s = [0] * M
    if coef_exp < M:
        s[coef_exp] = sign
    return ser_div_binom(s, -1, k)


def L_series(M):
    s = [0] * M
    r = 0
    while 4 * r * (r + 1) < M:
        s = [x + y for x, y in zip(s, geom(M, 1, 4 * r * (r + 1), 4 * r + 1))]
        r += 1
    r = 1
    while 4 * r * r - 1 < M:
        s = [x + y for x, y in zip(s, geom(M, -1, 4 * r * r - 1, 4 * r - 1))]
        r += 1
    return s


def L_series_direct_negative(M):
    """the r < 0 terms of sum q^{4r(r+1)}/(1-q^{4r+1}) rewritten term by term:
    for r=-s, q^{4s(s-1)}/(1-q^{1-4s}) = -q^{4s^2-1}/(1-q^{4s-1})."""
    s = [0] * M
    for ss in range(1, 40):
        e = 4 * ss * (ss - 1)
        # 1/(1-q^{-m}) = -q^m/(1-q^m)
        m = 4 * ss - 1
        if e + m < M:
            s = [x + y for x, y in zip(s, geom(M, -1, e + m, m))]
    return s


def U_series(M):
    s = [0] * M
    for r in range(1, 40):
        if r * r < M:
            s = [x + y for x, y in zip(s, geom(M, (-1) ** r, r * r, 2 * r - 1))]
    for t in range(0, 40):  # r = -t: (-1)^t q^{t^2}/(1-q^{-2t-1})
        e = t * t + 2 * t + 1
        if e < M:
            s = [x + y for x, y in zip(s, geom(M, -((-1) ** t), e, 2 * t + 1))]
    return s


def main():
    R = Report("check_series")

    R.section("The polynomial P(q,x) of (1.2)")
    R.check(sum(c for c, _, _ in P_TERMS) == 0, "sum of coefficients c_j = 0, so P(1,1) = 0")
    R.check(sum(abs(c) * u for c, u, _ in P_TERMS) == 112, "sum |c_j| u_j = 112")
    R.check(sum(abs(c) * 2 * v for c, _, v in P_TERMS) == 416,
            "sum |c_j| v_j = 416, where P_n(q) = sum c_j q^{u_j + v_j n} (so v_j = 2 deg_x)")
    su = sum(c * u for c, u, _ in P_TERMS)
    sv = sum(c * v for c, _, v in P_TERMS)
    R.check(su == -20 and sv == -84,
            f"P_n'(1) = sum c_j(u_j + 2n v_j) = {su} + {2*sv}n, i.e. -P_n'(1) = 4(42n+5)")

    R.section("The two product forms of B_n and the value a_0")
    ok = True
    for q in (Fr(1, 3), Fr(2, 5), Fr(-3, 7), Fr(5, 9)):
        for n in range(0, 8):
            alt = (poch(q, q * q, n) ** 4 * poch(q * q, q * q, 2 * n) * poch(q * q, q ** 4, n)
                   * poch(q ** 4, q ** 4, 4 * n)
                   / (poch(q * q, q * q, 4 * n) * poch(q ** 4, q ** 4, n)
                      * poch(q ** 4, q ** 4, 2 * n) ** 3))
            ok &= alt == B(q, n)
    R.check(ok, "B_n equals the even/odd separated form (4 rational q, 0 <= n < 8, exact)")
    ok = all(a_term(q, 0) == (1 - q ** 5) / ((1 + q) * (1 + q * q) ** 2)
             for q in (Fr(1, 3), Fr(2, 5), Fr(-3, 7), Fr(5, 9), Fr(7, 2)))
    R.check(ok, "a_0 = (1-q^5)/((1+q)(1+q^2)^2)  (1.4)")

    R.section(f"Theorem 1.1 as a power series identity mod q^{M}")
    S = S_series(M)
    R.info(f"S(q) = {S[:12]} ...")
    R.check(S[0] == 1, "S(0) = 1")
    rhs = ser_poch_mul(ser_one(M), 1, 1, 2, M, power=2)          # (q;q^2)_inf^2
    rhs = ser_mul(rhs, ser_J(4, M, power=4, divide=True))         # / J_4^4
    Lser = L_series(M)
    R.check(Lser[0] == 1, "L(0) = 1 (the apparent singularity at q = 0 is removable)")
    rhs = ser_mul(rhs, Lser)
    R.check(S == rhs, f"S(q) = (q;q^2)_inf^2 / J_4^4 * L(q)  holds mod q^{M}")

    R.section("Bilateral rearrangements used in the proof of Theorem 1.1")
    pos = [0] * M
    r = 0
    while 4 * r * (r + 1) < M:
        pos = [x + y for x, y in zip(pos, geom(M, 1, 4 * r * (r + 1), 4 * r + 1))]
        r += 1
    neg = L_series_direct_negative(M)
    R.check([x + y for x, y in zip(pos, neg)] == Lser,
            "L(q) = sum_{r>=0} q^{4r(r+1)}/(1-q^{4r+1}) - sum_{r>=1} q^{4r^2-1}/(1-q^{4r-1})")
    U = U_series(M)
    R.check(U == [-2 * c for c in ([0] + Lser[:-1])], "U(q) = -2 q L(q)")
    # sum_r b_r = -U(q)/(q J_2): both sides as series
    sb = [0] * M
    for r in range(0, 40):
        e = 2 * r * (r + 1)
        if e >= M:
            break
        t = [0] * M
        t[e] = 1
        ser_poch_mul(t, 1, 1, 2, r + 1, power=2, divide=True)
        sb = [x + y for x, y in zip(sb, t)]
    for j in range(0, 200):  # b_{-j-1} = q^{2j} (q;Q)_j^2
        if 2 * j >= M:
            break
        t = [0] * M
        t[2 * j] = 1
        ser_poch_mul(t, 1, 1, 2, j, power=2)
        sb = [x + y for x, y in zip(sb, t)]
    rhsb = ser_mul([-c for c in U[1:]] + [0], ser_J(2, M, divide=True))
    R.check(sb[:M - 1] == rhsb[:M - 1], "sum_{r in Z} b_r = -U(q)/(q J_2)   (3.4 after reindexing)")
    # (3.4) literally
    s1 = [0] * M
    for j in range(1, 200):
        if 2 * j >= M:
            break
        t = [0] * M
        t[2 * j] = 1
        ser_poch_mul(t, 1, 1, 2, j, power=2)
        s1 = [x + y for x, y in zip(s1, t)]
    s2 = [0] * M
    for j in range(0, 40):
        e = 2 * j * (j - 1)
        if e >= M:
            break
        t = [0] * M
        t[e] = 1
        ser_poch_mul(t, 1, 1, 2, j, power=2, divide=True)
        s2 = [x + y for x, y in zip(s2, t)]
    R.check([x + y for x, y in zip(s1, s2)][:M - 1] == rhsb[:M - 1],
            "(3.4): sum_{j>=1} q^{2j}(q;Q)_j^2 + sum_{j>=0} q^{2j(j-1)}/(q;Q)_j^2 = -U/(qJ_2)")

    R.section("Corollary 3.2: S = J_2^6 E_21 / (2 J_1^2 J_4^4)")
    E21 = [0] * M
    for n in range(0, M // 2 + 1):
        t = [0] * M
        for c, e in ((3, 0), (3, 2 * n + 1), (-1, 4 * n), (2, 4 * n + 2),
                     (-3, 6 * n + 1), (-3, 8 * n + 2), (-1, 10 * n + 3)):
            if 2 * n + e < M:
                t[2 * n + e] += c
        ser_poch_mul(t, 1, 1, 2, n, power=6)
        for _ in range(3):
            ser_div_binom(t, 1, 2 * n + 1)
        ser_poch_mul(t, 1, 2, 2, 2 * n, power=3, divide=True)
        E21 = [x + y for x, y in zip(E21, t)]
    lhs = ser_mul(S, ser_J(1, M, power=2))
    lhs = ser_mul(lhs, ser_J(4, M, power=4))
    lhs = [2 * c for c in lhs]
    rhs = ser_mul(E21, ser_J(2, M, power=6))
    R.check(lhs == rhs, f"2 J_1^2 J_4^4 S = J_2^6 E_21  mod q^{M}")

    R.section("Theorem 1.1 at complex q (direct summation vs closed form)")

    def S_direct(q):
        return sum(a_term(q, n) for n in range(400))

    def S_closed(q):
        p, j4 = 1, 1
        for j in range(4000):
            p *= 1 - q ** (2 * j + 1)
        for j in range(1, 4000):
            j4 *= 1 - q ** (4 * j)
        Lv = sum(q ** (4 * r * (r + 1)) / (1 - q ** (4 * r + 1)) for r in range(0, 60))
        Lv -= sum(q ** (4 * r * r - 1) / (1 - q ** (4 * r - 1)) for r in range(1, 60))
        return p * p / j4 ** 4 * Lv

    worst = 0
    for q in (0.3 + 0.2j, -0.5 + 0.1j, 0.55 * cmath.exp(2j), 0.6j, 0.45 * cmath.exp(-2.5j)):
        d, c = S_direct(q), S_closed(q)
        worst = max(worst, abs(d - c) / abs(c))
    R.check(worst < 1e-10, f"5 complex points with |q| <= 0.6: max relative difference {worst:.1e}")
    return R.done()


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
