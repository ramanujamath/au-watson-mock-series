"""Section 2 and Appendix A: the rational certificate, the telescoping lemma
(including the two removable singularities) and the Watson reduction.
All checks are exact (integers and Fractions)."""
import random
import sys
import time
from fractions import Fraction as Fr

from awlib import (Report, Poly, qP, xP, yP, UPoly, lin, ratio_limit, poch, poch_list,
                   P, P_TERMS, B, a_term, f0, F, T, C, H)

q_, x_, y_ = qP, xP, yP
one = Poly.const(1)

# ---------------------------------------------------------------- Appendix A
D = ((x_ - 1) ** 2 * (x_ + 1) ** 2 * (x_ ** 2 + 1) * (x_ ** 2 - y_) * (x_ ** 2 + y_) ** 2
     * (q_ * x_ + 1) * (q_ ** 2 * x_ ** 2 + 1) * (q_ ** 2 * x_ ** 2 - y_)
     * (q_ ** 2 * x_ ** 2 + y_) ** 2)
H_num = (x_ - q_ * y_) ** 2 * (x_ + q_ * y_) * (q_ ** 2 * y_ + x_ ** 4)
H_den = q_ * y_ * D * (q_ ** 3 * y_ ** 2 + x_ ** 3)


def m(c, a, b_):
    return Poly.mono(c, a, b_, 0)


v = [
    m(1, 9, 15),
    m(1, 9, 13) + m(1, 8, 14) - m(1, 8, 10) + m(1, 7, 13) + m(1, 7, 11) - m(1, 7, 9)
    + m(1, 6, 12) - m(1, 6, 8) + m(1, 5, 11) - m(1, 5, 9) - m(1, 4, 10),
    m(1, 8, 10) - m(1, 8, 8) - m(1, 7, 13) + m(1, 7, 11) + m(1, 6, 10) - m(1, 6, 6)
    + m(1, 5, 9) - m(1, 5, 5) + m(1, 4, 8) - m(2, 4, 6) - m(1, 2, 6),
    -m(1, 7, 9) - m(2, 5, 9) + m(1, 5, 7) - m(1, 4, 10) + m(1, 4, 6) - m(1, 3, 9)
    + m(1, 3, 5) + m(1, 2, 4) - m(1, 2, 2) - m(1, 1, 7) + m(1, 1, 5),
    -m(1, 5, 5) - m(1, 4, 6) + m(1, 4, 4) - m(1, 3, 7) + m(1, 3, 3) - m(1, 2, 6)
    + m(1, 2, 4) + m(1, 2, 2) - m(1, 1, 5) + m(1, 1, 1) + m(1, 0, 2),
    one,
]
V = sum((v[j] * y_ ** j for j in range(6)), Poly())
R_num = x_ * (y_ + 1) * V * H_num
R_den = H_den
NN = ((x_ ** 4 + y_) * (q_ * x_ - 1) * (q_ ** 2 * x_ ** 4 + y_) * (q_ ** 3 * x_ ** 3 + y_ ** 2)
      * (q_ ** 4 * x_ ** 4 + y_) * (q_ * x_ * y_ - 1))
rn_num, rn_den = NN * H_num, q_ * H_den
# r^(k) after multiplying numerator and denominator by x^10
rk_num = (q_ * (x_ ** 3 + q_ ** 7 * y_ ** 2) * (1 - q_ * x_ * y_) * (x_ ** 2 - q_ ** 2 * y_)
          * (x_ ** 2 + q_ ** 2 * y_) ** 2)
rk_den = ((x_ ** 3 + q_ ** 3 * y_ ** 2) * (x_ - q_ ** 3 * y_) ** 2 * (x_ + q_ ** 3 * y_)
          * (x_ ** 4 + q_ ** 4 * y_))


def rk_literal(q, x, y):
    """r^(k) exactly as displayed in Appendix A."""
    return (q / x * (1 + q ** 7 * y ** 2 / x ** 3) / (1 + q ** 3 * y ** 2 / x ** 3)
            * (1 - q * x * y) * (1 - q * q * y / x ** 2) * (1 + q * q * y / x ** 2) ** 2
            / ((1 - q ** 3 * y / x) ** 2 * (1 + q ** 3 * y / x) * (1 + q ** 4 * y / x ** 4)))


def ev(num, den, q, x, y):
    return num.ev(q, x, y) / den.ev(q, x, y)


def as_upoly(poly, q, y):
    """a Poly in q,x,y with q and y specialised -> UPoly in x."""
    deg = max(b_ for (_, b_, _) in poly.t)
    c = [Fr(0)] * (deg + 1)
    for (a, b_, e), val in poly.t.items():
        c[b_] += val * Fr(q) ** a * Fr(y) ** e
    return UPoly(c)


def phi_k(q, k):
    """F(n,k)/B_n as a rational function of X = x = q^{2n} (num, den UPolys)."""
    Q = q * q
    num = UPoly([q]) * (lin(1, 0) * UPoly([1, 0, 0, 0, -1])) * UPoly([1, 0, 0, 0, -1]) \
        * UPoly([1, 0, 0, 0, -1]) * UPoly([q * q, -q, 1])
    den = UPoly([2]) * UPoly([1, 0, 0, 0, 1]) * lin(-q, 1) * lin(-q, 1) * UPoly([q * q, 0, 0, 0, 1])
    num = num * UPoly([q ** k])                       # (q/X)^k: the X^-k is absorbed below
    num = num * UPoly([q ** 3 * Q ** (2 * k), 0, 0, 1])   # X^3 + q^3 Q^{2k}
    den = den * UPoly([q ** 3, 0, 0, 1])                  # X^3 + q^3
    for j in range(k):
        num = num * lin(1, -q * Q ** j)                   # 1 - q Q^j X
        num = num * UPoly([-Q ** (j + 1), 0, 1])          # X^2 - Q^{j+1}
        num = num * UPoly([Q ** (j + 1), 0, 1]) * UPoly([Q ** (j + 1), 0, 1])
        den = den * lin(-q ** 3 * Q ** j, 1) * lin(-q ** 3 * Q ** j, 1)
        den = den * lin(q ** 3 * Q ** j, 1)
        den = den * UPoly([q ** 4 * Q ** j, 0, 0, 0, 1])
    # net power of X: (q/X)^k gives X^-k, numerator products X^{-6k},
    # denominator products X^{-7k}; total X^0.
    return num, den


def watson_8phi7(a, sa, b_, c, d, e, Q, m_):
    z = a * a * Q ** (m_ + 2) / (b_ * c * d * e)
    tot = 0
    for k in range(m_ + 1):
        nu = poch_list([a, Q * sa, -Q * sa, b_, c, d, e, Q ** (-m_)], Q, k)
        de = poch_list([Q, sa, -sa, a * Q / b_, a * Q / c, a * Q / d, a * Q / e,
                        a * Q ** (m_ + 1)], Q, k)
        tot += nu / de * z ** k
    return tot


def watson_rhs(a, b_, c, d, e, Q, m_):
    pre = poch_list([a * Q, a * Q / (d * e)], Q, m_) / poch_list([a * Q / d, a * Q / e], Q, m_)
    tot = 0
    for k in range(m_ + 1):
        nu = poch_list([a * Q / (b_ * c), d, e, Q ** (-m_)], Q, k)
        de = poch_list([Q, a * Q / b_, a * Q / c, d * e * Q ** (-m_) / a], Q, k)
        tot += nu / de * Q ** k
    return pre * tot


def main():
    Rp = Report("check_wz")
    QS = (Fr(1, 3), Fr(2, 5), Fr(3, 4), Fr(-1, 2))

    Rp.section("Appendix A: polynomial identities (exact expansion in Z[q,x,y])")
    t0 = time.time()
    # lower boundary identity
    lhs = R_num.scale(0, 0).set_y(0) * q_ * (1 - x_ ** 4) ** 3 * (x_ ** 2 - q_ * x_ + q_ ** 2)
    lhs_den = R_den.set_y(0) * 2 * (1 + x_ ** 4) * (x_ - q_) ** 2 * (q_ ** 2 + x_ ** 4)
    Pp = sum((Poly.mono(c, u, w, 0) for c, u, w in P_TERMS), Poly())
    rhs = x_ * Pp
    rhs_den = (1 + x_ ** 4) * (1 + q_ * x_) * (1 + q_ ** 2 * x_ ** 2) ** 3
    Rp.check(lhs * rhs_den == rhs * lhs_den,
             "R(q,x,1) f0-factor = x P(q,x)/((1+x^4)(1+qx)(1+q^2x^2)^3)  (so G(n,0) = a_n)")
    # compatibility identity r^(n)(x,y) r^(k)(q^2 x, y) = r^(k)(x,y) r^(n)(x, q^2 y)
    L1 = rn_num * rk_num.scale(2, 0)
    L1d = rn_den * rk_den.scale(2, 0)
    R1 = rk_num * rn_num.scale(0, 2)
    R1d = rk_den * rn_den.scale(0, 2)
    Rp.check(L1 * R1d == R1 * L1d, "compatibility identity r^(n)(x,y) r^(k)(q^2x,y) = r^(k)(x,y) r^(n)(x,q^2y)")
    # (A.1): r^(k)(y) R(q^2 y) - R(y) = r^(n)(y) - 1
    a, bb = rk_num, rk_den
    c, d = R_num.scale(0, 2), R_den.scale(0, 2)
    e, f = R_num, R_den
    g = rn_num   # rn_den = q f
    left = q_ * a * c * f - q_ * e * bb * d
    right = (g - q_ * f) * bb * d
    Rp.check(left == right, f"certificate identity (A.1)  [{len(left.t)} monomials after clearing denominators]")
    e2 = (x_ * (y_ + 1) * (V + q_ * x_ * y_ ** 3) * H_num)        # perturbed certificate
    c2 = e2.scale(0, 2)
    Rp.check(q_ * a * c2 * f - q_ * e2 * bb * d != right,
             "negative control: changing one monomial of V breaks (A.1)")
    Rp.info(f"polynomial expansion took {time.time() - t0:.1f}s")
    # degrees of the k-quotient
    Rp.check(max(k3 for (_, _, k3) in rk_num.t) == 6 and max(k3 for (_, _, k3) in rk_den.t) == 6,
             "r^(k) has numerator and denominator degree 6 in y = Q^k")
    rnd = random.Random(1)
    ok = True
    for _ in range(200):
        qq, xx, yy = (Fr(rnd.randint(1, 30), rnd.randint(31, 60)) for _ in range(3))
        ok &= ev(rk_num, rk_den, qq, xx, yy) == rk_literal(qq, xx, yy)
    Rp.check(ok, "the cleared form of r^(k) equals the displayed form (200 random rational points)")
    # coprimality of the y-roots: numerator and denominator have no common factor
    ok = True
    for _ in range(20):
        qq, xx = Fr(rnd.randint(2, 9), rnd.randint(11, 19)), Fr(rnd.randint(2, 9), rnd.randint(11, 19))
        nroots = [x_r for x_r in (1 / (qq * xx), xx ** 2 / qq ** 2, -xx ** 2 / qq ** 2)]
        dden = [xx / qq ** 3, -xx / qq ** 3, -xx ** 4 / qq ** 4]
        ok &= all(rk_den.ev(qq, xx, r_) != 0 for r_ in nroots)
        ok &= all(rk_num.ev(qq, xx, r_) != 0 for r_ in dden)
    Rp.check(ok, "the k-quotient is reduced (no common rational root in y at random q, x)")

    Rp.section("The ratios r^(k), r^(n) match F(n,k) (exact, x = q^{2n}, y = q^{2k})")
    okk = okn = True
    for q in QS:
        for n in range(1, 5):
            x = q ** (2 * n)
            for k in range(0, 2 * n):
                y = q ** (2 * k)
                if k + 1 < 2 * n:
                    okk &= ev(rk_num, rk_den, q, x, y) == F(q, n, k + 1) / F(q, n, k)
                okn &= ev(rn_num, rn_den, q, x, y) == F(q, n + 1, k) / F(q, n, k)
            okn &= ev(rn_num, rn_den, q, x, 1) == f0(q, n + 1) / f0(q, n)
    Rp.check(okk, "r^(k) = F(n,k+1)/F(n,k) for 1<=n<=4, 0<=k<2n-1, four values of q")
    Rp.check(okn, "r^(n) = F(n+1,k)/F(n,k) for 1<=n<=4, 0<=k<2n, four values of q (incl. y=1)")

    Rp.section("Lemma 2.1 including the removable points k = 2n, 2n+1")
    ok_tel = ok_g0 = ok_rem = ok_zero = True
    for q in QS:
        for n in range(1, 5):
            X0 = q ** (2 * n)
            G = []
            for k in range(0, 2 * n + 4):
                yv = q ** (2 * k)
                pn, pd = phi_k(q, k)
                num = as_upoly(R_num, q, yv) * pn
                den = as_upoly(R_den, q, yv) * pd
                val, mult = ratio_limit(num, den, X0)
                G.append(val * B(q, n))
                if k in (2 * n, 2 * n + 1):
                    ok_rem &= mult == 1 and val != 0
                if k >= 2 * n + 2:
                    ok_zero &= val == 0
            ok_g0 &= G[0] == a_term(q, n)
            for k in range(0, 2 * n + 3):
                ok_tel &= (F(q, n + 1, k) - F(q, n, k)) == G[k + 1] - G[k]
    Rp.check(ok_g0, "G(n,0) = a_n  (1<=n<=4, four q)")
    Rp.check(ok_tel, "F(n+1,k) - F(n,k) = G(n,k+1) - G(n,k) for all 0 <= k <= 2n+2, limits taken at k=2n,2n+1")
    Rp.check(ok_rem, "at k = 2n and 2n+1 exactly one factor (x - q^{2n}) cancels and G is finite and nonzero")
    Rp.check(ok_zero, "G(n,k) = 0 for k >= 2n+2")
    ok = True
    for q in QS:
        ok &= T(q, 1) == -(1 - q ** 5) / ((1 + q) * (1 + q * q) ** 2)
        for N in range(1, 8):
            ok &= sum(a_term(q, n) for n in range(N)) == -T(q, N)
    Rp.check(ok, "T_1 = -a_0 and sum_{n<N} a_n = -T_N  (2.3) for 1<=N<=7, four q, exact")

    Rp.section("Watson's transformation (DLMF 17.9.15) and Lemma 2.2")
    ok, used = True, 0
    for _ in range(30):
        Q = Fr(rnd.randint(1, 5), rnd.randint(6, 11))
        sa = Fr(rnd.randint(1, 9), rnd.randint(2, 9)) * rnd.choice([1, -1])
        a = sa * sa
        b_, c, d, e = (Fr(rnd.randint(1, 20), rnd.randint(1, 20)) * rnd.choice([1, -1]) for _ in range(4))
        m_ = rnd.randint(0, 5)
        try:
            lhs, rhs = watson_8phi7(a, sa, b_, c, d, e, Q, m_), watson_rhs(a, b_, c, d, e, Q, m_)
        except ZeroDivisionError:
            continue
        ok &= lhs == rhs
        used += 1
    Rp.check(ok and used >= 25, f"the terminating 8phi7 -> 4phi3 transformation holds at {used} random rational points")
    okw = okp = okt = True
    for q in (Fr(-1, 4), Fr(-4, 9), Fr(-9, 25)):   # -q a square, so sqrt(a) is rational
        sq = {Fr(-1, 4): Fr(1, 2), Fr(-4, 9): Fr(2, 3), Fr(-9, 25): Fr(3, 5)}[q]
        Q = q * q
        for N in range(1, 5):
            x = q ** (2 * N)
            a = -q ** 3 / x ** 3
            sa = sq ** 3 / abs(q) ** (3 * N)
            assert sa * sa == a
            b_, c, d, e, m_ = q * x, -Q / x ** 2, Q, -Q / x ** 2, 2 * N - 1
            okw &= watson_8phi7(a, sa, b_, c, d, e, Q, m_) == T(q, N) / f0(q, N)
            pre = poch_list([a * Q, a * Q / (d * e)], Q, m_) / poch_list([a * Q / d, a * Q / e], Q, m_)
            okp &= pre == -q * x / (x * x - q * x + q * q)
            for k in range(2 * N):
                term = (poch_list([a * Q / (b_ * c), d, e, Q ** (-m_)], Q, k)
                        / poch_list([Q, a * Q / b_, a * Q / c, d * e * Q ** (-m_) / a], Q, k) * Q ** k)
                okt &= -f0(q, N) * pre * term == H(q, N, k)
    Rp.check(okw, "with a=-q^3/x^3, b=qx, c=e=-Q/x^2, d=Q, m=2N-1 the 8phi7 equals T_N/f0(N)")
    Rp.check(okp, "the product prefactor equals -qx/(x^2-qx+q^2)")
    Rp.check(okt, "-f0(N) * prefactor * (k-th term of the 4phi3) = H_N(k) term by term")
    ok = okpos = True
    for q in QS + (Fr(9, 10), Fr(1, 10)):
        for N in range(1, 8):
            hs = [H(q, N, k) for k in range(2 * N)]
            ok &= sum(a_term(q, n) for n in range(N)) == sum(hs)
            if 0 < q < 1:
                okpos &= all(h >= 0 for h in hs)
    Rp.check(ok, "(2.7): sum_{n<N} a_n = sum_{k<2N} H_N(k) for 1<=N<=7, six q, exact")
    Rp.check(okpos, "every H_N(k) is nonnegative for the five tested q in (0,1)")
    return Rp.done()


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
