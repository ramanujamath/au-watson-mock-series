"""Section 3 (Lemma 3.1, the centred limit) and Section 4 (the classical
limit 8S(q)/(1-q) -> 16/pi)."""
import math
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as Fr

from awlib import Report, P, P_TERMS, B, a_term, poch, C, H, b, E, rho


def S_closed_float(q, terms=6000):
    p, j2, j4 = 1.0, 1.0, 1.0
    for j in range(terms):
        p *= 1 - q ** (2 * j + 1)
    for j in range(1, terms):
        j2 *= 1 - q ** (2 * j)
        j4 *= 1 - q ** (4 * j)
    return p, j2, j4


def S_ratio(q, nmax=400):
    """sum a_n with B_n from the ratio B_{n+1}/B_n (safe near q = 1)."""
    tot, Bn = 0.0, 1.0
    for n in range(nmax):
        t = q ** (2 * n)
        x = t
        tot += (q ** (2 * n) * P(q, x) * Bn
                / ((1 + x ** 4) * (1 + q * x) * (1 + q * q * x * x) ** 3))
        Bn *= ((1 - q * t) / (1 - q * q * t)) ** 3 * (
            (1 + q ** 2 * t ** 4) * (1 + q ** 4 * t ** 4) * (1 + q ** 6 * t ** 4) * (1 + q ** 8 * t ** 4)
            / ((1 + q * t) * (1 + q * q * t) ** 3 * (1 + q * q * t * t) ** 3 * (1 + q ** 4 * t * t) ** 3))
        if Bn < 1e-300:
            break
    return tot


def main():
    R = Report("check_limits")
    QS = (Fr(1, 3), Fr(2, 5), Fr(3, 4), Fr(-1, 2), Fr(7, 3))

    R.section("Lemma 3.1: the centred form of H_N (exact)")
    ok = True
    for q in QS:
        for N in range(1, 7):
            for r in range(-N, N):
                ok &= H(q, N, N + r) == E(q, N) * b(q, r) * rho(q, N, r)
    R.check(ok, "H_N(N+r) = E_N b_r rho_{N,r} for 1<=N<=6, -N<=r<=N-1, five q")
    ok = True
    for q in QS:
        Q = q * q
        for N in range(1, 8):
            ok &= b(q, -N) == q ** (2 * N - 2) * poch(q, Q, N - 1) ** 2
            ok &= E(q, N) * b(q, -N) == C(q, N) and rho(q, N, -N) == 1
        for j in range(0, 8):
            ok &= b(q, -j - 1) == q ** (2 * j) * poch(q, Q, j) ** 2
    R.check(ok, "b_{-N} = q^{2N-2}(q;Q)_{N-1}^2, rho_{N,-N} = 1, b_{-j-1} = q^{2j}(q;Q)_j^2 (negative-index convention)")
    ok = True
    for q in QS:
        Q = q * q
        for N in range(1, 6):
            x = q ** (2 * N)
            for r in range(-N, N - 1):
                v = q ** (2 * r)
                left = (Q * (1 - Q * v / x) ** 2 * (1 + Q * v / x)
                        / ((1 + Q * Q * v / x ** 3) * (1 - q ** 3 * v) ** 2))
                right = (Q * Q * v * v / (1 - q ** 3 * v) ** 2 * (1 - x / (Q * v)) ** 2
                         * (1 + x / (Q * v)) / (1 + x ** 3 / (Q * Q * v)))
                ok &= left == H(q, N, N + r + 1) / H(q, N, N + r) == right
    R.check(ok, "the two displayed successive ratios agree with each other and with H_N")

    R.section("Uniform bounds and limits used for dominated convergence (0<q<1)")
    ok = True
    for qf in (0.2, 0.5, 0.8, 0.95):
        Q = qf * qf
        cap = 1.0
        for j in range(2000):
            cap *= 1 + qf ** 2 * Q ** j
        for N in range(1, 30):
            for r in range(-N, N):
                val = rho(qf, N, r)
                ok &= -1e-15 <= val <= cap * (1 + 1e-12)
    R.check(ok, "0 <= rho_{N,r} <= (-q^2;Q)_inf for N < 30 and q in {0.2,0.5,0.8,0.95}")
    ok = True
    for qf in (0.2, 0.5, 0.8):
        p, j2, j4 = S_closed_float(qf)
        Binf = p ** 4 * j2 / j4 ** 4
        ok &= abs(B(qf, 80) - Binf) < 1e-12 * Binf
        ok &= abs(E(qf, 80) - Binf / (2 * p * p)) < 1e-10 * Binf / (2 * p * p)
        for r in range(0, 20):
            ok &= b(qf, r) <= qf ** (2 * r * (r + 1)) / p ** 2 * (1 + 1e-12)
    R.check(ok, "B_N -> (q;Q)^4 J_2/J_4^4, E_N -> B_inf/(2(q;Q)^2), b_r <= q^{2r(r+1)}/(q;Q)^2 (floats)")
    worst = 0
    for qf in (0.1, 0.3, 0.6, 0.85, 0.93):
        p, j2, j4 = S_closed_float(qf)
        sb = (sum(b(qf, r) for r in range(0, 60))   # r < 0 via the exactly checked b_{-j-1} formula
              + sum(qf ** (2 * j) * poch(qf, qf * qf, j) ** 2 for j in range(0, 800)))
        direct = sum(a_term(qf, n) for n in range(3000))
        val = p * p * j2 / (2 * j4 ** 4) * sb
        worst = max(worst, abs(val - direct) / abs(direct))
    R.check(worst < 1e-9, f"(3.3): S(q) = (q;Q)^2 J_2/(2J_4^4) sum b_r at five real q, rel. diff {worst:.1e}")

    R.section("Proposition 4.1: termwise limits")
    ok = True
    for n in range(0, 12):
        # (1 - q^m)/(1 - q) -> m; orders of vanishing: 4n + 2n - 6n = 0
        num = Fr(1)
        for j in range(n):
            num *= Fr(2 * j + 1) ** 4 * Fr(4 * j + 2) ** 2
        num *= Fr(2) ** (4 * n)
        den = Fr(1)
        for j in range(2 * n):
            den *= Fr(4 * j + 4) ** 3
        lim_B = num / den
        half = Fr(1)
        for j in range(n):
            half *= Fr(1, 2) + j
        target = half ** 3 / (math.factorial(n) ** 3 * Fr(64) ** n)
        ok &= lim_B == target
        # 8 a_n/(1-q) -> 8 * 4(42n+5) * B / (2*2*8)
        dP = -sum(c * (u + 2 * n * w) for c, u, w in P_TERMS)
        ok &= Fr(8 * dP, 32) * lim_B == target * (42 * n + 5)
    R.check(ok, "lim B_n = (1/2)_n^3/(n!^3 64^n) and lim 8a_n/(1-q) = (1/2)_n^3 (42n+5)/(n!^3 64^n), n < 12")
    ok = True
    for qf in (Fr(1, 2), Fr(9, 10), Fr(99, 100)):
        for n in range(0, 10):
            t = qf ** (2 * n)
            q = qf
            ratio = ((1 - q * t) / (1 - q * q * t)) ** 3 * (
                (1 + q ** 2 * t ** 4) * (1 + q ** 4 * t ** 4) * (1 + q ** 6 * t ** 4) * (1 + q ** 8 * t ** 4)
                / ((1 + q * t) * (1 + q * q * t) ** 3 * (1 + q * q * t * t) ** 3 * (1 + q ** 4 * t * t) ** 3))
            ok &= ratio == B(q, n + 1) / B(q, n)
            ok &= ratio <= (1 + q ** (4 * n + 4)) ** -6
    R.check(ok, "B_{n+1}/B_n formula and the bound <= (1+q^{4n+4})^{-6} (exact, three q, n < 10)")
    ok = True
    for k in range(1, 400):
        qf = k / 400
        for n in range(0, 40):
            ok &= abs(P(qf, qf ** (2 * n))) / (1 - qf) <= 112 + 416 * n + 1e-9
    R.check(ok, "|P_n(q)|/(1-q) <= 112 + 416 n on a grid of 399 q in (0,1), n < 40")

    R.section("Proposition 4.1: the classical value and the limit itself")
    getcontext().prec = 60
    s, t = Decimal(0), Decimal(1)
    for n in range(0, 200):
        s += t * (42 * n + 5)
        t = t * (Decimal(n) + Decimal("0.5")) ** 3 / (Decimal(n + 1) ** 3 * 64)
    pi = Decimal("3.14159265358979323846264338327950288419716939937510582097494")
    R.check(abs(s - 16 / pi) < Decimal(10) ** -50, "(1.1): sum (1/2)_n^3 (42n+5)/(n!^3 64^n) = 16/pi to 50 digits")
    errs = []
    for h in (1e-1, 1e-2, 1e-3, 1e-4):
        qf = math.exp(-h)
        val = 8 * S_ratio(qf) / (1 - qf)
        errs.append(abs(val - 16 / math.pi))
        R.info(f"h = {h:g}: 8S(q)/(1-q) = {val:.12f}   (16/pi = {16 / math.pi:.12f})")
    R.check(errs[-1] < 1e-3 and all(errs[i + 1] < errs[i] for i in range(3)),
            f"8S(e^-h)/(1-e^-h) -> 16/pi as h -> 0 (errors {', '.join(f'{e:.1e}' for e in errs)})")
    return R.done()


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
