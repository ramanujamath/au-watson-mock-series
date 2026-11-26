"""Shared helpers for the Au Watson series checks (standard library only).

Contents
  Report            pass/fail bookkeeping
  poch, poch_list   exact q-shifted factorials (a;r)_n, including negative n
  B, a_term, P_coeffs, f0, F, T, C, H, b, E, rho   the objects of the paper
  Series helpers    truncated integer power series in q
  Poly              exact polynomials in q, x, y (dict of exponent triples)
  UPoly             univariate polynomials with Fraction coefficients
"""
from fractions import Fraction as Fr


class Report:
    def __init__(self, name):
        self.name, self.npass, self.nfail = name, 0, 0

    def check(self, cond, msg):
        if cond:
            self.npass += 1
            print(f"  [ok]   {msg}", flush=True)
        else:
            self.nfail += 1
            print(f"  [FAIL] {msg}", flush=True)
        return bool(cond)

    def info(self, msg):
        print(f"  [info] {msg}", flush=True)

    def section(self, msg):
        print(f"\n== {msg}", flush=True)

    def done(self):
        print(f"\n-- {self.name}: {self.npass} passed, {self.nfail} failed")
        return self.nfail == 0


# ---------------------------------------------------------------- products
def poch(a, r, n):
    """(a; r)_n = prod_{j<n} (1 - a r^j); for n<0, 1/(a r^n; r)_{-n}."""
    if n >= 0:
        p = a * 0 + 1  # keep the number type of a (Fraction, float, complex)
        for j in range(n):
            p *= 1 - a * r ** j
        return p
    return 1 / poch(a * r ** n, r, -n)


def poch_list(params, r, n):
    p = 1
    for a in params:
        p *= poch(a, r, n)
    return p


# coefficients of P(q, x) from (1.2): list of (c, u, v) meaning c q^u x^v
P_TERMS = [
    (1, 0, 0), (1, 1, 1), (1, 0, 2), (1, 2, 2), (2, 3, 3), (3, 2, 4), (2, 4, 4),
    (1, 3, 5), (-2, 2, 6), (1, 4, 6), (1, 6, 6), (-1, 1, 7), (-1, 3, 7),
    (2, 5, 7), (-1, 4, 8), (-2, 3, 9), (-3, 5, 9), (-2, 4, 10), (-1, 5, 11),
    (-1, 7, 11), (-1, 6, 12), (-1, 7, 13),
]


def P(q, x):
    return sum(c * q ** u * x ** v for c, u, v in P_TERMS)


def B(q, n):
    Q = q * q
    return (poch(q, Q, n) ** 4 * poch(Q, Q * Q, n) ** 2 * poch(-Q, Q, 4 * n)
            / poch(Q * Q, Q * Q, 2 * n) ** 3)


def a_term(q, n):
    x = q ** (2 * n)
    return (q ** (2 * n) * P(q, x) * B(q, n)
            / ((1 + x ** 4) * (1 + q * x) * (1 + q * q * x * x) ** 3))


def f0(q, n):
    x = q ** (2 * n)
    return (q * B(q, n) * (1 - x ** 4) ** 3 * (x * x - q * x + q * q)
            / (2 * (1 + x ** 4) * (x - q) ** 2 * (q * q + x ** 4)))


def F(q, n, k):
    Q, x = q * q, q ** (2 * n)
    num = poch_list([q * x, Q / x ** 2, -Q / x ** 2, -Q / x ** 2], Q, k)
    den = poch_list([q ** 3 / x, q ** 3 / x, -q ** 3 / x, -q ** 4 / x ** 4], Q, k)
    return (f0(q, n) * (q / x) ** k * (1 + q ** 3 * Q ** (2 * k) / x ** 3)
            / (1 + q ** 3 / x ** 3) * num / den)


def T(q, n):
    return sum(F(q, n, k) for k in range(2 * n))


def C(q, N):
    x = q ** (2 * N)
    return (q * q * x * B(q, N) * (1 - x ** 4) ** 3
            / (2 * (1 + x ** 4) * (x - q) ** 2 * (q * q + x ** 4)))


def H(q, N, k):
    Q, x = q * q, q ** (2 * N)
    return (C(q, N) * Q ** k * poch(Q / x ** 2, Q, k) ** 2 * poch(-Q / x ** 2, Q, k)
            / (poch(-Q * Q / x ** 4, Q, k) * poch(q ** 3 / x, Q, k) ** 2))


def b(q, r):
    return q ** (2 * r * (r + 1)) / poch(q, q * q, r + 1) ** 2


def E(q, N):
    return C(q, N) * q ** (2 - 2 * N) / poch(q, q * q, N - 1) ** 2


def rho(q, N, r):
    Q = q * q
    m = N + r
    return (poch(q ** (2 * N - 2 * r), Q, m) ** 2 * poch(-q ** (2 * N - 2 * r), Q, m)
            / poch(-q ** (6 * N - 2 * r - 2), Q, m))


# ---------------------------------------------------------- power series
def ser_one(M):
    s = [0] * M
    s[0] = 1
    return s


def ser_mul_binom(s, c, k):
    """s <- s * (1 + c q^k), k >= 1 (in place)."""
    M = len(s)
    for i in range(M - 1, k - 1, -1):
        s[i] += c * s[i - k]
    return s


def ser_div_binom(s, c, k):
    """s <- s / (1 + c q^k), k >= 1 (in place)."""
    M = len(s)
    for i in range(k, M):
        s[i] -= c * s[i - k]
    return s


def ser_mul(s, t):
    M = len(s)
    out = [0] * M
    for i, si in enumerate(s):
        if si:
            for j in range(M - i):
                if t[j]:
                    out[i + j] += si * t[j]
    return out


def ser_poch_mul(s, sign, e0, d, n, power=1, divide=False):
    """multiply (or divide) s by prod_{j<n} (1 - sign q^{e0 + j d})^power."""
    M = len(s)
    for j in range(n):
        e = e0 + j * d
        if e >= M:
            break
        for _ in range(power):
            if divide:
                ser_div_binom(s, -sign, e)
            else:
                ser_mul_binom(s, -sign, e)
    return s


def ser_J(d, M, power=1, divide=False):
    """(q^d; q^d)_infinity ^ power (or its inverse), truncated."""
    s = ser_one(M)
    return ser_poch_mul(s, 1, d, d, M // d + 1, power, divide)


# ------------------------------------------------- polynomials in q, x, y
class Poly:
    __slots__ = ("t",)

    def __init__(self, t=None):
        self.t = {k: v for k, v in (t or {}).items() if v}

    @staticmethod
    def const(c):
        return Poly({(0, 0, 0): c})

    @staticmethod
    def mono(c, a, b_, e):
        return Poly({(a, b_, e): c})

    def __add__(self, o):
        o = o if isinstance(o, Poly) else Poly.const(o)
        t = dict(self.t)
        for k, v in o.t.items():
            t[k] = t.get(k, 0) + v
        return Poly(t)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.t.items()})

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Poly) else Poly.const(o)))

    def __rsub__(self, o):
        return Poly.const(o) - self

    def __mul__(self, o):
        if not isinstance(o, Poly):
            return Poly({k: v * o for k, v in self.t.items()})
        t = {}
        for (a1, b1, c1), v1 in self.t.items():
            for (a2, b2, c2), v2 in o.t.items():
                k = (a1 + a2, b1 + b2, c1 + c2)
                t[k] = t.get(k, 0) + v1 * v2
        return Poly(t)

    __rmul__ = __mul__

    def __pow__(self, n):
        r = Poly.const(1)
        for _ in range(n):
            r = r * self
        return r

    def __eq__(self, o):
        return self.t == (o if isinstance(o, Poly) else Poly.const(o)).t

    def scale(self, qx=0, qy=0):
        """substitute x -> q^qx x and y -> q^qy y."""
        return Poly({(a + qx * b_ + qy * c, b_, c): v
                     for (a, b_, c), v in self.t.items()})

    def set_y(self, yval_q_exp):
        """substitute y -> q^e (a monomial in q)."""
        t = {}
        for (a, b_, c), v in self.t.items():
            k = (a + yval_q_exp * c, b_, 0)
            t[k] = t.get(k, 0) + v
        return Poly(t)

    def ev(self, q, x, y):
        return sum(v * q ** a * x ** b_ * y ** c for (a, b_, c), v in self.t.items())


qP, xP, yP = Poly.mono(1, 1, 0, 0), Poly.mono(1, 0, 1, 0), Poly.mono(1, 0, 0, 1)


# --------------------------------------- univariate polynomials (Fractions)
class UPoly:
    """coefficients c[0] + c[1] X + ... with Fraction entries."""

    def __init__(self, c):
        c = list(c)
        while len(c) > 1 and c[-1] == 0:
            c.pop()
        self.c = c

    def __mul__(self, o):
        out = [Fr(0)] * (len(self.c) + len(o.c) - 1)
        for i, a in enumerate(self.c):
            if a:
                for j, bb in enumerate(o.c):
                    out[i + j] += a * bb
        return UPoly(out)

    def ev(self, X):
        r = Fr(0)
        for a in reversed(self.c):
            r = r * X + a
        return r

    def div_root(self, X0):
        """synthetic division by (X - X0); requires ev(X0) == 0."""
        n = len(self.c) - 1
        out = [Fr(0)] * n
        acc = Fr(0)
        for i in range(n, 0, -1):
            acc = acc * X0 + self.c[i]
            out[i - 1] = acc
        return UPoly(out)


def lin(a, b_):
    """the univariate polynomial a + b X."""
    return UPoly([Fr(a), Fr(b_)])


def ratio_limit(num, den, X0):
    """value at X0 of the rational function num/den after cancelling
    common factors (X - X0)."""
    m = 0
    while num.ev(X0) == 0 and den.ev(X0) == 0:
        num, den = num.div_root(X0), den.div_root(X0)
        m += 1
    return num.ev(X0) / den.ev(X0), m
