"""Section 5: theta and mu identities, the transformation law of the completion
on Gamma_0(4), the shadow, and the cusp estimates.  Floating-point complex
arithmetic (double precision); Jacobi's identity is checked exactly."""
import cmath
import math
import sys

from awlib import Report, ser_one, ser_J

PI = math.pi
I = 1j


def e(z):
    return cmath.exp(2 * PI * I * z)


def qpow(tau, alpha):
    return e(alpha * tau)


def eta(tau, terms=400):
    q = e(tau)
    p = 1
    qj = 1
    for _ in range(1, terms):
        qj *= q
        p *= 1 - qj
        if abs(qj) < 1e-18:
            break
    return qpow(tau, 1 / 24) * p


def Jd(q, d, terms=4000):
    p = 1
    for j in range(1, terms):
        t = q ** (d * j)
        p *= 1 - t
        if abs(t) < 1e-18:
            break
    return p


def L(q):
    s = 0
    for r in range(0, 80):
        s += q ** (4 * r * (r + 1)) / (1 - q ** (4 * r + 1))
    for r in range(1, 80):
        s -= q ** (4 * r * r - 1) / (1 - q ** (4 * r - 1))
    return s


def U(q):
    return sum((-1) ** r * q ** (r * r) / (1 - q ** (2 * r - 1)) for r in range(-60, 61))


def Theta(q):
    return sum((-1) ** r * q ** (r * r) for r in range(-60, 61))


def S_of_tau(tau):
    q = e(tau)
    oddprod = Jd(q, 1) / Jd(q, 2)            # (q;q^2)_inf
    return oddprod ** 2 / Jd(q, 4) ** 4 * L(q)


def A(tau):
    return eta(tau) ** 4 / (eta(2 * tau) ** 3 * eta(4 * tau) ** 4)


def E0(tau):
    return eta(tau) ** 4 / (eta(2 * tau) ** 6 * eta(4 * tau) ** 4)


def vtheta(v, T, K=60):
    s = 0
    for m in range(-K, K):
        nu = m + 0.5
        s += cmath.exp(PI * I * nu * nu * T + 2 * PI * I * nu * (v + 0.5))
    return s


def mu(w, v, T, K=60):
    s = 0
    K = min(K, max(6, int(100 / T.imag)))   # keep exp() in range; dropped terms ~ e^{-pi r^2 Im T}
    for r in range(-K, K + 1):
        s += (-1) ** r * cmath.exp(PI * I * (r * r + r) * T + 2 * PI * I * r * v) / (
            1 - cmath.exp(2 * PI * I * (r * T + w)))
    return cmath.exp(PI * I * w) / vtheta(v, T) * s


def Ecal(t):
    return math.erf(math.sqrt(PI) * t)       # 2 * int_0^t exp(-pi s^2) ds


def R0T(T, K=60):
    """R(0;T) = sum_{nu in Z+1/2} (sgn nu - E(nu sqrt(2 Im T))) (-1)^{nu-1/2} e^{-pi i nu^2 T}."""
    y = T.imag
    s = 0
    for m in range(-K, K):
        nu = m + 0.5
        if PI * nu * nu * y > 600:   # |term| < erfc-tail * e^{pi nu^2 y} < e^{-600}
            continue
        sg = 1 if nu > 0 else -1
        coef = sg * math.erfc(math.sqrt(PI) * abs(nu) * math.sqrt(2 * y))  # = sgn - E(...)
        s += coef * (-1) ** m * cmath.exp(-PI * I * nu * nu * T)
    return s


def muhat(w, T):
    return mu(w, w, T) + 0.5 * I * R0T(T)


def Shat(tau):
    return A(tau) * (-0.5 * I) * muhat(tau, 2 * tau)


def K(tau):
    return eta(2 * tau) ** 3 * muhat(tau, 2 * tau)


def chi4(d):
    return 0 if d % 2 == 0 else (1 if d % 4 == 1 else -1)


def act(g, tau):
    a, b, c, d = g
    return (a * tau + b) / (c * tau + d)


def test_point(g, shift=0.031, s=1.0):
    a, b, c, d = g
    if c == 0:
        return complex(0.17, 0.6)
    return complex(-d / c + shift, s / abs(c))


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def main():
    R = Report("check_completion")

    R.section("Theta and mu at (tau; 2 tau), and S = A f")
    worst = [0] * 5
    for tau in (complex(0.1, 0.4), complex(-0.33, 0.25), complex(0.45, 0.7), complex(0.02, 0.15)):
        q = e(tau)
        worst[0] = max(worst[0], rel(vtheta(tau, 2 * tau), -I * qpow(tau, -0.25) * Theta(q)))
        worst[1] = max(worst[1], rel(Theta(q), Jd(q, 1) ** 2 / Jd(q, 2)))
        m_ = mu(tau, tau, 2 * tau)
        worst[2] = max(worst[2], rel(m_, -I * qpow(tau, -0.25) * U(q) / Theta(q)))
        worst[2] = max(worst[2], rel(m_, 2 * I * qpow(tau, 0.75) * L(q) / Theta(q)))
        worst[3] = max(worst[3], rel(A(tau), qpow(tau, -0.75) * Jd(q, 1) ** 4 / (Jd(q, 2) ** 3 * Jd(q, 4) ** 4)))
        worst[4] = max(worst[4], rel(S_of_tau(tau), A(tau) * (-0.5 * I) * m_))
    R.check(worst[0] < 1e-11, f"theta(tau;2tau) = -i q^(-1/4) Theta(q)   [{worst[0]:.1e}]")
    R.check(worst[1] < 1e-11, f"Theta(q) = J_1^2/J_2   [{worst[1]:.1e}]")
    R.check(worst[2] < 1e-10, f"mu(tau,tau;2tau) = -i q^(-1/4) U/Theta = 2i q^(3/4) L/Theta   [{worst[2]:.1e}]")
    R.check(worst[3] < 1e-11, f"A = q^(-3/4) J_1^4/(J_2^3 J_4^4)   [{worst[3]:.1e}]")
    R.check(worst[4] < 1e-10, f"S(e^(2 pi i tau)) = A(tau) f(tau)   [{worst[4]:.1e}]")
    tau = complex(0.2, 0.3)
    R.check(rel(Shat(tau), S_of_tau(tau) + 0.25 * A(tau) * R0T(2 * tau)) < 1e-11,
            "S^ = A f^ = S + A R(0;2tau)/4  (5.1)")

    R.section("Transformation laws")
    G04 = [(1, 1, 0, 1), (1, 0, 4, 1), (3, 1, 8, 3), (1, -1, 4, -3), (5, 3, 8, 5),
           (-1, 0, 4, -1), (3, -1, 4, -1), (-1, 0, 0, -1), (7, 2, 24, 7), (5, 2, 12, 5), (-5, 2, 12, -5)]
    for g in G04:
        assert g[0] * g[3] - g[1] * g[2] == 1 and g[2] % 4 == 0
    worst = 0
    for g in G04:
        tau = test_point(g)
        lhs = Shat(act(g, tau))
        rhs = chi4(g[3]) * (g[2] * tau + g[3]) ** -1 * Shat(tau)
        worst = max(worst, rel(lhs, rhs))
    R.check(worst < 1e-9, f"(5.2) S^(g tau) = chi_-4(d)(c tau+d)^-1 S^(tau) for {len(G04)} elements of Gamma_0(4)   [{worst:.1e}]")
    worst = 0
    for g in G04:
        tau = test_point(g)
        worst = max(worst, rel(E0(act(g, tau)), chi4(g[3]) * (g[2] * tau + g[3]) ** -3 * E0(tau)))
    R.check(worst < 1e-9, f"E_0 has weight -3 and character chi_-4 on Gamma_0(4)   [{worst:.1e}]")
    G02 = [(1, 0, 2, 1), (1, -1, 2, -1), (3, 1, 2, 1), (3, 2, 4, 3), (1, 0, 6, 1)]
    worst = 0
    for g in G02:
        tau = test_point(g)
        worst = max(worst, rel(K(act(g, tau)), (g[2] * tau + g[3]) ** 2 * K(tau)))
    R.check(worst < 1e-9, f"K = eta(2tau)^3 mu^(tau,tau;2tau) has weight 2, trivial character on Gamma_0(2)   [{worst:.1e}]")
    tau = complex(0.23, 0.37)
    R.check(rel(Shat(tau), -0.5 * I * E0(tau) * K(tau)) < 1e-12 and
            rel(A(tau), E0(tau) * eta(2 * tau) ** 3) < 1e-12, "S^ = -i E_0 K/2 and A = E_0 eta(2tau)^3")
    g = (1, 0, 4, 1)
    tau = test_point(g)
    bad = rel(S_of_tau(act(g, tau)), (4 * tau + 1) ** -1 * S_of_tau(tau))
    R.check(bad > 1e-3, f"control: the uncompleted S does not satisfy (5.2)   [discrepancy {bad:.2f}]")
    g = (1, 0, 2, 1)
    tau = test_point(g)
    bad = min(rel(Shat(act(g, tau)), s * (2 * tau + 1) ** -1 * Shat(tau)) for s in (1, -1, I, -I))
    R.check(bad > 1e-3, f"control: (5.2) fails for (1 0; 2 1), which lies in Gamma_0(2) but not Gamma_0(4)   [{bad:.2f}]")

    R.section("Shadow")
    # Jacobi: sum_{m>=0} (-1)^m (2m+1) q^{m(m+1)} = (q^2;q^2)^3  (exact series)
    M = 400
    lhs = [0] * M
    m = 0
    while m * (m + 1) < M:
        lhs[m * (m + 1)] += (-1) ** m * (2 * m + 1)
        m += 1
    rhs = ser_J(2, M, power=3)
    R.check(lhs == rhs, "Jacobi: sum_{l odd} (-1)^((l-1)/2) l q^(l^2/4) = 2 eta(2tau)^3  (exact mod q^400)")
    worst_r = worst_xi = worst_s = 0
    h = 1e-5
    for tau in (complex(0.1, 0.35), complex(-0.27, 0.8), complex(0.41, 0.22)):
        y = tau.imag
        dR = 0.5 * ((R0T(2 * (tau + h)) - R0T(2 * (tau - h))) / (2 * h)
                    + I * (R0T(2 * (tau + I * h)) - R0T(2 * (tau - I * h))) / (2 * h))
        target = -I / math.sqrt(y) * (eta(2 * tau) ** 3).conjugate()
        worst_r = max(worst_r, rel(dR, target))

        def fhat(t):
            return -0.5 * I * muhat(t, 2 * t)
        df = 0.5 * ((fhat(tau + h) - fhat(tau - h)) / (2 * h)
                    + I * (fhat(tau + I * h) - fhat(tau - I * h)) / (2 * h))
        xi = 2 * I * y ** 0.5 * df.conjugate()
        worst_xi = max(worst_xi, rel(xi, -0.5 * eta(2 * tau) ** 3))
        dS = 0.5 * ((Shat(tau + h) - Shat(tau - h)) / (2 * h)
                    + I * (Shat(tau + I * h) - Shat(tau - I * h)) / (2 * h))
        xim1 = 2 * I * y ** -1 * dS.conjugate()
        worst_s = max(worst_s, rel(xim1, -A(tau).conjugate() / (2 * y ** 1.5) * eta(2 * tau) ** 3))
    R.check(worst_r < 1e-6, f"d/d(tau-bar) R(0;2tau) = -(i/sqrt y) conj(eta(2tau)^3)  (5.4)   [{worst_r:.1e}]")
    R.check(worst_xi < 1e-6, f"xi_(1/2) f^ = -eta(2tau)^3/2  (5.3)   [{worst_xi:.1e}]")
    R.check(worst_s < 1e-6, f"xi_(-1) S^ = -conj(A) eta(2tau)^3/(2 y^(3/2)), not holomorphic   [{worst_s:.1e}]")

    R.section("Cusp estimates")
    worst = 0
    for z in (complex(0.13, 0.9), complex(-0.4, 0.5), complex(0.3, 1.3)):
        lhs = muhat(-1 / (2 * z), -1 / z)
        rhs = cmath.exp(3 * PI * I / 4) * cmath.sqrt(z) * (mu(0.5, 0.5, z) + 0.5 * I * R0T(z))
        worst = max(worst, rel(lhs, rhs))
    R.check(worst < 1e-9, f"mu^(-1/(2z),-1/(2z);-1/z) = e^(3 pi i/4) z^(1/2) mu^(1/2,1/2;z)   [{worst:.1e}]")
    ratios = []
    for v in (2.0, 4.0, 6.0):
        z = complex(0.1, v)
        p = e(z)
        ratios.append((abs(vtheta(0.5, z) / (-2 * qpow(z, 1 / 8)) - 1),
                       abs(mu(0.5, 0.5, z) / (-0.25 * I * qpow(z, -1 / 8)) - 1), abs(p)))
    R.check(all(a < 5 * pp + 1e-14 and bb < 5 * pp + 1e-14 for a, bb, pp in ratios),
            "theta(1/2;z) = -2p^(1/8)(1+O(p)) and mu(1/2,1/2;z) = -i p^(-1/8)(1+O(p))/4  (Im z = 2,4,6)")
    okR = True
    for v in (1.0, 2.0, 3.0, 5.0):
        okR &= abs(R0T(complex(0.2, v))) <= 2 * v ** -0.5 * math.exp(-PI * v / 4)
        okR &= abs(R0T(2 * complex(0.2, v))) <= 2 * v ** -0.5 * math.exp(-PI * v / 2)
    R.check(okR, "R(0;z) = O(v^(-1/2) e^(-pi v/4)) and R(0;2tau) = O(y^(-1/2) e^(-pi y/2))")
    g = (1, 0, 2, 1)
    R.check(abs(act(g, 1e12j) - 0.5) < 1e-9 and g[2] % 2 == 0,
            "the cusp 1/2 of Gamma_0(4) is the image of infinity under (1 0; 2 1) in Gamma_0(2)")
    R.check(sum(dd * r for dd, r in ((1, 4), (2, -6), (4, -4))) == -24 and
            sum((4 // dd) * r for dd, r in ((1, 4), (2, -6), (4, -4))) == 0 and (4 - 6 - 4) // 2 == -3,
            "eta-quotient conditions for E_0: sum d r_d = -24, sum (4/d) r_d = 0, weight -3")
    return R.done()


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
