# Audit: *Evaluation of Au's Watson series and its mixed mock modular completion*

Audited draft: `manuscript/original/au_watson.pdf` (11 pages, "Research draft").
Rewrite: `manuscript/main.pdf` (15 pages).
Verification: `verification/` (70 checks, all passing; about 20 s in total).

## Verdict

**No mathematical errors were found.** Every displayed identity in the draft was
checked, most of them in exact rational arithmetic:

- Theorem 1.1 and Corollary 3.2 as power series;
- the Wilf–Zeilberger certificate by full polynomial expansion;
- Lemmas 2.1, 2.2 and 3.1 at many rational values of q;
- the transformation law and shadow in Theorem 5.1 numerically.

The changes in the rewrite are expository: notation clashes removed, terms
defined, and terse steps written out. They are listed below.

## Numbering

The table below and the script messages use the **draft's** numbering. The rewrite renumbers as follows:

| Draft | Rewrite |
|---|---|
| §2, §3, §4, §5, §6 | §3, §4, §5, §6, §7 |
| Lemma 2.1, Lemma 2.2, Lemma 3.1 | Lemma 3.1, Lemma 3.2, Lemma 4.1 |
| Corollary 3.2 | Corollary 4.2 |
| Proposition 4.1 (classical limit) | Proposition 1.2 (proved in §5) |
| Theorem 5.1 (completion) | Theorem 1.3 (proved in §6) |
| Proposition 6.1 | Proposition 7.1 |
| Appendix A | Appendix A |

## What was checked

| Claim (draft numbering) | Script | Method | Result |
|---|---|---|---|
| (1.2) coefficient sums: Σc=0, Σ\|c\|u=112, Σ\|c\|v=416, −Pₙ′(1)=4(42n+5) | `check_series.py` | exact | ✓ |
| Two product forms of Bₙ; a₀ formula (1.4) | `check_series.py` | exact, 4–5 rational q | ✓ |
| **Theorem 1.1** S = (q;q²)²∞ L(q)/J₄⁴ | `check_series.py` | exact power series mod q²²⁰ | ✓ |
| Theorem 1.1 at complex q | `check_series.py` | direct sum vs closed form, 5 points | ✓ (4e-15) |
| L(q) rewritten with r ≥ 1; U = −2qL; Σb_r = −U/(qJ₂); identity (3.4) | `check_series.py` | exact power series | ✓ |
| **Corollary 3.2** S = J₂⁶ E₂₁/(2J₁²J₄⁴) | `check_series.py` | exact power series mod q²²⁰ | ✓ |
| **Certificate identity (A.1)** | `check_wz.py` | full expansion in ℤ[q,x,y] (7278 monomials) | ✓ |
| Negative control: perturbing one monomial of V breaks (A.1) | `check_wz.py` | exact | ✓ |
| Lower boundary identity R(q,x,1)·(…) = xP/(…); compatibility identity | `check_wz.py` | polynomial expansion | ✓ |
| r⁽ᵏ⁾ and r⁽ⁿ⁾ equal the ratios of F(n,k) | `check_wz.py` | exact, 4 q, n ≤ 4 | ✓ |
| r⁽ᵏ⁾ has degree 6/6 in y and is in lowest terms | `check_wz.py` | exact | ✓ |
| **Lemma 2.1** incl. the removable points k = 2n, 2n+1 | `check_wz.py` | univariate limits in x, exact | ✓ |
| Watson's transformation DLMF 17.9.15 | `check_wz.py` | 26 random rational parameter sets | ✓ |
| Specialization: ₈φ₇ = T_N/f₀(N), prefactor = −qx/(x²−qx+q²), termwise = H_N(k) | `check_wz.py` | exact (q = −1/4, −4/9, −9/25 so √a is rational) | ✓ |
| **Lemma 2.2** (2.7) and nonnegativity of H_N(k) | `check_wz.py` | exact, N ≤ 7, six q | ✓ |
| **Lemma 3.1** H_N(N+r) = E_N b_r ρ_{N,r}; b₋ₙ; both displayed ratios | `check_limits.py` | exact, N ≤ 6, five q | ✓ |
| Bounds 0 ≤ ρ ≤ (−q²;Q)∞, B∞, lim E_N, b_r bound, (3.3) | `check_limits.py` | floating point | ✓ |
| **Prop. 4.1**: termwise limits; ratio B_{n+1}/B_n and its bound; \|Pₙ\|/(1−q) ≤ 112+416n | `check_limits.py` | exact / grid | ✓ |
| (1.1) = 16/π; 8S(e^{−h})/(1−e^{−h}) → 16/π | `check_limits.py` | 50 digits; h down to 10⁻⁴ | ✓ |
| ϑ(τ;2τ), μ(τ,τ;2τ), A as eta quotient, S = Af, (5.1) | `check_completion.py` | complex floating point | ✓ (≤1e-15) |
| **(5.2)** transformation on Γ₀(4) | `check_completion.py` | 11 matrices | ✓ (4e-14) |
| E₀ weight −3, character χ₋₄ on Γ₀(4); K weight 2 on Γ₀(2); Ŝ = −iE₀K/2 | `check_completion.py` | numerical | ✓ |
| Controls: uncompleted S fails (5.2); (5.2) fails for (1 0; 2 1) ∉ Γ₀(4) | `check_completion.py` | numerical | ✓ |
| Jacobi's identity; (5.4); **(5.3)** ξ_{1/2} f̂ = −η(2τ)³/2; ξ₋₁Ŝ formula | `check_completion.py` | exact series / finite differences | ✓ (≤1e-8) |
| Inversion of μ̂ at the cusp 0; asymptotics of ϑ(½;z), μ(½,½;z); R bounds | `check_completion.py` | numerical | ✓ |
| **Prop. 6.1** cₙ = (1−q²)uₙ/q⁶ and the seven auxiliary product identities | `check_campbell.py` | exact, n < 12, six q | ✓ |
| Au's Example III has the classical limit (42n+5)·(½)ₙ³/(n!³64ⁿ) | `check_campbell.py` | exact termwise limit | ✓ |

## Expository changes in the rewrite

1. **Notation clashes removed.** The draft reused several letters for different objects:
   - **R:** the WZ certificate and Zwegers' correction R(0;T). The certificate is now 𝒲.
   - **H:** H_N(k) and H(q,x,y). These are now h_N(k) and Φ.
   - **E:** E_N, the error-type function E(t), the eta quotient E₀ and the series E₂₁. These are now κ_N, E(t), g(τ) and S_XXI.
   - **ρ:** ρ_{N,r} and the constant (1+e⁻⁴)⁻⁶ of §4. The constant is now λ.
   - **N:** the summation bound, ⌊1/h⌋ in §4, and 𝒩 in Appendix A. ⌊1/h⌋ is now m.
   - **T:** the finite sums Tₙ and Zwegers' modular variable. The modular variable is now σ.
   - **Section 6:** pₙ (vs Pₙ) and c_j/cₙ. These are now wₙ and α_j.
2. **Undefined or jargon terms explained.** This covers "Wilf–Zeilberger seed/pair", "certificate", "terminating" and "basic hypergeometric series". It also covers "normalized component", which is now "the weight-1/2 factor f̂". Harmonic Maass forms, mock modular forms, shadows, cusps and "mixed mock modular" are defined in a background section with numbered facts F1–F6.
3. **Lemma 2.2.** The draft's chain for the Watson prefactor skipped the intermediate step and did not list the ₄φ₃ parameters. Both are now written out, and the term-by-term match with h_N(k) is stated (and checked).
4. **Lemma 3.1.** "Multiplication by the powers of x, Q, v … shows the ratios are equal" is replaced by the common closed form of both ratios.
5. **Proof of Theorem 1.1.** The index shifts that turn (3.4) into Σb_r = −U/(qJ₂), and the even/odd split giving U = −2qL, are written out.
6. **Proposition 4.1.** The bound B_{n+1}/B_n ≤ λ is needed for n ≤ m−1 only; this is now explicit.
7. **Theorem 5.1.** The derivative computation behind (5.4) is written out. In the cusp-0 step, the draft's "simultaneous integer shifts" is spelled out: inversion with ν_η = e^{−πi/4}, then a shift by 1. The relation of the cusp ½ to ∞ under Γ₀(2) now names the matrix.
8. **Section 6.** The conclusion is stated plainly: Campbell's series is (1−q²)/q⁶ times Au's modular Example III, so it is not the mock modular Watson series. The degree argument now says why parameter specialization cannot raise the degree from 2 to 6.
9. **Front matter.** The rewrite adds author information, MSC codes, keywords, declarations and a link to the verification code.

## References

- **Confirmed online:**
  - Au's arXiv:2403.04555 (v3, Feb 2026);
  - Campbell, J. Symbolic Comput. 134 (2026) 102531 (also cited in Au's reference list);
  - Campbell arXiv:2602.14314;
  - Boylan–Swati, Res. Number Theory 12 (2026), doi 10.1007/s40993-026-00715-4.
- **Not confirmed:**
  - The Adv. Math. 481 (2026) 111181 publication data for Au's paper. Please check it before submission.
  - The exact statements of Au (3.13) and Lovejoy–Osburn (3.11). I did not compare these with the sources. Both consequences the paper uses are verified independently: identity (3.4) and Corollary 3.2 both hold as power series mod q²²⁰.

## Limits of the verification

- Section 5 is checked in double precision at finitely many points and matrices. Zwegers' Theorem 1.11 is used as a black box, though its consequences used here were all confirmed numerically.
- Cusp growth is checked through the stated asymptotics, not by a full analysis at every cusp.
- Lemmas 2.1–3.1 are rational-function identities. They are checked exactly at many rational q, while the certificate itself is proved by full polynomial expansion.

## Observations

- 8S(e^{−h})/(1−e^{−h}) − 16/π ≈ 2.5h for small h (2.5·10⁻⁴ at h = 10⁻⁴). This is recorded as a remark in §5.
- In Lemma 2.1, exactly one factor (x − q^{2n}) cancels at k = 2n and k = 2n+1, and G is finite and nonzero there.
- Floating-point evaluation of Tₙ suffers catastrophic cancellation already at q = 0.4, n = 5, because the terms Q/x² are huge. Exact arithmetic is needed for Lemma 2.1, and the scripts use it.
