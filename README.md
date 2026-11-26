# Au's Watson series: evaluation and mock modular completion

This repository holds the manuscript, audit and verification code for

> **Evaluation of Au's Watson series and its mixed mock modular completion**
> Sumit Kumar Jha (Greater Noida, India · sumitkumarjha75@gmail.com)

Ramanujan's series Σ (1/2)ₙ³(42n+5)/(n!³64ⁿ) = 16/π has many q-analogues. K. C. Au built one from a Watson-type Wilf–Zeilberger seed. He asked whether it is mock modular. The paper answers yes and evaluates the series in closed form:

S(q) = (q;q²)∞² / (q⁴;q⁴)∞⁴ · Σ_{r∈ℤ} q^{4r(r+1)} / (1 − q^{4r+1}).

The bilateral Lambert series is a value of Zwegers' μ-function, so S is an eta quotient times a mock theta function. The paper gives:

- an explicit completion of weight −1 with character χ₋₄ on Γ₀(4);
- the shadow −η(2τ)³/2 of its weight-1/2 factor;
- a direct proof that 8S(q)/(1−q) → 16/π;
- a comparison with Campbell's q-analogue of the same formula.

The evaluation has three steps:
1. **Finite identity.** A Wilf–Zeilberger pair with an explicit certificate turns each partial sum into a finite sum.
2. **Watson's transformation** rewrites that sum with nonnegative terms.
3. **Centred limit.** Re-indexing around the middle term and dominated convergence give a two-sided series. A classical identity then turns it into the Lambert series.

## Layout

```
manuscript/        main.tex, sections/, references.tex, main.pdf
  original/        the draft that was audited
verification/
  awlib.py           exact q-products, the paper's objects, power series, polynomials in q,x,y
  check_series.py    Theorem 1.1 and Corollary 4.2 (Example XXI) as power series mod q^220
  check_wz.py        Appendix certificate (full expansion), Lemma 3.1 incl. removable points,
                     Watson's transformation, Lemma 3.2
  check_limits.py    centred limit (Lemma 4.1) and the classical limit 16/pi (Prop. 1.2)
  check_completion.py  theta/mu identities, Gamma_0(4) transformation, shadow, cusps
  check_campbell.py  Campbell's summand vs Au's Example III
  run_all.py         runs everything; logs in verification/logs/
AUDIT.md           audit report
```

## Running the checks

```
cd verification && python3 run_all.py
```

This needs Python 3.9 or later and the standard library only. Script messages cite the draft's equation numbers; `AUDIT.md` maps them to the rewrite. All 70 checks pass in about 20 seconds.
