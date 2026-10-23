# C42 general-τ certificate: format and exact criterion (FROZEN)

This file is frozen. An independent verifier is written from this file, NOTES.md §4–5, and the papers. Changes after
the freeze commit need a new commit, a note here, and a message to the lead (5697e7f3).

## 1. File format

UTF-8 JSON object. Every number is a JSON **string** holding an exact rational, either a terminating decimal ("0.6883")
or "p/q" with integers p, q > 0. Floats, exponents and non-string numbers are not allowed. Fields:

| key | value |
|---|---|
| `"constant"` | `"C_42"` |
| `"C"` | the claimed bound C |
| `"alpha"` | `[Re α, Im α]` |
| `"breaks"` | `[t_0, t_1, …, t_M]`, M ≥ 1 |
| `"eta"` | `[[Re η_0, Im η_0], …, [Re η_{M−1}, Im η_{M−1}]]` |

Any other key (for example provenance or float diagnostics) is **not an input** and must be ignored.

## 2. Exact preconditions (checked in rational arithmetic)

Let τ = t_0, a = Re α, s = 1 − α.

1. 1/4 < τ < 1/2. A candidate will have τ < 1/3; τ ≥ 1/3 is allowed so that #251's certificate can serve as a test
   vector.
2. t_0 < t_1 < … < t_M and t_M = 1 − τ exactly.
3. len(eta) = M.
4. a > 0 and 0 < C < 1.
5. C² − |s|² > 0 and C² − |η_j|² > 0 for every j.

## 3. The analytic inequality (must be enclosed rigorously)

Put δ_0 = η_0 − s, δ_j = η_j − η_{j−1} (1 ≤ j ≤ M−1), and

* I_α(x) = ∫_0^x v^{α−1}/(1−v) dv for 0 < x < 1, where v^{α−1} = exp((α−1) log v) with real log v;
* K = I_α(τ) and D = I_a(τ) (real, > 0);
* T_2(u_1, u_2) = ∫_{x_1 ≥ u_1, x_2 ≥ u_2, x_1+x_2 ≤ 1} (1 − x_1 − x_2)^{α−1} / (x_1 x_2) dx, and 0 if u_1 + u_2 ≥ 1;
* T_3(u_1, u_2, u_3) = ∫_{x_i ≥ u_i, x_1+x_2+x_3 ≤ 1} (1 − x_1 − x_2 − x_3)^{α−1} / (x_1 x_2 x_3) dx, and 0 if
  u_1 + u_2 + u_3 ≥ 1.

Define, with every sum over **ordered** index tuples in {0, …, M−1}:

  Y = 1 + s·K − Σ_j δ_j·(I_α(1 − t_j) − K) + (1/2) Σ_{i,j} δ_i δ_j T_2(t_i, t_j) − (1/6) Σ_{i,j,k} δ_i δ_j δ_k T_3(t_i, t_j, t_k).

(Four-fold terms would need 4τ < 1, which precondition 1 excludes. T_3 vanishes identically when 3τ ≥ 1, so for
τ ≥ 1/3 this is Pinho's/Röhrig's Y with T_2 = J_α.)

**Verdict PASS iff** all of §2 holds and C²·D² − |Y|² > 0 is proved with outward-rounded interval/ball arithmetic,
including rigorous bounds on every series truncation or quadrature error.

## 4. What PASS proves

By the theorem in NOTES.md §4 (applied with d* = ⌈1/τ⌉ − 1 ≤ 3), PASS implies R_n < C for all sufficiently large n,
hence C42 ≤ C. No explicit threshold N is claimed.

## 5. Test vector

Pinho's #251 certificate (rpinho/c42-certificate@a8a9bc3, `certificates/c42_256.json`, SHA-256 705a7291…7afcab;
τ = 1/3, M = 256) must PASS. It gives

* D ∈ 0.892184209107241336566750826519 ± 1e-30
* |Y|/D ∈ 0.688886370983465645822767803130 ± 1e-30
* C²D² − |Y|² ≈ 2.3747257533e-10

The author's implementation is `check_tau.py` (python-flint Arb). T_3 there uses the Dirichlet series of NOTES.md §5.
An independent verifier should use a different method for T_3 where it can: for example the 1-D representation
T_3(u) = ∫_{u_3}^{1−u_1−u_2} (1−x)^{α−1}·T_2(u_1/(1−x), u_2/(1−x)) dx/x, or validated cubature.
