# C42 certificate for 1/5 < τ < 1/2: format and exact criterion (FROZEN)

This file is frozen, after the §7 review (7287498f, `REVIEW_tau_lt_quarter.md`) was applied. It extends FORMAT.md
(frozen, τ > 1/4) by one term, T_4. Nothing in FORMAT.md or the frozen τ = 13/50 candidate changes. An independent
verifier is written from this file, NOTES.md §4–5 and §7, and the papers. Changes after the freeze commit need a new
commit, a note here, and a message to the lead (5697e7f3).

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

1. 1/5 < τ < 1/2. A candidate will have τ < 1/4; larger τ is allowed so that #251's certificate and the frozen
   τ = 13/50 candidate serve as test vectors.
2. t_0 < t_1 < … < t_M and t_M = 1 − τ exactly.
3. len(eta) = M.
4. a > 0 and 0 < C < 1.
5. C² − |s|² > 0 and C² − |η_j|² > 0 for every j.

## 3. The analytic inequality (must be enclosed rigorously)

Put δ_0 = η_0 − s, δ_j = η_j − η_{j−1} (1 ≤ j ≤ M−1), and

* I_α(x) = ∫_0^x v^{α−1}/(1−v) dv for 0 < x < 1, where v^{α−1} = exp((α−1) log v) with real log v;
* K = I_α(τ) and D = I_a(τ) (real, > 0);
* T_2(u_1, u_2) = ∫_{x_1 ≥ u_1, x_2 ≥ u_2, x_1+x_2 ≤ 1} (1 − x_1 − x_2)^{α−1} / (x_1 x_2) dx, and 0 if u_1 + u_2 ≥ 1;
* for d = 3, 4: T_d(u) = ∫_{x_i ≥ u_i, Σx_i ≤ 1} (1 − Σx_i)^{α−1} / Πx_i dx (over R^d), and 0 if Σu_i ≥ 1.

Define, with every sum over **ordered** index tuples in {0, …, M−1}:

  Y = 1 + s·K − Σ_j δ_j·(I_α(1 − t_j) − K) + (1/2) Σ_{i,j} δ_i δ_j T_2(t_i, t_j) − (1/6) Σ_{i,j,k} δ_i δ_j δ_k T_3(t_i, t_j, t_k)
        + (1/24) Σ_{i,j,k,l} δ_i δ_j δ_k δ_l T_4(t_i, t_j, t_k, t_l).

(Five-fold terms would need 5τ < 1, which precondition 1 excludes. T_4 vanishes when 4τ ≥ 1 and T_3 when 3τ ≥ 1, so
for τ > 1/4 this is FORMAT.md's Y, and for τ ≥ 1/3 it is Pinho's/Röhrig's Y with T_2 = J_α.)

**Verdict PASS iff** all of §2 holds and C²·D² − |Y|² > 0 is proved with outward-rounded interval/ball arithmetic,
including rigorous bounds on every series truncation or quadrature error.

**Rigour points** (from the §7 review, `REVIEW_tau_lt_quarter.md`, 7287498f). Any implementation must handle these:

* d* = ⌈1/τ⌉ − 1 must be computed exactly. Precondition 1 gives d* ≤ 4; for τ ≤ 1/5 the T_5 terms would be lost
  silently.
* Whether Σu ≥ 1 (so T_d = 0) must be decided in exact rationals. Breakpoint sums can equal 1 exactly.
* Geometric tail bounds of the form 1/(1 − y) need y < 1 certified for the whole ball, or from exact rationals. In
  ball arithmetic 1/(1 − y) is finite but meaningless when y > 1.
* A tail bound "from N_0 on" must be added to the partial sum of terms 0..N_0−1.
* The integrand handed to a complex-ball integrator must agree with the real integrand on the path and be holomorphic
  on each ball it accepts. Do not evaluate fractional powers of complex x (such as (1−x)^{α−1}) inside it, and handle
  log's branch cut (for example with Arb's analytic flag).
* I_α(1 − t_j) is needed at arguments up to 1 − τ < 4/5. Tail bounds tuned for x ≤ 2/3 do not apply.

## 4. What PASS proves

By the theorem in NOTES.md §4 (applied with d* = ⌈1/τ⌉ − 1 ≤ 4), PASS implies R_n < C for all sufficiently large n,
hence C42 ≤ C. No explicit threshold N is claimed.

## 5. Test vector

Pinho's #251 certificate (rpinho/c42-certificate@a8a9bc3, `certificates/c42_256.json`, SHA-256 705a7291…7afcab;
τ = 1/3, M = 256) must PASS. It gives

* D ∈ 0.892184209107241336566750826519 ± 1e-30
* |Y|/D ∈ 0.688886370983465645822767803130 ± 1e-30
* C²D² − |Y|² ≈ 2.3747257533e-10

The frozen candidate `records/c42_turan/c42_tau26_M64_0.6878759240.json` (SHA-256 9cbc0019…1227) must also PASS,
with |Y|/D = 0.68787592366984836… as in FORMAT.md's checker.

The author's implementation is `check_tau4.py` (python-flint Arb). It uses the Dirichlet series of NOTES.md §7.1 for
every T_4 and for T_3 with q = (1 − Σu)/min u ≤ 4/5, and the 1-D integral of §7.2 for the other T_3. An independent
verifier should use a different method where it can, e.g. the 1-D scaling identity
T_d(u) = ∫_{u_d}^{1−Σ_{i<d}u_i} (1−x)^{α−1}·T_{d−1}(u_1/(1−x), …, u_{d−1}/(1−x)) dx/x, or validated cubature.
