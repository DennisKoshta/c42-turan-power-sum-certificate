# Lowering the first-block cut-off in the coefficient-cancellation bound for Turán's pure power-sum constant

Dennis Koshta, with Claude (Anthropic, Claude Opus 5.5). Unrefereed note, 2026.

The section numbers match the references in `spec/FORMAT.md` and `spec/FORMAT2.md`: the theorem is §4, the T_3
series is §5, and the T_4 series and the T_3 integral are §7.

**Summary.** The coefficient-cancellation construction of Griego [Gri26], refined by A. Röhrig (PR #184), K. Plata and
R. Pinho (PR #251), proves C42 ≤ C from a rational step profile with first-block cut-off τ = 1/3. We show that the
underlying criterion holds for every τ ∈ (0, 1/2), with additional d-fold correction terms for d < 1/τ. The resulting
bound improves as τ decreases. At τ = 41/200 we get C42 ≤ 0.6873358063, and at τ = 13/50 C42 ≤ 0.6878759240.

## 1. Definition

R_n = min over z_1..z_n ∈ C with max_j |z_j| = 1 of max_{1≤k≤n} |Σ_j z_j^k|, and C42 = limsup_n R_n.

## 2. Exact finite-n reformulation

For S = (S_1..S_n) ∈ C^n let

  B_S(z) = (1−z)^{-1} exp(−Σ_{m≤n} S_m z^m/m) = Σ_l b_l z^l,  so k·b_k = Σ_{m=1}^k (1−S_m) b_{k−m}.

Suppose b_n = 0. Let y_2..y_n be the roots of Z^{n−1} + b_1 Z^{n−2} + … + b_{n−1}, and y_1 = 1. Newton's identities
give Σ_j y_j^k = S_k for k ≤ n; at k = n this uses b_n = 0. Dividing all y_j by Λ = max|y_j| ≥ 1 puts the maximum
modulus at 1 and does not increase any |power sum|. So R_n ≤ max_k |S_k|.

Conversely, after a rotation a configuration has z_1 = 1, and then b_n = 0. Hence

  R_n = min { max_{k≤n} |S_k| : S ∈ C^n, b_n(S) = 0 }.

Only the inequality "≤" is used below; Step 1 of §4 proves it.

## 3. The construction and the criterion at τ = 1/3 (Griego, Röhrig, Plata, Pinho)

Parameters: α = a + ib with a > 0; s = 1 − α; a cut-off τ; breaks τ = t_0 < … < t_M = 1 − τ; and middle values η_j.
Put A = ⌊τn⌋ and split 1..n into E = {1..A}, Mid = {A+1..n−A−1} and F = {n−A..n}.

- On E, S_m = s. On Mid, S_m = η_j on the j-th step.
- On F, the values are chosen last, to force b_n = 0. F enters b_n only linearly, so this is possible, and it makes
  |S_m| on F equal to |Y_n|/D_n (see §4).

Expanding the exponential, [z^n]B is a sum over tuples of "correction" indices m_1..m_d > A with Σ m_i ≤ n. Each
correction index exceeds τn, so d < 1/τ.

At τ = 1/3, at most two corrections reach b_n. The limit Y then needs only the single integral I_α and the double
integral J_α:

  Y = 1 + sK − Σ_j δ_j (I_α(1−t_j) − K) + ½ Σ_{i,j} δ_i δ_j J_α(t_i, t_j).

The criterion |s| < C, |η_j| < C, |Y| < C·D gives C42 ≤ C (Pinho, Prop. 1). Pinho's paper names the bound of two
corrections as the reason for τ = 1/3. The general statement and its proof follow.

## 4. General cut-off τ ∈ (0, 1/2): theorem and proof

This section is self-contained.

**Data.** τ ∈ (0, 1/2); α = a + ib with a > 0; s = 1 − α; rationals τ = t_0 < t_1 < … < t_M = 1 − τ; and complex
η_0..η_{M−1}. Put δ_0 = η_0 − s and δ_j = η_j − η_{j−1}. Then the step function

  w(x) = Σ_j δ_j·1[x > t_j]   (τ < x < 1−τ)

equals η_j − s on (t_j, t_{j+1}]. Let d* = ⌈1/τ⌉ − 1, the largest integer d < 1/τ. Define

  I_α(x) = ∫_0^x v^{α−1}/(1−v) dv (0 < x < 1),   K = I_α(τ),   D = I_a(τ) > 0,

  T_d(u_1..u_d) = ∫_{x_i ≥ u_i, Σx_i ≤ 1} (1 − Σx_i)^{α−1} / (x_1⋯x_d) dx   (= 0 if Σu_i ≥ 1),

  Y = 1 + sK − Σ_j δ_j (I_α(1−t_j) − K) + Σ_{d=2}^{d*} ((−1)^d / d!) Σ_{j_1..j_d} δ_{j_1}⋯δ_{j_d} T_d(t_{j_1}, …, t_{j_d}),

where the inner sum runs over ordered d-tuples. Powers of positive reals use the real logarithm.

**Theorem.** If |s| < C, |η_j| < C for all j, and |Y| < C·D, then R_n < C for all sufficiently large n. Hence
C42 ≤ C.

For τ = 1/3 we have d* = 2 and T_2 = J_α, so Y is exactly the τ = 1/3 expression of §3.

### Proof

**Step 1 (power sums from coefficients).** Fix n and S ∈ C^n. Let B(z) = (1−z)^{−1} exp(−Σ_{m≤n} S_m z^m/m) =
Σ b_l z^l, and suppose b_n = 0.

- Factor Z^{n−1} + b_1 Z^{n−2} + … + b_{n−1} = Π_{j=2}^n (Z − y_j). This gives Π_{j≥2}(1 − y_j z) = Σ_{l<n} b_l z^l.
  Since b_n = 0, this agrees with B(z) mod z^{n+1}.
- Put y_1 = 1 and P(z) = Π_{j=1}^n (1 − y_j z). Then P(z) ≡ (1−z)B(z) = exp(−Σ_{m≤n} S_m z^m/m) mod z^{n+1}.
- Both sides have constant term 1, so taking logarithms gives −Σ_k p_k z^k/k ≡ −Σ_{m≤n} S_m z^m/m, with
  p_k = Σ_j y_j^k. Hence p_k = S_k for k ≤ n.
- Let Λ = max|y_j| ≥ |y_1| = 1. The points z_j = y_j/Λ have maximum modulus 1, and |Σ_j z_j^k| = |S_k|/Λ^k ≤ |S_k|.

So R_n ≤ max_{k≤n} |S_k|.

**Step 2 (the sequence S).** Let n be large enough that every step (t_j, t_{j+1}] contains an integer multiple of
1/n. Put A = ⌊τn⌋, E = {1..A}, Mid = {A+1..n−A−1} and F = {n−A..n}; these are nonempty because τ < 1/2.

Set S_m = s on E, and S_m = s + w_m on Mid with w_m = Σ_j δ_j·1[m > ⌊n t_j⌋], so S_m ∈ {η_j}. F is fixed below.

Let β_l = [z^l](1−z)^{−α} = α(α+1)⋯(α+l−1)/l!. This is nonzero because α is not an integer ≤ 0. Since
−s Σ_{m≤A} z^m/m = s·log(1−z) + s Σ_{m>A} z^m/m, we get

  B(z) ≡ (1−z)^{−α} exp(Σ_{A<m≤n} c_m z^m)  mod z^{n+1},   c_m = (s − S_m)/m.

The terms s z^m/m with m > n do not affect coefficients up to z^n. Expanding the exponential,

  b_n = Σ_{d≥0} (1/d!) Σ_{m_1..m_d ∈ (A, n], Σm_i ≤ n} c_{m_1}⋯c_{m_d} β_{n−Σm_i}.

- Each m_i ≥ A + 1 > τn, so a nonzero term needs d(A+1) ≤ n, hence d < 1/τ, i.e. d ≤ d*. This holds already at
  finite n, so no finite-n term is lost.
- If d ≥ 2 and some m_i ∈ F, then Σm_i ≥ (n−A) + (A+1) = n + 1. So F occurs only in d = 1, linearly. This holds for
  every τ.

Thus n·b_n = Y_n − Σ_{m∈F} S_m P_m, where P_m = n β_{n−m}/m ≠ 0 and

  Y_n = n β_n + n Σ_{d=1}^{d*} (1/d!) Σ_{m_i ∈ Mid, Σm_i ≤ n} Π_i (−w_{m_i}/m_i) β_{n−Σm_i} + s Σ_{m∈F} P_m.

For m ∈ F choose S_m = (Y_n / D_n)·conj(P_m)/|P_m|, with D_n = Σ_{m∈F} |P_m| (S_m = 0 if Y_n = 0). Then
Σ_F S_m P_m = Y_n, so b_n = 0, and |S_m| = |Y_n|/D_n on F. By Step 1,
R_n ≤ max(|s|, max_j |η_j|, |Y_n|/D_n).

**Step 3 (limits).** We show Γ(α) n^{−α} Y_n → Y and |Γ(α)| n^{−a} D_n → D. Then |Y_n|/D_n → |Y|/D < C, so R_n < C
for all large n. Two facts about β are used:

- (i) β_l = Γ(l+α)/(Γ(α) l!) = l^{α−1}(1 + ε_l)/Γ(α) with ε_l → 0;
- (ii) |β_l| ≤ B_0 (l+1)^{a−1} for all l ≥ 0, for some constant B_0.

Write k_n(l) = Γ(α) n^{1−α} β_l. By (i), k_n(l) − (l/n)^{α−1} → 0 uniformly on l ≥ εn, for each fixed ε > 0. By
(ii), Σ_{l<εn} |k_n(l)| ≤ B_1 n^{1−a} (εn+1)^a.

- **Constant term.** Γ(α) n^{−α}·n β_n = k_n(n) → 1.
- **D_n.** |Γ(α)| n^{−a} D_n = (1/n) Σ_{l=0}^{A} |k_n(l)|·n/(n−l), with l = n − m.
  - The part with l < εn is O(n^{−1}·n^{1−a}(εn+1)^a) = O(ε^a + n^{−a}).
  - The rest is a Riemann sum of the continuous v^{a−1}/(1−v) on [ε, τ], with the integrand converging uniformly.
  - The corresponding piece ∫_0^ε of D is O(ε^a). Let n → ∞, then ε → 0.
- **The F term s Σ_F P_m.** The same argument with k_n in place of |k_n| gives s ∫_0^τ v^{α−1}/(1−v) dv = sK.
- **The d-fold term (1 ≤ d ≤ d*).** With x_i = m_i/n and l = n − Σm_i,
  Γ(α) n^{−α}·n Π(−w_{m_i}/m_i) β_l = n^{−d}·Π_i(−w_{m_i}/x_i)·k_n(l). This is a Riemann sum with cell volume
  n^{−d} over the lattice points of the polytope {x_i > A/n, x_i < 1 − (A+1)/n, Σx_i ≤ 1}.
  - Tuples with l < εn: for each l there are at most n^{d−1} tuples, and |Π(−w/x_i)| ≤ (W/τ)^d with W = max|w|.
    Their total is at most (W/τ)^d n^{−1} Σ_{l<εn} |k_n(l)| = O(ε^a + n^{−a}).
  - The matching part of the integral, over 1−ε ≤ Σx_i ≤ 1, is O(ε^a): integrate (1−Σx)^{a−1} across the slices
    Σx_i = const, whose (d−1)-volume is bounded.
  - On Σx_i ≤ 1−ε, k_n(l) → (1−Σx_i)^{α−1} uniformly. The limit function is bounded and continuous off the finitely
    many hyperplanes x_i = t_j and x_i = τ. The discontinuity set and the polytope's boundary have measure 0, so the
    lattice Riemann sums converge to the integral.
  - Floors cause no error: for an integer m, m > ⌊y⌋ ⇔ m > y, so w_m = w(m/n) exactly.

  So the d-fold term tends to (1/d!) ∫_{x_i ∈ (τ,1−τ), Σx ≤ 1} Π_i(−w(x_i)/x_i)·(1−Σx_i)^{α−1} dx.
- **Identification.** Expand Π_i(−w(x_i)) = (−1)^d Σ_{j_1..j_d} Π δ_{j_i}·1[x_i > t_{j_i}].
  - For d = 1: −δ_j ∫_{t_j}^{1−τ} (1−x)^{α−1}/x dx = −δ_j (I_α(1−t_j) − K), by v = 1 − x.
  - For d ≥ 2: the cap x_i < 1−τ is automatic, because the other coordinates sum to more than τ. This gives
    ((−1)^d/d!) Σ Πδ·T_d(t_{j_1}..t_{j_d}).

  Together with 1 + sK, this is Y. ∎

Step 3 is Pinho's limit argument written out for each d ≤ d*. Nothing in Steps 1–3 uses τ = 1/3.

## 5. Evaluating T_3 by a Dirichlet series (1/4 < τ < 1/3)

Put x_i = t_i + p_i, with L = 1 − t_1 − t_2 − t_3 > 0. Since p_i ≤ L < τ ≤ t_i, we can expand
1/(t_i + p_i) = Σ_k (−p_i)^k / t_i^{k+1}. Each monomial integrates by the Dirichlet formula

  ∫_{Σp≤L} Π p_i^{k_i} (L − Σp)^{α−1} dp = L^{|k|+α+2} Γ(α) Π k_i! / Γ(|k| + α + 3).

Hence

  T_3 = L^{α+2} Σ_{N≥0} (−L)^N h_N/(α)_{N+3},   h_N = Σ_{k_1+k_2+k_3=N} Π_i k_i! / t_i^{k_i+1} ∈ Q.

The only transcendental factor is L^{α+2}.

**Interchange.** Σ_k Π_i p_i^{k_i}/t_i^{k_i+1} = Π_i 1/(t_i − p_i) ≤ Π_i 1/(t_i − L) is bounded on the simplex, and
|(L − Σp)^{α−1}| = (L − Σp)^{a−1} is integrable there. So Tonelli/Fubini allows termwise integration.

**Tail.** Put q = L/t_min < 1.

- Σ_{|k|=N} Π k_i! ≤ C(N+2,2)·N!, so h_N ≤ C(N+2,2)·N!·t_min^{−N}/(t_1 t_2 t_3).
- |(α)_{N+4}| = |α+N+3|·|(α)_{N+3}| ≥ (N+a+3)·|(α)_{N+3}|.
- So the bound B_N on the N-th term satisfies B_{N+1}/B_N ≤ q·(N+3)/(N+a+3) ≤ q, and
  Σ_{N≥N_0} |term_N| ≤ B_{N_0}/(1−q), where

  B_{N_0} = L^{a+2} q^{N_0} C(N_0+2,2) N_0! / (t_1 t_2 t_3 |(α)_{N_0+3}|).

**Scaling identity (any τ).** T_d(t) = ∫_{t_d}^{1−Σ_{i<d} t_i} (1−x)^{α−1}·T_{d−1}(t_1/(1−x), …, t_{d−1}/(1−x)) dx/x,
since Π dx_i/x_i is scale-invariant.

J_α = T_2 for pairs with t + r < 1 uses Pinho's log series. Its ratio L/(1−t) ≤ (1−2τ)/(1−τ) stays below 1 for every
τ < 1/2.

## 6. The certificates

### Primary: τ = 41/200, M = 64, d* = 4

File `certificates/c42_tau41_200_M64_0.6873358063.json` (SHA-256 be0c406f…cbef), with 64 uniform middle blocks on
[41/200, 159/200] and C = 6873358063/10000000000.

- `checker/check_tau4.py` (Arb, 384 and 640 bits): PASS, |Y|/D = 0.687335806186516616249961526259 ± 1e−31, and
  C²D² − |Y|² = 4.348178358e−11.
  - 1056 T_2 pairs, 2000 T_3 by the series of §7.1, 443 T_3 by the integral of §7.2, and 609 T_4 by the series.
- `verifier/turan_c42.py`, an independently written verifier (160 and 320 bits): VALID, with the same |Y|/D to 30
  digits. It uses its own series about the opposite simplex face, not the series of §5/§7.1. The four-fold term is
  Σ δδδδ T_4 = 1.214e−5 − 3.45e−6i, and it counts against the claim.
- Cross-check without any T_d: the exact finite-n recursion (`crosscheck/finite_n.py`) on the same rational profile
  gives 0.6875879, 0.6874327, 0.6873730 and 0.6873501 at n = 256k, 1.02M, 4.10M and 16.4M. Fitting c + b n^−p to the
  last three (p ≈ 0.69) gives 0.68733577, which is 3e−8 from the continuum value.
- M-convergence at τ = 41/200 (float optimum): 0.6873500, 0.6873358 and 0.6873323 at M = 32, 64 and 128. The error
  goes like M^−2, and the extrapolated limit is ≈ 0.6873311.

### Fallback: τ = 13/50, M = 64, d* = 3

File `certificates/c42_tau26_M64_0.6878759240.json` (SHA-256 9cbc0019…1227), with C = 171968981/250000000.

- `checker/check_tau.py` (384 and 640 bits): PASS, |Y|/D = 0.68787592366984836…, C²D² − |Y|² = 2.08e−10.
- The independent verifier (`--format 1`): VALID.
- Finite-n cross-check: 0.6878755.

## 7. 1/5 < τ < 1/4 (d* = 4)

By §4 with d* = ⌈1/τ⌉ − 1 = 4, Y gains the term +(1/24) Σ_{i,j,k,l} δ_i δ_j δ_k δ_l T_4(t_i, t_j, t_k, t_l) over
ordered 4-tuples. (τ = 1/4 exactly still gives d* = 3.)

**7.1 Dirichlet series for any d (when L < u_min).** For u ∈ (0,1)^d with L = 1 − Σu_i > 0 and L < u_min, the §5
argument gives, verbatim,

  T_d(u) = L^{α+d−1} Σ_{N≥0} (−L)^N h^{(d)}_N / (α)_{N+d},   h^{(d)}_N = Σ_{|k|=N} Π_i k_i!/u_i^{k_i+1} ∈ Q.

The tail bound is B_N = L^{a+d−1} q^N C(N+d−1, d−1) N!/(Πu_i·|(α)_{N+d}|) with q = L/u_min. Since
B_{N+1}/B_N = q(N+d)/|N+α+d| < q, we get Σ_{N≥N0} ≤ B_{N0}/(1−q).

- For d = 4 and τ > 1/5, every 4-tuple has L ≤ 1 − 4τ < τ ≤ u_min, so q ≤ (1−4τ)/τ < 1. All T_4 use 7.1.
- For d = 3 at τ < 1/4, L ≤ 1 − 3τ can exceed u_min. So 7.1 is used only when q ≤ 4/5, and 7.2 otherwise.

**7.2 T_3 when q is near or above 1: a 1-D integral.** By the scaling identity of §5,

  T_3(u) = ∫_{u_3}^{ℓ} (1−x)^{α−1} J(u_1/(1−x), u_2/(1−x)) dx/x,   ℓ = 1 − u_1 − u_2.

With t = u_1/(1−x) and r = u_2/(1−x), we have 1 − t − r = (ℓ−x)/(1−x) =: λ(x). J's series is
J(t,r) = λ^α Σ_k c_k λ^k/(α+k), with c_k = c_0 − Σ_{q≤k} (b^{−q} + d^{−q})/q, c_0 = log(bd/(tr)), b = 1−t, d = 1−r.
The factors (1−x)^{α−1}·(1−x)^{−α} combine to (1−x)^{−1}. Substitute x = ℓ − W e^{−v} with W = ℓ − u_3 = L and
v ∈ [0, ∞). Then

  T_3(u) = W^{α+1} ∫_0^∞ e^{−(α+1)v} · Φ(x(v)) dv,   Φ(x) = [Σ_k c_k(x) λ(x)^k/(α+k)] / (x(1−x)).

- **Analyticity.** On the real path we have x ∈ [u_3, ℓ) ⊂ (0, 1), t + r < 1, 0 < λ < b, d < 1, and bd/(tr) > 1
  (since bd − tr = 1 − t − r > 0). So Φ is analytic in a complex neighbourhood of the path. The only branch cut is
  from log in c_0, and it is handled with Arb's `log(analytic=True)`.
- **The k-series on complex balls.** It converges when |λ| < min(1, |b|, |d|). With ρ_b = max(1, 1/|b|) (and ρ_d
  likewise), |c_k| ≤ |c_0| + k(ρ_b^k + ρ_d^k). The tail beyond K terms is then at most
  (1/(a+K))·[|c_0||λ|^K/(1−|λ|) + Σ_{y ∈ {|λ|ρ_b, |λ|ρ_d}} y^K (K/(1−y) + y/(1−y)²)], provided y < 1. The checker
  certifies y < 1 on the whole ball.
- **Endpoint truncation.** Integrate v over [0, V] with `acb.integral`, and bound the rest by
  |∫_V^∞| ≤ |W^{α+1}| e^{−(a+1)V}/(a+1) · sup{|Φ(x)| : x ∈ [ℓ − We^{−V}, ℓ]}. The sup is one Arb evaluation of Φ on
  a real ball covering that interval.

Unit checks (Arb): the 7.2 integral and the 7.1 series agree to 1e−42 on single triples with q = 0.69, 0.80 and 0.85.
The T_4 series agrees to 1e−16 with quadrature of the 1-D scaling identity at q = 0.80 and 0.88.

## References

- [Gri26] S. Griego, *An improved asymptotic certificate for Turán's pure power sum constant C42*,
  https://github.com/sebastian-griego/turan-c42-certificate, tag v1.0.0, commit c8ddce14d9a5e898406d5dc6b8d08bb8a39507c7.
- A. Röhrig, teorth/optimizationproblems PR #184 (C42 certificate, 0.688983).
- R. Pinho, teorth/optimizationproblems PR #251, and https://github.com/rpinho/c42-certificate at commit
  a8a9bc309bd1992413fa07595bd60c5e3fa1b65b (paper `c42.tex`, certificate `certificates/c42_256.json`, 0.6888863712).
- K. Plata, six-level τ = 1/3 refinement (0.688970), https://github.com/krplatz/turan-c42-six-level-review at commit
  d3650caf944e8602e867e2e81a439e037a7b91a9 (as cited by Pinho).
