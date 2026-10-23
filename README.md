# C42: Turán's pure power-sum constant, upper bound 0.6873358063 (cut-off τ = 41/200)

**Claim.** C₄₂ ≤ 6873358063/10000000000 = 0.6873358063.

C₄₂ is constant 42 of [teorth/optimizationproblems](https://teorth.github.io/optimizationproblems/constants/42a.html):
C₄₂ = limsup_n R_n, with R_n = min over z_1..z_n ∈ ℂ, max|z_i| = 1, of max_{1≤k≤n} |Σ_i z_i^k|.

This is an unrefereed analytic asymptotic reduction with an exact rational certificate. The checker evaluates the
limiting inequality in rigorous ball arithmetic (Arb). We do not claim an explicit finite threshold N, optimality,
or formal verification.

**What is new.** Earlier certificates in this family use the first-block cut-off τ = 1/3, which keeps the limit to
two correction terms:

- Griego [Gri26], 0.6906538;
- Röhrig, PR #184, 0.688983;
- Plata, 0.688970;
- Pinho, PR #251, 0.6888863712.

`NOTES.md` §4 proves that their criterion holds for every τ ∈ (0, 1/2), with d-fold correction terms
((−1)^d/d!) Σ δ…δ T_d for d < 1/τ. The bound improves as τ decreases. At τ = 41/200 the terms
−(1/6) Σ δδδ T₃ and +(1/24) Σ δδδδ T₄ appear. They are evaluated by:

- exact-rational Dirichlet series with explicit tail bounds (`NOTES.md` §5, §7.1);
- a validated 1-D integral for the T₃ where the series does not converge fast enough (`NOTES.md` §7.2).

## Contents

| path | what |
|---|---|
| `certificates/c42_tau41_200_M64_0.6873358063.json` | **primary** certificate, τ = 41/200, M = 64 |
| `certificates/c42_tau26_M64_0.6878759240.json` | fallback certificate, τ = 13/50, M = 64 (C₄₂ ≤ 0.6878759240; needs only T₃) |
| `spec/FORMAT2.md` | certificate format and exact criterion for 1/5 < τ < 1/2 (governs the primary) |
| `spec/FORMAT.md` | the same for 1/4 < τ < 1/2 (governs the fallback) |
| `NOTES.md` | the mathematics: theorem and proof (§4), T_d series (§5, §7.1), T₃ integral (§7.2), results (§6) |
| `checker/check_tau4.py` | Arb checker for FORMAT2 (python-flint) |
| `checker/check_tau.py` | Arb checker for FORMAT |
| `checker/test_check_tau4.py`, `checker/test_check_tau.py` | rejection controls |
| `verifier/turan_c42.py` | independently written verifier (see Credits) |
| `crosscheck/finite_n.py`, `crosscheck/fn_profile.py` | exact finite-n recursion on the same profile, extrapolated in n (float; uses no T_d) |
| `expected/` | the checker outputs we obtained (JSON) |
| `requirements.txt` | python-flint 0.9.0 (checkers, verifier); numpy and scipy (cross-check only) |

The two spec files are reproduced verbatim as frozen during the work. Some of their lines refer to the original
repository layout:

- `NOTES.md §n` is this repository's `NOTES.md` (same section numbers).
- `records/c42_turan/…json` is `certificates/…json`.
- `REVIEW_tau_lt_quarter.md` is an internal review; its points are listed in FORMAT2.md itself.
- "message to the lead" is our internal process.

## Reproduce

Use Python 3.12 with `pip install -r requirements.txt`, and run from the repository root:

```
# primary: PASS, |Y|/D = 0.687335806186516616249961526259..., C^2D^2-|Y|^2 = 4.348178358e-11 (about 8 min, 1 core)
python -I checker/check_tau4.py certificates/c42_tau41_200_M64_0.6873358063.json --bits 384 --eps 1e-45
# same at 640 bits (slower)
python -I checker/check_tau4.py certificates/c42_tau41_200_M64_0.6873358063.json --bits 640 --terms 900 --eps 1e-70
# fallback: PASS, |Y|/D = 0.68787592366984836...
python -I checker/check_tau.py certificates/c42_tau26_M64_0.6878759240.json --bits 384
# independent verifier (exit 0 = VALID and BEATS)
python -I verifier/turan_c42.py certificates/c42_tau41_200_M64_0.6873358063.json --beat 0.6888863712
python -I verifier/turan_c42.py certificates/c42_tau26_M64_0.6878759240.json --format 1 --beat 0.6888863712
# rejection controls
python -I checker/test_check_tau4.py certificates/c42_tau41_200_M64_0.6873358063.json
python -I checker/test_check_tau.py certificates/c42_tau26_M64_0.6878759240.json
```

**Test vector (not redistributed).** Pinho's #251 certificate is `certificates/c42_256.json` from
https://github.com/rpinho/c42-certificate at commit a8a9bc309bd1992413fa07595bd60c5e3fa1b65b (SHA-256 705a7291…7afcab).
Both checkers must PASS it with |Y|/D = 0.688886370983465645822767803130:

```
python -I checker/check_tau4.py /path/to/c42_256.json --bits 384 --eps 1e-45
```

**Optional cross-check** (float, no T_d; the largest n takes hours):

```
python -I crosscheck/fn_profile.py certificates/c42_tau41_200_M64_0.6873358063.json 41/200 64 256000,1024000,4096000,16384000
```

This gives 0.6875879, 0.6874327, 0.6873730 and 0.6873501, which extrapolate to about 0.68733577.

## Credits

- **The construction and the τ = 1/3 criterion:**
  - S. Griego, [turan-c42-certificate v1.0.0](https://github.com/sebastian-griego/turan-c42-certificate/tree/c8ddce14d9a5e898406d5dc6b8d08bb8a39507c7);
  - A. Röhrig, [PR #184](https://github.com/teorth/optimizationproblems/pull/184);
  - K. Plata, [turan-c42-six-level-review](https://github.com/krplatz/turan-c42-six-level-review/tree/d3650caf944e8602e867e2e81a439e037a7b91a9);
  - R. Pinho, [PR #251](https://github.com/teorth/optimizationproblems/pull/251) and
    [c42-certificate](https://github.com/rpinho/c42-certificate/tree/a8a9bc309bd1992413fa07595bd60c5e3fa1b65b).

  Our contribution is the general-τ extension of their criterion and its proof, the T₃/T₄ evaluation, and the
  certificates. No third-party code is included in this repository.
- **AI use.** This result was found and prepared with Claude (Anthropic, Claude Opus 5.5), directed by Dennis
  Koshta, who takes responsibility for it. That covers the general-cut-off extension and its proof, the parameter
  search, the certificates, the Arb checkers and the text.
  - The proof was reviewed by two further Claude agents, working from the note and Pinho's paper only.
  - The certificates were checked by a separate Claude agent with its own verifier (`verifier/turan_c42.py`), written
    from the format specification without access to the search code.
  - No human mathematician has checked the argument.

## Licence

MIT; see `LICENSE`. Copyright (c) 2026 Dennis Koshta.
