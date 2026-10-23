#!/usr/bin/env python3
"""Independent verifier for asymptotic UPPER-bound certificates of Turan's pure power-sum
constant C_42 (optimizationproblems constant 42a), for cut-offs 1/5 < tau < 1/2 (FORMAT2.md,
frozen at faec3a5; default) or 1/4 < tau < 1/2 (FORMAT.md, frozen at 3031c2c; --format 1).
For tau < 1/4, d* = ceil(1/tau) - 1 = 4 and Y gains + 1/24 sum_{ijkl} d_i d_j d_k d_l T_4 (ordered).
T_4 uses the same face expansion as T_3, built recursively (see the comment above Evaluator.T4).

Definition (constants/42a.md): C_42 = limsup_n R_n,
    R_n = min_{z_1..z_n in C, max|z_i| = 1}  max_{1<=k<=n} |sum_i z_i^k|.

Certificate (problems/turan_c42/FORMAT.md, frozen at 3031c2c): a JSON object with keys
    "constant" = "C_42", "C", "alpha" = [Re, Im], "breaks" = [t_0..t_M], "eta" = [[Re, Im] x M];
every number a JSON string holding an exact rational: a terminating decimal "-?d+(.d+)?" or a
fraction "-?p/q" (q > 0; FORMAT says p, q > 0, the sign is accepted so negative fractions can be
written). Floats, exponents and bare JSON numbers are refused; other keys are ignored;
duplicate keys anywhere are refused.

Criterion (theorem of NOTES.md section 4, re-derived and checked independently; summary):
with tau = t_0, a = Re alpha, s = 1 - alpha, delta_0 = eta_0 - s, delta_j = eta_j - eta_{j-1},
    I(x) = int_0^x v^(alpha-1)/(1-v) dv,   K = I(tau),   D = int_0^tau v^(a-1)/(1-v) dv,
    T_d(u) = int_{x_i >= u_i, sum x <= 1} (1 - sum x)^(alpha-1) / prod x_i dx,
    Y = 1 + sK - sum_j delta_j (I(1-t_j) - K) + 1/2 sum_{ij} d_i d_j T_2 - 1/6 sum_{ijk} d_i d_j d_k T_3,
if |s| < C, |eta_j| < C and |Y| < C D then R_n < C for all large n, so C_42 <= C. For tau > 1/4
at most three correction indices fit below n, so T_4.. never enter; T_3 = 0 when 3 tau >= 1.

Evaluation (written from the definitions; no code of the certificate's author was read).
Polar coordinates on the simplex, p_i = x_i - u_i, rho = sum p_i, L = 1 - sum u_i:
    T_d(u) = int_0^L (L - rho)^(alpha-1) F_d(rho) d rho,
    F_d(rho) = int_{p >= 0, sum p = rho} prod_i 1/(u_i + p_i)    ((d-1)-dim, coordinates p_1..p_{d-1}).
  d = 1: F_1 = 1/(u + rho).
  d = 2: partial fractions (X + Y = U + rho):  F_2(rho) = [log(1 + rho/u_1) + log(1 + rho/u_2)] / (U + rho).
  d = 3: 1/(X1 X2 X3) = (sum of the three 1/(Xj Xk)) / (U + rho), and
         N_jk(rho) := int_{p_j + p_k <= rho} dp/((u_j+p_j)(u_k+p_k)) has N_jk' = F_2^{jk}, N_jk(0) = 0, and
         N_jk(rho) = log(c/u_k) log(1 + rho/u_j) - Li2((u_j+rho)/c) + Li2(u_j/c),  c = u_j+u_k+rho;
         F_3(rho) = sum_{pairs} N_jk(rho) / (U + rho).
Every F_d is analytic off the negative real axis; its singularities nearest to rho = L are at
rho = -u_min, so in v = L - rho its Taylor series at v = 0 has radius A = L + u_min > L. Expanding
there (the face sum x = 1 of the simplex, the opposite end to the vertex x = u) and integrating
v^(alpha-1) v^n exactly gives
    T_d = L^alpha sum_n G_n L^n / (n + alpha),
with REAL coefficients G_n (recurrences below) and an explicit geometric tail bound (ratio L/A,
which is < 1/2 whenever tau > 1/4). U + L = 1 makes 1/(U + rho) = 1/(1 - v).
I(x) uses x^alpha sum_k x^k/(k+alpha) with the tail bound x^(a+N) / ((N+a)(1-x)).

Arithmetic: preconditions in exact rationals (fractions.Fraction); everything analytic in Arb
balls (python-flint); the verdict is the ball test C^2 D^2 - |Y|^2 > 0 ("certainly positive").
No float decides anything (floats only choose truncation orders, whose tails are then bounded
rigorously).

Exit codes (repo convention): 0 VALID (and BEATS if --beat); 1 INVALID; 2 VALID but
DOES_NOT_BEAT; 3 usage error (bad arguments, python -O, flint missing).
"""

from __future__ import annotations

import argparse
import collections
import itertools
import json
import math
import re
import sys
import time
from fractions import Fraction

EXIT_VALID, EXIT_INVALID, EXIT_NOT_BEAT, EXIT_USAGE = 0, 1, 2, 3
MAX_BYTES = 64 << 20

_DEC = re.compile(r"-?[0-9]+(?:\.[0-9]+)?")
_FRAC = re.compile(r"(-?[0-9]+)/([0-9]+)")


class CertificateError(ValueError):
    """The certificate is malformed or violates a precondition."""


# ----------------------------------------------------------------------------- parsing

def parse_rational(x, what: str) -> Fraction:
    if not isinstance(x, str):
        raise CertificateError(f"{what}: expected a string holding an exact rational, got {type(x).__name__}")
    if _DEC.fullmatch(x):
        return Fraction(x)
    m = _FRAC.fullmatch(x)
    if m:
        q = int(m.group(2))
        if q == 0:
            raise CertificateError(f"{what}: zero denominator in {x!r}")
        return Fraction(int(m.group(1)), q)
    raise CertificateError(f"{what}: {x!r} is not a terminating decimal or p/q")


def _no_dupes(pairs):
    d = {}
    for k, v in pairs:
        if k in d:
            raise CertificateError(f"duplicate key {k!r}")
        d[k] = v
    return d


def _bad_constant(name):
    raise CertificateError(f"non-finite JSON constant {name}")


def load_certificate(path: str) -> dict:
    try:
        with open(path, "rb") as fh:
            raw = fh.read(MAX_BYTES + 1)
    except OSError as e:
        raise CertificateError(f"cannot read {path}: {e}") from None
    if len(raw) > MAX_BYTES:
        raise CertificateError("file too large")
    try:
        text = raw.decode("utf-8")
        obj = json.loads(text, object_pairs_hook=_no_dupes, parse_constant=_bad_constant)
    except CertificateError:
        raise
    except (UnicodeDecodeError, ValueError, RecursionError) as e:
        raise CertificateError(f"not valid UTF-8 JSON: {e}") from None
    if not isinstance(obj, dict):
        raise CertificateError("top level is not a JSON object")
    for key in ("constant", "C", "alpha", "breaks", "eta"):
        if key not in obj:
            raise CertificateError(f"missing key {key!r}")
    if obj["constant"] != "C_42":
        raise CertificateError(f"constant is {obj['constant']!r}, expected 'C_42'")

    def pair(v, what):
        if not isinstance(v, list) or len(v) != 2:
            raise CertificateError(f"{what}: expected [Re, Im]")
        return (parse_rational(v[0], what + "[0]"), parse_rational(v[1], what + "[1]"))

    C = parse_rational(obj["C"], "C")
    alpha = pair(obj["alpha"], "alpha")
    br = obj["breaks"]
    if not isinstance(br, list) or len(br) < 2:
        raise CertificateError("breaks: expected a list [t_0, ..., t_M] with M >= 1")
    breaks = [parse_rational(t, f"breaks[{i}]") for i, t in enumerate(br)]
    ev = obj["eta"]
    if not isinstance(ev, list):
        raise CertificateError("eta: expected a list")
    eta = [pair(e, f"eta[{j}]") for j, e in enumerate(ev)]
    return {"C": C, "alpha": alpha, "breaks": breaks, "eta": eta}


TAU_MIN = {1: Fraction(1, 4), 2: Fraction(1, 5)}   # FORMAT.md / FORMAT2.md precondition 1


def dstar(tau: Fraction) -> int:
    """d* = ceil(1/tau) - 1, the largest d with d*tau < 1, exactly."""
    inv = 1 / tau
    return -((-inv.numerator) // inv.denominator) - 1


def check_preconditions(cert: dict, fmt: int = 2) -> None:
    """FORMAT.md (fmt=1) / FORMAT2.md (fmt=2) section 2, exactly."""
    C, (a, b), t, eta = cert["C"], cert["alpha"], cert["breaks"], cert["eta"]
    tau, M = t[0], len(t) - 1
    lo = TAU_MIN[fmt]
    if not (lo < tau < Fraction(1, 2)):
        raise CertificateError(f"tau = t_0 = {tau} is not in ({lo}, 1/2)")
    if dstar(tau) > 4:  # cannot happen for tau > 1/5; kept so T_5 terms are never dropped silently
        raise CertificateError(f"d* = {dstar(tau)} > 4 is not supported")
    for i in range(M):
        if not t[i] < t[i + 1]:
            raise CertificateError(f"breaks not strictly increasing at index {i}")
    if t[M] != 1 - tau:
        raise CertificateError(f"t_M = {t[M]} != 1 - tau = {1 - tau}")
    if len(eta) != M:
        raise CertificateError(f"len(eta) = {len(eta)} != M = {M}")
    if not a > 0:
        raise CertificateError("Re alpha <= 0")
    if not (0 < C < 1):
        raise CertificateError("C not in (0, 1)")
    sr, si = 1 - a, -b
    if not C * C - (sr * sr + si * si) > 0:
        raise CertificateError("|s| = |1 - alpha| >= C")
    for j, (er, ei) in enumerate(eta):
        if not C * C - (er * er + ei * ei) > 0:
            raise CertificateError(f"|eta_{j}| >= C")


# ----------------------------------------------------------------------------- analysis

class Evaluator:
    """All analytic quantities, in Arb balls at the current flint precision."""

    def __init__(self, alpha, prec: int):
        from flint import acb, arb, ctx, fmpq
        self.acb, self.arb, self.fmpq = acb, arb, fmpq
        ctx.prec = prec
        ctx.threads = 1
        self.prec = prec
        self.target = 2.0 ** -(prec + 8)  # absolute tail target for each special value (float: chooses N only)
        self.a_q, self.b_q = alpha
        self.a = self.q(self.a_q)
        self.b = self.q(self.b_q)
        self.alpha = acb(self.a, self.b)
        self.af = float(self.a_q)
        self.unit = arb(0, 1)          # the ball [-1, 1]
        self._inv = []                 # (Re, Im) of 1/(n + alpha)

    def q(self, x: Fraction):
        return self.arb(self.fmpq(x.numerator, x.denominator))

    def cq(self, re: Fraction, im: Fraction):
        return self.acb(self.q(re), self.q(im))

    def inv_n_alpha(self, n):
        while len(self._inv) <= n:
            k = len(self._inv)
            re = k + self.a
            den = re * re + self.b * self.b
            self._inv.append((re / den, -self.b / den))
        return self._inv[n]

    def pow_alpha(self, x):
        """x^alpha for an arb ball x > 0 (real logarithm)."""
        return (self.alpha * self.acb(x.log())).exp()

    def _nterms(self, ratio: float, scale: float) -> int:
        """Smallest N with scale * ratio^N <= target (float, only picks the truncation order)."""
        if ratio <= 0.0 or scale <= 0.0:
            return 1
        need = math.log(self.target / scale) / math.log(ratio)
        return max(1, int(math.ceil(need)) + 2)

    # I(x) = int_0^x v^(alpha-1)/(1-v) dv = x^alpha sum_k x^k/(k+alpha), 0 < x < 1 rational.
    def I_alpha(self, x: Fraction, real_a: bool = False):
        xb = self.q(x)
        xf = float(x)
        N = self._nterms(xf, 1.0 / ((1 - xf) * self.af))
        re_s, im_s = self.arb(0), self.arb(0)
        p = self.arb(1)
        for k in range(N):
            if real_a:
                re_s += p / (k + self.a)
            else:
                ir, ii = self.inv_n_alpha(k)
                re_s += p * ir
                im_s += p * ii
            p *= xb
        # |tail| <= sum_{k>=N} x^k/(k+a) <= x^N / ((N+a)(1-x)); times |x^alpha| = x^a below
        tail = p / ((N + self.a) * self.q(1 - x))  # 1 - x > 0 exactly (rational)
        if real_a:
            xa = (self.a * xb.log()).exp()
            return xa * (re_s + self.unit * tail)
        return self.pow_alpha(xb) * self.acb(re_s + self.unit * tail, im_s + self.unit * tail)

    # T_2(u1, u2), u_i > 0 rational; 0 if u1 + u2 >= 1.
    def T2(self, u1: Fraction, u2: Fraction):
        Lq = 1 - u1 - u2
        if Lq <= 0:
            return self.acb(0)
        arb = self.arb
        L = self.q(Lq)
        U1, U2 = self.q(u1), self.q(u2)
        A1, A2 = L + U1, L + U2                  # = 1 - u2, 1 - u1
        lam = (A1 / U1).log() + (A2 / U2).log()  # real, > 0
        iA1, iA2 = 1 / A1, 1 / A2
        Af = float(min(1 - u1, 1 - u2))
        Lf = float(Lq)
        r = Lf / Af
        lamf = float(lam.mid())
        scale = Lf ** self.af / self.af * (lamf / (1 - Lf) + 2.0 / ((1 - Af) * (1 - r)))
        N = self._nterms(r, scale)
        # F_2(L - v) = [lam + log(1 - v/A1) + log(1 - v/A2)] / (1 - v) = sum_n e_n v^n,
        # e_0 = lam, e_n = e_{n-1} - (A1^-n + A2^-n)/n.
        e = lam
        p1, p2 = arb(1), arb(1)
        Ln = arb(1)
        re_s, im_s = arb(0), arb(0)
        for n in range(N):
            if n:
                p1 *= iA1
                p2 *= iA2
                e = e - (p1 + p2) / n
            ir, ii = self.inv_n_alpha(n)
            t = e * Ln
            re_s += t * ir
            im_s += t * ii
            Ln *= L
        # |e_n| <= lam + 2 A^-n / (1 - A) (A = min(A1, A2) < 1); |1/(n+alpha)| <= 1/(n+a), so
        # sum_{n>=N} |e_n| L^n / (n+a) <= [lam L^N/(1-L) + 2 (L/A)^N / ((1-A)(1-L/A))] / (N+a).
        Aq = min(1 - u1, 1 - u2)
        rq = Lq / Aq                             # < 1 exactly; 1 - L, 1 - A, 1 - r below are exact rationals
        tail = (lam * Ln / self.q(1 - Lq) + 2 * self.q(rq) ** N / self.q((1 - Aq) * (1 - rq))) / (N + self.a)
        return self.pow_alpha(L) * self.acb(re_s + self.unit * tail, im_s + self.unit * tail)

    @staticmethod
    def _li2(z):
        return z.polylog(2)

    def N_pair_at_L(self, uj, uk, L):
        """N_jk(L) = log(c/u_k) log(1 + L/u_j) - Li2((u_j+L)/c) + Li2(u_j/c), c = u_j+u_k+L (real)."""
        c = uj + uk + L
        return (c / uk).log() * ((uj + L) / uj).log() - self._li2((uj + L) / c) + self._li2(uj / c)

    # T_3(u1, u2, u3), u_i > 0 rational; 0 if u1 + u2 + u3 >= 1.
    def T3(self, u1: Fraction, u2: Fraction, u3: Fraction):
        Lq = 1 - u1 - u2 - u3
        if Lq <= 0:
            return self.acb(0)
        arb = self.arb
        L = self.q(Lq)
        U = [self.q(u1), self.q(u2), self.q(u3)]
        Aq = [Lq + u for u in (u1, u2, u3)]
        Amin_q = min(Aq)
        Amin = self.q(Amin_q)
        Lf, Af = float(Lq), float(Amin_q)
        r = Lf / Af
        pairs = ((0, 1, 2), (0, 2, 1), (1, 2, 0))
        pdata = []
        H0 = arb(0)
        Htot = arb(0)
        for j, k, l in pairs:
            Aj, Ak = L + U[j], L + U[k]
            B = 1 - U[l]                          # = u_j + u_k + L
            lam = (Aj / U[j]).log() + (Ak / U[k]).log()
            H0 += self.N_pair_at_L(U[j], U[k], L)
            # p_jk = lam + 2 >= max(lam, 2) bounds |P_m| A^m; |h_n| <= p_jk (A/B) A^-n for n >= 1.
            Htot += (lam + 2) * Amin / B
            pdata.append((1 / Aj, 1 / Ak, 1 / B, lam))
        H0f = abs(float(H0.mid()))
        Htf = float(Htot.mid()) * 1.01 + 1e-300
        scale = Lf ** self.af / self.af * (H0f / (1 - Lf) + Htf / ((1 - Af) * (1 - r)))
        N = self._nterms(r, scale)
        # Per pair: F_2^{jk}(L - s) = P(s) Q(s), P = lam - sum_{m>=1} (Aj^-m + Ak^-m) s^m / m,
        # Q = sum_q B^-(q+1) s^q; g_n = (g_{n-1} + P_n)/B; N_jk(L - v) = N_jk(L) - sum_n g_n v^{n+1}/(n+1).
        H = [arb(0)] * N
        H[0] = H0
        for iAj, iAk, iB, lam in pdata:
            pj, pk = arb(1), arb(1)
            g = lam * iB
            for n in range(N - 1):
                H[n + 1] = H[n + 1] - g / (n + 1)
                pj *= iAj
                pk *= iAk
                g = (g - (pj + pk) / (n + 1)) * iB
        # F_3(L - v) = (sum_pairs N_jk(L - v)) / (1 - v): G_n = H_0 + ... + H_n.
        G = arb(0)
        Ln = arb(1)
        re_s, im_s = arb(0), arb(0)
        for n in range(N):
            G += H[n]
            ir, ii = self.inv_n_alpha(n)
            t = G * Ln
            re_s += t * ir
            im_s += t * ii
            Ln *= L
        # |G_n| <= |H_0| + Htot A^-n/(1-A); tail <= [|H_0| L^N/(1-L) + Htot (L/A)^N/((1-A)(1-L/A))]/(N+a).
        rq = Lq / Amin_q                         # < 1 exactly
        tail = (abs(H0) * Ln / self.q(1 - Lq) + Htot * self.q(rq) ** N / self.q((1 - Amin_q) * (1 - rq))) / (N + self.a)
        return self.pow_alpha(L) * self.acb(re_s + self.unit * tail, im_s + self.unit * tail)


    # ---- T_4 by the same face expansion, recursively -------------------------------------------
    # For a subset S of the coordinates, put (all as functions of rho, expanded at rho = L in v = L - rho)
    #   F_S(rho) = int_{p >= 0, sum_{i in S} p_i = rho} prod_{i in S} 1/(u_i + p_i),
    #   N_S(rho) = int_{sum_{i in S} p_i <= rho} prod_{i in S} 1/(u_i + p_i) = int_0^rho F_S,
    # Then F_{j} = 1/(u_j + rho), N_{j} = log(1 + rho/u_j), and for |S| >= 2 (prod X_i with
    # sum X_i = U_S + rho, partial fractions, then integrate out the omitted p_i)
    #   F_S(rho) = (sum_{i in S} N_{S - i}(rho)) / (U_S + rho).
    # Every factor 1/(c - v) has c = U_S + L >= A := L + u_min, and log(1 - v/(u_j + L)) has radius
    # u_j + L >= A. Coefficient bounds |a_n| <= alpha_a (n+1)^p A^-n are carried along:
    #   1/(c - v): (1/c, 0);  N_{j}: (lam_j + 1, 0);  product with 1/(c - v): (alpha/c, p + 1);
    #   N = N(L) - int_0^v F: (|N(L)| + alpha_F A, max(p_F - 1, 0)), |N(L)| <= alpha_F L/(1 - r) if p_F = 1.
    # Here p never exceeds 1. r = L/A < 1 exactly.

    def _times_geom(self, a, ic):
        """Coefficients of a(v) / (c - v), given ic = 1/c: b_n = (a_n + b_{n-1}) / c."""
        out, prev = [], self.arb(0)
        for an in a:
            prev = (an + prev) * ic
            out.append(prev)
        return out

    @staticmethod
    def _antideriv(f, c0):
        """Coefficients of c0 - int_0^v f."""
        return [c0] + [-f[n] / (n + 1) for n in range(len(f) - 1)]

    def T4(self, *u: Fraction):
        if len(u) != 4:
            raise ValueError("T4 takes four arguments")
        Lq = 1 - sum(u)
        if Lq <= 0:
            return self.acb(0)
        arb, Q = self.arb, self.q
        Aq = Lq + min(u)
        rq = Lq / Aq                                   # < 1 exactly
        L, A, r = Q(Lq), Q(Aq), Q(rq)
        Lr = L / Q(1 - rq)                             # L / (1 - r)
        idx = range(4)
        # bounds first (independent of the truncation order K)
        lam1 = {j: Q((Lq + u[j]) / u[j]).log() for j in idx}
        b1 = {(j,): lam1[j] + 1 for j in idx}
        bF, bN = {}, {}
        for S in itertools.combinations(idx, 2):
            bF[S] = sum(b1[(i,)] for i in S) / Q(sum(u[i] for i in S) + Lq)
            bN[S] = bF[S] * Lr + bF[S] * A
        for S in itertools.combinations(idx, 3):
            bF[S] = sum(bN[P] for P in itertools.combinations(S, 2)) / Q(sum(u[i] for i in S) + Lq)
            bN[S] = bF[S] * Lr + bF[S] * A
        b4 = sum(bN[S] for S in itertools.combinations(idx, 3))   # U + L = 1
        # tail of sum_n f_n L^n/(n+alpha), n >= K: |f_n| <= b4 (n+1) A^-n, (n+1)/(n+a) <= m_K
        Lf, rf = float(Lq), float(rq)
        scale = float(b4.mid()) * 1.01 * max(1.0, 2.0 / self.af) * Lf ** self.af / (1 - rf)
        K = self._nterms(rf, scale)
        mK = arb(1) if self.a_q >= 1 else (K + 1) / (K + self.a)
        Ln = [arb(1)]
        for _ in range(K):
            Ln.append(Ln[-1] * L)
        # N_{j}(L - v) = lam_j + log(1 - v/A_j)
        N = {}
        for j in idx:
            iAj = 1 / Q(Lq + u[j])
            c, pw = [lam1[j]], arb(1)
            for n in range(1, K):
                pw *= iAj
                c.append(-pw / n)
            N[(j,)] = c
        for d in (2, 3):
            for S in itertools.combinations(idx, d):
                tot = [arb(0)] * K
                for P in itertools.combinations(S, d - 1):
                    tot = [x + y for x, y in zip(tot, N[P])]
                F = self._times_geom(tot, 1 / Q(sum(u[i] for i in S) + Lq))
                # N_S(L) = int_0^L F_S = sum_n F_n L^{n+1}/(n+1) + tail, tail <= bF L r^K / (1 - r)
                c0 = sum(F[n] * Ln[n + 1] / (n + 1) for n in range(K))
                c0 += self.unit * (bF[S] * Lr * r ** K)
                if d == 2:  # closed form as an internal consistency check
                    j, k = S
                    cf = self.N_pair_at_L(Q(u[j]), Q(u[k]), L)
                    if not cf.overlaps(c0):
                        raise ArithmeticError(f"N_pair(L) closed form vs series mismatch at {u}")
                    c0 = cf
                N[S] = self._antideriv(F, c0)
        tot = [arb(0)] * K
        for S in itertools.combinations(idx, 3):
            tot = [x + y for x, y in zip(tot, N[S])]
        F4 = self._times_geom(tot, arb(1))
        re_s, im_s = arb(0), arb(0)
        for n in range(K):
            ir, ii = self.inv_n_alpha(n)
            t = F4[n] * Ln[n]
            re_s += t * ir
            im_s += t * ii
        tail = b4 * mK * r ** K / Q(1 - rq)
        return self.pow_alpha(L) * self.acb(re_s + self.unit * tail, im_s + self.unit * tail)

    # ---- second, slower rigorous method (Arb's validated Gauss-Legendre integration) ----
    def F_closed(self, U, r, an: bool, deriv: bool):
        """F_d(r) (d = len(U) in {2, 3}) or F_d'(r) from the closed forms, on a complex ball r.
        With an=True returns a non-finite ball unless every log/Li2 is certainly off its cut."""
        acb = self.acb
        nan = acb("nan")

        def lg(z):
            return acb(z).log(analytic=an)

        def li2(z):
            z = acb(z)
            if an and not (z.real < 1):
                return None
            return z.polylog(2)

        Us = sum(U)
        if len(U) == 2:
            lsum = lg(1 + r / U[0]) + lg(1 + r / U[1])
            if not deriv:
                return lsum / (Us + r)
            return (1 / (U[0] + r) + 1 / (U[1] + r)) / (Us + r) - lsum / (Us + r) ** 2
        S, dS = acb(0), acb(0)
        for j, k in ((0, 1), (0, 2), (1, 2)):
            uj, uk = U[j], U[k]
            c = uj + uk + r
            l1, l2 = li2((uj + r) / c), li2(uj / c)
            if l1 is None or l2 is None:
                return nan
            S += lg(c / uk) * lg(1 + r / uj) - l1 + l2
            dS += (lg(1 + r / uj) + lg(1 + r / uk)) / c
        if not deriv:
            return S / (Us + r)
        return dS / (Us + r) - S / (Us + r) ** 2

    def T_integral(self, us):
        """T_d(u), d = len(us) in {2, 3}, as
        int_0^{L/2} (L-r)^(alpha-1) F_d + (L/2)^alpha F_d(L/2)/alpha + (1/alpha) int_{L/2}^L (L-r)^alpha F_d'
        (integration by parts on the half with the endpoint singularity). On balls touching r = L the
        integrand is bounded by |w|^a e^{|b| pi} |F_d'| (any branch value of w^alpha), w = L - r."""
        acb, arb = self.acb, self.arb
        Lq = 1 - sum(us)
        if Lq <= 0:
            return acb(0)
        al, a = self.alpha, self.a
        U = [self.q(u) for u in us]
        L = self.q(Lq)
        ebpi = (abs(self.b) * arb.pi()).exp()

        def f1(r, an):
            return ((al - 1) * acb(L - r).log(analytic=an)).exp() * self.F_closed(U, r, an, False)

        def f2(r, an):
            w = acb(L - r)
            if not an and w.abs_lower() == 0:
                R = w.abs_upper()
                X = (a * R.log()).exp() * ebpi if R > 0 else arb(0)
                return acb(self.unit * X, self.unit * X) * self.F_closed(U, r, an, True)
            return (al * w.log(analytic=an)).exp() * self.F_closed(U, r, an, True)

        tol = arb(2) ** -(self.prec - 16)
        h = L / 2
        I1 = acb.integral(f1, 0, h, abs_tol=tol)
        I2 = acb.integral(f2, h, L, abs_tol=tol)
        return I1 + self.pow_alpha(h) * self.F_closed(U, acb(h), False, False) / al + I2 / al


def evaluate(cert: dict, prec: int, log=None, t3: str = "series") -> dict:
    ev = Evaluator(cert["alpha"], prec)
    T3 = ev.T3 if t3 == "series" else (lambda *u: ev.T_integral(u))
    acb = ev.acb
    t = cert["breaks"]
    M = len(t) - 1
    tau = t[0]
    a_q, b_q = cert["alpha"]
    s = acb(1) - ev.alpha
    eta = [ev.cq(er, ei) for er, ei in cert["eta"]]
    dq = [(cert["eta"][0][0] - (1 - a_q), cert["eta"][0][1] + b_q)]
    for j in range(1, M):
        dq.append((cert["eta"][j][0] - cert["eta"][j - 1][0], cert["eta"][j][1] - cert["eta"][j - 1][1]))
    delta = [ev.cq(x, y) for x, y in dq]

    t0 = time.perf_counter()
    K = ev.I_alpha(tau)
    D = ev.I_alpha(tau, real_a=True)
    lin = acb(0)
    for j in range(M):
        lin += delta[j] * (ev.I_alpha(1 - t[j]) - K)
    # ordered-pair sum via symmetry of T_2
    S2 = acb(0)
    n2 = 0
    for i in range(M):
        if 2 * t[i] >= 1:
            break
        S2 += delta[i] * delta[i] * ev.T2(t[i], t[i])
        n2 += 1
        acc = acb(0)
        for j in range(i + 1, M):
            if t[i] + t[j] >= 1:
                break
            acc += delta[j] * ev.T2(t[i], t[j])
            n2 += 1
        S2 += 2 * delta[i] * acc
    t1 = time.perf_counter()
    # ordered-triple sum via symmetry of T_3 (multiplicity 1, 3 or 6)
    S3 = acb(0)
    n3 = 0
    for i in range(M):
        if 3 * t[i] >= 1:
            break
        for j in range(i, M):
            if t[i] + 2 * t[j] >= 1:
                break
            for k in range(j, M):
                if t[i] + t[j] + t[k] >= 1:
                    break
                mult = 6 if (i < j < k) else (1 if i == k else 3)
                S3 += mult * delta[i] * delta[j] * delta[k] * T3(t[i], t[j], t[k])
                n3 += 1
    t2 = time.perf_counter()
    # ordered-quadruple sum via symmetry of T_4 (multiplicity 4!/prod(repeat counts)!); only when d* = 4
    S4 = acb(0)
    n4 = 0
    if dstar(tau) >= 4:
        for i in range(M):
            if 4 * t[i] >= 1:
                break
            for j in range(i, M):
                if t[i] + 3 * t[j] >= 1:
                    break
                for k in range(j, M):
                    if t[i] + t[j] + 2 * t[k] >= 1:
                        break
                    for l in range(k, M):
                        if t[i] + t[j] + t[k] + t[l] >= 1:
                            break
                        mult = 24
                        for c in collections.Counter((i, j, k, l)).values():
                            mult //= math.factorial(c)
                        S4 += mult * delta[i] * delta[j] * delta[k] * delta[l] * ev.T4(t[i], t[j], t[k], t[l])
                        n4 += 1
    t4 = time.perf_counter()
    Y = 1 + s * K - lin + S2 / 2 - S3 / 6 + S4 / 24
    Cb = ev.q(cert["C"])
    absY2 = Y.real * Y.real + Y.imag * Y.imag
    gap = Cb * Cb * D * D - absY2
    if log:
        log(f"# prec={prec} M={M} tau={tau} T2 pairs={n2} T3 triples={n3} T4 quadruples={n4} d*={dstar(tau)} "
            f"time T2={t1 - t0:.1f}s T3({t3})={t2 - t1:.1f}s T4={t4 - t2:.1f}s")
    return {"D": D, "Y": Y, "gap": gap, "ratio": absY2.sqrt() / D, "S2": S2, "S3": S3, "S4": S4, "K": K,
            "n2": n2, "n3": n3, "n4": n4}


def round_up(x: Fraction, digits: int = 20) -> str:
    scale = 10 ** digits
    n = -((-x.numerator * scale) // x.denominator)  # ceil(x * scale)
    sign = "-" if n < 0 else ""
    n = abs(n)
    ip, fp = divmod(n, scale)
    fs = str(fp).rjust(digits, "0").rstrip("0")
    return f"{sign}{ip}" + (f".{fs}" if fs else "")


def main(argv=None) -> int:
    if sys.flags.optimize:
        print("refusing to run under python -O", file=sys.stderr)
        return EXIT_USAGE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("certificate")
    ap.add_argument("--beat", metavar="V", help="print BEATS only if the certified C < V exactly")
    ap.add_argument("--prec", type=int, default=160, help="Arb working precision in bits (default 160)")
    ap.add_argument("--t3", choices=("series", "integral"), default="series",
                    help="T_3 method: face series (default) or Arb rigorous integration (slower)")
    ap.add_argument("--format", type=int, choices=(1, 2), default=2,
                    help="precondition set: 2 = FORMAT2.md, 1/5 < tau < 1/2 (default); 1 = FORMAT.md, 1/4 < tau < 1/2")
    ap.add_argument("--quiet", action="store_true")
    try:
        args = ap.parse_args(argv)
    except SystemExit as e:
        return EXIT_USAGE if e.code else EXIT_VALID
    beat = None
    if args.beat is not None:
        try:
            beat = parse_rational(args.beat, "--beat")
        except CertificateError as e:
            print(f"bad --beat: {e}", file=sys.stderr)
            return EXIT_USAGE
    if not (64 <= args.prec <= 4096):
        print("--prec must be in [64, 4096]", file=sys.stderr)
        return EXIT_USAGE
    try:
        import flint  # noqa: F401
    except ImportError:
        print("python-flint is required", file=sys.stderr)
        return EXIT_USAGE

    log = None if args.quiet else (lambda m: print(m, file=sys.stderr, flush=True))
    try:
        cert = load_certificate(args.certificate)
        check_preconditions(cert, args.format)
    except CertificateError as e:
        print(f"INVALID {e}")
        return EXIT_INVALID

    prec = args.prec
    for attempt in range(2):
        res = evaluate(cert, prec, log, args.t3)
        if res["gap"] > 0:
            break
        if res["gap"] < 0:
            print(f"INVALID |Y| >= C*D: |Y|/D in {res['ratio'].str(20, radius=False)} vs C = {cert['C']}")
            return EXIT_INVALID
        prec *= 2  # ball straddles 0: retry once at double precision
    else:
        print(f"INVALID C^2 D^2 - |Y|^2 > 0 not proved: gap in {res['gap'].str(10)}")
        return EXIT_INVALID
    if log:
        log(f"# D       in {res['D'].str(30)}")
        log(f"# |Y|/D   in {res['ratio'].str(30)}")
        log(f"# C^2D^2-|Y|^2 in {res['gap'].str(12)}")
        log(f"# S3 = sum d_i d_j d_k T3 in {res['S3'].str(15)}")
        log(f"# S4 = sum d_i d_j d_k d_l T4 in {res['S4'].str(15)}")
    print(f"VALID {round_up(cert['C'])}")
    if beat is not None:
        if cert["C"] < beat:
            print(f"BEATS {args.beat}")
            return EXIT_VALID
        print(f"DOES_NOT_BEAT {args.beat}")
        return EXIT_NOT_BEAT
    return EXIT_VALID


if __name__ == "__main__":
    sys.exit(main())
