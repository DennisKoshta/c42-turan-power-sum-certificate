#!/usr/bin/env python3
"""Rigorous Arb checker for the general-tau C42 criterion with 1/5 < tau < 1/2 (d* <= 4), NOTES.md sections 4, 5, 7.

Same certificate format as FORMAT.md, except that tau may be in (1/5, 1/2) (FORMAT2.md).
    Y = 1 + sK - sum_j d_j (I(1-t_j) - K) + 1/2 sum d d T2 - 1/6 sum d d d T3 + 1/24 sum d d d d T4   (ordered tuples)
T2: log series (check_tau.J_series).  T_d with q = L/min u <= Q0, and every T4 (q < 1 as tau > 1/5): Dirichlet series (NOTES 7.1).
T3 with q > Q0: 1-D integral (NOTES 7.2), x = l - W e^{-v}, acb.integral on [0, V] plus an explicit tail at V.
PASS iff the exact preconditions hold and C^2 D^2 - |Y|^2 > 0 in ball arithmetic.  Exit 0 PASS, 1 FAIL.
Refuses python -O (no reliance on assert anyway).
"""
import sys
if sys.flags.optimize:
    sys.exit("refusing to run under python -O")
import argparse, hashlib, json, math, os
from fractions import Fraction
from itertools import combinations_with_replacement
from collections import Counter
from flint import arb, acb, fmpq, fmpq_poly, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_tau import Q, ball, err, I_series, J_series

Q0 = Fraction(4, 5)


def Td_series(alpha, a, us, eps):
    """NOTES 7.1: T_d = L^(alpha+d-1) sum_N (-L)^N h_N/(alpha)_{N+d}, h_N = sum_{|k|=N} prod k_i!/u_i^(k_i+1) (exact),
    tail sum_{N>=N0} <= B_N0/(1-q), B_N = L^(a+d-1) q^N C(N+d-1,d-1) N!/(prod u_i |(alpha)_{N+d}|), q = L/min u < 1."""
    d = len(us)
    Lq = 1 - sum(us)
    if Lq <= 0:
        return acb(0)
    umin = min(us)
    if not Lq < umin:
        raise ArithmeticError("Dirichlet series needs L < min u")
    N0 = max(20, int(math.log(eps) / math.log(float(Lq / umin))) + 10)
    prod = None
    for u in us:
        uu = fmpq(u.numerator, u.denominator)
        row, f, pw = [], fmpq(1), 1 / uu
        for k in range(N0):
            if k:
                f *= k
                pw /= uu
            row.append(f * pw)
        if prod is None:
            prod = fmpq_poly(row)
        else:
            full = prod * fmpq_poly(row)
            prod = fmpq_poly([full[n] for n in range(N0)])
    h = [prod[n] for n in range(N0)]
    L = ball(Lq)
    poch = acb(1)
    for i in range(d):
        poch *= alpha + i
    s, mL = acb(0), arb(1)
    for N in range(N0):
        s += mL * arb(h[N]) / poch
        mL *= -L
        poch *= alpha + N + d          # now (alpha)_{N+1+d}
    q = L / ball(umin)
    pu = arb(1)
    for u in us:
        pu *= ball(u)
    BN0 = L ** (a + d - 1) * q ** N0 * arb(math.comb(N0 + d - 1, d - 1)) * arb(math.factorial(N0)) / (pu * abs(poch))
    return acb(L) ** (alpha + d - 1) * s + err(BN0 / (1 - q))


def _phi(x, alpha, a, u1, u2, K, analytic):
    """Phi(x) = [sum_k c_k lam^k/(alpha+k)] / (x(1-x)), NOTES 7.2, with an explicit truncation bound; x an acb ball.
    Returns None if the convergence conditions cannot be certified on the ball."""
    t, r = u1 / (1 - x), u2 / (1 - x)
    lam, b, d = 1 - t - r, 1 - t, 1 - r
    c0 = (b * d / (t * r)).log(analytic=analytic)
    if not c0.is_finite():
        return None
    al, ab, ad = abs(lam), abs(b), abs(d)
    if not (al < 1 and ab > 0 and ad > 0):
        return None
    # rho = max(1, 1/|b|) as an upper bound over the ball (|b| > 0 certified above); review R1: y < 1 must be certain
    # for the whole ball, independently of `analytic` (1/(1-y) is finite but meaningless for y > 1)
    rb = arb((1 / ab).union(arb(1)).upper())
    rd = arb((1 / ad).union(arb(1)).upper())
    yb, yd = al * rb, al * rd
    if not (yb < 1 and yd < 1):
        return None
    c, s, p, ib, id_ = c0, c0 / alpha, acb(1), acb(1), acb(1)
    for k in range(1, K):              # terms 0..K-1 summed; the tail below bounds k >= K (review R5)
        ib /= b
        id_ /= d
        c -= (ib + id_) / k
        p *= lam
        s += c * p / (alpha + k)
    tail = abs(c0) * al ** K / (1 - al)
    for y in (yb, yd):
        tail += y ** K * (K / (1 - y) + y / (1 - y) ** 2)
    tail /= a + K
    return (s + err(tail)) / (x * (1 - x))


def T3_integral(alpha, a, us, eps):
    """NOTES 7.2: T3 = W^(alpha+1) int_0^inf e^{-(alpha+1)v} Phi(l - W e^{-v}) dv, W = L, l = 1 - u1 - u2.
    T3 is symmetric; the smallest u is the outer variable u3 (smaller series ratios on complex balls, review)."""
    u1q, u2q, u3q = sorted(us, reverse=True)
    lq = 1 - u1q - u2q
    Wq = lq - u3q
    if Wq <= 0:
        return acb(0)
    u1, u2, l, W = ball(u1q), ball(u2q), ball(lq), ball(Wq)
    # series length from the real path maximum of lam/b at v = 0: L/(L+u2) (and L/(L+u1)), plus slack
    y0 = float(max(Wq / (Wq + u2q), Wq / (Wq + u1q)))
    K = int(math.log(eps) / math.log(min(0.999, y0 + (1 - y0) / 3))) + 10

    def f(v, analytic):
        x = l - W * (-v).exp()
        ph = _phi(x, alpha, a, u1, u2, K, analytic)
        if ph is None:
            return acb(float("nan"))
        return (-(alpha + 1) * v).exp() * ph

    V = int(math.log(1 / eps) / (float(a.lower()) + 1)) + 8      # endpoint tail ~ e^{-(a+1)V}; bounded rigorously below
    val = acb.integral(f, 0, V, abs_tol=arb(eps) / 1000)
    if not val.is_finite():
        raise ArithmeticError("T3 integral not finite")
    # tail: x in [l - W e^{-V}, l]
    lo = l - W * arb(-V).exp()
    xb = acb(lo.union(l))
    ph = _phi(xb, alpha, a, u1, u2, K, False)
    if ph is None:
        raise ArithmeticError("T3 tail ball: series convergence not certified")
    tail = arb(-(a + 1) * V).exp() / (a + 1) * abs(ph)
    return acb(W) ** (alpha + 1) * (val + err(tail))


def check(data, terms, eps):
    if data.get("constant") != "C_42":
        raise ValueError("not a C_42 certificate")
    C = Q(data["C"])
    ar, ai = map(Q, data["alpha"])
    ts = [Q(x) for x in data["breaks"]]
    eta = [(Q(x), Q(y)) for x, y in data["eta"]]
    tau = ts[0]
    if not (Fraction(1, 5) < tau < Fraction(1, 2)):
        raise ValueError("need 1/5 < tau < 1/2")
    dstar = -((-1) // tau) - 1          # ceil(1/tau) - 1, exact (review R2: never drop d-fold terms silently)
    if not 2 <= dstar <= 4:
        raise ValueError(f"d* = {dstar} not supported")
    if ts[-1] != 1 - tau or any(not x < y for x, y in zip(ts, ts[1:])) or len(eta) != len(ts) - 1:
        raise ValueError("breaks must increase from tau to 1-tau, one eta per step")
    if not (0 < ar and 0 < C < 1):
        raise ValueError("need Re alpha > 0, 0 < C < 1")
    sr, si = 1 - ar, -ai
    margins = [C * C - sr * sr - si * si] + [C * C - x * x - y * y for x, y in eta]
    if not min(margins) > 0:
        raise ArithmeticError("profile radius exceeds C")
    a = ball(ar)
    alpha = acb(ball(ar), ball(ai))
    s = 1 - alpha
    K = I_series(alpha, a, tau, terms)
    Dc = I_series(acb(a), a, tau, terms)
    if not Dc.imag.contains(0):
        raise ArithmeticError("D not real")
    D = Dc.real
    if not D > 0:
        raise ArithmeticError("D > 0 not certified")
    E = [acb(ball(x), ball(y)) for x, y in eta]
    delta = [E[0] - s] + [E[j] - E[j - 1] for j in range(1, len(E))]
    M = len(E)
    t = ts[:-1]
    Y = 1 + s * K
    for j in range(M):
        Y -= delta[j] * (I_series(alpha, a, 1 - t[j], terms) - K)
    counts = {"J": 0, "T3_series": 0, "T3_integral": 0, "T4_series": 0}
    for i in range(M):
        for j in range(i + 1):
            if t[i] + t[j] < 1:
                w = arb(1) if i != j else arb(fmpq(1, 2))
                Y += w * delta[i] * delta[j] * J_series(alpha, a, t[i], t[j], terms)
                counts["J"] += 1
    for d in range(3, dstar + 1):
        sign = -1 if d == 3 else 1
        for tup in combinations_with_replacement(range(M), d):
            us = tuple(t[i] for i in tup)
            Lq = 1 - sum(us)
            if Lq <= 0:
                continue
            mult = math.factorial(d)
            for c in Counter(tup).values():
                mult //= math.factorial(c)
            if Lq / min(us) <= Q0 or d == 4:      # d = 4: q <= (1 - 4 tau)/tau < 1; Td_series checks L < min u exactly
                T = Td_series(alpha, a, us, eps)
                counts[f"T{d}_series"] += 1
            elif d == 3:
                T = T3_integral(alpha, a, us, eps)
                counts["T3_integral"] += 1
            pr = acb(1)
            for i in tup:
                pr *= delta[i]
            Y += arb(fmpq(sign * mult, math.factorial(d))) * pr * T
    # no (dstar+1)-fold terms: (dstar+1) tau >= 1, so every such tuple has L <= 0
    Cb = ball(C)
    margin = Cb * Cb * D * D - (Y.real ** 2 + Y.imag ** 2)
    if not margin > 0:
        raise ArithmeticError(f"|Y| < C D not certified: margin {margin}")
    return {"verdict": "PASS", "C": str(C), "tau": str(tau), "M": M, "D": D.str(30), "Y": Y.str(30),
            "ratio": (abs(Y) / D).str(30), "analytic_margin": margin.str(10),
            "min_exact_profile_margin": str(min(margins)), **counts}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("certificate")
    p.add_argument("--bits", type=int, default=256)
    p.add_argument("--terms", type=int, default=600)
    p.add_argument("--eps", type=float, default=1e-40)
    a = p.parse_args()
    ctx.prec = a.bits
    raw = open(a.certificate, "rb").read()
    try:
        res = check(json.loads(raw), a.terms, a.eps)
        code = 0
    except (ValueError, ArithmeticError, KeyError) as e:
        res, code = {"verdict": "FAIL", "reason": str(e)}, 1
    res["sha256"] = hashlib.sha256(raw).hexdigest()
    print(json.dumps(res, indent=1))
    return code


if __name__ == "__main__":
    sys.exit(main())
