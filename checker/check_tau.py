#!/usr/bin/env python3
"""Rigorous Arb checker for the general-tau C42 criterion (NOTES.md sections 3-4), 1/4 < tau < 1/2.

Certificate JSON (exact rational strings): {"constant": "C_42", "C": ..., "alpha": [re, im],
"breaks": [t_0 = tau, ..., t_M = 1 - tau], "eta": [[re, im] x M]}.
Checks, with outward-rounded Arb balls:
    |s| < C, |eta_j| < C          (exact rationals)
    |Y| < C * D                   (Arb), D = I_a(tau),
    Y = 1 + sK - sum_j d_j (I(1-t_j) - K) + 1/2 sum_{i,j} d_i d_j J(t_i,t_j) - 1/6 sum_{i,j,k} d_i d_j d_k T3(t_i,t_j,t_k).
Exit 0 = PASS, 1 = FAIL, 3 = usage error.  Refuses python -O (no reliance on assert anyway).
"""
import sys
if sys.flags.optimize:
    sys.exit("refusing to run under python -O")
import argparse, hashlib, json, math
from fractions import Fraction
from flint import arb, acb, fmpq, fmpq_poly, ctx


def Q(s):
    if not isinstance(s, (str, int)):
        raise ValueError("exact rational strings required")
    return Fraction(s)


def ball(x):
    return arb(fmpq(x.numerator, x.denominator))


def err(r):
    """acb ball [-r, r] + i[-r, r] for an arb upper bound r >= 0."""
    e = arb(0, r.upper())
    return acb(e, e)


def I_series(alpha, a, x, terms):
    """int_0^x v^(alpha-1)/(1-v) dv = sum_k x^(alpha+k)/(alpha+k), 0 < x < 1 rational."""
    X = ball(x)
    s, p = acb(0), arb(1)
    for k in range(terms):
        s += p / (alpha + k)
        p *= X
    tail = X ** (a + terms) / ((a + terms) * (1 - X))   # sum_{k>=terms} x^(a+k)/(a+k)
    return acb(X) ** alpha * s + err(tail)


def J_series(alpha, a, t, r, terms):
    """J(t,r) = int_0^L v^(alpha-1)/(1-v) log((b-v)(d-v)/(tr)) dv, L=1-t-r, b=1-t, d=1-r.
    1/(1-v) * log(...) = sum_k c_k v^k, c_k = c0 - sum_{q<=k} (b^-q + d^-q)/q, so J = L^alpha sum c_k L^k/(alpha+k).
    |c_k| <= c0 + sum_{q<=k}(b^-q+d^-q)/q <= c0 + k (b^-k + d^-k)  (b,d < 1).  Tail k >= terms:
    sum L^(a+k)/(a+k) |c_k| <= L^a/(a+terms) [c0 L^terms/(1-L) + sum_{y in L/b, L/d} sum_{k>=terms} k y^k]."""
    Lq = 1 - t - r
    if Lq <= 0:
        return acb(0)
    L, b, d = ball(Lq), ball(1 - t), ball(1 - r)
    c0 = (b * d / (ball(t) * ball(r))).log()
    if not c0 > 0:
        raise ArithmeticError("c0 > 0 not certified")
    c, s, p, ib, id_ = c0, acb(0), arb(1), arb(1), arb(1)
    for k in range(terms):
        if k:
            ib /= b
            id_ /= d
            c -= (ib + id_) / k
            p *= L
        s += c * p / (alpha + k)
    tail = c0 * L ** terms / (1 - L)
    for y in (L / b, L / d):
        if not y < 1:
            raise ArithmeticError("J series ratio not < 1")
        tail += y ** terms * (terms / (1 - y) + y / (1 - y) ** 2)
    tail *= L ** a / (a + terms)
    return acb(L) ** alpha * s + err(tail)


def T3_series(alpha, a, ts, prec_terms_eps):
    """T3 = L^(alpha+2) sum_N (-L)^N h_N / (alpha)_{N+3},  h_N = sum_{|k|=N} prod k_i!/t_i^(k_i+1)  (exact rationals),
    valid for L = 1 - sum t < min t.  Tail: B_N = L^(a+2) q^N C(N+2,2) N! / (t1 t2 t3 |(alpha)_{N+3}|), q = L/min t,
    B_{N+1} <= q B_N, so sum_{N>=N0} <= B_{N0}/(1-q)."""
    Lq = 1 - sum(ts)
    if Lq <= 0:
        return acb(0)
    tmin = min(ts)
    if not Lq < tmin:
        raise ArithmeticError("T3 Dirichlet series needs L < min t")
    qf = Lq / tmin
    # number of terms so that the geometric bound is below eps (float only chooses N0; the bound is rigorous)
    N0 = max(20, int(math.log(prec_terms_eps) / math.log(float(qf))) + 10)
    L = ball(Lq)
    # a_i[k] = k!/t_i^(k+1) as exact fmpq
    seqs = []
    for t in ts:
        tt = fmpq(t.numerator, t.denominator)
        row, f, pw = [], fmpq(1), 1 / tt
        for k in range(N0):
            if k:
                f *= k
                pw /= tt
            row.append(f * pw)
        seqs.append(row)
    # h = a1 * a2 * a3 (truncated product of exact rational polynomials, FLINT fmpq_poly)
    prod = fmpq_poly(seqs[0]) * fmpq_poly(seqs[1])
    prod = fmpq_poly([prod[n] for n in range(N0)]) * fmpq_poly(seqs[2])
    h = [prod[n] for n in range(N0)]
    s, poch, mL = acb(0), alpha * (alpha + 1) * (alpha + 2), arb(1)
    for N in range(N0):
        s += mL * arb(h[N]) / poch
        mL *= -L
        poch *= alpha + N + 3
    # tail bound at N0: poch now = (alpha)_{N0+3}
    q = L / ball(tmin)
    BN0 = L ** (a + 2) * q ** N0 * arb(math.comb(N0 + 2, 2)) * arb(math.factorial(N0)) / (
        ball(ts[0] * ts[1] * ts[2]) * abs(poch))
    tail = BN0 / (1 - q)
    return acb(L) ** (alpha + 2) * s + err(tail)


def check(data, terms, eps):
    if data.get("constant") != "C_42":
        raise ValueError("not a C_42 certificate")
    C = Q(data["C"])
    ar, ai = map(Q, data["alpha"])
    ts = [Q(x) for x in data["breaks"]]
    eta = [(Q(x), Q(y)) for x, y in data["eta"]]
    tau = ts[0]
    if not (Fraction(1, 4) < tau < Fraction(1, 2)):
        raise ValueError("need 1/4 < tau < 1/2")
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
    nJ = nT = 0
    for i in range(M):
        for j in range(i + 1):
            if t[i] + t[j] < 1:
                w = arb(1) if i != j else arb(fmpq(1, 2))
                Y += w * delta[i] * delta[j] * J_series(alpha, a, t[i], t[j], terms)
                nJ += 1
    for i in range(M):
        for j in range(i + 1):
            for k in range(j + 1):
                if t[i] + t[j] + t[k] < 1:
                    # ordered permutations of the unordered triple, times the 1/3! of the d = 3 term
                    mult = {3: 6, 2: 3, 1: 1}[len({i, j, k})]
                    Y -= arb(fmpq(mult, 6)) * delta[i] * delta[j] * delta[k] * T3_series(alpha, a, (t[i], t[j], t[k]), eps)
                    nT += 1
    if 4 * tau <= 1:
        raise ValueError("d = 4 terms not implemented")
    Cb = ball(C)
    margin = Cb * Cb * D * D - (Y.real ** 2 + Y.imag ** 2)
    if not margin > 0:
        raise ArithmeticError(f"|Y| < C D not certified: margin {margin}")
    return {"verdict": "PASS", "C": str(C), "tau": str(tau), "M": M, "D": D.str(30), "Y": Y.str(30),
            "ratio": (abs(Y) / D).str(30), "analytic_margin": margin.str(10),
            "min_exact_profile_margin": str(min(margins)), "J_pairs": nJ, "T3_triples": nT}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("certificate")
    p.add_argument("--bits", type=int, default=384)
    p.add_argument("--terms", type=int, default=400)
    p.add_argument("--eps", type=float, default=1e-45)
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
