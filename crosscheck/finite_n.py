"""Finite-n model of the Griego/Roehrig/Pinho coefficient-cancellation construction for C42.

Exact reduction (any n): for S_1..S_n in C put
    B(z) = (1-z)^{-1} exp(-sum_m S_m z^m / m) = sum_l b_l z^l,  k b_k = sum_{m<=k} (1-S_m) b_{k-m}.
If b_n = 0, the roots of Z^{n-1} + b_1 Z^{n-2} + ... + b_{n-1} together with y_1 = 1 have power sums
S_1..S_n, so R_n <= max_k |S_k| (divide by max|y_j| >= 1).  Conversely any z with z_1 = 1 gives b_n = 0,
so R_n = min{max|S_k| : b_n(S) = 0} exactly.

Construction (cut-off tau, M middle steps, alpha = a+ib, s = 1-alpha):
    A = floor(tau n);  S_m = s (m <= A);  S_m = eta_j on the middle (A < m < n-A), where step j
    covers floor(n t_j) < m <= floor(n t_{j+1}), t_j = tau + j(1-2tau)/M;  the final block
    F = {n-A..n} is chosen to force b_n = 0, which is possible with max_F |S_m| = |Y_n| / D_n,
    Y_n = n b_n(S_F = 0),  D_n = sum_F n|beta_{n-m}|/m,  beta_l = (alpha)_l / l!.
(A final index plus any other correction index is >= n+1, so b_n is affine in S_F for every tau.)
For tau < 1/3 the middle block may contribute 3 (or more) corrections; the exact series below keeps
them all, so the finite-n value is exact for every tau.

value(params) = max(|s|, max_j |eta_j|, |Y_n|/D_n) is a valid upper bound on R_n (float arithmetic).
"""
import numpy as np
from fractions import Fraction
from scipy.optimize import minimize


def _conv(a, b, L):
    """First L coefficients of a*b (complex), by FFT."""
    N = 1
    while N < len(a) + len(b):
        N *= 2
    r = np.fft.ifft(np.fft.fft(a, N) * np.fft.fft(b, N))[:L]
    return r


class Model:
    def __init__(self, n, tau, M):
        self.n, self.M = n, M
        tau = Fraction(tau)
        self.tau = tau
        A = (tau.numerator * n) // tau.denominator
        self.A = A
        cuts = [tau + Fraction(j, M) * (1 - 2 * tau) for j in range(M + 1)]
        fl = [(c.numerator * n) // c.denominator for c in cuts]
        assert fl[0] == A
        # middle indices A+1..n-A-1 ; step of m = max j with m > fl[j]
        mid = np.arange(A + 1, n - A)
        step = np.searchsorted(np.array(fl[:-1]), mid, side="left") - 1
        assert step.min() == 0 and step.max() <= M - 1
        self.mid, self.step = mid, step
        self.F = np.arange(n - A, n + 1)
        self.K = n // (A + 1)  # max number of corrections reaching index n

    def beta(self, alpha):
        n = self.n
        l = np.arange(1, n + 1)
        r = (alpha + l - 1) / l
        b = np.concatenate([[1.0 + 0j], np.cumprod(r)])
        # h_l = sum_{k<l} 1/(alpha+k) = d log beta_l / d alpha
        h = np.concatenate([[0j], np.cumsum(1.0 / (alpha + l - 1))])
        return b, h

    def evaluate(self, alpha, eta, grad=False):
        """Return (|s|, |eta|_max, ratio) and optionally gradients of ratio."""
        n, A = self.n, self.A
        s = 1 - alpha
        beta, h = self.beta(alpha)
        c = np.zeros(n + 1, complex)
        c[self.mid] = (s - eta[self.step]) / self.mid
        c[self.F] = s / self.F
        # e = exp(C) truncated at degree n; at most K factors reach degree n
        e = np.zeros(n + 1, complex)
        e[0] = 1.0
        term = np.zeros(n + 1, complex)
        term[0] = 1.0
        for k in range(1, self.K + 1):
            term = _conv(term, c, n + 1) / k
            e += term
        Y = n * np.dot(beta[::-1], e)
        absb = np.abs(beta)
        Fm = self.F
        D = n * np.sum(absb[n - Fm] / Fm)
        ratio = abs(Y) / D
        if not grad:
            return Y, D, ratio
        g = _conv(beta, e, n + 1)  # g_k = sum_i beta_{k-i} e_i ; dY/dc_m = n g_{n-m}
        dYdc = n * g[n - np.arange(n + 1)]
        # eta_j (holomorphic derivative)
        w = -dYdc[self.mid] / self.mid
        dYdeta = np.bincount(self.step, weights=w.real, minlength=self.M) + 1j * np.bincount(
            self.step, weights=w.imag, minlength=self.M)
        idx = np.arange(A + 1, n + 1)
        dYds = np.sum(dYdc[idx] / idx)
        dYdalpha = n * np.dot((beta * h)[::-1], e) - dYds
        # D derivatives wrt Re/Im alpha
        wD = n * absb[n - Fm] / Fm
        hF = h[n - Fm]
        dD_dar = np.sum(wD * hF.real)
        dD_dai = -np.sum(wD * hF.imag)
        u = np.conj(Y) / abs(Y)
        # d ratio / d(Re z), d(Im z) for holomorphic dependence Y_z: Re(u Y_z)/D, Re(u i Y_z)/D
        dr_dar = (u * dYdalpha).real / D - ratio / D * dD_dar
        dr_dai = (u * 1j * dYdalpha).real / D - ratio / D * dD_dai
        dr_deta_r = (u * dYdeta).real / D
        dr_deta_i = (u * 1j * dYdeta).real / D
        return Y, D, ratio, (dr_dar, dr_dai, dr_deta_r, dr_deta_i)

    def value(self, alpha, eta):
        Y, D, ratio = self.evaluate(alpha, eta)
        return max(abs(1 - alpha), np.abs(eta).max(), ratio), ratio

    def optimize(self, alpha0, eta0, maxiter=500, ftol=1e-15):
        M = self.M
        t0 = self.value(alpha0, eta0)[0]
        x0 = np.concatenate([[t0, alpha0.real, alpha0.imag], eta0.real, eta0.imag])
        cache = {}

        def ev(x):
            k = x.tobytes()
            if k not in cache:
                cache.clear()
                al = x[1] + 1j * x[2]
                et = x[3:3 + M] + 1j * x[3 + M:]
                cache[k] = self.evaluate(al, et, grad=True)
            return cache[k]

        def cons(x):
            Y, D, r, _ = ev(x)
            t = x[0]
            sr, si = 1 - x[1], -x[2]
            return np.concatenate([[t - r, t * t - sr * sr - si * si],
                                   t * t - x[3:3 + M] ** 2 - x[3 + M:] ** 2])

        def cjac(x):
            Y, D, r, (a1, a2, gr, gi) = ev(x)
            t = x[0]
            J = np.zeros((2 + M, 3 + 2 * M))
            J[0, 0] = 1
            J[0, 1], J[0, 2] = -a1, -a2
            J[0, 3:3 + M], J[0, 3 + M:] = -gr, -gi
            J[1, 0] = 2 * t
            J[1, 1] = 2 * (1 - x[1])
            J[1, 2] = -2 * x[2]
            J[2:, 0] = 2 * t
            J[2:, 3:3 + M] = np.diag(-2 * x[3:3 + M])
            J[2:, 3 + M:] = np.diag(-2 * x[3 + M:])
            return J

        res = minimize(lambda x: x[0], x0, jac=lambda x: np.eye(len(x))[0], method="SLSQP",
                       constraints=[{"type": "ineq", "fun": cons, "jac": cjac}],
                       options={"maxiter": maxiter, "ftol": ftol})
        x = res.x
        al = x[1] + 1j * x[2]
        et = x[3:3 + M] + 1j * x[3 + M:]
        return al, et, self.value(al, et)[0], res


def resample(eta, M):
    """Step profile on uniform relative positions -> M uniform steps (nearest/piecewise constant)."""
    L = len(eta)
    pos = (np.arange(M) + 0.5) / M
    return eta[np.minimum((pos * L).astype(int), L - 1)]
