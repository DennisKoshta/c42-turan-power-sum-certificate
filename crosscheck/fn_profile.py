"""Exact finite-n value of a fixed continuum profile (no T_d used), and its extrapolation in n.
usage: python -I fn_profile.py OPT.json TAU M n1,n2,...   (OPT.json from cont_scan/cont4_scan, or a certificate)
Fit ratio(n) = c + b n^-p through the last three n (geometric n): c is the n -> infinity estimate."""
import sys, json, time
import numpy as np
from fractions import Fraction
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from finite_n import Model

opt, tau_s, M, ns = sys.argv[1], sys.argv[2], int(sys.argv[3]), [int(x) for x in sys.argv[4].split(',')]
d = json.load(open(opt))
d = d.get(tau_s, d)
fl = lambda x: float(Fraction(x)) if isinstance(x, str) else float(x)
alpha = fl(d['alpha'][0]) + 1j * fl(d['alpha'][1])
eta = np.array([fl(x) + 1j * fl(y) for x, y in d['eta']])
assert len(eta) == M
rs = []
for n in ns:
    t = time.time()
    r = Model(n, Fraction(tau_s), M).evaluate(alpha, eta)[-1]
    rs.append(r)
    print(f"n={n} ratio={r:.12f} {time.time()-t:.0f}s", flush=True)
if len(rs) >= 3:
    (n1, n2, n3), (r1, r2, r3) = ns[-3:], rs[-3:]
    p = np.log((r1 - r2) / (r2 - r3)) / np.log(n2 / n1)
    b = (r2 - r3) / (n2 ** -p - n3 ** -p)
    print(f"fit p={p:.3f} limit={r3 - b * n3 ** -p:.10f}")
