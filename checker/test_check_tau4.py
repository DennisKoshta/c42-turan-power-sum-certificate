"""Rejection controls for check_tau4.py (review REVIEW_tau_lt_quarter.md R1, R2).
usage: python -I test_check_tau4.py CERT.json   (any valid FORMAT/FORMAT2 certificate; only its profile is reused)"""
import sys, json
from fractions import Fraction as F
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from flint import arb, acb, ctx
import check_tau4 as c
from check_tau import ball

ctx.prec = 256
ok = True
al, a = acb('0.67', '0.6'), arb('0.67')
u1, u2 = ball(F(21, 100)), ball(F(30, 100))
# R1: Phi must refuse a ball where y < 1 is not certain, whatever the analytic flag
for x in [acb('0.6', '0.3'), acb(arb('0.3', '0.29'), arb(0, '0.2'))]:
    for fl in (False, True):
        r = c._phi(x, al, a, u1, u2, 40, fl) is None
        print('R1 refuse', x.str(3), fl, r)
        ok &= r
# R1 converse: a point on the real path is accepted
r = c._phi(acb('0.25'), al, a, u1, u2, 40, False) is not None
print('R1 accept on path', r)
ok &= r
# R2: tau <= 1/5 refused
cert = json.load(open(sys.argv[1]))
for tau in ['1/5', '19/100']:
    d = dict(cert)
    d['breaks'] = [tau] + cert['breaks'][1:-1] + [str(1 - F(tau))]
    try:
        c.check(d, 50, 1e-10)
        print('R2', tau, 'ACCEPTED')
        ok = False
    except ValueError as e:
        print('R2', tau, 'refused:', e)
print('ALL OK' if ok else 'FAILURES')
sys.exit(0 if ok else 1)
