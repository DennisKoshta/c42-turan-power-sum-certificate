"""Rejection controls for check_tau.py (plain script; run with the flint venv: python -I test_check_tau.py CERT.json).
Every corruption must FAIL; the original must PASS."""
import sys, json, copy, subprocess, tempfile, os
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
cert = json.load(open(sys.argv[1]))


def run(d):
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
        json.dump(d, f)
    r = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'check_tau.py'), f.name, '--bits', '256',
                        '--terms', '250', '--eps', '1e-35'], capture_output=True, text=True)
    os.unlink(f.name)
    return r.returncode, json.loads(r.stdout).get('reason', 'PASS')


cases = []
c = copy.deepcopy(cert); cases.append(('original', c, 0))
c = copy.deepcopy(cert); r = Fraction(c['C']) - Fraction(1, 10**9); c['C'] = str(r); cases.append(('C lowered 1e-9 (below |Y|/D)', c, 1))
c = copy.deepcopy(cert); c['eta'][3] = [str(Fraction(c['C'])), '0']; cases.append(('|eta_3| = C', c, 1))
c = copy.deepcopy(cert); c['alpha'] = [c['alpha'][0], str(Fraction(c['alpha'][1]) + Fraction(1, 1000))]; cases.append(('alpha.im +1e-3', c, 1))
c = copy.deepcopy(cert); c['breaks'][-1] = str(Fraction(c['breaks'][-1]) - Fraction(1, 10**6)); cases.append(('last break != 1-tau', c, 1))
c = copy.deepcopy(cert); tau = Fraction(1, 4); M = len(c['eta']); c['breaks'] = [str(tau + Fraction(j, M) * (1 - 2 * tau)) for j in range(M + 1)]; cases.append(('tau = 1/4', c, 1))
c = copy.deepcopy(cert); c['eta'] = c['eta'][:-1]; cases.append(('missing eta', c, 1))
c = copy.deepcopy(cert); c['C'] = 0.7; cases.append(('float C', c, 1))
c = copy.deepcopy(cert); c['eta'][0], c['eta'][1] = c['eta'][1], c['eta'][0]; cases.append(('eta_0, eta_1 swapped', c, 1))
bad = 0
for name, d, want in cases:
    code, why = run(d)
    ok = (code == 0) == (want == 0)
    bad += not ok
    print(f"{'ok ' if ok else 'BAD'} {name:32s} exit={code} {why[:70]}")
print('ALL CONTROLS BEHAVE' if not bad else f'{bad} CONTROLS MISBEHAVE')
sys.exit(1 if bad else 0)
