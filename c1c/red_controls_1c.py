"""Red controls for the C_1c replay: each altered certificate must be refused or have its claim refuted.
Usage: python3 red_controls_1c.py certs/c1c/K24576_0.4103333.json"""
import json, subprocess, sys, os, tempfile
from fractions import Fraction
src = json.load(open(sys.argv[1])); here = os.path.dirname(os.path.abspath(__file__))
def run(cert):
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f: json.dump(cert, f); p = f.name
    r = subprocess.run([sys.executable, '-I', os.path.join(here, 'bigint_1c.py'), p], capture_output=True, text=True)
    os.unlink(p); return r.returncode, (r.stdout + r.stderr).strip().splitlines()[-1]
true = Fraction(src['score_num'], src['score_den'])
c = dict(src); up = true + Fraction(2, 10**15); c['score_num'], c['score_den'] = up.numerator, up.denominator
print('claim 2e-15 above the exact value :', run(c)[1])
c = dict(src); c['weights_int'] = list(src['weights_int']); c['weights_int'][100] = -1
print('a negative weight                 :', 'refused' if run(c)[0] != 0 else 'ACCEPTED (bad)')
c = dict(src); c['weights_int'] = list(src['weights_int']); c['weights_int'][100] = 0.5
print('a non-integer weight              :', 'refused' if run(c)[0] != 0 else 'ACCEPTED (bad)')
