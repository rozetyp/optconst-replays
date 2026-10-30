"""Red controls for the C86 replay: write three altered certificates; c86_replay.py must reject each.
Usage: python3 red_controls_86a.py certs/c86/certificate/c86_A19d7_robust_1.8084990.json outdir"""
import json, sys, os
from fractions import Fraction as F
src, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
def rebound(d):
    e = [F(x) for x in d['edges']]; den = [F(x) for x in d['densities']]
    m = sum(den[i] * (e[i + 1] - e[i]) for i in range(len(den))); T = sum(den[i] * (e[i + 1] ** 2 - e[i] ** 2) / 2 for i in range(len(den)))
    d['bound'] = str(T / m); return d
d = json.load(open(src)); nz = [i for i, x in enumerate(d['densities']) if F(x) > 0]; k = nz[len(nz) // 2]
d['densities'][k] = str(F(d['densities'][k]) * F(99, 100)); json.dump(rebound(d), open(f'{out}/density_cell{k}_minus1pct.json', 'w'))
d = json.load(open(src)); j = max(range(len(d['b'])), key=lambda i: F(d['b'][i])); d['b'][j] = str(F(d['b'][j]) * F(1001, 1000))
json.dump(d, open(f'{out}/b{j}_plus0.1pct.json', 'w'))
d = json.load(open(src)); d['bound'] = str(F(d['bound']) - F(1, 10 ** 15)); json.dump(d, open(f'{out}/bound_minus1e-15.json', 'w'))
print('written:', sorted(os.listdir(out)))
