"""Independent replay of a C86 (Schur-Siegel-Smyth trace) upper-bound certificate in the Pinho / Koshta format.
Written from the format description (edges, densities, A, b, bound) and the conditions of Smith's criterion
[Smi24, Prop. 5.15 and Thm 1.5; SO24, Cor. 1.6]; shares no code with the certificate's own checkers.

Conditions (numbering as in the certificate's FORMAT.md):
  E1 edges increasing, d_i >= 0, b_j >= 0, m > 0, min S > 0           exact rationals
  E2 Q_j primitive, irreducible over Z, pairwise non-associate          FLINT factorisation (exact)
  E3 sum b_j deg Q_j < m                                                exact
  E4 int log|Q_j| dmu = m log|lc Q_j| + sum_roots U_mu(r) >= 0          arb balls on certified root enclosures
  E5 f(x) = U_mu(x) - sum b_j log|Q_j(x)| >= 0 for every x in S          rigorous bisection (below)
  E6 T/m == bound                                                       exact
With mu = sum_i d_i 1[e_i, e_{i+1}], U_mu(x) = sum_k c_k G(x - e_k), c_k = d_k - d_{k-1} (d_{-1} = d_N = 0),
G(u) = u log|u| - u, G(0) = 0.

E5 method, on each support cell [e_k, e_{k+1}], bisected adaptively. On a piece I = [m - h, m + h]:
  * the two edge terms c G(x - e) with e in {e_k, e_{k+1}} have constant-sign argument on I, so c G is convex or
    concave there; replace it by its tangent at m (convex) or its chord over I (concave): a linear minorant;
  * the rest R (other edges, all at distance >= one cell, and the polynomial terms) is C^2 on I unless a root of
    some Q_j is in I: R(x) >= R(m) + R'(m)(x - m) - M (x - m)^2 / 2 with M >= sup_I |R''| from ball evaluation;
  * f >= V + D (x - m) - M (x - m)^2 / 2 on I, a concave function of x, so f >= V - |D| h - M h^2 / 2 on I.
  If a Q_j ball over I contains 0, fall back to a termwise bound: exact range of G on each argument interval
  (G is decreasing on [-1, 1], increasing outside) and -b_j log max_I |Q_j|.
Usage: python3 c86_replay.py certificate.json [target_margin] [workers]"""
import sys, json, time
from fractions import Fraction as Fr
from multiprocessing import Pool
import flint
from flint import arb, fmpq, fmpz_poly

flint.ctx.prec = 128

def A_(q): return arb(fmpq(q.numerator, q.denominator))

def load(path):
    d = json.load(open(path))
    e = [Fr(x) for x in d['edges']]; den = [Fr(x) for x in d['densities']]
    A = [[int(c) for c in q] for q in d['A']]; b = [Fr(x) for x in d['b']]
    return e, den, A, b, Fr(d['bound'])

def G_point(u):          # u: arb ball not containing 0, or exact zero
    if u.is_zero(): return arb(0)
    return u * abs(u).log() - u

def G_range(lo, hi):     # exact range of G on [lo, hi] (floats/arb endpoints), returned as an arb ball
    vals = [G_point(lo) if not lo.is_zero() else arb(0), G_point(hi) if not hi.is_zero() else arb(0)]
    if lo.lower() <= -1 and hi.upper() >= -1: vals.append(arb(1))
    if lo.lower() <= 1 and hi.upper() >= 1: vals.append(arb(-1))
    lo_v = min(v.lower() for v in vals); hi_v = max(v.upper() for v in vals)
    return arb((lo_v + hi_v) / 2, (hi_v - lo_v) / 2)

def G_ball(u):
    if not u.contains(0): return G_point(u)
    return G_range(u.lower(), u.upper())

# ---------------- worker globals ----------------
W = {}
def init(path):
    e, den, A, b, _ = load(path)
    N = len(den); c = [den[0]] + [den[k] - den[k - 1] for k in range(1, N)] + [-den[-1]]
    W['e'] = e; W['ea'] = [A_(x) for x in e]; W['c'] = c; W['ca'] = [A_(x) for x in c]
    W['nz'] = [k for k in range(len(c)) if c[k] != 0]
    W['polys'] = [(fmpz_poly(q), fmpz_poly(q).derivative(), fmpz_poly(q).derivative().derivative(), A_(bj)) for q, bj in zip(A, b) if bj > 0]

def peval(p, x):         # Horner with arb
    cs = [int(t) for t in p.coeffs()]; r = arb(0)
    for t in reversed(cs): r = r * x + t
    return r

def lb_taylor(k, a, bb):
    ea, ca, nz = W['ea'], W['ca'], W['nz']
    m = (a + bb) / 2; h = (bb - a) / 2; mA = A_(m); hA = A_(h); I = arb(mA, hA)
    V = arb(0); D = arb(0); M = arb(0)
    for j in nz:
        if j == k or j == k + 1:
            e = W['e'][j]; c = W['c'][j]; cA = ca[j]
            u_side = 1 if j == k else -1          # sign of x - e on the cell
            convex = (c > 0) == (u_side > 0)
            if convex:
                um = mA - ea[j]; V += cA * G_point(um); D += cA * abs(um).log()
            else:
                ga = cA * G_point(A_(a - e)) if a != e else arb(0)
                gb = cA * G_point(A_(bb - e)) if bb != e else arb(0)
                V += (ga + gb) / 2; D += (gb - ga) / (2 * hA)
            continue
        um = mA - ea[j]; uI = I - ea[j]
        V += ca[j] * G_point(um); D += ca[j] * abs(um).log(); M += abs(ca[j] / uI)
    for (P, P1, P2, bj) in W['polys']:
        qI = peval(P, I)
        if qI.contains(0): return None
        qm = peval(P, mA)
        V -= bj * abs(qm).log(); D -= bj * peval(P1, mA) / qm
        M += abs(bj * (peval(P2, I) * qI - peval(P1, I) ** 2) / (qI * qI))
    return (V - abs(D) * hA - M * hA * hA / 2).lower()

def lb_termwise(k, a, bb):
    I = arb(A_((a + bb) / 2), A_((bb - a) / 2)); tot = arb(0)
    for j in W['nz']:
        lo = A_(a - W['e'][j]); hi = A_(bb - W['e'][j])
        g = G_range(lo.lower() if not lo.is_zero() else arb(0), hi.upper() if not hi.is_zero() else arb(0))
        tot += (W['ca'][j] * g).lower()
    for (P, P1, P2, bj) in W['polys']:
        up = abs(peval(P, I)).upper()
        if up <= 0: return arb(10)  # cannot happen for a nonzero polynomial on a nondegenerate interval
        tot -= (bj * up.log()).upper()
    return tot.lower()

def check_cell(args):
    k, target = args; target = arb(target)
    a0, b0 = W['e'][k], W['e'][k + 1]
    stack = [(a0, b0)]; worst = None; pieces = 0; minw = b0 - a0
    while stack:
        a, bb = stack.pop(); pieces += 1
        lb = lb_taylor(k, a, bb)
        if lb is None: lb = lb_termwise(k, a, bb)
        if lb >= target:
            worst = lb if worst is None or lb < worst else worst; continue
        if bb - a < Fr(1, 10 ** 13):
            return (k, False, float(lb), pieces, float(a))
        mid = (a + bb) / 2; stack += [(a, mid), (mid, bb)]; minw = min(minw, (bb - a) / 2)
    return (k, True, float(worst.mid()), pieces, float(minw))

def main():
    path = sys.argv[1]; tval = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0; target = arb(tval)
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    t0 = time.time()
    e, den, A, b, bound = load(path); N = len(den)
    # E1
    assert all(e[i] < e[i + 1] for i in range(N)) and all(x >= 0 for x in den) and all(x >= 0 for x in b)
    m = sum(den[i] * (e[i + 1] - e[i]) for i in range(N)); T = sum(den[i] * (e[i + 1] ** 2 - e[i] ** 2) / 2 for i in range(N))
    S = [i for i in range(N) if den[i] > 0]; assert m > 0 and e[S[0]] > 0
    print(f'E1 ok: N = {N} cells, {len(S)} in support, S within [{e[S[0]]}, {e[S[-1] + 1]}], m = {m}')
    # E2
    polys = [fmpz_poly(q) for q in A]
    for q, P in zip(A, polys):
        assert P.degree() >= 1 and abs(P.content()) == 1, f'not primitive: {q}'
        cont, facs = P.factor(); assert len(facs) == 1 and facs[0][1] == 1, f'reducible: {q}'
    for i in range(len(polys)):
        for j in range(i + 1, len(polys)):
            assert polys[i] != polys[j] and polys[i] != -polys[j], 'associate pair'
    realrooted = all(all(r.imag.contains(0) for r, _ in P.complex_roots()) and len(P.complex_roots()) == P.degree() for P in polys)
    print(f'E2 ok: {len(A)} primitive irreducible pairwise non-associate polynomials (all real-rooted: {realrooted})')
    # E3
    sdeg = sum(bj * (len(q) - 1) for q, bj in zip(A, b)); assert sdeg < m
    print(f'E3 ok: sum b deg / m = {float(sdeg / m):.6f} < 1')
    # E6
    assert T / m == bound, f'E6 FAILED: bound {bound} differs from T/m = {T / m}'; print(f'E6 ok: T/m = bound = {bound} = {float(bound):.10f}')
    # E4
    init(path); minE4 = None
    for q, P, bj in zip(A, polys, b):
        if bj == 0: continue
        val = A_(m) * abs(arb(int(P.coeffs()[-1]))).log()
        for r, mult in P.complex_roots():
            rr = r.real
            if P.degree() == 1:                      # rational root: exact
                cs = [int(t) for t in P.coeffs()]; rq = Fr(-cs[0], cs[1]); U = arb(0)
                for k in W['nz']: U += W['ca'][k] * (A_(rq - W['e'][k]) * abs(A_(rq - W['e'][k])).log() - A_(rq - W['e'][k]) if rq != W['e'][k] else arb(0))
            else:
                U = arb(0)
                for k in W['nz']: U += W['ca'][k] * G_ball(rr - W['ea'][k])
            val += mult * U
        assert val > 0, f'E4 not certified for {q}: {val}'
        minE4 = val if minE4 is None or val.lower() < minE4.lower() else minE4
    print(f'E4 ok: every weighted Q has int log|Q| dmu > 0 (smallest {minE4.str(5)})')
    # E5
    with Pool(workers, initializer=init, initargs=(path,)) as pool:
        res = pool.map(check_cell, [(k, tval) for k in S], chunksize=4)
    bad = [r for r in res if not r[1]]
    pieces = sum(r[3] for r in res); worst = min(r[2] for r in res if r[1]) if len(bad) < len(res) else None
    if bad:
        print(f'E5 FAILED on {len(bad)} cells, e.g. cell {bad[0][0]} near x = {bad[0][4]}: lower bound {bad[0][2]}')
        sys.exit(1)
    print(f'E5 ok: f >= {target.str(3)} on all {len(S)} support cells ({pieces} pieces; smallest certified piece bound {worst:.3e})')
    print(f'CERTIFIED: C86 <= {bound} = {float(bound):.10f}   ({time.time() - t0:.0f} s)')

if __name__ == '__main__':
    main()
