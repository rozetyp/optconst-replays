"""Exact C_1c score of an equal-step certificate via one big-integer product (Kronecker substitution).
A_k = sum_i w_i w_{i+k} is the coefficient of x^(n-1+k) in W(x) * W_rev(x); with base 2^B larger than every
coefficient, the product of the two packed integers holds all A_k exactly. Usage: python3 -I bigint_1c.py cert.json"""
import json, sys
from fractions import Fraction
c = json.load(open(sys.argv[1])); K = c['K']; w = c['weights_int']
assert isinstance(K, int) and K > 0 and all(type(x) is int and x >= 0 for x in w) and any(w)
n = len(w); S = sum(w)
B = (max(w) ** 2 * n).bit_length() + 1
P = int(''.join(format(x, f'0{B}b') for x in w), 2)            # w_0 is the most significant digit
Q = int(''.join(format(x, f'0{B}b') for x in reversed(w)), 2)
R = P * Q; mask = (1 << B) - 1
coef = [(R >> (B * j)) & mask for j in range(2 * n - 1)]       # coef[j]: j counted from least significant
# digit position of w_i in P is (n-1-i), of w_l in Q is l; product digit (n-1-i)+l; A_k needs l = i+k -> digit n-1+k
A = [coef[n - 1 + k] if k < n else 0 for k in range(K + 1)]
kmin = min(range(K + 1), key=lambda k: A[k])
score = Fraction(K * A[kmin], S * S)
print(f'K={K} n={n} argmin k={kmin} score={score.numerator}/{score.denominator} = {float(score):.14f}')
if 'score_num' in c:
    claim = Fraction(c['score_num'], c['score_den'])
    print('claimed fraction:', 'EXACT MATCH' if claim == score else ('certified (claim below exact value)' if claim < score else 'REFUTED (claim above exact value)'))
