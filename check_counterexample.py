"""Exact interval check of the counterexample to Conjecture 3.8 of Shahverdi, Strom and Anden (arXiv:2409.04868).

Signal length 7, powers w = (1, 1, 2) in modes 1, 2, 3. For phases p = (p1, p2, p3) (integers, in units of pi/7)
z_p(t) = sum_k sqrt(2 w_k / 7) cos((2 k t + p_k) pi / 7), and the noiseless score is
F0(z_p) = max_j sum_k w_k cos((2 k j + p_k) pi / 7).
Every angle is a multiple of pi/7, so every cosine is +-1 or +-c_i, c_i = cos(2 pi i / 7), i = 1, 2, 3; these are
the three roots of 8t^3 + 4t^2 - 4t - 1. The script encloses each root in an interval of width 1e-6 by an exact
sign change (rational arithmetic), adds the enclosures, and prints for each p an interval containing F0(z_p).
A floating-point recomputation of every score from the time samples is printed beside it.
"""
from fractions import Fraction as Q
import math
import itertools

f = lambda t: 8 * t ** 3 + 4 * t ** 2 - 4 * t - 1
ENC = {}
for i, a in zip((1, 2, 3), (623489, -222521, -900969)):
    lo, hi = Q(a, 10 ** 6), Q(a + 1, 10 ** 6)
    assert f(lo) * f(hi) < 0, i
    ENC[i] = (lo, hi)
# the three intervals are disjoint and a cubic has three roots, so each holds exactly one root; order c1 > c2 > c3
assert ENC[1][0] > ENC[2][1] and ENC[2][0] > ENC[3][1]


def cos_enc(m):
    """interval for cos(m pi / 7), m an integer"""
    m %= 14
    if m % 2 == 1:  # cos(m pi/7) = -cos((m - 7) pi/7), and m - 7 is even
        lo, hi = cos_enc(m - 7)
        return -hi, -lo
    i = (m // 2) % 7
    i = min(i, 7 - i)
    if i == 0:
        return Q(1), Q(1)
    return ENC[i]


w = (1, 1, 2)


def score_enc(p):
    best_lo = best_hi = None
    for j in range(7):
        lo = hi = Q(0)
        for k in range(3):
            a, b = cos_enc(2 * (k + 1) * j + p[k])
            lo += w[k] * a
            hi += w[k] * b
        best_lo = lo if best_lo is None else max(best_lo, lo)
        best_hi = hi if best_hi is None else max(best_hi, hi)
    return best_lo, best_hi


def score_float(p):
    z = [sum(math.sqrt(2 * w[k] / 7) * math.cos((2 * (k + 1) * t + p[k]) * math.pi / 7) for k in range(3)) for t in range(7)]
    x = [sum(math.sqrt(2 * w[k] / 7) * math.cos(2 * (k + 1) * t * math.pi / 7) for k in range(3)) for t in range(7)]
    return max(sum(x[(t + j) % 7] * z[t] for t in range(7)) for j in range(7)), sum(v * v for v in z)


rows = [s for s in itertools.product((0, 7), repeat=3)] + [(0, 5, 10)]
res = {}
for p in rows:
    lo, hi = score_enc(p)
    fl, nrm = score_float(p)
    assert lo <= Q(fl) + Q(1, 10 ** 9) and Q(fl) - Q(1, 10 ** 9) <= hi
    res[p] = (lo, hi)
    print('p =', p, ' F0 in [%.6f, %.6f]  float %.9f  |z|^2 %.9f' % (lo, hi, fl, nrm))
gap = min(res[s][0] for s in rows[:8]) - res[(0, 5, 10)][1]
print('gap >= %s = %.6f ; 2*sqrt(7)/100 = %.6f' % (gap, gap, 2 * math.sqrt(7) / 100))
assert gap ** 2 > Q(28, 10 ** 4) and gap > 0  # gap > 2 sqrt(7) / 100  <=>  gap^2 > 28 / 10^4

# integer example of the remark
u = (0, 1, 0, 2, 0, 2, 1)
xi = [7 * a - 6 for a in u]
zi = [xi[(-t) % 7] for t in range(7)]
dot = lambda a, b, j: sum(a[(t + j) % 7] * b[t] for t in range(7))
assert sum(xi) == 0 and sum(a * a for a in xi) == 238
assert max(dot(xi, [-a for a in xi], j) for j in range(7)) == 154 and max(dot(xi, zi, j) for j in range(7)) == 140
assert 56 ** 2 > 7 * 238  # sqrt(1666) < 56
print('integer example: |x|^2 =', sum(a * a for a in xi), ' F0(-x) =', max(dot(xi, [-a for a in xi], j) for j in range(7)),
      ' F0(z) =', max(dot(xi, zi, j) for j in range(7)))
print('ALL CHECKS PASS')
