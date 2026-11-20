"""Monte Carlo computations for the remarks of the note on Conjecture 3.8 of Shahverdi, Strom and Anden
(arXiv:2409.04868v2). Seeded: numpy default_rng(1), 200000 noise samples shared by every z.
Loss (2.11): L(z;x,tau) = E 1/2 ||z - sigma(z; x+eps)||^2, sigma = the shift of x+eps best correlated with z.
=> L = 1/2||z||^2 + 1/2 E||x+eps||^2 - E max_j <R_j(x+eps), z>. On the manifold ||z|| is fixed, so
larger loss <=> smaller score F(z) = E max_j <R_j(x+eps), z>.
Every signal is built from its time samples; the FFT is used only to check the power spectrum.
"""
import itertools
import numpy as np

L = 7
t = np.arange(L)
w = np.array([1.0, 1.0, 2.0])


def sig(p):
    return sum(np.sqrt(2 * w[k] / L) * np.cos((2 * (k + 1) * t + p[k]) * np.pi / L) for k in range(3))


def shifts(v):
    return np.array([np.roll(v, -j) for j in range(L)])


x = sig((0, 0, 0))
P = np.abs(np.fft.fft(x)) ** 2


def score0(z):
    return max(shifts(x) @ z)


def score_mc(z, tau, G):
    X = x[None, :] + tau * G
    best = np.max(np.stack([np.roll(X, -j, axis=1) @ z for j in range(L)]), axis=0)
    return best.mean(), best.std() / np.sqrt(len(best))


# sign patterns via FFT sign flips of x (Prop 3.7 points), independent of the phase parametrisation
xf = np.fft.fft(x)
pats = []
for s in itertools.product([1, -1], repeat=3):
    zf = xf.copy()
    for k in range(3):
        zf[k + 1] *= s[k]; zf[L - k - 1] *= s[k]
    pats.append((s, np.real(np.fft.ifft(zf))))
zstar = sig((0, 5, 10))
print('power spectra equal:', all(np.allclose(np.abs(np.fft.fft(z)) ** 2, P) for _, z in pats),
      np.allclose(np.abs(np.fft.fft(zstar)) ** 2, P), 'means', round(x.sum(), 12), round(zstar.sum(), 12))
print('||x||^2', round(x @ x, 10))
for s, z in pats:
    print('sign', s, 'F0 = %.6f' % score0(z))
print('witness (0,5,10) F0 = %.6f' % score0(zstar))
gap = min(score0(z) for _, z in pats) - score0(zstar)
print('noiseless gap to every sign pattern %.6f ; bound 2*sqrt(7)*tau at tau=0.01: %.6f' % (gap, 2 * np.sqrt(7) * 0.01))

rng = np.random.default_rng(1)
G = rng.standard_normal((200000, L))
for tau in (0.01, 0.1, 0.3):
    m, e = score_mc(zstar, tau, G)
    lo = min(score_mc(z, tau, G)[0] for _, z in pats)
    print('tau %.2f  witness F=%.5f+-%.5f  min over sign patterns F=%.5f' % (tau, m, e, lo))

# integer witness
u = np.array([0, 1, 0, 2, 0, 2, 1]); xi = 7 * u - 6
zi = xi[(-t) % L]
print('int: ||x||^2', xi @ xi, 'same power', np.allclose(np.abs(np.fft.fft(zi)) ** 2, np.abs(np.fft.fft(xi)) ** 2),
      'F0(-x)', max(shifts(xi) @ (-xi)), 'F0(z)', max(shifts(xi) @ zi))
