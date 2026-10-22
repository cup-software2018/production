"""The single-photoelectron model of the WCMD dark-pulse charge, shared by 2_spe_fit.py and the check.

    n(q) = A0 G(q; m0, sp)                                         pedestal (baseline) excursions above 4σ
         + sum_n=1..3  An G(q; n mu, sqrt(n s1^2 + sp^2))         n photoelectrons

The standard single-photoelectron response without the Poisson weights (the An are free): the n p.e.
peak sits at n mu and is n single photoelectrons of width s1 on the pedestal width sp, the fitted
width of the pedestal peak.  mu (no offset) is the gain, s1 the SPE width.  Fitted with Poisson
weights over [LO, HI).
"""
import numpy as np
from scipy.optimize import curve_fit

LO, HI = -15.0, 160.0
NAMES = ["a0", "m0", "sp", "a1", "mu", "s1", "a2", "a3"]


def g(q, a, m, s, w):
    """Gaussian of area a (counts) per bin of width w."""
    return a * w / (np.sqrt(2 * np.pi) * s) * np.exp(-0.5 * ((q - m) / s) ** 2)


def pe(q, a, n, mu, s1, sp, w):
    return g(q, a, n * mu, np.sqrt(n * s1 ** 2 + sp ** 2), w)


def model(w):
    def f(q, a0, m0, sp, a1, mu, s1, a2, a3):
        return (g(q, a0, m0, sp, w) + pe(q, a1, 1, mu, s1, sp, w) + pe(q, a2, 2, mu, s1, sp, w)
                + pe(q, a3, 3, mu, s1, sp, w))
    return f


def fit(x, y):
    """x: bin centres, y: counts.  (params, errors, chi2, ndf), or None if the fit fails."""
    w = x[1] - x[0]
    s = (x >= LO) & (x < HI)
    big = y[s & (x >= 20)].sum()
    if big < 100:
        return None
    p0 = [y[s & (x < 20)].sum(), 5, 5, big, 37, 14, 0.05 * big, 0.005 * big]
    try:
        p, cov = curve_fit(model(w), x[s], y[s], p0=p0, sigma=np.sqrt(y[s] + 1), absolute_sigma=True,
                           bounds=([0, -10, 1, 0, 15, 3, 0, 0], [np.inf,
                                   20, 20, np.inf, 80, 40, np.inf, np.inf]),
                           maxfev=20000)
    except RuntimeError:
        return None
    chi2 = np.sum(((y[s] - model(w)(x[s], *p)) / np.sqrt(y[s] + 1)) ** 2)
    return p, np.sqrt(np.diag(cov)), chi2, int(s.sum() - len(p))
