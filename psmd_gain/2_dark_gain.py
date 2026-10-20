"""PSMD gain matching, step 2: per-channel SiPM gain from the dark pulses of the random triggers.

    python 2_dark_gain.py <run> [--dir DIR]       (needs ped_hist_<run>.npz)  -> gain_<run>.npz

g of a PID: slope of the 1, 2, 3 p.e. peak positions of the fitted dark-pulse area (pul_area),
in ADC*sample per p.e. Each peak: Gaussian fit within ±0.35 g of its expected position (two passes).
"""
import argparse
import os
import numpy as np
from scipy.optimize import curve_fit

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("run", type=int)
p.add_argument("--dir", default=".", help="where the input is and the output goes")
a = p.parse_args()
RUN = a.run
h = np.load(os.path.join(a.dir, f"ped_hist_{RUN}.npz"))


def gauss(x, a, m, s):
    return a * np.exp(-0.5 * ((x - m) / s) ** 2)


def peak(y, edges, m0, half):
    """Gaussian fit of the peak near m0; returns (mean, err, sigma) or None."""
    x = 0.5 * (edges[1:] + edges[:-1])
    for _ in range(2):
        s = (x > m0 - half) & (x < m0 + half) & (y > 0)
        if s.sum() < 8:
            return None
        try:
            p, c = curve_fit(gauss, x[s], y[s], p0=[y[s].max(), m0, half / 2],
                             sigma=np.sqrt(y[s] + 1), maxfev=4000)
        except RuntimeError:
            return None
        if not (m0 - half < p[1] < m0 + half) or p[2] <= 0:
            return None
        m0, half = p[1], min(half, 1.5 * abs(p[2]))
    return p[1], np.sqrt(c[1, 1]), abs(p[2])


def line(n, mu, err):
    """Weighted straight line mu = o + g n; returns (g, g_err, o)."""
    w = 1 / np.asarray(err) ** 2
    A = np.vstack([np.ones_like(n), n]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    o, g = cov @ (A.T @ (w * mu))
    return g, np.sqrt(cov[1, 1]), o


N = h["area"].shape[0]
out = {k: np.full(N, np.nan) for k in ["g", "g_err", "o", "w1"]}
ea = h["area_edges"]
xa = 0.5 * (ea[1:] + ea[:-1])
for i in range(N):
    ya = h["area"][i]
    if ya.sum() < 1e4:
        continue                                        # dead channel: no dark pulses
    s = (xa > 150) & (xa < 500)
    g0 = xa[s][np.argmax(ya[s])]                        # 1 p.e. pulse-area peak
    pk = [peak(ya, ea, n * g0, 0.35 * g0) for n in (1, 2, 3)]
    if all(pk):
        n = np.array([1., 2., 3.])
        g, ge, o = line(n, [p[0] for p in pk], [p[1] for p in pk])
        out["g"][i], out["g_err"][i], out["o"][i] = g, ge, o
        out["w1"][i] = pk[0][2] / pk[0][0]

out["pid"] = np.arange(1, N + 1)
np.savez(os.path.join(a.dir, f"gain_{RUN}.npz"), **out)

g = out["g"]
nlive = int(np.sum(h["area"].sum(1) >= 1e4))
print(f"live {nlive}, g fitted {np.sum(~np.isnan(g))}")
print(f"g median {np.nanmedian(g):.1f}, spread {np.nanstd(g) / np.nanmedian(g):.3%}, "
      f"stat err median {np.nanmedian(out['g_err']):.2f}, offset median {np.nanmedian(out['o']):.1f}")
print("failed fits (live):", [i + 1 for i in range(N) if h["area"][i].sum() >= 1e4 and np.isnan(g[i])])
