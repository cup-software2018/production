"""WCMD gain matching, check of step 2: the 47 SPE fits on one page, and g of every tank PMT.

    python check_spe_fit.py <run> [--dir DIR]     (needs dark_hist_<run>.root, spe_fit_<run>.npz)
        -> check_spe_fit_<run>.png, check_gain_<run>.png
"""
import argparse
import os
import numpy as np
import ROOT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spe

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("run", type=int)
p.add_argument("--dir", default=".", help="where the inputs are and the plots go")
a = p.parse_args()
RUN = a.run
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#1f1f1e", "#55544f", "#9a9990"
plt.rcParams.update({
    "font.size": 10, "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": MUTED,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": "#e4e3dd",
    "grid.linewidth": 0.6, "axes.axisbelow": True, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "axes.titlelocation": "left",
    "axes.titlesize": 11, "lines.linewidth": 1.4,
})
f = ROOT.TFile.Open(os.path.join(a.dir, f"dark_hist_{RUN}.root"))
z = np.load(os.path.join(a.dir, f"spe_fit_{RUN}.npz"))
P, E, chi2, ndf, live = z["p"], z["e"], z["chi2"], z["ndf"], z["live"]
i = {n: j for j, n in enumerate(z["names"])}
h0 = f.Get("h_q0")
ny = h0.GetNbinsY()
e = np.array([h0.GetYaxis().GetBinLowEdge(j) for j in range(1, ny + 2)])
x = 0.5 * (e[1:] + e[:-1])
w = e[1] - e[0]
xx = np.linspace(spe.LO, spe.HI, 400)

# ── the 47 fits ──────────────────────────────────────────────────────────────
fig, axs = plt.subplots(7, 7, figsize=(24, 21), sharex=True)
for k, ax in enumerate(axs.flat):
    if k >= 47:
        ax.axis("off")
        continue
    pid = k + 1
    y = sum(np.array([f.Get(f"h_q{r}").GetBinContent(pid, j) for j in range(1, ny + 1)]) for r in range(2))
    nz = y > 0
    ax.errorbar(x[nz], y[nz], yerr=np.sqrt(y[nz]), fmt="o", ms=1.8, lw=0.6, color=INK, capsize=0)
    q = P[k]
    if np.all(np.isfinite(q)):
        ax.plot(xx, spe.model(w)(xx, *q), color=INK2, lw=1.0, ls="--")
        ax.plot(xx, spe.g(xx, q[i["a0"]], q[i["m0"]], q[i["sp"]], w), color=MUTED, lw=1.0)
        for n, c in ((1, BLUE), (2, ORANGE), (3, GREEN)):
            ax.plot(xx, spe.pe(xx, q[i[f"a{n}"]], n, q[i["mu"]], q[i["s1"]], q[i["sp"]], w), color=c, lw=1.0)
        ax.set_title(f"PID {pid}: g {q[i['mu']]:.1f} ± {E[k, i['mu']]:.1f}, χ²/ndf {chi2[k]:.0f}/{ndf[k]}", fontsize=9)
    else:
        ax.set_title(f"PID {pid}: fit failed", fontsize=9)
    ax.set_yscale("log")
    ax.set_ylim(0.5, 2 * max(1, y.max()))
    ax.set_xlim(spe.LO, spe.HI)
    ax.tick_params(labelsize=8)
for ax in axs[-1]:
    ax.set_xlabel("dark_q  [ADC·sample]", fontsize=9)
fig.suptitle(f"WCMD tank PMTs, SPE fits of the dark pulses (grey pedestal, blue 1, orange 2, green 3 p.e.) — "
             f"run {RUN}", x=0.01, ha="left", color=INK, fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.975])
fig.savefig(os.path.join(a.dir, f"check_spe_fit_{RUN}.png"), dpi=90)

# ── g of every tank PMT ──────────────────────────────────────────────────────
g, ge = P[:, i["mu"]], E[:, i["mu"]]
ok = np.isfinite(g)
pid = np.arange(1, 48)
gm = np.median(g[ok])
fig, ax = plt.subplots(1, 3, figsize=(17, 4.6), gridspec_kw=dict(width_ratios=[2.2, 1, 1]))
ax[0].errorbar(pid[ok], g[ok], yerr=ge[ok], fmt="o", ms=4, color=BLUE, capsize=0)
ax[0].axhline(gm, color=MUTED, lw=1, ls="--")
ax[0].set_xlabel("PID")
ax[0].set_ylabel("g = 1 p.e. mean  [ADC·sample / p.e.]")
ax[0].set_title(f"gain per tank PMT, median {gm:.1f}, spread {np.std(g[ok]) / gm:.1%}")
ax[1].hist(g[ok], np.linspace(gm * 0.85, gm * 1.15, 31), histtype="step", color=BLUE, lw=1.6)
ax[1].set_xlabel("g  [ADC·sample / p.e.]")
ax[1].set_ylabel("PMTs")
ax[1].set_title(f"{ok.sum()} PMTs")
ax[2].plot(pid[ok], (chi2 / ndf)[ok], "o", ms=4, color=ORANGE)
ax[2].set_xlabel("PID")
ax[2].set_ylabel("χ² / ndf")
ax[2].set_title("fit quality")
fig.suptitle(f"WCMD PMT gain from the dark pulses — run {RUN}", x=0.01, ha="left", color=INK)
fig.tight_layout()
fig.savefig(os.path.join(a.dir, f"check_gain_{RUN}.png"), dpi=130)
print(f"-> check_spe_fit_{RUN}.png, check_gain_{RUN}.png")
