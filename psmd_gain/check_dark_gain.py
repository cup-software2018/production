"""PSMD gain matching, check of step 2: the dark-pulse area spectrum of one channel with the
fitted p.e. positions, and the gain of every channel.

    python check_dark_gain.py <run> [--pid N] [--dir DIR]   -> check_dark_gain_<pid>_<run>.png, check_dark_gain_<run>.png
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("run", type=int)
p.add_argument("--dir", default=".", help="where the inputs are and the plots go")
p.add_argument("--pid", type=int, default=1, help="channel whose spectra are shown")
a = p.parse_args()
RUN, PID = a.run, a.pid
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED = "#1f1f1e", "#55544f", "#9a9990"
plt.rcParams.update({
    "font.size": 10, "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": MUTED,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": "#e4e3dd",
    "grid.linewidth": 0.6, "axes.axisbelow": True, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "axes.titlelocation": "left",
    "axes.titlesize": 11, "lines.linewidth": 2,
})
h = np.load(os.path.join(a.dir, f"ped_hist_{RUN}.npz"))
g = np.load(os.path.join(a.dir, f"gain_{RUN}.npz"))
pid = g["pid"]

# ── one channel: the area spectrum with the fitted p.e. positions ────────────
i = PID - 1
fig, x = plt.subplots(figsize=(7, 4.4))
x.stairs(h["area"][i], h["area_edges"], color=BLUE, lw=1.5)
for n in (1, 2, 3):
    v = g["o"][i] + n * g["g"][i]
    x.axvline(v, color=MUTED, lw=1, ls="--")
    x.text(v, 1.0, f" {n} p.e.", transform=x.get_xaxis_transform(), va="top", color=INK2, fontsize=9)
x.set_yscale("log")
x.set_xlim(0, 1500)
x.set_xlabel("fitted dark-pulse area  [ADC·sample]")
x.set_ylabel("entries / bin")
x.set_title(f"gain {g['g'][i]:.1f} ± {g['g_err'][i]:.1f}")
fig.suptitle(f"PID {PID}, run {RUN} random triggers — dashed: fitted p.e. positions",
             x=0.01, ha="left", color=INK)
fig.tight_layout()
fig.savefig(os.path.join(a.dir, f"check_dark_gain_{PID}_{RUN}.png"), dpi=130)

# ── all channels ─────────────────────────────────────────────────────────────
live = ~np.isnan(g["g"])
gm = np.nanmedian(g["g"])
fig, ax = plt.subplots(1, 2, figsize=(13, 4.4), gridspec_kw=dict(width_ratios=[2.2, 1]))
ax[0].plot(pid[live], g["g"][live], "o", ms=3, color=BLUE)
ax[0].axhline(gm, color=MUTED, lw=1, ls="--")
for lo, hi in [(205, 260), (465, 520)]:
    ax[0].axvspan(lo - 0.5, hi + 0.5, color=MUTED, alpha=0.15, lw=0)
    ax[0].text(0.5 * (lo + hi), gm * 1.085, "no SiPM", ha="center", color=INK2, fontsize=9)
ax[0].set_ylim(gm * 0.88, gm * 1.1)
ax[0].set_xlim(0, 521)
ax[0].set_xlabel("PID")
ax[0].set_ylabel("gain  [ADC·sample / p.e.]")
ax[0].set_title(f"gain per channel, median {gm:.1f}")
b = np.linspace(gm * 0.88, gm * 1.12, 49)
ax[1].hist(g["g"][live], b, histtype="step", color=BLUE)
ax[1].set_xlabel("gain  [ADC·sample / p.e.]")
ax[1].set_ylabel("channels")
ax[1].set_title(f"{live.sum()} channels, rms {np.nanstd(g['g']) / gm:.1%}")
fig.suptitle(f"PSMD SiPM gain from dark pulses, run {RUN}", x=0.01, ha="left", color=INK)
fig.tight_layout()
fig.savefig(os.path.join(a.dir, f"check_dark_gain_{RUN}.png"), dpi=130)
