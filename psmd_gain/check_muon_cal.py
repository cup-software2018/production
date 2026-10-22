"""PSMD gain matching, check of step 3: the light-collection factor c of every channel, and the
matching per channel on data c was not taken from (c from the even subruns, applied to the odd ones).

    python check_muon_cal.py <run> [--dir DIR]     -> check_muon_cal_<run>.png, check_odd_match_<run>.png
"""
import matplotlib.pyplot as plt
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")


p = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("run", type=int)
p.add_argument("--dir", default=".",
               help="where the inputs are and the plots go")
a = p.parse_args()
RUN = a.run
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED = "#1f1f1e", "#55544f", "#9a9990"
plt.rcParams.update({
    "font.size": 10, "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": MUTED,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": "#e4e3dd",
    "grid.linewidth": 0.6, "axes.axisbelow": True, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "axes.titlelocation": "left",
    "axes.titlesize": 11, "lines.linewidth": 1.8,
})
z = np.load(os.path.join(a.dir, f"muon_cal_{RUN}.npz"))
c, ok, chk = z["c"], z["ok"], z["chk"]
pid = np.arange(1, 521)
E = z["edges"]

fig, ax = plt.subplots(figsize=(11, 4.6))
ax.plot(pid[ok], c[ok], "o", ms=3, color=BLUE,
        label=f"own factor ({ok.sum()})")
fb = np.isfinite(c) & ~ok
ax.plot(pid[fb], c[fb], "o", ms=4, mfc="none", color=ORANGE,
        label=f"too few muons: detector mean or 1 ({fb.sum()})")
ax.axhline(1, color=MUTED, lw=1, ls="--")
for lo, hi in [(205, 260), (465, 520)]:
    ax.axvspan(lo - 0.5, hi + 0.5, color=MUTED, alpha=0.15, lw=0)
ax.set_xlim(0, 521)
ax.set_xlabel("PID")
ax.set_ylabel("c = channel median / median of channel medians")
ax.set_title(
    f"light-collection factor per channel, rms {np.nanstd(c[ok]):.1%}")
ax.legend(loc="upper right")
fig.suptitle(
    f"PSMD per-channel muon response — run {RUN}", x=0.01, ha="left", color=INK)
fig.tight_layout()
fig.savefig(os.path.join(a.dir, f"check_muon_cal_{RUN}.png"), dpi=130)


# ── matching per channel on the odd subruns ──────────────────────────────────
def quantile(h, q):
    """q-quantile of a histogram on the edges E (linear inside the bin)."""
    cs = np.cumsum(h)
    t = q * cs[-1]
    k = np.searchsorted(cs, t)
    lo = cs[k - 1] if k > 0 else 0.0
    return E[k] + (t - lo) / max(h[k], 1e-12) * (E[k + 1] - E[k])


def response(H):
    """channel median / median of channel medians, and its statistical error, on the chk channels."""
    n = H.sum(1)
    med = np.full(520, np.nan)
    err = np.full(520, np.nan)
    for i in np.where(chk)[0]:
        d = 0.5 / np.sqrt(n[i])
        med[i] = quantile(H[i], 0.5)
        err[i] = 0.5 * (quantile(H[i], 0.5 + d) - quantile(H[i], 0.5 - d))
    m = np.median(med[chk])
    return med / m, err / m


r0, e0 = response(z["H_odd_gain"])
r1, e1 = response(z["H_odd_cal"])
fig, ax = plt.subplots(figsize=(11, 4.6))
ax.errorbar(pid[chk], r0[chk], yerr=e0[chk], fmt="o", ms=2.5, lw=0.6, color=BLUE, alpha=0.6,
            label=f"g only, rms {np.std(r0[chk]):.1%}")
ax.errorbar(pid[chk], r1[chk], yerr=e1[chk], fmt="o", ms=2.5, lw=0.6, color=ORANGE,
            label=f"g and c, rms {np.std(r1[chk]):.1%} (statistical error median {np.median(e1[chk]):.1%})")
ax.axhline(1, color=MUTED, lw=1, ls="--")
for lo, hi in [(205, 260), (465, 520)]:
    ax.axvspan(lo - 0.5, hi + 0.5, color=MUTED, alpha=0.15, lw=0)
ax.set_xlim(0, 521)
ax.set_xlabel("PID")
ax.set_ylabel("channel median / median of channel medians")
ax.set_title(f"odd subruns, c from the even subruns: {chk.sum()} channels")
ax.legend(loc="upper right")
fig.suptitle(
    f"PSMD gain matching per channel on independent data — run {RUN}", x=0.01, ha="left", color=INK)
fig.tight_layout()
fig.savefig(os.path.join(a.dir, f"check_odd_match_{RUN}.png"), dpi=130)
print(f"odd subruns, {chk.sum()} channels: rms g only {np.std(r0[chk]):.2%}, g and c {np.std(r1[chk]):.2%}; "
      f"median statistical error {np.median(e1[chk]):.2%}")
