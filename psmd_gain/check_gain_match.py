"""PSMD gain matching, check: does the SiPM gain g (step 2) equalise the channel response?
Clean muons are the reference; c is not applied here.

Per module in an event (all 4 channels): both slabs > 100 p.e. (gain-matched) and the largest
channel < 0.55 of the module total. For those modules, per-PID histograms of the channel charge
before (ADC·sample) and after (p.e.) gain matching; the per-PID medians are compared.

    python check_gain_match.py hist <run> [--prod DIR] [--dir DIR]   -> check_gain_match_<run>.npz
    python check_gain_match.py plot <run> [--dir DIR]                -> check_gain_match_<run>.png

Needs psmd_gain_<run>.txt (its g column).
"""
import argparse
import os
import numpy as np
from psmd_gain import prod_files

pa = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
pa.add_argument("what", choices=("hist", "plot"))
pa.add_argument("run", type=int)
pa.add_argument("--prod", default="/data/amore2/data/PROD", help="production directory")
pa.add_argument("--dir", default=".", help="where the inputs are and the outputs go")
A = pa.parse_args()
PROD = A.prod


def load_gain(run):
    """g by PID-1 from the table (nan for dead channels)."""
    pid, g = np.loadtxt(os.path.join(A.dir, f"psmd_gain_{run}.txt"), comments="#", usecols=(0, 1), unpack=True)
    return g[np.argsort(pid)]


def define_npe(df, run):
    import ROOT
    g = load_gain(run)
    if not hasattr(ROOT, "psmd_g_only"):
        vals = ", ".join("0" if not np.isfinite(x) else f"{x:.4f}" for x in g)
        ROOT.gInterpreter.Declare(f"const float psmd_g_only[{len(g)}] = {{{vals}}};")
    return df.Define("npe", "ROOT::RVec<float> o(nch); for (int i = 0; i < nch; i++) "
                            "{ float gg = psmd_g_only[pid[i] - 1]; o[i] = gg > 0 ? charge[i] / gg : 0.f; } return o;")


def hist(run):
    import ROOT
    ROOT.EnableImplicitMT(16)
    df = define_npe(ROOT.RDataFrame("normal", prod_files(PROD, run)), run)
    ROOT.gInterpreter.Declare("""
    ROOT::RVec<int> CleanMuonCh(const ROOT::RVec<short>& pid, const ROOT::RVec<float>& npe) {
      std::map<int, std::array<int, 4>> mods;
      std::map<int, int> cnt;
      for (size_t i = 0; i < pid.size(); i++) {
        int m = (pid[i] - 1) / 4, k = (pid[i] - 1) % 4;
        if (!cnt.count(m)) mods[m] = {-1, -1, -1, -1};
        mods[m][k] = i; cnt[m]++;
      }
      ROOT::RVec<int> out;
      for (auto& [m, idx] : mods) {
        if (cnt[m] != 4) continue;
        float q[4]; for (int k = 0; k < 4; k++) q[k] = npe[idx[k]];
        float a = q[0] + q[1], b = q[2] + q[3], t = a + b;
        float mx = std::max(std::max(q[0], q[1]), std::max(q[2], q[3]));
        if (a > 100 && b > 100 && mx < 0.55 * t) for (int k = 0; k < 4; k++) out.push_back(idx[k]);
      }
      return out;
    }""")
    d = (df.Define("sel", "CleanMuonCh(pid, npe)").Filter("sel.size() > 0")
           .Define("spid", "ROOT::VecOps::Take(pid, sel)")
           .Define("sq", "ROOT::VecOps::Take(charge, sel)")
           .Define("sn", "ROOT::VecOps::Take(npe, sel)"))
    ha = d.Histo2D(("ha", "", 520, 0.5, 520.5, 1000, 0, 100000), "spid", "sq")
    hn = d.Histo2D(("hn", "", 520, 0.5, 520.5, 1000, 0, 360), "spid", "sn")
    def arr(h):
        return np.array([[h.GetBinContent(i, j) for j in range(1, h.GetNbinsY() + 1)]
                         for i in range(1, h.GetNbinsX() + 1)])
    np.savez(os.path.join(A.dir, f"check_gain_match_{run}.npz"), adc=arr(ha), npe=arr(hn),
             adc_edges=np.linspace(0, 100000, 1001), npe_edges=np.linspace(0, 360, 1001))
    print("done")


def median(h, e):
    c = np.cumsum(h)
    if c[-1] == 0:
        return np.nan
    i = np.searchsorted(c, c[-1] / 2)
    return 0.5 * (e[i] + e[i + 1])


def plot(run):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BLUE, ORANGE = "#2a78d6", "#eb6834"
    INK, INK2, MUTED = "#1f1f1e", "#55544f", "#9a9990"
    plt.rcParams.update({
        "font.size": 10, "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": MUTED,
        "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": "#e4e3dd",
        "grid.linewidth": 0.6, "axes.axisbelow": True, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False, "axes.titlelocation": "left",
        "axes.titlesize": 11, "lines.linewidth": 1.8,
    })
    z = np.load(os.path.join(A.dir, f"check_gain_match_{run}.npz"))
    n = z["adc"].sum(1)
    g0 = load_gain(run)
    ok = (n >= 200) & np.isfinite(g0)
    ma = np.array([median(h, z["adc_edges"]) for h in z["adc"]])
    mn = np.array([median(h, z["npe_edges"]) for h in z["npe"]])
    ra, rn = ma / np.nanmedian(ma[ok]), mn / np.nanmedian(mn[ok])
    pid = np.arange(1, 521)
    g = load_gain(run)
    sa, sn = np.nanstd(ra[ok]), np.nanstd(rn[ok])
    print(f"channels with >= 200 muon entries: {ok.sum()} (median {np.median(n[ok]):.0f} entries)")
    print(f"muon median per channel, rms over channels: before {sa:.2%}, after {sn:.2%}")
    print(f"gain rms over the same channels: {np.nanstd(g[ok]) / np.nanmedian(g[ok]):.2%}")
    print(f"corr(before, gain) = {np.corrcoef(ra[ok], g[ok])[0, 1]:.2f}, "
          f"corr(after, gain) = {np.corrcoef(rn[ok], g[ok])[0, 1]:.2f}")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6), gridspec_kw=dict(width_ratios=[2.2, 1, 1]))
    ax[0].plot(pid[ok], ra[ok], "o", ms=3, color=BLUE, label=f"before (ADC·sample), rms {sa:.1%}")
    ax[0].plot(pid[ok], rn[ok], "o", ms=3, color=ORANGE, label=f"after (p.e.), rms {sn:.1%}")
    ax[0].axhline(1, color=MUTED, lw=1, ls="--")
    ax[0].set_xlim(0, 521)
    ax[0].set_xlabel("PID")
    ax[0].set_ylabel("muon median / median over channels")
    ax[0].set_title("per-channel muon response")
    ax[0].legend(loc="upper right")
    b = np.linspace(0.75, 1.25, 51)
    ax[1].hist(ra[ok], b, histtype="step", color=BLUE, label="before")
    ax[1].hist(rn[ok], b, histtype="step", color=ORANGE, label="after")
    ax[1].set_xlabel("muon median / median over channels")
    ax[1].set_ylabel("channels")
    ax[1].set_title(f"{ok.sum()} channels")
    ax[1].legend(loc="upper right")
    gr = g / np.nanmedian(g[ok])
    ax[2].plot(gr[ok], ra[ok], "o", ms=3, color=BLUE, label="before")
    ax[2].plot(gr[ok], rn[ok], "o", ms=3, color=ORANGE, label="after")
    ax[2].set_xlabel("channel gain / median gain")
    ax[2].set_ylabel("muon median / median over channels")
    ax[2].set_title("response vs gain")
    ax[2].legend(loc="upper left")
    fig.suptitle(f"Gain matching check with clean muons — run {run}", x=0.01, ha="left", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(A.dir, f"check_gain_match_{run}.png"), dpi=130)


if __name__ == "__main__":
    {"hist": hist, "plot": plot}[A.what](A.run)
