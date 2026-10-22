"""PSMD gain matching, step 3: per-channel light-collection factor c from clean muons, on top of
the SiPM gain g of step 2, and the final table.

    q = charge / (g[pid-1] * c[pid-1])

c is the channel's muon response relative to the median channel: per-PID median of the
gain-matched charge of clean muons, divided by the median over channels. Clean muon module
(all 4 channels in the event): both slabs > 100 p.e. and the largest channel < 0.55 of the
module total, always on the gain-matched charge -- the selection is NOT re-done with the
calibrated charge (doing so feeds back: a lowered factor inflates the channel, lets
lower-charge modules through, lowers the median, and runs away).

Check: c from the even subruns is applied to the odd ones (independent statistics). The
final c is then derived from all subruns.

    python 3_muon_cal.py <run> [--prod DIR] [--dir DIR]     (needs gain_<run>.npz)
        -> muon_cal_<run>.npz        intermediate: c, its checks, histograms
        -> psmd_gain_<run>.txt       the table analyses read: pid g c use

use = 1 when g and c exist, c comes from the channel's own muons, and all four channels of the
detector are alive (a detector with a dead channel biases the clean-muon selection).
"""
import argparse
import os
import numpy as np
from psmd_gain import prod_files

p = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("run", type=int)
p.add_argument("--prod", default="/data/amore2/data/PROD",
               help="production directory")
p.add_argument("--dir", default=".",
               help="where the input is and the output goes")
A = p.parse_args()
PROD = A.prod


def load_gain(run):
    return np.load(os.path.join(A.dir, f"gain_{run}.npz"))["g"]


NMIN = 200                      # muon entries a channel needs for its own factor
EDGES = np.linspace(0, 400, 2001)


def median(h, e):
    c = np.cumsum(h)
    if c[-1] == 0:
        return np.nan
    i = np.searchsorted(c, c[-1] / 2)
    return 0.5 * (e[i] + e[i + 1])


def one_pass(run, c, parity):
    """Clean-muon selection on the gain-matched charge; histogram charge / (g * c) per PID.
    parity 0/1: even/odd subruns."""
    import ROOT
    g = load_gain(run)
    good = np.isfinite(g) & (g > 0)
    name = f"psmd_cal_{run}_{one_pass.n}"
    one_pass.n += 1
    def arr(v): return "{" + ", ".join(f"{x:.5f}" for x in v) + "}"
    ROOT.gInterpreter.Declare(f"const float G_{name}[520] = {arr(np.where(good, g, 0))};\n"
                              f"const float C_{name}[520] = {arr(np.where(good, c, 1))};")
    ROOT.gInterpreter.Declare(f"""
    ROOT::RVec<float> Q_{name}(const ROOT::RVec<short>& pid, const ROOT::RVec<float>& q, bool cal) {{
      ROOT::RVec<float> o(q.size());
      for (size_t i = 0; i < q.size(); i++) {{
        float s = G_{name}[pid[i] - 1] * (cal ? C_{name}[pid[i] - 1] : 1.f);
        o[i] = s > 0 ? q[i] / s : 0.f;
      }}
      return o;
    }}
    ROOT::RVec<int> Sel_{name}(const ROOT::RVec<short>& pid, const ROOT::RVec<float>& q) {{
      std::map<int, std::array<int, 4>> mods; std::map<int, int> cnt;
      for (size_t i = 0; i < pid.size(); i++) {{
        int m = (pid[i] - 1) / 4, k = (pid[i] - 1) % 4;
        if (!cnt.count(m)) mods[m] = {{-1, -1, -1, -1}};
        mods[m][k] = i; cnt[m]++;
      }}
      ROOT::RVec<int> out;
      for (auto& [m, idx] : mods) {{
        if (cnt[m] != 4) continue;
        float x[4]; for (int k = 0; k < 4; k++) x[k] = q[idx[k]];
        float a = x[0] + x[1], b = x[2] + x[3], t = a + b;
        float mx = std::max(std::max(x[0], x[1]), std::max(x[2], x[3]));
        if (a > 100 && b > 100 && mx < 0.55 * t) for (int k = 0; k < 4; k++) out.push_back(idx[k]);
      }}
      return out;
    }}""")
    df = ROOT.RDataFrame("normal", prod_files(PROD, run))
    d = (df.Filter(f"subrun % 2 == {parity}")
           .Define("qg", f"Q_{name}(pid, charge, false)")
           .Define("qc", f"Q_{name}(pid, charge, true)")
           .Define("sel", f"Sel_{name}(pid, qg)").Filter("sel.size() > 0")
           .Define("spid", "ROOT::VecOps::Take(pid, sel)")
           .Define("sq", "ROOT::VecOps::Take(qc, sel)"))
    h = d.Histo2D(("h", "", 520, 0.5, 520.5, len(EDGES) -
                  1, EDGES[0], EDGES[-1]), "spid", "sq")
    return np.array([[h.GetBinContent(i, j) for j in range(1, h.GetNbinsY() + 1)]
                     for i in range(1, 521)])


one_pass.n = 0


def residual(H, ok):
    med = np.array([median(h, EDGES) for h in H])
    return med / np.median(med[ok]), med


def main(run):
    import ROOT
    ROOT.EnableImplicitMT(16)
    g = load_gain(run)
    live = np.isfinite(g)
    one = np.ones(520)
    He = one_pass(run, one, 0)                     # even: derive
    ne = He.sum(1)
    ok = live & (ne >= NMIN / 2)
    re, mede = residual(He, ok)
    c = np.where(ok, re, np.nan)
    # too few muons: mean of the module's other channels
    bad = live & ~ok
    for i in np.where(bad)[0]:
        oth = [j for j in range(4 * (i // 4), 4 * (i // 4) + 4) if ok[j]]
        c[i] = np.mean(c[oth]) if oth else 1.0
    Ho0 = one_pass(run, one, 1)                    # odd: gain matching only
    Ho1 = one_pass(run, np.where(np.isfinite(c), c, 1), 1)   # odd: with c
    no = Ho0.sum(1)
    chk = ok & (no >= NMIN / 2)
    r0, _ = residual(Ho0, chk)
    r1, _ = residual(Ho1, chk)
    print(
        f"even subruns: {ok.sum()} channels with >= {NMIN // 2} muons (median {np.median(ne[ok]):.0f})")
    print(f"odd subruns, {chk.sum()} channels: channel response rms  gain only {np.std(r0[chk]):.2%}  "
          f"-> with c {np.std(r1[chk]):.2%}")
    me, mo = residual(He, chk)[0], residual(Ho0, chk)[0]
    print(
        f"even vs odd reproducibility of the channel response (gain only): {np.std(me[chk] / mo[chk]):.2%}")
    # final factors: all subruns
    Ha = He + Ho0
    na = Ha.sum(1)
    ok = live & (na >= NMIN)
    ra, _ = residual(Ha, ok)
    c = np.where(ok, ra, np.nan)
    bad = live & ~ok
    for i in np.where(bad)[0]:
        oth = [j for j in range(4 * (i // 4), 4 * (i // 4) + 4) if ok[j]]
        c[i] = np.mean(c[oth]) if oth else 1.0
    print(f"final c (all subruns, {ok.sum()} channels, median {np.median(na[ok]):.0f} muons): "
          f"rms {np.nanstd(c[ok]):.2%}, range {np.nanmin(c[ok]):.3f}–{np.nanmax(c[ok]):.3f}; "
          f"set from the module or 1: {[i + 1 for i in np.where(bad)[0]]}")
    c[~live] = np.nan
    ne = na
    np.savez(os.path.join(A.dir, f"muon_cal_{run}.npz"), c=c, ok=ok, entries=ne, entries_odd=no, chk=chk,
             r_odd_gain=r0, r_odd_cal=r1, H_even=He, H_odd_gain=Ho0, H_odd_cal=Ho1, edges=EDGES)
    # the table: use only channels with their own c, in detectors with all four channels alive
    det_alive = np.repeat(live.reshape(-1, 4).all(axis=1), 4)
    use = live & ok & det_alive & np.isfinite(c)
    name = os.path.join(A.dir, f"psmd_gain_{run}.txt")
    with open(name, "w") as f:
        f.write(f"# PSMD gain matching, run {run}:  q = charge / (g * c)\n"
                f"# g   : SiPM gain, fitted area of the dark pulses of the random triggers [ADC*sample / p.e.]\n"
                f"# c   : muon response of the channel / median channel (clean muons)\n"
                f"# use : 1 = usable; 0 = dead, too few muons, or a detector with a dead channel\n"
                f"# pid        g        c  use\n")
        for i in range(520):
            f.write(f"{i + 1:5d} {g[i]:8.3f} {c[i]:8.4f} {int(use[i]):4d}\n")
    print(f"{name}: {use.sum()} channels with use = 1; use = 0 among the live ones: "
          f"{[i + 1 for i in np.where(live & ~use)[0]]}")


if __name__ == "__main__":
    main(A.run)
