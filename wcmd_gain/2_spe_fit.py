"""WCMD gain matching, step 2: the SPE fit of every tank PMT and the table.

    python 2_spe_fit.py <run> [--dir DIR]      (needs dark_hist_<run>.root)
        -> wcmd_gain_<run>.txt     the table analyses read
        -> spe_fit_<run>.npz       the fit parameters, for the check

The pedestal and normal-trigger spectra of a PMT are added and fitted with spe.py; the gain g is the
1 p.e. mean (no offset), charge / g is in photoelectrons.  The door PMT (PID 49, FADC2) is not
fitted: it triggers on its own and the front of its records holds correlated pulses.
"""
import argparse
import os
import numpy as np
import ROOT
import spe

p = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("run", type=int)
p.add_argument("--dir", default=".",
               help="where the input is and the output goes")
a = p.parse_args()
RUN = a.run
NTANK = 47

f = ROOT.TFile.Open(os.path.join(a.dir, f"dark_hist_{RUN}.root"))
h0 = f.Get("h_q0")
ny = h0.GetNbinsY()
e = np.array([h0.GetYaxis().GetBinLowEdge(j) for j in range(1, ny + 2)])
x = 0.5 * (e[1:] + e[:-1])
hb = f.Get("h_bias")

P = np.full((NTANK, len(spe.NAMES)), np.nan)
E = np.full((NTANK, len(spe.NAMES)), np.nan)
chi2 = np.full(NTANK, np.nan)
ndf = np.zeros(NTANK, int)
live = np.zeros(NTANK)
bias = np.array([hb.GetBinContent(k) for k in range(1, NTANK + 1)])
for k in range(NTANK):
    pid = k + 1
    y = sum(np.array([f.Get(f"h_q{r}").GetBinContent(pid, j)
            for j in range(1, ny + 1)]) for r in range(2))
    live[k] = sum(f.Get(f"h_n{r}").GetBinContent(pid) for r in range(2)) * 2e-9
    r = spe.fit(x, y)
    if r is None:
        print(f"PID {pid}: fit failed")
        continue
    P[k], E[k], chi2[k], ndf[k] = r
np.savez(os.path.join(a.dir, f"spe_fit_{RUN}.npz"), p=P, e=E, chi2=chi2, ndf=ndf, live=live, bias=bias,
         names=np.array(spe.NAMES), files=f.Get("files").GetTitle())

i = {n: j for j, n in enumerate(spe.NAMES)}
g, ge, s1, sp = P[:, i["mu"]], E[:, i["mu"]], P[:, i["s1"]], P[:, i["sp"]]
rate = P[:, i["a1"]] / live / 1e3
use = np.isfinite(g)
name = os.path.join(a.dir, f"wcmd_gain_{RUN}.txt")
with open(name, "w") as out:
    out.write(f"# WCMD gain matching, run {RUN}:  q = charge / g  [p.e.]\n"
              f"# g      : 1 p.e. mean of the dark pulses (no offset) [ADC*sample / p.e.]\n"
              f"# s_spe  : SPE width, s_p: pedestal width, rate: 1 p.e. dark rate [kHz]\n"
              f"# bias   : baseline bias of a 25-sample window [ADC*sample], the check on no offset\n"
              f"# use    : 1 = fitted; PID 49 (door, FADC2) is not fitted\n"
              f"# input  : {f.Get('files').GetTitle()}\n"
              f"# pid        g   g_err   s_spe     s_p  rate    bias    chi2  ndf  use\n")
    for k in range(NTANK):
        out.write(f"{k + 1:5d} {g[k]:8.3f} {ge[k]:7.3f} {s1[k]:7.2f} {sp[k]:7.2f} {rate[k]:5.2f} {bias[k]:+7.3f} "
                  f"{chi2[k]:7.1f} {ndf[k]:4d} {int(use[k]):4d}\n")
    out.write(f"{49:5d} {np.nan:8.3f} {np.nan:7.3f} {np.nan:7.2f} {np.nan:7.2f} {np.nan:5.2f} {np.nan:7.3f} "
              f"{np.nan:7.1f} {0:4d} {0:4d}\n")
ok = use
print(f"run {RUN}: {ok.sum()} of {NTANK} tank PMTs fitted; g median {np.median(g[ok]):.2f}, spread "
      f"{np.std(g[ok]) / np.median(g[ok]):.1%}, stat err median {np.median(ge[ok] / g[ok]):.2%}; chi2/ndf median "
      f"{np.median(chi2[ok] / ndf[ok]):.2f} (best PID {np.nanargmin(chi2 / ndf) + 1}, worst PID {np.nanargmax(chi2 / ndf) + 1}); "
      f"baseline bias median {np.median(bias):+.3f} ({np.median(bias) / np.median(g[ok]):+.2%} of g)")
print(f"-> {name}")
