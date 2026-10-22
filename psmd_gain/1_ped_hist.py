"""PSMD gain matching, step 1: per-channel histograms of the pedestal-trigger (random) records.

    python 1_ped_hist.py <run> [--prod DIR] [--dir DIR]   -> ped_hist_<run>.root, ped_hist_<run>.npz

Input: the pedestal tree of every prod_IADC_<run>.<sub>.root in <prod>/<run> (or <prod>).
  area     : pid vs pul_area (fitted dark-pulse area, ADC*sample)
  amp      : pid vs pul_h    (fitted dark-pulse height, ADC)
The .npz holds the same histograms as arrays [pid-1, bin] with their bin edges.
"""
import argparse
import os
import numpy as np
from psmd_gain import prod_files
import ROOT

p = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("run", type=int)
p.add_argument("--prod", default="/data/amore2/data/PROD",
               help="production directory")
p.add_argument("--dir", default=".", help="where the output goes")
a = p.parse_args()

ROOT.EnableImplicitMT(16)
ROOT.gInterpreter.Declare("""
ROOT::RVec<short> PulPid(const ROOT::RVec<short>& pid, const ROOT::RVec<short>& ch)
{ ROOT::RVec<short> o(ch.size()); for (size_t i = 0; i < ch.size(); i++) o[i] = pid[ch[i]]; return o; }
""")
df = ROOT.RDataFrame("pedestal", prod_files(a.prod, a.run))
d = df.Define("ppid", "PulPid(pid, pul_ch)")
E = {"area": (600, 0, 1500), "amp": (400, 0, 200)}
hs = {"area": d.Histo2D(("h_area", "dark-pulse area;PID;area [ADC*sample]", 520, 0.5, 520.5, *E["area"]), "ppid", "pul_area"),
      "amp": d.Histo2D(("h_amp", "dark-pulse height;PID;height [ADC]", 520, 0.5, 520.5, *E["amp"]), "ppid", "pul_h")}
n = d.Count()

f = ROOT.TFile(os.path.join(a.dir, f"ped_hist_{a.run}.root"), "recreate")
out = {}
for k, h in hs.items():
    h.GetValue().Write()
    out[k] = np.array([[h.GetBinContent(i, j) for j in range(
        1, h.GetNbinsY() + 1)] for i in range(1, 521)])
    nb, lo, hi = E[k]
    out[f"{k}_edges"] = np.linspace(lo, hi, nb + 1)
f.Close()
np.savez(os.path.join(a.dir, f"ped_hist_{a.run}.npz"), **out)
print(f"run {a.run}: {n.GetValue():,} pedestal events -> ped_hist_{a.run}.root, .npz")
