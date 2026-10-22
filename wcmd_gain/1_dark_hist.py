"""WCMD gain matching, step 1: the dark-pulse charge spectra of every tank PMT.

    python 1_dark_hist.py <run> [--prod DIR] [--dir DIR] [--subruns FIRST:LAST]   -> dark_hist_<run>.root

Input: dark_q of every prod_FADC1_<run>.<sub>.root in <prod>/<run> (or <prod>): the sum of the 4σ
pulses in a region with no light -- the front (window_normal) of the normal triggers, the whole
record of the pedestal ones -- 0 when there is none.
  h_q0 / h_q1 : PID x dark_q (dark_q != 0) of the pedestal / normal triggers
  h_n0 / h_n1 : PID: samples scanned (region length x records), for the rates
  h_bias      : PID: baseline bias of a 25-sample window -- the mean of the pedestal tree's pulse-free
                integral qp_8_24 scaled from its 33 samples -- the check on the no-offset gain
"""
import argparse
import glob
import os
import re
import ROOT
from wcmd_gain import prod_files

p = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("run", type=int)
p.add_argument("--prod", default="/data/amore2/data/PROD",
               help="production directory")
p.add_argument("--dir", default=".", help="where the output goes")
p.add_argument("--subruns", default=None,
               help="FIRST:LAST, only these subruns (for tests)")
a = p.parse_args()
ROOT.gROOT.SetBatch(True)

files = sorted(glob.glob(prod_files(a.prod, a.run)))
files = [f for f in files if re.search(r"\.\d{5}\.root$", f)]
if a.subruns:
    lo, hi = map(int, a.subruns.split(":"))
    files = [f for f in files if lo <= int(f[-10:-5]) <= hi]
if not files:
    raise SystemExit(f"no files {prod_files(a.prod, a.run)}")
f0 = ROOT.TFile.Open(files[0])
nlo, nhi = map(int, re.findall(r"\d+", f0.Get("window_normal").GetTitle()))
wped = f0.Get("window_pedestal").GetTitle()
nped = 1008 if wped == "whole record" else (
    lambda x, y: y - x)(*map(int, re.findall(r"\d+", wped)))
f0.Close()

out = ROOT.TFile(os.path.join(a.dir, f"dark_hist_{a.run}.root"), "recreate")
for r, (tree, length) in enumerate([("pedestal", nped), ("normal", nhi - nlo)]):
    ch = ROOT.TChain(tree)
    for f in files:
        ch.Add(f)
    hq = ROOT.TH2D(
        f"h_q{r}", f"{tree};PID;dark_q [ADC*sample]", 49, 0.5, 49.5, 300, -50, 550)
    hn = ROOT.TH1D(f"h_n{r}", f"{tree};PID;samples scanned", 49, 0.5, 49.5)
    ch.Draw(f"dark_q:pid>>h_q{r}", "dark_q != 0", "goff")
    ch.Draw(f"pid>>h_n{r}", f"{length}", "goff")
    out.cd()
    hq.Write()
    hn.Write()
    if tree == "pedestal":
        hb = ROOT.TProfile(
            "h_bias", "baseline bias of a 25-sample window;PID;ADC*sample", 49, 0.5, 49.5)
        ch.Draw("qp_8_24 * 25. / 33.:pid>>h_bias", "", "goff prof")
        out.cd()
        hb.Write()
    print(f"{tree}: {ch.GetEntries()} records, region {length} samples")
ROOT.TNamed("files", f"{len(files)} files, {files[0]} .. {files[-1]}").Write()
ROOT.TNamed("window_normal", f"[{nlo},{nhi})").Write()
out.Close()
print(f"run {a.run}: {len(files)} files -> dark_hist_{a.run}.root")
