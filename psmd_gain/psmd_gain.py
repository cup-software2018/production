"""Read the PSMD gain-matching table (psmd_gain_<run>.txt, from 3_muon_cal.py) and apply it.

    q = charge / (g[pid-1] * c[pid-1])

    import psmd_gain
    t = psmd_gain.load("psmd_gain_<run>.txt")    # dict of arrays by PID-1: pid, g, c, use
    q = psmd_gain.calibrate(charge, pid, t)    # numpy; nan where use == 0
    d = psmd_gain.define(df, t)                # RDataFrame: q[nch] next to charge[nch], 0 where use == 0
"""
import os
import numpy as np


def prod_files(prod, run):
    """Glob of the IADC production files of a run: <prod>/<run:06d>/ (as PROD), or <prod>/ directly."""
    sub = os.path.join(prod, f"{run:06d}")
    return os.path.join(sub if os.path.isdir(sub) else prod, f"prod_IADC_{run:06d}.*.root")


def load(path):
    pid, g, c, use = np.loadtxt(path, comments="#", unpack=True)
    order = np.argsort(pid)
    return {"pid": pid[order].astype(int), "g": g[order], "c": c[order], "use": use[order].astype(bool)}


def scale(t):
    """g * c by PID-1, nan where the channel is not usable."""
    return np.where(t["use"], t["g"] * t["c"], np.nan)


def calibrate(charge, pid, t):
    return np.asarray(charge) / scale(t)[np.asarray(pid) - 1]


def define(df, t, col="q", name="psmd_gain_scale"):
    import ROOT
    s = scale(t)
    if not hasattr(ROOT, name):
        vals = ", ".join("0" if not np.isfinite(x) else f"{x:.5f}" for x in s)
        ROOT.gInterpreter.Declare(
            f"const float {name}[{len(s)}] = {{{vals}}};")
    return df.Define(col, f"ROOT::RVec<float> o(nch); for (int i = 0; i < nch; i++) "
                          f"{{ float s = {name}[pid[i] - 1]; o[i] = s > 0 ? charge[i] / s : 0.f; }} return o;")
