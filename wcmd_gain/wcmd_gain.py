"""Read the WCMD gain table (wcmd_gain_<run>.txt, from 2_spe_fit.py) and apply it.

    q = charge / g[pid-1]                 in photoelectrons

    import wcmd_gain
    t = wcmd_gain.load("wcmd_gain_263.txt")    # dict of arrays by PID-1: pid, g, use
    q = wcmd_gain.calibrate(charge, pid, t)    # numpy; nan where use == 0
    d = wcmd_gain.define(df, t)                # RDataFrame: q[nch] next to charge[nch], 0 where use == 0
"""
import os
import numpy as np


def prod_files(prod, run, board=1):
    """Glob of the FADC production files of a run: <prod>/<run:06d>/ (as PROD), or <prod>/ directly."""
    sub = os.path.join(prod, f"{run:06d}")
    return os.path.join(sub if os.path.isdir(sub) else prod, f"prod_FADC{board}_{run:06d}.*.root")


def load(path):
    t = np.genfromtxt(path, comments="#", names=None)
    order = np.argsort(t[:, 0])
    t = t[order]
    return {"pid": t[:, 0].astype(int), "g": t[:, 1], "use": t[:, -1].astype(bool)}


def scale(t):
    """g by PID-1 (index pid - 1 up to the largest PID in the table), nan where not usable."""
    s = np.full(int(t["pid"].max()), np.nan)
    s[t["pid"] - 1] = np.where(t["use"], t["g"], np.nan)
    return s


def calibrate(charge, pid, t):
    return np.asarray(charge) / scale(t)[np.asarray(pid) - 1]


def define(df, t, col="q", name="wcmd_gain_scale"):
    import ROOT
    s = scale(t)
    if not hasattr(ROOT, name):
        vals = ", ".join("0" if not np.isfinite(x) else f"{x:.4f}" for x in s)
        ROOT.gInterpreter.Declare(f"const float {name}[{len(s)}] = {{{vals}}};")
    return df.Define(col, f"ROOT::RVec<float> o(nch); for (int i = 0; i < nch; i++) "
                          f"{{ int k = pid[i] - 1; float g = k >= 0 && k < {len(s)} ? {name}[k] : 0.f; "
                          f"o[i] = g > 0 ? charge[i] / g : 0.f; }} return o;")
