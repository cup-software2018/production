# PSMD gain matching

Offline calibration of the PSMD (IADC) channel charge, from the IADC production
(`prod_IADC_<run>.<sub>.root`).  Two steps, applied together:

    q = charge / (g[pid-1] * c[pid-1])

* **g** — the SiPM gain of the channel, from the dark pulses on the random
  (pedestal) triggers.  `charge / g` is in photoelectrons.
* **c** — the channel's light collection relative to the median channel, from
  clean muons.  It takes out what g leaves: scintillator–SiPM coupling,
  attenuation, photon detection efficiency.

Analyses read one file per run, `psmd_gain_<run>.txt`.

## Running

```
source ../../env.sh                          # ROOT + python; any environment with PyROOT, numpy, scipy, matplotlib
PY=$PY ./run_psmd_gain.sh <run> [prod]       # prod: /data/amore2/data/PROD by default
```

The production files are `<prod>/<run:06d>/prod_IADC_<run:06d>.*.root`, or
directly in `<prod>` when there is no run directory.  The driver runs steps
1–3 and the checks in a scratch directory under `<prod>/psmd_gain/`, then keeps

```
<prod>/psmd_gain/psmd_gain_<run>.txt            the table
<prod>/psmd_gain/check_dark_gain_<run>.png      g of every channel
<prod>/psmd_gain/check_dark_gain_1_<run>.png    dark-pulse spectra of PID 1 with the fitted p.e.
<prod>/psmd_gain/check_gain_match_<run>.png     clean-muon response per channel, before and after g
<prod>/psmd_gain/check_muon_cal_<run>.png       c of every channel
<prod>/psmd_gain/check_odd_match_<run>.png      matching per channel on the odd subruns, c from the even ones
<prod>/psmd_gain/psmd_gain_<run>.log            output of every step
```

and deletes the intermediate files.  If a step fails it stops and leaves the
scratch directory (`work_<run>_*`) for inspection.  A run of 811 subruns
takes about 9 minutes.

The steps can also be run one by one; every script takes `--dir DIR` (where its
inputs are and its outputs go, the current directory by default), and those
reading the production take `--prod DIR`:

```
$PY 1_ped_hist.py  <run>                     # pedestal-trigger histograms  -> ped_hist_<run>.root, .npz
$PY 2_dark_gain.py <run>                     # g                            -> gain_<run>.npz
$PY 3_muon_cal.py  <run>                     # c and the table              -> muon_cal_<run>.npz, psmd_gain_<run>.txt
$PY check_dark_gain.py  <run> [--pid N]      # needs ped_hist, gain npz
$PY check_muon_cal.py   <run>                # needs muon_cal npz
$PY check_gain_match.py hist <run>           # -> check_gain_match_<run>.npz ...
$PY check_gain_match.py plot <run>           # ... and the plot
```

## Step 1 — `1_ped_hist.py`

From the pedestal tree, two histograms against PID: the area `pul_area` and
height `pul_h` of the dark pulses the production fitted on the record.

## Step 2 — `2_dark_gain.py`: g

For every live PID, g is the slope through the 1, 2, 3 p.e. peaks of the
fitted dark-pulse area; each peak is a Gaussian fitted within ±0.35 g of its
expected place, twice.  The peaks are sharp, so this fixes each channel's gain
precisely (±0.15 %).  The p.e. unit is the same in every run, so cuts in p.e.
mean the same in every run.

## Step 3 — `3_muon_cal.py`: c

A **clean muon** is a detector (4 channels in the event) where both slabs have
more than 100 p.e. after g, and no channel carries more than 55 % of the
detector total.  For each PID, c = the median of its charge / g over its clean
muons, divided by the median over channels.

* The selection always uses `charge / g`, never `charge / (g c)`: selecting with
  the corrected charge feeds back (a low c inflates the channel, lets
  lower-charge detectors in, lowers the median further) and runs away.
* c is first taken from the even subruns and applied to the odd ones; then the
  final c comes from all subruns.
* A channel with fewer than 200 clean muons gets the mean c of the other
  channels of its detector (1 if none) and `use = 0`.

## The table — `psmd_gain_<run>.txt`

```
# pid        g        c  use
    1  295.567   0.9023    1
    2  296.159   0.9576    1
  ...
  205      nan      nan    0
```

| column | |
|---|---|
| `pid` | channel, 1..520 |
| `g` | SiPM gain, ADC·sample per p.e. of the fitted dark-pulse area; `nan` = dead (no dark pulses) |
| `c` | muon response / median channel; `nan` = dead |
| `use` | 1 when g and c exist, c comes from the channel's own muons, and all four channels of its detector are alive; 0 otherwise |

`use = 0` on a live channel means its c is not trustworthy: too few muons, or a
dead channel next to it in the detector (a slab read by one SiPM has half the
light, which biases the clean-muon selection).  Leave those channels out when
combining detectors.

In an analysis:

```python
import psmd_gain
t = psmd_gain.load("/data/amore2/data/PROD/psmd_gain/psmd_gain_<run>.txt")   # dict of arrays by PID-1: pid, g, c, use
q = psmd_gain.calibrate(charge, pid, t)    # numpy; nan where use == 0
d = psmd_gain.define(df, t)                # RDataFrame: q[nch] next to charge[nch], 0 where use == 0
```

or in C++ / a ROOT macro, with the header-only `psmd_gain.hh` (this directory on
the include path, e.g. `ROOT_INCLUDE_PATH`, or `-I`):

```cpp
#include "psmd_gain.hh"
psmd_gain::Table t = psmd_gain::Load("/data/amore2/data/PROD/psmd_gain/psmd_gain_<run>.txt");  // by PID-1: g, c, use
double q = t.Calibrate(charge, pid);       // NaN where use == 0
auto d = psmd_gain::Define(df, t);         // RDataFrame: q[nch] next to charge[nch], 0 where use == 0
```

Both also give the production file glob of a run, `prod_files(prod, run)` /
`psmd_gain::ProdFiles(prod, run)`.
