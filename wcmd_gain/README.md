# WCMD gain matching

Offline gain of the WCMD tank PMTs (FADC1, PID 1–47) from their own dark pulses, from the FADC
production (`prod_FADC1_<run>.<sub>.root`, branch `dark_q`):

    q = charge / g[pid-1]        in photoelectrons

g is the mean of the single-photoelectron peak of the dark pulses, with no offset. It is the same
reference as the PSMD's first step (`../psmd_gain`): every channel's photoelectron scale from its own
dark pulses, in every run, with no light source.

Analyses read one file per run, `wcmd_gain_<run>.txt`.

## Running

```
source ../../env.sh                          # ROOT + python; any environment with PyROOT, numpy, scipy, matplotlib
PY=$PY ./run_wcmd_gain.sh <run> [prod] [FIRST:LAST]     # prod: /data/amore2/data/PROD by default
```

The production files are `<prod>/<run:06d>/prod_FADC1_<run:06d>.*.root`, or directly in `<prod>`.
`FIRST:LAST` restricts to those subruns, for tests. The driver runs the steps and the check in a
scratch directory under `<prod>/wcmd_gain/`, then keeps

```
<prod>/wcmd_gain/wcmd_gain_<run>.txt          the table
<prod>/wcmd_gain/check_spe_fit_<run>.png      the 47 SPE fits
<prod>/wcmd_gain/check_gain_<run>.png         g of every tank PMT, fit quality
<prod>/wcmd_gain/wcmd_gain_<run>.log          output of every step
```

and deletes the intermediate files. If a step fails it stops and leaves the scratch directory
(`work_<run>_*`). On 384 subruns of run 263 it takes about 3 minutes.

The steps one by one (every script takes `--dir DIR`, those reading the production `--prod DIR`):

```
$PY 1_dark_hist.py <run> [--subruns F:L]     # dark_q spectra per PMT        -> dark_hist_<run>.root
$PY 2_spe_fit.py   <run>                     # SPE fits and the table        -> spe_fit_<run>.npz, wcmd_gain_<run>.txt
$PY check_spe_fit.py <run>                   # -> check_spe_fit_<run>.png, check_gain_<run>.png
```

## The dark pulses — `dark_q` of the production

`produce_fadc` runs the pulse finder of the signal charge (4σ excursions, widened by 3 samples
before and 21 after, merged, summed) over a region with no light and stores the sum as `dark_q`,
0 when nothing is found:

* **normal triggers:** the baseline region `window_normal` in front of the signal. The FADC delay
  puts the front of every record ahead of the trigger; on run 263 (DLY 1100) the first discriminator
  crossings stay at the dark level up to sample 208 and the trigger forms at 297. Every channel is
  used whatever its `tbit`: the discriminator watches the whole record, so a dark pulse above it in
  the front sets `tbit` by itself, and splitting on `tbit` would only move those pulses.
* **pedestal triggers:** the whole record.

## Step 1 — `1_dark_hist.py`

Per PID, `dark_q` (≠ 0) of the pedestal and of the normal triggers, the samples scanned (region
length × records, for the rates), and the baseline bias of a 25-sample window — the mean of the
pedestal tree's pulse-free integral `qp_8_24`, scaled from its 33 samples.

## Step 2 — `2_spe_fit.py`: g

The two spectra of a PMT are added and fitted (`spe.py`) over [−15, 160) ADC·sample:

    n(q) = A0 G(q; m0, σp) + Σ_{n=1..3} An G(q; n μ, sqrt(n σspe² + σp²))

the standard single-photoelectron response without Poisson weights (the An are free): the pedestal
peak — the baseline excursions the 4σ threshold lets through — of width σp, and n photoelectrons of
width σspe each on that width. g = μ. The baseline bias is −0.08 ADC·sample, 0.2 % of g, and the same
in runs 251, 263, 301; it is kept in the table as the check on the no-offset assumption.

## The table — `wcmd_gain_<run>.txt`

```
# pid        g   g_err   s_spe     s_p  rate    bias    chi2  ndf  use
    1   37.471   0.087   15.78    4.67  8.92  -0.086   525.2   80    1
  ...
   49      nan     nan     nan     nan   nan     nan     nan    0    0
```

| column | |
|---|---|
| `pid` | 1..47 tank, 49 door |
| `g`, `g_err` | 1 p.e. mean and its statistical error, ADC·sample per p.e. |
| `s_spe`, `s_p` | SPE width, pedestal width |
| `rate` | 1 p.e. dark rate, kHz |
| `bias` | baseline bias of a 25-sample window, ADC·sample |
| `chi2`, `ndf` | fit quality |
| `use` | 1 when fitted; 0 for the door PMT |

The door PMT (PID 49, FADC2) is not fitted: it triggers on its own (threshold 50 ADC, multiplicity 1),
and the front of its records holds pulses correlated with that trigger (SPE clusters ~500 ns before
it). On the pedestal triggers alone it looks like the tank PMTs; that is left for later.

In an analysis:

```python
import wcmd_gain
t = wcmd_gain.load("/data/amore2/data/PROD/wcmd_gain/wcmd_gain_263.txt")   # pid, g, use by PID-1
q = wcmd_gain.calibrate(charge, pid, t)    # numpy; nan where use == 0
d = wcmd_gain.define(df, t)                # RDataFrame: q[nch] next to charge[nch], 0 where use == 0
```

## Run 263 (subruns 0–383 of the local production)

| | |
|---|---|
| g | median 36.76 ADC·sample/p.e., PMT spread 3.2 % (PMT gains set to ~6×10⁶ by HV) |
| statistics | 0.23 % per PMT |
| fit | χ²/ndf median 5.7; best PID 24 (1.4), worst PID 41 (11.8, the lowest gain, 32.6) |
| stability | 10-subrun groups over 17.9 days: common drift 0.26 % rms, each PMT 1.5 % rms = its fit error |

## FADC delay

From DLY 1100 to DLY 1500 the record moves by 200 samples: the front grows from 208 to 408 samples
and doubles the dark pulses per record. For such runs the production needs `window_normal`
[0, 400] and `integration` [432, 616] (`WCMD/delay_opt/prod_fadc_dly1500.yaml`).
