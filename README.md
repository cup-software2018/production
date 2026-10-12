# Production — WCMD (FADC) and PSMD (IADC)

One directory produces both muon-veto detectors of a run.  Each has its own
producer and settings; the scripts, the build and the output location are
shared.

| detector | raw files | producer | settings |
|---|---|---|---|
| WCMD | `FADC<b>_<run>.root.<sub>` (b = 1, 2) | `produce_fadc.cc` | `prod_fadc.yaml` |
| PSMD | `IADC_<run>.root.<sub>` | `produce_iadc.cc` | `prod_iadc.yaml` |

Each producer reads a subrun once and fills two trees with the same names and
the same naming of branches and provenance:

| tree | trigger | |
|---|---|---|
| `normal` | physics triggers | the charge in the signal window |
| `pedestal` | `trgType == 1`, the periodic random trigger | the electronics with no light: baseline, noise, single photoelectrons |

The random trigger fires regardless of what the detector is doing.  It is the
only unbiased handle on the electronics, and on the PSMD, whose IADC drops
every channel that did not fire (zero suppression), the only record of every
channel.

---

## Running

A whole run, on the cluster, one job per subrun:

    ./submit.sh [run] [njobs] [from] [to]

    ./submit.sh <run> 8                      every subrun on disk, 8 at a time
    ./submit.sh <run> 8 0 600                subruns 0..600 only
    CFG_FADC=scan.yaml ./submit.sh <run> 8   another FADC config (CFG_IADC likewise)

One subrun of every detector, FADC1, FADC2 and then the IADC, wherever it is
called (each Slurm task is one of these):

    ./prod_all.sh <run> <subrun>

One subrun of one detector:

    ./prod.sh <fadc|iadc> <run> <subrun> [boards]

    OUT=/tmp/t ./prod.sh fadc <run> <subrun>     FADC1 and FADC2
    OUT=/tmp/t ./prod.sh fadc <run> <subrun> 2   FADC2 only
    OUT=/tmp/t ./prod.sh iadc <run> <subrun>     IADC

`prod.sh` decides every path: the raw file under `RAW`
(`/data/amore2/data/RAW`), the output and its log under `OUT` (`out/` here),
the settings `CFG` (`prod_<det>.yaml`) and the program `BIN`
(`build/produce_<det>`), each overridable from the environment.  A stream with
no raw file for the subrun is skipped with a note; the exit status is 1 if any
stream fails and 3 if none has input.  `prod_all.sh` calls it for `fadc` and
`iadc` with each one's `prod_<det>.yaml`, skips a detector with no input, and
fails if any producer fails.  Run by hand both need a build first, and an `OUT`
this machine can write: `out/` links to `/data/cupsoft`, which it cannot.

### The cluster

`submit.sh` runs one Slurm array task per subrun, every detector in it,
`njobs` at a time, on the cluster behind qserver (`cupsoft@192.168.0.160`,
nodes qnode1-2), submitted over passwordless ssh from this machine.  A slow
subrun then holds up one slot instead of a whole block, and a job that dies
costs one subrun to redo.  The script waits and prints a summary at the end;
Ctrl-C stops the waiting, not the jobs (`squeue -u cupsoft` / `scancel` on
qserver).

`submit.sh` builds both producers here with CMake before anything is
submitted: the nodes have the yaml-cpp runtime but not its headers.  The
programs, `prod.sh`, `prod_all.sh`, the two configs and the list of subruns
(those with raw data for any detector) then go to a fresh
`stage/<run>_<time>/` directory under `RDIR` on qserver, so a later
submission cannot change the code under jobs that are still running.  Task i
runs that `prod_all.sh` on line i of `subruns.txt`.

The nodes see qserver's `/home` and `/data`, not this machine's `/scratch`, and
cupsoft is uid 1000 here but 1005 there, so this side cannot write to
`/data/cupsoft`.  That is why everything is copied over through qserver.  The
programs are linked against ROOT and `CUPDATAMODEL/V1.1` under
`/home/cupsoft/Products`, the same path and the same files on both sides, and
their RPATH points there.

The producers of a task run one after another, so each task asks for 1 CPU
and 2 GB; a producer stays under 1 GB.  Arrays are split at 1000 tasks (the
cluster's `MaxArraySize` is 1001), and each piece starts when the previous one
finishes.  This side's NFS view of `/data` can lag the nodes by a few seconds,
so the end-of-run summary is computed on qserver.

`QHOST` and `RDIR` override the host and the remote directory:

    RDIR=/data/cupsoft/WCMD/test ./submit.sh <run> 4 0 3

### Building and running by hand

    source env.sh
    cmake -S . -B build && cmake --build build

    R=/data/amore2/data/RAW/<run:06d>
    build/produce_fadc -c prod_fadc.yaml -o t.root $R/FADC1_<run:06d>.root.<sub:05d>
    build/produce_fadc -c prod_fadc.yaml -o t.root $R/FADC2_<run:06d>.root.*      chained
    build/produce_fadc -c prod_fadc.yaml -o t.root -n 20000 <input>               first 20000 events
    build/produce_iadc -c prod_iadc.yaml -o t.root $R/IADC_<run:06d>.root.<sub:05d>

The producers take raw files and an output file and nothing else; where the
raw data lives and where the output goes is decided in `prod.sh`.  Run, board
and subrun are read from the DAQ's file names (`FADC<b>_<run>.root.<sub>`,
`IADC_<run>.root.<sub>`).  Several FADC inputs are chained into one output, in
the order given, and must be the same run and board.  The `input` provenance
is the file path, or `first .. last (N files)` for a chain.

`-h` lists the options.  Without `-c` the built-in defaults are used, which
are not those of the yaml files.  Both programs exit non-zero on an error.

The two producers were ROOT macros (ACLiC) until 2026-10-08; that version is
kept in `old/`, and the compiled programs reproduce its output value for value.

### Output location

`$RDIR/out/prod_FADC<board>_<run:06d>.<subrun:05d>.root` and
`$RDIR/out/prod_IADC_<run:06d>.<subrun:05d>.root`, with `RDIR` defaulting to
`/data/cupsoft/WCMD/production` (the name predates the PSMD); `out/` here is a
symlink to it, so `wcmd_io` finds the files there.  One file per subrun per
stream.  They are **not merged** — a `TChain` over them costs nothing to build
and there is then only one copy of the data to keep in step.  Per-subrun logs
are in `out/log/` (`FADC<b>_<sub>.log`, `IADC_<sub>.log`), Slurm's own in
`stage/*/slurm/`, and the build log is `build.log` here.

---

## Settings

Everything that shapes a number lives in the yaml files, with the reason for
each value written beside it.  Nothing is compiled in except the defaults, and
whatever was actually used is copied into every output file (see *Provenance*),
so a file says how it was made.  A key left out of the yaml keeps its default,
so a config can be a short list of overrides.

---

## FADC — WCMD

### Algorithm

One pass over the raw file.  Each event goes to one of the two trees by its
trigger type, and every channel in the `channels` PID range is processed as
below (a record too short for the windows is skipped and counted).

**Normal trigger** — per channel:

1. **Baseline.**  Over `window_normal`, take the mean of the even samples and
   of the odd samples separately, and the rms around them.  Drop the samples
   more than `clip_sigma` × rms *above* their level (a pulse only goes up) and
   recompute; repeat up to `iters` times or until nothing changes.
   → `base_even`, `base_odd`, `rms`
2. **Subtract** the even level from even samples and the odd level from odd
   samples.  Everything below uses this baseline-subtracted waveform.
3. **Window.**  In `integration`, sum every sample → `qwin`; find the largest
   sample and where it is → `hmax`, `hmaxx`.
4. **Pulses.**  Still in `integration`:
   1. a pulse is a run of samples above `nsigma` × `rms`; count them → `npul`;
   2. widen each pulse by `pre` samples before and `post` after;
   3. merge the widened pulses that touch or overlap;
   4. sum the samples they cover → `charge`.
   No pulse found means `charge` = 0 exactly.
5. Copy the DAQ trigger bit → `tbit`; count the channels with it set → `nhit`.

**Pedestal trigger** — per channel (no light, so the whole record is used):

1. **Baseline** as in step 1 above, over `window_pedestal` (the whole record
   by default).
2. **Peak.**  Find the largest sample in `hmax_range`.
3. **Baseline again**, this time leaving out `exclude_pedestal` around that
   peak, so the pulse does not lift it; find the peak again with the new
   baseline → `base_even`, `base_odd`, `rms`, `hmax`, `hmaxx`.
4. **Second peak.**  The largest sample at least `second_away` from the peak
   → `hmax2`, and its signed distance → `gap2`.
5. **Pulse-free position.**  Pick a random position in `hmax_range` at least
   `ped_away` from the peak (up to 8 tries) → `qp_pos`.
6. **Window charges.**  For each `(head, tail)` in `windows`, sum
   `[peak − head, peak + tail]` → `q_<head>_<tail>`, and the same window at the
   pulse-free position → `qp_<head>_<tail>`, the noise reference.

Several input files are read as one chain; with `--split` the output rolls
over to a new file past the given size.

### Settings (`prod_fadc.yaml`)

| key | |
|---|---|
| `baseline.window_normal` | baseline region of physics triggers, in front of the signal |
| `baseline.window_pedestal` | baseline region of pedestal triggers; `[0, 0]` for the whole record |
| `baseline.exclude_pedestal` | cut out of it around the peak, relative to the peak |
| `baseline.clip_sigma` / `iters` | one-sided rejection of pulses inside the region |
| `normal.integration` | the signal window |
| `normal.pulse` | `nsigma`, `pre`, `post`: pulse finding inside that window |
| `pedestal.windows` | `(head, tail)` charge windows around the peak, at most 8 |
| `pedestal.hmax_range` | where the peak may sit |
| `pedestal.ped_away` / `second_away` | how far `qp_*` and `hmax2` stay from it |
| `pedestal.seed` | random seed for the pulse-free position |
| `channels` | PID range kept per board |

#### Baselines

Both are parity-split: every channel carries a fixed even/odd sample offset,
and the two parities are measured and subtracted separately.  The clipping is
one-sided — a pulse only pushes samples up — and it is applied to the
parity-subtracted deviation, not to the raw sample against the averaged level,
which would bite far harder into one parity than the other.

The two triggers get different regions because they must.  A physics trigger
has light in the record, so only the part in front of the signal is usable; a
pedestal trigger has none, so the whole record is, apart from the pulse it just
found.  Past the signal the AC-coupling undershoot sits below the true level,
so a physics trigger cannot use the back of its record.

#### Why a pulse sum

The signal is a sparse train of narrow pulses, so most of a fixed window holds
only noise.  Summing only the pulses found keeps the signal and drops that
noise: a channel that saw nothing has `charge` exactly 0, while the plain
window integral `qwin`, kept beside it, wanders with the noise.

#### Channels

`channels` keeps a PID range per board, so an unconnected channel is dropped
here rather than written out and filtered later.  The door PMT reads out on
FADC2 on **its own trigger path**, so that board is produced separately and its
numbers are not mixed with the tank's.

### Output format

Each tree carries its own branches in full; nothing is stored once and shared
between them.  What they do share is naming: a quantity that means the same
thing in both is written under the same branch name with the same type, so a
reader that knows one tree can read the other.  `pedestal` is `normal` without
`run`, `nhit` and the branches that need a trigger, plus the branches that only
a signal-free record can give.

#### Tree `normal` — one entry per normal trigger

Per event:

| branch | type | unit | meaning |
|---|---|---|---|
| `run` | `Int_t` | | run number |
| `nch` | `Int_t` | | channels in this entry; the length of every array below |
| `subrun` | `Int_t` | | which subrun this event came from |
| `evtNum` | `UInt_t` | | event number from the DAQ |
| `trgNum` | `UInt_t` | | trigger number from the DAQ |
| `trgTime` | `ULong64_t` | ns | since the start of the run |
| `nhit` | `Int_t` | | how many channels have `tbit` set |

Per channel, all of length `nch` and indexed the same way as `pid`:

| branch | type | unit | meaning |
|---|---|---|---|
| `pid` | `Short_t` | | channel id as the DAQ writes it |
| `base_even` | `Float_t` | ADC | baseline of the even samples |
| `base_odd` | `Float_t` | ADC | baseline of the odd samples |
| `rms` | `Float_t` | ADC | its rms after parity subtraction |
| `hmax` | `Float_t` | ADC | largest sample over baseline, **inside the integration window** |
| `hmaxx` | `Short_t` | sample | where that sample is |
| `tbit` | `UShort_t` | | the DAQ's own per-channel discriminator, 0 or 1 |
| `charge` | `Float_t` | ADC×sample | **sum over the pulses found** in the integration window |
| `qwin` | `Float_t` | ADC×sample | plain integral of the whole window |
| `npul` | `Short_t` | | pulses the finder separated inside the window |

`nhit` is nothing but the count of channels the DAQ itself called hit, so it
can be compared with the trigger's own multiplicity requirement without
depending on anything this production decided.  Any offline hit definition —
a charge cut, a gain-corrected charge cut — belongs in the analysis, where it
can be changed without reproducing a run.

There is no `qtot`.  A total charge is a sum over a per-channel calibration,
and a sum of raw ADC over channels of different gain is not a quantity worth
storing.

#### Tree `pedestal` — one entry per pedestal trigger

Per event:

| branch | type | unit | meaning |
|---|---|---|---|
| `nch` | `Int_t` | | channels in this entry; the length of every array below |
| `subrun` | `Int_t` | | which subrun this event came from |
| `evtNum` | `UInt_t` | | event number from the DAQ |
| `trgNum` | `UInt_t` | | trigger number from the DAQ |
| `trgTime` | `ULong64_t` | ns | since the start of the run |

No `run`: the number is in the file's own `TNamed` provenance below, and in the
`normal` tree of the same file.  No `nhit`: there is no trigger bit to count.

Per channel, all of length `nch` and indexed the same way as `pid`:

| branch | type | unit | meaning |
|---|---|---|---|
| `pid` | `Short_t` | | channel id as the DAQ writes it |
| `base_even` | `Float_t` | ADC | baseline of the even samples |
| `base_odd` | `Float_t` | ADC | baseline of the odd samples |
| `rms` | `Float_t` | ADC | its rms after parity subtraction |
| `hmax` | `Float_t` | ADC | largest excursion over baseline, **anywhere in `hmax_range`** |
| `hmaxx` | `Short_t` | sample | where it is |
| `hmax2` | `Float_t` | ADC | next largest, at least `second_away` from the peak |
| `gap2` | `Short_t` | sample | signed distance to it |
| `q_<head>_<tail>` | `Float_t` | ADC×sample | charge over `[peak-head, peak+tail]`, one per window |
| `qp_<head>_<tail>` | `Float_t` | ADC×sample | the same integrator at a pulse-free position |
| `qp_pos` | `Short_t` | sample | where that position was |

`hmax` and `hmaxx` are the one pair whose shared name hides a difference: here
the peak is searched over the whole of `hmax_range`, because a pedestal record
has no signal window to look in, while in `normal` it is searched inside the
integration window.

The live time behind anything measured from this tree is the width of
`hmax_range` per entry per channel, not the record.

#### Provenance

Every file carries the settings it was made with as `TNamed` objects at the top
level, readable without the trees:

    run  board  subrun  input  config
    window_normal  window_pedestal  exclude_pedestal  clip_sigma  iters
    integration  pulse
    hmax_range  ped_away  second_away  seed
    channels  maxFileBytes

#### Units and conventions

A sample is 2 ns (NKFADC500 at 500 MS/s).  Charges are ADC × samples with the
baseline subtracted.  Trigger times are ns since the start of the run, counted
separately by each board, so FADC1 and FADC2 differ by a fixed offset for the
same physical event.

### Reading it

`../WCMD/wcmd_io.py` chains the per-subrun files and flattens the per-channel
arrays into numpy columns:

    import wcmd_io
    d = wcmd_io.load(["pid", "hmax", "q_4_8", "qp_4_8"], tree="pedestal", run=<run>)
    d["hmax"]          # one element per (record, channel)
    d["rec"]           # which record each element came from

The flattening runs inside ROOT through `TTree::Draw`, not an entry loop in
python; for millions of channels that is the difference between seconds and an
afternoon.  Ask for `pmt` instead of `pid` to get the door PMT numbered next to
the tank's, so the two boards can sit side by side.

The `normal` tree of a whole run is usually too large for that, so read it one
subrun file at a time and histogram as you go; `../WCMD/charge_dist.py` is a
worked example.

---

## IADC — PSMD

### The readout

* The IADC goes through the same DAQ writer as the FADC, so a file holds
  `FChannelData` waveforms: 12-bit, positive pulses, NDP = 8·RL − 2 samples per
  channel.  **16 ns per sample** is inferred from the register layout (one RL
  unit is 128 ns on every Notice board), not from a datasheet.  The record
  length is read from each channel, not compiled in.
* Zero suppression is on: a physics event holds only the channels over
  threshold, so the `normal` tree carries those and the `pedestal` tree carries
  every channel.
* **Module m = (pid − 1) / 4.**  The trigger is a coincidence over these
  groups.  Inside a module, `k = (pid − 1) % 4`: k 0,1 read one scintillator
  layer, k 2,3 the other; k 0 and 2 are at the same end, k 1 and 3 at the
  other.

### Algorithm

One pass over the raw file.  Each event goes to one of the two trees by its
trigger type; every channel in it is processed the same way, apart from the
baseline region in step 3 and the fit in step 5 (a record too short for the
windows is skipped and counted).

Per channel:

1. **Header pedestal.**  Take the pedestal the ADC writes in the channel
   header → `ped`, and subtract it from every sample.
2. **Level under the charge** (`normal.baseline: fit`, the default).  Fit the
   samples of `window_normal` with the dark-pulse fit below; its flat level is
   the true baseline, which the dark pulses had lifted the header pedestal
   above → `qbase` (0 if the fit fails, or with `normal.baseline: header`).
   This is done on both trees, so the pedestal tree stays the noise reference
   of the normal one.
3. **Baseline region.**  Mean and rms of the raw samples over `window_normal`
   on normal triggers, over `window_pedestal` (the whole record by default) on
   pedestal triggers → `base`, `rms`.  From here on every sample is taken
   against `ped + qbase`:
   * sum over `window_normal` → `qpre` (the tail of an earlier pulse shows up
     here);
   * sum over `integration` → `charge`; the largest sample and where →
     `hmax`, `hmaxx`;
   * the first point on the rising edge where the waveform crosses `cf` ×
     `hmax`, interpolated between samples → `tlead` (−1 if none).
4. Copy the DAQ trigger bit → `tbit`; count the channels with it set → `nhit`,
   and the modules they belong to → `nmod`.
5. **Dark pulses** → `fit_*`, `npul`, `pul_*`.  On pedestal triggers fit the
   whole record (if `pedestal.fit.enabled`); on normal triggers keep the
   `window_normal` fit of step 2.

**The dark-pulse fit** — a flat level plus a few pulses of one common shape:

1. **Level.**  Start from the lower 20% of the samples and average those
   within 2.5 `sigma_e` of it, a few times; dark pulses only push samples up.
2. **Seeds.**  Every local maximum at least `seed_height` above that level,
   at least 3 samples from a higher one and standing out by 0.8 ×
   `seed_height`, is a pulse to start from; so is a falling edge at the start
   of the record or a rising edge at its end (a pulse just outside it).
3. **Fit** the level, one common width (kept inside `sigma_range`) and a
   height and position per pulse, by Levenberg-Marquardt.
4. **Add pulses.**  Where the waveform is furthest above the fit, if by more
   than `add_residual`, add a pulse there and refit; keep it only if chi2
   drops by at least `add_dchi2`.  Repeat up to `max_pulses`.
5. **Keep** the pulses higher than `keep_height`; each one's area is
   4.1327 × height × width.

### Settings (`prod_iadc.yaml`)

| key | |
|---|---|
| `baseline.window_normal` | region in front of the signal on physics triggers |
| `baseline.window_pedestal` | baseline region of pedestal triggers; `[0, 0]` for the whole record |
| `normal.integration` | the signal window |
| `normal.cf` | constant fraction for `tlead` |
| `normal.baseline` | `fit` (dark-pulse fit level, the default) or `header` (header pedestal alone) |
| `pedestal.fit.*` | the dark-pulse fit: `enabled`, `sigma_e`, `seed_height`, `add_residual`, `add_dchi2`, `max_pulses` (≤ 12), `sigma_range`, `mu_min`, `keep_height` |

### Branches

Per event: `nch subrun evtNum trgNum trgTime(ns)`, plus `run nhit nmod` on
`normal` (`nhit` = channels with `tbit`, as in the WCMD; `nmod` = modules
among the channels).  Per channel (`[nch]`, indexed like `pid`), the same in
both trees:

| branch | unit | meaning |
|---|---|---|
| `pid` | | channel id as the DAQ writes it |
| `tbit` | | DAQ trigger bit |
| `ped` | ADC | header pedestal |
| `base`, `rms` | ADC | mean and rms of the baseline region: `window_normal` on `normal`, `window_pedestal` on `pedestal` |
| `qbase` | ADC over ped | the level charge, hmax, tlead and qpre are taken against |
| `qpre` | ADC·sample | Σ (w − ped − qbase) over `window_normal`: pile-up from an earlier pulse |
| `hmax`, `hmaxx` | ADC, sample | largest sample over ped + qbase inside the integration window and where |
| `tlead` | sample | constant-fraction crossing on the rising edge, −1 if none |
| `charge` | ADC·sample | as above |

There is no even/odd split: the IADC shows none, so `base` stands where the
WCMD has `base_even` / `base_odd`.

#### Dark-pulse fit

On the pedestal tree every channel is fitted over the whole record; on the
normal tree over `window_normal` (the fit that gives `qbase`).  The model is a
flat level plus n Moyal pulses of a common width (`DarkFit` in
`produce_iadc.cc`, settings under `pedestal.fit`, the same algorithm as the
python one in `../view_waveform_iadc.py`; steps in *Algorithm* above).  The point is
the gain: each dark pulse gives its own area, so the random triggers of a
single subrun already hold enough pulses per channel for a gain.  The same
branches on the normal tree hold the dark pulses found in front of the signal.

| branch | unit | meaning |
|---|---|---|
| `fit_base[nch]` | ADC over ped | fitted flat level (the pulse-free baseline) |
| `fit_sigma[nch]` | sample | common Moyal width, bounded to `sigma_range` |
| `fit_chi2[nch]` | | chi2/ndf with white noise `sigma_e` (the noise is correlated, so below 1 is normal) |
| `npul[nch]` | | pulses kept on the channel, −1 if the fit failed |
| `npulse` | | pulses in the event, length of the `pul_*` arrays |
| `pul_ch[npulse]` | | index into the `[nch]` arrays |
| `pul_t[npulse]` | sample | Moyal peak position (may lie outside the record for edge pulses) |
| `pul_h[npulse]` | ADC | peak height |
| `pul_area[npulse]` | ADC·sample | 4.1327 · h · sigma, the whole Moyal area |

#### Provenance

    run  subrun  input  config
    window_normal  window_pedestal  integration  cf  baseline  sample_ns
    pedestal_fit
