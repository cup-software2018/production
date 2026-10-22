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

    ./submit.sh <run> [njobs] [from] [to]

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
symlink to it.  One file per subrun per
stream.  They are **not merged** — a `TChain` over them costs nothing to build
and there is then only one copy of the data to keep in step.  Per-subrun logs
are in `out/log/` (`FADC<b>_<sub>.log`, `IADC_<sub>.log`), Slurm's own in
`stage/*/slurm/`, and the build log is `build.log` here.

---

## Settings

Everything that shapes a number lives in the yaml files, each key with a
one-line comment.  Nothing is compiled in except the defaults, and
whatever was actually used is copied into every output file (see *Provenance*),
so a file says how it was made.  A key left out of the yaml keeps its default,
so a config can be a short list of overrides.

---

## FADC — WCMD

### The readout

* NKFADC500: 2 ns per sample, 1008 samples (2.016 µs) per record, every
  channel in every event.  Charges are ADC × samples, baseline subtracted.
* Every channel has a fixed even/odd sample offset, so the baseline is taken
  per parity.
* FADC1 carries the tank PMTs (PID 1–47; 48 is unconnected), FADC2 the door
  PMT (PID 49) on its own trigger path.  `channels` keeps a PID range per
  board, so each board is produced separately.
* Trigger times are ns since the start of the run, counted by each board, so
  FADC1 and FADC2 differ by a fixed offset for the same event.

### Algorithm

One pass over the raw files (several are read as one chain; `--split` rolls
the output over past the given size).  Each event goes to one of the two trees
by its trigger type, and every channel in the `channels` range is processed as
below (a record too short for the windows is skipped and counted).

**Normal trigger** — per channel:

1. **Baseline.**  Over `window_normal`, in front of the signal, take the mean
   of the even and of the odd samples and the rms around them.  Drop the
   samples more than `clip_sigma` × rms *above* their level (a pulse only goes
   up) and recompute; repeat up to `iters` times or until nothing changes.
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
   No pulse found means `charge` = 0 exactly.  The signal is a sparse train of
   narrow pulses, so summing only the pulses drops the noise of the empty
   samples that `qwin` keeps.
5. Copy the DAQ trigger bit → `tbit`; count the channels with it set → `nhit`.
6. **Dark.**  Run the pulse finder of step 4 over `window_normal` and sum what
   it finds → `dark_q`, 0 if nothing.  The front of the record is delayed
   ahead of the trigger and holds no light, so these are dark pulses of every
   channel, whatever its `tbit`.

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
7. **Dark.**  The same pulse finder over `window_pedestal` → `dark_q`.

### Settings (`prod_fadc.yaml`)

| key | |
|---|---|
| `baseline.window_normal` | baseline region of normal triggers, in front of the signal; also the `dark_q` region |
| `baseline.window_pedestal` | baseline region of pedestal triggers; `[0, 0]` for the whole record |
| `baseline.exclude_pedestal` | cut out of it around the peak, relative to the peak |
| `baseline.clip_sigma` / `iters` | one-sided rejection of pulses inside the region |
| `normal.integration` | the signal window |
| `normal.pulse` | `nsigma`, `pre`, `post`: pulse finding |
| `pedestal.windows` | `(head, tail)` charge windows around the peak, at most 8 |
| `pedestal.hmax_range` | where the peak may sit |
| `pedestal.ped_away` / `second_away` | how far `qp_*` and `hmax2` stay from it |
| `pedestal.seed` | random seed for the pulse-free position |
| `channels` | PID range kept per board |

### Output

A quantity both trees have is written under the same name and type.

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
| `nhit` | `Int_t` | | channels with `tbit` set |

Per channel, all `[nch]`, indexed like `pid`:

| branch | type | unit | meaning |
|---|---|---|---|
| `pid` | `Short_t` | | channel id as the DAQ writes it |
| `tbit` | `UShort_t` | | DAQ per-channel discriminator, 0 or 1 |
| `base_even`, `base_odd` | `Float_t` | ADC | baseline of the even / odd samples |
| `rms` | `Float_t` | ADC | its rms after parity subtraction |
| `hmax`, `hmaxx` | `Float_t`, `Short_t` | ADC, sample | largest sample over baseline inside the integration window, and where |
| `charge` | `Float_t` | ADC·sample | **sum over the pulses found** in the integration window |
| `qwin` | `Float_t` | ADC·sample | plain integral of the window |
| `npul` | `Short_t` | | pulses found in the window |
| `dark_q` | `Float_t` | ADC·sample | sum of the pulses found in `window_normal`, 0 if none |

#### Tree `pedestal` — one entry per pedestal trigger

Per event: as `normal`, without `run` (it is in the provenance) and `nhit`.

Per channel, all `[nch]`, indexed like `pid`:

| branch | type | unit | meaning |
|---|---|---|---|
| `pid` | `Short_t` | | channel id as the DAQ writes it |
| `base_even`, `base_odd` | `Float_t` | ADC | baseline of the even / odd samples |
| `rms` | `Float_t` | ADC | its rms after parity subtraction |
| `hmax`, `hmaxx` | `Float_t`, `Short_t` | ADC, sample | largest excursion over baseline **anywhere in `hmax_range`**, and where |
| `hmax2`, `gap2` | `Float_t`, `Short_t` | ADC, sample | next largest, at least `second_away` from the peak, and its signed distance |
| `q_<head>_<tail>` | `Float_t` | ADC·sample | charge over `[peak − head, peak + tail]`, one per window |
| `qp_<head>_<tail>` | `Float_t` | ADC·sample | the same at a pulse-free position |
| `qp_pos` | `Short_t` | sample | where that position was |
| `dark_q` | `Float_t` | ADC·sample | sum of the pulses found in `window_pedestal`, 0 if none |

#### Provenance

The settings used, as `TNamed` objects at the top level of every file:

    run  board  subrun  input  config
    window_normal  window_pedestal  exclude_pedestal  clip_sigma  iters
    integration  pulse  dark
    hmax_range  ped_away  second_away  seed
    channels  maxFileBytes

---

## IADC — PSMD

### The readout

* NKIADC64: 16 ns per sample (inferred from the register layout: one RL unit
  is 128 ns on every Notice board), NDP = 8·RL − 2 samples per record,
  12-bit.  The record length is read from each channel.  Charges are
  ADC × samples over `ped + qbase`.
* No even/odd offset, so one baseline `base` where the FADC has two.
* Zero suppression: a normal trigger holds only the channels over threshold; a
  pedestal trigger holds all 520.
* Module m = (pid − 1) / 4.
* Trigger times are ns since the start of the run.

### Algorithm

One pass over the raw file.  Each event goes to one of the two trees by its
trigger type, and every channel in it is processed as below (a record too
short for the windows is skipped and counted).

**Per channel**, both triggers:

1. **Header pedestal.**  Take the pedestal the ADC writes in the channel
   header → `ped`, and subtract it from every sample.
2. **Level under the charge** (`normal.baseline: fit`, the default).  Run the
   dark-pulse fit below over `window_normal` on normal triggers, over the
   whole record on pedestal triggers; its flat level → `qbase` (0 if the fit
   fails, or with `normal.baseline: header`).
3. **Baseline region.**  Mean and rms of the raw samples over `window_normal`
   on normal triggers, over `window_pedestal` on pedestal triggers → `base`,
   `rms`.  From here on every sample is taken against `ped + qbase`:
   * sum over `window_normal` → `qpre` (the tail of an earlier pulse);
   * sum over `integration` → `charge`; the largest sample and where →
     `hmax`, `hmaxx`;
   * the first rising-edge crossing of `cf` × `hmax`, interpolated → `tlead`
     (−1 if none);
   * the same crossing found walking down from the peak, so an earlier,
     smaller pulse cannot take it, allowed in front of the integration window,
     interpolated → `tcf` in ns (−1 if none or saturated);
   * samples of the record at the 12-bit ceiling (4095) → `nsat`; there `hmax`
     is the ceiling, so no `tcf`.
4. Copy the DAQ trigger bit → `tbit`; count the channels with it set → `nhit`,
   and the modules they belong to → `nmod`.
5. **Dark pulses** from the fit of step 2 → `fit_*`, `npul`, `pul_*` (on
   pedestal triggers only if `pedestal.fit.enabled`).

**Dark-pulse fit** — a flat level plus n Moyal pulses of one common width
(`DarkFit`, the same algorithm as `../view_waveform_iadc.py`):

1. **Level.**  Start from the lower 20 % of the samples and average those
   within 2.5 `sigma_e` of it, a few times; dark pulses only push samples up.
2. **Seeds.**  Every local maximum at least `seed_height` above that level, at
   least 3 samples from a higher one and standing out by 0.8 × `seed_height`;
   also a falling edge at the start of the record or a rising edge at its end
   (a pulse just outside it).
3. **Fit** the level, the common width (inside `sigma_range`) and a height and
   position per pulse, by Levenberg-Marquardt.
4. **Add pulses.**  Where the waveform is furthest above the fit, if by more
   than `add_residual`, add a pulse and refit; keep it only if chi2 drops by
   at least `add_dchi2`.  Repeat up to `max_pulses`.
5. **Keep** the pulses higher than `keep_height`; area = 4.1327 × height ×
   width.

### Settings (`prod_iadc.yaml`)

| key | |
|---|---|
| `baseline.window_normal` | baseline region of normal triggers, in front of the signal |
| `baseline.window_pedestal` | baseline region of pedestal triggers; `[0, 0]` for the whole record |
| `normal.integration` | the signal window |
| `normal.cf` | constant fraction for `tlead` and `tcf` |
| `normal.baseline` | `fit` (dark-pulse fit level, the default) or `header` (header pedestal alone) |
| `pedestal.fit.*` | the dark-pulse fit: `enabled`, `sigma_e`, `seed_height`, `add_residual`, `add_dchi2`, `max_pulses` (≤ 12), `sigma_range`, `mu_min`, `keep_height` |

### Output

A quantity both trees have is written under the same name and type, and the
names follow the FADC.

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
| `nhit` | `Int_t` | | channels with `tbit` set |
| `nmod` | `Int_t` | | modules among the channels |
| `npulse` | `Int_t` | | dark pulses in the event, the length of the `pul_*` arrays |

Per channel, all `[nch]`, indexed like `pid`:

| branch | type | unit | meaning |
|---|---|---|---|
| `pid` | `Short_t` | | channel id as the DAQ writes it |
| `tbit` | `UShort_t` | | DAQ per-channel discriminator, 0 or 1 |
| `ped` | `UShort_t` | ADC | header pedestal |
| `qbase` | `Float_t` | ADC over ped | the level everything below is taken against |
| `base`, `rms` | `Float_t` | ADC | mean and rms of `window_normal` |
| `qpre` | `Float_t` | ADC·sample | sum over `window_normal` |
| `hmax`, `hmaxx` | `Float_t`, `Short_t` | ADC, sample | largest sample inside the integration window, and where |
| `tlead` | `Float_t` | sample | cf crossing on the rising edge, −1 if none |
| `tcf` | `Float_t` | ns | cf crossing walking down from the peak, −1 if none or saturated |
| `nsat` | `Short_t` | | samples at the ADC ceiling; > 0 = saturated |
| `charge` | `Float_t` | ADC·sample | integral of the integration window |
| `fit_base` | `Float_t` | ADC over ped | fitted flat level |
| `fit_sigma` | `Float_t` | sample | common Moyal width |
| `fit_chi2` | `Float_t` | | chi2/ndf with white noise `sigma_e` (the noise is correlated, so below 1 is normal) |
| `npul` | `Short_t` | | dark pulses kept, −1 if the fit failed |

Per dark pulse, all `[npulse]`:

| branch | type | unit | meaning |
|---|---|---|---|
| `pul_ch` | `Short_t` | | index into the `[nch]` arrays |
| `pul_t` | `Float_t` | sample | Moyal peak position (may lie outside the record) |
| `pul_h` | `Float_t` | ADC | peak height |
| `pul_area` | `Float_t` | ADC·sample | 4.1327 · h · sigma |

`fit_*`, `npul`, `npulse` and `pul_*` are written only when the fit runs
(`normal.baseline: fit` on `normal`, `pedestal.fit.enabled` on `pedestal`).

#### Tree `pedestal` — one entry per pedestal trigger

Per event: as `normal`, without `run` (it is in the provenance), `nhit` and
`nmod`.  Per channel and per dark pulse: the same branches as `normal`, with
`base`, `rms` over `window_pedestal` and the dark-pulse fit over the whole
record.

#### Provenance

The settings used, as `TNamed` objects at the top level of every file:

    run  subrun  input  config
    window_normal  window_pedestal  integration  cf  baseline  sample_ns
    adc_ceiling  pedestal_fit
